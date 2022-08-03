/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <cassert>
#include <functional>
#include <iostream>
#include <regex>

#include <cif++/v2/row.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------
// let's make life easier

iset get_category_fields(const category &cat);
uint16_t get_column_ix(const category &cat, std::string_view col);
bool is_column_type_uchar(const category &cat, std::string_view col);

// --------------------------------------------------------------------
// some more templates to be able to do querying

namespace detail
{
	struct condition_impl
	{
		virtual ~condition_impl() {}

		virtual void prepare(const category &c) {}
		virtual bool test(row_handle r) const = 0;
		virtual void str(std::ostream &os) const = 0;
	};

	struct all_condition_impl : public condition_impl
	{
		bool test(row_handle r) const override { return true; }
		void str(std::ostream &os) const override { os << "*"; }
	};

	struct or_condition_impl;
	struct and_condition_impl;
	struct not_condition_impl;
} // namespace detail

class condition
{
  public:
	using condition_impl = detail::condition_impl;

	condition()
		: m_impl(nullptr)
	{
	}
	condition(condition_impl *impl)
		: m_impl(impl)
	{
	}

	condition(const condition &) = delete;

	condition(condition &&rhs) noexcept
		: m_impl(nullptr)
	{
		std::swap(m_impl, rhs.m_impl);
	}

	condition &operator=(const condition &) = delete;

	condition &operator=(condition &&rhs) noexcept
	{
		std::swap(m_impl, rhs.m_impl);
		return *this;
	}

	~condition()
	{
		delete m_impl;
		m_impl = nullptr;
	}

	void prepare(const category &c)
	{
		if (m_impl)
			m_impl->prepare(c);
		m_prepared = true;
	}

	bool operator()(row_handle r) const
	{
		assert(this->m_impl != nullptr);
		assert(this->m_prepared);
		return m_impl ? m_impl->test(r) : false;
	}

	bool empty() const { return m_impl == nullptr; }

	friend condition operator||(condition &&a, condition &&b);
	friend condition operator&&(condition &&a, condition &&b);

	friend struct detail::or_condition_impl;
	friend struct detail::and_condition_impl;
	friend struct detail::not_condition_impl;

	void swap(condition &rhs)
	{
		std::swap(m_impl, rhs.m_impl);
		std::swap(m_prepared, rhs.m_prepared);
	}

	friend std::ostream &operator<<(std::ostream &os, const condition &cond)
	{
		if (cond.m_impl)
			cond.m_impl->str(os);
		return os;
	}

  private:
	condition_impl *m_impl;
	bool m_prepared = false;
};

namespace detail
{
	struct key_is_empty_condition_impl : public condition_impl
	{
		key_is_empty_condition_impl(const std::string &item_tag)
			: m_item_tag(item_tag)
		{
		}

		void prepare(const category &c) override
		{
			m_item_ix = get_column_ix(c, m_item_tag);
		}

		bool test(row_handle r) const override
		{
			return r[m_item_ix].empty();
		}

		void str(std::ostream &os) const override
		{
			os << m_item_tag << " IS NULL";
		}

		std::string m_item_tag;
		size_t m_item_ix = 0;
	};

	struct key_compare_condition_impl : public condition_impl
	{
		template <typename COMP>
		key_compare_condition_impl(const std::string &item_tag, COMP &&comp, const std::string &s)
			: m_item_tag(item_tag)
			, m_compare(std::move(comp))
			, m_str(s)
		{
		}

		void prepare(const category &c) override
		{
			m_item_ix = get_column_ix(c, m_item_tag);
			m_icase = is_column_type_uchar(c, m_item_tag);
		}

		bool test(row_handle r) const override
		{
			return m_compare(r, m_icase);
		}

		void str(std::ostream &os) const override
		{
			os << m_item_tag << (m_icase ? "^ " : " ") << m_str;
		}

		std::string m_item_tag;
		size_t m_item_ix = 0;
		bool m_icase = false;
		std::function<bool(row_handle, bool)> m_compare;
		std::string m_str;
	};

	struct key_matches_condition_impl : public condition_impl
	{
		key_matches_condition_impl(const std::string &item_tag, const std::regex &rx)
			: m_item_tag(item_tag)
			, m_item_ix(0)
			, mRx(rx)
		{
		}

		void prepare(const category &c) override
		{
			m_item_ix = get_column_ix(c, m_item_tag);
		}

		bool test(row_handle r) const override
		{
			std::string_view txt = r[m_item_ix].text();
			return std::regex_match(txt.begin(), txt.end(), mRx);
		}

		void str(std::ostream &os) const override
		{
			os << m_item_tag << " =~ expression";
		}

		std::string m_item_tag;
		size_t m_item_ix;
		std::regex mRx;
	};

	template <typename T>
	struct any_is_condition_impl : public condition_impl
	{
		typedef T valueType;

		any_is_condition_impl(const valueType &value)
			: mValue(value)
		{
		}

		bool test(row_handle r) const override
		{
			auto &c = r.cat();

			bool result = false;
			for (auto &f : get_category_fields(c))
			{
				try
				{
					if (r[f].compare(mValue) == 0)
					{
						result = true;
						break;
					}
				}
				catch (...)
				{
				}
			}

			return result;
		}

		void str(std::ostream &os) const override
		{
			os << "<any> == " << mValue;
		}

		valueType mValue;
	};

	struct any_matches_condition_impl : public condition_impl
	{
		any_matches_condition_impl(const std::regex &rx)
			: mRx(rx)
		{
		}

		bool test(row_handle r) const override
		{
			auto &c = r.cat();

			bool result = false;
			for (auto &f : get_category_fields(c))
			{
				try
				{
					std::string_view txt = r[f].text();
					if (std::regex_match(txt.begin(), txt.end(), mRx))
					{
						result = true;
						break;
					}
				}
				catch (...)
				{
				}
			}

			return result;
		}

		void str(std::ostream &os) const override
		{
			os << "<any> =~ expression";
		}

		std::regex mRx;
	};

	struct and_condition_impl : public condition_impl
	{
		and_condition_impl(condition &&a, condition &&b)
			: mA(nullptr)
			, mB(nullptr)
		{
			std::swap(mA, a.m_impl);
			std::swap(mB, b.m_impl);
		}

		~and_condition_impl()
		{
			delete mA;
			delete mB;
		}

		void prepare(const category &c) override
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		bool test(row_handle r) const override
		{
			return mA->test(r) and mB->test(r);
		}

		void str(std::ostream &os) const override
		{
			os << '(';
			mA->str(os);
			os << ") AND (";
			mB->str(os);
			os << ')';
		}

		condition_impl *mA;
		condition_impl *mB;
	};

	struct or_condition_impl : public condition_impl
	{
		or_condition_impl(condition &&a, condition &&b)
			: mA(nullptr)
			, mB(nullptr)
		{
			std::swap(mA, a.m_impl);
			std::swap(mB, b.m_impl);
		}

		~or_condition_impl()
		{
			delete mA;
			delete mB;
		}

		void prepare(const category &c) override
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		bool test(row_handle r) const override
		{
			return mA->test(r) or mB->test(r);
		}

		void str(std::ostream &os) const override
		{
			os << '(';
			mA->str(os);
			os << ") OR (";
			mB->str(os);
			os << ')';
		}

		condition_impl *mA;
		condition_impl *mB;
	};

	struct not_condition_impl : public condition_impl
	{
		not_condition_impl(condition &&a)
			: mA(nullptr)
		{
			std::swap(mA, a.m_impl);
		}

		~not_condition_impl()
		{
			delete mA;
		}

		void prepare(const category &c) override
		{
			mA->prepare(c);
		}

		bool test(row_handle r) const override
		{
			return not mA->test(r);
		}

		void str(std::ostream &os) const override
		{
			os << "NOT (";
			mA->str(os);
			os << ')';
		}

		condition_impl *mA;
	};

} // namespace detail

inline condition operator&&(condition &&a, condition &&b)
{
	if (a.m_impl and b.m_impl)
		return condition(new detail::and_condition_impl(std::move(a), std::move(b)));
	if (a.m_impl)
		return condition(std::move(a));
	return condition(std::move(b));
}

inline condition operator||(condition &&a, condition &&b)
{
	if (a.m_impl and b.m_impl)
		return condition(new detail::or_condition_impl(std::move(a), std::move(b)));
	if (a.m_impl)
		return condition(std::move(a));
	return condition(std::move(b));
}

struct empty_type
{
};

/// \brief A helper to make it possible to have conditions like ("id"_key == cif::null)

inline constexpr empty_type null = empty_type();

struct key
{
	explicit key(const std::string &itemTag)
		: m_item_tag(itemTag)
	{
	}

	explicit key(const char *itemTag)
		: m_item_tag(itemTag)
	{
	}

	key(const key &) = delete;
	key &operator=(const key &) = delete;

	std::string m_item_tag;
};

template <typename T>
condition operator==(const key &key, const T &v)
{
	std::ostringstream s;
	s << " == " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_tag, [tag = key.m_item_tag, v](row_handle r, bool icase)
		{ return r[tag].template compare<T>(v, icase) == 0; },
		s.str()));
}

inline condition operator==(const key &key, const char *value)
{
	if (value != nullptr and *value != 0)
	{
		std::ostringstream s;
		s << " == " << value;

		return condition(new detail::key_compare_condition_impl(
			key.m_item_tag, [tag = key.m_item_tag, value](row_handle r, bool icase)
			{ return r[tag].compare(value, icase) == 0; },
			s.str()));
	}
	else
		return condition(new detail::key_is_empty_condition_impl(key.m_item_tag));
}

// inline condition_t operator==(const key& key, const detail::ItemReference& v)
// {
// 	if (v.empty())
// 		return condition_t(new detail::key_is_empty_condition_impl(key.m_item_tag));
// 	else
// 		return condition_t(new detail::key_compare_condition_impl(key.m_item_tag, [tag = key.m_item_tag, v](const category& c, const row& r, bool icase)
// 			{ return r[tag].template compare<(v, icase) == 0; }));
// }

template <typename T>
condition operator!=(const key &key, const T &v)
{
	return condition(new detail::not_condition_impl(operator==(key, v)));
}

inline condition operator!=(const key &key, const char *v)
{
	std::string value(v ? v : "");
	return condition(new detail::not_condition_impl(operator==(key, value)));
}

template <typename T>
condition operator>(const key &key, const T &v)
{
	std::ostringstream s;
	s << " > " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_tag, [tag = key.m_item_tag, v](const category &c, row_handle r, bool icase)
		{ return r[tag].template compare<T>(v, icase) > 0; },
		s.str()));
}

template <typename T>
condition operator>=(const key &key, const T &v)
{
	std::ostringstream s;
	s << " >= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_tag, [tag = key.m_item_tag, v](const category &c, row_handle r, bool icase)
		{ return r[tag].template compare<T>(v, icase) >= 0; },
		s.str()));
}

template <typename T>
condition operator<(const key &key, const T &v)
{
	std::ostringstream s;
	s << " < " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_tag, [tag = key.m_item_tag, v](const category &c, row_handle r, bool icase)
		{ return r[tag].template compare<T>(v, icase) < 0; },
		s.str()));
}

template <typename T>
condition operator<=(const key &key, const T &v)
{
	std::ostringstream s;
	s << " <= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_tag, [tag = key.m_item_tag, v](const category &c, row_handle r, bool icase)
		{ return r[tag].template compare<T>(v, icase) <= 0; },
		s.str()));
}

inline condition operator==(const key &key, const std::regex &rx)
{
	return condition(new detail::key_matches_condition_impl(key.m_item_tag, rx));
}

inline condition operator==(const key &key, const empty_type &)
{
	return condition(new detail::key_is_empty_condition_impl(key.m_item_tag));
}

struct any_type
{
};

inline constexpr any_type any = any_type{};

template <typename T>
condition operator==(const any_type &, const T &v)
{
	return condition(new detail::any_is_condition_impl<T>(v));
}

inline condition operator==(const any_type &, const std::regex &rx)
{
	return condition(new detail::any_matches_condition_impl(rx));
}

inline condition all()
{
	return condition(new detail::all_condition_impl());
}

// inline condition_t<category> Not(condition_t<category> &&cond)
// {
// 	return condition_t<category>(new detail::not_condition_impl(std::move(cond)));
// }

namespace literals
{
	inline key operator""_key(const char *text, size_t length)
	{
		return key(std::string(text, length));
	}
} // namespace literals

} // namespace cif::v2