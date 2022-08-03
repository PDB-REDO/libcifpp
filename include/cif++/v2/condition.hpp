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
// some more templates to be able to do querying

namespace detail
{
	template <typename Category>
	struct condition_impl
	{
		using category_type = Category;
		using row_type = row_handle<category_type>;

		virtual ~condition_impl() {}

		virtual void prepare(const Category &c) {}
		virtual bool test(const Category &c, const row_type &r) const = 0;
		virtual void str(std::ostream &os) const = 0;
	};

	template <typename Category>
	struct all_condition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		virtual bool test(const Category &c, const row_type &r) const { return true; }
		virtual void str(std::ostream &os) const { os << "*"; }
	};

	template <typename>
	struct or_condition_impl;
	template <typename>
	struct and_condition_impl;
	template <typename>
	struct not_condition_impl;
} // namespace detail

template <typename Category>
class condition_t
{
  public:
	using category_type = Category;
	using condition_impl = detail::condition_impl<category_type>;
	using row_type = row_handle<category_type>;

	condition_t()
		: m_impl(nullptr)
	{
	}
	condition_t(condition_impl *impl)
		: m_impl(impl)
	{
	}

	condition_t(const condition_t &) = delete;

	condition_t(condition_t &&rhs) noexcept
		: m_impl(nullptr)
	{
		std::swap(m_impl, rhs.m_impl);
	}

	condition_t &operator=(const condition_t &) = delete;

	condition_t &operator=(condition_t &&rhs) noexcept
	{
		std::swap(m_impl, rhs.m_impl);
		return *this;
	}

	~condition_t()
	{
		delete m_impl;
		m_impl = nullptr;
	}

	void prepare(const category_type &c)
	{
		if (m_impl)
			m_impl->prepare(c);
		m_prepared = true;
	}

	bool operator()(const category_type &c, const row_type &r) const
	{
		assert(this->m_impl != nullptr);
		assert(this->m_prepared);
		return m_impl ? m_impl->test(c, r) : false;
	}

	bool empty() const { return m_impl == nullptr; }

	template<typename C> friend condition_t operator||(condition_t<C> &&a, condition_t<C> &&b);
	template<typename C> friend condition_t operator&&(condition_t<C> &&a, condition_t<C> &&b);

	template <typename>
	friend struct detail::or_condition_impl;
	template <typename>
	friend struct detail::and_condition_impl;
	template <typename>
	friend struct detail::not_condition_impl;

	void swap(condition_t &rhs)
	{
		std::swap(m_impl, rhs.m_impl);
		std::swap(m_prepared, rhs.m_prepared);
	}

	template<typename C> friend std::ostream &operator<<(std::ostream &os, const condition_t<C> &cond);

  private:
	condition_impl *m_impl;
	bool m_prepared = false;
};

template <typename Category>
inline std::ostream &operator<<(std::ostream &os, const condition_t<Category> &cond)
{
	if (cond.m_impl)
		cond.m_impl->str(os);
	return os;
}

namespace detail
{

	template <typename Category>
	struct keyIsemptycondition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		keyIsemptycondition_impl(const std::string &item_tag)
			: m_item_tag(item_tag)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const row_type &r) const
		{
			return r[m_item_ix].empty();
		}

		virtual void str(std::ostream &os) const
		{
			os << m_item_tag << " IS NULL";
		}

		std::string m_item_tag;
		size_t m_item_ix = 0;
	};

	template <typename Category>
	struct keyComparecondition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		template <typename COMP>
		keyComparecondition_impl(const std::string &item_tag, COMP &&comp, const std::string &s)
			: m_item_tag(item_tag)
			, m_compare(std::move(comp))
			, mStr(s)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const row_type &r) const
		{
			return m_compare(c, r, m_icase);
		}

		virtual void str(std::ostream &os) const
		{
			os << m_item_tag << (m_icase ? "^ " : " ") << mStr;
		}

		std::string m_item_tag;
		size_t m_item_ix = 0;
		bool m_icase = false;
		std::function<bool(const Category &, const row_type &, bool)> m_compare;
		std::string mStr;
	};

	template <typename Category>
	struct keyMatchescondition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		keyMatchescondition_impl(const std::string &item_tag, const std::regex &rx)
			: m_item_tag(item_tag)
			, m_item_ix(0)
			, mRx(rx)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const row_type &r) const
		{
			std::string_view txt = r[m_item_ix].text();
			return std::regex_match(txt.begin(), txt.end(), mRx);
		}

		virtual void str(std::ostream &os) const
		{
			os << m_item_tag << " =~ expression";
		}

		std::string m_item_tag;
		size_t m_item_ix;
		std::regex mRx;
	};

	template <typename Category, typename T>
	struct AnyIscondition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		typedef T valueType;

		AnyIscondition_impl(const valueType &value)
			: mValue(value)
		{
		}

		virtual bool test(const Category &c, const row_type &r) const;
		virtual void str(std::ostream &os) const
		{
			os << "<any> == " << mValue;
		}

		valueType mValue;
	};

	template <typename Category>
	struct AnyMatchescondition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		AnyMatchescondition_impl(const std::regex &rx)
			: mRx(rx)
		{
		}

		virtual bool test(const Category &c, const row_type &r) const;
		virtual void str(std::ostream &os) const
		{
			os << "<any> =~ expression";
		}

		std::regex mRx;
	};

	template <typename Category>
	struct and_condition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		and_condition_impl(condition_t<Category> &&a, condition_t<Category> &&b)
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

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		virtual bool test(const Category &c, const row_type &r) const
		{
			return mA->test(c, r) and mB->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << '(';
			mA->str(os);
			os << ") AND (";
			mB->str(os);
			os << ')';
		}

		base_type *mA;
		base_type *mB;
	};

	template <typename Category>
	struct or_condition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		or_condition_impl(condition_t<Category> &&a, condition_t<Category> &&b)
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

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		virtual bool test(const Category &c, const row_type &r) const
		{
			return mA->test(c, r) or mB->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << '(';
			mA->str(os);
			os << ") OR (";
			mB->str(os);
			os << ')';
		}

		base_type *mA;
		base_type *mB;
	};

	template <typename Category>
	struct not_condition_impl : public condition_impl<Category>
	{
		using base_type = condition_impl<Category>;
		using row_type = base_type::row_type;

		not_condition_impl(condition_t<Category> &&a)
			: mA(nullptr)
		{
			std::swap(mA, a.m_impl);
		}

		~not_condition_impl()
		{
			delete mA;
		}

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
		}

		virtual bool test(const Category &c, const row_type &r) const
		{
			return not mA->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << "NOT (";
			mA->str(os);
			os << ')';
		}

		base_type *mA;
	};

} // namespace detail

template <typename Category>
inline condition_t<Category> operator&&(condition_t<Category> &&a, condition_t<Category> &&b)
{
	if (a.m_impl and b.m_impl)
		return condition_t<Category>(new detail::and_condition_impl<Category>(std::move(a), std::move(b)));
	if (a.m_impl)
		return condition_t<Category>(std::move(a));
	return condition_t<Category>(std::move(b));
}

template <typename Category>
inline condition_t<Category> operator||(condition_t<Category> &&a, condition_t<Category> &&b)
{
	if (a.m_impl and b.m_impl)
		return condition_t<Category>(new detail::or_condition_impl<Category>(std::move(a), std::move(b)));
	if (a.m_impl)
		return condition_t<Category>(std::move(a));
	return condition_t<Category>(std::move(b));
}

struct empty
{
};

/// \brief A helper to make it possible to have conditions like ("id"_key == cif::null)

inline constexpr empty null = empty();

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

template <typename Category, typename T>
condition_t<Category> operator==(const key &key, const T &v)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	std::ostringstream s;
	s << " == " << v;

	return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
		key.m_item_tag, [tag = key.m_item_tag, v](const Category &c, const row_type &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) == 0; },
		s.str()));
}

template <typename Category>
inline condition_t<Category> operator==(const key &key, const char *value)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	if (value != nullptr and *value != 0)
	{
		std::ostringstream s;
		s << " == " << value;

		return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
			key.m_item_tag, [tag = key.m_item_tag, value](const Category &c, const row_type &r, bool icase)
			{ return r[tag].compare(value, icase) == 0; },
			s.str()));
	}
	else
		return condition_t(new detail::keyIsemptycondition_impl<Category>(key.m_item_tag));
}

// inline condition_t operator==(const key& key, const detail::ItemReference& v)
// {
// 	if (v.empty())
// 		return condition_t(new detail::keyIsemptycondition_impl(key.m_item_tag));
// 	else
// 		return condition_t(new detail::keyComparecondition_impl(key.m_item_tag, [tag = key.m_item_tag, v](const Category& c, const row_type& r, bool icase)
// 			{ return r[tag].template compare<(v, icase) == 0; }));
// }

template <typename Category, typename T>
condition_t<Category> operator!=(const key &key, const T &v)
{
	return condition_t<Category>(new detail::not_condition_impl<Category>(operator==(key, v)));
}

template <typename Category>
inline condition_t<Category> operator!=(const key &key, const char *v)
{
	std::string value(v ? v : "");
	return condition_t<Category>(new detail::not_condition_impl<Category>(operator==(key, value)));
}

template <typename Category, typename T>
condition_t<Category> operator>(const key &key, const T &v)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	std::ostringstream s;
	s << " > " << v;

	return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
		key.m_item_tag, [tag = key.m_item_tag, v](const Category &c, const row_type &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) > 0; },
		s.str()));
}

template <typename Category, typename T>
condition_t<Category> operator>=(const key &key, const T &v)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	std::ostringstream s;
	s << " >= " << v;

	return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
		key.m_item_tag, [tag = key.m_item_tag, v](const Category &c, const row_type &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) >= 0; },
		s.str()));
}

template <typename Category, typename T>
condition_t<Category> operator<(const key &key, const T &v)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	std::ostringstream s;
	s << " < " << v;

	return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
		key.m_item_tag, [tag = key.m_item_tag, v](const Category &c, const row_type &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) < 0; },
		s.str()));
}

template <typename Category, typename T>
condition_t<Category> operator<=(const key &key, const T &v)
{
	using category_type = Category;
	using row_type = row_handle<category_type>;

	std::ostringstream s;
	s << " <= " << v;

	return condition_t<Category>(new detail::keyComparecondition_impl<Category>(
		key.m_item_tag, [tag = key.m_item_tag, v](const Category &c, const row_type &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) <= 0; },
		s.str()));
}

template <typename Category>
inline condition_t<Category> operator==(const key &key, const std::regex &rx)
{
	return condition_t<Category>(new detail::keyMatchescondition_impl<Category>(key.m_item_tag, rx));
}

template <typename Category>
inline condition_t<Category> operator==(const key &key, const empty &)
{
	return condition_t<Category>(new detail::keyIsemptycondition_impl<Category>(key.m_item_tag));
}

template <typename Category>
struct any_t
{
};

template <typename Category, typename T>
condition_t<Category> operator==(const any_t<Category> &, const T &v)
{
	return condition_t<Category>(new detail::AnyIscondition_impl<Category, T>(v));
}

template <typename Category>
condition_t<Category> operator==(any_t<Category> &, const std::regex &rx)
{
	return condition_t<Category>(new detail::AnyMatchescondition_impl<Category>(rx));
}

template <typename Category>
inline condition_t<Category> all()
{
	return condition_t<Category>(new detail::all_condition_impl<Category>());
}

// inline condition_t<Category> Not(condition_t<Category> &&cond)
// {
// 	return condition_t<Category>(new detail::not_condition_impl<Category>(std::move(cond)));
// }

namespace literals
{
	inline key operator""_key(const char *text, size_t length)
	{
		return key(std::string(text, length));
	}

	inline constexpr empty null = empty();

} // namespace literals

} // namespace cif::v2