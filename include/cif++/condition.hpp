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

#include "cif++/row.hpp"

#include <cassert>
#include <concepts>
#include <functional>
#include <iostream>
#include <regex>
#include <utility>

/** \file condition.hpp
 * This file contains code to create conditions: object encapsulating a
 * query you can use to find rows in a @ref cif::category
 *
 * Conditions are created as standard C++ expressions. That means
 * you can use the standard comparison operators to compare item
 * contents with a value and boolean operators to chain everything
 * together.
 *
 * To create a query that simply compares one item with one value:
 *
 * @code {.cpp}
 * cif::condition c = cif::key("id") == 1;
 * @endcode
 * 
 * That will find rows where the ID item contains the number 1. If
 * using cif::key is a bit too much typing, you can also write:
 * 
 * @code{.cpp}
 * using namespace cif::literals;
 * 
 * cif::condition c2 = "id"_key == 1;
 * @endcode
 * 
 * Now if you want both ID = 1 and ID = 2 in the result:
 * 
 * @code{.cpp}
 * auto c3 = "id"_key == 1 or "id"_key == 2;
 * @endcode
 * 
 * There are some special values you can use. To find rows with item that
 * do not have a value:
 * 
 * @code{.cpp}
 * auto c4 = "type"_key == cif::null;
 * @endcode 
 * 
 * Of if it should not be NULL:
 * 
 * @code{.cpp}
 * auto c5 = "type"_key != cif::null;
 * @endcode 
 * 
 * There's even a way to find all records:
 * 
 * @code{.cpp}
 * auto c6 = cif::all;
 * @endcode
 * 
 * And when you want to search for any item containing the value 'foo':
 * 
 * @code{.cpp}
 * auto c7 = cif::any == "foo";
 * @endcode 
 * 
 * All these conditions can be chained together again:
 * 
 * @code{.cpp}
 * auto c8 = std::move(c3) and std::move(c5);
 * @endcode
 */

namespace cif
{

// --------------------------------------------------------------------
/// let's make life easier, since @ref cif::category is not known yet,
/// we declare a function to access its contents

/**
 * @brief Get the items that can be used as key in conditions for a category
 * 
 * @param cat The category whose items to return
 * @return iset The set of key item names
 */
[[deprecated("use get_category_items instead")]]
iset get_category_fields(const category &cat);

/**
 * @brief Get the items that can be used as key in conditions for a category
 * 
 * @param cat The category whose items to return
 * @return iset The set of key field names
 */
iset get_category_items(const category &cat);

/**
 * @brief Get the item index for item @a col in category @a cat
 * 
 * @param cat The category
 * @param col The name of the item
 * @return uint16_t The index
 */
uint16_t get_item_ix(const category &cat, std::string_view col);

/**
 * @brief Return whether the item @a col in category @a cat has a primitive type of *uchar*
 * 
 * @param cat The category
 * @param col The item name
 * @return true If the primitive type is of type *uchar*
 * @return false If the primitive type is not of type *uchar*
 */
bool is_item_type_uchar(const category &cat, std::string_view col);

// --------------------------------------------------------------------
// some more templates to be able to do querying

namespace detail
{
	struct condition_impl
	{
		virtual ~condition_impl() {}

		virtual condition_impl *prepare(const category &) { return this; }
		virtual bool test(row_handle) const = 0;
		virtual void str(std::ostream &) const = 0;
		virtual std::optional<row_handle> single() const { return {}; };

		virtual bool equals([[maybe_unused]] const condition_impl *rhs) const { return false; }
	};

	struct all_condition_impl : public condition_impl
	{
		bool test(row_handle) const override { return true; }
		void str(std::ostream &os) const override { os << "*"; }
	};

	struct or_condition_impl;
	struct and_condition_impl;
	struct not_condition_impl;
} // namespace detail

/**
 * @brief The interface class for conditions. This uses the bridge pattern,
 * which means the implementation is in the member m_impl
 */
class condition
{
  public:

	/** @cond */
	using condition_impl = detail::condition_impl;
	/** @endcond */

	/**
	 * @brief Construct a new, empty condition object
	 * 
	 */
	condition()
		: m_impl(nullptr)
	{
	}

	/**
	 * @brief Construct a new condition object with implementation @a impl
	 * 
	 * @param impl The implementation to use
	 */
	explicit condition(condition_impl *impl)
		: m_impl(impl)
	{
	}

	condition(const condition &) = delete;

	/**
	 * @brief Construct a new condition object moving the data from @a rhs
	 */
	condition(condition &&rhs) noexcept
		: m_impl(nullptr)
	{
		std::swap(m_impl, rhs.m_impl);
	}

	condition &operator=(const condition &) = delete;

	/**
	 * @brief Assignment operator moving the data from @a rhs
	 */
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

	/**
	 * @brief Prepare the condition to be used on category @a c. This will
	 * take care of setting the correct indices for items e.g.
	 * 
	 * @param c The category this query should act upon
	 */
	void prepare(const category &c);

	/**
	 * @brief This operator returns true if the row referenced by @a r is 
	 * a match for this condition.
	 * 
	 * @param r The reference to a row.
	 * @return true If there is a match
	 * @return false If there is no match
	 */
	bool operator()(row_handle r) const
	{
		assert(this->m_impl != nullptr);
		assert(this->m_prepared);
		return m_impl ? m_impl->test(r) : false;
	}

	/**
	 * @brief Return true if the condition is not empty
	 */
	explicit operator bool() { return not empty(); }

	/**
	 * @brief Return true if the condition is empty, has no condition
	 */
	bool empty() const { return m_impl == nullptr; }

	/**
	 * @brief If the prepare step found out there is only one hit
	 * this single hit can be returned by this method.
	 * 
	 * @return std::optional<row_handle> The result will contain
	 * a row reference if there is a single hit, it will be empty otherwise
	 */
	std::optional<row_handle> single() const
	{
		return m_impl ? m_impl->single() : std::optional<row_handle>();
	}

	friend condition operator||(condition &&a, condition &&b); /**< Return a condition which is the logical OR or condition @a and @b */
	friend condition operator&&(condition &&a, condition &&b); /**< Return a condition which is the logical AND or condition @a and @b */

	/// @cond
	friend struct detail::or_condition_impl;
	friend struct detail::and_condition_impl;
	friend struct detail::not_condition_impl;
	/// @endcond

	/**
	 * @brief Swap two conditions
	 */
	void swap(condition &rhs)
	{
		std::swap(m_impl, rhs.m_impl);
		std::swap(m_prepared, rhs.m_prepared);
	}

	/**
	 * @brief Operator to use to write out a condition to @a os, for debugging purposes
	 * 
	 * @param os The std::ostream to write to
	 * @param cond The condition to write
	 * @return std::ostream& The same as @a os
	 */
	friend std::ostream &operator<<(std::ostream &os, const condition &cond)
	{
		if (cond.m_impl)
			cond.m_impl->str(os);
		return os;
	}

  private:
	void optimise(condition_impl *&impl);

	condition_impl *m_impl;
	bool m_prepared = false;
};

namespace detail
{
	struct key_is_empty_condition_impl : public condition_impl
	{
		key_is_empty_condition_impl(const std::string &item_name)
			: m_item_name(item_name)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			return r[m_item_ix].empty();
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << " IS NULL";
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
	};

	struct key_is_not_empty_condition_impl : public condition_impl
	{
		key_is_not_empty_condition_impl(const std::string &item_name)
			: m_item_name(item_name)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			return not r[m_item_ix].empty();
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << " IS NOT NULL";
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
	};

	struct key_equals_condition_impl : public condition_impl
	{
		key_equals_condition_impl(item &&i)
			: m_item_name(i.name())
			, m_value(std::forward<item>(i).value())
		{
		}

		condition_impl *prepare(const category &c) override;

		bool test(row_handle r) const override
		{
			return m_single_hit.has_value() ? *m_single_hit == r : r[m_item_ix].compare(m_value, m_icase) == 0;
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << (m_icase ? "^ " : " ") << " == " << m_value;
		}

		virtual std::optional<row_handle> single() const override
		{
			return m_single_hit;
		}

		virtual bool equals(const condition_impl *rhs) const override
		{
			if (typeid(*rhs) == typeid(key_equals_condition_impl))
			{
				auto ri = static_cast<const key_equals_condition_impl *>(rhs);
				if (m_single_hit.has_value() or ri->m_single_hit.has_value())
					return m_single_hit == ri->m_single_hit;
				else
					// watch out, both m_item_ix might be the same while item_names might be diffent (in case they both do not exist in the category)
					return m_item_ix == ri->m_item_ix and m_value == ri->m_value and m_item_name == ri->m_item_name;
			}
			return this == rhs;
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
		bool m_icase = false;
		std::string m_value;
		std::optional<row_handle> m_single_hit;
	};

	struct key_equals_or_empty_condition_impl : public condition_impl
	{
		key_equals_or_empty_condition_impl(key_equals_condition_impl *equals)
			: m_item_name(equals->m_item_name)
			, m_value(equals->m_value)
			, m_icase(equals->m_icase)
			, m_single_hit(equals->m_single_hit)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			m_icase = is_item_type_uchar(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			bool result = false;
			if (m_single_hit.has_value())
				result = *m_single_hit == r;
			else
				result = r[m_item_ix].empty() or r[m_item_ix].compare(m_value, m_icase) == 0;
			return result;
		}

		void str(std::ostream &os) const override
		{
			os << '(' << m_item_name << (m_icase ? "^ " : " ") << " == " << m_value << " OR " << m_item_name << " IS NULL)";
		}

		virtual std::optional<row_handle> single() const override
		{
			return m_single_hit;
		}

		virtual bool equals(const condition_impl *rhs) const override
		{
			if (typeid(*rhs) == typeid(key_equals_or_empty_condition_impl))
			{
				auto ri = static_cast<const key_equals_or_empty_condition_impl *>(rhs);
				if (m_single_hit.has_value() or ri->m_single_hit.has_value())
					return m_single_hit == ri->m_single_hit;
				else
					// watch out, both m_item_ix might be the same while item_names might be diffent (in case they both do not exist in the category)
					return m_item_ix == ri->m_item_ix and m_value == ri->m_value and m_item_name == ri->m_item_name;
			}
			return this == rhs;
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
		std::string m_value;
		bool m_icase = false;
		std::optional<row_handle> m_single_hit;
	};

	struct key_equals_number_condition_impl : public condition_impl
	{
		key_equals_number_condition_impl(const std::string &name, double v)
			: m_item_name(name)
			, m_value(v)
		{
		}

		condition_impl *prepare(const category &c) override;

		bool test(row_handle r) const override
		{
			return m_single_hit.has_value() ? *m_single_hit == r : r[m_item_ix].compare(m_value) == 0;
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << " == " << m_value;
		}

		virtual std::optional<row_handle> single() const override
		{
			return m_single_hit;
		}

		virtual bool equals(const condition_impl *rhs) const override
		{
			if (typeid(*rhs) == typeid(key_equals_number_condition_impl))
			{
				auto ri = static_cast<const key_equals_number_condition_impl *>(rhs);
				if (m_single_hit.has_value() or ri->m_single_hit.has_value())
					return m_single_hit == ri->m_single_hit;
				else
					// watch out, both m_item_ix might be the same while item_names might be diffent (in case they both do not exist in the category)
					return m_item_ix == ri->m_item_ix and m_value == ri->m_value and m_item_name == ri->m_item_name;
			}
			return this == rhs;
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
		double m_value;
		std::optional<row_handle> m_single_hit;
	};

	struct key_equals_number_or_empty_condition_impl : public condition_impl
	{
		key_equals_number_or_empty_condition_impl(key_equals_number_condition_impl *equals)
			: m_item_name(equals->m_item_name)
			, m_value(equals->m_value)
			, m_single_hit(equals->m_single_hit)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			bool result = false;
			if (m_single_hit.has_value())
				result = *m_single_hit == r;
			else
				result = r[m_item_ix].empty() or r[m_item_ix].compare(m_value) == 0;
			return result;
		}

		void str(std::ostream &os) const override
		{
			os << '(' << m_item_name << " == " << m_value << " OR " << m_item_name << " IS NULL)";
		}

		virtual std::optional<row_handle> single() const override
		{
			return m_single_hit;
		}

		virtual bool equals(const condition_impl *rhs) const override
		{
			if (typeid(*rhs) == typeid(key_equals_number_or_empty_condition_impl))
			{
				auto ri = static_cast<const key_equals_number_or_empty_condition_impl *>(rhs);
				if (m_single_hit.has_value() or ri->m_single_hit.has_value())
					return m_single_hit == ri->m_single_hit;
				else
					// watch out, both m_item_ix might be the same while item_names might be diffent (in case they both do not exist in the category)
					return m_item_ix == ri->m_item_ix and m_value == ri->m_value and m_item_name == ri->m_item_name;
			}
			return this == rhs;
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
		double m_value;
		std::optional<row_handle> m_single_hit;
	};

	struct key_compare_condition_impl : public condition_impl
	{
		template <typename COMP>
		key_compare_condition_impl(const std::string &item_name, COMP &&comp, const std::string &s)
			: m_item_name(item_name)
			, m_compare(std::move(comp))
			, m_str(s)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			m_icase = is_item_type_uchar(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			return m_compare(r, m_icase);
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << (m_icase ? "^ " : " ") << m_str;
		}

		std::string m_item_name;
		uint16_t m_item_ix = 0;
		bool m_icase = false;
		std::function<bool(row_handle, bool)> m_compare;
		std::string m_str;
	};

	struct key_matches_condition_impl : public condition_impl
	{
		key_matches_condition_impl(const std::string &item_name, const std::regex &rx)
			: m_item_name(item_name)
			, m_item_ix(0)
			, mRx(rx)
		{
		}

		condition_impl *prepare(const category &c) override
		{
			m_item_ix = get_item_ix(c, m_item_name);
			return this;
		}

		bool test(row_handle r) const override
		{
			std::string_view txt = r[m_item_ix].text();
			return std::regex_match(txt.begin(), txt.end(), mRx);
		}

		void str(std::ostream &os) const override
		{
			os << m_item_name << " =~ expression";
		}

		std::string m_item_name;
		uint16_t m_item_ix;
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
			auto &c = r.get_category();

			bool result = false;
			for (auto &f : get_category_items(c))
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
			auto &c = r.get_category();

			bool result = false;
			for (auto &f : get_category_items(c))
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

	// TODO: Optimize and_condition by having a list of sub items.
	// That way you can also collapse multiple _is_ conditions in
	// case they make up an indexed tuple.
	struct and_condition_impl : public condition_impl
	{
		and_condition_impl() = default;

		and_condition_impl(condition &&a, condition &&b)
		{
			if (typeid(*a.m_impl) == typeid(*this))
			{
				and_condition_impl *ai = static_cast<and_condition_impl *>(a.m_impl);

				std::swap(m_sub, ai->m_sub);
				m_sub.emplace_back(std::exchange(b.m_impl, nullptr));
			}
			else if (typeid(*b.m_impl) == typeid(*this))
			{
				and_condition_impl *bi = static_cast<and_condition_impl *>(b.m_impl);

				std::swap(m_sub, bi->m_sub);
				m_sub.emplace_back(std::exchange(a.m_impl, nullptr));
			}
			else
			{
				m_sub.emplace_back(std::exchange(a.m_impl, nullptr));
				m_sub.emplace_back(std::exchange(b.m_impl, nullptr));
			}
		}

		~and_condition_impl()
		{
			for (auto sub : m_sub)
				delete sub;
		}

		condition_impl *prepare(const category &c) override
		{
			for (auto &sub : m_sub)
				sub = sub->prepare(c);
			return this;
		}

		bool test(row_handle r) const override
		{
			bool result = true;

			for (auto sub : m_sub)
			{
				if (sub->test(r))
					continue;

				result = false;
				break;
			}

			return result;
		}

		void str(std::ostream &os) const override
		{
			os << '(';

			bool first = true;
			for (auto sub : m_sub)
			{
				if (first)
					first = false;
				else
					os << " AND ";

				sub->str(os);
			}

			os << ')';
		}

		virtual std::optional<row_handle> single() const override
		{
			std::optional<row_handle> result;

			for (auto sub : m_sub)
			{
				auto s = sub->single();

				if (not result.has_value())
				{
					result = s;
					continue;
				}

				if (s == result)
					continue;

				result.reset();
				break;
			}

			return result;
		}

		static condition_impl *combine_equal(std::vector<and_condition_impl *> &subs, or_condition_impl *oc);

		std::vector<condition_impl *> m_sub;
	};

	struct or_condition_impl : public condition_impl
	{
		or_condition_impl(condition &&a, condition &&b)
		{
			if (typeid(*a.m_impl) == typeid(*this))
			{
				or_condition_impl *ai = static_cast<or_condition_impl *>(a.m_impl);

				std::swap(m_sub, ai->m_sub);
				m_sub.emplace_back(std::exchange(b.m_impl, nullptr));
			}
			else if (typeid(*b.m_impl) == typeid(*this))
			{
				or_condition_impl *bi = static_cast<or_condition_impl *>(b.m_impl);

				std::swap(m_sub, bi->m_sub);
				m_sub.emplace_back(std::exchange(a.m_impl, nullptr));
			}
			else
			{
				m_sub.emplace_back(std::exchange(a.m_impl, nullptr));
				m_sub.emplace_back(std::exchange(b.m_impl, nullptr));
			}
		}

		~or_condition_impl()
		{
			for (auto sub : m_sub)
				delete sub;
		}

		condition_impl *prepare(const category &c) override;

		bool test(row_handle r) const override
		{
			bool result = false;

			for (auto sub : m_sub)
			{
				if (not sub->test(r))
					continue;
				result = true;
				break;
			}

			return result;
		}

		void str(std::ostream &os) const override
		{
			bool first = true;

			os << '(';
			for (auto sub : m_sub)
			{
				if (first)
					first = false;
				else
					os << " OR ";
				sub->str(os);
			}
			os << ')';
		}

		virtual std::optional<row_handle> single() const override
		{
			std::optional<row_handle> result;

			for (auto sub : m_sub)
			{
				auto s = sub->single();

				if (not result.has_value())
				{
					result = s;
					continue;
				}

				if (s == result)
					continue;

				result.reset();
				break;
			}

			return result;
		}

		std::vector<condition_impl *> m_sub;
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

		condition_impl *prepare(const category &c) override
		{
			mA = mA->prepare(c);
			return this;
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

/**
 * @brief Create a condition containing the logical AND of conditions @a a and @a b
 */
inline condition operator and(condition &&a, condition &&b)
{
	if (a.m_impl and b.m_impl)
		return condition(new detail::and_condition_impl(std::move(a), std::move(b)));
	if (a.m_impl)
		return condition(std::move(a));
	return condition(std::move(b));
}

/**
 * @brief Create a condition containing the logical OR of conditions @a a and @a b
 */
inline condition operator or(condition &&a, condition &&b)
{
	if (a.m_impl and b.m_impl)
	{
		if (typeid(*a.m_impl) == typeid(detail::key_equals_condition_impl) and
			typeid(*b.m_impl) == typeid(detail::key_is_empty_condition_impl))
		{
			auto ci = static_cast<detail::key_equals_condition_impl *>(a.m_impl);
			auto ce = static_cast<detail::key_is_empty_condition_impl *>(b.m_impl);

			if (ci->m_item_name == ce->m_item_name)
				return condition(new detail::key_equals_or_empty_condition_impl(ci));
		}
		
		if (typeid(*b.m_impl) == typeid(detail::key_equals_condition_impl) and
				 typeid(*a.m_impl) == typeid(detail::key_is_empty_condition_impl))
		{
			auto ci = static_cast<detail::key_equals_condition_impl *>(b.m_impl);
			auto ce = static_cast<detail::key_is_empty_condition_impl *>(a.m_impl);

			if (ci->m_item_name == ce->m_item_name)
				return condition(new detail::key_equals_or_empty_condition_impl(ci));
		}

		if (typeid(*a.m_impl) == typeid(detail::key_equals_number_condition_impl) and
			typeid(*b.m_impl) == typeid(detail::key_is_empty_condition_impl))
		{
			auto ci = static_cast<detail::key_equals_number_condition_impl *>(a.m_impl);
			auto ce = static_cast<detail::key_is_empty_condition_impl *>(b.m_impl);

			if (ci->m_item_name == ce->m_item_name)
				return condition(new detail::key_equals_number_or_empty_condition_impl(ci));
		}
		
		if (typeid(*b.m_impl) == typeid(detail::key_equals_number_condition_impl) and
				 typeid(*a.m_impl) == typeid(detail::key_is_empty_condition_impl))
		{
			auto ci = static_cast<detail::key_equals_number_condition_impl *>(b.m_impl);
			auto ce = static_cast<detail::key_is_empty_condition_impl *>(a.m_impl);

			if (ci->m_item_name == ce->m_item_name)
				return condition(new detail::key_equals_number_or_empty_condition_impl(ci));
		}

		return condition(new detail::or_condition_impl(std::move(a), std::move(b)));
	}

	if (a.m_impl)
		return condition(std::move(a));

	return condition(std::move(b));
}

/**
 * @brief A helper class to make it possible to search for empty items (NULL)
 * 
 * @code{.cpp}
 * "id"_key == cif::empty_type();
 * @endcode
 */

struct empty_type
{
};

/**
 * @brief A helper to make it possible to have conditions like
 * 
 * @code{.cpp}
 * "id"_key == cif::null;
 * @endcode
 */

inline constexpr empty_type null = empty_type();

/**
 * @brief Class to use in creating conditions, creates a reference to a item or item
 * 
 */
struct key
{
	/**
	 * @brief Construct a new key object using @a item_name as name
	 * 
	 * @param item_name 
	 */
	explicit key(const std::string &item_name)
		: m_item_name(item_name)
	{
	}

	/**
	 * @brief Construct a new key object using @a item_name as name
	 * 
	 * @param item_name 
	 */
	explicit key(const char *item_name)
		: m_item_name(item_name)
	{
	}

	/**
	 * @brief Construct a new key object using @a item_name as name
	 * 
	 * @param item_name 
	 */
	explicit key(std::string_view item_name)
		: m_item_name(item_name)
	{
	}

	key(const key &) = delete;
	key &operator=(const key &) = delete;

	std::string m_item_name; ///< The item name
};

template <typename T>
concept Numeric = ((std::is_floating_point_v<T> or std::is_integral_v<T>) and not std::is_same_v<T, bool>);

/**
 * @brief Operator to create an equals condition based on a key @a key and a numeric value @a v
 */
template <Numeric T>
condition operator==(const key &key, const T &v)
{
	return condition(new detail::key_equals_number_condition_impl(key.m_item_name, v));
}

/**
 * @brief Operator to create an equals condition based on a key @a key and a value @a value
 */
inline condition operator==(const key &key, std::string_view value)
{
	if (not value.empty())
		return condition(new detail::key_equals_condition_impl({ key.m_item_name, value }));
	else
		return condition(new detail::key_is_empty_condition_impl(key.m_item_name));
}

/**
 * @brief Operator to create an equals condition based on a key @a key and a value @a value
 */
template <typename T>
	requires std::is_same_v<T, bool>
inline condition operator==(const key &key, T value)
{
	return condition(new detail::key_equals_condition_impl({ key.m_item_name, value ? "y" : "n" }));
}

/**
 * @brief Operator to create a not equals condition based on a key @a key and a value @a v
 */
template <typename T>
condition operator!=(const key &key, const T &v)
{
	return condition(new detail::not_condition_impl(operator==(key, v)));
}

/**
 * @brief Operator to create a not equals condition based on a key @a key and a value @a value
 */
inline condition operator!=(const key &key, std::string_view value)
{
	return condition(new detail::not_condition_impl(operator==(key, value)));
}

/**
 * @brief Operator to create a greater than condition based on a key @a key and a value @a v
 */
template <Numeric T>
condition operator>(const key &key, const T &v)
{
	std::ostringstream s;
	s << " > " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v) > 0; },
		s.str()));
}

/**
 * @brief Operator to create a greater than or equals condition based on a key @a key and a value @a v
 */
template <Numeric T>
condition operator>=(const key &key, const T &v)
{
	std::ostringstream s;
	s << " >= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v) >= 0; },
		s.str()));
}

/**
 * @brief Operator to create a less than condition based on a key @a key and a value @a v
 */
template <Numeric T>
condition operator<(const key &key, const T &v)
{
	std::ostringstream s;
	s << " < " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v) < 0; },
		s.str()));
}

/**
 * @brief Operator to create a less than or equals condition based on a key @a key and a value @a v
 */
template <Numeric T>
condition operator<=(const key &key, const T &v)
{
	std::ostringstream s;
	s << " <= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v) <= 0; },
		s.str()));
}

/**
 * @brief Operator to create a greater than condition based on a key @a key and a value @a v
 */
inline condition operator>(const key &key, std::string_view v)
{
	std::ostringstream s;
	s << " > " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v, icase) > 0; },
		s.str()));
}

/**
 * @brief Operator to create a greater than or equals condition based on a key @a key and a value @a v
 */
inline condition operator>=(const key &key, std::string_view v)
{
	std::ostringstream s;
	s << " >= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v, icase) >= 0; },
		s.str()));
}

/**
 * @brief Operator to create a less than condition based on a key @a key and a value @a v
 */
inline condition operator<(const key &key, std::string_view v)
{
	std::ostringstream s;
	s << " < " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v, icase) < 0; },
		s.str()));
}

/**
 * @brief Operator to create a less than or equals condition based on a key @a key and a value @a v
 */
inline condition operator<=(const key &key, std::string_view v)
{
	std::ostringstream s;
	s << " <= " << v;

	return condition(new detail::key_compare_condition_impl(
		key.m_item_name, [item_name = key.m_item_name, v](row_handle r, bool icase)
		{ return r[item_name].compare(v, icase) <= 0; },
		s.str()));
}

/**
 * @brief Operator to create a condition based on a key @a key and a regular expression @a rx
 */
inline condition operator==(const key &key, const std::regex &rx)
{
	return condition(new detail::key_matches_condition_impl(key.m_item_name, rx));
}

/**
 * @brief Operator to create a condition based on a key @a key which should be empty/null
 */
inline condition operator==(const key &key, const empty_type &)
{
	return condition(new detail::key_is_empty_condition_impl(key.m_item_name));
}

/**
 * @brief Operator to create a condition based on a key @a key which should be not empty/null
 */
inline condition operator!=(const key &key, const empty_type &)
{
	return condition(new detail::key_is_not_empty_condition_impl(key.m_item_name));
}

/**
 * @brief Create a condition to search any item for a value @a v if @a v contains a value
 * compare to null if not.
 */
template <typename T>
condition operator==(const key &key, const std::optional<T> &v)
{
	if (v.has_value())
		return condition(new detail::key_equals_condition_impl({ key.m_item_name, *v }));
	else
		return condition(new detail::key_is_empty_condition_impl(key.m_item_name));
}

/**
 * @brief Create a condition to search any item for a value @a v if @a v contains a value
 * compare to null if not.
 */
template <typename T>
condition operator!=(const key &key, const std::optional<T> &v)
{
	if (v.has_value())
		return condition(new detail::not_condition_impl(condition(new detail::key_equals_condition_impl({ key.m_item_name, *v }))));
	else
		return condition(new detail::not_condition_impl(condition(new detail::key_is_empty_condition_impl(key.m_item_name))));
}

/**
 * @brief Operator to create a boolean opposite of the condition in @a rhs
 */
inline condition operator not(condition &&rhs)
{
	return condition(new detail::not_condition_impl(std::move(rhs)));
}

/** @cond */
struct any_type
{
};
/** @endcond */

/**
 * @brief A helper for any item constructs
 */
inline constexpr any_type any = any_type{};

/**
 * @brief Create a condition to search any item for a value @a v
 */
template <typename T>
condition operator==(const any_type &, const T &v)
{
	return condition(new detail::any_is_condition_impl<T>(v));
}

/**
 * @brief Create a condition to search any item for a regular expression @a rx
 */
inline condition operator==(const any_type &, const std::regex &rx)
{
	return condition(new detail::any_matches_condition_impl(rx));
}

/**
 * @brief Create a condition to return all rows
 */
inline condition all()
{
	return condition(new detail::all_condition_impl());
}

namespace literals
{
	/**
	 * @brief Return a cif::key for the item name @a text
	 * 
	 * @param text The name of the item
	 * @param length The length of @a text
	 * @return key The cif::key created
	 */
	inline key operator""_key(const char *text, std::size_t length)
	{
		return key(std::string(text, length));
	}
} // namespace literals

} // namespace cif
