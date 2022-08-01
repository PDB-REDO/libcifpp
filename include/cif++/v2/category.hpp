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

#include "iterator.hpp"
#include "row.hpp"
#include "validate.hpp"

namespace cif::v2
{

// --------------------------------------------------------------------

template <
	typename Alloc = std::allocator<std::byte>>
class category_t
{
  private:
	// --------------------------------------------------------------------
	// Internal storage, strictly forward linked list with minimal space
	// requirements. Strings of size 7 or shorter are stored internally.
	// Typically, more than 99% of the strings in an mmCIF file are less
	// than 8 bytes in length.

	struct item_value
	{
		item_value(uint16_t column_ix, uint16_t length)
			: m_next(nullptr)
			, m_column_ix(column_ix)
			, m_length(length)
		{
		}

		item_value() = delete;
		item_value(const item_value &) = delete;
		item_value &operator=(const item_value &) = delete;

		item_value *m_next;
		uint16_t m_column_ix;
		uint16_t m_length;
		union
		{
			char m_local_data[8];
			char *m_data;
		};

		static constexpr size_t kBufferSize = sizeof(m_local_data);

		std::string_view text() const
		{
			return {m_length >= kBufferSize ? m_data : m_local_data, m_length};
		}

		const char *c_str() const
		{
			return m_length >= kBufferSize ? m_data : m_local_data;
		}
	};

	static_assert(sizeof(item_value) == 24, "sizeof(item_value) should be 24 bytes");

  public:
	using allocator_type = Alloc;

	allocator_type get_allocator() const
	{
		return m_allocator;
	}

	template <typename>
	friend class row_handle;

	template <typename, typename...>
	friend class iterator_impl;

	using value_type = row_handle<category_t>;
	using reference = value_type;
	using const_reference = const value_type;
	using iterator = iterator_impl<category_t>;
	using const_iterator = iterator_impl<const category_t>;

	class row
	{
	  public:
		row()
			: m_next(nullptr)
			, m_head(nullptr)
			, m_tail(nullptr)
		{
		}

		row(row *next, item_value *data)
			: m_next(next)
			, m_head(data)
			, m_tail(data)
		{
			auto n = m_tail ? m_tail->m_next : nullptr;
			while (n != nullptr)
			{
				m_tail = n;
				n = n->m_next;
			}
		}

	  private:
		template <typename>
		friend class item_handle;

		template <typename, typename...>
		friend class iterator_impl;

		friend class category_t;

		row *m_next = nullptr;
		item_value *m_head = nullptr, *m_tail = nullptr;
	};

	category_t() = default;

	category_t(std::string_view name, const allocator_type &alloc = allocator_type())
		: m_allocator(alloc)
		, m_name(name)
	{
	}

	category_t(const category_t &) = default;

	category_t(category_t &&) = default;

	template <typename Alloc2>
	category_t(const category_t &c, const Alloc2 &a)
		: m_allocator(a)
		, m_name(c.m_name)
	{
	}

	template <typename Alloc2>
	category_t(category_t &&c, const Alloc2 &a)
		: m_allocator(a)
		, m_name(std::move(c.m_name))
	{
	}

	category_t &operator=(const category_t &) = default;
	category_t &operator=(category_t &&) = default;

	~category_t()
	{
		auto r = m_head;
		while (r != nullptr)
		{
			auto i = r->m_head;
			while (i != nullptr)
			{
				auto ti = i->m_next;
				delete_item(i);
				i = ti;
			}

			auto t = r->m_next;
			delete_row(r);
			r = t;
		}
	}

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	reference front()
	{
		return {*this, *m_head};
	}

	const_reference front() const
	{
		return {*this, *m_head};
	}

	reference back()
	{
		return {*this, *m_tail};
	}

	const_reference back() const
	{
		return {*this, *m_tail};
	}

	iterator begin()
	{
		return {*this, m_head};
	}

	iterator end()
	{
		return {*this, nullptr};
	}

	const_iterator begin() const
	{
		return {*this, m_head};
	}

	const_iterator end() const
	{
		return {*this, nullptr};
	}

	const_iterator cbegin() const
	{
		return {*this, m_head};
	}

	const_iterator cend() const
	{
		return {*this, nullptr};
	}

	// --------------------------------------------------------------------

	template <typename... Ts, typename... Ns>
	iterator_proxy<const category_t, Ts...> rows(Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<const category_t, Ts...>(*this, begin(), {names...});
	}

	template <typename... Ts, typename... Ns>
	iterator_proxy<category_t, Ts...> rows(Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<category_t, Ts...>(*this, begin(), {names...});
	}

	// --------------------------------------------------------------------

	iterator emplace(std::initializer_list<item> items)
	{
		return this->emplace(items.begin(), items.end());
	}

	template <typename ItemIter>
	iterator emplace(ItemIter b, ItemIter e)
	{
		// First, make sure all mandatory fields are supplied
		// if (mCatValidator != nullptr and b != e)
		// {
		// 	for (auto &col : m_columns)
		// 	{
		// 		auto iv = mCatValidator->getValidatorForItem(col.mName);

		// 		if (iv == nullptr)
		// 			continue;

		// 		bool seen = false;

		// 		for (auto v = b; v != e; ++v)
		// 		{
		// 			if (iequals(v->name(), col.mName))
		// 			{
		// 				seen = true;
		// 				break;
		// 			}
		// 		}

		// 		if (not seen and iv->mMandatory)
		// 			throw std::runtime_error("missing mandatory field " + col.mName + " for Category " + mName);
		// 	}

		// 	if (mIndex != nullptr)
		// 	{
		// 		std::unique_ptr<ItemRow> nr(new ItemRow{nullptr, this, nullptr});
		// 		Row r(nr.get());
		// 		auto keys = keyFields();

		// 		for (auto v = b; v != e; ++v)
		// 		{
		// 			if (keys.count(v->name()))
		// 				r.assign(v->name(), v->value(), true);
		// 		}

		// 		auto test = mIndex->find(nr.get());
		// 		if (test != nullptr)
		// 		{
		// 			if (VERBOSE > 1)
		// 				std::cerr << "Not inserting new record in " << mName << " (duplicate Key)" << std::endl;
		// 			result = test;
		// 			isNew = false;
		// 		}
		// 	}
		// }

		row *r = this->create_row();

		try
		{
			for (auto i = b; i != e; ++i)
			{
				item_value *new_item = this->create_item(*i);

				if (r->m_head == nullptr)
					r->m_head = r->m_tail = new_item;
				else
					r->m_tail = r->m_tail->m_next = new_item;
			}
		}
		catch (...)
		{
			if (r != nullptr)
				this->delete_row(r);
			throw;
		}

		if (m_head == nullptr)
			m_head = m_tail = r;
		else
			m_tail = m_tail->m_next = r;

		return {*this, r};

		// result = r;

		// if (mIndex != nullptr)
		// 	mIndex->insert(nr);
	}

	// --------------------------------------------------------------------
	/// \brief Return the index number for \a column_name

	uint16_t get_column_ix(std::string_view column_name) const
	{
		uint16_t result;

		for (result = 0; result < m_columns.size(); ++result)
		{
			if (iequals(column_name, m_columns[result].m_name))
				break;
		}

		// if (VERBOSE > 0 and result == m_columns.size() and mCatValidator != nullptr) // validate the name, if it is known at all (since it was not found)
		// {
		// 	auto iv = mCatValidator->getValidatorForItem(name);
		// 	if (iv == nullptr)
		// 		std::cerr << "Invalid name used '" << name << "' is not a known column in " + mName << std::endl;
		// }

		return result;
	}

	uint16_t add_column(std::string_view column_name)
	{
		using namespace std::literals;

		size_t result = get_column_ix(column_name);

		if (result == m_columns.size())
		{
			const ValidateItem *itemValidator = nullptr;

			// if (mCatValidator != nullptr)
			// {
			// 	itemValidator = mCatValidator->getValidatorForItem(column_name);
			// 	if (itemValidator == nullptr)
			// 		mValidator->reportError("tag " + std::string(column_name) + " not allowed in Category " + mName, false);
			// }

			m_columns.emplace_back(column_name, itemValidator);
		}

		return result;
	}

  private:
	using char_allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<char>;
	using char_allocator_traits = std::allocator_traits<char_allocator_type>;

	using item_allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<item_value>;
	using item_allocator_traits = std::allocator_traits<item_allocator_type>;

	item_allocator_traits::pointer get_item()
	{
		item_allocator_type ia(get_allocator());
		return item_allocator_traits::allocate(ia, 1);
	}

	item_value *create_item(uint16_t column_ix, std::string_view text)
	{
		size_t text_length = text.length();

		if (text_length + 1 > std::numeric_limits<uint16_t>::max())
			throw std::runtime_error("libcifpp does not support string lengths longer than " + std::to_string(std::numeric_limits<uint16_t>::max()) + " bytes");

		auto p = this->get_item();
		item_allocator_type ia(get_allocator());
		item_allocator_traits::construct(ia, p, column_ix, static_cast<uint16_t>(text_length));

		char_allocator_type ca(get_allocator());

		char *data;
		if (text_length >= item_value::kBufferSize)
			data = p->m_data = char_allocator_traits::allocate(ca, text_length + 1);
		else
			data = p->m_local_data;

		std::copy(text.begin(), text.end(), data);
		data[text_length] = 0;

		return p;
	}

	item_value *create_item(const item &i)
	{
		uint16_t ix = add_column(i.name());
		return create_item(ix, i.value());
	}

	void delete_item(item_value *iv)
	{
		if (iv->m_length >= item_value::kBufferSize)
		{
			char_allocator_type ca(get_allocator());
			char_allocator_traits::deallocate(ca, iv->m_data, iv->m_length + 1);
		}

		item_allocator_type ia(get_allocator());
		item_allocator_traits::destroy(ia, iv);
		item_allocator_traits::deallocate(ia, iv, 1);
	}

	using row_allocator_type = typename std::allocator_traits<allocator_type>::template rebind_alloc<row>;
	using row_allocator_traits = std::allocator_traits<row_allocator_type>;

	row_allocator_traits::pointer get_row()
	{
		row_allocator_type ra(get_allocator());
		return row_allocator_traits::allocate(ra, 1);
	}

	template <typename... Args>
	row *create_row(Args... args)
	{
		auto p = this->get_row();
		row_allocator_type ra(get_allocator());
		row_allocator_traits::construct(ra, p, std::forward<Args>(args)...);
		return p;
	}

	void delete_row(row *r)
	{
		row_allocator_type ra(get_allocator());
		row_allocator_traits::destroy(ra, r);
		row_allocator_traits::deallocate(ra, r, 1);
	}

	struct item_column
	{
		std::string m_name;
		const ValidateItem *m_validator;

		item_column(std::string_view name, const ValidateItem *validator)
			: m_name(name)
			, m_validator(validator)
		{
		}
	};

	allocator_type m_allocator;
	std::string m_name;
	std::vector<item_column, typename std::allocator_traits<allocator_type>::template rebind_alloc<item_column>> m_columns;
	row *m_head = nullptr, *m_tail = nullptr;
};

using category = category_t<>;

} // namespace cif::v2