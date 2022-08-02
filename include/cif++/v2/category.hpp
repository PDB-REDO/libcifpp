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

#include <cif++/v2/condition.hpp>
#include <cif++/v2/iterator.hpp>
#include <cif++/v2/row.hpp>
#include <cif++/v2/validate.hpp>
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

		row() = default;

	  private:
		template <typename>
		friend class item_handle;

		template <typename, typename...>
		friend class iterator_impl;

		friend class category_t;

		void append(item_value *iv)
		{
			if (m_head == nullptr)
				m_head = m_tail = iv;
			else
				m_tail = m_tail->m_next = iv;
		}

		row *m_next = nullptr;
		item_value *m_head = nullptr, *m_tail = nullptr;
	};

	category_t() = default;

	category_t(std::string_view name, const allocator_type &alloc = allocator_type())
		: m_allocator(alloc)
		, m_name(name)
	{
	}

	category_t(const category_t &rhs)
		: m_allocator(std::allocator_traits<allocator_type>::select_on_container_copy_construction(rhs.get_allocator()))
		, m_name(rhs.m_name)
		, m_columns(rhs.m_columns)
	{
		for (auto r = rhs.m_head; r != nullptr; r = r->m_next)
			insert_impl(end(), clone_row(*r));
	}

	category_t(category_t &&rhs)
		: m_allocator(std::move(rhs.m_allocator))
		, m_name(std::move(rhs.m_name))
		, m_columns(std::move(rhs.m_columns))
		, m_head(rhs.m_head)
		, m_tail(rhs.m_tail)
	{
		rhs.m_head = nullptr;
		rhs.m_tail = nullptr;
	}

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

	category_t &operator=(const category_t &rhs)
	{
		if (this != &rhs)
		{
			if (not empty())
				clear();

			m_allocator = std::allocator_traits<allocator_type>::select_on_container_copy_construction(rhs.get_allocator());
			m_name = rhs.m_name;
			m_columns = rhs.m_columns;

			for (auto r = rhs.m_head; r != nullptr; r = r->m_next)
				insert_impl(cend(), clone_row(*r));
		}

		return *this;
	}

	category_t &operator=(category_t &&rhs)
	{
		if (this != &rhs)
		{
			if (not empty())
				clear();

			m_allocator = std::move(rhs.m_allocator);
			m_name = std::move(rhs.m_name);
			m_columns = std::move(rhs.m_columns);

			m_head = rhs.m_head;
			m_tail = rhs.m_tail;

			rhs.m_head = rhs.m_tail = nullptr;
		}

		return *this;
	}

	~category_t()
	{
		clear();
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

	size_t size() const
	{
		return std::distance(cbegin(), cend());
	}

	bool empty() const
	{
		return m_head == nullptr;
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

	void insert(const_iterator pos, const row &row)
	{
		insert_impl(pos, row);
	}

	void insert(const_iterator pos, row &&row)
	{
		insert_impl(pos, std::move(row));
	}

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
				r->append(new_item);
			}
		}
		catch (...)
		{
			if (r != nullptr)
				this->delete_row(r);
			throw;
		}

		return insert_impl(cend(), r);

		// result = r;

		// if (mIndex != nullptr)
		// 	mIndex->insert(nr);
	}

	void clear()
	{
		auto i = m_head;
		while (i != nullptr)
		{
			auto t = i;
			i = i->m_next;
			delete_row(t);
		}

		m_head = m_tail = nullptr;
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
			const ValidateItem *item_validator = nullptr;

			// if (mCatValidator != nullptr)
			// {
			// 	item_validator = mCatValidator->getValidatorForItem(column_name);
			// 	if (item_validator == nullptr)
			// 		m_validator->reportError("tag " + std::string(column_name) + " not allowed in Category " + mName, false);
			// }

			m_columns.emplace_back(column_name, item_validator);
		}

		return result;
	}

  private:
	void update_value(row *row, size_t column, std::string_view value, bool updateLinked, bool validate = true)
	{
		auto &col = m_columns[column];

		const char *oldValue = nullptr;
		for (auto iv = row->m_head; iv != nullptr; iv = iv->m_next)
		{
			assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));

			if (iv->m_column_ix == column)
			{
				oldValue = iv->c_str();
				break;
			}
		}

		if (oldValue != nullptr and value == oldValue) // no need to update
			return;

		std::string oldStrValue = oldValue ? oldValue : "";

		// // check the value
		// if (col.m_validator and validate)
		// 	(*col.m_validator)(value);

		// If the field is part of the Key for this Category, remove it from the index
		// before updating

		bool reinsert = false;

		// if (updateLinked and // an update of an Item's value
		// 	cat->mIndex != nullptr and cat->keyFieldsByIndex().count(column))
		// {
		// 	reinsert = cat->mIndex->find(mData);
		// 	if (reinsert)
		// 		cat->mIndex->erase(mData);
		// }

		// first remove old value with cix

		if (row->m_head == nullptr)
			; // nothing to do
		else if (row->m_head->m_column_ix == column)
		{
			auto iv = row->m_head;
			row->m_head = iv->m_next;
			iv->m_next = nullptr;
			delete_item(iv);
		}
		else
		{
			for (auto iv = row->m_head; iv->m_next != nullptr; iv = iv->m_next)
			{
				if (iv->m_next->m_column_ix != column)
					continue;

				auto nv = iv->m_next;
				iv->m_next = nv->m_next;
				nv->m_next = nullptr;
				delete_item(nv);

				break;
			}
		}

		if (not value.empty())
		{
			auto nv = create_item(column, value);

			if (row->m_head == nullptr)
				row->m_head = nv;
			else
			{
				auto iv = row->m_head;
				while (iv->m_next != nullptr)
					iv = iv->m_next;
				iv->m_next = nv;
			}
		}

		// if (reinsert)
		// 	cat->mIndex->insert(mData);

		// // see if we need to update any child categories that depend on this value
		// auto iv = col.m_validator;
		// if (not skipUpdateLinked and iv != nullptr and mCascade)
		// {
		// 	for (auto &&[childCat, linked] : cat->mChildLinks)
		// 	{
		// 		if (find(linked->mParentKeys.begin(), linked->mParentKeys.end(), iv->mTag) == linked->mParentKeys.end())
		// 			continue;

		// 		Condition cond;
		// 		std::string childTag;

		// 		for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
		// 		{
		// 			std::string pk = linked->mParentKeys[ix];
		// 			std::string ck = linked->mChildKeys[ix];

		// 			// TODO add code to *NOT* test mandatory fields for Empty

		// 			if (pk == iv->mTag)
		// 			{
		// 				childTag = ck;
		// 				cond = std::move(cond) && Key(ck) == oldStrValue;
		// 			}
		// 			else
		// 			{
		// 				const char *pk_value = (*this)[pk].c_str();
		// 				if (*pk_value == 0)
		// 					cond = std::move(cond) && Key(ck) == Empty();
		// 				else
		// 					cond = std::move(cond) && ((Key(ck) == pk_value) or Key(ck) == Empty());
		// 			}
		// 		}

		// 		auto rows = childCat->find(std::move(cond));
		// 		if (rows.empty())
		// 			continue;

		// 		// if (cif::VERBOSE > 2)
		// 		// {
		// 		// 	std::cerr << "Parent: " << linked->mParentCategory << " Child: " << linked->mChildCategory << std::endl
		// 		// 			  << cond << std::endl;
		// 		// }

		// 		// Now, suppose there are already rows in child that conform to the new value,
		// 		// we then skip this renam

		// 		Condition cond_n;

		// 		for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
		// 		{
		// 			std::string pk = linked->mParentKeys[ix];
		// 			std::string ck = linked->mChildKeys[ix];

		// 			// TODO add code to *NOT* test mandatory fields for Empty

		// 			if (pk == iv->mTag)
		// 				cond_n = std::move(cond_n) && Key(ck) == value;
		// 			else
		// 			{
		// 				const char *pk_value = (*this)[pk].c_str();
		// 				if (*pk_value == 0)
		// 					cond_n = std::move(cond_n) && Key(ck) == Empty();
		// 				else
		// 					cond_n = std::move(cond_n) && ((Key(ck) == pk_value) or Key(ck) == Empty());
		// 			}
		// 		}

		// 		auto rows_n = childCat->find(std::move(cond_n));
		// 		if (not rows_n.empty())
		// 		{
		// 			if (cif::VERBOSE > 0)
		// 				std::cerr << "Will not rename in child category since there are already rows that link to the parent" << std::endl;

		// 			continue;
		// 		}

		// 		for (auto &cr : rows)
		// 			cr.assign(childTag, value, false);
		// 	}
		// }		
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

	row *create_row()
	{
		auto p = this->get_row();
		row_allocator_type ra(get_allocator());
		row_allocator_traits::construct(ra, p);
		return p;
	}

	row *clone_row(const row &r)
	{
		row *result = create_row();

		try
		{
			for (auto i = r.m_head; i != nullptr; i = i->m_next)
			{
				item_value *v = create_item(i->m_column_ix, i->text());
				result->append(v);
			}
		}
		catch (...)
		{
			delete_row(result);
			throw;
		}
		
		return result;
	}

	void delete_row(row *r)
	{
		if (r != nullptr)
		{
			auto i = r->m_head;
			while (i != nullptr)
			{
				auto t = i;
				i = i->m_next;
				delete_item(t);
			}

			row_allocator_type ra(get_allocator());
			row_allocator_traits::destroy(ra, r);
			row_allocator_traits::deallocate(ra, r, 1);
		}
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

	// proxy methods for every insertion
	iterator insert_impl(const_iterator pos, row *n)
	{
		assert(n != nullptr);
		assert(n->m_next == nullptr);

		if (n == nullptr)
			throw std::runtime_error("Invalid pointer passed to insert");

		// insert at end, most often this is the case
		if (pos.m_current == nullptr)
		{
			if (m_head == nullptr)
				m_tail = m_head = n;
			else
				m_tail = m_tail->m_next = n;
		}
		else
		{
			assert(m_head != nullptr);

			if (pos.m_current == m_head)
				m_head = n->m_next = m_head;
			else
				n = n->m_next = m_head->m_next;
		}

		return iterator(*this, n);
	}

	iterator erase_impl(const_iterator pos)
	{
		if (pos == cend())
			return pos;

		row *n = const_cast<row *>(&*pos);
		row *cur;

		if (m_head == n)
		{
			m_head = static_cast<row *>(m_head->m_next);
			if (m_head != nullptr)
				m_head->m_prev = nullptr;
			else
				m_tail = nullptr;

			n->m_next = n->m_prev = n->m_parent = nullptr;
			delete_row(n);

			cur = m_head;
		}
		else
		{
			cur = static_cast<row *>(n->m_next);

			if (m_tail == n)
				m_tail = static_cast<row *>(n->m_prev);

			row *p = m_head;
			while (p != nullptr and p->m_next != n)
				p = p->m_next;

			if (p != nullptr and p->m_next == n)
			{
				p->m_next = n->m_next;
				if (p->m_next != nullptr)
					p->m_next->m_prev = p;
				n->m_next = nullptr;
			}
			else
				throw std::runtime_error("remove for a row not found in the list");

			delete_row(n);
		}

		return iterator(*this, cur);
	}

	allocator_type m_allocator;
	std::string m_name;
	std::vector<item_column, typename std::allocator_traits<allocator_type>::template rebind_alloc<item_column>> m_columns;
	row *m_head = nullptr, *m_tail = nullptr;
};

using category = category_t<>;

} // namespace cif::v2