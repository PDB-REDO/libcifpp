/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
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

#include <filesystem>
#include <forward_list>
#include <list>
#include <map>
#include <string>
#include <scoped_allocator>

#include "cif++/CifUtils.hpp"

#include "cif++/v2/item.hpp"

namespace cif::v2
{

// template <typename Alloc = std::allocator<void>>
// class item
// {
//   public:
// 	item() = default;

// 	item(std::string_view name, char value)
// 		: m_name(name)
// 		, m_value(value)
// 	{
// 	}

// #if defined(__cpp_lib_format)
// 	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
// 	item(std::string_view name, const T &value, int precision)
// 		: m_name(name)
// 		, m_value(std::format(".{}f", value, precision))
// 	{
// 	}
// #endif

// 	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
// 	item(const std::string_view name, const T &value)
// 		: m_name(name)
// 		, m_value(std::to_string(value))
// 	{
// 	}

// 	item(const std::string_view name, const std::string_view value)
// 		: m_name(name)
// 		, m_value(value)
// 	{
// 	}

// 	item(const item &rhs) = default;

// 	item(item &&rhs) noexcept = default;

// 	item &operator=(const item &rhs) = default;

// 	item &operator=(item &&rhs) noexcept = default;

// 	const std::string &name() const { return m_name; }
// 	const std::string &value() const { return m_value; }

// 	void value(const std::string &v) { m_value = v; }

// 	/// \brief empty means either null or unknown
// 	bool empty() const { return m_value.empty(); }

// 	/// \brief returns true if the field contains '.'
// 	bool is_null() const { return m_value == "."; }

// 	/// \brief returns true if the field contains '?'
// 	bool is_unknown() const { return m_value == "?"; }

// 	size_t length() const { return m_value.length(); }
// 	const char *c_str() const { return m_value.c_str(); }

//   private:
// 	std::string m_name;
// 	std::string m_value;
// };

// using item = item<>;

class item
{
  public:
	item() = default;

	item(std::string_view name, char value)
		: m_name(name)
		, m_value({ value })
	{
	}

	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	item(std::string_view name, const T &value, int precision)
		: m_name(name)
#if defined(__cpp_lib_format)
		, m_value(std::format(".{}f", value, precision))
#endif
	{
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
		, m_value(std::to_string(value))
	{
	}

	item(const std::string_view name, const std::string_view value)
		: m_name(name)
		, m_value(value)
	{
	}

	item(const item &rhs) = default;

	item(item &&rhs) noexcept = default;

	item &operator=(const item &rhs) = default;

	item &operator=(item &&rhs) noexcept = default;

	const std::string &name() const { return m_name; }
	const std::string &value() const { return m_value; }

	void value(const std::string &v) { m_value = v; }

	/// \brief empty means either null or unknown
	bool empty() const { return m_value.empty(); }

	/// \brief returns true if the field contains '.'
	bool is_null() const { return m_value == "."; }

	/// \brief returns true if the field contains '?'
	bool is_unknown() const { return m_value == "?"; }

	size_t length() const { return m_value.length(); }
	const char *c_str() const { return m_value.c_str(); }

  private:
	std::string m_name;
	std::string m_value;
};


// --------------------------------------------------------------------

template <typename Alloc = std::allocator<std::byte>>
class row_t
{
	using byte_allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<std::byte>;

	struct AllocTraitsImpl : std::allocator_traits<byte_allocator_type>
	{
		using base_type = std::allocator_traits<byte_allocator_type>;

		static constexpr typename base_type::pointer
		allocate(byte_allocator_type &a, typename base_type::size_type n)
		{
			return base_type::allocate(a, n);
		}
	};


	using AllocTraits = AllocTraitsImpl;

  public:

	using allocator_type = byte_allocator_type;

	row_t(const allocator_type &alloc = allocator_type())
		: m_allocator(alloc) {}

	row_t(const row_t &) = default;

	row_t(row_t &&) = default;

	// template<typename Alloc2>
	// row_t(const row_t &r, const Alloc2 &a)
	// 	: m_allocator(a)
	// {
	// }

	// template<typename Alloc2>
	// row_t(row_t &&r, const Alloc2 &a)
	// 	: m_allocator(a)
	// {
	// }

	row_t &operator=(const row_t &) = default;
	row_t &operator=(row_t &&) = default;

	row_t(std::initializer_list<item> items, const allocator_type &alloc = allocator_type())
		: m_allocator(alloc)
	{

	}

	item_handle<row_t> operator[](uint32_t column_ix)
	{
		return item_handle<row_t>(column_ix, *this);
	}

	const item_handle<const row_t> operator[](uint32_t column_ix) const
	{
		return item_handle<const row_t>(column_ix, *this);
	}

	item_handle<row_t> operator[](std::string_view column_name)
	{
		return item_handle<row_t>(column_name, get_column_ix(column_name), *this);
	}

	const item_handle<const row_t> operator[](std::string_view column_name) const
	{
		return item_handle<const row_t>(column_name, get_column_ix(column_name), *this);
	}




	// --------------------------------------------------------------------

  private:

	uint32_t get_column_ix(std::string_view name) const
	{
		return 0;
	}

	template<typename R> friend class item_handle;


	allocator_type m_allocator;

	item_value *m_head = nullptr, *m_tail = nullptr;
	size_t m_size = 0;
};

using row = row_t<>;

// --------------------------------------------------------------------

template <
	typename Row = row,
	typename Alloc = std::allocator<std::byte>>
class category_t : public std::forward_list<Row, std::scoped_allocator_adaptor<std::allocator<Row>>>
{
  public:
	using value_type = Row;
	using base_type = std::forward_list<Row, std::scoped_allocator_adaptor<std::allocator<Row>>>;
	using allocator_type = Alloc;

	category_t() = default;

	category_t(std::string_view name, const allocator_type &alloc = allocator_type())
		: base_type(alloc)
		, m_name(name)
	{
	}

	category_t(const category_t &) = default;

	category_t(category_t &&) = default;

	template <typename Alloc2>
	category_t(const category_t &c, const Alloc2 &a)
		: base_type(c, a)
		, m_name(c.m_name)
	{
	}

	template <typename Alloc2>
	category_t(category_t &&c, const Alloc2 &a)
		: base_type(std::move(c), a)
		, m_name(c.m_name)
	{
	}

	category_t &operator=(const category_t &) = default;
	category_t &operator=(category_t &&) = default;

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	void emplace(value_type &&row);

	void emplace(std::initializer_list<item> items)
	{
		this->emplace_back(value_type{items, this->get_allocator()});
	}

	// void write(std::ostream &os, const std::vector<size_t> &order, bool includeEmptyColumns)
	// {
	// 	if (empty())
	// 		return;

	// 	// If there are multiple rows in this category, we need a _loop
	// 	if (size() == 1)
	// 	{
	// 		os << "loop_" << std::endl;

	// 		std::vector<size_t> columnWidths;

	// 		for (auto cix : order)
	// 		{
	// 			auto &col = mColumns[cix];
	// 			os << '_' << mName << '.' << col.mName << ' ' << std::endl;
	// 			columnWidths.push_back(2);
	// 		}

	// 		for (auto Row = mHead; Row != nullptr; Row = Row->mNext)
	// 		{
	// 			for (auto v = Row->mValues; v != nullptr; v = v->mNext)
	// 			{
	// 				if (strchr(v->mText, '\n') == nullptr)
	// 				{
	// 					size_t l = strlen(v->mText);

	// 					if (not isUnquotedString(v->mText))
	// 						l += 2;

	// 					if (l > 132)
	// 						continue;

	// 					if (columnWidths[v->mColumnIndex] < l + 1)
	// 						columnWidths[v->mColumnIndex] = l + 1;
	// 				}
	// 			}
	// 		}

	// 		for (auto Row = mHead; Row != nullptr; Row = Row->mNext) // loop over rows
	// 		{
	// 			size_t offset = 0;

	// 			for (size_t cix : order)
	// 			{
	// 				size_t w = columnWidths[cix];

	// 				std::string s;
	// 				for (auto iv = Row->mValues; iv != nullptr; iv = iv->mNext)
	// 				{
	// 					if (iv->mColumnIndex == cix)
	// 					{
	// 						s = iv->mText;
	// 						break;
	// 					}
	// 				}

	// 				if (s.empty())
	// 					s = "?";

	// 				size_t l = s.length();
	// 				if (not isUnquotedString(s.c_str()))
	// 					l += 2;
	// 				if (l < w)
	// 					l = w;

	// 				if (offset + l > 132 and offset > 0)
	// 				{
	// 					os << std::endl;
	// 					offset = 0;
	// 				}

	// 				offset = detail::writeValue(os, s, offset, w);

	// 				if (offset > 132)
	// 				{
	// 					os << std::endl;
	// 					offset = 0;
	// 				}
	// 			}

	// 			if (offset > 0)
	// 				os << std::endl;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		// first find the indent level
	// 		size_t l = 0;

	// 		for (auto &col : mColumns)
	// 		{
	// 			std::string tag = '_' + mName + '.' + col.mName;

	// 			if (l < tag.length())
	// 				l = tag.length();
	// 		}

	// 		l += 3;

	// 		for (size_t cix : order)
	// 		{
	// 			auto &col = mColumns[cix];

	// 			os << '_' << mName << '.' << col.mName << std::string(l - col.mName.length() - mName.length() - 2, ' ');

	// 			std::string s;
	// 			for (auto iv = mHead->mValues; iv != nullptr; iv = iv->mNext)
	// 			{
	// 				if (iv->mColumnIndex == cix)
	// 				{
	// 					s = iv->mText;
	// 					break;
	// 				}
	// 			}

	// 			if (s.empty())
	// 				s = "?";

	// 			size_t offset = l;
	// 			if (s.length() + l >= kMaxLineLength)
	// 			{
	// 				os << std::endl;
	// 				offset = 0;
	// 			}

	// 			if (detail::writeValue(os, s, offset, 1) != 0)
	// 				os << std::endl;
	// 		}
	// 	}

	// 	os << "# " << std::endl;
	// }

	void write(std::ostream &os) const
	{
		// std::vector<size_t> order(mColumns.size());
		// iota(order.begin(), order.end(), 0);
		// write(os, order, false);

		os << '#' << m_name << std::endl;
		for (auto &r : *this)
		{
			for (auto &f : r)
				os << '_' << m_name << '.' << f.name() << ' ' << f.value() << std::endl;
		}
	}

	// void Category::write(std::ostream &os, const std::vector<std::string> &columns)
	// {
	// 	// make sure all columns are present
	// 	for (auto &c : columns)
	// 		addColumn(c);

	// 	std::vector<size_t> order;
	// 	order.reserve(mColumns.size());

	// 	for (auto &c : columns)
	// 		order.push_back(getColumnIndex(c));

	// 	for (size_t i = 0; i < mColumns.size(); ++i)
	// 	{
	// 		if (std::find(order.begin(), order.end(), i) == order.end())
	// 			order.push_back(i);
	// 	}

	// 	write(os, order, true);
	// }





  private:

	struct item_column
	{
		std::string m_name;
		// column_validator *validator;
	};

	struct item_row
	{

	};

	std::string m_name;
	std::vector<item_column, typename std::allocator_traits<allocator_type>::template rebind_alloc<item_column>> m_columns;
};

using category = category_t<>;

// --------------------------------------------------------------------

template <
	typename Category = category,
	typename Alloc = std::allocator<Category>>
class datablock_t : public std::list<Category, Alloc>
{
  public:
	using category_type = Category;
	using base_type = std::list<category_type, Alloc>;
	using allocator_type = Alloc;

	datablock_t(const std::string &name, const allocator_type &alloc = allocator_type())
		: base_type(alloc)
		, m_name(name)
	{
	}

	datablock_t(const datablock_t &) = default;

	datablock_t(datablock_t &&) = default;

	template <typename Alloc2>
	datablock_t(const datablock_t &db, const Alloc2 &a)
		: base_type(db, a)
		, m_name(db.m_name)
	{
	}

	template <typename Alloc2>
	datablock_t(datablock_t &&db, const Alloc2 &a)
		: base_type(std::move(db), a)
		, m_name(db.m_name)
	{
	}

	datablock_t &operator=(const datablock_t &) = default;
	datablock_t &operator=(datablock_t &&) = default;

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	// --------------------------------------------------------------------

	category_type &operator[](std::string_view name)
	{
		auto i = std::find_if(this->begin(), this->end(), [name](const category_type &c)
			{ return iequals(c.name(), name); });
		if (i == this->end())
			i = this->emplace(name);
		return *i;
	}

	const category_type &operator[](std::string_view name) const
	{
		static const category_type s_empty;
		auto i = std::find_if(this->begin(), this->end(), [name](const category_type &c)
			{ return iequals(c.name(), name); });
		return i == this->end() ? s_empty : *i;
	}

	void write(std::ostream &os) const
	{
		// std::shared_lock lock(mLock);

		os << "data_" << m_name << std::endl
		   << "# " << std::endl;

		// mmcif support, sort of. First write the 'entry' Category
		// and if it exists, _AND_ we have a Validator, write out the
		// audit_conform record.

		for (auto &cat : *this)
		{
			if (cat.name() != "entry")
				continue;

			cat.write(os);

			// if (mValidator != nullptr)
			// {
			// 	Category auditConform(*this, "audit_conform", nullptr);
			// 	auditConform.emplace({{"dict_name", mValidator->dictName()},
			// 		{"dict_version", mValidator->dictVersion()}});
			// 	auditConform.write(os);
			// }

			break;
		}

		for (auto &cat : *this)
		{
			if (cat.name() != "entry" and cat.name() != "audit_conform")
				cat.write(os);
		}
	}

	friend std::ostream &operator<<(std::ostream &os, const datablock_t &db)
	{
		db.write(os);
		return os;
	}

  private:
	std::string m_name;
};

using datablock = datablock_t<>;

// --------------------------------------------------------------------

template <
	typename Datablock = datablock,
	typename Alloc = std::allocator<Datablock>>
class file_t : public std::list<Datablock, Alloc>
{
  public:
	using value_type = Datablock;
	using base_type = std::list<value_type, Alloc>;
	using allocator_type = Alloc;

	file_t() = default;

	file_t(std::istream &is, const allocator_type &alloc = allocator_type())
	{
	}

	file_t(const file_t &) = default;
	file_t(file_t &&) = default;
	file_t &operator=(const file_t &) = default;
	file_t &operator=(file_t &&) = default;
};

using file = file_t<>;

} // namespace cif_v2
