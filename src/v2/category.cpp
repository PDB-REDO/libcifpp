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

#include <cif++/v2/category.hpp>

namespace cif::v2
{

void category::update_value(row *row, size_t column, std::string_view value, bool updateLinked, bool validate)
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

// proxy methods for every insertion
category::iterator category::insert_impl(const_iterator pos, row *n)
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

category::iterator category::erase_impl(const_iterator pos)
{
	if (pos == cend())
		return end();

	assert(false);
	// TODO: implement

	// row *n = const_cast<row *>(pos.row());
	// row *cur;

	// if (m_head == n)
	// {
	// 	m_head = static_cast<row *>(m_head->m_next);
	// 	if (m_head == nullptr)
	// 		m_tail = nullptr;

	// 	n->m_next = nullptr;
	// 	delete_row(n);

	// 	cur = m_head;
	// }
	// else
	// {
	// 	cur = static_cast<row *>(n->m_next);

	// 	if (m_tail == n)
	// 		m_tail = static_cast<row *>(n->m_prev);

	// 	row *p = m_head;
	// 	while (p != nullptr and p->m_next != n)
	// 		p = p->m_next;

	// 	if (p != nullptr and p->m_next == n)
	// 	{
	// 		p->m_next = n->m_next;
	// 		if (p->m_next != nullptr)
	// 			p->m_next->m_prev = p;
	// 		n->m_next = nullptr;
	// 	}
	// 	else
	// 		throw std::runtime_error("remove for a row not found in the list");

	// 	delete_row(n);
	// }

	// return iterator(*this, cur);
}

std::vector<std::string> get_category_fields(const category &cat)
{
	return {};
}

} // namespace cif::v2