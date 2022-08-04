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
#include <cif++/v2/datablock.hpp>

namespace cif::v2
{

template <typename V>
std::string join(const V &arr, std::string_view sep)
{
	std::ostringstream s;

	if (not arr.empty())
	{
		auto ai = arr.begin();
		auto ni = std::next(ai);

		for (;;)
		{
			s << *ai;
			ai = ni;
			ni = std::next(ai);

			if (ni == arr.end())
				break;

			s << sep;
		}
	}

	return s.str();
}

void category::set_validator(const validator *v, datablock &db)
{
	m_validator = v;

	// if (m_index != nullptr)
	// {
	// 	delete m_index;
	// 	m_index = nullptr;
	// }

	if (m_validator != nullptr)
	{
		m_cat_validator = m_validator->get_validator_for_category(m_name);

		// if (m_cat_validator != nullptr)
		// {
		// 	m_index = new CatIndex(this);
		// 	m_index->reconstruct();
		// 	//#if DEBUG
		// 	//			assert(m_index->size() == size());
		// 	//			m_index->validate();
		// 	//#endif
		// }
	}
	else
		m_cat_validator = nullptr;

	for (auto &&[column, cv] : m_columns)
		cv = m_cat_validator ? m_cat_validator->get_validator_for_item(column) : nullptr;

	update_links(db);
}

void category::update_links(datablock &db)
{
	m_child_links.clear();
	m_parent_links.clear();

	if (m_validator != nullptr)
	{
		for (auto link : m_validator->get_links_for_parent(m_name))
		{
			auto childCat = db.get(link->m_child_category);
			if (childCat == nullptr)
				continue;
			m_child_links.emplace_back(childCat, link);
		}

		for (auto link : m_validator->get_links_for_child(m_name))
		{
			auto parentCat = db.get(link->m_parent_category);
			if (parentCat == nullptr)
				continue;
			m_parent_links.emplace_back(parentCat, link);
		}
	}
}

bool category::is_valid() const
{
	bool result = true;

	if (m_validator == nullptr)
		throw std::runtime_error("no Validator specified");

	if (empty())
	{
		if (VERBOSE > 2)
			std::cerr << "Skipping validation of empty Category " << m_name << std::endl;
		return true;
	}

	if (m_cat_validator == nullptr)
	{
		m_validator->report_error("undefined Category " + m_name, false);
		return false;
	}

	auto mandatory = m_cat_validator->m_mandatory_fields;

	for (auto &col : m_columns)
	{
		auto iv = m_cat_validator->get_validator_for_item(col.m_name);
		if (iv == nullptr)
		{
			m_validator->report_error("Field " + col.m_name + " is not valid in Category " + m_name, false);
			result = false;
		}

		// col.m_validator = iv;
		if (col.m_validator != iv)
			m_validator->report_error("Column validator is not specified correctly", true);

		mandatory.erase(col.m_name);
	}

	if (not mandatory.empty())
	{
		m_validator->report_error("In Category " + m_name + " the following mandatory fields are missing: " + join(mandatory, ", "), false);
		result = false;
	}

	//#if not defined(NDEBUG)
	//	// check index?
	//	if (m_index)
	//	{
	//		m_index->validate();
	//		for (auto r: *this)
	//		{
	//			if (m_index->find(r.mData) != r.mData)
	//				m_validator->report_error("Key not found in index for Category " + m_name);
	//		}
	//	}
	//#endif

	// validate all values
	mandatory = m_cat_validator->m_mandatory_fields;

	for (auto ri = m_head; ri != nullptr; ri = ri->m_next)
	{
		for (size_t cix = 0; cix < m_columns.size(); ++cix)
		{
			bool seen = false;
			auto iv = m_columns[cix].m_validator;

			if (iv == nullptr)
			{
				m_validator->report_error("invalid field " + m_columns[cix].m_name + " for Category " + m_name, false);
				result = false;
				continue;
			}

			for (auto vi = ri->m_head; vi != nullptr; vi = vi->m_next)
			{
				if (vi->m_column_ix == cix)
				{
					seen = true;
					try
					{
						(*iv)(vi->text());
					}
					catch (const std::exception &e)
					{
						m_validator->report_error("Error validating " + m_columns[cix].m_name + ": " + e.what(), false);
						continue;
					}
				}
			}

			if (seen or ri != m_head)
				continue;

			if (iv != nullptr and iv->m_mandatory)
			{
				m_validator->report_error("missing mandatory field " + m_columns[cix].m_name + " for Category " + m_name, false);
				result = false;
			}
		}
	}

	return result;
}

category::iterator category::erase(iterator pos)
{
	row_handle rh = *pos;
	row *r = rh;
	iterator result = ++pos;

	iset keys;
	if (m_cat_validator)
		keys = iset(m_cat_validator->m_keys.begin(), m_cat_validator->m_keys.end());

	if (m_head == nullptr)
		throw std::runtime_error("erase");

	// if (mIndex != nullptr)
	// 	mIndex->erase(r.mData);

	if (r == m_head)
	{
		m_head = m_head->m_next;
		r->m_next = nullptr;
	}
	else
	{
		for (auto pi = m_head; pi != nullptr; pi = pi->m_next)
		{
			if (pi->m_next == r)
			{
				pi->m_next = r->m_next;
				r->m_next = nullptr;
				break;
			}
		}
	}

	// links are created based on the _pdbx_item_linked_group_list entries
	// in mmcif_pdbx_v50.dic dictionary.
	//
	// For each link group in _pdbx_item_linked_group_list
	// a std::set of keys from one category is mapped to another.
	// If all values in a child are the same as the specified parent ones
	// the child is removed as well, recursively of course.

	if (m_validator != nullptr)
	{
		for (auto &&[childCat, link] : m_child_links)
		{
			condition cond;

			for (size_t ix = 0; ix < link->m_parent_keys.size(); ++ix)
			{
				std::string_view value = rh[link->m_parent_keys[ix]].text();
				cond = std::move(cond) and (key(link->m_child_keys[ix]) == value);
			}

			childCat->erase_orphans(std::move(cond));
		}
	}

	delete_row(r);

	// reset mTail, if needed
	if (r == m_tail)
	{
		m_tail = m_head;
		if (m_tail != nullptr)
			while (m_tail->m_next != nullptr)
				m_tail = m_tail->m_next;
	}

	return result;	
}

size_t category::erase(condition &&cond)
{
	size_t result = 0;

	cond.prepare(*this);

	auto ri = begin();
	while (ri != end())
	{
		if (cond(*ri))
		{
			ri = erase(ri);
			++result;
		}
		else
			++ri;
	}

	return result;
}

size_t category::erase(condition &&cond, std::function<void(row_handle)> &&visit)
{
	size_t result = 0;

	cond.prepare(*this);

	auto ri = begin();
	while (ri != end())
	{
		if (cond(*ri))
		{
			visit(*ri);
			ri = erase(ri);
			++result;
		}
		else
			++ri;
	}

	return result;
}

bool category::is_orphan(row_handle r) const
{
	// be safe
	if (m_cat_validator == nullptr)
		return false;

	bool isOrphan = true;
	for (auto &&[parentCat, link] : m_parent_links)
	{
		condition cond;
		for (size_t ix = 0; ix < link->m_child_keys.size(); ++ix)
		{
			std::string_view value = r[link->m_child_keys[ix]].text();
			cond = std::move(cond) and (key(link->m_parent_keys[ix]) == value);
		}

		// if (VERBOSE > 2)
		// 	std::cerr << "Check condition '" << cond << "' in parent category " << link->mParentCategory << " for child cat " << mName << std::endl;

		if (parentCat->exists(std::move(cond)))
		{
			if (VERBOSE > 2)
				std::cerr << "Not removing because row has a parent in category " << link->m_parent_category << std::endl;

			isOrphan = false;
			break;
		}
	}

	return isOrphan;
}

void category::erase_orphans(condition &&cond)
{
	std::vector<row *> remove;

	cond.prepare(*this);

	for (auto r : *this)
	{
		if (cond(r) and is_orphan(r))
		{
			if (VERBOSE > 1)
				std::cerr << "Removing orphaned record: " << std::endl
						  << r << std::endl
						  << std::endl;

			remove.push_back(r);
		}
	}

	for (auto r : remove)
		erase(iterator(*this, r));
}


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

	// check the value
	if (col.m_validator and validate)
		col.m_validator->operator()(value);

	// If the field is part of the Key for this Category, remove it from the index
	// before updating

	bool reinsert = false;

	// if (updateLinked and // an update of an Item's value
	// 	cat->m_index != nullptr and cat->keyFieldsByIndex().count(column))
	// {
	// 	reinsert = cat->m_index->find(mData);
	// 	if (reinsert)
	// 		cat->m_index->erase(mData);
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

	// see if we need to update any child categories that depend on this value
	auto iv = col.m_validator;
	if (updateLinked and iv != nullptr /*and m_cascade*/)
	{
		row_handle rh(*this, *row);

		for (auto &&[childCat, linked] : m_child_links)
		{
			if (std::find(linked->m_parent_keys.begin(), linked->m_parent_keys.end(), iv->m_tag) == linked->m_parent_keys.end())
				continue;

			condition cond;
			std::string childTag;

			for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				// TODO add code to *NOT* test mandatory fields for Empty

				if (pk == iv->m_tag)
				{
					childTag = ck;
					cond = std::move(cond) and key(ck) == oldStrValue;
				}
				else
				{
					std::string_view pk_value = rh[pk].text();
					if (pk_value.empty())
						cond = std::move(cond) and key(ck) == null;
					else
						cond = std::move(cond) and ((key(ck) == pk_value) or key(ck) == null);
				}
			}

			auto rows = childCat->find(std::move(cond));
			if (rows.empty())
				continue;

			// if (cif::VERBOSE > 2)
			// {
			// 	std::cerr << "Parent: " << linked->mParentCategory << " Child: " << linked->mChildCategory << std::endl
			// 			  << cond << std::endl;
			// }

			// Now, suppose there are already rows in child that conform to the new value,
			// we then skip this renam

			condition cond_n;

			for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				// TODO add code to *NOT* test mandatory fields for Empty

				if (pk == iv->m_tag)
					cond_n = std::move(cond_n) and key(ck) == value;
				else
				{
					std::string_view pk_value = rh[pk].text();
					if (pk_value.empty())
						cond_n = std::move(cond_n) and key(ck) == null;
					else
						cond_n = std::move(cond_n) and ((key(ck) == pk_value) or key(ck) == null);
				}
			}

			auto rows_n = childCat->find(std::move(cond_n));
			if (not rows_n.empty())
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Will not rename in child category since there are already rows that link to the parent" << std::endl;

				continue;
			}

			for (auto cr : rows)
				cr.assign(childTag, value, false);
		}
	}
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

} // namespace cif::v2