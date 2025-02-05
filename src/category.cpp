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

#include "cif++/category.hpp"
#include "cif++/datablock.hpp"
#include "cif++/parser.hpp"
#include "cif++/utilities.hpp"

#include <numeric>
#include <stack>

// TODO: Find out what the rules are exactly for linked items, the current implementation
// is inconsistent. It all depends whether a link is satified if a item taking part in the
// set of linked items is null at one side and not null in the other.

namespace cif
{

const uint32_t kMaxLineLength = 132;

// --------------------------------------------------------------------

class row_comparator
{
  public:
	row_comparator(category &cat)
	{
		auto cv = cat.get_cat_validator();

		for (auto &k : cv->m_keys)
		{
			uint16_t ix = cat.add_item(k);

			auto iv = cv->get_validator_for_item(k);
			if (iv == nullptr)
				throw std::runtime_error("Incomplete dictionary, no Item Validator for Key " + k);

			auto tv = iv->m_type;
			if (tv == nullptr)
				throw std::runtime_error("Incomplete dictionary, no type Validator for Item " + k);

			using namespace std::placeholders;

			m_comparator.emplace_back(ix, std::bind(&type_validator::compare, tv, _1, _2));
		}
	}

	int operator()(const category &cat, const row *a, const row *b) const
	{
		assert(a);
		assert(b);

		row_handle rha(cat, *a);
		row_handle rhb(cat, *b);

		int d = 0;
		for (const auto &[k, f] : m_comparator)
		{
			std::string_view ka = rha[k].text();
			std::string_view kb = rhb[k].text();

			d = f(ka, kb);

			if (d != 0)
				break;
		}

		return d;
	}

	int operator()(const category &cat, const row_initializer &a, const row *b) const
	{
		assert(b);

		row_handle rhb(cat, *b);

		int d = 0;
		auto ai = a.begin();

		for (const auto &[k, f] : m_comparator)
		{
			assert(ai != a.end());

			std::string_view ka = ai->value();
			std::string_view kb = rhb[k].text();

			d = f(ka, kb);

			if (d != 0)
				break;

			++ai;
		}

		return d;
	}

  private:
	using compareFunc = std::function<int(std::string_view, std::string_view)>;
	using key_comparator = std::tuple<uint16_t, compareFunc>;

	std::vector<key_comparator> m_comparator;
};

// --------------------------------------------------------------------
//
//	class to keep an index on the keys of a category. This is a red/black
//	tree implementation.

class category_index
{
  public:
	category_index(category &cat);

	~category_index()
	{
		delete m_root;
	}

	row *find(const category &cat, row *k) const;
	row *find_by_value(const category &cat, row_initializer k) const;

	void insert(category &cat, row *r);
	void erase(category &cat, row *r);

	// reorder the row's and returns new head and tail
	std::tuple<row *, row *> reorder()
	{
		std::tuple<row *, row *> result = std::make_tuple(nullptr, nullptr);

		if (m_root != nullptr)
		{
			entry *head = find_min(m_root);
			entry *tail = reorder(m_root);

			tail->m_row->m_next = nullptr;

			result = std::make_tuple(head->m_row, tail->m_row);
		}

		return result;
	}

	std::size_t size() const;
	//	bool isValid() const;

  private:
	struct entry
	{
		entry(row *r)
			: m_row(r)
			, m_left(nullptr)
			, m_right(nullptr)
			, m_red(true)
		{
		}

		~entry()
		{
			delete m_left;
			delete m_right;
		}

		row *m_row;
		entry *m_left;
		entry *m_right;
		bool m_red;
	};

	entry *insert(category &cat, entry *h, row *v);
	entry *erase(category &cat, entry *h, row *k);

	//	void validate(entry* h, bool isParentRed, uint32_t blackDepth, uint32_t& minBlack, uint32_t& maxBlack) const;

	entry *rotateLeft(entry *h)
	{
		entry *x = h->m_right;
		h->m_right = x->m_left;
		x->m_left = h;
		x->m_red = h->m_red;
		h->m_red = true;
		return x;
	}

	entry *rotateRight(entry *h)
	{
		entry *x = h->m_left;
		h->m_left = x->m_right;
		x->m_right = h;
		x->m_red = h->m_red;
		h->m_red = true;
		return x;
	}

	void flipColour(entry *h)
	{
		h->m_red = not h->m_red;

		if (h->m_left != nullptr)
			h->m_left->m_red = not h->m_left->m_red;

		if (h->m_right != nullptr)
			h->m_right->m_red = not h->m_right->m_red;
	}

	constexpr bool is_red(entry *h) const
	{
		return h != nullptr and h->m_red;
	}

	entry *move_red_left(entry *h)
	{
		flipColour(h);

		if (h->m_right != nullptr and is_red(h->m_right->m_left))
		{
			h->m_right = rotateRight(h->m_right);
			h = rotateLeft(h);
			flipColour(h);
		}

		return h;
	}

	entry *move_red_right(entry *h)
	{
		flipColour(h);

		if (h->m_left != nullptr and is_red(h->m_left->m_left))
		{
			h = rotateRight(h);
			flipColour(h);
		}

		return h;
	}

	entry *fix_up(entry *h)
	{
		if (is_red(h->m_right))
			h = rotateLeft(h);

		if (is_red(h->m_left) and is_red(h->m_left->m_left))
			h = rotateRight(h);

		if (is_red(h->m_left) and is_red(h->m_right))
			flipColour(h);

		return h;
	}

	entry *find_min(entry *h)
	{
		while (h->m_left != nullptr)
			h = h->m_left;

		return h;
	}

	entry *erase_min(entry *h)
	{
		if (h->m_left == nullptr)
		{
			delete h;
			h = nullptr;
		}
		else
		{
			if (not is_red(h->m_left) and not is_red(h->m_left->m_left))
				h = move_red_left(h);

			h->m_left = erase_min(h->m_left);

			h = fix_up(h);
		}

		return h;
	}

	// Fix m_next items for rows in order of this index
	entry *reorder(entry *e)
	{
		auto result = e;

		if (e->m_left != nullptr)
		{
			auto l = reorder(e->m_left);
			l->m_row->m_next = e->m_row;
		}

		if (e->m_right != nullptr)
		{
			auto mr = find_min(e->m_right);
			e->m_row->m_next = mr->m_row;

			result = reorder(e->m_right);
		}

		return result;
	}

	row_comparator m_row_comparator;
	entry *m_root;
};

category_index::category_index(category &cat)
	: m_row_comparator(cat)
	, m_root(nullptr)
{
	for (auto r : cat)
		insert(cat, r.get_row());
}

row *category_index::find(const category &cat, row *k) const
{
	const entry *r = m_root;
	while (r != nullptr)
	{
		int d = m_row_comparator(cat, k, r->m_row);
		if (d < 0)
			r = r->m_left;
		else if (d > 0)
			r = r->m_right;
		else
			break;
	}

	return r ? r->m_row : nullptr;
}

row *category_index::find_by_value(const category &cat, row_initializer k) const
{
	// sort the values in k first

	row_initializer k2;
	for (auto &f : cat.key_item_indices())
	{
		auto fld = cat.get_item_name(f);

		auto ki = find_if(k.begin(), k.end(), [&fld](auto &i)
			{ return i.name() == fld; });
		if (ki == k.end())
			k2.emplace_back(fld, "");
		else
			k2.emplace_back(*ki);
	}

	const entry *r = m_root;
	while (r != nullptr)
	{
		int d = m_row_comparator(cat, k2, r->m_row);
		if (d < 0)
			r = r->m_left;
		else if (d > 0)
			r = r->m_right;
		else
			break;
	}

	return r ? r->m_row : nullptr;
}

void category_index::insert(category &cat, row *k)
{
	m_root = insert(cat, m_root, k);
	m_root->m_red = false;
}

category_index::entry *category_index::insert(category &cat, entry *h, row *v)
{
	if (h == nullptr)
		return new entry(v);

	int d = m_row_comparator(cat, v, h->m_row);
	if (d < 0)
		h->m_left = insert(cat, h->m_left, v);
	else if (d > 0)
		h->m_right = insert(cat, h->m_right, v);
	else
	{
		row_handle rh(cat, *v);

		std::ostringstream os;
		for (auto col : cat.key_items())
		{
			if (rh[col])
				os << col << ": " << std::quoted(rh[col].text()) << "; ";
		}

		throw duplicate_key_error("Duplicate Key violation, cat: " + cat.name() + " values: " + os.str());
	}

	if (is_red(h->m_right) and not is_red(h->m_left))
		h = rotateLeft(h);

	if (is_red(h->m_left) and is_red(h->m_left->m_left))
		h = rotateRight(h);

	if (is_red(h->m_left) and is_red(h->m_right))
		flipColour(h);

	return h;
}

void category_index::erase(category &cat, row *k)
{
	assert(find(cat, k) == k);

	m_root = erase(cat, m_root, k);
	if (m_root != nullptr)
		m_root->m_red = false;
}

category_index::entry *category_index::erase(category &cat, entry *h, row *k)
{
	if (m_row_comparator(cat, k, h->m_row) < 0)
	{
		if (h->m_left != nullptr)
		{
			if (not is_red(h->m_left) and not is_red(h->m_left->m_left))
				h = move_red_left(h);

			h->m_left = erase(cat, h->m_left, k);
		}
	}
	else
	{
		if (is_red(h->m_left))
			h = rotateRight(h);

		if (m_row_comparator(cat, k, h->m_row) == 0 and h->m_right == nullptr)
		{
			delete h;
			return nullptr;
		}

		if (h->m_right != nullptr)
		{
			if (not is_red(h->m_right) and not is_red(h->m_right->m_left))
				h = move_red_right(h);

			if (m_row_comparator(cat, k, h->m_row) == 0)
			{
				h->m_row = find_min(h->m_right)->m_row;
				h->m_right = erase_min(h->m_right);
			}
			else
				h->m_right = erase(cat, h->m_right, k);
		}
	}

	return fix_up(h);
}

std::size_t category_index::size() const
{
	std::stack<entry *> s;
	s.push(m_root);

	std::size_t result = 0;

	while (not s.empty())
	{
		entry *e = s.top();
		s.pop();

		if (e == nullptr)
			continue;

		++result;

		s.push(e->m_left);
		s.push(e->m_right);
	}

	return result;
}

// --------------------------------------------------------------------

category::category(std::string_view name)
	: m_name(name)
{
}

category::category(const category &rhs)
	: m_name(rhs.m_name)
	, m_items(rhs.m_items)
	, m_cascade(rhs.m_cascade)
{
	for (auto r = rhs.m_head; r != nullptr; r = r->m_next)
		insert_impl(end(), clone_row(*r));

	m_validator = rhs.m_validator;
	m_cat_validator = rhs.m_cat_validator;

	if (m_cat_validator != nullptr and m_index == nullptr)
		m_index = new category_index(*this);
}

void swap(category &a, category &b) noexcept
{
	std::swap(a.m_name, b.m_name);
	std::swap(a.m_items, b.m_items);
	std::swap(a.m_validator, b.m_validator);
	std::swap(a.m_cat_validator, b.m_cat_validator);
	std::swap(a.m_parent_links, b.m_parent_links);
	std::swap(a.m_child_links, b.m_child_links);
	std::swap(a.m_cascade, b.m_cascade);
	std::swap(a.m_index, b.m_index);
	std::swap(a.m_head, b.m_head);
	std::swap(a.m_tail, b.m_tail);
}

category::~category()
{
	clear();
}

// --------------------------------------------------------------------

void category::remove_item(std::string_view item_name)
{
	for (std::size_t ix = 0; ix < m_items.size(); ++ix)
	{
		if (not iequals(item_name, m_items[ix].m_name))
			continue;

		for (row *r = m_head; r != nullptr; r = r->m_next)
		{
			if (r->size() > ix)
				r->erase(r->begin() + ix);
		}

		m_items.erase(m_items.begin() + ix);

		break;
	}
}

void category::rename_item(std::string_view from_name, std::string_view to_name)
{
	for (std::size_t ix = 0; ix < m_items.size(); ++ix)
	{
		if (not iequals(from_name, m_items[ix].m_name))
			continue;

		m_items[ix].m_name = to_name;
		m_items[ix].m_validator = m_cat_validator ? m_cat_validator->get_validator_for_item(to_name) : nullptr;

		break;
	}
}

iset category::get_items() const
{
	iset result;

	for (auto &col : m_items)
		result.insert(col.m_name);

	return result;
}

iset category::key_items() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (m_cat_validator == nullptr)
		throw validation_exception(validation_error::undefined_category);

	iset result;
	for (auto &iv : m_cat_validator->m_item_validators)
		result.insert(iv.m_item_name);

	return result;
}

std::set<uint16_t> category::key_item_indices() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (m_cat_validator == nullptr)
		throw validation_exception(validation_error::undefined_category);

	std::set<uint16_t> result;
	for (auto &k : m_cat_validator->m_keys)
		result.insert(get_item_ix(k));

	return result;
}

// --------------------------------------------------------------------

void category::set_validator(const validator *v, datablock &db)
{
	m_validator = v;

	if (m_index != nullptr)
	{
		delete m_index;
		m_index = nullptr;
	}

	if (m_validator != nullptr)
	{
		m_cat_validator = m_validator->get_validator_for_category(m_name);

		if (m_cat_validator != nullptr)
		{
			std::set<std::string> missing;

			if (not empty())
			{
				std::vector<uint16_t> kix;
				for (auto k : m_cat_validator->m_keys)
				{
					kix.push_back(get_item_ix(k));
					if (kix.back() >= m_items.size())
						missing.insert(k);
				}
			}

			if (missing.empty())
				m_index = new category_index(*this);
			else
			{
				std::ostringstream msg;
				msg << "Cannot construct index since the key item" << (missing.size() > 1 ? "s" : "") << " "
					<< cif::join(missing, ", ") << " in " << m_name << " " << (missing.size() == 1 ? "is" : "are") << " missing\n";
				throw missing_key_error(msg.str(), *missing.begin());
			}
		}
	}
	else
		m_cat_validator = nullptr;

	for (auto &&[item, cv] : m_items)
		cv = m_cat_validator ? m_cat_validator->get_validator_for_item(item) : nullptr;

	update_links(db);
}

void category::update_links(const datablock &db)
{
	m_child_links.clear();
	m_parent_links.clear();

	if (m_validator != nullptr)
	{
		for (auto link : m_validator->get_links_for_parent(m_name))
		{
			auto childCat = const_cast<category *>(db.get(link->m_child_category));
			if (childCat == nullptr)
				continue;
			m_child_links.emplace_back(childCat, link);
		}

		for (auto link : m_validator->get_links_for_child(m_name))
		{
			auto parentCat = const_cast<category *>(db.get(link->m_parent_category));
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
			std::cerr << "Skipping validation of empty category " << m_name << '\n';
		return true;
	}

	if (m_cat_validator == nullptr)
	{
		m_validator->report_error(validation_error::undefined_category, m_name, {}, false);
		return false;
	}

	auto mandatory = m_cat_validator->m_mandatory_items;

	for (auto &col : m_items)
	{
		auto iv = m_cat_validator->get_validator_for_item(col.m_name);
		if (iv == nullptr)
		{
			m_validator->report_error(validation_error::unknown_item, col.m_name, m_name, false);
			result = false;
		}

		// col.m_validator = iv;
		if (col.m_validator != iv)
			m_validator->report_error(validation_error::incorrect_item_validator, true);

		mandatory.erase(col.m_name);
	}

	if (not mandatory.empty())
	{
		m_validator->report_error(validation_error::missing_mandatory_items, m_name, join(mandatory, ", "), false);
		result = false;
	}

	if (m_cat_validator->m_keys.empty() == false and m_index == nullptr)
	{
		std::set<std::string> missing;

		for (auto k : m_cat_validator->m_keys)
		{
			if (get_item_ix(k) >= m_items.size())
				missing.insert(k);
		}

		m_validator->report_error(validation_error::missing_key_items, m_name, join(missing, ", "), false);
		result = false;
	}

	// #if not defined(NDEBUG)
	// 	// check index?
	// 	if (m_index)
	// 	{
	// 		if (m_index->size() != size())
	// 			m_validator->report_error("size of index is not equal to size of category " + m_name, true);

	// 		// m_index->validate();
	// 		for (auto r : *this)
	// 		{
	// 			auto p = r.get_row();
	// 			if (m_index->find(*this, p) != p)
	// 				m_validator->report_error("Key not found in index for category " + m_name, true);
	// 		}
	// 	}
	// #endif

	// validate all values
	mandatory = m_cat_validator->m_mandatory_items;

	for (auto ri = m_head; ri != nullptr; ri = ri->m_next)
	{
		for (uint16_t cix = 0; cix < m_items.size(); ++cix)
		{
			bool seen = false;
			auto iv = m_items[cix].m_validator;

			if (iv == nullptr)
			{
				// no need to report, should have been reported already above
				result = false;
				continue;
			}

			auto vi = ri->get(cix);
			if (vi != nullptr)
			{
				seen = true;
				std::error_code ec;

				iv->validate_value(vi->text(), ec);

				if ((bool)ec)
				{
					m_validator->report_error(ec, m_name, m_items[cix].m_name, false);
					continue;
				}
			}

			if (seen or ri != m_head)
				continue;

			if (iv != nullptr and iv->m_mandatory)
			{
				m_validator->report_error(validation_error::missing_mandatory_items, m_name, m_items[cix].m_name, false);
				result = false;
			}
		}
	}

	return result;
}

bool category::validate_links() const
{
	if (not m_validator)
		return false;

	bool result = true;

	for (auto &link : m_parent_links)
	{
		auto parent = link.linked;

		if (parent == nullptr)
			continue;

		// this particular case should be skipped, that's because it is wrong:
		// there are atoms that are not part of a polymer, and thus will have no
		// parent in that category.
		if (name() == "atom_site" and (parent->name() == "pdbx_poly_seq_scheme" or parent->name() == "entity_poly_seq"))
			continue;

		std::size_t missing = 0;
		category first_missing_rows(name());

		for (auto r : *this)
		{
			auto cond = get_parents_condition(r, *parent);
			if (not cond)
				continue;
			if (not parent->contains(std::move(cond)))
			{
				++missing;
				if (VERBOSE and first_missing_rows.size() < 5)
					first_missing_rows.emplace(r);
			}
		}

		if (missing)
		{
			result = false;

			std::cerr << "Links for " << link.v->m_link_group_label << " are incomplete\n"
					  << "  There are " << missing << " items in " << m_name << " that don't have matching parent items in " << parent->m_name << '\n';

			if (VERBOSE)
			{
				std::cerr << "showing first " << first_missing_rows.size() << " rows\n"
						  << '\n';

				first_missing_rows.write(std::cerr, link.v->m_child_keys, false);

				std::cerr << '\n';
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------

row_handle category::operator[](const key_type &key)
{
	row_handle result{};

	if (not empty())
	{
		if (m_index == nullptr)
			throw std::logic_error("Category " + m_name + " does not have an index");

		auto row = m_index->find_by_value(*this, key);
		if (row != nullptr)
			result = { *this, *row };
	}

	return result;
}

// --------------------------------------------------------------------

condition category::get_parents_condition(row_handle rh, const category &parentCat) const
{
	if (m_validator == nullptr or m_cat_validator == nullptr)
		throw std::runtime_error("No validator known for category " + m_name);

	condition result;

	auto links = m_validator->get_links_for_child(m_name);
	links.erase(remove_if(links.begin(), links.end(), [n = parentCat.m_name](auto &l)
					{ return l->m_parent_category != n; }),
		links.end());

	if (not links.empty())
	{
		for (auto &link : links)
		{
			condition cond;

			for (std::size_t ix = 0; ix < link->m_child_keys.size(); ++ix)
			{
				auto childValue = rh[link->m_child_keys[ix]];

				if (childValue.empty())
					continue;

				cond = std::move(cond) and key(link->m_parent_keys[ix]) == childValue.text();
			}

			result = std::move(result) or std::move(cond);
		}
	}
	else if (cif::VERBOSE > 0)
		std::cerr << "warning: no child to parent links were found for child " << parentCat.name() << " and parent " << name() << '\n';

	return result;
}

condition category::get_children_condition(row_handle rh, const category &childCat) const
{
	if (m_validator == nullptr or m_cat_validator == nullptr)
		throw std::runtime_error("No validator known for category " + m_name);

	condition result;

	iset mandatoryChildItems;
	auto childCatValidator = m_validator->get_validator_for_category(childCat.name());
	if (childCatValidator != nullptr)
		mandatoryChildItems = childCatValidator->m_mandatory_items;

	auto links = m_validator->get_links_for_parent(m_name);
	links.erase(remove_if(links.begin(), links.end(), [n = childCat.m_name](auto &l)
					{ return l->m_child_category != n; }),
		links.end());

	if (not links.empty())
	{
		for (auto &link : links)
		{
			condition cond;

			for (std::size_t ix = 0; ix < link->m_parent_keys.size(); ++ix)
			{
				auto childKey = link->m_child_keys[ix];
				auto parentKey = link->m_parent_keys[ix];

				auto parentValue = rh[parentKey];

				if (parentValue.empty())
					cond = std::move(cond) and key(childKey) == null;
				else if (link->m_parent_keys.size() > 1 and not mandatoryChildItems.contains(childKey))
					cond = std::move(cond) and (key(childKey) == parentValue.text() or key(childKey) == null);
				else
					cond = std::move(cond) and key(childKey) == parentValue.text();
			}

			result = std::move(result) or std::move(cond);
		}
	}
	else if (cif::VERBOSE > 0)
		std::cerr << "warning: no parent to child links were found for parent " << name() << " and child " << childCat.name() << '\n';

	return result;
}

bool category::has_children(row_handle r) const
{
	bool result = false;

	for (auto &&[childCat, link] : m_child_links)
	{
		if (not childCat->contains(get_children_condition(r, *childCat)))
			continue;

		result = true;
		break;
	}

	return result;
}

bool category::has_parents(row_handle r) const
{
	bool result = false;

	for (auto &&[parentCat, link] : m_parent_links)
	{
		if (not parentCat->contains(get_parents_condition(r, *parentCat)))
			continue;

		result = true;
		break;
	}

	return result;
}

std::vector<row_handle> category::get_children(row_handle r, const category &childCat) const
{
	if (m_validator == nullptr or m_cat_validator == nullptr)
		throw std::runtime_error("No validator known for category " + m_name);

	std::vector<row_handle> result;

	for (auto child : childCat.find(get_children_condition(r, childCat)))
	{
		if (std::find(result.begin(), result.end(), child) == result.end())
			result.push_back(child);
	}

	return result;
}

std::vector<row_handle> category::get_parents(row_handle r, const category &parentCat) const
{
	assert(m_validator != nullptr);
	assert(m_cat_validator != nullptr);

	std::vector<row_handle> result;

	for (auto parent : parentCat.find(get_parents_condition(r, parentCat)))
	{
		if (std::find(result.begin(), result.end(), parent) == result.end())
			result.push_back(parent);
	}

	return result;
}

std::vector<row_handle> category::get_linked(row_handle r, const category &cat) const
{
	std::vector<row_handle> result = get_children(r, cat);
	if (result.empty())
		result = get_parents(r, cat);
	return result;
}

// --------------------------------------------------------------------

category::iterator category::erase(iterator pos)
{
	row_handle rh = *pos;
	row *r = rh.get_row();
	iterator result = ++pos;

	if (m_head == nullptr)
		throw std::runtime_error("erase");

	if (m_index != nullptr)
		m_index->erase(*this, r);

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
	// in mmcif_pdbx.dic dictionary.
	//
	// For each link group in _pdbx_item_linked_group_list
	// a set of keys from one category is mapped to another.
	// If all values in a child are the same as the specified parent ones
	// the child is removed as well, recursively of course.

	if (m_validator != nullptr)
	{
		for (auto &&[childCat, link] : m_child_links)
			childCat->erase_orphans(get_children_condition(rh, *childCat), *this);
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

template <typename T>
class save_value
{
  public:
	save_value(T &v, const T nv = {})
		: m_v(v)
		, m_sv(std::exchange(m_v, nv))
	{
	}

	~save_value()
	{
		m_v = m_sv;
	}

  private:
	T &m_v;
	const T m_sv;
};

std::size_t category::erase(condition &&cond)
{
	return erase(std::move(cond), {});
}

std::size_t category::erase(condition &&cond, std::function<void(row_handle)> &&visit)
{
	std::size_t result = 0;

	cond.prepare(*this);

	std::map<category *, condition> potential_orphans;

	auto ri = begin();
	while (ri != end())
	{
		if (cond(*ri))
		{
			if (visit)
				visit(*ri);

			for (auto &&[childCat, link] : m_child_links)
			{
				auto ccond = get_children_condition(*ri, *childCat);
				if (not ccond)
					continue;
				potential_orphans[childCat] = std::move(potential_orphans[childCat]) or std::move(ccond);
			}

			save_value sv(m_validator);

			ri = erase(ri);
			++result;
		}
		else
			++ri;
	}

	for (auto &&[childCat, condition] : potential_orphans)
		childCat->erase_orphans(std::move(condition), *this);

	return result;
}

void category::clear()
{
	auto i = m_head;
	while (i != nullptr)
	{
		auto t = i;
		i = i->m_next;
		delete_row(t);
	}

	m_head = m_tail = nullptr;

	delete m_index;
	m_index = nullptr;
}

void category::erase_orphans(condition &&cond, category &parent)
{
	std::vector<row *> remove;

	cond.prepare(*this);

	for (auto r : *this)
	{
		if (not cond(r))
			continue;

		if (parent.contains(get_parents_condition(r, parent)))
			continue;

		if (VERBOSE > 1)
		{
			category c(m_name);
			c.emplace(r);
			std::cerr << "Removing orphaned record: \n"
					  << c << '\n'
					  << '\n';
		}

		remove.emplace_back(r.m_row);
	}

	for (auto r : remove)
		erase(iterator(*this, r));
}

std::string category::get_unique_id(std::function<std::string(int)> generator)
{
	using namespace cif::literals;

	// calling size() often is a waste of resources
	if (m_last_unique_num == 0)
		m_last_unique_num = static_cast<uint32_t>(size());

	std::string result = generator(static_cast<int>(m_last_unique_num++));

	std::string id_name = "id";
	if (m_cat_validator != nullptr and m_cat_validator->m_keys.size() == 1)
	{
		id_name = m_cat_validator->m_keys.front();

		if (m_index == nullptr and m_cat_validator != nullptr)
			m_index = new category_index(*this);

		for (;;)
		{
			if (m_index->find_by_value(*this, { { id_name, result } }) == nullptr)
				break;
			result = generator(static_cast<int>(m_last_unique_num++));
		}
	}
	else
	{
		for (;;)
		{
			if (not contains(key(id_name) == result))
				break;

			result = generator(static_cast<int>(m_last_unique_num++));
		}
	}

	return result;
}

std::string category::get_unique_value(std::string_view item_name)
{
	std::string result;

	if (m_validator and m_cat_validator)
	{
		auto iv = m_cat_validator->get_validator_for_item(item_name);

		if (iv and iv->m_type and iv->m_type->m_primitive_type == DDL_PrimitiveType::Numb)
		{
			uint64_t v = find_max<uint64_t>(item_name);
			result = std::to_string(v + 1);
		}
	}

	if (result.empty())
	{
		// brain-dead implementation
		for (std::size_t ix = 0; ix < size(); ++ix)
		{
			// result = m_name + "-" + std::to_string(ix);
			result = cif_id_for_number(ix);
			if (not contains(key(item_name) == result))
				break;
		}
	}

	return result;
}

void category::update_value(const std::vector<row_handle> &rows, std::string_view item_name,
	value_provider_type &&value_provider)
{
	using namespace std::literals;

	if (rows.empty())
		return;

	auto colIx = get_item_ix(item_name);
	if (colIx >= m_items.size())
		throw validation_exception(validation_error::unknown_item, m_name, item_name);

	auto &col = m_items[colIx];

	// this is expensive, but better throw early on
	// check the value
	if (col.m_validator)
	{
		for (auto row : rows)
		{
			std::string value{ value_provider(row[item_name].text()) };

			std::error_code ec;
			col.m_validator->validate_value(value, ec);
			if (ec)
				throw validation_exception(ec, m_name, item_name);
		}
	}

	// update and see if we need to update any child categories that depend on this value
	for (auto parent : rows)
	{
		std::string oldValue{ parent[item_name].text() };
		std::string value{ value_provider(oldValue) };

		update_value(parent.get_row(), colIx, value, false, false);

		for (auto &&[childCat, linked] : m_child_links)
		{
			if (std::find(linked->m_parent_keys.begin(), linked->m_parent_keys.end(), item_name) == linked->m_parent_keys.end())
				continue;

			condition cond;
			std::string childItemName;

			for (std::size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				if (pk == item_name)
				{
					childItemName = ck;
					cond = std::move(cond) && key(ck) == oldValue;
				}
				else
					cond = std::move(cond) && key(ck) == parent[pk].text();
			}

			auto children = childCat->find(std::move(cond));
			if (children.empty())
				continue;

			std::vector<row_handle> child_rows;
			std::copy(children.begin(), children.end(), std::back_inserter(child_rows));

			// now be careful. If we search back from child to parent and still find a valid parent row
			// we cannot simply rename the child but will have to create a new child. Unless that new
			// child already exists of course.

			std::vector<row_handle> process;

			for (auto child : child_rows)
			{
				condition cond_c;

				for (std::size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
				{
					std::string pk = linked->m_parent_keys[ix];
					std::string ck = linked->m_child_keys[ix];

					cond_c = std::move(cond_c) && key(pk) == child[ck].text();
				}

				auto parents = find(std::move(cond_c));
				if (parents.empty())
				{
					process.push_back(child);
					continue;
				}

				// oops, we need to split this child, unless a row already exists for the new value
				condition check;

				for (std::size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
				{
					std::string pk = linked->m_parent_keys[ix];
					std::string ck = linked->m_child_keys[ix];

					if (pk == item_name)
						check = std::move(check) && key(ck) == value;
					else
						check = std::move(check) && key(ck) == parent[pk].text();
				}

				if (childCat->contains(std::move(check))) // phew..., narrow escape
					continue;

				// create the actual copy, if we can...
				if (childCat->m_cat_validator != nullptr and childCat->m_cat_validator->m_keys.size() == 1)
				{
					auto copy = childCat->create_copy(child);
					if (copy != child)
					{
						process.push_back(child);
						continue;
					}
				}

				// cannot update this...
				if (cif::VERBOSE > 0)
					std::cerr << "Cannot update child " << childCat->m_name << "." << childItemName << " with value " << value << '\n';
			}

			// finally, update the children
			if (not process.empty())
				childCat->update_value(process, childItemName, value);
		}
	}
}

void category::update_value(row *row, uint16_t item, std::string_view value, bool updateLinked, bool validate)
{
	// make sure we have an index, if possible
	if ((updateLinked or validate) and m_index == nullptr and m_cat_validator != nullptr)
		m_index = new category_index(*this);

	auto &col = m_items[item];

	std::string_view oldValue;

	auto ival = row->get(item);
	if (ival != nullptr)
		oldValue = ival->text();

	if (value == oldValue) // no need to update
		return;

	std::string oldStrValue{ oldValue };

	// check the value
	if (col.m_validator and validate)
		col.m_validator->operator()(value);

	// If the item is part of the Key for this category, remove it from the index
	// before updating

	bool reinsert = false;
	if (m_index != nullptr and key_item_indices().count(item))
	{
		reinsert = m_index->find(*this, row);
		if (reinsert)
			m_index->erase(*this, row);
	}

	// first remove old value with cix
	if (ival != nullptr)
		row->remove(item);

	if (not value.empty())
		row->append(item, { value });

	if (reinsert and m_index != nullptr)
		m_index->insert(*this, row);

	// see if we need to update any child categories that depend on this value
	auto iv = col.m_validator;
	if (updateLinked and iv != nullptr /*and m_cascade*/)
	{
		row_handle rh(*this, *row);

		for (auto &&[childCat, linked] : m_child_links)
		{
			if (std::find(linked->m_parent_keys.begin(), linked->m_parent_keys.end(), iv->m_item_name) == linked->m_parent_keys.end())
				continue;

			condition cond;
			std::string childItemName;

			for (std::size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				// TODO: add code to *NOT* test mandatory items for Empty

				if (pk == iv->m_item_name)
				{
					childItemName = ck;
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
			// 	std::cerr << "Parent: " << linked->mParentcategory << " Child: " << linked->m_child_category << '\n'
			// 			  << cond << '\n';
			// }

			// Now, suppose there are already rows in child that conform to the new value,
			// we then skip this rename

			condition cond_n;

			for (std::size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				if (pk == iv->m_item_name)
					cond_n = std::move(cond_n) and key(ck) == value;
				else
				{
					std::string_view pk_value = rh[pk].text();
					if (pk_value.empty())
						cond_n = std::move(cond_n) and key(ck) == null;
					else
						cond_n = std::move(cond_n) and key(ck) == pk_value;
				}
			}

			auto rows_n = childCat->find(std::move(cond_n));
			if (not rows_n.empty())
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Will not rename in child category since there are already rows that link to the parent\n";

				continue;
			}

			for (auto cr : rows)
				cr.assign(childItemName, value, false);
		}
	}
}

row *category::clone_row(const row &r)
{
	row *result = create_row();

	try
	{
		for (uint16_t ix = 0; ix < r.size(); ++ix)
		{
			auto &i = r[ix];
			if (not i)
				continue;

			result->append(ix, { i.text() });
		}
	}
	catch (...)
	{
		delete_row(result);
		throw;
	}

	return result;
}

void category::delete_row(row *r)
{
	if (r != nullptr)
	{
		row_allocator_type ra(get_allocator());
		row_allocator_traits::destroy(ra, r);
		row_allocator_traits::deallocate(ra, r, 1);
	}
}

row_handle category::create_copy(row_handle r)
{
	// copy the values
	std::vector<item> items;

	for (uint16_t ix = 0; ix < r.m_row->size(); ++ix)
	{
		auto i = r.m_row->get(ix);
		if (i != nullptr)
			items.emplace_back(m_items[ix].m_name, i->text());
	}

	if (m_cat_validator and m_cat_validator->m_keys.size() == 1)
	{
		auto key = m_cat_validator->m_keys.front();
		auto kv = m_cat_validator->get_validator_for_item(key);

		for (auto &item : items)
		{
			if (item.name() != key)
				continue;

			if (kv->m_type->m_primitive_type == DDL_PrimitiveType::Numb)
				item.value(get_unique_id(""));
			else
				item.value(get_unique_id(m_name + "_id_"));
			break;
		}
	}

	return emplace(items.begin(), items.end());
}

// proxy methods for every insertion
category::iterator category::insert_impl(const_iterator pos, row *n)
{
	if (m_index == nullptr and m_cat_validator != nullptr)
		m_index = new category_index(*this);

	assert(n != nullptr);
	assert(n->m_next == nullptr);

	if (n == nullptr)
		throw std::runtime_error("Invalid pointer passed to insert");

	// #ifndef NDEBUG
	// 	if (m_validator)
	// 		is_valid();
	// #endif

	try
	{
		// First, make sure all mandatory items are supplied
		if (m_cat_validator != nullptr)
		{
			for (uint16_t ix = 0; ix < static_cast<uint16_t>(m_items.size()); ++ix)
			{
				const auto &[item, iv] = m_items[ix];

				if (iv == nullptr)
					continue;

				bool seen = false;

				auto i = n->get(ix);
				if (i != nullptr)
				{
					iv->operator()(i->text());
					seen = true;
				}

				if (not seen and iv->m_mandatory)
					throw std::runtime_error("missing mandatory item " + item + " for category " + m_name);
			}
		}

		if (m_index != nullptr)
			m_index->insert(*this, n);

		// insert at end, most often this is the case
		if (pos.m_current.m_row == nullptr)
		{
			if (m_head == nullptr)
				m_tail = m_head = n;
			else
				m_tail = m_tail->m_next = n;
		}
		else
		{
			assert(m_head != nullptr);

			if (pos.m_current.m_row == m_head)
				m_head = n->m_next = m_head;
			else
				n = n->m_next = m_head->m_next;
		}

		return iterator(*this, n);
	}
	catch (const std::exception &e)
	{
		delete_row(n);
		throw;
	}

	// #ifndef NDEBUG
	// 	if (m_validator)
	// 		is_valid();
	// #endif
}

void category::swap_item(uint16_t item_ix, row_handle &a, row_handle &b)
{
	assert(this == a.m_category);
	assert(this == b.m_category);

	auto &ra = *a.m_row;
	auto &rb = *b.m_row;

	while (ra.size() <= item_ix)
		ra.emplace_back("");

	while (rb.size() <= item_ix)
		rb.emplace_back("");

	std::swap(ra.at(item_ix), rb.at(item_ix));
}

void category::sort(std::function<int(row_handle, row_handle)> f)
{
	if (m_head == nullptr)
		return;

	std::vector<row_handle> rows;
	for (auto itemRow = m_head; itemRow != nullptr; itemRow = itemRow->m_next)
		rows.emplace_back(*this, *itemRow);

	std::stable_sort(rows.begin(), rows.end(),
		[&f](row_handle ia, row_handle ib)
		{
			return f(ia, ib) < 0;
		});

	m_head = rows.front().get_row();
	m_tail = rows.back().get_row();

	auto r = m_head;
	for (std::size_t i = 1; i < rows.size(); ++i)
		r = r->m_next = rows[i].get_row();
	r->m_next = nullptr;

	assert(r == m_tail);
	assert(size() == rows.size());
}

void category::reorder_by_index()
{
	if (m_index)
		std::tie(m_head, m_tail) = m_index->reorder();
}

namespace detail
{
	std::size_t write_value(std::ostream &os, std::string_view value, std::size_t offset, std::size_t width, bool right_aligned)
	{
		if (value.find('\n') != std::string::npos or width == 0 or value.length() > 132) // write as text item
		{
			if (offset > 0)
				os << '\n';
			os << ';';

			char pc = 0;
			for (auto ch : value)
			{
				if (pc == '\n' and ch == ';')
					os << '\\';
				os << ch;
				pc = ch;
			}

			if (value.back() != '\n')
				os << '\n';
			os << ';' << '\n';
			offset = 0;
		}
		else if (sac_parser::is_unquoted_string(value))
		{
			if (right_aligned)
			{
				if (value.length() < width)
				{
					os << std::string(width - value.length() - 1, ' ');
					offset += width;
				}
				else
					offset += value.length() + 1;
			}

			os << value;

			if (right_aligned)
				os << ' ';
			else
			{
				if (value.length() < width)
				{
					os << std::string(width - value.length(), ' ');
					offset += width;
				}
				else
				{
					os << ' ';
					offset += value.length() + 1;
				}
			}
		}
		else
		{
			bool done = false;
			for (char q : { '\'', '"' })
			{
				auto p = value.find(q); // see if we can use the quote character
				while (p != std::string::npos and sac_parser::is_non_blank(value[p + 1]) and value[p + 1] != q)
					p = value.find(q, p + 1);

				if (p != std::string::npos)
					continue;

				os << q << value << q;

				if (value.length() + 2 < width)
				{
					os << std::string(width - value.length() - 2, ' ');
					offset += width;
				}
				else
				{
					os << ' ';
					offset += value.length() + 1;
				}

				done = true;
				break;
			}

			if (not done)
			{
				if (offset > 0)
					os << '\n';
				os << ';' << value << '\n'
				   << ';' << '\n';
				offset = 0;
			}
		}

		return offset;
	}

} // namespace detail

std::vector<std::string> category::get_item_order() const
{
	std::vector<std::string> result;
	for (auto &c : m_items)
		result.push_back("_" + m_name + "." + c.m_name);
	return result;
}

void category::write(std::ostream &os) const
{
	std::vector<uint16_t> order(m_items.size());
	iota(order.begin(), order.end(), static_cast<uint16_t>(0));
	write(os, order, false);
}

void category::write(std::ostream &os, const std::vector<std::string> &items, bool addMissingItems)
{
	// make sure all items are present
	for (auto &c : items)
		add_item(c);

	std::vector<uint16_t> order;
	order.reserve(m_items.size());

	for (auto &c : items)
		order.push_back(get_item_ix(c));

	if (addMissingItems)
	{
		for (uint16_t i = 0; i < m_items.size(); ++i)
		{
			if (std::find(order.begin(), order.end(), i) == order.end())
				order.push_back(i);
		}
	}

	write(os, order, true);
}

void category::write(std::ostream &os, const std::vector<uint16_t> &order, bool includeEmptyItems) const
{
	if (empty())
		return;

	// If the first Row has a next, we need a loop_
	bool needLoop = (m_head->m_next != nullptr);

	std::vector<bool> right_aligned(m_items.size(), false);

	if (m_cat_validator != nullptr)
	{
		for (auto cix : order)
		{
			auto &col = m_items[cix];
			right_aligned[cix] = col.m_validator != nullptr and
			                     col.m_validator->m_type != nullptr and
			                     col.m_validator->m_type->m_primitive_type == cif::DDL_PrimitiveType::Numb;
		}
	}

	if (needLoop)
	{
		os << "loop_\n";

		std::vector<std::size_t> itemWidths(m_items.size());

		for (auto cix : order)
		{
			auto &col = m_items[cix];
			os << '_';
			if (not m_name.empty())
				os << m_name << '.';
			os << col.m_name << ' ' << '\n';
			itemWidths[cix] = 2;
		}

		for (auto r = m_head; r != nullptr; r = r->m_next)
		{
			for (uint16_t ix = 0; ix < r->size(); ++ix)
			{
				auto v = r->get(ix);
				if (v == nullptr)
					continue;

				if (v->text().find('\n') == std::string_view::npos)
				{
					std::size_t l = v->text().length();

					if (not sac_parser::is_unquoted_string(v->text()))
						l += 2;

					if (l > 132)
						continue;

					if (itemWidths[ix] < l + 1)
						itemWidths[ix] = l + 1;
				}
			}
		}

		for (auto r = m_head; r != nullptr; r = r->m_next) // loop over rows
		{
			std::size_t offset = 0;

			for (uint16_t cix : order)
			{
				std::size_t w = itemWidths[cix];

				std::string_view s;
				auto iv = r->get(cix);
				if (iv != nullptr)
					s = iv->text();

				if (s.empty())
					s = "?";

				std::size_t l = s.length();
				if (not sac_parser::is_unquoted_string(s))
					l += 2;
				if (l < w)
					l = w;

				if (offset + l > 132 and offset > 0)
				{
					os << '\n';
					offset = 0;
				}

				offset = detail::write_value(os, s, offset, w, right_aligned[cix]);

				if (offset > 132)
				{
					os << '\n';
					offset = 0;
				}
			}

			if (offset > 0)
				os << '\n';
		}
	}
	else
	{
		// first find the indent level
		std::size_t l = 0;

		for (auto &col : m_items)
		{
			std::string item_name = '_' + m_name + '.' + col.m_name;

			if (l < item_name.length())
				l = item_name.length();
		}

		l += 3;

		std::size_t width = 1;

		for (auto cix : order)
		{
			if (not right_aligned[cix])
				continue;

			std::string_view s;
			auto iv = m_head->get(cix);
			if (iv != nullptr)
				s = iv->text();

			if (s.empty())
				s = "?";

			std::size_t l2 = s.length();

			if (not sac_parser::is_unquoted_string(s))
				l2 += 2;

			if (width < l2)
				width = l2;
		}

		for (uint16_t cix : order)
		{
			auto &col = m_items[cix];

			os << '_';
			if (not m_name.empty())
				os << m_name << '.';
			os << col.m_name << std::string(l - col.m_name.length() - m_name.length() - 2, ' ');

			std::string_view s;
			auto iv = m_head->get(cix);
			if (iv != nullptr)
				s = iv->text();

			if (s.empty())
				s = "?";

			std::size_t offset = l;
			if (s.length() + l >= kMaxLineLength)
			{
				os << '\n';
				offset = 0;
			}

			if (detail::write_value(os, s, offset, width, s.empty() or right_aligned[cix]) != 0)
				os << '\n';
		}
	}

	os << "# \n";
}

bool category::operator==(const category &rhs) const
{
	// shortcut
	if (this == &rhs)
		return true;

	auto &a = *this;
	auto &b = rhs;

	using namespace std::placeholders;

	//	set<std::string> item_namesA(a.items()), item_namesB(b.items());
	//
	//	if (item_namesA != item_namesB)
	//		std::cout << "Unequal number of items\n";

	const category_validator *catValidator = nullptr;

	auto validator = a.get_validator();
	if (validator != nullptr)
		catValidator = validator->get_validator_for_category(a.name());

	typedef std::function<int(std::string_view, std::string_view)> compType;
	std::vector<std::tuple<std::string, compType>> item_names;
	std::vector<std::string> keys;
	std::vector<std::size_t> keyIx;

	if (catValidator == nullptr)
	{
		for (auto &item_name : a.get_items())
		{
			item_names.push_back(std::make_tuple(item_name, [](std::string_view va, std::string_view vb)
				{ return va.compare(vb); }));
			keyIx.push_back(keys.size());
			keys.push_back(item_name);
		}
	}
	else
	{
		keys = catValidator->m_keys;

		for (auto &item_name : a.key_items())
		{
			auto iv = catValidator->get_validator_for_item(item_name);
			if (iv == nullptr)
				throw std::runtime_error("missing item validator");
			auto tv = iv->m_type;
			if (tv == nullptr)
				throw std::runtime_error("missing type validator");
			item_names.push_back(std::make_tuple(item_name, std::bind(&cif::type_validator::compare, tv, std::placeholders::_1, std::placeholders::_2)));

			auto pred = [item_name](const std::string &s) -> bool
			{
				return cif::iequals(item_name, s) == 0;
			};
			if (find_if(keys.begin(), keys.end(), pred) == keys.end())
				keyIx.push_back(item_names.size() - 1);
		}
	}

	// a.reorderByIndex();
	// b.reorderByIndex();

	auto rowEqual = [&](const row_handle &a, const row_handle &b)
	{
		int d = 0;

		for (auto kix : keyIx)
		{
			std::string item_name;
			compType compare;

			std::tie(item_name, compare) = item_names[kix];

			d = compare(a[item_name].text(), b[item_name].text());

			if (d != 0)
				break;
		}

		return d == 0;
	};

	auto ai = a.begin(), bi = b.begin();
	while (ai != a.end() or bi != b.end())
	{
		if (ai == a.end() or bi == b.end())
			return false;

		auto ra = *ai, rb = *bi;

		if (not rowEqual(ra, rb))
			return false;

		std::vector<std::string> missingA, missingB, different;

		for (auto &tt : item_names)
		{
			std::string item_name;
			compType compare;

			std::tie(item_name, compare) = tt;

			// make it an option to compare unapplicable to empty or something

			auto ta = ra[item_name].text();
			if (ta == "." or ta == "?")
				ta = "";
			auto tb = rb[item_name].text();
			if (tb == "." or tb == "?")
				tb = "";

			if (compare(ta, tb) != 0)
				return false;
		}

		++ai;
		++bi;
	}

	return true;
}

} // namespace cif