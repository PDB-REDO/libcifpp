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

#include <numeric>

#include <cif++/category.hpp>
#include <cif++/datablock.hpp>
#include <cif++/parser.hpp>
#include <cif++/utilities.hpp>

// TODO: Find out what the rules are exactly for linked items, the current implementation
// is inconsistent. It all depends whether a link is satified if a field taking part in the
// set of linked items is null at one side and not null in the other.

namespace cif
{

const uint32_t kMaxLineLength = 132;

// --------------------------------------------------------------------

class row_comparator
{
  public:
	row_comparator(category &cat)
		: m_category(cat)
	{
		auto cv = cat.get_cat_validator();

		for (auto k : cv->m_keys)
		{
			size_t ix = cat.get_column_ix(k);

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

	int operator()(const row *a, const row *b) const
	{
		assert(a);
		assert(b);

		row_handle rha(m_category, *a);
		row_handle rhb(m_category, *b);

		int d = 0;
		for (auto &c : m_comparator)
		{
			size_t k;
			compareFunc f;

			std::tie(k, f) = c;

			std::string_view ka = rha[k].text();
			std::string_view kb = rhb[k].text();

			d = f(ka, kb);

			if (d != 0)
				break;
		}

		return d;
	}

  private:
	typedef std::function<int(std::string_view, std::string_view)> compareFunc;
	typedef std::tuple<size_t, compareFunc> key_comparator;

	std::vector<key_comparator> m_comparator;
	category &m_category;
};

// --------------------------------------------------------------------
//
//	class to keep an index on the keys of a category. This is a red/black
//	tree implementation.

class category_index
{
  public:
	category_index(category *cat)
		: m_category(*cat)
		, m_row_comparator(m_category)
		, m_root(nullptr)
	{
		reconstruct();
	}

	~category_index()
	{
		delete m_root;
	}

	row *find(row *k) const;

	void insert(row *r);
	void erase(row *r);

	// batch create
	void reconstruct();

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

	size_t size() const;
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

	entry *insert(entry *h, row *v);
	entry *erase(entry *h, row *k);

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

	bool is_red(entry *h) const
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

	// Fix m_next fields for rows in order of this index
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

	category &m_category;
	row_comparator m_row_comparator;
	entry *m_root;
};

row *category_index::find(row *k) const
{
	const entry *r = m_root;
	while (r != nullptr)
	{
		int d = m_row_comparator(k, r->m_row);
		if (d < 0)
			r = r->m_left;
		else if (d > 0)
			r = r->m_right;
		else
			break;
	}

	return r ? r->m_row : nullptr;
}

void category_index::insert(row *k)
{
	m_root = insert(m_root, k);
	m_root->m_red = false;
}

category_index::entry *category_index::insert(entry *h, row *v)
{
	if (h == nullptr)
		return new entry(v);

	int d = m_row_comparator(v, h->m_row);
	if (d < 0)
		h->m_left = insert(h->m_left, v);
	else if (d > 0)
		h->m_right = insert(h->m_right, v);
	else
	{
		row_handle rh(m_category, *v);

		std::ostringstream os;
		for (auto col : m_category.fields())
			os << col << ": " << std::quoted(rh[col].text()) << "; ";

		throw std::runtime_error("Duplicate Key violation, cat: " + m_category.name() + " values: " + os.str());
	}

	if (is_red(h->m_right) and not is_red(h->m_left))
		h = rotateLeft(h);

	if (is_red(h->m_left) and is_red(h->m_left->m_left))
		h = rotateRight(h);

	if (is_red(h->m_left) and is_red(h->m_right))
		flipColour(h);

	return h;
}

void category_index::erase(row *k)
{
	assert(find(k) == k);

	m_root = erase(m_root, k);
	if (m_root != nullptr)
		m_root->m_red = false;
}

category_index::entry *category_index::erase(entry *h, row *k)
{
	if (m_row_comparator(k, h->m_row) < 0)
	{
		if (h->m_left != nullptr)
		{
			if (not is_red(h->m_left) and not is_red(h->m_left->m_left))
				h = move_red_left(h);

			h->m_left = erase(h->m_left, k);
		}
	}
	else
	{
		if (is_red(h->m_left))
			h = rotateRight(h);

		if (m_row_comparator(k, h->m_row) == 0 and h->m_right == nullptr)
		{
			delete h;
			return nullptr;
		}

		if (h->m_right != nullptr)
		{
			if (not is_red(h->m_right) and not is_red(h->m_right->m_left))
				h = move_red_right(h);

			if (m_row_comparator(k, h->m_row) == 0)
			{
				h->m_row = find_min(h->m_right)->m_row;
				h->m_right = erase_min(h->m_right);
			}
			else
				h->m_right = erase(h->m_right, k);
		}
	}

	return fix_up(h);
}

void category_index::reconstruct()
{
	delete m_root;
	m_root = nullptr;

	for (auto r : m_category)
		insert(r);

	// maybe reconstruction can be done quicker by using the following commented code.
	// however, I've not had the time to think of a way to set the red/black flag correctly in that case.

	//	std::vector<row*> rows;
	//	transform(mCat.begin(), mCat.end(), backInserter(rows),
	//		[](Row r) -> row* { assert(r.mData); return r.mData; });
	//
	//	assert(std::find(rows.begin(), rows.end(), nullptr) == rows.end());
	//
	//	// don't use sort here, it will run out of the stack of something.
	//	// quicksort is notorious for using excessive recursion.
	//	// Besides, most of the time, the data is ordered already anyway.
	//
	//	stable_sort(rows.begin(), rows.end(), [this](row* a, row* b) -> bool { return this->mComp(a, b) < 0; });
	//
	//	for (size_t i = 0; i < rows.size() - 1; ++i)
	//		assert(mComp(rows[i], rows[i + 1]) < 0);
	//
	//	deque<entry*> e;
	//	transform(rows.begin(), rows.end(), back_inserter(e),
	//		[](row* r) -> entry* { return new entry(r); });
	//
	//	while (e.size() > 1)
	//	{
	//		deque<entry*> ne;
	//
	//		while (not e.empty())
	//		{
	//			entry* a = e.front();
	//			e.pop_front();
	//
	//			if (e.empty())
	//				ne.push_back(a);
	//			else
	//			{
	//				entry* b = e.front();
	//				b->mLeft = a;
	//
	//				assert(mComp(a->mRow, b->mRow) < 0);
	//
	//				e.pop_front();
	//
	//				if (not e.empty())
	//				{
	//					entry* c = e.front();
	//					e.pop_front();
	//
	//					assert(mComp(b->mRow, c->mRow) < 0);
	//
	//					b->mRight = c;
	//				}
	//
	//				ne.push_back(b);
	//
	//				if (not e.empty())
	//				{
	//					ne.push_back(e.front());
	//					e.pop_front();
	//				}
	//			}
	//		}
	//
	//		swap (e, ne);
	//	}
	//
	//	assert(e.size() == 1);
	//	mRoot = e.front();
}

size_t category_index::size() const
{
	std::stack<entry *> s;
	s.push(m_root);

	size_t result = 0;

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
	, m_columns(rhs.m_columns)
	, m_validator(rhs.m_validator)
	, m_cat_validator(rhs.m_cat_validator)
	, m_parent_links(rhs.m_parent_links)
	, m_child_links(rhs.m_child_links)
	, m_cascade(rhs.m_cascade)
{
	for (auto r = rhs.m_head; r != nullptr; r = r->m_next)
		insert_impl(end(), clone_row(*r));

	if (m_cat_validator != nullptr)
		m_index = new category_index(this);
}

category::category(category &&rhs)
	: m_name(std::move(rhs.m_name))
	, m_columns(std::move(rhs.m_columns))
	, m_validator(rhs.m_validator)
	, m_cat_validator(rhs.m_cat_validator)
	, m_parent_links(std::move(rhs.m_parent_links))
	, m_child_links(std::move(rhs.m_child_links))
	, m_cascade(rhs.m_cascade)
	, m_index(rhs.m_index)
	, m_head(rhs.m_head)
	, m_tail(rhs.m_tail)
{
	rhs.m_head = nullptr;
	rhs.m_tail = nullptr;
	rhs.m_index = nullptr;
}

category &category::operator=(const category &rhs)
{
	if (this != &rhs)
	{
		if (not empty())
			clear();

		m_name = rhs.m_name;
		m_columns = rhs.m_columns;
		m_cascade = rhs.m_cascade;

		m_validator = nullptr;
		m_cat_validator = nullptr;

		delete m_index;
		m_index = nullptr;

		for (auto r = rhs.m_head; r != nullptr; r = r->m_next)
			insert_impl(cend(), clone_row(*r));

		m_validator = rhs.m_validator;
		m_cat_validator = rhs.m_cat_validator;
		m_parent_links = rhs.m_parent_links;
		m_child_links = rhs.m_child_links;

		if (m_cat_validator != nullptr)
			m_index = new category_index(this);
	}

	return *this;
}

category &category::operator=(category &&rhs)
{
	if (this != &rhs)
	{
		if (not empty())
			clear();

		m_name = std::move(rhs.m_name);
		m_columns = std::move(rhs.m_columns);
		m_cascade = rhs.m_cascade;
		m_validator = rhs.m_validator;
		m_cat_validator = rhs.m_cat_validator;
		m_parent_links = rhs.m_parent_links;
		m_child_links = rhs.m_child_links;
		m_index = rhs.m_index;
		m_head = rhs.m_head;
		m_tail = rhs.m_tail;

		rhs.m_head = rhs.m_tail = nullptr;
		rhs.m_index = nullptr;
	}

	return *this;
}

category::~category()
{
	clear();
}

// --------------------------------------------------------------------

iset category::fields() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (m_cat_validator == nullptr)
		m_validator->report_error("undefined Category", true);

	iset result;
	for (auto &iv : m_cat_validator->m_item_validators)
		result.insert(iv.m_tag);

	return result;
}

std::set<uint16_t> category::key_field_indices() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (m_cat_validator == nullptr)
		m_validator->report_error("undefined Category", true);

	std::set<uint16_t> result;
	for (auto &k : m_cat_validator->m_keys)
		result.insert(get_column_ix(k));

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
			m_index = new category_index(this);
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
			std::cerr << "Skipping validation of empty category " << m_name << std::endl;
		return true;
	}

	if (m_cat_validator == nullptr)
	{
		m_validator->report_error("undefined category " + m_name, false);
		return false;
	}

	auto mandatory = m_cat_validator->m_mandatory_fields;

	for (auto &col : m_columns)
	{
		auto iv = m_cat_validator->get_validator_for_item(col.m_name);
		if (iv == nullptr)
		{
			m_validator->report_error("Field " + col.m_name + " is not valid in category " + m_name, false);
			result = false;
		}

		// col.m_validator = iv;
		if (col.m_validator != iv)
			m_validator->report_error("Column validator is not specified correctly", true);

		mandatory.erase(col.m_name);
	}

	if (not mandatory.empty())
	{
		m_validator->report_error("In category " + m_name + " the following mandatory fields are missing: " + join(mandatory, ", "), false);
		result = false;
	}

#if not defined(NDEBUG)
	// check index?
	if (m_index)
	{
		// m_index->validate();
		for (auto r : *this)
		{
			if (m_index->find(r) != r)
				m_validator->report_error("Key not found in index for category " + m_name, true);
		}
	}
#endif

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
				m_validator->report_error("invalid field " + m_columns[cix].m_name + " for category " + m_name, false);
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
				m_validator->report_error("missing mandatory field " + m_columns[cix].m_name + " for category " + m_name, false);
				result = false;
			}
		}
	}

	return result;
}

void category::validate_links() const
{
	if (not m_validator)
		return;

	for (auto &link : m_parent_links)
	{
		auto parent = link.linked;

		if (parent == nullptr)
			continue;

		size_t missing = 0;
		category first_missing_rows(name());

		for (auto r : *this)
		{
			auto cond = get_parents_condition(r, *parent);
			if (not cond)
				continue;
			if (not parent->exists(std::move(cond)))
			{
				++missing;
				if (VERBOSE and first_missing_rows.size() < 5)
					first_missing_rows.emplace(r);
			}
		}

		if (missing)
		{
			std::cerr << "Links for " << link.v->m_link_group_label << " are incomplete" << std::endl
					  << "  There are " << missing << " items in " << m_name << " that don't have matching parent items in " << parent->m_name << std::endl;
			
			if (VERBOSE)
			{
				std::cerr << "showing first " << first_missing_rows.size() <<  " rows" << std::endl
						  << std::endl;

				first_missing_rows.write(std::cerr, link.v->m_child_keys, false);

				std::cerr << std::endl;
			}
		}
	}
}

// --------------------------------------------------------------------

condition category::get_parents_condition(row_handle rh, const category &parentCat) const
{
	if (m_validator == nullptr or m_cat_validator == nullptr)
		throw std::runtime_error("No validator known for category " + m_name);

	condition result;

	for (auto &link : m_validator->get_links_for_child(m_name))
	{
		if (link->m_parent_category != parentCat.m_name)
			continue;

		condition cond;

		for (size_t ix = 0; ix < link->m_child_keys.size(); ++ix)
		{
			auto childValue = rh[link->m_child_keys[ix]];

			if (childValue.empty())
				continue;

			cond = std::move(cond) and key(link->m_parent_keys[ix]) == childValue.text();
		}

		if (result)
			result = std::move(result) or std::move(cond);
		else
			result = std::move(cond);
	}

	return result;
}

condition category::get_children_condition(row_handle rh, const category &childCat) const
{
	if (m_validator == nullptr or m_cat_validator == nullptr)
		throw std::runtime_error("No validator known for category " + m_name);

	condition result;

	for (auto &link : m_validator->get_links_for_parent(m_name))
	{
		if (link->m_child_category != childCat.m_name)
			continue;

		condition cond;

		for (size_t ix = 0; ix < link->m_parent_keys.size(); ++ix)
		{
			auto childKey = link->m_child_keys[ix];
			auto parentKey = link->m_parent_keys[ix];

			auto parentValue = rh[parentKey];

			if (parentValue.empty())
				cond = std::move(cond) and key(childKey) == null;
			else
				cond = std::move(cond) and key(childKey) == parentValue.text();
		}

		if (result)
			result = std::move(result) or std::move(cond);
		else
			result = std::move(cond);
	}

	return result;
}

bool category::has_children(row_handle r) const
{
	bool result = false;

	for (auto &&[childCat, link] : m_child_links)
	{
		if (not childCat->exists(get_children_condition(r, *childCat)))
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
		if (not parentCat->exists(get_parents_condition(r, *parentCat)))
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
	row *r = rh;
	iterator result = ++pos;

	iset keys;
	if (m_cat_validator)
		keys = iset(m_cat_validator->m_keys.begin(), m_cat_validator->m_keys.end());

	if (m_head == nullptr)
		throw std::runtime_error("erase");

	if (m_index != nullptr)
		m_index->erase(r);

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
	// a std::set of keys from one category is mapped to another.
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
		
		if (parent.exists(get_parents_condition(r, parent)))
			continue;

		if (VERBOSE > 1)
			std::cerr << "Removing orphaned record: " << std::endl
						<< r << std::endl
						<< std::endl;
		
		remove.emplace_back(r.m_row);
	}

	for (auto r : remove)
		erase(iterator(*this, r));
}

std::string category::get_unique_id(std::function<std::string(int)> generator)
{
	using namespace cif::literals;

	std::string id_tag = "id";
	if (m_cat_validator != nullptr and m_cat_validator->m_keys.size() == 1)
		id_tag = m_cat_validator->m_keys.front();

	// calling size() often is a waste of resources
	if (m_last_unique_num == 0)
		m_last_unique_num = size();

	for (;;)
	{
		std::string result = generator(static_cast<int>(m_last_unique_num++));

		if (exists(key(id_tag) == result))
			continue;

		return result;
	}
}

void category::update_value(const std::vector<row_handle> &rows, std::string_view tag, std::string_view value)
{
	using namespace std::literals;

	if (rows.empty())
		return;

	auto colIx = get_column_ix(tag);
	if (colIx >= m_columns.size())
		throw std::runtime_error("Invalid column " + std::string{ value } + " for " + m_name);

	auto &col = m_columns[colIx];

	// check the value
	if (col.m_validator)
		(*col.m_validator)(value);

	// first some sanity checks, what was the old value and is it the same for all rows?
	std::string oldValue{ rows.front()[tag].text() };
	for (auto row : rows)
	{
		if (oldValue != row[tag].text())
			throw std::runtime_error("Inconsistent old values in update_value");
	}

	if (oldValue == value) // no need to do anything
		return;

	// update rows, but do not cascade
	for (auto row : rows)
		row.assign(colIx, value, false);

	// see if we need to update any child categories that depend on this value
	for (auto parent : rows)
	{
		for (auto &&[childCat, linked] : m_child_links)
		{
			if (std::find(linked->m_parent_keys.begin(), linked->m_parent_keys.end(), tag) == linked->m_parent_keys.end())
				continue;

			condition cond;
			std::string childTag;

			for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				if (pk == tag)
				{
					childTag = ck;
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

				for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
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

				for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
				{
					std::string pk = linked->m_parent_keys[ix];
					std::string ck = linked->m_child_keys[ix];

					if (pk == tag)
						check = std::move(check) && key(ck) == value;
					else
						check = std::move(check) && key(ck) == parent[pk].text();
				}

				if (childCat->exists(std::move(check))) // phew..., narrow escape
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
					std::cerr << "Cannot update child " << childCat->m_name << "." << childTag << " with value " << value << std::endl;
			}

			// finally, update the children
			if (not process.empty())
				childCat->update_value(process, childTag, value);
		}
	}
}

void category::update_value(row *row, size_t column, std::string_view value, bool updateLinked, bool validate)
{
	// make sure we have an index, if possible
	if (m_index == nullptr and m_cat_validator != nullptr)
		m_index = new category_index(this);

	auto &col = m_columns[column];

	std::string_view oldValue;
	for (auto iv = row->m_head; iv != nullptr; iv = iv->m_next)
	{
		assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));

		if (iv->m_column_ix == column)
		{
			oldValue = iv->text();
			break;
		}
	}

	if (value == oldValue) // no need to update
		return;

	std::string oldStrValue{ oldValue };

	// check the value
	if (col.m_validator and validate)
		col.m_validator->operator()(value);

	// If the field is part of the Key for this category, remove it from the index
	// before updating

	bool reinsert = false;
	if (updateLinked and // an update of an Item's value
		m_index != nullptr and key_field_indices().count(column))
	{
		reinsert = m_index->find(row);
		if (reinsert)
			m_index->erase(row);
	}

	// first remove old value with cix
	for (auto iv = row->m_head; iv != nullptr; iv = iv->m_next)
	{
		if (iv->m_column_ix != column)
			continue;

		row->remove(iv);
		delete_item(iv);
		break;
	}

	if (not value.empty())
		row->append(create_item(column, value));

	if (reinsert)
		m_index->insert(row);

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

				// TODO: add code to *NOT* test mandatory fields for Empty

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
			// 	std::cerr << "Parent: " << linked->mParentcategory << " Child: " << linked->m_child_category << std::endl
			// 			  << cond << std::endl;
			// }

			// Now, suppose there are already rows in child that conform to the new value,
			// we then skip this rename

			condition cond_n;

			for (size_t ix = 0; ix < linked->m_parent_keys.size(); ++ix)
			{
				std::string pk = linked->m_parent_keys[ix];
				std::string ck = linked->m_child_keys[ix];

				if (pk == iv->m_tag)
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
					std::cerr << "Will not rename in child category since there are already rows that link to the parent" << std::endl;

				continue;
			}

			for (auto cr : rows)
				cr.assign(childTag, value, false);
		}
	}
}

row *category::clone_row(const row &r)
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

void category::delete_row(row *r)
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

row_handle category::create_copy(row_handle r)
{
	// copy the values
	std::vector<item> items;

	for (item_value *iv = r.m_row->m_head; iv != nullptr; iv = iv->m_next)
		items.emplace_back(m_columns[iv->m_column_ix].m_name, iv->text());

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

	// auto &&[result, inserted] = emplace(items.begin(), items.end());
	// // assert(inserted);

	// return result;
}

// proxy methods for every insertion
category::iterator category::insert_impl(const_iterator pos, row *n)
{
	if (m_index == nullptr and m_cat_validator != nullptr)
		m_index = new category_index(this);

	assert(n != nullptr);
	assert(n->m_next == nullptr);

	if (n == nullptr)
		throw std::runtime_error("Invalid pointer passed to insert");

	try
	{
		// First, make sure all mandatory fields are supplied
		if (m_cat_validator != nullptr)
		{
			for (uint16_t ix = 0; ix < static_cast<uint16_t>(m_columns.size()); ++ix)
			{
				const auto &[column, iv] = m_columns[ix];

				if (iv == nullptr)
					continue;

				bool seen = false;

				for (auto i = n->m_head; i != nullptr; i = i->m_next)
				{
					if (i->m_column_ix == ix)
					{
						iv->operator()(i->text());

						seen = true;
						break;
					}
				}

				if (not seen and iv->m_mandatory)
					throw std::runtime_error("missing mandatory field " + column + " for category " + m_name);
			}

			// if (m_index != nullptr)
			// {
			// 	std::unique_ptr<ItemRow> nr(new ItemRow{nullptr, this, nullptr});
			// 	Row r(nr.get());
			// 	auto keys = keyFields();

			// 	for (auto v = b; v != e; ++v)
			// 	{
			// 		if (keys.count(v->name()))
			// 			r.assign(v->name(), v->value(), true);
			// 	}

			// 	auto test = m_index->find(nr.get());
			// 	if (test != nullptr)
			// 	{
			// 		if (VERBOSE > 1)
			// 			std::cerr << "Not inserting new record in " << m_name << " (duplicate Key)" << std::endl;
			// 		result = test;
			// 		isNew = false;
			// 	}
			// }
		}

		if (m_index != nullptr)
			m_index->insert(n);

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
	catch (const std::exception &e)
	{
		delete_row(n);
		throw;
	}
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

void category::swap_item(size_t column_ix, row_handle &a, row_handle &b)
{
	assert(this == a.m_category);
	assert(this == b.m_category);

	item_value *va = nullptr, *vb = nullptr;

	auto ra = a.m_row;
	auto rb = b.m_row;

	if (ra->m_head != nullptr and ra->m_head->m_column_ix == column_ix)
	{
		va = ra->m_head;
		ra->m_head = va->m_next;
		va->m_next = nullptr;

		if (ra->m_tail == va)
			ra->m_tail = ra->m_head;
	}
	else
	{
		for (auto v = ra->m_head; v->m_next != nullptr; v = v->m_next)
		{
			if (v->m_next->m_column_ix != column_ix)
				continue;
			
			va = v->m_next;
			v->m_next = va->m_next;
			va->m_next = nullptr;

			if (ra->m_tail == va)
				ra->m_tail = v;

			break;
		}
	}

	if (rb->m_head != nullptr and rb->m_head->m_column_ix == column_ix)
	{
		vb = rb->m_head;
		rb->m_head = vb->m_next;
		vb->m_next = nullptr;

		if (rb->m_tail == vb)
			rb->m_tail = rb->m_head;
	}
	else
	{
		for (auto v = rb->m_head; v->m_next != nullptr; v = v->m_next)
		{
			if (v->m_next->m_column_ix != column_ix)
				continue;
			
			vb = v->m_next;
			v->m_next = vb->m_next;
			vb->m_next = nullptr;

			if (rb->m_tail = vb)
				rb->m_tail = v;

			break;
		}
	}

	if (ra->m_head == nullptr)
		ra->m_head = ra->m_tail = vb;
	else
	{
		ra->m_tail->m_next = vb;
		ra->m_tail = vb;
	}

	if (rb->m_head == nullptr)
		rb->m_head = rb->m_tail = va;
	else
	{
		rb->m_tail->m_next = va;
		rb->m_tail = va;
	}
}

void category::reorder_by_index()
{
	if (m_index)
		std::tie(m_head, m_tail) = m_index->reorder();
}

namespace detail
{

	size_t write_value(std::ostream &os, std::string_view value, size_t offset, size_t width)
	{
		if (value.find('\n') != std::string::npos or width == 0 or value.length() > 132) // write as text field
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
			os << value;

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

std::vector<std::string> category::get_tag_order() const
{
	std::vector<std::string> result;
	for (auto &c : m_columns)
		result.push_back("_" + m_name + "." + c.m_name);
	return result;
}

void category::write(std::ostream &os) const
{
	std::vector<uint16_t> order(m_columns.size());
	iota(order.begin(), order.end(), 0);
	write(os, order, false);
}

void category::write(std::ostream &os, const std::vector<std::string> &columns, bool addMissingColumns)
{
	// make sure all columns are present
	for (auto &c : columns)
		add_column(c);

	std::vector<uint16_t> order;
	order.reserve(m_columns.size());

	for (auto &c : columns)
		order.push_back(get_column_ix(c));

	if (addMissingColumns)
	{
		for (size_t i = 0; i < m_columns.size(); ++i)
		{
			if (std::find(order.begin(), order.end(), i) == order.end())
				order.push_back(i);
		}
	}

	write(os, order, true);
}

void category::write(std::ostream &os, const std::vector<uint16_t> &order, bool includeEmptyColumns) const
{
	if (empty())
		return;

	// If the first Row has a next, we need a loop_
	bool needLoop = (m_head->m_next != nullptr);

	if (needLoop)
	{
		os << "loop_" << '\n';

		std::vector<size_t> columnWidths(m_columns.size());

		for (auto cix : order)
		{
			auto &col = m_columns[cix];
			os << '_' << m_name << '.' << col.m_name << ' ' << '\n';
			columnWidths[cix] = 2;
		}

		for (auto r = m_head; r != nullptr; r = r->m_next)
		{
			for (auto v = r->m_head; v != nullptr; v = v->m_next)
			{
				if (v->text().find('\n') == std::string_view::npos)
				{
					size_t l = v->text().length();

					if (not sac_parser::is_unquoted_string(v->text()))
						l += 2;

					if (l > 132)
						continue;

					if (columnWidths[v->m_column_ix] < l + 1)
						columnWidths[v->m_column_ix] = l + 1;
				}
			}
		}

		for (auto r = m_head; r != nullptr; r = r->m_next) // loop over rows
		{
			size_t offset = 0;

			for (size_t cix : order)
			{
				size_t w = columnWidths[cix];

				std::string_view s;
				for (auto iv = r->m_head; iv != nullptr; iv = iv->m_next)
				{
					if (iv->m_column_ix == cix)
					{
						s = iv->text();
						break;
					}
				}

				if (s.empty())
					s = "?";

				size_t l = s.length();
				if (not sac_parser::is_unquoted_string(s))
					l += 2;
				if (l < w)
					l = w;

				if (offset + l > 132 and offset > 0)
				{
					os << '\n';
					offset = 0;
				}

				offset = detail::write_value(os, s, offset, w);

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
		size_t l = 0;

		for (auto &col : m_columns)
		{
			std::string tag = '_' + m_name + '.' + col.m_name;

			if (l < tag.length())
				l = tag.length();
		}

		l += 3;

		for (size_t cix : order)
		{
			auto &col = m_columns[cix];

			os << '_' << m_name << '.' << col.m_name << std::string(l - col.m_name.length() - m_name.length() - 2, ' ');

			std::string_view s;
			for (auto iv = m_head->m_head; iv != nullptr; iv = iv->m_next)
			{
				if (iv->m_column_ix == cix)
				{
					s = iv->text();
					break;
				}
			}

			if (s.empty())
				s = "?";

			size_t offset = l;
			if (s.length() + l >= kMaxLineLength)
			{
				os << '\n';
				offset = 0;
			}

			if (detail::write_value(os, s, offset, 1) != 0)
				os << '\n';
		}
	}

	os << "# " << '\n';
}

bool category::operator==(const category &rhs) const
{
	auto &a = *this;
	auto &b = rhs;

	using namespace std::placeholders; 
	
//	set<std::string> tagsA(a.fields()), tagsB(b.fields());
//	
//	if (tagsA != tagsB)
//		std::cout << "Unequal number of fields" << std::endl;

	auto validator = a.get_validator();
	auto catValidator = validator->get_validator_for_category(a.name());
	if (catValidator == nullptr)
		throw std::runtime_error("missing cat validator");
	
	typedef std::function<int(std::string_view,std::string_view)> compType;
	std::vector<std::tuple<std::string,compType>> tags;
	auto keys = catValidator->m_keys;
	std::vector<size_t> keyIx;
	
	for (auto& tag: a.fields())
	{
		auto iv = catValidator->get_validator_for_item(tag);
		if (iv == nullptr)
			throw std::runtime_error("missing item validator");
		auto tv = iv->m_type;
		if (tv == nullptr)
			throw std::runtime_error("missing type validator");
		tags.push_back(std::make_tuple(tag, std::bind(&cif::type_validator::compare, tv, std::placeholders::_1, std::placeholders::_2)));
		
		auto pred = [tag](const std::string& s) -> bool { return cif::iequals(tag, s) == 0; };
		if (find_if(keys.begin(), keys.end(), pred) == keys.end())
			keyIx.push_back(tags.size() - 1);
	}
	
	// a.reorderByIndex();
	// b.reorderByIndex();
	
	auto rowEqual = [&](const row_handle& a, const row_handle& b)
	{
		int d = 0;

		for (auto kix: keyIx)
		{
			std::string tag;
			compType compare;
			
			std::tie(tag, compare) = tags[kix];

			d = compare(a[tag].text(), b[tag].text());

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
		
		for (auto& tt: tags)
		{
			std::string tag;
			compType compare;
			
			std::tie(tag, compare) = tt;
			
			// make it an option to compare unapplicable to empty or something
			
			auto ta = ra[tag].text();	if (ta == "." or ta == "?") ta = "";
			auto tb = rb[tag].text();	if (tb == "." or tb == "?") tb = "";
			
			if (compare(ta, tb) != 0)
				return false;
		}
		
		++ai;
		++bi;
	}

	return true;
}

} // namespace cif