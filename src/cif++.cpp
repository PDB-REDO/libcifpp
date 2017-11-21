// cif parsing library

#include <cassert>

#include <stack>
#include <tuple>
#include <regex>
#include <set>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#if defined(USE_RSRC)
#include "mrsrc.h"
#endif

#include "cif++.h"
#include "cif-parser.h"
#include "cif-validator.h"
#include "cif-utils.h"

using namespace std;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;

extern int VERBOSE;

namespace cif
{

static const char* kEmptyResult = "";
	
// --------------------------------------------------------------------
// most internal data structures are stored as linked lists
// item values are stored in a simple struct. They should be const anyway

struct item_value
{
	item_value*				m_next;
	uint32					m_column_index;
	char					m_text[0];
	
	item_value(const char* v, uint32 column_index);
	~item_value();
	
	void* operator new(size_t size, size_t data_size);
	void operator delete(void* p);
};

// --------------------------------------------------------------------

item_value::item_value(const char* value, uint32 column_index)
	: m_next(nullptr), m_column_index(column_index)
{
	strcpy(m_text, value);
}

item_value::~item_value()
{
	// remove recursion (and be paranoid)
	while (m_next != nullptr and m_next != this)
	{
		auto n = m_next;
		m_next = n->m_next;
		n->m_next = nullptr;
		delete n;
	}
}

void* item_value::operator new(size_t size, size_t data_size)
{
	return malloc(size + data_size + 1);
}

void item_value::operator delete(void* p)
{
	free(p);
}

// --------------------------------------------------------------------

// item_column contains info about a column or field in a category

struct item_column
{
	string					m_name;		// store lower-case, for optimization
	const validate_item*	m_validator;
};

// item_row contains the actual values for a row in a category

struct item_row
{
	~item_row();

	void drop(uint32 column_ix);
	const char* c_str(uint32 column_ix) const;
	
	string str() const
	{
		stringstream s;

		s << '{';
		for (auto v = m_values; v != nullptr; v = v->m_next)
		{
			s << m_category->get_column_name(v->m_column_index)
			  << ':'
			  << v->m_text;
			 if (v->m_next != nullptr)
			 	s << ", ";
		}
		s << '}';

		return s.str();
	}
	
	item_row*				m_next;
	category*				m_category;
	item_value*				m_values;
};

ostream& operator<<(ostream& os, const item_row& r)
{
	os << r.m_category->name() << '[';
	for (auto iv = r.m_values; iv != nullptr; iv = iv->m_next)
	{
		os << iv->m_text;
		if (iv->m_next)
			os << ',';
	}
	os << ']';
	
	return os;
}

// --------------------------------------------------------------------

item_row::~item_row()
{
	// remove recursive
	while (m_next != nullptr and m_next != this)
	{
		auto n = m_next;
		m_next = n->m_next;
		n->m_next = nullptr;
		delete n;
	}

	delete m_values;
}

void item_row::drop(uint32 column_ix)
{
	if (m_values != nullptr and m_values->m_column_index == column_ix)
	{
		auto v = m_values;
		m_values = m_values->m_next;
		v->m_next = nullptr;
		delete v;
	}
	else
	{
		for (auto v = m_values; v->m_next != nullptr; v = v->m_next)
		{
			if (v->m_next->m_column_index == column_ix)
			{
				auto vn = v->m_next;
				v->m_next = vn->m_next;
				vn->m_next = nullptr;
				delete vn;

				break;
			}
		}
	}

#if DEBUG
	for (auto iv = m_values; iv != nullptr; iv = iv->m_next)
		assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));
	
#endif
}

const char* item_row::c_str(uint32 column_ix) const
{
	const char* result = kEmptyResult;
	
	for (auto v = m_values; v != nullptr; v = v->m_next)
	{
		if (v->m_column_index == column_ix)
		{
			result = v->m_text;
			break;
		}
	}

	return result;
}

// --------------------------------------------------------------------

namespace detail
{

template<>
item_reference& item_reference::operator=(const string& value)
{
	row(m_row).assign(m_name, value, false);
	return *this;
}

const char*
item_reference::c_str() const
{
	const char* result = kEmptyResult;
	
	if (m_row != nullptr /* and m_row->m_category != nullptr*/)
	{
//		assert(m_row->m_category);
		
		auto cix = m_row->m_category->get_column_index(m_name);
		
		for (auto iv = m_row->m_values; iv != nullptr; iv = iv->m_next)
		{
			if (iv->m_column_index == cix)
			{
				if (iv->m_text[0] != '.' or iv->m_text[1] != 0)
					result = iv->m_text;
	
				break;
			}
		}
	}
	
	return result;
}

bool item_reference::empty() const
{
	return c_str() == kEmptyResult;
}

}

// --------------------------------------------------------------------
// datablock implementation

datablock::datablock(const string& name)
	: m_name(name), m_validator(nullptr), m_next(nullptr)
{
}

datablock::~datablock()
{
	delete m_next;
}

string datablock::first_item(const string& tag) const
{
	string result;

	string cat_name, item_name;
	std::tie(cat_name, item_name) = split_tag_name(tag);
	
	for (auto& cat: m_categories)
	{
		if (iequals(cat.name(), cat_name))
		{
			result = cat.get_first_item(item_name.c_str()).as<string>();
			break;
		}
	}

	return result;
}

auto datablock::emplace(const string& name) -> tuple<iterator,bool>
{
	bool isNew = false;
	iterator i = find_if(begin(), end(), [name](const category& cat) -> bool
		{ return iequals(cat.name(), name); });
	
	if (i == end())
	{
		isNew = true;
		i = m_categories.emplace(end(), *this, name, m_validator);
	}
	
	return make_tuple(i, isNew);
}

category& datablock::operator[](const string& name)
{
	iterator i;
	std::tie(i, ignore) = emplace(name);
	return *i;
}

category* datablock::get(const string& name)
{
	auto i = find_if(begin(), end(), [name](const category& cat) -> bool 
		{ return iequals(cat.name(), name); });
	
	return i == end() ? nullptr : &*i;
}

void datablock::validate()
{
	if (m_validator == nullptr)
		throw runtime_error("validator not specified");

	for (auto& cat: *this)
		cat.validate();
}

void datablock::set_validator(validator* v)
{
	m_validator = v;

	for (auto& cat: *this)
		cat.set_validator(v);
}

void datablock::get_tag_order(vector<string>& tags) const
{
	for (auto& cat: *this)
		cat.get_tag_order(tags);
}

void datablock::write(ostream& os)
{
	os << "data_" << m_name << endl
	   << "# " << endl;
	
	// mmcif support, sort of. First write the 'entry' category
	// and if it exists, _AND_ we have a validator, write out the
	// audit_conform record.

	for (auto& cat: m_categories)
	{
		if (cat.name() == "entry")
		{
			cat.write(os);
			
			if (m_validator != nullptr)
			{
				category audit_conform(*this, "audit_conform", nullptr);
				audit_conform.emplace({
					{ "dict_name", m_validator->dict_name() },
					{ "dict_version", m_validator->dict_version() }
				});
				audit_conform.write(os);
			}
			
			break;
		}
	}

	for (auto& cat: m_categories)
	{
		if (cat.name() != "entry" and cat.name() != "audit_conform")
			cat.write(os);
	}
}

void datablock::write(ostream& os, const vector<string>& order)
{
	os << "data_" << m_name << endl
	   << "# " << endl;
	
	vector<string> catOrder;
	for (auto& o: order)
	{
		string cat, item;
		std::tie(cat, item) = split_tag_name(o);
		if (find_if(catOrder.rbegin(), catOrder.rend(), [cat](const string& s) -> bool { return iequals(cat, s); }) == catOrder.rend())
			catOrder.push_back(cat);
	}

	for (auto& c: catOrder)
	{
		auto cat = get(c);
		if (cat == nullptr)
			continue;
		
		vector<string> items;
		for (auto& o: order)
		{
			string cat_name, item;
			std::tie(cat_name, item) = split_tag_name(o);
			
			if (cat_name == c)
				items.push_back(item);
		}
		
		cat->write(os, items);
	}
	
	// for any category we missed in the catOrder
	for (auto& cat: m_categories)
	{
		if (find_if(catOrder.begin(), catOrder.end(), [&](const string& s) -> bool { return iequals(cat.name(), s); }) != catOrder.end())
			continue;
		
		cat.write(os);
	}
	
	
//	// mmcif support, sort of. First write the 'entry' category
//	// and if it exists, _AND_ we have a validator, write out the
//	// audit_conform record.
//
//	for (auto& cat: m_categories)
//	{
//		if (cat.name() == "entry")
//		{
//			cat.write(os);
//			
//			if (m_validator != nullptr)
//			{
//				category audit_conform(*this, "audit_conform", nullptr);
//				audit_conform.emplace({
//					{ "dict_name", m_validator->dict_name() },
//					{ "dict_version", m_validator->dict_version() }
//				});
//				audit_conform.write(os);
//			}
//			
//			break;
//		}
//	}
//
//	for (auto& cat: m_categories)
//	{
//		if (cat.name() != "entry" and cat.name() != "audit_conform")
//			cat.write(os);
//	}
}

// --------------------------------------------------------------------
//
//	class to compare two rows based on their keys.

class row_comparator
{
  public:

	row_comparator(category* cat)
		: row_comparator(cat, cat->get_cat_validator()->m_keys.begin(), cat->get_cat_validator()->m_keys.end())
	{
	}
		
	template<typename KeyIter>
	row_comparator(category* cat, KeyIter b, KeyIter e);
	
	int operator()(const item_row* a, const item_row* b) const;

	int operator()(const row& a, const row& b) const
	{
		return operator()(a.m_data, b.m_data);
	}
	
  private:
	typedef function<int(const char*,const char*)>	compare_func;

	typedef tuple<size_t,compare_func>	key_comp;

	vector<key_comp>	m_comp;
};

template<typename KeyIter>
row_comparator::row_comparator(category* cat, KeyIter b, KeyIter e)
{
	auto cv = cat->get_cat_validator();
	
	for (auto ki = b; ki != e; ++ki)
	{
		string k = *ki;
		
		size_t ix = cat->get_column_index(k);

		auto iv = cv->get_validator_for_item(k);
		if (iv == nullptr)
			throw runtime_error("Incomplete dictionary, no item validator for key " + k);
		
		auto tv = iv->m_type;
		if (tv == nullptr)
			throw runtime_error("Incomplete dictionary, no type validator for item " + k);
		
		using namespace placeholders;
		
		m_comp.emplace_back(ix, bind(&validate_type::compare, tv, _1, _2));
	}
}

int row_comparator::operator()(const item_row* a, const item_row* b) const
{
	assert(a);
	assert(b);

	int d = 0;
	for (auto& c: m_comp)
	{
		size_t k;
		compare_func f;
		
		std::tie(k, f) = c;
		
		const char* ka = a->c_str(k);
		const char* kb = b->c_str(k);
		
		d = f(ka, kb);

		if (d != 0)
			break;
	}
	
	return d;
}

// --------------------------------------------------------------------
//
//	class to keep an index on the keys of a category. This is a red/black
//	tree implementation.

class cat_index
{
  public:
	cat_index(category* cat);
	~cat_index();
	
	item_row* find(item_row* k) const;

	void insert(item_row* r);
	void erase(item_row* r);

	// batch create
	void reconstruct();
	
	// reorder the item_row's and returns new head and tail
	tuple<item_row*,item_row*> reorder()
	{
		tuple<item_row*,item_row*> result = make_tuple(nullptr, nullptr);
		
		if (m_root != nullptr)
		{
			entry* head = findMin(m_root);
			entry* tail = reorder(m_root);
			
			tail->m_row->m_next = nullptr;
			
			result = make_tuple(head->m_row, tail->m_row);
		}
			
		return result;
	}

	size_t size() const;
	void validate() const;
	
  private:

	struct entry
	{
		entry(item_row* r)
			: m_row(r), m_left(nullptr), m_right(nullptr), m_red(true) {}
		
		~entry()
		{
			delete m_left;
			delete m_right;
		}

		item_row*		m_row;
		entry*			m_left;
		entry*			m_right;
		bool			m_red;
	};

	entry* insert(entry* h, item_row* v);
	entry* erase(entry* h, item_row* k);

	void validate(entry* h, bool isParentRed, uint32 blackDepth, uint32& minBlack, uint32& maxBlack) const;

	entry* rotateLeft(entry* h)
	{
		entry* x = h->m_right;
		h->m_right = x->m_left;
		x->m_left = h;
		x->m_red = h->m_red;
		h->m_red = true;
		return x;
	}
	
	entry* rotateRight(entry* h)
	{
		entry* x = h->m_left;
		h->m_left = x->m_right;
		x->m_right = h;
		x->m_red = h->m_red;
		h->m_red = true;
		return x;
	}
	
	void flipColour(entry* h)
	{
		h->m_red = not h->m_red;
		
		if (h->m_left != nullptr)
			h->m_left->m_red = not h->m_left->m_red;
	
		if (h->m_right != nullptr)
			h->m_right->m_red = not h->m_right->m_red;
	}
	
	bool isRed(entry* h) const
	{
		return h != nullptr and h->m_red;
	}
	
	entry* moveRedLeft(entry* h)
	{
		flipColour(h);
		
		if (h->m_right != nullptr and isRed(h->m_right->m_left))
		{
			h->m_right = rotateRight(h->m_right);
			h = rotateLeft(h);
			flipColour(h);
		}
		
		return h;
	}
	
	entry* moveRedRight(entry* h)
	{
		flipColour(h);
		
		if (h->m_left != nullptr and isRed(h->m_left->m_left))
		{
			h = rotateRight(h);
			flipColour(h);
		}
		
		return h;
	}
	
	entry* fixUp(entry* h)
	{
		if (isRed(h->m_right))
			h = rotateLeft(h);
		
		if (isRed(h->m_left) and isRed(h->m_left->m_left))
			h = rotateRight(h);
		
		if (isRed(h->m_left) and isRed(h->m_right))
			flipColour(h);
		
		return h;
	}
	
	entry* findMin(entry* h)
	{
		while (h->m_left != nullptr)
			h = h->m_left;

		return h;
	}
	
	entry* eraseMin(entry* h)
	{
		if (h->m_left == nullptr)
		{
			delete h;
			h = nullptr;
		}
		else
		{
			if (not isRed(h->m_left) and not isRed(h->m_left->m_left))
				h = moveRedLeft(h);
			
			h->m_left = eraseMin(h->m_left);
			
			h = fixUp(h);
		}
		
		return h;
	}
	
	// Fix m_next fields for rows in order of this index
	entry* reorder(entry* e)
	{
		auto result = e;
		
		if (e->m_left != nullptr)
		{
			auto l = reorder(e->m_left);
			l->m_row->m_next = e->m_row;
		}
		
		if (e->m_right != nullptr)
		{
			auto mr = findMin(e->m_right);
			e->m_row->m_next = mr->m_row;
			
			result = reorder(e->m_right);
		}
		
		return result;
	}
	
	category&			m_cat;
	row_comparator		m_comp;
	entry*				m_root;
};

cat_index::cat_index(category* cat)
	: m_cat(*cat), m_comp(cat), m_root(nullptr)
{
}

cat_index::~cat_index()
{
	delete m_root;
}

item_row* cat_index::find(item_row* k) const
{
	const entry* r = m_root;
	while (r != nullptr)
	{
		int d = m_comp(k, r->m_row);
		if (d < 0)
			r = r->m_left;
		else if (d > 0)
			r = r->m_right;
		else
			break;
	}
	
	return r ? r->m_row : nullptr;
}

void cat_index::insert(item_row* k)
{
	m_root = insert(m_root, k);
	m_root->m_red = false;
}

cat_index::entry* cat_index::insert(entry* h, item_row* v)
{
	if (h == nullptr)
		return new entry(v);
	
	int d = m_comp(v, h->m_row);
	if (d < 0)		h->m_left = insert(h->m_left, v);
	else if (d > 0)	h->m_right = insert(h->m_right, v);
	else
		throw runtime_error("Duplicate key violation, cat: " + m_cat.name() + " values: " + v->str());

	if (isRed(h->m_right) and not isRed(h->m_left))
		h = rotateLeft(h);

	if (isRed(h->m_left) and isRed(h->m_left->m_left))
		h = rotateRight(h);
	
	if (isRed(h->m_left) and isRed(h->m_right))	
		flipColour(h);
	
	return h;
}

void cat_index::erase(item_row* k)
{
	m_root = erase(m_root, k);
	if (m_root != nullptr)
		m_root->m_red = false;
}

cat_index::entry* cat_index::erase(entry* h, item_row* k)
{
	if (m_comp(k, h->m_row) < 0)
	{
		if (h->m_left != nullptr)
		{
			if (not isRed(h->m_left) and not isRed(h->m_left->m_left))
				h = moveRedLeft(h);

			h->m_left = erase(h->m_left, k);
		}
	}
	else
	{
		if (isRed(h->m_left))
			h = rotateRight(h);
			
		if (m_comp(k, h->m_row) == 0 and h->m_right == nullptr)
		{
			delete h;
			return nullptr;
		}
		
		if (h->m_right != nullptr)
		{
			if (not isRed(h->m_right) and not isRed(h->m_right->m_left))
				h = moveRedRight(h);
			
			if (m_comp(k, h->m_row) == 0)
			{
				h->m_row = findMin(h->m_right)->m_row;
				h->m_right = eraseMin(h->m_right);
			}
			else
				h->m_right = erase(h->m_right, k);
		}
	}
	
	return fixUp(h);
}

void cat_index::reconstruct()
{
	delete m_root;
	m_root = nullptr;
	
	for (auto r: m_cat)
		insert(r.m_data);

// maybe reconstruction can be done quicker by using the following commented code.
// however, I've not had the time to think of a way to set the red/black flag correctly in that case.
	
//	vector<item_row*> rows;
//	transform(m_cat.begin(), m_cat.end(), back_inserter(rows),
//		[](row r) -> item_row* { assert(r.m_data); return r.m_data; });
//	
//	assert(std::find(rows.begin(), rows.end(), nullptr) == rows.end());
//	
//	// don't use sort here, it will run out of the stack of something.
//	// quicksort is notorious for using excessive recursion.
//	// Besides, most of the time, the data is ordered already anyway.
//
//	stable_sort(rows.begin(), rows.end(), [this](item_row* a, item_row* b) -> bool { return this->m_comp(a, b) < 0; });
//	
//	for (size_t i = 0; i < rows.size() - 1; ++i)
//		assert(m_comp(rows[i], rows[i + 1]) < 0);
//	
//	deque<entry*> e;
//	transform(rows.begin(), rows.end(), back_inserter(e),
//		[](item_row* r) -> entry* { return new entry(r); });
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
//				b->m_left = a;
//				
//				assert(m_comp(a->m_row, b->m_row) < 0);
//
//				e.pop_front();
//				
//				if (not e.empty())
//				{
//					entry* c = e.front();
//					e.pop_front();
//
//					assert(m_comp(b->m_row, c->m_row) < 0);
//				
//					b->m_right = c;
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
//	m_root = e.front();
}

size_t cat_index::size() const
{
	stack<entry*> s;
	s.push(m_root);
	
	size_t result = 0;
	
	while (not s.empty())
	{
		entry* e = s.top();
		s.pop();
		
		if (e == nullptr)
			continue;
		
		++result;
		
		s.push(e->m_left);
		s.push(e->m_right);
	}
	
	return result;
}

void cat_index::validate() const
{
	if (m_root != nullptr)
	{
		uint32 minBlack = numeric_limits<uint32>::max();
		uint32 maxBlack = 0;
		
		assert(not m_root->m_red);
		
		validate(m_root, false, 0, minBlack, maxBlack);
		assert(minBlack == maxBlack);
	}
}

void cat_index::validate(entry* h, bool isParentRed, uint32 blackDepth, uint32& minBlack, uint32& maxBlack) const
{
	if (h->m_red)
		assert(not isParentRed);
	else
		++blackDepth;
	
	if (isParentRed)
		assert(not h->m_red);
	
	if (h->m_left != nullptr and h->m_right != nullptr)
	{
		if (isRed(h->m_left))
			assert(not isRed(h->m_right));
		if (isRed(h->m_right))
			assert(not isRed(h->m_left));
	}
	
	if (h->m_left != nullptr)
	{
		assert(m_comp(h->m_left->m_row, h->m_row) < 0);
		validate(h->m_left, h->m_red, blackDepth, minBlack, maxBlack);
	}
	else
	{
		if (minBlack > blackDepth)
			minBlack = blackDepth;
		if (maxBlack < blackDepth)
			maxBlack = blackDepth;
	}
	
	if (h->m_right != nullptr)
	{
		assert(m_comp(h->m_right->m_row, h->m_row) > 0);
		validate(h->m_right, h->m_right, blackDepth, minBlack, maxBlack);
	}
	else
	{
		if (minBlack > blackDepth)
			minBlack = blackDepth;
		if (maxBlack < blackDepth)
			maxBlack = blackDepth;
	}
}

// --------------------------------------------------------------------

rowset::rowset(category& cat)
	: m_cat(cat)
{
}

rowset& rowset::orderBy(initializer_list<string> items)
{
	row_comparator c(&m_cat, items.begin(), items.end());
	
	stable_sort(begin(), end(), c);
	
	return *this;
}

// --------------------------------------------------------------------

category::category(datablock& db, const string& name, validator* validator)
	: m_db(db), m_name(name), m_validator(validator)
	, m_head(nullptr), m_tail(nullptr), m_index(nullptr)
{
	if (m_name.empty())
		throw validation_error("invalid empty name for category");
	
	if (m_validator != nullptr)
	{
		m_cat_validator = m_validator->get_validator_for_category(m_name);
		if (m_cat_validator != nullptr)
		{
			// make sure all required columns are added
			
			for (auto& k: m_cat_validator->m_keys)
				add_column(k);

			for (auto& k: m_cat_validator->m_mandatory_fields)
				add_column(k);
			
			m_index = new cat_index(this);
		}
	}
}

category::~category()
{
	delete m_head;
	delete m_index;
}

void category::set_validator(validator* v)
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
			m_index = new cat_index(this);
			m_index->reconstruct();
#if DEBUG
			assert(m_index->size() == size());
			m_index->validate();
#endif
		}
	}
	else
		m_cat_validator = nullptr;
}

size_t category::get_column_index(const string& name) const
{
	size_t result;

	for (result = 0; result < m_columns.size(); ++result)
	{
		if (iequals(name, m_columns[result].m_name))
			break;
	}
	
	return result;
}

const string& category::get_column_name(size_t column_ix) const
{
	return m_columns.at(column_ix).m_name;
}

size_t category::add_column(const string& name)
{
	size_t result = get_column_index(name);
	
	if (result == m_columns.size())
	{
		const validate_item* item_validator = nullptr;
		
		if (m_cat_validator != nullptr)
		{
			item_validator = m_cat_validator->get_validator_for_item(name);
			if (item_validator == nullptr)
				m_validator->report_error("tag " + name + " not allowed in category " + m_name);
		}
		
		m_columns.push_back({name, item_validator});
	}
	
	return result;
}

void category::reorderByIndex()
{
	if (m_index != nullptr)
		std::tie(m_head, m_tail) = m_index->reorder();
}

size_t category::size() const
{
	size_t result = 0;
	
	for (auto pi = m_head; pi != nullptr; pi = pi->m_next)
		++result;
	
	return result;
}

bool category::empty() const
{
	return m_head == nullptr or m_head->m_values == nullptr;
}

void category::drop(const string& field)
{
	using namespace placeholders;
	auto ci = find_if(m_columns.begin(), m_columns.end(),
		[field](item_column& c) -> bool { return iequals(c.m_name, field); });

	if (ci != m_columns.end())
	{
		uint32 column_ix = ci - m_columns.begin();
		
		for (auto pi = m_head; pi != nullptr; pi = pi->m_next)
			pi->drop(column_ix);
		
		m_columns.erase(ci);
	}
}

row category::operator[](condition&& cond)
{
	row result;
	
	for (auto r: *this)
	{
		if (cond(*this, r))
		{
			result = r;
			break;
		}
	}

	return result;
}	

rowset category::find(condition&& cond)
{
	rowset result(*this);
	for (auto r: *this)
	{
		if (cond(*this, r))
			result.push_back(r);
	}
	return result;
}

bool category::exists(condition&& cond)
{
	bool result = false;
	
	for (auto r: *this)
	{
		if (cond(*this, r))
		{
			result = true;
			break;
		}
	}

	return result;
}

rowset category::orderBy(std::initializer_list<string> items)
{
	rowset result(*this);
	result.insert(result.begin(), begin(), end());
	
	return result.orderBy(items);
}

void category::clear()
{
	delete m_head;
	m_head = m_tail = nullptr;
	
	if (m_index != nullptr)
	{
		delete m_index;
		m_index = new cat_index(this);
	}
}

template<class Iter>
tuple<row,bool> category::emplace(Iter b, Iter e)
{
	// First, make sure all mandatory fields are supplied
	tuple<row,bool> result = make_tuple(row(), true);

	if (m_cat_validator != nullptr and b != e)
	{
		for (auto& col: m_columns)
		{
			auto iv = m_cat_validator->get_validator_for_item(col.m_name);
	
			if (iv == nullptr)
				continue;
			
			bool seen = false;
			
			for (auto v = b; v != e; ++v)
			{
				if (iequals(v->name(), col.m_name))
				{
					seen = true;
					break;
				}
			}
			
			if (not seen and iv->m_mandatory)
				throw runtime_error("missing mandatory field " + col.m_name + " for category " + m_name);
		}
		
		if (m_index != nullptr)
		{
			unique_ptr<item_row> nr(new item_row{nullptr, this, nullptr});
			row r(nr.get());
			auto keys = key_fields(); 
			
			for (auto v = b; v != e; ++v)
			{
				if (keys.count(v->name()))
					r.assign(v->name(), v->value(), true);
			}
			
			auto test = m_index->find(nr.get());
			if (test != nullptr)
			{
				if (VERBOSE > 1)
					cerr << "Not inserting new record in " << m_name << " (duplicate key)" << endl;
				result = make_tuple(row(test), false);
			}
		}
	}
	
	if (get<1>(result))
	{
		auto nr = new item_row{nullptr, this, nullptr};

		if (m_head == nullptr)
		{
			assert(m_tail == nullptr);
			m_head = m_tail = nr;
		}
		else
		{
			assert(m_tail != nullptr);
			assert(m_head != nullptr);
			m_tail->m_next = nr;
			m_tail = nr;
		}

		row r(nr);

		for (auto v = b; v != e; ++v)
			r.assign(*v, true);
		
		get<0>(result) = r;

		if (m_index != nullptr)
			m_index->insert(nr);
	}
	
	return result;
}

tuple<row,bool> category::emplace(row r)
{
	return emplace(r.begin(), r.end());
}

void category::erase(condition&& cond)
{
	rowset remove(*this);

	for (auto r: *this)
	{
		if (cond(*this, r))
			remove.push_back(r);
	}

	for (auto r: remove)
		erase(r);
}

void category::erase(iterator p)
{
	erase(*p);
}

void category::erase(row r)
{
	iset keys;
	if (m_cat_validator)
		keys = iset(m_cat_validator->m_keys.begin(), m_cat_validator->m_keys.end());
	
	for (auto& col: m_columns)
	{
		auto iv = col.m_validator;
		if (iv == nullptr or iv->m_children.empty())
			continue;
		
		if (not keys.count(col.m_name))
			continue;
		
		const char* value = r[col.m_name].c_str();
		
		for (auto child: iv->m_children)
		{
			if (child->m_category == nullptr)
				continue;
			
			auto child_cat = m_db.get(child->m_category->m_name);
			if (child_cat == nullptr)
				continue;
				
			auto rows = child_cat->find(key(child->m_tag) == value);
			for (auto& cr: rows)
				child_cat->erase(cr);
		}
	}

	if (m_head == nullptr)
		throw runtime_error("erase");

	if (m_index != nullptr)
		m_index->erase(r.m_data);
	
	if (r == m_head)
	{
		m_head = m_head->m_next;
		r.m_data->m_next = nullptr;
		delete r.m_data;
	}
	else
	{
		for (auto pi = m_head; pi != nullptr; pi = pi->m_next)
		{
			if (pi->m_next == r.m_data)
			{
				pi->m_next = r.m_data->m_next;
				r.m_data->m_next = nullptr;
				delete r.m_data;
				break;
			}
		}
	}
}

void category::get_tag_order(vector<string>& tags) const
{
	for (auto& c: m_columns)
		tags.push_back("_" + m_name + "." + c.m_name);
}

const detail::item_reference category::get_first_item(const char* item_name) const
{
	return detail::item_reference{item_name, m_head};
}

category::iterator category::begin()
{
	return iterator(m_head);
}

category::iterator category::end()
{
	return iterator(nullptr);
}

void category::validate()
{
	if (m_validator == nullptr)
		throw runtime_error("no validator specified");

	if (empty())
	{
		if (VERBOSE > 2)
			cerr << "Skipping validation of empty category " << m_name << endl;
		return;
	}
	
	if (m_cat_validator == nullptr)
	{
		m_validator->report_error("undefined category " + m_name);
		return;
	}
	
	auto mandatory = m_cat_validator->m_mandatory_fields;

	for (auto& col: m_columns)
	{
		auto iv = m_cat_validator->get_validator_for_item(col.m_name);
		if (iv == nullptr)
			m_validator->report_error("Field " + col.m_name + " is not valid in category " + m_name);
		
		col.m_validator = iv;
		
		mandatory.erase(col.m_name);
	}
	
	if (not mandatory.empty())
		m_validator->report_error("In category " + m_name + " the following mandatory fields are missing: " + ba::join(mandatory, ", "));
	
	// check index?
	if (m_index)
	{
#if not defined(NDEBUG)
		m_index->validate();
		for (auto r: *this)
		{
			if (m_index->find(r.m_data) != r.m_data)
				m_validator->report_error("Key not found in index for category " + m_name);
		}
#endif
	}
	
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
				m_validator->report_error("invalid field " + m_columns[cix].m_name + " for category " + m_name);
				continue;
			}
			
			for (auto vi = ri->m_values; vi != nullptr; vi = vi->m_next)
			{
				if (vi->m_column_index == cix)
				{
					seen = true;
 					(*iv)(vi->m_text);
				}
			}
			
			if (seen)
				continue;
			
			if (iv != nullptr and iv->m_mandatory)
				m_validator->report_error("missing mandatory field " + m_columns[cix].m_name + " for category " + m_name);
		}
	}
}

const validator& category::get_validator() const
{
	if (m_validator == nullptr)
		throw runtime_error("no validator defined yet");
	return *m_validator;
}

iset category::fields() const
{
	if (m_validator == nullptr)
		throw runtime_error("No validator specified");
	
	if (m_cat_validator == nullptr)
		m_validator->report_error("undefined category");
	
	iset result;
	for (auto& iv: m_cat_validator->m_item_validators)
		result.insert(iv.m_tag);
	return result;
}

iset category::mandatory_fields() const
{
	if (m_validator == nullptr)
		throw runtime_error("No validator specified");
	
	if (m_cat_validator == nullptr)
		m_validator->report_error("undefined category");
	
	return m_cat_validator->m_mandatory_fields;
}

iset category::key_fields() const
{
	if (m_validator == nullptr)
		throw runtime_error("No validator specified");
	
	if (m_cat_validator == nullptr)
		m_validator->report_error("undefined category");
	
	return iset{ m_cat_validator->m_keys.begin(), m_cat_validator->m_keys.end() };
}

auto category::iterator::operator++() -> iterator&
{
	m_current = row(m_current.data()->m_next);
	return *this;
}

namespace detail
{

size_t write_value(ostream& os, string value, size_t offset, size_t width)
{
	if (value.find('\n') != string::npos or width == 0 or value.length() >= 132)					// write as text field
	{
		ba::replace_all(value, "\n;", "\n\\;");

		if (offset > 0)
			os << endl;
		os << ';' << value;
		if (not ba::ends_with(value, "\n"))
			os << endl;
		os << ';' << endl;
		offset = 0;
	}
	else if (is_unquoted_string(value.c_str()))
	{
		os << value;

		if (value.length() < width)
		{
			os << string(width - value.length(), ' ');
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
		for (char q: { '\'', '"'})
		{
			auto p = value.find(q);	// see if we can use the quote character
			while (p != string::npos and is_non_blank(value[p + 1]) and value[p + 1] != q)
				p = value.find(q, p + 1);
			
			if (p != string::npos)
				continue;
			
			os << q << value << q;

			if (value.length() + 2 < width)
			{
				os << string(width - value.length() - 2, ' ');
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
				os << endl;
			os << ';' << value << endl
			   << ';' << endl;
			offset = 0;
		}
	}

	return offset;
}
	
}

void category::write(ostream& os, const vector<int>& order, bool includeEmptyColumns)
{
	if (empty())
		return;
	
	// If the first row has a next, we need a loop_
	bool need_loop = (m_head->m_next != nullptr);
	
	if (need_loop)
	{
		os << "loop_" << endl;
		
		vector<size_t> column_widths;
		
		for (auto cix: order)
		{
			auto& col = m_columns[cix];
			os << '_' << m_name << '.' << col.m_name << ' ' << endl;
			column_widths.push_back(2);
		}
		
		for (auto row = m_head; row != nullptr; row = row->m_next)
		{
			for (auto v = row->m_values; v != nullptr; v = v->m_next)
			{
				if (strchr(v->m_text, '\n') == nullptr)
				{
					size_t l = strlen(v->m_text);
					
					if (not is_unquoted_string(v->m_text))
						l += 2;

					if (l >= 132)
						continue;

					if (column_widths[v->m_column_index] < l + 1)
						column_widths[v->m_column_index] = l + 1;
				}
			}
		}
		
		for (auto row = m_head; row != nullptr; row = row->m_next)	// loop over rows
		{
			size_t offset = 0;
		
			for (size_t cix: order)
			{
				size_t w = column_widths[cix];
				
				string s;
				for (auto iv = row->m_values; iv != nullptr; iv = iv->m_next)
				{
					if (iv->m_column_index == cix)
					{
						s = iv->m_text;
						break;
					}
				}
				
				if (s.empty())
					s = "?";
				
				size_t l = s.length();
				if (not is_unquoted_string(s.c_str()))
					l += 2;
				if (l < w)
					l = w;

				if (offset + l >= 132 and offset > 0)
				{
					os << endl;
					offset = 0;
				}
				
				offset = detail::write_value(os, s, offset, w);
				
				if (offset >= 132)
				{
					os << endl;
					offset = 0;
				}
			}
			
			if (offset > 0)
				os << endl;
		}
	}
	else
	{
		// first find the indent level 
		size_t l = 0;
		
		for (auto& col: m_columns)
		{
			string tag = '_' + m_name + '.' + col.m_name;
			
			if (l < tag.length())
				l = tag.length();
		}
		
		l += 3;
		
		for (size_t cix: order)
		{
			auto& col = m_columns[cix];
			
			os << '_' << m_name << '.' << col.m_name << string(l - col.m_name.length() - m_name.length() - 2, ' ');
			
			string s;
			for (auto iv = m_head->m_values; iv != nullptr; iv = iv->m_next)
			{
				if (iv->m_column_index == cix)
				{
					s = iv->m_text;
					break;
				}
			}
			
			if (s.empty())
				s = "?";
			
			size_t offset = l;
			if (s.length() + l >= kMaxLineLength)
			{
				os << endl;
				offset = 0;
			}

			if (detail::write_value(os, s, offset, 1) != 0)
				os << endl;
		}
	}

	os << "# " << endl;
}

void category::write(ostream& os)
{
	vector<int> order(m_columns.size());
	iota(order.begin(), order.end(), 0);
	write(os, order, false);
}

void category::write(ostream& os, const vector<string>& columns)
{
	// make sure all columns are present
	for (auto& c: columns)
		add_column(c);
	
	vector<int> order;
	order.reserve(m_columns.size());

	for (auto& c: columns)
		order.push_back(get_column_index(c));

	for (size_t i = 0; i < m_columns.size(); ++i)
	{
		if (std::find(order.begin(), order.end(), i) == order.end())
			order.push_back(i);
	}

	write(os, order, true);
}

// --------------------------------------------------------------------

row::row(const row& rhs)
	: m_data(rhs.m_data)
{
}

row& row::operator=(const row& rhs)
{
	m_data = rhs.m_data;
	return *this;
}

void row::assign(const string& name, const string& value, bool emplacing)
{
	if (m_data == nullptr)
		throw logic_error("invalid row, no data");
	
	auto cat = m_data->m_category;
	auto cix = cat->add_column(name);
	auto& col = cat->m_columns[cix];
//	auto& db = cat->m_db;

	const char* oldValue = nullptr;
	for (auto iv = m_data->m_values; iv != nullptr; iv = iv->m_next)
	{
		assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));

		if (iv->m_column_index == cix)
		{
			oldValue = iv->m_text;
			break;
		}
	}
	
	if (oldValue != nullptr and value == oldValue)	// no need to update
		return;

	// check the value
	if (col.m_validator)
		(*col.m_validator)(value);

	// If the field is part of the key for this category, remove it from the index
	// before updating
	
	bool reinsert = false;
	
	if (not emplacing)	// an update of an item's value
	{
////#if DEBUG
////		if (VERBOSE)
////			cerr << "reassigning the value of key field _" << cat->m_name << '.' << name << endl;
////#endif
//		// see if we need to update any child categories that depend on this value
//		auto iv = col.m_validator;
//		if (iv != nullptr and not iv->m_children.empty())
//		{
//			for (auto child: iv->m_children)
//			{
//				if (child->m_category == nullptr)
//					continue;
//				
//				auto child_cat = db.get(child->m_category->m_name);
//				if (child_cat == nullptr)
//					continue;
//					
//				auto rows = child_cat->find(key(child->m_tag) == oldValue);
//				for (auto& cr: rows)
//					cr.assign(child->m_tag, value, false);
//			}
//		}

		if (cat->m_index != nullptr and cat->key_fields().count(name))
		{
			reinsert = cat->m_index->find(m_data);
			if (reinsert)
				cat->m_index->erase(m_data);
		}
	}

	// first remove old value with cix

	if (m_data->m_values == nullptr)
		;	// nothing to do
	else if (m_data->m_values->m_column_index == cix)
	{
		auto iv = m_data->m_values;
		m_data->m_values = iv->m_next;
		iv->m_next = nullptr;
		delete iv;
	}
	else
	{
		for (auto iv = m_data->m_values; iv->m_next != nullptr; iv = iv->m_next)
		{
			if (iv->m_next->m_column_index == cix)
			{
				auto nv = iv->m_next;
				iv->m_next = nv->m_next;
				nv->m_next = nullptr;
				delete nv;
				
				break;
			}
		}
	}

#if DEBUG
	for (auto iv = m_data->m_values; iv != nullptr; iv = iv->m_next)
		assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));
#endif

	if (not value.empty())
	{
		auto nv = new(value.length()) item_value(value.c_str(), cix);
	
		if (m_data->m_values == nullptr)
			m_data->m_values = nv;
		else
		{
			auto iv = m_data->m_values;
			while (iv->m_next != nullptr)
				iv = iv->m_next;
			iv->m_next = nv;
		}
	}

#if DEBUG
	for (auto iv = m_data->m_values; iv != nullptr; iv = iv->m_next)
		assert(iv != iv->m_next and (iv->m_next == nullptr or iv != iv->m_next->m_next));
#endif

	if (reinsert)
		cat->m_index->insert(m_data);
}

void row::assign(const item& value, bool emplacing)
{
	assign(value.name(), value.value(), emplacing);
}

bool row::empty() const
{
	return m_data == nullptr or m_data->m_values == nullptr;
}

auto row::begin() const -> const_iterator
{
	return const_iterator(m_data, m_data->m_values);
}

auto row::end() const -> const_iterator
{
	return const_iterator(m_data, nullptr);
}

row::const_iterator::const_iterator(item_row* data, item_value* ptr)
	: m_data(data), m_ptr(ptr)
{
	if (m_ptr != nullptr)
		fetch();
}

row::const_iterator& row::const_iterator::operator++()
{
	if (m_ptr != nullptr)
		m_ptr = m_ptr->m_next;

	if (m_ptr != nullptr)
		fetch();
	
	return *this;
}

void row::const_iterator::fetch()
{
	m_current = item(
		m_data->m_category->get_column_name(m_ptr->m_column_index),
		m_ptr->m_text);
}

// --------------------------------------------------------------------

file::file()
	: m_head(nullptr)
	, m_validator(nullptr)
{
}

file::file(istream& is, bool validate)
	: file()
{
//	parser p(is, *this);
//	p.parse_file();
	load(is);
}

file::file(file&& rhs)
	: m_head(nullptr), m_validator(nullptr)
{
	swap(m_head, rhs.m_head);
	swap(m_validator, rhs.m_validator);
}

file::~file()
{
	delete m_head;
	delete m_validator;
}

void file::append(datablock* e)
{
	e->set_validator(m_validator);
	
	if (m_head == nullptr)
		m_head = e;
	else
	{
		auto ie = m_head;
		for (;;)
		{
			if (iequals(ie->name(), e->name()))
				throw validation_error("datablock " + e->name() + " already defined in file");

			if (ie->m_next == nullptr)
			{
				ie->m_next = e;
				break;
			}
			
			ie = ie->m_next;
		}
	}
}

void file::load(istream& is)
{
	validator* saved = m_validator;
	set_validator(nullptr);

	parser p(is, *this);
	p.parse_file();
	
	if (saved != nullptr)
	{
		set_validator(saved);
		validate();
	}
}

void file::save(ostream& os)
{
	datablock* e = m_head;
	while (e != nullptr)
	{
		e->write(os);
		e = e->m_next;
	}
}

void file::write(ostream& os, const vector<string>& order)
{
	datablock* e = m_head;
	while (e != nullptr)
	{
		e->write(os, order);
		e = e->m_next;
	}
}

datablock& file::operator[](const string& name)
{
	datablock* result = m_head;
	while (result != nullptr and not iequals(result->m_name, name))
		result = result->m_next;
	if (result == nullptr)
		throw runtime_error("datablock " + name + " does not exist");
	return *result;
}

void file::validate()
{
	if (m_validator == nullptr)
	{
		if (VERBOSE)
			cerr << "No dictionary loaded explicitly, loading default" << endl;
		
		load_dictionary();
	}

	for (auto d = m_head; d != nullptr; d = d->m_next)
		d->validate();
}

const validator& file::get_validator() const
{
	if (m_validator == nullptr)
		throw runtime_error("no validator defined yet");
	return *m_validator;
}

void file::load_dictionary()
{
	load_dictionary("mmcif_ddl");
}

void file::load_dictionary(const char* dict)
{
	fs::path dict_file = string("dictionaries/") + dict + ".dic";
	
#if defined(USE_RSRC)
	mrsrc::rsrc dict_data(dict_file.string());

	if (not dict_data)
		throw invalid_argument("no such dictionary");
	
	struct membuf : public streambuf
	{
		membuf(char* dict, size_t length)
		{
			this->setg(dict, dict, dict + length);
		}
	} buffer(const_cast<char*>(dict_data.data()), dict_data.size());
	
	istream is(&buffer);
#else
	if (not fs::exists(dict_file))
		throw runtime_error("Dictionary not found (" + dict_file.string() + ")");
	fs::ifstream is(dict_file);
#endif

	load_dictionary(is);
}

void file::load_dictionary(istream& is)
{
	unique_ptr<validator> v(new validator());

	dict_parser p(*v, is);
	p.load_dictionary();

	set_validator(v.release());
}

void file::set_validator(validator* v)
{
	m_validator = v;

	for (auto d = m_head; d != nullptr; d = d->m_next)
		d->set_validator(m_validator);
}

void file::get_tag_order(vector<string>& tags) const
{
	for (auto d = m_head; d != nullptr; d = d->m_next)
		d->get_tag_order(tags);
}

auto file::iterator::operator++() -> iterator&
{
	m_current = m_current->m_next;
	return *this;
}

auto file::begin() const -> iterator
{
	return iterator(m_head);
}

auto file::end() const -> iterator
{
	return iterator(nullptr);
}

}
