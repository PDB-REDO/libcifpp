/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/cif++.hpp"
#include "cif++/item.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <exception>
#include <format>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ranges>
#include <regex>
#include <sqlite3.h>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

namespace cif::cql
{

struct result_impl
{
	category m_cat;
	std::string m_query;
};

// --------------------------------------------------------------------

size_t row_ref::size() const noexcept
{
	return m_result_impl->m_cat.get_item_count();
}

field_ref row_ref::operator[](std::string_view name) const
{
	for (uint16_t ix = 0; auto &item : m_result_impl->m_cat.get_items())
	{
		if (iequals(item, name))
			return { m_row, ix, m_result_impl };
		++ix;
	}
	throw std::runtime_error("Column not defined in query result");
}

// --------------------------------------------------------------------

result::result(category &&cat, const std::string &query)
	: m_impl(new result_impl{ .m_cat = std::forward<category>(cat), .m_query = query })
{
}

category &result::get_category() const
{
	return m_impl->m_cat;
}

size_t result::size() const noexcept
{
	return m_impl->m_cat.size();
}

size_t result::column_count() const
{
	return m_impl->m_cat.get_item_count();
}

result::iterator result::begin() const noexcept
{
	return { m_impl, m_impl->m_cat.begin() };
}

result::iterator result::cbegin() const noexcept
{
	return { m_impl, m_impl->m_cat.cbegin() };
}

result::iterator result::end() const noexcept
{
	return { m_impl, m_impl->m_cat.end() };
}

result::iterator result::cend() const noexcept
{
	return { m_impl, m_impl->m_cat.cend() };
}

row_ref result::front() const
{
	return { m_impl->m_cat.front(), m_impl };
}

row_ref result::back() const
{
	return { m_impl->m_cat.back(), m_impl };
}

// --------------------------------------------------------------------

struct virtual_table;

struct connection_impl
{
	datablock &m_db;
	sqlite3 *m_sqlite_db = nullptr;
	int m_next_result_nr = 1;
	std::vector<virtual_table *> m_vtabs;
	bool m_db_modified = false;

	connection_impl(datablock &db);

	~connection_impl()
	{
		sqlite3_close(m_sqlite_db);
	}

	int Connect(sqlite3 *db, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr);

	// The module interface

	static int Create(sqlite3 *db, void *pAux, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr);
	static int Connect(sqlite3 *db, void *pAux, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr);
	static int Destroy(sqlite3_vtab *pVtab);
	static int Disconnect(sqlite3_vtab *pVtab);
	static int Open(sqlite3_vtab *p, sqlite3_vtab_cursor **ppCursor);
	static int Close(sqlite3_vtab_cursor *cur);
	static int Next(sqlite3_vtab_cursor *cur);
	static int Column(sqlite3_vtab_cursor *cur, sqlite3_context *ctx, int i);
	static int Rowid(sqlite3_vtab_cursor *cur, sqlite_int64 *pRowid);
	static int Eof(sqlite3_vtab_cursor *cur);
	static int Filter(sqlite3_vtab_cursor *pVtabCursor, int idxNum, const char *idxStr, int argc, sqlite3_value **argv);
	static int BestIndex(sqlite3_vtab *tab, sqlite3_index_info *pIdxInfo);

	static int Update(sqlite3_vtab *pVTab, int argc, sqlite3_value **argv, sqlite_int64 *pRowid);

	static int Rename(sqlite3_vtab *pVtab, const char *zNew);

	// Transaction support
	static int Begin(sqlite3_vtab *pVTab);
	static int Commit(sqlite3_vtab *pVTab);
	static int Rollback(sqlite3_vtab *pVTab);

	static sqlite3_module s_module;
};

struct virtual_table
{
	sqlite3_vtab base;
	connection_impl &m_connection_impl;
	category &m_cat;
	datablock &m_db;
	std::stack<category> m_rollback_buffer;
	std::vector<std::string> m_items;
};

struct virtual_cursor
{
	sqlite3_vtab_cursor base;

	std::unique_ptr<conditional_iterator_proxy<>> m_result;
	conditional_iterator_proxy<>::iterator m_cur;
};

sqlite3_module connection_impl::s_module{
	/* iVersion    */ 0,
	/* xCreate     */ Create,
	/* xConnect    */ Connect,
	/* xBestIndex  */ BestIndex,
	/* xDisconnect */ Disconnect,
	/* xDestroy    */ Destroy,
	/* xOpen       */ Open,
	/* xClose      */ Close,
	/* xFilter     */ Filter,
	/* xNext       */ Next,
	/* xEof        */ Eof,
	/* xColumn     */ Column,
	/* xRowid      */ Rowid,
	/* xUpdate     */ Update,
	/* xBegin      */ Begin,
	/* xSync       */ nullptr,
	/* xCommit     */ Commit,
	/* xRollback   */ Rollback,
	/* xFindFunction */ nullptr,
	/* xRename     */ Rename,
	// /* xSavepoint  */ nullptr,
	// /* xRelease    */ nullptr,
	// /* xRollbackTo */ nullptr,
	// /* xShadowName */ nullptr,
	// /* xIntegrity  */ nullptr
};

/*
** The templatevtabConnect() method is invoked to create a new
** template virtual table.
**
** Think of this routine as the constructor for connection_impl objects.
**
** All this routine needs to do is:
**
**    (1) Allocate the connection_impl object and initialize all fields.
**
**    (2) Tell SQLite (via the sqlite3_declare_vtab() interface) what the
**        result set of queries against the virtual table will look like.
*/

int connection_impl::Connect(sqlite3 *db, void *pAux, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr)
{
	auto *impl = reinterpret_cast<connection_impl *>(pAux);
	try
	{
		return impl->Connect(db, argc, argv, ppVtab, pzErr);
	}
	catch (const std::exception &ex)
	{
		*pzErr = sqlite3_mprintf("%s", ex.what());
		return SQLITE_ERROR;
	}
}

int connection_impl::Create(sqlite3 *db, void *pAux, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr)
{
	return Connect(db, pAux, argc, argv, ppVtab, pzErr);
}

// --------------------------------------------------------------------

int connection_impl::Connect(sqlite3 *db, int argc, const char *const *argv, sqlite3_vtab **ppVtab, char **pzErr)
{
	if (argc < 3)
		throw std::runtime_error("Insufficient arguments to module connect");

	auto cat = m_db.get(argv[2]);
	if (cat == nullptr)
		throw std::runtime_error(std::format("Category {} is not known in this databank", argv[2]));

	auto vtab = std::make_unique<virtual_table>(sqlite3_vtab{}, *this, *cat, m_db);
	m_vtabs.emplace_back(vtab.get());

	std::vector<std::string> columns;

	if (auto cv = cat->get_cat_validator(); cv != nullptr)
	{
		// Try to keep the order the same as in the parsed category

		auto &items = vtab->m_items;

		for (auto item : cat->get_items())
			items.emplace_back(item);

		for (auto iv : cv->m_item_validators)
		{
			if (std::ranges::find_if(items, [&b = iv.m_item_name](const std::string &a)
					{ return iequals(a, b); }) == items.end())
				items.emplace_back(iv.m_item_name);
		}

		// Make sure all items are known in the category
		for (auto item : items)
		{
			auto iv = cv->get_validator_for_item(item);

			std::string primaryKey;
			if (cv->m_keys.size() == 1 and cv->m_keys.front() == item)
				primaryKey = " PRIMARY KEY";

			if (iv != nullptr and iv->m_type->m_primitive_type == DDL_PrimitiveType::Numb)
			{
				if (iequals(iv->m_type->m_name, "int"))
				{
					columns.emplace_back(std::format("'{}' {}", item, " INTEGER") + primaryKey);
					continue;
				}

				if (iequals(iv->m_type->m_name, "float"))
				{
					columns.emplace_back(std::format("'{}' {}", item, " REAL") + primaryKey);
					continue;
				}
			}

			columns.emplace_back(std::format("'{}' {}", item, " TEXT") + primaryKey);
		}
	}
	else
	{
		for (auto item : cat->get_items())
		{
			vtab->m_items.emplace_back(item);
			columns.emplace_back(std::format("'{}'", item));
		}
	}

	auto createStmt = std::format("CREATE TABLE {} ({})", cat->name(), join(columns, ", "));

	int rc = sqlite3_declare_vtab(db, createStmt.c_str());
	if (rc == SQLITE_OK)
		*ppVtab = reinterpret_cast<sqlite3_vtab *>(vtab.release());
	else
		std::clog << "statement:\n"
				  << createStmt << "\nresulted in error: " << sqlite3_errmsg(db) << '\n';

	return rc;
}

/*
** This method is used to drop tables.
*/
int connection_impl::Destroy(sqlite3_vtab *pVtab)
{
	auto *p = reinterpret_cast<virtual_table *>(pVtab);

	auto &conn = p->m_connection_impl;
	auto e = std::ranges::remove(conn.m_vtabs, p);
	conn.m_vtabs.erase(e.begin(), e.end());

	auto &db = p->m_db;
	db.remove(p->m_cat);
	delete p;

	return SQLITE_OK;
}

/*
** This method is the destructor for connection_impl objects.
*/
int connection_impl::Disconnect(sqlite3_vtab *pVtab)
{
	auto *p = reinterpret_cast<virtual_table *>(pVtab);

	auto &conn = p->m_connection_impl;
	auto e = std::ranges::remove(conn.m_vtabs, p);
	conn.m_vtabs.erase(e.begin(), e.end());

	delete p;

	return SQLITE_OK;
}

/*
** Constructor for a new templatevtab_cursor object.
*/
int connection_impl::Open(sqlite3_vtab *pVtab, sqlite3_vtab_cursor **ppCursor)
{
	auto cursor = std::make_unique<virtual_cursor>(sqlite3_vtab_cursor{});
	*ppCursor = reinterpret_cast<sqlite3_vtab_cursor *>(cursor.release());
	return SQLITE_OK;
}

/*
** Destructor for a templatevtab_cursor.
*/
int connection_impl::Close(sqlite3_vtab_cursor *cur)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(cur);
	delete pCur;
	return SQLITE_OK;
}

/*
** Advance a templatevtab_cursor to its next row of output.
*/
int connection_impl::Next(sqlite3_vtab_cursor *cur)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(cur);
	++pCur->m_cur;
	return SQLITE_OK;
}

/*
** Return values of columns for the row at which the templatevtab_cursor
** is currently pointing.
*/
int connection_impl::Column(sqlite3_vtab_cursor *cur, sqlite3_context *ctx, int i)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(cur);
	auto pVTab = reinterpret_cast<virtual_table *>(cur->pVtab);
	auto rh = *pCur->m_cur;
	auto &cat = pVTab->m_cat;
	auto ix = cat.get_item_ix(pVTab->m_items[i]);

	if (ix >= cat.get_item_count() or rh[ix].empty())
		sqlite3_result_null(ctx);
	else
	{
		auto item = rh[ix];

		switch (item.type())
		{
			using enum item_value_type;

			case FLOAT:
				sqlite3_result_double(ctx, item.get<double>());
				break;

			case INT:
				sqlite3_result_int64(ctx, item.get<int64_t>());
				break;

			case TEXT:
				sqlite3_result_text(ctx, item.sv().data(), static_cast<int>(item.sv().size()), SQLITE_STATIC);
				break;

			default:
				sqlite3_result_null(ctx);
				break;
		}
	}

	return SQLITE_OK;
}

/*
** Return the rowid for the current row.  In this implementation, the
** rowid is the same as the output value.
*/
int connection_impl::Rowid(sqlite3_vtab_cursor *cur, sqlite_int64 *pRowid)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(cur);
	row_handle rh = *pCur->m_cur;
	*pRowid = rh.row_id();
	return SQLITE_OK;
}

/*
** Return TRUE if the cursor has been moved off of the last
** row of output.
*/
int connection_impl::Eof(sqlite3_vtab_cursor *cur)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(cur);
	auto pVTab = reinterpret_cast<virtual_table *>(cur->pVtab);
	return pCur->m_cur == pVTab->m_cat.end();
}

/*
** This method is called to "rewind" the templatevtab_cursor object back
** to the first row of output.  This method is always called at least
** once prior to any call to templatevtabColumn() or templatevtabRowid() or
** templatevtabEof().
*/
int connection_impl::Filter(sqlite3_vtab_cursor *pVtabCursor, int idxNum, const char *idxStr, int argc, sqlite3_value **argv)
{
	auto pCur = reinterpret_cast<virtual_cursor *>(pVtabCursor);
	auto pVTab = reinterpret_cast<virtual_table *>(pCur->base.pVtab);
	auto &cat = pVTab->m_cat;

	pCur->m_result.reset();

	try
	{
		if (idxStr != nullptr)
		{
			struct membuf : public std::streambuf
			{
				membuf(char *text, std::size_t length)
				{
					this->setg(text, text, text + length);
				}
			} buffer(const_cast<char *>(idxStr), strlen(idxStr));

			std::istream is(&buffer);

			std::regex rx("^(.+?)( IS NULL| IS NOT NULL|(?: < | <= | == | >= | > ))(.+)?$");

			condition cond;
			std::string line;
			while (std::getline(is, line))
			{
				std::smatch m;
				if (not std::regex_match(line, m, rx))
					throw std::runtime_error("Internal error in cql, no match");

				if (m[2] == " IS NULL")
					cond = std::move(cond) and key(m[1]) == null;
				else if (m[2] == " IS NOT NULL")
					cond = std::move(cond) and key(m[1]) != null;
				else if (m[3].str().starts_with("\""))
				{
					std::istringstream isv(m[3]);
					std::string value;
					isv >> std::quoted(value);

					if (m[2] == " < ")
						cond = std::move(cond) and key(m[1]) < value;
					else if (m[2] == " <= ")
						cond = std::move(cond) and key(m[1]) <= value;
					else if (m[2] == " == ")
						cond = std::move(cond) and key(m[1]) == value;
					else if (m[2] == " >= ")
						cond = std::move(cond) and key(m[1]) >= value;
					else if (m[2] == " > ")
						cond = std::move(cond) and key(m[1]) > value;
				}
				else
				{
					double value;
					const auto &[ptr, ec] = from_chars(m[3].str().data(), m[3].str().data() + m[3].str().length(), value);
					if (ec != std::errc{})
						throw std::system_error(std::make_error_code(ec));

					if (m[2] == " < ")
						cond = std::move(cond) and key(m[1]) < value;
					else if (m[2] == " <= ")
						cond = std::move(cond) and key(m[1]) <= value;
					else if (m[2] == " == ")
						cond = std::move(cond) and key(m[1]) == value;
					else if (m[2] == " >= ")
						cond = std::move(cond) and key(m[1]) >= value;
					else if (m[2] == " > ")
						cond = std::move(cond) and key(m[1]) > value;
				}
			}

			pCur->m_result = std::make_unique<conditional_iterator_proxy<>>(cat.find(std::move(cond)));
			pCur->m_cur = pCur->m_result->begin();
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Internal error: " << ex.what() << "\n";
		pVtabCursor->pVtab->zErrMsg = sqlite3_mprintf("%s", ex.what());
	}

	if (not pCur->m_result)
	{
		condition cond = all();
		pCur->m_result = std::make_unique<conditional_iterator_proxy<>>(cat.find(std::move(cond)));
		pCur->m_cur = pCur->m_result->begin();
	}

	return SQLITE_OK;
}

/*
** SQLite will invoke this method one or more times while planning a query
** that uses the virtual table.  This routine needs to create
** a query plan for each invocation and compute an estimated cost for that
** plan.
*/
int connection_impl::BestIndex(sqlite3_vtab *pVtab, sqlite3_index_info *pIdxInfo)
{
	auto *p = reinterpret_cast<virtual_table *>(pVtab);

	try
	{
		std::ostringstream os;
		bool ok = true;

		if (pIdxInfo->nConstraint > 0)
		{
			auto constraint = [&os](std::string_view item, sqlite3_value *val, unsigned char op)
			{
				bool result = true;
				switch (op)
				{
					case SQLITE_INDEX_CONSTRAINT_EQ:
						os << item << " == ";
						break;
					case SQLITE_INDEX_CONSTRAINT_GT:
						os << item << " > ";
						break;
					case SQLITE_INDEX_CONSTRAINT_LE:
						os << item << " <= ";
						break;
					case SQLITE_INDEX_CONSTRAINT_LT:
						os << item << " < ";
						break;
					case SQLITE_INDEX_CONSTRAINT_GE:
						os << item << " >= ";
						break;
					default:
						result = false;
						break;
				}

				if (result)
				{
					switch (sqlite3_value_type(val))
					{
						case SQLITE_INTEGER:
							os << sqlite3_value_int64(val) << "\n";
							break;
						case SQLITE_FLOAT:
							os << sqlite3_value_double(val) << "\n";
							break;
						default:
						{
							std::string s = reinterpret_cast<const char *>(sqlite3_value_text(val));
							if (s.find("\n") == std::string::npos)
								os << std::quoted(s) << "\n";
							else
								result = false;
							break;
						}
					}
				}

				return result;
			};

			for (int i = 0; ok and i < pIdxInfo->nConstraint; ++i)
			{
				auto &info = pIdxInfo->aConstraint[i];
				auto item = p->m_items[info.iColumn];

				sqlite3_value *pVal;

				switch (info.op)
				{
					case SQLITE_INDEX_CONSTRAINT_EQ:
					case SQLITE_INDEX_CONSTRAINT_GT:
					case SQLITE_INDEX_CONSTRAINT_LE:
					case SQLITE_INDEX_CONSTRAINT_LT:
					case SQLITE_INDEX_CONSTRAINT_GE:
						if (sqlite3_vtab_rhs_value(pIdxInfo, i, &pVal) == SQLITE_OK and constraint(item, pVal, info.op))
						{
							pIdxInfo->aConstraintUsage[i].omit = 1;
							if (i < 63 and sqlite3_libversion_number() >= 3010000)
								pIdxInfo->colUsed |= 1ULL << i;
						}
						else
							ok = false;
						break;
					// case SQLITE_INDEX_CONSTRAINT_MATCH:
					// 	break;
					// case SQLITE_INDEX_CONSTRAINT_LIKE:
					// 	break;
					// case SQLITE_INDEX_CONSTRAINT_GLOB:
					// 	break;
					// case SQLITE_INDEX_CONSTRAINT_REGEXP:
					// 	break;
					case SQLITE_INDEX_CONSTRAINT_NE:
						break;
					// case SQLITE_INDEX_CONSTRAINT_ISNOT:
					// 	break;
					case SQLITE_INDEX_CONSTRAINT_ISNOTNULL:
						os << item << " IS NOT NULL\n";
						pIdxInfo->aConstraintUsage[i].omit = 1;
						if (i < 63 and sqlite3_libversion_number() >= 3010000)
							pIdxInfo->colUsed |= 1ULL << i;
						break;
					case SQLITE_INDEX_CONSTRAINT_ISNULL:
						os << item << " IS NULL\n";
						pIdxInfo->aConstraintUsage[i].omit = 1;
						if (i < 63 and sqlite3_libversion_number() >= 3010000)
							pIdxInfo->colUsed |= 1ULL << i;
						break;
						// case SQLITE_INDEX_CONSTRAINT_IS:
						// 	break;
						// case SQLITE_INDEX_CONSTRAINT_LIMIT:
						// 	break;
						// case SQLITE_INDEX_CONSTRAINT_OFFSET:
						// 	break;
						// case SQLITE_INDEX_CONSTRAINT_FUNCTION				:
						// 	break;

					default:
						ok = false;
						break;
				}
			}
		}

		if (auto cs = os.str(); ok and not cs.empty())
		{
			pIdxInfo->idxStr = sqlite3_mprintf("%s", cs.c_str());
			pIdxInfo->needToFreeIdxStr = 1;
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << ex.what() << "\n";
		if (sqlite3_libversion_number() >= 3010000)
			pIdxInfo->colUsed = 0;

		pVtab->zErrMsg = sqlite3_mprintf("%s", ex.what());
	}

	pIdxInfo->estimatedCost = static_cast<double>(p->m_cat.size());
	pIdxInfo->estimatedRows = static_cast<int64_t>(p->m_cat.size());
	return SQLITE_OK;
}

int bind_item_value(sqlite3_stmt *stmt, int ix, const item_value &value)
{
	switch (value.type())
	{
		using enum item_value_type;

		case FLOAT: return sqlite3_bind_double(stmt, ix, value.get<double>());
		case INT: return sqlite3_bind_int64(stmt, ix, value.get<int64_t>());
		case TEXT: return sqlite3_bind_text(stmt, ix, value.sv().data(), static_cast<int>(value.sv().size()), SQLITE_STATIC);
		default: return sqlite3_bind_null(stmt, ix);
	}
}

int connection_impl::Update(sqlite3_vtab *pVTab, int argc, sqlite3_value **argv, sqlite_int64 *pRowid)
{
	auto *p = reinterpret_cast<virtual_table *>(pVTab);

	int rc = SQLITE_ERROR;

	try
	{
		auto addr = sqlite3_value_int64(argv[0]);

		if (argc == 1) // DELETE
		{
			rc = SQLITE_OK;
			row_handle rh{ p->m_cat, *reinterpret_cast<row *>(addr) };

			if (auto v = p->m_cat.get_validator())
			{
				auto &db = p->m_db;

				for (auto link : v->get_links_for_parent(p->m_cat.name()))
				{
					auto childCat = const_cast<category *>(db.get(link->m_child_category));
					if (childCat == nullptr)
						continue;

					std::ostringstream sql;
					sql << "DELETE FROM " << link->m_child_category << " WHERE ";
					for (size_t i = 0; i < link->m_child_keys.size(); ++i)
					{
						if (i > 0)
							sql << " AND ";
						sql << link->m_child_keys[i] << " = ?" << (i + 1);
					}

					sqlite3_stmt *sub_stmt;
					auto sqls = sql.str();
					rc = sqlite3_prepare_v2(p->m_connection_impl.m_sqlite_db, sqls.c_str(), static_cast<int>(sqls.length()), &sub_stmt, nullptr);

					for (int i = 0; rc == SQLITE_OK and static_cast<size_t>(i) < link->m_parent_keys.size(); ++i)
						rc = bind_item_value(sub_stmt, i + 1, rh[link->m_parent_keys[i]].value());

					if (rc == SQLITE_OK)
						rc = sqlite3_step(sub_stmt);
					(void)sqlite3_finalize(sub_stmt);

					if (rc != SQLITE_DONE)
						break;

					rc = SQLITE_OK;
				}
			}

			if (rc == SQLITE_OK)
				p->m_cat.erase(rh);
		}
		else if (addr == 0) // INSERT
		{
			addr = sqlite3_value_int64(argv[1]);
			if (addr == 0) // We do not accept rowid's here
			{
				row_initializer data;
				for (int i = 2; i < argc; ++i)
				{
					switch (sqlite3_value_type(argv[i]))
					{
						case SQLITE_INTEGER:
							data.emplace_back(p->m_items[i - 2], sqlite3_value_int64(argv[i]));
							break;
						case SQLITE_FLOAT:
							data.emplace_back(p->m_items[i - 2], sqlite3_value_double(argv[i]));
							break;
						case SQLITE_NULL:
							data.emplace_back(p->m_items[i - 2], item_value_type::MISSING);
							break;
						default:
							data.emplace_back(p->m_items[i - 2], reinterpret_cast<const char *>(sqlite3_value_text(argv[i])));
							break;
					}
				}

				auto r = p->m_cat.emplace(std::move(data));
				*pRowid = r->row_id();
				rc = SQLITE_OK;
			}
		}
		else // UPDATE
		{
			row_handle rh{ p->m_cat, *reinterpret_cast<row *>(addr) };

			row_initializer data;
			for (int i = 2; i < argc; ++i)
			{
				switch (sqlite3_value_type(argv[i]))
				{
					case SQLITE_INTEGER:
						data.emplace_back(p->m_items[i - 2], sqlite3_value_int64(argv[i]));
						break;
					case SQLITE_FLOAT:
						data.emplace_back(p->m_items[i - 2], sqlite3_value_double(argv[i]));
						break;
					case SQLITE_NULL:
						data.emplace_back(p->m_items[i - 2], item_value_type::MISSING);
						break;
					default:
						data.emplace_back(p->m_items[i - 2], reinterpret_cast<const char *>(sqlite3_value_text(argv[i])));
						break;
				}
			}

			if (auto v = p->m_cat.get_validator())
			{
				auto &db = p->m_db;

				for (auto link : v->get_links_for_parent(p->m_cat.name()))
				{
					auto childCat = const_cast<category *>(db.get(link->m_child_category));
					if (childCat == nullptr)
						continue;

					std::vector<std::tuple<int, item_value>> ixs;
					for (auto &ri : data)
					{
						auto i = std::ranges::find(link->m_parent_keys, ri.name());
						if (i == link->m_parent_keys.end())
							continue; // no update needed
						ixs.emplace_back(i - link->m_parent_keys.begin(), ri.value());
					}

					if (ixs.empty())
						continue;

					std::ostringstream sql;
					sql << "UPDATE " << link->m_child_category;

					for (bool first = true; auto [i, txt] : ixs)
					{
						if (not std::exchange(first, false))
							sql << ",";
						sql << " SET " << link->m_child_keys[i] << " = ?" << (i + 1);
					}

					sql << " WHERE ";
					for (bool first = true; auto [i, txt] : ixs)
					{
						if (not std::exchange(first, false))
							sql << " AND ";
						sql << link->m_child_keys[i] << " = ?" << (i + ixs.size() + 1);
					}

					sqlite3_stmt *sub_stmt;
					auto sqls = sql.str();
					rc = sqlite3_prepare_v2(p->m_connection_impl.m_sqlite_db, sqls.c_str(), static_cast<int>(sqls.length()), &sub_stmt, nullptr);

					for (auto [i, txt] : ixs)
					{
						// set
						if (rc == SQLITE_OK)
							rc = bind_item_value(sub_stmt, i + 1, txt);

						// where
						rc = bind_item_value(sub_stmt, static_cast<int>(i + ixs.size() + 1), rh[link->m_parent_keys[i]].value());
					}

					if (rc == SQLITE_OK)
						rc = sqlite3_step(sub_stmt);
					(void)sqlite3_finalize(sub_stmt);

					if (rc != SQLITE_DONE)
						break;

					rc = SQLITE_OK;
				}
			}

			rh.assign(data);
			*pRowid = addr;
			rc = SQLITE_OK;
		}
	}
	catch (const std::exception &ex)
	{
		rc = SQLITE_ERROR;
		pVTab->zErrMsg = sqlite3_mprintf("%s", ex.what());
	}

	return rc;
}

int connection_impl::Rename(sqlite3_vtab *pVtab, const char *zNew)
{
	auto *p = reinterpret_cast<virtual_table *>(pVtab);
	p->m_cat.name(zNew);
	return SQLITE_OK;
}

int connection_impl::Begin(sqlite3_vtab *pVTab)
{
	auto *p = reinterpret_cast<virtual_table *>(pVTab);
	p->m_rollback_buffer.push(p->m_cat);
	return SQLITE_OK;
}

int connection_impl::Commit(sqlite3_vtab *pVTab)
{
	auto *p = reinterpret_cast<virtual_table *>(pVTab);
	if (not p->m_rollback_buffer.empty())
		p->m_rollback_buffer.pop();
	return SQLITE_OK;
}

int connection_impl::Rollback(sqlite3_vtab *pVTab)
{
	auto *p = reinterpret_cast<virtual_table *>(pVTab);
	if (not p->m_rollback_buffer.empty())
	{
		std::swap(p->m_cat, p->m_rollback_buffer.top());
		p->m_rollback_buffer.pop();
	}
	return SQLITE_OK;
}

// --------------------------------------------------------------------

connection_impl::connection_impl(datablock &db)
	: m_db(db)
{
	auto rc = sqlite3_open(":memory:", &m_sqlite_db);

	if (rc)
		throw std::runtime_error(std::format("Cannot open databank: {}", sqlite3_errmsg(m_sqlite_db)));

	rc = sqlite3_create_module_v2(m_sqlite_db, "CIFPP", &connection_impl::s_module, this, nullptr);

	if (rc)
		throw std::runtime_error(std::format("Cannot create module: {}", sqlite3_errmsg(m_sqlite_db)));

	// Now, create a table for all known categories in the datablock

	for (auto &cat : db)
	{
		char *errmsg;
		rc = sqlite3_exec(m_sqlite_db,
			("CREATE VIRTUAL TABLE " + cat.name() + " USING CIFPP;").c_str(),
			nullptr, nullptr, &errmsg);
		if (rc != SQLITE_OK)
		{
			if (errmsg != nullptr)
			{
				std::string err = errmsg;
				sqlite3_free(errmsg);

				throw std::runtime_error("Error creating virtual tables for the categories: " + err);
			}

			throw std::runtime_error("Error creating virtual tables for the categories");
		}
	}
}

connection::connection(datablock &db)
	: m_impl(new connection_impl(db))
{
}

connection::~connection()
{
	delete m_impl;
}

bool connection::is_complete_statement(const std::string &sql) const
{
	return sqlite3_complete(sql.c_str()) == 1;
}

bool connection::is_modified() const
{
	return m_impl->m_db_modified;
}

result connection::exec(std::string query)
{
	trim(query);

	for (;;)
	{
		auto r = exec(query, query);
		if (query.empty())
			return r;
	}
}

result connection::exec(std::string query, std::string &tail)
{
	category cat;

	sqlite3_stmt *stmt = nullptr;

	try
	{
		const char *tail_ptr = nullptr;
		int rc = sqlite3_prepare_v2(m_impl->m_sqlite_db, query.data(), static_cast<int>(query.length()),
			&stmt, &tail_ptr);

		if (rc != SQLITE_OK)
			throw std::runtime_error(std::format("Error preparing statement: {}", sqlite3_errmsg(m_impl->m_sqlite_db)));

		if (not sqlite3_stmt_readonly(stmt))
			m_impl->m_db_modified = true;

		auto used = tail_ptr - query.data();
		auto sql = query.substr(0, used);
		tail = trim_copy(query.substr(used));

		category ncat((std::format("result{}", m_impl->m_next_result_nr++)));
		std::swap(cat, ncat);

		for (;;)
		{
			rc = sqlite3_step(stmt);
			if (rc == SQLITE_ROW)
			{
				row_initializer data;

				for (int i = 0; i < sqlite3_column_count(stmt); ++i)
				{
					switch (sqlite3_column_type(stmt, i))
					{
						case SQLITE_INTEGER:
							data.emplace_back(sqlite3_column_name(stmt, i), sqlite3_column_int64(stmt, i));
							break;
						case SQLITE_FLOAT:
							data.emplace_back(sqlite3_column_name(stmt, i), sqlite3_column_double(stmt, i));
							break;
						case SQLITE_TEXT:
							data.emplace_back(sqlite3_column_name(stmt, i), reinterpret_cast<const char *>(sqlite3_column_text(stmt, i)));
							break;
						case SQLITE_BLOB:
							throw std::runtime_error("Unexpected: blob in result");
							break;
						case SQLITE_NULL:
						default:
							data.emplace_back(sqlite3_column_name(stmt, i), ".");
							break;
					}
				}

				cat.emplace(std::move(data));
				continue;
			}

			if (rc == SQLITE_BUSY)
				throw std::runtime_error("Oops, busy?");
			if (rc == SQLITE_DONE)
				break;
			if (rc == SQLITE_ERROR)
				throw std::runtime_error(std::format("Error in sqlite: {}", sqlite3_errmsg(m_impl->m_sqlite_db)));

			throw std::runtime_error("Unknown result from step");
		}

		sqlite3_finalize(stmt);

		return { std::move(cat), sql };
	}
	catch (const std::exception &ex)
	{
		if (stmt)
			sqlite3_finalize(stmt);
		throw;
	}
}

// --------------------------------------------------------------------

transaction::transaction(connection &conn)
	: m_conn(conn)
{
	char *errmsg = nullptr;
	std::string err;
	int rc = sqlite3_exec(m_conn.m_impl->m_sqlite_db, "BEGIN TRANSACTION;", nullptr, nullptr, &errmsg);
	if (errmsg)
	{
		err = errmsg;
		sqlite3_free(errmsg);
	}

	if (rc != SQLITE_OK)
		throw std::runtime_error("Error starting transaction: " + err);
	m_transaction_active = true;
}

transaction::~transaction()
{
	try
	{
		if (m_transaction_active)
			rollback();
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error in destructor of transaction: " << ex.what() << '\n';
	}
}

void transaction::commit()
{
	if (m_transaction_active)
	{
		char *errmsg = nullptr;
		std::string err;
		int rc = sqlite3_exec(m_conn.m_impl->m_sqlite_db, "COMMIT TRANSACTION;", nullptr, nullptr, &errmsg);
		if (errmsg)
		{
			err = errmsg;
			sqlite3_free(errmsg);
		}

		if (rc != SQLITE_OK)
			throw std::runtime_error("Error committing transaction: " + err);

		m_transaction_active = false;
	}
}

void transaction::rollback()
{
	if (m_transaction_active)
	{
		char *errmsg = nullptr;
		std::string err;
		int rc = sqlite3_exec(m_conn.m_impl->m_sqlite_db, "ROLLBACK TRANSACTION;", nullptr, nullptr, &errmsg);
		if (errmsg)
		{
			err = errmsg;
			sqlite3_free(errmsg);
		}

		if (rc != SQLITE_OK)
			throw std::runtime_error("Error rolling back transaction: " + err);

		m_transaction_active = false;
	}
}

result transaction::exec(std::string query)
{
	return m_conn.exec(query);
}

result transaction::exec(std::string query, std::string &tail)
{
	return m_conn.exec(query, tail);
}

} // namespace cif::cql