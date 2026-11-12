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

#include "test-main.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cif++/cif++.hpp>
#include <cif++/cql.hpp>
#include <cif++/utilities.hpp>

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	return cif::file(text, length);
}

// --------------------------------------------------------------------

const char *kAuthors[] = {
	"Kleywegt, G.J.",
	"Bergfors, T.",
	"Senn, H.",
	"Le Motte, P.",
	"Gsell, B.",
	"Shudo, K.",
	"Jones, T.A.",

	"Banaszak, L.",
	"Winter, N.",
	"Xu, Z.",
	"Bernlohr, D.A.",
	"Cowan, S.W.",
	"Jones, T.A.",
	"Bergfors, T.",
	"Kleywegt, G.J.",
	"Jones, T.A.",
	"Cowan, S.W.",
	"Newcomer, M.E.",
	"Jones, T.A.",
	"Jones, T.A.",
	"Bergfors, T.",
	"Sedzik, J.",
	"Unge, T."
};

// Test simple SELECT
TEST_CASE("cql-1")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT name, ordinal FROM citation_author WHERE citation_id = 'primary';");
	CHECK(r.size() == 7);

	for (size_t ix = 0; auto row : r)
	{
		REQUIRE(ix < (sizeof(kAuthors) / sizeof(char *)));

		CHECK(row[0].get<std::string>() == kAuthors[ix]);
		CHECK(row[1].get<size_t>() == ix + 1);

		CHECK(row["name"].get<std::string>() == kAuthors[ix]);
		CHECK(row["ordinal"].get<size_t>() == ix + 1);

		++ix;
	}

	r = tx.exec("SELECT ordinal, name FROM citation_author WHERE citation_id = 'primary';");
	CHECK(r.size() == 7);

	for (size_t ix = 0; auto row : r)
	{
		REQUIRE(ix < (sizeof(kAuthors) / sizeof(char *)));

		CHECK(row[1].get<std::string>() == kAuthors[ix]);
		CHECK(row[0].get<size_t>() == ix + 1);

		CHECK(row["name"].get<std::string>() == kAuthors[ix]);
		CHECK(row["ordinal"].get<size_t>() == ix + 1);

		++ix;
	}

	r = tx.exec("SELECT * FROM citation_author WHERE citation_id = 'primary';");
	CHECK(r.size() == 7);

	for (int ix = 0; auto row : r)
	{
		REQUIRE(static_cast<size_t>(ix) < (sizeof(kAuthors) / sizeof(char *)));

		for (auto fld : row)
		{
			switch (fld.num())
			{
				case 0:
					CHECK(fld.name() == "citation_id");
					CHECK(fld.get<std::string>() == "primary");
					break;
				case 1:
					CHECK(fld.name() == "name");
					CHECK(fld.get<std::string>() == kAuthors[ix]);
					break;
				case 2:
					CHECK(fld.name() == "ordinal");
					CHECK(fld.get<int>() == ix + 1);
					break;
				default:
					CHECK(fld.name() == "identifier_ORCID");
					CHECK(fld.is_null());
					break;
			}
		}

		CHECK(row["name"].get<std::string>() == kAuthors[ix]);
		CHECK(row["ordinal"].get<int>() == ix + 1);
		CHECK(row["citation_id"].get<std::string>() == "primary");

		++ix;
	}
}

// Test SELECT AS
TEST_CASE("cql-2")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT name AS v1, ordinal AS v2 FROM citation_author WHERE citation_id = 'primary';");
	CHECK(r.size() == 7);

	for (size_t ix = 0; auto row : r)
	{
		REQUIRE(ix < (sizeof(kAuthors) / sizeof(char *)));

		CHECK(row[0].get<std::string>() == kAuthors[ix]);
		CHECK(row[1].get<size_t>() == ix + 1);

		CHECK(row["v1"].get<std::string>() == kAuthors[ix]);
		CHECK(row["v2"].get<size_t>() == ix + 1);

		++ix;
	}
}

TEST_CASE("cql-3")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT name FROM citation_author WHERE ordinal = 10").one_field();
	CHECK(r.get<std::string>() == kAuthors[9]);
}

TEST_CASE("cql-4")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT name FROM citation_author WHERE ordinal BETWEEN 10 AND 15");
	REQUIRE(r.size() == 6);
}

TEST_CASE("cql-5")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT (SELECT year FROM citation WHERE id = citation_id) AS jaar FROM citation_author WHERE ordinal IS 23").one_field();
	CHECK(r.name() == "jaar");
	CHECK(r.get<int>() == 1988);
}

TEST_CASE("cql-6")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT COUNT(*) FROM citation WHERE page_last IS NULL").one_field();
	CHECK(r.get<int>() == 4);

	r = tx.exec("SELECT COUNT(*) FROM citation WHERE page_last IS NOT NULL").one_field();
	CHECK(r.get<int>() == 1);
}

TEST_CASE("cql-stream-1")
{
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	auto &db = f.front();
	db.load_dictionary("mmcif_pdbx.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	for (size_t ix = 0;
		const auto &[name, ordinal] : tx.stream<std::string, size_t>(
			"SELECT name, ordinal FROM citation_author WHERE citation_id = 'primary';"))
	{
		REQUIRE(ix < (sizeof(kAuthors) / sizeof(char *)));

		CHECK(name == kAuthors[ix]);
		CHECK(ordinal == (ix + 1));

		++ix;
	}
}

// --------------------------------------------------------------------

TEST_CASE("cql-insert-1")
{
	auto f1 = R"(
data_T1
loop_
_table1.id
_table1.name
1 aap
2 noot)"_cf;
	auto f0 = f1;

	auto &db = f1.front();

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto count = tx.exec("SELECT COUNT(*) FROM table1;").one_field().get<int>();
	CHECK(count == 2);

	auto r = tx.exec("INSERT INTO table1 (id, name) VALUES (3, 'mies')");

	count = tx.exec("SELECT COUNT(*) FROM table1").one_field().get<int>();
	CHECK(count == 3);

	(void)tx.exec("DELETE FROM table1 WHERE CAST(id AS INTEGER) = 1;");

	count = tx.exec("SELECT COUNT(*) FROM table1;").one_field().get<int>();
	CHECK(count == 2);

	(void)tx.exec("UPDATE table1 SET name = 'amandel' WHERE CAST(id AS INTEGER) = 2");

	auto f2 = R"(
data_T1
loop_
_table1.id
_table1.name
2 amandel
3 mies)"_cf;

	CHECK(f1 == f2);

	tx.rollback();

	CHECK(f1 == f0);
}

// --------------------------------------------------------------------

TEST_CASE("cql-rename")
{
	auto f1 = R"(
data_T1
loop_
_table1.id
_table1.name
1 aap
2 noot)"_cf;

	auto &db = f1.front();

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	(void)tx.exec("ALTER TABLE table1 RENAME TO 'table2'");

	auto f2 = R"(
data_T1
loop_
_table2.id
_table2.name
1 aap
2 noot)"_cf;

	CHECK(f1 == f2);
}

// --------------------------------------------------------------------

TEST_CASE("cql-foreign-keys-1")
{
	std::string_view dict = R"(
data_test_dict.dic
    _datablock.id	test_dict.dic
    _datablock.description
;
    A test dictionary
;
    _dictionary.title           test_dict.dic
    _dictionary.datablock_id    test_dict.dic
    _dictionary.version         1.0

     loop_
    _item_type_list.code
    _item_type_list.primitive_code
    _item_type_list.construct
    _item_type_list.detail
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'
;              code item types/single words ...
;
               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'
;              text item types / multi-line text ...
;
               int       numb
               '[+-]?[0-9]+'
;              int item types are the subset of numbers that are the negative
               or positive integers.
;

save_cat_1
    _category.description     'A simple test category'
    _category.id              cat_1
    _category.mandatory_code  no
    _category_key.name        '_cat_1.id'

    save_

save__cat_1.id
    _item.name                '_cat_1.id'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           code
    save_

save__cat_1.name
    _item.name                '_cat_1.name'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           text
    save_

save_cat_2
    _category.description     'A second simple test category'
    _category.id              cat_2
    _category.mandatory_code  no
    _category_key.name        '_cat_2.id'
    save_

save__cat_2.id
    _item.name                '_cat_2.id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           code
    save_

save__cat_2.desc
    _item.name                '_cat_2.desc'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           text
    save_
    )";

	cif::ispanstream is_dict(dict);
	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	std::string_view data = R"(
data_test
loop_
_cat_1.id
_cat_1.name
1 Aap
2 Noot
3 Mies

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.desc
1 1 'Een dier'
2 1 'Een andere aap'
3 2 'walnoot bijvoorbeeld'
    )";

	cif::ispanstream is_data(data);
	f.load(is_data);
	f.front().set_validator(&validator);

	auto &db = f.front();

	SECTION("stream")
	{
		cif::cql::connection connection(db);
		cif::cql::transaction tx(connection);

		for (const auto &desc : tx.stream<std::string>(R"(SELECT b.desc FROM cat_1 a, cat_2 b WHERE a.id = b.parent_id AND a.name = 'Noot')"))
		{
			CHECK(desc == "walnoot bijvoorbeeld");
		}
	}

	// Check cascading delete
	SECTION("delete")
	{
		cif::cql::connection connection(db);
		cif::cql::transaction tx(connection);

		tx.exec("DELETE FROM cat_1 WHERE id = 1");
		CHECK(db["cat_1"].size() == 2);
		CHECK(db["cat_2"].size() == 1);

		tx.rollback();
		CHECK(db["cat_1"].size() == 3);
		CHECK(db["cat_2"].size() == 3);
	}

	// Check cascading update
	SECTION("update")
	{
		cif::cql::connection connection(db);
		cif::cql::transaction tx(connection);

		tx.exec("UPDATE cat_1 SET id = '4' WHERE id = '1'");
		CHECK(db["cat_1"].size() == 3);
		CHECK(db["cat_2"].size() == 3);
		CHECK(db["cat_1"].count(cif::key("id") == 4) == 1);
		CHECK(db["cat_2"].count(cif::key("parent_id") == 4) == 2);

		std::cout << db;

		tx.rollback();
		CHECK(db["cat_1"].size() == 3);
		CHECK(db["cat_2"].size() == 3);
		CHECK_FALSE(db["cat_1"].contains(cif::key("id") == 4));
		CHECK_FALSE(db["cat_2"].contains(cif::key("parent_id") == 4));

		std::cout << db;
	}
}

// --------------------------------------------------------------------

TEST_CASE("drop-table")
{
	auto f1 = R"(
data_T1
loop_
_table1.id
_table1.name
1 aap
2 noot)"_cf;

	auto &db = f1.front();

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	SECTION("commit")
	{
		(void)tx.exec("DROP TABLE table1;");
		tx.commit();

		CHECK(db.empty());
	}

	// Ah, too bad: this doesn't work
	// SECTION("rollback")
	// {
	// 	(void)tx.exec("DROP TABLE table1;");
	// 	tx.rollback();

	// 	CHECK(not db.empty());
	// 	CHECK(db["table1"].size() == 2);
	// }
}

// Items without an _item_type.code have a null type validator; creating a
// virtual table for such a category used to dereference the null m_type.
TEST_CASE("cql-7")
{
	const char *dictText = R"(
data_test.dic

_datablock.id            test.dic

_dictionary.title            test.dic
_dictionary.datablock_id     test.dic
_dictionary.version          1.0

loop_
_item_type_list.code
_item_type_list.primitive_code
_item_type_list.construct
int   numb  '[-+]?[0-9]+'
float numb  '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
code  char  '.*'
text  char  '.*'

save_foo
   _category.id              foo
   _category.mandatory_code  no
   _category_key.name        '_foo.id'
   save_

save__foo.id
	_item.name                  '_foo.id'
	_item.category_id           foo
	_item.mandatory_code        no
	_item_type.code             int
	save_

save__foo.notes
	_item.name                  '_foo.notes'
	_item.category_id           foo
	_item.mandatory_code        no
	save_
)";

	std::istringstream dictIs(dictText);
	cif::validator v(dictIs);
	cif::validator_factory::instance().add(std::move(v));

	auto f = R"(data_TEST
#
loop_
_foo.id
_foo.notes
1 "first row"
2 "second row"
)"_cf;

	auto &db = f.front();
	db.load_dictionary("test.dic");

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT id, notes FROM foo;");
	REQUIRE(r.size() == 2);

	auto ri = r.begin();
	CHECK((*ri)[0].get<int>() == 1);
	CHECK((*ri)[1].get<std::string>() == "first row");
	++ri;
	CHECK((*ri)[0].get<int>() == 2);
}

// --------------------------------------------------------------------
// Regression tests for the CQL xBestIndex planning code.
//
// These cover constraints that used to be mishandled:
//  - rowid constraints have iColumn == -1 and used to index m_items[-1]
//  - unsupported constraints (LIKE) used to leave previously set omit flags
//    in place without producing an idxStr, so xFilter never applied them
//  - IS NULL / IS NOT NULL were added to the omit list but not to idxStr
// --------------------------------------------------------------------

TEST_CASE("cql-best-index-rowid")
{
	auto f = R"(
data_TEST
loop_
_t.id
_t.name
1 aap
2 noot
3 mies
)"_cf;

	cif::cql::connection connection(f.front());
	cif::cql::transaction tx(connection);

	auto rowid = tx.exec("SELECT rowid FROM t WHERE id = '2'").one_field().get<int64_t>();

	auto r = tx.exec("SELECT id FROM t WHERE rowid = " + std::to_string(rowid));
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "2");

	r = tx.exec("SELECT id FROM t WHERE rowid >= 1 AND rowid <= 2");
	CHECK(r.size() == 2);

	// An unknown rowid must yield an empty result rather than an error or a crash.
	r = tx.exec("SELECT id FROM t WHERE rowid = " + std::to_string(rowid + 100));
	CHECK(r.empty());
}

TEST_CASE("cql-best-index-mixed-constraints")
{
	auto f = R"(
data_TEST
loop_
_t.id
_t.name
1 aap
2 appel
3 noot
)"_cf;

	cif::cql::connection connection(f.front());
	cif::cql::transaction tx(connection);

	// id is handled by xBestIndex, LIKE is not. If the id constraint were
	// marked omitted but not added to idxStr, row 2 (appel) would leak through.
	auto r = tx.exec("SELECT id FROM t WHERE id = '1' AND name LIKE 'a%'");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "1");
}

TEST_CASE("cql-best-index-is-null")
{
	auto f = R"(
data_TEST
loop_
_t.id
_t.name
1 aap
2 .
3 noot
4 appel
)"_cf;

	cif::cql::connection connection(f.front());
	cif::cql::transaction tx(connection);

	auto r = tx.exec("SELECT id FROM t WHERE id >= '1' AND name IS NULL");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "2");

	r = tx.exec("SELECT id FROM t WHERE id >= '1' AND name IS NOT NULL");
	CHECK(r.size() == 3);

	// The NULL test must appear in the plan. Otherwise it is omitted (SQLite
	// may skip checking it) while xFilter never applies it.
	auto plan = tx.exec("EXPLAIN QUERY PLAN SELECT id FROM t WHERE id >= '1' AND name IS NULL");
	std::string detail;
	for (auto row : plan)
		detail += row[3].get<std::string>();
	CHECK(detail.find("IS NULL") != std::string::npos);
}

TEST_CASE("cql-best-index-large-integers")
{
	cif::datablock db("TEST");
	auto &t = db["t"];

	t.emplace(cif::row_initializer{ { "id", int64_t{ 9007199254740992LL } }, { "name", "small" } });
	t.emplace(cif::row_initializer{ { "id", int64_t{ 9007199254740993LL } }, { "name", "big" } });

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	// These values cannot both be represented exactly as a double.
	auto r = tx.exec("SELECT name FROM t WHERE id = 9007199254740993");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "big");

	r = tx.exec("SELECT name FROM t WHERE id = 9007199254740992");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "small");

	r = tx.exec("SELECT name FROM t WHERE id > 9007199254740992");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "big");
}

TEST_CASE("cql-rowid-not-reused")
{
	auto f = R"(
data_TEST
loop_
_t.id
_t.name
1 aap
2 noot
)"_cf;

	cif::cql::connection connection(f.front());
	cif::cql::transaction tx(connection);

	auto old_rowid = tx.exec("SELECT rowid FROM t WHERE id = '1'").one_field().get<int64_t>();

	(void)tx.exec("DELETE FROM t WHERE id = '1'");
	(void)tx.exec("INSERT INTO t (id, name) VALUES ('3', 'mies')");

	auto r = tx.exec("SELECT rowid FROM t WHERE id = '3'");
	REQUIRE(r.size() == 1);

	// Row ids must not be recycled: reusing the old id would make a stale
	// rowid silently refer to a different row.
	CHECK(r.front()[0].get<int64_t>() != old_rowid);
	CHECK(tx.exec("SELECT id FROM t WHERE rowid = " + std::to_string(old_rowid)).empty());
}

TEST_CASE("cql-best-index-non-finite")
{
	cif::datablock db("TEST");
	auto &t = db["t"];

	t.emplace(cif::row_initializer{ { "id", 1 }, { "x", 1.0 }, { "name", "a" } });
	t.emplace(cif::row_initializer{ { "id", 2 }, { "x", 2.0 }, { "name", "b" } });

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	// 1e999 is +Inf. Non-finite constraint values must not be serialized into
	// the plan; they have to be left for SQLite to evaluate.
	auto r = tx.exec("SELECT id FROM t WHERE x < 1e999");
	CHECK(r.size() == 2);

	r = tx.exec("SELECT id FROM t WHERE x = 1e999");
	CHECK(r.empty());

	r = tx.exec("SELECT id FROM t WHERE x > -1e999");
	CHECK(r.size() == 2);
}

TEST_CASE("cql-best-index-float-precision")
{
	cif::datablock db("TEST");
	auto &t = db["t"];

	t.emplace(cif::row_initializer{ { "id", 1 }, { "x", 123456.789 }, { "name", "a" } });
	t.emplace(cif::row_initializer{ { "id", 2 }, { "x", 1.23456789 }, { "name", "b" } });

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	// The constraint value must survive the round trip through idxStr without
	// losing precision (the stream's default precision is only 6 digits).
	auto r = tx.exec("SELECT name FROM t WHERE x = 123456.789");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "a");

	r = tx.exec("SELECT name FROM t WHERE x = 1.23456789");
	REQUIRE(r.size() == 1);
	CHECK(r.front()[0].get<std::string>() == "b");
}
TEST_CASE("cql-foreign-keys-composite")
{
	std::string_view dict = R"(
data_test_dict.dic
    _datablock.id	test_dict.dic
    _datablock.description
;
    A test dictionary
;
    _dictionary.title           test_dict.dic
    _dictionary.datablock_id    test_dict.dic
    _dictionary.version         1.0

     loop_
    _item_type_list.code
    _item_type_list.primitive_code
    _item_type_list.construct
    _item_type_list.detail
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'
;              code item types/single words ...
;
               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'
;              text item types / multi-line text ...
;
               int       numb
               '[+-]?[0-9]+'
;              int item types are the subset of numbers that are the negative
               or positive integers.
;

save_cat_1
    _category.description     'A simple test category'
    _category.id              cat_1
    _category.mandatory_code  no
    loop_
    _category_key.name
    '_cat_1.id'
    '_cat_1.id2'
    save_

save__cat_1.id
    _item.name                '_cat_1.id'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           code
    save_

save__cat_1.id2
    _item.name                '_cat_1.id2'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_linked.child_name   '_cat_2.parent_id2'
    _item_linked.parent_name  '_cat_1.id2'
    _item_type.code           code
    save_

save__cat_1.name
    _item.name                '_cat_1.name'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           text
    save_

save_cat_2
    _category.description     'A second simple test category'
    _category.id              cat_2
    _category.mandatory_code  no
    _category_key.name        '_cat_2.id'
    save_

save__cat_2.id
    _item.name                '_cat_2.id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           code
    save_

save__cat_2.parent_id2
    _item.name                '_cat_2.parent_id2'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           code
    save_

save__cat_2.desc
    _item.name                '_cat_2.desc'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_aliases.dictionary  cif_core.dic
    _item_aliases.version     2.0.1
    _item_type.code           text
    save_
    )";

	cif::ispanstream is_dict(dict);
	cif::validator validator(is_dict);

	cif::file f;

	std::string_view data = R"(
data_test
loop_
_cat_1.id
_cat_1.id2
_cat_1.name
1 10 Aap
2 20 Noot

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.parent_id2
_cat_2.desc
1 1 10 'Een dier'
2 1 10 'Een andere aap'
3 2 20 'walnoot bijvoorbeeld'
    )";

	cif::ispanstream is_data(data);
	f.load(is_data);
	f.front().set_validator(&validator);

	auto &db = f.front();

	cif::cql::connection connection(db);
	cif::cql::transaction tx(connection);

	// Updating a parent row must cascade both key items to the child. The
	// generated SQL used to contain a second SET keyword, which made the
	// statement invalid.
	tx.exec("UPDATE cat_1 SET id = '3' WHERE id = '1'");

	CHECK(db["cat_1"].count(cif::key("id") == 3) == 1);
	CHECK(db["cat_2"].count(cif::key("parent_id") == 3 and cif::key("parent_id2") == 10) == 2);
	CHECK(db["cat_2"].count(cif::key("parent_id") == 1) == 0);

	// The child rows are updated through the generated cascade statement, so
	// SQLite must also restore them when the transaction rolls back.
	tx.rollback();

	CHECK(db["cat_1"].count(cif::key("id") == 1) == 1);
	CHECK(db["cat_2"].count(cif::key("parent_id") == 1 and cif::key("parent_id2") == 10) == 2);
}
