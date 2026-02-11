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

#include "cif++/category.hpp"
#include "test-main.hpp"

// #include <cif++.hpp>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cif++/category.hpp>
#include <cif++/datablock.hpp>
#include <cif++/file.hpp>
#include <cif++/item.hpp>
#include <cif++/parser.hpp>
#include <cif++/row.hpp>
#include <cif++/validate.hpp>
#include <exception>
#include <stdexcept>

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

TEST_CASE("text_1")
{
	CHECK(cif::iequals("TEST", "test"));
	CHECK(cif::iequals(std::string_view{ "TEST" }, std::string_view{ "test" }));

	CHECK(cif::icompare("TEST", "test") == 0);
	CHECK(cif::icompare(std::string_view{ "TEST" }, std::string_view{ "test" }) == 0);

	CHECK(cif::icompare("TEST1", "test") > 0);
	CHECK(cif::icompare(std::string_view{ "TEST1" }, std::string_view{ "test" }) > 0);

	CHECK(cif::icompare("aap", "noot") < 0);
	CHECK(cif::icompare(std::string_view{ "aap" }, std::string_view{ "noot" }) < 0);
}

// // --------------------------------------------------------------------

// TEST_CASE("from_chars_1")
// {
// auto f = R"(data_TEST
// #
// loop_
// _test.v
// 616.487
// 616.487000
//     )"_cf;

// 	auto &db = f.front();
// 	auto &c = db.front();

// 	auto r1 = c.front();
// 	CHECK(r1.get<double>("v") == 616.487);
// 	CHECK(r1["v"].compare(616.487) == 0);

// 	auto r2 = c.back();
// 	CHECK(r2.get<double>("v") == 616.487);
// 	CHECK(r2["v"].compare(616.487) == 0);
// }

// --------------------------------------------------------------------

TEST_CASE("id_1")
{
	CHECK(cif::cif_id_for_number(0) == "A");
	CHECK(cif::cif_id_for_number(25) == "Z");
	CHECK(cif::cif_id_for_number(26) == "AA");
	CHECK(cif::cif_id_for_number(26 + 1) == "AB");

	CHECK(cif::cif_id_for_number(26 + 26 * 26 - 1) == "ZZ");
	CHECK(cif::cif_id_for_number(26 + 26 * 26) == "AAA");
	CHECK(cif::cif_id_for_number(26 + 26 * 26 + 1) == "AAB");

	std::set<std::string> testset;

	for (int i = 0; i < 100000; ++i)
	{
		std::string id = cif::cif_id_for_number(i);
		CHECK(testset.count(id) == 0);
		testset.insert(id);
	}
	CHECK(testset.size() == 100000);
}

// --------------------------------------------------------------------

TEST_CASE("cc_1")
{
	std::tuple<std::string_view, float, char> tests[] = {
		{ "1.0", 1.0f, 0 },
		{ "1.0e10", 1.0e10f, 0 },
		{ "-1.1e10", -1.1e10f, 0 },
		{ "-.2e11", -.2e11f, 0 },
		{ "1.3e-10", 1.3e-10f, 0 },

		{ "1.0 ", 1.0f, ' ' },
		{ "1.0e10 ", 1.0e10f, ' ' },
		{ "-1.1e10 ", -1.1e10f, ' ' },
		{ "-.2e11 ", -.2e11f, ' ' },
		{ "1.3e-10 ", 1.3e-10f, ' ' },

		{ "3.0", 3.0f, 0 },
		{ "3.0 ", 3.0f, ' ' },

		{ "3.000000", 3.0f, 0 },
		{ "3.000000 ", 3.0f, ' ' },
	};

	for (const auto &[txt, val, ch] : tests)
	{
		float tv;
		const auto &[ptr, ec] = cif::from_chars(txt.data(), txt.data() + txt.length(), tv);

		CHECK_FALSE((bool)ec);
		CHECK(tv == val);
		if (ch != 0)
			CHECK(*ptr == ch);
	}
}

TEST_CASE("cc_2")
{
	std::tuple<float, int, std::string_view> tests[] = {
		{ 1.1f, 1, "1.1" }
	};

	for (const auto &[val, prec, test] : tests)
	{
		char buffer[64] = {};
		const auto &[ptr, ec] = std::to_chars(buffer, buffer + sizeof(buffer), val, std::chars_format::fixed, prec);

		CHECK(ec == std::errc{});

		// This line generates a linker error on Windows:
		// CHECK(buffer == test);
		CHECK(std::string{ buffer } == std::string{ test });
	}
}

TEST_CASE("item_0")
{
	cif::item i1("v1", "tekst");
	cif::item i2("v2", 2);
	cif::item i3("v3", 3.0);
	cif::item i4("v4", true);
	cif::item i5("v5", nullptr);

	static_assert(cif::IntegralType<int>);
	static_assert(cif::FloatType<double>);

	CHECK(i1.value().get<std::string>() == "tekst");
	CHECK(i2.value().get<int>() == 2);
	CHECK(i3.value().get<double>() == 3.0);
	CHECK(i4.value().get<bool>() == true);
	CHECK(i5.value().is_null());
	CHECK(i5.value().is_missing());
	CHECK(i5.value().empty());

	i2.value() = false;
	CHECK(i2.value().type() == cif::item_value_type::INT);
	CHECK(i2.value().get<bool>() == false);

	cif::item i6 = std::move(i1);
	CHECK(i6.value().get<std::string>() == "tekst");
	CHECK(i1.value().is_null());
}

// TEST_CASE("row_1")
// {
// 	cif::row r;

// 	r.append(0, 1);
// 	r.append(1, "twee");

// 	cif::category cat("cat");

// 	cat.add_item("een");
// 	cat.add_item("twee");

// 	cif::row_handle rh(cat, r);

// 	const auto &[a, b] = rh.get<int, std::string>("een", "twee");
// 	CHECK(a == 1);
// 	CHECK(b == "twee");

// 	rh.assign("twee", 3.0, false);
// 	CHECK(rh[1].type() == cif::item_value_type::FLOAT);
// 	CHECK(rh[1].get<double>() == 3.0);
// }

TEST_CASE("cc_3")
{
	cif::category c("foo");
	c.emplace({
		{ "f-1", 1 },
		{ "f-2", "-1" },
		{ "f-3", "+1" },
		{ "f-4", " 1" },
		{ "f-5", " +1" },
		{ "f-6", "1 " },
	});

	auto row = c.front();
	CHECK(row["f-1"].get<int>() == 1);
	CHECK(row["f-2"].get<int>() == -1);

	CHECK_FALSE(row["f-3"].is_number());
	CHECK_FALSE(row["f-4"].is_number());
	CHECK_FALSE(row["f-5"].is_number());
	CHECK_FALSE(row["f-6"].is_number());
}

TEST_CASE("item_1")
{
	using namespace cif;

	item i1("1", "1");
	item i2("2", 2.0f);
	item i3("3", '3');

	item ci1(i1);
	item ci2(i2);
	item ci3(i3);

	CHECK(i1.value() == ci1.value());
	CHECK(i2.value() == ci2.value());
	CHECK(i3.value() == ci3.value());

	item mi1(std::move(ci1));
	item mi2(std::move(ci2));
	item mi3(std::move(ci3));

	CHECK(i1.value() == mi1.value());
	CHECK(i2.value() == mi2.value());
	CHECK(i3.value() == mi3.value());

	CHECK(ci1.empty());
	CHECK(ci2.empty());
	CHECK(ci3.empty());
}

TEST_CASE("item_2")
{
	using namespace cif;

	cif::item i0("test1");
	CHECK(i0.value().empty());
	CHECK(i0.value().type() == cif::item_value_type::MISSING);

	cif::item i1("test1", std::optional<float>());
	CHECK(i0.value().empty());
	// CHECK(i1.value() == "?");

	cif::item i2("test1", std::make_optional<float>(1.f));
	// CHECK(i2.value() == "1");

	// cif::item i3("test1", std::optional<float>(), 2);
	// CHECK(i3.value() == "?");

	// cif::item i4("test1", std::make_optional<float>(1.f), 2);
	// CHECK(i4.value() == "1.00");
}

// --------------------------------------------------------------------

TEST_CASE("r_1")
{
	cif::category c("foo");
	c.emplace({
		{ "f-1", 1 },
		{ "f-2", "two" },
		{ "f-3", 3.0f /* , 3 */ },
	});

	auto row = c.front();
	// CHECK(row["f-1"].compare(1) == 0);
	// CHECK(row["f-2"].compare("two") == 0);
	// CHECK(row["f-3"].compare(3.0f) == 0); // This fails when running in valgrind... sigh

	const auto &[f1, f2, f3] = row.get<int, std::string, float>("f-1", "f-2", "f-3");

	CHECK(f1 == 1);
	CHECK(f2 == "two");
	CHECK(f3 == 3.0f); // This fails when running in valgrind... sigh

	CHECK(row.get<int>("f-1") == 1);
	CHECK(row.get<std::string>("f-2") == "two");
	CHECK(row.get<float>("f-3") == 3.0f);

	int f_1;
	std::string f_2;
	float f_3;

	cif::tie(f_1, f_2, f_3) = row.get("f-1", "f-2", "f-3");

	CHECK(f_1 == 1);
	CHECK(f_2 == "two");
	CHECK(f_3 == 3.0f); // This fails when running in valgrind... sigh
}

TEST_CASE("r_2")
{
	cif::category c("foo");

	for (std::size_t i = 1; i < 256; ++i)
	{
		c.emplace({ { "id", i },
			{ "txt", std::string(i, 'x') } });
	}
}

TEST_CASE("c_1")
{
	cif::category c("foo");

	c.emplace({ { "id", 1 }, { "s", "aap" } });
	c.emplace({ { "id", 2 }, { "s", "noot" } });
	c.emplace({ { "id", 3 }, { "s", "mies" } });

	int n = 1;

	const char *ts[] = { "aap", "noot", "mies" };

	for (auto r : c)
	{
		CHECK(r["id"].get<int>() == n);
		CHECK(r["s"].compare(ts[n - 1]) == 0);
		++n;
	}

	n = 1;

	for (auto r : c)
	{
		int i;
		std::string s;

		cif::tie(i, s) = r.get("id", "s");

		CHECK(i == n);
		CHECK(s.compare(ts[n - 1]) == 0);
		++n;
	}

	n = 1;

	for (const auto &[i, s] : c.rows<int, std::string>("id", "s"))
	{
		CHECK(i == n);
		CHECK(s.compare(ts[n - 1]) == 0);
		++n;
	}
}

TEST_CASE("c_2")
{
	std::tuple<int, const char *> D[] = {
		{ 1, "aap" },
		{ 2, "noot" },
		{ 3, "mies" }
	};

	cif::category c("foo");

	for (const auto &[id, s] : D)
		c.emplace({ { "id", id }, { "s", s } });

	CHECK(not c.empty());
	CHECK(c.size() == 3);

	cif::category c2(c);

	CHECK(not c2.empty());
	CHECK(c2.size() == 3);

	cif::category c3(std::move(c));

	CHECK(not c3.empty());
	CHECK(c3.size() == 3);

	CHECK(c.empty());
	CHECK(c.size() == 0);

	c = c3;

	CHECK(not c.empty());
	CHECK(c.size() == 3);

	c = std::move(c2);

	CHECK(not c.empty());
	CHECK(c.size() == 3);
}

TEST_CASE("c_3")
{
	std::tuple<int, const char *> D[] = {
		{ 1, "aap" },
		{ 2, "noot" },
		{ 3, "mies" }
	};

	cif::category c("foo");

	for (const auto &[id, s] : D)
		c.emplace({ { "id", id }, { "s", s } });

	cif::category c2("bar");

	for (auto r : c)
		c2.emplace(r);

	// CHECK(c == c2);
}

TEST_CASE("ci_1")
{
	cif::category c("foo");

	c.emplace({ { "id", 1 }, { "s", "aap" } });
	c.emplace({ { "id", 2 }, { "s", "noot" } });
	c.emplace({ { "id", 3 }, { "s", "mies" } });

	cif::category::iterator i1 = c.begin();
	cif::category::const_iterator i2 = c.cbegin();
	cif::category::const_iterator i3 = c.begin();

	cif::category::const_iterator i4 = i2;
	cif::category::const_iterator i5 = i1;

	CHECK(i1 == i2);
	CHECK(i1 == i3);
	CHECK(i1 == i4);
	CHECK(i1 == i5);
}

// TEST_CASE("os_1")
// {
// 	using namespace cif::literals;
// 	using namespace std::literals;

// 	std::tuple<int, const char *> D[] = {
// 		{ 1, "aap" },
// 		{ 2, "noot" },
// 		{ 3, "mies" }
// 	};

// 	cif::category c("foo");

// 	for (const auto &[id, s] : D)
// 		c.emplace({ { "id", id }, { "s", s } });

// 	for (auto rh : c)
// 	{
// 		rh["o"].os(1, ',', 2, ": ", rh.get<std::string>("s"));
// 	}

// 	for (const auto &[id, s] : D)
// 	{
// 		auto rh = c.find1("id"_key == id);

// 		CHECK(rh.get<int>("id") == id);
// 		CHECK(rh.get<std::string>("s") == s);
// 		CHECK(rh.get<std::string>("o") == "1,2: "s + s);
// 	}
// }

// --------------------------------------------------------------------

TEST_CASE("get_1")
{
	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
4 ?
5 .
    )"_cf;

	REQUIRE(f.size() == 1);
	REQUIRE(f.contains("TEST"));
	REQUIRE(f.front().size() == 1);
	REQUIRE(f.front().contains("test"));
	REQUIRE(f.front()["test"].size() == 5);
	REQUIRE(f.front()["test"].front().size() == 2);

	for (auto r : f.front()["test"])
	{
		int id;
		std::optional<std::string> name;

		cif::tie(id, name) = r.get("id", "name");

		switch (id)
		{
			case 1: REQUIRE(name.has_value()); CHECK(*name == "aap"); break;
			case 2: REQUIRE(name.has_value()); CHECK(*name == "noot"); break;
			case 3: REQUIRE(name.has_value()); CHECK(*name == "mies"); break;
			case 4:
			case 5: CHECK_FALSE(name.has_value()); break;
			default: CHECK(false);
		}
	}
}

// --------------------------------------------------------------------

TEST_CASE("f_1")
{
	// using namespace mmcif;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
    )"_cf;

	CHECK(f.empty() == false);
	CHECK(f.size() == 1);

	auto &db = f.front();

	CHECK(db.name() == "TEST");

	auto &test = db["test"];
	CHECK(test.size() == 3);

	const char *ts[] = { "aap", "noot", "mies" };

	int n = 1;
	for (const auto &[i, s] : test.rows<int, std::string>("id", "name"))
	{
		CHECK(i == n);
		CHECK(s.compare(ts[n - 1]) == 0);
		++n;
	}

	auto n2 = test.erase(cif::key("id") == 1, [](cif::row_handle r)
		{
        CHECK(r["id"].get<int>() == 1);
        CHECK(r["name"].get<std::string>() == "aap"); });

	CHECK(n2 == 1);

	// for (auto r: test)
	// 	test.erase(r);

	test.clear();
	CHECK(test.empty());

	// fill again.

	test.emplace({ { "id", "1" }, { "name", "aap" } });
	test.emplace({ { "id", "2" }, { "name", "noot" } });
	test.emplace({ { "id", "3" }, { "name", "mies" } });

	n = 1;
	for (const auto &[i, s] : test.rows<int, std::string>("id", "name"))
	{
		CHECK(i == n);
		CHECK(s.compare(ts[n - 1]) == 0);
		++n;
	}
}

// --------------------------------------------------------------------

TEST_CASE("ut2")
{
	// using namespace mmcif;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
_test.value
1 aap   1.0
2 noot  1.1
3 mies  1.2
    )"_cf;

	auto &db = f.front();

	CHECK(db.name() == "TEST");

	auto &test = db["test"];
	CHECK(test.size() == 3);

	int n = 0;
	for (auto r : test.find(cif::key("name") == "aap"))
	{
		CHECK(++n == 1);
		CHECK(r["id"].get<int>() == 1);
		CHECK(r["name"].get<std::string>() == "aap");
		CHECK(r["value"].get<float>() == 1.0f);
	}

	auto t = test.find(cif::key("id") == 1);
	REQUIRE(not t.empty());
	CHECK(t.front()["name"].get<std::string>() == "aap");

	auto t2 = test.find(cif::key("value") == 1.2);
	REQUIRE(not t2.empty());
	CHECK(t2.front()["name"].get<std::string>() == "mies");
}

TEST_CASE("ut3")
{
	using namespace cif::literals;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
_test.value
1 aap   1.0
2 noot  1.1
3 mies  1.2
4 boom  .
5 roos  ?
    )"_cf;

	auto &db = f.front();

	CHECK(db.name() == "TEST");

	auto &test = db["test"];
	CHECK(test.size() == 5);

	CHECK(test.contains("value"_key == cif::null));
	CHECK(test.find("value"_key == cif::null).size() == 2);
}

// --------------------------------------------------------------------

TEST_CASE("compare-with-float")
{
	auto f = R"(data_TEST
#
loop_
_test.id
_test.value
1 1.0
2 2.0
3 3.0
4 ?
5 .
    )"_cf;

	auto &db = f.front();
	auto &cat = db["test"];

	CHECK(cat.find_first<int>(cif::key("value") == 1.0, "id") == 1);
}

// // --------------------------------------------------------------------

// TEST_CASE("sw_1")
// {
// 	using namespace cif::literals;

// 	auto f = R"(data_TEST
// #
// loop_
// _test.id
// _test.name
// _test.value
// 1 aap   1.0
// 2 noot  1.1
// 3 mies  1.2
//     )"_cf;

// 	auto &db = f.front();
// 	auto &test = db["test"];

// 	swap(test.front()["name"], test.back()["name"]);

// 	CHECK(test.find1<std::string>("id"_key == 1, "name") == "mies");
// 	CHECK(test.find1<std::string>("id"_key == 3, "name") == "aap");

// 	swap(test.front()["name"], test.back()["name"]);

// 	CHECK(test.find1<std::string>("id"_key == 1, "name") == "aap");
// 	CHECK(test.find1<std::string>("id"_key == 3, "name") == "mies");
// }

// --------------------------------------------------------------------

TEST_CASE("d0")
{
	auto v = cif::validator_factory::instance().get("mmcif_pdbx.dic");
	CHECK(v != nullptr);
}

// --------------------------------------------------------------------

TEST_CASE("d1")
{
	const char dict[] = R"(
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

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
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

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	SECTION("one")
	{
		auto &cat1 = f.front()["cat_1"];
		auto &cat2 = f.front()["cat_2"];

		CHECK(cat1.size() == 3);
		CHECK(cat2.size() == 3);

		cat1.erase(cif::key("id") == "1");

		CHECK(cat1.size() == 2);
		CHECK(cat2.size() == 1);

		// CHECK_THROWS_AS(cat2.emplace({
		//     { "id", 4 },
		//     { "parent_id", 4 },
		//     { "desc", "moet fout gaan" }
		// }), std::exception);

		CHECK_THROWS_AS(cat2.emplace({ { "id", "vijf" }, // <- invalid value
							{ "parent_id", 2 },
							{ "desc", "moet fout gaan" } }),
			std::exception);
	}

	// SECTION("two")
	// {
	// 	auto &cat1 = f.front()["cat_1"];
	// 	auto &cat2 = f.front()["cat_2"];

	// 	cat1.update_value(cif::all(), "id", [](std::string_view v) -> std::string
	// 	{
	// 		int vi;
	// 		auto [ec, ptr] = std::from_chars(v.data(), v.data() + v.length(), vi);
	// 		return std::to_string(vi + 1);
	// 	});

	// 	CHECK(cat1.find1<std::string>(cif::key("id") == 2, "name") == "Aap");
	// 	CHECK(cat1.find1<std::string>(cif::key("id") == 3, "name") == "Noot");
	// 	CHECK(cat1.find1<std::string>(cif::key("id") == 4, "name") == "Mies");

	// 	CHECK(cat2.find1<int>(cif::key("id") == 1, "parent_id") == 2);
	// 	CHECK(cat2.find1<int>(cif::key("id") == 2, "parent_id") == 2);
	// 	CHECK(cat2.find1<int>(cif::key("id") == 3, "parent_id") == 3);
	// }
}

// --------------------------------------------------------------------

TEST_CASE("d2")
{
	const char dict[] = R"(
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
               ucode     uchar
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'
;              code item types/single words, case insensitive
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
    _item_type.code           code
    save_

save__cat_1.c
    _item.name                '_cat_1.c'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           ucode
    save_
)";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.c
aap  Aap
noot Noot
mies Mies
)";

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];

	CHECK(cat1.size() == 3);

	cat1.erase(cif::key("id") == "AAP");

	CHECK(cat1.size() == 3);

	cat1.erase(cif::key("id") == "noot");

	CHECK(cat1.size() == 2);

	// should fail with duplicate key:
	CHECK_THROWS_AS(cat1.emplace({ { "id", "aap" },
						{ "c", "2e-aap" } }),
		std::exception);

	cat1.erase(cif::key("id") == "aap");

	CHECK(cat1.size() == 1);

	cat1.emplace({ { "id", "aap" },
		{ "c", "2e-aap" } });

	CHECK(cat1.size() == 2);
}

// --------------------------------------------------------------------

TEST_CASE("d3")
{
	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

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
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           code
    save_

save__cat_1.name1
    _item.name                '_cat_1.name1'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           text
    save_

save__cat_1.name2
    _item.name                '_cat_1.name2'
    _item.category_id         cat_1
    _item.mandatory_code      no
    _item_linked.child_name   '_cat_2.name2'
    _item_linked.parent_name  '_cat_1.name2'
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
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_2.name2
    _item.name                '_cat_2.name2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    _item_type.code           text
    save_

save__cat_2.desc
    _item.name                '_cat_2.desc'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           text
    save_
    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.name1
_cat_1.name2
1 Aap   aap
2 Noot  noot
3 Mies  mies

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.name2
_cat_2.desc
1 1 aap   'Een dier'
2 1 .     'Een andere aap'
3 2 noot  'walnoot bijvoorbeeld'
4 2 n2     hazelnoot
    )";

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];

	// check a rename in parent and child

	for (auto r : cat1.find(cif::key("id") == "1"))
	{
		r["id"] = "10";
		break;
	}

	CHECK(cat1.size() == 3);
	CHECK(cat2.size() == 4);

	CHECK(cat1.find(cif::key("id") == "1").size() == 0);
	CHECK(cat1.find(cif::key("id") == "10").size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == "1").size() == 0);
	CHECK(cat2.find(cif::key("parent_id") == "10").size() == 2);

	// check a rename in parent and child, this time only one child should be renamed

	for (auto r : cat1.find(cif::key("id") == "2"))
	{
		r["id"] = "20";
		break;
	}

	CHECK(cat1.size() == 3);
	CHECK(cat2.size() == 4);

	CHECK(cat1.find(cif::key("id") == "2").size() == 0);
	CHECK(cat1.find(cif::key("id") == "20").size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == "2").size() == 1);
	CHECK(cat2.find(cif::key("parent_id") == "20").size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == "2" and cif::key("name2") == "noot").size() == 0);
	CHECK(cat2.find(cif::key("parent_id") == "2" and cif::key("name2") == "n2").size() == 1);
	CHECK(cat2.find(cif::key("parent_id") == "20" and cif::key("name2") == "noot").size() == 1);
	CHECK(cat2.find(cif::key("parent_id") == "20" and cif::key("name2") == "n2").size() == 0);

	// --------------------------------------------------------------------

	cat1.erase(cif::key("id") == "10");

	CHECK(cat1.size() == 2);
	CHECK(cat2.size() == 2); // TODO: Is this really what we want?

	cat1.erase(cif::key("id") == "20");

	CHECK(cat1.size() == 1);
	CHECK(cat2.size() == 1); // TODO: Is this really what we want?
}

// --------------------------------------------------------------------

TEST_CASE("d4")
{
	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

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
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           int
    save_

save__cat_1.id2
    _item.name                '_cat_1.id2'
    _item.category_id         cat_1
    _item.mandatory_code      no
    _item_linked.child_name   '_cat_2.parent_id2'
    _item_linked.parent_name  '_cat_1.id2'
    _item_type.code           code
    save_

save__cat_1.id3
    _item.name                '_cat_1.id3'
    _item.category_id         cat_1
    _item.mandatory_code      no
    _item_linked.child_name   '_cat_2.parent_id3'
    _item_linked.parent_name  '_cat_1.id3'
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
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    save_

save__cat_2.parent_id2
    _item.name                '_cat_2.parent_id2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    save_

save__cat_2.parent_id3
    _item.name                '_cat_2.parent_id3'
    _item.category_id         cat_2
    _item.mandatory_code      no
    save_

    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.id2
_cat_1.id3
1 aap   aap
2 .     noot
3 mies  .
4 .     .

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.parent_id2
_cat_2.parent_id3
 1 1 aap   aap
 2 1 .     x
 3 1 aap   .
 4 2 noot  noot
 5 2 .     noot
 6 2 noot  .
 7 2 .     .
 8 3 mies  mies
 9 3 .     mies
10 3 mies  .
11 4 roos  roos
12 4 .     roos
13 4 roos  .
    )";

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];

	// check a rename in parent and child

	for (auto r : cat1.find(cif::key("id") == 1))
	{
		r["id"] = 10;
		break;
	}

	CHECK(cat1.size() == 4);
	CHECK(cat2.size() == 13);

	CHECK(cat1.find(cif::key("id") == 1).size() == 0);
	CHECK(cat1.find(cif::key("id") == 10).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 1).size() == 1);
	CHECK(cat2.find(cif::key("parent_id") == 10).size() == 2);

	for (auto r : cat1.find(cif::key("id") == 2))
	{
		r["id"] = 20;
		break;
	}

	CHECK(cat1.size() == 4);
	CHECK(cat2.size() == 13);

	CHECK(cat1.find(cif::key("id") == 2).size() == 0);
	CHECK(cat1.find(cif::key("id") == 20).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 2).size() == 2);
	CHECK(cat2.find(cif::key("parent_id") == 20).size() == 2);

	for (auto r : cat1.find(cif::key("id") == 3))
	{
		r["id"] = 30;
		break;
	}

	CHECK(cat1.size() == 4);
	CHECK(cat2.size() == 13);

	CHECK(cat1.find(cif::key("id") == 3).size() == 0);
	CHECK(cat1.find(cif::key("id") == 30).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 3).size() == 2);
	CHECK(cat2.find(cif::key("parent_id") == 30).size() == 1);

	for (auto r : cat1.find(cif::key("id") == 4))
	{
		r["id"] = 40;
		break;
	}

	CHECK(cat1.size() == 4);
	CHECK(cat2.size() == 13);

	CHECK(cat1.find(cif::key("id") == 4).size() == 0);
	CHECK(cat1.find(cif::key("id") == 10).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 4).size() == 3);
	CHECK(cat2.find(cif::key("parent_id") == 40).size() == 0);
}

// --------------------------------------------------------------------

TEST_CASE("d5")
{
	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

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
    _item_type.code           int
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
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    save_

save__cat_2.parent_id2
    _item.name                '_cat_2.parent_id2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    save_

save__cat_2.parent_id3
    _item.name                '_cat_2.parent_id3'
    _item.category_id         cat_2
    _item.mandatory_code      no
    save_

loop_
_pdbx_item_linked_group_list.child_category_id
_pdbx_item_linked_group_list.link_group_id
_pdbx_item_linked_group_list.child_name
_pdbx_item_linked_group_list.parent_name
_pdbx_item_linked_group_list.parent_category_id
cat_2 1 '_cat_2.parent_id'  '_cat_1.id' cat_1
cat_2 2 '_cat_2.parent_id2' '_cat_1.id' cat_1
cat_2 3 '_cat_2.parent_id3' '_cat_1.id' cat_1

loop_
_pdbx_item_linked_group.category_id
_pdbx_item_linked_group.link_group_id
_pdbx_item_linked_group.label
cat_2 1 cat_2:cat_1:1
cat_2 2 cat_2:cat_1:2
cat_2 3 cat_2:cat_1:3
    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
1
2
3

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.parent_id2
_cat_2.parent_id3
 1 1 ? ?
 2 ? 1 ?
 3 ? ? 1
 4 2 2 ?
 5 2 ? 2
 6 ? 2 2
 7 3 3 3
    )";

	// --------------------------------------------------------------------

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];

	// --------------------------------------------------------------------
	// check iterate children

	auto PR2set = cat1.find(cif::key("id") == 2);
	REQUIRE(PR2set.size() == 1);
	auto PR2 = PR2set.front();
	CHECK(PR2["id"].get<int>() == 2);

	auto CR2set = cat1.get_children(PR2, cat2);
	CHECK(CR2set.size() == 3);

	std::vector<int> CRids;
	std::ranges::transform(CR2set, std::back_inserter(CRids), [](cif::row_handle r)
		{ return r["id"].get<int>(); });
	std::ranges::sort(CRids);
	CHECK(CRids == std::vector<int>({ 4, 5, 6 }));

	// check a rename in parent and child

	for (auto r : cat1.find(cif::key("id") == 1))
	{
		r["id"] = 10;
		break;
	}

	CHECK(cat1.size() == 3);
	CHECK(cat2.size() == 7);

	CHECK(cat1.find(cif::key("id") == 1).size() == 0);
	CHECK(cat1.find(cif::key("id") == 10).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 1).size() == 0);
	CHECK(cat2.find(cif::key("parent_id2") == 1).size() == 0);
	CHECK(cat2.find(cif::key("parent_id3") == 1).size() == 0);
	CHECK(cat2.find(cif::key("parent_id") == 10).size() == 1);
	CHECK(cat2.find(cif::key("parent_id2") == 10).size() == 1);
	CHECK(cat2.find(cif::key("parent_id3") == 10).size() == 1);

	for (auto r : cat1.find(cif::key("id") == 2))
	{
		r["id"] = 20;
		break;
	}

	CHECK(cat1.size() == 3);
	CHECK(cat2.size() == 7);

	CHECK(cat1.find(cif::key("id") == 2).size() == 0);
	CHECK(cat1.find(cif::key("id") == 20).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 2).size() == 0);
	CHECK(cat2.find(cif::key("parent_id2") == 2).size() == 0);
	CHECK(cat2.find(cif::key("parent_id3") == 2).size() == 0);
	CHECK(cat2.find(cif::key("parent_id") == 20).size() == 2);
	CHECK(cat2.find(cif::key("parent_id2") == 20).size() == 2);
	CHECK(cat2.find(cif::key("parent_id3") == 20).size() == 2);

	for (auto r : cat1.find(cif::key("id") == 3))
	{
		r["id"] = 30;
		break;
	}

	CHECK(cat1.size() == 3);
	CHECK(cat2.size() == 7);

	CHECK(cat1.find(cif::key("id") == 3).size() == 0);
	CHECK(cat1.find(cif::key("id") == 30).size() == 1);

	CHECK(cat2.find(cif::key("parent_id") == 3).size() == 0);
	CHECK(cat2.find(cif::key("parent_id2") == 3).size() == 0);
	CHECK(cat2.find(cif::key("parent_id3") == 3).size() == 0);
	CHECK(cat2.find(cif::key("parent_id") == 30).size() == 1);
	CHECK(cat2.find(cif::key("parent_id2") == 30).size() == 1);
	CHECK(cat2.find(cif::key("parent_id3") == 30).size() == 1);

	// test delete

	cat1.erase(cif::key("id") == 10);
	CHECK(cat1.size() == 2);
	CHECK(cat2.size() == 4);

	cat1.erase(cif::key("id") == 20);
	CHECK(cat1.size() == 1);
	CHECK(cat2.size() == 1);

	cat1.erase(cif::key("id") == 30);
	CHECK(cat1.size() == 0);
	CHECK(cat2.size() == 0);
}

// --------------------------------------------------------------------

TEST_CASE("d6")
{
	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

save_cat_1
    _category.description     'A simple test category'
    _category.id              cat_1
    _category.mandatory_code  yes
    _category_key.name        '_cat_1.id'
    save_

save__cat_1.id
    _item.name                '_cat_1.id'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_1.id_2
    _item.name                '_cat_1.id_2'
    _item.category_id         cat_1
    _item.mandatory_code      no
    _item_type.code           int
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
    _item_type.code           int
    save_

save__cat_2.parent_id
    _item.name                '_cat_2.parent_id'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    save_

save__cat_2.parent_id_2
    _item.name                '_cat_2.parent_id_2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    save_

loop_
_pdbx_item_linked_group_list.child_category_id
_pdbx_item_linked_group_list.link_group_id
_pdbx_item_linked_group_list.child_name
_pdbx_item_linked_group_list.parent_name
_pdbx_item_linked_group_list.parent_category_id
cat_2 1 '_cat_2.parent_id'  '_cat_1.id' cat_1
cat_2 1 '_cat_2.parent_id_2' '_cat_1.id_2' cat_1

loop_
_pdbx_item_linked_group.category_id
_pdbx_item_linked_group.link_group_id
_pdbx_item_linked_group.label
cat_2 1 cat_2:cat_1:1
    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.id_2
1 1
2 2
3 ?

loop_
_cat_2.id
_cat_2.parent_id
_cat_2.parent_id_2
 0 1 1
 1 1 ?
 2 ? 1
 3 ? ?
 4 2 2
 5 3 1
 6 3 ?
    )";

	// --------------------------------------------------------------------

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	// auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];

	// f.front().validate_links();
	// CHECK(not );

	using namespace cif::literals;

	CHECK(cat2.has_parents(cat2.find1("id"_key == 0)));
	CHECK(cat2.has_parents(cat2.find1("id"_key == 1)));
	CHECK(cat2.has_parents(cat2.find1("id"_key == 2)));
	CHECK(not cat2.has_parents(cat2.find1("id"_key == 3)));
	CHECK(cat2.has_parents(cat2.find1("id"_key == 4)));
	CHECK(not cat2.has_parents(cat2.find1("id"_key == 5)));
	CHECK(cat2.has_parents(cat2.find1("id"_key == 6)));
}

// --------------------------------------------------------------------

TEST_CASE("c1")
{
	cif::VERBOSE = 1;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
4 .
5 ?
    )"_cf;

	auto &db = f.front();

	for (auto r : db["test"].find(cif::key("id") == 1))
	{
		const auto &[id, name] = r.get<int, std::string>("id", "name");
		CHECK(id == 1);
		CHECK(name == "aap");
	}

	for (auto r : db["test"].find(cif::key("id") == 4))
	{
		const auto &[id, name] = r.get<int, std::string>("id", "name");
		CHECK(id == 4);
		CHECK(name.empty());
	}

	for (auto r : db["test"].find(cif::key("id") == 5))
	{
		const auto &[id, name] = r.get<int, std::string>("id", "name");
		CHECK(id == 5);
		CHECK(name.empty());
	}

	// optional

	for (auto r : db["test"])
	{
		const auto &[id, name] = r.get<int, std::optional<std::string>>("id", "name");
		switch (id)
		{
			case 1: CHECK(name == "aap"); break;
			case 2: CHECK(name == "noot"); break;
			case 3: CHECK(name == "mies"); break;
			case 4:
			case 5: CHECK(not name); break;
			default:
				CHECK(false);
		}
	}
}

TEST_CASE("c2")
{
	cif::VERBOSE = 1;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
4 .
5 ?
    )"_cf;

	auto &db = f.front();

	// query tests

	for (const auto &[id, name] : db["test"].rows<int, std::optional<std::string>>("id", "name"))
	{
		switch (id)
		{
			case 1: CHECK(name == "aap"); break;
			case 2: CHECK(name == "noot"); break;
			case 3: CHECK(name == "mies"); break;
			case 4:
			case 5: CHECK(not name); break;
			default:
				CHECK(false);
		}
	}
}

TEST_CASE("c3")
{
	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
4 .
5 ?
    )"_cf;

	auto &db = f.front();

	// query tests
	for (const auto &[id, name] : db["test"].find<int, std::optional<std::string>>(cif::all(), "id", "name"))
	{
		switch (id)
		{
			case 1: CHECK(name == "aap"); break;
			case 2: CHECK(name == "noot"); break;
			case 3: CHECK(name == "mies"); break;
			case 4:
			case 5: CHECK(not name); break;
			default:
				CHECK(false);
		}
	}

	const auto &[id, name] = db["test"].find1<int, std::string>(cif::key("id") == 1, "id", "name");

	CHECK(id == 1);
	CHECK(name == "aap");
}

TEST_CASE("c4")
{
	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
4 .
5 ?
    )"_cf;

	auto &db = f.front();

	// query tests
	CHECK(db["test"].find_max<int>("id") == 5);
	CHECK(db["test"].find_max<int>("id", cif::key("name") != cif::null) == 3);

	CHECK(db["test"].find_min<int>("id") == 1);
	CHECK(db["test"].find_min<int>("id", cif::key("name") == cif::null) == 4);

	// count tests
	CHECK(db["test"].count(cif::all()) == 5);
	CHECK(db["test"].count(cif::key("name") != cif::null) == 3);
	CHECK(db["test"].count(cif::key("name") == cif::null) == 2);

	// find_first tests
	CHECK(db["test"].find_first<int>(cif::key("id") == 1, "id") == 1);
	CHECK(db["test"].find_first<int>(cif::all(), "id") == 1);

	std::optional<int> v;

	v = db["test"].find_first<std::optional<int>>(cif::key("id") == 1, "id");
	REQUIRE(v.has_value());
	CHECK(*v == 1);

	v = db["test"].find_first<std::optional<int>>(cif::key("id") == 6, "id");
	CHECK(not v.has_value());

	// find1 tests
	CHECK(db["test"].find1<int>(cif::key("id") == 1, "id") == 1);
	CHECK_THROWS_AS(db["test"].find1<int>(cif::all(), "id"), cif::multiple_results_error);
}

// --------------------------------------------------------------------
// rename test

TEST_CASE("r1")
{
	/*
	    Rationale:

	    The pdbx_mmcif dictionary contains inconsistent child-parent relations. E.g. atom_site is parent
	    of pdbx_nonpoly_scheme which itself is a parent of pdbx_entity_nonpoly. If I want to rename a residue
	    I cannot update pdbx_nonpoly_scheme since changing a parent changes children, but not vice versa.

	    But if I change the comp_id in atom_site, the pdbx_nonpoly_scheme is updated, that's good, and then
	    pdbx_entity_nonpoly is updated and that's bad.

	    The idea is now that if we update a parent and a child that must change as well, we first check
	    if there are more parents of this child that will not change. In that case we have to split the
	    child into two, one with the new value and one with the old. We then of course have to split all
	    children of this split row that are direct children.
	*/

	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

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
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           int
    save_

save__cat_1.name
    _item.name                '_cat_1.name'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_1.desc
    _item.name                '_cat_1.desc'
    _item.category_id         cat_1
    _item.mandatory_code      yes
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
    _item_type.code           int
    save_

save__cat_2.name
    _item.name                '_cat_2.name'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_2.num
    _item.name                '_cat_2.num'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_2.desc
    _item.name                '_cat_2.desc'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           text
    save_

save_cat_3
    _category.description     'A third simple test category'
    _category.id              cat_3
    _category.mandatory_code  no
    _category_key.name        '_cat_3.id'
    save_

save__cat_3.id
    _item.name                '_cat_3.id'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_3.name
    _item.name                '_cat_3.name'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_3.num
    _item.name                '_cat_3.num'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           int
    save_

loop_
_pdbx_item_linked_group_list.child_category_id
_pdbx_item_linked_group_list.link_group_id
_pdbx_item_linked_group_list.child_name
_pdbx_item_linked_group_list.parent_name
_pdbx_item_linked_group_list.parent_category_id
cat_1 1 '_cat_1.name' '_cat_2.name' cat_2
cat_2 1 '_cat_2.name' '_cat_3.name' cat_3
cat_2 1 '_cat_2.num'  '_cat_3.num'  cat_3

    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.name
_cat_1.desc
1 aap  Aap
2 noot Noot
3 mies Mies

loop_
_cat_2.id
_cat_2.name
_cat_2.num
_cat_2.desc
1 aap  1 'Een dier'
2 aap  2 'Een andere aap'
3 noot 1 'walnoot bijvoorbeeld'

loop_
_cat_3.id
_cat_3.name
_cat_3.num
1 aap 1
2 aap 2
    )";

	using namespace cif::literals;

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];
	auto &cat3 = f.front()["cat_3"];

	cat3.update_value("name"_key == "aap" and "num"_key == 1, "name", "aapje");

	CHECK(cat3.size() == 2);

	{
		int id, num;
		std::string name;
		cif::tie(id, name, num) = cat3.front().get("id", "name", "num");
		CHECK(id == 1);
		CHECK(num == 1);
		CHECK(name == "aapje");

		cif::tie(id, name, num) = cat3.back().get("id", "name", "num");
		CHECK(id == 2);
		CHECK(num == 2);
		CHECK(name == "aap");
	}

	int i = 0;
	for (const auto &[id, name, num, desc] : cat2.rows<int, std::string, int, std::string>("id", "name", "num", "desc"))
	{
		switch (++i)
		{
			case 1:
				CHECK(id == 1);
				CHECK(num == 1);
				CHECK(name == "aapje");
				CHECK(desc == "Een dier");
				break;

			case 2:
				CHECK(id == 2);
				CHECK(num == 2);
				CHECK(name == "aap");
				CHECK(desc == "Een andere aap");
				break;

			case 3:
				CHECK(id == 3);
				CHECK(num == 1);
				CHECK(name == "noot");
				CHECK(desc == "walnoot bijvoorbeeld");
				break;

			default:
				CHECK(false /*"Unexpected record"*/);
		}
	}

	CHECK(cat1.size() == 4);
	i = 0;
	for (const auto &[id, name, desc] : cat1.rows<int, std::string, std::string>("id", "name", "desc"))
	{
		switch (++i)
		{
			case 1:
				CHECK(id == 1);
				CHECK(name == "aapje");
				CHECK(desc == "Aap");
				break;

			case 2:
				CHECK(id == 2);
				CHECK(name == "noot");
				CHECK(desc == "Noot");
				break;

			case 3:
				CHECK(id == 3);
				CHECK(name == "mies");
				CHECK(desc == "Mies");
				break;

			case 4:
				CHECK(id == 4);
				CHECK(name == "aap");
				CHECK(desc == "Aap");
				break;

			default:
				CHECK(false /* "Unexpected record" */);
		}
	}

	// f.save(std::cout);
}

TEST_CASE("pc_1")
{
	/*
	    Parent/child tests

	    Note that the dictionary is different than the one in test r1
	*/

	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

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
    _item_linked.child_name   '_cat_2.parent_id'
    _item_linked.parent_name  '_cat_1.id'
    _item_type.code           int
    save_

save__cat_1.name
    _item.name                '_cat_1.name'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_1.desc
    _item.name                '_cat_1.desc'
    _item.category_id         cat_1
    _item.mandatory_code      yes
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
    _item_type.code           int
    save_

save__cat_2.name
    _item.name                '_cat_2.name'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_2.num
    _item.name                '_cat_2.num'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_2.desc
    _item.name                '_cat_2.desc'
    _item.category_id         cat_2
    _item.mandatory_code      yes
    _item_type.code           text
    save_

save_cat_3
    _category.description     'A third simple test category'
    _category.id              cat_3
    _category.mandatory_code  no
    _category_key.name        '_cat_3.id'
    save_

save__cat_3.id
    _item.name                '_cat_3.id'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_3.name
    _item.name                '_cat_3.name'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           code
    save_

save__cat_3.num
    _item.name                '_cat_3.num'
    _item.category_id         cat_3
    _item.mandatory_code      yes
    _item_type.code           int
    save_

loop_
_pdbx_item_linked_group_list.parent_category_id
_pdbx_item_linked_group_list.link_group_id
_pdbx_item_linked_group_list.parent_name
_pdbx_item_linked_group_list.child_name
_pdbx_item_linked_group_list.child_category_id
cat_1 1 '_cat_1.name' '_cat_2.name' cat_2
cat_2 1 '_cat_2.name' '_cat_3.name' cat_3
cat_2 1 '_cat_2.num'  '_cat_3.num'  cat_3

    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.name
_cat_1.desc
1 aap  Aap
2 noot Noot
3 mies Mies

loop_
_cat_2.id
_cat_2.name
_cat_2.num
_cat_2.desc
1 aap  1 'Een dier'
2 aap  2 'Een andere aap'
3 noot 1 'walnoot bijvoorbeeld'

loop_
_cat_3.id
_cat_3.name
_cat_3.num
1 aap 1
2 aap 2
    )";

	using namespace cif::literals;

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];
	auto &cat2 = f.front()["cat_2"];
	auto &cat3 = f.front()["cat_3"];

	// some parent/child tests

	// find all children in cat2 for the row with id == 1 in cat1
	auto rs1 = cat1.get_children(cat1.find1("id"_key == 1), cat2);
	CHECK(rs1.size() == 2);

	auto rs2 = cat1.get_children(cat1.find1("id"_key == 2), cat2);
	CHECK(rs2.size() == 1);

	auto rs3 = cat1.get_children(cat1.find1("id"_key == 3), cat2);
	CHECK(rs3.size() == 0);

	// finding parents
	auto rs4 = cat2.get_parents(cat2.find1("id"_key == 1), cat1);
	CHECK(rs4.size() == 1);

	auto rs5 = cat3.get_parents(cat3.find1("id"_key == 1), cat2);
	CHECK(rs5.size() == 1);

	// This link is not defined:
	auto rs6 = cat3.get_parents(cat3.find1("id"_key == 1), cat1);
	CHECK(rs6.size() == 0);
}

// --------------------------------------------------------------------

// TEST_CASE("bondmap_1")
// {
// 	cif::VERBOSE = 2;

// 	// sections taken from CCD compounds.cif
// 	auto components = R"(
// data_ASN
// loop_
// _chem_comp_bond.comp_id
// _chem_comp_bond.atom_id_1
// _chem_comp_bond.atom_id_2
// _chem_comp_bond.value_order
// _chem_comp_bond.pdbx_aromatic_flag
// _chem_comp_bond.pdbx_stereo_config
// _chem_comp_bond.pdbx_ordinal
// ASN N   CA   SING N N 1
// ASN N   H    SING N N 2
// ASN N   H2   SING N N 3
// ASN CA  C    SING N N 4
// ASN CA  CB   SING N N 5
// ASN CA  HA   SING N N 6
// ASN C   O    DOUB N N 7
// ASN C   OXT  SING N N 8
// ASN CB  CG   SING N N 9
// ASN CB  HB2  SING N N 10
// ASN CB  HB3  SING N N 11
// ASN CG  OD1  DOUB N N 12
// ASN CG  ND2  SING N N 13
// ASN ND2 HD21 SING N N 14
// ASN ND2 HD22 SING N N 15
// ASN OXT HXT  SING N N 16
// data_PHE
// loop_
// _chem_comp_bond.comp_id
// _chem_comp_bond.atom_id_1
// _chem_comp_bond.atom_id_2
// _chem_comp_bond.value_order
// _chem_comp_bond.pdbx_aromatic_flag
// _chem_comp_bond.pdbx_stereo_config
// _chem_comp_bond.pdbx_ordinal
// PHE N   CA  SING N N 1
// PHE N   H   SING N N 2
// PHE N   H2  SING N N 3
// PHE CA  C   SING N N 4
// PHE CA  CB  SING N N 5
// PHE CA  HA  SING N N 6
// PHE C   O   DOUB N N 7
// PHE C   OXT SING N N 8
// PHE CB  CG  SING N N 9
// PHE CB  HB2 SING N N 10
// PHE CB  HB3 SING N N 11
// PHE CG  CD1 DOUB Y N 12
// PHE CG  CD2 SING Y N 13
// PHE CD1 CE1 SING Y N 14
// PHE CD1 HD1 SING N N 15
// PHE CD2 CE2 DOUB Y N 16
// PHE CD2 HD2 SING N N 17
// PHE CE1 CZ  DOUB Y N 18
// PHE CE1 HE1 SING N N 19
// PHE CE2 CZ  SING Y N 20
// PHE CE2 HE2 SING N N 21
// PHE CZ  HZ  SING N N 22
// PHE OXT HXT SING N N 23
// data_PRO
// loop_
// _chem_comp_bond.comp_id
// _chem_comp_bond.atom_id_1
// _chem_comp_bond.atom_id_2
// _chem_comp_bond.value_order
// _chem_comp_bond.pdbx_aromatic_flag
// _chem_comp_bond.pdbx_stereo_config
// _chem_comp_bond.pdbx_ordinal
// PRO N   CA  SING N N 1
// PRO N   CD  SING N N 2
// PRO N   H   SING N N 3
// PRO CA  C   SING N N 4
// PRO CA  CB  SING N N 5
// PRO CA  HA  SING N N 6
// PRO C   O   DOUB N N 7
// PRO C   OXT SING N N 8
// PRO CB  CG  SING N N 9
// PRO CB  HB2 SING N N 10
// PRO CB  HB3 SING N N 11
// PRO CG  CD  SING N N 12
// PRO CG  HG2 SING N N 13
// PRO CG  HG3 SING N N 14
// PRO CD  HD2 SING N N 15
// PRO CD  HD3 SING N N 16
// PRO OXT HXT SING N N 17
// )"_cf;

// 	const std::filesystem::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
// 	mmcif::File file(example.string());
// 	mmcif::Structure structure(file);

// 	(void)file.isValid();

// 	mmcif::BondMap bm(structure);

// 	// Test the bonds of the first three residues, that's PRO A 1, ASN A 2, PHE A 3

// 	for (const auto &[compound, seqnr] : std::initializer_list<std::tuple<std::string, int>>{{"PRO", 1}, {"ASN", 2}, {"PHE", 3}})
// 	{
// 		auto &res = structure.get_residue("A", compound, seqnr, "");
// 		auto atoms = res.atoms();

// 		auto dc = components.get(compound);
// 		CHECK(dc != nullptr);

// 		auto cc = dc->get("chem_comp_bond");
// 		CHECK(cc != nullptr);

// 		std::set<std::tuple<std::string, std::string>> bonded;

// 		for (const auto &[atom_id_1, atom_id_2] : cc->rows<std::string, std::string>("atom_id_1", "atom_id_2"))
// 		{
// 			if (atom_id_1 > atom_id_2)
// 				bonded.insert({atom_id_2, atom_id_1});
// 			else
// 				bonded.insert({atom_id_1, atom_id_2});
// 		}

// 		for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
// 		{
// 			auto label_i = atoms[i].labelAtomID();

// 			for (std::size_t j = i + 1; j < atoms.size(); ++j)
// 			{
// 				auto label_j = atoms[j].labelAtomID();

// 				bool bonded_1 = bm(atoms[i], atoms[j]);
// 				bool bonded_1_i = bm(atoms[j], atoms[i]);

// 				bool bonded_t = label_i > label_j
// 				                    ? bonded.count({label_j, label_i})
// 				                    : bonded.count({label_i, label_j});

// 				CHECK(bonded_1 == bonded_t);
// 				CHECK(bonded_1_i == bonded_t);
// 			}
// 		}
// 	}

// 	// And check the inter-aminoacid links

// 	auto &poly = structure.polymers().front();

// 	for (std::size_t i = 0; i + 1 < poly.size(); ++i)
// 	{
// 		auto C = poly[i].atomByID("C");
// 		auto N = poly[i + 1].atomByID("N");

// 		CHECK(bm(C, N));
// 		CHECK(bm(N, C));
// 	}
// }

// TEST_CASE("bondmap_2")
// {
// 	CHECK_THROWS_AS(mmcif::BondMap::atomIDsForCompound("UN_"), mmcif::BondMapException);

// 	mmcif::CompoundFactory::instance().pushDictionary(gTestDir / "UN_.cif");

// 	CHECK(mmcif::BondMap::atomIDsForCompound("UN_").empty() == false);
// }

TEST_CASE("reading_file_1")
{
	std::istringstream is("Hello, world!");

	cif::file file;
	CHECK_THROWS_AS(file.load(is), std::runtime_error);
}

TEST_CASE("parser_test_1")
{
	auto data1 = R"(
data_QM
_test.text ??
)"_cf;

	auto &db1 = data1.front();
	auto &test1 = db1["test"];

	CHECK(test1.size() == 1);

	for (auto r : test1)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == "??");
	}

	std::stringstream ss;
	data1.save(ss);

	auto data2 = cif::file(ss);

	auto &db2 = data2.front();
	auto &test2 = db2["test"];

	CHECK(test2.size() == 1);

	for (auto r : test2)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == "??");
	}
}

TEST_CASE("output_test_1")
{
	auto data1 = R"(
data_Q
loop_
_test.text
stop_the_crap
'and stop_ this too'
'data_dinges'
blablaglobal_bla
boo.data_.whatever
'data_.whatever'
'stop_'
'loop_'
'global_'
'_with.leading_underscore'
)"_cf;

	auto &db1 = data1.front();
	auto &test1 = db1["test"];

	struct T
	{
		const char *s;
		bool q;
	} kS[] = {
		{ "stop_the_crap", true },
		{ "and stop_ this too", false },
		{ "data_dinges", false },
		{ "blablaglobal_bla", true },
		{ "boo.data_.whatever", true },
		{ "data_.whatever", false },
		{ "stop_", false },
		{ "loop_", false },
		{ "global_", false },
		{ "_with.leading_underscore", false }
	};

	CHECK(test1.size() == sizeof(kS) / sizeof(T));

	std::size_t i = 0;
	for (auto r : test1)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == kS[i].s);
		CHECK(cif::sac_parser::is_unquoted_string(kS[i].s) == kS[i].q);
		++i;
	}

	std::stringstream ss;
	data1.save(ss);

	auto data2 = cif::file(ss);

	auto &db2 = data2.front();
	auto &test2 = db2["test"];

	CHECK(test2.size() == sizeof(kS) / sizeof(T));

	i = 0;
	for (auto r : test2)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == kS[i++].s);
	}
}

TEST_CASE("output_test_2")
{
	auto data1 = R"(
data_Q
loop_
_test.text
;A very, very loooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooong line
;
;A line with a newline, look:
There it was!
;
)"_cf;

	auto &db1 = data1.front();
	auto &test1 = db1["test"];

	struct T
	{
		const char *s;
		bool q;
	} kS[] = {
		{ "A very, very loooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooong line", false },
		{ R"(A line with a newline, look:
There it was!)",
			false }
	};

	CHECK(test1.size() == sizeof(kS) / sizeof(T));

	std::size_t i = 0;
	for (auto r : test1)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == kS[i].s);
		CHECK(cif::sac_parser::is_unquoted_string(kS[i].s) == kS[i].q);
		++i;
	}

	std::stringstream ss;
	data1.save(ss);

	auto data2 = cif::file(ss);

	auto &db2 = data2.front();
	auto &test2 = db2["test"];

	CHECK(test2.size() == sizeof(kS) / sizeof(T));

	i = 0;
	for (auto r : test2)
	{
		auto text = r.get<std::string>("text");
		CHECK(text == kS[i++].s);
	}
}

TEST_CASE("trim_test")
{
	CHECK(cif::trim_copy("aap") == "aap");
	CHECK(cif::trim_copy(" aap") == "aap");
	CHECK(cif::trim_copy(" aap ") == "aap");
	CHECK(cif::trim_copy("aap ") == "aap");
	CHECK(cif::trim_copy("	 aap	") == "aap");

	CHECK(cif::trim_left_copy("aap") == "aap");
	CHECK(cif::trim_left_copy(" aap") == "aap");
	CHECK(cif::trim_left_copy(" aap ") == "aap ");
	CHECK(cif::trim_left_copy("aap ") == "aap ");
	CHECK(cif::trim_left_copy("aap	") == "aap	");

	CHECK(cif::trim_right_copy("aap") == "aap");
	CHECK(cif::trim_right_copy(" aap") == " aap");
	CHECK(cif::trim_right_copy(" aap ") == " aap");
	CHECK(cif::trim_right_copy("aap ") == "aap");
	CHECK(cif::trim_right_copy("	 aap	") == "	 aap");

	std::string s;

	s = "aap";
	cif::trim(s);
	CHECK(s == "aap");
	s = " aap";
	cif::trim(s);
	CHECK(s == "aap");
	s = " aap ";
	cif::trim(s);
	CHECK(s == "aap");
	s = "aap ";
	cif::trim(s);
	CHECK(s == "aap");
	s = "	 aap	";
	cif::trim(s);
	CHECK(s == "aap");

	s = "aap";
	cif::trim_left(s);
	CHECK(s == "aap");
	s = " aap";
	cif::trim_left(s);
	CHECK(s == "aap");
	s = " aap ";
	cif::trim_left(s);
	CHECK(s == "aap ");
	s = "aap ";
	cif::trim_left(s);
	CHECK(s == "aap ");
	s = "aap	";
	cif::trim_left(s);
	CHECK(s == "aap	");

	s = "aap";
	cif::trim_right(s);
	CHECK(s == "aap");
	s = " aap";
	cif::trim_right(s);
	CHECK(s == " aap");
	s = " aap ";
	cif::trim_right(s);
	CHECK(s == " aap");
	s = "aap ";
	cif::trim_right(s);
	CHECK(s == "aap");
	s = "	 aap	";
	cif::trim_right(s);
	CHECK(s == "	 aap");
}

TEST_CASE("split_test")
{
	std::vector<std::string_view> v, t;

	v = cif::split<>("aap;noot;mies", ";");
	t = std::vector<std::string_view>{ "aap", "noot", "mies" };

	CHECK(v == t);

	v = cif::split("aap;noot,mies", ";,");
	// t = std::vector<std::string>{ "aap", "noot", "mies" };

	CHECK(v == t);

	v = cif::split(";aap;noot,mies;", ";,");
	t = std::vector<std::string_view>{ "", "aap", "noot", "mies", "" };

	CHECK(v == t);

	v = cif::split(";aap;noot,mies;", ";,", true);
	t = std::vector<std::string_view>{ "aap", "noot", "mies" };

	CHECK(v == t);
}

TEST_CASE("join_test")
{
	CHECK(cif::join(std::vector<std::string>{ "aap" }, ", ") == "aap");
	CHECK(cif::join(std::vector<std::string>{ "aap", "noot" }, ", ") == "aap, noot");
	CHECK(cif::join(std::vector<std::string>{ "aap", "noot", "mies" }, ", ") == "aap, noot, mies");
}

TEST_CASE("replace_all_test")
{
	std::string s("aap, noot, mies");
	cif::replace_all(s, ", ", ",");
	CHECK(s == "aap,noot,mies");

	cif::replace_all(s, ",", ", ");
	CHECK(s == "aap, noot, mies");

	cif::replace_all(s, ", ", ", ");
	CHECK(s == "aap, noot, mies");
}

// --------------------------------------------------------------------

TEST_CASE("reorder_test")
{

	const char dict[] = R"(
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
    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.name
2 Noot
1 Aap
3 Mies
    )";

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	CHECK(f.is_valid());

	auto &cat1 = f.front()["cat_1"];
	cat1.reorder_by_index();

	int n = 1;

	const char *ts[] = { "Aap", "Noot", "Mies" };

	for (const auto &[id, name] : cat1.rows<int, std::string>("id", "name"))
	{
		CHECK(id == n);
		CHECK(name == ts[n - 1]);
		++n;
	}
}

// --------------------------------------------------------------------

TEST_CASE("audit_conform_test")
{

	const char dict[] = R"(
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

###################
## AUDIT_CONFORM ##
###################

save_audit_conform
    _category.description
;              Data items in the AUDIT_CONFORM category describe the
               dictionary versions against which the data names appearing in
               the current data block are conformant.
;
    _category.id                  audit_conform
    _category.mandatory_code      no
    loop_
    _category_key.name          '_audit_conform.dict_name'
                                '_audit_conform.dict_version'
    loop_
    _category_group.id           'inclusive_group'
                                 'audit_group'
    loop_
    _category_examples.detail
    _category_examples.case
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;
    Example 1 - any file conforming to the current CIF core dictionary.
;
;
    _audit_conform.dict_name         cif_core.dic
    _audit_conform.dict_version      2.3.1
    _audit_conform.dict_location
                         ftp://ftp.iucr.org/pub/cif_core.2.3.1.dic
;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     save_

save__audit_conform.dict_location
    _item_description.description
;              A file name or uniform resource locator (URL) for the
               dictionary to which the current data block conforms.
;
    _item.name                  '_audit_conform.dict_location'
    _item.category_id             audit_conform
    _item.mandatory_code          no
    _item_aliases.alias_name    '_audit_conform_dict_location'
    _item_aliases.dictionary      cif_core.dic
    _item_aliases.version         2.0.1
    _item_type.code               text
     save_

save__audit_conform.dict_name
    _item_description.description
;              The string identifying the highest-level dictionary defining
               data names used in this file.
;
    _item.name                  '_audit_conform.dict_name'
    _item.category_id             audit_conform
    _item.mandatory_code          yes
    _item_aliases.alias_name    '_audit_conform_dict_name'
    _item_aliases.dictionary      cif_core.dic
    _item_aliases.version         2.0.1
    _item_type.code               text
     save_

save__audit_conform.dict_version
    _item_description.description
;              The version number of the dictionary to which the current
               data block conforms.
;
    _item.name                  '_audit_conform.dict_version'
    _item.category_id             audit_conform
    _item.mandatory_code          yes
    _item_aliases.alias_name    '_audit_conform_dict_version'
    _item_aliases.dictionary      cif_core.dic
    _item_aliases.version         2.0.1
    _item_type.code               text
     save_

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
    )";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	auto &validator = cif::validator_factory::instance().add(cif::validator(is_dict));

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
#
_audit_conform.dict_name test_dict.dic
_audit_conform.dict_version 1.0
#
loop_
_cat_1.id
_cat_1.name
2 Noot
1 Aap
3 Mies
    )";

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	CHECK(f.is_valid());

	std::stringstream ss;
	ss << f;

	cif::file f2(ss);
	CHECK(f2.empty() == false);
	f2.front().load_dictionary();
	CHECK(f2.is_valid());

	auto &audit_conform = f2.front()["audit_conform"];
	CHECK(audit_conform.front()["dict_name"].get<std::string>() == "test_dict.dic");
	CHECK(audit_conform.front()["dict_version"].get<float>() == 1.0);
}

// --------------------------------------------------------------------

TEST_CASE("ix_op_1")
{
	const char dict[] = R"(
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
               code      char
               '[][_,.;:"&<>()/\{}'`~!@#$%A-Za-z0-9*|+-]*'

               text      char
               '[][ \n\t()_,.;:"&<>/\{}'`~!@#$%?+=*A-Za-z0-9|^-]*'

               int       numb
               '[+-]?[0-9]+'

save_cat_1
    _category.description     'A simple test category'
    _category.id              cat_1
    _category.mandatory_code  yes
	loop_
	_category_key.name        '_cat_1.id'
	                          '_cat_1.id_2'
    save_

save__cat_1.id
    _item.name                '_cat_1.id'
    _item.category_id         cat_1
    _item.mandatory_code      yes
    _item_type.code           int
    save_

save__cat_1.id_2
    _item.name                '_cat_1.id_2'
    _item.category_id         cat_1
    _item.mandatory_code      no
    _item_type.code           int
    save_
)";

	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(dict), sizeof(dict) - 1);

	std::istream is_dict(&buffer);

	cif::validator validator(is_dict);

	cif::file f;

	// --------------------------------------------------------------------

	const char data[] = R"(
data_test
loop_
_cat_1.id
_cat_1.id_2
1 10
2 20
3 ?
    )";

	// --------------------------------------------------------------------

	struct data_membuf : public std::streambuf
	{
		data_membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} data_buffer(const_cast<char *>(data), sizeof(data) - 1);

	std::istream is_data(&data_buffer);
	f.load(is_data, validator);

	auto &cat1 = f.front()["cat_1"];

	using key_type = cif::category::key_type;
	using test_tuple_type = std::tuple<key_type, bool>;

	test_tuple_type TESTS[] = {
		{ { { "id", 1 }, { "id_2", 10 } }, true },
		{ { { "id_2", 10 }, { "id", 1 } }, true },
		{ { { "id", 1 }, { "id_2", 20 } }, false },
		{ { { "id", 3 }, { "id_2", nullptr } }, true },
	};

	for (const auto &[key, test] : TESTS)
		CHECK((bool)cat1[key] == test);
}

// --------------------------------------------------------------------

TEST_CASE("cifv1_0_1")
{
	auto f = R"(data_TEST
#
loop_
_id
_name
1 aap
2 noot
3 mies
4 ?
5 .
    )"_cf;

	auto &db = f.front();

	auto &cat = db[""];

	for (auto r : cat)
	{
		int id;
		std::optional<std::string> name;

		cif::tie(id, name) = r.get("id", "name");

		switch (id)
		{
			case 1: REQUIRE(name.has_value()); CHECK(*name == "aap"); break;
			case 2: REQUIRE(name.has_value()); CHECK(*name == "noot"); break;
			case 3: REQUIRE(name.has_value()); CHECK(*name == "mies"); break;
			default: CHECK(name.has_value() == false);
		}
	}

	std::stringstream ss;
	ss << db;

	auto f2 = cif::file(ss);
	auto &db2 = f2.front();

	CHECK(db == db2);
}

// TEST_CASE("cifv1_0_2")
// {
// 	CHECK_THROWS_AS(R"(data_TEST
// #
// _version 1.0
// loop_
// _id
// _name
// 1 aap
// 2 noot
// 3 mies
// 4 ?
// 5 .
//     )"_cf, cif::parse_error);
// }

TEST_CASE("cifv1_0_3")
{
	auto f = R"(data_TEST
#
_version 1.0
_date    today
)"_cf;

	auto &db = f.front();

	auto &cat = db[""];
	CHECK(not cat.empty());

	auto r = cat.front();
	CHECK(r["version"].get<float>() == 1.0);
	CHECK(r["date"].get<std::string>() == "today");

	std::stringstream ss;
	ss << db;

	auto f2 = cif::file(ss);
	auto &db2 = f2.front();

	CHECK(db == db2);
}

TEST_CASE("find1_opt_1")
{
	using namespace cif::literals;
	using namespace std::literals;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
_test.value
1 aap   1.0
2 noot  1.1
3 mies  1.2
    )"_cf;

	auto &db = f.front();
	auto &test = db["test"];

	auto v = test.find1<std::optional<float>>("id"_key == 1, "value");
	CHECK(v.has_value());
	CHECK(*v == 1.0f);

	v = test.find1<std::optional<float>>("id"_key == 4, "value");
	CHECK(v.has_value() == false);
}

// // --------------------------------------------------------------------

// TEST_CASE("compound_test_1")
// {
// 	cif::compound_factory::instance().push_dictionary(gTestDir / "REA_v2.cif");
// 	auto compound = cif::compound_factory::instance().create("REA_v2");
// 	CHECK(compound != nullptr);
// 	CHECK(cif::iequals(compound->id(), "REA_v2"));
// }

// // --------------------------------------------------------------------

// TEST_CASE("pdb_parser_test_1")
// {
// 	char k1CBS[] = R"(HEADER    RETINOIC-ACID TRANSPORT                 28-SEP-94   1CBS
// TITLE     CRYSTAL STRUCTURE OF CELLULAR RETINOIC-ACID-BINDING
// TITLE    2 PROTEINS I AND II IN COMPLEX WITH ALL-TRANS-RETINOIC ACID
// TITLE    3 AND A SYNTHETIC RETINOID
// COMPND    MOL_ID: 1;
// COMPND   2 MOLECULE: CELLULAR RETINOIC ACID BINDING PROTEIN TYPE II;
// COMPND   3 CHAIN: A;
// COMPND   4 ENGINEERED: YES
// SOURCE    MOL_ID: 1;
// SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;
// SOURCE   3 ORGANISM_COMMON: HUMAN;
// SOURCE   4 ORGANISM_TAXID: 9606;
// SOURCE   5 CELL_LINE: BL21;
// SOURCE   6 GENE: HUMAN CRABP-II;
// SOURCE   7 EXPRESSION_SYSTEM: ESCHERICHIA COLI BL21(DE3);
// SOURCE   8 EXPRESSION_SYSTEM_TAXID: 469008;
// SOURCE   9 EXPRESSION_SYSTEM_STRAIN: BL21 (DE3);
// SOURCE  10 EXPRESSION_SYSTEM_PLASMID: PET-3A
// KEYWDS    RETINOIC-ACID TRANSPORT
// EXPDTA    X-RAY DIFFRACTION
// AUTHOR    G.J.KLEYWEGT,T.BERGFORS,T.A.JONES
// ATOM      1  N   PRO A   1      16.979  13.301  44.555  1.00 30.05           N
// ATOM      2  CA  PRO A   1      18.150  13.525  43.680  1.00 28.82           C
// ATOM      3  C   PRO A   1      18.656  14.966  43.784  1.00 26.59           C
// ATOM      4  O   PRO A   1      17.890  15.889  44.078  1.00 26.84           O
// ATOM      5  CB  PRO A   1      17.678  13.270  42.255  1.00 29.24           C
// ATOM      6  CG  PRO A   1      16.248  13.734  42.347  1.00 29.29           C
// ATOM      7  CD  PRO A   1      15.762  13.216  43.724  1.00 30.71           C)";

// 	struct membuf : public std::streambuf
// 	{
// 		membuf(char *text, std::size_t length)
// 		{
// 			this->setg(text, text, text + length);
// 		}
// 	} buffer(k1CBS, sizeof(k1CBS) - 1);

// 	std::istream is(&buffer);

// 	auto f = cif::pdb::read(is);
// }

// // --------------------------------------------------------------------

// TEST_CASE("compound_not_found_test_1")
// {
// 	auto cmp = cif::compound_factory::instance().create("&&&");
// 	CHECK(cmp == nullptr);
// }

// // --------------------------------------------------------------------
// // PDB2CIF tests

// TEST_CASE("pdb2cif_formula_weight")
// {
// 	cif::compound_factory::instance().push_dictionary(gTestDir / "REA.cif");

// 	cif::file a = cif::pdb::read(gTestDir / "pdb1cbs.ent.gz");

// 	auto fw = a.front()["entity"].find1<float>(cif::key("id") == 1, "formula_weight");
// 	CHECK(std::abs(fw - 15581.802f) < 0.1f);

// 	fw = a.front()["entity"].find1<float>(cif::key("id") == 2, "formula_weight");
// 	CHECK(fw == 300.435f);

// 	fw = a.front()["entity"].find1<float>(cif::key("id") == 3, "formula_weight");
// 	CHECK(fw == 18.015f);
// }

// // --------------------------------------------------------------------

// TEST_CASE("update_values_with_provider")
// {

// }
