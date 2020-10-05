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

#define BOOST_TEST_MODULE LibCifPP_Test
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>

// #include "cif++/DistanceMap.hpp"
#include "cif++/Cif++.hpp"

// --------------------------------------------------------------------

cif::File operator""_cf(const char* text, size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::File(is);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ut1)
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

    auto& db = f.firstDatablock();

    BOOST_CHECK(db.getName() == "TEST");
    
    auto& test = db["test"];
    BOOST_CHECK(test.size() == 3);

    // wrong! the next lines will crash. And that's OK, don't do that
    // for (auto r: test)
    // 	test.erase(r);
    
    // BOOST_CHECK(test.empty());

    // test.purge();

    auto n = test.erase(cif::Key("id") == 1, [](const cif::Row& r) {
        BOOST_CHECK_EQUAL(r["id"].as<int>(), 1);
        BOOST_CHECK_EQUAL(r["name"].as<std::string>(), "aap");
      });

    BOOST_CHECK_EQUAL(n, 1);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ut2)
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

    auto& db = f.firstDatablock();

    BOOST_CHECK(db.getName() == "TEST");
    
    auto& test = db["test"];
    BOOST_CHECK(test.size() == 3);

    int n = 0;
    for (auto r: test.find(cif::Key("name") == "aap"))
    {
        BOOST_CHECK(++n == 1);
        BOOST_CHECK(r["id"].as<int>() == 1);
        BOOST_CHECK(r["name"].as<std::string>() == "aap");
        BOOST_CHECK(r["value"].as<float>() == 1.0);
    }

    auto t = test.find(cif::Key("id") == 1);
    BOOST_CHECK(not t.empty());
    BOOST_CHECK(t.front()["name"].as<std::string>() == "aap");

    auto t2 = test.find(cif::Key("value") == 1.2);
    BOOST_CHECK(not t2.empty());
    BOOST_CHECK(t2.front()["name"].as<std::string>() == "mies");
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(d1)
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
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(dict), sizeof(dict) - 1);

    std::istream is_dict(&buffer);

    cif::File f;
    f.loadDictionary(is_dict);

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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];
    auto& cat2 = f.firstDatablock()["cat_2"];

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 3);

    cat1.erase(cif::Key("id") == 1);

    BOOST_CHECK(cat1.size() == 2);
    BOOST_CHECK(cat2.size() == 1);

    // BOOST_CHECK_THROW(cat2.emplace({
    //     { "id", 4 },
    //     { "parent_id", 4 },
    //     { "desc", "moet fout gaan" }
    // }), std::exception);

    BOOST_CHECK_THROW(cat2.emplace({
        { "id", "vijf" },   // <- invalid value
        { "parent_id", 2 },
        { "desc", "moet fout gaan" }
    }), std::exception);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(d2)
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
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(dict), sizeof(dict) - 1);

    std::istream is_dict(&buffer);

    cif::File f;
    f.loadDictionary(is_dict);

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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];

    BOOST_CHECK(cat1.size() == 3);

    cat1.erase(cif::Key("id") == "AAP");

    BOOST_CHECK(cat1.size() == 3);

    cat1.erase(cif::Key("id") == "noot");

    BOOST_CHECK(cat1.size() == 2);


}