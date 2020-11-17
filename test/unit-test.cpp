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
    cif::VERBOSE = 1;

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

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(d3)
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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];
    auto& cat2 = f.firstDatablock()["cat_2"];

    // check a rename in parent and child

    for (auto r: cat1.find(cif::Key("id") == 1))
    {
        r["id"] = 10;
        break;
    }

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 4);

    BOOST_CHECK(cat1.find(cif::Key("id") == 1).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 10).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 1).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 10).size() == 2);

    // check a rename in parent and child, this time only one child should be renamed

    for (auto r: cat1.find(cif::Key("id") == 2))
    {
        r["id"] = 20;
        break;
    }

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 4);

    BOOST_CHECK(cat1.find(cif::Key("id") == 2).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 20).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 2).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 20).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 2 and cif::Key("name2") == "noot").size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 2 and cif::Key("name2") == "n2").size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 20 and cif::Key("name2") == "noot").size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 20 and cif::Key("name2") == "n2").size() == 0);



    // // --------------------------------------------------------------------
    
    // cat1.erase(cif::Key("id") == 10);

    // BOOST_CHECK(cat1.size() == 2);
    // BOOST_CHECK(cat2.size() == 2);

    // cat1.erase(cif::Key("id") == 20);

    // BOOST_CHECK(cat1.size() == 1);
    // BOOST_CHECK(cat2.size() == 1);



}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(d4)
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
    _item_type.code           int
    save_

save__cat_2.parent_id2
    _item.name                '_cat_2.parent_id2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    _item_type.code           code
    save_

save__cat_2.parent_id3
    _item.name                '_cat_2.parent_id3'
    _item.category_id         cat_2
    _item.mandatory_code      no
    _item_type.code           code
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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];
    auto& cat2 = f.firstDatablock()["cat_2"];

    // check a rename in parent and child

    for (auto r: cat1.find(cif::Key("id") == 1))
    {
        r["id"] = 10;
        break;
    }

    BOOST_CHECK(cat1.size() == 4);
    BOOST_CHECK(cat2.size() == 13);

    BOOST_CHECK(cat1.find(cif::Key("id") == 1).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 10).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 1).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 10).size() == 2);


    for (auto r: cat1.find(cif::Key("id") == 2))
    {
        r["id"] = 20;
        break;
    }

    BOOST_CHECK(cat1.size() == 4);
    BOOST_CHECK(cat2.size() == 13);

    BOOST_CHECK(cat1.find(cif::Key("id") == 2).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 20).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 2).size() == 2);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 20).size() == 2);


    for (auto r: cat1.find(cif::Key("id") == 3))
    {
        r["id"] = 30;
        break;
    }

    BOOST_CHECK(cat1.size() == 4);
    BOOST_CHECK(cat2.size() == 13);

    BOOST_CHECK(cat1.find(cif::Key("id") == 3).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 30).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 3).size() == 2);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 30).size() == 1);


    for (auto r: cat1.find(cif::Key("id") == 4))
    {
        r["id"] = 40;
        break;
    }

    BOOST_CHECK(cat1.size() == 4);
    BOOST_CHECK(cat2.size() == 13);

    BOOST_CHECK(cat1.find(cif::Key("id") == 4).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 10).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 4).size() == 3);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 40).size() == 0);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(d5)
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
    _item_type.code           int
    save_

save__cat_2.parent_id2
    _item.name                '_cat_2.parent_id2'
    _item.category_id         cat_2
    _item.mandatory_code      no
    _item_type.code           code
    save_

save__cat_2.parent_id3
    _item.name                '_cat_2.parent_id3'
    _item.category_id         cat_2
    _item.mandatory_code      no
    _item_type.code           code
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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];
    auto& cat2 = f.firstDatablock()["cat_2"];

    // --------------------------------------------------------------------
    // check iterate children

    auto PR2set = cat1.find(cif::Key("id") == 2);
    BOOST_ASSERT(PR2set.size() == 1);
    auto PR2 = PR2set.front();
    BOOST_CHECK(PR2["id"].as<int>() == 2);

    auto CR2set = cat1.getChildren(PR2, "cat_2");
    BOOST_ASSERT(CR2set.size() == 3);

    std::vector<int> CRids;
    std::transform(CR2set.begin(), CR2set.end(), std::back_inserter(CRids), [](cif::Row r) { return r["id"].as<int>(); });
    std::sort(CRids.begin(), CRids.end());
    BOOST_CHECK(CRids == std::vector<int>({ 4, 5, 6}));

    // check a rename in parent and child

    for (auto r: cat1.find(cif::Key("id") == 1))
    {
        r["id"] = 10;
        break;
    }

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 7);

    BOOST_CHECK(cat1.find(cif::Key("id") == 1).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 10).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 1).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 1).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 1).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 10).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 10).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 10).size() == 1);

    for (auto r: cat1.find(cif::Key("id") == 2))
    {
        r["id"] = 20;
        break;
    }

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 7);

    BOOST_CHECK(cat1.find(cif::Key("id") == 2).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 20).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 2).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 2).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 2).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 20).size() == 2);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 20).size() == 2);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 20).size() == 2);

    for (auto r: cat1.find(cif::Key("id") == 3))
    {
        r["id"] = 30;
        break;
    }

    BOOST_CHECK(cat1.size() == 3);
    BOOST_CHECK(cat2.size() == 7);

    BOOST_CHECK(cat1.find(cif::Key("id") == 3).size() == 0);
    BOOST_CHECK(cat1.find(cif::Key("id") == 30).size() == 1);

    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 3).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 3).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 3).size() == 0);
    BOOST_CHECK(cat2.find(cif::Key("parent_id") == 30).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id2") == 30).size() == 1);
    BOOST_CHECK(cat2.find(cif::Key("parent_id3") == 30).size() == 1);

    // test delete

    cat1.erase(cif::Key("id") == 10);
    BOOST_CHECK(cat1.size() == 2);
    BOOST_CHECK(cat2.size() == 4);

    cat1.erase(cif::Key("id") == 20);
    BOOST_CHECK(cat1.size() == 1);
    BOOST_CHECK(cat2.size() == 1);

    cat1.erase(cif::Key("id") == 30);
    BOOST_CHECK(cat1.size() == 0);
    BOOST_CHECK(cat2.size() == 0);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(c1)
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

    auto& db = f.firstDatablock();

    for (auto r: db["test"].find(cif::Key("id") == 1))
    {
        const auto& [id, name] = r.get<int, std::string>({"id", "name"});
        BOOST_CHECK(id == 1);
        BOOST_CHECK(name == "aap");
    }

    for (auto r: db["test"].find(cif::Key("id") == 4))
    {
        const auto& [id, name] = r.get<int, std::string>({"id", "name"});
        BOOST_CHECK(id == 4);
        BOOST_CHECK(name.empty());
    }

    for (auto r: db["test"].find(cif::Key("id") == 5))
    {
        const auto& [id, name] = r.get<int, std::string>({"id", "name"});
        BOOST_CHECK(id == 5);
        BOOST_CHECK(name.empty());
    }

    // optional

    for (auto r: db["test"])
    {
        const auto& [id, name] = r.get<int, std::optional<std::string>>({"id", "name"});
        switch (id)
        {
            case 1: BOOST_CHECK(name == "aap"); break;
            case 2: BOOST_CHECK(name == "noot"); break;
            case 3: BOOST_CHECK(name == "mies"); break;
            case 4:
            case 5: BOOST_CHECK(not name); break;
            default: 
                    BOOST_CHECK(false);
        }
    }
}

BOOST_AUTO_TEST_CASE(c2)
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

    auto& db = f.firstDatablock();

    // query tests

    for (const auto& [id, name]: db["test"].rows<int, std::optional<std::string>>({ "id", "name" }))
    {
        switch (id)
        {
            case 1: BOOST_CHECK(name == "aap"); break;
            case 2: BOOST_CHECK(name == "noot"); break;
            case 3: BOOST_CHECK(name == "mies"); break;
            case 4:
            case 5: BOOST_CHECK(not name); break;
            default: 
                    BOOST_CHECK(false);
        }
    }
}


BOOST_AUTO_TEST_CASE(c3)
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

    auto& db = f.firstDatablock();

    // query tests
    for (const auto& [id, name]: db["test"].find<int, std::optional<std::string>>(cif::All(), { "id", "name" }))
    {
        switch (id)
        {
            case 1: BOOST_CHECK(name == "aap"); break;
            case 2: BOOST_CHECK(name == "noot"); break;
            case 3: BOOST_CHECK(name == "mies"); break;
            case 4:
            case 5: BOOST_CHECK(not name); break;
            default: 
                    BOOST_CHECK(false);
        }
    }

    const auto& [id, name] = db["test"].find1<int, std::string>(cif::Key("id") == 1, { "id", "name" });

    BOOST_CHECK(id == 1);
    BOOST_CHECK(name == "aap");
}
