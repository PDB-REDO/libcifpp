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
#include "cif++/BondMap.hpp"

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

	// do this now, avoids the need for installing
	cif::addFileResource("mmcif_pdbx_v50.dic", "../rsrc/mmcif_pdbx_v50.dic");

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

// --------------------------------------------------------------------
// rename test

BOOST_AUTO_TEST_CASE(r1)
{
	/*
		Rationale:

		The pdbx_mmcif dictionary contains inconsistent child-parent relations. E.g. atom_site is parent
		of pdbx_nonpoly_scheme which itself is a parent of pdbx_entity_nonpoly. If I want to rename a residue
		I cannot update pdbx_nonpoly_scheme since changing a parent changes children, but not vice versa.

		But if I change the comp_id in atom_site, the pdbx_nonpoly_scheme is update, that's good, and then
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
        data_membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } data_buffer(const_cast<char*>(data), sizeof(data) - 1);

    std::istream is_data(&data_buffer);
    f.load(is_data);

    auto& cat1 = f.firstDatablock()["cat_1"];
    auto& cat2 = f.firstDatablock()["cat_2"];
	auto& cat3 = f.firstDatablock()["cat_3"];

	cat3.update_value("name"_key == "aap" and "num"_key == 1, "name", "aapje");

	BOOST_CHECK(cat3.size() == 2);
	
	int id, num;
	std::string name;
	cif::tie(id, name, num) = cat3.front().get("id", "name", "num");
	BOOST_CHECK(id == 1);
	BOOST_CHECK(num == 1);
	BOOST_CHECK(name == "aapje");

	cif::tie(id, name, num) = cat3.back().get("id", "name", "num");
	BOOST_CHECK(id == 2);
	BOOST_CHECK(num == 2);
	BOOST_CHECK(name == "aap");
	
	int i = 0;	
	for (const auto &[id, name, num, desc]: cat2.rows<int,std::string,int,std::string>({"id", "name", "num", "desc"}))
	{
		switch (++i)
		{
			case 1:
				BOOST_CHECK(id == 1);
				BOOST_CHECK(num == 1);
				BOOST_CHECK(name == "aapje");
				BOOST_CHECK(desc == "Een dier");
				break;

			case 2:
				BOOST_CHECK(id == 2);
				BOOST_CHECK(num == 2);
				BOOST_CHECK(name == "aap");
				BOOST_CHECK(desc == "Een andere aap");
				break;

			case 3:
				BOOST_CHECK(id == 3);
				BOOST_CHECK(num == 1);
				BOOST_CHECK(name == "noot");
				BOOST_CHECK(desc == "walnoot bijvoorbeeld");
				break;
			
			default:
				BOOST_FAIL("Unexpected record");
		}
	}

	BOOST_CHECK(cat1.size() == 4);
	i = 0;
	for (const auto &[id, name, desc]: cat1.rows<int,std::string,std::string>({"id", "name", "desc"}))
	{
		switch (++i)
		{
			case 1:
				BOOST_CHECK(id == 1);
				BOOST_CHECK(name == "aapje");
				BOOST_CHECK(desc == "Aap");
				break;

			case 2:
				BOOST_CHECK(id == 2);
				BOOST_CHECK(name == "noot");
				BOOST_CHECK(desc == "Noot");
				break;

			case 3:
				BOOST_CHECK(id == 3);
				BOOST_CHECK(name == "mies");
				BOOST_CHECK(desc == "Mies");
				break;
			
			case 4:
				BOOST_CHECK(id == 4);
				BOOST_CHECK(name == "aap");
				BOOST_CHECK(desc == "Aap");
				break;
			
			default:
				BOOST_FAIL("Unexpected record");
		}
	}

	f.save(std::cout);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bondmap_1)
{
    cif::VERBOSE = 2;

	cif::addFileResource("components.cif", "../data/components.cif");

	// sections taken from CCD compounds.cif
	auto components = R"(
data_ASN
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
ASN N   CA   SING N N 1
ASN N   H    SING N N 2
ASN N   H2   SING N N 3
ASN CA  C    SING N N 4
ASN CA  CB   SING N N 5
ASN CA  HA   SING N N 6
ASN C   O    DOUB N N 7
ASN C   OXT  SING N N 8
ASN CB  CG   SING N N 9
ASN CB  HB2  SING N N 10
ASN CB  HB3  SING N N 11
ASN CG  OD1  DOUB N N 12
ASN CG  ND2  SING N N 13
ASN ND2 HD21 SING N N 14
ASN ND2 HD22 SING N N 15
ASN OXT HXT  SING N N 16
data_PHE
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
PHE N   CA  SING N N 1
PHE N   H   SING N N 2
PHE N   H2  SING N N 3
PHE CA  C   SING N N 4
PHE CA  CB  SING N N 5
PHE CA  HA  SING N N 6
PHE C   O   DOUB N N 7
PHE C   OXT SING N N 8
PHE CB  CG  SING N N 9
PHE CB  HB2 SING N N 10
PHE CB  HB3 SING N N 11
PHE CG  CD1 DOUB Y N 12
PHE CG  CD2 SING Y N 13
PHE CD1 CE1 SING Y N 14
PHE CD1 HD1 SING N N 15
PHE CD2 CE2 DOUB Y N 16
PHE CD2 HD2 SING N N 17
PHE CE1 CZ  DOUB Y N 18
PHE CE1 HE1 SING N N 19
PHE CE2 CZ  SING Y N 20
PHE CE2 HE2 SING N N 21
PHE CZ  HZ  SING N N 22
PHE OXT HXT SING N N 23
data_PRO
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
PRO N   CA  SING N N 1
PRO N   CD  SING N N 2
PRO N   H   SING N N 3
PRO CA  C   SING N N 4
PRO CA  CB  SING N N 5
PRO CA  HA  SING N N 6
PRO C   O   DOUB N N 7
PRO C   OXT SING N N 8
PRO CB  CG  SING N N 9
PRO CB  HB2 SING N N 10
PRO CB  HB3 SING N N 11
PRO CG  CD  SING N N 12
PRO CG  HG2 SING N N 13
PRO CG  HG3 SING N N 14
PRO CD  HD2 SING N N 15
PRO CD  HD3 SING N N 16
PRO OXT HXT SING N N 17
)"_cf;

	const std::filesystem::path example("../examples/1cbs.cif.gz");
	mmcif::File file(example.string());
	mmcif::Structure structure(file);

	mmcif::BondMap bm(structure);

	// Test the bonds of the first three residues, that's PRO A 1, ASN A 2, PHE A 3	

	for (const auto& [compound, seqnr]: std::initializer_list<std::tuple<std::string,int>>{ { "PRO", 1 }, { "ASN", 2 }, { "PHE", 3 } })
	{
		auto& res = structure.getResidue("A", compound, seqnr);
		auto atoms = res.atoms();

		auto dc = components.get(compound);
		BOOST_ASSERT(dc != nullptr);

		auto cc = dc->get("chem_comp_bond");
		BOOST_ASSERT(cc != nullptr);

		std::set<std::tuple<std::string,std::string>> bonded;

		for (const auto& [atom_id_1, atom_id_2]: cc->rows<std::string,std::string>({ "atom_id_1", "atom_id_2" }))
		{
			if (atom_id_1 > atom_id_2)
				bonded.insert({ atom_id_2, atom_id_1 });
			else
				bonded.insert({ atom_id_1, atom_id_2 });
		}

		for (size_t i = 0; i + 1 < atoms.size(); ++i)
		{
			auto label_i = atoms[i].labelAtomID();

			for (size_t j = i + 1; j < atoms.size(); ++j)
			{
				auto label_j = atoms[j].labelAtomID();

				bool bonded_1 = bm(atoms[i], atoms[j]);
				bool bonded_1_i = bm(atoms[j], atoms[i]);

				bool bonded_t = label_i > label_j
					? bonded.count({ label_j, label_i })
					: bonded.count({ label_i, label_j });

				BOOST_CHECK(bonded_1 == bonded_t);
				BOOST_CHECK(bonded_1_i == bonded_t);
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(bondmap_2)
{
	BOOST_CHECK_THROW(mmcif::BondMap::atomIDsForCompound("UN_"), mmcif::BondMapException);

	mmcif::CompoundFactory::instance().pushDictionary("./UN_.cif");

	BOOST_CHECK(mmcif::BondMap::atomIDsForCompound("UN_").empty() == false);
}