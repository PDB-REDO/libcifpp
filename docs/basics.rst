Basic usage
===========

This library, *libcifpp*, is a generic *CIF* library with some specific additions to work with *mmCIF* files. The main focus of this library is to make sure that files read or written are valid. That is, they are syntactically valid *and* their content is valid with respect to a CIF dictionary, if such a dictionary is available and specified.

Reading a file is as simple as:

.. code-block:: cpp

    #include <cif++.hpp>

    cif::file f("/path/to/file.cif");

The file may also be compressed using *gzip* which is detected automatically.

Writing out the file again is also simple, to write out the terminal you can do:

.. code-block:: cpp

    std::cout << f;

    // or
    f.save(std::cout);

    // or write a compressed file using gzip compression:
    f.save("/tmp/f.cif.gz");

CIF files contain one or more datablocks. To print out the names of all datablocks in our file:

.. code-block:: cpp

    for (auto &db : f)
        std::cout << db.name() << '\n';

Most often *libcifpp* is used to read in structure files in mmCIF format. These files only contain one datablock and so you can safely use code like this:

.. code-block:: cpp

    // get a reference to the first datablock in f
    auto &db = f.front();

But if you know the name of the datablock, this also works:

.. code-block:: cpp

    // get a reference to the datablock name '1CBS'
    auto &db = f["1CBS"];

Now, each datablock contains categories. To print out all their names:

.. code-block:: cpp

    for (auto &cat : db)
        std::cout << cat.name() << '\n';

But you probably know what category you need to use, so lets fetch it by name:

.. _atom_site-label:
.. code-block:: cpp

    // get a reference to the atom_site category in db
    auto &atom_site = db["atom_site"];

    // and make sure there's some data in it:
    assert(not atom_site.empty());

.. note::
    
    Note that we omit the leading underscore in the name of the category here.

Categories contain rows of data and each row has fields or items. Referencing a row in a category results in a :cpp:class:`cif::row_handle` object which you can use to request or manipulate item data.

.. code-block:: cpp

    // Get the first row in atom_site
    auto rh = atom_site.front();

    // Get the label_atom_id value from this row handle as a std::string
    std::string atom_id = rh["label_atom_id"].as<std::string>();

    // Get the x, y and z coordinates using structered binding
    const auto &[x, y, z] = rh.get<float,float,float>("Cartn_x", "Cartn_y", "Cartn_z");

    // Assign a new value to the x coordinate or our atom
    rh["Cartn_x"] = x + 1;

Querying
--------

Walking over the rows in a category is often not very useful. More often you are interested in specific rows in a category. The function :cpp:func:`cif::category::find` and friends are here to help.

What these functions have in common is that they return data based on a query implemented by :cpp:class:`cif::condition`. These condition objects are built in code using regular C++ syntax. The most basic example of a query is:

.. code-block:: cpp

    cif::condition c = cif::key("id") == 1;

Here the condition is that all rows returned should have a value of 1 in there item named *id*. Likewise you can use other data types and even combine those. Oh, and I said we use regular C++ syntax for conditions, so you may as well use other operators to compare values:

.. code-block:: cpp

    // condition for C-alpha atoms having an occupancy less than 1.0
    cif::condition c = cif::key("occupancy") < 1.0f and cif::key("label_atom_id") == "CA";

Using the namespace *cif::literals* that code becomes a little less verbose:

.. code-block:: cpp

    using namespace cif::literals;
    cif::condition c = "occupancy"_key < 1.0f and "label_atom_id"_key == "CA";

Conditions can also be combined:

.. code-block:: cpp

    cif::condition c = "occupancy"_key < 1.0f and "label_atom_id"_key == "CA";

    // extend the condition by requiring the compound ID to be unequal to PRO
    c = std::move(c) and "label_comp_id"_key != "PRO";

.. note::

    Note the use of std::move here. 

Using queries constructed in this way is simple:

.. code-block:: cpp

    cif::condition c = ...
    auto result = atom_site.find(std::move(c));

    // or construct a condition inline:
    auto result = atom_site.find("label_atom_id"_key == "CA");

In the example above the result is a range of :cpp:class:`cif::row_handle` objects. Often, using individual field values is more useful:

.. code-block:: cpp

    // Requesting a single item:
    for (auto id : atom_site.find<std::string>("label_atom_id"_key == "CA", "id"))
        std::cout << "ID for CA: " << id << '\n';

    // Requesting multiple items:
    for (const auto &[id, x, y, z] : atom_site.find<std::string,float,float,float>("label_atom_id"_key == "CA",
            "id", "Cartn_x", "Cartn_y", "Cartn_z"))
    {
        std::cout << "Atom " << id << " is at [" << x << ", " << y << ", " z << "]\n";
    }

Returning a complete set if often not required, if you only want to have the first you can use :cpp:func:`cif::category::find_first` as shown here:

.. code-block:: cpp

    // return the ID item for the first C-alpha atom
    std::string v1 = atom_site.find_first<std::string>("label_atom_id"_key == "CA", "id");

    // If you're not sure the row exists, use std::optional
    auto v2 = atom_site.find_first<std::optional<std::string>>("label_atom_id"_key == "CA", "id");
    if (v2.has_value())
        ...

There are cases when you really need exactly one result. The :cpp:func:`cif::category::find1` can be used in that case, it will throw an exception if the query does not result in exactly one row.

NULL and ANY
------------

Sometimes items may be empty. The trouble is a bit that empty comes in two flavors: unknown and null. Null in *CIF* parlance means the item should not contain a value since it makes no sense in this case, the value stored in the file is a single dot character: ``'.'``. E.g. *atom_site* records may have a NULL value for label_seq_id for atoms that are part of a *non-polymer*.

The other empty value is indicated by a question mark character: ``'?'``. This means the value is simply unknown.

Both these are NULL in *libcifpp* conditions and can be searched for using :cpp:var:`cif::null`.

So you can search for:

.. code-block:: cpp

    cif::condition c = "label_seq_id"_key == cif::null;

You might also want to look for a certain value and don't care in which item it is stored, in that case you can use :cpp:var:`cif::any`.

.. code-block:: cpp

    cif::condition c = cif::any == "foo";

And in linked record you might have the items that have a value in both parent and child or both should be NULL. For that, you can request the value to return by find to be of type std::optional and then use that value to build the query. An example to explain this, let's find the location of the atom that is referenced as the first atom in a struct_conn record:

.. code-block:: cpp

    // Take references to the two categories we need
    auto struct_conn = db["struct_conn"];
    auto atom_site = db["atom_site"];

    // Loop over all rows in struct_conn taking only the values we need
    // Note that the label_seq_id is returned as a std::optional<int>
    // That means it may contain an integer or may be empty
    for (const auto &[asym1, seqid1, authseqid1, atomid1] :
        struct_conn.rows<std::string,std::optional<int>,std::string,std::string,std::string>(
            "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id"
        ))
    {
        // Find the location of the first atom
        cif::point p1 = atom_site.find1<float,float,float>(
            "label_asym_id"_key == asym1 and "label_seq_id"_key == seqid1 and "auth_seq_id"_key == authseqid1 and "label_atom_id"_key == atomid1,
            "cartn_x", "cartn_y", "cartn_z");
    }
    

Validation
----------

CIF files can have a dictionary attached. And based on such a dictionary a :cpp:class:`cif::validator` object can be constructed which in turn can be used to validate the content of the file.

A simple case:

.. code-block:: cpp

    #include <cif++.hpp>

    cif::file f("1cbs.cif.gz");
    f.load_dictionary("mmcif_pdbx");

    if (not f.is_valid())
        std::cout << "This file is not valid\n";

If you want to know why it is not valid, you should set the global variable :cpp:var:`cif::VERBOSE` to something higer than zero. Depending on the value more or less diagnostic output is sent to std::cerr.

In the case above we load a dictionary based on its name. You can of course also load dictionaries based on a specific file, that's a bit more work:

.. code-block:: cpp

    std::filesystem::ifstream dictFile("/tmp/my-dictionary.dic");
    auto &validator = cif::parse_dictionary("my-dictionary", dictFile);

    cif::file f("1cbs.cif.gz");

    // assign the validator
    f.set_validator(&validator);

    // alternatively, load it by name
    f.load_dictionary("my-dictionary");

    if (not f.is_valid())
        std::cout << "This file is not valid\n";

Creating your own dictionary is a lot of work, especially if you are only extending an existing dictionary with a couple of new categories or items. So, what you can do is extend a loaded validator like this (code taken from DSSP):

.. code-block:: cpp

    // db is a cif::datablock reference containing an mmCIF file with DSSP annotations
    auto &validator = const_cast<cif::validator &>(*db.get_validator());
    if (validator.get_validator_for_category("dssp_struct_summary") == nullptr)
    {
        auto dssp_extension = cif::load_resource("dssp-extension.dic");
        if (dssp_extension)
            cif::extend_dictionary(validator, *dssp_extension);
    }

.. note::

    In the example above we're loading the data using :doc:`/resources`. See the documentation on that for more information.

If a validator has been assigned to a file, assignments to items are checked for valid data. So the following code will throw an exception (see: :ref:`_atom_site-label`):

.. code-block:: cpp
    
    auto rh = atom_site.front();
    rh["Cartn_x"] = "foo";

Linking
-------

Based on information recorded in dictionary files (see :ref:`Validation`) you can locate linked records in parent or child categories.

To make this example not too complex, lets assume the following example file:

.. code-block:: cif

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

And we have a dictionary containing the following link definition:

.. code-block:: cif

    loop_
    _pdbx_item_linked_group_list.parent_category_id
    _pdbx_item_linked_group_list.link_group_id
    _pdbx_item_linked_group_list.parent_name
    _pdbx_item_linked_group_list.child_name
    _pdbx_item_linked_group_list.child_category_id
    cat_1 1 '_cat_1.name' '_cat_2.name' cat_2

So, there are links between *cat_1* and *cat_2* based on the value in items named *name*. Using this information, we can now locate children and parents:

.. code-block:: cpp

    // Assuming the file was loaded in f:
    auto &cat1 = f.front()["cat_1"];
    auto &cat2 = f.front()["cat_2"];
    auto &cat3 = f.front()["cat_3"];

    // Loop over all ape's in cat2
    for (auto r : cat1.get_children(cat1.find1("name"_key == "aap"), cat2))
        std::cout << r.get<std::string>("desc") << '\n';

Updating a value in an item in a parent category will update the corresponding value in all related children:

.. code-block:: cpp

    auto r1 = cat1.find1("id"_key == 1);
    r1["name"] = "aapje";

    auto rs1 = cat2.find("name"_key == "aapje");
    assert(rs1.size() == 2);

However, changing a value in a child record will not update the parent. This may result in an invalid file since you may then have a child that has no parent:

.. code-block:: cpp

    auto r2 = cat2.find1("id"_key == 3);
    r2["name"] = "wim";

    assert(f.is_valid() == false);

So you have to fix this yourself by inserting a new item in cat1 with the new value.

.. _splitting-rows:
Another situation is when you change a value in a parent and updating children might introduce a situation where you need to split a child. To give an example, consider this:

.. code-block:: cif

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

And we have a dictionary containing the following link definition (reversed compared to the previous example):

.. code-block:: cif

    loop_
    _pdbx_item_linked_group_list.parent_category_id
    _pdbx_item_linked_group_list.link_group_id
    _pdbx_item_linked_group_list.parent_name
    _pdbx_item_linked_group_list.child_name
    _pdbx_item_linked_group_list.child_category_id
    cat_2 1 '_cat_2.name' '_cat_1.name' cat_1
    cat_3 1 '_cat_3.name' '_cat_2.name' cat_2
    cat_3 1 '_cat_3.num'  '_cat_2.num'  cat_2

So *cat3* is a parent of *cat2* and *cat2* is a parent of *cat1*. Now, if you change the *name* value of the first row of *cat3* to 'aapje', the corresponding row in *cat2* is updated as well. But when you update *cat2* you have to update *cat1* too. And simply changing the name field in row 1 of *cat1* is wrong. The default behaviour in libcifpp is to split the record in *cat1* and have a new child with the new name whereas the other remains as is.

The new *cat1* will thus be like:

.. code-block:: cif

    loop_
    _cat_1.id
    _cat_1.name
    _cat_1.desc
    1 aapje Aap
    2 noot  Noot
    3 mies  Mies
    5 aap   Aap

