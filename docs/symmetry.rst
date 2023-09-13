Symmetry & Geometry
===================

Although not really a core *CIF* functionality, when working with *mmCIF* files you often need to work with symmetry information. And symmetry works on points in a certain space and thus geometry calculations are also something you need often. Former versions of *libcifpp* used to use `clipper <http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/index.html>`_ to do many of these calculations, but that introduces a dependency and besides, the way clipper numbers symmetry operations is not completely compatible with the way this is done in the PDB.

Points
------

The most basic type in use is :cpp:type:`cif::point`. It can be thought of as a point in space with three coordinates, but it is also often used as a vector in 3d space. To keep the interface simple there's no separate vector type.

Many functions are available in :ref:`file_cif++_point.hpp` that work on points. There are functions to calculate the :cpp:func:`cif::distance` between two points and also function to calculate dot products, cross products and dihedral angles between sets of points.

Quaternions
-----------

All operations inside *libcifpp* that perform some kind of rotation use :cpp:type:`cif::quaternion`. The reason to use Quaternions is not only that they are cool, they are faster than multiplying with a matrix and the results also suffer less from numerical instability.

Matrix
------

Although Quaternions are the preferred way of doing rotations, not every manipulation is a rotation and thus we need a matrix class as well. Matrices and their operations are encoded as matrix_expressions in *libcifpp* allowing the compiler to generate very fast code. See the :ref:`file_cif++_matrix.hpp` for what is on offer.

Crystals
--------

The *CIF* and *mmCIF* were initially developed to store crystallographic information on structures. Apart from coordinates and the chemical information the crystallographic information is important. This information can be split into two parts, a unit cell and a set of  :ref:`symmetry-ops` making up a spacegroup. The spacegroup number and name are stored in the *symmetry* category. The corresponding symmetry operations can be obtained in *libcifpp* by using the :cpp:class:`cif::spacegroup`. The cell is stored in the category *cell* and likewise can be loaded using the :cpp:class:`cif::cell`. Together these two classes make up a crystal and so we have a :cpp:class:`cif::crystal` which contains both. You can easily create such a crystal object by passing the datablock containing the data to the constructor. As in:

.. code:: cpp

    // Load the file
    cif::file f("1cbs.cif.gz");

    auto &db = f.front();
    cif::crystal c(db);

.. _symmetry-ops:
Symmetry operations
-------------------

Each basic symmetry operation in the crystallographic world consists of a matrix multiplication followed by a translation. To apply such an operation on a carthesian coordinate you first have to convert the point into a fractional coordinate with respect to the unit cell of the crystal, then apply the matrix and translation operations and then convert the result back into carthesian coordinates. This is all done by the proper routines in *libcifpp*.

Symmetry operations are encoded as a string in *mmCIF* PDBx files. The format is a string with the rotational number followed by an underscore and then the encoded translation in each direction where 5 means no translation. So, the identity operator is ``1_555`` meaning that we have rotational number 1 (which is always the identity rotation, point multiplied with the identity matrix) and a translation of zero in each direction.

To give an idea how this works, here's a piece of code copied from one of the unit tests in *libcifpp*. It takes the *struct_conn* records in a certain PDB file and checks wether the distances in each row correspond to what we can calculate.

.. code:: cpp

    using namespace cif::literals;

    // Load the file
    cif::file f(gTestDir / "2bi3.cif.gz");

    // Simply assume we can use the first datablock
    auto &db = f.front();

    // Load the crystal information
    cif::crystal c(db);

    // Take references to the two categories we need
    auto struct_conn = db["struct_conn"];
    auto atom_site = db["atom_site"];

    // Loop over all rows in struct_conn taking only the values we need
    for (const auto &[
            asym1, seqid1, authseqid1, atomid1, symm1,
            asym2, seqid2, authseqid2, atomid2, symm2,
            dist] : struct_conn.find<
                std::string,std::optional<int>,std::string,std::string,std::string,
                std::string,std::optional<int>,std::string,std::string,std::string,
                float>(
            cif::key("ptnr1_symmetry") != "1_555" or cif::key("ptnr2_symmetry") != "1_555",
            "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry", 
            "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry", 
            "pdbx_dist_value"
        ))
    {
        // Find the location of the first atom
        cif::point p1 = atom_site.find1<float,float,float>(
            "label_asym_id"_key == asym1 and "label_seq_id"_key == seqid1 and "auth_seq_id"_key == authseqid1 and "label_atom_id"_key == atomid1,
            "cartn_x", "cartn_y", "cartn_z");

        // Find the location of the second atom
        cif::point p2 = atom_site.find1<float,float,float>(
            "label_asym_id"_key == asym2 and "label_seq_id"_key == seqid2 and "auth_seq_id"_key == authseqid2 and "label_atom_id"_key == atomid2,
            "cartn_x", "cartn_y", "cartn_z");

        // Calculate the position of the first atom using the symmetry operator defined in struct_conn
        auto sa1 = c.symmetry_copy(p1, cif::sym_op(symm1));

        // Calculate the position of the second atom using the symmetry operator defined in struct_conn
        auto sa2 = c.symmetry_copy(p2, cif::sym_op(symm2));

        // The distance between these symmetry atoms should be equal to the distance in the struct_conn record
        assert(cif::distance(sa1, sa2) == dist);

        // And to show how you can obtain the closest symmetry copy of an atom near another one:
        // here we request the symmetry copy of p2 that lies closest to p1
        const auto &[d, p, so] = c.closest_symmetry_copy(p1, p2);

        // And that should of course be equal to the location in struct_conn for p2
        assert(p.m_x == sa2.m_x);
        assert(p.m_y == sa2.m_y);
        assert(p.m_z == sa2.m_z);

        // Distance and symmetry operator string should also be the same
        assert(d == dist);
        assert(so.string() == symm2);
    }