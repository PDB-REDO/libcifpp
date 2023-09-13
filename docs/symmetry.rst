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

Symmetry operations
-------------------

Each basic symmetry operation in the crystallographic world consists of a matrix multiplication followed by a translation. To apply such an operation on a carthesian coordinate you first have to convert the point into a fractional coordinate with respect to the unit cell of the crystal, then apply the matrix and translation operations and then convert the result back into carthesian coordinates. This is all done by the proper routines in *libcifpp*.

Symmetry operations are encoded as a string in *mmCIF* PDBx files. The format is a string with the rotational number followed by an underscore and then the encoded translation in each direction where 5 means no translation. So, the identity operator is ``1_555`` meaning that we have rotational number 1 (which is always the identity rotation, point multiplied with the identity matrix) and a translation of zero in each direction.

To give an idea how this works, here's a piece of code copied from one of the unit tests in *libcifpp*. It takes the *struct_conn* records in a certain PDB file and checks wether the distances in each row correspond to what we can calculate.

.. code:: cpp

    // Load the file
    cif::file f(gTestDir / "2bi3.cif.gz");

    auto &db = f.front();
    cif::mm::structure s(db);

    cif::crystal c(db);

    auto struct_conn = db["struct_conn"];
    for (const auto &[
            asym1, seqid1, authseqid1, atomid1, symm1,
            asym2, seqid2, authseqid2, atomid2, symm2,
            dist] : struct_conn.find<
                std::string,int,std::string,std::string,std::string,
                std::string,int,std::string,std::string,std::string,
                float>(
            cif::key("ptnr1_symmetry") != "1_555" or cif::key("ptnr2_symmetry") != "1_555",
            "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry", 
            "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry", 
            "pdbx_dist_value"
        ))
    {
        auto &r1 = s.get_residue(asym1, seqid1, authseqid1);
        auto &r2 = s.get_residue(asym2, seqid2, authseqid2);

        auto a1 = r1.get_atom_by_atom_id(atomid1);
        auto a2 = r2.get_atom_by_atom_id(atomid2);

        auto sa1 = c.symmetry_copy(a1.get_location(), cif::sym_op(symm1));
        auto sa2 = c.symmetry_copy(a2.get_location(), cif::sym_op(symm2));

        BOOST_TEST(cif::distance(sa1, sa2) == dist);

        auto pa1 = a1.get_location();

        const auto &[d, p, so] = c.closest_symmetry_copy(pa1, a2.get_location());

        BOOST_TEST(p.m_x == sa2.m_x);
        BOOST_TEST(p.m_y == sa2.m_y);
        BOOST_TEST(p.m_z == sa2.m_z);

        BOOST_TEST(d == dist);
        BOOST_TEST(so.string() == symm2);
    }