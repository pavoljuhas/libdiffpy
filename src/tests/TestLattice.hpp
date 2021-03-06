/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class TestLattice -- unit tests for the Lattice class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/Lattice.hpp>
#include <diffpy/srreal/R3linalg.hpp>

using namespace std;
using namespace diffpy::srreal;
using diffpy::mathutils::EpsilonEqual;

class TestLattice : public CxxTest::TestSuite
{

private:

    Lattice* lattice;
    double precision;
    EpsilonEqual allclose;

public:

    void setUp()
    {
        lattice = new Lattice();
        precision = 1.0e-8;
        allclose = EpsilonEqual(precision);
    }

    void tearDown()
    {
        delete lattice;
    }

    void test_Lattice()
    {
        // default lattice should be cartesian
        TS_ASSERT(allclose(lattice->base(), R3::identity()));
        // lattice parameters constructor
        Lattice lattice1(1.0, 2.0, 3.0, 90, 90, 120);
        R3::Vector va = lattice1.va();
        R3::Vector vb = lattice1.vb();
        R3::Vector vc = lattice1.vc();
        double adotb = R3::dot(va, vb);
        double adotc = R3::dot(va, vc);
        double bdotc = R3::dot(vb, vc);
        TS_ASSERT_DELTA(1.0, R3::norm(va), precision);
        TS_ASSERT_DELTA(2.0, R3::norm(vb), precision);
        TS_ASSERT_DELTA(3.0, R3::norm(vc), precision);
        TS_ASSERT_DELTA(-0.5*1.0*2.0, adotb, precision);
        TS_ASSERT_DELTA(0.0, adotc, precision);
        TS_ASSERT_DELTA(0.0, bdotc, precision);
        TS_ASSERT_DELTA(sqrt(4.0/3), lattice1.ar(), precision);
        TS_ASSERT_DELTA(sqrt(1.0/3), lattice1.br(), precision);
        TS_ASSERT_DELTA(1.0/3, lattice1.cr(), precision);
        TS_ASSERT_DELTA(90.0, lattice1.alphar(), precision);
        TS_ASSERT_DELTA(90.0, lattice1.betar(), precision);
        TS_ASSERT_DELTA(60.0, lattice1.gammar(), precision);
        // lattice vectors constructor
        va = R3::Vector(1.0, 1.0, 0.0);
        vb = R3::Vector(0.0, 1.0, 1.0);
        vc = R3::Vector(1.0, 0.0, 1.0);
        Lattice lattice2(va, vb, vc);
        TS_ASSERT_DELTA(sqrt(2.0), lattice2.a(), precision);
        TS_ASSERT_DELTA(60.0, lattice2.alpha(), precision);
    }

    void test_setLatPar()
    {
        lattice->setLatPar(1.0, 2.0, 3.0, 90, 90, 120);
        R3::Matrix base_check(
            sqrt(0.75),    -0.5,      0.0,
            0.0,           +2.0,      0.0,
            0.0,           +0.0,      3.0);
        TS_ASSERT(allclose(base_check, lattice->base()));
        R3::Matrix recbase_check(
            sqrt(4.0/3),    sqrt(1.0/12),   0.0,
            0.0,            0.5,            0.0,
            0.0,            0.0,            1.0/3.0);
        TS_ASSERT(allclose(recbase_check, lattice->recbase()));
    }

    void test_setLatBase()
    {
        R3::Vector va(1.0,  1.0,  0.0);
        R3::Vector vb(0.0,  1.0,  1.0);
        R3::Vector vc(1.0,  0.0,  1.0);
        lattice->setLatBase(va, vb, vc);
        TS_ASSERT_DELTA(sqrt(2.0), lattice->a(), precision);
        TS_ASSERT_DELTA(sqrt(2.0), lattice->b(), precision);
        TS_ASSERT_DELTA(sqrt(2.0), lattice->c(), precision);
        TS_ASSERT_DELTA(60.0, lattice->alpha(), precision);
        TS_ASSERT_DELTA(60.0, lattice->beta(), precision);
        TS_ASSERT_DELTA(60.0, lattice->gamma(), precision);
        // check determinant of rotation matrix
        double detRot0 = R3::determinant(lattice->baserot());
        TS_ASSERT_DELTA(1.0, detRot0, precision);
        // check if rotation matrix works
        R3::Matrix base_check(
                va[0], va[1], va[2],
                vb[0], vb[1], vb[2],
                vc[0], vc[1], vc[2]);
        TS_ASSERT(allclose(base_check, lattice->base()));
        R3::Matrix recbase_check(
                0.5,   -0.5,    0.5,
                0.5,    0.5,   -0.5,
                -0.5,   0.5,    0.5);
        TS_ASSERT(allclose(recbase_check, lattice->recbase()));
        lattice->setLatPar( lattice->a(), lattice->b(), lattice->c(),
                44.0, 66.0, 88.0 );
        TS_ASSERT(!allclose(base_check, lattice->base()));
        TS_ASSERT(!allclose(recbase_check, lattice->recbase()));
        lattice->setLatPar( lattice->a(), lattice->b(), lattice->c(),
                60.0, 60.0, 60.0 );
        TS_ASSERT(allclose(base_check, lattice->base()));
        TS_ASSERT(allclose(recbase_check, lattice->recbase()));
    }

    void test_dist()
    {
        R3::Vector va(1.0, 2.0, 2.0);
        R3::Vector vb(0.0, 0.0, 0.0);
        TS_ASSERT_DELTA(3.0,
                lattice->distance(va, vb), precision);
        lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 90.0);
        TS_ASSERT_DELTA(6.0,
                lattice->distance(va, vb), precision);
    }

    void test_angle()
    {
        R3::Vector va(1.0, 0.0, 0.0);
        R3::Vector vb(0.0, 1.0, 0.0);
        TS_ASSERT_DELTA(90.0,
                lattice->angledeg(va, vb), precision);
        lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 120.0);
        TS_ASSERT_DELTA(120.0,
                lattice->angledeg(va, vb), precision);
    }

    void test_ucvCartesian()
    {
        R3::Vector ucv;
        ucv = lattice->ucvCartesian(R3::Vector(0.1, 0.2, 0.3));
        R3::Vector ucv_check(0.1, 0.2, 0.3);
        TS_ASSERT(allclose(ucv_check, ucv));
        ucv = lattice->ucvCartesian(R3::Vector(1.1, 13.2, -0.7));
        TS_ASSERT(allclose(ucv_check, ucv));
        ucv = lattice->ucvCartesian(R3::Vector(0.5, 0.5, 0.5));
        ucv_check = R3::Vector(0.5, 0.5, 0.5);
        TS_ASSERT(allclose(ucv_check, ucv));
        lattice->setLatPar(13, 17, 19, 37, 41, 47);
        ucv = lattice->ucvCartesian(R3::Vector(100.0, 100.0, 100.0));
        ucv_check = R3::Vector(8.09338442077, 9.55056747225, 30.0389043325);
        TS_ASSERT(allclose(ucv_check, ucv));
    }

};  // class TestLattice

// End of file
