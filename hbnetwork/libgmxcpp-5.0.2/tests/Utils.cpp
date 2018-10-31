#include <assert.h>
#include "tests.h"
#include "gmxcpp/coordinates.h"
#include "gmxcpp/coordinates4.h"
#include "gmxcpp/triclinicbox.h"
#include "gmxcpp/Utils.h"
#include <iomanip>

int main()
{
    /* Test some utils*/
    srand(0);
    triclinicbox b1(3.5, 4.5, 4.0);

    coordinates c6(3.6, 4.7, 5.0);
    coordinates c7 = pbc(c6, b1);
    assert(test_equal(c7[X], 0.1));
    assert(test_equal(c7[Y], 0.2));
    assert(test_equal(c7[Z], 1.0));

    triclinicbox b2(5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 2.5, 2.5, 3.5);
    coordinates c8(5.5, 5.5, 3.5);

    coordinates c9 = pbc(c8, b2);
    assert(test_equal(c9[X], -2.0));
    assert(test_equal(c9[Y], -2.0));
    assert(test_equal(c9[Z], 0.0));

    assert(test_equal(dot(c6, c7), 6.3));
    assert(test_equal(magnitude(c6), 7.749193506));

    assert(test_equal(distance(c8, c9, b2), 0.0));
    assert(test_equal(distance(c6, c8, b1), 2.33452));

/*
    coordinates4 c_4(c6, c7, c8, c9);
    vector <float> d2_4 = distance2(c_4, c9, b2);
    assert(test_equal(distance2(c6, c9, b2), d2_4[0]));
    assert(test_equal(distance2(c7, c9, b2), d2_4[1]));
    assert(test_equal(distance2(c8, c9, b2), d2_4[2]));
    assert(test_equal(distance2(c9, c9, b2), d2_4[3]));
*/

    vector <coordinates> sites;
    triclinicbox b3(10.0, 10.0, 10.0);
    coordinates site1(1.0, 1.0, 1.0);
    sites.push_back(site1);
    assert(test_equal(get_surf_area(sites, 1.0, 10000, b3), 12.56637061));

    coordinates site2(1.0, 1.0, 4.0);
    sites.push_back(site2);
    assert(test_equal(get_surf_area(sites, 1.0, 10000, b3), 25.13274123));

	coordinates a(0.0,0.0,0.0);
	coordinates b(1.0,0.0,0.0);
	coordinates c(1.0,1.0,0.0);
	assert(test_equal(bond_angle(a,b,c,b1),M_PI/2.0));

	coordinates d(0.5,0.0,0.86602540378443864);
	assert(test_equal(bond_angle(a,b,d,b1),M_PI/3.0));

	coordinates e(1.0,1.0,-1.0);
	assert(test_equal(dihedral_angle(a,b,c,e,b1),M_PI/2.0));
	coordinates f(1.0,1.0,1.0);
	assert(test_equal(dihedral_angle(a,b,c,f,b1),-M_PI/2.0));

    coordinates com1(0.0,0.0,0.0);
    coordinates com2(0.0,0.0,0.0);
    coordinates com3(0.0,0.0,0.0);

    vector <coordinates> com_test;
    vector <double> mass;
    mass.push_back(1.0);
    mass.push_back(2.0);
    mass.push_back(3.0);

    com_test.push_back(com1);
    com_test.push_back(com2);
    com_test.push_back(com3);

    cubicbox combox(4.0,4.0,4.0);

    coordinates com = center_of_mass(com_test,mass,combox);
    assert(test_equal(com[X],0.0));
    assert(test_equal(com[Y],0.0));
    assert(test_equal(com[Z],0.0));

    coordinates com4(1.0,1.0,1.0);
    coordinates com5(3.0,3.0,3.0);
    coordinates com6(0.0,0.0,0.0);
    com_test[0] = com4;
    com_test[1] = com5;
    com_test[2] = com6;

    com = center_of_mass(com_test,mass,combox);
    assert(test_equal(com[X],-0.16667));
    com = center_of_mass(com_test,mass);
    assert(test_equal(com[X],1.16667));

    return 0;
}
