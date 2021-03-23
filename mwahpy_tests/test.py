'''
The contents of this file are focused on testing the verious classes and routines
that are written out in the other files. All unit tests should be implemented
in this file so only one file needs to be run in order to complete all tests.
'''

import unittest
import numpy as np
import mwahpy.output_handler as oh
from mwahpy.timestep import Timestep
import mwahpy.coords as co

#===============================================================================
# UNIT TESTING
#===============================================================================
#Most of the tests in this file read in the test_bodies.out MW@h output file that
#is packaged in the tests folder along with this file. Without the test_bodies file,
#these tests will fail. If a different file is used, these tests will also fail.

#TODO: Complete tests for all functions in class

#NOTE: test_bodies.out does not have correct COM information, since it is a
#truncated version of a much larger MW@h output file.

#to run tests, run this file like you would run a regular python script

#to add more tests,
# 1) Find a test class that fits the description of the test you want to add,
#    then add that test to the class as a function beginning with "test", OR
# 2) Create a new class beginning with "test" and composed of just (unittest.TestCase)
#    and then add new functions that begin with "test" to that class.
#    This should be used if new tests do not fit the descriptions of other
#    test classes

#after changing any mwahpy file, this file should be run in order to make sure
#that functionality is still as is intended
#if tests don't pass, congratulations you broke something
#please fix it

prec = 8 #number of digits to round to when comparing floats
#WARNING: Tests may (but are not expected to) fail at high levels of prec

class TestTimestepClass(unittest.TestCase):

    def test_timestep_initialize(self):
        t = Timestep() #check default initialization
        t = oh.read_output('test_bodies.out') #check loading from a file

    def test_timestep_dict(self):
        t = oh.read_output('test_bodies.out')
        self.assertTrue(t.x[0] == t['x'][0])

        t.x[0] = 1
        self.assertTrue(t['x'][0] == 1)

    def test_timestep_iter(self):
        t = oh.read_output('test_bodies.out')

        for key, k in zip(t, t.index_list):
            self.assertTrue(t[key][0] == t[k][0])
            t[key][0] = 1

        self.assertTrue(t.x[0] == t['x'][0] == 1)

        for key in t:
            t[key] = np.append(t[key], 0)

        self.assertTrue(t.x[-1] == t['x'][-1] == 0)

    def test_copy(self):
        t = oh.read_output('test_bodies.out')

        t2 = t.copy()
        test = t.x[0]
        t.x[0] = 0

        self.assertTrue(t.x[0] != t2.x[0])
        self.assertTrue(t2.x[0] == test)

    def test_append_point(self):
        t = oh.read_output('test_bodies.out')

        t2 = t.copy()

        t2.append_point(t, 5)

        self.assertTrue(t.x[5] == t2.x[-1])
        self.assertTrue(len(t2) == len(t)+1)

    def test_split(self):
        t = oh.read_output('test_bodies.out')

        t1, t2 = t.split(5)
        self.assertTrue(t1.x[0] == t.x[0])
        self.assertTrue(t2.x[0] == t.x[5])

    def test_subset_rect(self):
        t = oh.read_output('test_bodies.out')
        tc = t.copy()

        tc.subset_rect(('x',), ((-1,1),))

        self.assertTrue(len(tc) == 7)

    def test_calcs(self): #makes sure that the values that
        #are not initially calculated on instantiation then run when the user
        #tries to get those values

        t = oh.read_output('test_bodies.out')

        self.assertTrue(len(t.pmra) == len(t))
        self.assertTrue(len(t.energy) == len(t))

class TestTimestepMethods(unittest.TestCase):

    def test_update(self):
        t = oh.read_output('test_bodies.out')
        t.update(force=True)

        #these values for the COM's are hard coded, unfortunately. But they are correct.
        self.assertTrue(len(t.center_of_mass) == 3)
        self.assertTrue(round(abs(t.center_of_mass[0] - 0.7835947518080879) + abs(t.center_of_mass[1] - 0.2947546230649471) + abs(t.center_of_mass[2] + 0.1318053758650839), prec) == 0)

        self.assertTrue(len(t.center_of_momentum) == 3)
        self.assertTrue(round(abs(t.center_of_momentum[0] - 6.001115213580641) + abs(t.center_of_momentum[1] + 65.29652414026405) + abs(t.center_of_momentum[2] + 26.462554427407223), prec) == 0)

        self.assertTrue(len(t.distFromCOM) == len(t))
        self.assertTrue(round(abs(t.distFromCOM[0] - ((t.center_of_mass[0] - t.x[0])**2 + (t.center_of_mass[1] - t.y[0])**2 + (t.center_of_mass[2] - t.z[0])**2)**0.5), prec) == 0)

class TestCoordsInverses(unittest.TestCase):

    #First order coordinate (position) transformations

    def test_cart_to_gal_rh(self):
        #right-handed
        x, y, z = 1., 2., 3.
        l, b, r = co.cart_to_gal(x, y, z)
        new_x, new_y, new_z = co.gal_to_cart(l, b, r)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_gal_lh(self):
        #left-handed
        x, y, z = 1., 2., 3.
        l, b, r = co.cart_to_gal(x, y, z, left_handed=True)
        new_x, new_y, new_z = co.gal_to_cart(l, b, r, left_handed=True)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_cyl(self):
        x, y, z = 1., 2., 3.
        R, z, phi = co.cart_to_cyl(x, y, z)
        new_x, new_y, new_z = co.cyl_to_cart(R, z, phi)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_sph(self):
        x, y, z = 1., 2., 3.
        phi, theta, r = co.cart_to_sph(x, y, z)
        new_x, new_y, new_z = co.sph_to_cart(phi, theta, r)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    #Second order coordinate (position) transformations

    def test_cyl_to_gal(self):
        R, z, phi = 10., 5., 120.
        l, b, r = co.cyl_to_gal(R, z, phi)
        new_R, new_z, new_phi = co.gal_to_cyl(l, b, r)
        self.assertTrue((round(new_R, prec) == round(R, prec)) and (round(new_z, prec) == round(z, prec)) and (round(new_phi, prec) == round(phi, prec)))

    #higher order coordinate transformations

    def test_sky_to_pole(self):
        tol = 1e-6 #tolerance for this test
        ra = np.linspace(10, 350, 10) #avoid poles, which do not test well but work
        #the ra rotates strangely but the dec = +/-90 so it doesn't matter
        dec = np.linspace(-80, 80, 10)
        ra_new, dec_new = co.sky_to_pole(ra, dec, (0, 90), (0, 0)) #null transformation
        self.assertTrue((np.sum(np.abs(ra_new - ra)) <= tol) and (np.sum(np.abs(dec_new - dec)) <= tol))

        ra = np.linspace(10, 350, 10) #avoid poles, which do not test well but work
        dec = np.linspace(-80, 80, 10)
        L, B = co.sky_to_pole(ra, dec, (0, 0), (90, 0))
        ra_new, dec_new = co.sky_to_pole(L, B, (90, 0), (0, 90)) #this *should* be the inverse transformation
        self.assertTrue((np.sum(np.abs(ra_new - ra)) <= tol) and (np.sum(np.abs(dec_new - dec)) <= tol))

#===============================================================================
# RUNTIME
#===============================================================================

if __name__ == '__main__':
    unittest.main()
