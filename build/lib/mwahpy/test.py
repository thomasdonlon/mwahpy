'''
The contents of this file are focused on testing the verious classes and routines
that are written out in the other files. All unit tests should be implemented
in this file so only one file needs to be run in order to complete all tests.
'''

import unittest
import numpy as np
import output_handler
from timestep import Timestep

#===============================================================================
# UNIT TESTING
#===============================================================================
#Most of the tests in this file read in the test.out MW@h output file that is
#packaged along with the developer version of mwahpy. Without this file, these
#tests will fail. If a different file is used, these tests will also fail.

#TODO: Complete tests for all functions in class

#NOTE: test.out does not have correct COM information, since it is a
#truncated version of a much larger MW@h output file.

#to run tests, run this file like you would run a regular python script

#to add more tests,
# 1) Find a test class that fits the description of the test you want to add,
#    then add that test to the class as a function beginning with "test", OR
# 2) Create a new class beginning with "test" and composed of just (unittest.TestCase)
#    and then add new functions that begin with "test" to that class.
#    This should be used if new tests do not fit the descriptions of other
#    test classes

#after changing this file, it should be run in order to make sure all tests pass
#if tests don't pass, congratulations you broke something
#please fix it

prec = 8 #number of digits to round to when comparing floats
#WARNING: Tests may (but are not expected to) fail at high levels of prec

class TestTimestepClass(unittest.TestCase):

    import output_handler

    def testTimestepInitialize(self):
        t = Timestep() #check default initialization
        t = output_handler.readOutput('../test/test.out') #check loading from a file

    def testTimestepDict(self):
        t = output_handler.readOutput('../test/test.out')
        self.assertTrue(t.x[0] == t['x'][0])

        t.x[0] = 1
        self.assertTrue(t['x'][0] == 1)

    def testTimestepIter(self):
        t = output_handler.readOutput('../test/test.out')

        for key, k in zip(t, t.indexList):
            self.assertTrue(t[key][0] == t[k][0])
            t[key][0] = 1

        self.assertTrue(t.x[0] == t['x'][0] == 1)

        for key in t:
            t[key] = np.append(t[key], 0)

        self.assertTrue(t.x[-1] == t['x'][-1] == 0)

    def testCopy(self):
        t = output_handler.readOutput('../test/test.out')

        t2 = t.copy()
        test = t.x[0]
        t.x[0] = 0

        self.assertTrue(t.x[0] != t2.x[0])
        self.assertTrue(t2.x[0] == test)

    def testAppendPoint(self):
        t = output_handler.readOutput('../test/test.out')

        t2 = t.copy()

        t2.appendPoint(t, 5)

        self.assertTrue(t.x[5] == t2.x[-1])
        self.assertTrue(len(t2) == len(t)+1)

    def testSplit(self):
        t = output_handler.readOutput('../test/test.out')

        t1, t2 = t.split(5)
        self.assertTrue(t1.x[0] == t.x[0])
        self.assertTrue(t2.x[0] == t.x[5])

    def testSubsetRect(self):
        t = output_handler.readOutput('../test/test.out')
        tc = t.copy()

        tc.subsetRect(('x',), ((-1,1),))

        self.assertTrue(len(tc) == 7)

    def testCalcs(self): #this just makes sure that the values that
        #are not initially calculated on instantiation then run when the user
        #tries to get those values

        t = output_handler.readOutput('../test/test.out')

        print(len(t.pmra))

        self.assertTrue(len(t.pmra) == len(t))
        self.assertTrue(len(t.energy) == len(t))

class TestTimestepMethods(unittest.TestCase):

    def testUpdate(self):
        t = output_handler.readOutput('../test/test.out')
        t.update(force=True)

        #these values for the COM's are hard coded, unfortunately. But they are correct.
        self.assertTrue(len(t.centerOfMass) == 3)
        print(t.centerOfMass)
        self.assertTrue(round(abs(t.centerOfMass[0] - 0.7835947518080879) + abs(t.centerOfMass[1] - 0.2947546230649471) + abs(t.centerOfMass[2] + 0.1318053758650839), prec) == 0)

        self.assertTrue(len(t.centerOfMomentum) == 3)
        self.assertTrue(round(abs(t.centerOfMomentum[0] - 6.001115213580641) + abs(t.centerOfMomentum[1] + 65.29652414026405) + abs(t.centerOfMomentum[2] + 26.462554427407223), prec) == 0)

        self.assertTrue(len(t.distFromCOM) == len(t))
        self.assertTrue(round(abs(t.distFromCOM[0] - ((t.centerOfMass[0] - t.x[0])**2 + (t.centerOfMass[1] - t.y[0])**2 + (t.centerOfMass[2] - t.z[0])**2)**0.5), prec) == 0)

#===============================================================================
# RUNTIME
#===============================================================================

if __name__ == '__main__':
    unittest.main()
