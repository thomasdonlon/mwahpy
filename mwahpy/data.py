'''
The contents of this file are focused on the Data class, which is used for storage of
imported data from N-body output files.
'''

#TODO: split functions should fix id values? Maybe not, depends on what the behavior is supposed to do
#   this could probably be added into update
#TODO: UNIT TESTS

#===============================================================================
# IMPORTS
#===============================================================================

#external imports
import numpy as np
import coords as co
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import random
import galpy
import galpy.potential
import unittest
import os

#mwahpy imports
import mwahpy_glob
import flags
import output_handler
import pot

#===============================================================================
# DATA CLASS
#===============================================================================

#AttrDict is used as a helper class in Data to allow referencing attributes
#as dict keys and vice-versa.
#this is probably a bad way to implement this but it works, and it's better than
#making Data inherit from dict, which was the other solution I was able to strum up
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class Data():

    def __init__(self, id_val=[], x=[], y=[], z=[], l=[], b=[], r=[], vx=[], vy=[], vz=[], mass=[], vlos=[], centerOfMass=[0, 0, 0], centerOfMomentum=[0, 0, 0], pot_offset=0, *args, **kwargs):
        #all this is typically specified by the readOutput function in output_handler
        #readOutput is the preferred way to input data to this data structure
        #but if you're really feeling adventurous you can always do it yourself

        #-----------------------------------------------------------------------
        # HOUSEKEEPING
        #-----------------------------------------------------------------------

        #this data structure allows access to a dictionary like an attribute, i.e.
        #   d = AttrDict({'a':1, 'b':2})
        #   d['a'] == d.a //returns True
        #it seems like magic but it works by accessing the way python handles attributes
        #within each class. There are pros and cons, the most notable thing is
        #that it creates a memory leak in python <3.2.3
        ad = AttrDict() #this is a private helper class that allows for the
                        #desired indexing behavior
        self.__dict__ = ad.__dict__

        #-----------------------------------------------------------------------

        self.id = np.array(id_val)
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.l = np.array(l)
        self.b = np.array(b)
        self.dist = np.array(r)
        self.vx = np.array(vx)
        self.vy = np.array(vy)
        self.vz = np.array(vz)
        self.mass = np.array(mass)
        self.vlos = np.array(vlos)

        #NOTE: Any time you update the position data, you have to update
        #   the center of mass (and same for velocity and COMomentum)
        #   This is done automatically if flags.updateData is on
        self.centerOfMass = centerOfMass
        self.centerOfMomentum = centerOfMomentum

        self.msol = self.mass * mwahpy_glob.structToSol

        #ICRS information
        c = SkyCoord(l=self.l*u.degree, b=self.b*u.degree, frame='galactic')
        c_trans = c.transform_to('icrs')
        self.ra = c_trans.ra.degree
        self.dec = c_trans.dec.degree

        #angular momentum information
        self.lx = self.y * self.vz - self.z * self.vy
        self.ly = self.x * self.vz - self.z * self.vx
        self.lz = self.x * self.vy - self.y * self.vx
        self.lperp = (self.lx**2 + self.ly**2)**0.5
        self.ltot = (self.lx**2 + self.ly**2 + self.lz**2)**0.5

        #galactocentric information
        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5
        self.R = (self.x**2 + self.y**2)**0.5
        self.vgsr = self.vlos + 10.1*np.cos(self.b*np.pi/180)*np.cos(self.l*np.pi/180) + 224*np.cos(self.b*np.pi/180)*np.sin(self.l*np.pi/180) + 6.7*np.sin(self.b*np.pi/180)
        self.rad = (self.x*self.vx + self.y*self.vy + self.z*self.vz)/self.r
        self.rot = self.lz/(self.x**2 + self.y**2)**0.5

        #relative information

        #this shouldn't have *too* much overhead, but it could be added to the
        #   things that aren't computed until they're called
        #   If that's the case, adding 'xFromCOM' etc. is probably also a good idea
        self.distFromCOM = ((self.x - self.centerOfMass[0])**2 + (self.y - self.centerOfMass[1])**2 + (self.z - self.centerOfMass[2])**2)**0.5

        #-----------------------------------------------------------------------
        # HOUSEKEEPING
        #-----------------------------------------------------------------------

        #this has to be manually updated any time a new iterable quantity is added
        #to the Data class. This allows us to control what values are iterated over.
        self.indexList = ['id', 'x', 'y', 'z', 'l', 'b', 'dist', 'vx', 'vy', 'vz', \
                          'mass', 'vlos', 'msol', 'ra', 'dec', \
                          'lx', 'ly', 'lz', 'lperp', 'ltot', 'r', 'R', 'vgsr', 'rad', 'rot', \
                          'distFromCOM']
        self.index = 0

        #these should initially be set to false. If the user tries to get a value
        #that isn't yet calculated, then the getter calculates it. There's a
        #gigantic overhead on the rv/pm calculation and a pretty big overhead on
        #the energy calculation, so we avoid it until the user asks for those values.
        self.have_rvpm = False
        self.have_energy = False

    #---------------------------------------------------------------------------
    # METHODS
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # ITERATOR CONTROL
    #self.indexList allows us to ignore things that are not identically iterable in self.__dict__(),
    #such as the centers of mass and momentum

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index == len(self.indexList):
            raise StopIteration
        else:
            self.index += 1
            return self.indexList[self.index - 1]

    #---------------------------------------------------------------------------
    # These allow for the behavior we want from the class attributes, i.e. self.id == self['id']
    #__getattr__ and __getitem__ have qualifiers to make sure that loading rv/pm and energy
    #related values works as intended

    def __getattr__(self, i):
        if (not(self.have_rvpm) and i in ['rv', 'pmra', 'pmdec',  'pmtot', 'vtan']):
            self.calcrvpm()
        if (not(self.have_energy) and i in ['PE', 'KE', 'energy']):
            self.calcEnergy()
        return self.__dict__[i]

    def __getitem__(self, i):
        if (not(self.have_rvpm) and i in ['rv', 'pmra', 'pmdec',  'pmtot', 'vtan']):
            self.calcrvpm()
        if (not(self.have_energy) and i in ['PE', 'KE', 'energy']):
            self.calcEnergy()
        return self.__dict__[i]

    def __setitem__(self, i, val):
        self.__dict__[i] = val

    #---------------------------------------------------------------------------

    def __len__(self):
        return len(self.id)

    def update(self):
        self.centerOfMass = [np.sum(self.x*self.mass/sum(self.mass)), np.sum(self.y*self.mass/sum(self.mass)), np.sum(self.z*self.mass/sum(self.mass))]
        self.centerOfMomentum = [np.sum(self.vx*self.mass/sum(self.mass)), np.sum(self.vy*self.mass/sum(self.mass)), np.sum(self.vz*self.mass/sum(self.mass))]
        self.distFromCOM = ((self.x - self.centerOfMass[0])**2 + (self.y - self.centerOfMass[1])**2 + (self.z - self.centerOfMass[2])**2)**0.5

    #creates a deep copy of the Data object
    #this can't be done by iterating over the object, since comass etc. have to be copied as well
    def copy(self):
        out = Data()
        for key in self.__dict__.keys():
            if type(out[str(key)]) == type(np.array([])) or type(out[str(key)]) == type([]):
                out[str(key)] = self[str(key)].copy()
            else:
                out[str(key)] = self[str(key)]

        return out

    #---------------------------------------------------------------------------

    def calcrvpm(self):
        #the biggest source of overhead in this class is co.getrvpm
        #so, don't run it unless you have to.

        self.have_rvpm = True #make sure the getter doesn't try to run this again

        self.rv, self.pmra, self.pmdec = co.getrvpm(self.ra, self.dec, self.dist, self.vx, self.vy, self.vz)
        self.pmtot = (self.pmra**2 + self.pmdec**2)**0.5
        #4.848e-6 is arcsec->rad, 3.086e16 is kpc->km, and 3.156e7 is sidereal yr -> seconds
        self.vtan = 4.74*self.dist*self.pmtot #eq. to self.r*np.tan(self.pmtot*4.848e-6) * 3.086e16 / 3.156e7

        self.indexList = self.indexList + ['rv', 'pmra', 'pmdec',  'pmtot', 'vtan']

    def calcEnergy(self):
        #calculating the energy of every particle can generate some overhead,
        #so I've quarantined it with a flag.

        self.have_energy = True #make sure the getter doesn't try to run this again

        #in a logarithmic halo, the magnitude of the potential doesn't impact the result,
        #just the difference in potentials. So, you can specify a potential offset
        #to keep bound objects' total energy negative.
        PE = galpy.potential.evaluatePotentials(pot.pot, self.R * u.kpc, self.z*u.kpc, ro=8., vo=220.) + pot.energy_offset
        KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)

        #set attributes
        self.PE = PE
        self.KE = KE
        self.energy = PE + KE

        #allow iteration over these attributes
        self.indexList = self.indexList + ['PE', 'KE', 'energy']

    #---------------------------------------------------------------------------

    #cuts the first n entries from every attribute in the data structure
    def cutFirstN(self, n):
        for key in self:
            self[key] = self[key][n:]
        if flags.updateData:
            self.update()

    #cuts the last n entries from every attribute in the data structure
    def cutLastN(self, n):
        l = len(self) #length changes during this, so have to save it
        for key in self:
            self[key] = self[key][:l-n]
        if flags.updateData:
            self.update()

    #cut the data to only include the values at the given indices
    def take(self, indices):
        #indices: the indices you want to take, must be array-like
        for key in self:
            self[key] = np.take(self[key], indices)
        if flags.updateData:
            self.update()

    #resets the IDs of the data instance
    def resetIds(self):
        self.id = np.arange(len(self))

    #splits the Data into two new Data structures,
    #the first has the data points up to entry n and
    #the second has the data points after entry n
    def split(self, n):
        Data1 = self.copy()
        Data2 = self.copy()

        Data1.cutLastN(len(self) - n)
        Data2.cutFirstN(n)

        if flags.updateData:
            Data1.update()
            Data2.update()

        return Data1, Data2

    #splits the Data into a list of new Data structures,
    #where the Data is split every time ID wraps back to zero
    def splitAtIdWrap(self):

        outlist = []
        indices = np.where(self.id==0)[0] #1D list of arrays

        Data2 = self.copy()
        i = 1
        while i < len(indices): #pseudo recursive
            Data1, Data2 = Data2.split(indices[i] - indices[i-1])
            outlist.append(Data1)
            i += 1
        outlist.append(Data2)

        if flags.updateData:
            for d in outlist:
                d.update()

        return outlist

    #append a data object onto this one
    def appendData(self, d):
        #d: the data to append to this object
        for key, dkey in zip(self, d):
            self[key] = np.append(self[key], d[dkey])

        if flags.updateData:
            self.update()

    #append the nth item of another data object onto this one
    #WARNING: Extremely time intensive if used repeatedly, np.append is not a speedy function
    def appendPoint(self, d, n=0, id=None):
        #d: the data object with the item being appended
        #n: the index of the item in d to be appended
        #id: if not None, uses finds the first item with matching id and appends that
        if id:
            n = np.where(d['id'] == id)[0]
        for key, dkey in zip(self, d):
            self[key] = np.append(self[key], d[dkey][n])

        if flags.updateData:
            self.update()

    #---------------------------------------------------------------------------
    # OUTPUT
    #---------------------------------------------------------------------------

    #prints out a '.csv' file of the Data class structure
    def makeCSV(self, f_name):
        #f: the path for the output csv file

        f = open(f_name, 'w')

        #make the header
        #TODO: Print out COM's
        if flags.verbose:
            print('Writing header...')
        header = ''
        for key in self:
            header += (key + ',')
        header += '\n'
        f.write(header)

        #iterate through the data and print each line
        i = 0
        if flags.verbose:
            print('Printing data...')
        while i < len(self): #i avoided zip() here because it's messy to zip
        #                        like a billion different arrays, although zip()
        #                        is "better programming"
            if flags.progressBars:
                mwahpy_glob.progressBar(i, len(self))
            line = ''
            for key in self:
                line += (str(self[key][i]) + ',')
            line += '\n'
            f.write(line)
            i += 1

        print('Data output to ' + f_name)

    #---------------------------------------------------------------------------
    # MISCELLANEOUS
    #---------------------------------------------------------------------------

    #NOTE: Does not correctly adjust all parameters
    #might be better to rotate x, y, vx, vy and then re-initialize the data object
    def rotateAroundZAxis(self, theta, rad=False): #rotates counter-clockwise
        #rad: True if theta is in radians

        if not rad:
            theta *= np.pi/180

        cos = np.cos(theta)
        sin = np.sin(theta)

        R = np.array([[cos, -1*sin],
                      [sin, cos]])

        #update position
        xy = np.array([self.x, self.y])
        new_xy = np.matmul(R, xy)

        self.x = new_xy[0]
        self.y = new_xy[1]

        #update velocities
        vxvy = np.array([self.vx, self.vy])
        new_vxvy = np.matmul(R, vxvy)

        self.vx = new_vxvy[0]
        self.vy = new_vxvy[1]

        #update angular momenta
        self.lx = self.y * self.vz - self.z * self.vy
        self.ly = self.x * self.vz - self.z * self.vx
        self.lz = self.x * self.vy - self.y * self.vx

        #TODO: update other parameters that will be impacted by this (vlos, pm's)

        if flags.updateData:
            self.update()

    #TODO: subsets should return the same instance
    #   the user should copy the data instance themselves if they want
    #   that functionality. This way, we aren't wasting time copying the
    #   object every time you subset it

    #make an n-dimensional rectangular cut on the data
    #TODO: make an inverted method, i.e. cut out things within the bounds
    def subsetRect(self, axes, bounds):
        #axes ([str, ...]): the parameters/values that you are cutting on. Input as a list of strings
        #   The strings must be in self.indexList
        #bounds ([(float, float), ...]): The boundary conditions that you are cutting on
        #   given as a list of tuples, each with two floats.
        #   If you wish to only declare a minimum or maximum, then leave the other value as None
        #
        #Both axes and bounds can be a single input of their respective type
        #
        #Both axes and bounds can be given as a tuple of tuples, instead of a list
        #
        #This is substantially faster than the previous subset function, which performed like 30
        #np.appends for every index that fit the criteria in each axis.

        d = self.copy() #return a new array instead of altering the old one
        #probably slows things down a little but you should only have to make cuts once anyways

        if type(axes) not in [type([]), type((1,1))]: #make the input compatible if it is only 1 axis
            axes = [axes]
            bounds = [bounds]

        if len(axes) != len(bounds): #need the same number of axes and bounds
            raise Exception('Number of axes is ' + str(len(axes)) + ', but number of bounds is ' + str(len(bounds)))

        #make the first cut, gives us some place to start on for the indices
        if bounds[0][0] > bounds[0][1]: #make sure the bounds are input correctly
            raise Exception('First value in bound was larger than second value')

        indices = np.intersect1d(np.where(d[axes[0]] > bounds[0][0]), np.where(d[axes[0]] < bounds[0][1]))

        for a, b in zip(axes[1:], bounds[1:]): #already performed first cut
            if b[0] > b[1]: #make sure the bounds are input correctly
                raise Exception('First value in bound was larger than second value')

            #slowly wittle down the number of indices that fit the criteria
            indices = np.intersect1d(indices, np.intersect1d(np.where(d[a] > b[0]), np.where(d[a] < b[1])))

        #cut the sample down to the indices that work
        d.take(indices)
        d.update()

        return d

    #make a circular cut of the data in 2 dimensions
    def subsetCirc(self, ax, zy, rad, center):
        #ax, ay (str): the axes to cut on
        #rad (float): The radius of the circular cut
        #center (tuple(float, float)): The center around which to make the circular cut
        #
        #This is substantially faster than the previous subset function, which performed like 30
        #np.appends for every index that fit the criteria in each axis.

        d = self.copy() #return a new array instead of altering the old one
        #probably slows things down a little but you should only have to make cuts once anyways

        #get the indices that lie within the cut
        dist = ((d[ax] - center[0])**2 + (d[ay] - center[1])**2)**0.5
        indices = np.where(dist < rad)[0]

        #cut the sample down to the indices that work
        d.take(indices)
        d.update()

        return d

    #cut the data instance to n random stars from the data instance
    #uses the reservoir algorithm for single-pass random sampling
    def randSample(self, n):
        #n: the number of stars to sample

        reservoir = self.id[:n]
        for i in self.id[n:]:
            r = random.randint(0,len(self))
            if r <= n:
                reservoir[r] = self.id[r]

        self.take(reservoir)
        self.update()

    #cut the data instance to every nth star from the data instance
    def subsample(self, n):
        #n: take the star every n rows (n should be integer >= 1)

        i = 0
        j = 1
        indices = []
        while i < len(self):
            if j == n:
                indices.append(i)
                j = 1
            else:
                j += 1

        self.take(indices)
        self.update()


#===============================================================================
# FUNCTIONS INVOLVING DATA CLASSES
#===============================================================================



#===============================================================================
# UNIT TESTING
#===============================================================================
#Most of the tests in this file read in the test.out MW@h output file that is
#packaged along with the developer version of mwahpy. Without this file, these
#tests will fail. If a different file is used, these tests will also fail.

#NOTE: test.out does not have correct COM information, since it is a
#truncated version of a much larger MW@h output file.

#to run tests, run this file like you would run a regular python script

#to add more tests,
# 1) Find a test class that fits the description of the test you want to add,
#    then add that test to the class as a function beginning with "test", OR
# 2) Create a new class beginning with "test" and composed of just (unittest.TestCase)
#    and then add new functions that begin with "test_" to that class.
#    This should be used if new tests do not fit the descriptions of other
#    test classes

#after changing this file, it should be run in order to make sure all tests pass
#if tests don't pass, congratulations you broke something
#please fix it

prec = 8 #number of digits to round to when comparing values
#WARNING: Tests may (but are not expected to) fail at high levels of prec

class TestDataClass(unittest.TestCase):

    def testDataInitialize(self):
        d = Data() #check default initialization
        d = output_handler.readOutput('../test/test.out') #check loading from a file

    def testDataDict(self):
        d = output_handler.readOutput('../test/test.out')
        self.assertTrue(d.x[0] == d['x'][0])

        d.x[0] = 1
        self.assertTrue(d['x'][0] == 1)

    def testDataIter(self):
        d = output_handler.readOutput('../test/test.out')

        for key, k in zip(d, d.indexList):
            self.assertTrue(d[key][0] == d[k][0])
            d[key][0] = 1

        self.assertTrue(d.x[0] == d['x'][0] == 1)

        for key in d:
            d[key] = np.append(d[key], 0)

        self.assertTrue(d.x[-1] == d['x'][-1] == 0)

    def testCopy(self):
        d = output_handler.readOutput('../test/test.out')

        d2 = d.copy()
        test = d.x[0]
        d.x[0] = 0

        self.assertTrue(d.x[0] != d2.x[0])
        self.assertTrue(d2.x[0] == test)

    def testAppendPoint(self):
        d = output_handler.readOutput('../test/test.out')

        d2 = d.copy()

        d2.appendPoint(d, 5)

        self.assertTrue(d.x[5] == d2.x[-1])
        self.assertTrue(len(d2) == len(d)+1)

    def testSplit(self):
        d = output_handler.readOutput('../test/test.out')

        d1, d2 = d.split(5)
        self.assertTrue(d1.x[0] == d.x[0])
        self.assertTrue(d2.x[0] == d.x[5])

    def testSubsetRect(self):
        d = output_handler.readOutput('../test/test.out')
        dc = d.copy()

        dc.subsetRect('x', (-1,1))

        self.assertTrue(len(dc) == 7)

    def testCalcs(self): #this just makes sure that the values that
        #are not initially calculated on instantiation then run when the user
        #tries to get those values

        d = output_handler.readOutput('../test/test.out')

        print(len(d.rv))

        self.assertTrue(len(d.rv) == len(d))
        self.assertTrue(len(d.energy) == len(d))

class TestDataMethods(unittest.TestCase):

    def testUpdate(self):
        d = output_handler.readOutput('../test/test.out')
        d.update()

        #these values for the COM's are hard coded, unfortunately. But they are correct.
        self.assertTrue(len(d.centerOfMass) == 3)
        self.assertTrue(round(abs(d.centerOfMass[0] - 0.7835947518080879) + abs(d.centerOfMass[1] - 0.2947546230649471) + abs(d.centerOfMass[2] + 0.1318053758650839), prec) == 0)

        self.assertTrue(len(d.centerOfMomentum) == 3)
        self.assertTrue(round(abs(d.centerOfMomentum[0] - 6.001115213580641) + abs(d.centerOfMomentum[1] + 65.29652414026405) + abs(d.centerOfMomentum[2] + 26.462554427407223), prec) == 0)

        self.assertTrue(len(d.distFromCOM) == len(d))
        self.assertTrue(round(abs(d.distFromCOM[0] - ((d.centerOfMass[0] - d.x[0])**2 + (d.centerOfMass[1] - d.y[0])**2 + (d.centerOfMass[2] - d.z[0])**2)**0.5), prec) == 0)

#===============================================================================
# RUNTIME
#===============================================================================

if __name__ == '__main__':
    unittest.main()
