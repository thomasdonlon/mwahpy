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

import numpy as np
import coords as co
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import random
import galpy
import unittest
import os

import mwahpy_glob
import flags
import output_handler

#===============================================================================
# DATA CLASS
#===============================================================================

#AttrDict is used as a helper class in Data to allow referencing attributes
#as dict keys and vice-versa.
#this is probably a bad way to implement this but it works, and it's better than
#making Data inherit from dict, which was the other solution I was able to
#strum up
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
        self.rv, self.pmra, self.pmdec = co.getrvpm(self.ra, self.dec, self.dist, self.vx, self.vy, self.vz)
        self.pmtot = (self.pmra**2 + self.pmdec**2)**0.5
        #4.848e-6 is arcsec->rad, 3.086e16 is kpc->km, and 3.156e7 is sidereal yr -> seconds
        self.vtan = 4.74*self.dist*self.pmtot #eq. to self.r*np.tan(self.pmtot*4.848e-6) * 3.086e16 / 3.156e7

        #angular momentum information
        self.lx = self.y * self.vz - self.z * self.vy
        self.ly = self.x * self.vz - self.z * self.vx
        self.lz = self.x * self.vy - self.y * self.vx
        self.lperp = (self.lx**2 + self.ly**2)**0.5
        self.ltot = (self.lx**2 + self.ly**2 + self.lz**2)**0.5

        #galactocentric information
        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5
        self.vgsr = self.vlos + 10.1*np.cos(self.b*np.pi/180)*np.cos(self.l*np.pi/180) + 224*np.cos(self.b*np.pi/180)*np.sin(self.l*np.pi/180) + 6.7*np.sin(self.b*np.pi/180)
        self.rad = (self.x*self.vx + self.y*self.vy + self.z*self.vz)/self.r
        self.rot = self.lz/(self.x**2 + self.y**2)**0.5

        #relative information
        self.distFromCOM = ((self.x - self.centerOfMass[0])**2 + (self.y - self.centerOfMass[1])**2 + (self.z - self.centerOfMass[2])**2)**0.5

        #-----------------------------------------------------------------------
        if flags.calcEnergy:
            #get the energy info
            #calculating the energy of every particle can generate some overhead,
            #so I've quarantined it with a flag.

            #in a logarithmic halo, the magnitude of the potential doesn;t impact the result,
            #just the difference in potentials. So, you can specify a potential offset
            #to keep bound objects' total energy negative.
            PE = galpy.potential.evaluatePotentials(mwahpy_glob.pot, (self.x**2 + self.y**2)**0.5 * u.kpc, self.z*u.kpc, ro=8., vo=220.) - pot_offset
            KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)

            #set attributes
            self.PE = PE
            self.KE = KE
            self.energy = PE + KE

            #allow iteration over these attributes
            self.indexList = self.indexList + ['PE', 'KE', 'energy']

        #-----------------------------------------------------------------------
        # HOUSEKEEPING
        #-----------------------------------------------------------------------

        #this has to be manually updated any time a new iterable quantity is added
        #to the Data class. This allows us to control what values are iterated over.
        self.indexList = ['id', 'x', 'y', 'z', 'l', 'b', 'dist', 'vx', 'vy', 'vz', \
                          'mass', 'vlos', 'msol', 'ra', 'dec', 'rv', 'pmra', 'pmdec',  'pmtot', 'vtan', \
                          'lx', 'ly', 'lz', 'lperp', 'ltot', 'r', 'vgsr', 'rad', 'rad', \
                          'distFromCOM']
        self.index = 0

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

    def __getitem__(self, i):
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

    #cuts the first n entries from every attribute in the data structure
    def cutFirstN(self, n):
        for key in self:
            self[key] = self[key][n:]
        if flags.updateData:
            self.update()

    #cuts the last n entries from every attribute in the data structure
    def cutLastN(self, n):
        for key in self:
            self[key] = self[key][:n]
        if flags.updateData:
            self.update()

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
        indices = np.where(self.id==0)[0] #1D list of arrays,

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

#===============================================================================
# FUNCTIONS INVOLVING DATA CLASSES
#===============================================================================
#TODO: Split a Data class based on when id loops back to 0 (meaning that you used
#   manual body input as well as generating a dwarf, or multiple manual body inputs)

#subset: data object -> data object
#takes in a data object and cuts it according to the specified parameters
#USAGE:
#   You MUST specify values for radius and center, or set rect=True
#   if rect=True, y, xbounds, and ybounds must be specified
#   if rect=False and y is given, then do 2D radial cut
#TODO: Allow an arbitrary length array of axes and bounds, instead of this
#   manually setting rectangular or not nonsense. Should make a range subset
#   function as well as a rectangular one
#TODO: Calculate the indices that are correct, then append all the points with
#   those indices to each attribute. This appendPoint stuff takes way too long
#   because numpy.append sucks
def subset(data, x, y=None, rect=False, center=None, radius=None, xbounds=None, ybounds=None):
    #data (Data): the data object being cut
    #x (str): the x-axis parameter
    #y (str): the y-axis parameter
    #   any Data array_dict value can be used for x and y, e.g. "ra", "dec", "x", "vx", etc.
    #rect (bool): if True, do a 2D rectangular subset cut instead
    #----------------------
    #1D cut:
    #----------------------
    #   center (float):
    #   radius (float): data points under <radius> away from <center> on the <x> axis
    #                   are added to a new Data object
    #----------------------
    #2D cut: (rect=False)
    #----------------------
    #   center (tuple(float, float)):
    #   radius (float): data points under <radius> away from <center> on both axes
    #                   are added to a new Data object
    #----------------------
    #2D cut: (rect=True)
    #----------------------
    #   xbounds (tuple(float, float)): x range of allowed subsection
    #   ybounds (tuple(float, float)): y range of allowed subsection
    #   if xbounds[0] !< xbounds[1] etc, then the routine will throw an exception

    data_out = Data()

    if rect: #Rectangular 2D cut
        if not(y) or not(xbounds) or not(ybounds): #check for necessary inputs
            raise Exception('Must provide <y>, <xbounds>, and <ybounds> if <rect>==True')
        elif (xbounds[0] >= xbounds[1]):
            raise Exception('First value in <xbounds> was greater than the second value')
        elif (ybounds[0] >= ybounds[1]):
            raise Exception('First value in <ybounds> was greater than the second value')
        else: #do rectangular cut
            i = 0
            while i < len(data):
                if flags.progressBars:
                    mwahpy_glob.progressBar(i, len(data[x]))
                #check if within bounds
                if xbounds[0] < data[x][i] < xbounds[1] and ybounds[0] < data[y][i] < ybounds[1]:
                    data_out.appendPoint(data, i)
                i+=1
            if flags.verbose:
                print('\n' + str(len(data_out)) + ' objects found in bounds')
            if flags.updateData:
                data_out.update()
            return data_out

    else: #Not rectangular cut
        if radius==None or center==None: #radius and/or center weren't provided
            raise Exception('Must provide <center> and <radius> for a non-rectangular cut')
        else:
            if y: #do 2D cut
                if not(type(center) is tuple): #make sure center is a tuple
                    raise Exception('Must provide tuple for <center> if <y> is provided')
                else:
                    i = 0
                    while i < len(data[x]):
                        if flags.progressBars:
                            mwahpy_glob.progressBar(i, len(data[x]))
                        if radius >= ((array_dict[x][i] - center[0])**2 + (array_dict[y][i] - center[1])**2)**0.5:
                            data_out.appendPoint(data, i)
                        i+=1
                    if flags.verbose:
                        print('\n' + str(len(data_out)) + ' objects found in bounds')
                    if flags.updateData:
                        data_out.update()
                    return data_out
            else: #do 1D cut
                i = 0
                while i < len(data[x]):
                    if flags.progressBars:
                        mwahpy_glob.progressBar(i, len(data[x]))
                    if radius >= abs(data[x][i] - center):
                        data_out.appendPoint(data, i)
                    i+=1
                if flags.verbose:
                    print('\n' + str(len(data_out)) + ' objects found in bounds')
                if flags.updateData:
                    data_out.update()
                return data_out

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
        d.x[0] = 0

        self.assertTrue(d.x[0] != d2.x[0])

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
