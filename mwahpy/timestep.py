'''
The contents of this file are focused on the Timestep class, which is used for storage of
imported data from N-body output files.
'''

#TODO: UNIT TESTS
#TODO: Create sort function that sorts based on id (not really importamt)

#===============================================================================
# IMPORTS
#===============================================================================

#external imports
import numpy as np
import coords as co
from astropy.coordinates import SkyCoord
import astropy.units as u
import random
import galpy.potential
import unittest

#mwahpy imports
import mwahpy_glob
import flags
import plot
import pot

#===============================================================================
# TIMESTEP CLASS
#===============================================================================

#AttrDict is used as a helper class in Timestep to allow referencing attributes
#as dict keys and vice-versa.
#this is probably a bad way to implement this but it works, and it's better than
#making Timestep inherit from dict, which was the other solution I was able to strum up
class AttrDict(dict):
    def __init__(self):
        super(AttrDict, self).__init__()
        self.__dict__ = self

class Timestep():

    def __init__(self, id_val=[], x=[], y=[], z=[], vx=[], vy=[], vz=[], mass=[], centerOfMass=[0, 0, 0], centerOfMomentum=[0, 0, 0], potential=None, time=None, nbody=None):
        #all this is typically specified by the readOutput function in output_handler
        #readOutput is the preferred way to input data to this data structure
        #but if you're really feeling adventurous you can always do it yourself

        #-----------------------------------------------------------------------
        # HOUSEKEEPING
        #-----------------------------------------------------------------------

        #this Timestep structure allows access to a dictionary like an attribute, i.e.
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
        self.vx = np.array(vx)
        self.vy = np.array(vy)
        self.vz = np.array(vz)
        self.mass = np.array(mass)

        #NOTE: Any time you update the position data, you have to update
        #   the center of mass (and same for velocity and COMomentum)
        #   This is done automatically if flags.autoUpdate is on
        #   If you are manually screwing around with the provided values,
        #   you will need to update it manually
        self.centerOfMass = centerOfMass
        self.centerOfMomentum = centerOfMomentum

        #-----------------------------------------------------------------------
        # HOUSEKEEPING
        #-----------------------------------------------------------------------

        #if the time in Gyr is known for this simulation (i.e. is provided),
        #it is saved here.
        self.time = time

        #if the Timestep is part of a Nbody object, it saves that information here
        self.nbody = nbody

        #just in case the user wants to specify a different potential than the
        #mwahpy default, this is included
        self.potential = potential

        #this has to be manually updated any time a new iterable quantity is added
        #to the Timestep class. This allows us to control what values are iterated over.
        self.indexList = ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass']
        self.index = 0

        #these should initially be set to false. If the user tries to get a value
        #that isn't yet calculated, then the getter calculates it. There's a
        #gigantic overhead on the rv/pm calculation and a pretty big overhead on
        #the energy calculation, so we avoid it until the user asks for those values.
        #The basic values are included here so that updating these values are easy
        self.have_basic = False
        self.have_rvpm = False
        self.have_energy = False

        #these flags should also initially be set to false. If the position or
        #velocity of one or more particles is changed at any point, the
        #corresponding flag should be swapped. This tells update() how much it
        #needs to update instead of running the whole routine every time, which
        #would be very time consuming.
        #these flags allow all calculated attributes to be updated whenever
        #changes are made, and at a reasonable speed.
        self.changed_pos = False
        self.changed_vel = False

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
        if (not(self.have_basic) and i in ['msol', 'l', 'b', 'ra', 'dec', 'dist', 'lx', 'ly', 'lz', 'lperp', 'ltot', 'r', 'R', 'vlos', 'vgsr', 'rad', 'rot', 'distFromCOM']):
            self.calcBasic()
        if (not(self.have_rvpm) and i in ['rv', 'pmra', 'pmdec',  'pmtot', 'vtan']):
            self.calcrvpm()
        if (not(self.have_energy) and i in ['PE', 'KE', 'energy']):
            self.calcEnergy()
        return self.__dict__[i]

    def __getitem__(self, i):
        if (not(self.have_basic) and i in ['msol', 'l', 'b', 'ra', 'dec', 'dist', 'lx', 'ly', 'lz', 'lperp', 'ltot', 'r', 'R', 'vlos', 'vgsr', 'rad', 'rot', 'distFromCOM']):
            self.calcBasic()
        if (not(self.have_rvpm) and i in ['rv', 'pmra', 'pmdec',  'pmtot', 'vtan']):
            self.calcrvpm()
        if (not(self.have_energy) and i in ['PE', 'KE', 'energy']):
            self.calcEnergy(potential=self.potential)
        return self.__dict__[i]

    def __setitem__(self, i, val):
        self.__dict__[i] = val

        #set the flags that say we need to update the corresponding attributes
        #based on what was changed
        if i in ['x', 'y', 'z']:
            self.changed_pos = True
        elif i in ['vx', 'vy', 'vz']:
            self.changed_vel = True

    #---------------------------------------------------------------------------
    #a few routines that are important for the functionality of Timestep, but
    #don't really fit in elsewhere

    def __len__(self):
        return len(self.id)

    #replaces this Timestep with a blank Timestep
    def reset(self):
        self.id = np.array([])
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.vx = np.array([])
        self.vy = np.array([])
        self.vz = np.array([])
        self.mass = np.array([])

        self.centerOfMass = [0,0,0]
        self.centerOfMomentum = [0,0,0]

        self.time = None
        self.nbody = None
        self.potential = None

        self.indexList = ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass']
        self.index = 0

        self.have_basic = False
        self.have_rvpm = False
        self.have_energy = False

        self.changed_pos = False
        self.changed_vel = False

    def update(self, force=False):
        #only update the necessary attributes based on what changed
        #if force==True, then force updating everything
        if force:
            self.changed_pos = True
            self.changed_vel = True

        if self.changed_pos:
            self.centerOfMass = [np.sum(self.x*self.mass/sum(self.mass)), np.sum(self.y*self.mass/sum(self.mass)), np.sum(self.z*self.mass/sum(self.mass))]
        if self.changed_vel:
            self.centerOfMomentum = [np.sum(self.vx*self.mass/sum(self.mass)), np.sum(self.vy*self.mass/sum(self.mass)), np.sum(self.vz*self.mass/sum(self.mass))]

        if self.changed_pos or self.changed_vel:
            if self.have_basic:
                self.calcBasic()
            if self.have_rvpm:
                self.calcrvpm()
            if self.have_energy:
                self.calcEnergy(potential=self.potential)

        #the positions and velocities have not changed since updating
        #(since we just updated)
        changed_pos = False
        changed_vel = False

    #creates a deep copy of the Timestep object
    #this can't be done by iterating over the object, since comass etc. have to be copied as well
    def copy(self):
        out = Timestep()
        if self.have_basic:
            out.calcBasic()
        if self.have_rvpm:
            out.calcrvpm()
        if self.have_energy:
            out.calcEnergy()
        for key in self.__dict__.keys():
            if type(out[str(key)]) == type(np.array([])) or type(out[str(key)]) == type([]):
                out[str(key)] = self[str(key)].copy()
            else:
                out[str(key)] = self[str(key)]

        return out

    #---------------------------------------------------------------------------

    #the following methods are used to calculate values only if the user asks for them
    #this dramatically increases the speed of the package if the user is only using
    #a few values that need to be calculated

    def calcBasic(self):
        #none of these calculations take very long, so it's fair to group them up
        #this is only taken out of the original initialization because then you can update them
        #if you end up changing the provided values

        #this could be split up even more if any of the component value calculations
        #ends up taking a significant amount of time, but I doubt that will ever be the case

        if flags.verbose:
            print('Calculating basic values...')

        self.msol = self.mass * mwahpy_glob.structToSol

        #position information
        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5
        self.dist = ((self.x + 8)**2 + self.y**2 + self.z**2)**0.5
        self.R = (self.x**2 + self.y**2)**0.5

        #galactic coordinate information
        self.l = np.arctan2(self.y, self.x)*180/np.pi
        self.b = np.arcsin(self.z/self.r)*180/np.pi

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

        #velocity information
        self.vgsr = ((self.x+8)*self.vx + self.y*self.vy + self.x*self.vz)/self.dist
        self.vlos = self.vgsr - 10.1*np.cos(self.b*np.pi/180)*np.cos(self.l*np.pi/180) - 224*np.cos(self.b*np.pi/180)*np.sin(self.l*np.pi/180) - 6.7*np.sin(self.b*np.pi/180)
        self.rad = (self.x*self.vx + self.y*self.vy + self.z*self.vz)/self.r
        self.rot = self.lz/(self.x**2 + self.y**2)**0.5
        self.vR = (self.x*self.vx + self.y*self.vy)/self.R

        #relative information
        self.distFromCOM = ((self.x - self.centerOfMass[0])**2 + (self.y - self.centerOfMass[1])**2 + (self.z - self.centerOfMass[2])**2)**0.5

        if not(self.have_basic): #only run if first time running method, not if updating
            self.indexList = self.indexList + ['msol', 'l', 'b', 'ra', 'dec', 'dist', 'lx', 'ly', 'lz', 'lperp', 'ltot', 'r', 'R', 'vlos', 'vgsr', 'rad', 'rot', 'distFromCOM']

        self.have_basic = True #make sure the getter doesn't try to run this again

    def calcrvpm(self):
        #the biggest source of overhead in this class is co.getrvpm
        #so, don't run it unless you have to.

        if flags.verbose:
            print('Calculating proper motion values...')

        self.pmra, self.pmdec = co.getrvpm(self.ra, self.dec, self.dist, self.vx, self.vy, self.vz)[1:]
        #already get rv from vlos, don't need to save it as something else
        self.pmtot = (self.pmra**2 + self.pmdec**2)**0.5
        #4.848e-6 is arcsec->rad, 3.086e16 is kpc->km, and 3.156e7 is sidereal yr -> seconds
        self.vtan = 4.74*self.dist*self.pmtot #eq. to self.r*np.tan(self.pmtot*4.848e-6) * 3.086e16 / 3.156e7

        if not(self.have_rvpm): #only run if first time running method, not if updating
            self.indexList = self.indexList + ['pmra', 'pmdec',  'pmtot', 'vtan']

        self.have_rvpm = True #make sure the getter doesn't try to run this again

    def calcEnergy(self, potential=None):
        #calculating the energy of every particle can generate some overhead,
        #so I've quarantined it with a flag.

        if flags.verbose:
            print('Calculating energy values...')

        if potential == None:
            potential = pot.pot

        #in a logarithmic halo, the magnitude of the potential doesn't impact the result,
        #just the difference in potentials. So, you can specify a potential offset
        #to keep bound objects' total energy negative.
        PE = galpy.potential.evaluatePotentials(potential, self.R * u.kpc, self.z*u.kpc, ro=8., vo=220.) + pot.energy_offset
        KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)

        #set attributes
        self.PE = PE
        self.KE = KE
        self.energy = PE + KE

        #allow iteration over these attributes
        if not(self.have_energy):  #only run if first time running method, not if updating
            self.indexList = self.indexList + ['PE', 'KE', 'energy']

        self.have_energy = True #make sure the getter doesn't try to run this again

    #---------------------------------------------------------------------------

    #cuts the first n entries from every attribute in the Timestep structure
    def cutFirstN(self, n):
        for key in self:
            self[key] = self[key][n:]
        if flags.autoUpdate:
            self.update()

    #cuts the last n entries from every attribute in the Timestep structure
    def cutLastN(self, n):
        l = len(self) #length changes during this, so have to save it
        for key in self:
            self[key] = self[key][:l-n]
        if flags.autoUpdate:
            self.update()

    #cut the data to only include the values at the given indices
    def take(self, indices):
        #indices: the indices you want to take, must be array-like
        for key in self:
            self[key] = np.take(self[key], indices)
        if flags.autoUpdate:
            self.update()

    #resets the IDs of the Timestep instance
    def resetIds(self):
        self.id = np.arange(len(self))

    #splits the Timestep into two new Timestep structures,
    #the first has the data points up to entry n and
    #the second has the data points after entry n
    def split(self, n):
        Timestep1 = self.copy()
        Timestep2 = self.copy()

        Timestep1.cutLastN(len(self) - n)
        Timestep2.cutFirstN(n)

        if flags.autoUpdate:
            Timestep1.update()
            Timestep2.update()

        return Timestep1, Timestep2

    #splits the Timestep into a list of new Timestep structures,
    #where the Timestep is split every time ID wraps back to zero
    #TODO: the dwarf id start at 1, not 0
    def splitAtIdWrap(self):

        outlist = []
        indices = np.where(self.id==0)[0] #1D list of arrays

        Timestep2 = self.copy()
        i = 1
        while i < len(indices): #pseudo-recursive
            Timestep1, Timestep2 = Timestep2.split(indices[i] - indices[i-1])
            outlist.append(Timestep1)
            i += 1
        outlist.append(Timestep2)

        if flags.autoUpdate:
            for t in outlist:
                t.update()

        return outlist

    #append a Timestep object onto this one
    def appendTimestep(self, t):
        #t: the Timestep to append to this object
        for key, dkey in zip(self, t):
            self[key] = np.append(self[key], t[dkey])

        if flags.autoUpdate:
            self.update()

    #append the nth item of another Timestep object onto this one
    #WARNING: Extremely time intensive if used repeatedly, np.append is not a speedy function
    #   you should always try to append a timestep instead if possible
    def appendPoint(self, t, n=0, id=None):
        #t: the Timestep object with the item being appended
        #n (optional): the index of the item in t to be appended
        #id (optional): if not None, uses finds the first item with matching id and appends that
        if id:
            n = np.where(t['id'] == id)[0]
        for key, tkey in zip(self, t):
            self[key] = np.append(self[key], t[tkey][n])

        if flags.autoUpdate:
            self.update()

    #centers the timestep so that the COMs are now both zero vectors
    #this can be useful if running live simulations, which can pick up nonzero
    #   overall velocities when they relax from unstable virial equilibrium ICs
    #Conserves the physics of the system under a change of reference frame
    def recenter(self):
        #center the positions of the particles
        self.x = self.x - self.centerOfMass[0]
        self.y = self.y - self.centerOfMass[1]
        self.z = self.z - self.centerOfMass[2]

        #center the velocities of the particles
        self.vx = self.vx - self.centerOfMomentum[0]
        self.vy = self.vy - self.centerOfMomentum[1]
        self.vz = self.vz - self.centerOfMomentum[2]

        if flags.autoUpdate:
            self.update()

    #centers the timestep so that the COMs of all individual components are now zero vectors
    #this can be useful if running live simulations, which can pick up nonzero
    #   overall velocities when they relax from unstable virial equilibrium ICs
    #WARNING: Does not conserve the physics of the system. This is more for testing
    #   or debugging purposes, or if you really know what you're doing.
    def recenterEachComponent(self):
        #center the positions of the particles
        data_list = self.splitAtIdWrap()

        for i in data_list:
            i.recenter()

        self.reset()
        for i in data_list:
            self.appendTimestep(i)

        if flags.autoUpdate:
            self.update()

    #---------------------------------------------------------------------------
    # OUTPUT
    #---------------------------------------------------------------------------

    def printParticle(self, n, dec=8):
        #n (int): the location of the particle that you want the information for
        #TODO: allow for using the id to find the particle
        print('Printing data for Particle '+str(n)+':')

        outstr = '('
        for key in self:
            outstr = outstr + key + ':' + str(round(self[key][n],dec)) + ', '
        outstr = outstr + ')'
        print(outstr)

    #---------------------------------------------------------------------------
    # PLOTTING
    #---------------------------------------------------------------------------
    #NOTE: the actual plotting routines should be implemented in plot.py
    #this is just to allow you to access the routines as timestep.<method>()

    def scatter(self, x, y, **kwargs):
        plot.scatter(self, x, y, **kwargs)

    def hist2d(self, x, y, **kwargs):
        plot.hist2d(self, x, y, **kwargs)

    #---------------------------------------------------------------------------
    # MISCELLANEOUS
    #---------------------------------------------------------------------------

    #make an n-dimensional rectangular cut on the data
    #TODO: make an inverted method, i.e. cut out things within the bounds
    def subsetRect(self, axes, bounds):
        #axes ([str, ...]): the parameters/values that you are cutting on. Input as a list of strings
        #   The strings must be in self.indexList
        #bounds ([(float, float), ...]): The boundary conditions that you are cutting on
        #   given as a list of tuples, each with two floats.
        #   If you wish to only declare a minimum or maximum, then leave the other value as None
        #
        #Both axes and bounds can be given as a tuple instead of a list
        #
        #This is substantially faster than the previous subset function, which performed like 30
        #np.appends for every index that fit the criteria in each axis.

        if type(axes) not in [type([]), type((1,1))]: #make the input compatible if it is only 1 axis
            raise Exception('Axes must be a list or a tuple, but was of type ' + str(type(axes)))

        if len(axes) != len(bounds): #need the same number of axes and bounds
            raise Exception('Number of axes is ' + str(len(axes)) + ', but number of bounds is ' + str(len(bounds)))

        #make the first cut, gives us some place to start on for the indices
        if bounds[0][0] > bounds[0][1]: #make sure the bounds are input correctly
            raise Exception('First value in bound was larger than second value')

        indices = np.intersect1d(np.where(self[axes[0]] > bounds[0][0]), np.where(self[axes[0]] < bounds[0][1]))

        for a, b in zip(axes[1:], bounds[1:]): #already performed first cut
            if b[0] > b[1]: #make sure the bounds are input correctly
                raise Exception('First value in bound was larger than second value')

            #slowly wittle down the number of indices that fit the criteria
            indices = np.intersect1d(indices, np.intersect1d(np.where(self[a] > b[0]), np.where(self[a] < b[1])))

        #cut the sample down to the indices that work
        self.take(indices)
        self.update()

    #TODO: make n-dimensional like subsetRect
    #make a circular cut of the data in 2 dimensions
    def subsetCirc(self, axes, rads, centers):
        #axs ([str]): the axes to cut on
        #rads ([float]): The radiii of the circular cut (along each axis)
        #center ([float], len=#axes): The center around which to make the circular cut
        #
        #This is substantially faster than the previous subset function, which performed like 30
        #np.appends for every index that fit the criteria in each axis.

        #get the indices that lie within the cut
        dist = 0
        for a, r, c in zip(axes, rads, centers):
            dist += (self[a] - c)**2/r**2

        indices = np.where(dist < 1)[0]

        #cut the sample down to the indices that work
        self.take(indices)
        self.update()

    #cut the Timestep instance to n random stars from the Timestep instance
    #uses the reservoir algorithm for single-pass random sampling
    def randSample(self, n):
        #n: the number of stars to sample

        reservoir = np.arange(0,n)
        for i in np.arange(n,len(self)-1):
            r = random.randint(0,len(self)-1)
            if r < n:
                reservoir[r] = i

        self.take(reservoir.astype('int'))
        self.update()

    #cut the Timestep instance to every nth star from the Timestep instance
    def subsample(self, n, offset=0):
        #n: take the star every n rows (n should be integer >= 1)

        i = offset
        indices = []
        while i < len(self):
            indices.append(i)
            i += n

        self.take(indices)
        self.update()


#===============================================================================
# FUNCTIONS INVOLVING TIMESTEP CLASSES
#===============================================================================

    #Timestep -> np array of floats
    #computes the energy of each particle based only on the self-gravity
    #of the particles and their velocities w.r.t. the COM's of the particles.
    #slow because I am just implementing straight N^2 gravity calcs, but
    #not sure if there's a better way to do this
    #NOTE: This is energy per unit mass, not energy
    #NOTE: Only returns initial energy of the particles in a dwarf when run
    #      on a Timestep where the entire dwarf is still bound (e.g. the
    #      beginning of a simulation), otherwise output is ~meaningless
    def getSelfEnergies(t):
        #t: the Timestep

        gc = 4.30091e-3 * mwahpy_glob.structToSol / 1000 #Newton's gravitational constant in units of

        t = t.copy() #don't hurt the poor innocent Timestep
        t.recenter() #recenters both positions and velocities

        pot_energies = []

        #compute the self gravitational potential energy of each particle
        for x1, y1, z1 in zip(t.x, t.y, t.z):
            pe = 0
            for x2, y2, z2, m in zip(t.x, t.y, t.z, t.mass):
                r = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
                pe -= gc * m / r #potential energy is negative
            pot_energies.append(pe)
        pot_energies = np.array(pot_energies) #do it this way because np.append is horrendously slow

        #compute the kinetic energy of each particle (within reference frame of
        #the COM's of the particles)
        kin_energies = t.vx**2 + t.vy**2 + t.vz**2

        energies = pot_energies + kin_energies

        return energies
