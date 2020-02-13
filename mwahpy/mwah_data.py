'''
MilkyWay@home Python Package (mwahpy) v1.0
Copyright Tom Donlon 2020, RPI

See provided readme and documentation for information
or reach out to me on github
user: thomasdonlon

The contents of this file are focused on the Data class, which is used for storage of
imported data from N-body output files.
'''

import numpy as np
import coord_trans as ct
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import random

import galpy
from galpy.orbit import Orbit
from galpy.potential import HernquistPotential
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import PlummerPotential

#===========================================
# CONSTANTS
#===========================================

m_bulge = 3.4e10*u.solMass #solar masses
m_disk = 1.0e11*u.solMass
v_halo = 74.61*u.km/u.s #km/s

G = 6.67e-11*u.m**3/(u.kg*u.s**2)

pot_bulge = HernquistPotential(amp=2*m_bulge, a=0.7*u.kpc, ro=8., vo=220.)
pot_disk = MiyamotoNagaiPotential(amp=G*m_disk, a=6.5*u.kpc, b=0.26*u.kpc, ro=8., vo=220.)
pot_halo = LogarithmicHaloPotential(amp=2*v_halo**2, q=1., core=12.0*u.kpc, ro=8., vo=220.)
pot = [pot_bulge, pot_disk, pot_halo]

m_plummer = 1e9*u.solMass
r_scale_plummer = 3*u.kpc
plummer_pot = PlummerPotential(amp=G*m_plummer, b=r_scale_plummer, ro=10*u.kpc, vo=20*u.km/u.s)

struct_to_sol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

#===========================================

class Data():

    def __init__(self, id_val=[], x=[], y=[], z=[], l=[], b=[], r=[], vx=[], vy=[], vz=[], mass=[], vlos=[], centerOfMass=[], centerOfMomentum=[], *args, **kwargs):
        #all this is typically specified by the readOutput function in mwahpy.py
        #readOutput is the preferred way to input data to this data structure
        #but if you're really feeling adventurous you can always do it yourself

        #this data structure allows access to a dictionary like an attribute, i.e.
        #   d = AttrDict({'a':1, 'b':2})
        #   d['a'] == d.a //returns True
        #it seems like magic but it works by accessing the way python handles attributes
        #within each class. There are pros and cons, the most notable thing is
        #that it creates a memory leak in python <3.2.3
        super(Data, self).__init__(*args, **kwargs)
        self.__dict__ = self

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

        self.centerOfMass = centerOfMass
        self.centerOfMomentum = centerOfMomentum

        self.msol = self.mass * struct_to_sol

        #ICRS information
        c = SkyCoord(l=self.l*u.degree, b=self.b*u.degree, frame='galactic')
        c_trans = c.transform_to('icrs')
        self.ra = c_trans.ra.degree
        self.dec = c_trans.dec.degree
        self.rv, self.pmra, self.pmdec = ct.getrvpm(self.ra, self.dec, self.dist, self.vx, self.vy, self.vz)
        self.pmtot = (self.pmra**2 + self.pmdec**2)**0.5
        #4.848e-6 is arcsec->rad, 3.086e16 is kpc->km, and 3.156e7 is sidereal yr -> seconds
        self.vtan = 4.74*self.dist*self.pmtot #eq. to self.r*np.tan(self.pmtot*4.848e-6) * 3.086e16 / 3.156e7

        #galactocentric information
        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5
        self.vgsr = self.vlos + 10.1*np.cos(self.b*np.pi/180)*np.cos(self.l*np.pi/180) + 224*np.cos(self.b*np.pi/180)*np.sin(self.l*np.pi/180) + 6.7*np.sin(self.b*np.pi/180)
        self.rad = (self.x*self.vx + self.y*self.vy + self.z*self.vz)/self.r
        self.rot = self.lz/(self.x**2 + self.y**2)**0.5

        #angular momentum information
        self.lx = self.y * self.vz - self.z * self.vy
        self.ly = self.x * self.vz - self.z * self.vx
        self.lz = self.x * self.vy - self.y * self.vx
        self.lperp = (self.lx**2 + self.ly**2)**0.5
        self.ltot = (self.lx**2 + self.ly**2 + self.lz**2)**0.5

    #calculating the energy of every particle can generate some overhead,
    #so I've quarantined it to its own function.
    #Running this function will set attributes for the calculated energies, so
    #you don't have to re-calculate every time you want to get an energy
    def energy(self, potential_offset = 0):
        #get the energy info

        #in a logarithmic halo, the magnitude of the potential doesn;t impact the result,
        #just the difference in potentials. So, you can specify a potential offset
        #to keep bound objects' total energy negative.
        PE = galpy.potential.evaluatePotentials(pot, (self.x**2 + self.y**2)**0.5 * u.kpc, self.z*u.kpc, ro=8., vo=220.) - pot_offset
        KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)

        #set attributes
        self.PE = PE
        self.KE = KE
        self.energy = PE + KE

        #update array_dict
        self.array_dict['PE'] = self.PE
        self.array_dict['KE'] = self.KE
        self.array_dict['energy'] = self.energy

        #returns the energy so that the function semantically makes sense
        return self.energy

    #cuts the first n entries from every attribute in the data structure
    def cutFirstN(self, n):
        for key in self.keys():
            self.key = self.key[n:]

    #cuts the last n entries from every attribute in the data structure
    def cutLastN(self, n):
        for key in self.keys():
            self.key = self.key[:len(self.key) - n]

    #splits the Data into two new Data structures,
    #the first has the data points up to entry n and
    #the second has the data points after entry n
    def split(self, n):
        Data1 = self.copy()
        Data2 = self.copy()

        Data1.cutLastN(n)
        Data2.cutFirstN(n)

        return Data1, Data2
