#purpose: to take output from a MW@h .out file and produce workable data/plots to look at the resulting output in meaningful ways
#this is still too hard-coded for my liking, but it'll have to do for now
#i.e. if you want to add new attributes to the data class then you manually have to go through and fix the appending functions

import matplotlib.pyplot as plt
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

class Data():

    def __init__(self, id_val=[], x=[], y=[], z=[], l=[], b=[], r=[], vx=[], vy=[], vz=[], mass=[], vlos=[], centerOfMass=[], centerOfMomentum=[], pot_offset=0):
        #these should all be lists of floats
        self.id = np.array(id_val)
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.l = np.array(l)
        self.b = np.array(b)
        self.r = np.array(r)
        self.vx = np.array(vx)
        self.vy = np.array(vy)
        self.vz = np.array(vz)
        self.mass = np.array(mass)
        self.vlos = np.array(vlos)
        self.centerOfMass = centerOfMass
        self.centerOfMomentum = centerOfMomentum

        self.msol = self.mass * struct_to_sol
        c = SkyCoord(l=self.l*u.degree, b=self.b*u.degree, frame='galactic')
        c_trans = c.transform_to('icrs')
        self.ra = c_trans.ra.degree
        self.dec = c_trans.dec.degree
        self.d = (self.x**2 + self.y**2 + self.z**2)**0.5 #TODO: correct d vs r
        self.vgsr = self.vlos + 10.1*np.cos(self.b*np.pi/180)*np.cos(self.l*np.pi/180) + 224*np.cos(self.b*np.pi/180)*np.sin(self.l*np.pi/180) + 6.7*np.sin(self.b*np.pi/180)

        self.rv, self.pmra, self.pmdec = ct.getrvpm(self.ra, self.dec, self.r, self.vx, self.vy, self.vz)
        self.pmtot = (self.pmra**2 + self.pmdec**2)**0.5
        #4.848e-6 is arcsec->rad, 3.086e16 is kpc->km, and 3.156e7 is sidereal yr -> seconds
        self.vtan = 4.74*self.r*self.pmtot #self.r*np.tan(self.pmtot*4.848e-6) * 3.086e16 / 3.156e7

        #get the angular momentum info
        self.lx = self.y * self.vz - self.z * self.vy
        self.ly = self.x * self.vz - self.z * self.vx
        self.lz = self.x * self.vy - self.y * self.vx
        self.lperp = (self.lx**2 + self.ly**2)**0.5
        self.ltot = (self.lx**2 + self.ly**2 + self.lz**2)**0.5

        #galactocentric velocity information
        self.rad = (self.x*self.vx + self.y*self.vy + self.z*self.vz)/self.r
        self.rot = self.lz/(self.x**2 + self.y**2)**0.5

        #get the energy info
        PE = galpy.potential.evaluatePotentials(pot, (self.x**2 + self.y**2)**0.5 * u.kpc, self.z*u.kpc, ro=8., vo=220.) - pot_offset
        KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)
        self.energy = PE + KE

        self.array_dict = {'id':self.id, 'x':self.x, 'y':self.y, 'z':self.z, 'l':self.l, 'b':self.b, 'r':self.r, 'vx':self.vx, \
                            'vy':self.vy, 'vz':self.vz, 'mass':self.mass, 'vtan':self.vtan, 'vlos':self.vlos, \
                            'msol':self.msol, 'ra':self.ra, 'dec':self.dec, 'd':self.d, 'vgsr':self.vgsr, \
                            'energy':self.energy, 'ltot':self.ltot, 'lz':self.lz, 'lperp':self.lperp, 'pmra':self.pmra, \
                            'pmdec':self.pmdec, 'rad':self.rad}

    def cutFirstN(self, n):
        self.id = self.id[n:]
        self.x = self.x[n:]
        self.y = self.y[n:]
        self.z = self.z[n:]
        self.l = self.l[n:]
        self.b = self.b[n:]
        self.r = self.r[n:]
        self.vx = self.vx[n:]
        self.vy = self.vy[n:]
        self.vz = self.vz[n:]
        self.mass = self.mass[n:]
        self.vlos = self.vlos[n:]
        self.msol = self.msol[n:]
        self.ra = self.ra[n:]
        self.dec = self.dec[n:]
        self.d = self.d[n:]
        self.vgsr = self.vgsr[n:]
        self.vtan = self.vtan[n:]
        self.lz = self.lz[n:]
        self.lperp = self.lperp[n:]
        self.energy = self.energy[n:]
        self.pmra = self.pmra[n:]
        self.pmdec = self.pmdec[n:]
        self.rad = self.rad[n:]
        self.updateArrayDict()

    def cutLastN(self, n):
        self.id = self.id[:len(self.rad) - n]
        self.x = self.x[:len(self.rad) - n]
        self.y = self.y[:len(self.rad) - n]
        self.z = self.z[:len(self.rad) - n]
        self.l = self.l[:len(self.rad) - n]
        self.b = self.b[:len(self.rad) - n]
        self.r = self.r[:len(self.rad) - n]
        self.vx = self.vx[:len(self.rad) - n]
        self.vy = self.vy[:len(self.rad) - n]
        self.vz = self.vz[:len(self.rad) - n]
        self.mass = self.mass[:len(self.rad) - n]
        self.vlos = self.vlos[:len(self.rad) - n]
        self.msol = self.msol[:len(self.rad) - n]
        self.ra = self.ra[:len(self.rad) - n]
        self.dec = self.dec[:len(self.rad) - n]
        self.d = self.d[:len(self.rad) - n]
        self.vgsr = self.vgsr[:len(self.rad) - n]
        self.vtan = self.vtan[:len(self.rad) - n]
        self.lz = self.lz[:len(self.rad) - n]
        self.lperp = self.lperp[:len(self.rad) - n]
        self.energy = self.energy[:len(self.rad) - n]
        self.pmra = self.pmra[:len(self.rad) - n]
        self.pmdec = self.pmdec[:len(self.rad) - n]
        self.rad = self.rad[:len(self.rad) - n]
        self.updateArrayDict()

    def initial_energy(self):
        self.x = self.x - self.centerOfMass[0]
        self.y = self.y - self.centerOfMass[1]
        self.z = self.z - self.centerOfMass[2]

        self.vx = self.vx - self.centerOfMomentum[0]
        self.vy = self.vy - self.centerOfMomentum[1]
        self.vz = self.vz - self.centerOfMomentum[2]

        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5

        PE = galpy.potential.evaluatePotentials(plummer_pot, (self.x**2 + self.y**2)**0.5 * u.kpc, self.z*u.kpc, ro=10., vo=20.)
        KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)
        self.energy = PE + KE

    #data.plot(d1='var1', d2='var2'): data, str, str -> plot
    #takes in the 2 coordinates of the data you want to plot and plots them in a 2d scatter plot
    def plot(self, d1='r', d2='z', overplot=False, s=5.0, color='k', marker='o', **kwargs):
        #d1: the x-axis variable
        #d2: the y-axis variable
        #these can be x, y, z, vlos, l, b, etc. any attribute of the data class
        array_dict = self.array_dict
        x_array = array_dict[d1]
        y_array = array_dict[d2]
        plt.scatter(x_array, y_array, s=s, c=color, marker=marker, **kwargs)

        #plt.xlim([-60000, 10000])
        #plt.ylim([-1000, 1000])
        if not overplot:
            plt.xlabel(d1)
            plt.ylabel(d2)
            plt.show()

    #sticks a big fat red dot wherever the specific star is, given an id
    def trace_particle(self, id, d1='r', d2='z', overplot=False, s=50.0, color='r', marker='o', vel=False, **kwargs):
        #right now, id is the index of the star in the data structure
        #TODO: allow id matching to the id array in data structure as well
        array_dict = self.array_dict
        x_array = array_dict[d1][id]
        y_array = array_dict[d2][id]
        plt.scatter(x_array, y_array, s=s, c=color, marker=marker, **kwargs)
        if vel:
            vx = array_dict['vx'][id]
            vy = array_dict['vz'][id]
            #TODO: should alter to properly allow other velocities than what I'm hard coding
            #TODO: allow to change scaling of arrow length
            plt.arrow(x_array, y_array, vx/50, vy/50, color=color, head_width=1, **kwargs)

        if not overplot:
            plt.xlabel(d1)
            plt.ylabel(d2)
            plt.show()

    #data.hist(d='r'): data, str -> histogram plot
    #takes in the coordinate of the data you want in your histogram and then produces the relevant plot
    def hist(self, d='r', overplot=False, hist_range=None, hist_bins=10, *args, **kwargs):
        #d: the variable being binned
        #again, can be any attribute of the data class
        array_dict = self.array_dict
        x_array = array_dict[d]

        h = plt.hist(x_array, range=hist_range, bins=hist_bins, *args, **kwargs)
        if not overplot:
            plt.xlabel(d)
            plt.show()
        return h

    def hist2d(self, d1='r', d2='lz', bins=10, range=None, density=False, weights=None, cmin=None, cmax=None, data=None, overplot=False, **kwargs):
        array_dict = self.array_dict
        x = array_dict[d1]
        y = array_dict[d2]
        h = plt.hist2d(x, y, bins=bins, range=range, weights=weights, cmin=cmin, cmax=cmax, data=data, **kwargs)
        if not overplot:
            plt.xlabel(d1)
            plt.ylabel(d2)
            plt.show()
        return h

    def makeCSVFile(self, f_name):
        #f: the filename for the output csv file
        array_dict = self.array_dict
        array_order = []
        f = open(f_name, 'w')

        #make the header
        header = ''
        for key in array_dict:
            header += (key + ',')
        header += '\n'
        f.write(header)

        #iterate through the data and print each line
        i = 0
        while i < len(self.id):
            line = ''
            for key in array_dict:
                line += (str(array_dict[key][i]) + ',')
            line += '\n'
            f.write(line)
            i += 1

        print('Data output to ' + f_name)

    def updateArrayDict(self):
        self.array_dict = {'id':self.id, 'x':self.x, 'y':self.y, 'z':self.z, 'l':self.l, 'b':self.b, 'r':self.r, 'vx':self.vx, \
                            'vy':self.vy, 'vz':self.vz, 'mass':self.mass, 'vtan':self.vtan, 'vlos':self.vlos, \
                            'msol':self.msol, 'ra':self.ra, 'dec':self.dec, 'd':self.d, 'vgsr':self.vgsr, \
                            'energy':self.energy, 'ltot':self.ltot, 'lz':self.lz, 'lperp':self.lperp, 'pmra':self.pmra, \
                            'pmdec':self.pmdec, 'rad':self.rad}

    def append_data(self, append_data):
        self.id = np.append(self.id, append_data.id)
        self.x = np.append(self.x, append_data.x)
        self.y = np.append(self.y, append_data.y)
        self.z = np.append(self.z, append_data.z)
        self.l = np.append(self.l, append_data.l)
        self.b = np.append(self.b, append_data.b)
        self.r = np.append(self.r, append_data.r)
        self.vx = np.append(self.vx, append_data.vx)
        self.vy = np.append(self.vy, append_data.vy)
        self.vz = np.append(self.vz, append_data.vz)
        self.mass = np.append(self.mass, append_data.mass)
        self.vlos = np.append(self.vlos, append_data.vlos)
        self.msol = np.append(self.msol, append_data.msol)
        self.ra = np.append(self.ra, append_data.ra)
        self.dec = np.append(self.dec, append_data.dec)
        self.d = np.append(self.d, append_data.d)
        self.vgsr = np.append(self.vgsr, append_data.vgsr)
        self.vtan = np.append(self.vtan, append_data.vtan)
        self.lz = np.append(self.lz, append_data.lz)
        self.lperp = np.append(self.lperp, append_data.lperp)
        self.energy = np.append(self.energy, append_data.energy)
        self.pmra = np.append(self.pmra, append_data.pmra)
        self.pmdec = np.append(self.pmdec, append_data.pmdec)
        self.rad = np.append(self.rad, append_data.rad)
        self.updateArrayDict()

    def append_point(self, append_data, i):
        self.id = np.append(self.id, append_data.id[i])
        self.x = np.append(self.x, append_data.x[i])
        self.y = np.append(self.y, append_data.y[i])
        self.z = np.append(self.z, append_data.z[i])
        self.l = np.append(self.l, append_data.l[i])
        self.b = np.append(self.b, append_data.b[i])
        self.r = np.append(self.r, append_data.r[i])
        self.vx = np.append(self.vx, append_data.vx[i])
        self.vy = np.append(self.vy, append_data.vy[i])
        self.vz = np.append(self.vz, append_data.vz[i])
        self.mass = np.append(self.mass, append_data.mass[i])
        self.vlos = np.append(self.vlos, append_data.vlos[i])
        self.msol = np.append(self.msol, append_data.msol[i])
        self.ra = np.append(self.ra, append_data.ra[i])
        self.dec = np.append(self.dec, append_data.dec[i])
        self.d = np.append(self.d, append_data.d[i])
        self.vgsr = np.append(self.vgsr, append_data.vgsr[i])
        self.vtan = np.append(self.vtan, append_data.vtan[i])
        self.lz = np.append(self.lz, append_data.lz[i])
        self.lperp = np.append(self.lperp, append_data.lperp[i])
        self.energy = np.append(self.energy, append_data.energy[i])
        self.pmra = np.append(self.pmra, append_data.pmra[i])
        self.pmdec = np.append(self.pmdec, append_data.pmdec[i])
        self.rad = np.append(self.rad, append_data.rad[i])
        self.updateArrayDict()

#read_output(f): filename -> data class
#reads a milky way at home output file and turns it into a data class
def readOutput(f, init_energy=False, subsample=1.0, pot_offset=0):
    #f: the filename, formatted ('~/milkywayathome_client/nbody/...')
    #subsample: the percentage [0.0, 1.0] of the sample to use
    f = open(f)
    #remove the header: this is info we don't need
    comass = []
    comom = []
    if init_energy:
        for i in range(0, 4):
            f.readline()
        #read in data for initial_energy work
        line = f.readline()
        line = line.split(',')
        line[0] = line[0].strip('centerOfMass = ')
        line[3] = line[3].strip('centerOfMomentum = ')
        comass = [float(line[0]), float(line[1]), float(line[2])]
        comom = [float(line[3]), float(line[4]), float(line[5])]

        f.readline()
    else:
        for i in range(0, 5):
            f.readline()

    #store the data here temporarily
    #indexed this way to avoid the 'ignore' column
    array_dict = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[], 12:[]}

    #place all the data from the file into the dictionary
    for line in f:
        m = random.random()
        if m <= subsample:
            line = line.strip().split(',')
            i = 1
            while i < len(line):
                array_dict[i].append(float(line[i]))
                i += 1

    #return the data class using the array dictionary we built
    return Data(array_dict[1], array_dict[2], array_dict[3], array_dict[4], array_dict[5], array_dict[6], array_dict[7], array_dict[8], array_dict[9], array_dict[10], array_dict[11], array_dict[12], centerOfMass=comass, centerOfMomentum=comom, pot_offset=pot_offset)

#subset(data): data_object -> data_object
#takes in a data object and outputs a cut data object. Can cut within some radius or a rectangle cut. Can specify the axes, or if there is only 1 axis.
def subset(data, ax1='ra', ax2=None, rect=False, radius=None, center=None, corner1=None, corner2=None):
    #data: the data object being cut
    #ax1: the x-axis for the subset procedure
    #ax2: the y-axis for the procedure. If not specified, function assumes 1-dimensional cut on ax1
    #rect: If True, use a rectangular cut (corner1 & corner2) instead of radius to cut
    #radius: the radius for a radius cut
    #center: the center of the radius cut tuple(x,y) if 2d or a number if 1d
    #corner1: Bottom left corner for rectangular cut tuple(x, y)
    #corner2: Top right corner for rectangular cut tuple(x, y)

    array_dict = data.array_dict
    data_out = Data()

    if rect:
        if (not ax2) or (not corner1) or (not corner2):
            print('Must provide ax2, corner1 and corner2 for a rectangular cut')
            return None
        else: #do rectangular cut                    data_out.id = np.append(data_out.id, Data.id[i])
            i = 0
            while i < len(data.ra):
                if corner1[0] < array_dict[ax1][i] < corner2[0] and corner1[1] < array_dict[ax2][i] < corner2[1]:
                    data_out.append_point(data, i)
                i+=1
            return data_out

    else: #do radius cut
        if (not radius) or center == None: #radius and/or center weren't provided
            print('Must provide radius and center for a radial cut')
            return None
        else:
            if ax2: #do 2d cut
                i = 0
                while i < len(data.ra):
                    if radius >= ((array_dict[ax1][i] - center[0])**2 + (array_dict[ax2][i] - center[1])**2)**0.5:
                        data_out.append_point(data, i)
                    i+=1
                return data_out
            else: #do 1d cut
                i = 0
                while i < len(data.ra):
                    if radius >= abs(array_dict[ax1][i] - center):
                        data_out.append_point(data, i)
                    i+=1
                return data_out
