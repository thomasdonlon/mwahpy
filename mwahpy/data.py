'''
The contents of this file are focused on the Data class, which is used for storage of
imported data from N-body output files.
'''

#TODO: split functions should fix id values? Maybe not, depends on what the behavior is supposed to do

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

import mwahpy_glob
import flags

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

    def __init__(self, id_val=[], x=[], y=[], z=[], l=[], b=[], r=[], vx=[], vy=[], vz=[], mass=[], vlos=[], centerOfMass=[], centerOfMomentum=[], pot_offset=0, *args, **kwargs):
        #all this is typically specified by the readOutput function in output_handler
        #readOutput is the preferred way to input data to this data structure
        #but if you're really feeling adventurous you can always do it yourself

        #this data structure allows access to a dictionary like an attribute, i.e.
        #   d = AttrDict({'a':1, 'b':2})
        #   d['a'] == d.a //returns True
        #it seems like magic but it works by accessing the way python handles attributes
        #within each class. There are pros and cons, the most notable thing is
        #that it creates a memory leak in python <3.2.3
        ad = AttrDict() #this is a private helper class that allows for the
                        #desired indexing behavior
        self.__dict__ = ad.__dict__

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
        self.centerOfMass = centerOfMass
        self.centerOfMomentum = centerOfMomentum

        self.msol = self.mass * mwahpy_glob.struct_to_sol

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

        #-----------------------------------------------------------------------
        if flags.calcEnergy:
            #get the energy info
            #calculating the energy of every particle can generate some overhead,
            #so I've quarantined it with a flag.

            #in a logarithmic halo, the magnitude of the potential doesn;t impact the result,
            #just the difference in potentials. So, you can specify a potential offset
            #to keep bound objects' total energy negative.
            PE = galpy.potential.evaluatePotentials(mg.pot, (self.x**2 + self.y**2)**0.5 * u.kpc, self.z*u.kpc, ro=8., vo=220.) - pot_offset
            KE = 0.5*(self.vx**2 + self.vy**2 + self.vz**2)

            #set attributes
            self.PE = PE
            self.KE = KE
            self.energy = PE + KE
        #-----------------------------------------------------------------------

    def __getitem__(self, i):
        return self.__dict__[i]

    def __setitem__(self, i, val):
        self.__dict__[i] = val

    def __len__(self):
        return len(self.id)

    #creates a deep copy of the Data object
    def copy(self):
        out = Data()
        for key in self.__dict__.keys():
            out[str(key)] = self[str(key)].copy()

        return out

    #cuts the first n entries from every attribute in the data structure
    def cutFirstN(self, n):
        for key in self.__dict__.keys():
            self[str(key)] = self[str(key)][n:]

    #cuts the last n entries from every attribute in the data structure
    def cutLastN(self, n):
        for key in self.__dict__.keys():
            self[str(key)] = self[str(key)][:len(self[str(key)]) - n]

    #splits the Data into two new Data structures,
    #the first has the data points up to entry n and
    #the second has the data points after entry n
    def split(self, n):
        Data1 = self.copy()
        Data2 = self.copy()

        Data1.cutLastN(len(Data1) - n)
        Data2.cutFirstN(n)

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

        return outlist

    #append a data object onto this one
    #TODO: fix breakage if this object has energy but the appended object doesn't
    def appendData(self, d):
        #d: the data to append to this object
        for key in self.keys():
            self[str(key)] = np.append(self[str(key)], d[str(key)])

    #append the nth item of another data object onto this one
    #TODO: fix breakage if this object has energy but the appended object doesn't
    def appendPoint(self, d, n=0, id=None):
        #d: the data object with the item being appended
        #n: the index of the item in d to be appended
        #id: if not None, uses finds the first item with matching id and appends that
        if id:
            n = np.where(d['id'] == id)[0]
        for key in self.keys():
            self[str(key)] = np.append(self[str(key)], d[str(key)][n])

    #---------------------------------------------------------------------------
    # OUTPUT
    #---------------------------------------------------------------------------

    #prints out a '.csv' file of the Data class structure
    def makeCSV(self, f_name):
        #f: the path for the output csv file

        f = open(f_name, 'w')

        #make the header
        if flags.verbose:
            print('Writing header...')
        header = ''
        for key in self.array_dict:
            header += (key + ',')
        header += '\n'
        f.write(header)

        #iterate through the data and print each line
        i = 0
        if flags.verbose:
            print('Printing data...')
        while i < len(self.id): #i avoided zip() here because it's messy to zip
        #                        like a billion different arrays, although zip()
        #                        is "better programming"
            if flags.progressBars:
                mg.progressBar(i, len(self.id))
            line = ''
            for key in array_dict:
                line += (str(array_dict[key][i]) + ',')
            line += '\n'
            f.write(line)
            i += 1

        print('Data output to ' + f_name)

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
            while i < len(data[x]):
                if flags.progressBars:
                    mg.progressBar(i, len(data[x]))
                #check if within bounds
                if xbounds[0] < data[x][i] < xbounds[1] and ybounds[0] < data[y][i] < ybounds[1]:
                    data_out.append_point(data, i)
                i+=1
            if flags.verbose:
                print(str(len(data_out[x])) + ' objects found in bounds')
            return data_out

    else: #Not rectangular cut
        if not(radius) or not(center): #radius and/or center weren't provided
            raise Exception('Must provide <center> and <radius> for a radial cut')
        else:
            if ax2: #do 2D cut
                if not(type(center) is tuple): #make sure center is a tuple
                    raise Exception('Must provide tuple for <center> if <y> is provided')
                else:
                    i = 0
                    while i < len(data[x]):
                        if flags.progressBars:
                            mg.progressBar(i, len(data[x]))
                        if radius >= ((array_dict[x][i] - center[0])**2 + (array_dict[y][i] - center[1])**2)**0.5:
                            data_out.append_point(data, i)
                        i+=1
                    if flags.verbose:
                        print(str(len(data_out[x])) + ' objects found in bounds')
                    return data_out
            else: #do 1D cut
                i = 0
                while i < len(data[x]):
                    if flags.progressBars:
                        mg.progressBar(i, len(data[x]))
                    if radius >= abs(array_dict[x][i] - center):
                        data_out.append_point(data, i)
                    i+=1
                if flags.verbose:
                    print(str(len(data_out[x])) + ' objects found in bounds')
                return data_out
