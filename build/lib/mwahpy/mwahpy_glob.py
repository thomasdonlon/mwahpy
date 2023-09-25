'''
This file is a place to store constants and other information
that are widely used throughout mwahpy.
'''

#===============================================================================
# IMPORTS
#===============================================================================

import astropy.units as u

#===============================================================================
# CLASSES
#===============================================================================

#currently just used for the position of the Sun, but could be more useful in other places
class PhaseSpace():

    def __init__(self, ps):
        self.x = ps[0]
        self.y = ps[1]
        self.z = ps[2]
        self.vx = ps[3]
        self.vy = ps[4]
        self.vz = ps[5]

    #this way you can change the phase space used for the Sun within the entire package
    def update(self, ps):
        self.x = ps[0]
        self.y = ps[1]
        self.z = ps[2]
        self.vx = ps[3]
        self.vy = ps[4]
        self.vz = ps[5]

#===============================================================================
# CONSTANTS
#===============================================================================

G = 6.674e-11*u.m**3/(u.kg*u.s**2)

struct_to_sol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

kms_to_kpcgyr = 1.023 #1 km/s is 1.023 kpc/Gyrb

kpcgyr_to_kms = 0.978 #1 kpc/Gyr is 0.978 km/s

#define transformations from heliocentric -> galactocentric frame
#by defining heliocentric frame here

#right-handed coordinate frame(!)
#Hogg et al. (2005)
solar_ps = PhaseSpace([-8., 0., 0., 10.1, 224.0, 6.7])

#===============================================================================
# FUNCTIONS
#===============================================================================

#prints out a progress bar in the terminal
#"borrowed" from https://stackoverflow.com/questions/6169217/replace-console-output-in-python
def progress_bar(value, endvalue, bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    print('\r[{0}] {1}%'.format(arrow + spaces, int(round(percent * 100))), end='')

#get length of a file (in lines), given the filename
def file_len(f):
    with open(f) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
