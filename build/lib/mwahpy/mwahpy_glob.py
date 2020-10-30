'''
This file is a place to store constants and other information
that is widely used throughout mwahpy.
'''

#===============================================================================
# IMPORTS
#===============================================================================

import sys
import astropy.units as u

#===============================================================================
# CONSTANTS
#===============================================================================

G = 6.674e-11*u.m**3/(u.kg*u.s**2)

structToSol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

kmsToKpcgyr = 1.023 #1 km/s is 1.023 kpc/Gyr
kpcgyrToKms = 0.978 #1 kpc/Gyr is 0.978 km/s

#===============================================================================
# FUNCTIONS
#===============================================================================

#prints out a progress bar in the terminal
#"borrowed" from https://stackoverflow.com/questions/6169217/replace-console-output-in-python
def progressBar(value, endvalue, bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

#get length of a file (in lines), given the filename
def fileLen(f):
    with open(f) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
