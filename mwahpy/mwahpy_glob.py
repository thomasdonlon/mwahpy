'''
This file is a place to store constants and other information
that is widely used throughout mwahpy.
'''

#TODO: Get more accurate solar velocity and position of Sun from GC
#TODO: Make more tests and put them all in a single test.py file
#TODO: Ensure all imports are actually used (and necessary) for each file

#===============================================================================
# IMPORTS
#===============================================================================

import sys
import astropy.units as u

#===============================================================================
# CONSTANTS
#===============================================================================

G = 6.674e-11*u.m**3/(u.kg*u.s**2)

struct_to_sol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

kms_to_kpcgyr = 1.023 #1 km/s is 1.023 kpc/Gyr
kpcgyr_to_kms = 0.978 #1 kpc/Gyr is 0.978 km/s

#===============================================================================
# FUNCTIONS
#===============================================================================

#prints out a progress bar in the terminal
#"borrowed" from https://stackoverflow.com/questions/6169217/replace-console-output-in-python
def progress_bar(value, endvalue, bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

#get length of a file (in lines), given the filename
def file_len(f):
    with open(f) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
