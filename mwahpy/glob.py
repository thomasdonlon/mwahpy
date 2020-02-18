'''
This file is a place to store constants and other information
that is widely used throughout mwahpy.
'''

#===============================================================================
# IMPORTS
#===============================================================================

import astropy.units as u

#===============================================================================
# CONSTANTS
#===============================================================================

G = 6.67e-11*u.m**3/(u.kg*u.s**2)

struct_to_sol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

#===============================================================================
# FUNCTIONS
#===============================================================================

#prints out a progress bar in the terminal
#"borrowed" from https://stackoverflow.com/questions/6169217/replace-console-output-in-python
def progressBar(value, endvalue, bar_length=20):

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

#get length of a file (in lines), given the filename
def fileLen(f):
    with open(f) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
