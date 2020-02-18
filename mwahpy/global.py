'''
This file is a place to store constants and other information
that is widely used throughout mwahpy.
'''

#===============================================================================
# IMPORTS
#===============================================================================

import astropy
import astropy.units as u

import galpy
from galpy.potential import HernquistPotential
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import MiyamotoNagaiPotential

#===============================================================================
# CONSTANTS
#===============================================================================

G = 6.67e-11*u.m**3/(u.kg*u.s**2)

struct_to_sol = 222288.47 #this many solar masses make up one structural nass unit (the output of mwah)

#-------------------------------------------------------------------------------
# GALPY DEFINITIONS
#-------------------------------------------------------------------------------

m_bulge = 3.4e10*u.solMass #solar masses
m_disk = 1.0e11*u.solMass
v_halo = 74.61*u.km/u.s #km/s

pot_bulge = HernquistPotential(amp=2*m_bulge, a=0.7*u.kpc, ro=8., vo=220.)
pot_disk = MiyamotoNagaiPotential(amp=G*m_disk, a=6.5*u.kpc, b=0.26*u.kpc, ro=8., vo=220.)
pot_halo = LogarithmicHaloPotential(amp=2*v_halo**2, q=1., core=12.0*u.kpc, ro=8., vo=220.)

#this potential is from Newberg et al. 2010, Orphan Stream Model 5. It's basically Law 2005
pot = [pot_bulge, pot_disk, pot_halo]

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
