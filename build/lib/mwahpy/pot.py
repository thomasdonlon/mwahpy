'''
This file is a place to store constants and other information
that is widely used throughout mwahpy.
'''

#===============================================================================
# IMPORTS
#===============================================================================

import astropy
import astropy.units as u
import matplotlib.pyplot as plt

import galpy
from galpy.potential import HernquistPotential
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import MiyamotoNagaiPotential

from .mwahpy_glob import G

#===============================================================================
# CONSTANTS
#===============================================================================

m_bulge = 3.4e10*u.solMass #solar masses
m_disk = 1.0e11*u.solMass
v_halo = 74.61*u.km/u.s #km/s

pot_bulge = HernquistPotential(amp=2*m_bulge, a=0.7*u.kpc, ro=8., vo=220.)
pot_disk = MiyamotoNagaiPotential(amp=G*m_disk, a=6.5*u.kpc, b=0.26*u.kpc, ro=8., vo=220.)
pot_halo = LogarithmicHaloPotential(amp=2*v_halo**2, q=1., core=12.0*u.kpc, ro=8., vo=220.)

#this potential is from Newberg et al. 2010, Orphan Stream Model 5. It's basically Law 2005
mwahpy_default_pot = [pot_bulge, pot_disk, pot_halo]

energy_offset = -60000 #adjusts the energy to be consistent with Donlon et al. 2019

#===============================================================================
# FUNCTIONS
#===============================================================================

#produces a rotation curve for the given potential
def plot_potential(potential, Rrange=[0.01,10.]):
    fig = plt.figure(figsize=(12,8))
    potential.plotRotcurve(Rrange=Rrange, overplot=True)
    plt.show()
