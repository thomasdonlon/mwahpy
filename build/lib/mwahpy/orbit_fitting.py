'''
This is a self-contained orbit fitting routine.

This orbit fitter is unique compared to other common orbit fitters in that it
uses a galactocentric generalized plane coordinate system when fitting data
'''

#TODO: Allow fitting on b, ra, or dec instead of l

import numpy as np
import scipy as sc
import scipy.optimize as scopt
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

import galpy
from galpy.orbit import Orbit

from .flags import verbose
from .pot import mwahpy_default_pot

'''
================================================================================
 FLAGS
================================================================================
'''

#-------------------------------------------------------------------------------
#DO NOT TOUCH
#Any routine meant to be used by an end user will configure this from input data
vx_flag = 0
vy_flag = 0
vz_flag = 0
vgsr_flag = 0
#-------------------------------------------------------------------------------

'''
================================================================================
 PARAMETERS FOR OPTIMIZATION
================================================================================
'''
t_length = 0.5 #Gyr
resolution = 1000
ts = np.linspace(0, t_length, num=resolution)*u.Gyr

punishment = 1000 #multiplier for every degree that the lambda fitting function is off
#this can be tweaked based on how many degrees of freedom you have - fewer DOF
#means punishment can/should be smaller

'''
================================================================================
 HELPER FUNCTIONS FOR OPTIMIZATION
================================================================================
'''

class OrbitData():
    #all of these parameters can be np.arrays of floats
    def __init__(self, l, b, d, vx, vy, vz, vgsr, b_err, d_err, vx_err, vy_err, vz_err, vgsr_err):
        self.l = l
        self.b = b
        self.d = d #heliocentric distance
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.vgsr = vgsr

        self.b_err = b_err
        self.d_err = d_err
        self.vx_err = vx_err
        self.vy_err = vy_err
        self.vz_err = vz_err
        self.vgsr_err = vgsr_err

        self.x = self.d*np.cos(np.pi/180*self.l)*np.cos(np.pi/180*self.b) - 8
        self.y = self.d*np.sin(np.pi/180*self.l)*np.cos(np.pi/180*self.b)
        self.z = self.d*np.sin(np.pi/180*self.b)

    #add the icrs converted coordinates to this data instance
    def icrs(self):
        s = SkyCoord(self.l, self.b, frame='galactic', unit=(u.deg, u.deg))
        s = s.transform_to('icrs')
        self.ra = s.ra
        self.dec = s.dec

    #gets the orbit with the correct vgsr values
    #since you can't just multiply vgsr by -1 to get the correct value along an orbit
    #this needs to be run to compare the vgsr of reverse orbits
    def correctVgsr(self):
        self.vgsr = ((self.x + 8) * self.vx + self.y * self.vy + self.z * self.vz)/self.d

#this function just makes a few places in the code cleaner
def getOrbitDataFromOrbit(o, o_rev):
    #o: the orbit
    #o: the reversed orbit
    #both of these should be integrated prior to calling this function

    data_orbit = OrbitData(np.array(o.ll(ts)), np.array(o.bb(ts)), np.array(o.dist(ts)), np.array(o.vx(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o.vy(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o.vz(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o.vlos(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    data_orbit_rev = OrbitData(np.array(o_rev.ll(ts)), np.array(o_rev.bb(ts)), np.array(o_rev.dist(ts)), np.array(o_rev.vx(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o_rev.vy(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o_rev.vz(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    data_orbit_rev.correctVgsr()

    return data_orbit, data_orbit_rev

def getClosestIndex(val, Lam):
    L = np.abs(Lam - val)
    m = np.min(L)
    ind = np.where(L == m)

    return ind[0], m*punishment

#getPointList: np.array([]), np.array([]), int, int -> np.array([])
#given the full list of Lambdas, outputs the indices of the points within that list closest to our data's Lambdas
#(while keeping in mind that it may wrap from 0 to 360 degrees and vice versa)
def getPointList(vals, Lam):
    #vals: the Lambda values that you want to find the closest indices to
    #Lam: the array of Lambda values that you are searching through

    point_list = []
    Lam_list = []
    costs = 0

    for val in vals:

        #within that segment, find the index which produces the value closest to val
        point, c = getClosestIndex(val, Lam)
        costs += c

        #toss it in the list
        point_list.append(point)
        Lam_list.append(Lam[point])

    return point_list, costs

#getModelFromOrbit: data, orbit, vector, vector -> list(int) x3
#take in data, orbit, and plane info: Output model data corresponding to each data point
def getModelFromOrbit(data, o):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    #initialize the orbit we are fitting --
    #we flip it around so that we are fitting both the forwards and the backwards orbit
    ts = np.linspace(0, t_length, num=resolution)*u.Gyr
    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    #sign swap on vx because galpy is left-handed, and we are inputting data in a right-handed coordinate system
    data_orbit, data_orbit_rev = getOrbitDataFromOrbit(o, o_rev)

    #grab full lists so that we can select the closest points once we get a list
    Lam = np.append(np.flip(data_orbit_rev.l), data_orbit.l)

    #get the list of points closest to each data point in Lambda
    point_list, costs = getPointList(data.l, Lam)

    #grab the model points from the point list we grabbed
    Bet = np.append(np.flip(data_orbit_rev.b), data_orbit.b)
    B_model = np.array([Bet[p] for p in point_list]).flatten()

    D = np.append(np.flip(data_orbit_rev.d), data_orbit.d)
    D_model = np.array([D[p] for p in point_list]).flatten()

    if vx_flag:
        vx = np.append(np.flip(data_orbit_rev.vx), data_orbit.vx)
        vx_model = np.array([vx[p] for p in point_list]).flatten()
    else:
        vx_model = np.zeros(len(B_model))

    if vy_flag:
        vy = np.append(np.flip(data_orbit_rev.vy), data_orbit.vy)
        vy_model = np.array([vy[p] for p in point_list]).flatten()
    else:
        vy_model = np.zeros(len(B_model))

    if vz_flag:
        vz = np.append(np.flip(data_orbit_rev.vz), data_orbit.vz)
        vz_model = np.array([vz[p] for p in point_list]).flatten()
    else:
        vz_model = np.zeros(len(B_model))

    if vgsr_flag:
        vgsr = np.append(np.flip(data_orbit_rev.vgsr), data_orbit.vgsr)
        vgsr_model = np.array([vgsr[p] for p in point_list]).flatten()
    else:
        vgsr_model = np.zeros(len(B_model))

    return B_model, D_model, vx_model, vy_model, vz_model, vgsr_model, costs

#chi_squared: data, galpy.Orbit() --> float
#takes in the observed data and a test orbit and calculates the goodness-of-fit using a chi-squared method
def chiSquared(params, data=[]):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    o = Orbit(vxvv=[params[0], params[1], params[2], params[3], params[4]-220, params[5]], uvw=True, lb=True, ro=8., vo=220., zo=0.) #generate the orbit
    o.integrate(ts, mwahpy_default_pot) #integrate the orbit

    B_model, d_model, vx_model, vy_model, vz_model, vgsr_model, costs = getModelFromOrbit(data, o) #get model data from orbit

    #B_model sometimes has different length than data.b, no idea why
    #I think it might be a race condition
    #this keeps the script running and tells the optimizer that the parameters are bad
    if len(B_model) != len(data.b):
        return 1e10

    x2_B = sum(((B_model - data.b)/data.b_err)**2)
    x2_d = sum(((d_model - data.d)/data.d_err)**2)

    if vx_flag:
        x2_vx = sum(((vx_model - data.vx)/data.vx_err)**2)
    else:
        x2_vx = 0

    if vy_flag:
        x2_vy = sum(((vy_model - data.vy)/data.vy_err)**2)
    else:
        x2_vy = 0

    if vz_flag:
        x2_vz = sum(((vz_model - data.vz)/data.vz_err)**2)
    else:
        x2_vz = 0

    if vgsr_flag:
        x2_vgsr = sum(((vgsr_model - data.vgsr)/data.vgsr_err)**2)
    else:
        x2_vgsr = 0

    #get normalization factor
    N = len(data.l) #number of data points
    n = 5 #number of parameters
    eta = N - n - 1 #normalizing parameter
    if eta <= 0:
        eta = 1 #if you use fewer data points than needed to constrain the problem, then this will still work but it won't be normalized correctly

    x2 = (1/eta) * (x2_B + x2_d + x2_vx + x2_vy + x2_vz + x2_vgsr) + costs #Willett et al. 2009, give or take

    #there's a weird edge case where occasionally x2 is a short array of floats
    #this bit prevents scipy from throwing an error
    #no idea what causes that
    if type(x2) == type(np.array([])):
        x2 = x2[0]

    if flags.verbose:
        print('X^2: ' + str(x2))

    return x2

#optimize: data -> [float, float, float, float, float], (float, float, float), (float, float, float)
#takes in data, then fits a Great Circle to that data and minimizes the chi_squared to fit an orbit to the data
def optimize(data_opt, max_it, bounds, **kwargs):

    '''
    ============================================================================
    DIFFERENTIAL EVOLUTION CONSTANTS
    ============================================================================
    '''
    #DO NOT TOUCH
    pop_size = 50 #10 times number of parameters
    diff_scaling_factor = 0.8
    crossover_rate = 0.9
    '''
    ============================================================================
    '''

    params = scopt.differential_evolution(chiSquared, bounds, args=(data_opt,), strategy='rand1bin', maxiter=max_it, popsize=pop_size, mutation=diff_scaling_factor, recombination=crossover_rate, workers=-1, disp=not(flags.verbose), **kwargs).x
    #'''

    x2 = chiSquared(params, data_opt)

    return params, x2

'''
================================================================================
 FUNCTIONS
================================================================================
'''

def fit_orbit(l, b, b_err, d, d_err, vx=None, vy=None, vz=None, vgsr=None, \
              vx_err=None, vy_err=None, vz_err=None, vgsr_err=None, max_it=20, \
              bounds=[(0, 360), (-90, 90), (0, 100), (-500, 500), (-500, 500), (-500, 500)], \
              t_len=None, **kwargs):

    #construct data
    #set proper flags based on input data
    if type(vx) == type(np.array([])):
        global vx_flag
        vx_flag = 1
    if type(vy) == type(np.array([])):
        global vy_flag
        vy_flag = 1
    if type(vz) == type(np.array([])):
        global vz_flag
        vz_flag = 1
    if type(vgsr) == type(np.array([])):
        global vgsr_flag
        vgsr_flag = 1

    #update t_length if necessary
    if t_len != None:
        global t_length
        t_length = t_len
        global ts
        ts = np.linspace(0, t_length, num=resolution)*u.Gyr

    if flags.verbose:
        print('===================================')
        print('Optimizing:')
        print('===================================')

    data_opt = OrbitData(l, b, d, vx, vy, vz, vgsr, b_err, d_err, vx_err, vy_err, vz_err, vgsr_err)

    #optimization
    params, x2 = optimize(data_opt, max_it, bounds, **kwargs)

    print('===================================')
    print('Params: l, b, d, vx, vy, vz')
    print(params)
    print()
    print('Chi Squared:')
    print(x2)
    print('===================================')

    return params, x2

'''
================================================================================
 PLOTTING
================================================================================
'''

#TODO: Implement unwrap in a way that actually makes sense
#splits every array in the list of arrays, a, every time that the position wraps
#from 0 to 360 or vice versa in a[0]
#therefore, a[0] should be the parameter you are trying to unwrap (longitude)
#returns a list of lists of the unwrapped arrays
def unwrap(a, threshold=10):
    #t: difference in position needed to trigger a split
    split = np.nonzero(np.abs(a[0][:-1] - a[0][1:]) > threshold)[0] + 1

    out = []
    for arr in a:
        if len(split) > 0:
            out.append(np.split(arr, split))
        else:
            out.append(np.array([arr])) #didn't find a place to split on

    return out

#TODO: Expand this and plotOrbiticrs to allow other velocities
#possibly make them the same function with a switch
#TODO: Split values so wrapping lines don't happen

def plotOrbitgal(l, b, d, params, vgsr=None):
    o = Orbit(vxvv=[params[0], params[1], params[2], params[3], params[4] - 220, params[5]], uvw=True, lb=True, ro=8., vo=220.) #generate the orbit
    o.integrate(ts, mwahpy_default_pot) #integrate the orbit

    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    #sign swap on vx because galpy is left-handed, and we are inputting data in a right-handed coordinate system
    data_orbit, data_orbit_rev = getOrbitDataFromOrbit(o, o_rev)

    fig = plt.figure(figsize=(24, 6))
    nplots = 2
    if type(vgsr) == type(np.array([])):
        nplots += 1
    ax1 = fig.add_subplot(1, nplots, 1)
    ax2 = fig.add_subplot(1, nplots, 2)
    if type(vgsr) == type(np.array([])):
        ax3 = fig.add_subplot(1, nplots, 3)

    ax1.plot(data_orbit.l, data_orbit.b, c='b')
    ax1.plot(data_orbit_rev.l, data_orbit_rev.b, c='r')
    ax1.scatter(l, b, c='k')

    ax1.set_xlim(0, 360)
    ax1.set_ylim(-90, 90)
    ax1.set_xlabel('l')
    ax1.set_ylabel('b')

    ax2.plot(data_orbit.l, data_orbit.d, c='b')
    ax2.plot(data_orbit_rev.l, data_orbit_rev.d, c='r')
    ax2.scatter(l, d, c='k')

    ax2.set_xlim(0, 360)
    ax2.set_xlabel('l')
    ax2.set_ylabel('d (helio)')

    if type(vgsr) == type(np.array([])):
        ax3.plot(data_orbit.l, data_orbit.vgsr, c='b')
        ax3.plot(data_orbit_rev.l, data_orbit_rev.vgsr, c='r')
        ax3.scatter(l, vgsr, c='k')

        ax3.set_xlim(0, 360)
        ax3.set_xlabel('l')
        ax3.set_ylabel('vgsr (km/s)')

    plt.show()

def plotOrbiticrs(l, b, d, params, vgsr=None):

    s = SkyCoord(l, b, frame='galactic', unit=(u.deg, u.deg))
    s = s.transform_to('icrs')
    ra = s.ra
    dec = s.dec

    o = Orbit(vxvv=[params[0], params[1], params[2], params[3], params[4] - 220, params[5]], uvw=True, lb=True, ro=8., vo=220.) #generate the orbit
    o.integrate(ts, mwahpy_default_pot) #integrate the orbit

    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    #sign swap on vx because galpy is left-handed, and we are inputting data in a right-handed coordinate system
    data_orbit, data_orbit_rev = getOrbitDataFromOrbit(o, o_rev)

    fig = plt.figure(figsize=(24, 6))
    nplots=2
    if type(vgsr) == type(np.array([])):
        nplots += 1
    ax1 = fig.add_subplot(1,nplots,1)
    ax2 = fig.add_subplot(1,nplots,2)
    if type(vgsr) == type(np.array([])):
        ax3 = fig.add_subplot(1,nplots,3)

    data_orbit.icrs()
    data_orbit_rev.icrs()
    #TODO: This will break if vgsr isn't used
    #TODO: Unwrap should really be a orbit data method
    o_unwrapped = unwrap([data_orbit.ra, data_orbit.dec, data_orbit.d, data_orbit.vgsr])
    o_rev_unwrapped = unwrap([data_orbit_rev.ra, data_orbit_rev.dec, data_orbit_rev.d, data_orbit_rev.vgsr])

    for o_ra, o_dec in zip(o_unwrapped[0], o_unwrapped[1]):
        ax1.plot(o_ra, o_dec, c='b')
    for o_ra, o_dec in zip(o_rev_unwrapped[0], o_rev_unwrapped[1]):
        ax1.plot(o_ra, o_dec, c='r')
    ax1.scatter(ra, dec, c='k')

    ax1.set_xlim(360, 0)
    ax1.set_ylim(-90, 90)
    ax1.set_xlabel('ra')
    ax1.set_ylabel('dec')

    for o_ra, o_d in zip(o_unwrapped[0], o_unwrapped[2]):
        ax2.plot(o_ra, o_d, c='b')
    for o_ra, o_d in zip(o_rev_unwrapped[0], o_rev_unwrapped[2]):
        ax2.plot(o_ra, o_d, c='r')
    ax2.scatter(ra, d, c='k')

    ax2.set_xlim(360, 0)
    ax2.set_xlabel('ra')
    ax2.set_ylabel('d (helio)')

    if type(vgsr) == type(np.array([])):
        for o_ra, o_vgsr in zip(o_unwrapped[0], o_unwrapped[-1]):
            ax1.plot(o_ra, o_vgsr, c='b')
        for o_ra, o_vgsr in zip(o_rev_unwrapped[0], o_rev_unwrapped[-1]):
            ax1.plot(o_ra, o_vgsr, c='r')
        ax3.scatter(ra, d, c='k')

        ax3.set_xlim(360, 0)
        ax3.set_xlabel('ra')
        ax3.set_ylabel('vgsr (km/s)')

    plt.show()

'''
================================================================================
 TESTING
================================================================================
'''
def test():
    l = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180])
    b = np.array([0, 10, 20, 10, 0, -10, -20, -10, 0, 10])
    b_err = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    d = np.array([20, 18, 15, 12, 10, 12, 15, 18, 20, 23])
    d_err = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    params, x2 = fit_orbit(l, b, b_err, d, d_err)
    plotOrbitgal(l, b, d, params)

'''
================================================================================
 RUNTIME
================================================================================
'''

if __name__ == "__main__":
    test()
