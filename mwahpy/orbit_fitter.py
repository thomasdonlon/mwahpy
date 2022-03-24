'''
This is a self-contained orbit fitting routine.
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
from .coords import cart_to_lambet

'''
================================================================================
 FLAGS
================================================================================
'''

#-------------------------------------------------------------------------------
#DO NOT TOUCH
#Any routine meant to be used by an end user will configure this from input data
vx_flag = False
vy_flag = False
vz_flag = False
vgsr_flag = False
#-------------------------------------------------------------------------------

'''
================================================================================
 PARAMETERS FOR OPTIMIZATION
================================================================================
'''
t_length = 0.5 #Gyr
resolution = 1000
ts = np.linspace(0, t_length, num=resolution)*u.Gyr

deg_sep_pun = 1000 #multiplier for every degree that the lambda fitting function is off
#this can be tweaked based on how many degrees of freedom you have - fewer DOF
#means punishment can/should be smaller

'''
================================================================================
 CLASSES
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

        self.x = self.d*np.cos(self.l*np.pi/180)*np.cos(self.b*np.pi/180) - 8
        self.y = self.d*np.sin(self.l*np.pi/180)*np.cos(self.b*np.pi/180)
        self.z = self.d*np.sin(self.b*np.pi/180)

        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5

    #add the icrs converted coordinates to this data instance
    def icrs(self):
        s = SkyCoord(self.l, self.b, frame='galactic', unit=(u.deg, u.deg))
        s = s.transform_to('icrs')
        self.ra = s.ra
        self.dec = s.dec

    #compute the lat/long for an arbitrary plane, given by its normal vector and rough origin point. \
    #used in orbit_fitter_gc
    def calc_lam_bet(self, normal, point):
        #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
        #point: parameter for the axis generation of the Great Circle coordinates
        self.L, self.B = cart_to_lambet(self.x, self.y, self.z, normal, point)
        self.B_err = self.b_err #should be error propogated (Newby et al. 2013) but is not yet

'''
================================================================================
 HELPER FUNCTIONS FOR OPTIMIZATION
================================================================================
'''

#do this in one place so it's uniform across the file and easily changeable
def make_orbit(params, int_ts=None):

    #negate the x velocity because galpy is in left-handed frame and we aren't barbarians
    o = Orbit(vxvv=[params[0]*u.deg, params[1]*u.deg, params[2]*u.kpc, -1*params[3]*u.km/u.s, params[4]*u.km/u.s, params[5]*u.km/u.s], uvw=True, lb=True, ro=8., vo=220., zo=0., solarmotion=[0, -220, 0])

    if int_ts: #use the global ts unless explicitly provided with a different ts
        o.integrate(int_ts, mwahpy_default_pot)
    else:
        o.integrate(ts, mwahpy_default_pot)

    return o

#-------------------------------------------------------------------------------

#this function just makes a few places in the code cleaner
def get_OrbitData_from_orbit(o, o_rev):
    #o: the orbit
    #o: the reversed orbit
    #both of these should be integrated prior to calling this function

    #don't think that the velocities need obs=[8., 0., etc] except the vlos (to make it vgsr)
    data_orbit = OrbitData(np.array(o.ll(ts)), np.array(o.bb(ts)), np.array(o.dist(ts)), np.array(o.vx(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o.vy(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o.vz(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o.vlos(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    data_orbit_rev = OrbitData(np.array(o_rev.ll(ts)), np.array(o_rev.bb(ts)), np.array(o_rev.dist(ts)), np.array(o_rev.vx(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o_rev.vy(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o_rev.vz(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o.vlos(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))

    return data_orbit, data_orbit_rev

#-------------------------------------------------------------------------------

#given the full list of Lambdas, outputs the indices of the points within that list closest to our data's Lambdas
#(while keeping in mind that it may wrap from 0 to 360 degrees and vice versa)
def get_point_list(vals, Lam):
    #vals: the Lambda values that you want to find the closest indices to
    #Lam: the array of Lambda values that you are searching through

    #---------------------------------------------------------------------------
    #grabs the index of the value in Lam with the closest value to val
    def get_closest_index(val, Lam):
        Lam = np.abs(Lam - val)
        m = np.min(Lam)

        try: #another race condition?
            ind = np.where(Lam == m)[0][0]
        except IndexError:
            ind = -1

        return ind, m*deg_sep_pun
    #---------------------------------------------------------------------------

    point_list = []
    Lam_list = []
    costs = 0

    for val in vals:

        #within that segment, find the index which produces the value closest to val
        point, c = get_closest_index(val, Lam)
        costs += c

        #toss it in the list
        point_list.append(point)
        Lam_list.append(Lam[point])

    return point_list, costs

#get_model_from_orbit: data, orbit, vector, vector -> list(int) x3
#take in data, orbit, and plane info: Output model data corresponding to each data point
def get_model_from_orbit(data, params):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    #initialize the orbit we are fitting
    o = make_orbit(params)

    #we flip it around so that we are fitting both the forwards and the backwards orbit
    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    data_orbit, data_orbit_rev = get_OrbitData_from_orbit(o, o_rev)

    #grab full lists so that we can select the closest points once we get a list
    Lam = np.append(np.flip(data_orbit_rev.l), data_orbit.l)

    #get the list of points closest to each data point in Lambda
    point_list, costs = get_point_list(data.l, Lam)

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
def chi_squared(params, data=[]):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    B_model, d_model, vx_model, vy_model, vz_model, vgsr_model, costs = get_model_from_orbit(data, params) #get model data from orbit

    #B_model sometimes has different length than data.b, no idea why
    #I think it might be a race condition
    #this keeps the script running and tells the optimizer that the parameters are bad
    if len(B_model) != len(data.b):
        return 1e10

    x2_B = sum(((B_model - data.b)/data.b_err)**2)
    x2_d = sum(((d_model - data.d)/data.d_err)**2)

    #this is the main chi squared that we are minimizing
    #it's not a "chi-squared" anymore, but multiplying by a constant seems dumb
    x2 = x2_B + x2_d + costs #Willett et al. 2009, give or take

    if vx_flag:
        x2_vx = sum(((vx_model - data.vx)/data.vx_err)**2)
        x2 += x2_vx

    if vy_flag:
        x2_vy = sum(((vy_model - data.vy)/data.vy_err)**2)
        x2 += x2_vy

    if vz_flag:
        x2_vz = sum(((vz_model - data.vz)/data.vz_err)**2)
        x2 += x2_vz

    if vgsr_flag:
        x2_vgsr = sum(((vgsr_model - data.vgsr)/data.vgsr_err)**2)
        x2 += x2_vgsr

    '''
    #get normalization factor
    N = len(data.l) #number of data points
    n = 5 #number of parameters
    eta = N - n - 1 #normalizing parameter
    if eta <= 0:
        eta = 1 #if you use fewer data points than needed to constrain the problem, then this will still work but it won't be normalized correctly

    x2 = (1/eta) * (x2_B + x2_d + x2_vx + x2_vy + x2_vz + x2_vgsr) + costs #Willett et al. 2009, give or take
    '''

    #there's a weird edge case where occasionally x2 is a short array of floats
    #this bit prevents scipy from throwing an error
    #no idea what causes that
    if type(x2) == type(np.array([])):
        x2 = x2[0]

    if verbose:
        print('GoF: ' + str(x2))

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

    params = scopt.differential_evolution(chi_squared, bounds, args=(data_opt,), maxiter=max_it, popsize=pop_size, mutation=diff_scaling_factor, recombination=crossover_rate, workers=-1, disp=not(verbose), **kwargs).x
    #'''

    x2 = chi_squared(params, data_opt)

    return params, x2

'''
================================================================================
 FUNCTIONS
================================================================================
'''

def fit_orbit(l, b, b_err, d, d_err, vx=None, vy=None, vz=None, vgsr=None, \
              vx_err=None, vy_err=None, vz_err=None, vgsr_err=None, max_it=100, \
              bounds=[(0, 360), (-90, 90), (0, 100), (-500, 500), (-500, 500), (-500, 500)], \
              t_len=None, **kwargs):

    #construct data
    #set proper flags based on input data
    if vx is not None:
        global vx_flag
        vx_flag = True
    if vy is not None:
        global vy_flag
        vy_flag = True
    if vz is not None:
        global vz_flag
        vz_flag = True
    if vgsr is not None:
        global vgsr_flag
        vgsr_flag = True

    #update t_length if necessary
    if t_len is not None:
        global t_length
        t_length = t_len
        global ts
        ts = np.linspace(0, t_length, num=resolution)*u.Gyr

    if verbose:
        print('Building orbit data...')

    data_opt = OrbitData(l, b, d, vx, vy, vz, vgsr, b_err, d_err, vx_err, vy_err, vz_err, vgsr_err)

    if verbose:
        print('===================================')
        print('Optimizing:')
        print('===================================')

    #optimization
    params, x2 = optimize(data_opt, max_it, bounds, **kwargs)

    print('===================================')
    print('Params: l, b, d, vx, vy, vz')
    print(params)
    print()
    print('Goodness-of-Fit:')
    print(x2)
    print('===================================')

    return params, x2

'''
================================================================================
 PLOTTING
================================================================================
'''

def plot_orbit_gal(l, b, d, params):

    o = make_orbit(params)

    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    data_orbit, data_orbit_rev = get_OrbitData_from_orbit(o, o_rev)

    fig = plt.figure(figsize=(18, 8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

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

    plt.show()

def plot_orbit_icrs(l, b, d, params):

    s = SkyCoord(l, b, frame='galactic', unit=(u.deg, u.deg))
    s = s.transform_to('icrs')
    ra = s.ra
    dec = s.dec

    o = make_orbit(params)

    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    data_orbit, data_orbit_rev = get_OrbitData_from_orbit(o, o_rev)
    data_orbit.icrs()
    data_orbit_rev.icrs()

    fig = plt.figure(figsize=(18, 8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot(data_orbit.ra, data_orbit.dec, c='b')
    ax1.plot(data_orbit_rev.ra, data_orbit_rev.dec, c='r')
    ax1.scatter(ra, dec, c='k')

    ax1.set_xlim(0, 360)
    ax1.set_ylim(-90, 90)
    ax1.set_xlabel('l')
    ax1.set_ylabel('b')

    ax2.plot(data_orbit.ra, data_orbit.d, c='b')
    ax2.plot(data_orbit_rev.ra, data_orbit_rev.d, c='r')
    ax2.scatter(ra, d, c='k')

    ax2.set_xlim(0, 360)
    ax2.set_xlabel('l')
    ax2.set_ylabel('d (helio)')

    plt.show()

'''
================================================================================
 TESTING
================================================================================
'''

def test():
    test_o_l = 90
    test_o_b = 20
    test_o_d = 25
    test_o_vx = 150
    test_o_vy = -200
    test_o_vz = 100

    ts = np.linspace(0, 0.25, 1000)*u.Gyr
    sample = np.array([100, 250, 400, 500, 600, 750, 850])

    o = make_orbit([test_o_l, test_o_b, test_o_d, test_o_vx, test_o_vy, test_o_vz])

    l = np.take(o.ll(ts), sample)
    b = np.take(o.bb(ts), sample)
    d = np.take(o.dist(ts), sample)

    plt.scatter(l, b, s=5)
    plt.scatter(o.ll(ts), o.bb(ts), s=1)
    plt.show()

    plt.scatter(l, d, s=5)
    plt.scatter(o.ll(ts), o.dist(ts), s=1)
    plt.show()

    b_err = np.zeros(len(b)) + 0.05
    d_err = np.zeros(len(d)) + 1

    test_orbit_data = OrbitData(l, b, d, None, None, None, None, b_err, d_err, None, None, None, None)

    print('Goodness-of-Fit of actual values:')
    print(chi_squared([test_o_l, test_o_b, test_o_d, test_o_vx, test_o_vy, test_o_vz], data=test_orbit_data))

    params, x2 = fit_orbit(l, b, b_err, d, d_err)
    plot_orbit_gal(l, b, d, params)
    plot_orbit_icrs(l, b, d, params)

'''
================================================================================
 RUNTIME
================================================================================
'''

if __name__ == "__main__":
    test()
