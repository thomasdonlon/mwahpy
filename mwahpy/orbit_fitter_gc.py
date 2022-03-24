'''
This is a self-contained orbit fitting routine.

This orbit fitter is unique compared to other common orbit fitters in that it
uses a galactocentric generalized plane coordinate system when fitting data
'''

#TODO: Allow different coordinate systems input
#TODO: more settings put in through fit_orbit()

import numpy as np
import scipy.optimize as scopt
import matplotlib.pyplot as plt

from astropy import units as u
from galpy.orbit import Orbit

from .coords import get_plane_normal, plane_OLS, gal_to_lambet_galcentric
from .flags import verbose
from .pot import mwahpy_default_pot
from .orbit_fitter import make_orbit, get_OrbitData_from_orbit, get_point_list, OrbitData

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

#i think these can get mismatched between orbit_fitter.py and this file, so maybe
#we fix that if it's a problem

t_length = 0.5 #Gyr
resolution = 1000
ts = np.linspace(0, t_length, num=resolution)*u.Gyr

'''
================================================================================
 HELPER FUNCTIONS FOR OPTIMIZATION
================================================================================
'''

#take in data, orbit, and plane info: Output model data corresponding to each data point
def get_model_from_params(data, params, normal, point):

    #initialize the orbit we are fitting
    o = make_orbit(params)

    #we flip it around so that we are fitting both the forwards and the backwards orbit
    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    data_orbit, data_orbit_rev = get_OrbitData_from_orbit(o, o_rev)
    data_orbit.calc_lam_bet(normal, point)
    data_orbit_rev.calc_lam_bet(normal, point)

    #grab full lists so that we can select the closest points once we get a list
    Lam = np.append(np.flip(data_orbit_rev.L), data_orbit.L)

    #get the list of points closest to each data point in Lambda
    point_list, costs = get_point_list(data.L, Lam)

    #grab the model points from the point list we grabbed
    Bet = np.append(np.flip(data_orbit_rev.B), data_orbit.B)
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
        vgsr = np.append(np.flip(data_orbit_rev.vlos), data_orbit.vlos)
        vgsr_model = np.array([vgsr[p] for p in point_list]).flatten()
    else:
        vgsr_model = np.zeros(len(B_model))

    return B_model, D_model, vx_model, vy_model, vz_model, vgsr_model, costs

#-------------------------------------------------------------------------------

#takes in the observed data and a test orbit and calculates the goodness-of-fit using a chi-squared method
def chi_squared(params, data=[], normal=(0, 0, 0), point=(1, 0, 0)):

    B_model, d_model, vx_model, vy_model, vz_model, vgsr_model, costs = get_model_from_params(data, params, normal, point) #get model data from orbit

    x2_B = sum(((B_model - data.B)/data.B_err)**2)
    x2_d = sum(((d_model - data.d)/data.d_err)**2)

    #this is the main chi squared that we are minimizing
    x2 = x2_B + x2_d + costs

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

    #there's a weird edge case where occasionally x2 is a short array of floats
    #this bit prevents scipy from throwing an error
    #no idea what causes that
    if type(x2) == type(np.array([])):
        x2 = x2[0]

    if verbose:
        print('GoF: ' + str(x2))

    return x2

#-------------------------------------------------------------------------------

#takes in data, then fits a Great Circle to that data and
#minimizes the chi_squared to fit an orbit to the data
def optimize(data_opt, max_it, bounds, **kwargs):

    #compute the preferred galactocentric spherical frame to fit this orbit in
    #(should be close to the orbital plane of the observed data)
    normal = get_plane_normal(plane_OLS(data_opt.x, data_opt.y, data_opt.z)) #get normal vector for fit Great Circle

    #construct the point vector (must be orthogonal to normal and be in the x-z plane)
    px = 1 / ((normal[0]/normal[2])**2 + 1)**0.5
    pz = -1 * px * (normal[0] / normal[2])
    point = (px, 0, pz)

    data_opt.calc_lam_bet(normal, point)

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

    #the actual bulk of the optimization
    params = scopt.differential_evolution(chi_squared, bounds, args=(data_opt, normal, point), maxiter=max_it, popsize=pop_size, mutation=diff_scaling_factor, recombination=crossover_rate, workers=-1, disp=not(verbose), **kwargs).x

    x2 = chi_squared(params, data_opt, normal, point)

    return params, normal, point, x2

'''
================================================================================
 FUNCTIONS
================================================================================
'''

#pass in l, b, d values (and optional vx, vy, vz, vlos velocities)
#computes the best orbit fit to those parameters
#outputs the parameters for that best orbit fit
#also need to pass in errors for each value used (except l)

#additional kwargs are passed to scipy.optimize.differential_evolution()
def fit_orbit(l, b, b_err, d, d_err, vx=None, vy=None, vz=None, vgsr=None, \
              vx_err=None, vy_err=None, vz_err=None, vgsr_err=None, max_it=100, \
              bounds=[(0, 360), (-90, 90), (0, 100), (-1000, 1000), (-1000, 1000), (-1000, 1000)], \
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
    params, normal, point, x2 = optimize(data_opt, max_it, bounds, **kwargs)

    if verbose:
        print('===================================')
        print('Params: l, b, d, vx, vy, vz')
        print(params)
        print()
        print('Normal Vector:')
        print(normal)
        print()
        print('Point Vector:')
        print(point)
        print()
        print('Goodness-of-Fit:')
        print(x2)
        print('===================================')

    return params, normal, point, x2

'''
================================================================================
 PLOTTING
================================================================================
'''

def plot_orbit_lam_bet(L, B, params, normal, point):

    o = make_orbit(params)

    o_rev = o.flip()
    o_rev.integrate(ts, mwahpy_default_pot)

    data_orbit, data_orbit_rev = get_OrbitData_from_orbit(o, o_rev)
    data_orbit.calc_lam_bet(normal, point)
    data_orbit_rev.calc_lam_bet(normal, point)

    fig = plt.figure(figsize=(9, 6))

    plt.plot(data_orbit.L, data_orbit.B, c='b')
    plt.plot(data_orbit_rev.L, data_orbit_rev.B, c='r')
    plt.scatter(L, B, c='k')

    plt.xlim(0, 360)
    plt.ylim(-90, 90)
    plt.xlabel('$\\Lambda$')
    plt.ylabel('$\\beta$')

    plt.show()

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
    test_orbit_data.calc_lam_bet((0, 0, 1), (1, 0, 0))

    print('Goodness-of-Fit of actual values:')
    print(chi_squared([test_o_l, test_o_b, test_o_d, test_o_vx, test_o_vy, test_o_vz], data=test_orbit_data, normal=(0, 0, 1), point=(1, 0, 0)))

    params, normal, point, x2 = fit_orbit(l, b, b_err, d, d_err)
    L, B = gal_to_lambet_galcentric(l, b, d, normal, point)
    plot_orbit_lam_bet(L, B, params, normal, point)
    plot_orbit_gal(l, b, d, params)

'''
================================================================================
 RUNTIME
================================================================================
'''

if __name__ == "__main__":
    test()
