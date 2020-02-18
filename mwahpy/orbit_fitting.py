#============================
#Tom Donlon, 2018
#============================

#TODO: Allow different coordinate systems input

import numpy as np
import scipy as sc
from astropy import units as u
import scipy.optimize as scopt
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import galpy
import coords as co
import flags
import global as glob

from glob import pot

'''
================================================================================
FLAGS
================================================================================
'''

#-------------------------------------------------------------------------------
#DO NOT TOUCH
#Any routine meant to be used by an end user will configure this from input data
d_flag = 0
vx_flag = 0
vy_flag = 0
vz_flag = 0
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

'''
================================================================================
'''

class OrbitData():
    #all of these parameters can be np.arrays of floats
    def __init__(self, l, b, d, vx, vy, vz, b_err, d_err, vx_err, vy_err, vz_err):
        self.l = l
        self.b = b
        if d_flag:
            self.d = d #heliocentric distance
        else:
            self.d = np.zeros(len(l)) + 1 #populate with ones so that Lam/Bet calculation still functions
        self.vx = vx
        self.vy = vy
        self.vz = vz

        self.l_err = l_err
        self.b_err = b_err
        self.d_err = d_err
        self.vx_err = vx_err
        self.vy_err = vy_err
        self.vz_err = vz_err

        self.x = d*np.cos(self.l)*np.cos(self.b) - 8
        self.y = d*np.sin(self.l)*np.cos(self.b)
        self.z = d*np.sin(self.b)

    def LamBet(self, normal, point):
        #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
        #point: parameter for the axis generation of the Great Circle coordinates
        self.L, self.B = co.gal2LamBet(self.x, self.y, self.z, normal=normal, point=point)
        self.r = (self.x**2 + self.y**2 + self.z**2)**0.5
        self.L_err = self.l_err #should be error propogated (Newby et al. 2013) but is not yet
        self.B_err = self.b_err #should be error propogated (Newby et al. 2013) but is not yet
        self.D_err = self.d_err #should be error propogated but is not yet
        return self

def getClosestIndex(val, Lam):
    Lam = np.abs(Lam - val)
    m = np.min(Lam)
    ind = np.where(Lam == m)

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
def getModelFromOrbit(data, o, normal, point):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    #initialize the orbit we are fitting --
    #we flip it around so that we are fitting both the forwards and the backwards orbit
    ts = np.linspace(0, t_length, num=resolution)*u.Gyr
    o_rev = o.flip()
    o_rev.integrate(ts, pot)

    #sign swap on vx because galpy is left-handed, and we are inputting data in a right-handed coordinate system
    data_orbit = OrbitData(np.array(o.ll(ts)), np.array(o.bb(ts)), np.array(o.dist(ts)), np.array(o.vx(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o.vy(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o.vz(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array([]), np.array([]), np.array([]), np.array([]), , np.array([]))
    data_orbit = data_orbit.LamBet(normal, point)
    data_orbit_rev = OrbitData(np.array(o_rev.ll(ts)), np.array(o_rev.bb(ts)), np.array(o_rev.dist(ts)), np.array(o_rev.vx(ts, obs=[8., 0., 0., 0., 0., 0.])), np.array(o_rev.vy(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array(o_rev.vz(ts, obs=[8., 0., 0., 0., 0., 0.]))*-1, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    data_orbit_rev = data_orbit_rev.LamBet(normal, point)

    #grab full lists so that we can select the closest points once we get a list
    Lam = np.append(np.flip(data_orbit_rev.L), data_orbit.L)

    #get the list of points closest to each data point in Lambda
    point_list, costs = getPointList(data.L, Lam)

    #grab the model points from the point list we grabbed
    Bet = np.append(np.flip(data_orbit_rev.B), data_orbit.B)
    B_model = np.array([Bet[p] for p in point_list])
    if d_flag:
        D = np.append(np.flip(data_orbit_rev.D), data_orbit.D)
        D_model = np.array([D[p] for p in point_list])
    else:
        D_model = np.zeros(len(B_model))

    if vx_flag:
        vx = np.append(np.flip(data_orbit_rev.vx), data_orbit.vx)
        vx_model = np.array([vx[p] for p in point_list])
    else:
        vx_model = np.zeros(len(B_model))

    if vy_flag:
        vy = np.append(np.flip(data_orbit_rev.vy), data_orbit.vy)
        vy_model = np.array([vy[p] for p in point_list])
    else:
        vy_model = np.zeros(len(B_model))

    if vz_flag:
        vz = np.append(np.flip(data_orbit_rev.vz), data_orbit.vz)
        vz_model = np.array([vz[p] for p in point_list])
    else:
        vz_model = np.zeros(len(B_model))

    return B_model, D_model, vx_model, vy_model, vz_model, costs

#chi_squared: data, galpy.Orbit() --> float
#takes in the observed data and a test orbit and calculates the goodness-of-fit using a chi-squared method
def chiSquared(params, data=[], normal=(0, 0, 0), point=(1, 0, 0)):
    #data: the data that the orbit is being fit to
    #o: the test orbit that we are calculating the goodness-of-fit of
    #normal: the normal vector to the plane of the Great Circle we are estimating for the orbit
    #point: parameter for the axis generation of the Great Circle coordinates

    o = Orbit(vxvv=[params[0], params[1], params[2], params[3], params[4] - 220, params[5]], uvw=True, lb=True, ro=8., vo=220.) #generate the orbit
    o.integrate(ts, pot) #integrate the orbit

    B_model, d_model, vx_model, vy_model, vz_model, costs = getModelFromOrbit(data, o, normal, point) #get model data from orbit

    x2_B = sum(((B_model - data.B)/data.B_err)**2)
    if d_flag:
        x2_d = sum(((d_model - data.d)/data.d_err)**2)
    else:
        x2_d = 0

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

    #get normalization factor
    N = len(data.L) #number of data points
    n = 5 #number of parameters
    eta = N - n - 1 #normalizing parameter
    if eta <= 0:
        eta = 1 #if you use fewer data points than needed to constrain the problem, then this will still work but it won't be normalized correctly

    x2 = (1/eta) * (x2_B + x2_D + x2_vx + x2_vy + x2_vz) + costs #Willett et al. 2009

    if flags.verbose:
        print('=================================')
        print('Optimizing:')
        print('Chi Squared: ' + str(x2))
        print('=================================')

    return x2

#optimize: data -> [float, float, float, float, float], (float, float, float), (float, float, float)
#takes in data, then fits a Great Circle to that data and minimizes the chi_squared to fit an orbit to the data
def optimize(data_opt, max_it, bounds):
    normal = co.getPlaneNormal(gc.plane_OLS(data_opt.x, data_opt.y, data_opt.z)) #get normal vector for fit Great Circle
    point = (1, 0, 0) #not currently fitting this, but this can be changed or fit at a later date
    #this way it makes 0 deg. in the Great Circle ~0 deg. in galactic longitude

    data_opt = data_opt.LamBet(normal, point)

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

    params = scopt.differential_evolution(chiSquared, bounds, args=(data_opt, normal, point), strategy='rand1bin', maxiter=max_it, popsize=pop_size, mutation=diff_scaling_factor, recombination=crossover_rate).x
    #'''

    x2 = chiSquared(params, data_opt, normal, point)

    return params, normal, point, x2

def fit_orbit(l, b, b_err, d=None, vx=None, vy=None, vz=None, d_err=None, \
              vx_err=None, vy_err=None, vz_err=None, max_it=20, \
              bounds=[(0, 360), (-90, 90), (0, 100), (-1000, 1000), (-1000, 1000), (-1000, 1000)], \
              t_len=None):

    #construct data
    #set proper flags based on input data
    if d:
        global d_flag
        d_flag = 1
    if vx:
        global vx_flag
        vx_flag = 1
    if vy:
        global vy_flag
        vy_flag = 1
    if vz:
        global vz_flag
        vz_flag = 1

    #update t_length if necessary
    if t_len:
        global t_length:
        t_length = t_len
        global ts
        ts = np.linspace(0, t_length, num=resolution)*u.Gyr

    data_opt = Data(l, b, d, vx, vy, vz, b_err, vx_err, vy_err, vz_err)

    #optimization
    params, normal, point, x2 = optimize(data_opt, max_it, bounds)

    if flags.verbose:
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
        print('Chi Squared:')
        print(x2)
        print('===================================')

    return params, normal, point, x2
