import numpy as np
from numpy import linalg
import scipy as sc
import math as ma
import unittest

#TODO: Make sure all transformations have inverse transformations
#TODO: Make sure all transormations have unit tests
#TODO: Comment out the code
#TODO: Standardize function naming conventions
#TODO: (low priority) Allow certain quantities (like galactocentric phi/theta) to be expressed/input in radians
#TODO: (low priority) add left-handed versions to all galactic cartesian transformations

#===============================================================================
#HELPER FUNCTIONS
#===============================================================================

def wrapLong(long):

    use_array = True
    if type(long) != type(np.array([])):
        use_array = False
        long = np.array([long])

    for i in range(len(long)):
        if long[i] < 0:
            long[i] += 360
        elif long[i] > 360:
            long[i] -= 360

    if not use_array:
        long = long[0]

    return long

#rotate the given data around the provided axis
def rotAroundArbAxis(x, y, z, ux, uy, uz, theta):
    #x, y, z: 3D cartesian coordinates of the data
    #ux, uy, uz: 3d cartesian coordinates of the axis vector u = (ux, uy, uz)
    #theta: angle to rotate data counter-clockwise around axis (rad)

    #make sure that u is normalized
    norm_u = (ux**2 + uy**2 + uz**2)**0.5
    if round(norm_u, 8 != 1.):
        ux /= norm_u
        uy /= norm_u
        uz /= norm_u

    xyz = np.array([x, y, z])

    cos = np.cos(theta)
    sin = np.sin(theta)
    R = np.array([[cos+ux**2*(1-cos),    ux*uy*(1-cos)-uz*sin, ux*uz*(1-cos)+uy*sin],
                  [uy*ux*(1-cos)+uz*sin, cos+uy**2*(1-cos),    uy*uz*(1-cos)-ux*sin],
                  [uz*ux*(1-cos)-uy*sin, uz*uy*(1-cos)+ux*sin, cos+uz**2*(1-cos)   ]])

    xyz = np.matmul(R, xyz)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    return x, y, z

#===============================================================================
#FIRST ORDER COORDINATE/POSITION TRANSFORMATIONS
#===============================================================================

def cart_to_gal(x, y, z, left_handed=False):
    #get l, b, r (helio) from galactocentric X, Y, Z coords
    use_array = True
    if type(x) != type(np.array([])):
        use_array = False
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])

    if left_handed:
        r = ((x-8)**2 + y**2 + z**2)**0.5
        l = np.arctan2(y,-1*(x-8))*180/np.pi
    else:
        r = ((x+8)**2 + y**2 + z**2)**0.5
        l = np.arctan2(y,(x+8))*180/np.pi
    b = np.arcsin(z/r)*180/np.pi

    if not use_array:
        l = l[0]
        b = b[0]
        r = r[0]

    return l, b, r

#-------------------------------------------------------------------------------

def gal_to_cart(l, b, r, left_handed=False):

    use_array = True
    if type(l) != type(np.array([])):
        use_array = False
        l = np.array([l])
        b = np.array([b])
        r = np.array([r])

    l_rad = l*np.pi/180
    b_rad = b*np.pi/180

    if left_handed:
        #left-handed
        x = 8 - r*np.cos(l_rad)*np.cos(b_rad)
    else:
        #right-handed
        x = r*np.cos(l_rad)*np.cos(b_rad) - 8

    y = r*np.sin(l_rad)*np.cos(b_rad)
    z = r*np.sin(b_rad)

    if not use_array:
        x = x[0]
        y = y[0]
        z = z[0]

    return x, y, z

#-------------------------------------------------------------------------------

#input: x[kpc], y[kpc], z[kpc]
#output: Cylindrical coordinates (R, [kpc], Z [kpc], Phi [deg])
#phi = 0 when oriented along the positive x-axis
def cart_to_cyl(x, y, z):

    use_array = True
    if type(x) != type(np.array([])):
        use_array = False
        x = np.array([x])
        y = np.array([y])

    R = (x**2 + y**2)**0.5
    phi = np.arctan2(y, x)*180/np.pi

    if not use_array:
        R = R[0]
        phi = phi[0]

    return R, z, phi

#-------------------------------------------------------------------------------

def cyl_to_cart(R, z, phi):

    use_array = True
    if type(R) != type(np.array([])):
        use_array = False
        R = np.array([R])
        phi = np.array([phi])

    phi = phi*np.pi/180

    x = R*np.cos(phi)
    y = R*np.sin(phi)

    if not use_array:
        x = x[0]
        y = y[0]

    return x, y, z

#-------------------------------------------------------------------------------

def cart_to_sph(x, y, z):

    use_array = True
    if type(x) != type(np.array([])):
        use_array = False
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])

    r = (x**2 + y**2 + z**2)**0.5

    phi = np.arctan2(y, x)*180/np.pi
    theta = np.arcsin(z/r)*180/np.pi

    if not use_array:
        phi = phi[0]
        theta = theta[0]
        r = r[0]

    return phi, theta, r

#-------------------------------------------------------------------------------

def sph_to_cart(phi, theta, r):

    use_array = True
    if type(phi) != type(np.array([])):
        use_array = False
        phi = np.array([phi])
        theta = np.array([theta])
        r = np.array([r])

    phi = phi*np.pi/180
    theta = theta*np.pi/180

    x = r*np.cos(phi)*np.cos(theta)
    y = r*np.sin(phi)*np.cos(theta)
    z = r*np.sin(theta)

    if not use_array:
        x = x[0]
        y = y[0]
        z = z[0]

    return x, y, z

#===============================================================================
#SECOND ORDER COORDINATE/POSITION TRANSFORMATIONS
#===============================================================================

def cyl_to_gal(R, z, phi):
    x, y, z = cyl_to_cart(R, z, phi)
    l, b, r = cart_to_gal(x, y, z)
    return l, b, r

#-------------------------------------------------------------------------------

def gal_to_cyl(l, b, r):
    x, y, z = gal_to_cart(l, b, r)
    R, z, phi = cart_to_cyl(x, y, z)
    return R, z, phi

#===============================================================================
#VELOCITY TRANSFORMATIONS
#===============================================================================

# Input: distance [kpc], radial velocity [km/s], RA/DEC [degrees], and pmRA/pmDEC [mas/yr]
# Returns: Galactic vx, vy, vz velocities [km/s]
# NOTE: pmRA = d/dt(RA) * cos(DEC)
# Arguments should be numpy arrays for most efficient usage
def getUVW(dist, rv, ra, dec, pmra, pmde):

    # Conversion from Equatorial (J2000) Cartesian to Galactic Cartesian
    EQ2GC = np.array( [[-0.05487572, -0.87343729, -0.48383453],
                      [  0.49410871, -0.44482923,  0.7469821 ],
                      [ -0.86766654, -0.19807649,  0.45598456]], dtype=np.float32 )

    ra_rad = ra * (np.pi/180.0)
    dec_rad = dec * (np.pi/180.0)

    sina = np.sin(ra_rad)
    cosa = np.cos(ra_rad)
    sind = np.sin(dec_rad)
    cosd = np.cos(dec_rad)

    vra =  4.741067035842384 * pmra * dist
    vdec = 4.741067035842384 * pmde * dist

    vx_Eq = rv * cosd * cosa   -   vdec * sind * cosa   -   vra * sina
    vy_Eq = rv * cosd * sina   -   vdec * sind * sina   +   vra * cosa
    vz_Eq = rv * sind          +   vdec * cosd

    vel_Eq = np.array([ vx_Eq, vy_Eq, vz_Eq ])
    U, V, W = np.dot(EQ2GC, vel_Eq)

    return U, V, W

#-------------------------------------------------------------------------------

def getvxvyvz(dist, rv, ra, dec, pmra, pmde):

    U, V, W = getUVW(dist, rv, ra, dec, pmra, pmde)

    # Sun's velocity is (10.1, 224.0, 6.7)_GSR
    vx = U + 10.1
    vy = V + 224.0
    vz = W + 6.7

    return vx, vy, vz

#-------------------------------------------------------------------------------

#solar reflex motion will be already removed if UVW are galactocentric
#inputs must be arrays, even if just of length 1
def getrvpm(ra, dec, dist, U, V, W):
    k = 4.74057
    rv = np.array([0.0]*len(ra))
    pmra = np.array([0.0]*len(ra))
    pmdec = np.array([0.0]*len(ra))

    i = 0
    while i < len(ra):

        ra_rad = ra[i] * np.pi/180
        dec_rad = dec[i] * np.pi/180

        T = np.array( [[-0.05487572, -0.87343729, -0.48383453],
                      [  0.49410871, -0.44482923,  0.7469821 ],
                      [ -0.86766654, -0.19807649,  0.45598456]], dtype=np.float32 )
        A = np.array([[np.cos(ra_rad)*np.cos(dec_rad), -1*np.sin(ra_rad), -1*np.cos(ra_rad)*np.sin(dec_rad)], [np.sin(ra_rad)*np.cos(dec_rad), np.cos(ra_rad), -1*np.sin(ra_rad)*np.sin(dec_rad)], [np.sin(dec_rad), 0, np.cos(dec_rad)]])
        B = np.matmul(T,A)
        B_inv = np.linalg.inv(B)

        uvw = np.array([[U[i]], [V[i]], [W[i]]])
        rvpm = np.matmul(B_inv, uvw)

        rv[i] = rvpm[0][0]
        pmra[i] = rvpm[1][0]/dist[i]/k
        pmdec[i] = rvpm[2][0]/dist[i]/k

        i += 1

    return rv, pmra, pmdec

#-------------------------------------------------------------------------------

#inputs must be arrays, even if just of length 1
def removeSolarMotionFromPM(ra, dec, dist, pmra, pmdec):
    vx = np.array([-10] * len(ra))
    vy = np.array([-224] * len(ra))
    vz = np.array([-7] * len(ra))

    rv, mura, mudec = getrvpm(ra, dec, dist, vx, vy, vz)

    pmra -= mura
    pmdec -= mudec

    return pmra, pmdec

#-------------------------------------------------------------------------------

#does not work for arrays!
#TODO: make work with arrays
def getUVWerrors(dist, ra, dec, pmra, pmdec, err_pmra, err_pmdec, err_rv, err_dist):
    #distance in pc
    k = 4.74057

    ra_rad = ra * np.pi/180
    dec_rad = dec * np.pi/180

    T = np.array([[-0.05487572,-0.87343729, -0.48383453], [0.49410871, -0.44482923,  0.7469821], [-0.86766654, -0.19807649,  0.45598456]])
    A = np.array([[np.cos(ra_rad)*np.cos(dec_rad), -1*np.sin(ra_rad), -1*np.cos(ra_rad)*np.sin(dec_rad)], [np.sin(ra_rad)*np.cos(dec_rad), np.cos(ra_rad), -1*np.sin(ra_rad)*np.sin(dec_rad)], [np.sin(dec_rad), 0, np.cos(dec_rad)]])
    B = np.matmul(T,A)
    C = B**2

    M = np.array([[err_rv**2],
                  [(k*dist)**2 * (err_pmra**2 + (pmra*err_dist/dist)**2)],
                  [(k*dist)**2 * (err_pmdec**2 + (pmdec*err_dist/dist)**2)]])
    N = 2*pmra*pmdec*k**2*err_dist**2*np.array([[B[0][1]*B[0][2]],
                                               [B[1][1]*B[1][2]],
                                               [B[2][1]*B[2][2]]])

    uvw_var = np.matmul(C,M) + N

    err_u = (uvw_var[0][0])**0.5
    err_v = (uvw_var[1][0])**0.5
    err_w = (uvw_var[2][0])**0.5

    return err_u, err_v, err_w

#=====================================
#MISC TOOLS
#=====================================

#calculate approximate local standard of rest (LSR) velocity at a point in the galactic plane
def approx_lsr(x, y):
    v_lsr = 220
    vx = v_lsr / (1 + (x/y)**2)**0.5
    vy = v_lsr / (1 + (y/x)**2)**0.5
    if x > 0:
        vy = -1 * vy
    if y < 0:
        vx = -1 * vx
    return vx, vy

#-------------------------------------------------------------------------------

#borrowed from NewbyTools
def plane_dist(x,y,z, params):
    a,b,c = params
    return (a*x + b*y + c*z)

#-------------------------------------------------------------------------------

#x, y, z should be numpy arrays
#takes in galactocentric x, y, z data and outputs parameters for the best fit plane to those points incident the galactic center
#adapted from NewbyTools
def plane_OLS(x,y,z, print_distances=False):
    """ Solves for the best-fit plane to a set of x,y,z data using ordinary least-squares.
        Equation is of form z = Ax + By.
        DIFFERENT FROM NEWBY TOOLS plane_OLS() IN THAT WE CONSTRAIN THE PLANE THROUGH GALACTIC CENTER
        Output is normalized a,b,c of a plane of the form ax+by+cz=0"""
    A = np.array([x, y]).T
    B = z.T

    #solve Ax=B
    p = np.matmul(np.matmul(linalg.inv(np.matmul(A.T, A)), A.T), B.T) #uses left pseudo-inverse {(A^T * A)^-1 * A^T} due to system being overconstrained (A doesn't have a true inverse)
    params = [-float(p[0]), -float(p[1]), 1.0]  #c=1.0 by default
    bottom = np.sqrt(params[0]*params[0] + params[1]*params[1] + params[2]*params[2])
    for i in range(len(params)):  params[i] = params[i]/bottom
    print("# - Normalized best-fit plane parameters: {0}".format(params))
    if print_distances:
        for i in range(len(x)):
            print(plane_dist(x[i], y[i], z[i], params))
    return params

#-------------------------------------------------------------------------------

#getPlaneNormal: [float, float, float, float] --> np.array([float, float, float])
#takes in parameters that define a plane in 3D and returns a normalized normal vector to that plane
def getPlaneNormal(params):
    #params ([a, b, c, d]) corresponding to the equation for a plane ax + by + cz + d = 0
    #comes from the plan fitting method above

    #definition of a normal vector, given equation of a plane
    normal = np.array([params[0], params[1], params[2]])

    #normalize the normal vector
    len_normal = (normal[0]**2 + normal[1]**2 + normal[2]**2)**0.5
    normal = (normal[0]/len_normal, normal[1]/len_normal, normal[2]/len_normal)

    return normal

#===============================================================================
#NONSTANDARD/UNIQUE COORDINATE TRANFORMATIONS
#===============================================================================

#gal2plane: np.array(floats), np.array(floats), np.array(floats), (float, float, float), (float, float, float) --> np.array(floats), np.array(floats), np.array(floats)
#takes in galactic coordinates for a star(s) and returns their x,y,z coordinates with respect to a rotated plane with the normal vector provided
#Newby 2013 et al, appendix
def cart_to_plane(x, y, z, normal, point):

    #ensure that normal and point are normalized
    len_normal = (normal[0]**2 + normal[1]**2 + normal[2]**2)**0.5
    if len_normal != 1.0:
        normal = (normal[0]/len_normal, normal[1]/len_normal, normal[2]/len_normal)
    len_point = (point[0]**2 + point[1]**2 + point[2]**2)**0.5
    if len_point != 1.0:
        point = (point[0]/len_point, point[1]/len_point, point[2]/len_point)

    #define new axes along the plane
    z_plane = np.array(normal)
    y_plane = np.cross(z_plane, np.array(point))
    y_plane = y_plane / (y_plane[0]**2 + y_plane[1]**2 + y_plane[2]**2)**0.5 #normalize y_plane vector to prevent skewing
    x_plane = np.cross(y_plane, z_plane)
    x_plane = x_plane / (x_plane[0]**2 + x_plane[1]**2 + x_plane[2]**2)**0.5 #normalize x_plane vector to prevent skewing

    #get new x, y, z through change of basis
    xyz = np.array([x, y, z])

    COB = linalg.inv(np.array([x_plane, y_plane, z_plane]).T)
    new_xyz = np.matmul(COB, xyz)

    return new_xyz[0], new_xyz[1], new_xyz[2]

#gal2plane: np.array(floats), np.array(floats), np.array(floats), (float, float, float), (float, float, float) --> np.array(floats), np.array(floats), np.array(floats)
#takes in galactic coordinates for a star(s) and returns their x,y,z coordinates with respect to a rotated plane with the normal vector provided
#Newby 2013 et al, appendix
def gal_to_plane(l, b, d, normal, point):

    x, y, z = gal_to_cart(l, b, d)
    x, y, z = cart_to_plane(x, y, z, normal, point)

    return x, y, z

#-------------------------------------------------------------------------------

def gal_to_LamBet(l, b, d, normal, point):

    x_prime, y_prime, z_prime = gal_to_plane(l, b, d, normal, point)
    Lam = np.arctan2(y_prime, x_prime)*180/np.pi #convert to degrees
    #correct Lam to be between 0 and 360 instead of -180 to 180
    i = 0
    while i < len(Lam):
        if Lam[i] < 0:
            Lam[i] += 360
        i += 1

    Bet = np.arcsin(z_prime/(x_prime**2 + y_prime**2 + z_prime**2)**0.5)*180/np.pi #convert to degrees

    return Lam, Bet

#-------------------------------------------------------------------------------

def cart_to_LamBet(x,y,z, normal, point):

    x_prime, y_prime, z_prime = cart_to_plane(x,y,z, normal=normal, point=point)
    Lam = np.arctan2(y_prime, x_prime)*180/np.pi #convert to degrees
    #correct Lam to be between 0 and 360 instead of -180 to 180
    i = 0
    while i < len(Lam):
        if Lam[i] < 0:
            Lam[i] += 360
        i += 1

    Bet = np.arcsin(z_prime/(x_prime**2 + y_prime**2 + z_prime**2)**0.5)*180/np.pi #convert to degrees
    
    return Lam, Bet

#-------------------------------------------------------------------------------

def xyz2Sgr(x,y,z):

    x_prime = x_prime + 8
    xyz = np.array([x, y, z])

    #Euler rotation matrix from Solar frame into Sgr
    phi = 183.8 * np.pi/180
    theta = 76.5 * np.pi/180
    psi = 194.1 * np.pi/180
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi= np.cos(psi)
    spsi= np.sin(psi)

    R = np.array([[cpsi*cphi-ctheta*sphi*spsi, cpsi*sphi+ctheta*cphi*spsi, spsi*stheta],
                   [-1*spsi*cphi-ctheta*sphi*cpsi, -1*spsi*sphi+ctheta*cphi*cpsi, cpsi*stheta],
                   [stheta*sphi, -1*stheta*cphi, ctheta]])

    xyz = np.matmul(R, xyz)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    Lam = np.arctan2(y, x)*180/np.pi
    Bet = np.arcsin(z/(x**2 + y**2 + z**2)**0.5)*180/np.pi

    return Lam, Bet

#-------------------------------------------------------------------------------

#gal2plane: np.array(floats), np.array(floats), np.array(floats), (float, float, float), (float, float, float) --> np.array(floats), np.array(floats)
#takes in galactic coordinates for a star(s) and returns their Lamba, Beta coordinates with respect to a rotated plane with the normal vector provided
def gal2Sgr(l, b):
    #x, y, z: galactocentric x, y, z coordinates
    #normal: 3-vector providing the orientation of the plane to rotate into
    #point: 3-vector 'suggesting' the direction of the new x-axis

    x = np.cos(l*np.pi/180)*np.cos(b*np.pi/180)
    y = np.sin(l*np.pi/180)*np.cos(b*np.pi/180)
    z = np.sin(b*np.pi/180)
    xyz = np.array([x, y, z])

    #Euler rotation matrix from Solar frame into Sgr
    phi = 183.8 * np.pi/180
    theta = 76.5 * np.pi/180
    psi = 194.1 * np.pi/180
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi= np.cos(psi)
    spsi= np.sin(psi)

    R = np.array([[   cpsi*cphi-ctheta*sphi*spsi,    cpsi*sphi+ctheta*cphi*spsi, spsi*stheta],
                  [-1*spsi*cphi-ctheta*sphi*cpsi, -1*spsi*sphi+ctheta*cphi*cpsi, cpsi*stheta],
                  [   stheta*sphi,                -1*stheta*cphi,                ctheta     ]])

    xyz = np.matmul(R, xyz)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    Lam = wrapLong(np.arctan2(y, x)*180/np.pi)
    Bet = np.arcsin(z/(x**2 + y**2 + z**2)**0.5)*180/np.pi

    return Lam, Bet

#-------------------------------------------------------------------------------

#gal2plane: np.array(floats), np.array(floats), np.array(floats), (float, float, float), (float, float, float) --> np.array(floats), np.array(floats)
#takes in galactic coordinates for a star(s) and returns their Lamba, Beta coordinates with respect to a rotated plane with the normal vector provided
def Sgr2gal(Lam, Bet):
    #x, y, z: galactocentric x, y, z coordinates
    #normal: 3-vector providing the orientation of the plane to rotate into
    #point: 3-vector 'suggesting' the direction of the new x-axis

    x = np.cos(Lam*np.pi/180)*np.cos(Bet*np.pi/180)
    y = np.sin(Lam*np.pi/180)*np.cos(Bet*np.pi/180)
    z = np.sin(Bet*np.pi/180)
    xyz = np.array([x, y, z])

    #Euler rotation matrix from Solar frame into Sgr
    phi = 183.8 * np.pi/180
    theta = 76.5 * np.pi/180
    psi = 194.1 * np.pi/180
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi= np.cos(psi)
    spsi= np.sin(psi)

    R = np.array([[   cpsi*cphi-ctheta*sphi*spsi,    cpsi*sphi+ctheta*cphi*spsi, spsi*stheta],
                  [-1*spsi*cphi-ctheta*sphi*cpsi, -1*spsi*sphi+ctheta*cphi*cpsi, cpsi*stheta],
                  [   stheta*sphi,                -1*stheta*cphi,                ctheta     ]])
    R_inv = np.linalg.inv(R)

    xyz = np.matmul(R_inv, xyz)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    l = wrapLong(np.arctan2(y, x)*180/np.pi)
    b = np.arcsin(z/(x**2 + y**2 + z**2)**0.5)*180/np.pi

    return l, b

#-------------------------------------------------------------------------------

#sky2pole: array(float), array(float), tuple(float, float), tuple(float, float) -> array(float), array(float)
#rotate the positions on the sky (sky1, sky2) into a new frame determined by
#the pole of the new frame (pole1, pole2) and the origin of the new frame (origin1, origin2)
#The output is the longitude and latitude in the new coordinate frame

#NOTE: Inputs can be any spherical geometry, as long as the pole & origin arguments
#are in the same coordinate system as the (sky1, sky2) coordinates.
#e.g. if sky1 is a list of RAs and sky2 is a list of Decs, then
#the pole and origin arguments must also be specified in RA Dec.
def sky2pole(sky1, sky2, pole, origin, wrap=False, rad=False):
    #sky1, sky2: positions of the data on the sky (e.g. sky1 = array(RA), sky2 = array(Dec), etc.)
    #pole1, pole2: position of the pole of the new coordinate system
    #origin1, origin2: position of the origin of the new coordinate system
    #wrap: if True, Lam is constrained to only positive values. Otherwise, [-180,180] is the typical Lam output
    #rad: is True, inputs are in radians. ALL inputs should be in degrees if rad=False

    pole1 = pole[0]
    pole2 = pole[1]
    origin1 = origin[0]
    origin2 = origin[1]

    if not rad:
        sky1 *= np.pi/180
        sky2 *= np.pi/180
        pole1 *= np.pi/180
        pole2 *= np.pi/180
        origin1 *= np.pi/180
        origin2 *= np.pi/180

    #construct a 3D frame in the current long/lat frame
    #all points have unit distance
    x = np.cos(sky1)*np.cos(sky2)
    y = np.sin(sky1)*np.cos(sky2)
    z = np.sin(sky2)
    xyz = np.array([x, y, z])

    #rotate old pole into correct lat by rotating around y-axis
    theta = np.pi/2 - pole2 #angle to rotate
    cos = np.cos(theta)
    sin = np.sin(theta)
    Ry = np.array([[cos,    0, sin],
                  [0,      1, 0  ],
                  [-1*sin, 0, cos]])

    #rotate old pole into correct long by rotating around z-axis
    theta = pole1
    cos = np.cos(theta)
    sin = np.sin(theta)
    Rz = np.array([[cos, -1*sin, 0],
                   [sin, cos,    0],
                   [0,   0,      1]])

    #rotate data into new frame
    R = np.matmul(Rz, Ry)
    xyz = np.matmul(R, xyz)
    print('xyz =',xyz)

    #rotate z=0 plane of new pole around new pole to align origin
    r1 = np.array([1, 0, 0]) #rotate old origin into new rotated plane
    r1 = np.matmul(R, r1)

    #construct vector of desired origin
    r2 = np.array([np.cos(origin1)*np.cos(origin2), np.sin(origin1)*np.cos(origin2), np.sin(origin2)])

    #angle to rotate
    psi = np.arccos(np.dot(r1, r2))
    print('psi =', psi)

    #rotate data by psi around the new pole axis
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    x, y, z = rotAroundArbAxis(x, y, z, np.cos(pole1)*np.cos(pole2), np.sin(pole1)*np.cos(pole2), np.sin(pole2), psi)

    #calculate Lam, Bet in new frame
    Lam = np.arctan2(y, x)
    Bet = np.arcsin(z/(x**2 + y**2 + z**2)**0.5)

    if not rad:
        Lam *= 180/np.pi
        Bet *= 180/np.pi

    if wrap:
        Lam = wrapLong(Lam) #TODO: wrapLong should allow for radians

    return Lam, Bet

#===============================================================================
# UNIT TESTING
#===============================================================================
#to run tests, run this file like you would run a regular python script

#to add more tests,
# 1) Find a test class that fits the description of the test you want to add,
#    then add that test to the class as a function beginning with "test_"
# 2) Create a new class beginning with "Test" and composed of just (unittest.TestCase)
#    and then add new functions that begin with "test_" to that class.
#    This should be used if new tests do not fit the descriptions of other
#    test classes

#after changing this file, it should be run in order to make sure all tests pass
#if tests don't pass, congratulations you broke something
#please fix it

prec = 8 #number of digits to round to when comparing values
#WARNING: Tests may fail at high levels of prec or after repeated conversions
#due to computer precision problems

class TestInverses(unittest.TestCase):

    #First order coordinate (position) transformations

    def test_cart_to_gal_rh(self):
        #right-handed
        x, y, z = 1., 2., 3.
        l, b, r = cart_to_gal(x, y, z)
        new_x, new_y, new_z = gal_to_cart(l, b, r)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_gal_lh(self):
        #left-handed
        x, y, z = 1., 2., 3.
        l, b, r = cart_to_gal(x, y, z, left_handed=True)
        new_x, new_y, new_z = gal_to_cart(l, b, r, left_handed=True)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_cyl(self):
        x, y, z = 1., 2., 3.
        R, z, phi = cart_to_cyl(x, y, z)
        new_x, new_y, new_z = cyl_to_cart(R, z, phi)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    def test_cart_to_sph(self):
        x, y, z = 1., 2., 3.
        phi, theta, r = cart_to_sph(x, y, z)
        new_x, new_y, new_z = sph_to_cart(phi, theta, r)
        self.assertTrue((round(new_x, prec) == round(x, prec)) and (round(new_y, prec) == round(y, prec)) and (round(new_z, prec) == round(z, prec)))

    #Second order coordinate (position) transformations

    def test_cyl_to_gal(self):
        R, z, phi = 10., 5., 120.
        l, b, r = cyl_to_gal(R, z, phi)
        new_R, new_z, new_phi = gal_to_cyl(l, b, r)
        self.assertTrue((round(new_R, prec) == round(R, prec)) and (round(new_z, prec) == round(z, prec)) and (round(new_phi, prec) == round(phi, prec)))

#===============================================================================
# RUNTIME
#===============================================================================

if __name__ == '__main__':
    unittest.main()
