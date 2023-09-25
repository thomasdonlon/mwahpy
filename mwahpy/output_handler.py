'''
The contents of this file are used to import data from milkyway .out files,
as well as write data out in a variety of formats
'''

#===============================================================================
# IMPORTS
#===============================================================================

import numpy as np
from pathlib import Path
import pandas as pd

from .flags import verbose, progress_bars
from .mwahpy_glob import file_len, progress_bar
from .timestep import Timestep
from .nbody import Nbody
from .histogram import histData

#===============================================================================
# FUNCTIONS
#===============================================================================

#-------------------------------------------------------------------------------
# INPUT
#-------------------------------------------------------------------------------

def remove_header(f):
    #f (open file): the milkyway ".out" file

    has_bodies_tag = False #tracks if there are leading and trailing <bodies>/</bodies> tags in the .out file

    comass = []
    comom = []

    #first few lines are junk
    line = f.readline()
    if line.strip() == '<bodies>': #exactly how many lines to skip depends on if we have the tags or not
        n_skip = 3
        has_bodies_tag = True
    else:
        n_skip = 2
    for i in range(n_skip):
        f.readline()

    #next line has COM info
    line = f.readline()
    line = line.split(',')
    line[0] = line[0].strip('centerOfMass = ')
    line[3] = line[3].strip('centerOfMomentum = ')
    comass = [float(line[0]), float(line[1]), float(line[2])]
    comom = [float(line[3]), float(line[4]), float(line[5])]

    return comass, comom, has_bodies_tag

#parses a milkyway ".out" file and returns a Timestep class structure
def read_output(f, start=None, stop=None):
    #f (str): the path of the milkyway ".out" file
    #start (optional): the line to actually start reading in data
    #stop (optional): the line to stop reading in data

    flen = 0 #allows this to be used outside of the next if statement's scope
    if progress_bars:
        flen = file_len(f)
        #properly adjust for length change based on start and stop points
        if start is not None:
            flen -= start
        if stop is not None:
            flen -= stop
    if verbose:
        print(f"\nReading in data from {f}...")

    f = open(f, 'r')

    #remove the header, get relevant info from header
    comass, comom, has_bodies_tag = remove_header(f)

    if has_bodies_tag: #have to account for the </bodies> tag at the end of the file
        if stop is not None: #if stop is given, then you assume it cuts off before the end of the file
            stop = flen - 1
            flen -= 1

    #store the data here temporarily
    array_dict = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[]}

    #place all the data from the file into the dictionary
    j = 0 #line num tracker

    #if start is specified,
    #skip lines until we're at the starting line
    if start is not None:
        while j < start:
            f.readline()
            j += 1

    j = 0 #reset line counter for progress bar
    for line in f:

        line = line.strip().split(',')
        i = 0
        while i < len(line) - 1: #this grabs l, b, r data even though that is calculated from x, y, z in the Timestep class implementation
                                 #it's mostly just for simplicity when reading in, although I guess it probably slows down the code somewhat
            array_dict[i].append(float(line[i]))
            i += 1
        j += 1

        if progress_bars:
            progress_bar(j, flen)

        #stop reading in lines if we're past the stop setting
        if stop is not None:
            if j >= stop:
                break
        else:
            continue  # only executed if the inner loop did NOT break
        break #only executed if the inner loop DID break

    #return the Timestep class using the array dictionary we built
    if verbose:
        print(f"\n{len(array_dict[1])} objects read in")
        print('\rConverting data...', end='')
    d = Timestep(typ=array_dict[0], id_val=array_dict[1], x=array_dict[2], y=array_dict[3], z=array_dict[4], vx=array_dict[8], vy=array_dict[9], vz=array_dict[10], mass=array_dict[11], center_of_mass=comass, center_of_momentum=comom)
    if verbose:
        print('done')

    f.close()

    return d

#parses a milkyway ".in" file and returns a Timestep class structure
def read_input(f):
    #f (str): the path of the milkyway ".in" file

    if progress_bars:
        flen = file_len(f)
    if verbose:
        print('\nReading in data from ' + str(f) + '...')

    f = open(f, 'r')

    #remove the header
    f.readline()

    #store the data here temporarily
    array_dict = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[]}

    #place all the data from the file into the dictionary
    if progress_bars:
        j = 0
    for line in f:
        line = line.strip().split('\t')
        i = 0
        while i < len(line):
            array_dict[i].append(float(line[i]))
            i += 1
        if progress_bars:
            j += 1
            progress_bar(j, flen)

    #return the Timestep class using the array dictionary we built
    if verbose:
        print('\n'+ str(len(array_dict[1])) + ' objects read in')
        print('\rConverting data...', end='')
    d = Timestep(typ=array_dict[0], id_val=array_dict[1], x=array_dict[2], y=array_dict[3], z=array_dict[4], vx=array_dict[5], vy=array_dict[6], vz=array_dict[7], mass=array_dict[8], center_of_mass=[0,0,0], center_of_momentum=[0,0,0])
    if verbose:
        print('done')

    f.close()

    d.update(force=True) #this generates Comass and Comomentum, as well as the other bits
    #that you expect a Timestep to have initially when read in

    return d

def read_folder(f, ts_scale=None):
    #f: the path to the folder that you want to create an Nbody structure out of
    #ts_scale: the scale of a single timestep in the Nbody sim

    if verbose:
        print('\nReading in data from directory ' + str(f) + '...')

    n = Nbody(ts_scale=ts_scale) #create the Nbody instance

    #iterate through the folder
    for i in Path(f).iterdir():

        t = read_output(str(i))
        time = int(str(i).split('/')[-1])
        n[time] = t #assign each timestep to the correct place in the Nbody

    return n

#Reads in histogram file path as a string and creates histData class
def read_histogram(hist_file_path):
    
    #initialising histData class instance
    hist = histData()

    if (verbose):
        print('Creating python class for given histogram data...')

    with open(str(hist_file_path)) as f:
        line = f.readline()
        lineNumber = 1
        while (line):
            line = line.strip()
            line = line.split(' ')
            if ('' in line):
                line = list(filter(None, line))
            #print(line)
            if (len(line)> 1):
                #Date histogram was generated
                if (line[1]=='Generated'):
                    generationDate = ''
                    for i in range(len(line) - 2):
                        generationDate = generationDate + line[i+2]+ ' '
                    hist.generationDate = generationDate
                    #print(hist.generationDate)
                #Euler Angles in degrees
                if (line[1]=='(phi,'):
                    eulerPhi = float(line[5].strip('(').strip(','))
                    eulerTheta = float(line[6].strip(','))
                    eulerPsi = float(line[7].strip(')').strip(','))
                    hist.eulerAngles = [eulerPhi, eulerTheta, eulerPsi]
                    #print(hist.eulerAngles)
                #Lambda Range
                if (line[1]=='(lambdaStart,'):
                    lambdaStart = float(line[5].strip('(').strip(','))
                    lambdaCenter = float(line[6].strip(','))
                    lambdaEnd = float(line[7].strip(')').strip(','))
                    hist.lambdaRange = [lambdaStart, lambdaCenter, lambdaEnd]
                    #print(hist.lambdaRange)
                #Beta Range
                if (line[1]=='(betaStart,'):
                    betaStart = float(line[5].strip('(').strip(','))
                    betaCenter = float(line[6].strip(','))
                    betaEnd = float(line[7].strip(')').strip(','))
                    hist.betaRange = [betaStart, betaCenter, betaEnd]
                    #print(hist.betaRange)
                #Lambda Bin Size
                if (line[1]=='Lambda'):
                    hist.lambdaBinSize = float(line[5])
                    #print(hist.lambdaBinSize)
                #Beta Bin Size
                if (line[1]=='Beta'):
                    hist.betaBinSize = float(line[5])
                    #print(hist.betaBinSize)
                #Nbody (total numbers of body in Nbody simulation)
                if (line[1]=='Nbody'):
                    hist.nbody = int(line[3])
                    #print(hist.nbody)
                #Evolve backward time and evolve foward time in Gyr
                if (line[1]=='Evolve'):
                    if (line[2]=='backward'):
                        hist.evolveBackwardTime = float(line[5])
                        #print(hist.evolveBackwardTime)
                    elif (line[2]=='forward'):
                        hist.evolveFowardTime = float(line[5])
                        #print(hist.evolveFowardTime)
                #Best Likelihood (absolute value of the likelihood score)
                if (line[1]=='Best'):
                    hist.bestLikelihood = float(line[4])
                    #print(hist.bestLikelihood)
                #Timestep in Gyr
                if (line[1]=='Timestep'):
                    hist.timestep = float(line[3])
                    #print(hist.timestep)
                #Sun GC Distance in kpc
                if (line[1]=='Sun'):
                    hist.sunGCDist = float(line[5])
                    #print(hist.sunGCDist)
                #Criterion (method used to caclulate nbody interections)
                if (line[1]=='Criterion'):
                    hist.criterion = line[3]
                    #print(hist.criterion)
                #Theta (ratio used by treecode to determine whne to calculate using a single body vs using a collection of bodies)
                if (line[1]=='Theta'): 
                    hist.theta = float(line[3])
                    #print(hist.theta)
                #Quadrupole Moments (used by tree code)
                if (line[1]=='Quadrupole'):
                    hist.quadrupoleMoments = line[4]
                    #print(hist.quadrupoleMoments)
                #Eps (softening length in kpc)
                if (line[1]=='Eps'):
                    hist.eps = float(line[3])
                    #print(hist.eps)
                #Bar time error in Gyr and bar angle error in radians
                if (line[1]=='Bar'):
                    if (line[2]=='Time'):
                        hist.barTimeError = float(line[5])
                        #print(hist.barTimeError)
                    if (line[2]=='Angle'):
                        hist.barAngleError = float(line[5])
                        #print(hist.barAngleError)
                #Spherical Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Spherical:'):
                    hist.spherical = line[2]
                    #print(hist.spherical)
                #Primary Disk Potential Type  (just saves the type of potential and not the mass and such)
                if (line[1]=='Primary'):
                    hist.primaryDisk = line[3]
                    #print(hist.primaryDisk)
                #Secondary Disk Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Secondary'):
                    hist.secondaryDisk = line[3]
                    #print(hist.secondaryDisk)
                #Halo Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Halo:'):
                    hist.halo = line [2]
                    #print(hist.halo)
                #Column headers
                if (line[1]=='UseBin,'):
                    #print(len(line))
                    if (len(line) > 22):
                        hist.orbitFitting = True
                        if (verbose):
                            print('Data includes orbit fitting parameters')
                    else:
                        hist.orbitFitting = False
                        if (verbose):
                            print('Data does not include orbit fitting parameters')
                    columnHeaders = ''
                    for i in range(len(line) - 1):
                        columnHeaders += line[i + 1] + ' '
                    hist.columnHeaders = columnHeaders
                    #print(hist.columnHeaders)
                #n (total number of bodies within ranges of histogram)
                if (line[0]=='n'):
                    hist.n = int(line[2])
                    #print(hist.n)
                #Mass per particle (structure mass units)
                if (line[0]=='massPerParticle'):
                    hist.massPerParticle = float(line[2])
                    #print(hist.massPerParticle)
                if (line[0]=='totalSimulated'):
                    hist.totalSimulated = int(line[2])
                    #print(hist.totalSimulated)
                #Number of lambda bins
                if (line[0]=='lambdaBins'):
                    hist.lambdaBins = int(line[2])
                    #print(hist.lambdaBins)
                #number of Beta Bins
                if (line[0]=='betaBins'):
                    hist.betaBins = int(line[2])
                    #print(hist.betaBins)
                #Getting Data (Not sure if useBin can have any other value other than 1)
                if (line[0]=='1'):
                    if (hist.orbitFitting): #checks to see how many columns there are since different in different versions
                        try:
                            hist.useBin.append(float(line[0]))
                            hist.lamb.append(float(line[1]))
                            hist.beta.append(float(line[2]))
                            hist.normalizedCounts.append(float(line[3]))
                            hist.countError.append(float(line[4]))
                            hist.betaDispersion.append(float(line[5]))
                            hist.betaDispersionError.append(float(line[6]))
                            hist.losVelocityDispersion.append(float(line[7]))
                            hist.velocityDispersionError.append(float(line[8]))
                            hist.losVelocity.append(float(line[9]))
                            hist.losVelocityError.append(float(line[1]))
                            hist.betaAverage.append(float(line[11]))
                            hist.betaAverageError.append(float(line[12]))
                            hist.distance.append(float(line[13]))
                            hist.distanceError.append(float(line[12]))
                        except ValueError:
                            print('Invalid histagram data entry: Non-numerical value in histogram data on line ' + str(lineNumber))
                    else:
                        try:
                            hist.useBin.append(float(line[0]))
                            hist.lamb.append(float(line[1]))
                            hist.beta.append(float(line[2]))
                            hist.normalizedCounts.append(float(line[3]))
                            hist.countError.append(float(line[4]))
                            hist.betaDispersion.append(float(line[5]))
                            hist.betaDispersionError.append(float(line[6]))
                            hist.losVelocityDispersion.append(float(line[7]))
                            hist.velocityDispersionError.append(float(line[8]))
                        except ValueError:
                            print('Invalid histagram data entry: Non-numerical value in histogram data on line ' + str(lineNumber))
                
            line = f.readline()
            lineNumber += 1
    f.close()

    if (verbose):
        print('Done')
    return hist

def hist_to_df(hist):
    if (verbose):
        print('Creating data frame for given histogram data...')
    #Making Pandas Data Frame
    df = pd.DataFrame()
    if (verbose):
        if (hist.orbitFitting):
            print('Orbit fitting data included')
        else:
            print('Orbit fitting data not included')
    df['Use Bin'] = hist.useBin
    df['Lambda'] = hist.lamb
    df['Beta'] = hist.beta 
    df['Normalized Counts'] = hist.normalizedCounts
    df['Count Error'] = hist.countError
    df['Beta Dispersion'] = hist.betaDispersion
    df['Beta Dispersion Error'] = hist.betaDispersionError
    df['LOS Velocity Dispersion'] = hist.losVelocityDispersion
    df['Velocity Dispersion Error'] = hist.velocityDispersionError
    if (hist.orbitFitting):
        df['LOS Velocity'] = hist.losVelocity
        df['LOS Velocity Error'] = hist.losVelocityError
        df['Beta Average'] = hist.betaAverage
        df['Beta Average Error'] = hist.betaAverageError
        df['Distance'] = hist.distance
        df['Distance Error'] = hist.distanceError
    if (verbose):
        print('Done')
    return df

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------

#parses a Timestep class object and outputs a file that can be read into a
#MilkyWay@home N-body simulation as the 'manual bodies' parameter
def make_nbody_input(t, f, recenter=True):
    #t (Timestep): the Timestep object that will be printed out
    #f (str): the path of the file that will be printed to
    if verbose:
        print('Writing Timestep as N-body input to '+f+'...')

    if recenter:
        t.recenter()

    f = open(f, 'w')
    f.write('#ignore\tid\tx\ty\tz\tvx\tvy\tvz\tm')

    i = 0

    while i < len(t):
        f.write('\n{}\t{}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{} '.format(1, t.id[i], t.x[i], t.y[i], t.z[i], t.vx[i], t.vy[i], t.vz[i], t.mass[i]))
        #f.write('\n'+str(1)+' '+str(t.id[i])+' '+str(t.x[i])+' '+str(t.y[i])+' '+str(t.z[i])+' '+\
        #        str(t.vx[i])+' '+str(t.vy[i])+' '+str(t.vz[i])+' '+str(t.mass[i]))
        if progress_bars:
            progress_bar(i, len(t))
        i += 1

    f.write('\n')

    if verbose:
        print('\ndone')

#prints out a '.csv' file of the Timestep class structure
def make_csv(t, f_name):
    #t: the Timestep object being written out
    #f: the path for the output csv file

    f = open(f_name, 'w')

    #make the header
    if verbose:
        print('Writing header...')
    header = ''
    for key in t:
        header += (key + ',')
    header += '\n'
    f.write(header)

    #iterate through the data and print each line
    i = 0
    if verbose:
        print('Printing data...')
    while i < len(t): #i avoided zip() here because it's messy to zip
    #                        like a billion different arrays, although zip()
    #                        is "better programming"
        if progress_bars:
            progress_bar(i, len(t))
        line = ''
        for key in t:
            line += (str(t[key][i]) + ',')
        line += '\n'
        f.write(line)
        i += 1

    print('Timestep output to ' + f_name)
