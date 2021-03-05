'''
The contents of this file are used to import data from milkyway .out files,
as well as write data out in a variety of formats
'''

#===============================================================================
# IMPORTS
#===============================================================================

import numpy as np
import sys
from pathlib import Path

import flags
import mwahpy_glob
from timestep import Timestep
from nbody import Nbody

#===============================================================================
# FUNCTIONS
#===============================================================================

#-------------------------------------------------------------------------------
# INPUT
#-------------------------------------------------------------------------------

def removeHeader(f):
    #f (open file): the milkyway ".out" file

    comass = []
    comom = []

    #first 3 lines are junk
    for i in range(0, 3):
        f.readline()

    #next line has COM info
    line = f.readline()
    line = line.split(',')
    line[0] = line[0].strip('centerOfMass = ')
    line[3] = line[3].strip('centerOfMomentum = ')
    comass = [float(line[0]), float(line[1]), float(line[2])]
    comom = [float(line[3]), float(line[4]), float(line[5])]

    return comass, comom

#parses a milkyway ".out" file and returns a Timestep class structure
def readOutput(f):
    #f (str): the path of the milkyway ".out" file

    if flags.progressBars:
        flen = mwahpy_glob.fileLen(f)
    if flags.verbose:
        print('\nReading in data from ' + str(f) + '...')

    f = open(f, 'r')

    #remove the header, get relevant info from header
    comass, comom = removeHeader(f)

    #store the data here temporarily
    #indexed this way to avoid the 'ignore' column
    array_dict = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[]}

    #place all the data from the file into the dictionary
    if flags.progressBars:
        j = 0
    for line in f:
        line = line.strip().split(',')
        i = 1
        while i < len(line) - 1: #this grabs l, b, r data even though that is calculated from x, y, z in the Timestep class implementation
                                 #  this is mostly just for simplicity when reading in, although I guess it probably slows down the code somewhat
            array_dict[i].append(float(line[i]))
            i += 1
        if flags.progressBars:
            j += 1
            mwahpy_glob.progressBar(j, flen)

    #return the Timestep class using the array dictionary we built
    if flags.verbose:
        print('\n'+ str(len(array_dict[1])) + ' objects read in')
        sys.stdout.write('\rConverting data...')
    d = Timestep(id_val=array_dict[1], x=array_dict[2], y=array_dict[3], z=array_dict[4], vx=array_dict[8], vy=array_dict[9], vz=array_dict[10], mass=array_dict[11], centerOfMass=comass, centerOfMomentum=comom)
    if flags.verbose:
        sys.stdout.write('done\n')

    f.close()

    return d

#parses a milkyway ".in" file and returns a Timestep class structure
def readInput(f):
    #f (str): the path of the milkyway ".in" file

    if flags.progressBars:
        flen = mwahpy_glob.fileLen(f)
    if flags.verbose:
        print('\nReading in data from ' + str(f) + '...')

    f = open(f, 'r')

    #remove the header
    f.readline()

    #store the data here temporarily
    #indexed this way to avoid the 'ignore' column
    array_dict = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[]}

    #place all the data from the file into the dictionary
    if flags.progressBars:
        j = 0
    for line in f:
        line = line.strip().split('\t')
        i = 1
        while i < len(line):
            array_dict[i].append(float(line[i]))
            i += 1
        if flags.progressBars:
            j += 1
            mwahpy_glob.progressBar(j, flen)

    #return the Timestep class using the array dictionary we built
    if flags.verbose:
        print('\n'+ str(len(array_dict[1])) + ' objects read in')
        sys.stdout.write('\rConverting data...')
    d = Timestep(id_val=array_dict[1], x=array_dict[2], y=array_dict[3], z=array_dict[4], vx=array_dict[5], vy=array_dict[6], vz=array_dict[7], mass=array_dict[8], centerOfMass=[0,0,0], centerOfMomentum=[0,0,0])
    if flags.verbose:
        sys.stdout.write('done\n')

    f.close()

    d.update(force=True) #this generates Comass and Comomentum, as well as the other bits
    #that you expect a Timestep to have initially when read in

    return d

#TODO: add nested progress bars
def readFolder(f, ts_scale=None):
    #f: the path to the folder that you want to create an Nbody structure out of
    #ts_scale: the scale of a single timestep in the Nbody sim

    if flags.verbose:
        print('\nReading in data from directory ' + str(f) + '...')

    n = Nbody(ts_scale=ts_scale) #create the Nbody instance

    #iterate through the folder
    for i in Path(f).iterdir():

        t = readOutput(str(i))
        time = int(str(i).split('/')[-1])
        n[time] = t #assign each timestep to the correct place in the Nbody


    return n

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------

#parses a Timestep class object and outputs a file that can be read into a
#MilkyWay@home N-body simulation as the 'manual bodies' parameter
def makeNbodyInput(t, f):
    #t (Timestep): the Timestep object that will be printed out
    #f (str): the path of the file that will be printed to
    if flags.verbose:
        print('Writing Timestep as N-body input to '+f+'...')

    f = open(f, 'w')
    f.write('#ignore\tid\tx\ty\tz\tvx\tvy\tvz\tm')

    i = 0

    #recenter on center of mass for ease of inserting a dwarf
    xnew = t.x - t.centerOfMass[0]
    ynew = t.y - t.centerOfMass[1]
    znew = t.z - t.centerOfMass[2]

    while i < len(t):
        f.write('\n1\t'+str(int(t.id[i]))+'\t'+str(xnew[i])+'\t'+str(ynew[i])+'\t'+str(znew[i])+'\t'+\
                str(t.vx[i])+'\t'+str(t.vy[i])+'\t'+str(t.vz[i])+'\t'+str(t.mass[i]))
        if flags.progressBars:
            mwahpy_glob.progressBar(i, len(t))
        i += 1

    if flags.verbose:
        print('\ndone')

#prints out a '.csv' file of the Timestep class structure
def makeCSV(t, f_name):
    #t: the Timestep object being written out
    #f: the path for the output csv file

    f = open(f_name, 'w')

    #make the header
    #TODO: Print out COM's
    if flags.verbose:
        print('Writing header...')
    header = ''
    for key in t:
        header += (key + ',')
    header += '\n'
    f.write(header)

    #iterate through the data and print each line
    i = 0
    if flags.verbose:
        print('Printing data...')
    while i < len(self): #i avoided zip() here because it's messy to zip
    #                        like a billion different arrays, although zip()
    #                        is "better programming"
        if flags.progressBars:
            mwahpy_glob.progressBar(i, len(t))
        line = ''
        for key in t:
            line += (str(t[key][i]) + ',')
        line += '\n'
        f.write(line)
        i += 1

    print('Timestep output to ' + f_name)
