'''
The contents of this file are used to import data from milkyway .out files,
as well as write data out in a variety of formats
'''

#===============================================================================
# IMPORTS
#===============================================================================

import numpy as np
import random
from data import Data

import flags
import glob

#===============================================================================
# FUNCTIONS
#===============================================================================

#removes the header of a milkyway ".out" file
#WARNING: if the file already had a header (or any starting lines) removed,
#   this will delete data from your file and return junk for the COM's
#Returns the center of mass and center of momentum from the header as well
def removeHeader(f):
    #f (open file): the milkyway ".out" file

    comass = []
    comom = []

    #first 4 lines are junk
    for i in range(0, 4):
        f.readline()

    #next line has COM info
    line = f.readline()
    line = line.split(',')
    line[0] = line[0].strip('centerOfMass = ')
    line[3] = line[3].strip('centerOfMomentum = ')
    comass = [float(line[0]), float(line[1]), float(line[2])]
    comom = [float(line[3]), float(line[4]), float(line[5])]

    return comass, comom

#parses a milkyway ".out" file and returns a Data class structure
def readOutput(f, subsample=1.0):
    #f (str): the path of the milkyway ".out" file
    #subsample (float): the percentage [0.0, 1.0] of the data to use
    #   This is done randomly: for more specific subsampling, use the built-in
    #   cutFirstN etc. in the Data class, or build your own and
    #   push it to github

    if flags.progressBars:
        flen = glob.fileLen(f)

    f = open(f, 'r')

    #remove the header, get relevant info from header
    comass, comom = removeHeader(f)

    #store the data here temporarily
    #indexed this way to avoid the 'ignore' column
    array_dict = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[], 8:[], 9:[], 10:[], 11:[], 12:[]}

    if flags.verbose:
        print('Reading in data...')
    #place all the data from the file into the dictionary
    if flags.progressBars:
        j = 0
    for line in f:
        m = random.random()
        if m <= subsample:
            line = line.strip().split(',')
            i = 1
            while i < len(line):
                array_dict[i].append(float(line[i]))
                i += 1
            if flags.progressBars:
                j += 1
                glob.progressBar(j, flen)

    #return the data class using the array dictionary we built
    if flags.verbose:
        print(str(len(array_dict[1])) + ' objects read in')
    return Data(array_dict[1], array_dict[2], array_dict[3], array_dict[4], array_dict[5], array_dict[6], array_dict[7], array_dict[8], array_dict[9], array_dict[10], array_dict[11], array_dict[12], centerOfMass=comass, centerOfMomentum=comom)
