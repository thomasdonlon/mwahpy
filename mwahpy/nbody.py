'''
The contents of this file are focused on the Nbody class, which is used for storage of
imported data from N-body output files. Specifically, an Nbody instance is a collection
of Timestep classes, with overarching capability for plotting entire simulations and
cutting simulations as a whole, instead of manually going through each timestep
'''

#===============================================================================
# IMPORTS
#===============================================================================

#external imports
import numpy as np
import random
import unittest

#mwahpy imports
import mwahpy_glob
import flags
import output_handler
import timestep

#TODO: make a Nbody class that reads in from a folder full of timesteps, and acts
#   as a collection of Timestep instances. Each Timestep instance is one of the timesteps,
#   and can be called as that timestep.

#TODO: write a function that reads in Nbody data from a folder of Nbody .out files

#===============================================================================
# NBODY CLASS
#===============================================================================

class Nbody():

    def __init__(self, ts={}, ts_scale=None):

        self.ts = ts #dictionary of Timestep instances with the numerical timestep
            #(the given filename, e.g. 99, 2099, etc.) as the key

        if ts_scale != None:
            self.scaleTimes(ts_len)

    #---------------------------------------------------------------------------
    # INDEXING
    #---------------------------------------------------------------------------

    def __getitem__(self, i):
        return self.ts[i]

    def __setitem__(self, i, val):
        val.time = i
        val.nbody = self
        self.ts[i] = val

    #---------------------------------------------------------------------------
    # ITERATOR CONTROL
    #---------------------------------------------------------------------------

    def __iter__(self):
        self.sortedts = sorted(self.ts.items()) #implementing this here prevents the program from
            #having to sort the dictionary every iteration
        self.index = 0
        return self

    def __next__(self):
        if self.index == len(self.sortedts):
            raise StopIteration
        else:
            self.index += 1
            return self.sortedts[self.index - 1][1] #sortedts is actual a list of tuples,
                                                #but we don't care about the keys

    #---------------------------------------------------------------------------
    # METHODS
    #---------------------------------------------------------------------------

    #adds the actual times of each timestep in the specified unit to a list
    #this can be provided in the initialization of the Nbody class, or adjusted later
    def scaleTimes(self, l):
        #l (float): the length of each timestep

        #save the times in the nbody structure
        self.times = [l * x for x in sorted(self.ts.keys())]

        new_ts = {} #have to build a new dict in order to not delete information
                   #while renaming the keys of the old dict

        #save the time in each Timestep so that it can be accessed from there
        for key in sorted(self.ts.keys()):
            self.ts[key].time = l*key
            new_ts[l*key] = self.ts[key] #build a new dict

        self.ts = new_ts

#===============================================================================
# FUNCTIONS INVOLVING NBODY CLASSES
#===============================================================================


#===============================================================================
# UNIT TESTING
#===============================================================================


#===============================================================================
# RUNTIME
#===============================================================================
