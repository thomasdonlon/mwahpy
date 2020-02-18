'''
This file is a place to store flags that affect the operation of mwahpy
'''

verbose=1 #If on, functions are a lot noisier
#will certainly slow things down, but it may be helpful for debugging

progressBars=1 #Display progress bars
#slows down a program, but is useful for making sure that
#a long program is running, as well as for debugging

calcEnergy=1 #calculate energies in a Data object
#this has some decent overhead associated with it, so if speed is super
#important to you, you can turn off this flag
