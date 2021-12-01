'''
This file is a place to store flags that affect the operation of mwahpy
'''

#===============================================================================
# FUNCTIONS
#===============================================================================

#checks if the user is running mwahpy in a jupyter notebook
#borrowed from https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
#(thanks Gustavo Bezerra!)
def isnotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return 1   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return 0  # Terminal running IPython
        else:
            return 0  # Other type (?)
    except NameError:
        return 0      # Probably standard Python interpreter

#===============================================================================
# CONSTANTS
#===============================================================================

verbose=1 #If on, functions are a lot noisier
#will slow things down slightly, but it may be helpful for debugging

progress_bars = int(not(isnotebook())) #Display progress bars, on by default
#slows down a program somewhat, but is useful for making sure that
#a long program is running, as well as for debugging

#automatically turns progress bars off if the user is running a jupyter notebook,
#since the current implementation of progres bars clog up the client
#and lag out jupyter notebook

auto_update=1 #update all Timestep objects immediately when necessary
#this will keep things like center of mass, center of momentum, etc.
#accurate for the Timestep class, but will reduce performance somewhat.
#Should only be turned off if you REALLY know what you're doing, and you're sure
#that not updating Timestep objects won't impact what you're doing.

#Timestep classes can still be updated manually when necessary with <Timestep>.update()
