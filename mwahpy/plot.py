'''
The contents of this file are focused on the plotting of the Data structure
in various projections and formats

These functions do this in whatever matplotlib instance you've got going on,
unless you toggle <show>
e.g. consider that below each function, I've added
    #show: if True, opens a window and shows the plot. Otherwise, adds it to
    #whatever matplotlib figure instance you have open
'''

#===============================================================================
# IMPORTS
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt
from data import Data
import data

#===============================================================================
# FUNCTIONS
#===============================================================================

#generates a scatter plot of your Data object
def scatter(d, x, y, show=False, s=5.0, color='k', marker='o', **kwargs):
    #d: the Data object to plot
    #x: the x-axis parameter
    #y: the y-axis parameter

    plt.scatter(d[x], d[y], s=s, c=color, marker=marker, **kwargs)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

#sticks a big fat red dot wherever the specific star(s) is, given an id(s)
#if vx and vy are specified, an arrow is drawn
def traceParticle(d, id, x, y, vx=None, vy=None, vscale=0.02, show=False, s=50., color='r', marker='o', **kwargs):
    #d (Data): the Data object that your particle is plotted over
    #id (int, can be array-like): the (list of) id(s) for your particle(s)
    #x (str): the x-axis parameter
    #y (str): the y-axis parameter
    #vx (str, optional): the x-axis arrow parameter
    #vy (str, optional): the y-axis arrow parameter
    #vscale (float, >0, optional): scales the arrow size
    #TODO: more customization of the arrow, probably

    #make sure id is always a numpy array of id values
    if not hasattr(id, "__len__"): #if true, id is not array-like:
        id = [id]
    else:
        if type(id)!=type(np.array([])):
            id = np.array(id)

    #get the indices of the particles from the id's
    #gotta be a 1-line way to do this with numpy, but I couldn't figure it out
    n = []
    for i in id:
        n.append(d['id'].index(i))

    #plot the particles
    for i in n:
        plt.scatter(d[x][i], d[y][i], s=s, c=color, marker=marker, **kwargs)
        if vx and vy:
            #TODO: figure out how to pass kwargs to the arrow and scatter separately
            plt.arrow(d[x][i], d[y][i], d[vx][i]*vscale, d[vy][i]*vscale, color=color, head_width=1)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

#plots a histogram of the Data
#see np.hist for usage
def hist(d, x, show=False, *args, **kwargs):
    #d (Data): the data object being plotted
    #x (str): the axis parameter

    h = plt.hist(d[x], range=range, bins=bins, *args, **kwargs)

    if show:
        plt.xlabel(x)
        plt.show()

    return h

def hist2d(d, x, y, show=false, *args, **kwargs):
    #d (Data): the data object being plotted
    #x (str): the x-axis parameter
    #y (str): the y-axis parameter

    h = plt.hist2d(d[x], d[y], *args, **kwargs)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

    return h
