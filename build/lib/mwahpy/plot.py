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

#===============================================================================
# CONSTANTS & FLAGS
#===============================================================================

#used as the default background for fancy plots
default_fancy_background_url = 'https://upload.wikimedia.org/wikipedia/commons/thumb/4/4f/NGC_5170_HST_9766_R814B435.png/2560px-NGC_5170_HST_9766_R814B435.png'

#the licensing warning given when the user uses fancy plot
licensing_warning = "WARNING: The fancy_plot() routine uses an image that has been licensed by Rensselaer Polytechnic Institute (RPI). If you are not affiliated with RPI, then public usage of the default background image is not allowed. Please provide another URL for the background image with the `bg_url` keyword."

#===============================================================================
# FUNCTIONS
#===============================================================================

#generates a scatter plot of your Timestep object
def scatter(t, x, y, cbar=None, show=False, s=5.0, color=None, marker='o', **kwargs):
    #t: the Timestep object to plot
    #x: the x-axis parameter
    #y: the y-axis parameter

    if (cbar is not None) and (color is not None):
        print(UserWarning("Provided both 'cbar' and 'color' argument. Using cbar argument for this plot."))

    if color is None:
        color = 'k' #default set to black, but this can't be assigned in the def line because we need to check if a user set it

    if cbar is not None:
        plt.scatter(t[x], t[y], s=s, c=t[cbar], marker=marker, **kwargs)
        plt.colorbar()
    else:
        plt.scatter(t[x], t[y], s=s, c=color, marker=marker, **kwargs)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

#sticks a big fat red dot wherever the specific star(s) is, given an id(s)
#if vx and vy are specified, an arrow is drawn
def trace_particle(t, id, x, y, vx=None, vy=None, vscale=0.02, show=False, s=50., color='r', marker='o', **kwargs):
    #t (Timestep): the Timestep object that your particle is plotted over
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
        n.append(t['id'].index(i))

    #plot the particles
    for i in n:
        plt.scatter(t[x][i], t[y][i], s=s, c=color, marker=marker, **kwargs)
        if vx and vy:
            #TODO: figure out how to pass kwargs to the arrow and scatter separately
            plt.arrow(t[x][i], t[y][i], t[vx][i]*vscale, t[vy][i]*vscale, color=color, head_width=1)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

#plots a histogram of the Timestep
#see np.hist for usage
def hist(t, x, show=False, range=None, bins=10, *args, **kwargs):
    #t (Timestep): the Timestep object being plotted
    #x (str): the axis parameter

    if range is None: #default range
        range = (min(t[x]), max(t[x]))

    h = plt.hist(t[x], range=range, bins=bins, *args, **kwargs)

    if show:
        plt.xlabel(x)
        plt.show()

    return h

def hist2d(t, x, y, show=False, *args, **kwargs):
    #t (Timestep): the Timestep object being plotted
    #x (str): the x-axis parameter
    #y (str): the y-axis parameter

    h = plt.hist2d(t[x], t[y], *args, **kwargs)

    if show:
        plt.xlabel(x)
        plt.ylabel(y)
        plt.show()

    return h

#-------------------------------------------------------------------------------
# fancy plots
#-------------------------------------------------------------------------------

#fancy plot is technically built off of plt.plot, not plt.scatter
#so the kwargs need to be for plot, not scatter.
def fancy_plot(t, ax1, ax2, show=False, ax=None, bg_url=default_fancy_background_url, no_axes=True, extent=(-100, 100, -100, 100), **kwargs):
    #t (Timestep): the Timestep that contains the data
    #ax1 (str): the x-axis, i.e. 'x' or 'z'
    #ax2 (str): the y-axis
    #show (bool): if True, shows the plot
    #ax (pyplot axis object): can pass an axis to place the fancy image in
    #    if None, then this routine just makes its own figure
    #bg_url (str): the url of an image online to use for the background image
    #    (by default this is an image that RPI has licensing for)
    #no_axes (bool): if True, removes the axes labels
    #extent (tuple of floats): (xmin, xmax, ymin, ymax) of plot bounds

    #out: returns the active pyplot axis, so that you can alter the image as needed

    #print the licensing warning if this is the first time that the user has made a fancy plot
    if bg_url == default_fancy_background_url:
        print(UserWarning(licensing_warning))

    # First set up the figure (if we have to)
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    #set current axis
    plt.sca(ax)

    # read the image file in a numpy array
    a = plt.imread(bg_url)
    ax.imshow(a, extent=extent)

    #bottom right text
    if bg_url == default_fancy_background_url:
        ax.text(extent[0] + 0.05*(extent[1] - extent[0]), extent[2] + 0.05*(extent[3] - extent[2]), 'Background image credit: ESA/NASA', color='w',
                fontsize=20, zorder=1000, horizontalalignment='left', verticalalignment='bottom')

    #split up into dark and baryonic matter
    bm_indx = (t.typ == 0)
    dm_indx = (t.typ == 1)

    x_bm = t[ax1][bm_indx]
    y_bm = t[ax2][bm_indx]
    x_dm = t[ax1][dm_indx]
    y_dm = t[ax2][dm_indx]

    #plot baryonic matter
    if len(x_bm) > 0:
        bm_points, = ax.plot(x_bm, y_bm, color=(1, 1, 1, 1), marker='o', linestyle='', ms=1., zorder=100, **kwargs)
        bm_points2, = ax.plot(x_bm, y_bm, color=(0.8, 0.8, 0.8, 0.8), marker='o', linestyle='', ms=2., zorder=99, **kwargs)
        bm_points3, = ax.plot(x_bm, y_bm, color=(0.6, 0.6, 0.6, 0.6), marker='o', linestyle='', ms=4., zorder=98, **kwargs)
        bm_points4, = ax.plot(x_bm, y_bm, color=(0.4, 0.4, 0.4, 0.4), marker='o', linestyle='', ms=8., zorder=97, **kwargs)
        bm_points5, = ax.plot(x_bm, y_bm, color=(0.2, 0.2, 0.2, 0.2), marker='o', linestyle='', ms=16., zorder=96, **kwargs)

    #plot dark matter
    if len(x_dm) > 0:
        size_mult = 1.2
        dm_points, = ax.plot(x_dm, y_dm, color=(238/256,130/256,238/256, 1), marker='o', linestyle='', ms=1.*size_mult, zorder=95, **kwargs)
        dm_points2, = ax.plot(x_dm, y_dm, color=(221/256,160/256,221/256, 0.8), marker='o', linestyle='', ms=2.*size_mult, zorder=94, **kwargs)
        dm_points3, = ax.plot(x_dm, y_dm, color=(186/256,85/256,211/256, 0.6), marker='o', linestyle='', ms=4.*size_mult, zorder=93, **kwargs)
        dm_points4, = ax.plot(x_dm, y_dm, color=(148/256,0/256,211/256, 0.4), marker='o', linestyle='', ms=8.*size_mult, zorder=92, **kwargs)
        dm_points5, = ax.plot(x_dm, y_dm, color=(139/256,0/256,139/256, 0.2), marker='o', linestyle='', ms=16.*size_mult, zorder=91, **kwargs)

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    if no_axes:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    if show:
        plt.show()

    #return axis so you can manipulate it outside this function
    return ax

#-------------------------------------------------------------------------------
