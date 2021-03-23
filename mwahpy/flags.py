'''
This file is a place to store flags that affect the operation of mwahpy
'''

verbose=1 #If on, functions are a lot noisier
#will slow things down slightly, but it may be helpful for debugging

progress_bars=1 #Display progress bars
#slows down a program somewhat, but is useful for making sure that
#a long program is running, as well as for debugging

auto_update=1 #update all Timestep objects immediately when necessary
#this will keep things like center of mass, center of momentum, etc.
#accurate for the Timestep class, but will reduce performance somewhat. Should only be
#turned off if you REALLY know what you're doing, and you're sure that not updating
#Timestep objects won't impact what you're doing.

#Timestep classes can still be updated manually when necessary with <Timestep>.update()
