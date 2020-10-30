 MILKYWAY@HOME PYTHON PACKAGE (MWAHPY)
========================================

Copyright Tom Donlon, 2020 RPI

github user: thomasdonlon

Requires Python v.>3.2.3

-----------------------------------------

MilkyWay@home is a computational astrophysics project located at RPI and
run by Dr. Heidi Newberg. The project contains two main components,
(i) the Separation Application and (ii) the N-body Application. I have
worked with both, and over a few years have developed several useful bits of
python code to supplement analysis and usage of MilkyWay@home. This code has
been cleaned and organized here, so that others can benefit from the work I did
and not have to spend time rewriting thing's I've already done.

The purpose of the mwahpy package is to provide a collection of useful tools
for people working with/on the MilkyWay@home project at RPI. Ideally, the
package will contain tools for both Separation and N-body, although that might
not initially be the case.

In practice it would be good for this package to be updated along with the
MilkyWay@home project and/or python updates after I'm gone. No clue if that's
going to happen, but I plan on maintaining support at least through the spring
of 2022 (while I'm in graduate school).

Issues with the package can be directed to my github profile or on the
MilkyWay@home forums at milkyway.cs.rpi.edu. Improvements or additions are
welcome, just send a pull request to the mwahpy github repository.


INSTALLATION
========================================

FOR USERS:

1. Open your terminal, and run

> python3 -m pip install mwahpy

2. Insert import statements for the subpackages that you want to use in your .py files:

> import mwahpy.{mwahpy subpackage}

> import mwahpy.{other mwahpy subpackage}

> ...

3. Do science

FOR DEVELOPERS:

1. Download the mwahpy github repository

2. Place the files in the repository in your sys-path directory
(this can be found by running "python3 -m site" in your terminal
and is given after "USER-SITE")
(also see https://stackoverflow.com/questions/122327/how-do-i-find-the-location-of-my-python-site-packages-directory)

3. Insert import statements for the subpackages that you want to use:

> import mwahpy.{mwahpy subpackage}

> import mwahpy.{other mwahpy subpackage}

> ...

in your code

4. Do science

Alternatively,

1. Download the mwahpy github repository

2. At the beginning of your .py file, put the path to the /mwahpy/mwahpy directory,
and then import the relevant mwahpy .py files that you need to access

> sys.path.insert(1, '.../mwahpy/mwahpy')

> import {mwahpy subpackage}

> import {other mwahpy subpackage}

> ...

3. Do science

CONTENTS
========================================

A non-exhaustive list of contents of the package is given below:

 - easy importing of data from N-body output
 - easy manipulation of data after reading in
 - A variety of coordinate transformations for typical coordinate systems used in Galactic astronomy
 - easy visualization of the imported data through plotting functionality

TODO
========================================

 - Finish refactoring coords.py
 - Finish unit testing of coordinate transformations
 - Refactor and polish orbit_fitting_gc.py
 - Write API
