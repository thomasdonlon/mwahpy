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

CONTENTS
========================================

A non-exhaustive list of contents of the package is given below:

 - easy importing of data from N-body output
 - easy manipulation of data after reading in
 - A variety of coordinate transformations for typical coordinate systems used in Galactic astronomy
 - easy visualization of the imported data through plotting functionality
 - a tutorial .pdf document along with documentation for each function in the package

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

1. Clone the mwahpy github repository

2. Make the desired changes in the source code

3. Navigate to the directory where you cloned the repo, and then run

> python3 setup.py develop --user

(note, you will probably need to uninstall any previous versions of mwahpy you had on your machine before running this)

4. To test your changes, insert import statements for the subpackages that you want to use in your .py files as you normally would:

> import mwahpy.{mwahpy subpackage}

> import mwahpy.{other mwahpy subpackage}

> ...

Any time you make a change to the source code, you'll have to repeat step 3

5. Once you are done making changes to the source code, put in a pull request to master

6. Navigate to the directory where you cloned the repo, and then run

> python3 setup.py develop --uninstall

> pip3 install mwahpy

Your changes will not be available in the main mwahpy build until a new release comes out.

TODO
========================================

MAJOR:
 - Expand plot capabilities
 - Finish refactoring coords.py
 - Finish unit testing of coordinate transformations
 - Refactor and polish orbit_fitting_gc.py

MINOR:
- Apparently there's a MW@h readout method with a different number of values or
  different values than normally used? Should probably support that.
- Inverse (cut out within range) subset_circ() & subset_rect() options
- Implement better linear algebra to reduce the computation time of coords.get_rvpm()
- Play around with turning off mwahpy_glob.verbose flag for some things to quiet unnecessary output

BUGS
========================================

 No known bugs.
