 MILKYWAY@HOME PYTHON PACKAGE (MWAHPY)
========================================

COPYRIGHT TOM DONLON, 2020 RPI

github user: thomasdonlon

REQUIRES PYTHON >3.2.3

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

NOT TESTED AND NOT GUARANTEED TO WORK

1. Download the mwahpy github repository

2. Place the files in the repository in your sys-path directory
(this can be found by running "python3 -m site" in your terminal
and is given after "USER-SITE")
(also see https://stackoverflow.com/questions/122327/how-do-i-find-the-location-of-my-python-site-packages-directory)

3. Insert "import mwahpy" in your code

4. Do science

Alternatively, if you don't want to dig around in your python directory then
you can always place the files from the repository in whatever directory you
are using for your python code that will use mwahpy, and then import it as you
would any other file. This may get gross though, as I expect that the package
will expand into several dependency files over time.

I plan on eventually figuring out how to make it installable with pip, but
that's not a huge concern at the moment, as this is not meant to be a
widely used package, and is oriented mainly at a specific group.

CONTENTS
========================================

A non-exhaustive list of contents of the package are listed below:

A non-exhaustive list of contents to be added to the package at a later
date are listed below:

 - easy importing of data from N-body output
 - easy manipulation of data after reading in
 - easy visualization of the imported data through plotting functionality
