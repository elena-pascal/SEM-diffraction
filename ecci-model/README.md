# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* The Python scripts collected here can be used to compute and plot the strain profile of a threading dislocation.
* A Fortran module containg the beta as a numerical function can also be obtained.
* The beta module can then be called by the Fortran Howie-Whelan implementation to calculate intensity profiles as observed in ECCI around a dislocation. 
* Version 1.0 released in March 2017
* For the physics applied see E. Pascal, Dislocation contrast in electron channelling contrast images as projections of strain-like components, Materials Today: Proceedings, 2017

### Set up? ###

* Requirement: Fortran compiler and Python 2.7 or later.
* geometry.py contains crystalographic calculus for a cubic or hexagonal crystal.
* coordinates.py contains the relationship between the different reference frames.
* calculateBeta2.py calculates the beta function as an object dependent on the position from the dislocation, the tilt and rotation of the sample, the diffraction condition.
* generateBetaModules.py uses previous funciton to create Fortran callable modules containg the beta function. 
* plotField.py is used to the displacement field.
* plotStrain.py is used for a number of plots for the strain profile.
* the fortran modules should be build using the Makefile.


### Future work ###
* A version of screw dislocation will soon be available.
* Mixed dislocations.
* Multiple dislocations.
* Non Laue symmetric version (we don't expect this to make significant quantitative difference)

### Contact ###

* Elena Pascal
* elena.pascal@strath.ac.uk

### Licence ###
	# Copyright (c) Elena Pascal Physics, University of Strathclyde
    #  All rights reserved.
	#
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation version 3.
	#
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
	#
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <http://www.gnu.org/licenses/>.