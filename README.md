# Senior Design Heliostat
Ben Procknow senior design project - Heliostat design

This part of the project is a visualization tool used to simulate a multi-row heliostat system.

A heliostat is part of a Concentrated Solar Power (CSP) system that is used to redirect sunlight to a collection site, called a power tower, using mirrors, where the energy can be stored and used to generate electricity.  The mirrors are controlled by two motors, which control two axises of the mirror rotation.  This project is different in that instead of controlling one mirror with two motors, there are four mirrors linked together controlled by two motors.  In order for this system to work, there needs to be an extra mirror for each rotated mirror (i.e for every extra mirror linked together that is controlled by the motors, there needs to be an extra stationary mirror).  In the case of four mirrors linked together and controlled by two motors, there are four stationary mirrors.  To reach the collection site, the sunlight is bounced off the rotated mirror set, to the stationary mirror set, then to the power tower.  In a typical Concentrated Solar Power system there are hundreds to thousands of heliostats in a field around a power tower. 

The goal of the senior design project was to reduce the cost of the heliostat in order to reduce the cost of the entire CSP system cost.  Right now the heliostat accounts for around 50% of the cost of a CSP system.  By linking mirrors together, the idea is that by reducing the number of expensive motors used to control the mirrors the cost per square meter of light reflected to the power tower would decrease despite doubling the number of mirrors needed.

The visualization tool is used to optimize different parts of the heliostat system including mirror sizes, find optimal locations and orientations for the heliostats, and to find where error in the system can be tolerated that doesn't have a large impact on the efficiency of each heliostat.  

Usage:
Needs Python 2.7 to run because of the module "pygame" that is used as the screen to display the heliostat. 
Modules:  pygame, openGL (pyOpenGL), numpy

Press "p" to see a simulation of the heliostat moving throughout the day, tracking the sun and moving such that the light hits the power tower.
