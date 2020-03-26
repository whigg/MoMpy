# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Fri Nov  8 20:57:12 2019

version: Base

Description:
    1. Generic canvas to test different modules
"""

from Modules.BasicAstrodynamics import convertCartesianToKepler
from Modules.BasicAstrodynamics import convertKeplerToCartesian

## following method to import all the functions from a module is not recommended 
## in the cheat sheets
#from Modules.BasicAstrodynamics import * 

import math # to use value of Pi
import numpy as np #NumPy library

import spiceypy as spice

import pygmo



# Part(Problem) 1 of Basics-I assignment 
S_bar = np.array([8751268.4691, -7041314.6869, 4846546.9938, 332.2601039, -2977.0815768, -4869.8462227]) 


# Extracting mu from SPICE kernels
print (spice.tkvrsn('TOOLKIT'))

spice.furnsh("./External_files/Spice_kernels/kernel_load.txt")

muE = spice.bodvrd( 'EARTH', 'GM', 1 );

#mu = 398600.441E9 # Gravitational parameter for Earth [m^3/s^2]
mu = muE[1][0] * 1e9

spice.kclear

CovertedKep = convertCartesianToKepler(S_bar,mu,True,True)    # Position arguments are passed

# Part(Problem) 2 of Basics-I assignment 
a = 12158817.9615 # Semi-major axis[m]
e = 0.014074320051 # Eccentricity
i = 52.666016957 # Inclination [deg]
RAAN = 323.089150643 # Right ascention of ascending node [deg]
omega = 148.382589129 # Argument of pericenter [deg]
M = 112.192638384 # Mean anomaly[deg] 

Kepler = np.array([a,e,i,RAAN,omega,M])

ConvertedCarte = convertKeplerToCartesian(Kepler,mu,1,isInputInDegree = True, isPrint=True)
ConvertedCarte2 = convertKeplerToCartesian(Kepler,mu,7,isInputInDegree = True, isPrint=True)

