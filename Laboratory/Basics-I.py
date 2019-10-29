# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Tue Oct 29 20:43:25 2019

Description:
    Solving BASICS-I assignment of MGOD course 
    1. Conversion from Kepler Elements to Cartesian coordinates
    2. Conversion from Cartesian coordinates to Kepler elements
    

"""
import math # to use value of Pi
import numpy as np #NumPy library

mu = 398600.441E9 # Gravitational parameter for Earth [m^3/s^2]
r2d = 180/math.pi #bconversion factor for radians to degrees

# State vector in cartesian coordinates:
S_bar = np.array([8751268.4691, -7041314.6869, 4846546.9938, 332.2601039, -2977.0815768, -4869.8462227]) 


####################################################
######### Cartesian to Kepler Conversion ###########
####################################################
## 1. Computation of required Vectors and their norms
r_bar = S_bar[:3] # position vector [m]
V_bar = S_bar[3:6] # velocity vector [m/s]
    
h_bar = np.cross(r_bar,V_bar) # angular momentum vector
h = np.linalg.norm(h_bar) # magnitude of angular momentum vector
r = np.linalg.norm(r_bar) # magnitude of position vector [m]
V = np.linalg.norm(V_bar) # magnitude of velocity vector [m/s]

N_bar = np.cross(np.array([0,0,1]),h_bar) # N vector
Nxy = np.linalg.norm(N_bar) # squareroot of N_bar(1,1)^2 + N_bar(2,1)^2;
e_bar = np.cross(V_bar,h_bar)/mu - r_bar/r # eccentricity vector

####################################################