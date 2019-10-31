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
# S_bar = np.array([8751268.4691, -7041314.6869, 4846546.9938, 332.2601039, -2977.0815768, -4869.8462227]) 
# S_bar = np.array([-2700816.14,-3314092.80, 5266346.42, 5168.606550, -5597.546618, -868.878445]) 
S_bar = np.array([3126974.99, -6374445.74, 28673.59, -254.91197, -83.30107, 7485.70674]) 


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


## 2. Computation of orbital elements
a = 1/(2/r - V**2/mu) # Semi-major axis [m]
e = np.linalg.norm(e_bar) # Eccentricity
i = math.acos(h_bar[2]/h) # Inclination [radian]
RAAN = math.atan2(N_bar[1]/Nxy,N_bar[0]/Nxy) # RAAN [radian]

## 3. Settting up the sign for argument of peri-center
if np.dot(np.cross(N_bar,e_bar),h_bar) > 0:
    sign1 = 1
else:
    sign1 = -1
    
## 4. Setting up the sign for true anomaly 
if np.dot(np.cross(e_bar,r_bar),h_bar) > 0:
    sign2 = 1
else:
    sign2 = -1
    
## 5. Computation of orbital elements - continued
omega = sign1 * math.acos(np.dot((e_bar/e),(N_bar/Nxy))) # Argument of periapsis [radian]   
theta = sign2 * math.acos(np.dot((r_bar/r),(e_bar/e))) # True anomaly [radian]


## 6. Computation of mean and eccentric anomalies
E = 0
M = 0
if e > 0 and e < 1:
    E = 2 * math.atan(math.tan(theta/2)*math.sqrt( (1-e) / (1+e) )) # Eccentric anomaly [radian]
    M = E - e * math.sin(E) # Mean anomaly [radian]

## 7. Correcting for negative values of angles if any   
if RAAN < 0:
    RAAN = RAAN + 2 * math.pi
if omega < 0:
    omega = omega + 2 * math.pi
if theta < 0:    
    theta = theta+ 2 * math.pi
if E < 0:
    E = E + 2 * math.pi
if M < 0:
    M = M + 2 * math.pi

## 8. Displaying results
print("Conversion from cartesian coordinates to Kepler elements is successful.")
print ("Converted Kepler elements are:")    
print ("Semi-major axis [m]: a = ",a)
print ("Eccentricity: e = ",e)
print ("Inclination [deg]: i = ",i*r2d)
print ("RAAN [deg]: \u03A9 = ",RAAN*r2d)
print ("Argument of periapsis [deg]: \u03C9 = ",omega*r2d)
print ("True anomaly [deg]: \u03B8 = ",theta*r2d)
print ("Eccentric anomaly [deg]: E = ",E*r2d)
print ("Mean anomaly [deg]: M = ",M*r2d)

####################################################
######### Kepler to Cartesian Conversion ###########
####################################################


## 1. Setting up Kepler elements
a = 12158817.9615 # Semi-major axis[m]
e = 0.014074320051 # Eccentricity
i = 52.666016957 # Inclination [deg]
RAAN = 323.089150643 # Right ascention of ascending node [deg]
omega = 148.382589129 # Argument of pericenter [deg]
M = 112.192638384 # Mean anomaly[deg] 

## 1. Converting angles from degree to radian
i = i/r2d # inclination [radian]
RAAN = RAAN/r2d # RAAN [radian]
omega = omega/r2d # Argument of periapsis [radian]
M = M/r2d # Mean anomaly[radias] 

## 1. Computation of Eccentric and True anomalies  
E0 = math.pi
E1 = 2 * math.pi
while True: 
    E1 = E0 + (M-E0+e*math.sin(E0))/(1-e*math.cos(E0))
    if abs(E1 - E0) > 1e-20:
        E0 = E1
    else:
        E = E1
        break
theta = 2 * math.atan( math.tan(E/2) * math.sqrt((1+e)/(1-e)))

## 3. Computation of r
r = a * (1 - e * math.cos(E))

## 4. Computation of required variables
xi = r * math.cos(theta) 
eta = r * math.sin(theta)
l1 = math.cos(RAAN)*math.cos(omega) - math.sin(RAAN)*math.sin(omega)*math.cos(i)
l2 = -math.cos(RAAN)*math.sin(omega) - math.sin(RAAN)*math.cos(omega)*math.cos(i)
m1 = math.sin(RAAN)*math.cos(omega) + math.cos(RAAN)*math.sin(omega)*math.cos(i)
m2 = -math.sin(RAAN)*math.sin(omega) + math.cos(RAAN)*math.cos(omega)*math.cos(i)
n1 = math.sin(omega)*math.sin(i)
n2 = math.cos(omega)* math.sin(i)

## 5. Computation of position and velocity coordinates
S = np.array([[l1,l2],[m1,m2],[n1,n2]])
S2 = np.array([[xi],[eta]]) 
S = np.dot(S,S2)

    
        
####################################################