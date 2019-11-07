# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Tue Nov  5 21:21:57 2019

Version: 1

Description:
    Solving BASICS-I assignment of MGOD course 
    1. Conversion from Kepler Elements to Cartesian coordinates 
            - through a function
    2. Conversion from Cartesian coordinates to Kepler elements
        - through a function
       
    In follow-ups:
    1. Generilse the functions as part of a module, called from another file. 
       Extract mu (manually/SPICE) in this file and pass it to the function.
    2. Use switch cases to decide what kind of anomaly is used when Kepler to Cartesian.

    

"""
import math # to use value of Pi
import numpy as np #NumPy library





def convertCartesianToKepler(S_bar,mu,isOutputInDegree = False,isPrint=False):
    """Converts cartesian coordinates into equivalent Kepler elements."""

    r2d = 180/math.pi #bconversion factor for radians to degrees
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
        
    ############### Conversion Complete #################    
       
    ## 8. Saving the results in one array and displaying results, if required   
    if isOutputInDegree == True:
        Kepler = np.array([a, e, i*r2d, RAAN*r2d, omega*r2d, theta*r2d, E*r2d, M*r2d])
        if  isPrint == True:
            print("\nConversion from cartesian coordinates to Kepler elements is successful.")
            print ("Converted Kepler elements are:")    
            print ("Semi-major axis [m]: a = ",a)
            print ("Eccentricity: e = ",e)
            print ("Inclination [deg]: i = ",i*r2d)
            print ("RAAN [deg]: \u03A9 = ",RAAN*r2d)
            print ("Argument of periapsis [deg]: \u03C9 = ",omega*r2d)
            print ("True anomaly [deg]: \u03B8 = ",theta*r2d)
            print ("Eccentric anomaly [deg]: E = ",E*r2d)
            print ("Mean anomaly [deg]: M = ",M*r2d)
    else:
        Kepler = np.array([a, e, i, RAAN, omega, theta, E, M])
        if  isPrint == True:
            print("\nConversion from cartesian coordinates to Kepler elements is successful.")
            print ("Converted Kepler elements are:")    
            print ("Semi-major axis [m]: a = \t",a)
            print ("Eccentricity: e = \t",e)
            print ("Inclination [rad]: i = \t",i)
            print ("RAAN [rad]: \u03A9 = \t",RAAN)
            print ("Argument of periapsis [rad]: \u03C9 = \t",omega)
            print ("True anomaly [rad]: \u03B8 = \t",theta)
            print ("Eccentric anomaly [rad]: E = \t",E)
            print ("Mean anomaly [rad]: M = \t",M)
    
    return Kepler    

def convertKeplerToCartesian(KeplerElements,mu,typeOfAnomaly =1,isInputInDegree = False,isPrint=False):
    
    """Converts Kepler elements into equivalent Cartesian coordinates."""

    r2d = 180/math.pi #bconversion factor for radians to degrees
    ####################################################
    ######### Kepler to Cartesian Conversion ###########
    ####################################################
    ## 1. Setting up Kepler elements
    a = KeplerElements[0] # Semi-major axis[m]
    e = KeplerElements[1] # Eccentricity
    i = KeplerElements[2] # Inclination [radian]
    RAAN = KeplerElements[3] # Right ascention of ascending node [radian]
    omega = KeplerElements[4] # Argument of pericenter [radian]
    
    ## 1. Converting angles from degree to radian, if required
    if isInputInDegree==True:
        i = i/r2d # inclination [radian]
        RAAN = RAAN/r2d # RAAN [radian]
        omega = omega/r2d # Argument of periapsis [radian]
        KeplerElements[5] = KeplerElements[5]/r2d # Anomaly[radian] 
    
    ## 2. Deciding the anomaly type passed        
    if typeOfAnomaly == 1: # Mean anolmaly is passed
        M = KeplerElements[5] # Mean anomaly[radian] 
        E = 0
        theta = 0
    elif typeOfAnomaly == 2:  # Eccentric anolmaly is passed    
        E = KeplerElements[5] # Eccentric anomaly[radian] 
        theta = 0
    else:
        typeOfAnomaly = 1 # Assuming that Mean anomaly was passed
        M = KeplerElements[5] # Mean anomaly[radian] 
            
    ## 3. Computation of Eccentric and True anomalies, if required  
    if typeOfAnomaly == 1: 
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
    elif typeOfAnomaly == 2: 
        theta = 2 * math.atan( math.tan(E/2) * math.sqrt((1+e)/(1-e)))
        
    ## 4. Computation of r
    r = a * (1 - e * math.cos(E))
    
    ## 5. Computation of required variables
    xi = r * math.cos(theta) 
    eta = r * math.sin(theta)
    l1 = math.cos(RAAN)*math.cos(omega) - math.sin(RAAN)*math.sin(omega)*math.cos(i)
    l2 = -math.cos(RAAN)*math.sin(omega) - math.sin(RAAN)*math.cos(omega)*math.cos(i)
    m1 = math.sin(RAAN)*math.cos(omega) + math.cos(RAAN)*math.sin(omega)*math.cos(i)
    m2 = -math.sin(RAAN)*math.sin(omega) + math.cos(RAAN)*math.cos(omega)*math.cos(i)
    n1 = math.sin(omega)*math.sin(i)
    n2 = math.cos(omega)* math.sin(i)
    
    ## 6. Computation of position and velocity coordinates
    S = np.array([[l1,l2],[m1,m2],[n1,n2]])
    S2 = np.array([[xi],[eta]]) 
    S = np.dot(S,S2)
    
    ## 7. Magnitutde of specific angular momentum vector
    h = math.sqrt(mu*a*(1-e**2))
    
    ## 8. Computation of velocity components
    x_dot = mu * (-l1*math.sin(theta)+l2*(e+math.cos(theta)))/h # x velocity [m/s]
    y_dot = mu * (-m1*math.sin(theta)+m2*(e+math.cos(theta)))/h # y velocity [m/s]
    z_dot = mu * (-n1*math.sin(theta)+n2*(e+math.cos(theta)))/h # z velocity [m/s]
    
    ############### Conversion Complete #################    
       
    ## 9. Saving the results in one array and displaying results, if required   
    Cartesian = np.array([S[0][0], S[1][0], S[2][0], x_dot, y_dot, z_dot])
    if  isPrint == True:
        print("\nConversion from Kepler elements to Cartesian coordinates is successful.")
        print ("Converted Cartesian coordinates are:")    
        print ("X-position [m]: x = \t",S[0][0])
        print ("Y-position [m]: y = \t",S[1][0])
        print ("Z-position [m]: z = \t",S[2][0])   
        print ("X-velocity [m/s]: x_dot = \t", x_dot)
        print ("Y-velocity [m/s]: y_dot = \t", y_dot)
        print ("Z-velocity [m/s]: z_dot = \t", z_dot) 
    
    return Cartesian


#########################################################
############    Testing the functions ####################
#########################################################
    
# State vector in cartesian coordinates:
# S_bar = np.array([8751268.4691, -7041314.6869, 4846546.9938, 332.2601039, -2977.0815768, -4869.8462227]) 
# S_bar = np.array([-2700816.14,-3314092.80, 5266346.42, 5168.606550, -5597.546618, -868.878445]) 
S_bar = np.array([3126974.99, -6374445.74, 28673.59, -254.91197, -83.30107, 7485.70674]) 
mu = 398600.441E9 # Gravitational parameter for Earth [m^3/s^2]
    
# Kepler elements   
#a = 12158817.9615 # Semi-major axis[m]
#e = 0.014074320051 # Eccentricity
#i = 52.666016957 # Inclination [deg]
#RAAN = 323.089150643 # Right ascention of ascending node [deg]
#omega = 148.382589129 # Argument of pericenter [deg]
#M = 112.192638384 # Mean anomaly[deg] 

#a = 6787746.891  # Semi-major axis[m]
#e = 0.000731104 # Eccentricity
#i = 51.68714486 # Inclination [deg]
#RAAN = 127.5486706 # Right ascention of ascending node [deg]
#omega = 74.21987137 # Argument of pericenter [deg]
#M = 24.06608426 # Mean anomaly[deg] 

a = 7096137  # Semi-major axis[m]
e = 0.0011219 # Eccentricity
i = 92.0316 # Inclination [deg]
RAAN = 296.1384 # Right ascention of ascending node [deg]
omega = 120.6878 # Argument of pericenter [deg]
M = 239.6546 # Mean anomaly[deg]   
theta = 239.5437 # True Anomaly [deg]  
E = 239.5991 # Eccentric Anomaly [deg]

Kepler = np.array([a,e,i,RAAN,omega,M])
Kepler2 = np.array([a,e,i,RAAN,omega,theta])
Kepler3 = np.array([a,e,i,RAAN,omega,E])

CovertedKep = convertCartesianToKepler(S_bar,mu,False,True)    # Position arguments are passed

ConvertedCarte = convertKeplerToCartesian(Kepler,mu,1,isInputInDegree = True, isPrint=True)

ConvertedCarte4 = convertKeplerToCartesian(Kepler,mu,isInputInDegree = True, isPrint=True)

ConvertedCarte2 = convertKeplerToCartesian(Kepler3,mu,typeOfAnomaly =2,isInputInDegree = True, isPrint=True)
