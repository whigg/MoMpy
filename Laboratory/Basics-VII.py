# -*- coding: utf-8 -*-
"""

Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Sat Mar  7 14:42:52 2020

version: Base

Description: (Original assignment)
    1. Earth repeat orbits - without functions, approach 1: Question 2A 
    2. Earth repeat orbits - without functions, approach 2: Question 1, 2B 
"""

import math

import spiceypy as spice

import numpy as np

import matplotlib.pyplot as plt

import sys

spice.furnsh("./External_files/Spice_kernels/kernel_load.txt")

Question = 3 # 1 - Question 1, 2 - Question 2A, 3 - Question 2B, 
           

if Question ==1:
############ Question 1 - validation, solved using approach 2 
    
    
    i = 28.0 # %deg
    e = 0.0 # assumption
    D = 86164.1004 #s - import from spice
    muE = spice.bodvrd( 'EARTH', 'GM', 1 ); #km3/s2 for earth 
    mu = muE[1][0]
    J2 = 1082.63E-6 #J2 for earth - import from spice
    RE = spice.bodvrd('EARTH','RADII',3)
    Re = RE[1][0] # Radius of Earth, from SPICE - pck00010.tpc, 6378.1366 km
    L_dot = 360.0
    r2d = 1/spice.rpd() # Radians per degree
    print('r2d2 is: ',r2d)
    k2 = 0.75 * J2 * math.sqrt(mu)* Re*Re
    
    
    
    # Method 1 # READ 10.2.1. of the course handbook
    jkStore = np.matrix('14, 1; 43, 3; 29, 2; 59, 4; 74, 5 ; 15,1')
    #jk = jkStore[0][0]
    #print(jk)
    
    
    ## Method 2 # READ 10.2.1. of the course handbook
    #jkStore2 = np.array([14, 1, 43, 3, 29, 2, 59, 4, 74, 5 , 15, 1])
    #jkStore2.shape = [6,2] # READ 10.2.1. of the course handbook
    
    h_store = np.zeros(int(jkStore.size /2))
    a_store = np.zeros(int(jkStore.size /2))
    T_store = np.zeros(int(jkStore.size /2))
    
    
    # for more on 'range', READ 7.2.2. of the course handbook
    for row in range(1,1+int(jkStore.size /2)): 
       j = jkStore[row-1,0]
       k = jkStore[row-1,1]
       
       a0 = math.pow((mu*D**2*k**2/(4*math.pi**2*j**2)),1/3)
       
       iterations = 0
       
       while True:
           iterations = iterations+1
           
           #Computations
           RAAN_dot = - 2 * k2 * math.pow(a0,-3.5) * math.cos(i/r2d) * math.pow(1-e**2,-2) * D * r2d
    #       print(RAAN_dot)
           omega_dot = k2 * math.pow(a0,-3.5) * (5*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-2) * D * r2d
    #       print(omega_dot)
           M_dot = k2 * math.pow(a0,-3.5) * (3*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-1.5) * D * r2d
    #       print(M_dot)
           n = (j/k) * (L_dot - RAAN_dot) -(omega_dot + M_dot)
    #       print(n)
           a1 = (mu/(n/(D*r2d))**2)**(1/3)
    #       print('a1 is:',a1)
           
    #       print('At iteration:',iterations,' a1 was ',a1, 'a0 was', a0)
           if iterations > 1000:
               print('Maximum number of iterations executed')
               a0 = a1
               break
           elif (abs(a1-a0) < 10**(-10)):  
               print('Convergence reached after ',iterations, 'iterations')    
               a0 = a1
               break
           else:
               a0 = a1
           
       a_store[row-1] = a0
       h_store[row-1] = a0-Re # Altitude in km 
       T_store[row-1] = 2 * math.pi *math.sqrt(a0**3/mu)/60; #Period in minutes
            
    print('Altitudes are:',h_store)
    print('Corrosponding periods are:',T_store)

elif Question ==2 or Question == 3 or Question == 4:
############ Question 2, solved using approach 1 (Q2A)
    ## input parameters
    e = 0
    De = 86164.1004 # s, to be extracted from spice
    muE = spice.bodvrd('EARTH','GM',1)
    mu = muE[1][0]
    J2 = 1082.63E-6 # %J2 for earth, to be extracted from spice
    RE = spice.bodvrd('EARTH','RADII',3)
    Re = RE[1][0]
    r2d = 1/spice.rpd()
    k2 = 0.75 * J2 * math.sqrt(mu)* Re*Re
        
    k = 3 
    j_store = np.arange(39,49)  
    hMax = 1200
    hMin = 200
     
    # Storage
    hi_store = np.zeros((int(j_store.size)) * 3*((hMax-hMin)+1))
    hi_store.shape = [int(j_store.size), ((hMax-hMin)+1),3]
    
    # Solution tracker
    Solindex = np.zeros(int(j_store.size)*3)
    Solindex.shape = [int(j_store.size),3]
    
    
    # computations
    for count in range(0,int(j_store.size)):
        j = j_store[count]
        Solindex[count,0] = j
        minIndex = (hMax-hMin)+1
        maxIndex = 0
#        hiStore2 = np.zeros(6) # attempt to append row to a matrix - failure
#        hiStore2.shape = [2,3] # attempt to append row to a matrix - failure

        for h in range(hMin,hMax+1): # altitude variation between two limits
            if e == 0:
               a = Re + h
            else:
               sys.exit('Problem not defined yet for non-circular orbits.')
           
            if Question ==2: 
                T = 2 * math.pi * math.sqrt(a**3/mu)
                #           print(T)
                DeltaL1 = -2 * math.pi * (T/De)          
                # DeltaL2 = C1 * cosi
                C1 = (-3*math.pi * J2 * Re**2)/(a**2 * (1 - e**2)**2)
                C2 = -2 * math.pi* k/j  - DeltaL1
                C3 = 2 * math.pi* k/j  - DeltaL1
    #           print('DeltaL1: ', DeltaL1,', C1:', C1, ', C2:',C2, ', C3:',C3)
              
                # Inverse cosine
                if (C2 <0) & (abs(C2/C1) <=1):
                    DeltaL2 = C2
                    i = math.acos(DeltaL2/C1) * r2d # for i = [0, 90] deg 
                    hi_store[count,h-hMin,0] = h
                    hi_store[count,h-hMin,1] = i
                    if h-hMin < minIndex:
                        minIndex = h-hMin
                    if h-hMin > maxIndex:
                        maxIndex = h-hMin
                    # temp = np.array([h,i,0]) # attempt to append row to a matrix - failure
                    # np.append(hiStore2,temp) # attempt to append row to a matrix - failure
                    # print(hi_store)
                elif (C2 > 0) & (abs(C2/C1) <=1):
                    DeltaL2 = C2
                    i = math.acos(DeltaL2/C1) * r2d # for i = (90, 180] deg 
                    hi_store[count,h-hMin,0] = h
                    hi_store[count,h-hMin,1] = i
                    if h-hMin < minIndex:
                        minIndex = h-hMin
                    if h-hMin > maxIndex:
                        maxIndex = h-hMin    
                    # temp = np.array([h,i,0]) # attempt to append row to a matrix - failure
                    # np.append(hiStore2,temp) # attempt to append row to a matrix - failure
                elif (C3>0) & (abs(C3/C1) <=1):
                    DeltaL2 = C3
                    i = math.acos(DeltaL2/C1) * r2d # for i = (90, 180] deg 
                    hi_store[count,h-hMin,0] = h
                    hi_store[count,h-hMin,2] = i
                    if h-hMin < minIndex:
                        minIndex = h-hMin
                    if h-hMin > maxIndex:
                        maxIndex = h-hMin    
#                else:
#                    print('No solution at h:', h, ' and j:' , j)
                    
            elif Question == 3:
                C1 = -2 * k2 * (a**-3.5) * ((1-e**2)**-2) * De*r2d  # RAAN_dot = C1 * cos(i)
                C2 = k2 * (a**-3.5) *((1-e**2)**-2) * De*r2d # omega_dot = C2 *(5*(cosd(i))^2 -1)
                C3 = k2 * (a**-3.5) *((1-e**2)**-1.5) * De * r2d # M_dot = C3 *(3*(cosd(i))^2 -1)
                n = De * r2d * math.sqrt(mu/a**3);

                # For the quadratic equation A*x^2 + B*x +C = 0 with x = cos(i) 
                A = 5* C2 + 3 * C3
                B = C1 * j/k
                C = n - 360 *j/k - (C2 + C3)
                D = B**2-4*A*C    
                
                if D >=0:
                    if (D == 0) & (abs(-B/(2*A)) <=1):
                        i1 = math.acos(-B/(2*A))
                        hi_store[count,h-hMin,0] = h
                        hi_store[count,h-hMin,1] = i1*r2d
                        if h-hMin < minIndex:
                            minIndex = h-hMin
                        if h-hMin > maxIndex:
                            maxIndex = h-hMin 
                    elif (D > 0): 
                        x1 = (-B + math.sqrt(D))/ (2 * A)
                        x2 = (-B - math.sqrt(D))/ (2 * A)
                        if abs(x1) <=1:
                            i1 = math.acos(x1)
                            hi_store[count,h-hMin,0] = h
                            hi_store[count,h-hMin,1] = i1*r2d
                            if h-hMin < minIndex:
                                minIndex = h-hMin
                            if h-hMin > maxIndex:
                                maxIndex = h-hMin 
                               
                        if abs(x2) <=1:
                            i2 = math.acos(x2) 
                            hi_store[count,h-hMin,0] = h
                            hi_store[count,h-hMin,2] = i2*r2d
                            if h-hMin < minIndex:
                                minIndex = h-hMin
                            if h-hMin > maxIndex:
                                maxIndex = h-hMin 

        Solindex[count,1] = minIndex         
        Solindex[count,2] = maxIndex   
    
    # Plotting
    if Question ==2: 
        for count in range(0,j_store.size):
            plt.plot(hi_store[count,int(Solindex[count,1]):int(Solindex[count,2])+1,0], hi_store[count,int(Solindex[count,1]):int(Solindex[count,2])+1,1]) 
            plt.title("Earth-repeat orbits cosindering the effect of $J_2$ perturbation on $\Omega$ only")
    elif Question == 3:
        for count in range(0,j_store.size):
            plt.plot(hi_store[count,int(Solindex[count,1]):int(Solindex[count,2])+1,0], hi_store[count,int(Solindex[count,1]):int(Solindex[count,2])+1,2]) 
            plt.title("Earth-repeat orbits cosindering the effect of $J_2$ perturbation on $\Omega$, $\omega$, and $M$")
    plt.xlabel("Altitude [km]")
    plt.ylabel("Inclination [degree]")       
    plt.axis([200, 1200, 0, 180])        
       

spice.kclear #spice kernel clear   