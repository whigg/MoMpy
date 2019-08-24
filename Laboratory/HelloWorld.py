# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: 24 August 2019

Description:
    First script in Python. Printing "Hello World". 
    Includes definition of a function and calling it.
    Includes experiments with user input and basic error handling, 
    various variable types and simple algebra. 
    

"""

# this is a comment line
## how different is this from using single '#'?


print ("Hello World !") # this is adding a comment on the side of a statement

###################################
# defining a function

def cohle():
    print ("Death created time to grow the things it would kill...")
   
# function ends when indentation is stopped

cohle()
##########################################################
##### simple algebric calculations #######################
##########################################################

import math

print ("To solve the quadratic eqn: ax^2 + bx + c = 0. a, b, c should have integer values.")

########################################
##### with simple error handling #######
########################################

S = 1

try:

    a = int(input("Enter the value of a:"))
except ValueError:
    S = 0;
    print("Invalid input")

try:

    b = int(input("Enter the value of b:"))
except ValueError:
    S = 0;
    print("Invalid input")
    
try:

    c = int(input("Enter the value of c:"))
except ValueError:
    S = 0;
    print("Invalid input")    


if S==0:
    print ("At least one of the user input is invalid")
else:
########################################
##### without simple error handling #######
########################################
    
    
#a = int(input("Enter the value of a:"))
#b = int(input("Enter the value of b:"))
#c = int(input("Enter the value of c:"))
       
    D = b**2 - 4*a*c 
    
    if D <0:
    
        print ("This equation has no solutions")
        
    else:    
        x1 = (-b - math.sqrt(D)) / (2*a)
        x2 = (-b + math.sqrt(D)) / (2*a)
        
        print("The solution to above quadractic equation is:")
        print ("x1 = ",x1)
        print ("x2 = ",x2)









   