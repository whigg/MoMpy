# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: 24 August 2019

Description:
    First script in Python. Printing "Hello World". 
    Includes definition of a function and calling it.
    Includes experiments with user input, 
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
###################################
### simple algebric calculations

import math

print ("To solve the quadratic eqn: ax^2 + bx + c = 0:")

a = input("Enter the value of a:")
b = input("Enter the value of b:")
c = input("Enter the value of c:")
 
D = b**2 - 4*a*c #problem with this statement

if D <0:

    print ("This equation has no solutions")
    
else:    
    x1 = (-b - math.sqrt(D)) / (2*a)
    x2 = (-b + math.sqrt(D)) / (2*a)
    
    print ("x1 = "),x1
    print ("x2 = "),x2









   