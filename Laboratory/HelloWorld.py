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
##### variable assignments and operations ################
##########################################################
print ("\n\n%%%%%Variable types and their assignments:\n")

#Integers
i = 200
k = -5
divider = 4

print ("Value of k is",k, ". And the division -5/4 is: ", k/divider, " (even when integer type variables are used.")

print ("Division -5/4 when the variables are explictly converted into integers: ", int(k)/int(divider))

print ("Modulo function: -5%4 is:",k%divider)

print ("Modulo function: 3/2 is:", 3%2)

print ("Modulo function: -6/2 is:", -6%2)

print ("Modulo function: -7/2 is:", -7%2, "\n")

#Floats
pi = 3.142
trillion = 1e12
milli = 1e-3
half = .5

#Logicals
isEarthflat = False
isthisPython = True

#Strings
name = "MoMPy"
purpose = 'fun'
number = "100"
floatN = "200.0"

print ("Converting string \"100\" into integer: ", int(number))
print ("Converting string \"200.0\" into float: ",float(floatN))

print ("The string - name - is:", name, " and its length is ",len(name))
print("The string - floatN is: ",floatN, " and its length is ",len(floatN))
print("The concentated string is: ", name + "@" + floatN)
print ("Converting integer 100 into a string and then converting float 200.0 into a string and concenating them: ", str(number) + " and " + str(floatN))
#print ("Evaluating string \"fun\" we get: ", eval("x+2=0"))
print ("Converting ascii code 22 into character: ",chr(22))
print ("ASCII code for 't' is", ord('t'))

print("%%%%%%%%%%%Accessing parts of a string \n\n")
Original = name + number + " " + purpose
print("Original string: ", Original )
print ("Original[3]: ",Original[3])
print ("Original[2:5]: ",Original[2:5])
print("Original[:3]: ",Original[:3])
print ("Original[6:2:-1]: ",Original[6:2:-1])
print ("Original[::3]: ",Original[::3])
print("Original[-1]: ",Original[-1])
print("Original[-2]: ",Original[-2])
print ("Original[:4]+Original[-2]: ",Original[:4]+Original[-2])


##########################################################
##### simple algebric calculations #######################
##########################################################

import math

print ("\n\n%%%%%%To solve the quadratic eqn: ax^2 + bx + c = 0.\n")
print ("a, b, c can have float values.")
########################################
##### with simple error handling #######
########################################

S = 1

try:

    a = float(input("Enter the value of a:"))
except ValueError:
    S = 0;
    print("Invalid input")

try:

    b = float(input("Enter the value of b:"))
except ValueError:
    S = 0;
    print("Invalid input")
    
try:

    c = float(input("Enter the value of c:"))
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









   