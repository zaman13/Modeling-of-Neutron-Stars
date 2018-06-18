# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 21:40:45 2014

@author: Mohammad Asif Zaman
"""


# Using Plank system units (hcut = c = 1)
# Ref.: http://newfi.narod.ru/Newfi/Constants.htm

hc = 197.327                # Conversion factor in MeV fm (hut * c)
G = hc * 6.67259e-45        # Gravitational constant
Ms = 1.1157467e60
rho_s = 1665.3              # Central density (density at r = 0)
M0 = (4*3.14159265*(G**3)*rho_s)**(-0.5)
R0 = G*M0
mn = 938.926                # Mass of neutron in MeV c^-2


# Defining the functions

# Function for determining initial value of n(r=0)
def initial_n():          
    n = 1
    err = 1
    tol = 1e-15
    count = 0
    # Newton-Raphson method
    while err > tol : 
        count += 1
        fn = n*mn + 236*n**(2.54) - rho_s
        dfn = mn + 236*2.54*n**(1.54)
        temp = n - fn/dfn
        err = np.abs(n-temp)
        n = temp
    print "Newton-Raphson Converged after ", count, "iterations"
    return n
    
    
def rho(p):
    n = (p*rho_s/363.44)**(1./2.54)
    return (236. * n**2.54 + n *mn)/rho_s 
    

def dp_dr(r,m,p,flag):
    if flag == 0:                              # classical model
        y = -m*rho(p)/(r**2 + 1e-20)
    else:                                      # relativistic model
        rh = rho(p)                            
        y = -(p+rh)*(p*r**3 + m)/(r**2 - 2*m*r + 1e-20)
    return y

def dm_dr(r,m,p):
    return rho(p)*r**2
    

def EulerSolver(r,m,p,h,flag):
    y = np.zeros(2)
    y[0] = m + dm_dr(r,m,p)*h
    y[1] = p + dp_dr(r,m,p,flag)*h
    return y
    
def RK4Solver(r,m,p,h,flag):
    y = np.zeros(2)
    k11 = dm_dr(r,m,p)
    k21 = dp_dr(r,m,p,flag)
    
    k12 = dm_dr(r+0.5*h,m+0.5*k11*h,p+0.5*k21*h)
    k22 = dp_dr(r+0.5*h,m+0.5*k11*h,p+0.5*k21*h,flag)
    
    k13 = dm_dr(r+0.5*h,m+0.5*k12*h,p+0.5*k22*h)
    k23 = dp_dr(r+0.5*h,m+0.5*k12*h,p+0.5*k22*h,flag)
    
    k14 = dm_dr(r+h,m+h*k13,p+h*k23)    
    k24 = dp_dr(r+h,m+h*k13,p+h*k23,flag)    
    
    y[0] = m + h*(k11 + 2.*k12 + 2.*k13 + k14)/6.
    y[1] = p + h*(k21 + 2.*k22 + 2.*k23 + k24)/6.
    return y
    
    
def mplot(fign,x,y,xl,yl,clr,lbl):
    py.figure(fign)
    py.xlabel(xl)    
    py.ylabel(yl)
    return py.plot(x,y,clr, linewidth =2.0,label = lbl)
 




import time
import math
import matplotlib as mp
import numpy as np
import scipy as sp
import pylab as py




Np = 100
rho_s_set = np.logspace(1,4.5,Np)
global rho_s

N = 4000
r = np.linspace(0,15,N)
h = r[1]-r[0]



rf = np.zeros([2, Np])
mf = np.zeros([2, Np])

flag_set = [0,1]

print "Simulation range, R = 0 to ", r[-1]*R0*1e-18, "km" 
tol = 1e-4

for glb in range(1,Np):
    m = np.zeros(N)
    p = np.zeros(N)
    rh = np.zeros(N)
       
    r[0] = 0
    m[0] = 0
    rh[0] = 1
    
    rho_s = rho_s_set[glb]
    M0 = (4*np.pi*(G**3)*rho_s)**(-0.5)
    R0 = G*M0
    ni = initial_n()

    p[0] = 363.44 * (ni**2.54)/rho_s
    
    for k in range(0,2):
        flag = flag_set[k]
        for i in range(0,N-1):
            [m[i+1], p[i+1]] = EulerSolver(r[i],m[i],p[i],h,flag)
            rh[i+1] = rho(r[i])
            if p[i+1] < tol:
                rf[flag, glb] = r[i] - p[i]*(r[i]-r[i-1])/(p[i]-p[i-1])
                mf[flag, glb] = m[i]  -p[i]*(m[i]-m[i-1])/(p[i]-p[i-1])
                mf[flag, glb] = mf[flag, glb]*M0/Ms
                rf[flag, glb] = rf[flag, glb]*R0*1e-18
                break
        
        print 
        print "Running RK4"
        if i == N-2:
            print "Program didn't converge to P = 0, extend the maximum value of r"
        else:
            print "P <", tol, "found after", i, "runs"
        print 
        
            
               
        if flag == 0:
            lbl = "Classical Model"
            clr = "blue"
        else:
            lbl = "Relativistic Model"
            clr = "red"
        
        print 
        print "=============================================="
        print lbl, "Results"
        print "=============================================="
        print "Initial density, rho_s = ", rho_s, "MeV/fm3"
        print "Total mass = ", mf[flag, glb], "times Solar mass"
        print "Radius of the Neutron star = ", rf[flag, glb], "km"
    

mplot(1,rf[0,:],mf[0,:],'Radius, $R$ (km)',r'Mass, $M/M_{sun}$','blue','Classical Model')
q = py.legend(loc = 0)
q.draw_frame(False)
mplot(2,rf[1,:],mf[1,:],'Radius, $R$ (km)',r'Mass, $M/M_{sun}$','red','Relativistic Model')
q = py.legend(loc = 0)
q.draw_frame(False)



