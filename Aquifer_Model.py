#!/usr/bin/env python
import numpy as np
from scipy.stats import qmc
from functools import partial
from scipy.optimize import fsolve

#Latin Hypercube Sampling
sampler = qmc.LatinHypercube(d=5)
sample = sampler.random(n=5)
l_bounds = [39.5, 0.39, 30.0, 0.0001, 10.0]
u_bounds = [40.5, 1.38, 50.0, 0.0002, 30.0]
sample_scaled = qmc.scale(sample, l_bounds, u_bounds)
print(sample_scaled)

#a; Density Ratio
#q0; uniform discharge to the sea per unit length of coastline
#z0; depth below mean sea level of the aquifer bottom (m)
#W; Uniform Net Recharge (m/d)
#K; Hydraulic Conductivity (m/d)

# Basic Model of Aquifer System
def aquifer_model(a,W,K,x_length,x_i,q0,z0):

    
    print("Aquifer Length: " + str(x_length) + "m")
    print("a: " + str(a))
    print("q0: " + str(q0) + "m^2/day")
    print("z0: " + str(z0) + "m")
    print("W: " + str(W) + "m/day")
    print("K: " + str(K) + "m/day")
        
    a = float(a)
    q0 = float(q0)
    z0 = float(z0)
    W = float(W)
    K = float(K)

    #Position of Toe of Sea Water Wedge
    xt = (q0/W)-np.sqrt(((q0/W)**2)-((K*(1+a)*(z0)**2)/(W*(a)**2)))
    xt = round(xt,2)
    print("xt Value: " + str(xt) + "m")
           
    #Elevation of Water Table at Xt
    ht = np.sqrt(((2*q0*xt)-(W*xt**2))/(K*(1+a)))
    ht = round(ht,3)
    ht = ht
    print("ht Value: " + str(ht) + "m")

    #Elevation of Water Table at Inland Head
    h_inland = np.sqrt(((2/K)*(x_i-xt)*(q0-(W/2)*(x_i+xt))+(ht+z0)**2)-z0)
    h_inland = round(h_inland,3)
    print("h_inland value at x=2000: " + str(h_inland) + "m")

def new_q_formula(a,W,K,x_i,h_i,z0):

    def q_equation(q0):
        xt = (q0/W)-np.sqrt(((q0/W)**2)-((K*(1+a)*(z0)**2)/(W*(a)**2)))
        ht = np.sqrt(((2*q0*xt)-(W*xt**2))/(K*(1+a)))
        eq3 = h_i - np.sqrt(((2/K)*(x_i-xt)*(q0-(W/2)*(x_i+xt))+(ht+z0)**2)-z0)
        return eq3
    
    initial_guess = 0.9

    q_solution = fsolve(q_equation, x0=initial_guess)
    q_solution = float(q_solution)
    q_solution = round(q_solution,5)
    return(q_solution)
    

count = 0

for i in sample_scaled:

    print(" ")

    a = i[0]
    q0 = i[1]
    z0 = i[2]
    W = i[3]
    K = i[4]
    
    a = round(a,3)
    q0 = round(q0,3)
    z0 = round(z0,3)
    W = round(W,4)
    K = round(K,3)

    #Number of Runs
    count = int(count)
    count = count + 1
    count = str(count)
    print("Run Number " + count)

    #Partial Function of Flux-Controlled System 
    flux_controlled = partial(aquifer_model, a,W,K,3000,2000,q0)

    #Partial Function of Head-Controlled System
    head_controlled = partial(aquifer_model, a,W,K,3000,2000)

    #Inital Eval of Both Systems at Base SL
    print(" ")
    print("SLR: 0m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0))
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0))
    print(" ")
    h_i = input("Input h_inland Value for SLR = 0m: ")
    h_i = float(h_i)

    #Partial Function of q0 with Respect to SLR
    new_q = partial(new_q_formula,40,0.0002,10,2000,h_i)

    #Eval of Eval of Both Systems at 0.25m SLR
    print(" ")
    print("SLR: 0.25m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 0.25))
    q0 = new_q(z0 + 0.25)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 0.25))

    #Eval of Eval of Both Systems at 0.5m SLR
    print(" ")
    print("SLR: 0.5m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 0.5))
    q0 = new_q(z0 + 0.5)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 0.5))

    #Eval of Eval of Both Systems at 0.75m SLR
    print(" ")
    print("SLR: 0.75m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 0.75))
    q0 = new_q(z0 + 0.75)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 0.75))

    #Eval of Eval of Both Systems at 1m SLR
    print(" ")
    print("SLR: 1m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 1))
    q0 = new_q(z0 + 1)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 1))

    #Eval of Eval of Both Systems at 1.25m SLR
    print(" ")
    print("SLR: 1.25m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 1.25))
    q0 = new_q(z0 + 1.25)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 1.25))

    #Eval of Eval of Both Systems at 1.5m SLR
    print(" ")
    print("SLR: 1.5m")
    print("FLUX-CONTROLLED SYSTEM:")
    print(flux_controlled(z0 + 1.5))
    q0 = new_q(z0 + 1.5)
    print("HEAD-CONTROLLED SYSTEM:")
    print(head_controlled(q0,z0 + 1.5))
    print(" ")



        

    



