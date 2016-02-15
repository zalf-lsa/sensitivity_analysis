#!/usr/bin/python
# -*- coding: UTF-8-*-


import sys
import numpy
sys.path.append('../sensitivity_analysis')
import math
import random
from saparameter import SAParameter

max_harmonic = 4



"""

"""
def get_transformation(omega, s, param=None):

    phi = random.uniform(0, 2*math.pi)
    val = 0.5 + (1/math.pi) * math.asin(math.sin(omega*s))    

    # Normierung zum Wertebereich des Parameters
    if (param!=None):
        val = val*(float(param.getMax()) - float(param.getMin())) + float(param.getMin())
    
    return val

##############################################################################


def get_spectrum(result_fft):
    real = numpy.real(result_fft)
    imag = numpy.imag(result_fft)
    A = (numpy.add(numpy.multiply(real, real), numpy.multiply(imag,imag)))
    return A


def get_all_variance(result_fft, omegas):
    # overall variance used for D_y 
    real = numpy.real(result_fft)
    imag = numpy.imag(result_fft)
    A = (numpy.add(numpy.multiply(real, real), numpy.multiply(imag,imag)))

    D_y = sum(A[1:max(omegas)*4])
    #print "\nget_all_variance", D_y, max(omegas),max(omegas)*4, A[0]
    
    return D_y



"""
Calculation of first order indices of FAST either
for overall output or sinlge input
"""
def get_first_order_indices(result_fft, index=None, omegas=None):    

    #print "get_first_order for parameter ", index, omegas

    harmonics_input = get_harmonics_of_input(omegas, index)
    real = numpy.real(result_fft)
    imag = numpy.imag(result_fft)    
    A = (numpy.add(numpy.multiply(real, real), numpy.multiply(imag,imag)))/2.0

    #print "A: ", A

    D_input = sum(A[harmonics_input])    
    
    return D_input

##############################################################################

"""
Calculation of total indices of FAST
"""
def get_total_indices(result_fft, index=None, omegas=None):    
    

    max_omega_o = max(omegas) / (2 * max_harmonic)
   
    real = numpy.real(result_fft)
    imag = numpy.imag(result_fft)
    
    A = numpy.add(numpy.multiply(real, real), numpy.multiply(imag,imag))
    
    om = max(omegas)-1
    
    D_not_i = 2*(sum(A[1:om]))
    D_y = get_all_variance(result_fft, omegas)
     


    #print "get_total_indices for parameter", index
    #D_y = get_all_variance(result_fft, omegas)
    #print "Output variance", D_y
    #D_not_i = 0
    #D_total = 0
    
    #for o_index, o in enumerate(omegas):
        # add main effects of all inputs that are not x_i
        #first_order_o = get_first_order_indices(result_fft, o_index, omegas)           
        #D_total += first_order_o

        #if (o_index != index):        
            #D_not_i += first_order_o    
            #print "TOTAL", o_index, first_order_o, D_not_i,"\n"            

                
            
    TS_i =  1-(D_not_i / D_y)
    return (TS_i, D_y)

##############################################################################

"""
"""
def get_harmonics_of_input(omegas, index):
    harmonics = []
    for i in range(1,max_harmonic+1):
        harmonics.append(omegas[index]*i)
    return harmonics

##############################################################################

def get_ts_frequencies(omega_i, index, max_param):
    
    M=max_harmonic
    
    max_omega_o = omega_i / (2 * M)
    
    omegas = []
    omega_o = 1
    step = 16

    for i in range(max_param):
        
        if (omega_o>max_omega_o):
            omega_o=1
            
        if (i==index):
            # frequency of parameter of interest
            omegas.append(omega_i)
        else:
            # frequency of other parameters
            omegas.append(omega_o)
            omega_o+=step
    print "Omegas: ", omegas, "max_param", max_param
    return omegas
