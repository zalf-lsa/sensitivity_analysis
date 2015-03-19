#!/usr/bin/python
# -*- coding: UTF-8-*-

import numpy as numpy


""" 
Calculates the root mean square error of two lists    
"""
def calc_rmse(p, o):
    
    # p = predicted or modelled data
    # o = observed, measured data
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    mse = calc_mse(p,o)
    rmse = numpy.sqrt(mse)
    
    return (rmse)


def calc_mse(p,o):
    # p = predicted or modelled data
    # o = observed, measured data
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    n = len(p)
    
    diff = numpy.subtract(p,o)
    squared = numpy.multiply(diff,diff)
    #print "squared: ", squared, "\tcumsum: ", numpy.cumsum(squared)
    mse = (numpy.cumsum(squared)[n-1]) / n
    
    return mse
""" 
Calculates the normalised mean absolute error of two lists    
"""
def calc_nrmse(p,o):
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    # p = predicted or modelled data
    # o = observed, measured data
    mse = calc_mse(p,o)
    var_o = numpy.var(o)
    
    if (var_o == 0):
        var_o = numpy.mean(o)
        
    if (var_o == 0):
        return -1
    
    nrmse = numpy.sqrt(mse/var_o) 
    return (nrmse)

""" 
Calculates the mean absolute error of two lists    
"""
def calc_mae(p,o):
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    # p = predicted or modelled data
    # o = observed, measured data
    n = len(p)
    diff = numpy.subtract(p,o)
    abs_diff = numpy.absolute(diff)
    mae = (numpy.cumsum(abs_diff)[n-1]) / n
    
    return (mae)


""" 
Calculates the normalised mean absolute error of two lists    
"""
def calc_nmae(p,o):
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    # p = predicted or modelled data
    # o = observed, measured data  
    #n = len(o)
    #diff = abs(numpy.subtract(p,o))
    #mae = (numpy.cumsum(diff)[n-1]) / n
    
    mae = calc_mae(p,o)
    mean_o = numpy.mean(o)
    mean_p = numpy.mean(p)
    if (mean_o==0):
        return(-1)
   
    
    nmae = mae/mean_p

    return (nmae)


""" 
Calculates the modelling efficiency of two lists    
"""
def calc_ef(p,o):
    # p = predicted or modelled data
    # o = observed, measured data
    if (len(p) == 0 or len(o) == 0):
        return -1
  
    n = len(p)
    x = numpy.cumsum((numpy.subtract(p,o))**2)[n-1]
    y = numpy.cumsum((o-numpy.mean(o))**2)[n-1]
    ef = 1-(x/y)

    return (ef)

