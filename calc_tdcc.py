#! /usr/bin/python
# -*- coding: UTF-8-*-

import os
import csv
import re
import numpy
import sys
import pandas
import scipy.stats as ss

ranking_file = "D:/daten_specka/ZALF/publications/2015-Time-dependent-SA/data/morris_ranking_primyield.csv"

######################################################
######################################################
######################################################

"""
Main method
"""
def main():

    parameter_ids, parameter_names, rank_array = get_ranking_map(ranking_file)
    print(rank_array)
    
    # test array
    #rank_array=numpy.array([1,3,2,2,1,3,3,2,1,4,4,4]).reshape(4,3)
    
    savage_scores = get_savage_score_array (rank_array)
    
    calc_TDCC(savage_scores)

######################################################
######################################################
######################################################

"""
Calculates savage scores the passed array
"""
def get_savage_score_array(r_array):

    print("Calc savage scores of parameter ranks")
    
    # number of different rankings
    different_ranking_count = len(r_array[0])
    
    # initialise result array with zeroes
    savage_score_array = numpy.zeros((r_array.shape))
    
    print ("Found different rankings: ", different_ranking_count)

    for rank_index in range(different_ranking_count):                
        savage_score_array[:,rank_index] = calc_savage_scores(r_array[:,rank_index])        
    
    print("\nSavage scores:\n", savage_score_array,"\n")
    
    return savage_score_array

######################################################
######################################################
######################################################

"""
Calculates savage scores of a ranked list
"""
def calc_savage_scores(ranks):
    
    n=len(ranks)
    savage_scores = []
    
    for rank in ranks:
        savage = 0
        
        for j in range(rank, n+1):
            savage += 1/j
        
        savage_scores.append(savage) 
    
    return savage_scores

######################################################
######################################################
######################################################

def calc_TDCC(savage_score_array):
    print("Calculation of TDCC")

    # formula extracted from "Iman, R.L., Conover, W.J. 1987: A Measure of Top-Down-Correlation. Technometrics Vol.29 (3)"
    
    n = len(savage_score_array[:, 0])
    b = len(savage_score_array[0])
        
    print("Number of parameters:\t\t", n)
    print("Number of different rankings:\t", b)

    sum_zaehler = 0
    S1 = 0

    for i in range(0,n):
        Si = numpy.sum(savage_score_array[i])
        sum_zaehler += Si**2         
    
        S1 += 1/(i+1)

    zaehler = sum_zaehler - b**2*(n)
    nenner = b**2 * (n-S1)
    
    TDCC= zaehler/nenner

    print("\nTDCC:", TDCC)

    return

    
######################################################
######################################################
######################################################

"""
Reads in MONICA files and saves the results into a list of maps
"""
def get_ranking_map(filename):
         
    parameter_ids = []
    parameter_names = []
    rank_array = []


    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter=";")
        
        # skip header
        next(reader) 

        # store parameter id, name and ranks into separate lists
        for row in reader:      
               
            parameter_ids.append(row[0])
            parameter_names.append(row[1])

            ranks = []
            for rank in row[2:]:
                ranks.append(int(rank))
            rank_array.append(ranks)

    
    rank_array = numpy.array(rank_array)
    
    
    return parameter_ids, parameter_names, rank_array

######################################################
######################################################
######################################################

main()