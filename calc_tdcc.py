#! /usr/bin/python
# -*- coding: UTF-8-*-

import os
import csv
import re
import numpy
import sys


ranking_file = "D:/daten_specka/ZALF/publications/2015-Time-dependent-SA/data/morris_ranking_primyield.csv"

######################################################
######################################################
######################################################

"""
Main method
"""
def main():

    #parameter_ids, parameter_names, rank_array = get_ranking_map(ranking_file)
    #print(rank_array)
    
    #array = rank_array[:,[0,1]]

    array=numpy.array([1,1,1,2,4,3,3,2,4,4,3,2]).reshape(4,3)
    print("Array:", array)

    savage_scores = get_savage_score_array (array)
    print("Savage scores:\t", savage_scores)
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
         
    print (savage_score_array)

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
        
        start = n-rank+1
        x=list(range(start, n+1))
        
        for j in x:
            savage += 1/j
        savage_scores.append(1-savage)
    
    return savage_scores

######################################################
######################################################
######################################################

def calc_TDCC(savage_score_array):
    print("Calc TDCC")
    
    k = len(savage_score_array[:, 0])
    ranking_count = len(savage_score_array[0])
        
    
    sum_nenner = 0
    for i in range(1,k+1):
        sum_nenner += (1/(i+1))
        

    sum_savage_scores = numpy.sum(savage_score_array)
    print("Sum:", sum_savage_scores, sum_nenner)
    
    zaehler = sum_savage_scores - ((ranking_count**2)*k)
    nenner = ranking_count**2 * (k - sum_nenner)

    TDCC = zaehler / nenner
    print("ns_a^2*k", (ranking_count**2)*k)
    print("zaehler:", zaehler)
    print("nenner:", nenner)
    print("TDCC:", TDCC)

    


        


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