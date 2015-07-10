#! /usr/bin/python
# -*- coding: UTF-8-*-

import os
import csv
import re
import numpy
import sys
import pandas
import scipy.stats as ss
import matplotlib.pyplot as plot
from matplotlib.lines import Line2D
import sa_functions as saf
from mpl_toolkits.axes_grid1 import make_axes_locatable

path = "D:/daten_specka/ZALF/publications/2015-Time-dependent-SA/data/morris/"
output_path = "D:/daten_specka/ZALF/publications/2015-Time-dependent-SA/tex/"


input_configs = [
               #["morris_ranking_primyield.csv", "TDCC_matrix-primyield.png"],
               #["morris_ranking_biomNContent.csv","TDCC_matrix-biomNContent.png"],
               #["morris_ranking_ETaCrop.csv","TDCC_matrix-Eta_Crop.png"],
               #["morris_ranking_TranspCrop.csv","TDCC_matrix-Transp_Crop.png"],
               #["morris_ranking_anthesis.csv","TDCC_matrix-Anthesis.png"],
               #["morris_ranking_moist90.csv","TDCC_matrix-moist90.png"],
               #["morris_ranking_nmin.csv","TDCC_matrix-nmin90.png"]
               ["morris_ranking_winter_wheat.csv","TDCC_matrix_winter_wheat.png"],
               ["morris_ranking_winter_rape.csv","TDCC_matrix_winter_rape.png"],
               ["morris_ranking_clover.csv","TDCC_matrix_clover.png"],
               ["morris_ranking_winter_wheat.csv","TDCC_matrix_winter_wheat.png"],
               ["morris_ranking_winter_wheat.csv","TDCC_matrix_winter_wheat.png"],
               ["morris_ranking_winter_wheat.csv","TDCC_matrix_winter_wheat.png"],
               
               
               ]


######################################################
######################################################
######################################################

"""
Main method
"""
def main():

    for input in input_configs:
        print("#######################################################")
        print("Analysing file", input[0])
        ranking_file = input[0]
        matrix_figure_filename = input[1]

        parameter_ids, parameter_names, rank_array, header = get_ranking_map(path + ranking_file)
        
        # test array
        #rank_array=numpy.array([1,3,2,2,1,3,3,2,1,4,4,4]).reshape(4,3)

        # Remove rows with lowest equal ranks
        corrected_rank_array = saf.correct_array(rank_array)
        
        # calculates TDCC value for all rankings
        saf.calc_TDCC(corrected_rank_array)
        
        # calculate several TDCCs for a pair of rankings
        print("Calc TDCC matrix")
        
        size = len(rank_array[0])        
        tdcc_array = numpy.zeros((size,size))

        

        for x in range(size):
            for y in range(0,size):            
                #if (x<y):
                #    continue
        
                tdcc= None
                

                if (x==y):
                    tdcc = 1
                else:
                    s_array = numpy.zeros((len(rank_array[:,1]),2))              
                    s_array[:,0] = rank_array[:,x]
                    s_array[:,1] = rank_array[:,y]

                    corrected_s_array = saf.correct_array(s_array)
                    tdcc = saf.calc_TDCC(corrected_s_array)  

                print("x:",x,"\ty:",y,"\tTDCC:", tdcc)          
                
                tdcc_array[y,x] = tdcc

        # save the matrix of TDCCs into an image
        save_tdcc_matrix_figure(tdcc_array, header, output_path + matrix_figure_filename)

        





    
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
    header = None

    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter=";")
        
        # skip header
        header = next(reader) 

        # store parameter id, name and ranks into separate lists
        for row in reader:      
               
            parameter_ids.append(row[0])
            parameter_names.append(row[1])

            ranks = []
            for rank in row[2:]:
                ranks.append(float(rank))
            rank_array.append(ranks)

    
    rank_array = numpy.array(rank_array)    
    
    return parameter_ids, parameter_names, rank_array, header[2:]

######################################################
######################################################
######################################################

main()