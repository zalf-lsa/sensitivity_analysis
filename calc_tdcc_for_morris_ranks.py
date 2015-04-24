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
               ["morris_ranking_primyield.csv", "TDCC_matrix-primyield.png"],
               ["morris_ranking_biomNContent.csv","TDCC_matrix-biomNContent.png"],
               ["morris_ranking_ETaCrop.csv","TDCC_matrix-Eta_Crop.png"],
               ["morris_ranking_TranspCrop.csv","TDCC_matrix-Transp_Crop.png"],
               ["morris_ranking_anthesis.csv","TDCC_matrix-Anthesis.png"],
               ["morris_ranking_moist90.csv","TDCC_matrix-moist90.png"],
               ["morris_ranking_nmin.csv","TDCC_matrix-nmin90.png"]
               #["morris_ranking_abbiom.csv","TDCC_matrix-ABG.png"],
               #["morris_ranking_abbiomNcontent.csv","TDCC_matrix-ABG_N_content.png"],                              
               #["morris_ranking_yearlyGWrecharge.csv","TDCC_matrix-yearly_GW_Recharge.png"],
               #["morris_ranking_yearlyLeachN.csv","TDCC_matrix-yearly_N_Leaching.png"],
               #["morris_ranking_corg.csv","TDCC_matrix-Corg.png"]
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
Creates an image with the visualisation of
the TDCC matrix.
"""    
def save_tdcc_matrix_figure(tdcc_array, header, filename):

    plot.rcParams['figure.subplot.left'] = 0.2
    plot.rcParams['figure.subplot.right'] = 0.9
    plot.rcParams['figure.subplot.bottom'] = 0.1
    plot.rcParams['figure.subplot.top'] = 0.99
    plot.rcParams['savefig.dpi'] = 200      # figure dots per inch
    #plot.rcParams['figure.subplot.wspace'] = 0.6    # the amount of width reserved for blank space between subplots
    #plot.rcParams['figure.subplot.hspace'] = 0.2    # the amount of height reserved for white space between subplots
    plot.rcParams['legend.fontsize'] = 11
    plot.rcParams['xtick.major.pad'] = 10
    plot.rcParams['ytick.major.pad'] = 10
    plot.rc('text', usetex=True)

    width=5
    height=4.5
    font_size = 14

    size =len(header)

    fig = plot.figure(figsize=(width,height)) #,frameon=False
    ax = fig.add_subplot(111)
    #extent = -delta,len(crop_names)-delta, -delta, len(parameter_names)-delta
    im1 = plot.imshow(tdcc_array, cmap=plot.cm.Greens,vmin=0.0, vmax=1.0 ) # gist_yarg    #plot.cm.Paired,  extent=extent,
    im1.set_interpolation('none')

    x_pos = range(0,size,1)
    y_pos = range(0,size,1)
    delta = 0.5

    plot.xticks(x_pos, header, fontsize=font_size, rotation=45, ha='center')
    plot.yticks(y_pos, header, fontsize=font_size, va='center')

    ax.set_xlim(-delta, size-delta)
    ax.set_ylim(-delta, size-delta)

    for x in x_pos:
        ax.axvline(x=x-delta, ls='solid', color='#ffffff', lw=1.5)
    for y in y_pos:
        ax.axhline(y=y-delta, ls='solid', color='#ffffff', lw=1.5)

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = plot.colorbar(im1, cax=cax)
    #cbar = plot.colorbar()
    #cbar.set_clim(0, 1.0)
    #cbar.ax.set_title('TDCC', fontsize=10)
    
    cbar.ax.tick_params(labelsize=10) 
    fig.savefig(filename, dpi_value=200)

    del fig
    del im1




    
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