#!/usr/bin/python
# -*- coding: UTF-8-*-

import numpy
import matplotlib.pyplot as plot
import sys
import os
import matplotlib.pyplot as plot
import matplotlib as mpl
from matplotlib import cm
import re
import csv
import datetime

sys.path.append('..')
sys.path.append('../../../../user-libs/lib')
sys.path.append('/media/san1_data1/data1/specka/devel/master/monica/python')

# import monica
#import station


#######################################################
# Methods:
#
# processResult(result):
# calc_output_year(start_date, end_date):
# get_year(date_string):
#
# get_individual_range(rank, nodes, grid_size):
#
# debug(rank=0, message=""):
#
# plot_array(result_grid, filename):
# plot_histogramm(result_array, filename):
# plot_grid(grid, filename, grid_type=""):
#
#######################################################

"""
Calculates mean and standard deviation for all outputs for 
CCGermany simulation. Results are sequentially stored in a list
for further processing
"""
def processResult(result):
    
    result_ids = monica.CCGermanyResultIds()
     
    processed_results=[]
    
    for id in result_ids:
        
        tlist = result.getResultsById(id)
        mean = numpy.mean(tlist)
        std = numpy.std(tlist)
        processed_results.append(mean)
        processed_results.append(std)     
        processed_results.append(list(tlist))   
    
    return processed_results

#######################################################

"""
Calculates the middle year of given range
"""
def calc_output_year(start_date, end_date):
    output_year = (get_year(start_date) + get_year(end_date))/2
    return output_year

#######################################################

"""
Returns year of passed date string of format "yyyy-mm-dd"
"""
def get_year(date_string):
    return int(date_string[0:4])


#######################################################

"""
Calculates start and end index of a two dimensional grid,
that can be distributed among different nodes
"""
def get_individual_range(rank, nodes, grid_size):

    grid_size = 50
    start = 1500000
    #start = 0
    max_cells_per_node = grid_size
    if (nodes>1):
        # divide grid into separate parts for cluster computing
        # calculate number of cells that must be calculated per node
        max_cells_per_node = grid_size // (nodes-1)        

    # get start and end index of range, that must be calculated on a node
    start_index = rank*max_cells_per_node
    end_index = ((rank+1) * max_cells_per_node)-1
    if (rank == nodes-1 or nodes == 0):
        end_index = grid_size-1

    return (start + start_index, start + end_index, max_cells_per_node)


#######################################################

"""
Simple output message with equal formatting of rank and following message
"""
def debug(rank=0, message=""):
    print "Rank:", rank, "\t" + message



#####################################################################

"""
Plots histogramm of result
"""
def plot_histogramm(array, filename):

    y=numpy.array(array)
    new_y = y.compress((y!=0).flat)
    if (len(new_y)>0):
        fig = plot.figure() 
        bx = fig.add_subplot(111)
        n, bins, patches = bx.hist(new_y, 100, normed=0, facecolor='green', alpha=0.75)
        plot.xlabel('Tavg')
        plot.ylabel('Frequency')
        plot.grid(True)
        plot.savefig(filename)

#####################################################################


"""
Plots resulting grid
"""
def plot_array(result_grid, filename,img_norm=None, colormap = plot.cm.jet):

    fig = plot.figure() 
    ax = fig.add_subplot(111)

    if (img_norm != None):
        normalize_range = mpl.colors.Normalize(vmin=img_norm[0], vmax=img_norm[1])
    else:
        normalize_range = None
    
    im = ax.imshow(result_grid, cmap = colormap, norm=normalize_range)
    cbar = plot.colorbar(im)
    plot.savefig(filename)

#####################################################################
"""
Plot a grid map
"""
def plot_grid(grid, filename, img_norm=None, colormap = plot.cm.jet):

    grid_rows = grid.nrows
    grid_cols = grid.ncols

    # saving original maps into a png for later analysis
    # saving height and buek soil types
    array = numpy.zeros((grid_rows,grid_cols))

    hist_array = []

    index = 0
    for r in range(grid_rows):
        for c in range(grid_cols):

            val = grid.get_xy(r,c)

            # set -9999 values to None, so they won't be plotted
            if (val==-9999):
                val = None
            array[r,c] = val

            # reset value to zero because of plotting of histogram
            if (val != None):
                hist_array.append(val)

            index += 1

    plot_array(array, filename, img_norm, colormap)
    
    plot_histogramm(hist_array, filename[:-4] + "-hist" + filename[-4:])




""" 
Returns a list with all files that are located in the 
directory specified by 'path'
"""
def getFilesInDirectory(path):
    directory_list = os.listdir(path)    
    files = []
    
    for item in directory_list:
        if os.path.isfile(path + '/' + item):
            files.append(item)
    
    return (files)

########################################################

def test_location(filename, location):
    # separately analyse for each location
    regex = re.compile(location)
    result = regex.search(filename)
    return result

########################################################

def test_anlage(filename, anlage):
    # only calculate statistics for files of anlage 1 und 2
    regex = re.compile("anlage-"+str(anlage))
    result = regex.search(filename)
    return result

#############################################################
#############################################################
#############################################################

def get_list(smout, key):

    values = []

    for index, map in enumerate(smout):
        if (key == "Datum"):
            d = map["Datum"].split("/")
            date = datetime.date(int(d[2]), int(d[1]), int(d[0]))               
            values.append(date)
        elif (key == "Yield"):
            y = float(map[key])
            y = y/100.0     # Umwandlung von kg / ha in dt / ha
            values.append(y)
        else:
            values.append(float(map[key]))
        
    return values

#############################################################
#############################################################
#############################################################

"""
Reads in MONICA files and saves the results into a list of maps
"""
def read_monica_file(file):
         
    data = csv.reader(open(file), delimiter="\t")  
    # Read the column names from the first line of the file  
    fields = data.next()
    
    # skip the line with units
    data.next()
    
    map_list = []  
    for row in data:  
        # Zip together the field names and values  
        items = zip(fields, row)  
        item = {}  
        # Add the value to our dictionary  
        
        for (name, value) in items:
            #print name, value  
            item[name.strip()] = value.strip()
        map_list.append(item)
    
    return map_list
