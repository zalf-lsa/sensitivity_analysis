#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys

sys.path.append('../..')

import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D
import analysing_lib

treshold_mu=0.02

width=4
height=4

path = "../morris/runs/TDSA-2015-03-20/Ascha/"


infos = ["crop-Winter wheat-min-step6_range20_startvector20", "Winterweizen"]


def main():

    parameter_name_map = analysing_lib.get_parameter_map()

    parameter_map = {}
    dtype = [('no', 'a15'), ('name', 'a40')]
    header = ['No', 'Parameter']

    directory = infos[0]
    crop_name = infos[1]

    print
    print crop_name
    
  
    mean_file = open(path + directory + "/morris_mean.txt", "rb")
    mean_csv = csv.reader(mean_file, delimiter="\t")

    row_no = 0
    name = []
    no = []
    mu = []


    for mean in mean_csv:
        # first row contains output names
        if (row_no==0):
            for item_index, item in enumerate(mean):
                # ignore first column, only store output names of the other columns
                if (item_index>0):
                    header.append(item)

        else:
            fname = None
            pname = None
            pno = None
            for item_index, item in enumerate(mean):

                if (item_index == 0):
                    # get parameter name that is stored in first column
                    fname = item
                    pname = parameter_name_map[fname].name
                    pno = parameter_name_map[fname].no
                    if (pname not in parameter_map):        

                        tlist = [pno, pname]
                        tlist.extend(list(np.zeros(len(mean)-1)))
                        print "Adding new list to map", pname, tlist
                        parameter_map[pname] = tlist
                        name.append(str(parameter_name_map[fname].name))
                        no.append(str(parameter_name_map[fname].no))
                else:            
                    print item
                    # get SA values for each output
   
                    tlist = parameter_map[pname]
                    print tlist, len(tlist), item_index+1
                    tlist[item_index+1] = float(item)
                    parameter_map[pname] = tlist
            

            
        row_no += 1



    scaled_parameter_map = scale_parameter_effects(parameter_map, len(mean)+1)

    csv_file = open(path + directory + "/" + "morris_ranking.csv", 'wb')
    csv_writer = csv.writer(csv_file, delimiter='&')
    csv_writer.writerow(header)

    for pname, plist in scaled_parameter_map.iteritems():
        csv_writer.writerow(plist)

    csv_file.close()



def scale_parameter_effects(parameter_map, max_index):

    scaled_parameter_map = {}

    for item_index in range(2,max_index):
        names = []
        mu_list = []
        for pname, plist in parameter_map.iteritems():

            if (pname not in scaled_parameter_map):        
                scaled_list = plist[0:2]
                scaled_list.extend(list(np.zeros(max_index-1)))
                scaled_parameter_map[pname] = scaled_list                

            names.append(pname)
            mu = float(plist[item_index])
            mu_list.append(mu)

        mu_list= np.divide(mu_list, max(mu_list))   

        for name_index, name in enumerate(names):
            scaled_list = scaled_parameter_map[name]
            mu = mu_list[name_index]
            if (mu<treshold_mu):
                mu = 0.0
            scaled_list[item_index] = mu

            scaled_parameter_map[name] = scaled_list
            
        print
    return scaled_parameter_map

########################################################
# start
########################################################




main()


