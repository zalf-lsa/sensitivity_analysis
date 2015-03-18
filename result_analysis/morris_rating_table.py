#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys

sys.path.append('..')

import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D
import analysing_lib

xlim=0.02

width=4
height=4

path = "2013-01-25-all/"


folders = [["crop-Winterweizen-min-step5_range20_startvector40", "Winterweizen"],
            ["crop-Zuckerrube-min-step5_range20_startvector40",  "Zuckerrube"],
            ["crop-Sommergerste-min-step5_range20_startvector40", "Sommergerste"],
            ["crop-Mais-min-step5_range20_startvector40", "Mais"],
            ["crop-Kleegras-min-step5_range20_startvector40", "Kleegras"],
            ["crop-Winterraps-min-step5_range20_startvector40", "Winterraps"],
            ["crop-Hafer-min-step5_range20_startvector40", "Hafer"],
            ["crop-Luzernegras-min-step5_range20_startvector40", "Luzernegras"],
            ["crop-Phacelia-min-step5_range20_startvector40", "Phacelia"],
            ["crop-Sommertriticale-min-step5_range20_startvector40", "Sommertriticale"],
            ["crop-Sudangras-min-step5_range20_startvector40", "Sudangras"],
            ["crop-Weidelgras-min-step5_range20_startvector40", "Weidelgras"],
            ["crop-Wintergerste-min-step5_range20_startvector40", "Wintergerste"],
            ["crop-Wintertriticale-min-step5_range20_startvector40", "Wintertriticale"]
            ]

m_color = "#c85a2e"

def main():

    parameter_name_map = analysing_lib.get_parameter_map()

    parameter_map = {}
    dtype = [('no', 'a15'), ('name', 'a40')]
    header = ['No', 'Parameter']
    for crop_index,infos in enumerate(folders):

        directory = infos[0]
        crop_name = infos[1]

        print
        print crop_name
        
        dtype.append((crop_name, 'f4'))
        dtype.append((crop_name + '-rank', 'i8'))
        header.append(crop_name)
        header.append(crop_name+ '-rank')

        print dir, crop_name


        mean_file = open(path + directory + "/morris_mean.txt", "rb")
        mean_csv = csv.reader(mean_file, delimiter="\t")

        row_no = 0
        name = []
        no = []
        mu = []

        for mean in mean_csv:

            if (row_no > 0):
                m = 0
                if (len(mean[0])>0):
                    fname = mean[0]
                    m = float(mean[1])
                name.append(str(parameter_name_map[fname].name))
                no.append(str(parameter_name_map[fname].no))
                mu.append(m)
                if (crop_index == 0):
                    print fname, m
            row_no += 1


        # scale the mu's
        mu= np.divide(mu, max(mu))

        # remove elements from list that are below the treshold value
        mu2 = []
        name2 = []
        no2 = []
        for (pno, pname, pmu) in zip(no, name, mu):
            if (pmu>xlim):
                mu2.append(pmu)
                name2.append(pname)
                no2.append(pno)

        mu = mu2
        name = name2
        no = no2

        # save values into the parameter map
        for (pname,pno, pmu) in zip(name, no, mu):
            print pno, pname, pmu, crop_index
            if (pname not in parameter_map):        
                print "Adding new list to map"
                parameter_map[pname] = [[pno],list(np.zeros(len(folders)*2))]

            tlist = parameter_map[pname]
            pno = tlist[0]
            values = tlist[1]
            values[crop_index*2] = pmu
            parameter_map[pname] = [pno,values]
                


    morris_array = np.zeros((len(parameter_map.keys())),dtype=dtype)

    # analyse parameter map
    y_index = 0
    for pname, parameter in parameter_map.iteritems():
        print pname, parameter
        pno = parameter[0]
        values = parameter[1] 
        print pno, values

        array_list = pno        
        array_list.extend([pname])
        array_list.extend(values)
        array_list = tuple(array_list)
        morris_array[y_index] = array_list
        y_index += 1
        
    # create an array with parameter number, name and mu
    sorted_morris_array =  np.sort(morris_array, order='no')
    #print sorted_morris_array

    csv_file = open('morris_ranking.csv', 'wb')
    csv_writer = csv.writer(csv_file, delimiter='&')

    csv_writer.writerow(header)
    csv_writer.writerows(sorted_morris_array)

    csv_file.close()




########################################################
# start
########################################################




main()


