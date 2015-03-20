#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys


import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D
import analysing_lib

xlim=0.0
ylim=0.0

width=4
height=4

path = "../morris/runs/TDSA-2015-03-20/Ascha/"


folders = [["crop-Winter wheat-min-step6_range20_startvector20", "Winterweizen"]]

m_color = "#c85a2e"

output_anzahl = 11

def main(infos):
    #print infos

    directory = infos[0]
    crop_name = infos[1]

    print dir, crop_name
    parameter_map = analysing_lib.get_parameter_map()


    for output_index in range(1,output_anzahl+1):

        print output_index, "/", output_anzahl
        plot.rcParams['figure.subplot.left'] = 0.12
        plot.rcParams['figure.subplot.right'] = 0.98
        plot.rcParams['figure.subplot.bottom'] = 0.12
        plot.rcParams['figure.subplot.top'] = 0.98
        plot.rcParams['savefig.dpi'] = 200      # figure dots per inch
        plot.rcParams['figure.subplot.wspace'] = 0.2    # the amount of width reserved for blank space between subplots
        plot.rcParams['figure.subplot.hspace'] = 0.3    # the amount of height reserved for white space between subplots



        fig = plot.figure(figsize=(width,height))



        mean_file = open(path + directory + "/morris_mean.txt", "rb")
        std_file = open(path + directory + "/morris_std.txt", "rb")

        mean_csv = csv.reader(mean_file, delimiter="\t")
        std_csv = csv.reader(std_file, delimiter="\t")



        row_no = 0
        name = []
        no = []
        x = []
        y = []
        output_name = None
        for (mean, std) in zip(mean_csv, std_csv):

            if (row_no  == 0):
                output_name = mean[output_index]
            else:
                m = 0
                s = 0
                
                if (len(mean[0])>0):
                    fname = mean[0]
                    m = float(mean[output_index])
                if (len(std[0])>0):
                    s = float(std[output_index])

     
                if (m>xlim or s>ylim):
                    name.append(str(parameter_map[fname].name))
                    no.append(str(parameter_map[fname].no))
                    x.append(m)
                    y.append(s)


            row_no += 1

        if (len(x)>0):
            x= np.divide(x, max(x))

        if (len(y)>0):
            y= np.divide(y, max(y))

        ax = fig.add_subplot(1,1,1)

        ax.plot(x, y, 's', color = m_color)
        # ax.set_title(crop_name)
        ax.set_xlim(0, 1.1)
        ax.set_ylim(0, 1.1)
        plot.xlabel(r'$\mu^*$')
        plot.ylabel(r'$\sigma$')

        for (x_pos, y_pos, n) in zip(x,y,no):
            if (x_pos>0.2 or y_pos>0.2):
                #print x_pos, y_pos, n
                plot.text(x_pos, y_pos+0.04, r'' + str(n), horizontalalignment='center', verticalalignment='center', fontsize=11)

        fig.savefig(path + directory + "/" + output_name + ".png", dpi_value=200)






########################################################
# start
########################################################




for folder in folders:
    main(folder)


