#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys
import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D
sys.path.append('..')
import analysing_lib

xlim=0.00
ylim=0.00

width=4.5
height=4.5

path = "2013-01-29-neu/"

folders = [["crop-1", "Winterweizen", 1, "winter_wheat_c.png"],
           ["crop-19", "Wintertriticale", 19, "winter_triticale_c.png"],
            ["crop-2",  "Wintergerste",2, "winter_barley_c.png"],
            ["crop-9", "Winterraps",9, "winter_rape_c.png"],
            ["crop-4", "Sommergerste",4, "spring_barley_c.png"],
            ["crop-23", "Sommertriticale",23, "spring_triticale_c.png"],
            ["crop-10", "Zuckerrübe",10, "sugar_beet_c.png"],
            ["crop-7", "Mais",7, "maize_c.png"],
            ["crop-18", "Sudangras",18, "sudan_gras_c.png"],
            ["crop-12", "Phacelia",12, "phacelia_c.png"],
            ["crop-13", "Kleegras",13, "clover_c.png"],
            ["crop-14", "Luzernegras",14, "alfalfa_c.png"],
            ["crop-16", "Weidelgras",16, "rye_gras_c.png"],
            ["crop-22", "Hafer",22, "oat_c.png"]
            ]


def main(info_list):

    directory = info_list[0]
    crop_name = info_list[1]
    crop_id = info_list[2]
    png_file = info_list[3]
    
    parameter_map = analysing_lib.get_parameter_map()

    result_file = open(path + directory + "/tsi-abBiom.csv", "rb")
    result_csv = csv.reader(result_file, delimiter=",")

    parameter_no = []
    parameter_name = []
    si_list = []
    tsi_list = []
    row_no = 0
    for result in result_csv:
        if (row_no>1):

            parameter = parameter_map[str(result[0]) + ".txt"]
            si = float(result[1])
            tsi = float(result[2])

            parameter_no.append(parameter.no)
            parameter_name.append(parameter.name)
            si_list.append(si)
            tsi_list.append(tsi)
        row_no += 1


    plot.rcParams['figure.subplot.left'] = 0.3
    plot.rcParams['figure.subplot.right'] = 0.95
    plot.rcParams['figure.subplot.bottom'] = 0.12
    plot.rcParams['figure.subplot.top'] = 0.98
    plot.rcParams['savefig.dpi'] = 200      # figure dots per inch
    plot.rcParams['figure.subplot.wspace'] = 0.6    # the amount of width reserved for blank space between subplots
    plot.rcParams['figure.subplot.hspace'] = 0.2    # the amount of height reserved for white space between subplots
    plot.rcParams['legend.fontsize'] = 9

    fig = plot.figure(figsize=(width,height))
    ax = fig.add_subplot(111)
    #ax.set_title(crop_name)

    bars_pos = range(len(parameter_no),0,-1)

    print crop_name, len(parameter_no),len(tsi_list), len(si_list), bars_pos
    tsi = ax.barh(bars_pos, tsi_list, align='center', color='#000000', height=0.5) # #ffac0d
    si = ax.barh(bars_pos, si_list, align='center', color='#a0ff54', height=0.5) # #58822e
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0,len(tsi_list)+1)
    plot.yticks(bars_pos, parameter_name, fontsize='9')
    plot.xlabel(r'${TS_i}$')
    
    
    print crop_id
    if (crop_id in (23, 14,22, 23, 7) ):
        #l = ax.legend([tsi[0], si[0]],[r' Main effect', r' Interactions'], loc='lower right')
        l = ax.legend([si[0], tsi[0]],[r' Haupteffekt', r' Interaktionseffekt'], loc='lower right')
    
    fig.savefig(png_file, dpi_value=200)
    

########################################################
# start
########################################################


for folder in folders:
    main(folder)


