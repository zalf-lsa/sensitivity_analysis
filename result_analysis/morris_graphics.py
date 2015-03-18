#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys
import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D

xlim=0.00
ylim=0.00

width=12
height=18

path = "2012-12-19_15-59-20/"

crop_indices = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]

def main(crop_index):
    print crop_index
    mean_file = open(path + "morris_mean.csv", "rb")
    std_file = open(path + "morris_std.csv", "rb")

    mean_csv = csv.reader(mean_file, delimiter=";")
    std_csv = csv.reader(std_file, delimiter=";")


    row_no = 0
    name = []
    no = []
    x = []
    y = []

    crop_name = ""   
    for (mean, std) in zip(mean_csv, std_csv):

        if (row_no > 0):
            m = 0
            s = 0

            if (len(mean[crop_index])>0):
                m = float(mean[crop_index])
            if (len(std[crop_index])>0):
                s = float(std[crop_index])

            if (m>xlim or s>ylim):
                name.append(mean[1])
                no.append(str(mean[0]))
                x.append(m)
                y.append(s)
        else:
            crop_name = mean[crop_index]

        row_no += 1

    x= np.divide(x, max(x))
    y= np.divide(y, max(y))

    ax = None
    if (crop_index == 2):
        ax = fig.add_subplot(5,3,1)
    elif (crop_index == 3):
        ax = fig.add_subplot(5,3,2)
    elif (crop_index == 4):
        ax = fig.add_subplot(5,3,3)
    elif (crop_index == 5):
        ax = fig.add_subplot(5,3,4)
    elif (crop_index == 6):
        ax = fig.add_subplot(5,3,5)
    elif (crop_index == 7):
        ax = fig.add_subplot(5,3,6)
    elif (crop_index == 8):
        ax = fig.add_subplot(5,3,7)
    elif (crop_index == 9):
        ax = fig.add_subplot(5,3,8)
    elif (crop_index == 10):
        ax = fig.add_subplot(5,3,9)
    elif (crop_index == 11):
        ax = fig.add_subplot(5,3,10)
    elif (crop_index == 12):
        ax = fig.add_subplot(5,3,11)
    elif (crop_index == 13):
        ax = fig.add_subplot(5,3,12)
    elif (crop_index == 14):
        ax = fig.add_subplot(5,3,13)
    elif (crop_index == 15):
        ax = fig.add_subplot(5,3,14)
    elif (crop_index == 16):
        ax = fig.add_subplot(5,3,15)

    ax.plot(x, y, 's', color = "#c85a2e")
    ax.set_title(crop_name)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    plot.xlabel(r'$\mu^*$')
    plot.ylabel(r'$\sigma$')

    for (x_pos, y_pos, n) in zip(x,y,no):
        if (x_pos>0.18 or y_pos>0.18):
            print x_pos, y_pos, n
            plot.text(x_pos, y_pos+0.04, r'' + str(n), horizontalalignment='center', verticalalignment='center', fontsize=8)

    

########################################################
# start
########################################################
plot.rcParams['figure.subplot.left'] = 0.05
plot.rcParams['figure.subplot.right'] = 0.98
plot.rcParams['figure.subplot.bottom'] = 0.06
plot.rcParams['figure.subplot.top'] = 0.95
plot.rcParams['savefig.dpi'] = 200      # figure dots per inch
plot.rcParams['figure.subplot.wspace'] = 0.2    # the amount of width reserved for blank space between subplots
plot.rcParams['figure.subplot.hspace'] = 0.3    # the amount of height reserved for white space between subplots

fig = plot.figure(figsize=(width,height))

for crop_index in crop_indices:
    main(crop_index)

fig.savefig("morris.png", dpi_value=200)
