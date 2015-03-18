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

width=16
height=11

path = "2013-01-29-neu/"

folders = [["crop-1", "Winterweizen", 1, "winter_wheat.png"],
           ["crop-19", "Wintertriticale", 19, "winter_triticale.png"],
            ["crop-2",  "Wintergerste",2, "winter_barley.png"],
            ["crop-9", "Winterraps",9, "winter_rape.png"],
            ["crop-4", "Sommergerste",4, "spring_barley.png"],
            ["crop-23", "Sommertriticale",23, "spring_triticale.png"],
            ["crop-22", "Hafer",22, "oat.png"],
            ["crop-10", u"Zuckerr\u00fcbe",10, "sugar_beet.png"],
            ["crop-7", "Mais",7, "maize.png"],
            ["crop-18", "Sudangras",18, "sudan_gras.png"],
            ["crop-12", "Phacelia",12, "phacelia.png"],
            ["crop-13", "Kleegras",13, "clover.png"],
            ["crop-14", "Luzernegras",14, "alfalfa.png"],
            ["crop-16", "Weidelgras",16, "rye_gras.png"]

            ]

font_size = '20'

def main():

    plot.rcParams['figure.subplot.left'] = 0.25
    plot.rcParams['figure.subplot.right'] = 0.99
    plot.rcParams['figure.subplot.bottom'] = 0.17
    plot.rcParams['figure.subplot.top'] = 0.95
    plot.rcParams['savefig.dpi'] = 200      # figure dots per inch
    plot.rcParams['figure.subplot.wspace'] = 0.6    # the amount of width reserved for blank space between subplots
    plot.rcParams['figure.subplot.hspace'] = 0.2    # the amount of height reserved for white space between subplots
    plot.rcParams['legend.fontsize'] = 11
    plot.rcParams['xtick.major.pad'] = 10
    plot.rcParams['ytick.major.pad'] = 10

    tsi_map = {}
    crop_names = []
    parameter_names = []
    parameter_pos = []
    
    for num_crop, info_list in enumerate(folders):


        directory = info_list[0]
        crop_name = unicode(info_list[1])
        crop_id = info_list[2]
        png_file = info_list[3]

        print "Analysing crop", num_crop, crop_name

        crop_names.append(crop_name)
        
        parameter_map = analysing_lib.get_parameter_map()

        result_file = open(path + directory + "/tsi-abBiom.csv", "rb")
        result_csv = csv.reader(result_file, delimiter=",")

        row_no = 0
        for result in result_csv:
            if (row_no>1):

                parameter = parameter_map[str(result[0]) + ".txt"]
                tsi = float(result[2])

                if (tsi>0.05):

                    if (parameter.name not in tsi_map):
                        tmp_list= []
                        while (len(tmp_list)<(num_crop)):
                            tmp_list.append(0.0)
                        tsi_map[parameter.name] = tmp_list

                        position = get_param_pos(parameter_pos, int(parameter.pos))
                        print parameter_pos, "ppos", parameter.pos, "pos",position
                        parameter_names.insert(position,parameter.name)
                        parameter_pos.insert(position,int(parameter.pos))

                    # add the values of tsi to the according parameter list

                    tmp_list = tsi_map[parameter.name]
                    tmp_list.append(tsi)
                    tsi_map[parameter.name] = tmp_list
                    #print crop_name, parameter.name, tsi

            row_no += 1
        
        # add a zero into the parameter list where no value was added
        # for this crop
        for pname  in parameter_names:
            tmp_list = tsi_map[pname]
            #print "Vergleiche", len(tmp_list), (num_crop+1)
            # add a zero because this parameter was not analysed for
            # the crop
            while (len(tmp_list)<(num_crop+1)):
                tmp_list.append(0.0)
                #print "Add a zero", len(tmp_list), num_crop+1
            tsi_map[pname] = tmp_list
            



    tsi_array = np.zeros([len(parameter_names),len(crop_names)])
    print tsi_array.shape
    print "Construct array"
    for crop_index, cname in enumerate(crop_names):

        print cname, crop_index        
        for p_index, pname in enumerate(parameter_names):
            
            plist = tsi_map[pname]
            tsi = plist[crop_index]
            tsi_array[p_index, crop_index] = tsi
            #print cname, pname, tsi

    x_pos = range(0,len(crop_names),1)
    y_pos = range(len(parameter_names)-1,-1,-1)

    delta = 0.5

    fig = plot.figure(figsize=(width,height)) #,frameon=False
    ax = fig.add_subplot(111)
    extent = -delta,len(crop_names)-delta, -delta, len(parameter_names)-delta
    im1 = plot.imshow(tsi_array, cmap=plot.cm.gist_heat_r, extent=extent, aspect='auto' ) # gist_yarg    #plot.cm.Paired
    im1.set_interpolation('nearest')

    ax.set_xlim(-delta, len(crop_names)-delta)
    ax.set_ylim(-delta, len(parameter_names)-delta)

    plot.xticks(x_pos, crop_names, fontsize=font_size, rotation=45, ha='right')
    plot.yticks(y_pos, parameter_names, fontsize=font_size)

    for x in x_pos:
        ax.axvline(x=x-delta, ls='solid', color='#333333', lw=1.5)
    for y in y_pos:
        ax.axhline(y=y-delta+0.03, ls='solid', color='#333333', lw=1.5)
    
    cbar = plot.colorbar()
    cbar.ax.tick_params(labelsize=font_size) 
    cbar.ax.set_title('Totaleffekt', fontsize=font_size)
    fig.savefig("compare_crop_tsi_vortrag.png", dpi_value=200)
        

def get_param_pos(liste,value):

    for pos,item in enumerate(liste):
        #print "Compare", item, value
        if (value<item):
            #print "Return ", pos+1
            return pos
        
    return len(liste)

main()




