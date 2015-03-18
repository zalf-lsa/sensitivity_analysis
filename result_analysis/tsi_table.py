#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys
import csv
import datetime
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.lines import Line2D
sys.path.append('../..')
import analysing_lib

xlim=0.00
ylim=0.00

width=14
height=20

path = "2013-02-23/"

output_list = ["abBiom","primYield","rootBiomass","shootBiomass","leafBiomass","lai","GPP","Ra","abiomNContent","cropHeight","Eta","soilMoist","soilNmin","soilTemp","Corg"]
crop_names = ['AGB', 'Yield', 'Root', 'Shoot', 'Leaf', 'LAI', 'GPP', r'R$_a$', 'NConc', 'Height', r'ET$_a$', 'Moist', r'N$_{min}$', 'Temp', r'C$_{org}$']

font_size = '16'

def main():

    tsi_map = {}

    parameter_names = []
    parameter_pos = []
    parameter_no = []
    for num_crop, output in enumerate(output_list):


        directory = output
        print "Analysing output", num_crop, output

        parameter_map = analysing_lib.get_parameter_map()

        result_file = open(path + directory + "/tsi-" + output + ".csv", "rb")
        result_csv = csv.reader(result_file, delimiter=",")

        row_no = 0
        for result in result_csv:
            if (row_no>1):

                parameter = parameter_map[str(result[0]) + ".txt"]
                tsi = float(result[2])

                if (parameter.name not in tsi_map):
                    tmp_list= []
                    while (len(tmp_list)<(num_crop)):
                        tmp_list.append(0.0)
                    tsi_map[parameter.name] = tmp_list

                    position = get_param_pos(parameter_pos, int(parameter.pos))
                    print parameter_pos, "ppos", parameter.pos, "pos",position, "no", parameter.no
                    parameter_names.insert(position,parameter.name)
                    parameter_pos.insert(position,int(parameter.pos))
                    parameter_no.insert(position,str(parameter.no))

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
            




    output_file = open("tsi_latex_table.txt", 'wb')
    csv_output = csv.writer(output_file, delimiter='&')

    for p_index, pname in enumerate(parameter_names):
    
        csv_row = [parameter_no[p_index], pname ]

        for crop_index, cname in enumerate(crop_names):

            #print cname, crop_index        
            
            plist = tsi_map[pname]
            tsi = float(plist[crop_index])
            tsi = str(round(tsi,2))
            if (tsi=="0.0"):
                tsi=""
            csv_row.append(tsi)
            #print cname, pname, tsi
        
        print csv_row
        csv_output.writerow(csv_row)

    output_file.close()

        

def get_param_pos(liste,value):

    for pos,item in enumerate(liste):
        #print "Compare", item, value
        if (value<item):
            #print "Return ", pos+1
            return pos
        
    return len(liste)

main()




