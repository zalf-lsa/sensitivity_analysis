#! /usr/bin/python


import os
import csv
import re
import numpy
import sys
import math

"""
Main method
"""
def main():
    
    outpath = ""
    dir_file = open("output_dir.r", "r")
    for line in dir_file:
        mo = re.match("directory=\"(.*)\"", line)
        if (mo != None):        
            outpath = str(mo.group(1))
    
    path = outpath + "/outputs/"
    print "Analysing files in", path
    
    
    files = getFilesInDirectory(path)
    
    mean_filename = "morris_mean.txt"
    mean_file = open(outpath + "/" + mean_filename, "w")
    mean_csv_writer = csv.writer(mean_file,delimiter='\t')

    std_filename = "morris_std.txt"
    std_file = open(outpath + "/" + std_filename, "w")
    std_csv_writer = csv.writer(std_file,delimiter='\t')
    
    output_names = ["Parametername"]
    mean_data = []
    std_data = []
    for file_index, file in enumerate(files):
        print "Analysing file \"", file, "\""
        data = csv.reader(open(path + "/" + file, "rb"),delimiter='\t')

        # swap rows and columns from input file
        new_data = []
        for d in data:
            new_data.append(d)
        data = zip(*new_data)

        

        mean_row = [file]
        std_row = [file]
        for index, row in enumerate(data):
            if (file_index==0):
                # initialise header only when reading first file         
                output_names.append(row[0])                
            print row
            data_row = row[1:len(row)]
            new_mean = []
            new_std = []
            for x in data_row:
                if(x!='' and x!='nan ' and x!='nan'):
                    print x, type(x)
                    new_mean.append(math.fabs(float(x)))
                    new_std.append(float(x))

            print getMean(new_mean), getSTD(new_std)
            mean_row.append(getMean(new_mean))
            std_row.append(getSTD(new_std))
        mean_data.append(mean_row)
        std_data.append(std_row)

    mean_csv_writer.writerow(output_names)
    std_csv_writer.writerow(output_names)

    mean_csv_writer.writerows(mean_data)
    std_csv_writer.writerows(std_data)
        
    mean_file.close()
    std_file.close()

""" 
Analyses file and returns a list with
mean values
"""    
def getMean(list):
    new_list = []
    for item in list:
        if(item=='' or item=='nan ' or item=='nan'):
            new_list.append(0.0)
        else:
            new_list.append(float(item))

    mean = numpy.mean(new_list)
    return (round(mean,2))
   
   
    
""" 
Analyses file and returns a list with
standard mean deviation values
"""     
def getSTD(list):
    new_list = []
    for item in list: 
        if(item=='' or item=="nan "):
            new_list.append(0.0)
        else:
            new_list.append(float(item))

    std = numpy.std(new_list)
    return (round(std,2))
    


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

main()
