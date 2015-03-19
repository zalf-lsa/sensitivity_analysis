#! /usr/bin/python
import os
import csv
import re
import numpy


max_rank = 30

"""
Main method
"""
def main():
    
    outpath = "runs/crop-Winterweizen-min-step1_range5_startvector10"
    dir_file = open("output_dir.r", "r")
    for line in dir_file:
        mo = re.match("directory=\"(.*)\"", line)
        if (mo != None):        
            path = str(mo.group(1))
    
    print "Ranking files in", path
    
    mean_filename = path + "/morris_mean.txt"
    
    csv_data = csv.reader(open(mean_filename, "rb"),delimiter='\t')

    mean_data=[]
    index = 0
    for d in csv_data:        
        if (index==0):
          mean_data.append(list(d))          
          print "d",d
        else:           
          print d
          tmp_list = [d[0]]
          for item_index in range(1,len(d)): 
            tmp_list.append(float(d[item_index]))
          mean_data.append(tmp_list)    
          

        index += 1

    if (not os.path.isdir(path + "/ranking")):
        os.makedirs(path + "/ranking")     

    array = []
    array.append(["Parameter","Rang","Wert"])
    number_rows = len(mean_data)

    row_number = len(mean_data)
    column_number = len(mean_data[0])
    print "Col: ", column_number, "    row: ", row_number
    for index in range(1,column_number):
       	print "SORT: ", index, column_number
        for i,r in enumerate(mean_data):        
            i,r
            if (len(r)<column_number):
                mean_data[i].append(0.0)
        sorted_array = sorted(mean_data, key=lambda columns: columns[index], reverse=True) # sort whole array by column index
	
        swap_data = zip(*sorted_array)

        array = []
        parameter_names = swap_data[0]
        sorted_row = swap_data[index]

        
        ranking_file = open(path+"/ranking/"+mean_data[0][index] + ".txt", "w")
        ranking_csv = csv.writer(ranking_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)        
        rank = 0
        print sorted_row
        parameter_row = []
        rank_row = []
        value_row = []
        for item_index, item in enumerate(sorted_row):
            
            if (is_number(item) and (float(item)>=0.001) and (rank<max_rank) and item!="nan"): # 
                parameter_row.append(parameter_names[item_index])
                rank_row.append(rank + 1)
                value_row.append(item)
                rank += 1
        if (len(value_row)>0 and value_row[0]!=0.0):
          value_row = numpy.true_divide(value_row, value_row[0])
        array.append(parameter_row)
        array.append(rank_row)
        array.append(value_row)              
      
        ranking_csv.writerows(array)    
        ranking_file.close()



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False



main()
    
