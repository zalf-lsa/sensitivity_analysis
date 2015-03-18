#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys
import csv

output_path = "fast-diss/"
path_morris_files = "./"
morris_ranking_filename = "morris_ranking_final2.csv"

infos = ["parameter_definitions_winter_wheat.csv", "Winterweizen", 2]

output_map = {}
output_map[2] = "abBiom"
output_map[3] = "primYield"
output_map[4] = "rootBiomass"
output_map[5] = "shootBiomass"
output_map[6] = "leafBiomass"
output_map[7] = "lai"
output_map[8] = "GPP"
output_map[9] = "Ra"
output_map[10] = "abiomNContent"
output_map[11] = "cropHeight"
output_map[12] = "Eta"
output_map[13] = "soilMoist"
output_map[14] = "soilNmin"
output_map[15] = "soilTemp"
output_map[16] = "Corg"

def main():
    
  print infos
  crop_parameter_file = infos[0]
  row_index = infos[2]

  for key, value in output_map.iteritems():
    print key, value

    morris_ranking_file = open(morris_ranking_filename, "rb")
    ranking_csv = csv.reader(morris_ranking_file, delimiter="&")

    ranking_csv.next() # skip first row

    fast_parameter_file = open(output_path + value + ".csv", 'wb')
    fast_parameter_csv = csv.writer(fast_parameter_file, delimiter=';')
    fast_parameter_csv.writerow(["","","min","max","nom"])

    for ranking_row in ranking_csv:

        pno = ranking_row[0].strip()    
        rank = ranking_row[int(key)]
        print pno, rank, key, value
        if (rank==''):
          rank = 0
        else:
          rank = int(rank)
        if (pno!="" and rank>0):
          print ranking_row
          print pno, rank
          morris_parameter_file = open(path_morris_files + crop_parameter_file, 'rb')          
          morris_csv = csv.reader(morris_parameter_file, delimiter=';')
          morris_csv.next()

          for morris_p_row in morris_csv:
            if (pno == morris_p_row[0].strip()):
              print pno, morris_p_row[0]
              print
              fast_parameter_csv.writerow(morris_p_row)
              break
          morris_parameter_file.close()
    fast_parameter_file.close()
   
    morris_ranking_file.close()

       
  



########################################################
# start
########################################################




main()


