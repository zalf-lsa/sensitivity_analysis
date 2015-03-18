#!/usr/bin/python
# -*- coding: UTF-8-*-

import sys
import csv

output_path = "dissertation/optimization/"
path_morris_files = "dissertation/"
morris_ranking_filename = "optimization_helper.csv"

folders = [["parameter_definitions_winter_wheat.csv", "Winterweizen", 2],
            ["parameter_definitions_winter_triticale.csv", "Wintertriticale", 3],
            ["parameter_definitions_winter_barley.csv", "Winter barley", 4],
            ["parameter_definitions_winter_rape.csv", "Winter rape", 5],
            ["parameter_definitions_spring_barley.csv", "Spring barley", 6],
            ["parameter_definitions_spring_triticale.csv", "Spring triticale", 7],
            ["parameter_definitions_sugarbeet.csv", "Sugar beet", 8],
            ["parameter_definitions_maize.csv", "Maize", 9],
            ["parameter_definitions_sudangras.csv", "Sudangras", 10],
            ["parameter_definitions_phacelia.csv", "Phacelia", 11],
            ["parameter_definitions_clover.csv", "Clover", 12],
            ["parameter_definitions_alfalfa.csv", "Alfalfa", 13],
            ["parameter_definitions_rye_grass.csv", "Rye grass", 14],            
            ["parameter_definitions_oat.csv", "Oat", 15]
          ]


def main():
    
  for infos in folders:
    print infos
    crop_parameter_file = infos[0]
    row_index = infos[2]

    morris_ranking_file = open(morris_ranking_filename, "rb")
    ranking_csv = csv.reader(morris_ranking_file, delimiter="&")

    ranking_csv.next() # skip first row

    fast_parameter_file = open(output_path + crop_parameter_file, 'wb')
    fast_parameter_csv = csv.writer(fast_parameter_file, delimiter=';')
    fast_parameter_csv.writerow(["","","min","max","nom"])

    for ranking_row in ranking_csv:
        print ranking_row
        pno = ranking_row[0].strip()    
        rank = int(ranking_row[row_index])

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


