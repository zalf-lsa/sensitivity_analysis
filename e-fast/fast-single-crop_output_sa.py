#!/usr/bin/python
# -*- coding: UTF-8-*-

from mpi4py import MPI
import numpy
import math
import matplotlib.pyplot as plot
import sys
import datetime
import csv
import os

# add path to monica module to PATH

root_path = ".."

sys.path.append("..")
sys.path.append("../monica-src")



import monica
from monica_simulation import initializeMonicaSimulation
from monica_simulation import initializeEVA2MonicaSimulation
from saparameter import SAParameter
from apply_sa_values import applySAValues
import sa_functions 
import fast_lib 
import mpi_helper
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

#crops = [1, 2, 4, 7, 9, 10, 13, 14,18,19,12, 22,23,16]
#crops = [12, 19, 22,23,16]
crops = [1,4,7,9,10,13]

#output_list = ["abBiom","rootBiomass","shootBiomass","leafBiomass","lai","GPP","Ra","abiomNContent","cropHeight","Eta","soilMoist","soilNmin","soilTemp","Corg"]
#output_list = ["primYield", "abiomNContent","abBiom","rootBiomass","shootBiomass","leafBiomass","lai","GPP","Ra","Eta","soilMoist","soilNmin","soilTemp","Corg"]
#output_list = ["soilNmin","soilTemp","Corg"]
output_list = ["aboveGroundBiomass"]

sites = ["Ascha","Guelzow","Werlte"]


max_omega = 4096
sample_size = 20000

parameters_path = "../configs/2015-07-sapaper_agronomy/parameters/fast_parameters"

""" 
####################################################################
Main routine for parallel MPI execution 
####################################################################
"""
def mpi_main(crop):
    
    
    for output_index, output in enumerate(output_list):
      for site in sites:

          crop_map = sa_functions.getCropsForSA()
          crop_info = crop_map[crop]
          
          directory = datetime.datetime.today().strftime("runs/2015-07-13/crop-" + str(crop) + "/" + site)
          output_path = "runs/"
       
          # HERMES configuration
          simulation_path = getHermesSimulationPath(crop, site)
          hermes_config = monica.getHermesConfigFromIni(simulation_path)
          env = monica.getHermesEnvFromConfiguration(hermes_config)
          env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)
          
          if (rank == 0):
            if (not os.path.exists(directory)):  
              os.makedirs(directory)
            if (not os.path.exists(directory + "/parameters")):
              os.makedirs(directory + "/parameters") 

          output_path = directory  
          filename = parameters_path + "/" + crop_info.parameter_file 

          if (rank == 0):
            if (not os.path.exists(output_path)):  
              os.makedirs(output_path) 
          
          # reads parameter specification from a file
          complete_list = sa_functions.readParameterFile(filename)
          nominal_list = sa_functions.getNominalList(complete_list)

          max_param = len(complete_list)
          result_ids = monica.sensitivityAnalysisResultIds() 
          env = applySAValues(complete_list, nominal_list, env, crop)     

          if (rank==0):
              tsi_filename = output_path + "/tsi-" + output + ".csv"
              tsi_filehandle = open(tsi_filename, "w")
              tsi_file = csv.writer(tsi_filehandle)
              tsi_file.writerow(["Parameter", "FI", "TSI"])

          # individual analysis for each parameter as proposed in Saltelli (2000)
          for parameter_index, parameter in enumerate(complete_list):

              print 2, parameter_index, parameter
              all_node_list = []
              if (rank == 0):
                  (all_node_list, sample_list) = initialise_parameter_value_list(complete_list, parameter_index)

              t_morris_start = datetime.datetime.now()

              print rank, "Start Parallel", datetime.datetime.now()
              ###################################################
              # parallel part
              ##################################################
              sample_list = comm.scatter(all_node_list, root=0)      

              print rank, "Received local list ", parameter.getName() #, sample_list , 

              #print rank, monica.resultIdInfo(output_id[output_index]).shortName

              local_result_map = {}
              for sample_index, sample in enumerate(sample_list):
                  print rank, sample_index, "/", len(sample_list), "\t",sample
                  new_env = applySAValues(complete_list, sample, env, crop)        
                  monica.activateDebugOutput(0)
                  result = monica.runMonica(new_env)       
                  #print result.getResultsById(output_id[output_index])       
                  
                  value = sa_functions.getMaxOfList(result.getResultsById(output_index)) # original getMeanOfList
                  print value
                  local_result_map[str(sample)] = value

              ###################################################
              # end of parallel part
              ##################################################
              result_list = comm.gather(local_result_map, root=0)   

              t_parallel_end = datetime.datetime.now()
              print rank, "end_parallel part", t_parallel_end, "parallel_part1:", t_parallel_end - t_morris_start

              if (rank == 0):


                  output_filename = output_path + "/parameters/" + parameter.getName() + ".csv"
                  csv_filehandle = open(output_filename, "wb")
                  csv_file = csv.writer(csv_filehandle)
                  names = []
                  for n in range(max_param):
                      names.append(complete_list[n].getName())
                  names.append(output)
                  csv_file.writerow(names)

                  # iterate sequentially through sample list
                  v_input = []
                  for n in range(max_param):
                      v_input.append([])
                  v_output = []    
                  for node_index, node in enumerate(all_node_list):
                      for item in node:
                          value = result_list[node_index][str(item)]     
                          
                          item.append(value)
                          for n in range(max_param):
                              v_input[n].append(item[n])
                          if (value != None):
                              #print item, value
                              v_output.append(value)
                              csv_file.writerow(item)
                                             
                  csv_filehandle.close()

                  omegas = fast_lib.get_ts_frequencies(max_omega, parameter_index, max_param)
                  print omegas

                  max_omega_o = int(max_omega/8)

                  # fast fourier transformation
                  fft = numpy.fft.fft(v_output)
                  real=numpy.real(fft)
                  imag= numpy.imag(fft)
                  A=numpy.add(numpy.multiply(real,real),numpy.multiply(imag,imag))
                  D_all = sum(A[1:])

                  first_order = 0
                  tsi = 0
                  analysed_omegas = []
                  print "Analysing TSI for", names[parameter_index]
                  print "D_all", D_all
                  for index in range(max_param):
                     
                      harmonics = []
                      for i in range(1,5):
                          harmonics.append(omegas[index]*i)
                      
                      D_i = 2*sum(A[harmonics])
                      if (parameter_index == index):
                          first_order = (D_i/D_all)
                      else:
                          if (omegas[index] not in analysed_omegas):
                              tsi +=D_i
                              print "Add first order of", names[index], "(", D_i,")",tsi,  omegas[index]
                      analysed_omegas.append(omegas[index])

                  tmp_tsi = 2*sum(A[1:4*max_omega_o])
                  tsi = 1-(tmp_tsi/D_all)

                  if (parameter_index == 0):
                      tsi_file.writerow([output, D_all, None])

                  name = complete_list[parameter_index].getName()
                  tsi_file.writerow([name, first_order, tsi,tmp_tsi])
                  print "--------------------"
                  print "Seq time2", datetime.datetime.now() - t_parallel_end 


          if (rank==0):            
              tsi_filehandle.close()                   

                        
        
def getHermesSimulationPath(crop_id, site):
  

    path = "../configs/2015-07-sapaper_agronomy/" + site + "/"

    crop_map = sa_functions.getCropsForSA()
    crop_info = crop_map[crop_id]
    path += crop_info.simulation_files_dir
    
    return path


"""
Initialisation routine for root node
"""
def initialise_parameter_value_list(complete_list, parameter_index):

    # initialisation of integer frequencies for TS analysis
    omegas = []
    max_param = len(complete_list)
    for i in range(max_param): 
        omegas.append(fast_lib.get_ts_frequencies(max_omega, i, max_param))

    print omegas
    
    # initialisation of input signal range from -pi to pi
    step = (2*math.pi) / sample_size
    input_signal = numpy.arange(-math.pi, math.pi, step)  

    sample_list = []    
    for s in input_signal:                            
        sample = []
        for i in range(max_param):            
            sample.append(fast_lib.get_transformation(omegas[parameter_index][i], s, complete_list[i]))
        sample_list.append(sample)

    all_node_list = mpi_helper.splitListForNodes(sample_list, size)
    return (all_node_list,input_signal)


t_start = datetime.datetime.now()
for crop in crops:
  print rank, "New crop", crop
  mpi_main(crop)

t_end = datetime.datetime.now()
time_simulation = t_end - t_start 
print "Node: ", rank, "\tSimulationszeit: ", time_simulation

