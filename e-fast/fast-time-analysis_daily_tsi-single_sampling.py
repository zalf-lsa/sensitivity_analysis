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

crops = [1]


output_list = [ ["primaryYield", "parameter_definitions_winter_wheat-primyield.csv"]
                ]                #,"dailyAGB", "dailyAGB_N","ETa","soilMoist0_90cm","nmin0_90cm"]

#sites = ["Ascha","Dornburg","Ettlingen""Guelzow","Werlte"]
sites = ["Ascha"]

max_omega = 2048
sample_size = 10000

parameters_path = "../configs/2015-03-time-dependent-SA/parameters/winter-wheat/tdsa_parameters"

""" 
####################################################################
Main routine for parallel MPI execution 
####################################################################
"""
def mpi_main(crop):
    

    for output_index, output_info in enumerate(output_list):
        output = output_info[0]
        parameter_filename = output_info[1]
        for site in sites:

          crop_map = sa_functions.getCropsForSA()
          crop_info = crop_map[crop]
          
          directory = datetime.datetime.today().strftime("runs/time/2015-09-17/" + str(output) + "/" + site)
          output_path = "runs/"
       
          # HERMES configuration
          simulation_path = getHermesSimulationPath(crop,site)
          hermes_config = monica.getHermesConfigFromIni(simulation_path)
          env = monica.getHermesEnvFromConfiguration(hermes_config)
          env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)

          if (rank == 0):
            if (not os.path.exists(directory)):  
              os.makedirs(directory) 

          output_path = directory + "/" + crop_info.simulation_files_dir
          filename = parameters_path + "/" + parameter_filename

          if (rank == 0):
            if (not os.path.exists(directory)):  
              os.makedirs(directory) 
          
          # reads parameter specification from a file
          complete_list = sa_functions.readParameterFile(filename) # original mc_parameter.txt    
          nominal_list = sa_functions.getNominalList(complete_list)

          max_param = len(complete_list)
          result_ids = monica.sensitivityAnalysisResultIds() 
          env = applySAValues(complete_list, nominal_list, env, crop)     

          # individual analysis for each parameter as proposed in Saltelli (2000)
          for parameter_index, parameter in enumerate(complete_list):

              print 2, parameter_index, parameter

              dates = ["Date"]
              dev_stages = ["Stage"]

              all_node_list = []
              if (rank == 0):
                  (all_node_list, sample_list) = initialise_parameter_value_list(complete_list, parameter_index)

              ###################################################
              # parallel part
              ##################################################
              sample_list = comm.scatter(all_node_list, root=0)      

              print rank, "Received local list ", parameter.getName() #, sample_list , 

              #print rank, monica.resultIdInfo(output_id[output_index]).shortName

              local_result_map = {}
              for sample_index, sample in enumerate(sample_list):
                  #print rank, sample_index, "/", len(sample_list), "\t",sample
                  new_env = applySAValues(complete_list, sample, env, crop)        
                  monica.activateDebugOutput(0)
                  result = monica.runMonica(new_env)            
                  values = result.getResultsById(output_id[output_index])
                  local_result_map[str(sample)] = values

                  if (sample_index==0):
                    dev_stage_id = (sa_functions.getOutputId(["devStage"]))[0]
                    dates.extend(list(result.dates))
                    dev_stages.extend(list(result.getResultsById(dev_stage_id)))
                    
              ###################################################
              # end of parallel part
              ##################################################
              result_list = comm.gather(local_result_map, root=0)   
              print "end_parallel part"

              if (rank == 0):

                  result_array = []
                  for node_index, node in enumerate(all_node_list):
                      for item in node:
                          values = result_list[node_index][str(item)]     
                          if (values != None):
                              #print item, value
                              result_array.append(values)
                              

                  sample_count = len(result_array)
                  days = len(result_array[0])

                  minimum_values = ["Min"]
                  maximum_values = ["Max"]
                  mean_values = ["Mean"]
                  tsi_values = ["TSi"]
                  si_values = ["Si"]

                  for day in range(0,days,1):
                    values_for_day = []
                    for sample in range(0,sample_count):
                      values_for_day.append(result_array[sample][day])

                    tsi = 0.0
                    first_order= 0.0
                    if ( min(values_for_day) == 0 and max(values_for_day) == 0):
                      pass
                    else:
              
                      
                      omegas = fast_lib.get_ts_frequencies(max_omega, parameter_index, max_param)
        
                      max_omega_o = int(max_omega/8)

                      fft = numpy.fft.fft(values_for_day)
                      real=numpy.real(fft)
                      imag= numpy.imag(fft)
                      A=numpy.add(numpy.multiply(real,real),numpy.multiply(imag,imag))
                      
                      D_all = sum(A[1:])
                      

                      first_order = 0
                      tsi = 0
                      analysed_omegas = []
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

                          analysed_omegas.append(omegas[index])

                      tmp_tsi = 2*sum(A[1:4*max_omega_o])
                      tsi = 1-(tmp_tsi/D_all)

                    minimum_values.append(numpy.min(values_for_day))
                    maximum_values.append(numpy.max(values_for_day))
                    mean_values.append(numpy.mean(values_for_day))
                    tsi_values.append(tsi)
                    si_values.append(first_order)  

                  time_filename = directory + "/" + parameter.getName() + ".csv"
                  time_filehandle = open(time_filename, "wb")
                  time_file = csv.writer(time_filehandle)
                  time_file.writerow(dates)
                  time_file.writerow(dev_stages)
                  time_file.writerow(minimum_values)
                  time_file.writerow(maximum_values)
                  time_file.writerow(mean_values)
                  time_file.writerow(tsi_values)
                  time_file.writerow(si_values)
                
                  time_filehandle.close()
                
        
def getHermesSimulationPath(crop_id, site):
  

    path = "../configs/2015-03-time-dependent-SA/sites/" + site + "/"

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
            if (i == parameter_index):
              sample.append(fast_lib.get_transformation(omegas[parameter_index][i], s, complete_list[i]))
            else:
              #print complete_list[i].getName(), complete_list[i].getNominal()
              sample.append(complete_list[i].getNominal())

        sample_list.append(sample)

    all_node_list = mpi_helper.splitListForNodes(sample_list, size)
    return (all_node_list,input_signal)


for crop in crops:
  print rank, "New crop", crop
  mpi_main(crop)

