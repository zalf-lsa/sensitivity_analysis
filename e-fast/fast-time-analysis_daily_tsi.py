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


###########################################################
# adding the path of monica's python module to search path
# that is defined by the environment variable PYTHONPATH
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

##############################################################    
##############################################################
##############################################################

crops = [1]

output_list = [ #[ "dailyPrimYield",    "parameter_definitions_winter_wheat-primyield.csv"],
                [ "dailyAGB",        "parameter_definitions_winter_wheat-agb.csv"]
                #[ "dailyAGB_N",      "parameter_definitions_winter_wheat-nagb.csv"],
                #[ "ETa",             "parameter_definitions_winter_wheat-eta.csv"],
                #[ "soilMoist0_90cm", "parameter_definitions_winter_wheat-moist.csv"],
                #[ "nmin0_90cm",      "parameter_definitions_winter_wheat-nmin.csv"]
                ]


output_names = [o[0] for o in output_list]
output_id = sa_functions.getOutputId(output_names)
print "Output_id:", output_id

max_omega = 4096    
sample_size = 32769

sites = ["Ascha", "Dornburg", "Ettlingen", "Guelzow","Werlte"]

parameters_path = "../configs/2016-01-time-dependent-SA/parameters/winter-wheat/tdsa_parameters"


##############################################################    
##############################################################
##############################################################

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
          
            directory = datetime.datetime.today().strftime("runs/time/2016-02-16a/" + str(output) + "/" + site)
            output_path = "runs/"
       
            # HERMES configuration
            simulation_path = getHermesSimulationPath(crop,site)
            hermes_config = monica.getHermesConfigFromIni(simulation_path)
            env = monica.getHermesEnvFromConfiguration(hermes_config)
            env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)

            if (rank == 0):
                if (not os.path.exists(directory)):  
                    os.makedirs(directory) 

            output_path = root_path + "/python/sensitivity_analysis/fast/output/" + directory
            filename = parameters_path + "/" + parameter_filename

                   
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
                    #print output_index, output_id            
                    values = result.getResultsById(output_id[output_index])
                    local_result_map[str(sample)] = values

                    if (sample_index==0):
                        dev_stage_id = (sa_functions.getOutputId(["dev_stage"]))[0]
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
                    std_values = ["STD"]
                    #mean_si_values = ["Mean-Si"]
                    #mean_tsi_values = ["Mean-TSi"]

                    for day in range(0,days,1):
                        values_for_day = []
                        for sample in range(0,sample_count):
                            values_for_day.append(result_array[sample][day])

                        #print
                        #print "Tag:", day
                        #print len(values_for_day)
                        #print len(mean_values_for_day)
    
                        tsi = 0.0
                        first_order= 0.0
                        if ( min(values_for_day) == 0 and max(values_for_day) == 0):
                            pass
                        else:          

                            ######################################
                            tsi, first_order = get_TSi_Si(max_omega, parameter_index, max_param, values_for_day)
                            ######################################
                            
                        #^print tsi, first_order

                        minimum_values.append(numpy.min(values_for_day))
                        maximum_values.append(numpy.max(values_for_day))
                        mean_values.append(numpy.mean(values_for_day))

                        std_values.append(numpy.std(values_for_day))
                        tsi_values.append(tsi)
                        si_values.append(first_order)  
                        #mean_tsi_values.append(tsi_mean)
                        #mean_si_values.append(si_mean)

                    time_filehandle = open(time_filename, "w")
                    time_file = csv.writer(time_filehandle)
                    time_file.writerow(dates)
                    time_file.writerow(dev_stages)
                    time_file.writerow(minimum_values)
                    time_file.writerow(maximum_values)
                    time_file.writerow(mean_values)
                    time_file.writerow(std_values)
                    time_file.writerow(tsi_values)
                    time_file.writerow(si_values)
                    #time_file.writerow(mean_tsi_values)
                    #time_file.writerow(mean_si_values)            
                    time_filehandle.close()
                
##############################################################    
##############################################################
##############################################################
        
def get_TSi_Si(max_omega, parameter_index, max_param, values):
    omegas = fast_lib.get_ts_frequencies(max_omega, parameter_index, max_param)    
    max_omega_o = int(max_omega/8)
    fft = numpy.fft.fft(values)
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

    return tsi, first_order
    
##############################################################    
##############################################################
##############################################################

def getHermesSimulationPath(crop_id, site):  

    path = "../configs/2016-01-time-dependent-SA/sites/" + site + "/"

    crop_map = sa_functions.getCropsForSA()
    crop_info = crop_map[crop_id]
    path += crop_info.simulation_files_dir
    
    return path

##############################################################    
##############################################################
##############################################################

"""
Initialisation routine for root node
"""
def initialise_parameter_value_list(complete_list, parameter_index):

    # initialisation of integer frequencies for TS analysis
    omegas = []
    max_param = len(complete_list)
    for i in range(max_param): 
        omegas.append(fast_lib.get_ts_frequencies(max_omega, i, max_param))

    #print omegas
    
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


##############################################################    
##############################################################
##############################################################

# Call of main routine
for crop in crops:
  print rank, "New crop", crop
  mpi_main(crop)

