#!/usr/bin/python
# -*- coding: UTF-8-*-

# adding the path of monica's python module to search path
# that is defined by the environment variable PYTHONPATH

import sys


# add path to monica module to PATH
sys.path.append("..")
sys.path.append("../monica-src")


import math
import random
import datetime
import copy
import sa_functions 
import mpi_helper
import numpy
import os
import csv
import monica

from scipy.stats import rankdata

from monica_simulation import initializeMonicaSimulation
from saparameter import SAParameter
from apply_sa_values import applySAValues 
from string import lower
from mpi4py import MPI


###############################
size = 0
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
###############################

#############################################################################
# Configuration
#############################################################################
ranges =  20
start_vector_count = 40
random_start_vector_count = 500
schrittweite = 5

dx_fix = float(schrittweite) /  (float(ranges)-1)

add_analyse_after_calculation = 1

###########################################################################
# configuration for SA paper rev 1 (European Journal of Agronomy)
###########################################################################
#parameter_files_directory = "../configs/2015-07-sapaper_agronomy/parameters"
#simulation_files_path = "../configs/2015-07-sapaper_agronomy/"
#sites = ["Ascha", "Guelzow", "Werlte"]
#crops_to_analyse = [1,4,7,9,10,13] # [1,4,7,9,10,13]
###########################################################################


###########################################################################
# configuration for time-dependent SA paper
###########################################################################
parameter_files_directory = "../configs/2016-01-time-dependent-SA/parameters"
simulation_files_path = "../configs/2016-01-time-dependent-SA/sites/"
sites = ["Ascha" , "Dornburg", "Ettlingen", "Guelzow", "Werlte"]
crops_to_analyse=[1]
###########################################################################




#############################################################################
# End of Configuration
#############################################################################
    
""" 
####################################################################
Main routine for parallel MPI execution 
####################################################################
"""
def mpi_main():

  basis_output_dir = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
  #basis_output_dir = "2016-02-04"
  crop_map = getCropsForSA()

  for site in sites:
    
    output_dir = basis_output_dir + "/" + site + "/"

    for crop_id in crops_to_analyse:

      crop_info = crop_map[crop_id]
      
      # get information of analysed crop
      crop_name = crop_info.name
      crop_dir = crop_info.simulation_files_dir
      crop_parameter_file = parameter_files_directory + "/" + crop_dir + crop_info.parameter_file
     
      crop_simulation_files_dir = simulation_files_path + site + "/" + crop_info.simulation_files_dir
      print crop_id, crop_name, crop_parameter_file, crop_simulation_files_dir

      # create output director
     
      output_path = "runs/" + output_dir + "/crop-" + crop_name + "-min-step" + str(schrittweite) + "_range" + str(ranges) + "_startvector"+ str(start_vector_count)

      monica.activateDebugOutput(0)
	   
      node_list = []
      
      full_parameter_list = sa_functions.readParameterFile(crop_parameter_file)
      nominal_list = sa_functions.getNominalList(full_parameter_list)

      run_parameter_list= full_parameter_list
      run_parameter_count = len(run_parameter_list)
      
      # reads parameter specification from a file 
      parameter_count = len(full_parameter_list)      
      parameter_grid = generateParameterGrid(run_parameter_list,ranges)

      # HERMES configuration
      hermes_config = monica.getHermesConfigFromIni(crop_simulation_files_dir)
      env = monica.getHermesEnvFromConfiguration(hermes_config)
      env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)
      env = applySAValues(full_parameter_list, nominal_list, env, crop_id)    


      if (rank==0):
          if (not os.path.exists(output_path)):        
              os.makedirs(output_path)        

      trajectories = sa_functions.get_optimised_trajectories(full_parameter_list, parameter_grid, start_vector_count, ranges, schrittweite, random_start_vector_count)
          
      if (rank == 0):    
        trajectory_list = mpi_helper.splitListForNodes(trajectories, size)

        for n in range(size):
            #print n
            node_list.append([trajectory_list[n], parameter_grid, full_parameter_list])

  #        print node_list, len(node_list)

      ###################################################
      # parallel part
      ##################################################
      local_list = comm.scatter(node_list, root=0)        

	   
      local_trajectory_list = local_list[0]
      local_parameter_grid = local_list[1]
      local_parameter_list = local_list[2]
      print rank, "Received", len(local_trajectory_list), "elements" #, local_sv_list
      #print rank, "Received", local_parameter_grid
      #print rank, "Received", local_parameter_list

      effect_list = []
      if (len(local_trajectory_list) > 0):
        effect_list = runMorrisSA(local_parameter_list, local_parameter_grid, local_trajectory_list, env, output_path, crop_id)
      	
      global_effect_list = comm.gather(effect_list, root=0)
	
      ##################################################
      # end of parallel part
      ##################################################
       
	
      if (rank == 0):

          print ("Analyse daily parameter effects")
          output_names = sa_functions.getOutputNames()


          for output_index, output in enumerate(output_names):

              print "Analyse results for output:", output

              # create separate directory for each output
              output_directory = output_path + "/" + output
              if (not os.path.exists(output_directory)):
                  os.makedirs(output_directory)

              output_mu_array = []
              for parameter_index, parameter in enumerate(full_parameter_list):

                   print "Analyse parameter:", parameter.getName(), parameter_index              

                   file_handle = open(output_directory + "/" + parameter.getName() + ".csv", "w")
                   parameter_csv = csv.writer(file_handle, delimiter=";")
                   
                   daily_parameter_effects = []
                   for sample_index, sample_list in enumerate(global_effect_list):
              
		       daily_effect_list = sample_list[parameter_index][0][output_index]
                       parameter_csv.writerow(daily_effect_list)
                       daily_parameter_effects.append(daily_effect_list)

                   # calculate mu for each day
                   effect_array = numpy.array(daily_parameter_effects)
                   daily_mu_values = numpy.mean(effect_array, axis=0)
                   output_mu_array.append(daily_mu_values)

                   # close file of daily effect values
                   file_handle.close()

              # write scaled mu values to file
              mu_array = numpy.array(output_mu_array)

              # scale effects for one output   
              max_of_array = numpy.max(mu_array, axis=0)
              max_of_array[max_of_array == 0] = 1                                    
              mu_array /= (max_of_array*1.2) # cheating factor !!!
              
              rank_array = numpy.zeros(mu_array.shape)
              # calculate rank for each parameter at each simulation day
              for day_index, day in enumerate(max_of_array): 
                  rank_array[:, day_index] =  (len(mu_array[:, day_index])+1) - rankdata(mu_array[:,day_index], method='min')

              # save ranks to a csv file
              output_filehandle = open(output_path + "/" + output + "-mus.csv", "w")
              rank_filehandle = open(output_path + "/" + output + "-rank.csv", "w")
              output_csv = csv.writer(output_filehandle, delimiter=";")
              rank_csv = csv.writer(rank_filehandle, delimiter=";")

              for parameter_index, parameter in enumerate(full_parameter_list):
                   row = () 
                   rank_row = () 

                   row += (parameter.getName(), )
                   rank_row += (parameter.getName(), )

                   row += tuple(mu_array[parameter_index])
                   rank_row += tuple(rank_array[parameter_index])

                   output_csv.writerow(row)
                   rank_csv.writerow(rank_row)

              # close output csv file
              output_filehandle.close()
              rank_filehandle.close()             

          t_end = datetime.datetime.now()
          time_simulation = t_end - t_start 
          print "Node: ", rank, "\tSimulationszeit: ", time_simulation
   

##########################################################    
##########################################################    
##########################################################    


"""
generates the two dimensional parameter grid that is used by
morris sensitivity analysis
"""
def generateParameterGrid(parameter_list,range_number):
    
    parameter_number = len(parameter_list)       
    
    dx = []
    parameter_grid = []
    for param in range(parameter_number):
        
        max = float(parameter_list[param].getMax())
        min = float(parameter_list[param].getMin())
        
        dx.append(math.fabs(max-min)/(range_number-1.0))
        
        # generate row for parameter grid
        row = []
        # inner elements of the parameter grid
        for r in range(range_number):            
            row.append(min + dx[param]*(r))
            
        parameter_grid.append(row)

    return parameter_grid


##########################################################    
##########################################################    
##########################################################    

"""  
Method that performs a sensitivity analysis according to morris
with the given parameters. The parameters that should be analysed
are passed through the parameter list.
"""
def runMorrisSA(parameter_list, parameter_grid, local_trajectory_list, env, output_path, crop_id):
    
    t_morris_start = datetime.datetime.now()
     
    parameter_number = len(parameter_list)
    traj_size = len(local_trajectory_list)
    result_map = {}
     
     
    # initialisation of effect list and intput value list, needed for scaling the EE
    parameter_effects = []
    inputs = []    
    for i in range(len(parameter_list)):
        parameter_effects.append([])
        inputs.append([])

    outputs=[]
    for i in range(len(monica.sensitivityAnalysisResultIds())):
        outputs.append([])
         
    #print rank, parameter_grid
    for traj in local_trajectory_list:        

      result_old = None
      point_old = None
      for point_number, point in enumerate(traj):
        
        # first model evaluation
        if (point_number==0):
          start_vector_index = point
          print rank, "\t",  point_number, "/", len(traj) #, "\t", point 
          result_old =getResult(result_map, start_vector_index, env, parameter_list, parameter_grid, crop_id)            
        else:
          print rank, "\t",  point_number, "/", len(traj) #, "\t", point  
          # next model evaluation
          result_new = getResult(result_map, point, env, parameter_list, parameter_grid, crop_id)
                         
          parameter_index = 0
          for o,n in zip(point_old, point):
            if (o == n):
              parameter_index += 1
            else:
              break
          
          #dx_new = (parameter_list[parameter_index].getMax() - parameter_list[parameter_index].getMin())/dx_fix
          #dx = math.fabs(parameter_grid[parameter_index][0] - parameter_grid[parameter_index][1])
          dx = dx_fix
          parameter_value = parameter_grid[parameter_index][point[parameter_index]]
          parameter_name = parameter_list[parameter_index].getName()
          inputs[parameter_index].append(parameter_value)
          
          effects = analyseResults(result_old, result_new, parameter_index, dx, outputs, parameter_value, parameter_name)
          #print parameter_name, parameter_index, effects
          parameter_effects[parameter_index].append(effects)                    
          result_old = result_new    
        point_old = point
     
    t_morris_end = datetime.datetime.now()
    time_morris_step = t_morris_end - t_morris_start
    print rank, "Morris_Endtime: ", time_morris_step
    
    return parameter_effects
        
##########################################################    
##########################################################    
##########################################################    
    
def getResult(result_map, start_vector_index, env, parameter_list, parameter_grid, crop_id):
    #print start_vector_index
    if (tuple(start_vector_index) in result_map):
        result = result_map[tuple(start_vector_index)]
    else:        
        # this point has not been calculated by model yet             
        start_vector = [parameter_grid[p][start_vector_index[p]] for p in range(len(parameter_list)) ]
        #print start_vector
        new_env = applySAValues(parameter_list, start_vector, env, crop_id)
        monica.activateDebugOutput(0)
        result = monica.runMonica(new_env)
        result_map[tuple(start_vector_index)] = copy.copy(result)
    return result   
        
    
##########################################################    
##########################################################    
##########################################################    


""" 
Analyses result of different model runs and calculates
elementary effect for a parameter specified by parameter index
"""
def analyseResults(result_old, result_new, parameter_index, dx, outputs, parameter_value, name=None):
    
    #dx = 0.0333
    result_ids = monica.sensitivityAnalysisResultIds()    
    effects_per_day = []
    
    output_names = sa_functions.getOutputNames()
    for output_index, output_id in enumerate(result_ids):
        
        #print monica.resultIdInfo(id).shortName
        old = result_old.getResultsById(output_id)
        new = result_new.getResultsById(output_id)
      
        simulation_day = 0
        
        effects = []
        daily_outputs = []
        for (old_value, new_value) in zip (old, new):
            
            if (sa_functions.is_nan(old_value) or 
                sa_functions.is_inf(old_value) or 
                sa_functions.is_nan(new_value) or 
                sa_functions.is_inf(new_value)):            
                print ("INF or NAN at day: ", simulation_day, "\t",  monica.resultIdInfo(output_id).shortName,"\tOLD: ", old_value, "\tNEW: ", new_value, "DX: ", dx, "param: ", parameter_index, parameter_value)
              #quit()
            else:
                result_per_day = math.fabs(old_value-new_value) / dx                          
                effects.append(result_per_day)                        
        
            if (not(sa_functions.is_nan(new_value)) or not(sa_functions.is_inf(new_value))):
                daily_outputs.append(new_value)

            simulation_day +=1
            
        outputs[output_index].append(daily_outputs)
	    
        effects_per_day.append(effects)
    return effects_per_day


    
            
        
##########################################################
##########################################################
##########################################################

class CropInfo:
  def __init__(self, crop_id, german_name, name, parameter_file, simulation_files_dir=None):
    self.crop_id = crop_id
    self.german_name = german_name
    self.name = name
    self.parameter_file = parameter_file
    self.simulation_files_dir = simulation_files_dir

        
#####################################################
#####################################################
#####################################################

def getCropsForSA():

  crop_map = {}
  crop_map[1]   = CropInfo(1, "Winterweizen", "Winter wheat",  "parameter_definitions_winter_wheat.csv", "winter-wheat/")
  crop_map[19]  = CropInfo(19,"Wintertriticale", "Winter triticale","parameter_definitions_winter_triticale.csv", "winter-triticale/")
  crop_map[2]   = CropInfo(2, "Wintergerste", "Winter barley", "parameter_definitions_winter_barley.csv", "winter-barley/")
  crop_map[9]   = CropInfo(9, "Winterraps", "Winter rapeseed", "parameter_definitions_winter_rape.csv", "winter-rape/")
  crop_map[4]   = CropInfo(4, "Sommergerste",  "Spring barley", "parameter_definitions_spring_barley.csv", "spring-barley/")
  crop_map[23]  = CropInfo(23,"Sommertriticale", "Spring triticale","parameter_definitions_spring_triticale.csv", "spring-triticale/")
  crop_map[10]  = CropInfo(10,"Zuckerrübe", "Sugar beet", "parameter_definitions_sugarbeet.csv",  "sugar-beet/")
  crop_map[7]   = CropInfo(7, "Silomais", "Silage Maize", "parameter_definitions_maize.csv", "maize/")
  crop_map[18]  = CropInfo(18,"Sudangras", "Sorghum", "parameter_definitions_sudangras.csv", "sudangras/")
  crop_map[12]  = CropInfo(12,"Phacelia", "Phacelia", "parameter_definitions_phacelia.csv", "phacelia/")
  crop_map[13]  = CropInfo(13,"Kleegras", "Clover grass ley", "parameter_definitions_clover.csv", "clover/")
  crop_map[14]  = CropInfo(14,"Luzernegras", "Alfalfa", "parameter_definitions_alfalfa.csv", "alfalfa/")
  crop_map[16]  = CropInfo(16,"Weidelgras", "Rye grass", "parameter_definitions_rye_grass.csv", "rye_grass/")
  crop_map[22]  = CropInfo(22,"Hafer", "Oat", "parameter_definitions_oat.csv", "oat/")
  return crop_map


#####################################################
#####################################################
#####################################################

""" 
Hauptprogramm
"""
t_start = datetime.datetime.now()
mpi_main()
