#!/usr/bin/python
# -*- coding: UTF-8-*-

from saparameter import SAParameter
import monica
import sys
import numpy
import csv
import math
import random
from mpi4py import MPI
import mpi_helper

# constants definition
ADDITION = 0
SUBTRACTION = 1

Infinity = 1e10000
NaN = Infinity / Infinity
 
##########################################################
##########################################################
##########################################################

def is_nan(x):
    return type(x) is float and x != x

##########################################################
##########################################################
##########################################################

def is_inf(x):
    return x in (Infinity, -Infinity)

##########################################################
##########################################################
##########################################################

"""
"""
def removeWeatherItems(list):
    
    new_list = []
    for p in list:
        if (p.getName() not in ("precip", "globrad", "tmin", "tmax", "tavg", "wind", "relhumid", "vs_SoilMoisture")):
            new_list.append(p)    
            
    return new_list

##########################################################
##########################################################
##########################################################

""" 
"""
def getNominalList(list):
    
    new_list = []
    for p in list:        
        new_list.append(float(p.getNominal()))    
        
    return new_list

##########################################################
##########################################################
##########################################################

        
"""
Writes output effects into several files. For each analysed 
input parameter a new file is created that contains effects
for each simulation run. 
"""     
def saveEffects(parameter_effects, parameter_list, output_path):
    print "save parameter_effects", output_path
    #print parameter_effects
    output_names = getOutputNames()
    parameter_names = []
    file_list = []
    for p in parameter_list:
        file = csv.writer(open(output_path+"/outputs/"+p.getName() +".txt", "w"), delimiter='\t')
        file.writerow(output_names)                  
        file_list.append(file)
        parameter_names.append(p.getName())
    
    
    for proc_list in parameter_effects:
        index =0    
        for parameter_effects in proc_list:
            file = file_list[index]            
            
            for effects in parameter_effects:
                # write effect values of all outputs to file
                #print parameter_names[index], effects
                file.writerow(effects)                
            index = index+1       

##########################################################
##########################################################
##########################################################
        
      
"""
Returns a list of names for the analysed outputs. 
"""  
def getOutputNames():
    
    names = []
    result_ids = monica.sensitivityAnalysisResultIds()
    for id  in result_ids:
        names.append(monica.resultIdInfo(id).shortName)
    return names

##########################################################
##########################################################
##########################################################

"""
Reads parameter configuration that are specified in a file
to a parameter list. For each parameter a name, the minimum
and maximum value of range ist given.
According to these parameter specification the SA will be done.
"""
def readParameterFile(filename):
    parameter_list = []

    print "Reading parameters from file: ", filename
    parameter_file = csv.reader(open(filename,'r'), delimiter=';')

    line_nr = 0
    for line in parameter_file: 
        #print line, len(line)
        if (line_nr==0):
            line_nr = line_nr + 1
            continue

        if (len(line)>0):            
            parameter = SAParameter(line[1],float(line[2]), float(line[3]),float(line[4]))
            parameter_list.append(parameter)
            #parameter.display()
        else:
            print "Error: Cannot read parameter file\"", filename, "\"\n", len(line), "\n", "Line:", line
            sys.exit(-1)
    
    return parameter_list

##########################################################
##########################################################
##########################################################

"""
Reads parameter configuration that are specified in a file
to a parameter list. For each parameter a name, the minimum
and maximum value of range ist given.
According to these parameter specification the SA will be done.
"""
def readNewParameterFile(filename, crop_id):

    print "Reading new parameters from file: ", filename
    
    p_file = open(filename,'r')
    parameter_file = csv.reader(p_file, delimiter=';')

    line_nr = 0
    min_index = 2
    max_index = 3
    crop_index = -1
    parameter_list = []
    for line in parameter_file: 
        #print line
        if (line_nr==0):

            # get column index of the analysed crop
            for i, r in enumerate(line):
                if (i>3 and len(r)>0 and int(r) == crop_id):
                    crop_index = i
                
            line_nr += 1
            continue

        #print line
        if (len(line[crop_index])>0):
            parameter = SAParameter(line[1],float(line[min_index]), float(line[max_index]),float(line[crop_index]))
            parameter_list.append(parameter)
            #parameter.display()
    
    p_file.close()
    print "Have read ", len(parameter_list), "parameters"
    return parameter_list

##########################################################
##########################################################
##########################################################

"""
Returns last element of a list 
"""    
def lastElement(list):
    
    length = len(list)
    if (length>0):
        elem = list[length-1]  
        return elem
    
    return 0 

##########################################################
##########################################################
##########################################################


"""
Returns last element of a list 
"""    
def getMeanOfList(list): 
    if (len(list)>0):   
        return numpy.mean(list)
    else:
        return None

##########################################################
##########################################################
##########################################################

def getOutputId(output_list):
    
    result_ids = monica.sensitivityAnalysisResultIds() 

    name_map = {}
    name_map["Corg"]        = "dailyCorg30"
    name_map["primYield"]   = "primYield"
    #name_map["yield"]       = "primYield"
    name_map["abBiom"]      = "abBiom"
    name_map["earBiomass"]  = "fruitBiomass"
    name_map["rootBiomass"] = "rootBiomass"
    name_map["shootBiomass"]= "shootBiomass"
    name_map["leafBiomass"] = "leafBiomass"
    name_map["lai"]         = "lai"
    name_map["GPP"]         = "dailyGPP"
    name_map["Ra"]          = "dailyRa"
    name_map["abiomNContent"] = "aboveBiomassNContent"
    name_map["cropHeight"]  = "cropHeight"
    name_map["Eta"]         = "dailyETa"
    name_map["soilMoist"]   = "dailySoilMoist90"
    name_map["soilNmin"]    = "dailyNmin90"
    name_map["soilTemp"]    = "dailySoilTemp30"
    name_map["abiomNContent"] = "aboveBiomassNContent"
    name_map["dailyPrimYield"] = "dailyPrimYield"
    name_map["dailyAGB"]    = "dailyAGB"
    name_map["dailyRoot"]   = "dailyRoot"
    name_map["devStage"] = "devStage"



    output_id = []

    for o in output_list:
        short_name = name_map[o]

        for r_id in result_ids:
            if (monica.resultIdInfo(r_id).shortName == short_name):
                output_id.append(r_id)
                break
    return output_id


##########################################################
##########################################################
##########################################################

def get_optimised_trajectories(parameter_list, parameter_grid, start_vector_count, ranges, schrittweite, random_traj_count):

  parameter_number = len(parameter_list)

  print "Parameter_number:", parameter_number

  trajectories = []
  for traj_count in range(random_traj_count):
  
    trajectory = []
    # calculate starting point of the trajectory
    start_vector = []
    index = 0;
    for param in range(parameter_number) : 
      index = int(round(random.uniform(0,ranges-1)))
      start_vector.append(index)
      # values.append(parameter_grid[param][index])
 
    trajectory.append(list(start_vector))
    #print "Start vector", start_vector
    # calculate the rest of the trajectory by changing only one
    # parameter at the time
   
    initialized = [False for i in range(parameter_number)]
    
    # increment every parameter only once
    for index in range(parameter_number):            

      parameter_index = int(round(random.uniform(0,parameter_number-1)))

      # get parameter index randomly but be aware
      # that every parameter can be incremented
      # only for once
      while (initialized[parameter_index] == True):
      
        # determine randomly the next parameter
        parameter_index = int(random.uniform(0,parameter_number))

      # set status, so this parameter will be marked as incremented
      initialized[parameter_index] = True

      # increase0or decrease 1 one element of start vector
      operation = getOperation(start_vector[parameter_index], ranges,schrittweite)
           
      # apply operation to startvector
      if (operation == ADDITION) :
          start_vector[parameter_index]+=schrittweite
      else :
          start_vector[parameter_index]-=schrittweite
      next_point = list(start_vector)

      #print next_point
      trajectory.append(next_point)

    #print trajectory

    trajectories.append(trajectory)

  opt_trajectories = find_trajectories_with_maximised_distance(trajectories, parameter_grid, parameter_number, start_vector_count)

  return opt_trajectories

##########################################################
##########################################################
##########################################################

"""
"""
def find_trajectories_with_maximised_distance(trajectories, parameter_grid, parameter_number, start_vector_count):
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  name = MPI.Get_processor_name()

  print "Find trajectories with maximum distance"
  traj_count = len(trajectories)

  parameter_distance_array = numpy.zeros((traj_count, traj_count))

  # generate list for scattering to nodes  
  trajectory_list =  []  

  trajectory_list = mpi_helper.splitListForNodes(range(traj_count), size)


  local_list = comm.scatter(trajectory_list, root=0)    
    
  print rank, "received: ", local_list

  for m in range(local_list[0],traj_count):
    for l in local_list:##########################################################
      #print rank, "\t", m, l, parameter_distance_array[m,l]
      if ((m != l)):
        # calculate distrance of two trajectories


        dist_ml = 0  

        for point in range(parameter_number+1):
          m_traj = trajectories[m][point]
          l_traj = trajectories[l][point]

          sum_sqrt = 0
          for z in range(parameter_number):  
            dist = (m_traj[z]-l_traj[z])**2
            sum_sqrt += dist

          sum_sqrt = numpy.sqrt(sum_sqrt)
          dist_ml += sum_sqrt
        parameter_distance_array[m,l] = dist_ml
      
        #print rank, m,l, parameter_distance_array 

  distance_array = comm.gather(numpy.array(parameter_distance_array), root=0)

  parameter_distance_array = numpy.zeros((traj_count, traj_count))
  opt_trajectories = []
  if (rank == 0):
    #print rank, len(distance_array)

    
    for d_a in distance_array:
      #print d_a, "\n"
      parameter_distance_array = numpy.add(d_a, numpy.array(parameter_distance_array))
      #print parameter_distance_array, "\n", "\n"

    #print parameter_distance_array

    # pick the trajectories with the most distance in between
    max_traj_indizes = []
    print "search for max"
    while (len(max_traj_indizes)<start_vector_count):
      
      i,j = numpy.unravel_index(parameter_distance_array.argmax(), parameter_distance_array.shape)
      print i,j
      if (parameter_distance_array[i][j] == 0):
        break

      parameter_distance_array[i][j] = 0
      parameter_distance_array[j][i] = 0

      if (i not in max_traj_indizes):
        max_traj_indizes.append(i)
        print len(max_traj_indizes)
      if (j not in max_traj_indizes):
        max_traj_indizes.append(j)
        print len(max_traj_indizes)
    
    #print "MAX", max_traj_indizes
    
    
    for t in max_traj_indizes:
      opt_trajectories.append(trajectories[t])

  return opt_trajectories
  
##########################################################
##########################################################
##########################################################

"""
Calculates number of steps based on the parameter number
and range number; Algorithm for this calculation can be
find in Saltelli 2000
"""
def getNumberOfMorrisSteps(parameter_number, range_number, schrittweite):
    
    delta = schrittweite/(range_number-1.0);
    steps = int(math.pow(range_number,parameter_number-1.0)*(range_number-(delta * (range_number-1.0))));
    #steps = 1
    return int(steps)

##########################################################
##########################################################
##########################################################

""" 
Randomly generates operation for manipulating the
startvector. Returns index of operation that should be
applied.
"""
def getOperation(parameter_index, range_number, schrittweite):
    
    x = random.uniform(0,1)
    operation = int(round(x))
    
    # test if operation is possible
    if (operation == ADDITION) :
        if (parameter_index+schrittweite >= range_number) :
            operation = SUBTRACTION
    else :
        if (parameter_index == 0) :
            operation = ADDITION

    return operation
