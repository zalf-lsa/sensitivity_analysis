# adding the path of monica's python module to search path
# that is defined by the environment variable PYTHONPATH
import sys
sys.path.append('/media/san1_data1/data1/specka/devel/models/monica/python')

import monica


def initializeMonicaSimulation(input_path, output_path):
    print "initializeMonicaSimulation"
    central_parameters = monica.readUserParameterFromDatabase(monica.Env.MODE_HERMES)
    # soil data
    soil_parameters = monica.soilParametersFromHermesFile(1, input_path + "SOIL.txt")
    # climate data
    climateData = monica.climateDataFromHermesFiles(input_path + "MET_BS.0", 1999, 2001)
    # monica.testClimateData(climateData);
    # fruchtfolge
    ff = monica.cropRotationFromHermesFile(input_path + "ROTATION.txt")
    
    # irrigation
    #monica.attachIrrigationApplicationsToCropRotation(ff, input_path + "Irrig.TXT")
    ff = monica.attachFertiliserSA(ff, input_path + "SLAGDUNG.txt")
    
    # build up the monica environment
    env = monica.Env(soil_parameters, central_parameters)
    env.pathToOutputDir = output_path
    siteParams = monica.SiteParameters()
    env.site = siteParams;
    env.da = climateData
    env.setCropRotation(ff)
    env.autoIrrigationParams = monica.AutomaticIrrigationParameters(20, 0.1, 0.01)
    env.nMinUserParams = monica.NMinUserParameters(10.0, 100.0, 30)
    env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)

    # mineral fertiliser
    mineral_fert_id = 1
    env.nMinFertiliserPartition = monica.getMineralFertiliserParametersFromMonicaDB(mineral_fert_id);
    return env    


def initializeEVA2MonicaSimulation(output_path):
  
  config = monica.Eva2SimulationConfiguration()
  config.setOutputPath(output_path)

  config.setFruchtFolge("03")
  config.setVariante(1)
  config.setLocation(11)
  config.setClassification(1)
  config.setStartDate(1,1,2005,1)
  config.setEndDate(31,12,2008,1)   

  standort="ascha"
  config.setLocationName(standort)


  config.setPseudoSimulation(0)
  config.setProfil_number(27) 
    
  config.setFruchtfolgeGlied(1)
  config.setFruchtArt("141")          # Silomais 2005
  config.setFruchtfolgeYear("2005")

  config.setFruchtfolgeGlied(2)
  config.setFruchtArt("172")          # Winterroggen 2005/2006
  config.setFruchtfolgeYear("2006")

  config.setFruchtfolgeGlied(3)
  config.setFruchtArt("160")          # Sudangras 2006
  config.setFruchtfolgeYear("2006")

  config.setFruchtfolgeGlied(4)
  config.setFruchtArt("175")          # Wintertriticale 2006/2007
  config.setFruchtfolgeYear("2007")

  config.setFruchtfolgeGlied(5)
  config.setFruchtArt("182")          # Einj. Weidelgras 2007 
  config.setFruchtfolgeYear("2007")

  config.setFruchtfolgeGlied(6)
  config.setFruchtArt("176")          # Winterweizen 2008
  config.setFruchtfolgeYear("2008")

  env = monica.getEVA2Env(config)
  env.setMode(monica.Env.MODE_SENSITIVITY_ANALYSIS)

  return env      

