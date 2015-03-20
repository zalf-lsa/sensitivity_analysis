# adding the path of monica's python module to search path
# that is defined by the environment variable PYTHONPATH
import sys
sys.path.append('../..')

import monica
import re

"""
Changes values that are tested in the SA directly in 
the datastructure objects that are used in monica simulation.
Returns copy of manipulated env object for running monica simulation. 
"""
def applySAValues(parameter_list, start_vector, env, crop_id=None):
    
    index = 0
    #print env
    new_env = monica.Env(env)
    if (crop_id!=None):    
        new_env.centralParameterProvider.sensitivityAnalysisParameters.sa_crop_id = crop_id

    # initialisation of dev stage vectors
    size = 7 # developmental stages
    organ_count = 4 # number of organs



    daylength_req_vector = monica.DoubleVector(size)
    stage_kc_vector = monica.DoubleVector(size)
    stage_temp_sum_vector = monica.DoubleVector(size)
    vern_req_vector = monica.DoubleVector(size)
    spec_leaf_area_vector = monica.DoubleVector(size)
    base_daylength_vector = monica.DoubleVector(size)
    drought_stress_threshold_vector = monica.DoubleVector(size)
    critical_oxygen_content_vector = monica.DoubleVector(size)
    stage_max_root_nconcentration_vector = monica.DoubleVector(size)
    base_temperature_vector = monica.DoubleVector(size)
    optimum_temperature_vector = monica.DoubleVector(size)

    organ_growth_respiration_vector = monica.DoubleVector(organ_count)
    organ_maintenance_respiration_vector = monica.DoubleVector(organ_count)

    assimilate_partitioning_vector = []
    for n in range(size):
      assimilate_partitioning_vector.append(monica.DoubleVector(organ_count))

    for stage_index in range(size):
      for organ_index in range(organ_count):
        assimilate_partitioning_vector[stage_index][organ_index] = -9999

    for i in range(size):
        daylength_req_vector[i] = -9999
        stage_kc_vector[i] = -9999
        stage_temp_sum_vector[i] = -9999
        vern_req_vector[i] = -9999
        spec_leaf_area_vector[i] = -9999
        base_daylength_vector[i] = -9999
        drought_stress_threshold_vector[i] = -9999
        critical_oxygen_content_vector[i] = -9999
        stage_max_root_nconcentration_vector[i] = -9999
        base_temperature_vector[i] = -9999
        optimum_temperature_vector[i] = -9999

    for i in range(organ_count):
        organ_growth_respiration_vector[i] = -9999
        organ_maintenance_respiration_vector[i] = -9999

    daylength_req_changed = False
    stage_kc_changed = False
    stage_temp_sum_changed = False
    vern_req_changed = False
    spec_leaf_area_changed = False
    base_daylength_changed = False
    drought_stress_threshold_changed = False
    critical_oxygen_content_changed = False
    stage_max_root_nconcentration_changed = False
    base_temperature_changed = False
    optimum_temperature_changed = False
    organ_growth_respiration_changed = False
    organ_maintenance_respiration_changed = False        
    assimilate_partitioning_changed = False

    # go through all parameters in list
    for p in parameter_list:
        
        name = p.getName()
        #print "Apply : ", name, start_vector[index], len(start_vector), index
        
        # soil temperature parameters
        if (name in ("pt_BaseTemperature", 
                     "pt_InitialSurfaceTemperature",
                     "pt_NTau", 
                     "pt_SoilMoisture")):
            # replace values in user soil temperature parameters
            new_env.centralParameterProvider.userSoilTemperatureParameters.__setattr__(name, start_vector[index])
        
        # climate data for the data accessor object    
#        elif (name in ("tmin", 
#                       "tmax", 
#                       "tavg", 
#                       "precip", 
#                       "globrad", 
#                       "wind", 
#                       "sunhours", 
#                       "relhumid")):
#            
#            # replace climate data
#            size = new_env.numberOfPossibleSteps()
#            new_vector = monica.DoubleVector(size)
#            for i in range(0,size):
#                new_vector[i] = start_vector[index]
#
#            new_env.addOrReplaceClimateData(name,new_vector)
            
        # soil moisture parameters
        elif (name in ("pm_SnowMeltTemperature",
                       "pm_SnowAccumulationThresholdTemperature",           
                       "pm_TemperatureLimitForLiquidWater",
                       "pm_CorrectionRain", 
                       "pm_CorrectionSnow",
                       "pm_RefreezeTemperature",
                       "pm_NewSnowDensityMin",
                       "pm_SnowMaxAdditionalDensity",
                       "pm_SnowPacking", 
                       "pm_SnowRetentionCapacityMin",
                       "pm_SnowRetentionCapacityMax",
                       "pm_SurfaceRoughness",
                       "pm_MaxPercolationRate",
                       "pm_HydraulicConductivityRedux",
                       "pm_GroundwaterDischarge")):
            
            # replace values of user soil moisture parameters
            new_env.centralParameterProvider.userSoilMoistureParameters.__setattr__(name, start_vector[index])
         
     
        # normal crop parameters            
        elif (name in ("pc_InitialKcFactor",
                       "pc_StageAtMaxHeight",
                       "pc_CropHeightP1",
                       "pc_CropHeightP2",
                       "pc_LuxuryNCoeff",
                       "pc_ResidueNRatio",
                       "pc_CropSpecificMaxRootingDepth",
                       "pc_RootPenetrationRate",
                       "pc_RootGrowthLag",
                       "pc_InitialRootingDepth",
                       "pc_RootFormFactor",
                       "pc_MaxNUptakeParam",
                       "pc_CarboxylationPathway",  
                       "pc_DefaultRadiationUseEfficiency",
                       "pc_MaxAssimilationRate",
                       "pc_MaxCropDiameter",
                       "pc_MinimumNConcentration",
                       "pc_NConcentrationB0",
                       "pc_NConcentrationPN",
                       "pc_NConcentrationRoot",
                       "pc_PlantDensity",
                       "pc_ResidueNRatio",
                       "pc_MinimumTemperatureForAssimilation",
                       "pc_NConcentrationAbovegroundBiomass",
                       "pc_MaxCropHeight",
                       "pc_SamplingDepth",
                       "pc_TargetNSamplingDepth",
                       "pc_TargetN30",
                       "pc_StageAtMaxDiameter",
                       "pc_HeatSumIrrigationStart",
                       "pc_HeatSumIrrigationEnd",
                       "pc_RootDistributionParam",
                       "pc_MinimumTemperatureRootGrowth",
                       "pc_SpecificRootLength",
                       "pc_CriticalTemperatureHeatStress",
                       "pc_LimitingTemperatureHeatStress",
                       "pc_BeginSensitivePhaseHeatStress",
                       "pc_EndSensitivePhaseHeatStress",
                       "pc_DroughtImpactOnFertilityFactor",
                       "pc_AssimilateReallocation",
                       "pc_LT50cultivar",
                       "pc_FrostHardening;",
                       "pc_FrostDehardening",
                       "pc_LowTemperatureExposure",
                       "pc_RespiratoryStress"                       
                        )):

            crop_parameters = monica.getCropParameters(crop_id, new_env.cropRotation)
            crop_parameters.__setattr__(name, start_vector[index])
            new_crop_rotation = monica.applyNewCropParameters(crop_id, env.cropRotation, crop_parameters); 
            new_env.setCropRotation(new_crop_rotation)



            # iterate through all production processes to changes special parameter
            #new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__(name, start_vector[index])
            
            #new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
            #new_env.setCropRotation(new_crop_rotation)
            
        # user crop parameters
        elif (name in ("pc_MaintenanceRespirationParameter1",
                       "pc_MaintenanceRespirationParameter2",
                       "pc_GrowthRespirationParameter1",
                       "pc_GrowthRespirationParameter2",
                       "pc_CanopyReflectionCoefficient",
                       "pc_SaturationBeta",
                       "pc_StomataConductanceAlpha",
                       "pc_MinimumAvailableN",
                       "pc_MinimumNConcentrationRoot",
                       "pc_Tortuosity",
                       "pc_MaxCropNDemand",
                       "pc_ReferenceAlbedo", 
                       "pc_ReferenceLeafAreaIndex",
                       "pc_ReferenceMaxAssimilationRate"
                       )):

            # replace values in user crop parameters
            new_env.centralParameterProvider.userCropParameters.__setattr__(name, start_vector[index])
        
       
        
        
        # organic parameters parameters
        elif (name in ("po_SOM_SlowDecCoeffStandard",
                       "po_SOM_FastDecCoeffStandard",
                       "po_SMB_SlowMaintRateStandard",
                       "po_SMB_FastMaintRateStandard",
                       "po_SMB_SlowDeathRateStandard",
                       "po_SMB_FastDeathRateStandard",
                       "po_SMB_UtilizationEfficiency",
                       "po_SOM_SlowUtilizationEfficiency",
                       "po_SOM_FastUtilizationEfficiency",
                       "po_AOM_SlowUtilizationEfficiency",
                       "po_AOM_FastUtilizationEfficiency",
                       "po_AOM_FastMaxC_to_N",
                       "po_PartSOM_Fast_to_SOM_Slow",
                       "po_PartSMB_Slow_to_SOM_Fast",
                       "po_PartSMB_Fast_to_SOM_Fast",
                       "po_PartSOM_to_SMB_Slow",
                       "po_PartSOM_to_SMB_Fast",
                       "po_CN_Ratio_SMB",
                       "po_LimitClayEffect",
                       "po_NitrificationRateCoeffStandard",
                       "po_TransportRateCoeff",
                       "po_SpecAnaerobDenitrification",
                       "po_ImmobilisationRateCoeffNO3",
                       "po_ImmobilisationRateCoeffNH4",
                       "po_Denit1",
                       "po_Denit2",
                       "po_Denit3",
                       "po_UreaMolecularWeight",
                       "po_Urea_to_N",
                       "po_HydrolysisKM",
                       "po_ActivationEnergy",
                       "po_HydrolysisP1",
                       "po_HydrolysisP2",
                       "po_NH4MolecularWeight",
                       "po_H2OIonConcentration",
                       "po_AtmosphericResistance"
                       "po_AmmoniaOxidationRateCoeffStandard",
                       "po_NitriteOxidationRateCoeffStandard"
                       )):



            # replace values in user crop parameters
            new_env.centralParameterProvider.userSoilOrganicParameters.__setattr__(name, start_vector[index])
            
        # user soil transport parameters
        elif (name in ("pq_DiffusionCoefficientStandard",
                       "pq_AD",
                       "pq_DispersionLength")):
            
            # replace values in user crop parameters
            new_env.centralParameterProvider.userSoilTransportParameters.__setattr__(name, start_vector[index])
            


        # daylength requirement -----------------------------
        m = re.match(r"""pc_DaylengthRequirement(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                daylength_req_vector[int(stage)-1] = start_vector[index]
                daylength_req_changed = True



        # stage_kc_vector -----------------------------
        m = re.match(r"""pc_StageKcFactor(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                stage_kc_vector[int(stage)-1] = start_vector[index]
                stage_kc_changed = True



        # stage_temp_sum -----------------------------
        m = re.match(r"""pc_StageTemperatureSum(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                stage_temp_sum_vector[int(stage)-1] = start_vector[index]
                stage_temp_sum_changed = True


        # pc_VernalisationRequirement -----------------------------
        m = re.match(r"""pc_VernalisationRequirement(\d)""", name)
        if (m and (crop_id==1 or crop_id == 4 or crop_id ==9)):
            stage = m.group(1)
            if (stage != None):
                vern_req_vector[int(stage)-1] = start_vector[index]
                vern_req_changed = True

        # pc_SpecificLeafArea -----------------------------
        m = re.match(r"""pc_SpecificLeafArea(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                spec_leaf_area_vector[int(stage)-1] = start_vector[index]
                spec_leaf_area_changed = True

        # pc_BaseDaylength -----------------------------
        m = re.match(r"""pc_BaseDaylength(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                base_daylength_vector[int(stage)-1] = start_vector[index]
                base_daylength_changed = True


        # pc_DroughtStressThreshold -----------------------------
        m = re.match(r"""pc_DroughtStressThreshold(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                drought_stress_threshold_vector[int(stage)-1] = start_vector[index]
                drought_stress_threshold_changed = True


        # pc_CriticalOxygenContent -----------------------------
        m = re.match(r"""pc_CriticalOxygenContent(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                critical_oxygen_content_vector[int(stage)-1] = start_vector[index]
                critical_oxygen_content_changed = True


        # pc_StageMaxRootNConcentration -----------------------------
        m = re.match(r"""pc_StageMaxRootNConcentration(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                stage_max_root_nconcentration_vector[int(stage)-1] = start_vector[index]
                stage_max_root_nconcentration_changed = True

        # pc_BaseTemperature -----------------------------
        m = re.match(r"""pc_BaseTemperature(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                base_temperature_vector[int(stage)-1] = start_vector[index]
                base_temperature_changed = True

        m = re.match(r"""pc_OptimumTemperature(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                optimum_temperature_vector[int(stage)-1] = start_vector[index]
                optimum_temperature_changed = True

        # pc_OrganGrowthRespiration
        m = re.match(r"""pc_OrganGrowthRespiration(\d)""", name)
        if (m):
            organ_id = m.group(1)
            if (organ_id != None):
                organ_growth_respiration_vector[int(organ_id)-1] = start_vector[index]
                organ_growth_respiration_changed = True

        # pc_OrganMaintenanceRespiration
        m = re.match(r"""pc_OrganMaintenanceRespiration(\d)""", name)
        if (m):
            organ_id = m.group(1)
            if (organ_id != None):
                organ_maintenance_respiration_vector[int(organ_id)-1] = start_vector[index]
                organ_maintenance_respiration_changed = True

        # pc_AssimilatePartitioningDevStage1Organ1
        m = re.match(r"""pc_AssimilatePartitioningDevStage(\d)Organ(\d)""", name)
        if (m):
            dev_stage = m.group(1)
            organ_id = m.group(2)
            if (organ_id != None and dev_stage != None):
                assimilate_partitioning_vector[int(dev_stage)-1][int(organ_id)-1]  = start_vector[index]

                assimilate_partitioning_changed = True

        index = index+1 # index for the start_vector values
    

    # end of for p in parameters

    if (daylength_req_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_DaylengthRequirement", daylength_req_vector)            

    if (stage_kc_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_StageKcFactor", stage_kc_vector)            

    if (stage_temp_sum_changed):           
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_StageTemperatureSum", stage_temp_sum_vector)        

    if (vern_req_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_VernalisationRequirement", vern_req_vector)       

    if (spec_leaf_area_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_SpecificLeafArea", spec_leaf_area_vector)       

    if (base_daylength_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_BaseDaylength", base_daylength_vector)

    if (drought_stress_threshold_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_DroughtStressThreshold", drought_stress_threshold_vector)              

    if (critical_oxygen_content_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_BaseDaylength", critical_oxygen_content_vector)       

    if (stage_max_root_nconcentration_changed):            
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_StageMaxRootNConcentration", stage_max_root_nconcentration_vector)      

    if (base_temperature_changed):  
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_BaseTemperature", base_temperature_vector)   

    if (optimum_temperature_changed):  
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_OptimumTemperature", optimum_temperature_vector)    

    if (organ_growth_respiration_changed):  
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_OrganGrowthRespiration", organ_growth_respiration_vector)    
 
    if (organ_maintenance_respiration_changed):  
        new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_OrganMaintenanceRespiration", organ_maintenance_respiration_vector)   



    if (daylength_req_changed or 
        stage_kc_changed or 
        stage_temp_sum_changed or 
        vern_req_changed or 
        spec_leaf_area_changed or 
        base_daylength_changed or
        drought_stress_threshold_changed or
        critical_oxygen_content_changed or
        stage_max_root_nconcentration_changed or
        base_temperature_changed or
        optimum_temperature_changed or
        organ_maintenance_respiration_changed or
        organ_growth_respiration_changed or
        assimilate_partitioning_changed):

        new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
        new_env.setCropRotation(new_crop_rotation)
        
    #print "\n"

    if (assimilate_partitioning_changed):  
        for i in range(size):
          new_crop_rotation = monica.setAssimilatePartitioningCoefficient(i, assimilate_partitioning_vector[i], env.cropRotation, new_env.centralParameterProvider)
          new_env.setCropRotation(new_crop_rotation)


    return new_env

#"""
#Returns default value of a parameter
#"""
#def getNominalValue(parameter_list, parameter_name):
#    print "GetNominalValue: ", parameter_name
#    for p in parameter_list:
#        if (p.getName() == parameter_name):
#            print p.getNominal()
#            return p.getNominal()

#    return None
