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

    for  i in range(size):
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
    
    # go through all parameters in list
    for p in parameter_list:
        
        name = p.getName()
        #print "Apply : ", name, start_vector[index], index #Prozessor", rank,"\t
        #print "Apply : ", len(start_vector), index, type(start_vector[index])
        #print "Apply1 : ", name, start_vector[index] #Prozessor", rank,"\t
        
        # soil temperature parameters
        if (name in ("pt_BaseTemperature", 
                     "pt_InitialSurfaceTemperature",
                     "pt_NTau", 
                     "pt_SoilMoisture")):
            # replace values in user soil temperature parameters
            new_env.centralParameterProvider.userSoilTemperatureParameters.__setattr__(name, start_vector[index])
        
        # climate data for the data accessor object    
        elif (name in ("tmin", 
                       "tmax", 
                       "tavg", 
                       "precip", 
                       "globrad", 
                       "wind", 
                       "sunhours", 
                       "relhumid")):
            
            # replace climate data
            size = new_env.numberOfPossibleSteps()
            new_vector = monica.DoubleVector(size)
            for i in range(0,size):
                new_vector[i] = start_vector[index]

            new_env.addOrReplaceClimateData(name,new_vector)
            
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
         
        # sensitivity analysis parameters
        elif (name in ("p_MeanBulkDensity",
                       "p_MeanFieldCapacity",
                       "p_HeatConductivityFrozen",
                       "p_HeatConductivityUnfrozen",
                       "p_LatentHeatTransfer",
                       "p_ReducedHydraulicConductivity",
                       "vs_FieldCapacity",
                       "vs_Saturation",
                       "vs_PermanentWiltingPoint",
                       "vc_SoilCoverage",
                       "vc_MaxRootingDepth",
                       "vc_RootDiameter",
                       "vs_SoilMoisture", 
                       "vs_SoilTemperature")):   
            
            new_env.centralParameterProvider.sensitivityAnalysisParameters.__setattr__(name, start_vector[index])
            
        # env parameters    
        elif (name in ("vs_Slope")):
            new_env.site.__setattr__(name, start_vector[index])
            
        # crop parameters that have different values for developmental stages            


     
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
            
            # round values, because pathway can be either 1 or 2 (for a method selection)
            if (name == "pc_CarboxylationPathway") :
                val = int(round(start_vector[index] % 2)+1)
                start_vector[index] = val
            
            # iterate through all production processes to changes special parameter
            new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__(name, start_vector[index])
            
            new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
            new_env.setCropRotation(new_crop_rotation)
            
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
                       "pc_ReferenceAlbedo")):
            
            # replace values in user crop parameters
            new_env.centralParameterProvider.userCropParameters.__setattr__(name, start_vector[index])
        
        # organ dependent crop parameters
        elif (name in ("pc_OrganGrowthRespiration",
                       "pc_OrganMaintenanceRespiration")):
            
            size = 4 # assuming there are 4 different organs
            param_vector = monica.DoubleVector(size)
            
            for i in range(0,size):                
                
                # changes all values of vector
                param_vector[i] = start_vector[index]
            
            # iterate through all production processes to changes special parameter
            new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__(name, param_vector)
                
            new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
            new_env.setCropRotation(new_crop_rotation)
        
        
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
            new_env.organic.__setattr__(name, start_vector[index])
            
        # user soil transport parameters
        elif (name in ("pq_DiffusionCoefficientStandard",
                       "pq_AD",
                       "pq_DispersionLength")):
            
            # replace values in user crop parameters
            new_env.centralParameterProvider.userSoilTransportParameters.__setattr__(name, start_vector[index])
            
        # user environment parameters
        elif (name in ("p_LeachingDepth"),
                       "p_MinGroundwaterDepth"):
            

            # replace values in user crop parameters
            new_env.centralParameterProvider.userEnvironmentParameters.__setattr__(name, start_vector[index])
            if (name == "p_MinGroundwaterDepth"):
                new_env.centralParameterProvider.userEnvironmentParameters.__setattr__("p_MaxGroundwaterDepth", start_vector[index]+2)
            
        # organic matter parameters
        elif (name in ("vo_AOM_DryMatterContent",
                       "vo_AOM_NH4Content",
                       "vo_AOM_NO3Content",
                       "vo_AOM_CarbamidContent",
                       "vo_PartAOM_to_AOM_Slow",
                       "vo_PartAOM_to_AOM_Fast",
                       "vo_CN_Ratio_AOM_Slow",
                       "vo_CN_Ratio_AOM_Fast")):
            new_env.centralParameterProvider.sensitivityAnalysisParameters.organic_matter_parameters.__setattr__(name, start_vector[index])                      
            new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
            new_env.setCropRotation(new_crop_rotation)
            


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

        # pc_OptimumTemperature -----------------------------
        m = re.match(r"""pc_OptimumTemperature(\d)""", name)
        if (m):
            stage = m.group(1)
            if (stage != None):
                optimum_temperature_vector[int(stage)-1] = start_vector[index]
                optimum_temperature_changed = True





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
       pc_BaseTemperature new_env.centralParameterProvider.sensitivityAnalysisParameters.crop_parameters.__setattr__("pc_BaseDaylength", base_daylength_vector)

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
        optimum_temperature_changed
        ):

        new_crop_rotation = monica.applySAChanges(env.cropRotation, new_env.centralParameterProvider)
        new_env.setCropRotation(new_crop_rotation)
        
    #print "\n"

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
