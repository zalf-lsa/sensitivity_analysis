/**
Authors: 
Xenia Specka <xenia.specka@zalf.de>

Maintainers: 
Currently maintained by the authors.

This file is part of the MONICA model. 
Copyright (C) 2007-2013, Leibniz Centre for Agricultural Landscape Research (ZALF)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "sensitivity_analysis_interface.h"


using namespace std;
using namespace Monica;




/**
 * Check, if some parameters should be changed according to sensitivity analysis simulation.
 * Replace old parameters with new values from sensitivity analysis object
 * @param ff
 * @param centralParameterProvider
 * @return New crop rotation vector
 */
std::vector<ProductionProcess>
Monica::applySAChanges(std::vector<ProductionProcess> ff,
                           const CentralParameterProvider &centralParameterProvider) 
{
//  cout << "Apply SA values method" << endl;
  std::vector<ProductionProcess> new_ff;

  
  for (ProductionProcess pp : ff)
  {
    CropPtr crop = pp.crop();

    const SensitivityAnalysisParameters& saps =  centralParameterProvider.sensitivityAnalysisParameters;

    if (saps.sa_crop_id != crop->id() && saps.sa_crop_id>0){
        //cout << "Do not apply SA values" << endl;
      continue;
    } else {
        //cout << "CropIds: SA:\t"<< saps.sa_crop_id << "\tMONICA:\t" << crop->id() << endl;
    }

    CropParameters* cps = new CropParameters((*crop->cropParameters()));

    
    
    
  
    // pc_DaylengthRequirement
    if (saps.crop_parameters.pc_DaylengthRequirement.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_DaylengthRequirement.size(); i++) {
        double sa_value = saps.crop_parameters.pc_DaylengthRequirement.at(i);
        double default_value = cps->pc_DaylengthRequirement.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_DaylengthRequirement = new_values;
    }

    // pc_VernalisationRequirement
    if (saps.crop_parameters.pc_VernalisationRequirement.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_VernalisationRequirement.size(); i++) {
        double sa_value = saps.crop_parameters.pc_VernalisationRequirement.at(i);
        double default_value = cps->pc_VernalisationRequirement.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }


      cps->pc_VernalisationRequirement = new_values;
    }

    // CriticalOxygenContent
    if (saps.crop_parameters.pc_CriticalOxygenContent.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_CriticalOxygenContent.size(); i++) {
        double sa_value = saps.crop_parameters.pc_CriticalOxygenContent.at(i);
        double default_value = cps->pc_CriticalOxygenContent.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_CriticalOxygenContent = new_values;
    }

    
    // pc_StageKcFactor
    if (saps.crop_parameters.pc_StageKcFactor.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_StageKcFactor.size(); i++) {
        double sa_value = saps.crop_parameters.pc_StageKcFactor.at(i);
        double default_value = cps->pc_StageKcFactor.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_StageKcFactor = new_values;
    }

    // pc_SpecificLeafArea
    if (saps.crop_parameters.pc_SpecificLeafArea.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_SpecificLeafArea.size(); i++) {
        double sa_value = saps.crop_parameters.pc_SpecificLeafArea.at(i);
        double default_value = cps->pc_SpecificLeafArea.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_SpecificLeafArea = new_values;
    }

    // pc_StageTemperatureSum
    if (saps.crop_parameters.pc_StageTemperatureSum.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_StageTemperatureSum.size(); i++) {
          //cout << "C++ Apply: i " << i << "\t" <<  saps.crop_parameters.pc_StageTemperatureSum.at(i) << endl;
        double sa_value = saps.crop_parameters.pc_StageTemperatureSum.at(i);
        double default_value = cps->pc_StageTemperatureSum.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_StageTemperatureSum = new_values;
      
      //for (int sum : new_values) {
      //    cout << "TempSum: " << sum << endl;
          
      //  }
    }

    // pc_BaseTemperature
    if (saps.crop_parameters.pc_BaseTemperature.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_BaseTemperature.size(); i++) {
        double sa_value = saps.crop_parameters.pc_BaseTemperature.at(i);
        double default_value = cps->pc_BaseTemperature.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_BaseTemperature = new_values;
    }
    
    // pc_OptimumTemperature
    if (saps.crop_parameters.pc_OptimumTemperature.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_OptimumTemperature.size(); i++) {
        double sa_value = saps.crop_parameters.pc_OptimumTemperature.at(i);
        double default_value = cps->pc_OptimumTemperature.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_OptimumTemperature = new_values;
    }
    
    // pc_OrganGrowthRespiration
    if (saps.crop_parameters.pc_OrganGrowthRespiration.size() > 0) {
      std::vector<double> new_values;
      
      for (unsigned int i=0; i<cps->pc_OrganGrowthRespiration.size(); i++) {
        double sa_value = saps.crop_parameters.pc_OrganGrowthRespiration.at(i);
        double default_value = cps->pc_OrganGrowthRespiration.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_OrganGrowthRespiration = new_values;
    }

    // pc_OrganMaintenanceRespiration
    if (saps.crop_parameters.pc_OrganMaintenanceRespiration.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_OrganMaintenanceRespiration.size(); i++) {
        double sa_value = saps.crop_parameters.pc_OrganMaintenanceRespiration.at(i);
        double default_value = cps->pc_OrganMaintenanceRespiration.at(i);
        
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_OrganMaintenanceRespiration = new_values;
    }

    
    // pc_StageMaxRootNConcentration
    if (saps.crop_parameters.pc_StageMaxRootNConcentration.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_StageMaxRootNConcentration.size(); i++) {
        double sa_value = saps.crop_parameters.pc_StageMaxRootNConcentration.at(i);
        double default_value = cps->pc_StageMaxRootNConcentration.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_StageMaxRootNConcentration = new_values;
    }

   
    // pc_BaseDaylength
    if (saps.crop_parameters.pc_BaseDaylength.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_BaseDaylength.size(); i++) {
        double sa_value = saps.crop_parameters.pc_BaseDaylength.at(i);
        double default_value = cps->pc_BaseDaylength.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_BaseDaylength = new_values;
    }

    // pc_DroughtStressThreshold {
    if (saps.crop_parameters.pc_DroughtStressThreshold.size() > 0) {
      std::vector<double> new_values;
      for (unsigned int i=0; i<cps->pc_DroughtStressThreshold.size(); i++) {
        double sa_value = saps.crop_parameters.pc_DroughtStressThreshold.at(i);
        double default_value = cps->pc_DroughtStressThreshold.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_DroughtStressThreshold = new_values;
    }
    

    //cout << cps->toString().c_str() << endl;
    crop->setCropParameters(cps);
    new_ff.push_back(pp);
  }

  return ff;
}

/**
 * Interface method to apply sensitivity analysis values to
 * assimilate partitioning coefficient. There have been problems using
 * a vector<vector<double>> from python scripts as it was not wrapped
 * properly by SWIG. Swig returned a tuple which could not be changed/access
 * by the python scripts. Therefore the need of this interface method.
 * 
 * Assigns the values of assimilate partitioning separated for each developmental
 * stage.
 */
std::vector<ProductionProcess> 
Monica::setAssimilatePartitioningCoefficient(
    int dev_stage, 
    std::vector<double> partitioning_coefficient,
    std::vector<ProductionProcess> ff,
    const CentralParameterProvider &centralParameterProvider)
{
    std::vector<ProductionProcess> new_ff;

  for(ProductionProcess pp : ff)
  {
    CropPtr crop = pp.crop();

    const SensitivityAnalysisParameters& saps =  centralParameterProvider.sensitivityAnalysisParameters;

    if (saps.sa_crop_id != crop->id() && saps.sa_crop_id>0){
        //cout << "Do not apply SA values" << endl;
      continue;
    } else {
        //cout << "CropIds: SA:\t"<< saps.sa_crop_id << "\tMONICA:\t" << crop->id() << endl;
    }

    CropParameters* cps = new CropParameters((*crop->cropParameters()));
    
    std::vector<std::vector<double> > new_pc_AssimilatePartitioningCoeff;
    new_pc_AssimilatePartitioningCoeff.resize(cps->pc_NumberOfDevelopmentalStages,
                                        std::vector<double>(cps->pc_NumberOfOrgans));
    
    for (unsigned int stage_index = 0; stage_index < cps->pc_AssimilatePartitioningCoeff.size(); stage_index ++) {
      
        if (stage_index == dev_stage) {
            
            std::vector<double> old_organ_coefficients = cps->pc_AssimilatePartitioningCoeff[stage_index];    
            std::vector<double> new_values;
            
            for (unsigned int organ_i=0; organ_i<old_organ_coefficients.size(); organ_i++) {
                
                double sa_value = partitioning_coefficient.at(organ_i);
                double default_value = old_organ_coefficients.at(organ_i);
                if (sa_value == -9999) {
                    new_values.push_back(default_value);
                } else {
                    new_values.push_back(sa_value);
                }
             }

            new_pc_AssimilatePartitioningCoeff[stage_index] = new_values;
            
        } else  {
            // use old values
            new_pc_AssimilatePartitioningCoeff[stage_index] = cps->pc_AssimilatePartitioningCoeff[stage_index];
        }
    }
    
    cps->pc_AssimilatePartitioningCoeff = new_pc_AssimilatePartitioningCoeff;
  
    //cout << cps->toString().c_str() << endl;
    crop->setCropParameters(cps);
    new_ff.push_back(pp);
  
  } // for

  return ff;

}


CropParameters* Monica::getCropParameters(int id, std::vector<ProductionProcess> ff)
{
  for (ProductionProcess pp : ff)
  {
    CropPtr crop = pp.crop();

    if (id == crop->id()){
      //cout << "return crop parameters" << endl;
      return new CropParameters((*crop->cropParameters()));
    }
  }
  return NULL;   

}


std::vector<ProductionProcess> 
Monica::applyNewCropParameters(int crop_id, std::vector<ProductionProcess> ff, CropParameters *cps) 
{
  std::vector<ProductionProcess> new_ff;

  for (ProductionProcess pp : ff)
  {
    CropPtr crop = pp.crop();

    
    if (crop_id == crop->id()){
      //cout << "Apply new crop parameters" << endl;
      crop->setCropParameters(cps);
    }
    new_ff.push_back(pp);      
  }
  return new_ff;
}


//------------------------------------------------------------------------------
const vector<int>& Monica::sensitivityAnalysisResultIds()
{
  static ResultId ids[] =
  {
    //primaryYield,                   // done
    //biomassNContent,
    aboveGroundBiomass
    //aboveBiomassNContent,
    //avg0_90cmSoilMoisture,         // done
    //mean90cmMonthlyAvgWaterContent, // done
    //yearlySumGroundWaterRecharge,    
    //yearlySumNLeaching,
    //sum90cmYearlyNatDay,
    //sumETaPerCrop,                  // done
    //avg30cmMonthlyAvgCorg
    
  };

  //static vector<int> v(ids, ids+2);
  static vector<int> v(ids, ids+1);

  return v;
}

