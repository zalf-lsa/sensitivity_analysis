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

#include "boost/foreach.hpp"

#define LOKI_OBJECT_LEVEL_THREADING
#include "loki/Threads.h"

using namespace std;
using namespace Monica;


/**
 * @brief Lockable object
 */
struct L: public Loki::ObjectLevelLockable<L> {};

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

  BOOST_FOREACH(ProductionProcess pp, ff)
  {
    CropPtr crop = pp.crop();

    const SensitivityAnalysisParameters& saps =  centralParameterProvider.sensitivityAnalysisParameters;

    if (saps.sa_crop_id != crop->id() && saps.sa_crop_id>0){
//        cout << "Do not apply SA values" << endl;
      continue;
    } else {
//        cout << "CropIds: SA:\t"<< saps.sa_crop_id << "\tMONICA:\t" << crop->id() << endl;
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

    if (saps.crop_parameters.pc_InitialKcFactor != UNDEFINED) {
      cps->pc_InitialKcFactor = saps.crop_parameters.pc_InitialKcFactor;
    }

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

    // pc_StageAtMaxHeight
    if (saps.crop_parameters.pc_StageAtMaxHeight != UNDEFINED) {
      cps->pc_StageAtMaxHeight = saps.crop_parameters.pc_StageAtMaxHeight;
    }

    if (saps.crop_parameters.pc_CropHeightP1 != UNDEFINED) {
      cps->pc_CropHeightP1 = saps.crop_parameters.pc_CropHeightP1;
    }

    // pc_CropHeightP2
    if (saps.crop_parameters.pc_CropHeightP2 != UNDEFINED) {
      cps->pc_CropHeightP2 = saps.crop_parameters.pc_CropHeightP2;
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
        double sa_value = saps.crop_parameters.pc_StageTemperatureSum.at(i);
        double default_value = cps->pc_StageTemperatureSum.at(i);
        if (sa_value == -9999) {
            new_values.push_back(default_value);
        } else {
            new_values.push_back(sa_value);
        }
      }

      cps->pc_StageTemperatureSum = new_values;
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




    // pc_LuxuryNCoeff
    if (saps.crop_parameters.pc_LuxuryNCoeff != UNDEFINED) {
      cps->pc_LuxuryNCoeff = saps.crop_parameters.pc_LuxuryNCoeff;
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

    // pc_ResidueNRatio
    if (saps.crop_parameters.pc_ResidueNRatio != UNDEFINED) {
      cps->pc_ResidueNRatio = saps.crop_parameters.pc_ResidueNRatio;
    }

    // pc_CropSpecificMaxRootingDepth
    if (saps.crop_parameters.pc_CropSpecificMaxRootingDepth != UNDEFINED) {
      cps->pc_CropSpecificMaxRootingDepth = saps.crop_parameters.pc_CropSpecificMaxRootingDepth;
    }

    // pc_RootPenetrationRate
    if (saps.crop_parameters.pc_RootPenetrationRate != UNDEFINED) {
      cps->pc_RootPenetrationRate = saps.crop_parameters.pc_RootPenetrationRate;
    }

    // pc_RootGrowthLag
    if (saps.crop_parameters.pc_RootGrowthLag != UNDEFINED) {
      cps->pc_RootGrowthLag = saps.crop_parameters.pc_RootGrowthLag;
    }

    // pc_InitialRootingDepth
    if (saps.crop_parameters.pc_InitialRootingDepth != UNDEFINED) {
      cps->pc_InitialRootingDepth = saps.crop_parameters.pc_InitialRootingDepth;
    }

    // pc_RootFormFactor
    if (saps.crop_parameters.pc_RootFormFactor != UNDEFINED) {
      cps->pc_RootFormFactor = saps.crop_parameters.pc_RootFormFactor;
    }

    // pc_MaxNUptakeParam
    if (saps.crop_parameters.pc_MaxNUptakeParam != UNDEFINED) {
      cps->pc_MaxNUptakeParam = saps.crop_parameters.pc_MaxNUptakeParam;
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

    // pc_CarboxylationPathway
    if (saps.crop_parameters.pc_CarboxylationPathway > -9999) { // UNDEFINED
      cps->pc_CarboxylationPathway = saps.crop_parameters.pc_CarboxylationPathway;
    }

    // pc_DefaultRadiationUseEfficiency
    if (saps.crop_parameters.pc_DefaultRadiationUseEfficiency != UNDEFINED) {
      cps->pc_DefaultRadiationUseEfficiency = saps.crop_parameters.pc_DefaultRadiationUseEfficiency;
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

    // pc_MaxAssimilationRate
    if (saps.crop_parameters.pc_MaxAssimilationRate != UNDEFINED) {
      cps->pc_MaxAssimilationRate = saps.crop_parameters.pc_MaxAssimilationRate;
    }

    // pc_MaxCropDiameter
    if (saps.crop_parameters.pc_MaxCropDiameter != UNDEFINED) {
      cps->pc_MaxCropDiameter = saps.crop_parameters.pc_MaxCropDiameter;
    }

    // pc_MinimumNConcentration
    if (saps.crop_parameters.pc_MinimumNConcentration != UNDEFINED) {
      cps->pc_MinimumNConcentration = saps.crop_parameters.pc_MinimumNConcentration;
    }

    // pc_NConcentrationB0
    if (saps.crop_parameters.pc_NConcentrationB0 != UNDEFINED) {
      cps->pc_NConcentrationB0 = saps.crop_parameters.pc_NConcentrationB0;
    }

    // pc_NConcentrationPN
    if (saps.crop_parameters.pc_NConcentrationPN != UNDEFINED) {
      cps->pc_NConcentrationPN = saps.crop_parameters.pc_NConcentrationPN;
    }

    // pc_NConcentrationRoot
    if (saps.crop_parameters.pc_NConcentrationRoot != UNDEFINED) {
      cps->pc_NConcentrationRoot = saps.crop_parameters.pc_NConcentrationRoot;
    }

    // pc_OrganGrowthRespiration
    if (saps.crop_parameters.pc_OrganGrowthRespiration.size() > 0) {
      cps->pc_OrganGrowthRespiration = saps.crop_parameters.pc_OrganGrowthRespiration;
    }

    // pc_OrganMaintenanceRespiration
    if (saps.crop_parameters.pc_OrganMaintenanceRespiration.size() > 0) {
      cps->pc_OrganMaintenanceRespiration = saps.crop_parameters.pc_OrganMaintenanceRespiration;
    }

    // pc_PlantDensity
    if (saps.crop_parameters.pc_PlantDensity != UNDEFINED) {
      cps->pc_PlantDensity = saps.crop_parameters.pc_PlantDensity;
    }

    // pc_ResidueNRatio
    if (saps.crop_parameters.pc_ResidueNRatio != UNDEFINED) {
      cps->pc_ResidueNRatio = saps.crop_parameters.pc_ResidueNRatio;
    }


    //cout << cps->toString().c_str() << endl;
    crop->setCropParameters(cps);
    new_ff.push_back(pp);
  }

  return ff;
}
