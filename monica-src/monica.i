%module monica
%include "std_string.i"
%include "std_vector.i"
%include "boost_shared_ptr.i"
%include "std_map.i"

%shared_ptr(Monica::WorkStep)
%shared_ptr(Monica::Seed)
%shared_ptr(Monica::Harvest)
%shared_ptr(Monica::Cutting)
%shared_ptr(Monica::MineralFertiliserApplication)
%shared_ptr(Monica::OrganicFertiliserApplication)
%shared_ptr(Monica::TillageApplication)
%shared_ptr(Monica::IrrigationApplication)

%{
//#define SWIG_COMPILATION
/* Includes the header in the wrapper code */
#include "../../monica/src/simulation.h"
#include "../../monica/src/monica-parameters.h"
#include "../../monica/src/monica.h"
#include "../../util/soil/conversion.h"
#include "sensitivity_analysis_interface.h"
%}


// Instantiate templates used by example
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<std::string>;
   %template(PPVector) vector<Monica::ProductionProcess>;   
   %nodefaultdtor Climate::DataAccessor;
   %nodefaultdtor Monica::CropPtr;
   %template(IntDoubleMap) map<int, double>;
   
}


/* Parse the header file to generate wrappers */
%include "../../monica/src/simulation.h"
%include "../../monica/src/monica-typedefs.h"
%include "../../monica/src/monica-parameters.h"
%include "../../monica/src/soiltemperature.h"
%include "../../monica/src/soilmoisture.h"
%include "../../monica/src/soilorganic.h"
%include "../../monica/src/soiltransport.h"
%include "../../monica/src/monica.h"
%include "../../util/soil/conversion.h"
%include "../../util/tools/date.h"
%include "sensitivity_analysis_interface.h"
