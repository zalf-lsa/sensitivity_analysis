%module monica
%include "std_string.i"
%include "std_vector.i"

%{
#define SWIG_COMPILATION
/* Includes the header in the wrapper code */
#include <src/simulation.h>
#include <src/monica-parameters.h>
#include <src/monica.h>
#include <src/eva_methods.h>
#include <src/conversion.h>
#include <tools/date.h>
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
   
}


/* Parse the header file to generate wrappers */
%include "src/simulation.h"
%include <src/monica-typedefs.h>
%include <src/monica-parameters.h>
%include <src/soiltemperature.h>
%include <src/soilmoisture.h>
%include <src/soilorganic.h>
%include <src/soiltransport.h>
%include <src/monica.h>
%include <src/eva_methods.h>
%include <src/conversion.h>
%include <tools/date.h>
%include "sensitivity_analysis_interface.h"
