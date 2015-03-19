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


#ifndef SENSITIVITY_ANALYIS_INTERFFACE_H_
#define SENSITIVITY_ANALYIS_INTERFFACE_H_


#include "tools/date.h"
#include <src/monica.h>
#include <src/monica-parameters.h>
#include <src/debug.h>

namespace Monica 
{

	std::vector<ProductionProcess>
	applySAChanges(std::vector<ProductionProcess> ff,
					 const CentralParameterProvider &centralParameterProvider);
                                 
    std::vector<ProductionProcess> setAssimilatePartitioningCoefficient( int dev_stage, 
                    std::vector<double> partitioning_coefficient,
                    std::vector<ProductionProcess> ff,
                    const CentralParameterProvider &centralParameterProvider);

}

#endif /*SENSITIVITY_ANALYIS_H_*/

