
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 * 
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 * 
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 * 
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 */


#include "boundarysystem.h"

DLLSOL double params::regional_resize_factor = 1.5;
DLLSOL double params::regional_resize_time = 2.5;
DLLSOL double params::regional_resize_cutoff_eps = 1e-6;


bool add_solution_params(param_map_type& param_map)
{
	using namespace params;

	param_map["regional-resize-factor"] = { &regional_resize_factor, new param_assign<double>, 'z' };
	param_map["regional-factor"] = { &regional_resize_factor, new param_assign<double>, 'z' };
	param_map["regional-resize-time"] = { &regional_resize_time, new param_assign<double>, 'Z' };
	param_map["regional-delta"] = { &regional_resize_time, new param_assign<double>, 'Z' };
	param_map["regional-resize-eps"] = { &regional_resize_cutoff_eps, new param_assign<double>, 'E' };
	param_map["regional-eps"] = { &regional_resize_cutoff_eps, new param_assign<double>, 'E' };

	return true;
}

