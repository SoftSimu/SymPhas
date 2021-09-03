
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
 *
 * This file is part of the SymPhas API. 
 *
 * ***************************************************************************
 */

#pragma once

#include "prereq.h"


namespace symphas
{


	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Makes checkpoints and performs
	 * any postprocessing which is necessary (and if the module is enabled).
	 *
	 * \param model The model which is solved.
	 * \param dt The time stepping distance.
	 * \param n The number of solution intervals to perform.
	 */
	template<typename M>
	void find_solution(M& model, double dt, iter_type n)
	{
		double run_time = 0;
		model.print_info(SYMPHAS_LOG);
		{
			symphas::Time t;
			run_model(model, n, dt, TIME_INIT + model.index() * dt);
			run_time += t.current_duration();
		}
		fprintf(SYMPHAS_LOG, "completed %d iterations in %lf seconds.\n", n, run_time);
		fprintf(SYMPHAS_LOG, OUTPUT_BANNER);
	}

}


