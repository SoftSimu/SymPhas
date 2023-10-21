
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
 * MODULE:  lib
 * PURPOSE: Defines functions for MPI functionality.
 *
 * ***************************************************************************
 */

#pragma once

#ifdef USING_MPI

#define SYMPHAS_MPI_HOST_RANK 0

#include "definitions.h"
#include <mpi.h>

namespace symphas
{
	namespace mpi {}
}


namespace symphas::mpi
{

	inline int get_current_rank()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}

    inline bool is_host_rank()
    {
        return get_current_rank() == SYMPHAS_MPI_HOST_RANK;
    }

	inline bool is_host_rank(int rank)
	{
		return rank == SYMPHAS_MPI_HOST_RANK;
	}

	inline int get_num_nodes()
	{
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		return size;
	}

	inline std::pair<int, int> get_lower_upper(int num_fields, int rank)
	{
		int N = get_num_nodes();
		int size = num_fields / N;
		int R = num_fields % N;

		int lower = size * rank + std::min(R, rank);
		int upper = size * (rank + 1) + std::min(R, rank + 1);

		return { lower, upper };
	}


	inline std::pair<int, int> get_lower_upper(int num_fields)
	{
		int rank = get_current_rank();
		return get_lower_upper(num_fields, rank);
	}

	//! Node manages the index when 
	inline bool index_in_node(int index, int num_fields)
	{
		auto [lower, upper] = get_lower_upper(num_fields);
		return (index >= lower && index < upper);
	}

	struct info_type
	{
		int rank;
		int index;
		int num_fields;

		info_type(int index, int num_fields) : rank{ get_current_rank() }, index{ index }, num_fields{ num_fields }
		{

		}

		info_type() : rank{ 0 }, index{ 0 }, num_fields{ 1 } {}

		//! Node manages the index when 
		inline bool index_in_node() const
		{
			auto [lower, upper] = get_lower_upper(num_fields);
			return (index >= lower && index < upper);
		}
	};
}

#endif

namespace symphas
{


#ifdef USING_MPI
	using multi_thr_info_type = symphas::mpi::info_type;
#else
	using multi_thr_info_type = size_t;
#endif
}