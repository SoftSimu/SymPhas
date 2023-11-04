
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


#ifdef _MSC_VER
#include <process.h>
#else
#include <unistd.h> 
#endif


namespace symphas
{
	namespace parallel {}
}


namespace symphas::parallel
{

	inline int get_node_pid()
	{
#ifdef _MSC_VER
		return _getpid();
#else
		return getpid();
#endif
	}

	inline int get_node_rank()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}

    inline bool is_host_node()
    {
        return get_node_rank() == SYMPHAS_MPI_HOST_RANK;
    }

	inline bool is_host_node(int rank)
	{
		return rank == SYMPHAS_MPI_HOST_RANK;
	}

	inline int get_num_nodes()
	{
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		return size;
	}

	inline std::pair<iter_type, iter_type> get_index_range(len_type num_fields, iter_type rank)
	{
		int N = get_num_nodes();
		int size = num_fields / N;
		int R = num_fields % N;

		int lower = size * rank + std::min(R, rank);
		int upper = size * (rank + 1) + std::min(R, rank + 1);

		return { lower, upper };
	}


	inline std::pair<iter_type, iter_type> get_index_range(len_type num_fields)
	{
		int rank = get_node_rank();
		return get_index_range(num_fields, rank);
	}

	//! Node manages the index when 
	inline bool index_in_node(int index, int num_fields)
	{
		auto [lower, upper] = get_index_range(num_fields);
		return (index >= lower && index < upper);
	}

	struct info_type
	{
		int rank;
		int index;
		int num_fields;

		info_type(int index, int num_fields) : rank{ get_node_rank() }, index{ index }, num_fields{ num_fields }
		{

		}

		info_type() : rank{ 0 }, index{ 0 }, num_fields{ 1 } {}

		//! Node manages the index when 
		inline bool index_in_node() const
		{
			auto [lower, upper] = get_index_range(num_fields);
			return (index >= lower && index < upper);
		}
	};
}

#else

namespace symphas::parallel
{

	inline int get_node_rank()
	{
		return SYMPHAS_MPI_HOST_RANK;
	}

	inline bool is_host_node()
	{
		return true;
	}

	inline bool is_host_node(int rank)
	{
		return true;
	}

	inline int get_num_nodes()
	{
		return 1;
	}

	inline std::pair<int, int> get_index_range(int num_fields, int rank)
	{
		return { 0, num_fields };
	}


	inline std::pair<int, int> get_index_range(int num_fields)
	{
		return get_index_range(num_fields, get_node_rank());
s	}

	//! Node manages the index when 
	inline bool index_in_node(int index, int num_fields)
	{
		return true;
	}

	using info_type = size_t;
	{
		int rank;
		int index;
		int num_fields;

		info_type(int index, int num_fields) : rank{ get_node_rank() }, index{ index }, num_fields{ num_fields }
		{

		}

		info_type() : rank{ 0 }, index{ 0 }, num_fields{ 1 } {}

		//! Node manages the index when 
		inline bool index_in_node() const
		{
			auto [lower, upper] = get_index_range(num_fields);
			return (index >= lower && index < upper);
		}
	};
}

#endif

namespace symphas
{
	using multi_thr_info_type = symphas::parallel::info_type;
}