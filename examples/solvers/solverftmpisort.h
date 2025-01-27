
/* ***************************************************************************
 * This file is part of the SymPhas package, containing a framework for
 * implementing solvers for phase-field problems with compile-time symbolic 
 * algebra.
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
 * This file is part of the distributed solvers.
 * 
 * FORWARD-TIME CENTRAL-SPACE METHOD (FORWARD EULER)
 *
 * ***************************************************************************
 */

#pragma once

#include "solver.h"
#include "expressions.h"



template<size_t D>
void naive_domain_splitting(len_type const (&dimensions)[D], len_type (&boundaries)[D][2])
{
	int num_nodes = symphas::parallel::get_num_nodes();
	int rank = symphas::parallel::get_node_rank();

	int per_dimension = int(std::pow(double(num_nodes), 1. / D));
	int num_domains = std::pow(per_dimension, D);
	if (num_domains != num_nodes)
	{
		throw;
	}

	int pos[D]{};
	int dims[D]{};
	for (iter_type i = 0; i < D; ++i) dims[i] = per_dimension;

	int dims2[D](dims);

	grid::get_grid_position(pos, dims, rank);
	for (iter_type n = 0; n < D; ++n)
	{
		len_type len = dimensions[n] / per_dimension;
		boundaries[n][0] = pos[n] * len;
		boundaries[n][1] = (pos[n] + i) * len;
	}
	
}


// Get the domain that the given position is in
template<size_t D>
void naive_domain_indexing(len_type const (&dimensions)[D], len_type(&position)[D])
{
	len_type boundaries[D][2]{};
	naive_domain_splitting(dimensions, boundaries);

	int pos[D]{};
	int dims[D]{};
	for (iter_type i = 0; i < D; ++i) dims[i] = per_dimension;

	int dims2[D](dims);

	grid::get_grid_position(pos, dims, rank);
	for (iter_type n = 0; n < D; ++n)
	{
		len_type len = dimensions[n] / per_dimension;
		boundaries[n][0] = pos[n] * len;
		boundaries[n][1] = (pos[n] + i) * len;
	}

}

template<typename T, size_t D>
void reorganize_grids(SolverSystemFDwSDMPI<T, D>* systems, len_type len, len_type (&dimensions)[D])
{
	len_type boundaries[D][2]{};
	naive_domain_splitting(dimensions, boundaries);

	bool* in_boundaries_list = new bool[len] {};


	int num_nodes = symphas::parallel::get_num_nodes();
	int node0 = symphas::parallel::get_node_rank();

	for (iter_type i = 0; i < len; ++i)
	{
		grid::select_region<D> &region = systems[i].region;
		
		in_boundaries_list[i] = true;
		for (iter_type n = 0; n < D; ++n)
		{
			len_type center_n = region.origin[n] + region.dims[n] / 2;
			if (center_n < boundaries[n][0] || center_n > boundaries[n][1])
			{
				in_boundaries_list[i] = false;
			}
		}

	}

	//MPI_Request* req_recv = &req[index_count * N];
	for (iter_type i = 0; i < len; ++i)
	{
		bool owns = systems[i].thr_info.owning_node == node0;
		if (in_boundaries_list[i])
		{
			if (!owns)
			{
				//MPI_Irecv(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node, i + len * node + len * N * node0, MPI_COMM_WORLD, &req_recv[i]);
				MPI_Recv(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node0, i, MPI_COMM_WORLD);
			}
		}
		else
		{
			if (owns)
			{
				MPI_Recv(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node0, i, MPI_COMM_WORLD);
			}
		}
	}
}
