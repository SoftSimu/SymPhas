
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



//! Finite difference solver.
/*!
 * The boundaries are updated after the provisional variables are computed
 * from the corresponding equations. The boundary data is typically the same
 * as the first field that is given, but in practise it does not matter
 * for the numerical results unless finite difference approximations are
 * applied extensively to the provisional variables. This is not recommended.
 */
NEW_SOLVER_WITH_STENCIL(SolverFT)



	/*
	 * forward euler method for actually updating the grid
	 */
	template<typename S>
	void step(S &sys) const
	{
		expr::result(expr::make_term(sys.as_grid()) + expr::make_term(dt, sys.dframe), sys.as_grid(), expr::iterable_domain(sys.as_grid()));
	}


	/*
	 * the parameter is a tuple with the first parameter being a Variable object with ref
	 * base type, and the ref is a reference to the system; like this, the variable index is
	 * packaged and the equality operator is deferred to the oplvariable
	 */
	template<typename S, typename E>
	inline void equation(std::pair<S, E>& r) const
	{
		TIME_THIS_CONTEXT_LIFETIME(solverft_equation);
		expr::prune::update<expr::not_<expr::matches_series>>(r.second);
		expr::result_by_term<expr::matches_series>(r.second, r.first.get().dframe);
	}

	/*
	 * some solvers require the equation that is provided
	 * to be in a slightly different form
	 * the forward euler solver has no such requirement and so
	 * it returns its argument directly
	 */

	template<size_t En, typename SS, typename S, typename E>
	auto form_expr_one(SS&&, std::pair<S, E> const& e) const
	{
		auto [sys, equation] = e;
		auto eq_ft = expr::transform::optimize(expr::apply_operators(equation));
		expr::prune::update(eq_ft);
		expr::printe(eq_ft, "scheme");
		return std::make_pair(sys, eq_ft);
	}


	static auto make_solver(symphas::problem_parameters_type const& parameters)
	{
		if (parameters.length())
		{
			double h = parameters.get_interval_data()[0].at(Axis::X).width();
			size_t dim = parameters.get_dimension();

			// integrity check: all grid widths must be the same
			for (iter_type i = 0; i < parameters.length(); ++i)
			{
				for (iter_type n = 0; n < dim; ++n)
				{
					if (h != parameters.get_interval_data()[i].at(symphas::index_to_axis(n)).width())
					{
						char axis =
							(symphas::index_to_axis(n) == Axis::X) ? 'x' :
							(symphas::index_to_axis(n) == Axis::Y) ? 'y' :
							(symphas::index_to_axis(n) == Axis::Z) ? 'z' : '?';
						fprintf(SYMPHAS_WARN, "the grid spacing of system %d for axis '%c' is "
							"not consistent, results will not reflect the given problem!\n",
							i, axis);
					}
				}
			}

			/* the dimensions of the problem are taken from the first system
			 * since the solver assumes they are homogeneous
			 */
			len_type* dims = new len_type[dim];
			for (iter_type i = 0; i < dim; ++i)
			{
				Axis side = symphas::index_to_axis(i);
				dims[i] = parameters.get_interval_data()[0].at(side).get_count() + 2 * BOUNDARY_DEPTH;
			}

			auto s = this_type{ dims, h };
			delete[] dims;
			return s;
		}
		else
		{
			return this_type{ grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3) };
		}

	}





	// given the grid/equation pair, evaluate the equation into the grid
	// element by element
	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>& r) const
	{
		auto& [grid, equation] = r;
		if constexpr (expr::has_state<E>::value)
		{
			expr::prune::update(equation);
		}

		expr::result(equation, expr::BaseData<G>::get(grid));
	}

	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>&& r) const
	{
		auto& [grid, equation] = r;
		if constexpr (expr::has_state<E>::value)
		{
			expr::prune::update(equation);
		}

		expr::result(equation, expr::BaseData<G>::get(grid));
	}

};


//ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFD)
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDwSD)
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDwSDMPI)
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverFT, ProvisionalSystemFD)
SYMPHAS_SOLVER_ALL_SUPPORTED(SolverFT)


template<size_t D>
bool check_overlapping_domain(grid::region_interval<D> const& region0, const iter_type(&dims)[D], const iter_type(&origin)[D])
{
	grid::region_interval<D> region1(origin, dims);
	return grid::is_overlapping(region0, region1);
}

template<size_t D>
bool check_overlapping_domain(grid::region_interval_multiple<D> region, const iter_type (&dims)[D], const iter_type (&origin)[D])
{
	grid::region_interval<D> region0(origin, dims);
	region /= region0;
	if (region.regions.empty())
	{
		return false;
	}
	else
	{
		return true;
	}
}

template<typename T, size_t D>
void synchronize_regional_pos(SolverSystemFDwSDMPI<T, D>* systems, len_type len, int* region_info = nullptr, int* len_info = nullptr)
{
	auto [start, end] = symphas::parallel::get_index_range(len);
	len_type index_count = end - start;

	int node0 = symphas::parallel::get_node_rank();

	int N = symphas::parallel::get_num_nodes();
	int region_info_entries = ((len + N - 1) / N) * N;
	int info_count = D * 2;
	
	bool clear_region_info = false;
	if (!region_info)
	{
		region_info = new int[region_info_entries * info_count] {};
		clear_region_info = true;
	}

	int offset = (region_info_entries / N) * symphas::parallel::get_node_rank();
	for (iter_type i = start; i < end; ++i)
	{
		iter_type ind0 = offset + (i - start);
		iter_type ind1 = ind0 * info_count;
		for (iter_type n = 0; n < D; ++n)
		{
			region_info[ind1 + n] = systems[i].region.dims[n];
			region_info[ind1 + D + n] = systems[i].region.origin[n];
		}

		if (len_info)
		{
			len_info[ind0] = systems[i].region.len;
		}
	}

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
		&region_info[0], (region_info_entries / N) * info_count, MPI_INT, MPI_COMM_WORLD);

	if (len_info)
	{
		MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
			&len_info[0], len / N, MPI_INT, MPI_COMM_WORLD);
	}

	for (iter_type node = 0; node < N; ++node)
	{
		if (node != node0)
		{
			auto [start0, end0] = symphas::parallel::get_index_range(len, node);

			int offset = (region_info_entries / N) * node;
			for (iter_type i = start0; i < end0; ++i)
			{
				iter_type ind0 = (offset + (i - start0)) * info_count;

				iter_type dims[D]{};
				iter_type origin[D]{};
				for (iter_type n = 0; n < D; ++n)
				{
					dims[n] = region_info[ind0 + n];
					origin[n] = region_info[ind0 + D + n];
				}
				systems[i].region.update(origin, dims);
			}
		}
	}

	if (clear_region_info)
	{
		delete[] region_info;
	}
}

template<typename V, typename E, typename T, size_t D>
void synchronize_regional_grids(std::pair<V, E>* es, SolverSystemFDwSDMPI<T, D>* systems, len_type len, int *working_data)
{
	auto [start, end] = symphas::parallel::get_index_range(len);
	len_type index_count = end - start;

	int node0 = symphas::parallel::get_node_rank();

	int N = symphas::parallel::get_num_nodes();
	int region_info_entries = ((len + N - 1) / N) * N;
	int info_count = D * 2;

	int* region_info = &working_data[0];
	int* request_info = &working_data[region_info_entries * info_count];
	int* request_info_buffer = &working_data[region_info_entries * info_count + len];
	int* len_info = &working_data[region_info_entries * info_count + 2 * len];

	synchronize_regional_pos(systems, len, region_info, len_info);

	for (iter_type node = 0; node < N; ++node)
	{
		if (node != node0)
		{
			auto [start0, end0] = symphas::parallel::get_index_range(len, node);
			for (iter_type i = start0; i < end0; ++i)
			{
				for (iter_type j = start; j < end; ++j)
				{
					const auto& eq = es[j].second;
					auto region = expr::iterable_domain(eq);
 					if (check_overlapping_domain(region, systems[i].region.dims, systems[i].region.origin))
					{
						request_info_buffer[i] = (1ul << (int)(node0));
					}
				}
			}
		}
	}

	MPI_Win window;
	MPI_Win_create(&request_info[0], sizeof(int) * len, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);

	MPI_Win_fence(0, window);
	for (iter_type node = 0; node < N; ++node)
	{
		MPI_Accumulate(&request_info_buffer[0], len, MPI_INT, node, 0, len, MPI_INT, MPI_SUM, window);
		MPI_Win_fence(0, window);
	}

	len_type req_count = index_count * N + len;
	MPI_Request* req = new MPI_Request[req_count];
	MPI_Request** req_send = new MPI_Request*[N];
	MPI_Request* req_recv = &req[index_count * N];

	for (iter_type i = 0; i < N; ++i)
	{
		req_send[i] = &req[index_count * i];
	}

	for (iter_type i = 0; i < req_count; ++i)
	{
		req[i] = MPI_REQUEST_NULL;
	}

	// send to nodes
	for (iter_type i = start; i < end; ++i)
	{
		for (iter_type node = 0; node < N; ++node)
		{
			if ((request_info[i] & (1ul << node)) != 0)
			{
				MPI_Isend(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node, i + len * node0 + len * N * node, MPI_COMM_WORLD, &req_send[node][i - start]);
			}
		}
	}

	for (iter_type node = 0; node < symphas::parallel::get_num_nodes(); ++node)
	{
		// receive from node
		if (node != node0)
		{
			auto [start0, end0] = symphas::parallel::get_index_range(len, node);

			for (iter_type i = start0; i < end0; ++i)
			{
				// check if the current system needs to be received from the iteration node
				if ((request_info[i] & (1ul << node0)) != 0)
				{
					// allocate space
					if (systems[i].values == nullptr)
					{
						systems[i].values = new T[systems[i].region.len];
					}
					else
					{
						iter_type ind0 = (region_info_entries / N) * node + (i - start0);
						if (len_info[ind0] != systems[i].region.len)
						{
							delete[] systems[i].values;
							systems[i].values = new T[systems[i].region.len];
						}
					}
					// receive data
					MPI_Irecv(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node, i + len * node + len * N * node0, MPI_COMM_WORLD, &req_recv[i]);
				}
				
			}
		}
	}

	for (iter_type i = 0; i < req_count; ++i)
	{
		MPI_Wait(&req[i], MPI_STATUS_IGNORE);
	}

	delete[] req_send;
	delete[] req;
}


template<typename T, size_t D>
void regional_grids_to_host(SolverSystemFDwSDMPI<T, D>* systems, len_type len)
{
	synchronize_regional_pos(systems, len);

	MPI_Request* rqst = new MPI_Request[len];
	for (iter_type i = 0; i < len; ++i)
	{
		rqst[i] = MPI_REQUEST_NULL;
	}

	int node0 = symphas::parallel::get_node_rank();
	if (symphas::parallel::is_host_node(node0))
	{
		for (iter_type node = 0; node < symphas::parallel::get_num_nodes(); ++node)
		{
			// receive from nodes
			if (node != node0)
			{
				auto [start0, end0] = symphas::parallel::get_index_range(len, node);
				for (iter_type i = start0; i < end0; ++i)
				{
					delete[] systems[i].values;
					systems[i].values = new T[systems[i].region.len];
					MPI_Irecv(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, node, i, MPI_COMM_WORLD, &rqst[i]);
				}
			}
		}
	}
	// send to host node
	else
	{
		auto [start, end] = symphas::parallel::get_index_range(len, node0);
		for (iter_type i = start; i < end; ++i)
		{
			MPI_Isend(&(systems[i].values[0]), systems[i].region.len, MPI_DOUBLE, SYMPHAS_MPI_HOST_RANK, i, MPI_COMM_WORLD, &rqst[i]);
		}
	}

	for (iter_type i = 0; i < len; ++i)
	{
		if (rqst[i] != MPI_REQUEST_NULL)
		{
			MPI_Wait(&rqst[i], MPI_STATUS_IGNORE);
		}
	}

	delete[] rqst;

}


template<size_t D, typename Sp, typename... S>
void synchronize(void* es, ArrayModel<D, Sp, S...>* model)
{
	regional_grids_to_host(model->systems(), model->len);
}

template<size_t D>
void new_synchronize_working_data(int* &working_data, len_type len)
{
	int N = symphas::parallel::get_num_nodes();
	int region_info_entries = ((len + N - 1) / N) * N;
	int info_count = D * 2;

	len_type region_info_len = region_info_entries * info_count;
	len_type request_info_len = len;
	len_type len_info_len = region_info_entries * info_count;

	len_type data_len = region_info_len + request_info_len * 2 + len_info_len;
	if (working_data == nullptr)
	{
		working_data = new int[data_len] {};
	}
	else
	{
		std::fill(working_data, working_data + data_len, 0);
	}

}

template<typename V, typename E, size_t D, typename Sp, typename... S>
void synchronize(std::pair<V, E>* es, ArrayModel<D, Sp, S...>* model)
{
	static int* working_data = nullptr;
	new_synchronize_working_data<D>(working_data, model->len);
	synchronize_regional_grids(es, model->systems(), model->len, working_data);
}




template<typename T, size_t D>
template<typename... Ts>
void PhaseFieldSystem<RegionalGridMPI, T, D>::synchronize(Ts&&... args)
{
	::synchronize(std::forward<Ts>(args)...);
}




