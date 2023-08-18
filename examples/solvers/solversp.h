
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
 * SEMI-IMPLICIT FOURIER SPECTRAL SOLVER
 *
 * ***************************************************************************
 */

#include "solver.h"

#include "stencilincludes.h"
#include "spectrallib.h"
#include "gridfunctions.h"


inline len_type SP_DIMS[] = { BOUNDARY_DEPTH * 2 + 1, BOUNDARY_DEPTH * 2 + 1 };


//! Semi-implicit Fourier Spectral solver.
/*!
 * The spectral method approximates the derivatives in Fourier space by
 * solving the resulting ODE with nonhomogeneous term and then approximating
 * the operators.
 * 
 * Nonlinear terms are evaluated in real space and then the transform is
 * taken to evaluate the full equation.
 * 
 * This solver does not support a changing `dt` (time step) value.
 */
NEW_SOLVER(SolverSP)

	double h[3];

	SolverSP(const double *h = nullptr, double dt = 1, size_t dim = 0) : parent_type(dt), h{ 0 }
	{
		std::copy(h, h + dim, this->h);
	}


	/* ptr is the pointer to the specific index in the grid
	 * n is the numerical index of the given pointer
	 */
	template<typename T, size_t D>
	decltype(auto) subgrid(Grid<T, D> const& src, iter_type n) const
	{
		len_type extent[D]{ 0 };
		std::fill(extent, extent + D, BOUNDARY_DEPTH);
		return grid::get_subgrid(src, extent, n);
	}


	template<typename S>
	void step(S&&) const {}

	template<size_t Z, typename S, typename T>
	void equation(std::pair<Variable<Z, symphas::ref<S>>, T>& r) const
	{
		auto& [sys, data] = r;
		data.update();
		expr::result(data.scheme, sys.get().dframe, sys.get().transformed_len);
	}

	/*
	 * the spectral solver does some work on the given equation in order
	 * to return two separate equations; one containing linear terms
	 * and the other containing nonlinear terms
	 */
	template<size_t En, typename... Ss, size_t Z, typename S, typename E, 
		typename T_src = typename grid::value_type_of<S>::type, size_t D = grid::dimension_of<S>::value>
	decltype(auto) form_expr_one(std::tuple<Ss...> const& systems, std::pair<Variable<Z, symphas::ref<S>>, E> const& e) const
	{
		auto&& [sys, equation] = e;
		auto&& [linear, nonlinear] = expr::split::by_linear(expr::apply_operators(equation));
		//nonlinear.print(stdout);

		auto&& [linear_in_Z, linear_in_nonZ] = expr::split::separate_var<Z>(linear);
		auto&& [l_op, non_op] = solver_sp::get_l_op<Z>(linear_in_Z, h);

		auto&& A_expression = solver_sp::form_A_op<D>(l_op, dt, sys.get().dims);
		auto&& B_expression = solver_sp::form_B_op<D>(l_op, dt, sys.get().dims);

		auto&& nonlinear_scheme = solver_sp::construct_nonlinear<Z, D>(systems, B_expression, linear_in_nonZ + non_op + nonlinear, h, sys.get().dims);


		auto&& name = expr::get_fourier_name(expr::get_op_name(expr::get_variable<Z>(equation)));

		len_type dims[D];
		expr::fill_data_dimensions(equation, dims);
		auto term_ft = expr::as_grid_data<D>(sys.get().frame_t, dims);

		auto&& A_term = solver_sp::get_A_term<T_src>(A_expression);
		auto&& scheme = expr::make_add(A_term * expr::make_term<Z + sizeof...(Ss)>(NamedData(std::move(term_ft), name)), nonlinear_scheme);
		expr::printe(scheme, "spectral scheme");



		//scheme.print(stdout);
		//fprintf(stdout, "\n");

		auto data = SpectralData(scheme);

		return std::make_pair(sys, data);
	}

	template<size_t En, typename S, typename E,
		typename T_src = typename grid::value_type_of<S>::type, size_t D = grid::dimension_of<S>::value>
	decltype(auto) form_expr_one(std::pair<S*, len_type> const& systems, std::pair<DynamicVariable<S>, E> const& e) const
	{
		//auto&& [sys, equation] = e;
		//auto&& [linear, nonlinear] = expr::split::by_linear(expr::apply_operators(equation));
		////nonlinear.print(stdout);

		//auto&& [linear_in_Z, linear_in_nonZ] = expr::split::separate_var<Z>(linear);
		//auto&& [l_op, non_op] = solver_sp::get_l_op<Z>(linear_in_Z, h);

		//auto&& A_expression = solver_sp::form_A_op<D>(l_op, dt, sys.get().dims);
		//auto&& B_expression = solver_sp::form_B_op<D>(l_op, dt, sys.get().dims);

		//auto&& nonlinear_scheme = solver_sp::construct_nonlinear<Z, D>(systems, B_expression, linear_in_nonZ + non_op + nonlinear, h, sys.get().dims);


		//auto&& name = expr::get_fourier_name(expr::get_op_name(expr::get_variable<Z>(equation)));

		//len_type dims[D];
		//expr::fill_data_dimensions(equation, dims);
		//auto term_ft = expr::as_grid_data<D>(sys.get().frame_t, dims);

		//auto&& A_term = solver_sp::get_A_term<T_src>(A_expression);
		//auto&& scheme = expr::make_add(A_term * expr::make_term<Z + sizeof...(Ss)>(NamedData(std::move(term_ft), name)), nonlinear_scheme);
		//expr::printe(scheme, "spectral scheme");



		////scheme.print(stdout);
		////fprintf(stdout, "\n");

		//auto data = SpectralData(scheme);

		//return std::make_pair(sys, data);
		return std::make_pair(0, 0);
	}



	static auto make_solver(symphas::problem_parameters_type const& parameters)
	{
		size_t dim = parameters.get_dimension();
		double* h = new double[dim];

		for (iter_type i = 0; i < dim; ++i)
		{
			h[i] = parameters.get_interval_data()[0].at(symphas::index_to_axis(i)).width();
		}

		auto s = SolverSP{ h, parameters.get_time_step(), dim};
		delete[] h;
		return s;

	}



	//! Evaluate the result of the equation into the given data.
	/*!
	 * Evaluate the equation into the given data.
	 *
	 * \param r The pair consisting of the data and equation.
	 */
	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>& r) const
	{
		if constexpr (!expr::is_nonlinear<E>::value)
		{
			static bool warned = false;
			if (!warned)
			{
				fprintf(SYMPHAS_WARN, "A provisional equation which is used contains linear terms. "
					"Linear terms may cause incorrect construction of the spectral operators!\n");
				warned = true;
			}
		}

		parent_type::evaluate_one(r);
	}

	//! Evaluate the result of the equation into the given data.
	/*!
	 * Evaluate the equation into the given data.
	 *
	 * \param r The pair consisting of the data and equation.
	 */
	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>&& r) const
	{
		if constexpr (!expr::is_nonlinear<E>::value)
		{
			static bool warned = false;
			if (!warned)
			{
				fprintf(SYMPHAS_WARN, "A provisional equation which is used contains linear terms. "
					"Linear terms may cause incorrect construction of the spectral operators!\n");
				warned = true;
			}
		}
		parent_type::evaluate_one(r);

	}


};



ASSOCIATE_SOLVER_SYSTEM_TYPE(SolverSP, SolverSystemSpectral)
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverSP, ProvisionalSystemSpectral)
SYMPHAS_SOLVER_ALL_SUPPORTED(SolverSP)



