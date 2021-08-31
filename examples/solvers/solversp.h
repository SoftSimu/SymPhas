
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


inline len_type SP_DIMS[] = { THICKNESS * 2 + 1, THICKNESS * 2 + 1 };


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

	double dt;
	double h[3];

	SolverSP(const double *h, double dt, size_t dim) : St(SP_DIMS, h[0]), dt{ dt }, h{ 0 } 
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
		std::fill(extent, extent + D, THICKNESS);
		return grid::get_subgrid(src, extent, n);
	}

	template<size_t O, typename T, size_t D>
	auto applied_generalized_derivative(Grid<T, D> const& e, iter_type n) const
	{
		auto g = subgrid(e, n);
		return St::template apply<O>(g.values + g.len / 2);
	}

	template<typename T, size_t D>
	auto applied_laplacian(Grid<T, D> const& e, iter_type n) const
	{
		auto g = subgrid(e, n);
		return St::laplacian(g.values + g.len / 2);
	}

	template<typename T, size_t D>
	auto applied_bilaplacian(Grid<T, D> const& e, iter_type n) const
	{
		auto g = subgrid(e, n);
		return St::bilaplacian(g.values + g.len / 2);
	}

	template<typename T, size_t D>
	auto applied_gradlaplacian(Grid<T, D> const& e, iter_type n) const
	{
		auto g = subgrid(e, n);
		return St::gradlaplacian(g.values + g.len / 2);
	}

	template<typename T, size_t D>
	auto applied_gradient(Grid<T, D> const& e, iter_type n) const
	{
		auto g = subgrid(e, n);
		return St::gradient(g.values + g.len / 2);
	}


	template<size_t O, typename T>
	auto applied_generalized_derivative(Block<T> const& e, iter_type n) const
	{
		auto& v = e.values[n];
		return (*static_cast<St const*>(this)).template apply<O>(&v);
	}

	template<typename T>
	auto applied_laplacian(Block<T> const& e, iter_type n) const
	{
		auto& v = e.values[n];
		return (*static_cast<St const*>(this)).laplacian(&v);
	}

	template<typename T>
	auto applied_bilaplacian(Block<T> const& e, iter_type n) const
	{
		auto& v = e.values[n];
		return (*static_cast<St const*>(this)).bilaplacian(&v);
	}

	template<typename T>
	auto applied_gradlaplacian(Block<T> const& e, iter_type n) const
	{
		auto& v = e.values[n];
		return (*static_cast<St const*>(this)).gradlaplacian(&v);
	}

	template<typename T>
	auto applied_gradient(Block<T> const& e, iter_type n) const
	{
		auto& v = e.values[n];
		return (*static_cast<St const*>(this)).gradient(&v);
	}




	template<typename S>
	void step(S&&, double) {}

	template<size_t Z, typename S, typename T>
	void equation(std::pair<Variable<Z, symphas::ref<S>>, T>& r)
	{
		auto& [sys, data] = r;
		data.update();
		expr::result(data.evolution_equation, sys.get().dframe, sys.get().transformed_len);
	}


	/*
	 * the spectral solver does some work on the given equation in order
	 * to return two separate equations; one containing linear terms
	 * and the other containing nonlinear terms
	 */
	template<size_t En, typename... Ss, size_t Z, typename S, typename E, typename T_src = typename grid::value_type_of<S>::type, size_t D = grid::dimension_of<S>::value>
	decltype(auto) form_expr_one(std::tuple<Ss...> const& systems, std::pair<Variable<Z, symphas::ref<S>>, E>&& e)
	{
		auto&& [sys, equation] = e;
		auto&& [linear, nonlinear] = expr::split::by_linear(equation);

		auto [linear_in_Z, nonlinear_in_Z] = expr::split::separate_var<Z>(linear);
		auto [l_op, non_op] = solver_sp::get_l_op<Z>(linear_in_Z, h);


		auto linear_swapped = solver_sp::swap_var(systems, nonlinear_in_Z + non_op);

		auto&& nl_split_Z = expr::split::separate_deriv(linear_swapped);
		auto&& nl_split_R = expr::split::separate_deriv(nonlinear);
		auto&& [differential_ops, nonlinear_pair_list] = 
			symphas::lib::unfurl_tuple(
				symphas::lib::unzip_pot<OpVoid>(nl_split_Z, nl_split_R));



		auto&& A_expression = solver_sp::form_A_op<D>(l_op, dt, sys.get().dims);
		auto&& B_expressions = solver_sp::form_B_ops<D>(differential_ops, l_op, dt, h, sys.get().dims);

		// the A operator for the linear term
		auto&& A = expr::transform::to_fftw_grid(A_expression);
		
		// a tuple of B operators associated with all the nonlinear terms and their respective derivatives
		auto&& Bs = solver_sp::pair_B_with_working<T_src, D>(B_expressions, sys.get().dims);

		// a list of equations of the nonlinear terms, their list index corresponding to the B operators
		auto&& nls = solver_sp::get_nonlinear_exprs<D>(nonlinear_pair_list, h, sys.get().dims);

		// the nonlinear equation which is solved, append it to the tuple of equation data
		auto&& data = SpectralData(
			expr::make_op(NamedData(A, A_expression)),
			Bs, 
			nls, 
			sys.get().frame_t, 
			expr::get_op_name(expr::property::get_data_variable<Z>(equation))
		);

		expr::printf(A_expression, "A operator");
		expr::printf(data.evolution_equation, "evolution equation");

		return std::make_pair(sys, data);
	}




	static auto make_solver(symphas::problem_parameters_type const& parameters)
	{
		size_t dim = parameters.get_dimension();
		double* h = new double[dim];

		for (iter_type i = 0; i < dim; ++i)
		{
			h[i] = parameters.get_interval_data()[0].at(symphas::index_to_axis(i)).width();
		}

		auto s = SolverSP<St>{ h, parameters.dt, dim };
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
	void evaluate_one(std::pair<G, E>& r)
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
	void evaluate_one(std::pair<G, E>&& r)
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
SYMPHAS_SOLVER_SUPPORTED_TYPE(SolverSP, scalar_t)
SYMPHAS_SOLVER_SUPPORTED_TYPE(SolverSP, complex_t)



