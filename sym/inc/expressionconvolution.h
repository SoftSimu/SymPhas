
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
 * MODULE:  expr
 * PURPOSE: Implements the convolution operation. Represents taking the
 * convolution between two expressions, with specializations based on
 * when one of the terms is a variable or smoothing kernel.
 *
 * ***************************************************************************
 */

#pragma once


#include "expressionaggregates.h"
#include "convolutionlib.h"

//! \cond

#ifdef LATEX_PLOT
#define SYEX_CONVOLUTION_FMT_A "\\left("
#define SYEX_CONVOLUTION_FMT_B "\\right)"
#define SYEX_CONVOLUTION_FMT_SEP "*"
#else
#define SYEX_CONVOLUTION_FMT_A "("
#define SYEX_CONVOLUTION_FMT_B ")"
#define SYEX_CONVOLUTION_FMT_SEP " x "
#endif

#define SYEX_CONVOLUTION_FMT SYEX_CONVOLUTION_FMT_A "%s" SYEX_CONVOLUTION_FMT_SEP "%s" SYEX_CONVOLUTION_FMT_B
#define SYEX_CONVOLUTION_FMT_LEN (STR_ARR_LEN(SYEX_CONVOLUTION_FMT_A SYEX_CONVOLUTION_FMT_B SYEX_CONVOLUTION_FMT_SEP) - 1)

//! \endcond

namespace symphas::internal
{
	//! Implementation of convolution expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct convolution expressions. Wraps the template deduction necessary
	 * to initialize a convolution expression.
	 */
	struct make_convolution
	{
		//! Constructs the convolution with the identity coefficient.
		template<typename A, typename B>
		static auto get(A&&, B&&);

		//! Constructs the convolution applied to the given expression.
		template<typename V, typename E1, typename E2>
		static auto get(V, OpExpression<E1> const&, OpExpression<E2> const&);



		//! Constructs the convolution applied to the given expression.
		template<typename V, size_t D, typename E>
		static auto get(V, GaussianSmoothing<D> const&, OpExpression<E> const&);

		//! Constructs the convolution applied to the given expression.
		template<typename V, size_t D, typename E>
		static auto get(V, OpExpression<E> const&, GaussianSmoothing<D> const&);



		//! Constructs the convolution between a variable and Gaussian kernel.
		template<typename V, size_t D, typename S, typename G>
		static auto get(V, OpLVariable<S, G> const&, GaussianSmoothing<D> const&);

		//! Constructs the convolution between a variable and Gaussian kernel.
		template<typename V, size_t D, typename S, typename G>
		static auto get(V, GaussianSmoothing<D> const&, OpLVariable<S, G> const&);



		//! Constructs the convolution using a grid instead of an expression.
		/*!
		 * Used for the convolution specialization for the OpLVariable.
		 */
		template<typename V, size_t D, typename G>
		static auto get_g(V v, GaussianSmoothing<D> const&, G g);

	};
}


namespace expr
{
	//! Create a convolution expression.
	/*!
	 * Create a convolution expression.
	 * 
	 * \param a The left hand side of the convolution operation.
	 * \param b The right hand side of the convolution operation.
	 */
	template<typename A, typename B>
	auto make_convolution(A&& a, B&& b)
	{
		return symphas::internal::make_convolution::template get(std::forward<A>(a), std::forward<B>(b));
	}

	//! Create a convolution expression.
	/*!
	 * Create a convolution expression.
	 *
	 * \param v The coefficient of the convolution expression.
	 * \param a The left hand side of the convolution operation.
	 * \param b The right hand side of the convolution operation.
	 */
	template<typename V, typename A, typename B>
	auto make_convolution(V&& v, A&& a, B&& b)
	{
		return symphas::internal::make_convolution::template get(std::forward<V>(v), std::forward<A>(a), std::forward<B>(b));
	}
}



//! Convolution of two arbitrary expressions.
/*!
 * An implementation of the convolution operator; the convolution takes two
 * functions and upon evaluating, computes the convolution.
 * In order to evaluate the convolution, the update function must be called in 
 * the pruning step.
 *
 * In order to compute the convolution, the convolution operator uses the
 * convolution theorem, evaluating the terms in Fourier space before performing
 * an inverse Fourier transform back into real space.
 * 
 * \tparam V Type of the coefficient.
 * \tparam E1 Left hand expression type.
 * \tparam E2 Right hand expression type.
 */
template<typename V, typename E1, typename E2>
struct OpFuncConvolution : OpExpression<OpFuncConvolution<V, E1, E2>>
{
	static const int e1_U = expr::grid_dim<E1>::dimension;
	static const int e2_U = expr::grid_dim<E2>::dimension;
	static_assert(e1_U == e2_U);

	static const int D = e1_U;

	// identifying the type that comes out as a result of the derivative
	// get the type of the system holding the intermediate results; based on 
	// the multiplication of the 'real' space terms
	using e1_T = typename expr::eval_type<E1>::type;
	using e2_T = typename expr::eval_type<E2>::type;
	using G_T = mul_result_t<e1_T, e2_T>;
	
	//! Generate the convolution expression.
	/*!
	 * The convolution expression is generated using the given value as the
	 * coefficient and the two terms, with the order representing the
	 * side they appear on relative to the convolution operator.
	 * 
	 * \param value The coefficient value of the convolution term.
	 * \param a The left hand side of the convolution operator.
	 * \param b The right hand side of the convolution operator.
	 */
	OpFuncConvolution(V value, E1 const& a, E2 const& b) :
		a{ a }, b{ b }, value{ value }, g0{ expr::property::data_dimensions(a, b) }, 
		data_a{ expr::property::data_len(a, b) }, 
		data_b{ expr::property::data_len(a, b) },
		compute{ data_a.values, data_b.values, g0 } { /*update();*/ }

	OpFuncConvolution(OpFuncConvolution<V, E1, E2> const& other) :
		a{ other.a }, b{ other.b }, value{ other.value }, g0{ expr::property::data_dimensions(a, b) }, 
		data_a{ expr::property::data_len(a, b) }, 
		data_b{ expr::property::data_len(a, b) },
		compute{ data_a.values, data_b.values, g0 } { /*update();*/; }

	OpFuncConvolution(OpFuncConvolution<V, E1, E2>&& other) noexcept :
		a{ other.a }, b{ other.b }, value{ other.value }, g0{ std::move(other.g0) }, 
		data_a{ std::move(other.data_a) }, 
		data_b{ std::move(other.data_b) },
		compute{ data_a.values, data_b.values, g0 } { /*update();*/ }


	inline auto eval(iter_type n) const
	{
		return value * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get(-value, a, b);
	}


	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
		n += a.print(out);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
		n += b.print(out);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
		n += b.print(out + n);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN
			+ a.print_length() + b.print_length();
	}

	E1 a;		//!< First expression in the convolution.
	E2 b;		//!< Second expression in the convolution.
	V value;	//!< Coefficient of this convolution expression.


	friend struct expr::compound_get;


	//! Update the convolution by computing the result into the stored grid. 
	/*!
	 * Update the convolution by computing the result into the stored grid.
	 */
	void update()
	{
		// compute the result of the expressions and update the grids
		expr::result(a, data_a.values, data_a.len);
		expr::result(b, data_b.values, data_b.len);

		symphas::dft::fftw_execute(compute.p_in_out_0);
		symphas::dft::fftw_execute(compute.p_in_out_1);

		for (iter_type i = 0; i < g0.len; ++i)
		{
			compute.in_2[i][0] = compute.out_0[i][0] * compute.out_1[i][0] - compute.out_0[i][1] * compute.out_1[i][1];
			compute.in_2[i][1] = compute.out_0[i][0] * compute.out_1[i][1] + compute.out_0[i][1] * compute.out_1[i][0];
		}

		symphas::dft::fftw_execute(compute.p_out_in);
		symphas::dft::scale(g0);
	}

protected:

	Block<e1_T> data_a;			//!< Data of the result of the first expression.
	Block<e2_T> data_b;			//!< Data of the result of the second expression.
	Grid<G_T, D> g0;			//!< Grid storing the final result of the convolution.

	expr::ConvolutionDataPair<D> compute;

};


template<typename S1, typename V2, typename E1, typename E2, typename T2>
auto operator*(OpLiteral<S1> const& a, OpFuncConvolution<V2, E1, E2> const& b)
{
	return OpFuncConvolution(a.value * b.value, b.a, b.b);
}

// ******************************************************************************************************************


//! Convolution of an expression with a Gaussian smoothing kernel.
/*!
 * Specialization of the convolution of two arbitrary expressions.
 * 
 * \tparam V Type of the coefficient.
 * \tparam D the dimension of the Gaussian smoothing kernel.
 * \tparam E The type of the expression convoluted with the Gaussian smoother.
 */
template<typename V, size_t D, typename E>
struct OpFuncConvolution<V, GaussianSmoothing<D>, E> : OpExpression<OpFuncConvolution<V, GaussianSmoothing<D>, E>>
{
	using G_T = typename expr::eval_type<E>::type;


	//! Generate the convolution expression.
	/*!
	 * The convolution expression is generated using the given value as the
	 * coefficient and the two terms, with the order representing the
	 * side they appear on relative to the convolution operator.
	 *
	 * \param value The coefficient value of the convolution term.
	 * \param smoother The Gaussian smoothing kernel.
	 * \param e The expression that is smoothed with this expression.
	 */
	OpFuncConvolution(V value, GaussianSmoothing<D> const& smoother, E const& e) :
		e{ e }, 
		value{ value }, 
		smoother{ smoother }, 
		g0{ expr::property::data_dimensions(smoother) },
		data{ expr::property::data_len(smoother) }, 
		compute{ data.values, g0 } 
	{ update(); }

	OpFuncConvolution(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& other) :
		e{ other.e }, value{ other.value }, smoother{ other.smoother }, g0{ expr::property::data_dimensions(other.smoother) },
		data{ expr::property::data_len(e) }, compute{ data.values, g0 } { /*update();*/ }
	OpFuncConvolution(OpFuncConvolution<V, GaussianSmoothing<D>, E>&& other) noexcept :
		e{ other.e }, value{ other.value }, smoother{ other.smoother }, g0{ std::move(other.g0) },
		data{ std::move(other.data) }, compute{ data.values, g0 } { /*update();*/ }


	inline auto eval(iter_type n) const
	{
		return value * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get(-value, smoother, e);
	}


	/* print statements
	 */

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
		n += smoother.print(out);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
		n += e.print(out);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
		n += smoother.print(out + n);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN
			+ e.print_length() + smoother.print_length();
	}


	E e;							//!< The expression to populate the first grid.
	V value;						//!< Value multiplying the result of this convolution.
	GaussianSmoothing<D> smoother;	//!< The smoothing kernel.


	friend struct expr::compound_get;

	//! Update the convolution by computing the result into the stored grid. 
	/*!
	 * Update the convolution by computing the result into the stored grid.
	 */
	void update()
	{
		expr::result(e, data.values, data.len);
		symphas::dft::fftw_execute(compute.p_in_out);

		auto f = [&](iter_type i, iter_type ft_i)
		{
			compute.in_1[ft_i][0] = smoother.eval(i) * compute.out_0[ft_i][0];
			compute.in_1[ft_i][1] = smoother.eval(i) * compute.out_0[ft_i][1];
		};

		symphas::dft::iterate_rc<G_T, D>(f, g0.dims);
		symphas::dft::fftw_execute(compute.p_out_in);
		symphas::dft::scale(g0);
	}


protected:

	Grid<G_T, D> g0;			//!< The grid storing the result of this convolution.
	Block<G_T> data;			//!< Grid storing the result of the second expression.

	expr::ConvolutionData<D> compute;
};



template<typename S1, typename V2, size_t U2, typename E2>
auto operator*(OpLiteral<S1> const& a, OpFuncConvolution<V2, GaussianSmoothing<U2>, E2> const& b)
{
	return symphas::internal::make_convolution::get(a.value * b.value, b.smoother, b.e);
}

// ******************************************************************************************************************

//! Convolution of an OpLVariable with a Gaussian smoothing kernel.
/*!
 * Specialization of the convolution of two arbitrary expressions.
 * 
 * \tparam V Type of the coefficient.
 * \tparam D the dimension of the Gaussian smoothing kernel.
 * \tparam G The type of variable data.
 */
template<typename V, size_t D, typename G>
struct OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>> : 
	OpExpression<OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>>
{
	using E = OpLVariable<V, G>;
	using G_T = typename expr::eval_type<E>::type;
	using G0 = Grid<G_T, D>;


	//! Generate the convolution expression.
	/*!
	 * The convolution expression is generated using the given value as the
	 * coefficient and the two terms, with the order representing the
	 * side they appear on relative to the convolution operator.
	 *
	 * \param value The coefficient value of the convolution term.
	 * \param smoother The Gaussian smoothing kernel.
	 * \param a The variable which is smoothed.
	 */
	template<typename V0, typename V1, typename std::enable_if<std::is_convertible<mul_result_t<V0, V1>, V>::value, int>::type = 0>
	OpFuncConvolution(V0 value, GaussianSmoothing<D> const& smoother, OpLVariable<V1, G> const& a) :
		data{ a.data }, value{ value * a.value }, smoother{ smoother }, g0{ expr::property::data_dimensions(smoother) },
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/; }
	OpFuncConvolution(V value, GaussianSmoothing<D> const& smoother, G grid) :
		data{ grid }, value{ value }, smoother{ smoother }, g0{ expr::property::data_dimensions(smoother) },
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/ }


	OpFuncConvolution(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>> const& other) :
		data{ other.data }, value{ other.value }, smoother{ other.smoother }, g0{ expr::property::data_dimensions(other.smoother) },
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/ }
	OpFuncConvolution(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>&& other) noexcept :
		data{ other.data }, value{ other.value }, smoother{ other.smoother }, g0{ std::move(other.g0) },
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/ }

	inline auto eval(iter_type n) const
	{
		return value * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get_g(-value, smoother, data);
	}


	/* print statements
	 */


	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
		n += smoother.print(out);
		n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
		n += fprintf(out, "%s", expr::get_op_name(data));
		n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
		return n;;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
		n += smoother.print(out + n);
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
		n += sprintf(out + n, "%s", expr::get_op_name(data));
		n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
		return n;;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN
			+ smoother.print_length() + std::strlen(expr::get_op_name(data));
	}


	V value;						//!< Value multiplying the result of this convolution.
	GaussianSmoothing<D> smoother;	//!< The smoothing kernel.

	friend struct expr::compound_get;


	void update()
	{
		symphas::dft::fftw_execute(compute.p_in_out);

		auto f = [&](iter_type ft_i, iter_type i)
		{
			compute.in_1[i][0] = smoother.eval(ft_i) * compute.out_0[i][0];
			compute.in_1[i][1] = smoother.eval(ft_i) * compute.out_0[i][1];

		};

		symphas::dft::iterate_rc<G_T, D>(f, g0.dims);
		symphas::dft::fftw_execute(compute.p_out_in);
		symphas::dft::scale(g0);
	}

protected:

	G data;							//!< The data from the OpLVariable.
	G0 g0;							//!< Grid storing the result of the convolution.
	expr::ConvolutionData<D> compute;

};


template<typename V0, size_t D, typename T, typename G>
OpFuncConvolution(V0, GaussianSmoothing<D>, OpLVariable<T, G>)->OpFuncConvolution<mul_result_t<V0, T>, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>;
template<typename V, size_t D, typename E>
OpFuncConvolution(V, GaussianSmoothing<D>, OpExpression<E>)->OpFuncConvolution<V, GaussianSmoothing<D>, E>;
template<typename V, typename E1, typename E2>
OpFuncConvolution(V, OpExpression<E1>, OpExpression<E2>)->OpFuncConvolution<V, E1, E2>;



template<typename S1, typename V2, size_t U2, typename G2>
auto operator*(OpLiteral<S1> const& a, OpFuncConvolution<V2, GaussianSmoothing<U2>, OpLVariable<OpIdentity, G2>> const& b)
{
	return symphas::internal::make_convolution::get_g(a.value * b.value, b.smoother, b.op_g);
}

// ************************************************************************************************

namespace symphas::internal
{

	template<typename A, typename B>
	auto make_convolution::get(A&& a, B&& b)
	{
		return get(OpIdentity{}, std::forward<A>(a), std::forward<B>(b));
	}

	template<typename V, typename E1, typename E2>
	auto make_convolution::get(V value, OpExpression<E1> const& e1, OpExpression<E2> const& e2)
	{
		return OpFuncConvolution<V, E1, E2>(value, *static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
	}

	template<typename V, size_t D, typename E>
	auto make_convolution::get(V value, GaussianSmoothing<D> const& smoother, OpExpression<E> const& e)
	{
		return OpFuncConvolution<V, GaussianSmoothing<D>, E>(value, smoother, *static_cast<const E*>(&e));
	}

	template<typename V, size_t D, typename E>
	auto make_convolution::get(V value, OpExpression<E> const& e, GaussianSmoothing<D> const& smoother)
	{
		return get(value, smoother, *static_cast<E const*>(&e));
	}

	template<typename V, size_t D, typename T, typename G>
	auto make_convolution::get(V value, GaussianSmoothing<D> const& smoother, OpLVariable<T, G> const& e)
	{
		return OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>(value, smoother, e);
	}

	template<typename V, size_t D, typename T, typename G>
	auto make_convolution::get(V value, OpLVariable<T, G> const& e, GaussianSmoothing<D> const& smoother)
	{
		return get(value, smoother, e);
	}


	template<typename V, size_t D, typename G>
	auto make_convolution::get_g(V value, GaussianSmoothing<D> const& smoother, G grid)
	{
		return OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>(value, smoother, grid);
	}

}




