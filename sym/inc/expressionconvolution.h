
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


//#include "expressionaggregates.h"
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
		template<typename V, typename E2>
		static auto get(V const&, OpVoid, OpExpression<E2> const&)
		{
			return OpVoid{};
		}

		//! Constructs the convolution applied to the given expression.
		template<typename V, typename E1>
		static auto get(V const&, OpExpression<E1> const&, OpVoid)
		{
			return OpVoid{};
		}

		//! Constructs the convolution applied to the given expression.
		template<typename V, typename E1, typename E2>
		static auto get(V const&, OpExpression<E1> const&, OpExpression<E2> const&);

		//! Constructs the convolution applied to the given expression.
		template<typename V, size_t D, typename E>
		static auto get(V const&, GaussianSmoothing<D> const&, OpExpression<E> const&);

		//! Constructs the convolution applied to the given expression.
		template<typename V, size_t D, typename E>
		static auto get(V const& value, OpExpression<E> const& e, GaussianSmoothing<D> const& smoother)
		{
			return get(value, smoother, *static_cast<E const*>(&e));
		}

		//! Constructs the convolution between a variable and Gaussian kernel.
		template<typename V, size_t D, typename S, typename G>
		static auto get(V const&, GaussianSmoothing<D> const&, OpTerm<S, G> const&);

		//! Constructs the convolution between a variable and Gaussian kernel.
		template<typename V, size_t D, typename S, typename G>
		static auto get(V const& value, OpTerm<S, G> const& e, GaussianSmoothing<D> const& smoother)
		{
			return get(value, smoother, e);
		}



		//! Constructs the convolution using a grid instead of an expression.
		/*!
		 * Used for the convolution specialization for the OpTerm.
		 */
		template<typename V, size_t D, typename G>
		static auto get_g(V const&, GaussianSmoothing<D> const&, G g);

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
struct OpConvolution : OpExpression<OpConvolution<V, E1, E2>>
{
	static const size_t e1_U = expr::grid_dim<E1>::dimension;
	static const size_t e2_U = expr::grid_dim<E2>::dimension;

	static const size_t D = fixed_max<e1_U, e1_U>;

	// identifying the type that comes out as a result of the derivative
	// get the type of the system holding the intermediate results; based on 
	// the multiplication of the 'real' space terms
	using e1_T = typename expr::eval_type<E1>::type;
	using e2_T = typename expr::eval_type<E2>::type;
	using G_T = mul_result_t<e1_T, e2_T>;

	OpConvolution() : data_a{ 0 }, data_b{ 0 }, g0{ 0 }, value{ V{} }, a{}, b{} {}

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
	OpConvolution(V value, E1 const& a, E2 const& b) :
		data_a{ expr::data_length(a, b) },
		data_b{ expr::data_length(a, b) },
		g0{ expr::data_dimensions(a, b) },
		value{ value }, a{ a }, b{ b },
		compute{ data_a.values, data_b.values, g0 } { /*update();*/ }

	OpConvolution(OpConvolution<V, E1, E2> const& other) :
		data_a{ expr::data_length(other.a, other.b) },
		data_b{ expr::data_length(other.a, other.b) },
		g0{ expr::data_dimensions(other.a, other.b) },
		value{ other.value }, a{ other.a }, b{ other.b },
		compute{ data_a.values, data_b.values, g0 } { /*update();*/; }

	OpConvolution(OpConvolution<V, E1, E2>&& other) noexcept
		: OpConvolution() 
	{
		swap(*this, other);
	}


	friend void swap(OpConvolution<V, E1, E2>& first,
		OpConvolution<V, E1, E2>& second)
	{
		using std::swap;
		swap(first.data_a, second.data_a);
		swap(first.data_b, second.data_b);
		swap(first.g0, second.g0);
		swap(first.value, second.value);
		swap(first.a, second.a);
		swap(first.b, second.b);
		swap(first.compute, second.compute);
	}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get(-value, a, b);
	}

	//template<typename V1,
	//	typename std::enable_if_t<expr::is_combinable<E1>, int> = 0>
	//auto operator+(OpConvolution<V1, E1, E2> const& other) const
	//{
	//	return symphas::internal::make_convolution::get(value + other.value, a, b);
	//}

	//template<typename V1>
	//auto operator+(OpConvolution<V1, E2, E1> const& other) const
	//{
	//	return symphas::internal::make_convolution::get(value + other.value, a, b);
	//}

	//template<typename V1>
	//auto operator-(OpConvolution<V1, E1, E2> const& other) const
	//{
	//	return symphas::internal::make_convolution::get(value - other.value, a, b);
	//}

	//template<typename V1>
	//auto operator-(OpConvolution<V1, E2, E1> const& other) const
	//{
	//	return symphas::internal::make_convolution::get(value - other.value, a, b);
	//}

#ifdef PRINTABLE_EQUATIONS

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

#endif

protected:

	Block<e1_T> data_a;			//!< Data of the result of the first expression.
	Block<e2_T> data_b;			//!< Data of the result of the second expression.
	Grid<G_T, D> g0;			//!< Grid storing the final result of the convolution.

public:

	V value;					//!< Coefficient of this convolution expression.
	E1 a;						//!< First expression in the convolution.
	E2 b;						//!< Second expression in the convolution.


    template<typename V0, typename E10, typename E20>
	friend auto const& expr::get_result_data(OpConvolution<V0, E10, E20> const&);
    template<typename V0, typename E10, typename E20>
	friend auto& expr::get_result_data(OpConvolution<V0, E10, E20>&);


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

		len_type len = symphas::dft::length<G_T, D>(g0.dims);
		if constexpr (std::is_same<G_T, complex_t>::value)
		{
			grid::scale(compute.out_0, len);
			grid::scale(compute.out_1, len);
		}

#		pragma omp parallel for
		for (iter_type i = 0; i < len; ++i)
		{
			compute.in_2[i][0] = compute.out_0[i][0] * compute.out_1[i][0] - compute.out_0[i][1] * compute.out_1[i][1];
			compute.in_2[i][1] = compute.out_0[i][0] * compute.out_1[i][1] + compute.out_0[i][1] * compute.out_1[i][0];
		}

		symphas::dft::fftw_execute(compute.p_out_in);

		if constexpr (std::is_same<G_T, scalar_t>::value)
		{
			grid::scale(g0);
		}
	}

protected:

	expr::ConvolutionDataPair<D> compute;


};


template<typename coeff_t, typename V2, typename E1, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<V2, E1, E2> const& b)
{
	return symphas::internal::make_convolution::get(value * b.value, b.a, b.b);
}

template<typename coeff_t, typename tensor_t, typename E1, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<tensor_t, E1, E2> const& b)
{
	return (value * b.value) * symphas::internal::make_convolution::get(OpIdentity{}, b.a, b.b);
}

template<typename tensor_t, typename V, typename E1, typename E2,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V>
		&& expr::eval_type<E1>::rank == 0 && expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b)
{
	return symphas::internal::make_convolution::get(tensor * b.value, b.a, b.b);
}

template<typename tensor_t, typename V, typename E1, typename E2,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V>
		&& expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b)
{
	return symphas::internal::make_convolution::get(b.value, tensor * b.a, b.b);
}

template<typename tensor_t, typename V, typename E1, typename E2,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V>
		&& expr::eval_type<E1>::rank == 0 && expr::eval_type<E2>::rank > 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b)
{
	return symphas::internal::make_convolution::get(b.value, b.a, tensor * b.b);
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
struct OpConvolution<V, GaussianSmoothing<D>, E> : OpExpression<OpConvolution<V, GaussianSmoothing<D>, E>>
{
	using G_T = typename expr::eval_type<E>::type;

	OpConvolution() : g0{ 0 }, data{ 0 }, value{ V{} }, smoother{ GaussianSmoothing<D>() } {}

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
	OpConvolution(V value, GaussianSmoothing<D> const& smoother, E const& e) :
		g0{ expr::data_dimensions(smoother) }, data{ expr::data_length(smoother) },
		value{ value }, e{ e }, smoother{ smoother }, 
		compute{ data.values, g0 } 
	{ /*update();*/ }

	OpConvolution(OpConvolution<V, GaussianSmoothing<D>, E> const& other) :
		g0{ expr::data_dimensions(other.smoother) }, data{ expr::data_length(other.e) },
		value{ other.value }, e{ other.e }, smoother{ other.smoother },
		compute{ data.values, g0 } { /*update();*/ }

	OpConvolution(OpConvolution<V, GaussianSmoothing<D>, E>&& other) noexcept
		: OpConvolution()
	{
		swap(*this, other);
	}

	friend void swap(OpConvolution<V, GaussianSmoothing<D>, E>& first,
		OpConvolution<V, GaussianSmoothing<D>, E>& second)
	{
		using std::swap;
		swap(first.g0, second.g0);
		swap(first.data, second.data);
		swap(first.value, second.value);
		swap(first.e, second.e);
		swap(first.smoother, second.smoother);
		swap(first.compute, second.compute);
	}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get(-value, smoother, e);
	}

	template<typename V1>
	auto operator+(OpConvolution<V1, GaussianSmoothing<D>, E> const& other) const
	{
		return symphas::internal::make_convolution::get(value + other.value, smoother, e);
	}

	template<typename V1>
	auto operator-(OpConvolution<V1, GaussianSmoothing<D>, E> const& other) const
	{
		return symphas::internal::make_convolution::get(value - other.value, smoother, e);
	}

#ifdef PRINTABLE_EQUATIONS

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

#endif

protected:

	Grid<G_T, D> g0;				//!< The grid storing the result of this convolution.
	Block<G_T> data;				//!< Grid storing the result of the second expression.

public:

	V value;						//!< Value multiplying the result of this convolution.
	E e;							//!< The expression to populate the first grid.
	GaussianSmoothing<D> smoother;	//!< The smoothing kernel.

    template<typename V0, size_t D0, typename E0>
	friend auto const& expr::get_enclosed_expression(OpConvolution<V0, GaussianSmoothing<D0>, E0> const&);
    template<typename V0, size_t D0, typename E0>
	friend auto& expr::get_enclosed_expression(OpConvolution<V0, GaussianSmoothing<D0>, E0>&);
    template<typename V0, size_t D0, typename E0>
	friend auto const& expr::get_result_data(OpConvolution<V0, GaussianSmoothing<D0>, E0> const&);
    template<typename V0, size_t D0, typename E0>
	friend auto& expr::get_result_data(OpConvolution<V0, GaussianSmoothing<D0>, E0>&);


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
		grid::scale(g0);
	}


protected:
	
	expr::ConvolutionData<D> compute;

};



template<typename coeff_t, typename V, size_t D, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<V, GaussianSmoothing<D>, E> const& b)
{
	return symphas::internal::make_convolution::get(value * b.value, b.smoother, expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename tensor_t, size_t D, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<tensor_t, GaussianSmoothing<D>, E> const& b)
{
	return (value * b.value) * symphas::internal::make_convolution::get(OpIdentity{}, b.smoother, expr::get_enclosed_expression(b));
}

template<typename tensor_t, typename V, size_t D, typename E,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V> && expr::eval_type<E>::rank == 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, GaussianSmoothing<D>, E> const& b)
{
	return symphas::internal::make_convolution::get(tensor * b.value, b.smoother, expr::get_enclosed_expression(b));
}

template<typename tensor_t, typename V, size_t D, typename E,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V> && expr::eval_type<E>::rank > 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, GaussianSmoothing<D>, E> const& b)
{
	return symphas::internal::make_convolution::get(b.value, b.smoother, tensor * expr::get_enclosed_expression(b));
}

// ******************************************************************************************************************

//! Convolution of an OpTerm with a Gaussian smoothing kernel.
/*!
 * Specialization of the convolution of two arbitrary expressions.
 * 
 * \tparam V Type of the coefficient.
 * \tparam D the dimension of the Gaussian smoothing kernel.
 * \tparam G The type of variable data.
 */
template<typename V, size_t D, typename G>
struct OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> : 
	OpExpression<OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>>
{
	using E = OpTerm<V, G>;
	using G_T = typename expr::eval_type<E>::type;
	using G0 = Grid<G_T, D>;

	OpConvolution() : g0{ 0 }, data{ 0 }, value{ V{} }, smoother{ GaussianSmoothing<D>() } {}

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
	template<typename V0, typename V1, typename std::enable_if_t<std::is_convertible<mul_result_t<V0, V1>, V>::value, int> = 0>
	OpConvolution(V0 value, GaussianSmoothing<D> const& smoother, OpTerm<V1, G> const& a) :
		g0{ expr::data_dimensions(smoother) }, data{ expr::data(a) }, 
		value{ value * expr::coeff(a) }, smoother{smoother},
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/; }
	OpConvolution(V value, GaussianSmoothing<D> const& smoother, G grid) :
		g0{ expr::data_dimensions(smoother) }, data{ grid }, 
		value{ value }, smoother{ smoother }, 
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/ }


	OpConvolution(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& other) :
		g0{ expr::data_dimensions(other.smoother) }, data{ other.data }, 
		value{ other.value }, smoother{ other.smoother },
		compute{ expr::BaseData<G>::get(data), g0 } { /*update();*/ }
	
	OpConvolution(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>&& other) noexcept :
		OpConvolution()
	{
		swap(*this, other);
	}

	friend void swap(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& first, 
		OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& second)
	{
		using std::swap;
		swap(first.g0, second.g0);
		swap(first.data, second.data);
		swap(first.value, second.value);
		swap(first.smoother, second.smoother);
		swap(first.compute, second.compute);
	}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * g0[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_convolution::get_g(-value, smoother, data);
	}

	template<typename V1>
	auto operator+(OpConvolution<V1, GaussianSmoothing<D>, E> const& other) const
	{
		return symphas::internal::make_convolution::get_g(value + other.value, smoother, data);
	}

	template<typename V1>
	auto operator-(OpConvolution<V1, GaussianSmoothing<D>, E> const& other) const
	{
		return symphas::internal::make_convolution::get_g(value - other.value, smoother, data);
	}

#ifdef PRINTABLE_EQUATIONS

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

#endif

protected:

	G0 g0;							//!< Grid storing the result of the convolution.
	G data;							//!< The data from the OpTerm.

public:

	V value;						//!< Value multiplying the result of this convolution.
	GaussianSmoothing<D> smoother;	//!< The smoothing kernel.

    template<typename V0, size_t D0, typename G0>
	friend auto expr::get_enclosed_expression(OpConvolution<V0, GaussianSmoothing<D0>, OpTerm<OpIdentity, G0>> const&);
    template<typename V0, size_t D0, typename G0>
	friend auto expr::get_enclosed_expression(OpConvolution<V0, GaussianSmoothing<D0>, OpTerm<OpIdentity, G0>>&);
    template<typename V0, size_t D0, typename G0>
	friend auto const& expr::get_result_data(OpConvolution<V0, GaussianSmoothing<D0>, OpTerm<OpIdentity, G0>> const&);
    template<typename V0, size_t D0, typename G0>
	friend auto& expr::get_result_data(OpConvolution<V0, GaussianSmoothing<D0>, OpTerm<OpIdentity, G0>>&);


	void update()
	{
		expr::prune::update(data);
		symphas::dft::fftw_execute(compute.p_in_out);

		auto f = [&](iter_type ft_i, iter_type i)
		{
			compute.in_1[i][0] = smoother.eval(ft_i) * compute.out_0[i][0];
			compute.in_1[i][1] = smoother.eval(ft_i) * compute.out_0[i][1];

		};

		symphas::dft::iterate_rc<G_T, D>(f, g0.dims);
		symphas::dft::fftw_execute(compute.p_out_in);
		grid::scale(g0);
	}

protected:

	expr::ConvolutionData<D> compute;

};


template<typename V0, size_t D, typename T, typename G>
OpConvolution(V0, GaussianSmoothing<D>, OpTerm<T, G>)->OpConvolution<mul_result_t<V0, T>, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>;
template<typename V, size_t D, typename E>
OpConvolution(V, GaussianSmoothing<D>, OpExpression<E>)->OpConvolution<V, GaussianSmoothing<D>, E>;
template<typename V, typename E1, typename E2>
OpConvolution(V, OpExpression<E1>, OpExpression<E2>)->OpConvolution<V, E1, E2>;



template<typename coeff_t, typename V2, size_t D, typename G2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<V2, GaussianSmoothing<D>, OpTerm<OpIdentity, G2>> const& b)
{
	return symphas::internal::make_convolution::get_g(value * b.value, b.smoother, b.op_g);
}

template<typename coeff_t, typename tensor_t, size_t D, typename G2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<tensor_t, GaussianSmoothing<D>, OpTerm<OpIdentity, G2>> const& b)
{
	return (value * b.value) * symphas::internal::make_convolution::get_g(OpIdentity{}, b.smoother, b.op_g);
}

template<typename tensor_t, typename V, size_t D, typename G2,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V>&& expr::eval_type<OpTerm<OpIdentity, G2>>::rank == 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G2>> const& b)
{
	return symphas::internal::make_convolution::get(tensor * b.value, b.smoother, b.op_g);
}

template<typename tensor_t, typename V, size_t D, typename G2,
	typename std::enable_if_t<(expr::is_tensor<tensor_t> && !expr::is_tensor<V> && expr::eval_type<OpTerm<OpIdentity, G2>>::rank > 0), int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G2>> const& b)
{
	return symphas::internal::make_convolution::get(b.value, b.smoother, tensor * b.op_g);
}


// ************************************************************************************************

namespace symphas::internal
{

	template<typename V, typename E1, typename E2>
	auto construct_convolution(V const& value, OpExpression<E1> const& e1, OpExpression<E2> const& e2)
	{
		return OpConvolution<V, E1, E2>(value, *static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
	}

	template<typename V, typename E2>
	auto construct_convolution(V const& value, OpVoid, OpExpression<E2> const&)
	{
		return OpVoid{};
	}

	template<typename V, typename E1>
	auto construct_convolution(V const& value, OpExpression<E1> const&, OpVoid)
	{
		return OpVoid{};
	}

	template<typename V, typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs),
		typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank == 0), int> = 0>
	auto make_convolution_tensor(V const& value, OpExpression<E1> const& e1, OpExpression<E2> const& e2, std::index_sequence<Rs...>)
	{
		return (construct_convolution(
			expr::make_column_vector<Rs, R>() * value,
			expr::make_row_vector<Rs, R>() * (*static_cast<const E1*>(&e1)), 
			*static_cast<const E2*>(&e2))
		+ ...);
	}

	template<typename V, typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs),
		typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 && expr::eval_type<E2>::rank > 0), int> = 0>
		auto make_convolution_tensor(V const& value, OpExpression<E1> const& e1, OpExpression<E2> const& e2, std::index_sequence<Rs...>)
	{
		return (construct_convolution(
			expr::make_column_vector<Rs, R>() * value,
			*static_cast<const E1*>(&e1),
			expr::make_row_vector<Rs, R>() * (*static_cast<const E2*>(&e2)))
			+ ...);
	}

	template<typename A, typename B>
	auto make_convolution::get(A&& a, B&& b)
	{
		return get(OpIdentity{}, std::forward<A>(a), std::forward<B>(b));
	}

	template<typename V, typename E1, typename E2>
	auto make_convolution::get(V const& value, OpExpression<E1> const& e1, OpExpression<E2> const& e2)
	{
		constexpr size_t R = fixed_max<expr::eval_type<E1>::rank, expr::eval_type<E2>::rank>;
		if constexpr (R > 0)
		{
			return make_convolution_tensor(value, *static_cast<E1 const*>(&e1), *static_cast<E2 const*>(&e2), std::make_index_sequence<R>{});
		}
		else
		{
			return OpConvolution<V, E1, E2>(value, *static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
		}
	}

	template<typename V, size_t D, typename E>
	auto make_convolution::get(V const& value, GaussianSmoothing<D> const& smoother, OpExpression<E> const& e)
	{
		constexpr size_t R = expr::eval_type<E>::rank;
		if constexpr (R > 0)
		{
			return make_convolution_tensor(value, smoother, *static_cast<E const*>(&e), std::make_index_sequence<R>{});
		}
		else
		{
			return OpConvolution<V, GaussianSmoothing<D>, E>(value, smoother, *static_cast<E const*>(&e));
		}
	}

	template<typename V, size_t D, typename T, typename G>
	auto make_convolution::get(V const&  value, GaussianSmoothing<D> const& smoother, OpTerm<T, G> const& e)
	{
		constexpr size_t R = expr::eval_type<OpTerm<T, G>>::rank;
		if constexpr (R > 0)
		{
			return make_convolution_tensor(value, smoother, e, std::make_index_sequence<R>{});
		}
		else
		{
			return OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>(value, smoother, e);
		}
	}

	template<typename V, size_t D, typename G>
	auto make_convolution::get_g(V const&  value, GaussianSmoothing<D> const& smoother, G grid)
	{
		return OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>(value, smoother, grid);
	}

}




