
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
 * PURPOSE: Defines the foundational symbolic algebra elements, such as the
 * base expression object and the binary operations.
 *
 * ***************************************************************************
 */

#pragma once

#include <iostream>

#include "expressionlib.h"



// ******************************************************************************************

//! \cond


#define SYEX_BINARY_FMT "%s %s %s"
#define SYEX_ADD_SEP " + "
#define SYEX_SUB_SEP " - "
#define SYEX_BINARY_FMT_LEN (sizeof(SYEX_ADD_SEP) / sizeof(char) - 1)



#ifdef LATEX_PLOT

#define SYEX_MUL_SEP_OP " "
#define SYEX_DIV_SEP_OP "}{"
#define SYEX_NLV_SEP_OP " "

//#define SYEX_MUL_FMT_AA "\\left["
#define SYEX_MUL_FMT_AA ""
#define SYEX_MUL_FMT_AB "\\left("
#define SYEX_MUL_FMT_BA "\\right)"
//#define SYEX_MUL_FMT_BB "\\right]"
#define SYEX_MUL_FMT_BB ""

#define SYEX_DIV_SEP SYEX_DIV_SEP_OP
#define SYEX_DIV_FMT_A "\\frac{"
#define SYEX_DIV_FMT_B "}" 

#else

#define SYEX_MUL_SEP_OP "*"
#define SYEX_DIV_SEP_OP "/"
#define SYEX_NLV_SEP_OP "*"

//#define SYEX_MUL_FMT_AA "["
#define SYEX_MUL_FMT_AA ""
#define SYEX_MUL_FMT_AB "("
#define SYEX_MUL_FMT_BA ")"
//#define SYEX_MUL_FMT_BB "]"
#define SYEX_MUL_FMT_BB ""

#define SYEX_DIV_SEP  ") " SYEX_DIV_SEP_OP " ("
#define SYEX_DIV_FMT_A SYEX_MUL_FMT_AB
#define SYEX_DIV_FMT_B SYEX_MUL_FMT_BA



#endif

#define SYEX_MUL_SEP  SYEX_MUL_FMT_BA SYEX_MUL_SEP_OP SYEX_MUL_FMT_AB
#define SYEX_MUL_SEP_A  SYEX_MUL_FMT_BA SYEX_MUL_SEP_OP
#define SYEX_MUL_SEP_B  SYEX_MUL_SEP_OP SYEX_MUL_FMT_AB
#define SYEX_MUL_SEP_AB  SYEX_MUL_SEP_OP

#define SYEX_MUL_FMT_A SYEX_MUL_FMT_AA SYEX_MUL_FMT_AB
#define SYEX_MUL_FMT_B SYEX_MUL_FMT_BA SYEX_MUL_FMT_BB

#define SYEX_MUL_FMT SYEX_MUL_FMT_A "%s" SYEX_MUL_SEP "%s" SYEX_MUL_FMT_B
#define SYEX_MUL_FMT_LEN (STR_ARR_LEN(SYEX_MUL_FMT_A SYEX_MUL_SEP SYEX_MUL_FMT_B) - 1)

#define SYEX_DIV_FMT SYEX_DIV_FMT_A "%s" SYEX_DIV_SEP "%s" SYEX_DIV_FMT_B
#define SYEX_DIV_FMT_LEN (STR_ARR_LEN(SYEX_DIV_FMT_A SYEX_DIV_SEP SYEX_DIV_FMT_B) - 1)


#define SYEX_FRA_SEP SYEX_DIV_SEP_OP
#define SYEX_FRA_FMT_A SYEX_DIV_FMT_A
#define SYEX_FRA_FMT_B SYEX_MUL_FMT_B 

#define SYEX_FRA_FMT SYEX_DIV_FMT_A "%zd" SYEX_DIV_SEP_OP "%zd" SYEX_DIV_FMT_B
#define SYEX_FRA_FMT_LEN (STR_ARR_LEN(SYEX_DIV_FMT_A SYEX_DIV_SEP_OP SYEX_DIV_FMT_B) - 1)

#define SYEX_LAMBDA_FUNC_FMT_A "%s("
#define SYEX_LAMBDA_FUNC_FMT_B ")"
#define SYEX_LAMBDA_FUNC_FMT SYEX_DIV_FMT_A SYEX_LAMBDA_FUNC_FMT_A "%s" SYEX_LAMBDA_FUNC_FMT_B
#define SYEX_LAMBDA_FUNC_FMT_LEN (STR_ARR_LEN(SYEX_LAMBDA_FUNC_FMT_A SYEX_LAMBDA_FUNC_FMT_B) - 3)

#ifdef LATEX_PLOT
#define SYEX_POW_SEP_A "^{"
#define SYEX_POW_SEP_B "}"
#define SYEX_POW_DIV_SEP "/"
#else
#define SYEX_POW_SEP_A "^"
#define SYEX_POW_SEP_B ""
#define SYEX_POW_DIV_SEP "/"
#endif



//! \endcond


//! The maximum supported power that an expression can be evaluated at.
#define MAX_EXPONENT 20


//! The display precision of values which appear in expressions.
/*!
 * Defines the number of digits which appear in floating point numbers
 * for values from expressions that are printed to the screen.
 */
#define EXPR_VALUE_DISPLAY_PRECISION 3

namespace symphas::internal
{
	struct tensor_cast
	{
		template<typename T, size_t... Ns>
		static auto const& cast(OpTensor<T, Ns...> const& tensor)
		{
			return tensor.value;
		}

		template<typename T, size_t... Ns>
		static auto& cast(OpTensor<T, Ns...>& tensor)
		{
			return tensor.value;
		}
	};

	template<typename T, size_t N0, size_t N1, size_t NA, size_t NB, size_t... Ns, 
		typename S, size_t M0, size_t M1, size_t MA, size_t MB, size_t... Ms>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB, Ns...>, OpTensor<S, M0, M1, MA, MB, Ms...>) = delete;

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t N2, size_t N3, size_t NA, size_t NB, size_t NC>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB> const& a, OpTensor<S, N2, N3, NB, NC> const& b)
	{
		if constexpr (N1 == N2)
		{
			T at = tensor_cast::cast(a);
			S bt = tensor_cast::cast(b);
			return expr::make_tensor<N0, N3, NA, NC>(at * bt);
		}
		else
		{
			return OpVoid{};
		}
	}

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t NA, size_t NB, 
		typename std::enable_if_t<(NA > 0), int> = 0>
	auto tensor_multiply(OpTensor<T, N0, NA> const& a, OpTensor<S, 0, N1, 1, NB> const& b)
	{
		T at = tensor_cast::cast(a);
		S bt = tensor_cast::cast(b);
		return expr::make_tensor<N0, N1, NA, NB>(at * bt);
	}

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t N2, size_t NA, size_t NB>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB> const& a, OpTensor<S, N2, NB> const& b)
	{
		if constexpr (N1 == N2)
		{
			T at = tensor_cast::cast(a);
			S bt = tensor_cast::cast(b);
			
			if constexpr (NA == 1)
			{
				return expr::make_literal(at * bt);
			}
			else
			{
				return expr::make_tensor<N0, 0, NA, 1>(at * bt);
			}
		}
		else
		{
			return OpVoid{};
		}
	}

	template<typename T, size_t... Ms>
	auto tensor_multiply(OpTensor<OpIdentity, 0, 0> const& a, OpTensor<T, Ms...> const& b)
	{
		return OpIdentity{};
	}

	template<typename T, size_t... Ms>
	auto tensor_multiply(OpTensor<T, Ms...> const& a, OpTensor<OpIdentity, 0, 0> const& b)
	{
		return OpIdentity{};
	}

	template<typename T, typename S, size_t... Ms>
	auto tensor_multiply(OpTensor<T, Ms...> const& a, OpTensor<S, 0, 1> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<Ms...>(v);
	}

	template<typename T, typename S, size_t... Ms>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, Ms...> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<Ms...>(v);
	}

	template<typename T, typename S, size_t N0, size_t N>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, 0, N0, 1, N> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<0, N0, 1, N>(v);
	}


	template<typename S>
	auto tensor_multiply(OpTensor<OpIdentity, 0, 0> const& a, OpTensor<S, 0, 1> const& b)
	{
		return OpIdentity{};
	}

	template<typename S>
	auto tensor_multiply(OpTensor<S, 0, 0, 1, 1> const& a, OpTensor<OpIdentity, 0, 0> const& b)
	{
		return OpIdentity{};
	}




	template<typename T, typename S>
	auto tensor_multiply(OpTensor<T, 0, 1> const& a, OpTensor<S, 0, 1> const& b)
	{
		return symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
	}

	template<typename T, typename S>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, 0, 1> const& b)
	{
		return symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
	}

	template<typename T, typename S>
	auto tensor_multiply(OpTensor<S, 0, 1> const& a, OpTensor<T, 0, 0, 1, 1> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<0, 0, 1, 1>(v);
	}

	template<size_t N0, size_t N1, size_t NA, size_t NB, typename T>
	auto tensor_as_coeff(T const& value)
	{
		//if constexpr (N0 == 0 && N1 == 0 && NA == 1 && NB == 1)
		//{
		//	return value.eval();
		//}
		//else
		{
			using elem_type = std::invoke_result_t<decltype(&T::eval), T, iter_type>;
			any_matrix_t<elem_type, NA, NB> matrix;
			matrix[N0][N1] = value.eval();
			return matrix;
		}
	}

	template<size_t N0, size_t NA, typename T>
	auto tensor_as_coeff(T const& value)
	{
		if constexpr (N0 == 0 && NA == 0)
		{
			return value;
		}
		else
		{
			using elem_type = std::invoke_result_t<decltype(&T::eval), T, iter_type>;
			any_vector_t<elem_type, NA> vector;
			vector[N0] = value.eval();
			return vector;
		}
	}

	template<typename>
	struct tensor_to_vector;

	template<typename T, size_t N1, size_t NA>
	struct tensor_to_vector<OpTensor<T, N1, NA>>
	{
		using type = any_vector_t<std::invoke_result_t<decltype(&T::eval), T, iter_type>, NA>;
	};

	template<typename T, size_t N1, size_t N2, size_t NA, size_t NB>
	struct tensor_to_vector<OpTensor<T, N1, N2, NA, NB>>
	{
		using type = any_vector_t<any_vector_t<std::invoke_result_t<decltype(&T::eval), T, iter_type>, NB>, NA>;
	};


	template<typename T>
	using tensor_to_vector_t = typename tensor_to_vector<T>::type;

	template<size_t N1, size_t NA, typename T, typename S>
	void set_vector(any_vector_t<T, NA>& vector, S const& value)
	{
		vector[N1] = value;
	}

	template<size_t N1, size_t N2, size_t NA, size_t NB, typename T, typename S>
	void set_vector(any_vector_t<any_vector_t<T, NB>, NA>& vector, S const& value)
	{
		vector[N1][N2] = value;
	}
}


//! A value of a tensor, acts as a coefficient.
/*! 
 * A tensor turns the expression for which it is a coefficient into an entry of a
 * tensor. 
 * 
 * A tensor of rank 1 is a column vector, and of rank 2 is a matrix. A row vector is a
 * tensor of rank 2 and has only 1 row. Multiplication only for rank 1 and 2 vectors is
 * supported, and there are no tensor products.
 */
template<typename T, size_t... Ns>
struct OpTensor : OpExpression<OpTensor<T, Ns...>>
{

	explicit OpTensor(T const& entry) : value{ entry } {}
	explicit OpTensor(T&& entry) noexcept : value{ std::move(entry) } {}
	constexpr OpTensor() : value{ T{} } {}

	auto eval(iter_type n = 0) const
	{
		return symphas::internal::tensor_as_coeff<Ns...>(cast());
	}

	template<typename S, size_t... Ms>
	auto operator*(OpTensor<S, Ms...> const& a) const
	{
		return symphas::internal::tensor_multiply(*this, a);
	}

	template<typename S>
	auto operator-(OpTensor<S, Ns...> const& a) const
	{
		return expr::make_tensor<Ns...>(cast() - symphas::internal::tensor_cast::cast(a));
	}

	template<typename S>
	auto operator+(OpTensor<S, Ns...> const& a) const
	{
		return expr::make_tensor<Ns...>(cast() + symphas::internal::tensor_cast::cast(a));
	}

	auto operator-() const
	{
		return expr::make_tensor<Ns...>(-cast());
	}

	operator symphas::internal::tensor_to_vector_t<OpTensor<T, Ns...>>() const
	{
		symphas::internal::tensor_to_vector_t<OpTensor<T, Ns...>> vector;
		symphas::internal::set_vector<Ns...>(vector, symphas::internal::tensor_cast::cast(*this));
		return vector;
	}

	operator T const& () const = delete;
	operator T& () = delete;


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, "", *this);
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, "", *this);
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(*this);
	}

#endif

	friend struct symphas::internal::tensor_cast;

protected:

	const T& cast() const
	{
		return symphas::internal::tensor_cast::cast(*this);
	}

	T value;
};

namespace expr
{
	template<typename T, size_t... Ns>
	auto coeff(OpTensor<T, Ns...> const& tensor)
	{
		return symphas::internal::tensor_cast::cast(tensor);
	}
}

namespace symphas::internal
{
	using tensor_cancel = OpTensor<OpIdentity, 0, 0>;
}


template<typename coeff_t, typename T, size_t... Ns,
	typename std::enable_if_t<(expr::is_identity<coeff_t> && expr::is_fraction<coeff_t>), int> = 0>
auto operator*(coeff_t const& value, OpTensor<T, Ns...> const& tensor)
{
	return expr::make_tensor<Ns...>(value * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t... Ns>
auto operator*(OpLiteral<T1> const& value, OpTensor<T2, Ns...> const& tensor)
{
	return expr::make_tensor<Ns...>(value * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t... Ns>
auto operator*(OpTensor<T1, Ns...> const& tensor, OpLiteral<T2> const& value)
{
	return expr::make_tensor<Ns...>(symphas::internal::tensor_cast::cast(tensor) * value);
}


template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 1> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	return expr::make_tensor<0, N1, 1, NB>(value[0] * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 2> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	auto v = (*static_cast<T2 const*>(&tensor));
	return expr::make_tensor<0, N1, 2, NB>(value[0] * v)
		+ expr::make_tensor<1, N1, 2, NB>(value[1] * v);
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 3> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	auto v = (*static_cast<T2 const*>(&tensor));
	return expr::make_tensor<0, N1, 3, NB>(value[0] * v)
		+ expr::make_tensor<1, N1, 3, NB>(value[1] * v)
		+ expr::make_tensor<2, N1, 3, NB>(value[2] * v);
}


template<typename T1, typename T2, size_t N0, size_t NA>
auto operator*(OpTensor<T1, N0, NA> const& tensor, any_vector_t<T2, 1> const& value)
{
	return expr::make_tensor<N0, NA>(value[0] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NA>
auto operator*(OpTensor<T1, N0, N1, NA, 2> const& tensor, any_vector_t<T2, 2> const& value)
{
	return expr::make_tensor<N0, NA>(value[N1] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NA>
auto operator*(OpTensor<T1, N0, N1, NA, 3> const& tensor, any_vector_t<T2, 3> const& value)
{
	return expr::make_tensor<N0, NA>(value[N1] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T, size_t D>
auto operator*(OpTensor<OpIdentity, 0, 0> const&, OpLiteral<any_vector_t<T, D>> const&)
{
	return OpIdentity{};
}

template<typename T, size_t D>
auto operator*(OpTensor<OpIdentity, 0, 0> const&, any_vector_t<T, D> const&)
{
	return OpIdentity{};
}


//! Representation of a constant.
/*!
 * Stores a single value of any type. 
 * 
 * \tparam T The type of constant being stored.
 */
template<typename T>
struct OpLiteral : OpExpression<OpLiteral<T>>
{
	T value;

	OpLiteral(T value) : value{ value } {}
	constexpr OpLiteral(OpIdentity);
	constexpr OpLiteral(OpNegIdentity);
	constexpr OpLiteral(OpVoid);
	constexpr OpLiteral() : value{ T{} } {}

	inline T eval(iter_type = 0) const
	{
		return value;
	}

	operator const T() const
	{
		return value;
	}

	operator T&()
	{
		return value;
	}


	auto operator^(size_t exp) const
	{
		return expr::make_literal(std::pow(value, exp));
	}

	template<typename S>
	auto operator*(OpLiteral<S> const& other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}
	auto operator-() const
	{
		return expr::make_literal(-value);
	}

	template<typename S>
	auto operator+(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value + a.value);
		return v;
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const;
	size_t print(char* out) const;
	size_t print_length() const;

#endif

};

#ifdef PRINTABLE_EQUATIONS

template<>
inline size_t OpLiteral<double>::print(FILE* out) const
{
	return fprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", value);
}

template<>
inline size_t OpLiteral<int>::print(FILE* out) const
{
	return fprintf(out, "%d", value);
}

template<>
inline size_t OpLiteral<complex_t>::print(FILE* out) const
{
	return fprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf + %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), imag(value));
}

template<>
inline size_t OpLiteral<double>::print(char* out) const
{
	return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", value);
}

template<>
inline size_t OpLiteral<int>::print(char* out) const
{
	return sprintf(out, "%d", value);
}

template<>
inline size_t OpLiteral<complex_t>::print(char* out) const
{
	if (real(value) == 0)
	{
		if (imag(value) < 0)
		{
			return sprintf(out, "-%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", std::abs(imag(value)));
		}
		else if (imag(value) == 0)
		{
			return sprintf(out, "0");
		}
		else
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", imag(value));
		}
	}
	else
	{
		if (imag(value) < 0)
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf - %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), std::abs(imag(value)));
		}
		else if (imag(value) == 0)
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", real(value)); //1
		}
		else
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf + %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), imag(value));
		}
	}
}

template<>
inline size_t OpLiteral<double>::print_length() const
{
	size_t len = symphas::lib::num_digits(static_cast<int>(std::abs(value)));
	if (value < 0)
	{
		return len + 2 + EXPR_VALUE_DISPLAY_PRECISION;
	}
	else
	{
		return len + 1 + EXPR_VALUE_DISPLAY_PRECISION;
	}
}

template<>
inline size_t OpLiteral<int>::print_length() const
{
	size_t len = symphas::lib::num_digits(std::abs(value));
	if (value < 0)
	{
		return len + 1;
	}
	else
	{
		return len;
	}
}

template<>
inline size_t OpLiteral<complex_t>::print_length() const
{
	size_t len = 0;
	len += symphas::lib::num_digits(static_cast<int>(std::abs(real(value))));
	len += symphas::lib::num_digits(static_cast<int>(std::abs(imag(value))));

	if (real(value) == 0)
	{
		if (imag(value) < 0)
		{
			len += 3 + EXPR_VALUE_DISPLAY_PRECISION;
		}
		else if (imag(value) == 0)
		{
			len += 1;
		}
		else
		{
			len += 2 + EXPR_VALUE_DISPLAY_PRECISION;
		}
	}
	else
	{
		if (imag(value) == 0)
		{
			len += 1 + EXPR_VALUE_DISPLAY_PRECISION;
		}
		else
		{
			len += 6 + EXPR_VALUE_DISPLAY_PRECISION * 2;
		}
	}


	if (real(value) < 0)
	{
		return len + 1;
	}
	else
	{
		return len;
	}
}

#endif


template<typename T, typename I>
struct OpCoeff : OpLiteral<T>
{
	using parent_type = OpLiteral<T>;
	using parent_type::parent_type;

	OpCoeff(T* data, len_type len, len_type stride = 1) :
		parent_type(T{}), data{ (len > 0) ? new T[len] : nullptr }, len{ len }
	{
		for (iter_type i = 0; i < len; ++i)
		{
			this->data[i * stride] = data[i * stride];
		}
	}

	constexpr OpCoeff() : OpCoeff(nullptr, 0) {}

	OpCoeff(OpCoeff<T, I> const& other) : OpCoeff(other.data, other.len) {}
	OpCoeff(OpCoeff<T, I>&& other) : OpCoeff()
	{
		swap(*this, other);
	}

	auto operator=(OpCoeff<T, I> other)
	{
		swap(*this, other);
	}

	friend void swap(OpCoeff<T, I>& first, OpCoeff<T, I>& second)
	{
		using std::swap;
		swap(first.data, second.data);
		swap(first.len, second.len);
	}

	template<int N, int P>
	auto operator()(expr::symbols::i_<N, P>) const
	{
		return OpCoeff<T, expr::symbols::i_<N, P>>(data, len);
	}

	auto operator[](iter_type i) const
	{
		return expr::make_literal(data[i]);
	}

	T* data;
	len_type len;


	~OpCoeff()
	{
		delete[] data;
	}

};

template<typename T>
struct OpCoeff<T, int> : OpCoeff<T, void>
{
	using parent_type = OpCoeff<T, void>;
	using parent_type::parent_type;

	auto eval(iter_type = 0) const
	{
		return parent_type::operator[](*ptr);
	}

	void set_pointer(int* ptr)
	{
		this->ptr = ptr;
	}

	size_t print(FILE* out) const
	{
		OpLiteral<T>::value = eval();
		return OpLiteral<T>::print(out);
	}

	size_t print(char* out) const
	{
		OpLiteral<T>::value = eval();
		return OpLiteral<T>::print(out);
	}

	size_t print_length() const
	{
		OpLiteral<T>::value = eval();
		return OpLiteral<T>::print_length();
	}

	T* data;
	len_type stride;
	int* ptr;
};


template<typename I>
struct OpCoeff<void, I> {};

template<>
struct OpCoeff<void, void> 
{
	OpCoeff(iter_type i) : i{ i } {}
	int i;
};

template<typename T>
OpCoeff(T*, len_type) -> OpCoeff<T, void>;
OpCoeff(iter_type) -> OpCoeff<void, void>;


template<size_t N, size_t D>
struct OpNegFractionLiteral;
template<size_t N, size_t D>
struct OpFractionLiteral;

template<size_t N, size_t D>
struct OpFractionLiteral : OpExpression<OpFractionLiteral<N, D>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return static_cast<scalar_t>(N) / D;
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpNegFractionLiteral<N, D>{};
	}



#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return fprintf(out, SYEX_FRA_FMT, N, D);
	}

	size_t print(char* out) const
	{
		return sprintf(out, SYEX_FRA_FMT, N, D);
	}

	size_t print_length() const
	{
		size_t len1 = symphas::lib::num_digits<N>();
		size_t len2 = symphas::lib::num_digits<N>();
		return len1 + len2 + SYEX_FRA_FMT_LEN;
	}

#endif

};

template<size_t N>
struct OpFractionLiteral<N, 1> : OpExpression<OpFractionLiteral<N, 1>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return static_cast<scalar_t>(N);
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpNegFractionLiteral<N, 1>{};
	}



#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return fprintf(out, "%zd", N);
	}

	size_t print(char* out) const
	{
		return sprintf(out, "%zd", N);
	}

	size_t print_length() const
	{
		return symphas::lib::num_digits<N>();
	}

#endif

};


template<size_t N, size_t D>
struct OpNegFractionLiteral : OpExpression<OpNegFractionLiteral<N, D>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return -static_cast<scalar_t>(N) / D;
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpFractionLiteral<N, D>{};
	}



#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, "-");
		n += OpFractionLiteral<N, D>{}.print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out, "-");
		n += OpFractionLiteral<N, D>{}.print(out + n);
		return n;
	}

	size_t print_length() const
	{
		return OpFractionLiteral<N, D>{}.print_length() + 1;
	}

#endif

};


template<typename T, size_t N, size_t D>
auto operator*(OpLiteral<T> const& a, OpFractionLiteral<N, D>)
{
	return expr::make_literal(a.value * OpFractionLiteral<N, D>{}.eval());
}

template<typename T, size_t N, size_t D>
auto operator*(OpLiteral<T> const& a, OpNegFractionLiteral<N, D>)
{
	return expr::make_literal(a.value * OpNegFractionLiteral<N, D>{}.eval());
}

// ******************************************************************************************


//! Specialized version which returns multiplicative identity.
template<>
inline auto expr::make_literal<OpIdentity>(OpIdentity const&)
{
	return OpIdentity{};
}

//! Specialized version which returns negative multiplicative identity.
template<>
inline auto expr::make_literal<OpNegIdentity>(OpNegIdentity const&)
{
	return OpNegIdentity{};
}

//! Specialized version which returns additive identity.
template<>
inline auto expr::make_literal<OpVoid>(OpVoid const&)
{
	return OpVoid{};
}


template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpIdentity) : value{ symphas::lib::get_identity<T>() } {}
template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpNegIdentity) : value{ -symphas::lib::get_identity<T>() } {}
template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpVoid) : value{ 0 } {}


namespace expr
{
	namespace symbols
	{
		//! Common fractions.
		constexpr OpFractionLiteral<1, 2> _2{};
		constexpr OpFractionLiteral<1, 4> _4{};

		//! Pi to 7th order
		constexpr OpFractionLiteral<355, 113> Pi{};
	}
}

// ******************************************************************************************

namespace symphas::internal
{
	using namespace symphas::lib;


	template<int... Is>
	constexpr static bool all_positive(std::integer_sequence<int, Is...>)
	{
		return ((Is > 0) && ...) || (sizeof...(Is) == 1 && ((Is == 0) && ...));
	}

	template<size_t... Ns>
	struct construct_tensor
	{

		static constexpr size_t N = sizeof...(Ns);
		using pos_t = seq_join_t<types_before_index<N / 2, std::index_sequence<Ns>...>>;
		using dim_t = seq_join_t<types_after_at_index<N / 2, std::index_sequence<Ns>...>>;

		template<typename T>
		auto operator()(T const& v)
		{
			static_assert(symphas::internal::all_positive(seq_sub_t<dim_t, pos_t>{}),
				"incorrect tensor arguments");
			return OpTensor<T, Ns...>{ v };
		}
	};

	template<size_t N0, size_t NA>
	struct construct_tensor<N0, 0, NA, 1> : construct_tensor<N0, NA> {};


	template<>
	struct construct_tensor<0, 0, 1, 1>
	{
		template<typename T>
		auto operator()(T const& v)
		{
			return OpTensor<T, 0, 0, 1, 1>{ v };
		}
	};
}


template<size_t N0, size_t N1, size_t... Ns, typename T>
constexpr auto expr::make_tensor(T const& v)
{
	return symphas::internal::construct_tensor<N0, N1, Ns...>{}(expr::make_literal(v));
}

template<size_t N0, size_t N1, size_t... Ns>
constexpr auto expr::make_tensor()
{
	return expr::make_tensor<N0, N1, Ns...>(OpIdentity{});
}

template<size_t I, size_t N, typename T>
constexpr auto expr::make_column_vector(T&& v)
{
	return make_tensor<I, N>(std::forward<T>(v));
}

template<size_t I, size_t N>
constexpr auto expr::make_column_vector()
{
	return make_column_vector<I, N>(OpIdentity{});
}

template<size_t I, size_t N, typename T>
constexpr auto expr::make_row_vector(T&& v)
{
	return make_tensor<0, I, 1, N>(std::forward<T>(v));
}

template<size_t I, size_t N>
constexpr auto expr::make_row_vector()
{
	return make_row_vector<I, N>(OpIdentity{});
}


template<size_t R, size_t R0>
auto expr::make_filled_column_vector()
{
	if constexpr (R0 >= R)
	{
		return OpIdentity{};
	}
	else
	{
		auto c = expr::make_column_vector<R0, R>();
		if constexpr (R0 == R - 1)
		{
			return c;
		}
		else
		{
			return c + make_filled_column_vector<R, R0 + 1>();
		}
	}
}

template<size_t R, size_t R0>
auto expr::make_filled_row_vector()
{
	if constexpr (R0 >= R)
	{
		return OpIdentity{};
	}
	else
	{
		auto c = expr::make_row_vector<R0, R>();
		if constexpr (R0 == R - 1)
		{
			return c;
		}
		else
		{
			return c + make_filled_row_vector<R, R0 + 1>();
		}
	}
}

template<typename T>
constexpr auto expr::make_literal(T const& v)
{
	if constexpr (expr::is_fraction<T>)
	{
		return T{};
	}
	else if constexpr (expr::is_tensor<T>)
	{
		return v;
	}
	else if constexpr (expr::is_arr_coeff<T>)
	{
		return v;
	}
	else
	{
		return OpLiteral{ v };
	}
}

template<size_t N, size_t D>
constexpr auto expr::make_fraction()
{
	static_assert(D != 0, "dividing by zero");
	if constexpr (N == 0)
	{
		return OpVoid{};
	}
	else if constexpr (N == D)
	{
		return OpIdentity{};
	}
	else
	{
		constexpr size_t GCD = GCD_of<N, D>;
		return OpFractionLiteral<N / GCD, D / GCD>{};
	}
}

template<int I>
constexpr auto expr::make_integer()
{
	if constexpr (I < 0)
	{
		return -expr::make_fraction<static_cast<size_t>(-I), 1>();
	}
	else
	{
		return expr::make_fraction<static_cast<size_t>(I), 1>();
	}
}

// ******************************************************************************************


//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::is_identity<E1> && expr::is_identity<E2>), int> = 0>
//constexpr bool operator==(E1 const& a, E2 const& b)
//{
//	return a.eval() == b.eval();
//}
//
//template<typename A, typename B,
//	typename std::enable_if_t<(expr::is_identity<A> && expr::is_identity<B>), int> = 0>
//constexpr bool operator<(A const& a, B const& b)
//{
//	return a.eval() < b.eval();
//}
//
//template<typename E1, typename E2, 
//	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
//inline bool operator==(E1 const& a, E2 const& b)
//{
//	return OpLiteral<identity_eval_t>(a).eval() == OpLiteral<identity_eval_t>(b).eval();
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
//inline bool operator<(E1 const& a, E2 const& b)
//{
//	return OpLiteral<identity_eval_t>(a).eval() < OpLiteral<identity_eval_t>(b).eval();
//}
//
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator!=(E1 const& a, E2 const& b)
//{
//	return !(a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator<=(E1 const& a, E2 const& b)
//{
//	return (a < b || a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator>(E1 const& a, E2 const& b)
//{
//	return !(a < b || a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator>=(E1 const& a, E2 const& b)
//{
//	return !(a < b);
//}

namespace expr
{
	//! Constructs the addition expression.
	/*!
	 * Directly constructs the addition expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The first term in the addition.
	 * \param b The second term in the addition.
	 * \param es The rest of the terms in the addition.
	 */
	OpVoid make_add();
	OpVoid make_add(OpVoid, OpVoid);
	template<typename E>
	auto make_add(E&& a);
	template<typename E>
	auto make_add(E&& e, OpVoid);
	template<typename E>
	auto make_add(OpVoid, E&& e);
	template<typename... E0s, typename E0>
	auto make_add(OpAdd<E0s...> const& add, E0 const& e0);
	template<typename... E0s>
	auto make_add(OpAdd<E0s...> const& add, OpVoid);
	template<typename... E0s, typename E0>
	auto make_add(E0 const& e0, OpAdd<E0s...> const& add);
	template<typename... E0s>
	auto make_add(OpVoid, OpAdd<E0s...> const& add);
	template<typename A, typename B>
	auto make_add(OpExpression<A> const& a, OpExpression<B> const& b);
	template<typename E1, typename E2, typename E3, typename... Es>
	auto make_add(E1&& a, E2&& b, E3&& c, Es&& ...es);


	template<typename... As, typename... Bs, size_t... Is>
	auto make_add(OpAdd<As...> const& a, OpAdd<Bs...> const& b, std::index_sequence<Is...>)
	{
		return make_add(a, expr::get<Is>(b)...);
	}

	template<typename... As, typename... Bs>
	auto make_add(OpAdd<As...> const& a, OpAdd<Bs...> const& b)
	{
		if constexpr (sizeof...(As) > sizeof...(Bs))
		{
			return make_add(b, a, std::make_index_sequence<sizeof...(As)>{});
		}
		else
		{
			return make_add(a, b, std::make_index_sequence<sizeof...(Bs)>{});
		}
	}
}


template<typename T, size_t D, typename S>
bool operator>(any_vector_t<T, D> const& value, S&&)
{
	return true;
}

template<typename T, size_t D, typename S>
bool operator<(any_vector_t<T, D> const& value, S&&)
{
	return false;
}


namespace symphas::internal
{
	using expr::has_coeff;

#ifdef PRINTABLE_EQUATIONS


	inline size_t print_sep(char* out, const char* sep)
	{
		return sprintf(out, "%s", sep);
	}

	inline size_t print_sep(FILE* out, const char* sep)
	{
		return fprintf(out, "%s", sep);
	}

	template<typename E>
	void print_one(E const& e, FILE* out, size_t& n)
	{
		if (expr::eval(expr::coeff(e)) > 0)
		{
			n += print_sep(out, SYEX_ADD_SEP);
			n += e.print(out);
		}
		else
		{
			n += print_sep(out, SYEX_SUB_SEP);
			n += (-e).print(out);
		}
	}

	template<typename E>
	void print_one(E const& e, char* out, size_t& n)
	{
		if (expr::eval(expr::coeff(e)) > 0)
		{
			n += print_sep(out + n, SYEX_ADD_SEP);
			n += e.print(out + n);
		}
		else
		{
			n += print_sep(out + n, SYEX_SUB_SEP);
			n += (-e).print(out + n);
		}
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename... Es, typename os_type, size_t... Is>
	size_t binary_print(OpAdd<Es...> const& e, os_type* out, std::index_sequence<Is...>)
	{
		size_t n = 0;
		(print_one(expr::get<Is + 1>(e), out, n), ...);
		return n;
	}


#endif

	template<size_t N, typename E>
	struct Nth_type_of_add;

	template<size_t N, typename E0, typename... Es>
	struct Nth_type_of_add<N, OpAdd<E0, Es...>>
	{
		using type = typename Nth_type_of_add<N - 1, OpAdd<Es...>>::type;
	};

	template<size_t N>
	struct Nth_type_of_add<N, OpAdd<>>
	{
		using type = OpAdd<>;
	};

	template<typename E0, typename... Es>
	struct Nth_type_of_add<0, OpAdd<E0, Es...>>
	{
		using type = OpAdd<E0, Es...>;
	};

	template<typename cast_type>
	struct cast_add
	{
		template<typename... Es>
		static cast_type const& cast(OpAdd<Es...> const& adds)
		{
			return *static_cast<cast_type const*>(&adds);
		}

		template<typename... Es>
		static cast_type& cast(OpAdd<Es...>& adds)
		{
			return *static_cast<cast_type*>(&adds);
		}
	};

	template<typename... E0s>
	struct cast_add<symphas::lib::types_list<E0s...>>
	{
		template<typename... Es>
		static OpAdd<E0s...> const& cast(OpAdd<Es...> const& adds)
		{
			return *static_cast<OpAdd<E0s...> const*>(&adds);
		}

		template<typename... Es>
		static OpAdd<E0s...>& cast(OpAdd<Es...>& adds)
		{
			return *static_cast<OpAdd<E0s...>*>(&adds);
		}
	};

	template<>
	struct cast_add<void>
	{
		template<typename... Es, size_t... Is>
		static decltype(auto) cast(OpAdd<Es...> const& adds, std::index_sequence<Is...>)
		{
			return expr::make_add(expr::get<Is>(adds)...);
		}

		template<typename... Es, size_t... Is>
		static decltype(auto) cast(OpAdd<Es...>& adds, std::index_sequence<Is...>)
		{
			return expr::make_add(expr::get<Is>(adds)...);
		}
	};
}


namespace expr
{

	template<size_t N, typename... Es>
	const auto& get(OpAdd<Es...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).data;
	}

	template<size_t N, typename... Es>
	auto& get(OpAdd<Es...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).data;
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...> const& e)
	{
		if constexpr (N == 0)
		{
			return OpVoid{};
		}
		else if constexpr (N == 1)
		{
			return expr::get<0>(e);
		}
		else
		{
			return symphas::internal::cast_add<void>::cast(e, std::make_index_sequence<N>{});
		}
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...>& e)
	{
		if constexpr (N == 0)
		{
			return OpVoid{};
		}
		else if constexpr (N == 1)
		{
			return expr::get<0>(e);
		}
		else
		{
			return symphas::internal::cast_add<void>::cast(e, std::make_index_sequence<N>{});
		}
	}


	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...> const& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, Es...>;
		if constexpr (N + 1 == sizeof...(Es) - 1)
		{
			return symphas::internal::cast_add<ts>::cast(e).data;
		}
		else if constexpr (N + 1 < sizeof...(Es))
		{
			return symphas::internal::cast_add<ts>::cast(e);
		}
		else
		{
			return OpVoid{};
		}
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...>& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, Es...>;
		if constexpr (N + 1 == sizeof...(Es) - 1)
		{
			return symphas::internal::cast_add<ts>::cast(e).data;
		}
		else if constexpr (N + 1 < sizeof...(Es))
		{
			return symphas::internal::cast_add<ts>::cast(e);
		}
		else
		{
			return OpVoid{};
		}
	}

	template<typename E0, typename E1, typename E2, typename... Es>
	const auto& terms_after_first(OpAdd<E0, E1, E2, Es...> const& e)
	{
		return symphas::internal::cast_add<OpAdd<E1, E2, Es...>>::cast(e);
	}

	template<typename E0, typename E1>
	const auto& terms_after_first(OpAdd<E0, E1> const& e)
	{
		return symphas::internal::cast_add<OpAdd<E1>>::cast(e).data;
	}

	template<typename E0, typename E1, typename E2, typename... Es>
	auto& terms_after_first(OpAdd<E0, E1, E2, Es...>& e)
	{
		return symphas::internal::cast_add<OpAdd<E1, E2, Es...>>::cast(e);
	}

	template<typename E0, typename E1>
	auto& terms_after_first(OpAdd<E0, E1>& e)
	{
		return symphas::internal::cast_add<OpAdd<E1>>::cast(e).data;
	}
}

template<>
struct OpAdd<> : OpExpression<OpAdd<>>
{
	inline auto eval(iter_type) const
	{
		return OpVoid{};
	}

#ifdef PRINTABLE_EQUATIONS

	inline int print_length() const
	{
		return -int(SYEX_BINARY_FMT_LEN);
	}

	bool coeff_sign() const
	{
		return false;
	}

#endif
	

};

template<typename E0, typename... Es>
struct OpExpression<OpAdd<E0, Es...>> : OpAdd<Es...>
{
	using parent_type = OpAdd<Es...>;
	using parent_type::parent_type;
	using E = OpAdd<E0, Es...>;

	explicit OpExpression(OpAdd<Es...> const& rest) : parent_type(rest) {}
	explicit OpExpression(OpAdd<Es...>&& rest) noexcept : parent_type(std::move(rest)) {}

	auto operator()(iter_type n) const
	{
		return cast().eval(n);
	}

	template<typename EE>
	auto operator()(OpExpression<EE> const& e) const
	{
		return cast() * (*static_cast<EE const*>(&e));
	}

	symphas::internal::expression_iterator<E> begin() const
	{
		return symphas::internal::expression_iterator<E>(cast());
	}

	symphas::internal::expression_iterator<E> end(len_type len) const
	{
		return symphas::internal::expression_iterator<E>(cast(), len);
	}

	auto& cast() const
	{
		return *static_cast<E const*>(this);
	}

};

// ******************************************************************************************

//! Binary expression, the addition of two terms.
/*!
 * Binary addition between two expressions.
 * 
 * \tparam E1 The type of the left hand side expression.
 * \tparam E2 The type of the right hand side expression.
 */
template<typename E0, typename... Es>
struct OpAdd<E0, Es...> : OpExpression<OpAdd<E0, Es...>>
{
protected:
	using parent_type = OpExpression<OpAdd<E0, Es...>>;

public:

	//! Create the binary addition expression between two expressions.
	/*!
	 * Create the binary addition expression between two expressions.
	 * 
	 * \param a The expression on the left hand side of the addition operator.
	 * \param b The expression on the right hand side of the addition operator.
	 */
	OpAdd(E0 const& e, Es const&... es) : parent_type(es...), data{ e } {}
	OpAdd(E0 const& e, OpAdd<Es...> const& rest) : parent_type(rest), data{ e } {}
	OpAdd(OpAdd<Es...> const& rest, E0 const& e) : parent_type(rest), data{ e } {}

	//template<typename data_type = E0, typename nested_type = OpAdd<Es...>,
	//	typename = std::enable_if_t<(
	//		std::is_default_constructible<data_type>::value &&
	//		std::is_default_constructible<nested_type>::value), int>>
	OpAdd() : parent_type(), data{} {}

	inline auto eval(iter_type n = 0) const
	{
		return (data.eval(n) + parent_type::eval(n));
	}

	auto operator-() const;


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = data.print(out);
		return n + symphas::internal::binary_print(*this, out, std::make_index_sequence<sizeof...(Es)>{});
	}

	size_t print(char* out) const
	{
		size_t n = data.print(out);
		return n + symphas::internal::binary_print(*this, out + n, std::make_index_sequence<sizeof...(Es)>{});
	}

	size_t print_length() const
	{
		size_t n = data.print_length() + SYEX_BINARY_FMT_LEN;
		if (parent_type::coeff_sign())
		{
			n -= 1;
		}
		return n + parent_type::print_length();
	}

protected:

	bool coeff_sign() const
	{
		return (expr::coeff(data).eval() < 0);
	}

public:

#endif

	E0 data;

	//template<typename S = symphas::lib::type_at_index<sizeof...(Es), E0, Es...>, 
	//	typename T = add_result_t<expr::eval_type_t<E0>, expr::eval_type_t<Es>...>,
	//	typename std::enable_if_t<(expr::is_coeff<S> && expr::is_coeff<E0> && (expr::is_coeff<Es> && ...)), int> = 0>
	//explicit operator decltype(T{})() const
	//{
	//	return eval(0);
	//}

	template<typename cast_type>
	friend struct symphas::internal::cast_add;

};

template<typename E0, typename... Es>
OpAdd(E0, OpAdd<Es...>)->OpAdd<E0, Es...>;
template<typename E0, typename... Es>
OpAdd(E0, Es...)->OpAdd<E0, Es...>;

inline OpVoid expr::make_add()
{
	return OpVoid{};
}

inline OpVoid expr::make_add(OpVoid, OpVoid)
{
	return OpVoid{};
}

template<typename E>
auto expr::make_add(E&& a)
{
	return std::forward<E>(a);
}

template<typename E>
auto expr::make_add(OpVoid, E&& e)
{
	return std::forward<E>(e);
}

template<typename E>
auto expr::make_add(E&& e, OpVoid)
{
	return std::forward<E>(e);
}

template<typename... E0s, typename E0>
auto expr::make_add(OpAdd<E0s...> const& add, E0 const& e0)
{
	return OpAdd<E0, E0s...>(add, e0);
}

template<typename... E0s>
auto expr::make_add(OpAdd<E0s...> const& add, OpVoid)
{
	return add;
}

template<typename... E0s, typename E0>
auto expr::make_add(E0 const& e0, OpAdd<E0s...> const& add)
{
	return OpAdd<E0, E0s...>(add, e0);
}

template<typename... E0s>
auto expr::make_add(OpVoid, OpAdd<E0s...> const& add)
{
	return add;
}

template<typename A, typename B>
auto expr::make_add(OpExpression<A> const& a, OpExpression<B> const& b)
{
	return OpAdd<A, B>(*static_cast<A const*>(&a), *static_cast<B const*>(&b));
}

template<typename E1, typename E2, typename E3, typename... Es>
auto expr::make_add(E1&& a, E2&& b, E3&& c, Es&& ...es)
{
	return make_add(make_add(std::forward<E1>(a), std::forward<E2>(b)), std::forward<E3>(c), std::forward<Es>(es)...);
}


template<typename E1, typename E2>
auto operator+(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator-(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), -*static_cast<E2 const*>(&b));
}


namespace symphas::internal
{
	using expr::has_coeff;
	using expr::has_pmi_coeff;
	using expr::has_nmi_coeff;

#ifdef PRINTABLE_EQUATIONS

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, SYEX_MUL_SEP_B);
		n += static_cast<E2 const*>(&b)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_B);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, SYEX_MUL_SEP_AB);
		n += static_cast<E2 const*>(&b)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_AB);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_pmi_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, SYEX_MUL_SEP_B);
		n += static_cast<E2 const*>(&b)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_pmi_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_B);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_nmi_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_A);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, SYEX_MUL_SEP);
		n += static_cast<E2 const*>(&b)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_nmi_coeff<E1> && has_nmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_A);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, SYEX_MUL_SEP_OP);
		n += static_cast<E2 const*>(&b)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_OP);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_BB);
		return n;
	}


#endif

}

//! Binary expression, the multiplication of two terms.
template<typename E1, typename E2>
struct OpBinaryMul : OpExpression<OpBinaryMul<E1, E2>>
{
	OpBinaryMul(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	//template<typename AA = E1, typename BB = E2, 
	//	typename = std::enable_if_t<(std::is_default_constructible<AA>::value && std::is_default_constructible<BB>::value), int>>
	OpBinaryMul() : a{ E1{} }, b{ E2{} } {}

	inline auto eval(iter_type n) const
	{
		return a.eval(n) * b.eval(n);
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return symphas::internal::mul_print(out, a, b);
	}

	size_t print(char* out) const
	{
		return symphas::internal::mul_print(out, a, b);
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_MUL_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.

};

namespace expr
{
	//! Constructs the binary multiplication expression.
	/*!
	 * Directly constructs the binary multiplication expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2>
	auto make_mul(E1 const& a, E2 const& b)
	{
		return OpBinaryMul<E1, E2>(a, b);
	}

	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b);
}

template<typename E1, typename E2, 
	typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_mul(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return OpOperatorChain(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

//! Distributing an RHS expression between operands in addition expression.
template<typename... As, typename B2>
auto operator*(OpAdd<As...> const& a, OpOperator<B2> const& b)
{
	return OpOperatorChain(a, *static_cast<const B2*>(&b));
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank > 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::dot(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename A, typename B, typename E>
auto operator*(OpBinaryMul<A, B> const& a, OpOperator<E> const& b)
{
	return OpOperatorChain(a, *static_cast<const E*>(&b));
}

template<typename E1, typename E2>
auto operator*(OpOperator<E1> const& a, OpOperator<E2> const& b)
{
	return (*static_cast<const E1*>(&a)).operator*(*static_cast<const E2*>(&b));
}

//! Binary expression, the division of two terms.
template<typename E1, typename E2>
struct OpBinaryDiv : OpExpression<OpBinaryDiv<E1, E2>>
{
	OpBinaryDiv(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	//template<typename AA = E1, typename BB = E2,
	//	typename = std::enable_if_t<(std::is_default_constructible<AA>::value&& std::is_default_constructible<BB>::value), int>>
	OpBinaryDiv() : a{ E1{} }, b{ E2{} } {}

	inline auto eval(iter_type n) const
	{
		return a.eval(n) / b.eval(n);
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, SYEX_DIV_FMT_A);
		n += a.print(out);
		n += fprintf(out, SYEX_DIV_SEP);
		n += b.print(out);
		n += fprintf(out, SYEX_DIV_FMT_B);
		return n;

	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_DIV_FMT_A);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_DIV_SEP);
		n += b.print(out + n);
		n += sprintf(out + n, SYEX_DIV_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_DIV_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.
};

namespace symphas::internal
{

	template<typename>
	struct coeff_factor_numerator
	{
		static const size_t value = 1;
	};

	template<size_t N, size_t D>
	struct coeff_factor_numerator<OpFractionLiteral<N, D>>
	{
		static const size_t value = N;
	};

	template<typename>
	struct coeff_factor_denominator
	{
		static const size_t value = 1;
	};

	template<size_t N, size_t D>
	struct coeff_factor_denominator<OpFractionLiteral<N, D>>
	{
		static const size_t value = D;
	};

	template<typename E>
	struct coeff_factor_impl
	{
		using type = typename coeff_factor_impl<expr::coeff_t<E>>::type;
	};

	template<typename T, size_t... Ns>
	struct coeff_factor_impl<OpTensor<T, Ns...>>
	{
		using type = typename coeff_factor_impl<T>::type;
	};

	template<typename T>
	struct coeff_factor_impl<OpLiteral<T>>
	{
		using type = OpIdentity;
	};

	template<>
	struct coeff_factor_impl<OpNegIdentity>
	{
		using type = OpIdentity;
	};

	template<>
	struct coeff_factor_impl<OpIdentity>
	{
		using type = OpIdentity;
	};

	template<size_t N, size_t D>
	struct coeff_factor_impl<OpNegFractionLiteral<N, D>>
	{
		using type = OpFractionLiteral<N, D>;
	};

	template<size_t N, size_t D>
	struct coeff_factor_impl<OpFractionLiteral<N, D>>
	{
		using type = OpFractionLiteral<N, D>;
	};

	template<typename... As>
	struct coeff_factor_impl<OpAdd<As...>>
	{
		using type = decltype(expr::make_fraction<
			GCD_of<coeff_factor_numerator<typename coeff_factor_impl<As>::type>::value...>,
			GCD_of<coeff_factor_denominator<typename coeff_factor_impl<As>::type>::value...>>());
	};

	template<>
	struct coeff_factor_impl<OpVoid>
	{
		using type = OpVoid;
	};


	template<typename E>
	using coeff_factor = typename coeff_factor_impl<E>::type;
}

namespace expr
{

	//! Constructs the binary division expression.
	/*!
	 * Directly constructs the binary division expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2, size_t R1 = expr::eval_type<E1>::rank, size_t R2 = expr::eval_type<E2>::rank,
		typename std::enable_if_t<(R1 == 0 && R2 == 0), int> = 0>
	auto make_div(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		auto fr1 = symphas::internal::coeff_factor<E1>{};
		auto fr2 = symphas::internal::coeff_factor<E2>{};

		auto fr = (fr1 / fr2);
		constexpr int N = symphas::internal::coeff_factor_numerator<decltype(fr)>::value;
		constexpr int D = symphas::internal::coeff_factor_denominator<decltype(fr)>::value;

		auto _fr1 = expr::make_integer<N>() / fr1;
		auto _fr2 = expr::make_integer<D>() / fr2;


		return OpBinaryDiv(_fr1 * (*static_cast<const E1*>(&a)), _fr2 * (*static_cast<const E2*>(&b)));
	}

	template<typename E1, typename E2, size_t R1 = expr::eval_type<E1>::rank, size_t R2 = expr::eval_type<E2>::rank,
		typename std::enable_if_t<(R1 > 0 || R2 > 0), int> = 0>
	auto make_div(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return OpBinaryDiv((*static_cast<const E1*>(&a)), (*static_cast<const E2*>(&b)));
	}

	template<typename E2>
	auto make_div(OpVoid, OpExpression<E2> const& b)
	{
		return OpVoid{};
	}
}




namespace expr
{
	//! Get the expression that the OpMap applies to.
	template<typename G, typename V, typename E>
	auto const& get_enclosed_expression(OpMap<G, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto const& get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
	{
		return e.e;
	}
	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto const& get_enclosed_expression(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename V, typename E>
	auto const& get_enclosed_expression(OpMap<void, V, E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename G, typename Sp>
	decltype(auto) get_enclosed_expression(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename G, typename G0>
	auto const& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G>, SymbolicDerivative<G0>> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename E, typename G0>
	auto const& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename E, typename Sp>
	auto const& get_enclosed_expression(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpChain applies to.
	template<typename A1, typename A2, typename E>
	auto const& get_enclosed_expression(OpChain<A1, A2, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpCombination applies to.
	template<typename A1, typename A2, typename E>
	auto const& get_enclosed_expression(OpCombination<A1, A2, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpFunction applies to.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto const& get_enclosed_expression(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<auto f, typename V, typename E>
	auto const& get_enclosed_expression(OpFunctionApply<f, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename sub_t, typename E, typename... Ts>
	const auto& get_enclosed_expression(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		return e.f.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename G>
	decltype(auto) get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename E>
	auto const& get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return e.e;
	}


	//! Get the expression that the OpExponential applies to.
	template<typename V, typename E>
	auto const& get_enclosed_expression(OpExponential<V, E> const& e)
	{
		return e.e;
	}

	template<expr::exp_key_t X, typename V, typename E>
	auto const& get_enclosed_expression(OpPow<X, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename G, typename V, typename E>
	auto& get_enclosed_expression(OpMap<G, V, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto& get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto& get_enclosed_expression(OpMap<symphas::internal::STHC, OpIdentity, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename V, typename E>
	auto& get_enclosed_expression(OpMap<void, V, E>& e)
	{
		return *static_cast<E*>(&e);
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename G, typename Sp>
	decltype(auto) get_enclosed_expression(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename G, typename G0>
	auto& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G>, SymbolicDerivative<G0>>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename E, typename Sp>
	auto& get_enclosed_expression(OpDerivative<Dd, V, E, Sp>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename E, typename G0>
	auto& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpChain applies to.
	template<typename A1, typename A2, typename E>
	auto& get_enclosed_expression(OpChain<A1, A2, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpCombination applies to.
	template<typename A1, typename A2, typename E>
	auto& get_enclosed_expression(OpCombination<A1, A2, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto& get_enclosed_expression(OpFunction<V, E, F, Arg0, Args...>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<auto f, typename V, typename E>
	auto& get_enclosed_expression(OpFunctionApply<f, V, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename sub_t, typename E, typename... Ts>
	auto& get_enclosed_expression(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>& e)
	{
		return e.f.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename G>
	decltype(auto) get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename E>
	auto& get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, E>& e)
	{
		return e.e;
	}

	template<expr::exp_key_t X, typename V, typename E>
	auto& get_enclosed_expression(OpPow<X, V, E> &e)
	{
		return e.e;
	}

	//! Get the expression that the OpExponential applies to.
	template<typename V, typename E>
	auto& get_enclosed_expression(OpExponential<V, E>& e)
	{
		return e.e;
	}


	template<typename E>
	auto& get_result_data(OpExpression<E>& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename G, typename Sp>
	auto& get_result_data(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename E, typename Sp>
	auto& get_result_data(OpDerivative<Dd, V, E, Sp>& e)
	{
		return e.grid;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<size_t O, typename V, typename E, typename G0>
	auto& get_result_data(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>>&) = delete;


	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, typename E1, typename E2>
	auto& get_result_data(OpConvolution<V, E1, E2>& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename G>
	auto& get_result_data(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename E>
	auto& get_result_data(OpConvolution<V, GaussianSmoothing<D>, E>& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename E, typename Sp>
	auto const& get_result_data(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return e.grid;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename G, typename Sp>
	auto& get_result_data(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<size_t O, typename V, typename E, typename G0>
	auto const& get_result_data(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>> const&) = delete;

	//! Get the grid storing the underlying data of the OpExponential.
	template<typename V, typename E>
	auto const& get_result_data(OpExponential<V, E> const& e)
	{
		return e.data;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, typename E1, typename E2>
	auto const& get_result_data(OpConvolution<V, E1, E2> const& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename G>
	auto const& get_result_data(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename E>
	auto const& get_result_data(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return e.g0;
	}

}

namespace expr
{

	template<int I>
	constexpr auto sym_N = expr::make_integer<I>();
	template<size_t N0, size_t N1, size_t... Ns>
	constexpr auto sym_T = expr::make_tensor<N0, N1, Ns...>();
	template<size_t I, size_t N>
	constexpr auto sym_R = expr::make_row_vector<I, N>(OpIdentity{});
	template<size_t I, size_t N>
	constexpr auto sym_C = expr::make_column_vector<I, N>(OpIdentity{});
}

#undef SYEX_BINARY_FMT
#undef SYEX_BINARY_FMT_LEN
#undef SYEX_ADD_SEP
#undef SYEX_SUB_SEP
#undef SYEX_MUL_FMT
#undef SYEX_MUL_FMT_LATEX
#undef SYEX_MUL_FMT_LEN
#undef SYEX_MUL_FMT_LATEX_LEN
#undef SYEX_DIV_FMT
#undef SYEX_DIV_FMT_LATEX
#undef SYEX_DIV_FMT_LEN
#undef SYEX_DIV_FMT_LATEX_LEN



