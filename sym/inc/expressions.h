
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

#include <tuple>
#include <iostream>


#include "expressioniterator.h"



// ******************************************************************************************

//! \cond


#define SYEX_BINARY_FMT "%s %s %s"
#define SYEX_ADD_SEP " + "
#define SYEX_SUB_SEP " - "
#define SYEX_BINARY_FMT_LEN (sizeof(SYEX_ADD_SEP) / sizeof(char) - 1)



#ifdef LATEX_PLOT

#define SYEX_MUL_SEP_OP "\\cdot "
#define SYEX_DIV_SEP_OP "}{"
#define SYEX_NLV_SEP_OP " "

#define SYEX_MUL_FMT_AA "\\left["
#define SYEX_MUL_FMT_AB "\\left("
#define SYEX_MUL_FMT_BA "\\right)"
#define SYEX_MUL_FMT_BB "\\right]"

#define SYEX_DIV_SEP SYEX_DIV_SEP_OP
#define SYEX_DIV_FMT_A "\\frac{"
#define SYEX_DIV_FMT_B "}" 

#else

#define SYEX_MUL_SEP_OP "*"
#define SYEX_DIV_SEP_OP "/"
#define SYEX_NLV_SEP_OP "*"

#define SYEX_MUL_FMT_AA "["
#define SYEX_MUL_FMT_AB "("
#define SYEX_MUL_FMT_BA ")"
#define SYEX_MUL_FMT_BB "]"


#define SYEX_DIV_SEP  ") " SYEX_DIV_SEP_OP " ("
#define SYEX_DIV_FMT_A SYEX_MUL_FMT_A
#define SYEX_DIV_FMT_B SYEX_MUL_FMT_B 

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


#define SYEX_LAMBDA_FUNC_FMT_A "%s("
#define SYEX_LAMBDA_FUNC_FMT_B ")"
#define SYEX_LAMBDA_FUNC_FMT SYEX_DIV_FMT_A SYEX_LAMBDA_FUNC_FMT_A "%s" SYEX_LAMBDA_FUNC_FMT_B
#define SYEX_LAMBDA_FUNC_FMT_LEN (STR_ARR_LEN(SYEX_LAMBDA_FUNC_FMT_A SYEX_LAMBDA_FUNC_FMT_B) - 3)

#ifdef LATEX_PLOT
#define SYEX_POW_SEP_A "^{"
#define SYEX_POW_SEP_B "}"
#else
#define SYEX_POW_SEP_A "^"
#define SYEX_POW_SEP_B
#endif

//! \endcond


//! The display precision of values which appear in expressions.
/*!
 * Defines the number of digits which appear in floating point numbers
 * for values from expressions that are printed to the screen.
 */
#define EXPR_VALUE_DISPLAY_PRECISION 3


//! Base expression object which is inherited from with the CRTP technique.
/*
 * applying Expression Templates to create the expression tree for the
 * evaluation of the equations of motion
 */
template<typename E>
struct OpExpression
{
	//! Return the value of this expression at the given index.
	/*!
	 * Evaluate the expression, which typically refers to a list of data points
	 * that are indexed sequentially.
	 * 
	 * \param n The index at which to evaluate the data.
	 */
	inline auto eval(iter_type n) const
	{
		return cast().eval(n);
	}

	//! Return the value of this expression at the given index.
	/*!
	 * Evaluate the expression, which typically refers to a list of data points
	 * that are indexed sequentially.
	 *
	 * \param n The index at which to evaluate the data.
	 */
	auto operator()(iter_type n) const
	{
		return cast().eval(n);
	}


#ifdef PRINTABLE_EQUATIONS

	//! Print the string representation of this expression to the file.
	/*!
	 * The string representation of this expression is printed to the given
	 * file. The string is assumed to already have enough memory allocated.
	 *
	 * \param out The file to which the expression is printed.
	 */
	size_t print(FILE *out) const 
	{
		return cast().print(out);
	}

	//! Print the string representation of this expression to the string.
	/*!
	 * The string representation of this expression is printed to the given
	 * string. The string is assumed to already have enough memory allocated.
	 * 
	 * \param out The string array into which to write.
	 */
	size_t print(char *out) const
	{
		return cast().print(out);
	}

	//! Returns the number of characters in the string representation.
	/*!
	 * Returns the number of characters in the string that represents this
	 * expression. This value does not include a terminating character.
	 */
	size_t print_length() const
	{
		return cast().print_length();
	}

#else

	size_t print(...) const
	{
		return 0;
	}

#endif

	auto& cast() const
	{
		return *static_cast<E const*>(this);
	}

	//! Return an iterator the beginning of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the beginning of the data, used when evaluating
	 * the expression.
	 */
	symphas::internal::expression_iterator<E> begin() const 
	{ 
		return symphas::internal::expression_iterator<E>(cast());
	}


	//! Return an iterator the end of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the end of the data, used when evaluating
	 * the expression. The end point has to be provided, as the length
	 * of the data is not known directly by the expression.
	 * 
	 * \param len The end point of the data, for the end iterator to
	 * point to.
	 */
	symphas::internal::expression_iterator<E> end(len_type len) const 
	{ 
		return symphas::internal::expression_iterator<E>(cast(), len);
	}
};

//! Contains all elements constituting the symbolic algebra functionality.
/*!
 * Defines elements which support the symbolic algebra functionality, such
 * as creating new expression terms, manipulating and transforming expressions,
 * and printing expressions.
 */
namespace expr
{

	//! Constructs a constant of the given value.
	/*!
	 * Constructs a constant of the given value.
	 * 
	 * \param v The value to give to the literal.
	 */
	template<typename T>
	auto make_literal(T v);
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
	OpLiteral(OpIdentity);
	OpLiteral(OpNegIdentity);
	OpLiteral(OpVoid);

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
	auto operator*(OpLiteral<S> const other) const
	{
		return expr::make_literal(value * other.value);
	}

	auto operator*(OpLiteral<T> const other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S> const a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}
	auto operator-() const
	{
		return expr::make_literal(-value);
	}

	template<typename S>
	auto operator+(OpLiteral<S> const a) const
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

// ******************************************************************************************

//! Initialize an expression for a constant.
/*!
 * Create a literal expression using the given value.
 */
template<typename T>
inline auto expr::make_literal(T v)
{
	return OpLiteral{ v };
}

// ******************************************************************************************

//! Additive identity.
/*!
 * Expression representing the additive identity, a special type of operator 
 * used primarily in the simplification of expressions.
 */
struct OpVoid : OpExpression<OpVoid>
{
	constexpr int eval(iter_type = 0) const
	{
		return 0;
	}

	size_t print(FILE* out) const
	{
		return fprintf(out, "0");
	}

	size_t print(char* out) const
	{
		return sprintf(out, "0");
	}

	auto operator-() const
	{
		return OpVoid{};
	}

	size_t print_length() const
	{
		return 1;
	}
};

// ******************************************************************************************


//! Multiplicative identity.
/*!
 * Expression representing the multiplicative identity, a special type of 
 * expression used primarily in the simplification of expressions. It is also
 * used as the coefficient to other expressions and can be substituted directly
 * as a value.
 */
struct OpIdentity : OpExpression<OpIdentity>
{
	constexpr auto eval(iter_type = 0) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return fprintf(out, "1");
	}

	size_t print(char *out) const
	{
		return sprintf(out, "1");
	}

	size_t print_length() const
	{
		return 1;
	}

#endif

	auto operator-() const;

};

//! Negative of the multiplicative identity.
/*!
 * Expression representing the negative of the multiplicative identity, a special type of 
 * expression used primarily in the simplification of expressions. It is also
 * used as the coefficient to other expressions and can be substituted directly
 * as a value.
 */
struct OpNegIdentity : OpExpression<OpNegIdentity>
{
	constexpr auto eval(iter_type = 0) const
	{
		return -OpIdentity{}.eval();
	}

	size_t print(FILE* out) const
	{
		return fprintf(out, "-1");
	}

	size_t print(char* out) const
	{
		return sprintf(out, "-1");
	}


	size_t print_length() const
	{
		return 2;
	}

	auto operator-() const;

};

inline auto OpIdentity::operator-() const
{
	return OpNegIdentity{};
}

inline auto OpNegIdentity::operator-() const
{
	return OpIdentity{};
}


//! Specialized version which returns multiplicative identity.
template<>
inline auto expr::make_literal<OpIdentity>(OpIdentity)
{
	return OpIdentity{};
}

//! Specialized version which returns negative multiplicative identity.
template<>
inline auto expr::make_literal<OpNegIdentity>(OpNegIdentity)
{
	return OpNegIdentity{};
}

//! Specialized version which returns additive identity.
template<>
inline auto expr::make_literal<OpVoid>(OpVoid)
{
	return OpVoid{};
}


template<typename T>
OpLiteral<T>::OpLiteral(OpIdentity) : value{ symphas::lib::get_identity<T>() } {}
template<typename T>
OpLiteral<T>::OpLiteral(OpNegIdentity) : value{ -symphas::lib::get_identity<T>() } {}
template<typename T>
OpLiteral<T>::OpLiteral(OpVoid) : value{ 0 } {}




//! Return type of evaluating an expression identity.
/*!
 * Returns the type of the multiplicative identity. The same
 * type would be returned for the negative of the multiplicative identity.
 */
using identity_eval_t = decltype(std::declval<OpIdentity>().eval());



namespace symphas::internal
{
	/* testing to check if the type has the "value" attribute, which
	* is the coefficient to an expression.
	 */

	template <typename E> static char test_coeff_attribute(decltype(&E::value));
	template <typename E> static long test_coeff_attribute(...);

	template <typename E> static auto test_coeff_neg(decltype(E::value))
		->std::conditional_t<std::is_same<decltype(E::value), OpNegIdentity>::value, char, long>;
	template <typename E> static long test_coeff_neg(...);

	template <typename E> static auto test_coeff_pos(decltype(E::value))
		->std::conditional_t<std::is_same<decltype(E::value), OpIdentity>::value, char, long>;
	template <typename E> static long test_coeff_pos(...);
}

namespace expr
{
	//! True if the given type is expression identity.
	/*!
	 * This value is true if the given type is any of OpIdentity, 
	 * OpNegIdentity or OpVoid, meaning all the expression identities.
	 * 
	 * \tparam E The given type that is checked.
	 */
	template<typename E>
	constexpr bool is_identity =
		(std::is_same<E, OpIdentity>::value ||
			std::is_same<E, OpNegIdentity>::value ||
			std::is_same<E, OpVoid>::value);


	//! True when the types are compatible with identity types.
	/*!
	 * This is true only when one of the given types is an identity and
	 * the other is a type which is convertible to the evaluate type of
	 * an identity.
	 * 
	 * \tparam The first given type, either an identity or a value type.
	 * \tparam The second given type, either an identity or a value type.
	 */
	template<typename E1, typename E2>
	constexpr bool identity_comparator_enabled =
		(std::is_same<E1, OpIdentity>::value && std::is_convertible<E2, identity_eval_t>::value ||
			std::is_same<E2, OpIdentity>::value && std::is_convertible<E1, identity_eval_t>::value ||
			std::is_same<E1, OpNegIdentity>::value && std::is_convertible<E2, identity_eval_t>::value ||
			std::is_same<E2, OpNegIdentity>::value && std::is_convertible<E1, identity_eval_t>::value);


	using symphas::internal::test_coeff_attribute;
	using symphas::internal::test_coeff_neg;
	using symphas::internal::test_coeff_pos;

	//! Tests whether the given expression has a coefficient. 
	/*!
	 * Tests whether the given expression has a member called `value`. This
	 * member is always the coefficient of the expression. For example,
	 * during printing, a different format may be used depending on the
	 * sign of the coefficient.
	 * 
	 * \tparam E The expression type to check the coefficient.
	 */
	template<typename E>
	constexpr bool has_coeff = sizeof(test_coeff_attribute<E>(0)) == sizeof(char);

	//! Tests if the coefficient of the expression is minus of the ::OpIdentity.
	/*!
	 * Tests whether the given expression has a member called `value`. This
	 * member is always the coefficient of the expression. If there is the
	 * `value` member, returns true when `value` is an ::OpNegIdentity.
	 * (nmi means negative multiplicative identity).
	 *
	 * See expr::has_coeff.
	 * 
	 * \tparam E The expression type to check the coefficient.
	 */
	template<typename E>
	constexpr bool has_nmi_coeff = sizeof(test_coeff_neg<E>(0)) == sizeof(char);

	//! Tests if the coefficient of the expression is ::OpIdentity.
	/*!
	 * Tests whether the given expression has a member called `value`. This
	 * member is always the coefficient of the expression. If there is the
	 * `value` member, returns true when `value` is an ::OpIdentity.
	 * (pmi means positive multiplicative identity).
	 *
	 * See expr::has_coeff.
	 *
	 * \tparam E The expression type to check the coefficient.
	 */
	template<typename E>
	constexpr bool has_pmi_coeff = sizeof(test_coeff_pos<E>(0)) == sizeof(char);


	//! Tests if the coefficient of the expression is not an identity.
	/*!
	 * Tests whether the given expression has a member called `value`. This
	 * member is always the coefficient of the expression. If there is the
	 * `value` member, returns true when `value` is not ::OpIdentity and
	 * is not ::OpNegIdentity. This means that the coefficient
	 * has a particular value.
	 *
	 * See expr::has_coeff.
	 *
	 * \tparam E The expression type to check the coefficient.
	 */
	template<typename E>
	constexpr bool has_val_coeff = has_coeff<E> && !(has_nmi_coeff<E> || has_pmi_coeff<E>);
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::is_identity<E1> && expr::is_identity<E2>), int> = 0>
constexpr bool operator==(E1 const& a, E2 const& b)
{
	return a.eval() == b.eval();
}

template<typename A, typename B,
	typename std::enable_if_t<(expr::is_identity<A> && expr::is_identity<B>), int> = 0>
constexpr bool operator<(A const& a, B const& b)
{
	return a.eval() < b.eval();
}

template<typename E1, typename E2, 
	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
inline bool operator==(E1 const& a, E2 const& b)
{
	return OpLiteral<identity_eval_t>(a).eval() == OpLiteral<identity_eval_t>(b).eval();
}

template<typename E1, typename E2,
	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
inline bool operator<(E1 const& a, E2 const& b)
{
	return OpLiteral<identity_eval_t>(a).eval() < OpLiteral<identity_eval_t>(b).eval();
}


template<typename E1, typename E2,
	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
inline bool operator!=(E1 const& a, E2 const& b)
{
	return !(a == b);
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
inline bool operator<=(E1 const& a, E2 const& b)
{
	return (a < b || a == b);
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
inline bool operator>(E1 const& a, E2 const& b)
{
	return !(a < b || a == b);
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
inline bool operator>=(E1 const& a, E2 const& b)
{
	return !(a < b);
}

namespace expr
{
	//! Constructs the binary addition expression.
	/*!
	 * Directly constructs the binary addition expression between two
	 * expressions without applying any rules.
	 * 
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2>
	auto make_add(OpExpression<E1> const& a, OpExpression<E2> const& b);

	//! Constructs the binary subtraction expression.
	/*!
	 * Directly constructs the binary subtraction expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2>
	auto make_sub(OpExpression<E1> const& a, OpExpression<E2> const& b);
}

namespace symphas::internal
{
	using expr::has_coeff;

#ifdef PRINTABLE_EQUATIONS

	inline std::tuple<size_t, size_t> print_sep(char* out, const char* sep)
	{
		size_t n = sprintf(out, "%s", sep);
		return { n, n };
	}

	inline std::tuple<size_t, size_t> print_sep(FILE* out, const char* sep)
	{
		return { fprintf(out, "%s", sep), 0 };
	}

	//! Pretty print for compound expressions where a double sign may appear.
	/*!
	 * In order for a subtraction to be printed correctly if an OpBinaryAdd
	 * or OpBinarySub is itself subtracted from another expression, the minus
	 * sign has to be distributed inside.
	 */
	template<typename E1, typename E2, typename os_type, typename std::enable_if_t<!has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, E2> const& b, os_type* out)
	{
		auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
		n += b.b.print(out + offset);
		return n;
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename os_type, typename std::enable_if_t<!has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinaryAdd<E1, E2> const& b, os_type* out)
	{
		auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
		n += b.b.print(out + offset);
		return n;
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinaryAdd<E1, OpBinarySub<E2, E3>> const& b, os_type* out)
	{
		if (b.b.a.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += expr::make_sub(-b.b.a, b.b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += expr::make_sub(b.b.a, b.b.b).print(out + offset);
			return n;
		}
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinaryAdd<E1, OpBinaryAdd<E2, E3>> const& b, os_type* out)
	{
		if (b.b.a.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += expr::make_add(-b.b.a, b.b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += expr::make_add(b.b.a, b.b.b).print(out + offset);
			return n;
		}
	}


	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, OpBinarySub<E2, E3>> const& b, os_type* out)
	{
		if (b.b.a.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += expr::make_add(-b.b.a, b.b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += expr::make_add(b.b.a, b.b.b).print(out + offset);
			return n;
		}
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<!has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, OpBinarySub<E2, E3>> const& b, os_type* out)
	{
		auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
		n += expr::make_add(b.b.a, b.b.b).print(out + offset);
		return n;
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, OpBinaryAdd<E2, E3>> const& b, os_type* out)
	{
		if (b.b.a.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += expr::make_sub(-b.b.a, b.b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += expr::make_sub(b.b.a, b.b.b).print(out + offset);
			return n;
		}
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename E3, typename os_type, typename std::enable_if_t<!has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, OpBinaryAdd<E2, E3>> const& b, os_type* out)
	{
		auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
		n += expr::make_sub(b.b.a, b.b.b).print(out + offset);
		return n;
	}

	

	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinarySub<E1, E2> const& b, os_type* out)
	{
		if (b.b.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += (-b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += b.b.print(out + offset);
			return n;
		}
	}


	//! Specialization based on symphas::internal::binary_print().
	template<typename E1, typename E2, typename os_type, typename std::enable_if_t<has_coeff<E2>, int> = 0>
	size_t binary_print(OpBinaryAdd<E1, E2> const& b, os_type* out)
	{
		if (b.b.value < 0)
		{
			auto&& [n, offset] = print_sep(out, SYEX_SUB_SEP);
			n += (-b.b).print(out + offset);
			return n;
		}
		else
		{
			auto&& [n, offset] = print_sep(out, SYEX_ADD_SEP);
			n += b.b.print(out + offset);
			return n;
		}
	}

#endif
}



// ******************************************************************************************

//! Binary expression, the addition of two terms.
/*!
 * Binary addition between two expressions.
 * 
 * \tparam E1 The type of the left hand side expression.
 * \tparam E2 The type of the right hand side expression.
 */
template<typename E1, typename E2>
struct OpBinaryAdd : OpExpression<OpBinaryAdd<E1, E2>>
{
	//! Create the binary addition expression between two expressions.
	/*!
	 * Create the binary addition expression between two expressions.
	 * 
	 * \param a The expression on the left hand side of the addition operator.
	 * \param b The expression on the right hand side of the addition operator.
	 */
	OpBinaryAdd(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	inline auto eval(iter_type n) const
	{
		return a.eval(n) + b.eval(n);
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = a.print(out);
		n += symphas::internal::binary_print(*this, out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = a.print(out);
		n += symphas::internal::binary_print(*this, out + n);
		return n;
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_BINARY_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.
};

template<typename E1, typename E2>
auto expr::make_add(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return OpBinaryAdd<E1, E2>(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

//! Binary expression, the subtraction of two terms.
/*!
 * Binary subtraction between two expressions.
 *
 * \tparam E1 The type of the left hand side expression.
 * \tparam E2 The type of the right hand side expression.
 */
template<typename E1, typename E2>
struct OpBinarySub : OpExpression<OpBinarySub<E1, E2>>
{
	//! Create the binary subtraction expression between two expressions.
	/*!
	 * Create the binary subtraction expression between two expressions.
	 *
	 * \param a The expression on the left hand side of the addition operator.
	 * \param b The expression on the right hand side of the addition operator.
	 */
	OpBinarySub(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	inline auto eval(iter_type n) const
	{
		return a.eval(n) - b.eval(n);
	}

	auto operator-() const;


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = a.print(out);
		n += symphas::internal::binary_print(*this, out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = a.print(out);
		n += symphas::internal::binary_print(*this, out + n);
		return n;
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_BINARY_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.
};

template<typename E1, typename E2>
auto expr::make_sub(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return OpBinarySub<E1, E2>(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}



template<typename E1, typename E2>
auto operator-(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_sub(a, b);
}

template<typename E1, typename E2>
auto operator+(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(a, b);
}


namespace symphas::internal
{
	using expr::has_coeff;
	using expr::has_pmi_coeff;

#ifdef PRINTABLE_EQUATIONS

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && !has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_A);
		n += a.print(out);
		n += fprintf(out, SYEX_MUL_SEP);
		n += b.print(out);
		n += fprintf(out, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && !has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_A);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP);
		n += b.print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_A);
		n += a.print(out);
		n += fprintf(out, SYEX_MUL_SEP_A);
		n += b.print(out);
		n += fprintf(out, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(!has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_A);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_A);
		n += b.print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && !has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += a.print(out);
		n += fprintf(out, SYEX_MUL_SEP_B);
		n += b.print(out);
		n += fprintf(out, SYEX_MUL_FMT_B);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && !has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_B);
		n += b.print(out + n);
		n += sprintf(out + n, SYEX_MUL_FMT_B);
		return n;
	}


	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AA);
		n += a.print(out);
		n += fprintf(out, SYEX_MUL_SEP_OP);
		n += b.print(out);
		n += fprintf(out, SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename E1, typename E2, std::enable_if_t<(has_coeff<E1> && has_pmi_coeff<E2>), int> = 0>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AA);
		n += a.print(out + n);
		n += sprintf(out + n, SYEX_MUL_SEP_OP);
		n += b.print(out + n);
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
	auto make_mul(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return OpBinaryMul<E1, E2>(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
	}
}

template<typename E1, typename E2>
auto operator*(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_mul(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}


 //! Binary expression, the division of two terms.
template<typename E1, typename E2>
struct OpBinaryDiv : OpExpression<OpBinaryDiv<E1, E2>>
{
	OpBinaryDiv(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

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
	template<typename E1, typename E2>
	auto make_div(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return OpBinaryDiv<E1, E2>(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
	}
}


namespace expr
{
	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 * 
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename E>
	auto pow(OpExpression<E> const& e)
	{
		if constexpr (N == 1)
		{
			return *static_cast<E const*>(&e);
		}
		else
		{
			return *static_cast<E const*>(&e)
				* pow<N - 1>(*static_cast<E const*>(&e));
		}
	}
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



