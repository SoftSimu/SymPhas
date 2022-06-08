
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
 * PURPOSE: Defines elements which assist in printing expressions such as
 * derivatives and functions.
 *
 * ***************************************************************************
 */

#pragma once


#include <stdarg.h>

#include "expressions.h"
#include "expressionproperties.h"


template<typename Sp>
struct Solver;

#ifdef PRINTABLE_EQUATIONS

//! \cond

#ifdef LATEX_PLOT
#define SYEX_DERIV_APPLIED_EXPR_FMT_A "\\left("
#define SYEX_DERIV_APPLIED_EXPR_FMT_B "\\right)"
#else
#define SYEX_DERIV_APPLIED_EXPR_FMT_A "("
#define SYEX_DERIV_APPLIED_EXPR_FMT_B ")"
#endif
#define SYEX_DERIV_APPLIED_EXPR_FMT SYEX_DERIV_APPLIED_EXPR_FMT_A "%s" SYEX_DERIV_APPLIED_EXPR_FMT_B


#ifdef LATEX_PLOT
//! Display string of the derivative.
#define SYEX_DERIV_STR(VALUE) "\\frac{d" VALUE "}{d\\vec{x}" VALUE "}"
#define SYEX_DERIV_STR_1 SYEX_DERIV_STR("")
#else
//! Display string of the derivative.
#define SYEX_DERIV_STR(VALUE) "(d" VALUE "/dx" VALUE ")"
#define SYEX_DERIV_STR_1 SYEX_DERIV_STR("")
#endif

#define SYEX_DERIV_STR_SUB(ORDER) SYEX_DERIV_STR("%zd"), ORDER, ORDER
#define SYEX_DERIV_STR_LEN(ORDER) (STR_ARR_LEN(SYEX_DERIV_STR("")) + 2 * symphas::lib::num_digits<ORDER>() - 1)

//! Character used in place of an unknown identifier or derivative.
#define SYEX_UNKNOWN_FUNCTION_TOKEN "?"

//! The separation string used to separate a description from an expression.


#ifdef LATEX_PLOT
#define SYEX_EQN_SEP "\\rightarrow"
#else
#define SYEX_EQN_SEP "->"
#endif

/* Definitions concerning printing new functions
 */


#define SYEX_PRINT_DERIV_STRUCT(ORDER, FMT) \
	template<>																		\
	struct print_deriv<ORDER>														\
	{																				\
		static size_t print(FILE* out) { return fprintf(out, FMT); }				\
		static size_t print(char* out) { return sprintf(out, FMT); }				\
		static constexpr size_t print_length() { return STR_ARR_LEN(FMT) - 1; }		\
	};


//! Generates code which will associate a function with a display name.
#define SYEX_PRINT_FUNCTION_STRUCT_IMPL(F, FMT_STR) \
using f_arg_ ## F = typename expr::eval_type<E>::type; \
static size_t print(wrap_f<&symphas::math:: F <typename expr::eval_type<E>::type>>, FILE* out, const char* str) { return fprintf(out, FMT_STR, str); } \
static size_t print(wrap_f<&symphas::math:: F <typename expr::eval_type<E>::type>>, char* out, const char* str) { return sprintf(out, FMT_STR, str); } \
static constexpr size_t print_length(wrap_f<&symphas::math:: F <typename expr::eval_type<E>::type>>) { return sizeof(FMT_STR) / sizeof(char) - 3; }

#ifdef LATEX_PLOT
#define SYEX_DEFINE_FUNCTION_STRUCT_DELEGATE(F, FMT_STR, _) SYEX_PRINT_FUNCTION_STRUCT_IMPL(F, FMT_STR)
#else
#define SYEX_DEFINE_FUNCTION_STRUCT_DELEGATE(F, _, FMT_STR) SYEX_PRINT_FUNCTION_STRUCT_IMPL(F, FMT_STR)
#endif

//! Generates code which will associate a function with a display name.
#define SYEX_PRINT_FUNCTION_STRUCT(F, FMT_STR_LATEX, FMT_STR) SYEX_DEFINE_FUNCTION_STRUCT_DELEGATE(F, FMT_STR_LATEX, FMT_STR)


/* Definitions concerning the coefficients.
 */

#ifdef LATEX_PLOT
#define SYEX_COEFF_PRINT_FMT_SEP " "
#else
#define SYEX_COEFF_PRINT_FMT_SEP "."
#endif
//! Appearance of a coefficient in front of another string.
#define SYEX_COEFF_PRINT_FMT "%s" SYEX_COEFF_PRINT_FMT_SEP "%s"
#define SYEX_COEFF_SEP_LEN (sizeof(SYEX_COEFF_PRINT_FMT_SEP) / sizeof(char) - 1)

#define SYEX_ZERO_COEFF "0"

//! \endcond


// ************************************************************************************

namespace symphas::internal
{

	//! Generates an expression display string that use functions.
	/*!
	 * Uses expression template to avoid restating
	 * the print function for each specialization. Generates a display string
	 * that is used when printing the expression. 
	 */
	template<typename... Ts>
	struct print_f_op;


	//! Generates an expression display string that use functions.
	/*!
	 * See symphas::internal::print_f_op.
	 */
	template<typename E>
	struct print_f_op<E>
	{
	protected:

		static inline const char default_fmt[] = SYEX_UNKNOWN_FUNCTION_TOKEN SYEX_DERIV_APPLIED_EXPR_FMT;

		struct wrap_base {};
		template<auto f>
		struct wrap_f : wrap_base {};

		static size_t print(wrap_base, FILE* out) { return fprintf(out, default_fmt); }
		static size_t print(wrap_base, char* out) { return sprintf(out, default_fmt); }
		static constexpr size_t print_length(wrap_base, char* out, const char* str) { return sizeof(default_fmt) / sizeof(char) - 3; }

		SYEX_PRINT_FUNCTION_STRUCT(conj, "\\overline{%s}", "conj" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(modulus, "\\left|{%s}\\right|", "|%s|");
		SYEX_PRINT_FUNCTION_STRUCT(real, "\\Re" SYEX_DERIV_APPLIED_EXPR_FMT, "Re" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(imag, "\\Im" SYEX_DERIV_APPLIED_EXPR_FMT, "Im" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(cos, "\\cos" SYEX_DERIV_APPLIED_EXPR_FMT, "cos" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(sin, "\\sin" SYEX_DERIV_APPLIED_EXPR_FMT, "sin" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(tan, "\\tan" SYEX_DERIV_APPLIED_EXPR_FMT, "tan" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(cosh, "\\cosh" SYEX_DERIV_APPLIED_EXPR_FMT, "cosh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(sinh, "\\sinh" SYEX_DERIV_APPLIED_EXPR_FMT, "sinh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(tanh, "\\tanh" SYEX_DERIV_APPLIED_EXPR_FMT, "tanh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(acos, "\\acos" SYEX_DERIV_APPLIED_EXPR_FMT, "acos" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(asin, "\\asin" SYEX_DERIV_APPLIED_EXPR_FMT, "asin" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(atan, "\\atan" SYEX_DERIV_APPLIED_EXPR_FMT, "atan" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(acosh, "\\acosh" SYEX_DERIV_APPLIED_EXPR_FMT, "acosh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(asinh, "\\asinh" SYEX_DERIV_APPLIED_EXPR_FMT, "asinh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(atanh, "\\atanh" SYEX_DERIV_APPLIED_EXPR_FMT, "atanh" SYEX_DERIV_APPLIED_EXPR_FMT);
		SYEX_PRINT_FUNCTION_STRUCT(sqrt, "\\sqrt" SYEX_DERIV_APPLIED_EXPR_FMT, "sqrt" SYEX_DERIV_APPLIED_EXPR_FMT);

	public:

		template<auto f>
		static size_t print(FILE* out, const char* str)
		{
			return print(wrap_f<f>{}, out, str);
		}

		template<auto f>
		static size_t print(char* out, const char* str)
		{
			return print(wrap_f<f>{}, out, str);
		}


		template<auto f>
		static constexpr size_t print_length()
		{
			return print_length(wrap_f<f>{});
		}
	};

	namespace deriv_formats
	{

#ifdef LATEX_PLOT
		static inline const char laplacian_fmt[] = "\\nabla^2";
		static inline const char bilaplacian_fmt[] = "\\nabla^4";
		static inline const char gradlaplacian_fmt[] = "\\nabla\\nabla^2";
		static inline const char gradient_fmt[] = "\\nabla";

#else
		static inline const char laplacian_fmt[] = SYEX_DERIV_STR("2");
		static inline const char bilaplacian_fmt[] = SYEX_DERIV_STR("4");
		static inline const char gradlaplacian_fmt[] = SYEX_DERIV_STR_1 SYEX_DERIV_STR("2");
		static inline const char gradient_fmt[] = SYEX_DERIV_STR_1;

#endif
	}




	template<size_t O>
	struct print_deriv
	{

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output file, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		static size_t print(FILE* out)
		{
			return fprintf(out, SYEX_DERIV_STR_SUB(O));
		}

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output string, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		static size_t print(char* out)
		{
			return sprintf(out, SYEX_DERIV_STR_SUB(O));
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 *
		 * \tparam f The function which applies the derivative and from
		 * which to infer the order of the derivative. This function must
		 * be a member of the given solver type.
		 */
		static size_t print_length()
		{
			return SYEX_DERIV_STR_LEN(O);
		}
	};

	SYEX_PRINT_DERIV_STRUCT(1, deriv_formats::gradient_fmt);
	SYEX_PRINT_DERIV_STRUCT(2, deriv_formats::laplacian_fmt);
	SYEX_PRINT_DERIV_STRUCT(3, deriv_formats::gradlaplacian_fmt);
	SYEX_PRINT_DERIV_STRUCT(4, deriv_formats::bilaplacian_fmt);


	template<typename Sp, typename G>
	struct print_f_op<Sp, G>
	{
	protected:

		template<auto f>
		static constexpr size_t get_order()
		{
			return Sp::template order_of<f, G>;
		}

		template<size_t O>
		struct wrap_base {};
		template<auto f>
		struct wrap_f : wrap_base<get_order<f>()> {};

		using wrap_f_laplacian = wrap_f<&Solver<Sp>::template laplacian<G>>;
		using wrap_f_bilaplacian = wrap_f<&Solver<Sp>::template bilaplacian<G>>;
		using wrap_f_gradlaplacian = wrap_f<&Solver<Sp>::template gradlaplacian<G>>;
		using wrap_f_gradient = wrap_f<&Solver<Sp>::template gradient<G>>;

		template<size_t O, typename std::enable_if_t<(O > DERIV_MAX_ORDER), int> = 0>
		static size_t print(wrap_base<O>, FILE* out)
		{
			return fprintf(out, SYEX_UNKNOWN_FUNCTION_TOKEN);
		}

		template<size_t O, typename std::enable_if_t<(O > DERIV_MAX_ORDER), int> = 0>
		static size_t print(wrap_base<O>, char* out)
		{
			return sprintf(out, SYEX_UNKNOWN_FUNCTION_TOKEN);
		}


		template<size_t O, typename std::enable_if_t<(O <= DERIV_MAX_ORDER), int> = 0>
		static size_t print(wrap_base<O>, FILE* out)
		{
			return print_deriv<O>::print(out);
		}

		template<size_t O, typename std::enable_if_t<(O <= DERIV_MAX_ORDER), int> = 0>
		static size_t print(wrap_base<O>, char* out)
		{
			return print_deriv<O>::print(out);
		}


		template<size_t O, typename std::enable_if_t<(O > DERIV_MAX_ORDER), int> = 0>
		static constexpr size_t print_length(wrap_base<O>)
		{
			return STR_ARR_LEN(SYEX_UNKNOWN_FUNCTION_TOKEN) - 1;
		}

		template<size_t O, typename std::enable_if_t<(O <= DERIV_MAX_ORDER), int> = 0>
		static constexpr size_t print_length(wrap_base<O>)
		{
			return print_deriv<O>::print_length();
		}


	public:
		

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output file, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 * 
		 * \param out The file to which the derivative is printed.
		 * \param str The string representing the expression the derivative
		 * is applied to.
		 * 
		 * \tparam f The function which applies the derivative and from
		 * which to infer the order of the derivative. This function must
		 * be a member of the given solver type.
		 */
		template<auto f>
		static size_t print(FILE* out)
		{
			return print(wrap_f<f>{}, out);
		}

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output string, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param str The string representing the expression the derivative
		 * is applied to.
		 *
		 * \tparam f The function which applies the derivative and from
		 * which to infer the order of the derivative. This function must
		 * be a member of the given solver type.
		 */
		template<auto f>
		static size_t print(char* out)
		{
			return print(wrap_f<f>{}, out);
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 *
		 * \tparam f The function which applies the derivative and from
		 * which to infer the order of the derivative. This function must
		 * be a member of the given solver type.
		 */
		template<auto f>
		static constexpr size_t print_length()
		{
			return print_length(wrap_f<f>{});
		}

	};


}

// ************************************************************************************


namespace expr
{
	//! Print a formatted expression.
	/*!
	 * The expression is printed to the log with the given description, which
	 * can be formatted.
	 * 
	 * \param expr The expression to be printed.
	 * \param fmt The string format used in the printing.
	 * \param ... The arguments to the print format.
	 */
	template<typename E>
	size_t printf(E const& expr, const char* fmt = "", ...)
	{
		size_t n = 0;
		if (*fmt)
		{
			va_list list;
			va_start(list, fmt);

			char buffer[LINE_READ_BUFFER];
			vsnprintf(buffer, LINE_READ_BUFFER, fmt, list);
			n += fprintf(SYMPHAS_LOG, "%s " SYEX_EQN_SEP " ", buffer);

			va_end(list);
		}

		n += expr.print(SYMPHAS_LOG);
		n += fprintf(SYMPHAS_LOG, "\n");
		return n;
	}


	//! Print a formatted expression list.
	/*!
	 * The list of expressions is printed to the log, sequentially with respect
	 * to the order given. Applies the individual expr::printf() function to
	 * each of the expressions using the same format and argument.
	 * 
	 * Returns a tuple of the length of each string that was written.
	 *
	 * \param exprs The list of expressions to be printed.
	 * \param fmt The string format used in the printing of all expressions.
	 * \param args... The arguments to the print format, for each expression.
	 */
	template<typename E0, typename... Es, size_t... Is, typename... Args>
	auto printf(std::tuple<E0, Es...> const& exprs, std::index_sequence<Is...>, const char* fmt = "", Args&&... args)
	{
		return std::make_tuple(printf(std::get<Is>(exprs), fmt, std::forward<Args>(args)...)...);
	}

	//! Print a formatted expression list.
	/*!
	 * The list of expressions is printed to the log, sequentially with respect
	 * to the order given. Applies the individual expr::printf() function to
	 * each of the expressions using the same format and argument.
	 *
	 * \param exprs The list of expressions to be printed.
	 * \param fmt The string format used in the printing of all expressions.
	 * \param args... The arguments to the print format, for each expression.
	 */
	template<typename E0, typename... Es, typename... Args>
	size_t printf(std::tuple<E0, Es...> const& exprs, const char* fmt = "", Args&&... args)
	{
		return printf(exprs, std::make_index_sequence<1 + sizeof...(Es)>{}, fmt, std::forward<Args>(args)...);
	}

	//! Print a list of expressions in a format.
	/*!
	 * If the list is empty, then nothing is printed.
	 * 
	 * See expr::printf().
	 */
	template<typename... Args>
	auto printf(std::tuple<> const&, const char* = "", Args&&...) { return std::tuple<size_t>{}; }



	//! Prepends coefficient to an expression.
	/*!
	 * The given value is prepended to the expression as strings. The format
	 * of value changes depending on what it is.
	 * 
	 * \param out The string to put the printed expression.
	 * \param expr The expression which is printed.
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t print_with_coeff(char* out, const char* expr, V value);

	//! Prepends coefficient to an expression.
	/*!
	 * Specialization based on expr::print_with_coeff() when the identity
	 * constant is passed.
	 *
	 * \param out The string to put the printed expression.
	 * \param expr The expression which is printed.
	 */
	template<>
	inline size_t print_with_coeff(char* out, const char* expr, OpIdentity)
	{
		return sprintf(out, "%s", expr);
	}

	//! Prepends coefficient to an expression.
	/*!
	 * Specialization based on expr::print_with_coeff() when the negative
	 * identity constant is passed.
	 *
	 * \param out The string to put the printed expression.
	 * \param expr The expression which is printed.
	 */
	template<>
	inline size_t print_with_coeff(char* out, const char* expr, OpNegIdentity)
	{
		return sprintf(out, "-%s", expr);
	}

	//! Prints only the coefficient and coefficient separator.
	/*!
	 * Specialization which only prints the coefficient, without an
	 * appended expression. The appended expression is considered to be
	 * an empty string. The format
	 * of value changes depending on what it is.
	 *
	 * \param out The string to put the printed expression.
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t print_with_coeff(char* out, V value)
	{
		return print_with_coeff(out, "", value);
	}


	template<typename V>
	size_t print_with_coeff(char* out, const char* expr, V value)
	{
		if (value == V{})
		{
			return sprintf(out, SYEX_COEFF_PRINT_FMT, SYEX_ZERO_COEFF, expr);
		}
		else
		{
			if (value == value * value)
			{
				return print_with_coeff(out, expr, OpIdentity{});
			}
			else if (value == value * value * value)
			{
				return print_with_coeff(out, expr, OpNegIdentity{});
			}
			else
			{
				size_t n = make_literal(value).print(out);
				n += sprintf(out + n, SYEX_COEFF_PRINT_FMT_SEP "%s", expr);
				return n;
			}
		}
	}

	//! Prepends coefficient to an expression.
	/*!
	 * Overload based on expr::print_with_coeff() when printing
	 * to a file.
	 *
	 * \param out The file to put the printed expression.
	 * \param expr The expression which is printed.
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t print_with_coeff(FILE* out, const char* expr, V value);

	//! Prepends coefficient to an expression.
	/*!
	 * Overload based on expr::print_with_coeff() when printing
	 * to a file and the value is the identity.
	 *
	 * \param out The file to put the printed expression.
	 * \param expr The expression which is printed.
	 */
	template<>
	inline size_t print_with_coeff(FILE* out, const char* expr, OpIdentity)
	{
		return fprintf(out, "%s", expr);
	}

	//! Prepends coefficient to an expression.
	/*!
	 * Overload based on expr::print_with_coeff() when printing
	 * to a file and the value is the negative identity.
	 *
	 * \param out The file to put the printed expression.
	 * \param expr The expression which is printed.
	 */
	template<>
	inline size_t print_with_coeff(FILE* out, const char* expr, OpNegIdentity)
	{
		return fprintf(out, "-%s", expr);
	}


	//! Prints only the coefficient and coefficient separator.
	/*!
	 * Specialization which only prints the coefficient, without an
	 * appended expression. 
	 * Overload based on expr::print_with_coeff() when printing
	 * to a file and the value is the negative identity.
	 *
	 * \param out The file to put the printed expression.
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t print_with_coeff(FILE* out, V value)
	{
		return print_with_coeff(out, "", value);
	}


	template<typename V>
	size_t print_with_coeff(FILE* out, const char* expr, V value)
	{
		if (value == V{})
		{
			return fprintf(out, SYEX_COEFF_PRINT_FMT, SYEX_ZERO_COEFF, expr);
		}
		else
		{
			if (value == value * value)
			{
				return print_with_coeff(out, expr, OpIdentity{});
			}
			else if (value == value * value * value)
			{
				return print_with_coeff(out, expr, OpNegIdentity{});
			}
			else
			{
				size_t n = make_literal(value).print(out);
				n += fprintf(out, SYEX_COEFF_PRINT_FMT_SEP "%s", expr);
				return n;
			}
		}
	}

	//! Print the length of the coefficient string that will be preprended.
	/*! 
	 * Prints the number of characters that will be used to print the 
	 * coefficient which is prepended to an expression.
	 * 
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t coeff_print_length(V value);

	//! Print the length of the zero coefficient string.
	/*!
	 * Overload for the length of a coefficient equal to 0 (the additive
	 * identity).
	 * Prints the number of characters that will be used to print the
	 * coefficient which is prepended to an expression.
	 */
	template<>
	inline size_t coeff_print_length(OpVoid)
	{
		return STR_ARR_LEN(SYEX_COEFF_PRINT_FMT_SEP SYEX_ZERO_COEFF) - 1;
	}

	//! Print the length of the string for the coefficient equal to 1.
	/*!
	 * Overload for the length of a coefficient equal to 1 (the multiplicative
	 * identity).
	 * Prints the number of characters that will be used to print the
	 * coefficient which is prepended to an expression.
	 */
	template<>
	inline size_t coeff_print_length(OpIdentity)
	{
		return 0;
	}

	//! Print the length of the string for the coefficient equal to -1.
	/*!
	 * Overload for the length of a coefficient equal to 1 (the negative of the
	 * multiplicative identity).
	 * Prints the number of characters that will be used to print the
	 * coefficient which is prepended to an expression.
	 */
	template<>
	inline size_t coeff_print_length(OpNegIdentity)
	{
		return 1;
	}

	template<typename V>
	size_t coeff_print_length(V value)
	{

		if (value == V{})
		{
			return coeff_print_length(OpVoid{});
		}
		else
		{
			if (value == value * value)
			{
				return coeff_print_length(OpIdentity{});
			}
			else if (value == value * value * value)
			{
				return coeff_print_length(OpNegIdentity{});
			}
			else
			{
				size_t n = make_literal(value).print_length();
				n += STR_ARR_LEN(SYEX_COEFF_PRINT_FMT_SEP) - 1;
				return n;
			}
		}
	}

}


#undef SYEX_EQN_SEP
#undef SYEX_PRINT_FUNCTION_STRUCT
#undef SYEX_DEFINE_FUNCTION_STRUCT_DELEGATE
#undef SYEX_PRINT_FUNCTION_STRUCT

#else

namespace expr
{
	//! Empty implementation of printing an expression.
	/*!
	 * Empty implementation for the print function when #PRINTABLE_EQUATIONS is
	 * not used. Therefore, the print function can still be invoked but there is no effect.
	 */ 
	template<typename... T>
	auto printf(T&&...) {}


	//! Empty implementation for obtaining a variable's string represenation.
	/*!
	 * Empty implementation for the obtaining the name of a variable when #PRINTABLE_EQUATIONS is
	 * not used. Therefore, this function can still be invoked but there is no effect.
	 */
	template<typename... T>
	auto get_op_name(T&&...)
	{
		return "";
	}

	template<typename... T>
	size_t print_with_coeff(T&&...)
	{
		return 0;
	}

	template<typename... T>
	size_t coeff_print_length(T&&...)
	{
		return 0;
	}
}

#endif

