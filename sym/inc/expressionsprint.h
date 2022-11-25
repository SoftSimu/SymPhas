
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

#include "expressionlib.h"


template<typename Sp>
struct Solver;


namespace symphas::internal
{

	struct wrap_base {};
	template<auto f>
	struct wrap_f : wrap_base {};
}

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
#define SYEX_DERIV_APPLIED_EXPR_LEN STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B)


#define ENABLE_UNICODE
#undef ENABLE_UNICODE

#ifdef LATEX_PLOT
//! Display string of the derivative.
#define SYEX_DIRECTIONAL_DERIV_VAR_STR(NAME, VALUE, AXIS) "\\frac{\\partial^{" VALUE "}" NAME "}{\\partial " AXIS "^{" VALUE "}}"
#define SYEX_DIRECTIONAL_DERIV_STR(VALUE, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR("", VALUE, AXIS)
#define SYEX_DIRECTIONAL_DERIV_1_VAR_STR(NAME, AXIS) "\\frac{\\partial" NAME "}{\\partial " AXIS "}"
#define SYEX_DIRECTIONAL_DERIV_1_STR(AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", AXIS)
#define SYEX_DERIV_STR(VALUE) "\\nabla^{" VALUE "}"
#define SYEX_DERIV_STR_1 "\\vec{\\nabla}"
#else
//! Display string of the derivative.
#define SYEX_DIRECTIONAL_DERIV_VAR_STR(NAME, VALUE, AXIS) SYEX_DERIV_APPLIED_EXPR_FMT_A "d" VALUE NAME "/d" AXIS VALUE SYEX_DERIV_APPLIED_EXPR_FMT_B
#define SYEX_DIRECTIONAL_DERIV_STR(VALUE, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR("", VALUE, AXIS)
#define SYEX_DIRECTIONAL_DERIV_1_VAR_STR(NAME, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR(NAME, "", AXIS)
#define SYEX_DIRECTIONAL_DERIV_1_STR(AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", AXIS)

#ifdef ENABLE_UNICODE
#define SYEX_DERIV_STR(VALUE) "▼^" VALUE ""
#else
#define SYEX_DERIV_STR(VALUE) "V^" VALUE ""
#endif
#define SYEX_DERIV_STR_1 SYEX_DERIV_STR("")
#endif


#define SYEX_DIRECTIONAL_DERIV_VAR_FMT(NAME, ORDER, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR("%s", "%zd", "%c"), ORDER, NAME, \
	((AXIS == Axis::X) ? 'x' : (AXIS == Axis::Y) ? 'y' : (AXIS == Axis::Y) ? 'z' : '?'), ORDER
#define SYEX_DIRECTIONAL_DERIV_VAR_LEN(NAME, ORDER) \
	(STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_VAR_STR("", "", "")) + 2 * symphas::lib::num_digits<ORDER>() + std::strlen(NAME) + 1)

#define SYEX_DIRECTIONAL_DERIV_FMT(ORDER, AXIS) SYEX_DIRECTIONAL_DERIV_STR("%zd", "%c"), ORDER, \
	((AXIS == Axis::X) ? 'x' : (AXIS == Axis::Y) ? 'y' : (AXIS == Axis::Y) ? 'z' : '?'), ORDER
#define SYEX_DIRECTIONAL_DERIV_LEN(ORDER) (STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_STR("", "")) + 2 * symphas::lib::num_digits<ORDER>() + 1)

#define SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(NAME, AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("%s", "%c"), \
	NAME, ((AXIS == Axis::X) ? 'x' : (AXIS == Axis::Y) ? 'y' : (AXIS == Axis::Y) ? 'z' : '?')
#define SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(NAME) (STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", "")) + std::strlen(NAME) + 1)

#define SYEX_DIRECTIONAL_DERIV_1_FMT(AXIS) SYEX_DIRECTIONAL_DERIV_1_STR("%c"), \
	((AXIS == Axis::X) ? 'x' : (AXIS == Axis::Y) ? 'y' : (AXIS == Axis::Y) ? 'z' : '?')
#define SYEX_DIRECTIONAL_DERIV_1_LEN STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_1_STR("") + 1)



#define SYEX_DERIV_STR_FMT(ORDER) SYEX_DERIV_STR("%zd"), ORDER
#define SYEX_DERIV_STR_LEN(ORDER) (STR_ARR_LEN(SYEX_DERIV_STR("")) + symphas::lib::num_digits<ORDER>())

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
//
//
//#define SYEX_PRINT_DERIV_STRUCT(ORDER, FMT) \
//	template<Axis ax>																\
//	struct print_deriv<ax, ORDER>													\
//	{																				\
//		static const char c = (ax == Axis::X) ? 'x' : (ax == Axis::Y) ? 'y' : (ax == Axis::Y) ? 'z' : '?' \
//		\
//		static size_t print(FILE* out) { return fprintf(out, FMT, c); }				\
//		static size_t print(char* out) { return sprintf(out, FMT, c); }				\
//		static constexpr size_t print_length() { return STR_ARR_LEN(FMT) - 1; }		\
//	};


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

#define SYEX_TENSOR_FMT_A "\\begin{bmatrix} "
#define SYEX_TENSOR_FMT_B "\\end{bmatrix} "
#define SYEX_TENSOR_COLUMN_SEP_FMT " \\\\"
#define SYEX_TENSOR_ROW_SEP_FMT " & "
#define SYEX_TENSOR_EMPTY_FMT "0"

#define SYEX_VECTOR_FMT_A SYEX_TENSOR_FMT_A
#define SYEX_VECTOR_FMT_B SYEX_TENSOR_FMT_B
#define SYEX_VECTOR_FMT_SEP SYEX_TENSOR_COLUMN_SEP_FMT
#else

#define SYEX_COEFF_PRINT_FMT_SEP "."


#define SYEX_TENSOR_FMT_A "["
#define SYEX_TENSOR_FMT_B "]"
#define SYEX_TENSOR_COLUMN_SEP_FMT ";"
#define SYEX_TENSOR_ROW_SEP_FMT ","
#define SYEX_TENSOR_EMPTY_FMT "0"

#define SYEX_VECTOR_FMT_A SYEX_TENSOR_FMT_A
#define SYEX_VECTOR_FMT_B SYEX_TENSOR_FMT_B
#define SYEX_VECTOR_FMT_SEP SYEX_TENSOR_COLUMN_SEP_FMT
#endif

// N is the number of rows (entries in a column)
// M is the number of columns (entries in a row)
#define SYEX_TENSOR_FMT_LEN(N, M) \
STR_ARR_LEN(SYEX_TENSOR_FMT_A SYEX_TENSOR_FMT_B) + \
((N - 1) * STR_ARR_LEN(SYEX_TENSOR_COLUMN_SEP_FMT)) + \
((M - 1) * STR_ARR_LEN(SYEX_TENSOR_ROW_SEP_FMT)) + \
((M * N - 1) * STR_ARR_LEN(SYEX_TENSOR_EMPTY_FMT))

#define SYEX_VECTOR_FMT_LEN STR_ARR_LEN(SYEX_VECTOR_FMT_A SYEX_VECTOR_FMT_SEP SYEX_VECTOR_FMT_B)

//! Appearance of a coefficient in front of another string.
#define SYEX_COEFF_PRINT_FMT "%s" SYEX_COEFF_PRINT_FMT_SEP "%s"
#define SYEX_COEFF_SEP_LEN (sizeof(SYEX_COEFF_PRINT_FMT_SEP) / sizeof(char) - 1)

#define SYEX_ZERO_COEFF "0"




#ifdef LATEX_PLOT
#define SYEX_FT_OF_EXPR_FMT_A "\\hat{\\mathcal{F}}\\left\\{"
#define SYEX_FT_OF_EXPR_FMT_B "\\right\\}"
#define SYEX_FT_OF_OP_FMT_A "\\hat{"
#define SYEX_FT_OF_OP_FMT_B "}"

#define SYEX_IFT_OF_EXPR_FMT_A "\\hat{\\mathcal{F}}^{-1}\\left\\{"
#define SYEX_IFT_OF_EXPR_FMT_B SYEX_FT_OF_EXPR_FMT_B
#define SYEX_IFT_OF_OP_FMT_A SYEX_FT_OF_OP_FMT_A
#define SYEX_IFT_OF_OP_FMT_B SYEX_FT_OF_OP_FMT_B

#else
#define SYEX_FT_OF_EXPR_FMT_A "F{"
#define SYEX_FT_OF_EXPR_FMT_B "}"
#define SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_EXPR_FMT_A
#define SYEX_FT_OF_OP_FMT_B SYEX_FT_OF_EXPR_FMT_B

#define SYEX_IFT_OF_EXPR_FMT_A "F-1{"
#define SYEX_IFT_OF_EXPR_FMT_B SYEX_FT_OF_EXPR_FMT_B
#define SYEX_IFT_OF_OP_FMT_A SYEX_IFT_OF_EXPR_FMT_A
#define SYEX_IFT_OF_OP_FMT_B SYEX_FT_OF_EXPR_FMT_B

#endif

#define SYEX_FT_OF_EXPR_FMT SYEX_FT_OF_EXPR_FMT_A "%s" SYEX_FT_OF_EXPR_FMT_B
#define SYEX_FT_OF_EXPR_FMT_LEN (STR_ARR_LEN(SYEX_FT_OF_EXPR_FMT_A SYEX_FT_OF_EXPR_FMT_B) - 1)
#define SYEX_FT_OF_OP_FMT SYEX_FT_OF_OP_FMT_A "%s" SYEX_FT_OF_OP_FMT_B
#define SYEX_FT_OF_OP_FMT_LEN (STR_ARR_LEN(SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_OP_FMT_B) - 1)

#define SYEX_IFT_OF_EXPR_FMT SYEX_IFT_OF_EXPR_FMT_A "%s" SYEX_IFT_OF_EXPR_FMT_B
#define SYEX_IFT_OF_EXPR_FMT_LEN (STR_ARR_LEN(SYEX_IFT_OF_EXPR_FMT_A SYEX_IFT_OF_EXPR_FMT_B) - 1)
#define SYEX_IFT_OF_OP_FMT SYEX_IFT_OF_OP_FMT_A "%s" SYEX_IFT_OF_OP_FMT_B
#define SYEX_IFT_OF_OP_FMT_LEN (STR_ARR_LEN(SYEX_IFT_OF_OP_FMT_A SYEX_IFT_OF_OP_FMT_B) - 1)


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
		static inline const char gradlaplacian_fmt[] = "\\nabla^2";
		static inline const char gradient_fmt[] = "";

#else
		static inline const char laplacian_fmt[] = SYEX_DERIV_STR("2");
		static inline const char bilaplacian_fmt[] = SYEX_DERIV_STR("4");
		static inline const char gradlaplacian_fmt[] = SYEX_DERIV_STR("2");
		static inline const char gradient_fmt[] = "";

#endif
	}




	template<size_t O, Axis ax = Axis::X, bool is_directional = false>
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
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					size_t n = 0;
					n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					return n;
				}
				else if constexpr (O % 2 == 1)
				{
					size_t n = 0;
					n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					n += fprintf(out, SYEX_DERIV_STR_FMT(O - 1));
					return n;
				}
				else
				{
					return fprintf(out, SYEX_DERIV_STR_FMT(O));
				}
			}
			else
			{
				return fprintf(out, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
			}
		}

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output file, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 *
		 * \param out The file to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		static size_t print(FILE* out, const char* name)
		{
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					size_t n = 0;
					n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(name, ax));
					return n;
				}
				else if constexpr (O % 2 == 1)
				{
					size_t n = 0;
					n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					n += fprintf(out, SYEX_DERIV_STR_FMT(O - 1));
					n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT, name);
					return n;
				}
				else
				{
					size_t n = 0;
					n += fprintf(out, SYEX_DERIV_STR_FMT(O));
					n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT, name);
					return n;
				}
			}
			else
			{
				size_t n = 0;
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
				n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT, name);
				return n;
			}
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
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					size_t n = 0;
					n += sprintf(out + n, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					return n;
				}
				else if constexpr (O % 2 == 1)
				{
					size_t n = 0;
					n += sprintf(out + n, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					n += sprintf(out + n, SYEX_DERIV_STR_FMT(O - 1));
					return n;
				}
				else
				{
					return sprintf(out, SYEX_DERIV_STR_FMT(O));
				}
			}
			else
			{
				return sprintf(out, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
			}
		}

		//! Print the derivative of an order inferred from the function.
		/*!
		 * Print the derivative of an order inferred from the function.
		 * The derivative output is printed the output string, and
		 * formatted using the given string representing the
		 * expression which the derivative is applied to.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		static size_t print(char* out, const char* name)
		{
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					size_t n = 0;
					n += sprintf(out + n, SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(name, ax));
					return n;
				}
				else if constexpr (O % 2 == 1)
				{
					size_t n = 0;
					n += sprintf(out + n, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
					n += sprintf(out + n, SYEX_DERIV_STR_FMT(O - 1));
					n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT, name);
					return n;
				}
				else
				{
					size_t n = 0;
					n += sprintf(out + n, SYEX_DERIV_STR_FMT(O));
					n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT, name);
					return n;
				}
			}
			else
			{
				size_t n = 0;
				n += sprintf(out + n, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
				n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT, name);
				return n;
			}
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
			if constexpr (O == 1)
			{
				return SYEX_DIRECTIONAL_DERIV_1_STR_LEN;
			}
			else if constexpr (O % 2 == 1)
			{
				return SYEX_DIRECTIONAL_DERIV_1_STR_LEN + SYEX_DERIV_STR_LEN(O - 1);
			}
			else
			{
				return SYEX_DERIV_STR_LEN(O);
			}
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
		static size_t print_length(const char* name)
		{
			if constexpr (O == 1)
			{
				return SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(NAME);
			}
			else if constexpr (O % 2 == 1)
			{
				return std::strlen(name) + SYEX_DIRECTIONAL_DERIV_1_LEN
					+ SYEX_DERIV_STR_LEN(O - 1) + SYEX_DERIV_APPLIED_EXPR_LEN;
			}
			else
			{
				return std::strlen(name) + SYEX_DERIV_STR_LEN(O) + SYEX_DERIV_APPLIED_EXPR_LEN;
			}
		}
	};

	//SYEX_PRINT_DERIV_STRUCT(1, deriv_formats::gradient_fmt);
	//SYEX_PRINT_DERIV_STRUCT(2, deriv_formats::laplacian_fmt);
	//SYEX_PRINT_DERIV_STRUCT(3, deriv_formats::gradlaplacian_fmt);
	//SYEX_PRINT_DERIV_STRUCT(4, deriv_formats::bilaplacian_fmt);


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



	//! Print the length of the coefficient string that will be preprended.
	/*!
	 * Prints the number of characters that will be used to print the
	 * coefficient which is prepended to an expression.
	 *
	 * \param value The value of the coefficient.
	 */
	template<typename V>
	size_t coeff_print_length(V const& value);

	template<typename T>
	size_t print_tensor_entries_1(char* out, T const& value, size_t P0, size_t P1, size_t N, size_t M)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_TENSOR_FMT_A);

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < M; ++j)
			{
				if (i == P0 && j == P1)
				{
					n += make_literal(expr::eval(value)).print(out + n);
				}
				else
				{
					n += sprintf(out + n, SYEX_TENSOR_EMPTY_FMT);
				}

				if (j < M - 1)
				{
					n += sprintf(out + n, SYEX_TENSOR_ROW_SEP_FMT);
				}
			}
			if (i < N - 1)
			{
				n += sprintf(out + n, SYEX_TENSOR_COLUMN_SEP_FMT);
			}
		}
		n += sprintf(out + n, SYEX_TENSOR_FMT_B);
		return n;
	}

	template<typename T, size_t P0, size_t P1, size_t N, size_t M>
	size_t print_tensor(char* out, T const& value, std::index_sequence<P0, P1>, std::index_sequence<N, M>)
	{
		return print_tensor_entries_1(out, value, P0, P1, N, M);
	}

	template<typename T, size_t P0, size_t N>
	size_t print_tensor(char* out, T const& value, std::index_sequence<P0>, std::index_sequence<N>)
	{
		return print_tensor_entries_1(out, value, P0, 0, N, 1);
	}

	template<typename T, size_t N, size_t M>
	size_t tensor_print_length(T const& value, std::index_sequence<N, M>)
	{
		size_t n = coeff_print_length(value);
		n += SYEX_TENSOR_FMT_LEN(N, M);
		return n;
	}

	template<typename T, size_t N>
	size_t tensor_print_length(T const& value, std::index_sequence<N>)
	{
		size_t n = coeff_print_length(value);
		n += SYEX_TENSOR_FMT_LEN(N, 1);
		return n;
	}

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
	inline size_t print_with_coeff(char* out, const char* expr, OpNegIdentity)
	{
		return sprintf(out, "-%s", expr);
	}

	template<size_t N, size_t D>
	size_t print_with_coeff(char* out, const char* expr, OpFractionLiteral<N, D>)
	{
		size_t n = 0;
		n += OpFractionLiteral<N, D>{}.print(out);
		n += sprintf(out + n, "%s", expr);
		return n;
	}

	template<size_t N, size_t D>
	size_t print_with_coeff(char* out, const char* expr, OpNegFractionLiteral<N, D>)
	{
		size_t n = 0;
		n += OpNegFractionLiteral<N, D>{}.print(out);
		n += sprintf(out + n, "%s", expr);
		return n;
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

	template<typename T, size_t... Ns>
	size_t print_with_coeff(char* out, const char* expr, OpTensor<T, Ns...> const& value)
	{
		using p_seq = symphas::lib::seq_join_t<symphas::lib::types_before_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;
		using n_seq = symphas::lib::seq_join_t<symphas::lib::types_after_at_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;

		size_t n = print_tensor(out, symphas::internal::tensor_cast::cast(value), p_seq{}, n_seq{});
		n += sprintf(out + n, "%s", expr);
		return n;
	}

	template<typename T, size_t D>
	size_t print_with_coeff(char* out, const char* expr, any_vector_t<T, D> const& value)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_VECTOR_FMT_A);
		for (iter_type i = 0; i < D; ++i)
		{
			n += make_literal(value[i]).print(out);
			if (i < D - 1)
			{
				n += sprintf(out + n, SYEX_VECTOR_FMT_SEP);
			}
		}
		n += sprintf(out + n, SYEX_VECTOR_FMT_B "%s", expr);
		return n;
	}

	template<typename... Ts>
	size_t print_with_coeff(char* out, const char* expr, OpAdd<Ts...> const& value)
	{
		return print_with_coeff(out, expr, value.eval());
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
	size_t print_with_coeff(FILE* out, const char* expr, V const& value);

	//! Prepends coefficient to an expression.
	/*!
	 * Overload based on expr::print_with_coeff() when printing
	 * to a file and the value is the identity.
	 *
	 * \param out The file to put the printed expression.
	 * \param expr The expression which is printed.
	 */
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
	inline size_t print_with_coeff(FILE* out, const char* expr, OpNegIdentity)
	{
		return fprintf(out, "-%s", expr);
	}

	template<size_t N, size_t D>
	size_t print_with_coeff(FILE* out, const char* expr, OpFractionLiteral<N, D>)
	{
		size_t n = 0;
		n += OpFractionLiteral<N, D>{}.print(out);
		n += fprintf(out, "%s", expr);
		return n;
	}

	template<size_t N, size_t D>
	size_t print_with_coeff(FILE* out, const char* expr, OpNegFractionLiteral<N, D>)
	{
		size_t n = 0;
		n += OpNegFractionLiteral<N, D>{}.print(out);
		n += fprintf(out, "%s", expr);
		return n;
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
	size_t print_with_coeff(FILE* out, V const& value)
	{
		return print_with_coeff(out, "", value);
	}

	template<typename V>
	size_t print_with_coeff(FILE* out, const char* expr, V const& value)
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

	template<typename T, size_t... Ns>
	size_t print_with_coeff(FILE* out, const char* expr, OpTensor<T, Ns...> const& value)
	{
		size_t len = coeff_print_length(value) + 1;
		char* buffer = new char[len];
		print_with_coeff(buffer, "", value);
		return fprintf(out, "%s%s", buffer, expr);
	}

	template<typename... Ts>
	size_t print_with_coeff(FILE* out, const char* expr, OpAdd<Ts...> const& value)
	{
		return print_with_coeff(out, expr, value.eval());
	}

	template<typename T, size_t D>
	size_t print_with_coeff(FILE* out, const char* expr, any_vector_t<T, D> const& value)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_VECTOR_FMT_A);
		for (iter_type i = 0; i < D; ++i)
		{
			n += make_literal(value[i]).print(out);
			if (i < D - 1)
			{
				n += fprintf(out, SYEX_VECTOR_FMT_SEP);
			}
		}
		n += fprintf(out, SYEX_VECTOR_FMT_B "%s", expr);
		return n;
	}

	//! Print the length of the zero coefficient string.
	/*!
	 * Overload for the length of a coefficient equal to 0 (the additive
	 * identity).
	 * Prints the number of characters that will be used to print the
	 * coefficient which is prepended to an expression.
	 */
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
	inline size_t coeff_print_length(OpNegIdentity)
	{
		return 1;
	}

	template<typename V>
	size_t coeff_print_length(V const& value)
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

	template<typename T, size_t... Ns>
	size_t coeff_print_length(OpTensor<T, Ns...> const& value)
	{
		using n_seq = symphas::lib::seq_join_t<symphas::lib::types_after_at_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;

		return tensor_print_length(symphas::internal::tensor_cast::cast(value), n_seq{});
	}

	template<typename... Ts>
	size_t coeff_print_length(OpAdd<Ts...> const& value)
	{
		return coeff_print_length(value.eval());
	}

	template<typename T, size_t D>
	size_t coeff_print_length(any_vector_t<T, D> const& value)
	{
		size_t n = 0;
		n += 2;
		for (iter_type i = 0; i < D; ++i)
		{
			n += coeff_print_length(value[i]);
			if (i < D - 1)
			{
				n += 1;
			}
		}
		return n;
	}

	template<size_t N>
	struct expr_name_arr
	{
		char value[N];
		char* new_str()
		{
			char* name = new char[N];
			std::strcpy(name, value);
			return name;
		}
	};

	namespace
	{


		template<size_t N, size_t LN, size_t I, size_t L0, size_t L = symphas::lib::num_digits<N>()>
		void print_number(expr_name_arr<L0>& ptr)
		{
			constexpr size_t N0 = N / symphas::math::fixed_pow<10, L - I - 1>;
			ptr.value[LN + I] =
				(N0 == 0) ? '0' :
				(N0 == 1) ? '1' :
				(N0 == 2) ? '2' :
				(N0 == 3) ? '3' :
				(N0 == 4) ? '4' :
				(N0 == 5) ? '5' :
				(N0 == 6) ? '6' :
				(N0 == 7) ? '7' :
				(N0 == 8) ? '8' :
				(N0 == 9) ? '9' :
				' ';

			if constexpr (L > 1)
			{
				print_number<N0* symphas::math::fixed_pow<10, L - I - 1> -N, LN, I + 1>(ptr);
			}
		}
	}


	template<int N, size_t L, size_t L0 = L - 1, size_t NN = (N < 0) ? static_cast<size_t>(-N) : static_cast<size_t>(N),
#ifdef LATEX_PLOT
		size_t L1 = 1
#else
		size_t L1 = 0
#endif
		, size_t L2 = ((N < 0) ? 1 : 0), size_t L3 = symphas::lib::num_digits<NN>()>
		expr_name_arr<L0 + L1 + L2 + L3 + 3> print_with_subscript(const char(&term)[L])
	{
		expr_name_arr<L0 + L1 + L2 + L3 + 3> out;
		for (iter_type i = 0; i < L0; ++i)
		{
			out.value[i] = term[i];
		}

		if constexpr (L1)
		{
			out.value[L0] = '_';
		}
		out.value[L0 + L1] = '{';
		if constexpr (L2 > 0)
		{
			out.value[L0 + L1 + 1] = '-';
		}

		constexpr size_t LN = L0 + L1 + 1 + L2;
		print_number<NN, LN, 0>(out);

		out.value[LN + L3] = '}';
		out.value[LN + L3 + 1] = '\0';

		return out;
	}

	template<int N0, int... Ns, size_t L>
	auto print_with_subscripts(const char(&term)[L])
	{
		if constexpr (sizeof...(Ns) == 0)
		{
			return print_with_subscript<N0>(term);
		}
		else
		{
			return print_with_subscripts<Ns...>(print_with_subscript<N0>(term).value);
		}
	}

	inline auto get_fourier_name(const char* name)
	{
		return std::string(SYEX_FT_OF_OP_FMT_A) + std::string(name) + std::string(SYEX_FT_OF_OP_FMT_B);
	}

	template<typename E>
	auto get_fourier_name(OpExpression<E> const& e)
	{
		char* name = new char[e.print_length() + 1];
		e.print(name);
		return std::string(SYEX_FT_OF_EXPR_FMT_A) + std::string(name) + std::string(SYEX_FT_OF_EXPR_FMT_B);
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

