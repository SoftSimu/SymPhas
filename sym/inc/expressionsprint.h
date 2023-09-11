
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

#include "expressionlogic.h"

 //! \cond
#ifdef EXPR_EXPORTS
#define DLLEXPR DLLEXPORT
#else
#define DLLEXPR DLLIMPORT
#endif
//! \endcond

namespace symphas::internal
{

	struct wrap_base {};
	template<auto f>
	struct wrap_f : wrap_base {};
}

namespace expr
{

	template<size_t N = 0>
	struct expr_name_arr
	{
        expr_name_arr(const char* arr) : value{}
        {
            std::strcpy(arr, value);
        }

        expr_name_arr() : value{} {}

		char value[N];
		char* new_str()
		{
			char* name = new char[N];
			std::strcpy(name, value);
			return name;
		}
	};

	template<>
	struct expr_name_arr<0>
	{
		char* value;

		expr_name_arr(len_type len = 0) :
			value{ (len > 0) ? new char[len] {} : nullptr } {}
		expr_name_arr(const char* value) :
			value{ (std::strlen(value) > 0) ? new char[std::strlen(value) + 1] {} : nullptr} 
		{
			if (this->value != nullptr)
			{
				std::strcpy(this->value, value);
			}
		}
		expr_name_arr(expr_name_arr<0> const& other) :
			expr_name_arr((other.value != nullptr) ? len_type(std::strlen(other.value) + 1) : 0)
		{
			if (value != nullptr)
			{
				std::strcpy(value, other.value);
			}
		}

		expr_name_arr(expr_name_arr<0>&& other) : expr_name_arr()
		{
			swap(*this, other);
		}

		expr_name_arr<0> operator=(expr_name_arr<0> other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(expr_name_arr<0>& first, expr_name_arr<0>& second)
		{
			using std::swap;
			swap(first.value, second.value);
		}

		operator const char* () const
		{
			return value;
		}

		operator char* ()
		{
			return value;
		}

		char* new_str() const
		{
			char* name = new char[std::strlen(value) + 1];
			std::strcpy(name, value);
			return name;
		}

		~expr_name_arr()
		{
			delete[] value;
		}
	};

	expr_name_arr(len_type)->expr_name_arr<0>;
	expr_name_arr(const char*)->expr_name_arr<0>;

	namespace
	{

		template<typename E = OpVoid>
		struct limit_data
		{
			limit_data(E const& e = E{}, iter_type index = 0, iter_type offset = 0, bool dynamic_index = false) :
				e{ e }, index{ index }, offset{ offset }, dynamic_index{ dynamic_index } {}
			limit_data(iter_type index = 0, iter_type offset = 0, bool dynamic_index = false) :
				e{}, index{ index }, offset{ offset }, dynamic_index{ dynamic_index } {}

			template<typename E0>
			limit_data<add_result_t<E, E0>> operator+(limit_data<E0> const& other)
			{
				return { e + other.e, index + other.index, offset + other.offset };
			}

			auto operator!() const
			{
				return !std::is_same<E, OpVoid>::value;
			}

			template<typename E0>
			auto operator>(limit_data<E0> const& other) const
			{
				return index > other.index;
			}

			template<typename E0>
			auto operator<(limit_data<E0> const& other) const
			{
				return index < other.index;
			}

			E e;
			iter_type index;
			iter_type offset;
			bool dynamic_index;
		};

		limit_data(iter_type, iter_type)->limit_data<OpVoid>;


		auto get_limit_data(DynamicIndex const& index);
		template<typename V>
		auto get_limit_data(OpBinaryMul<V, DynamicIndex> const& index);

		inline auto get_limit_data(int index)
		{
			return limit_data<OpVoid>(index);
		}

		template<typename V, int N, int P>
		auto get_limit_data(OpTerm<V, expr::symbols::i_<N, P>> const& index)
		{
			return limit_data<OpTerm<V, expr::symbols::i_<N, 0>>>(expr::coeff(index) * expr::symbols::i_<N, 0>{}, 0, P);
		}

		template<size_t N>
		auto get_limit_data(expr::symbols::placeholder_N_symbol_<N>)
		{
			return limit_data<expr::symbols::placeholder_N_symbol_<N>>();
		}

		template<typename E>
		auto get_limit_data(OpExpression<E> const& e)
		{
			return limit_data<OpVoid>(static_cast<E const*>(&e)->eval());
		}

		template<typename... Es, size_t... Is>
		auto get_limit_data(OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return (get_limit_data(expr::get<Is>(e)) + ...);
		}

		template<typename... Es, size_t... Is>
		auto get_limit_data(OpAdd<Es...> const& e)
		{
			return get_limit_data(e, std::make_index_sequence<sizeof...(Es)>{});
		}

	}
}



#ifdef PRINTABLE_EQUATIONS

//! \cond

#define SYEX_BINARY_FMT "%s %s %s"
#define SYEX_ADD_SEP " + "
#define SYEX_SUB_SEP " - "
#define SYEX_BINARY_FMT_LEN (sizeof(SYEX_ADD_SEP) / sizeof(char) - 1)


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
#define SYEX_DIRECTIONAL_DERIV_1_VAR_STR(NAME, AXIS) "\\frac{\\partial " NAME "}{\\partial " AXIS "}"
#define SYEX_DIRECTIONAL_DERIV_1_STR(AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", AXIS)
#define SYEX_DERIV_STR(VALUE) "\\nabla^{" VALUE "}"
#define SYEX_DERIV_STR_1 "\\vec{\\nabla}"
#define SYEX_FUNCTIONAL_DERIV_VAR_STR(NAME, VAR) "\\frac{\\delta " NAME "}{\\delta " VAR "}"
#define SYEX_FUNCTIONAL_DERIV_STR(VAR) SYEX_FUNCTIONAL_DERIV_VAR_STR("", VAR)
#else
//! Display string of the derivative.
#define SYEX_DIRECTIONAL_DERIV_VAR_STR(NAME, VALUE, AXIS) "d" VALUE NAME "/d" AXIS VALUE
#define SYEX_DIRECTIONAL_DERIV_STR(VALUE, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR("", VALUE, AXIS)
#define SYEX_DIRECTIONAL_DERIV_1_VAR_STR(NAME, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_STR(NAME, "", AXIS)
#define SYEX_DIRECTIONAL_DERIV_1_STR(AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", AXIS)
#define SYEX_FUNCTIONAL_DERIV_VAR_STR(NAME, VAR) "{d}" NAME "/{d}" VAR
#define SYEX_FUNCTIONAL_DERIV_STR(VAR) SYEX_FUNCTIONAL_DERIV_VAR_STR("", VAR)

#ifdef ENABLE_UNICODE
#define SYEX_DERIV_STR(VALUE) "▼^" VALUE ""
#else
#define SYEX_DERIV_STR(VALUE) "V^" VALUE ""
#endif
#define SYEX_DERIV_STR_1 SYEX_DERIV_STR("")
#endif


#define SYEX_DIRECTIONAL_DERIV_VAR_WITH_OTHER_FMT(NAME, ORDER, AXIS, OTHER) SYEX_DIRECTIONAL_DERIV_VAR_STR("%s", "%zd", "%s"), ORDER, NAME, \
	((AXIS == Axis::X) ? "x" : (AXIS == Axis::Y) ? "y" : (AXIS == Axis::Z) ? "z" : OTHER), ORDER
#define SYEX_DIRECTIONAL_DERIV_VAR_FMT(NAME, ORDER, AXIS) SYEX_DIRECTIONAL_DERIV_VAR_WITH_OTHER_FMT(NAME, ORDER, AXIS, "?")
#define SYEX_DIRECTIONAL_DERIV_VAR_LEN(NAME, ORDER) \
	(STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_VAR_STR("", "", "")) + 2 * symphas::lib::num_digits<ORDER>() + std::strlen(NAME) + 1)

#define SYEX_DIRECTIONAL_DERIV_WITH_OTHER_FMT(ORDER, AXIS, OTHER) SYEX_DIRECTIONAL_DERIV_STR("%zd", "%s"), ORDER, \
	((AXIS == Axis::X) ? "x" : (AXIS == Axis::Y) ? "y" : (AXIS == Axis::Z) ? "z" : OTHER), ORDER
#define SYEX_DIRECTIONAL_DERIV_FMT(ORDER, AXIS) SYEX_DIRECTIONAL_DERIV_WITH_OTHER_FMT(ORDER, AXIS, "?")
#define SYEX_DIRECTIONAL_DERIV_LEN(ORDER) (STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_STR("", "")) + 2 * symphas::lib::num_digits<ORDER>() + 1)

#define SYEX_DIRECTIONAL_DERIV_1_VAR_WITH_OTHER_FMT(NAME, AXIS, OTHER) SYEX_DIRECTIONAL_DERIV_1_VAR_STR("%s", "%s"), \
	NAME, ((AXIS == Axis::X) ? "x" : (AXIS == Axis::Y) ? "y" : (AXIS == Axis::Z) ? "z" : OTHER)
#define SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(NAME, AXIS) SYEX_DIRECTIONAL_DERIV_1_VAR_WITH_OTHER_FMT(NAME, AXIS, "?")
#define SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(NAME) (STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_1_VAR_STR("", "")) + std::strlen(NAME) + 1)

#define SYEX_DIRECTIONAL_DERIV_1_WITH_OTHER_FMT(AXIS, OTHER) SYEX_DIRECTIONAL_DERIV_1_STR("%s"), \
	((AXIS == Axis::X) ? "x" : (AXIS == Axis::Y) ? "y" : (AXIS == Axis::Z) ? "z" : OTHER)
#define SYEX_DIRECTIONAL_DERIV_1_FMT(AXIS) SYEX_DIRECTIONAL_DERIV_1_WITH_OTHER_FMT(AXIS, "?")
#define SYEX_DIRECTIONAL_DERIV_1_LEN STR_ARR_LEN(SYEX_DIRECTIONAL_DERIV_1_STR("") + 1)

#define SYEX_FUNCTIONAL_DERIV_VAR_FMT(NAME, VAR) SYEX_FUNCTIONAL_DERIV_VAR_STR("%s", "%s"), NAME, VAR
#define SYEX_FUNCTIONAL_DERIV_VAR_LEN(NAME, VAR) (STR_ARR_LEN(SYEX_FUNCTIONAL_DERIV_VAR_STR("", "")) + std::strlen(NAME) + std::strlen(VAR))

#define SYEX_FUNCTIONAL_DERIV_FMT(VAR) SYEX_FUNCTIONAL_DERIV_STR("%s"), VAR
#define SYEX_FUNCTIONAL_DERIV_LEN(VAR) (STR_ARR_LEN(SYEX_FUNCTIONAL_DERIV_STR("")) + std::strlen(VAR))


#define SYEX_DERIV_STR_FMT(ORDER) SYEX_DERIV_STR("%zd"), ORDER
#define SYEX_DERIV_STR_LEN(ORDER) (STR_ARR_LEN(SYEX_DERIV_STR("")) + symphas::lib::num_digits<ORDER>())

//! Character used in place of an unknown identifier or derivative.
#define SYEX_UNKNOWN_FUNCTION_TOKEN "?"

//! The separation string used to separate a description from an expression.


#ifndef LATEX_PLOT
#define SYEX_SUM_SYMBOL "Sum"
#define SYEX_SUM_A "{"
#define SYEX_SUM_B "}"
#define SYEX_SUM_LIM_A "["
#define SYEX_SUM_LIM_SEP ","
#define SYEX_SUM_LIM_B "]"

#define SYEX_SUM_LIM_COMPARE_A "("
#define SYEX_SUM_LIM_COMPARE_B ")"
#define SYEX_SUM_LIM_COMPARE_SEP ","

#define SYEX_LIMIT_OFFSET_A "("
#define SYEX_LIMIT_OFFSET_B ")"

#else
#define SYEX_SUM_SYMBOL ""
#define SYEX_SUM_A "\\left("
#define SYEX_SUM_B "\\right)"
#define SYEX_SUM_LIM_A "\\sum_{"
#define SYEX_SUM_LIM_SEP "}^{"
#define SYEX_SUM_LIM_B "}"

#define SYEX_SUM_LIM_COMPARE_A SYEX_SUM_A
#define SYEX_SUM_LIM_COMPARE_B SYEX_SUM_B
#define SYEX_SUM_LIM_COMPARE_SEP ","

#define SYEX_LIMIT_OFFSET_A "_{"
#define SYEX_LIMIT_OFFSET_B "}"

#endif



#ifndef LATEX_PLOT
#define SYEX_INTEGRAL_SYMBOL "Int"
#define SYEX_INTEGRAL_A "{"
#define SYEX_INTEGRAL_B "}"
#define SYEX_INTEGRAL_LIM_A "["
#define SYEX_INTEGRAL_LIM_SEP ","
#define SYEX_INTEGRAL_LIM_B "]"
#define SYEX_INTEGRAL_DOMAIN_SYM "RxR"
#define SYEX_INTEGRAL_INTEGRATION_SYM "dx"
#else
#define SYEX_INTEGRAL_SYMBOL ""
#define SYEX_INTEGRAL_A "\\left["
#define SYEX_INTEGRAL_B "\\right]"
#define SYEX_INTEGRAL_LIM_A "\\int_{"
#define SYEX_INTEGRAL_LIM_SEP "}^{"
#define SYEX_INTEGRAL_LIM_B "}"
#define SYEX_INTEGRAL_DOMAIN_SYM "\\Omega"
#define SYEX_INTEGRAL_INTEGRATION_SYM "d\\vec{x}"
#endif


#ifdef LATEX_PLOT
#define SYEX_EQN_SEP "\\rightarrow"
#else
#define SYEX_EQN_SEP "->"
#endif

#ifdef LATEX_PLOT


#else


#endif

#ifdef LATEX_PLOT

#define SYEX_NOISE_TOKEN_POISSON "\\mathbf{P}"
#define SYEX_NOISE_TOKEN_WHITE "\\eta"
#define SYEX_NOISE_TOKEN_NONE "\\hat{" SYEX_NOISE_TOKEN_WHITE "}"
#define SYEX_NOISE_A SYEX_SUM_A
#define SYEX_NOISE_B SYEX_SUM_B

#define SYEX_ARRAY_A SYEX_INTEGRAL_A
#define SYEX_ARRAY_B SYEX_INTEGRAL_B
#define SYEX_ARRAY_SUBSCRIPT_A "_{"
#define SYEX_ARRAY_SUBSCRIPT_B "}"

#else

#define SYEX_NOISE_TOKEN_POISSON "Pois"
#define SYEX_NOISE_TOKEN_WHITE "n"
#define SYEX_NOISE_TOKEN_NONE SYEX_NOISE_TOKEN_WHITE "'"
#define SYEX_NOISE_A SYEX_SUM_A
#define SYEX_NOISE_B SYEX_SUM_B

#define SYEX_ARRAY_A SYEX_INTEGRAL_A
#define SYEX_ARRAY_B SYEX_INTEGRAL_B
#define SYEX_ARRAY_SUBSCRIPT_A "["
#define SYEX_ARRAY_SUBSCRIPT_B "]"

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
(STR_ARR_LEN(SYEX_TENSOR_FMT_A SYEX_TENSOR_FMT_B) + \
((N - 1) * STR_ARR_LEN(SYEX_TENSOR_COLUMN_SEP_FMT)) + \
((M - 1) * STR_ARR_LEN(SYEX_TENSOR_ROW_SEP_FMT)) + \
((M * N - 1) * STR_ARR_LEN(SYEX_TENSOR_EMPTY_FMT)))

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


namespace expr
{

	DLLEXPR extern int NAME_PTR_POS;					//!< Current position in selecting name for arbitrary data pointers.
	DLLEXPR extern std::vector<const void*> NAME_PTRS;	//!< List of all data pointers with names associated with them.
	DLLEXPR extern std::vector<char*> MORE_NAMES;		//!< List of overflow names for data pointers.


	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(symphas::ref<A> const& a);

	//! Gets the string name associated with the data.
	template<typename T>
	const char* get_op_name(T* ptr);


	//! Gets the string name associated with the data.
	template<>
	inline const char* get_op_name(char* a)
	{
		return a;
	}

	//! Gets the string name associated with the data.
	template<>
	inline const char* get_op_name(const char* a)
	{
		return a;
	}

	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(A const& a)
	{
		return get_op_name(expr::SymbolID<A>::get(a));
	}

	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(symphas::ref<A> const& a)
	{
		return get_op_name(expr::SymbolID<A>::get(a));
	}

	//! Specialization based on SymbolID.
	template<Axis ax, typename G>
	const char* get_op_name(VectorComponent<ax, G> const& a)
	{
		try
		{
			assert(static_cast<G const*>(&a));
			static std::map<std::string, char*> map;
			const char* name = get_op_name(SymbolID<G>::get(*static_cast<G const*>(&a)));

			auto with_component = map.find(name);
			if (with_component == map.end())
			{
				map[name] = new char[std::strlen(name) + 3];

				std::strcpy(map[name], name);
				map[name][std::strlen(name)] = '_';
				map[name][std::strlen(name) + 1] = (ax == Axis::X) ? 'x' : (ax == Axis::Y) ? 'y' : (ax == Axis::Z) ? 'z' : '?';
				map[name][std::strlen(name) + 2] = '\0';
			}

			return map[name];
		}
		catch (...)
		{
			return "???";
		}
	};

	//! Specialization based on SymbolID.
	template<size_t Z, Axis ax, typename G>
	const char* get_op_name(Variable<Z, VectorComponent<ax, G>> const& a)
	{
		return get_op_name(*static_cast<VectorComponent<ax, G> const*>(&a));
	};

	//! Gets the string name associated with the data.
	template<size_t N>
	const char* get_op_name(Variable<N> const& a)
	{
		static size_t NN = 0;
		static char** names;
		const char prefix[] = "var";
		if (N >= NN)
		{
			char** new_names = new char* [N + 1];
			for (iter_type i = 0; i < NN; ++i)
			{
				new_names[i] = new char[std::strlen(names[i]) + 1];
				std::strcpy(new_names[i], names[i]);
			}
			for (size_t i = NN; i <= N; ++i)
			{
				new_names[i] = new char[STR_ARR_LEN(prefix) + symphas::lib::num_digits(N)];
				sprintf(new_names[i], "%s%zd", prefix, N);
			}
			delete[] names;
			names = new_names;
			return names[N];
		}
		else
		{
			return names[N];
		}
	}


	//! Gets the string name associated with the data.
	inline char*& get_op_name_reference(const void* ptr)
	{
		if (!ptr)
		{
			auto get_none_name = [&] (const char* name)
			{
				char* none_name = new char[std::strlen(name) + 1];
				std::strcpy(none_name, name);
				return none_name;
			};

			static char* none = get_none_name("?");
			return none;
		}
		else
		{
			static size_t NAME_PTR_POS = 0;
			static std::vector<const void*> name_ptrs;
			static std::vector<char*> more_names;
			
			constexpr size_t MAX_NAME_COUNT = sizeof(VARIABLE_NAMES) / sizeof(*VARIABLE_NAMES);

			auto copy_var_names = [&] () 
			{
				char** variable_names = new char*[MAX_NAME_COUNT];
				for (iter_type i = 0; i < MAX_NAME_COUNT; ++i)
				{
					variable_names[i] = new char[std::strlen(VARIABLE_NAMES[i]) + 1];
					std::strcpy(variable_names[i], VARIABLE_NAMES[i]);
				}
				return variable_names;
			};

			static char** variable_names = copy_var_names();

			for (iter_type i = 0; i < NAME_PTR_POS; ++i)
			{
				if (name_ptrs[i] == ptr)
				{
					if (i < MAX_NAME_COUNT)
					{
						return variable_names[i];
					}
					else
					{
						return more_names[i - MAX_NAME_COUNT];
					}
				}
			}
			name_ptrs.push_back(ptr);

			if (NAME_PTR_POS < MAX_NAME_COUNT)
			{
				return variable_names[NAME_PTR_POS++];
			}
			else
			{
				char* name = new char[BUFFER_LENGTH_R4];
				snprintf(name, BUFFER_LENGTH_R4, VARIABLE_NAME_EXTRA_FMT, NAME_PTR_POS);
				more_names.push_back(name);
				return more_names[NAME_PTR_POS++ - MAX_NAME_COUNT];
			}
		}
	}

	//! Gets the string name associated with the data.
	template<typename T>
	const char* get_op_name(T* ptr)
	{
		return get_op_name_reference(static_cast<const void*>(ptr));
	}

	//! Gets the string name associated with the data.
	template<typename G>
	const char* get_op_name(DynamicVariable<NamedData<G*>> const& a);

}


namespace symphas::internal
{


	inline size_t print_sep(char* out, const char* sep)
	{
		return sprintf(out, "%s", sep);
	}

	inline size_t print_sep(FILE* out, const char* sep)
	{
		return fprintf(out, "%s", sep);
	}

	template<typename E>
	auto set_var_string(OpExpression<E> const& var);
	template<typename E>
	auto set_var_string(OpOperator<E> const& var);
	template<typename G>
	auto set_var_string(G const& var);

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




	template<size_t O, Axis ax = Axis::NONE, bool is_directional = false>
	struct print_deriv;

	template<size_t O, Axis ax, bool is_directional>
	struct print_deriv
	{		
		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		template<typename Sp = int>
		static size_t print(FILE* out, Sp const& = 0)
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
				if constexpr (O == 1)
				{
					return fprintf(out, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
				}
				else
				{
					return fprintf(out, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
				}
			}
		}

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The file to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename Sp = int>
		static size_t print(FILE* out, const char* name, Sp const& = 0)
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
				if constexpr (O == 1)
				{
					return fprintf(out, SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(name, ax));
				}
				else
				{
					return fprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_FMT(name, O, ax));
				}
			}
		}


		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		template<typename Sp = int>
		static size_t print(char* out, Sp const& = 0)
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
				if constexpr (O == 1)
				{
					return sprintf(out, SYEX_DIRECTIONAL_DERIV_1_FMT(ax));
				}
				else
				{
					return sprintf(out, SYEX_DIRECTIONAL_DERIV_FMT(O, ax));
				}
			}
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename Sp = int>
		static size_t print(char* out, const char* name, Sp const& = 0)
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
				if constexpr (O == 1)
				{
					return sprintf(out, SYEX_DIRECTIONAL_DERIV_1_VAR_FMT(name, ax));
				}
				else
				{
					return sprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_FMT(name, O, ax));
				}
			}
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename Sp = int>
		static size_t print_length(Sp const& = 0)
		{
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					return SYEX_DIRECTIONAL_DERIV_1_LEN;
				}
				else if constexpr (O % 2 == 1)
				{
					return SYEX_DIRECTIONAL_DERIV_1_LEN + SYEX_DERIV_STR_LEN(O - 1);
				}
				else
				{
					return SYEX_DERIV_STR_LEN(O);
				}
			}
			else
			{
				if constexpr (O == 1)
				{
					return SYEX_DIRECTIONAL_DERIV_1_LEN;
				}
				else
				{
					return SYEX_DIRECTIONAL_DERIV_LEN(O);
				}
			}
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename Sp = int>
		static size_t print_length(const char* name, Sp const& = 0)
		{
			if constexpr (!is_directional)
			{
				if constexpr (O == 1)
				{
					return SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(name);
				}
				else if constexpr (O % 2 == 1)
				{
					return std::strlen(name) + SYEX_DIRECTIONAL_DERIV_1_LEN
						+ SYEX_DERIV_STR_LEN(O - 1) + SYEX_DERIV_APPLIED_EXPR_LEN;
				}
				else
				{
					return std::strlen(name) + SYEX_DERIV_STR_LEN(O - 1) + SYEX_DERIV_APPLIED_EXPR_LEN;
				}
			}
			else
			{
				if constexpr (O == 1)
				{
					return SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(name);
				}
				else
				{
					return SYEX_DIRECTIONAL_DERIV_VAR_LEN(name, O);
				}
			}
		}

	};


	template<size_t O>
	struct print_deriv<O, Axis::NONE, false>
	{
		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		template<typename Sp = int>
		static size_t print(FILE* out, Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				size_t n = 0;
				n += fprintf(out, SYEX_DERIV_STR_1);
				return n;
			}
			else if constexpr (O % 2 == 1)
			{
				size_t n = 0;
				n += fprintf(out, SYEX_DERIV_STR_1);
				n += fprintf(out, SYEX_DERIV_STR_FMT(O - 1));
				return n;
			}
			else
			{
				return fprintf(out, SYEX_DERIV_STR_FMT(O));
			}
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename Sp = int>
		static size_t print(FILE* out, const char* name, Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				size_t n = 0;
                n += fprintf(out, SYEX_DERIV_STR_1);
				n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT, name);
				return n;
			}
			else if constexpr (O % 2 == 1)
			{
				size_t n = 0;
                n += fprintf(out, SYEX_DERIV_STR_1);
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

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		template<typename Sp = int>
		static size_t print(char* out, Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				size_t n = 0;
				n += sprintf(out + n, SYEX_DERIV_STR_1);
				return n;
			}
			else if constexpr (O % 2 == 1)
			{
				size_t n = 0;
				n += sprintf(out + n, SYEX_DERIV_STR_1);
				n += sprintf(out + n, SYEX_DERIV_STR_FMT(O - 1));
				return n;
			}
			else
			{
				return sprintf(out, SYEX_DERIV_STR_FMT(O));
			}
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename Sp = int>
		static size_t print(char* out, const char* name, Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				size_t n = 0;
				n += sprintf(out + n, SYEX_DERIV_STR_1);
				n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT, name);
				return n;
			}
			else if constexpr (O % 2 == 1)
			{
				size_t n = 0;
				n += sprintf(out + n, SYEX_DERIV_STR_1);
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

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename Sp = int>
		static size_t print_length(Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				return SYEX_DERIV_STR_1;
			}
			else if constexpr (O % 2 == 1)
			{
				return SYEX_DERIV_STR_1 + SYEX_DERIV_STR_LEN(O - 1);
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
		 */
		template<typename Sp = int>
		static size_t print_length(const char* name, Sp const& = 0)
		{
			if constexpr (O == 1)
			{
				return SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(name);
			}
			else if constexpr (O % 2 == 1)
			{
				return std::strlen(name) + SYEX_DIRECTIONAL_DERIV_1_LEN
					+ SYEX_DERIV_STR_LEN(O - 1) + SYEX_DERIV_APPLIED_EXPR_LEN;
			}
			else
			{
				return std::strlen(name) + SYEX_DERIV_STR_LEN(O - 1) + SYEX_DERIV_APPLIED_EXPR_LEN;
			}
		}


		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(FILE* out, SymbolicDerivative<G>)
		{
			return fprintf(out, SYEX_DIRECTIONAL_DERIV_STR("%zd", "%s"), O, expr::get_op_name(G{}), O);
		}

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The file to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(FILE* out, const char* name, SymbolicDerivative<G>)
		{
			return fprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_STR("%s", "%zd", "%s"), O, name, expr::get_op_name(G{}), O);
		}


		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(char* out, SymbolicDerivative<G>)
		{
			return sprintf(out, SYEX_DIRECTIONAL_DERIV_STR("%zd", "%s"), O, expr::get_op_name(G{}), O);
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(char* out, const char* name, SymbolicDerivative<G>)
		{
			return sprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_STR("%s", "%zd", "%s"), O, name, expr::get_op_name(G{}), O);
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(SymbolicDerivative<G>)
		{
			return SYEX_DIRECTIONAL_DERIV_LEN(O) + std::strlen(expr::get_op_name(G{}));
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(const char* name, SymbolicDerivative<G>)
		{
			return SYEX_DIRECTIONAL_DERIV_LEN(O) + std::strlen(expr::get_op_name(name))
				+ std::strlen(expr::get_op_name(G{}));
		}
	};

	template<size_t O>
	struct print_deriv<O, Axis::NONE, true>
	{

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(FILE* out, G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_WITH_OTHER_FMT(Axis::NONE, buffer));
			}
			else
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_WITH_OTHER_FMT(O, Axis::NONE, buffer));
			}
			delete[] buffer;
			return n;
		}

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The file to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(FILE* out, const char* name, G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_VAR_WITH_OTHER_FMT(name, Axis::NONE, buffer));
			}
			else
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_WITH_OTHER_FMT(name, O, Axis::NONE, buffer));
			}
			delete[] buffer;
			return n;
		}


		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(char* out, G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_1_WITH_OTHER_FMT(Axis::NONE, buffer));
			}
			else
			{
				n += fprintf(out, SYEX_DIRECTIONAL_DERIV_WITH_OTHER_FMT(O, Axis::NONE, buffer));
			}
			delete[] buffer;
			return n;
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(char* out, const char* name, G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += sprintf(out, SYEX_DIRECTIONAL_DERIV_1_VAR_WITH_OTHER_FMT(name, Axis::NONE, buffer));
			}
			else
			{
				n += sprintf(out, SYEX_DIRECTIONAL_DERIV_VAR_WITH_OTHER_FMT(name, O, Axis::NONE, buffer));
			}
			delete[] buffer;
			return n;
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += SYEX_DIRECTIONAL_DERIV_1_LEN + std::strlen(buffer);
			}
			else
			{
				n += SYEX_DIRECTIONAL_DERIV_LEN(O) + std::strlen(buffer);
			}
			delete[] buffer;
			return n;
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(const char* name, G const& var)
		{
			size_t n = 0;
			auto buffer = set_var_string(var);
			if constexpr (O == 1)
			{
				n += SYEX_DIRECTIONAL_DERIV_1_VAR_LEN(name) + std::strlen(buffer);
			}
			else
			{
				n += SYEX_DIRECTIONAL_DERIV_VAR_LEN(name, O) + std::strlen(buffer);
			}
			delete[] buffer;
			return n;
		}

	};



	template<Axis ax, bool is_directional>
	struct print_deriv<0, ax, is_directional>
	{
		static size_t print(...)
		{
			return 0;
		}

		static size_t print_length(...)
		{
			return 0;
		}
	};



	template<Axis ax, size_t... Os>
	struct print_mixed_deriv;

	template<Axis ax>
	struct print_mixed_deriv<ax>
	{
		static size_t print(...)
		{
			return 0;
		}

		static size_t print_length(...)
		{
			return 0;
		}
	};


	template<Axis ax, size_t O, size_t... Os>
	struct print_mixed_deriv<ax, O, Os...>
	{
		static const Axis next_ax = (ax == Axis::X) ? Axis::Y : (ax == Axis::Y) ? Axis::Z : Axis::NONE;

		using next_type = std::conditional_t<
			(sizeof...(Os) > 1),
			print_mixed_deriv<next_ax, Os...>,
			print_deriv<Os..., next_ax, true>>;

		static size_t print(FILE* out)
		{
			size_t n = print_deriv<O, ax, true>::print(out);
			n += next_type::print(out);
			return n;
		}

		static size_t print(FILE* out, const char* name)
		{
			size_t n = print_deriv<O, ax, true>::print(out);
			n += next_type::print(out, name);
			return n;
		}

		static size_t print(char* out)
		{
			size_t n = print_deriv<O, ax, true>::print(out);
			n += next_type::print(out + n);
			return n;
		}

		static size_t print(char* out, const char* name)
		{
			size_t n = print_deriv<O, ax, true>::print(out);
			n += next_type::print(out + n, name);
			return n;
		}

		static size_t print_length()
		{
			size_t n = print_deriv<O, ax, true>::print_length();
			n += next_type::print_length();
			return n;
		}

		static size_t print_length(const char* name)
		{
			size_t n = print_deriv<O, ax, true>::print_length();
			n += next_type::print_length(name);
			return n;
		}
	};

	template<typename G>
	auto var_str(SymbolicFunctionalDerivative<G> const& solver)
	{
		return get_op_name(G{});
	}

	template<typename G>
	auto var_str(SymbolicFunctionalDerivative<DynamicVariable<NamedData<G>>> const& solver)
	{
		return solver.get_name();
	}

	template<size_t Z, typename G>
	auto var_str(SymbolicFunctionalDerivative<Variable<Z, NamedData<G>>> const& solver)
	{
		return solver.get_name();
	}


	template<typename G0>
	struct print_functional_deriv
	{
	/*	print_functional_deriv(const char* name) : var{ (std::strlen(name) > 0) ? new char[std::strlen(name) + 1] : nullptr } 
		{
			std::strcpy(var, name);
		}
		print_functional_deriv() : print_functional_deriv(get_op_name(G0{})) {}

		print_functional_deriv(print_functional_deriv<G0> const& other) : print_functional_deriv(other.var) {}
		print_functional_deriv(print_functional_deriv<G0>&& other) : print_functional_deriv() 
		{
			swap(*this, other);
		}
		print_functional_deriv<G0>& operator=(print_functional_deriv<G0> other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(print_functional_deriv<G0>& first, print_functional_deriv<G0>& second)
		{
			using std::swap;
			swap(first.var, second.var);
		}*/

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The file to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(FILE* out, SymbolicFunctionalDerivative<G> const& solver)
		{
			return fprintf(out, SYEX_FUNCTIONAL_DERIV_FMT(var_str(solver)));
		}

		//! Print the derivative the given order to a file.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The file to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(FILE* out, SymbolicFunctionalDerivative<G> const& solver, const char* name)
		{
			return fprintf(out, SYEX_FUNCTIONAL_DERIV_VAR_FMT(name, var_str(solver)));
		}


		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order.
		 *
		 * \param out The string to which the derivative is printed.
		 */
		template<typename G>
		static size_t print(char* out, SymbolicFunctionalDerivative<G> const& solver)
		{
			return sprintf(out, SYEX_FUNCTIONAL_DERIV_FMT(var_str(solver)));
		}

		//! Print the derivative the given order to a string.
		/*!
		 * Print the derivative of the given order, and formatted using the order. A
		 * name is also provided, which is the name of the variable to which the derivative
		 * is applied.
		 *
		 * \param out The string to which the derivative is printed.
		 * \param name A string that appears in the numerator after the partial symbol.
		 */
		template<typename G>
		static size_t print(char* out, SymbolicFunctionalDerivative<G> const& solver, const char* name)
		{
			return sprintf(out, SYEX_FUNCTIONAL_DERIV_VAR_FMT(name, var_str(solver)));
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(SymbolicFunctionalDerivative<G> const& solver)
		{
			return SYEX_FUNCTIONAL_DERIV_LEN(var_str(solver));
		}

		//! Get the print length of the derivative output string.
		/*!
		 * Returns the number of characters in the format string that is
		 * printed to display the derivative. Only includes characters that
		 * are printed as part of the format, and not substituted expression
		 * strings.
		 */
		template<typename G>
		static size_t print_length(SymbolicFunctionalDerivative<G> const& solver, const char* name)
		{
			return SYEX_FUNCTIONAL_DERIV_VAR_LEN(name, var_str(solver));
		}
	};

	template<size_t O>
	auto select_print_deriv(std::index_sequence<O>)
	{
		return print_deriv<O, Axis::NONE, true>{};
	}

	template<typename Sp, Axis ax, size_t O>
	auto select_print_deriv(typename Solver<Sp>::template derivative<ax, O>)
	{
		using Dd = typename Solver<Sp>::template derivative<ax, O>;
		return print_deriv<O, ax, Dd::is_directional>{};
	}

	template<typename Sp, Axis ax, size_t O>
	auto select_print_deriv(typename Solver<Sp>::template directional_derivative<ax, O>)
	{
		return print_deriv<O, ax, true>{};
	}

	template<typename Sp, size_t O1, size_t O2>
	auto select_print_deriv(typename Solver<Sp>::template mixed_derivative<O1, O2>)
	{
		return print_mixed_deriv<Axis::X, O1, O2>{};
	}

	template<typename Sp, size_t O1, size_t O2, size_t O3>
	auto select_print_deriv(typename Solver<Sp>::template mixed_derivative<O1, O2, O3>)
	{
		return print_mixed_deriv<Axis::X, O1, O2, O3>{};
	}

	template<typename G>
	auto select_print_deriv(SymbolicFunctionalDerivative<G> const& solver)
	{
		return print_functional_deriv<G>{};
	}
		
}

// ************************************************************************************

namespace expr
{

	template<typename T>
	auto get_fourier_name(T const& t)
	{
		char* name = expr::get_op_name(std::forward<T>(t));
		return std::string(SYEX_FT_OF_OP_FMT_A) + std::string(name) + std::string(SYEX_FT_OF_OP_FMT_B);
	}



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
	size_t printe(E const& expr, const char* fmt = "", ...)
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
	size_t fprinte(E const& expr, FILE* out, const char* fmt = "", ...)
	{
		size_t n = 0;
		if (*fmt)
		{
			va_list list;
			va_start(list, fmt);

			char buffer[LINE_READ_BUFFER];
			vsnprintf(buffer, LINE_READ_BUFFER, fmt, list);
			n += fprintf(out, "%s " SYEX_EQN_SEP " ", buffer);

			va_end(list);
		}

		n += expr.print(out);
		n += fprintf(out, "\n");
		return n;
	}

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
	size_t sprinte(E const& expr, char* out, const char* fmt = "", ...)
	{
		size_t n = 0;
		if (*fmt)
		{
			va_list list;
			va_start(list, fmt);

			char buffer[LINE_READ_BUFFER];
			vsnprintf(buffer, LINE_READ_BUFFER, fmt, list);
			n += sprintf(out + n, "%s " SYEX_EQN_SEP " ", buffer);

			va_end(list);
		}

		n += expr.print(out + n);
		n += sprintf(out + n, "\n");
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
	template<typename T, size_t D>
	size_t coeff_print_length(any_vector_t<T, D> const& value);
	template<typename T, typename I>
	size_t coeff_print_length(OpCoeff<T, I> const& value);

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

	template<size_t N, size_t D>
	inline size_t coeff_print_length(OpFractionLiteral<N, D>)
	{
		return OpFractionLiteral<N, D>{}.print_length();
	}

	template<size_t N, size_t D>
	inline size_t coeff_print_length(OpNegFractionLiteral<N, D>)
	{
		return OpNegFractionLiteral<N, D>{}.print_length();
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
		n += STR_ARR_LEN(SYEX_VECTOR_FMT_A SYEX_VECTOR_FMT_B)
			+ (D - 1) * STR_ARR_LEN(SYEX_VECTOR_FMT_SEP);
		for (iter_type i = 0; i < D; ++i)
		{
			n += coeff_print_length(value[i]);
		}
		return n;
	}


	template<typename T, typename I>
	size_t coeff_print_length(OpCoeff<T, I> const& value)
	{
		return coeff_print_length(expr::make_literal(value.eval()));
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

	template<typename T>
	size_t print_tensor_entries_1(FILE* out, T const& value, size_t P0, size_t P1, size_t N, size_t M)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_TENSOR_FMT_A);

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < M; ++j)
			{
				if (i == P0 && j == P1)
				{
					n += make_literal(expr::eval(value)).print(out);
				}
				else
				{
					n += fprintf(out, SYEX_TENSOR_EMPTY_FMT);
				}

				if (j < M - 1)
				{
					n += fprintf(out, SYEX_TENSOR_ROW_SEP_FMT);
				}
			}
			if (i < N - 1)
			{
				n += fprintf(out, SYEX_TENSOR_COLUMN_SEP_FMT);
			}
		}
		n += fprintf(out, SYEX_TENSOR_FMT_B);
		return n;
	}


	template<typename ostream_t, typename T, size_t P0, size_t P1, size_t N, size_t M>
	size_t print_tensor(ostream_t* out, T const& value, std::index_sequence<P0, P1>, std::index_sequence<N, M>)
	{
		return print_tensor_entries_1(out, value, P0, P1, N, M);
	}

	template<typename ostream_t, typename T, size_t P0, size_t N>
	size_t print_tensor(ostream_t* out, T const& value, std::index_sequence<P0>, std::index_sequence<N>)
	{
		return print_tensor_entries_1(out, value, P0, 0, N, 1);
	}

	template<typename T>
	size_t print_tensor_length_1(T const& value, size_t P0, size_t P1, size_t N, size_t M)
	{
		size_t n = 0;
		n += STR_ARR_LEN(SYEX_TENSOR_FMT_A);

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < M; ++j)
			{
				if (i == P0 && j == P1)
				{
					n += make_literal(expr::eval(value)).print_length();
				}
				else
				{
					n += STR_ARR_LEN(SYEX_TENSOR_EMPTY_FMT);
				}

				if (j < M - 1)
				{
					n += STR_ARR_LEN(SYEX_TENSOR_ROW_SEP_FMT);
				}
			}
			if (i < N - 1)
			{
				n += STR_ARR_LEN(SYEX_TENSOR_COLUMN_SEP_FMT);
			}
		}
		n += STR_ARR_LEN(SYEX_TENSOR_FMT_B);
		return n;
	}


	template<typename T, size_t P0, size_t P1, size_t N, size_t M>
	size_t print_tensor_length(T const& value, std::index_sequence<P0, P1>, std::index_sequence<N, M>)
	{
		return print_tensor_length_1(value, P0, P1, N, M);
	}

	template<typename T, size_t P0, size_t N>
	size_t print_tensor_length(T const& value, std::index_sequence<P0>, std::index_sequence<N>)
	{
		return print_tensor_length_1(value, P0, 0, N, 1);
	}


	template<typename T, size_t... Ns>
	size_t coeff_print_length(OpTensor<T, Ns...> const& value)
	{
		using p_seq = symphas::lib::seq_join_t<symphas::lib::types_before_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;
		using n_seq = symphas::lib::seq_join_t<symphas::lib::types_after_at_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;
		return print_tensor_length(T(value), p_seq{}, n_seq{});
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
	template<typename T, size_t D>
	size_t print_with_coeff(char* out, const char* expr, any_vector_t<T, D> const& value);
	template<typename T, typename I>
	size_t print_with_coeff(char* out, const char* expr, OpCoeff<T, I> const& value);

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

	template<typename T, size_t... Ns>
	size_t print_with_coeff(char* out, const char* expr, OpTensor<T, Ns...> const& value)
	{
		using p_seq = symphas::lib::seq_join_t<symphas::lib::types_before_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;
		using n_seq = symphas::lib::seq_join_t<symphas::lib::types_after_at_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;

		size_t n = print_tensor(out, T(value), p_seq{}, n_seq{});
		n += sprintf(out + n, "%s", expr);
		return n;
	}

	template<typename... Ts>
	size_t print_with_coeff(char* out, const char* expr, OpAdd<Ts...> const& value)
	{
		return print_with_coeff(out, expr, value.eval());
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

	template<typename T, size_t D>
	size_t print_with_coeff(char* out, const char* expr, any_vector_t<T, D> const& value)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_VECTOR_FMT_A);
		for (iter_type i = 0; i < D; ++i)
		{
			n += make_literal(value[i]).print(out + n);
			if (i < D - 1)
			{
				n += sprintf(out + n, SYEX_VECTOR_FMT_SEP);
			}
		}
		n += sprintf(out + n, SYEX_VECTOR_FMT_B "%s", expr);
		return n;
	}

	template<typename T, typename I>
	size_t print_with_coeff(char* out, const char* expr, OpCoeff<T, I> const& value)
	{
		return print_with_coeff(out, expr, expr::make_literal(value.eval()));
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
	template<typename T, size_t D>
	size_t print_with_coeff(FILE* out, const char* expr, any_vector_t<T, D> const& value);
	template<typename T, typename I>
	size_t print_with_coeff(FILE* out, const char* expr, OpCoeff<T, I> const& value);

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
	 * \param out The file to put the printed expression.  * \param expr The expression which is printed.
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

	template<typename T, size_t... Ns>
	size_t print_with_coeff(FILE* out, const char* expr, OpTensor<T, Ns...> const& value)
	{
		using p_seq = symphas::lib::seq_join_t<symphas::lib::types_before_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;
		using n_seq = symphas::lib::seq_join_t<symphas::lib::types_after_at_index<sizeof...(Ns) / 2, std::index_sequence<Ns>...>>;

		size_t n = print_tensor(out, T(value), p_seq{}, n_seq{});
		n += fprintf(out, "%s", expr);
		return n;
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

	template<typename T, typename I>
	size_t print_with_coeff(FILE* out, const char* expr, OpCoeff<T, I> const& value)
	{
		return print_with_coeff(out, expr, expr::make_literal(value.eval()));
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

	template<size_t L, size_t L0 = L - 1,
#ifdef LATEX_PLOT
		size_t L1 = 1
#else
		size_t L1 = 0
#endif
	>
	expr_name_arr<0> print_with_subscript(const char(&term)[L], const char* subscript)
	{

		size_t LL = std::strlen(subscript);
		expr_name_arr out(len_type(L0 + L1 + LL + 3));
		for (iter_type i = 0; i < L0; ++i)
		{
			out.value[i] = term[i];
		}

		if constexpr (L1)
		{
			out.value[L0] = '_';
		}
		out.value[L0 + L1] = '{';

		std::strcat(out.value, subscript);

		size_t LN = L0 + L1 + 1;

		out.value[LN + LL] = '}';
		out.value[LN + LL + 1] = '\0';

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

	template<typename E0, typename... Es>
	auto print_list(E0&& e0, Es&&... es)
	{
		const char* names[]{ expr::get_op_name(std::forward<E0>(e0)), expr::get_op_name(std::forward<Es>(es))... };
		len_type len = 0;
		for (const char* name : names)
		{
			len += (len_type)std::strlen(name) + 1;
		}

		expr_name_arr list_str(len);
		for (const char* name : names)
		{
			std::strcpy(list_str, name);
			if (name != names[sizeof...(Es)])
			{
				std::strcat(list_str, ",");
			}
		}

		return list_str;
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


	template<typename T>
	struct symbolic_eval_print
	{
		template<typename E>
		size_t operator()(FILE* out, T const&, E const& e)
		{
			return e.print(out);
		}

		template<typename E>
		size_t operator()(char* out, T const&, E const& e)
		{
			return e.print(out);
		}

		template<typename E>
		size_t operator()(T const&, E const& e)
		{
			return e.print_length();
		}
	};

	template<expr::NoiseType nt>
	struct noise_name_print
	{
		auto operator()(FILE* out)
		{
			return fprintf(out, "%s", SYEX_NOISE_TOKEN_NONE);
		}

		auto operator()(char* out)
		{
			return sprintf(out, "%s", SYEX_NOISE_TOKEN_NONE);
		}

		auto operator()()
		{
			return STR_ARR_LEN(SYEX_NOISE_TOKEN_NONE);
		}
	};

	template<>
	struct noise_name_print<expr::NoiseType::POISSON>
	{
		auto operator()(FILE* out)
		{
			return fprintf(out, "%s", SYEX_NOISE_TOKEN_POISSON);
		}

		auto operator()(char* out)
		{
			return sprintf(out, "%s", SYEX_NOISE_TOKEN_POISSON);
		}

		auto operator()()
		{
			return STR_ARR_LEN(SYEX_NOISE_TOKEN_POISSON);
		}
	};

	template<>
	struct noise_name_print<expr::NoiseType::WHITE>
	{
		auto operator()(FILE* out)
		{
			return fprintf(out, "%s", SYEX_NOISE_TOKEN_WHITE);
		}

		auto operator()(char* out)
		{
			return sprintf(out, "%s", SYEX_NOISE_TOKEN_WHITE);
		}

		auto operator()()
		{
			return STR_ARR_LEN(SYEX_NOISE_TOKEN_WHITE);
		}
	};

	template<expr::NoiseType nt, typename T, size_t D>
	struct symbolic_eval_print<NoiseData<nt, T, D>>
	{
		template<typename E>
		size_t operator()(FILE* out, NoiseData<nt, T, D> const&, E const& e)
		{
			size_t n = 0;
			n += noise_name_print<nt>{}(out);
			n += fprintf(out, "%s", SYEX_NOISE_A);
			n += e.print(out);
			n += fprintf(out, "%s", SYEX_NOISE_B);
			return n;
		}

		template<typename E>
		size_t operator()(char* out, NoiseData<nt, T, D> const&, E const& e)
		{
			size_t n = 0;
			n += noise_name_print<nt>{}(out + n);
			n += sprintf(out + n, "%s", SYEX_NOISE_A);
			n += e.print(out + n);
			n += sprintf(out + n, "%s", SYEX_NOISE_B);
			return n;
		}

		template<typename E>
		size_t operator()(NoiseData<nt, T, D> const&, E const& e)
		{
			return noise_name_print<nt>{}() + STR_ARR_LEN(SYEX_NOISE_A SYEX_NOISE_B)
				+ e.print_length();
		}
	};


	template<typename E0>
	struct symbolic_eval_print<SymbolicListIndex<E0>>
	{
		template<typename E>
		size_t operator()(FILE* out, SymbolicListIndex<E0> const& data, E const& e)
		{
			size_t n = 0;
			n += fprintf(out, "%s", SYEX_ARRAY_A);
			n += e.print(out);
			n += fprintf(out, "%s", SYEX_ARRAY_B SYEX_ARRAY_SUBSCRIPT_A);
			n += data.e.print(out);
			n += fprintf(out, "%s", SYEX_ARRAY_SUBSCRIPT_B);
			return n;
		}

		template<typename E>
		size_t operator()(char* out, SymbolicListIndex<E0> const& data, E const& e)
		{
			size_t n = 0;
			n += sprintf(out + n, "%s", SYEX_ARRAY_A);
			n += e.print(out + n);
			n += sprintf(out + n, "%s", SYEX_ARRAY_B SYEX_ARRAY_SUBSCRIPT_A);
			n += data.e.print(out + n);
			n += sprintf(out + n, "%s", SYEX_ARRAY_SUBSCRIPT_B);
			return n;
		}

		template<typename E>
		size_t operator()(SymbolicListIndex<E0> const& data, E const& e)
		{
			return data.e.print_length() + e.print_length()
				+ STR_ARR_LEN(SYEX_ARRAY_A SYEX_ARRAY_B SYEX_ARRAY_SUBSCRIPT_A SYEX_ARRAY_SUBSCRIPT_B);
		}
	};

	namespace
	{


		template<typename E>
		size_t print_limit(FILE* out, limit_data<E> const& e, bool upper = false);
		inline size_t print_limit(FILE* out, limit_data<OpVoid> const& e, bool upper = false);
		template<typename E>
		size_t print_limit(char* out, limit_data<E> const& e, bool upper = false);
		inline size_t print_limit(char* out, limit_data<OpVoid> const& e, bool upper = false);
		template<typename E>
		size_t print_limit_length(limit_data<E> const& e, bool upper = false);
		inline size_t print_limit_length(limit_data<OpVoid> const& e, bool upper = false);


		size_t print_limit(FILE* out, int l, bool upper = false)
		{
			return fprintf(out, "%d", l);
		}

		template<typename E0>
		size_t print_limit(FILE* out, OpExpression<E0> const& e, bool upper = false)
		{
			return static_cast<E0 const*>(&e)->print(out);
		}

		template<typename V, int N, int P>
		size_t print_limit(FILE* out, OpTerm<V, expr::symbols::i_<N, P>>, bool upper = false)
		{
			if constexpr (P > 0)
			{
				return fprintf(out, "%s+%d", expr::get_op_name(expr::symbols::i_<N, 0>{}), P);
			}
			else if constexpr (P < 0)
			{
				return fprintf(out, "%s-%d", expr::get_op_name(expr::symbols::i_<N, 0>{}), -P);
			}
			else
			{
				return fprintf(out, "%s", expr::get_op_name(expr::symbols::i_<N, 0>{}));
			}
		}

		template<typename... Es, size_t I0, size_t... Is>
		size_t print_limit(FILE* out, OpAdd<Es...> const& e, std::index_sequence<I0, Is...>, bool upper = false)
		{
			size_t n = 0;
			n += print_limit(out, expr::get<I0>(e));
			
			if constexpr (sizeof...(Is) > 0)
			{
				if (expr::eval(expr::coeff(expr::get<I0 + 1>(e))) >= 0)
				{
					n += symphas::internal::print_sep(out, SYEX_ADD_SEP);
				}
				n += print_limit(out, e, std::index_sequence<Is...>{});
			}
			
			return n;
		}

		template<typename... Es>
		size_t print_limit(FILE* out, OpAdd<Es...> const& e, bool upper = false)
		{
			return print_limit(out, get_limit_data(e));
		}

		template<typename L, typename R>
		size_t print_limit(FILE* out, std::pair<L, R> const& e, bool upper = false)
		{
			size_t n = 0;
			auto _a = get_limit_data(e.first);
			auto _b = get_limit_data(e.second);

			if (!_a || !_b)
			{
				n += fprintf(out, "%s", SYEX_SUM_LIM_COMPARE_A);
				n += print_limit(out, e.first);
				n += fprintf(out, "%s", SYEX_SUM_LIM_COMPARE_SEP);
				n += print_limit(out, e.second);
				n += fprintf(out, "%s", SYEX_SUM_LIM_COMPARE_B);
			}
			else
			{
				if (upper)
				{
					if (_a > _b)
					{
						n += print_limit(out, _b);
					}
					else
					{
						n += print_limit(out, _a);
					}
				}
				else
				{
					if (_a > _b)
					{
						n += print_limit(out, _a);
					}
					else
					{
						n += print_limit(out, _b);
					}
				}
			}
			return n;
		}


		size_t print_limit(char* out, int l, bool upper = false)
		{
			return sprintf(out, "%d", l);
		}

		template<typename E0>
		size_t print_limit(char* out, OpExpression<E0> const& e, bool upper = false)
		{
			return static_cast<E0 const*>(&e)->print(out);
		}

		template<typename V, int N, int P>
		size_t print_limit(char* out, OpTerm<V, expr::symbols::i_<N, P>>, bool upper = false)
		{
			if constexpr (P > 0)
			{
				return sprintf(out, "%s+%d", expr::get_op_name(expr::symbols::i_<N, 0>{}), P);
			}
			else if constexpr (P < 0)
			{
				return sprintf(out, "%s-%d", expr::get_op_name(expr::symbols::i_<N, 0>{}), -P);
			}
			else
			{
				return sprintf(out, "%s", expr::get_op_name(expr::symbols::i_<N, 0>{}));
			}
		}

		template<typename... Es, size_t I0, size_t... Is>
		size_t print_limit(char* out, OpAdd<Es...> const& e, std::index_sequence<I0, Is...>, bool upper = false)
		{
			size_t n = 0;
			n += print_limit(out, expr::get<I0>(e));

			if constexpr (sizeof...(Is) > 0)
			{
				if (expr::eval(expr::coeff(expr::get<I0 + 1>(e))) >= 0)
				{
					n += symphas::internal::print_sep(out + n, SYEX_ADD_SEP);
				}
				n += print_limit(out + n, e, std::index_sequence<Is...>{});
			}

			return n;
		}

		template<typename... Es>
		size_t print_limit(char* out, OpAdd<Es...> const& e, bool upper = false)
		{
			return print_limit(out, get_limit_data(e));
		}

		template<typename L, typename R>
		size_t print_limit(char* out, std::pair<L, R> const& e, bool upper = false)
		{
			size_t n = 0;
			n += sprintf(out + n, "%s", SYEX_SUM_LIM_COMPARE_A);
			n += print_limit(out + n, e.first);
			n += sprintf(out + n, "%s", SYEX_SUM_LIM_COMPARE_SEP);
			n += print_limit(out + n, e.second);
			n += sprintf(out + n, "%s", SYEX_SUM_LIM_COMPARE_B);
			return n;
		}

		template<int N, int P, typename T1, typename T2>
		size_t print_limit(expr::symbols::i_<N, P>, FILE* out, expr::series_limits<T1, T2> const& limit, bool upper = false)
		{
			size_t n = 0;
			n += fprintf(out, "%s%s=", SYEX_SUM_LIM_A, expr::get_op_name(expr::symbols::i_<N, P>{}));
			n += print_limit(out, expr::limit_0(limit));
			n += fprintf(out, "%s", SYEX_SUM_LIM_SEP);
			n += print_limit(out, expr::limit_1(limit));
			n += fprintf(out, "%s", SYEX_SUM_LIM_B);
			return n;
		}

		template<int N, int P, typename T1, typename T2>
		size_t print_limit(expr::symbols::i_<N, P>, char* out, expr::series_limits<T1, T2> const& limit, bool upper = false)
		{
			size_t n = 0;
			n += sprintf(out, "%s%s=", SYEX_SUM_LIM_A, expr::get_op_name(expr::symbols::i_<N, P>{}));
			n += print_limit(out + n, expr::limit_0(limit));
			n += sprintf(out + n, "%s", SYEX_SUM_LIM_SEP);
			n += print_limit(out + n, expr::limit_1(limit));
			n += sprintf(out + n, "%s", SYEX_SUM_LIM_B);
			return n;
		}

		size_t print_limit_length(int l, bool upper = false)
		{
			return symphas::lib::num_digits(l) + ((l < 0) ? 1 : 0);
		}

		template<typename E0>
		size_t print_limit_length(OpExpression<E0> const& e, bool upper = false)
		{
			return static_cast<E0 const*>(&e)->print_length();
		}

		template<typename V, int N, int P>
		size_t print_limit_length(OpTerm<V, expr::symbols::i_<N, P>>, bool upper = false)
		{
			if constexpr (P > 0)
			{
				return 1 + std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{})) + symphas::lib::num_digits(P);
			}
			else if constexpr (P < 0)
			{
				return 2 + std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{})) + symphas::lib::num_digits(P);
			}
			else
			{
				return std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{}));
			}
		}

		template<typename... Es, size_t I0, size_t... Is>
		size_t print_limit_length(OpAdd<Es...> const& e, std::index_sequence<I0, Is...>, bool upper = false)
		{
			size_t n = 0;
			n += print_limit_length(expr::get<I0>(e));

			if constexpr (sizeof...(Is) > 0)
			{
				if (expr::eval(expr::coeff(expr::get<I0 + 1>(e))) >= 0)
				{
					n += STR_ARR_LEN(SYEX_ADD_SEP);
				}
				n += print_limit_length(e, std::index_sequence<Is...>{});
			}

			return n;
		}

		template<typename... Es>
		size_t print_limit_length(OpAdd<Es...> const& e, bool upper = false)
		{
			return print_limit_length(get_limit_data(e));
		}

		template<typename L, typename R>
		size_t print_limit_length(std::pair<L, R> const& e, bool upper = false)
		{
			size_t n = 0;
			n += STR_ARR_LEN(SYEX_SUM_LIM_COMPARE_A SYEX_SUM_LIM_COMPARE_SEP SYEX_SUM_LIM_COMPARE_B);
			n += print_limit_length(e.first);
			n += print_limit_length(e.second);
			return n;
		}

		template<int N, int P>
		size_t print_limit_length(expr::symbols::i_<N, P>, bool upper = false)
		{
			size_t n = 0;
			if constexpr (P > 0)
			{
				n += std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{}));
				n += 1 + symphas::lib::num_digits(P);
			}
			else if constexpr (P < 0)
			{
				n += std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{}));
				n += 1 + symphas::lib::num_digits(-P);
			}
			else
			{
				n += std::strlen(expr::get_op_name(expr::symbols::i_<N, 0>{}));
			}
			return n;
		}

		template<int N, int P, typename T1, typename T2>
		size_t print_limit_length(expr::symbols::i_<N, P>, expr::series_limits<T1, T2> const& limit)
		{
			size_t n = 0;
			n += STR_ARR_LEN(SYEX_SUM_LIM_A SYEX_SUM_LIM_SEP SYEX_SUM_LIM_B);
			n += 1 + std::strlen(expr::get_op_name(expr::symbols::i_<N, P>{}));
			n += print_limit_length(limit._0, true);
			n += print_limit_length(limit._1, false);
			return n;
		}


		template<typename E>
		size_t print_limit(FILE* out, limit_data<E> const& e, bool upper)
		{
			size_t n = 0;

			n += (e.e + expr::make_literal(e.index)).print(out);
			if (e.offset > 0)
			{
				n += fprintf(out, SYEX_LIMIT_OFFSET_A);
				n += print_limit(out, e.offset);
				n += fprintf(out, SYEX_LIMIT_OFFSET_B);
			}

			return n;
		}

		inline size_t print_limit(FILE* out, limit_data<OpVoid> const& e, bool upper)
		{
			size_t n = 0;

			n += print_limit(out, e.index);
			if (e.offset > 0)
			{
				n += fprintf(out, SYEX_LIMIT_OFFSET_A);
				n += print_limit(out, e.offset);
				n += fprintf(out, SYEX_LIMIT_OFFSET_B);
			}

			return n;
		}

		template<typename E>
		size_t print_limit(char* out, limit_data<E> const& e, bool upper)
		{
			size_t n = 0;

			n += (e.e + expr::make_literal(e.index)).print(out + n);
			if (e.offset > 0)
			{
				n += sprintf(out + n, SYEX_LIMIT_OFFSET_A);
				n += print_limit(out + n, e.offset);
				n += sprintf(out + n, SYEX_LIMIT_OFFSET_B);
			}

			return n;
		}

		inline size_t print_limit(char* out, limit_data<OpVoid> const& e, bool upper)
		{
			size_t n = 0;

			n += print_limit(out + n, e.index);
			if (e.offset > 0)
			{
				n += sprintf(out + n, SYEX_LIMIT_OFFSET_A);
				n += print_limit(out + n, e.offset);
				n += sprintf(out + n, SYEX_LIMIT_OFFSET_B);
			}

			return n;
		}


		template<typename E>
		size_t print_limit_length(limit_data<E> const& e, bool upper)
		{
			size_t n = 0;

			n += (e.e + expr::make_literal(e.index)).print_length();
			if (e.offset > 0)
			{
				n += STR_ARR_LEN(SYEX_LIMIT_OFFSET_A SYEX_LIMIT_OFFSET_B);
				n += print_limit_length(e.offset);
			}

			return n;
		}

		inline size_t print_limit_length(limit_data<OpVoid> const& e, bool upper)
		{
			size_t n = 0;

			n += print_limit_length(e.index);
			if (e.offset > 0)
			{
				n += STR_ARR_LEN(SYEX_LIMIT_OFFSET_A SYEX_LIMIT_OFFSET_B);
				n += print_limit_length(e.offset);
			}

			return n;
		}
	}

	template<typename sub_t, typename E, int... I0s, int... P0s, typename... T1s, typename... T2s, typename B>
	struct symbolic_eval_print<
		SymbolicSeries<expr::sum_op, sub_t,
		symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		B>>>
	{
		using sum_t = SymbolicSeries<expr::sum_op, sub_t,
			symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			B>>;


		template<size_t... Ns>
		size_t print_limits(FILE* out, std::tuple<expr::series_limits<T1s, T2s>...> const& limits, std::index_sequence<Ns...>)
		{
			return (print_limit(symphas::lib::type_at_index<Ns, expr::symbols::i_<I0s, P0s>...>{}, out, std::get<Ns>(limits)) + ...);
		}

		template<size_t... Ns>
		size_t print_limits(char* out, std::tuple<expr::series_limits<T1s, T2s>...> const& limits, std::index_sequence<Ns...>)
		{
			size_t n = 0;
			((n += print_limit(symphas::lib::type_at_index<Ns, expr::symbols::i_<I0s, P0s>...>{}, out + n, std::get<Ns>(limits))), ...);
			return n;
		}

		template<size_t... Ns>
		size_t print_limits_length(std::tuple<expr::series_limits<T1s, T2s>...> const& limits, std::index_sequence<Ns...>)
		{
			size_t n = 0;
			((n += print_limit_length(symphas::lib::type_at_index<Ns, expr::symbols::i_<I0s, P0s>...>{}, std::get<Ns>(limits))), ...);
			return n;
		}

		template<typename E0>
		size_t operator()(FILE* out, sum_t const& sum, E0 const& e)
		{
			size_t n = 0;
			n += fprintf(out, SYEX_SUM_SYMBOL);
			n += print_limits(out, sum.limits, std::make_index_sequence<sizeof...(I0s)>{});
			n += fprintf(out, SYEX_SUM_A);
			n += e.print(out);
			n += fprintf(out, SYEX_SUM_B);
			return n;
		}

		template<typename E0>
		size_t operator()(char* out, sum_t const& sum, E0 const& e)
		{
			size_t n = 0;
			n += sprintf(out, SYEX_SUM_SYMBOL);
			n += print_limits(out + n, sum.limits, std::make_index_sequence<sizeof...(I0s)>{});
			n += sprintf(out + n, SYEX_SUM_A);
			n += e.print(out + n);
			n += sprintf(out + n, SYEX_SUM_B);
			return n;
		}

		template<typename E0>
		size_t operator()(sum_t const& sum, E0 const& e)
		{
			size_t n = 0;
			n += STR_ARR_LEN(SYEX_SUM_SYMBOL SYEX_SUM_A SYEX_SUM_B);
			n += print_limits_length(sum.limits, std::make_index_sequence<sizeof...(I0s)>{});
			n += e.print_length();
			return n;
		}
	};

	symbolic_eval_print()->symbolic_eval_print<void>;


	template<typename T>
	struct integral_print
	{
		template<typename E, typename T0>
		size_t operator()(FILE* out, T const& domain, E const& e)
		{
			size_t n = 0;
			n += fprintf(out, SYEX_INTEGRAL_SYMBOL);
			n += fprintf(out, SYEX_INTEGRAL_LIM_A "%s" SYEX_INTEGRAL_LIM_SEP "%s" SYEX_INTEGRAL_LIM_B, "", "");
			n += fprintf(out, expr::get_op_name(domain));
			n += fprintf(out, SYEX_INTEGRAL_A);
			n += e.print(out);
			n += fprintf(out, SYEX_INTEGRAL_B);
			return n;
		}

		template<typename E>
		size_t operator()(char* out, T const& domain, E const& e)
		{
			size_t n = 0;
			n += sprintf(out + n, SYEX_INTEGRAL_SYMBOL);
			n += sprintf(out + n, SYEX_INTEGRAL_LIM_A "%s" SYEX_INTEGRAL_LIM_SEP "%s" SYEX_INTEGRAL_LIM_B, "", "");
			n += fprintf(out, expr::get_op_name(domain));
			n += sprintf(out + n, SYEX_INTEGRAL_A);
			n += e.print(out + n);
			n += sprintf(out + n, SYEX_INTEGRAL_B);
			return n;
		}

		template<typename E>
		size_t operator()(symphas::grid_info const&, E const& e)
		{
			size_t n = 0;
			n += STR_ARR_LEN(SYEX_INTEGRAL_SYMBOL SYEX_INTEGRAL_LIM_A SYEX_INTEGRAL_LIM_SEP SYEX_INTEGRAL_LIM_B
				SYEX_INTEGRAL_DOMAIN_SYM SYEX_INTEGRAL_INTEGRATION_SYM SYEX_INTEGRAL_A SYEX_INTEGRAL_B);
			n += e.print_length();
			return n;
		}
	};

	template<typename T>
	struct integral_print<expr::variational_t<T>>
	{
		template<typename E>
		size_t operator()(FILE* out, symphas::grid_info const&, E const& e)
		{
			size_t n = 0;
			n += fprintf(out, SYEX_INTEGRAL_SYMBOL);
			n += fprintf(out, SYEX_INTEGRAL_LIM_A "%s" SYEX_INTEGRAL_LIM_SEP "%s" SYEX_INTEGRAL_LIM_B, SYEX_INTEGRAL_DOMAIN_SYM, "");
			n += fprintf(out, SYEX_INTEGRAL_INTEGRATION_SYM);
			n += fprintf(out, SYEX_INTEGRAL_A);
			n += e.print(out);
			n += fprintf(out, SYEX_INTEGRAL_B);
			return n;
		}

		template<typename E>
		size_t operator()(char* out, symphas::grid_info const&, E const& e)
		{
			size_t n = 0;
			n += sprintf(out + n, SYEX_INTEGRAL_SYMBOL);
			n += sprintf(out + n, SYEX_INTEGRAL_LIM_A "%s" SYEX_INTEGRAL_LIM_SEP "%s" SYEX_INTEGRAL_LIM_B, SYEX_INTEGRAL_DOMAIN_SYM, "");
			n += sprintf(out + n, SYEX_INTEGRAL_INTEGRATION_SYM);
			n += sprintf(out + n, SYEX_INTEGRAL_A);
			n += e.print(out + n);
			n += sprintf(out + n, SYEX_INTEGRAL_B);
			return n;
		}

		template<typename E>
		size_t operator()(symphas::grid_info const&, E const& e)
		{
			size_t n = 0;
			n += STR_ARR_LEN(SYEX_INTEGRAL_SYMBOL SYEX_INTEGRAL_LIM_A SYEX_INTEGRAL_LIM_SEP SYEX_INTEGRAL_LIM_B
				SYEX_INTEGRAL_DOMAIN_SYM SYEX_INTEGRAL_INTEGRATION_SYM SYEX_INTEGRAL_A SYEX_INTEGRAL_B);
			n += e.print_length();
			return n;
		}
	};

}

#undef SYEX_EQN_SEP
#undef SYEX_PRINT_FUNCTION_STRUCT
#undef SYEX_DEFINE_FUNCTION_STRUCT_DELEGATE
#undef SYEX_PRINT_FUNCTION_STRUCT

#else

namespace expr
{
	template<size_t... Os, typename... T>
	expr_name_arr<1> print_with_subscript(T&&...)
	{
		return "";
	}

	template<size_t... Os, typename... T>
	expr_name_arr<1> print_with_subscripts(T&&...)
	{
		return "";
	}

	//! Empty implementation of printing an expression.
	/*!
	 * Empty implementation for the print function when #PRINTABLE_EQUATIONS is
	 * not used. Therefore, the print function can still be invoked but there is no effect.
	 */
	template<typename... T>
	auto printe(T&&...) {}


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

	template<typename... T>
	auto get_fourier_name(T&&...)
	{
		return "";
	}

	template<typename E0, typename... Es>
	auto print_list(E0&& e0, Es&&... es)
	{
		expr_name_arr list_str(0);
		return list_str;
	}

}

#endif

