
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
 * PURPOSE: Defines elements that can allow definition of additional symbols
 * that can be used in the symbolic algebra framework, as well as type traits
 * that define rules and factoring of expressions. 
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"





/*! 
 * \defgroup expraddsyms Adding Symbols to the Expression Library
 * \ingroup Op
 *
 * Definitions used in adding a new object to structs used in accessing an objects data.
 * Adds a new object to the expression algorithm; allows the object to be given a name and its data accessed.
 *
 * - The `TEMPLATES` argument is the argument to the template parameters; the first argument may
 * be empty if the object is not parameterized. It is always surrounded in parentheses.
 * - The `TYPES` argument second argument is the object type, it is also surrounded in parentheses.
 * - The `DATA_NAME` parameter is the name of the data member within the object.
 * - The `RETURN` parameter depends on the context of the macro, but is typically
 * the data element corresponding to the new data object being added.
 * 
 * Symbols can also have the rules of symbolic algebra modified for them, 
 * such as multiplication rules and whether or not they should be combined
 * as like terms.
 * @{
 */

//! Define a new Symbol in the symbolic algebra functionality.
/*!
 * The symbol is added so as to be printed by the symbolic algebra
 * functionality. 
 * 
 * The macro wraps a specialization which defines a function to return the
 * name of the data element (the new symbol being added). The parameter
 * name to the function is `data`. The function body is defined in the 
 * parameter \p RETURN and must return something that can be interpreted by
 * expr::get_op_name.
 * 
 * For example, the following will add a new symbol using the existing class
 * `new_symbol`, templated on type `T`:
 * 
 * ```cpp
 * DEFINE_SYMBOL_ID((typename T), (new_symbol<T>), return "S")
 * ```
 * 
 * \param TEMPLATES The templates which are passed to the new
 * symbol type in the specialization. The argument is always surrounded
 * by parentheses.
 * \param TYPE The new type, which is provided with the template
 * arguments that are given in \p TEMPLATES. The argument is always surrounded
 * by parentheses.
 * \param ... The body of the function which returns the symbol identifier.
 */
#define DEFINE_SYMBOL_ID(TEMPLATES, TYPE, ...) \
template<SINGLE_ARG TEMPLATES> \
struct expr::SymbolID<SINGLE_ARG TYPE> \
{ \
	static auto get([[maybe_unused]] SINGLE_ARG TYPE const& data) \
	{ \
		__VA_ARGS__; \
	} \
};

//! Define a new Symbol in the symbolic algebra functionality.
/*!
 * The symbol is added so as to be printed by the symbolic algebra
 * functionality. 
 * Defines a new symbol using the condition to instantiate this
 * symbol specialization.
 * 
 * See #DEFINE_SYMBOL_ID.
 * 
 * \param TEMPLATES The templates which are passed to the new
 * symbol type in the specialization. The argument is always surrounded
 * by parentheses.
 * \param TYPE The new type, which is provided with the template
 * arguments that are given in \p TEMPLATES. The argument is always surrounded
 * by parentheses.
 * \param CONDITION The condition which conditionally compiles the
 * specialization, typically based on some feature of the symbol.
 * The argument is always surrounded by parentheses.
 * \param ... The body of the function which returns the symbol identifier.
 */
#define DEFINE_SYMBOL_ID_CONDITION(TEMPLATES, TYPE, CONDITION, ...) \
template<SINGLE_ARG TEMPLATES> \
struct expr::SymbolID<SINGLE_ARG TYPE, typename std::enable_if<SINGLE_ARG CONDITION>::type> \
{ \
	static auto get([[maybe_unused]] SINGLE_ARG TYPE const& data) \
	{ \
		__VA_ARGS__; \
	} \
};

//! \cond

 /* definition helpers for creating a new data type for use in expression templates
 */
#define SYEX_IMPL_DEFINE_BASE_DATA(TEMPLATES, TYPE, CONST_RETURN_TYPE, RETURN_TYPE, RETURN_EVAL, RETURN_DATA) \
template<SINGLE_ARG TEMPLATES> \
struct expr::BaseData<SINGLE_ARG TYPE> \
{ \
	static CONST_RETURN_TYPE get(SINGLE_ARG TYPE const& data, [[maybe_unused]] iter_type n) \
	{ \
	  return SINGLE_ARG RETURN_EVAL; \
	} \
	static CONST_RETURN_TYPE get(SINGLE_ARG TYPE const& data) \
	{ \
	  return SINGLE_ARG RETURN_DATA; \
	} \
	static RETURN_TYPE get(SINGLE_ARG TYPE & data, [[maybe_unused]] iter_type n) \
	{ \
	  return SINGLE_ARG RETURN_EVAL; \
	} \
	static RETURN_TYPE get(SINGLE_ARG TYPE & data) \
	{ \
	  return SINGLE_ARG RETURN_DATA; \
	} \
};
 //! \endcond

//! Applies an implicit cast of the base data access.
/*!
 * The given data is redirected to the base data of its template arguments.
 * 
 * \param TYPENAME_TEMPLATES The templates which are passed to the 
 * specialization of the base data object.
 * \param REDIRECT_TEMPLATE The template which is passed to the
 * redirected base object, meaning that \p TYPE should be implicitly
 * convertible to this type.
 * \param TYPE The type which is redirected to another base object
 * through the cast to \p REDIRECT_TEMPLATE.
 */
#define DEFINE_BASE_DATA_REDIRECT(TYPENAME_TEMPLATES, REDIRECT_TEMPLATE, TYPE) \
SYEX_IMPL_DEFINE_BASE_DATA((SINGLE_ARG TYPENAME_TEMPLATES), (SINGLE_ARG TYPE), decltype(auto), decltype(auto), \
	(expr::BaseData<SINGLE_ARG REDIRECT_TEMPLATE>::get(data, n)), (expr::BaseData<SINGLE_ARG REDIRECT_TEMPLATE>::get(data)))

//! Define a new base data access method for an array.
/*!
 * Use this on a new object which contains data and which should be used
 * by the symbolic algebra functionality. Provide the name (this can also be
 * a member function) of the object that will be returned. This is done by
 * specializing the expr::BaseData class. The underlying data is
 * assumed (and to use this macro, must be) to be an array type object.
 * 
 * \param TEMPLATES Template arguments used by the type specialization.
 * \param TYPE The types which specialize the expr::BaseData class.
 * \param DATA_NAME The variable name or member function which is used to 
 * get the data. This should return a reference to the underlying data
 * of the data object.
 */
#define DEFINE_BASE_DATA(TEMPLATES, TYPE, DATA_NAME) \
SYEX_IMPL_DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), const auto&, auto&, (data.DATA_NAME[n]), (data.DATA_NAME))

//! Define a new base data access method for a point.
/*!
 * See #DEFINE_BASE_DATA. Applies the access method for a data
 * which has only a point as the underlying data.
 */
#define DEFINE_BASE_DATA_POINT(TEMPLATES, TYPE, DATA_NAME) \
((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), &, (data.DATA_NAME), (data.DATA_NAME))



/* definition helpers that combine data type creation with symbol definition to entirely
* define a new struct for use in expressions
*/

//! Defines a new type of object to be used as data in expressions.
/*!
 * Creates expression specific definitions around an object of block- (array-) 
 * type, which will satisfy the prerequisites of using it in the symbolic 
 * algebra. The given object has must have its template parameters specified in 
 * brackets. The object type specification must also be surrounded in brackets. 
 * The last parameter is the name of the objects underlying data member.
 * 
 * See #DEFINE_BASE_DATA and #DEFINE_SYMBOL_ID
 */
#define ADD_EXPR_TYPE(TEMPLATES, TYPE, DATA_NAME) \
DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), DATA_NAME) \
DEFINE_SYMBOL_ID((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), return data.DATA_NAME)

 //! Defines a new type of object to be used as data in expressions.
 /*!
  * Creates expression specific definitions around an object of block- (array-)
  * type, which will satisfy the prerequisites of using it in the symbolic
  * algebra. The given object has must have its template parameters specified in
  * brackets. The object type specification must also be surrounded in brackets.
  * The last parameter is the name of the objects underlying data member.
  */
#define ADD_EXPR_TYPE_POINT(TEMPLATE, TYPE, DATA_NAME) \
DEFINE_BASE_DATA_POINT(TEMPLATE, TYPE, DATA_NAME) \
DEFINE_SYMBOL_ID(TEMPLATE, TYPE, return &data.DATA_NAME)



// **************************************************************************************






/*
 * Allowing combination means that the type passed to the argument will be
 * collected as like terms by the expression logic
 *
 * NOTE: this is for the original types; this means that only the top level
 * (i.e. children) are check against this rule, and not the parents
 * in other words, the child will not satisfy this rule unless explicitly added
 * moreover, child types that are different will not be combined
 */
#define ALLOW_COMBINATION(TEMPLATES, TYPES) \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_can_combine<SINGLE_ARG TYPES> \
{ \
	static const bool value = true; \
};

#define ALLOW_COMBINATION_CONDITION(TEMPLATES, TYPES, CONDITION) ALLOW_COMBINATION((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPES) COMMA typename std::enable_if<(CONDITION)>::type)

//! Restrict adding a type to nonlinear variables when multiplying.
/*!
 * Application of #RESTRICT_MULTIPLICATION_CONDITION without condition.
 */
#define RESTRICT_MULTIPLICATION(TEMPLATES, TYPES) \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_can_multiply<SINGLE_ARG TYPES> \
{ \
	static const bool value = false; \
};

//! Restrict adding a type to nonlinear variables when multiplying.
/*!
 * Restricting multiplication means that the type passed to the argument will
 * not be allowed to be combined with other objects in the same nonlinear variable,
 * instead it will be accumulated into a nonlinear variable of only its own types.
 * Specializes the behaviour of expr::grid_can_multiply.
 *
 * **NOTE**: this is for the base types; this means that all children will be combined
 * if the parent is added to this rule.
 * Adding a child type that has its parent registered as a base will not add
 * the child to this expression logic.
 *
 * \param TEMPLATES The templates used in specializing the types.
 * \param TYPES The types which specialize the behaviour of expr::grid_can_multiply.
 * \param CONDITION The condition under which to apply the multiple restriction.
 */
#define RESTRICT_MULTIPLICATION_CONDITION(TEMPLATES, TYPES, CONDITION) \
RESTRICT_MULTIPLICATION((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPES COMMA typename std::enable_if<(CONDITION)>::type))





//! Add a rule which is executed in a nonlinear variable.
/*!
 * Creates a new rule which is used when the types are combined in
 * a nonlinear variable, OpNLVariable. The `STATEMENT` argument is
 * code which uses the list of types `datas` and computes a new object
 * which is the result of the rule being applied to the list. This 
 * specializes expr::nl_applied_rule.
 * 
 * A condition is provided to the specialization to evaluate
 * whether or not the rule should be applied.
 * 
 * \param CONDITION A condition for using the rule.
 * \param TEMPLATES the template defining the specialization of expr::nl_applied_rule,
 * and the templates that are used by `TYPES`.
 * \param TYPES A list of types in parentheses, which must directly specialize
 * expr::nl_applied_rule.
 * \param STATEMENT The code which computes the result of the rule.
 */
#define ADD_RULE_NL_CONDITION(CONDITION, TEMPLATES, TYPES, STATEMENT) \
template<SINGLE_ARG TEMPLATES> \
struct expr::nl_applied_rule<SINGLE_ARG TYPES> \
{ \
	static const bool enable = CONDITION; \
	template<typename... Gs> \
	auto operator()(std::tuple<Gs...> const& datas) \
	{ \
		SINGLE_ARG STATEMENT; \
	} \
};

//! Add a rule which is executed in a nonlinear variable.
/*!
 * Creates a new rule which is used when the types are combined in
 * a nonlinear variable, OpNLVariable. The `STATEMENT` argument is
 * code which uses the list of types `datas` and computes a new object
 * which is the result of the rule being applied to the list. This
 * specializes expr::nl_applied_rule.
 *
 * \param TYPES A list of types in parentheses, which must directly specialize
 * expr::nl_applied_rule.
 * \param TEMPLATES The templates which specialize the rule.
 * \param STATEMENT The code which computes the result of the rule.
 */
#define ADD_RULE_NL(TYPES, TEMPLATES, STATEMENT) \
ADD_RULE_NL_CONDITION(true, (SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPES), (SINGLE_ARG STATEMENT))




/*
 * in the case of more complicated object inheritance heirarchies, a base type needs to be defined
 * in order for the expressions to detect the correct type in their evaluations, as the type
 * heirarchy will be descended into to determine relationships
 * (for instance, it may be useful to define
 * a child of an object if it has certain identites, that way it is parameterized as a special type to
 * be accepted by the appropriate function)
 *
 * the first argument are the template arguments to the to the defining objects
 * the second argument is the base type
 * the third argument is the parent of the base type
 */

#define DEFINE_BASE_TYPE(TEMPLATE, TYPE, PARENT) \
template<SINGLE_ARG TEMPLATE> \
struct expr::base_data_type<SINGLE_ARG TYPE> \
{ \
	using type = PARENT; \
};


#define DEFINE_BASE_TYPE_CONDITION(TEMPLATE, TYPE, CONDITION, PARENT) \
DEFINE_BASE_TYPE(TEMPLATE, TYPE COMMA typename std::enable_if<(CONDITION)>::type, PARENT)





/* in some specializations of grids, variables should still be able to divide each other
 * in particular, if the specialization has a power associated with it, then the expression logic
 * should be able to manually divide the grids even if their types don't match exactly
 */

#define ALLOW_DIVISION(TEMPLATES, TYPES) \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_divides<SINGLE_ARG TYPES> \
{ \
	static const bool value = true; \
}; 

#define ALLOW_DIVISION_CONDITION(TEMPLATES, TYPES, CONDITION) \
ALLOW_DIVISION(TEMPLATES, TYPES COMMA typename std::enable_if<(CONDITION)>::type)

 //! @}


namespace expr
{

	//! Mechanism to get the ID of a symbol in the expression algebra.
	/*!
	 * Every kind of base symbol that is defined (such as a grid or a variable)
	 * needs some way to be uniquely represented and printed to output. expr::SymbolID
	 * class needs to be specialized and a get function implemented for each new
	 * object that could be used in an expression, in order to return the address
	 * of some part of the unique object. i.e. the address of a data, or simply
	 * the pointer of an array, is the unique element.
	 * 
	 * \ingroup Op
	 */
	template<typename A, typename Enable = void>
	struct SymbolID
	{
		static A* get(A& a)
		{
			return &a;
		}
		static const A* get(A const& a)
		{
			return &a;
		}
	};
}


/* implementations for grid, variable and an arbitrary pointer class
 * additional grids need to implement how they will return the symbolid struct
 */


DEFINE_SYMBOL_ID((typename A), (A*), return data)
DEFINE_SYMBOL_ID((typename T, size_t D), (Grid<T, D>), return data.values)



// **************************************************************************************

namespace expr
{

	//! Helper for accessing data encapsulated by different types.
	/*!
	 * For the different types of objects that be data to an OpLVariable, a
	 * method is implemented to correctly obtain the value (the data point) at a
	 * given index. Additionally, the entire data itself can be obtained.
	 */
	template<typename A>
	struct BaseData
	{
		static const auto& get(A const& data, iter_type)
		{
			return data;
		}
		static const auto& get(A const& data)
		{
			return data;
		}
		static auto& get(A& data, iter_type)
		{
			return data;
		}
		static auto& get(A& data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<typename A>
	struct BaseData<A*>
	{
		static auto& get(A* data, iter_type n)
		{
			return data[n];
		}
		static auto get(A* data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<typename T>
	struct BaseData<Block<T>>
	{
		static auto const& get(Block<T> const& data, iter_type n)
		{
			return data[n];
		}
		static auto const& get(Block<T> const& data)
		{
			return data;
		}
		static auto& get(Block<T>& data, iter_type n)
		{
			return data[n];
		}
		static auto& get(Block<T>& data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<template<typename, size_t> typename G, typename T, size_t D>
	struct BaseData<G<T, D>>
	{
		static auto const& get(G<T, D> const& data, iter_type n)
		{
			return data[n];
		}
		static auto const& get(G<T, D> const& data)
		{
			return data;
		}
		static auto& get(G<T, D>& data, iter_type n)
		{
			return data[n];
		}
		static auto& get(G<T, D>& data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<typename A>
	struct BaseData<A&> : BaseData<A> {};
}

// *****************************************************************************



DEFINE_BASE_DATA_REDIRECT((typename G), (G), (NamedData<G>))
DEFINE_SYMBOL_ID((typename G), (NamedData<G>), return data.name.c_str())


namespace expr
{

	//! Used to identify the representative base type of a data type.
	/*!
	 * A struct based on type traits that prunes, from a given type, the underlying
	 * data object for all the used op-type wrappers. Typically used in identifying
	 * whether something satisfies an identity or getting the original object type
	 * of an expression data.
	 * 
	 * By default, it returns the type of the object that it is given, but for
	 * the data of OpLVariables, including data using types where there is the
	 * underlying or inherited data object which should represent it as the base
	 * type, it is returned instead.
	 *
	 * This type trait is used for identities in multiplication and division
	 * because it is specialized for those types.
	 */
	template<typename A, typename Enable = void>
	struct base_data_type
	{
		using type = A;
	};

	////! Specialization based on expr::base_data_type.
	//template<typename T, typename G>
	//struct base_data_type<OpLVariable<T, G>>
	//{
	//	using type = typename expr::base_data_type<G>::type;
	//};

	////! Specialization based on expr::base_data_type.
	//template<typename E>
	//struct base_data_type<OpExpression<E>>
	//{
	//	using type = typename expr::base_data_type<E>::type;
	//};

	// *************************************************************************


	//! Indicates whether a data may be multiplied. 
	/*!
	 * Always assumed true for itself. Specialized for types which are the
	 * exception to this behaviour. Specialize this for desired data (grid)
	 * types. If this is false for any grids, they will only be able to combine
	 * into each other and when multiplied by other variables, will result in
	 * OpBinaryMul type.
	 */
	template<typename E, typename Enable = void>
	struct grid_can_multiply
	{
		static const bool value = true;
	};

	//! Indicates whether a data may be multiplied. 
	/*!
	 * Always assumed true for itself. Specialized for types that must act
	 * as the exception to this behaviour
	 */
	template<typename E, typename Enable = void>
	struct grid_can_combine
	{
		static const bool value = false;
	};


	//! Provides the number of times that `C` can divide `G`. 
	/*!
	 * In the specialized case, if `C` is X^2 and G is `X^6` then the value 
	 * should be 3. Used in particular for rules in specialized grids, but in 
	 * general this is 0 as the division is performed primarily between data 
	 * which can be linearly combined (the default behaviour is for Variable)
	 *
	 * for other types that won't follow exponent laws, but can combine, 
	 * the predicate is automatically satisfied.
	 */
	template<typename C, typename G, typename Enable = void>
	struct grid_divides
	{
		static const size_t value = (grid_can_combine<C>::value && std::is_same<C, G>::value) ? 1 : 0;
	};


	// *************************************************************************

	//! Values representing the possible identity types.
	enum class IdentityType { ADD, SUB, MUL, DIV };

	//! A class to be specialized in order to implement identities.
	/*!
	 * The identity needs to be enabled for any specializations. The `apply`
	 * function is implemented to return the result of the identity.
	 * 
	 * \tparam type The type of operation which applies this identity.
	 * \tparam Es The types of expressions involved in the identity.
	 */
	template<IdentityType type, typename... Es>
	struct identity_rule
	{
		static bool get_value() { return false; }
		static const bool enable = false;

		template<typename E1, typename E2>
		static auto apply(E1 const& a, E2 const& b) {}
	};


	//! Type trait standardizing the type for an enumerated factor.
	template<size_t Z>
	using variable_data_factor = std::index_sequence<Z>;


	//! Extends factor checking to enumerated types.
	template<size_t Y, typename G>
	struct factor_count_index;

	//! Return the total number of times that `C` is a factor of `G`.
	template<typename C, typename G>
	struct factor_count
	{
		static const size_t value = expr::grid_divides<
			typename expr::base_data_type<C>::type, 
			typename expr::base_data_type<G>::type>::value;
	};

	//! Gives the sum of factor_count applied with each `G`.
	template<typename C, typename... Gs>
	struct factor_count_list;

	//! Number of times an enumerated factor `A` appears in `E1` and `E2`.
	template<typename E1, typename E2, typename A>
	struct factor_list_variables
	{
		using type = typename std::tuple<>;
	};

	//! Number of times that factor `A` appears in `E1` and `E2`.
	template<typename E1, typename E2, typename A>
	struct factor_list_data
	{
		using type = typename std::tuple<>;
	};

	//! Provides all shared factors between `E1` and `E2`.
	/*!
	 * The factors between the two expressions is computed.
	 * 
	 * \tparam E1 The first expression type.
	 * \tparam E2 The second expression type.
	 */
	template<typename E1, typename E2>
	struct factor_list_all;

}

namespace expr
{

	//! Wrapper for multiplication identity checking.
	/*!
	 * An alias for an identity of the given list of types in the
	 * context of multiplication.
	 * 
	 * \tparam The given types which form the identity.
	 */
	template<typename G0, typename... Gs>
	using m_idy = expr::identity_rule<
		expr::IdentityType::MUL,
		typename expr::base_data_type<G0>::type,
		typename expr::base_data_type<Gs>::type...>;

	//! Wrapper for division identity checking.
	/*!
	 * An alias for an identity of the given list of types in the
	 * context of division.
	 *
	 * \tparam The given types which form the identity.
	 */
	template<typename G0, typename... Gs>
	using d_idy = expr::identity_rule<
		expr::IdentityType::DIV,
		typename expr::base_data_type<G0>::type,
		typename expr::base_data_type<Gs>::type...>;

}

// ******************************************************************************************




/* returns the smallest number of variable data that matches the type C
 * this is the same as expr::property::derivative_index except its for counting grids
 */


//! Specialization based on expr::factor_count;
template<typename C, typename G>
struct expr::factor_count<C, G&>
{
	static const size_t value = expr::factor_count<C, G>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename G>
struct expr::factor_count<C, const G>
{
	static const size_t value = expr::factor_count<C, G>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename G>
struct expr::factor_count<C, symphas::ref<G>>
{
	static const size_t value = expr::factor_count<C, G>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename T, typename G>
struct expr::factor_count<C, OpLVariable<T, G>>
{
	static const size_t value = expr::factor_count<C, G>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename T, typename... Gs>
struct expr::factor_count<C, OpNLVariable<T, Gs...>>
{
	static const size_t value = ((expr::factor_count<C, Gs>::value + ...));
};

//! Specialization based on expr::factor_count;
template<typename C, typename E1, typename E2>
struct expr::factor_count<C, OpBinaryAdd<E1, E2>>
{
	static const size_t value = fixed_min<expr::factor_count<C, E1>::value, expr::factor_count<C, E2>::value>;
};

//! Specialization based on expr::factor_count;
template<typename C, typename E1, typename E2>
struct expr::factor_count<C, OpBinarySub<E1, E2>>
{
	static const size_t value = fixed_min<expr::factor_count<C, E1>::value, expr::factor_count<C, E2>::value>;
};

//! Specialization based on expr::factor_count;
template<typename C, typename E1, typename E2>
struct expr::factor_count<C, OpBinaryMul<E1, E2>>
{
	static const size_t value = expr::factor_count<C, E1>::value + expr::factor_count<C, E2>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename E1, typename E2>
struct expr::factor_count<C, OpBinaryDiv<E1, E2>>
{
	static const size_t value = expr::factor_count<C, E1>::value;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, typename G>
struct expr::factor_count<Variable<Y, G>, Variable<Z, G>>
{
	static const size_t value = (Y == Z) ? 1 : 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, typename G>
struct expr::factor_count<expr::variable_data_factor<Y>, Variable<Z, G>>
{
	static const size_t value = expr::factor_count<Variable<Y, G>, Variable<Z, G>>::value;
};




template<typename C, typename... Gs>
struct expr::factor_count_list
{
	static const size_t value = ((expr::factor_count<C, Gs>::value + ...));
};

//! Specialization based on expr::factor_count_list;
template<typename C, typename... Gs>
struct expr::factor_count_list<C, std::tuple<Gs...>>
{
	static const size_t value = expr::factor_count_list<C, Gs...>::value;
};

template<size_t Y, typename G>
struct expr::factor_count_index
{
	static const size_t value = expr::factor_count<expr::variable_data_factor<Y>, G>::value;
};




/*
 * type traits for evaluating an expression and collecting all the datas which can be factored out
 */



 //! Specialization based on expr::factor_list_variables;
template<typename E1, typename E2, size_t I0>
struct expr::factor_list_variables<E1, E2, std::index_sequence<I0>>
{
protected:

	static const size_t N1 = expr::factor_count_index<I0, E1>::value;
	static const size_t N2 = expr::factor_count_index<I0, E2>::value;
	static const size_t N = fixed_min<N1, N2>;

public:
	using type = typename std::conditional_t<(N > 0), std::tuple<std::pair<std::index_sequence<N>, std::index_sequence<I0>>>, std::tuple<>>;
};

//! Specialization based on expr::factor_list_variables;
template<typename E1, typename E2, size_t I0, size_t I1, size_t... Is>
struct expr::factor_list_variables<E1, E2, std::index_sequence<I0, I1, Is...>>
{
protected:

	using types_1 = typename expr::factor_list_variables<E1, E2, std::index_sequence<I0>>::type;
	using types_2 = typename expr::factor_list_variables<E1, E2, std::index_sequence<I1, Is...>>::type;

public:
	using type = typename symphas::lib::combine_types<types_1, types_2>::type;
};


//! Specialization based on expr::factor_list_data;
template<typename E1, typename E2, typename T0>
struct expr::factor_list_data<E1, E2, std::tuple<T0>>
{
protected:

	static const size_t N1 = expr::factor_count<T0, E1>::value;
	static const size_t N2 = expr::factor_count<T0, E2>::value;
	static const size_t N = fixed_min<N1, N2>;

public:
	using type = typename std::conditional_t<(N > 0), std::tuple<std::pair<std::index_sequence<N>, T0>>, std::tuple<>>;
};

//! Specialization based on expr::factor_list_data;
template<typename E1, typename E2, typename T0, typename T1, typename... Ts>
struct expr::factor_list_data<E1, E2, std::tuple<T0, T1, Ts...>>
{
protected:

	using types_1 = typename expr::factor_list_data<E1, E2, std::tuple<T0>>::type;
	using types_2 = typename expr::factor_list_data<E1, E2, std::tuple<T1, Ts...>>::type;

public:
	using type = typename symphas::lib::combine_types<types_1, types_2>::type;
};


template<typename E1, typename E2>
struct expr::factor_list_all
{
protected:
	using variable_ids = decltype(expr::property::vars<E1>::get_ids());
	using grid_types = typename expr::grid_types<E1>::type;

public:

	using type = typename symphas::lib::combine_types<
		typename expr::factor_list_variables<E1, E2, variable_ids>::type, 
		typename expr::factor_list_data<E1, E2, grid_types>::type>::type;

	//! The number of factors between the two expressions.
	static const bool value = std::tuple_size<type>::value;

};









// ******************************************************************************************

namespace expr
{

	//! Introduces order between expression data. 
	/*!
	 * The variables in an OpNLVariable object can be sorted according to a
	 * greater than inequality between types.
	 *
	 * If one type is greater than another, it is placed afterwards in the tuple.
	 * The algorithm is used for new types by through specialization.
	 *
	 * In general, reordering of types means the types have to be commutative,
	 * although order between matching types will remain unchanged.
	 */
	template<typename A, typename B>
	struct nl_greater_than
	{
		static const bool value = false;
	};

	//! The order between equivalent types is always the same.
	template<typename A>
	struct nl_greater_than<A, A>
	{
		static const bool value = false;
	};

	//! Alias for applying order, applying the base types.
	template<typename A, typename B>
	using nl_gt = nl_greater_than<typename expr::base_data_type<A>::type, typename expr::base_data_type<B>::type>;




	//! An identity applied to an OpNLVariable.
	/*!
	 * Specialized in order to apply specific identities, but default behaviour
	 * is simply returning the OpNLVariable formed directly with the given 
	 * list of data.
	 */
	template<typename... As>
	struct nl_applied_rule
	{
		static const bool enable = false;
		template<typename... Gs>
		auto operator()(std::tuple<Gs...> const& datas)
		{
			return OpNLVariable(OpIdentity{}, datas);
		}
	};

}


namespace symphas::internal
{

	template<typename G0, typename G1, typename... Gs>
	auto NL_sort(std::tuple<G0, G1, Gs...> const& datas);
	template<typename G0>
	auto NL_sort(std::tuple<G0> const& datas);


	template<typename G0, typename G1, typename... Gs, typename std::enable_if_t<expr::nl_gt<G0, G1>::value, int> = 0>
	auto NL_sort_terminate(std::tuple<G0, G1, Gs...> const& datas)
	{
		return symphas::internal::NL_sort(datas);
	}

	template<typename G0, typename G1, typename... Gs, typename std::enable_if_t<!expr::nl_gt<G0, G1>::value, int> = 0>
	auto NL_sort_terminate(std::tuple<G0, G1, Gs...> const& datas)
	{
		return datas;
	}

	template<typename... Gs>
	auto NL_rules(std::tuple<Gs...> const& datas)
	{
		return expr::nl_applied_rule<typename expr::base_data_type<Gs>::type...>{}(datas);
	}

	template<typename T, typename... Gs, typename std::enable_if<expr::nl_applied_rule<typename expr::base_data_type<Gs>::type...>::enable, int>::type = 0>
	auto NL_rules(OpNLVariable<T, Gs...> const& e)
	{
		return expr::make_literal(e.value) * symphas::internal::NL_rules(e.datas);
	}


	template<typename T, typename... Gs, typename std::enable_if<!expr::nl_applied_rule<typename expr::base_data_type<Gs>::type...>::enable, int>::type = 0>
	auto NL_rules(OpNLVariable<T, Gs...> const& e)
	{
		return e;
	}


	struct NL_sort_recurse
	{
		template<typename G0, typename G1, typename... Gs, typename std::enable_if_t<expr::nl_gt<G0, G1>::value, int> = 0>
		auto operator()(std::tuple<G0, G1, Gs...> const& datas)
		{
			return symphas::internal::NL_sort_terminate(
				std::tuple_cat(
					std::tuple<G1>(std::get<1>(datas)),
					symphas::internal::NL_sort(std::tuple_cat(std::tuple<G0>(std::get<0>(datas)), symphas::lib::get_tuple_ge<2>(datas)))
				));
		}

		template<typename G0, typename G1, typename... Gs, typename std::enable_if_t<!expr::nl_gt<G0, G1>::value, int> = 0>
		auto operator()(std::tuple<G0, G1, Gs...> const& datas)
		{
			return symphas::internal::NL_sort_terminate(
				std::tuple_cat(
					std::tuple<G0>(std::get<0>(datas)),
					symphas::internal::NL_sort(symphas::lib::get_tuple_ge<1>(datas))
				));
		}

	};

}


template<typename G0>
auto symphas::internal::NL_sort(std::tuple<G0> const& datas)
{
	return datas;
}


template<typename G0, typename G1, typename... Gs>
auto symphas::internal::NL_sort(std::tuple<G0, G1, Gs...> const& datas)
{
	return NL_sort_recurse{}(datas);
}














