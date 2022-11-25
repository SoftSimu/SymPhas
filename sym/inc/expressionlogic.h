
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

#include "expressionlib.h"


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
	static decltype(auto) get([[maybe_unused]] SINGLE_ARG TYPE const& data) \
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
struct expr::SymbolID<SINGLE_ARG TYPE, typename std::enable_if_t<SINGLE_ARG CONDITION>> \
{ \
	static decltype(auto) get([[maybe_unused]] SINGLE_ARG TYPE const& data) \
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
 * \param REDIRECT_TEMPLATE The template which is passed to the static cast
 * to cast the object to its base type. This also implies that \p TYPE should 
 * be implicitly convertible to this type.
 * \param TYPE The type which is redirected to another base object
 * through the cast to \p REDIRECT_TEMPLATE.
 */
#define DEFINE_BASE_DATA_REDIRECT(TYPENAME_TEMPLATES, REDIRECT_TEMPLATE, TYPE) \
SYEX_IMPL_DEFINE_BASE_DATA((SINGLE_ARG TYPENAME_TEMPLATES), (SINGLE_ARG TYPE), decltype(auto), decltype(auto), \
	(expr::BaseData<SINGLE_ARG REDIRECT_TEMPLATE>::get(data, n)), (expr::BaseData<SINGLE_ARG REDIRECT_TEMPLATE>::get(data)))

//! Define a new base data access method for an array.
/*!
 * Use this for a new object which contains data and which should be used
 * by the symbolic algebra functionality. By defining a expr::BaseData, a variable
 * defined with this type can be evaluated and return a value.
 * 
 * When writing a parameter for the \p RETURN_ELEMENT or \p RETURN_DATA arguments,
 * the instance of type \p TYPE is named `data`. As can example, one can write:
```
	DEFINE_BASE_DATA((typename T, size_t D), (Grid<T, D>), data[n], data)
```
 * which would define a new ::BaseData specialization for the Grid class (one already
 * exists though), and get the element of the grid, and the grid itself.
 *
 * \param TEMPLATES Template arguments used by the type specialization. Wrapped
 * in parentheses.
 * \param TYPE The type that the expr::BaseData is specialized for.
 * \param RETURN_ELEMENT A statement which will access an element 
 * of the specialized type instance, and placed after a `return` keyword
 * in order for it to be returned.
 * \param RETURN_DATA A statement which will return a data object, which can be
 * the instance itself, or the instance cast to a different type.
 */
#define DEFINE_BASE_DATA(TEMPLATES, TYPE, RETURN_ELEMENT, RETURN_DATA) \
SYEX_IMPL_DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), decltype(auto), decltype(auto), (RETURN_ELEMENT), (RETURN_DATA))

//! Define a new base data access method for an array type.
/*!
 * See #DEFINE_BASE_DATA. Applies the access method for a data
 * which has only a point as the underlying data.
 *
 * The type must have its elements able to be accessed using `operator[]`.
 */
#define DEFINE_BASE_DATA_ARRAY(TEMPLATES, TYPE) \
DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), (data[n]), (data))

//! Define a new base data access method for a point type.
/*!
 * See #DEFINE_BASE_DATA. Applies the access method for a data
 * which has only a point as the underlying data.
 * 
 * The name of the value of the type instance is given by \p DATA_NAME.
 */
#define DEFINE_BASE_DATA_POINT(TEMPLATES, TYPE, DATA_NAME) \
DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), (data.DATA_NAME), (data))



/* definition helpers that combine data type creation with symbol definition to entirely
* define a new struct for use in expressions
*/

//! Defines a new array-like type to be used as data in expressions.
/*!
 * Creates expression specific definitions around an object of block- (array-) 
 * type, which will satisfy the prerequisites of using it in the symbolic 
 * algebra. The given object has must have its template parameters specified in 
 * brackets. The object type specification must also be surrounded in brackets. 
 * The last parameter is the name of the objects underlying data member.
 * 
 * The type must have its elements able to be accessed using `operator[]`.
 *
 * See #DEFINE_BASE_DATA and #DEFINE_SYMBOL_ID
 */
#define ADD_EXPR_TYPE(TEMPLATES, TYPE, DATA_NAME) \
DEFINE_BASE_DATA((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), data[n], data.data_NAME) \
DEFINE_SYMBOL_ID((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPE), return data.DATA_NAME)

 //! Defines a new type of object to be used as data in expressions.
 /*!
  * Creates expression specific definitions around an object of block- (array-)
  * type, which will satisfy the prerequisites of using it in the symbolic
  * algebra. The given object has must have its template parameters specified in
  * brackets. The object type specification must also be surrounded in brackets.
  * The last parameter is the name of the objects underlying data member.
  * 
  * Additionally, the type must have its elements able to be accessed using `operator[]`. 
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

#define ALLOW_COMBINATION_CONDITION(TEMPLATES, TYPES, CONDITION) \
ALLOW_COMBINATION((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPES) COMMA typename std::enable_if_t<(CONDITION) COMMA int>)

//! Restrict adding a type to nonlinear variables when multiplying.
/*!
 * Application of #RESTRICT_MULTIPLICATION_CONDITION without condition.
 */
#define RESTRICT_COMMUTATIVITY(TEMPLATES, TYPES) \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_can_commute<SINGLE_ARG TYPES> \
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
#define RESTRICT_COMMUTATIVITY_CONDITION(TEMPLATES, TYPES, CONDITION) \
RESTRICT_COMMUTATIVITY((SINGLE_ARG TEMPLATES), (SINGLE_ARG TYPES COMMA typename std::enable_if_t<(CONDITION) COMMA int>))


//! Define the base type that the new data object corresponds to.
/*!
 * In the case of more complicated object inheritance heirarchies, a base type needs to be defined
 * in order for the expressions to detect the correct type in their evaluations, as the type
 * heirarchy will be descended into to determine relationships.
 * 
 * (for instance, it may be useful to define
 * a child of an object if it has certain identites, that way it is parameterized as a special type to
 * be accepted by the appropriate function)
 *
 * \param TEMPLATE Template arguments to the to the defining objects.
 * \param TYPE The base type itself.
 * \param PARENT Parent of the base type that the given TYPE is translated to when using this object.
 */
#define DEFINE_BASE_TYPE(TEMPLATE, TYPE, PARENT) \
template<SINGLE_ARG TEMPLATE> \
struct expr::base_data_type<SINGLE_ARG TYPE> \
{ \
	using type = PARENT; \
};


 //! Define the base type that the new data object corresponds to.
  /*!
  * In the case of more complicated object inheritance heirarchies, a base type needs to be defined
  * in order for the expressions to detect the correct type in their evaluations, as the type
  * heirarchy will be descended into to determine relationships.
  * 
  * (for instance, it may be useful to define
  * a child of an object if it has certain identites, that way it is parameterized as a special type to
  * be accepted by the appropriate function)
  *
  * \param TEMPLATE Template arguments to the to the defining objects.
  * \param TYPE The base type itself.
 * \param PARENT Parent of the base type that the given TYPE is translated to when using this object.
  * \param CONDITION The compile-time condition under which the base type will be valid.
  */
#define DEFINE_BASE_TYPE_CONDITION(TEMPLATE, TYPE, PARENT, CONDITION) \
DEFINE_BASE_TYPE(TEMPLATE, TYPE COMMA typename std::enable_if_t<(CONDITION) COMMA int>, PARENT)





/* in some specializations of grids, variables should still be able to divide each other
 * in particular, if the specialization has a power associated with it, then the expression logic
 * should be able to manually divide the grids even if their types don't match exactly
 */


//! Allows division between objects of the specified type.
/*!
 * Defines a rule to allow division between variables which define this type. That is, the given
 * type \p TYPES
 */
#define DEFINE_FACTOR_COUNT(TEMPLATES, TYPE1, TYPE2, COUNT) \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_divides_count<SINGLE_ARG TYPE1, SINGLE_ARG TYPE2> \
{ \
	static const size_t value = COUNT; \
}; \
template<SINGLE_ARG TEMPLATES> \
struct expr::grid_has_identity<SINGLE_ARG TYPE1, SINGLE_ARG TYPE2> \
{ \
	static const bool value = true; \
};

#define DEFINE_FACTOR_COUNT_CONDITION(TEMPLATES, TYPE1, TYPE2, COUNT, CONDITION) \
DEFINE_FACTOR_COUNT(TEMPLATES, TYPE1, TYPE2 COMMA typename std::enable_if_t<(CONDITION) COMMA int>, COUNT)


//! Add a 'dataless' symbol to the symbolic algebra.
/*! 
 * Adds a dataless symbol to the symbolic algebra library. Dataless means that the
 * symbol doesn't store or define any data, it is simply a symbol for algebra manipulations,
 * and being distinguished from other symbols.
 * 
 * The name is defined by the provided parameter, \p TYPE, and a new class of that name is 
 * added to the namespace expr::symbols. 
 * 
 * \param TYPE The name of the new symbol that is added. 
 */
#define ADD_EXPR_TYPE_SYMBOL(TYPE) \
namespace expr::symbols { \
struct TYPE ## _symbol : Symbol {}; \
using TYPE = OpTerm<OpIdentity, TYPE ## _symbol>; } \
DEFINE_SYMBOL_ID((), (expr::symbols::TYPE ## _symbol), return #TYPE) \
DEFINE_BASE_TYPE((), (expr::symbols:: TYPE), expr::symbols::TYPE ## _symbol) \
ALLOW_COMBINATION((), (expr::symbols::TYPE ## _symbol)) \


 // **************************************************************************************



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

	template<typename A, typename Enable>
	struct SymbolID<const A, Enable> : SymbolID<A, Enable> {};

	template<typename A, typename Enable>
	struct SymbolID<A&, Enable> : SymbolID<A, Enable> {};
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
	 * For the different types of objects that be data to an OpTerm, a
	 * method is implemented to correctly obtain the value (the data point) at a
	 * given index. Additionally, the entire data itself can be obtained.
	 */
	template<typename A>
	struct BaseData
	{
		template<typename T>
		static const auto& get(Block<T> const& data, iter_type n)
		{
			return data[n];
		}

		template<typename T>
		static const auto& get(Block<T> const& data)
		{
			return data;
		}

		template<typename T>
		static auto& get(Block<T>& data, iter_type n)
		{
			return data[n];
		}

		template<typename T>
		static auto& get(Block<T>& data)
		{
			return data;
		}

		template<typename T>
		static const auto& get(T const& data, iter_type n)
		{
			return *(&data + n);
		}

		template<typename T>
		static const auto& get(T const& data)
		{
			return data;
		}

		template<typename T>
		static auto& get(T& data, iter_type n)
		{
			return *(&data + n);
		}

	};

	//! Helper for accessing data encapsulated by different types.
	/*!
	 * For the different types of objects that be data to an OpTerm, a
	 * method is implemented to correctly obtain the value (the data point) at a
	 * given index. Additionally, the entire data itself can be obtained.
	 */
	template<typename A>
	struct BaseData<A const> : BaseData<A> {};

	template<typename A>
	struct BaseData<A const*> : BaseData<A*> {};

	//! Specialization based on expr::BaseData.
	template<typename A>
	struct BaseData<A*>
	{
		static A& get(A* data, iter_type n)
		{
			return data[n];
		}

		static A* get(A* data)
		{
			return data;
		}

		static A const& get(const A* data, iter_type n)
		{
			return data[n];
		}

		static const A* get(const A* data)
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
	template<size_t N, typename T>
	struct BaseData<MultiBlock<N, T>>
	{
		static decltype(auto) get(MultiBlock<N, T> const& data, iter_type n)
		{
			return data[n];
		}
		static auto const& get(MultiBlock<N, T> const& data)
		{
			return data;
		}
		static decltype(auto) get(MultiBlock<N, T>& data, iter_type n)
		{
			return data[n];
		}
		static auto& get(MultiBlock<N, T>& data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<size_t N>
	struct BaseData<Variable<N>>
	{
		static auto get(Variable<N> const& data, iter_type n)
		{
			return OpVoid{}.eval();
		}
		static auto get(Variable<N> const& data)
		{
			return OpVoid{};
		}
		static auto get(Variable<N>& data, iter_type n)
		{
			return OpVoid{}.eval();
		}
		static auto get(Variable<N>& data)
		{
			return OpVoid{};
		}
	};

	//! Specialization based on expr::BaseData.
	template<Axis ax, typename T, size_t D>
	struct BaseData<VectorComponentData<ax, T*, D>>
	{
		static auto const& get(VectorComponentData<ax, T*, D> const& data, iter_type n)
		{
			return data.values[n];
		}
		static auto get(VectorComponentData<ax, T*, D> const& data)
		{
			return data;
		}
		static auto& get(VectorComponentData<ax, T*, D>& data, iter_type n)
		{
			return data.values[n];
		}
		static auto get(VectorComponentData<ax, T*, D>& data)
		{
			return data;
		}
	};

	//! Specialization based on expr::BaseData.
	template<typename T, size_t D>
	struct BaseData<GridData<T, D>>
	{
		static auto const& get(GridData<T, D> const& data, iter_type n)
		{
			return BaseData<T*>::get(data, n);
		}
		static auto get(GridData<T, D> const& data)
		{
			return BaseData<T*>::get(data);
		}
		static auto& get(GridData<T, D>& data, iter_type n)
		{
			return BaseData<T*>::get(data, n);
		}
		static auto get(GridData<T, D>& data)
		{
			return BaseData<T*>::get(data);
		}
	};

	//! Specialization based on expr::BaseData.
	template<typename T, size_t N>
	struct BaseData<GridData<MultiBlock<N, T>, N>>
	{
		static decltype(auto) get(GridData<MultiBlock<N, T>, N> const& data, iter_type n)
		{
			return BaseData<MultiBlock<N, T>>::get(data, n);
		}
		static decltype(auto) get(GridData<MultiBlock<N, T>, N> const& data)
		{
			return BaseData<MultiBlock<N, T>>::get(data);
		}
		static decltype(auto) get(GridData<MultiBlock<N, T>, N>& data, iter_type n)
		{
			return BaseData<MultiBlock<N, T>>::get(data, n);
		}
		static decltype(auto) get(GridData<MultiBlock<N, T>, N>& data)
		{
			return BaseData<MultiBlock<N, T>>::get(data);
		}
	};

	//! Specialization based on expr::BaseData.
	template<Axis ax, typename T, size_t D>
	struct BaseData<VectorComponentData<ax, T, D>>
	{
		static T const& get(VectorComponentData<ax, T, D> const& data, iter_type n)
		{
			return data.values;
		}
		static auto const& get(VectorComponentData<ax, T, D> const& data)
		{
			return data;
		}
		static T& get(VectorComponentData<ax, T, D>& data, iter_type n)
		{
			return data.values;
		}
		static auto& get(VectorComponentData<ax, T, D>& data)
		{
			return data;
		}
	};


	//! Specialization based on expr::BaseData.
	template<Axis ax, typename G>
	struct BaseData<VectorComponent<ax, G>>
	{
		template<size_t N, typename T>
		static VectorComponentData<ax, T*, N> resolve_axis_component(MultiBlock<N, T> const& data)
		{
			return { data.values[symphas::axis_to_index(ax)] };
		}

		template<size_t N, typename T>
		static VectorComponentData<ax, T*, N> resolve_axis_component(MultiBlock<N, T>& data)
		{
			return { data.values[symphas::axis_to_index(ax)] };
		}

		template<typename T, size_t D>
		static VectorComponentData<ax, T, D> resolve_axis_component(any_vector_t<T, D> const& data)
		{
			return { data[symphas::axis_to_index(ax)] };
		}

		template<typename T, size_t D>
		static VectorComponentData<ax, T, D> resolve_axis_component(any_vector_t<T, D>& data)
		{
			return { data[symphas::axis_to_index(ax)] };
		}

		template<typename T>
		static auto const& get_index(T const& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}

		template<typename T>
		static auto& get_index(T& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}

		static auto const& get(VectorComponent<ax, G> const& data, iter_type n)
		{
			return get_index(resolve_axis_component(BaseData<G>::get(*static_cast<G const*>(&data))), n);
		}
		static auto get(VectorComponent<ax, G> const& data)
		{
			return resolve_axis_component(BaseData<G>::get(*static_cast<G const*>(&data)));
		}
		static auto& get(VectorComponent<ax, G>& data, iter_type n)
		{
			return get_index(resolve_axis_component(BaseData<G>::get(*static_cast<G*>(&data))), n);
		}
		static auto get(VectorComponent<ax, G>& data)
		{
			return resolve_axis_component(BaseData<G>::get(*static_cast<G*>(&data)));
		}
	};


	//! Specialization based on expr::BaseData.
	template<typename A>
	struct BaseData<A&> : BaseData<A> {};
	template<typename A>
	struct BaseData<NamedData<A>> : BaseData<A> {};
	template<typename A>
	struct BaseData<OpLiteral<A>> : BaseData<A> {};
}

DEFINE_BASE_DATA((typename T, size_t D), (Grid<T, D>), (T)data[n], data)
DEFINE_BASE_DATA((typename T, size_t D), (BoundaryGrid<T, D>), (T)data[n], data)
//DEFINE_BASE_DATA((template<typename, size_t> typename F, typename T, size_t D), (F<T, D>), data.as_grid()[n], data.as_grid())

// *****************************************************************************



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
	 * the data of OpTerms, including data using types where there is the
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

	template<typename A, typename Enable = void>
	using base_data_t = typename base_data_type<A, Enable>::type;

	// *************************************************************************


	//! Indicates whether a data can commute. 
	/*!
	 * Specialized for types which are the exception to this behaviour. 
	 * 
	 * Specialize this for desired data (grid)
	 * types. This affects how terms are combined during multiplication of
	 * OpTerms.
	 */
	template<typename A, typename B, typename Enable = int>
	struct grid_can_commute
	{
		static const bool value = true;
	};

	template<template<typename, size_t> typename G, typename T, typename S, size_t D>
	struct grid_can_commute<G<any_vector_t<T, D>, D>, G<any_vector_t<S, D>, D>, int>
	{
		static const bool value = false;
	};

	template<size_t N, typename T, typename S>
	struct grid_can_commute<MultiBlock<N, T>, MultiBlock<N, S>, int>
	{
		static const bool value = false;
	};

	template<typename T, typename S, size_t D>
	struct grid_can_commute<any_vector_t<T, D>*, any_vector_t<S, D>*, int>
	{
		static const bool value = false;
	};



	//! Indicates whether a data may be multiplied. 
	/*!
	 * Always assumed true for itself. Specialized for types that must act
	 * as the exception to this behaviour
	 */
	template<typename E, typename Enable = int>
	struct grid_can_combine
	{
		static const bool value = false;
	};

	template<typename G1, typename G2, typename Enable = int>
	struct grid_has_identity
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
	template<typename C, typename G, typename Enable = int>
	struct grid_divides_count
	{
		static const size_t value = (grid_can_combine<C>::value && std::is_same<C, G>::value) ? 1 : 0;
	};


	// *************************************************************************


	//! Extends factor checking to enumerated types.
	/*!
	 * Returns the number of factors of `Variable<Y, G>` (where `G` is any type): it just acts 
	 * as a wrapper so that the enumerated type `Y` can be processed by expr::factor_count.
	 * 
	 * \tparam Y The index of the variable that is searched as the factor.
	 * \tparam G The expression type or grid to compare.
	 */
	template<size_t Y, typename G>
	struct factor_count_index;

	//! Return the total number of times that `C` is a factor of `G`.
	template<typename C, typename G>
	struct factor_count
	{
		static const size_t value = expr::grid_divides_count<
			typename expr::base_data_type<C>::type, 
			typename expr::base_data_type<G>::type>::value;
	};

	//! Gives the sum of factor_count applied with each `G`.
	template<typename C, typename G0, typename... Gs>
	struct factor_count_list;

	//! Number of times an enumerated factor `A` appears in `E1` and `E2`.
	template<typename A, typename E1, typename E2>
	struct factor_list_variables
	{
		using type = typename symphas::lib::types_list<>;
	};

	//! Number of times that factor `A` appears in `E1` and `E2`.
	template<typename A, typename E1, typename E2>
	struct factor_list_data
	{
		using type = typename symphas::lib::types_list<>;
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



	namespace
	{

		template<typename E>
		struct check_is_expression
		{
			template<typename EE>
			static constexpr std::true_type _get_value(OpExpression<EE>)
			{
				return {};
			}

			static constexpr std::false_type _get_value(...)
			{
				return {};
			}


			static const bool value = decltype(_get_value(std::declval<E>()))::value;
		};
	}

	template<typename E>
	constexpr bool is_expression = check_is_expression<E>::value;

	template<typename E>
	constexpr bool is_expression<E&> = is_expression<E>;
	template<typename E>
	constexpr bool is_expression<E&&> = is_expression<E>;

	//! Indicates whether a term can be collected into like terms.
	/*!
	 * Predicate indicating whether a data type can combine with another data type.
	 * Combining means that two data are considered "equivalent" and can be put
	 * together into like terms.
	 *
	 * Grid combination requires a level of indirection in order
	 * to allow variables of all types to be combined. This means that some kind
	 * of subtype needs to be compared, as the encapsulating types can vary.
	 * By default, variables are allowed to be combined into like terms.
	 *
	 * \tparam E The type which is checked if it satisfies the predicate.
	 */
	template<typename E, typename Enable = int>
	constexpr bool is_combinable = grid_can_combine<typename base_data_type<E>::type>::value;

	//! Specialization based on expr::is_combinable.
	template<size_t Z, typename G>
	constexpr bool is_combinable<Variable<Z, G>> = true;


	//! Indicates whether a term is commutable (in multiplication).
	/*!
	 * It may be desirable to ensure that some terms are not commutable if there
	 * are special rules of algebra specific to those terms.
	 *
	 * \tparam E The type which is checked if it satisfies the predicate.
	 */
	template<typename G1, typename G2>
	constexpr bool is_commutable = expr::grid_can_commute<
		typename expr::base_data_type<G1>::type, 
		typename expr::base_data_type<G2>::type, int>::value;

	//! Indicates whether or not there is a multiplication rule between two grids.
	template<typename G1, typename G2>
	constexpr bool has_identity = grid_has_identity<
		typename base_data_type<G1>::type, 
		typename base_data_type<G2>::type, int>::value;



	template<typename T>
	auto eval(T const& value)
	{
		if constexpr (expr::is_expression<T>)
		{
			if constexpr (expr::is_identity<T> || expr::is_fraction<T>)
			{
				return value.eval();
			}
			else
			{
				return eval(value.eval());
			}
		}
		else
		{
			return value;
		}
	}

}

// ******************************************************************************************
// Specializations of SymbolID and BaseData.
// ******************************************************************************************

namespace expr
{
	//! Specialization based on SymbolID.
	template<typename A>
	struct SymbolID<symphas::ref<A>>
	{
		static decltype(auto) get(symphas::ref<A> const& data)
		{
			return SymbolID<A>::get(data);
		}
	};

	//! Specialization based on SymbolID.
	template<size_t Z, typename G>
	struct SymbolID<Variable<Z, G>>
	{
		static decltype(auto) get(Variable<Z, G> const& data)
		{
			return SymbolID<G>::get(data);
		}
	};


	//! Specialization based on expr::BaseData.
	template<typename T>
	struct BaseData<symphas::ref<T>>
	{
		static decltype(auto) get(symphas::ref<T> data, iter_type n)
		{
			return BaseData<T>::get(data.get(), n);
		}

		static decltype(auto) get(symphas::ref<T> data)
		{
			return BaseData<T>::get(data.get());
		}
	};

	//! Specialization based on BaseData.
	template<size_t Z, typename T>
	struct BaseData<Variable<Z, T>>
	{
		static decltype(auto) get(Variable<Z, T> const& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}
		static decltype(auto) get(Variable<Z, T> const& data)
		{
			return BaseData<T>::get(data);
		}
		static decltype(auto) get(Variable<Z, T>& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}
		static decltype(auto) get(Variable<Z, T>& data)
		{
			return BaseData<T>::get(data);
		}
	};
}


//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<symphas::ref<A>>
{
	using type = typename expr::base_data_type<A>::type;
};

//! Specialization based on expr::base_data_type.
template<size_t Z, typename T>
struct expr::base_data_type<Variable<Z, T>>
{
	using type = typename expr::base_data_type<T>::type;
};

//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<const A>
{
	using type = typename expr::base_data_type<A>::type;
};

//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<A&>
{
	using type = typename expr::base_data_type<A>::type;
};





// ******************************************************************************************




/* returns the smallest number of variable data that matches the type C
 * this is the same as expr::derivative_index except its for counting grids
 */


 //! Specialization based on expr::factor_count;
template<typename C, typename G>
struct expr::factor_count<OpTerm<OpIdentity, C>, G>
{
	static const size_t value = expr::factor_count<C, G>::value;
};

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
template<typename C, typename G, expr::exp_key_t X>
struct expr::factor_count<C, Term<G, X>>
{
	static const size_t value = 
		((expr::_Xk_t<X>::N - expr::_Xk_t<X>::N % expr::_Xk_t<X>::D) / expr::_Xk_t<X>::D) 
		* expr::factor_count<C, G>::value;
};

//! Specialization based on expr::factor_count;
template<typename C, typename V, typename... Gs, expr::exp_key_t... Xs>
struct expr::factor_count<C, OpTerms<V, Term<Gs, Xs>...>>
{
	static const size_t value = ((expr::factor_count<C, Term<Gs, Xs>>::value + ...));
};

//! Specialization based on expr::factor_count;
template<typename C, typename... Es>
struct expr::factor_count<C, OpAdd<Es...>>
{
	static const size_t value = fixed_min<expr::factor_count<C, Es>::value...>;
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
template<size_t Y, size_t Z>
struct expr::factor_count<Variable<Y>, Variable<Z>>
{
	static const size_t value = (Y == Z) ? 1 : 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, typename G>
struct expr::factor_count<Variable<Y, G>, Variable<Z, G>>
{
	static const size_t value = (Y == Z) ? 1 : 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, typename G>
struct expr::factor_count<Variable<Y>, Variable<Z, G>>
{
	static const size_t value = expr::factor_count<Variable<Y>, Variable<Z>>::value;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, Axis ax0, Axis ax1, typename G, typename F>
struct expr::factor_count<Variable<Y, VectorComponent<ax0, G>>, Variable<Z, VectorComponent<ax1, F>>>
{
	static const size_t value = (ax0 == ax1) ? expr::factor_count<Variable<Y, G>, Variable<Z, F>>::value : 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, Axis ax0, Axis ax1, typename F>
struct expr::factor_count<Variable<Y, VectorComponent<ax0>>, Variable<Z, VectorComponent<ax1, F>>>
{
	static const size_t value = (ax0 == ax1) ? expr::factor_count<Variable<Y>, Variable<Z, F>>::value : 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, Axis ax1, typename G, typename F>
struct expr::factor_count<Variable<Y, G>, Variable<Z, VectorComponent<ax1, F>>>
{
	static const size_t value = 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, Axis ax1, typename F>
struct expr::factor_count<Variable<Y>, Variable<Z, VectorComponent<ax1, F>>>
{
	static const size_t value = 0;
};

//! Specialization based on expr::factor_count;
template<size_t Y, size_t Z, Axis ax0, typename G, typename F>
struct expr::factor_count<Variable<Y, VectorComponent<ax0, G>>, Variable<Z, F>>
{
	static const size_t value = 0;
};



template<typename C, typename G0, typename... Gs>
struct expr::factor_count_list
{
	static const size_t value = (expr::factor_count<C, G0>::value + ... + expr::factor_count<C, Gs>::value);
};

//! Specialization based on expr::factor_count_list;
template<typename C, typename G0, typename... Gs>
struct expr::factor_count_list<C, symphas::lib::types_list<G0, Gs...>>
{
	static const size_t value = expr::factor_count_list<C, G0, Gs...>::value;
};

//! Specialization based on expr::factor_count_list;
template<typename C>
struct expr::factor_count_list<C, symphas::lib::types_list<>>
{
	static const size_t value = 0;
};

template<size_t Y, typename G>
struct expr::factor_count_index
{
	static const size_t value = expr::factor_count<Variable<Y>, G>::value;
};




/*
 * type traits for evaluating an expression and collecting all the datas which can be factored out
 */



 //! Specialization based on expr::factor_list_variables;
template<size_t I0, typename E1, typename E2>
struct expr::factor_list_variables<std::index_sequence<I0>, E1, E2>
{
protected:

	static const size_t N1 = expr::factor_count_index<I0, E1>::value;
	static const size_t N2 = expr::factor_count_index<I0, E2>::value;
	static const size_t N = fixed_min<N1, N2>;

public:
	using type = typename std::conditional_t<
		(N > 0), 
		symphas::lib::types_list<
			std::pair<std::index_sequence<N>, 
			std::index_sequence<I0>>>, 
		symphas::lib::types_list<>>;
};

//! Specialization based on expr::factor_list_variables;
template<size_t I0, size_t I1, size_t... Is, typename E1, typename E2>
struct expr::factor_list_variables<std::index_sequence<I0, I1, Is...>, E1, E2>
{
protected:

	using types_1 = typename expr::factor_list_variables<std::index_sequence<I0>, E1, E2>::type;
	using types_2 = typename expr::factor_list_variables<std::index_sequence<I1, Is...>, E1, E2>::type;

public:
	using type = symphas::lib::expand_types_list<types_1, types_2>;
};


//! Specialization based on expr::factor_list_data;
template<typename E1, typename E2, typename T0>
struct expr::factor_list_data<symphas::lib::types_list<T0>, E1, E2>
{
protected:

	static const size_t N1 = expr::factor_count<T0, E1>::value;
	static const size_t N2 = expr::factor_count<T0, E2>::value;
	static const size_t N = fixed_min<N1, N2>;

public:
	using type = typename std::conditional_t<
		(N > 0), 
		symphas::lib::types_list<std::pair<std::index_sequence<N>, T0>>, 
		symphas::lib::types_list<>>;
};

//! Specialization based on expr::factor_list_data;
template<typename E1, typename E2, typename T0, typename T1, typename... Ts>
struct expr::factor_list_data<symphas::lib::types_list<T0, T1, Ts...>, E1, E2>
{
protected:

	using types_1 = typename expr::factor_list_data<symphas::lib::types_list<T0>, E1, E2>::type;
	using types_2 = typename expr::factor_list_data<symphas::lib::types_list<T1, Ts...>, E1, E2>::type;

public:
	using type = symphas::lib::expand_types_list<types_1, types_2>;
};


template<typename E1, typename E2>
struct expr::factor_list_all
{
protected:
	using variable_ids = decltype(expr::vars<E1>::get_ids());
	using term_types = typename expr::term_types<E1>::type;

public:

	using type = symphas::lib::expand_types_list<
		typename expr::factor_list_variables<variable_ids, E1, E2>::type,
		typename expr::factor_list_data<term_types, E1, E2>::type>;

	//! The number of factors between the two expressions.
	static const bool value = symphas::lib::types_list_size<type>::value;

};


namespace symphas::internal
{
	//! Helper in order to construct an OpTerms.
	/*!
	 * Allows the type to be delegated to the correct function, so that
	 * the correct OpTerms object is generated.
	 * In particular, it is important that if a reference is made, that
	 * an additional reference is not nested.
	 *
	 * The usual template will assume that an rvalue is given, and
	 * therefore the data can be moved into the OpTerms.
	 *
	 * \tparam A The type being given to create an OpTerms.
	 */
	template<typename A>
	struct construct_term;

}















