
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
 * PURPOSE: Defines the basic functionality of using the symbolic algebra
 * framework.
 *
 * ***************************************************************************
 */

#pragma once

#include <iterator>


#include "grid.h"
#include "gridfunctions.h"
#include "expressionprototypes.h"
#include "indexseqhelpers.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

#define DERIV_MAX_ORDER 24				//!< Maximum order a derivative can be.



namespace expr
{

	//! Type trait used to identify the original type of a data.
	/*!
	 * Applied to the data object of a variable expression in order to obtain
	 * the original data type after certain wrapping types are taken away.
	 * Applied to any other data type, it will simply directly return it.
	 */
	template<typename A>
	struct original_data_type
	{
		using type = A;
	};

	//! Specialization for a reference. See expr::original_data_type.
	template<typename A>
	struct original_data_type<symphas::ref<A>>
	{
		using type = typename original_data_type<A>::type;
	};

	//! Specialization for NamedData type. See expr::original_data_type.
	template<typename G>
	struct original_data_type<NamedData<G>>
	{
		using type = typename original_data_type<G>::type;
	};

	//! Specialization for Variable type. See expr::original_data_type.
	template<size_t Z, typename T>
	struct original_data_type<Variable<Z, T>>
	{
		using type = typename original_data_type<T>::type;
	};

	//! Specialization for a constant type. See expr::original_data_type.
	template<typename A>
	struct original_data_type<const A>
	{
		using type = typename original_data_type<A>::type;
	};

	//! Specialization for a reference type. See expr::original_data_type.
	template<typename A>
	struct original_data_type<A&>
	{
		using type = typename original_data_type<A>::type;
	};


}



// *******************************************************************************

namespace expr
{

	//! Determines the dimension of the data in the expression.
	/*!
	 * Type trait that will prune the grid dimension from an expression based
	 * on its type. If there is no dimension, the evaluated dimension is 0.
	 */
	template<typename E>
	struct grid_dim
	{
		static const int dimension = 0;
		static const int value = dimension;
	};

	//! Constructs the list of all data types.
	/*!
	 * Constructs a tuple type where the tuple template arguments are (uniquely
	 * listed) template arguments of the OpLVariables from the given
	 * expression.
	 */
	template<typename E>
	struct op_types
	{
		using type = std::tuple<>;
	};

	//! Determines the underlying data in the expression.
	/*!
	 * Type trait returns as its type a tuple of data that is found in the
	 * expression.
	 */
	template<typename E>
	struct grid_types;

	//! Determines the type that the given expression evaluates to.
	/*!
	 * Determines the type that the given expression evaluates to.
	 */
	template<typename E>
	struct eval_type;


	//! Determines the grid which would enclose the data of the expression.
	/*!
	 * Determines the evaluation type of the expression, and then builds the
	 * type based on the type of the grid that first appears in the given 
	 * expression type. This is distinct from expr::eval_type, in that
	 * the resulting type here is the full grid.
	 */
	template<typename E>
	struct grid_type;

	//! Checks whether the expression manages data to be updated.
	/*!
	 * Typically for optimization purposes (and correctness purposes also),
	 * parse the expression hierarchy and look for one of the expressions 
	 * having a type that stores a grid that must be updated. For instance,
	 * the prune function must be run on expressions with a state, but this
	 * check isn't strictly required.
	 */
	template<typename E>
	struct has_state
	{
		static const bool value = false;
	};


	//! Determines if an expression contains only constants and derivatives.
	/*!
	 * Parses the expression to determine if it composed of only constants
	 * and derivative operators. This would constitute an expression which is 
	 * only an operator. I.e. this predicate checks if the expression is an
	 * operator type expression.
	 */
	template<typename E>
	struct is_operator_like
	{
		static const bool value = false;
	};

	//! Determines whether the expression is linear in complex space.
	/*!
	 * A type trait that identifies the expressions as either
	 * linear or not. By default, expressions are nonlinear.
	 * 
	 * Expressions are considered linear if there are no OpNLVariables. If an 
	 * expression is a function of more than one variable, they must all appear 
	 * linearly.
	 * 
	 * Linear expressions to which an operator satisfying is_operator_like
	 * is applied are considered linear.
	 */
	template<typename E>
	struct is_linear
	{
		static const bool value = false;
	};


	//! Determines whether the expression is entirely nonlinear.
	/*!
	 * A type trait that identifies the expressions as either
	 * nonlinear or not. This checks that each individual term does not satisfy
	 * linear. If so, then it is nonlinear, but if any terms are linear, then
	 * this predicate is not satisfied.
	 */
	template<typename E>
	struct is_nonlinear
	{
		static const bool value = !is_linear<E>::value;
	};

	//! Determines if the expression is a linear combination.
	/*
	 * As opposed to expr::is_linear, this has the added value of returning 
	 * whether the expression is a linear combination type; only add and 
	 * subtract expressions satisfy this.
	 */
	template<typename E>
	struct is_combination
	{
		static const bool value = false;
	};

	template<typename E>
	struct is_directional_derivative
	{

	private:
		template <typename T, typename = int>
		struct HasAx : std::false_type {};

		template <typename T>
		struct HasAx<T, decltype((void)T::ax, 0)> : std::true_type { };

	public:

		static const bool value = HasAx<E>::value;
	};


	//! Combines type traits is_linear and is_combination.
	template<typename E>
	struct expression_predicate;


	//! Reevaluate an entire expression.
	/*!
	 * Reevaluation of an expression, used for instance, after it had 
	 * identities substituted.
	 */
	template<typename E>
	auto reevaluate(OpExpression<E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinaryAdd<E1, E2> const& e)
	{
		return e.a + e.b;
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinarySub<E1, E2> const& e)
	{
		return e.a - e.b;
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinaryMul<E1, E2> const& e)
	{
		return e.a * e.b;
	}


	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinaryDiv<E1, E2> const& e)
	{
		return e.a / e.b;
	}



}


namespace expr::property
{


	//! Get the index corresponding to the order of a derivative.
	/*!
	 * Return the index for a derivative with the given order value. The index
	 * is generated by signing the OD-th order bit; consequently linear
	 * combinations of derivatives are formed by each term signing the
	 * corresponding bit, where combinations are derivatives applied to
	 * other derivatives.
	 */
	template<size_t OD>
	struct derivative_index_raw
	{
		static const size_t value = 1ull << OD;
	};

	//! Determine the derivative index order in the expression.
	/*!
	 * Returns the smallest index (associated with the derivative) for the
	 * derivative with index greater than the given index `L`. If there are no
	 * greater derivatives, it will return the largest derivative. If there
	 * are terms with no derivative, the index is 1.
	 *
	 * \tparam L The index which the derivative order must be greater.
	 * \tparam E The expression type to parse.
	 */
	template<size_t L, typename E>
	struct derivative_index
	{
		static const size_t value = 1;
	};


}


//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<const E>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<E&>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<OpExpression<E>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<template<typename, size_t> typename G, typename T, size_t D>
struct expr::grid_dim<G<T, D>>
{
	static const int dimension = D;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A>
struct expr::grid_dim<symphas::ref<A>>
{
	static const int dimension = expr::grid_dim<A>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G>
struct expr::grid_dim<NamedData<G>>
{
	static const int dimension = expr::grid_dim<G>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<size_t Z, typename G>
struct expr::grid_dim<Variable<Z, G>>
{
	static const int dimension = expr::grid_dim<G>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename T, typename G>
struct expr::grid_dim<OpLVariable<T, G>>
{
	static const int dimension = expr::grid_dim<G>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename T, typename G0, typename G1>
struct expr::grid_dim<OpNLVariable<T, G0, G1>>
{
	static const int dimension = fixed_max<
		expr::grid_dim<G0>::dimension, 
		expr::grid_dim<G1>::dimension>;

	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename T, typename G0, typename G1, typename G2, typename... Gs>
struct expr::grid_dim<OpNLVariable<T, G0, G1, G2, Gs...>>
{
protected:
	static const int dimensionN = expr::grid_dim<OpNLVariable<T, G1, G2, Gs...>>::dimension;

public:
	static const int dimension = fixed_max<expr::grid_dim<G0>::dimension, dimensionN>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinaryAdd<A, B>>
{
protected:
	static const int Ad = expr::grid_dim<A>::dimension;
	static const int Bd = expr::grid_dim<B>::dimension;

public:
	static const int dimension = fixed_max<Ad, Bd>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinaryMul<A, B>>
{
protected:
	static const int Ad = expr::grid_dim<A>::dimension;
	static const int Bd = expr::grid_dim<B>::dimension;

public:
	static const int dimension = fixed_max<Ad, Bd>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinaryDiv<A, B>>
{
protected:
	static const int Ad = expr::grid_dim<A>::dimension;
	static const int Bd = expr::grid_dim<B>::dimension;

public:
	static const int dimension = fixed_max<Ad, Bd>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinarySub<A, B>>
{
protected:
	static const int Ad = expr::grid_dim<A>::dimension;
	static const int Bd = expr::grid_dim<B>::dimension;

public:
	static const int dimension = fixed_max<Ad, Bd>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E1, typename E2>
struct expr::grid_dim<OpFuncConvolution<V, E1, E2>>
{
protected:
	static const int Ad = expr::grid_dim<E1>::dimension;
	static const int Bd = expr::grid_dim<E2>::dimension;

public:
	static const int dimension = fixed_max<Ad, Bd>;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<size_t D>
struct expr::grid_dim<GaussianSmoothing<D>>
{
	static const int dimension = D;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G, typename V, typename E>
struct expr::grid_dim<OpMap<G, V, E>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::grid_dim<OpFuncDerivative<Dd, V, E, Sp>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A1, typename A2, typename E>
struct expr::grid_dim<OpCombination<A1, A2, E>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A1, typename A2, typename E>
struct expr::grid_dim<OpChain<A1, A2, E>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E>
struct expr::grid_dim<OpExponential<V, E>>
{
	static const int dimension = expr::grid_dim<E>::dimension;
	static const int value = dimension;
};


// *******************************************************************************


//! Specialization based on expr::op_types.
template<typename E>
struct expr::op_types<const E>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename E>
struct expr::op_types<E&>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename E>
struct expr::op_types<OpExpression<E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename T, typename G>
struct expr::op_types<OpLVariable<T, G>>
{
	using type = std::tuple<G>;
};

//! Specialization based on expr::op_types.
template<typename T, typename G0>
struct expr::op_types<OpNLVariable<T, G0>>
{
	using type = typename expr::op_types<OpLVariable<T, G0>>::type;
};

//! Specialization based on expr::op_types.
template<typename T, typename G0, typename G1, typename... Gs>
struct expr::op_types<OpNLVariable<T, G0, G1, Gs...>>
{
	using type = typename symphas::lib::combine_types_unique<
		typename expr::op_types<OpLVariable<T, G0>>::type, 
		typename expr::op_types<OpNLVariable<T, G1, Gs...>>::type>::type;
};

//! Specialization based on expr::op_types.
template<typename A, typename B>
struct expr::op_types<OpBinaryAdd<A, B>>
{
protected:
	using At = typename expr::op_types<A>::type;
	using Bt = typename expr::op_types<B>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A, typename B>
struct expr::op_types<OpBinaryMul<A, B>>
{
protected:
	using At = typename expr::op_types<A>::type;
	using Bt = typename expr::op_types<B>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A, typename B>
struct expr::op_types<OpBinaryDiv<A, B>>
{
protected:
	using At = typename expr::op_types<A>::type;
	using Bt = typename expr::op_types<B>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A, typename B>
struct expr::op_types<OpBinarySub<A, B>>
{
protected:
	using At = typename expr::op_types<A>::type;
	using Bt = typename expr::op_types<B>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename V, typename E1, typename E2>
struct expr::op_types<OpFuncConvolution<V, E1, E2>>
{
protected:
	using At = typename expr::op_types<E1>::type;
	using Bt = typename expr::op_types<E2>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<size_t D>
struct expr::op_types<GaussianSmoothing<D>>
{
	using type = typename expr::op_types<OpLVariable<OpIdentity, decltype(GaussianSmoothing<D>::data)>>::type;
};

//! Specialization based on expr::op_types.
template<typename G, typename V, typename E>
struct expr::op_types<OpMap<G, V, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::op_types<OpFuncDerivative<Dd, V, E, Sp>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2, typename E>
struct expr::op_types<OpCombination<A1, A2, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2, typename E>
struct expr::op_types<OpChain<A1, A2, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename V, typename E>
struct expr::op_types<OpExponential<V, E>>
{
	using type = typename expr::op_types<E>::type;
};

// *******************************************************************************


template<typename... Gs>
struct expr::grid_types<std::tuple<Gs...>>
{
	using type = std::tuple<typename original_data_type<Gs>::type...>;
};

template<typename E>
struct expr::grid_types
{
	using type = typename grid_types<typename expr::op_types<E>::type>::type;
};




// *******************************************************************************


template<typename E>
struct expr::eval_type
{

protected:

	static auto _test_eval(E e, iter_type n)
	{
		return e.eval(n);
	}

public:

	using type = typename std::invoke_result_t<decltype(&expr::eval_type<E>::_test_eval), E, iter_type>;
};

// *******************************************************************************


namespace symphas::lib
{

	//! Get the first type of types that are in a tuple.
	/*!
	 * If the given type is not a tuple of types, then type trait always maps
	 * to void.
	 */
	template<typename G>
	struct get_first_type
	{
		using type = void;
	};

	//! Specialization based on get_first_type.
	template<typename G0, typename... Gs>
	struct get_first_type<std::tuple<G0, Gs...>>
	{
		using type = G0;
	};
}


//! Return the type of the grid in the expression as a complete grid.
template<typename E>
struct expr::grid_type
{
protected:

	using grid_t = typename symphas::lib::get_first_type<typename expr::grid_types<E>::type>::type;

	template<template<typename, size_t> typename enc_type>
	struct grid_class_wrap 
	{
		using type = enc_type<
			typename expr::eval_type<E>::type,
			expr::grid_dim<E>::value>;
	};

	template<typename und_type>
	struct block_class_wrap
	{
		using type = Block<und_type>;
	};

	template<template<typename, size_t> typename G, typename T, size_t D>
	static grid_class_wrap<G> _pack_grid(G<T, D>)
	{
		return {};
	}

	template<typename T>
	static block_class_wrap<T> _pack_grid(T*)
	{
		return {};
	}

	static int _pack_grid()
	{
		return 0;
	}


	static auto pack_grid(grid_t g)
	{
		return _pack_grid(g);
	}


public:
	using type = typename std::invoke_result_t<decltype(&expr::grid_type<E>::pack_grid), grid_t>::type;
};



//! Specialization based on expr::has_state.
template<typename V, typename E1, typename E2>
struct expr::has_state<OpFuncConvolution<V, E1, E2>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename Dd, typename V, typename G, typename Sp>
struct expr::has_state<OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>>
{
	static const bool value = false;
};

//! Specialization based on expr::has_state.
template<typename V, typename S, typename G>
struct expr::has_state<OpExponential<V, OpLVariable<S, G>>>
{
	static const bool value = false;
};

//! Specialization based on expr::has_state.
template<typename V, typename E>
struct expr::has_state<OpExponential<V, E>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename G, typename V, typename E>
struct expr::has_state<OpMap<G, V, E>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::has_state<OpFuncDerivative<Dd, V, E, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename A1, typename A2, typename E>
struct expr::has_state<OpChain<A1, A2, E>>
{
	static const bool value = expr::has_state<E>::value;
};

//! Specialization based on expr::has_state.
template<typename A1, typename A2, typename E>
struct expr::has_state<OpCombination<A1, A2, E>>
{
	static const bool value = expr::has_state<E>::value;
};

//! Specialization based on expr::has_state.
template<typename E1, typename E2>
struct expr::has_state<OpBinaryAdd<E1, E2>>
{
	static const bool value = expr::has_state<E1>::value || expr::has_state<E2>::value;
};

//! Specialization based on expr::has_state.
template<typename E1, typename E2>
struct expr::has_state<OpBinaryMul<E1, E2>>
{
	static const bool value = expr::has_state<E1>::value || expr::has_state<E2>::value;
};

//! Specialization based on expr::has_state.
template<typename E1, typename E2>
struct expr::has_state<OpBinaryDiv<E1, E2>>
{
	static const bool value = expr::has_state<E1>::value || expr::has_state<E2>::value;
};

//! Specialization based on expr::has_state.
template<typename E1, typename E2>
struct expr::has_state<OpBinarySub<E1, E2>>
{
	static const bool value = expr::has_state<E1>::value || expr::has_state<E2>::value;
};

//! Specialization based on expr::has_state.
template<typename... Rs>
struct expr::has_state<std::tuple<Rs...>>
{
	static const bool value = ((expr::has_state<Rs>::value || ...));
};

//! Specialization based on expr::has_state.
template<typename G, typename E>
struct expr::has_state<std::pair<G&, E>>
{
	static const bool value = expr::has_state<E>::value;
};


// *******************************************************************************

template<typename E>
struct expr::is_operator_like<OpExpression<E>>
{
	static const bool value = expr::is_operator_like<E>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename T>
struct expr::is_operator_like<OpLiteral<T>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<>
struct expr::is_operator_like<OpIdentity>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<>
struct expr::is_operator_like<OpNegIdentity>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<size_t D>
struct expr::is_operator_like<GaussianSmoothing<D>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<size_t O, typename V, typename Sp>
struct expr::is_operator_like<OpOperatorDerivative<O, V, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_operator_like<OpFuncDerivative<Dd, V, E, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<typename A1, typename A2>
struct expr::is_operator_like<OpOperatorCombination<A1, A2>>
{
	static const bool value = expr::is_operator_like<A1>::value && expr::is_operator_like<A2>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename A1, typename A2, typename E>
struct expr::is_operator_like<OpCombination<A1, A2, E>>
{
	static const bool value = expr::is_operator_like<OpOperatorCombination<A1, A2>>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename A1, typename A2>
struct expr::is_operator_like<OpOperatorChain<A1, A2>>
{
	static const bool value = expr::is_operator_like<A1>::value && expr::is_operator_like<A2>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename A1, typename A2, typename E>
struct expr::is_operator_like<OpChain<A1, A2, E>>
{
	static const bool value = expr::is_operator_like<OpOperatorChain<A1, A2>>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename E1, typename E2>
struct expr::is_operator_like<OpBinaryAdd<E1, E2>>
{
	static const bool value = expr::is_operator_like<E1>::value && expr::is_operator_like<E2>::value;
};

//! Specialization based on expr::is_operator_like.
template<typename E1, typename E2>
struct expr::is_operator_like<OpBinarySub<E1, E2>>
{
	static const bool value = expr::is_operator_like<E1>::value && expr::is_operator_like<E2>::value;
};

// *******************************************************************************


template<typename E>
struct expr::is_linear<OpExpression<E>>
{
	static const bool value = expr::is_linear<E>::value;
};

//! Specialization based on expr::is_linear.
template<typename T>
struct expr::is_linear<OpLiteral<T>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<>
struct expr::is_linear<OpIdentity>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<>
struct expr::is_linear<OpNegIdentity>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<size_t D>
struct expr::is_linear<GaussianSmoothing<D>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename Dd, typename V, typename G, typename Sp>
struct expr::is_linear<OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_linear<OpFuncDerivative<Dd, V, E, Sp>>
{
	static const bool value = expr::is_linear<E>::value;
};

//! Specialization based on expr::is_linear.
template<size_t O, typename V, typename Sp>
struct expr::is_linear<OpOperatorDerivative<O, V, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename V, typename E1, typename E2>
struct expr::is_linear<OpFuncConvolution<V, E1, E2>>
{
	static const bool value = expr::is_linear<E1>::value && expr::is_linear<E2>::value;
};

//! Specialization based on expr::is_linear.
template<typename T, typename G>
struct expr::is_linear<OpLVariable<T, G>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename A1, typename A2>
struct expr::is_linear<OpOperatorCombination<A1, A2>>
{
	static const bool value = expr::is_operator_like<OpOperatorCombination<A1, A2>>::value;
};

//! Specialization based on expr::is_linear.
template<typename A1, typename A2, typename E>
struct expr::is_linear<OpCombination<A1, A2, E>>
{
	static const bool value = expr::is_linear<OpOperatorCombination<A1, A2>>::value && expr::is_linear<E>::value;
};

//! Specialization based on expr::is_linear.
template<typename A1, typename A2>
struct expr::is_linear<OpOperatorChain<A1, A2>>
{
	static const bool value = expr::is_operator_like<OpOperatorChain<A1, A2>>::value;
};

//! Specialization based on expr::is_linear.
template<typename A1, typename A2, typename E>
struct expr::is_linear<OpChain<A1, A2, E>>
{
	static const bool value = expr::is_linear<OpOperatorChain<A1, A2>>::value && expr::is_linear<E>::value;
};

//! Specialization based on expr::is_linear.
template<typename E1, typename E2>
struct expr::is_linear<OpBinaryAdd<E1, E2>>
{
	static const bool value = expr::is_linear<E1>::value && expr::is_linear<E2>::value;
};

//! Specialization based on expr::is_linear.
template<typename E1, typename E2>
struct expr::is_linear<OpBinarySub<E1, E2>>
{
	static const bool value = expr::is_linear<E1>::value && expr::is_linear<E2>::value;
};


// *******************************************************************************


//! Specialization based on expr::is_nonlinear.
template<typename E>
struct expr::is_nonlinear<OpExpression<E>>
{
	static const bool value = expr::is_nonlinear<E>::value;
};

//! Specialization based on expr::is_nonlinear.
template<typename E1, typename E2>
struct expr::is_nonlinear<OpBinaryAdd<E1, E2>>
{
	static const bool value = expr::is_nonlinear<E1>::value && expr::is_nonlinear<E2>::value;
};

//! Specialization based on expr::is_nonlinear.
template<typename E1, typename E2>
struct expr::is_nonlinear<OpBinarySub<E1, E2>>
{
	static const bool value = expr::is_nonlinear<E1>::value && expr::is_nonlinear<E2>::value;
};





//! Specialization based on expr::is_combination.
template<typename E1, typename E2>
struct expr::is_combination<OpBinaryAdd<E1, E2>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_combination.
template<typename E1, typename E2>
struct expr::is_combination<OpBinarySub<E1, E2>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_combination.
template<typename A1, typename A2, typename E>
struct expr::is_combination<OpCombination<A1, A2, E>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_combination.
template<typename A1, typename A2, typename E>
struct expr::is_combination<OpChain<A1, A2, E>>
{
	static const bool value = expr::is_combination<A1>::value || expr::is_combination<A2>::value || expr::is_combination<E>::value;
};


//! Specialization based on expr::expression_predicate.
template<typename E>
struct expr::expression_predicate : expr::is_linear<E>, expr::is_combination<E>
{
	static const bool linear = expr::is_linear<E>::value;
	static const bool nonlinear = expr::is_nonlinear<E>::value;
	static const bool combination = expr::is_combination<E>::value;
	static const bool operator_like = expr::is_operator_like<E>::value;
};

template<Axis ax, size_t O, typename V, typename Sp>
struct expr::is_directional_derivative<OpOperatorDirectionalDerivative<ax, O, V, Sp>>
{
	static const bool value = true;
};


template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_directional_derivative<OpFuncDerivative<Dd, V, E, Sp>>
{
	static const bool value = expr::is_directional_derivative<Dd>::value;
};


namespace expr
{


	//! Returns underlying data of the expression.
	/*!
	 * Get underlying data related to the expression passed to the member
	 * function. Useful for getting the data that the expression is applied
	 * to in a standardized way. Typically not used for modifying the
	 * expression.
	 */
	struct compound_get
	{

		//! Get the expression that the OpMap applies to.
		template<typename G, typename V, typename E>
		static auto expr(const OpMap<G, V, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncDerivative applies to.
		template<typename Dd, typename V, typename G, typename Sp>
		static auto expr(const OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>& e)
		{
			return OpLVariable<OpIdentity, G>(e.data);
		}

		//! Get the expression that the OpFuncDerivative applies to.
		template<typename Dd, typename V, typename E, typename Sp>
		static auto& expr(const OpFuncDerivative<Dd, V, E, Sp>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpChain applies to.
		template<typename A1, typename A2, typename E>
		static auto& expr(const OpChain<A1, A2, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpCombination applies to.
		template<typename A1, typename A2, typename E>
		static auto& expr(const OpCombination<A1, A2, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, typename E, typename F, typename... Args>
		static auto expr(const OpFunc<V, E, F, Args...>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, size_t D, typename G>
		static auto expr(const OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e)
		{
			return OpLVariable<OpIdentity, G>(e.data);
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, size_t D, typename E>
		static auto& expr(const OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			return e.e;
		}


		//! Get the expression that the OpExponential applies to.
		template<typename V, typename E>
		static auto& expr(const OpExponential<V, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpExponential applies to.
		template<typename V, typename S, typename G>
		static auto expr(const OpExponential<V, OpLVariable<S, G>>& e)
		{
			return OpLVariable<S, G>(e.pow, e.data);
		}





		//! Get the expression that the OpMap applies to.
		template<typename G, typename V, typename E>
		static auto expr(OpMap<G, V, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncDerivative applies to.
		template<typename Dd, typename V, typename G, typename Sp>
		static auto expr(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>& e)
		{
			return OpLVariable<OpIdentity, G>(e.data);
		}

		//! Get the expression that the OpFuncDerivative applies to.
		template<typename Dd, typename V, typename E, typename Sp>
		static auto& expr(OpFuncDerivative<Dd, V, E, Sp>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpChain applies to.
		template<typename A1, typename A2, typename E>
		static auto& expr(OpChain<A1, A2, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpCombination applies to.
		template<typename A1, typename A2, typename E>
		static auto& expr(OpCombination<A1, A2, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, typename E, typename F, typename... Args>
		static auto expr(OpFunc<V, E, F, Args...>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, size_t D, typename G>
		static auto expr(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e)
		{
			return OpLVariable<OpIdentity, G>(e.data);
		}

		//! Get the expression that the OpFuncConvolution applies to.
		template<typename V, size_t D, typename E>
		static auto& expr(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpExponential applies to.
		template<typename V, typename E>
		static auto& expr(OpExponential<V, E>& e)
		{
			return e.e;
		}

		//! Get the expression that the OpExponential applies to.
		template<typename V, typename S, typename G>
		static auto expr(OpExponential<V, OpLVariable<S, G>>& e)
		{
			return OpLVariable<S, G>(e.pow, e.data);
		}




		//! Get the grid storing the underlying data of the OpMap.
		template<typename G, typename V, typename E>
		static auto& grid(OpMap<G, V, E>& e)
		{
			return e.result;
		}

		//! Get the grid storing the underlying data of the OpFuncDerivative.
		template<typename Dd, typename V, typename G, typename Sp>
		static G& grid(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>& e)
		{
			return e.grid;
		}

		//! Get the grid storing the underlying data of the OpFuncDerivative.
		template<typename Dd, typename V, typename E, typename Sp>
		static auto& grid(OpFuncDerivative<Dd, V, E, Sp>& e)
		{
			return e.grid;
		}

		//! Get the grid storing the underlying data of the OpCombination.
		template<typename A1, typename A2, typename E>
		static auto& grid(OpCombination<A1, A2, E>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpCombination.
		template<typename A1, typename A2, typename T, typename G>
		static auto& grid(OpCombination<A1, A2, OpLVariable<T, G>>& e)
		{
			return e.e.data;
		}

		//! Get the grid storing the underlying data of the OpExponential.
		template<typename V, typename E>
		static auto& grid(OpExponential<V, E>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpExponential.
		template<typename V, typename S, typename G>
		static auto& grid(OpExponential<V, OpLVariable<S, G>>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpFuncConvolution.
		template<typename V, size_t D, typename G>
		static auto& grid(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e)
		{
			return e.g0;
		}

		//! Get the grid storing the underlying data of the OpFuncConvolution.
		template<typename V, size_t D, typename E>
		static auto& grid(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			return e.g0;
		}

		//! Get the grid storing the underlying data of the OpChain.
		template<typename A1, typename A2, typename E>
		static auto& grid(OpChain<A1, A2, E>& e)
		{
			return grid(e.outer);
		}



		//! Get the grid storing the underlying data of the OpMap.
		template<typename G, typename V, typename E>
		static auto& grid(const OpMap<G, V, E>& e)
		{
			return e.result;
		}

		//! Get the grid storing the underlying data of the OpFuncDerivative.
		template<typename Dd, typename V, typename E, typename Sp>
		static auto& grid(const OpFuncDerivative<Dd, V, E, Sp>& e)
		{
			return e.grid;
		}

		//! Get the grid storing the underlying data of the OpFuncDerivative.
		template<typename Dd, typename V, typename G, typename Sp>
		static auto& grid(const OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpCombination.
		template<typename A1, typename A2, typename E>
		static auto& grid(const OpCombination<A1, A2, E>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpCombination.
		template<typename A1, typename A2, typename T, typename G>
		static auto& grid(const OpCombination<A1, A2, OpLVariable<T, G>>& e)
		{
			return e.e.data;
		}

		//! Get the grid storing the underlying data of the OpExponential.
		template<typename V, typename E>
		static auto& grid(const OpExponential<V, E>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpExponential.
		template<typename V, typename S, typename G>
		static auto& grid(const OpExponential<V, OpLVariable<S, G>>& e)
		{
			return e.data;
		}

		//! Get the grid storing the underlying data of the OpFuncConvolution.
		template<typename V, size_t D, typename G>
		static auto& grid(const OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e)
		{
			return e.g0;
		}

		//! Get the grid storing the underlying data of the OpFuncConvolution.
		template<typename V, size_t D, typename E>
		static auto& grid(const OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			return e.g0;
		}

		//! Get the grid storing the underlying data of the OpChain.
		template<typename A1, typename A2, typename E>
		static auto& grid(const OpChain<A1, A2, E>& e)
		{
			return grid(e.outer);
		}

	};
}

namespace expr::property
{

	//! Get the unique list of variables appearing in the given expression.
	/*!
	 * Get the unique list of variables appearing in the given expression.
	 * 
	 * \tparam E The expression type that is parsed.
	 */
	template<typename E>
	struct vars
	{
		//! Returns a unique index list of the set of IDs.
		static constexpr auto get_ids()
		{
			return std::index_sequence<>{};
		}

		//! Returns whether the expression has at least one instance of an ID.
		/*!
		 * \tparam Y The ID that is checked if there is one instance of.
		 */
		template<size_t Y>
		static constexpr bool has_id()
		{
			return false;
		}

		//! Returns whether the given ID is the only one.
		/*!
		 * \tparam Y The ID that is checked if it is the only one.
		 */
		template<size_t Y>
		static constexpr auto only_id()
		{
			return has_id<Y>();
		}

		//! Returns whether the ID appears in every combination.
		/*!
		 * A combination is anything that satisfies the predicate evaluated
		 * by expr::is_combination.
		 * 
		 * \tparam Y The ID that is checked if it is the only one.
		 */
		template<size_t Y>
		static constexpr auto each_id()
		{
			return has_id<Y>();
		}
	};
}


//! Specialization based on expr::property::vars.
template<size_t Z, typename G>
struct expr::property::vars<Variable<Z, G>>
{
	static constexpr auto get_ids()
	{
		return std::index_sequence<Z>{};
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return Y == Z;
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return has_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename T, typename G>
struct expr::property::vars<OpLVariable<T, G>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<G>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<G>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<G>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return has_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename T, typename... Gs>
struct expr::property::vars<OpNLVariable<T, Gs...>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<Gs>::get_ids()...));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return ((expr::property::vars<Gs>::template has_id<Y>() || ...));
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return ((expr::property::vars<Gs>::template has_id<Y>() && ...));
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return has_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<auto f, typename V, typename E>
struct expr::property::vars<OpFuncApply<f, V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::property::vars.
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct expr::property::vars<OpFunc<V, E, F, Arg0, Args...>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename G, typename V, typename E>
struct expr::property::vars<OpMap<G, V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::property::vars<OpFuncDerivative<Dd, V, E, Sp>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::property::vars.
template<typename A1, typename A2, typename E>
struct expr::property::vars<OpChain<A1, A2, E>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::property::vars.
template<typename A1, typename A2, typename E>
struct expr::property::vars<OpCombination<A1, A2, E>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::property::vars.
template<typename V, typename E1, typename E2>
struct expr::property::vars<OpFuncConvolution<V, E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<E1>::get_ids(), expr::property::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E1>::template has_id<Y>() || expr::property::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E1>::template only_id<Y>() && expr::property::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E1>::template each_id<Y>() && expr::property::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename V, size_t D, typename E>
struct expr::property::vars<OpFuncConvolution<V, GaussianSmoothing<D>, E>>
{
	static constexpr auto get_ids()
	{
		return expr::property::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename E1, typename E2>
struct expr::property::vars<OpBinaryAdd<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<E1>::get_ids(), expr::property::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E1>::template has_id<Y>() || expr::property::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E1>::template only_id<Y>() && expr::property::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E1>::template each_id<Y>() && expr::property::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename E1, typename E2>
struct expr::property::vars<OpBinarySub<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<E1>::get_ids(), expr::property::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E1>::template has_id<Y>() || expr::property::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E1>::template only_id<Y>() && expr::property::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E1>::template each_id<Y>() && expr::property::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename E1, typename E2>
struct expr::property::vars<OpBinaryMul<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<E1>::get_ids(), expr::property::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E1>::template has_id<Y>() || expr::property::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E1>::template only_id<Y>() && expr::property::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E1>::template each_id<Y>() || expr::property::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::property::vars.
template<typename E1, typename E2>
struct expr::property::vars<OpBinaryDiv<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::property::vars<E1>::get_ids(), expr::property::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::property::vars<E1>::template has_id<Y>() || expr::property::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::property::vars<E1>::template only_id<Y>() && expr::property::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::property::vars<E1>::template each_id<Y>() || expr::property::vars<E2>::template each_id<Y>();
	}
};





//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2>
struct expr::property::derivative_index<L, OpOperatorCombination<A1, A2>>
{
protected:
	static const size_t raw_value = 
		(expr::property::derivative_index<L, A1>::value 
			| expr::property::derivative_index<L, A2>::value);

public:
	static const size_t value = raw_value;
};

//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2>
struct expr::property::derivative_index<L, OpOperatorChain<A1, A2>>
{
protected:

	static const size_t raw_value_1 = expr::property::derivative_index<L, A1>::value;
	static const size_t raw_value_2 = expr::property::derivative_index<L, A2>::value;

	static constexpr size_t get_value()
	{
		size_t i = 0;
		size_t shift_value = raw_value_2;
		size_t final_value = 0;
		do
		{
			final_value |= (shift_value & 1ull) ? (raw_value_1 << i) : 0;
			++i;
		} while (shift_value >>= 1ull);
		return final_value;
	}

public:
	static const size_t value = get_value();
};


//! Specialization based on expr::property::derivative_index.
template<size_t L, typename A1, typename A2, typename E>
struct expr::property::derivative_index<L, OpCombination<A1, A2, E>>
{
public:
	static const size_t value = expr::property::derivative_index<L, OpOperatorCombination<A1, A2>>::value;
};

//! Specialization based on expr::property::derivative_index.
template<size_t L, typename A1, typename A2, typename E>
struct expr::property::derivative_index<L, OpChain<A1, A2, E>>
{
public:
	static const size_t value = expr::property::derivative_index<L, OpOperatorChain<A1, A2>>::value;
};



//! Specialization based on expr::property::derivative_index.
template<size_t L, size_t O, typename V, typename Sp>
struct expr::property::derivative_index<L, OpOperatorDerivative<O, V, Sp>>
{
	static const size_t value = expr::property::derivative_index_raw<O>::value;
};


//! Specialization based on expr::property::derivative_index.
template<size_t L, typename Dd, typename V, typename E, typename Sp>
struct expr::property::derivative_index<L, OpFuncDerivative<Dd, V, E, Sp>>
{
protected:
	static const size_t raw_value = expr::property::derivative_index_raw<Dd::order>::value;

public:
	static const size_t value = raw_value;
};

//! Specialization based on expr::property::derivative_index.
template<size_t L, typename E1, typename E2>
struct expr::property::derivative_index<L, OpBinaryAdd<E1, E2>>
{
protected:
	static const size_t value_1 = fixed_min<expr::property::derivative_index<L, E1>::value, expr::property::derivative_index<L, E2>::value>;
	static const size_t value_2 = fixed_max<expr::property::derivative_index<L, E1>::value, expr::property::derivative_index<L, E2>::value>;

public:
	static const size_t value = (value_1 < L) ? value_2 : value_1;
};

template<size_t L, typename E1, typename E2>
struct expr::property::derivative_index<L, OpBinarySub<E1, E2>>
{
	static const size_t value = expr::property::derivative_index<L, OpBinaryAdd<E1, E2>>::value;
};




// ******************************************************************************************

/*
 * evaluate the expression in every index in the array
 */

namespace expr
{

#ifdef MULTITHREAD

	template<typename E, typename T>
	void result_interior(OpExpression<E> const& e, T* values, iter_type* inners, len_type len)
	{
		std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par, 
#endif
			inners, inners + len, [&](iter_type index) { values[index] = static_cast<const E*>(&e)->eval(index); });
	}

	template<typename T, size_t D, typename E>
	void result_interior(OpExpression<E> const& e, Grid<T, D>& grid)
	{
		expr::result_interior(*static_cast<const E*>(&e), grid.values, grid::interior_indices_list<D>(grid.dims), grid::length_interior<D>(grid.dims));
	}

#else

	//! Evaluate the expression into the interior of the array.
	/*!
	 * Evaluate and store the result of the expression on the interior values
	 * of a grid logically represented by the given dimensions. The direct
	 * array pointer is given. The expression must be of the same
	 * dimensionality as the array.
	 * 
	 * \param e The expression which is evaluated into the array.
	 * \param values The array which contains the result of the expression on
	 * its interior values.
	 * \param dim The logical dimensions of the grid that the expression is
	 * evaluated over.
	 */
	template<typename E, typename T>
	void result_interior(OpExpression<E> const& e, T* values, len_type(&dim)[3])
	{
		ITER_GRID3(values[INDEX] = static_cast<const E*>(&e)->eval(INDEX), dim[0], dim[1], dim[2])
	}

	//! Specialization based on result_interior(OpExpression<E> const&, T*, len_type(&)[3]).
	template<typename E, typename T>
	void result_interior(OpExpression<E> const& e, T* values, len_type(&dim)[2])
	{
		ITER_GRID2(values[INDEX] = static_cast<const E*>(&e)->eval(INDEX), dim[0], dim[1])
	}

	//! Specialization based on result_interior(OpExpression<E> const&, T*, len_type(&)[3]).
	template<typename E, typename T>
	void result_interior(OpExpression<E> const& e, T* values, len_type(&dim)[1])
	{
		ITER_GRID1(values[INDEX] = static_cast<const E*>(&e)->eval(INDEX), dim[0])
	}

	//! Evaluate the expression into the interior of the array.
	/*!
	 * See result_interior(OpExpression<E> const&, T*, len_type(&)[3]).
	 * Implementation based on parameter input of a Grid.
	 * Given the grid, the interior values are initialized to the result of
	 * the expression at each respective point.
	 * 
	 * \param e The expression which is evaluated into the array.
	 * \param grid The grid that contains the result of the expression on its
	 * interior values.
	 */
	template<typename E, typename T, size_t D>
	void result_interior(OpExpression<E> const& e, Grid<T, D>& grid)
	{
		expr::result_interior(*static_cast<const E*>(&e), grid.values, grid.dims);
	}

#endif



	//! Evaluate the expression into the underlying data member.
	/*!
	 * The expression must be iterable over the entire given length.
	 * 
	 * \param e Expression that is evaluated.
	 * \param data The array containing the result of the expression.
	 * \param len The number of elements in the array.
	 */
	template<typename E, typename T>
	void result(OpExpression<E> const& e, T* data, len_type len)
	{
#if defined(EXECUTION_HEADER_AVAILABLE)
		std::copy(std::execution::par, e.begin(), e.end(len), data);
#else
		for (iter_type i = 0; i < len; i++)
		{
			data[i] = static_cast<const E*>(&e)->eval(i);
		}
#endif
	}
	
	//! See result(OpExpression<E> const&, T*, len_type).
	/*!
	 * The given expression is evaluated at every point of the grid and the
	 * result stored in grid array.
	 * 
	 * \param e The expression which is evaluated.
	 * \param grid The grid that contains the result of the expression.
	 */
	template<typename E, size_t D, typename T>
	void result(OpExpression<E> const& e, Grid<T, D>& grid)
	{
		expr::result(*static_cast<const E*>(&e), grid.values, grid.len);
	}

	//! Specialization based on result(OpExpression<E> const&, Grid<T, D>&).
	template<typename E, typename T>
	void result(OpExpression<E> const& e, Block<T>& grid)
	{
		expr::result(*static_cast<const E*>(&e), grid.values, grid.len);
	}

	//! Specialization based on result(OpExpression<E> const&, Grid<T, D>&).
	template<typename G, typename E>
	void result(std::pair<G, E> const& r)
	{
		auto&& [grid, expression] = r;
		expr::result(expression, grid);
	}

	//! Add the result of the expression into the underlying data member.
	/*!
	 * The expression is evaluated and the result is added to the existing
	 * values in the data array.
	 * 
	 * \param e Expression that is evaluated.
	 * \param data The array of data.
	 * \param len The length of the array.
	 */
	template<typename E, typename T>
	void result_accumulate(OpExpression<E> const& e, T* data, len_type len)
	{
#if defined(MULTITHREAD) && defined(EXECUTION_HEADER_AVAILABLE)
		std::transform(std::execution::par, e.begin(), e.end(len), data, data, [](auto expr_value, auto data_value) { return data_value + expr_value; });
#else

		for (iter_type i = 0; i < len; i++)
		{
			data[i] += static_cast<const E*>(&e)->eval(i);
		}
#endif
	}
}


