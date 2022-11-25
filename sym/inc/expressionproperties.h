
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
 * PURPOSE: Defines methods of obtaining properties of expressions, in
 * particular about the underlying data of the expressions.
 *
 * ***************************************************************************
 */

#pragma once

#include <tuple>

#include "expressionlib.h"


namespace expr
{
	namespace
	{
		//! Obtains the data_len from the Block compatible instance.
		template<typename T>
		len_type data_len_cast(Block<T> const* data)
		{
			return data->len;
		}

		//! Obtains the data_len from the Block compatible instance.
		template<typename T, size_t D>
		len_type data_len_cast(GridData<T, D> const* data)
		{
			return grid::length<D>(data->dims);
		}

		//! The data_len of a typical data object is 1.
		inline len_type data_len_cast(...)
		{
			return 1;
		}

		//! Determines the data_len of the data.
		/*!
		 * Determines the data_len by attempting to 
		 * implicitly cast the given type to a one compatible with getting
		 * the data_len.
		 */
		template<typename T>
		len_type data_len_data(T const& data)
		{
			return data_len_cast(&data);
		}

		//! Specialization based on expr::data_len_data().
		template<typename G>
		len_type data_len_data(symphas::ref<G> const& data)
		{
			return data_len_data(data.get());
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		len_type data_len_data(NamedData<G> const& data)
		{
			return data_len_data(static_cast<G const&>(data));
		}

		//! Specialization based on expr::data_len_data().
		template<size_t Z, typename G>
		len_type data_len_data(Variable<Z, G> const& data)
		{
			return data_len_data(*static_cast<G const*>(&data));
		}

		//! Specialization based on expr::data_len_data().
		template<Axis ax, typename G>
		len_type data_len_data(VectorComponent<ax, G> const& data)
		{
			return data_len_data(*static_cast<G const*>(&data));
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T, size_t D>
		const len_type* data_dimensions_cast(Grid<T, D> const* data)
		{
			return data->dims;
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T, size_t D>
		const len_type* data_dimensions_cast(GridData<T, D> const* data)
		{
			return data->dims;
		}

		//! The dimension of a typical data object is undefined.
		const len_type* data_dimensions_cast(...)
		{
			return nullptr;
		}

		//! Determines the dimensions of the data.
		/*!
		 * Determines the data_len by attempting to
		 * implicitly cast the given type to a one compatible with getting
		 * the data_len.
		 */
		template<typename T>
		const len_type* data_dimensions_data(T const& data)
		{
			return data_dimensions_cast(&data);
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		const len_type* data_dimensions_data(symphas::ref<G> const& data)
		{
			return data_dimensions_data(data.get());
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		const len_type* data_dimensions_data(NamedData<G> const& data)
		{
			return data_dimensions_data(static_cast<G const&>(data));
		}

		//! Specialization based on expr::data_dimensions_data().
		template<size_t Z, typename G>
		const len_type* data_dimensions_data(Variable<Z, G> const& data)
		{
			return data_dimensions_data(*static_cast<G const*>(&data));
		}

		//! Specialization based on expr::data_dimensions_data().
		template<Axis ax, typename G>
		const len_type* data_dimensions_data(VectorComponent<ax, G> const& data)
		{
			return data_dimensions_data(*static_cast<G const*>(&data));
		}
	}




	//! Obtain the list (which may include repeats) of datas in an expression.
	/*!
	 * Returns the list of all data elements that are used in an expression,
	 * which can be any of Variable, NamedData, Grid, and anything else which
	 * is managed by the OpTerm and OpNLVariable terms.
	 */
	template<typename E>
	auto data_list(E const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename E>
	auto data_list(OpExpression<E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	auto data_list(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	auto data_list(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	auto data_list(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename A1, typename A2, typename E>
	auto data_list(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename A1, typename A2, typename E>
	auto data_list(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename E>
	auto data_list(OpExponential<V, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<auto f, typename V, typename E>
	auto data_list(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto data_list(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename E1, typename E2>
	auto data_list(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<size_t D>
	auto data_list(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, size_t D, typename E>
	auto data_list(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename G, typename V, typename E>
	auto data_list(OpMap<G, V, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename... Es>
	auto data_list(OpAdd<Es...> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename E1, typename E2>
	auto data_list(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename E1, typename E2>
	auto data_list(OpBinaryDiv<E1, E2> const& e);
	
	template<typename E1, typename E2>
	auto data_list(OpExpression<E1> const& a, OpExpression<E2> const& b);

	template<typename... Es, size_t... Is>
	auto data_list(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return std::tuple_cat(data_list(expr::get<Is>(e))...);
	}

	template<typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	auto data_list(OpTerms<Term<Gs, Xs>...> const& e, std::index_sequence<Is...>)
	{
		return std::make_tuple(expr::get<Is>(e).data()...);
	}

	template<typename E>
	auto data_list(E const& e)
	{
		return std::tuple<>{};
	}

	template<typename E>
	auto data_list(OpExpression<E> const& e)
	{
		return data_list(*static_cast<E const*>(&e));
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	auto data_list(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return data_list(expr::terms_after_first(e), std::make_index_sequence<sizeof...(Gs)>{});
	}

	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	auto data_list(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e)
	{
		return data_list(e, std::make_index_sequence<sizeof...(Gs)>{});
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto data_list(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	auto data_list(OpCombination<A1, A2, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	auto data_list(OpChain<A1, A2, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<typename V, typename E>
	auto data_list(OpExponential<V, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<auto f, typename V, typename E>
	auto data_list(OpFuncApply<f, V, E> const& e)
	{
		return data_list(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto data_list(OpFunc<V, E, F, Arg0, Args...> const& e)
	{
		return data_list(e.e);
	}

	template<typename V, typename E1, typename E2>
	auto data_list(OpFuncConvolution<V, E1, E2> const& e)
	{
		return data_list(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	auto data_list(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e), e.smoother);
	}

	template<typename G, typename V, typename E>
	auto data_list(OpMap<G, V, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<size_t D>
	auto data_list(GaussianSmoothing<D> const& e)
	{
		return std::make_tuple(e.data);
	}

	template<typename... Es>
	auto data_list(OpAdd<Es...> const& e)
	{
		return data_list(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename E1, typename E2>
	auto data_list(OpBinaryMul<E1, E2> const& e)
	{
		return data_list(e.a, e.b);
	}

	template<typename E1, typename E2>
	auto data_list(OpBinaryDiv<E1, E2> const& e)
	{
		return data_list(e.a, e.b);
	}

	//! Obtain all datas used in the expression.
	/*!
	 * Concatenate the list of data elements obtained from two expressions.
	 */
	template<typename E1, typename E2>
	auto data_list(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return std::tuple_cat(
			data_list(*static_cast<const E1*>(&a)),
			data_list(*static_cast<const E2*>(&b))
		);
	}


	namespace
	{
		template<size_t Z>
		auto get_variable_apply(std::tuple<> const&)
		{
			return nullptr;
		}

		template<size_t Z, typename G, typename... Ds>
		auto get_variable_apply(std::tuple<Variable<Z, G>, Ds...> const& datas)
		{
			return std::get<0>(datas);
		}

		template<size_t Z, Axis ax, typename G, typename... Ds>
		auto get_variable_apply(std::tuple<Variable<Z, VectorComponent<ax, G>>, Ds...> const& datas)
		{
			return Variable<Z, G>(*static_cast<G const*>(&std::get<0>(datas)));
		}

		/*template<size_t Z, Axis ax, typename G, typename... Ds>
		Variable<Z, G> get_variable_apply(std::tuple<VectorComponent<ax, Variable<Z, G>>, Ds...> const& datas)
		{
			return std::get<0>(datas);
		}*/

		template<size_t Z, typename D0, typename... Ds>
		auto get_variable_apply(std::tuple<D0, Ds...> const& datas)
		{
			return get_variable_apply<Z>(symphas::lib::get_tuple_ge<1>(datas));
		}
	}

	//! Return the variable of the specified index from the expression. 
	/*!
	 * Returns an instance of the variable from the expression.
	 */
	template<size_t Z, typename E>
	auto get_variable(OpExpression<E> const& e)
	{
		auto datas = data_list(*static_cast<E const*>(&e));
		return get_variable_apply<Z>(datas);
	}

	//! Obtain the dimensions of the data in the expression.
	/*!
	 * if the underlying expresions contain a data object with dimensions, 
	 * such as a grid, then it will return the dimensions as a pointer.
	 * The pointer to the dimensions cannot be deallocated before the 
	 * dimensions are used.
	 * This is particularly important for data sets which contain grids.
	 * If the data set does not contain a grid, then a `nullptr` is
	 * returned.
	 */
	template<typename E>
	const len_type* data_dimensions(E const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E>
	const len_type* data_dimensions(OpExpression<E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	const len_type* data_dimensions(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	const len_type* data_dimensions(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	const len_type* data_dimensions(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E>
	const len_type* data_dimensions(OpExponential<V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<auto f, typename V, typename E>
	const len_type* data_dimensions(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	const len_type* data_dimensions(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E1, typename E2>
	const len_type* data_dimensions(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, size_t D, typename E>
	const len_type* data_dimensions(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename G, typename V, typename E>
	const len_type* data_dimensions(OpMap<G, V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<size_t D>
	const len_type* data_dimensions(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E0, typename... Es>
	const len_type* data_dimensions(OpAdd<E0, Es...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the dimensions of the data in the expression.
	/*!
	 * Obtain the dimensions from either of the two expressions. 
	 * Chooses the appropriate expr::data_dimensions() between two expressions.
	 * This standardizes the choosing procedure for binary expressions by
	 * checking if one expression is `nullptr`, and subsequently getting
	 * the dimensions from the other.
	 */
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpExpression<E1> const& a, OpExpression<E2> const& b);


	template<typename E>
	const len_type* data_dimensions(E const& e)
	{
		return data_dimensions_data(e);
	}

	template<typename E>
	const len_type* data_dimensions(OpExpression<E> const& e)
	{
		return data_dimensions(*static_cast<E const*>(&e));
	}


	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	const len_type* data_dimensions(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e)
	{
		auto data_dimensions1 = data_dimensions(e.term.data());
		return (data_dimensions1 != nullptr) ? data_dimensions1 : data_dimensions(expr::terms_after_first(e));

	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	const len_type* data_dimensions(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return data_dimensions(expr::terms_after_first(e));
	}

	template<typename Dd, typename V, typename E, typename Sp>
	const len_type* data_dimensions(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpCombination<A1, A2, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpChain<A1, A2, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename E>
	const len_type* data_dimensions(OpExponential<V, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<auto f, typename V, typename E>
	const len_type* data_dimensions(OpFuncApply<f, V, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	const len_type* data_dimensions(OpFunc<V, E, F, Arg0, Args...> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename E1, typename E2>
	const len_type* data_dimensions(OpFuncConvolution<V, E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	const len_type* data_dimensions(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e), e.smoother);
	}

	template<typename G, typename V, typename E>
	const len_type* data_dimensions(OpMap<G, V, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<size_t D>
	const len_type* data_dimensions(GaussianSmoothing<D> const& e)
	{
		return data_dimensions_data(e.data);
	}

	template<typename E0, typename... Es>
	const len_type* data_dimensions(OpAdd<E0, Es...> const& e)
	{
		auto data_dimensions1 = data_dimensions(e.data);
		return (data_dimensions1 != nullptr) ? data_dimensions1 : data_dimensions(expr::terms_after_first(e));
	}

	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryMul<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryDiv<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename E1, typename E2>
	const len_type* data_dimensions(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		auto data_dimensions1 = data_dimensions(*static_cast<const E1*>(&a));
		return (data_dimensions1 != nullptr) ? data_dimensions1 : data_dimensions(*static_cast<const E2*>(&b));
	}



	//! Obtain the length of the data in the expression.
	/*!
	 * Return the data_len of the underlying data set.
	 * This is particularly important for data sets which contain grids.
	 * If the data set does not contain a grid, then a default data_len is
	 * returned.
	 */
	template<typename E>
	len_type data_len(E const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename E>
	len_type data_len(OpExpression<E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	len_type data_len(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	len_type data_len(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename Dd, typename V, typename G, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_len(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_len(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename V, typename E>
	len_type data_len(OpExponential<V, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<auto f, typename V, typename E>
	len_type data_len(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	len_type data_len(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename V, typename E1, typename E2>
	len_type data_len(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<size_t D>
	len_type data_len(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename V, size_t D, typename E>
	len_type data_len(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename G, typename V, typename E>
	len_type data_len(OpMap<G, V, E> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename... Es>
	len_type data_len(OpAdd<Es...> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the length of the data in the expression.
	/*!
	 * Return the data_len of the underlying data set.
	 * Chooses the appropriate expr::data_len() between two expressions.
	 * This standardizes the choosing procedure for binary expressions.
	 */
	template<typename E1, typename E2>
	len_type data_len(OpExpression<E1> const& a, OpExpression<E2> const& b);

	template<typename... Es, size_t... Is>
	len_type data_len(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return std::max({ data_len(expr::get<Is>(e))... });
	}

	template<typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	len_type data_len(OpTerms<Term<Gs, Xs>...> const& e, std::index_sequence<Is...>)
	{
		return std::max({ 1, data_len_data(expr::get<Is>(e).data())... });
	}

	template<typename E>
	len_type data_len(E const& e)
	{
		return data_len_data(e);
	}

	template<typename E>
	len_type data_len(OpExpression<E> const& e)
	{
		return data_len(*static_cast<E const*>(&e));
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	len_type data_len(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return data_len(expr::terms_after_first(e), std::make_index_sequence<sizeof...(Gs)>{});
	}

	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	len_type data_len(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e)
	{
		return data_len(e, std::make_index_sequence<sizeof...(Gs) + 1>{});
	}

	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return data_len_data(expr::get_result_data(e));
	}

	template<typename Dd, typename V, typename G, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
	{
		return data_len_data(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	len_type data_len(OpCombination<A1, A2, E> const& e)
	{
		return data_len(e.e);
	}

	template<typename A1, typename A2, typename E>
	len_type data_len(OpChain<A1, A2, E> const& e)
	{
		return data_len(e.e);
	}

	template<typename V, typename E>
	len_type data_len(OpExponential<V, E> const& e)
	{
		return data_len(e.e);
	}

	template<auto f, typename V, typename E>
	len_type data_len(OpFuncApply<f, V, E> const& e)
	{
		return data_len(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	len_type data_len(OpFunc<V, E, F, Arg0, Args...> const& e)
	{
		return data_len(e.e);
	}

	template<typename V, typename E1, typename E2>
	len_type data_len(OpFuncConvolution<V, E1, E2> const& e)
	{
		return data_len(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	len_type data_len(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return data_len(expr::get_enclosed_expression(e), e.smoother);
	}

	template<typename G, typename V, typename E>
	len_type data_len(OpMap<G, V, E> const& e)
	{
		return data_len(expr::get_enclosed_expression(e));
	}

	//template<typename E>
	//len_type data_len(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
	//{
	//	const len_type* dims = data_dimensions(expr::get_enclosed_expression(e));
	//	len_type len = dims[0] / 2 + 1;
	//	for (iter_type i = 1; i < expr::grid_dim<E>::value; ++i)
	//	{
	//		len *= dims[i];
	//	}
	//}

	template<size_t D>
	len_type data_len(GaussianSmoothing<D> const& e)
	{
		return data_len_data(e.data);
	}

	template<typename... Es>
	len_type data_len(OpAdd<Es...> const& e)
	{
		return data_len(e, std::make_index_sequence<sizeof...(Es)>{});
	}


	template<typename E1, typename E2>
	len_type data_len(OpBinaryMul<E1, E2> const& e)
	{
		return data_len(e.a, e.b);
	}

	template<typename E1, typename E2>
	len_type data_len(OpBinaryDiv<E1, E2> const& e)
	{
		return data_len(e.a, e.b);
	}

	template<typename E1, typename E2>
	len_type data_len(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return std::max(
			data_len(*static_cast<const E1*>(&a)),
			data_len(*static_cast<const E2*>(&b)));
	}



}

