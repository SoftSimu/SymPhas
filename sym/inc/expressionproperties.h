
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


		//! Specialization based on expr::data_len_data().
		template<typename G>
		len_type data_len_data(symphas::ref<G> const& data);
		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		len_type data_len_data(NamedData<G> const& data);
		//! Specialization based on expr::data_len_data().
		template<size_t Z, typename G>
		len_type data_len_data(Variable<Z, G> const& data);
		//! Specialization based on expr::data_len_data().
		template<Axis ax, typename G>
		len_type data_len_data(VectorComponent<ax, G> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename T>
		len_type data_len_data(SymbolicData<T> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename T>
		len_type data_len_data(SymbolicDataArray<T> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename... Ts>
		len_type data_len_data(SymbolicTermArray<Ts...> const& data);


		//! Specialization based on expr::data_len_data().
		template<typename... Ts>
		len_type data_len_data(SymbolicTermArray<Ts...> const& data, std::index_sequence<>)
		{
			return 0;
		}

		//! Specialization based on expr::data_len_data().
		template<typename... Ts, size_t I0, size_t... Is>
		len_type data_len_data(SymbolicTermArray<Ts...> const& data, std::index_sequence<I0, Is...>)
		{
			len_type len = data_len_data(data.data[I0]);
			return (len > 0) ? len : data_len_data(data, std::index_sequence<Is...>{});
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

		//! Specialization based on expr::data_len_data().
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

		//! Specialization based on expr::data_len_data().
		template<typename T>
		len_type data_len_data(SymbolicData<T> const& data)
		{
			return data_len_data(*data.data);
		}

		//! Specialization based on expr::data_len_data().
		template<typename T>
		len_type data_len_data(SymbolicDataArray<T> const& data)
		{
			return data_len_data(*data.data);
		}

		//! Specialization based on expr::data_len_data().
		template<typename... Ts>
		len_type data_len_data(SymbolicTermArray<Ts...> const& data)
		{
			return data_len_data(data, std::make_index_sequence<sizeof...(Ts)>{});
		}


		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(GridData<T, 1> const* data)
		{
			return { data->dims[0] };
		}
		
		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(GridData<T, 2> const* data)
		{
			return { data->dims[0], data->dims[1] };
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(GridData<T, 3> const* data)
		{
			return { data->dims[0], data->dims[1], data->dims[2] };
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(Grid<T, 1> const* data)
		{
			return { data->dims[0] };
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(Grid<T, 2> const* data)
		{
			return { data->dims[0], data->dims[1] };
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(Grid<T, 3> const* data)
		{
			return { data->dims[0], data->dims[1], data->dims[2] };
		}

		//! The dimension of a typical data object is undefined.
		grid::dim_list data_dimensions_cast(...)
		{
			return {};
		}

		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T>
		grid::dim_list data_dimensions_cast(SymbolicData<T> const* data)
		{
			return data_dimensions_cast(data->data);
		}

		//! Specialization based on expr::data_len_data().
		template<typename G>
		grid::dim_list data_dimensions_data(symphas::ref<G> const& data);
		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		grid::dim_list data_dimensions_data(NamedData<G> const& data);
		//! Specialization based on expr::data_len_data().
		template<size_t Z, typename G>
		grid::dim_list data_dimensions_data(Variable<Z, G> const& data);
		//! Specialization based on expr::data_len_data().
		template<Axis ax, typename G>
		grid::dim_list data_dimensions_data(VectorComponent<ax, G> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename T>
		grid::dim_list data_dimensions_data(SymbolicData<T> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename T>
		grid::dim_list data_dimensions_data(SymbolicDataArray<T> const& data);
		//! Specialization based on expr::data_len_data().
		template<typename... Ts>
		grid::dim_list data_dimensions_data(SymbolicTermArray<Ts...> const& data);




		//! Determines the dimensions of the data.
		/*!
		 * Determines the data_len by attempting to
		 * implicitly cast the given type to a one compatible with getting
		 * the data_len.
		 */
		template<typename T>
		grid::dim_list data_dimensions_data(T const& data)
		{
			return data_dimensions_cast(&data);
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		grid::dim_list data_dimensions_data(symphas::ref<G> const& data)
		{
			return data_dimensions_data(data.get());
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename G>
		grid::dim_list data_dimensions_data(NamedData<G> const& data)
		{
			return data_dimensions_data(static_cast<G const&>(data));
		}

		//! Specialization based on expr::data_dimensions_data().
		template<size_t Z, typename G>
		grid::dim_list data_dimensions_data(Variable<Z, G> const& data)
		{
			return data_dimensions_data(*static_cast<G const*>(&data));
		}

		//! Specialization based on expr::data_dimensions_data().
		template<Axis ax, typename G>
		grid::dim_list data_dimensions_data(VectorComponent<ax, G> const& data)
		{
			return data_dimensions_data(*static_cast<G const*>(&data));
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename T>
		grid::dim_list data_dimensions_data(SymbolicData<T> const& data)
		{
			return data_dimensions_data(*data.data);
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename T>
		grid::dim_list data_dimensions_data(SymbolicDataArray<T> const& data)
		{
			return data_dimensions_data(*data.data);
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename... Ts>
		grid::dim_list data_dimensions_data(SymbolicTermArray<Ts...> const& data, std::index_sequence<>)
		{
			return {};
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename... Ts, size_t I0, size_t... Is>
		grid::dim_list data_dimensions_data(SymbolicTermArray<Ts...> const& data, std::index_sequence<I0, Is...>)
		{
			auto dims = data_dimensions_data(data.data[I0]);
			return (dims.n > 0) ? dims : data_dimensions_data(data, std::index_sequence<Is...>{});
		}

		//! Specialization based on expr::data_dimensions_data().
		template<typename... Ts>
		grid::dim_list data_dimensions_data(SymbolicTermArray<Ts...> const& data)
		{
			return data_dimensions_data(data, std::make_index_sequence<sizeof...(Ts)>{});
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
	template<typename Dd, typename V, typename E, typename Sp>
	auto data_list(OpDerivative<Dd, V, E, Sp> const& e);
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
	template<expr::exp_key_t X, typename V, typename E>
	auto data_list(OpPow<X, V, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<auto f, typename V, typename E>
	auto data_list(OpFunctionApply<f, V, E> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto data_list(OpFunction<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::vars.
	template<typename V, typename sub_t, typename E, typename... Ts>
	auto data_list(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e);
	//! Specialization based on expr::vars.
	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	auto data_list(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, typename E1, typename E2>
	auto data_list(OpConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<size_t D>
	auto data_list(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_list(E const&).
	template<typename V, size_t D, typename E>
	auto data_list(OpConvolution<V, GaussianSmoothing<D>, E> const& e);
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


	template<typename T>
	auto data_list_check(SymbolicData<T> const& e)
	{
		return data_list_check(*e.data);
	}

	template<typename T>
	auto data_list_check(T const& e)
	{
		return e;
	}

	inline auto data_list(SymbolicDataArray<expr::symbols::Symbol> const& e)
	{
		return std::make_tuple();
	}

	template<typename T>
	auto data_list(SymbolicDataArray<T> const& e)
	{
		return std::make_tuple(*e.data);
	}

	template<typename... Ts, size_t... Is>
	auto data_list(SymbolicTermArray<Ts...> const& e, std::index_sequence<Is...>)
	{
		return std::make_tuple(data_list_check(e.data[Is])...);
	}

	template<typename... Ts>
	auto data_list(SymbolicTermArray<Ts...> const& e)
	{
		return data_list(e, std::make_index_sequence<sizeof...(Ts)>{});
	}

	template<typename... Ts, size_t... Is>
	auto data_list(Substitution<SymbolicDataArray<Ts>...> const& e, std::index_sequence<Is...>)
	{
		return std::tuple_cat(expr::get<Is>(e)...);
	}


	template<typename... Es, size_t... Is>
	auto data_list(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return std::tuple_cat(data_list(expr::get<Is>(e))...);
	}

	template<typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	auto data_list(OpTerms<Term<Gs, Xs>...> const& e, std::index_sequence<Is...>)
	{
		return std::make_tuple(data_list_check(expr::get<Is>(e).data())...);
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
		if constexpr (expr::has_coeff<OpTerms<V, Term<Gs, Xs>...>>)
		{
			return data_list(expr::terms_after_first(e), std::make_index_sequence<sizeof...(Gs)>{});
		}
		else
		{
			return data_list(e, std::make_index_sequence<sizeof...(Gs)>{});
		}
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto data_list(OpDerivative<Dd, V, E, Sp> const& e)
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

	template<expr::exp_key_t X, typename V, typename E>
	auto data_list(OpPow<X, V, E> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<auto f, typename V, typename E>
	auto data_list(OpFunctionApply<f, V, E> const& e)
	{
		return data_list(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto data_list(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return data_list(e.e);
	}

	template<typename V, typename sub_t, typename E, typename... Ts>
	auto data_list(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		return data_list(expr::get_enclosed_expression(e));
	}

	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	auto data_list(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e)
	{
		return std::tuple_cat(data_list(e.data.e), data_list(e.data.substitution, std::make_index_sequence<sizeof...(Ts)>{}));
	}

	template<typename V, typename E1, typename E2>
	auto data_list(OpConvolution<V, E1, E2> const& e)
	{
		return data_list(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	auto data_list(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
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
	grid::dim_list data_dimensions(E const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E>
	grid::dim_list data_dimensions(OpExpression<E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	grid::dim_list data_dimensions(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	grid::dim_list data_dimensions(OpDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	grid::dim_list data_dimensions(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	grid::dim_list data_dimensions(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E>
	grid::dim_list data_dimensions(OpExponential<V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<expr::exp_key_t X, typename V, typename E>
	grid::dim_list data_dimensions(OpPow<X, V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<auto f, typename V, typename E>
	grid::dim_list data_dimensions(OpFunctionApply<f, V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	grid::dim_list data_dimensions(OpFunction<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename sub_t, typename E, typename... Ts>
	grid::dim_list data_dimensions(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	grid::dim_list data_dimensions(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, typename E1, typename E2>
	grid::dim_list data_dimensions(OpConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename V, size_t D, typename E>
	grid::dim_list data_dimensions(OpConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename G, typename V, typename E>
	grid::dim_list data_dimensions(OpMap<G, V, E> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<size_t D>
	grid::dim_list data_dimensions(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E0, typename... Es>
	grid::dim_list data_dimensions(OpAdd<E0, Es...> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::data_dimensions(E const&).
	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the dimensions of the data in the expression.
	/*!
	 * Obtain the dimensions from either of the two expressions. 
	 * Chooses the appropriate expr::data_dimensions() between two expressions.
	 * This standardizes the choosing procedure for binary expressions by
	 * checking if one expression is `nullptr`, and subsequently getting
	 * the dimensions from the other.
	 */
	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpExpression<E1> const& a, OpExpression<E2> const& b);


	template<typename E>
	grid::dim_list data_dimensions(E const& e)
	{
		return data_dimensions_data(e);
	}

	template<typename E>
	grid::dim_list data_dimensions(OpExpression<E> const& e)
	{
		return data_dimensions(*static_cast<E const*>(&e));
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	grid::dim_list data_dimensions(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		if constexpr (expr::has_coeff<OpTerms<V, Term<Gs, Xs>...>>)
		{
			return data_dimensions(expr::terms_after_first(e));
		}
		else
		{
			auto data_dimensions1 = data_dimensions(e.term.data());
			return (data_dimensions1.n > 0) ? data_dimensions1 : data_dimensions(expr::terms_after_first(e));
		}
	}

	template<typename Dd, typename V, typename E, typename Sp>
	grid::dim_list data_dimensions(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	grid::dim_list data_dimensions(OpCombination<A1, A2, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	grid::dim_list data_dimensions(OpChain<A1, A2, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename E>
	grid::dim_list data_dimensions(OpExponential<V, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<expr::exp_key_t X, typename V, typename E>
	grid::dim_list data_dimensions(OpPow<X, V, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<auto f, typename V, typename E>
	grid::dim_list data_dimensions(OpFunctionApply<f, V, E> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	grid::dim_list data_dimensions(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return data_dimensions(e.e);
	}

	template<typename V, typename sub_t, typename E, typename... Ts>
	grid::dim_list data_dimensions(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	grid::dim_list data_dimensions(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e)
	{
		auto dims = data_dimensions(expr::get_enclosed_expression(e));
		return (dims.n > 0) ? dims : data_dimensions(e.data.substitution);
	}

	template<typename V, typename E1, typename E2>
	grid::dim_list data_dimensions(OpConvolution<V, E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	grid::dim_list data_dimensions(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e), e.smoother);
	}

	template<typename G, typename V, typename E>
	grid::dim_list data_dimensions(OpMap<G, V, E> const& e)
	{
		return data_dimensions(expr::get_enclosed_expression(e));
	}

	template<size_t D>
	grid::dim_list data_dimensions(GaussianSmoothing<D> const& e)
	{
		return data_dimensions_data(e.data);
	}

	template<typename E0, typename... Es>
	grid::dim_list data_dimensions(OpAdd<E0, Es...> const& e)
	{
		auto dims = data_dimensions(expr::get<0>(e));
		return (dims.n > 0) ? dims : data_dimensions(expr::terms_after_first(e));
	}

	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpBinaryMul<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpBinaryDiv<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename E1, typename E2>
	grid::dim_list data_dimensions(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		auto dims = data_dimensions(*static_cast<const E1*>(&a));
		return (dims.n > 0) ? dims : data_dimensions(*static_cast<const E2*>(&b));
	}

	template<typename E, size_t D>
	void fill_data_dimensions(OpExpression<E> const& a, len_type (&dims)[D])
	{
		auto data_dims = data_dimensions(*static_cast<const E*>(&a));
		for (iter_type i = 0; i < D; ++i)
		{
			dims[i] = data_dims[i];
		}
	}


	//! Obtain the length of the data in the expression.
	/*!
	 * Return the data_len of the underlying data set.
	 * This is particularly important for data sets which contain grids.
	 * If the data set does not contain a grid, then a default data_len is
	 * returned.
	 */
	template<typename E>
	len_type data_length(E const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename E>
	len_type data_length(OpExpression<E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	len_type data_length(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_length(OpDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename Dd, typename V, typename G, typename Sp>
	len_type data_length(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_length(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_length(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename E>
	len_type data_length(OpExponential<V, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<expr::exp_key_t X, typename V, typename E>
	len_type data_length(OpPow<X, V, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<auto f, typename V, typename E>
	len_type data_length(OpFunctionApply<f, V, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	len_type data_length(OpFunction<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename sub_t, typename E, typename... Ts>
	len_type data_length(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	len_type data_length(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, typename E1, typename E2>
	len_type data_length(OpConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<size_t D>
	len_type data_length(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename V, size_t D, typename E>
	len_type data_length(OpConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename G, typename V, typename E>
	len_type data_length(OpMap<G, V, E> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename... Es>
	len_type data_length(OpAdd<Es...> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename E1, typename E2>
	len_type data_length(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::data_length(E const&).
	template<typename E1, typename E2>
	len_type data_length(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the length of the data in the expression.
	/*!
	 * Return the data_len of the underlying data set.
	 * Chooses the appropriate expr::data_length() between two expressions.
	 * This standardizes the choosing procedure for binary expressions.
	 */
	template<typename E1, typename E2>
	len_type data_length(OpExpression<E1> const& a, OpExpression<E2> const& b);

	template<typename... Es, size_t... Is>
	len_type data_length(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return std::max({ data_length(expr::get<Is>(e))... });
	}

	template<typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	len_type data_length(OpTerms<Term<Gs, Xs>...> const& e, std::index_sequence<Is...>)
	{
		return std::max({ 1, data_len_data(expr::get<Is>(e).data())... });
	}

	template<typename E>
	len_type data_length(E const& e)
	{
		return data_len_data(e);
	}

	template<typename E>
	len_type data_length(OpExpression<E> const& e)
	{
		return data_length(*static_cast<E const*>(&e));
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	len_type data_length(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		if constexpr (expr::has_coeff<OpTerms<V, Term<Gs, Xs>...>>)
		{
			return data_length(expr::terms_after_first(e), std::make_index_sequence<sizeof...(Gs)>{});
		}
		else
		{
			return data_length(e, std::make_index_sequence<sizeof...(Gs) + 1>{});
		}
	}

	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_length(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return data_len_data(expr::get_result_data(e));
	}

	template<typename Dd, typename V, typename G, typename Sp>
	len_type data_length(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
	{
		return data_len_data(expr::get_enclosed_expression(e));
	}

	template<typename A1, typename A2, typename E>
	len_type data_length(OpCombination<A1, A2, E> const& e)
	{
		return data_length(e.e);
	}

	template<typename A1, typename A2, typename E>
	len_type data_length(OpChain<A1, A2, E> const& e)
	{
		return data_length(e.e);
	}

	template<typename V, typename E>
	len_type data_length(OpExponential<V, E> const& e)
	{
		return data_length(e.e);
	}

	template<expr::exp_key_t X, typename V, typename E>
	len_type data_length(OpPow<X, V, E> const& e)
	{
		return data_length(expr::get_enclosed_expression(e));
	}

	template<auto f, typename V, typename E>
	len_type data_length(OpFunctionApply<f, V, E> const& e)
	{
		return data_length(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	len_type data_length(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return data_length(e.e);
	}

	template<typename V, typename sub_t, typename E, typename... Ts>
	len_type data_length(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		return data_length(expr::get_enclosed_expression(e));
	}

	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename E0, typename... T0s, typename B, typename C>
	len_type data_length(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, SymbolicFunction<E0, T0s...>, B, C> const& e)
	{
		auto len = data_length(expr::get_enclosed_expression(e));
		return (len > 0) ? len : data_length(e.data.substitution);
	}

	template<typename V, typename E1, typename E2>
	len_type data_length(OpConvolution<V, E1, E2> const& e)
	{
		return data_length(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	len_type data_length(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return data_length(expr::get_enclosed_expression(e), e.smoother);
	}

	template<typename G, typename V, typename E>
	len_type data_length(OpMap<G, V, E> const& e)
	{
		return data_length(expr::get_enclosed_expression(e));
	}

	template<size_t D>
	len_type data_length(GaussianSmoothing<D> const& e)
	{
		return data_len_data(e.data);
	}

	template<typename... Es>
	len_type data_length(OpAdd<Es...> const& e)
	{
		return data_length(e, std::make_index_sequence<sizeof...(Es)>{});
	}


	template<typename E1, typename E2>
	len_type data_length(OpBinaryMul<E1, E2> const& e)
	{
		return data_length(e.a, e.b);
	}

	template<typename E1, typename E2>
	len_type data_length(OpBinaryDiv<E1, E2> const& e)
	{
		return data_length(e.a, e.b);
	}

	template<typename E1, typename E2>
	len_type data_length(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return std::max(
			data_length(*static_cast<const E1*>(&a)),
			data_length(*static_cast<const E2*>(&b)));
	}


	template<typename T>
	struct vars_in_ops;

	template<>
	struct vars_in_ops<symphas::lib::types_list<>>
	{
		using type = std::index_sequence<>;
	};

	template<typename G0, typename... Gs>
	struct vars_in_ops<symphas::lib::types_list<G0, Gs...>>
	{
		using type = typename vars_in_ops<symphas::lib::types_list<Gs...>>::type;
	};

	template<size_t N, typename G0, typename... Gs>
	struct vars_in_ops<symphas::lib::types_list<Variable<N, G0>, Gs...>>
	{
		using type = symphas::lib::seq_join_t<
			std::index_sequence<N>,
			typename vars_in_ops<symphas::lib::types_list<Gs...>>::type>;
	};

	template<typename T>
	using vars_in_ops_t = typename vars_in_ops<T>::type;

	//! Returns an index sequence with all Variable IDs.
	/*!
	 * Returns an index sequence containing the ID of all Variable types within the
	 * given expression.
	 */
	template<typename E>
	auto get_independent_variables(E const& e)
	{
		return vars_in_ops_t<expr::op_types_t<E>>{};
	}

	//! Returns an index sequence with all Variable IDs.
	/*!
	 * Returns an index sequence containing the ID of all Variable types within the
	 * given expression.
	 */
	template<typename E, size_t... ArgNs>
	auto get_independent_variables(SymbolicTemplate<E, ArgNs...> const& e)
	{
		return std::index_sequence<ArgNs...>{};
	}

	template<typename E>
	constexpr auto independent_variables_of = decltype(get_independent_variables(std::declval<E>())){};

}

