
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

#include "expressionlib.h"


namespace expr
{
	//! Identifies properties about expressions.
	/*!
	 * Functions which identify properties about the given expression, such as
	 * the size of the dimensions.
	 */
	namespace property {}
}

namespace expr::property
{
	namespace
	{
		//! Obtains the data_len from the Block compatible instance.
		template<typename T>
		len_type data_len_cast(Block<T> const* data)
		{
			return data->len;
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

		//! Specialization based on expr::property::data_len_data().
		template<typename G>
		len_type data_len_data(symphas::ref<G> const& data)
		{
			return data_len_data(data.get());
		}

		//! Specialization based on expr::property::data_dimensions_data().
		template<typename G>
		len_type data_len_data(NamedData<G> const& data)
		{
			return data_len_data(static_cast<G const&>(data));
		}

		//! Specialization based on expr::property::data_len_data().
		template<size_t Z, typename G>
		len_type data_len_data(Variable<Z, G> const& data)
		{
			return data_len_data(*static_cast<G const*>(&data));
		}

		//! Specialization based on expr::property::data_len_data().
		template<typename T0>
		len_type data_len_data(T0 const& data0, std::tuple<> const&)
		{
			return data_len_data(data0);
		}

		//! Specialization based on expr::property::data_len_data().
		template<typename T0, typename T1, typename... Ts>
		len_type data_len_data(T0 const& data0, std::tuple<T1, Ts...> const& data)
		{
			return std::max(data_len_data(data0), data_len_data(std::get<0>(data), symphas::lib::get_tuple_ge<1>(data)));
		}


		//! Obtains the dimensions from the Grid compatible instance.
		template<typename T, size_t D>
		const len_type* data_dimensions_cast(Grid<T, D> const* data)
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

		//! Specialization based on expr::property::data_dimensions_data().
		template<typename G>
		const len_type* data_dimensions_data(symphas::ref<G> const& data)
		{
			return data_dimensions_data(data.get());
		}

		//! Specialization based on expr::property::data_dimensions_data().
		template<typename G>
		const len_type* data_dimensions_data(NamedData<G> const& data)
		{
			return data_dimensions_data(static_cast<G const&>(data));
		}

		//! Specialization based on expr::property::data_dimensions_data().
		template<size_t Z, typename G>
		const len_type* data_dimensions_data(Variable<Z, G> const& data)
		{
			return data_dimensions_data(*static_cast<G const*>(&data));
		}


		//! Specialization based on expr::property::data_dimensions_data().
		template<typename T0>
		const len_type* data_dimensions_data(T0 const& data0, std::tuple<> const&)
		{
			return data_dimensions_data(data0);
		}

		//! Specialization based on expr::property::data_dimensions_data().
		template<typename T0, typename T1, typename... Ts>
		const len_type* data_dimensions_data(T0 const& data0, std::tuple<T1, Ts...> const& data)
		{
			const len_type* d = data_dimensions_data(data0);
			return (d != nullptr) ? d : data_dimensions_data(std::get<0>(data), symphas::lib::get_tuple_ge<1>(data));
		}
	}




	//! Obtain the list (which may include repeats) of datas in an expression.
	/*!
	 * Returns the list of all data elements that are used in an expression,
	 * which can be any of Variable, NamedData, Grid, and anything else which
	 * is managed by the OpLVariable and OpNLVariable terms.
	 */
	template<typename E>
	auto get_data(E const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename E>
	auto get_data(OpExpression<E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename T, typename G>
	auto get_data(OpLVariable<T, G> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename T, typename... Gs>
	auto get_data(OpNLVariable<T, Gs...> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	auto get_data(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename A1, typename A2, typename E>
	auto get_data(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename A1, typename A2, typename E>
	auto get_data(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename V, typename E>
	auto get_data(OpExponential<V, E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<auto f, typename V, typename E>
	auto get_data(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto get_data(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename V, typename E1, typename E2>
	auto get_data(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<size_t D>
	auto get_data(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename V, size_t D, typename E>
	auto get_data(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename E1, typename E2>
	auto get_data(OpBinaryAdd<E1, E2> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename E1, typename E2>
	auto get_data(OpBinarySub<E1, E2> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename E1, typename E2>
	auto get_data(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::property::get_data(E const&).
	template<typename E1, typename E2>
	auto get_data(OpBinaryDiv<E1, E2> const& e);
	template<typename E1, typename E2>
	auto get_data(OpExpression<E1> const& a, OpExpression<E2> const& b);


	template<typename E>
	auto get_data(E const& e)
	{
		return std::tuple<>{};
	}

	template<typename E>
	auto get_data(OpExpression<E> const& e)
	{
		return get_data(*static_cast<E const*>(&e));
	}

	template<typename T, typename G>
	auto get_data(OpLVariable<T, G> const& e)
	{
		return std::make_tuple(e.data);
	}

	template<typename T, typename... Gs>
	auto get_data(OpNLVariable<T, Gs...> const& e)
	{
		return e.datas;
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto get_data(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return get_data(expr::compound_get::template expr(e));
	}

	template<typename A1, typename A2, typename E>
	auto get_data(OpCombination<A1, A2, E> const& e)
	{
		return get_data(expr::compound_get::template expr(e));
	}

	template<typename A1, typename A2, typename E>
	auto get_data(OpChain<A1, A2, E> const& e)
	{
		return get_data(expr::compound_get::template expr(e));
	}

	template<typename V, typename E>
	auto get_data(OpExponential<V, E> const& e)
	{
		return get_data(expr::compound_get::template expr(e));
	}

	template<auto f, typename V, typename E>
	auto get_data(OpFuncApply<f, V, E> const& e)
	{
		return get_data(e.e);
	}

	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto get_data(OpFunc<V, E, F, Arg0, Args...> const& e)
	{
		return get_data(e.e);
	}

	template<typename V, typename E1, typename E2>
	auto get_data(OpFuncConvolution<V, E1, E2> const& e)
	{
		return get_data(e.a, e.b);
	}

	template<typename V, size_t D, typename E>
	auto get_data(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return get_data(expr::compound_get::template expr(e), e.smoother);
	}

	template<size_t D>
	auto get_data(GaussianSmoothing<D> const& e)
	{
		return std::make_tuple(e.data);
	}

	template<typename E1, typename E2>
	auto get_data(OpBinaryAdd<E1, E2> const& e)
	{
		return get_data(e.a, e.b);
	}

	template<typename E1, typename E2>
	auto get_data(OpBinarySub<E1, E2> const& e)
	{
		return get_data(e.a, e.b);
	}

	template<typename E1, typename E2>
	auto get_data(OpBinaryMul<E1, E2> const& e)
	{
		return get_data(e.a, e.b);
	}

	template<typename E1, typename E2>
	auto get_data(OpBinaryDiv<E1, E2> const& e)
	{
		return get_data(e.a, e.b);
	}

	//! Obtain all datas used in the expression.
	/*!
	 * Concatenate the list of data elements obtained from two expressions.
	 */
	template<typename E1, typename E2>
	auto get_data(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return std::tuple_cat(
			get_data(*static_cast<const E1*>(&a)),
			get_data(*static_cast<const E2*>(&b))
		);
	}


	namespace
	{
		template<size_t Z>
		auto get_data_variable_apply(std::tuple<> const&)
		{
			return nullptr;
		}

		template<size_t Z, typename G, typename... Ds>
		auto get_data_variable_apply(std::tuple<Variable<Z, G>, Ds...> const& datas)
		{
			return std::get<0>(datas);
		}

		template<size_t Z, typename D0, typename... Ds>
		auto get_data_variable_apply(std::tuple<D0, Ds...> const& datas)
		{
			return get_data_variable_apply<Z>(symphas::lib::get_tuple_ge<1>(datas));
		}
	}

	//! Return the variable of the specified index from the expression. 
	/*!
	 * Returns an instance of the variable from the expression.
	 */
	template<size_t Z, typename E>
	auto get_data_variable(OpExpression<E> const& e)
	{
		auto datas = get_data(*static_cast<E const*>(&e));
		return get_data_variable_apply<Z>(datas);
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
	//! Specialization based on expr::property::data_len(E const&).
	template<typename E>
	len_type data_len(OpExpression<E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename T, typename G>
	len_type data_len(OpLVariable<T, G> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename T, typename... Gs>
	len_type data_len(OpNLVariable<T, Gs...> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_len(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename A1, typename A2, typename E>
	len_type data_len(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename V, typename E>
	len_type data_len(OpExponential<V, E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename V, typename T, typename G>
	len_type data_len(OpExponential<V, OpLVariable<T, G>> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<auto f, typename V, typename E>
	len_type data_len(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	len_type data_len(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename V, typename E1, typename E2>
	len_type data_len(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<size_t D>
	len_type data_len(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename V, size_t D, typename E>
	len_type data_len(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinaryAdd<E1, E2> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinarySub<E1, E2> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::property::data_len(E const&).
	template<typename E1, typename E2>
	len_type data_len(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the length of the data in the expression.
	/*!
	 * Return the data_len of the underlying data set.
	 * Chooses the appropriate expr::property::data_len() between two expressions.
	 * This standardizes the choosing procedure for binary expressions.
	 */
	template<typename E1, typename E2>
	len_type data_len(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return std::max(
			data_len(*static_cast<const E1*>(&a)),
			data_len(*static_cast<const E2*>(&b)));
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

	template<typename T, typename G>
	len_type data_len(OpLVariable<T, G> const& e)
	{
		return data_len_data(e.data);
	}
	
	template<typename T, typename... Gs>
	len_type data_len(OpNLVariable<T, Gs...> const& e)
	{
		return data_len_data(std::get<0>(e.datas), symphas::lib::get_tuple_ge<1>(e.datas));
	}

	template<typename Dd, typename V, typename E, typename Sp>
	len_type data_len(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return data_len_data(expr::compound_get::template grid(e));
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

	template<typename V, typename T, typename G>
	len_type data_len(OpExponential<V, OpLVariable<T, G>> const& e)
	{
		return data_len_data(expr::compound_get::template grid(e));
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
		return data_len(expr::compound_get::template expr(e), e.smoother);
	}

	template<size_t D>
	len_type data_len(GaussianSmoothing<D> const& e)
	{
		return data_len_data(e.data);
	}

	template<typename E1, typename E2>
	len_type data_len(OpBinaryAdd<E1, E2> const& e)
	{
		return data_len(e.a, e.b);
	}

	template<typename E1, typename E2>
	len_type data_len(OpBinarySub<E1, E2> const& e)
	{
		return data_len(e.a, e.b);
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
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename E>
	const len_type* data_dimensions(OpExpression<E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename T, typename G>
	const len_type* data_dimensions(OpLVariable<T, G> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename T, typename... Gs>
	const len_type* data_dimensions(OpNLVariable<T, Gs...> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename Dd, typename V, typename E, typename Sp>
	const len_type* data_dimensions(OpFuncDerivative<Dd, V, E, Sp> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpCombination<A1, A2, E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpChain<A1, A2, E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename V, typename E>
	const len_type* data_dimensions(OpExponential<V, E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename V, typename T, typename G>
	const len_type* data_dimensions(OpExponential<V, OpLVariable<T, G>> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<auto f, typename V, typename E>
	const len_type* data_dimensions(OpFuncApply<f, V, E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	const len_type* data_dimensions(OpFunc<V, E, F, Arg0, Args...> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename V, typename E1, typename E2>
	const len_type* data_dimensions(OpFuncConvolution<V, E1, E2> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename V, size_t D, typename E>
	const len_type* data_dimensions(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<size_t D>
	const len_type* data_dimensions(GaussianSmoothing<D> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryAdd<E1, E2> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinarySub<E1, E2> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryMul<E1, E2> const& e);
	//! Specialization based on expr::property::data_dimensions(E const&).
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryDiv<E1, E2> const& e);


	//! Obtain the dimensions of the data in the expression.
	/*!
	 * Obtain the dimensions from either of the two expressions. 
	 * Chooses the appropriate expr::property::data_dimensions() between two expressions.
	 * This standardizes the choosing procedure for binary expressions by
	 * checking if one expression is `nullptr`, and subsequently getting
	 * the dimensions from the other.
	 */
	template<typename E1, typename E2>
	const len_type* data_dimensions(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		auto data_dimensions1 = data_dimensions(*static_cast<const E1*>(&a));
		return (data_dimensions1 != nullptr) ? data_dimensions1 : data_dimensions(*static_cast<const E2*>(&b));
	}

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

	template<typename T, typename G>
	const len_type* data_dimensions(OpLVariable<T, G> const& e)
	{
		return data_dimensions_data(e.data);
	}

	template<typename T, typename... Gs>
	const len_type* data_dimensions(OpNLVariable<T, Gs...> const& e)
	{
		return data_dimensions_data(std::get<0>(e.datas), symphas::lib::get_tuple_ge<1>(e.datas));
	}

	template<typename Dd, typename V, typename E, typename Sp>
	const len_type* data_dimensions(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return data_dimensions_data(expr::compound_get::template grid(e));
	}

	template<typename A1, typename A2, typename E>
	const len_type* data_dimensions(OpCombination<A1, A2, E> const& e)
	{
		return data_dimensions(e.e);
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

	template<typename V, typename T, typename G>
	const len_type* data_dimensions(OpExponential<V, OpLVariable<T, G>> const& e)
	{
		return data_dimensions_data(expr::compound_get::template grid(e));
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
		return data_dimensions(expr::compound_get::template expr(e), e.smoother);
	}

	template<size_t D>
	const len_type* data_dimensions(GaussianSmoothing<D> const& e)
	{
		return data_dimensions_data(e.data);
	}

	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinaryAdd<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
	}

	template<typename E1, typename E2>
	const len_type* data_dimensions(OpBinarySub<E1, E2> const& e)
	{
		return data_dimensions(e.a, e.b);
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





}

