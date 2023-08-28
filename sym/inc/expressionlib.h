
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
#include "gridpair.h"
#include "expressioniterator.h"
#include "indexseqhelpers.h"
#include "symbolicprototypes.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

#define DERIV_MAX_ORDER 12				//!< Maximum order a derivative can be.


namespace expr
{
	template<typename A, typename B>
	constexpr auto numeric_range(A, B);
	template<typename A>
	constexpr auto numeric_range(A);

	template<typename E>
	len_type data_length_for_iterator(OpEvaluable<E> const& e);

	template<typename T>
	struct iterator_policy_expression;

	template<>
	struct iterator_policy_expression<symphas::sequential_iterator_policy>
	{
		template<typename E>
		using type = expression_iterator<E>;

		template<typename E, size_t D>
		static auto begin(OpEvaluable<E> const& e, len_type begin_pos = 0)
		{
			return expression_iterator(*static_cast<E const*>(&e), begin_pos);
		}

		template<typename E, size_t D>
		static auto end(OpEvaluable<E> const& e, len_type end_pos = 0)
		{
			return expression_iterator(*static_cast<E const*>(&e), end_pos);
		}

		template<typename E>
		static auto begin(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e));
		}

		template<typename E>
		static auto end(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e), interval.len);
		}
	};

	template<>
	struct iterator_policy_expression<symphas::selection_iterator_policy>
	{
		template<typename E>
		using type = expression_iterator_selection<E>;

		template<typename E, size_t D>
		static auto begin(OpEvaluable<E> const& e, grid::region_index_list<D> const& iters)
		{
			return expression_iterator_selection(*static_cast<E const*>(&e), iters.iters);
		}

		template<typename E, size_t D>
		static auto end(OpEvaluable<E> const& e, grid::region_index_list<D> const& iters)
		{
			return expression_iterator_selection(*static_cast<E const*>(&e), iters.iters, iters.len);
		}
	};

	template<>
	struct iterator_policy_expression<symphas::region_iterator_policy>
	{
		template<typename E, size_t D>
		using type = expression_iterator_region<E, D>;

		template<typename E, size_t D>
		static auto begin(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
		{
			return expression_iterator_region(*static_cast<E const*>(&e), interval);
		}

		template<typename E, size_t D>
		static auto end(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
		{
			return expression_iterator_region(*static_cast<E const*>(&e), interval, grid::length<D>(interval));
		}

		template<typename E>
		static auto begin(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e));
		}

		template<typename E>
		static auto end(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e), 1);
		}

		template<typename E>
		static auto begin(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e));
		}

		template<typename E>
		static auto end(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e), interval.len);
		}
	};

	template<>
	struct iterator_policy_expression<symphas::group_iterator_policy>
	{
		template<typename E, size_t D>
		using type = expression_iterator_group<E, D>;

		template<typename E, size_t D>
		static auto begin(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
		{
			return expression_iterator_group(*static_cast<E const*>(&e), interval);
		}

		template<typename E, size_t D>
		static auto end(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
		{
			len_type len = grid::length<D>(interval);
			if (len > 0)
			{
				return expression_iterator_group(*static_cast<E const*>(&e), interval, len / (interval[0][1] - interval[0][0]));
			}
			else
			{
				return expression_iterator_group(*static_cast<E const*>(&e), interval, len);
			}
		}

		template<typename E>
		static auto begin(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e));
		}

		template<typename E>
		static auto end(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e), 1);
		}

		template<typename E>
		static auto begin(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e));
		}

		template<typename E>
		static auto end(OpEvaluable<E> const& e, grid::region_size const& interval)
		{
			return expression_iterator(*static_cast<E const*>(&e), interval.len);
		}
	};


	template<typename T>
	using iterator_policy_t = typename iterator_policy_expression<T>::type;
}

template<typename E>
struct OpEvaluable 
{


	//! Return an iterator the beginning of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the beginning of the data, used when evaluating
	 * the expression.
	 */
	expr::expression_iterator<E> begin() const
	{
		return expr::expression_iterator<E>(cast());
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
	expr::expression_iterator<E> end(len_type len) const
	{
		return expr::expression_iterator<E>(cast(), len);
	}


	//! Return an iterator the beginning of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the beginning of the data, used when evaluating
	 * the expression.
	 */
	template<size_t D>
	expr::expression_iterator_region<E, D> begin(grid::region_interval<D> const& interval) const
	{
		return expr::expression_iterator_region<E, D>(cast(), interval);
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
	template<size_t D>
	expr::expression_iterator_region<E, D> end(grid::region_interval<D> const& interval) const
	{
		return expr::expression_iterator_region<E, D>(cast(), interval, grid::length<D>(interval));
	}

	//! Return an iterator the beginning of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the beginning of the data, used when evaluating
	 * the expression.
	 */
	expr::expression_iterator_selection<E> begin(iter_type* iters) const
	{
		return expr::expression_iterator_selection<E>(cast(), iters);
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
	expr::expression_iterator_selection<E> end(iter_type* iters, len_type len) const
	{
		return expr::expression_iterator_selection<E>(cast(), iters, len);
	}


	//! Return an iterator the beginning of the data.
	/*!
	 * For the data related to the expression, return an iterator
	 * representing the beginning of the data, used when evaluating
	 * the expression.
	 */
	template<typename picked_iterator_t, typename arg_t>
	auto begin(picked_iterator_t, arg_t&& data = arg_t{}) const
	{
		return expr::iterator_policy_expression<picked_iterator_t>::begin(cast(), std::forward<arg_t>(data));
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
	template<typename picked_iterator_t, typename arg_t>
	auto end(picked_iterator_t, arg_t&& data = arg_t{}) const
	{
		return expr::iterator_policy_expression<picked_iterator_t>::end(cast(), std::forward<arg_t>(data));
	}

	auto& cast() const
	{
		return *static_cast<E const*>(this);
	}
};

 //! Base expression object which is inherited from with the CRTP technique.
 /*
  * applying Expression Templates to create the expression tree for the
  * evaluation of the equations of motion
  */
template<typename E>
struct OpExpression : OpEvaluable<E>
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
	auto operator[](iter_type n) const
	{
		return cast().eval(n);
	}

	template<typename E0>
	auto operator()(OpExpression<E0> const& e) const
	{
		return cast() * (*static_cast<E0 const*>(&e));
	}

	template<typename E0>
	auto operator()(OpOperator<E0> const& e) const
	{
		return cast() * (*static_cast<E0 const*>(&e));
	}

#ifdef PRINTABLE_EQUATIONS

	//! Print the string representation of this expression to the file.
	/*!
	 * The string representation of this expression is printed to the given
	 * file. The string is assumed to already have enough memory allocated.
	 *
	 * \param out The file to which the expression is printed.
	 */
	size_t print(FILE* out) const
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
	size_t print(char* out) const
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

	size_t print_length() const
	{
		return 0;
	}

#endif

	auto& cast() const
	{
		return *static_cast<E const*>(this);
	}

};

// **************************************************************************************


template<Axis ax, typename G = expr::symbols::Symbol>
struct VectorComponent : G
{
	using G::G;

	VectorComponent(G const& data) : G(data) {}
	VectorComponent(G&& data) noexcept : G(std::forward<G>(data)) {}
	constexpr VectorComponent() : VectorComponent(((G*)0)[0]) {}
};

template<Axis ax, typename G>
struct VectorComponent<ax, symphas::ref<G>> : symphas::ref<G>
{
	using parent_type = symphas::ref<G>;
	using parent_type::parent_type;

	VectorComponent(symphas::ref<G> const& data) : parent_type(data) {}
	VectorComponent(symphas::ref<G>&& data) noexcept : parent_type(std::move(data)) {}
	VectorComponent(G&& data) noexcept : parent_type(std::move(data)) {}
	constexpr VectorComponent() : VectorComponent(((G*)0)[0]) {}
};


template<Axis ax>
struct VectorComponent<ax, void> {};

template<Axis ax, typename G>
struct VectorComponent<ax, const G> : VectorComponent<ax, G>
{
	using parent_type = VectorComponent<ax, G>;
};

template<Axis ax, typename T, size_t D>
struct VectorComponentData
{
	T values;
};

//! Wraps a pointer in order to represent it as a grid.
/*!
 * Wraps a pointer as a grid type, so that the symbolic algebra can interpret
 * the data as having a dimension and the size of the dimensions.
 * 
 * The given pointer is not managed by GridData, it is simply a wrapper.
 */
template<typename T, size_t D>
struct GridData
{
	GridData(T* data, len_type const* dims) : data{ data }, dims{ 0 }, len{ 0 }
	{
		for (iter_type i = 0; i < D; ++i)
		{
			this->dims[i] = dims[i];
		}
		len = grid::length<D>(dims);
	}

	GridData() : data{ nullptr }, dims{ 0 }, len{ 0 } {}

	const T& operator[](iter_type i) const
	{
		return data[i];
	}
	
	T& operator[](iter_type i)
	{
		return data[i];
	}

	const T& operator*() const
	{
		return *data;
	}

	operator const T* () const
	{
		return data;
	}

	operator T* ()
	{
		return data;
	}

protected:

	T* data;

public:

	len_type dims[D];
	len_type len;

};

template<size_t D>
struct GridData<void, D> 
{
	GridData(len_type const* dims) : dims{ 0 }, len{ 0 }
	{
		for (iter_type i = 0; i < D; ++i)
		{
			this->dims[i] = dims[i];
		}
		len = grid::length<D>(dims);
	}

	GridData() : dims{ 0 }, len{ 0 } {}

	len_type dims[D];
	len_type len;
};


template<typename T, size_t D>
struct GridData<expr::symbols::SymbolType<T>, D>
{
	GridData(len_type const* dims) : dims{ 0 }, len{ 0 }
	{
		for (iter_type i = 0; i < D; ++i)
		{
			this->dims[i] = dims[i];
		}
		len = grid::length<D>(dims);
	}

	GridData() : dims{ 0 }, len{ 0 } {}
	
	operator T() const
	{
		return T{};
	}

	len_type dims[D];
	len_type len;
};

template<typename T, size_t N>
struct GridData<MultiBlock<N, T>, N>
{
	GridData(MultiBlock<N, T> const& data = MultiBlock<N, T>(0), len_type const* dims = nullptr) 
		: data{ std::ref(const_cast<MultiBlock<N, T>&>(data)) }, dims{ 0 }, len{ 0 }
	{
		if (dims != nullptr)
		{
			for (iter_type i = 0; i < N; ++i)
			{
				this->dims[i] = dims[i];
			}
			len = grid::length<N>(dims);
		}
	}

	const auto& operator[](iter_type i) const
	{
		return data.get()[i];
	}

	auto& operator[](iter_type i)
	{
		return data.get()[i];
	}

	operator const MultiBlock<N, T>& () const
	{
		return data.get();
	}

	operator MultiBlock<N, T>& ()
	{
		return data.get();
	}

protected:

	symphas::ref<MultiBlock<N, T>> data;

public:

	len_type dims[N];
	len_type len;

};



// **************************************************************************************


template<>
inline constexpr bool is_simple_data<expr::symbols::Symbol> = true;

namespace expr
{
	template<Axis ax, size_t D, typename T>
	auto as_component_data(T* data)
	{
		return VectorComponentData<ax, T*, D>{ data };
	}

	template<Axis ax, size_t N, typename T>
	VectorComponentData<ax, T*, N> resolve_axis_component(MultiBlock<N, T> const& data)
	{
		return { data.values[symphas::axis_to_index(ax)] };
	}

	template<Axis ax, size_t N, typename T>
	VectorComponentData<ax, T*, N> resolve_axis_component(MultiBlock<N, T>& data)
	{
		return { data.values[symphas::axis_to_index(ax)] };
	}

	template<Axis ax, typename T, size_t D>
	VectorComponentData<ax, T, D> resolve_axis_component(any_vector_t<T, D> const& data)
	{
		return { data[symphas::axis_to_index(ax)] };
	}

	template<Axis ax, typename T, size_t D>
	VectorComponentData<ax, T, D> resolve_axis_component(any_vector_t<T, D>& data)
	{
		return { data[symphas::axis_to_index(ax)] };
	}


	//! Constructs a constant of the given value.
	/*!
	 * Constructs a constant of the given value.
	 *
	 * \param v The value to give to the literal.
	 */
	template<typename T>
	constexpr decltype(auto) make_literal(T const& v);

	//! Constructs a fraction with the given values.
	/*!
	 * Constructs a compile time constant fraction with positive whole numbers.
	 *
	 * \tparam N The value of the numerator, a positive whole number.
	 * \tparam D The value of the denominator, a positive whole number.
	 */
	template<size_t N, size_t D>
	constexpr auto make_fraction();

	template<size_t N, size_t D>
	constexpr auto frac = make_fraction<N, D>();

	//! Constructs an integer as a fraction type.
	/*!
	 * Constructs a compile time constant integer as a fraction type. This is a wrapper
	 * for the expr::make_fraction() function, particularly useful when the given integer
	 * should be negative.
	 */
	template<int I>
	constexpr auto make_integer();

	template<int I>
	constexpr auto val = make_integer<I>();

	//! Construct a tensor entry, which acts as a coefficient.
	/*!
	 * The first half of the given values represent the position of the
	 * value in the tensor, and the last half represent the rank of the tensor.
	 */
	template<size_t N0, size_t N1, size_t... Ns, typename T>
	constexpr auto make_tensor(T const& v);

	//! Construct a tensor entry with an identity element.
	/*!
	 * The first half of the given values represent the position of the
	 * value in the tensor, and the last half represent the rank of the tensor.
	 */
	template<size_t N0, size_t N1, size_t... Ns>
	constexpr auto make_tensor();

	//! Construct an entry of a column vector, behaves like a coefficient.
	/*!
	 * A vector is a specialization of a tensor. It is a column vector. Taking
	 * the transpose of a column vector converts it into a row vector.
	 */
	template<size_t I, size_t N, typename T>
	constexpr auto make_column_vector(T&& v);

	template<size_t I, size_t N>
	constexpr auto make_column_vector();

	//! Construct selectable coefficients using a coefficient array.
	/*!
	 * For coefficients which should be selected from a list, construct a coefficient which
	 * can select them based on their index. 
	 */
	template<typename T>
	auto make_coeff(T* data, len_type len, len_type stride = 1);

	template<typename T, typename I, typename T0, typename I0>
	auto init_coeff_from(OpCoeff<T0, I0> const& coeff);

	//! Construct an entry of a row vector, behaves like a coefficient.
	/*!
	 * A vector is a specialization of a tensor. It is a column vector. Taking
	 * the transpose of a row vector converts it into a column vector. A row
	 * vector is treated as a matrix of one row.
	 */
	template<size_t I, size_t N, typename T>
	constexpr auto make_row_vector(T&& v);

	template<size_t I, size_t N>
	constexpr auto make_row_vector();

	template<size_t R, size_t R0 = 0>
	auto make_filled_column_vector();

	template<size_t R, size_t R0 = 0>
	auto make_filled_row_vector();

	//! Constructs the expression representing an expression to a power.
	/*!
	 * Directly constructs the exponent expression of an
	 * expression without applying any rules.
	 *
	 * \param value The coefficient.
	 * \param e The expression.
	 */
	template<expr::exp_key_t X, typename V, typename E>
	auto make_pow(V const& value, OpExpression<E> const& e);

	//! Constructs the expression representing an expression to a power.
	/*!
	 * Directly constructs the exponent expression of an
	 * expression without applying any rules.
	 *
	 * \param value The coefficient.
	 * \param e The expression.
	 */
	template<expr::exp_key_t X, typename V, typename V0, typename... Gs, expr::exp_key_t... Xs>
	auto make_pow(V const& value, OpTerms<V0, Term<Gs, Xs>...> const& e);

	//! Constructs the expression representing an expression to a power.
	/*!
	 * Directly constructs the exponent expression of an
	 * expression without applying any rules.
	 *
	 * \param value The coefficient.
	 * \param e The expression.
	 */
	template<expr::exp_key_t X, typename E>
	auto make_pow(E const& e);


	//! Returns the given data in a wrapper imitating a Grid with specified dimensions.
	/*!
	 * The given data is wrapped in an object which imitates a Grid in the symbolic
	 * algebra functionality, so that dimensions can be interpreted. This is used when
	 * the given data is a non-grid type, such as a raw pointer.
	 * 
	 * \param values The data.
	 * \param dims The dimensions of the data.
	 */
	template<size_t D, typename T>
	auto as_grid_data(T* values, const len_type* dims)
	{
		return GridData<T, D>(values, dims);
	}

	//! Returns the given data in a wrapper imitating a Grid with specified dimensions.
	/*!
	 * The given data is wrapped in an object which imitates a Grid in the symbolic
	 * algebra functionality, so that dimensions can be interpreted. This is used when
	 * the given data is a non-grid type, such as a raw pointer.
	 *
	 * \param values The data.
	 * \param dims The dimensions of the data.
	 */
	template<size_t N, typename T>
	auto as_grid_data(MultiBlock<N, T> const& data, const len_type* dims)
	{
		return GridData<MultiBlock<N, T>, N>(data, dims);
	}

	//! Returns the given data in a wrapper imitating a Grid with specified dimensions.
	/*!
	 * The given data is wrapped in an object which imitates a Grid in the symbolic
	 * algebra functionality, so that dimensions can be interpreted. This is used when
	 * the given data is a non-grid type, such as a raw pointer.
	 *
	 * \param values The data.
	 * \param dims The dimensions of the data.
	 */
	template<size_t D, typename T>
	auto as_grid_data(T const& data, const len_type* dims)
	{
		return data;
	}

	//! The unit vector, which can be defined with one or two (for 3D) angles.
	/*!
	 * The correct unit vector will be chosen according to the dimension. In two dimensions,
	 * only the first direction will be used, and in three, the second direction will be used.
	 * If only one direction is provided, then in 3D, the second direction is assumed to be 0.
	 */
	template<size_t D, typename T0, typename T1>
	auto make_unit_vector(T0 const& direction0, T1 const& direction1);

	//! The unit vector, which can be defined with one or two (for 3D) angles.
	/*!
	 * The correct unit vector will be chosen according to the dimension. In 3D, the second 
	 * direction is assumed to be 0.
	 */
	template<size_t D, typename T0>
	auto make_unit_vector(T0 const& direction0);

	template<typename T, size_t... Ns>
	auto inverse(OpTensor<T, Ns...> const& tensor) = delete;

	//! Apply an inverse to a scalar value.
	inline auto inverse(scalar_t e);
	//! Apply an inverse to an integer value.
	inline auto inverse(int e);
	//! Apply an inverse to a complex value.
	inline auto inverse(complex_t const& e);
	//! Apply an inverse to an expression.
	template<typename E>
	auto inverse(OpExpression<E> const& e);
	//! Apply an inverse to an expression.
	template<typename A, typename B>
	auto inverse(OpBinaryDiv<A, B> const& e);
	//! Apply an inverse to an expression literal.
	template<typename T>
	auto inverse(OpLiteral<T> const& e);
	//! Apply an inverse to an expression literal.
	template<size_t N, size_t D>
	constexpr auto inverse(OpFractionLiteral<N, D>);
	//! Apply an inverse to an expression literal.
	template<size_t N, size_t D>
	constexpr auto inverse(OpNegFractionLiteral<N, D>);
	//! Apply an inverse to an expression literal.
	constexpr inline auto inverse(OpIdentity);
	//! Apply an inverse to an expression literal.
	constexpr inline auto inverse(OpNegIdentity);
	//! Apply an inverse to an expression.
	template<typename V, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	auto inverse(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e);
	//! Apply an inverse to an expression.
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	auto inverse(OpTerms<V, Term<Gs, Xs>...> const& e);
	//! Apply an inverse to an expression.
	template<typename V, typename E>
	auto inverse(OpExponential<V, E> const& e);
	//! Apply an inverse to an expression.
	template<expr::exp_key_t X, typename V, typename E>
	auto inverse(OpPow<X, V, E> const& e);

}

namespace symphas::internal 
{


	template<int N0, typename V>
	struct search_index_in_v
	{
		static const bool value = 0;
	};

	template<int N0, int N1, int P, typename... Is>
	struct search_index_in_v<N0, expr::symbols::v_id_type<expr::symbols::i_<N1, P>, Is...>>
	{
		static const bool value = (N0 == N1) ? P : search_index_in_v<N0, expr::symbols::v_id_type<Is...>>::value;
	};

	template<typename I, typename T>
	constexpr bool has_matching_i = false;

	template<int N0, int P0, int... Ns, int... Ps>
	constexpr bool has_matching_i<expr::symbols::i_<N0, P0>, expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>
		= (symphas::lib::index_of_value<int, N0, Ns...> >= 0);
	
	template<int N0, int P0, int... Ns, int... Ps, size_t D>
	constexpr bool has_matching_i<expr::symbols::i_<N0, P0>, GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>>
		= has_matching_i<expr::symbols::i_<N0, P0>, expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>;

	template<bool flag, typename T0>
	struct test_select_v_i_impl;


	template<typename T0>
	struct test_select_v_i_impl<true, T0>
	{
		using type = symphas::lib::types_list<T0>;
	};

	template<typename T0>
	struct test_select_v_i_impl<false, T0>
	{
		using type = symphas::lib::types_list<>;
	};

	template<typename I, typename T0, typename List>
	struct select_v_i_impl;

	template<int N0, int P0, typename... Vs>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int N0, int P0, typename... Vs, typename T0, typename... Ts>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T0, Ts...>>
	{
		using type = typename select_v_i_impl<expr::symbols::i_<N0, P0>, 
			symphas::lib::expand_types_list<symphas::lib::types_list<Vs...>, typename test_select_v_i_impl<has_matching_i<expr::symbols::i_<N0, P0>, T0>, T0>::type>,
			symphas::lib::types_list<Ts...>>::type;
	};

	//template<int N0, int P0, typename... Vs, int N1, int P1, int... Ns, int... Ps, typename... Rest>
	//struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N1, P1>, expr::symbols::i_<Ns, Ps>...>, Rest...>>
	//{
	//	//static const bool flag = symphas::lib::index_of_value<int, N0, Ns...> >= 0;
	//	//
	//	//using type = typename select_v_i_impl<
	//	//	expr::symbols::i_<N0, P0>,
	//	//	std::conditional_t<flag,
	//	//		symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>,
	//	//		symphas::lib::types_list<Vs...>>,
	//	//	symphas::lib::types_list<Rest...>>::type;
	//	//using type = typename select_v_i_impl < expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<N1, P1>, expr::symbols::i_<Ns, Ps>...>>, symphas::lib::types_list<Rest...>>::type;
	//};

	//template<int N0, int P0, typename... Vs, int... Ns, int... Ps, size_t D, typename... Rest>
	//struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>, Rest...>>
	//{
	//	using type = typename select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>::type;
	//};

	//template<int N0, int P0, typename... Vs, typename T, typename... Rest>
	//struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	//{
	//	using type = typename select_v_i_impl<
	//		expr::symbols::i_<N0, P0>,
	//		symphas::lib::types_list<Vs...>,
	//		symphas::lib::types_list<Rest...>>::type;
	//};

	//template<int N0, int P0, typename... Vs>
	//struct select_v_i_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>>
	//{
	//	using type = typename select_v_i_<
	//		expr::symbols::i_<N0, P0>,
	//		symphas::lib::types_list<>,
	//		symphas::lib::types_list<Vs...>>::type;
	//};

	template<typename I, typename List>
	using select_v_i_ = typename select_v_i_impl<I, symphas::lib::types_list<>, List>::type;

    template<typename T>
    constexpr bool is_v_type = false;

    template<int... Ns, int... Ps>
    constexpr bool is_v_type<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>> = true;
    template<int... Ns, int... Ps, size_t D>
    constexpr bool is_v_type<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>> = true;


	template<size_t flag, typename Is, typename T0, typename List>
	struct select_v_nested_impl;

	template<int... N0s, int... P0s, typename... Vs>
	struct select_v_nested_impl<2, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int... N0s, int... P0s, typename... Vs, int... Ns, int... Ps, typename... Rest>
	struct select_v_nested_impl<1, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>
	{
		static const bool flag = symphas::lib::filter_seq_t<std::integer_sequence<int, Ns...>, std::integer_sequence<int, N0s...>>::size() == 0;

		using type = typename select_v_nested_impl<2,
			symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>,
			std::conditional_t<false,
				symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>,
				symphas::lib::types_list<Vs...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int... N0s, int... P0s, typename... Vs, int... Ns, int... Ps, size_t D, typename... Rest>
	struct select_v_nested_impl<1, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>, Rest...>>
	{
		using type = typename select_v_nested_impl<1, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>::type;
	};

	template<int... N0s, int... P0s, typename... Vs, typename T, typename... Rest>
	struct select_v_nested_impl<2, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_nested_impl<
            size_t(is_v_type<T>),
			symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>,
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<T, Rest...>>::type;
	};

	template<int... N0s, int... P0s, typename... Vs, typename T, typename... Rest>
	struct select_v_nested_impl<0, symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_nested_impl<2,
			symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>,
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename Is, typename List>
	using select_v_nested_ = typename select_v_nested_impl<2, Is, symphas::lib::types_list<>, List>::type;


	template<size_t flag, typename... Ts>
	struct select_v_impl;
	
	template<typename... Vs, typename I0, typename... Is, typename... Rest>
	struct select_v_impl<1, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<I0, Is...>, Rest...>>
	{
		using type = typename select_v_impl<2,
			symphas::lib::types_list<Vs..., expr::symbols::v_id_type<I0, Is...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Vs, typename I0, typename... Is, size_t D, typename... Rest>
	struct select_v_impl<1, symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<I0, Is...>, D>, Rest...>>
	{
		using type = typename select_v_impl<1, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<I0, Is...>, Rest...>>::type;
	};

	template<typename... Vs, typename T, typename... Rest>
	struct select_v_impl<2, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_impl<
            size_t(is_v_type<T>),
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<T, Rest...>>::type;
	};

	template<typename... Vs, typename T, typename... Rest>
	struct select_v_impl<0, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_impl<2, 
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};
    
	template<typename... Vs>
	struct select_v_impl<2, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<typename List>
	using select_v_ = typename select_v_impl<2, symphas::lib::types_list<>, List>::type;



	template<typename... Ts>
	struct select_all_i_impl_;

	template<int N, typename... Is>
	struct select_all_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Is...>;
	};

	template<int N, typename... Is, int P0, typename... Rest>
	struct select_all_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N, P0>>, Rest...>>
	{
		using type = typename select_all_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is..., expr::symbols::i_<N, P0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, int N0, int P0, typename... Rest>
	struct select_all_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::i_<N0, P0>, Rest...>>
	{
		using type = typename select_all_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, int P0, typename... Rest>
	struct select_all_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::i_<N, P0>, Rest...>>
	{
		using type = typename select_all_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is..., expr::symbols::i_<N, P0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, typename Other, typename... Rest>
	struct select_all_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<Other, Rest...>>
	{
		using type = typename select_all_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename... Vs>
	struct select_all_i_impl_<expr::symbols::i_<N, P>, symphas::lib::types_list<Vs...>>
	{
		using type = typename select_all_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type;
	};

	template<int N0, int P0, int... Ns, int... Ps, typename... Vs>
	struct select_all_i_impl_<symphas::lib::types_list<
			expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>, 
		symphas::lib::types_list<Vs...>>
	{
		using type = symphas::lib::expand_types_list<
			typename select_all_i_impl_<expr::symbols::i_<N0, 0>, symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type,
			typename select_all_i_impl_<expr::symbols::i_<Ns, 0>, symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type...
			>;
	};

	template<typename I, typename T>
	using select_all_i_ = typename select_all_i_impl_<I, T>::type;




	template<typename... Ts>
	struct select_only_i_impl_;

	template<int N, typename... Is>
	struct select_only_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Is...>;
	};

	template<int N, typename... Is, int P0, typename... Rest>
	struct select_only_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::i_<N, P0>, Rest...>>
	{
		using type = typename select_only_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is..., expr::symbols::i_<N, P0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, typename Other, typename... Rest>
	struct select_only_i_impl_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<Other, Rest...>>
	{
		using type = typename select_only_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename... Vs>
	struct select_only_i_impl_<expr::symbols::i_<N, P>, symphas::lib::types_list<Vs...>>
	{
		using type = typename select_only_i_impl_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type;
	};

	template<int N0, int P0, int... Ns, int... Ps, typename... Vs>
	struct select_only_i_impl_<symphas::lib::types_list<
		expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>,
		symphas::lib::types_list<Vs...>>
	{
		using type = symphas::lib::expand_types_list<
			typename select_only_i_impl_<expr::symbols::i_<N0, 0>, symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type,
			typename select_only_i_impl_<expr::symbols::i_<Ns, 0>, symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type...
		>;
	};

	template<typename I, typename T>
	using select_only_i_ = typename select_only_i_impl_<I, T>::type;


	template<typename... Ts>
	struct select_unique_i_impl_;

	template<int... Ns>
	struct select_unique_i_impl_<std::integer_sequence<int, Ns...>>
	{
		using type = symphas::lib::types_list<expr::symbols::i_<Ns, 0>...>;
	};

	template<int... Ns>
	struct select_unique_i_impl_<symphas::lib::types_list<expr::symbols::i_<Ns, 0>...>, symphas::lib::types_list<>>
	{
		using type = typename select_unique_i_impl_<symphas::lib::sorted_seq<std::integer_sequence<int, Ns...>>>::type;
	};

	template<int N, int P, int M0, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::index_neq_N<expr::symbols::i_<N, P>, M0>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename B, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::index_neq<expr::symbols::i_<N, P>, B>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<B, Rest...>>::type;
	};

	template<int N, int P, int M0, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<N, P>, M0>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename B, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::index_eq<expr::symbols::i_<N, P>, B>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<B, Rest...>>::type;
	};

	template<typename A, typename B, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<symphas::lib::types_list<A, B>, Rest...>>
	{
		using type = symphas::lib::expand_types_list<
			typename select_unique_i_impl_<symphas::lib::types_list<>, symphas::lib::types_list<A>>::type,
			typename select_unique_i_impl_<
			symphas::lib::types_list<Is...>, symphas::lib::types_list<B, Rest...>>::type>;
	};

	template<int N, int P, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::i_<N, P>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename... Is, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N, P>>, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is..., expr::symbols::i_<N, 0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Is, typename Other, typename... Rest>
	struct select_unique_i_impl_<symphas::lib::types_list<Is...>, symphas::lib::types_list<Other, Rest...>>
	{
		using type = typename select_unique_i_impl_<
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Es>
	struct select_unique_i_impl_<symphas::lib::types_list<Es...>>
	{
		using type = typename symphas::lib::combine_types_unique<typename select_unique_i_impl_<
			symphas::lib::types_list<>, symphas::lib::types_list<Es...>>::type>::type;
	};

	template<typename T>
	using select_unique_i_ = typename select_unique_i_impl_<T>::type;


	template<typename... Is>
	struct filter_for_unique_vs
	{
		using type = symphas::lib::types_list<>;
	};

	template<int I0, int P0, int... I0s, int... P0s>
	struct filter_for_unique_vs<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>
	{
		using type = typename symphas::lib::combine_types_unique<
			expr::symbols::i_<I0, 0>,
			expr::symbols::i_<I0s, 0>...>::type;
	};

	template<int I0, int P0, int... I0s, int... P0s>
	struct filter_for_unique_vs<expr::symbols::v_<expr::symbols::i_<I0, P0>>, expr::symbols::v_<expr::symbols::i_<I0s, P0s>>...>
	{
		using type = typename filter_for_unique_vs<
			expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
			expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>::type;
	};


}



namespace expr
{

	//! Type trait used to identify the original type of a data.
	/*!
	 * Applied to the data object of a variable expression in order to obtain
	 * the original data type after certain wrapping types are taken away.
	 * Applied to any other data type, it will simply directly return it.
	 * 
	 * This is different from ::base_data_type because this will return the data type
	 * that the variable represents, whereas ::base_data_type has fewer specializations,
	 * and returns only the variable itself.
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

	//! Specialization for Variable type. See expr::original_data_type.
	template<size_t Z, typename T>
	struct original_data_type<Variable<Z, T*>>
	{
		using type = typename original_data_type<T>::type;
	};

	//! Specialization for Variable type. See expr::original_data_type.
	template<typename T>
	struct original_data_type<DynamicVariable<T>>
	{
		using type = typename original_data_type<T>::type;
	};

	//! Specialization for Variable type. See expr::original_data_type.
	template<typename T>
	struct original_data_type<DynamicVariable<NamedData<T*>>>
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

	//! Specialization for Variable type. See expr::original_data_type.
	template<Axis ax, typename G>
	struct original_data_type<VectorComponent<ax, G>>
	{
		using type = VectorComponent<ax, typename original_data_type<G>::type>;
	};

	//! Specialization for Variable type. See expr::original_data_type.
	template<typename G>
	using original_data_t = typename original_data_type<G>::type;
}



// *******************************************************************************

namespace symphas::internal
{

	/* testing to check if the type has the "value" attribute, which
	* is the coefficient to an expression.
	 */

	template<typename E>
	static char test_coeff_attribute(decltype(E::value)*);
	template<typename E>
	static long test_coeff_attribute(...);

	template<typename E>
	struct coeff_neg_trait
	{
		using type = long;
	};

	template<>
	struct coeff_neg_trait<OpNegIdentity>
	{
		using type = char;
	};

	template<size_t N, size_t D>
	struct coeff_neg_trait<OpNegFractionLiteral<N, D>>
	{
		using type = char;
	};

	template<>
	struct coeff_neg_trait<scalar_t>
	{
		using type = char;
	};

	template<>
	struct coeff_neg_trait<complex_t>
	{
		using type = char;
	};

	template<typename T>
	struct coeff_neg_trait<OpLiteral<T>>
	{
		using type = char;
	};

	template<typename T, typename I>
	struct coeff_neg_trait<OpCoeff<T, I>>
	{
		using type = char;
	};


	template <typename E> static auto test_coeff_neg(decltype(E::value)*)
		-> typename coeff_neg_trait<decltype(E::value)>::type;
	template <typename E> static long test_coeff_neg(...);


	template<typename E>
	struct test_is_op_fraction
	{
		static const bool value = false;
	};

	template<size_t N, size_t D>
	struct test_is_op_fraction<OpFractionLiteral<N, D>>
	{
		static const bool value = true;
	};

	template<size_t N, size_t D>
	struct test_is_op_fraction<OpNegFractionLiteral<N, D>>
	{
		static const bool value = true;
	};


	template<typename E>
	struct test_is_symbol
	{
		static const bool value = std::is_convertible_v<E, expr::symbols::Symbol>;
	};

	template<typename T, size_t D>
	struct test_is_symbol<GridData<T, D>>
	{
		static const bool value = test_is_symbol<T>::value;
	};

	template<typename T, size_t D>
	struct test_is_symbol<GridSymbol<T, D>>
	{
		static const bool value = true;
	};

	template<typename E>
	struct test_is_op_tensor
	{
		static const bool value = false;
	};

	template<typename T, size_t... Ns>
	struct test_is_op_tensor<OpTensor<T, Ns...>>
	{
		static const bool value = true;
	};

	template<typename... Ts>
	struct test_is_op_tensor<OpAdd<Ts...>>
	{
		static const bool value = (test_is_op_tensor<Ts>::value && ...);
	};

}


namespace expr
{

	//! Return type of evaluating an expression identity.
	/*!
		* Returns the type of the multiplicative identity. The same
		* type would be returned for the negative of the multiplicative identity.
		*/
		//using identity_eval_t = decltype(std::declval<OpIdentity>().eval());

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
	 //template<typename E1, typename E2>
	 //constexpr bool identity_comparator_enabled =
	 //	((std::is_same<E1, OpIdentity>::value && std::is_convertible<E2, identity_eval_t>::value) ||
	 //		(std::is_same<E2, OpIdentity>::value && std::is_convertible<E1, identity_eval_t>::value) ||
	 //		(std::is_same<E1, OpNegIdentity>::value && std::is_convertible<E2, identity_eval_t>::value) ||
	 //		(std::is_same<E2, OpNegIdentity>::value && std::is_convertible<E1, identity_eval_t>::value));
	 
	 
	 //! True if the given type is a fraction.
	 /*!
	 * This value is true if the given type is any of OpFractionLiteral, or its negative.
	 *
	 * \tparam E The given type that is checked.
	 */
	template<typename E>
	constexpr bool is_fraction = symphas::internal::test_is_op_fraction<E>::value;

	template<typename E>
	constexpr bool is_literal = false;

	template<typename T>
	constexpr bool is_literal<OpLiteral<T>> = true;

	template<typename E>
	constexpr bool is_arr_coeff = false;

	template<typename T, typename I>
	constexpr bool is_arr_coeff<OpCoeff<T, I>> = true;

	//! True if the given type has the base type expr::symbols::Symbol.
	template<typename E>
	constexpr bool is_symbol = symphas::internal::test_is_symbol<E>::value;

	//! True if the given type is a tensor element.
	/*!
		* This value is true if the given type is any of OpTensor.
		*
		* \tparam E The given type that is checked.
		*/
	template<typename E>
	constexpr bool is_tensor = symphas::internal::test_is_op_tensor<E>::value;

	//template<typename... Es>
	//constexpr bool is_tensor<OpAdd<Es...>> = (is_tensor<Es> && ...);

	template<typename E>
	constexpr bool is_coeff =
		(is_fraction<E> || is_identity<E> || is_tensor<E> || is_literal<E> || is_arr_coeff<E>) && !std::is_same<E, OpVoid>::value;

	template<typename... Es>
	constexpr bool is_coeff<OpAdd<Es...>> = (is_coeff<Es> && ...);

	template<typename E>
	constexpr bool is_add = false;
	template<typename... Es>
	constexpr bool is_add<OpAdd<Es...>> = true;
	
    
    template<typename T>
    constexpr bool is_id_variable = false;

    template<size_t Z, typename G>
    constexpr bool is_id_variable<Variable<Z, G>> = true;

    template<typename G>
    constexpr bool is_id_variable<DynamicVariable<G>> = true;

	using symphas::internal::test_coeff_attribute;
	using symphas::internal::test_coeff_neg;

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

	template<typename V, typename... Gs>
	constexpr bool has_coeff<OpTerms<V, Gs...>> = expr::is_coeff<V> || is_simple_data<V>;

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

	template<typename V, typename... Gs>
	constexpr bool has_nmi_coeff<OpTerms<V, Gs...>> = sizeof(typename symphas::internal::coeff_neg_trait<V>::type) == sizeof(char);

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
	constexpr bool has_pmi_coeff = has_coeff<E> && !has_nmi_coeff<E>;


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

	template<typename T>
	struct clear_named_data;

	//! Determines the dimension of the data in the expression.
	/*!
	 * Type trait that will prune the grid dimension from an expression based
	 * on its type. If there is no dimension, the evaluated dimension is 0.
	 */
	template<typename E>
	struct grid_dim
	{
	protected:

		template<typename E0>
		static constexpr auto _infer_dimension(OpExpression<E0>)
		{
			if constexpr (is_coeff<E0>)
			{
				return std::index_sequence<0>{};
			}
			else
			{
				return std::index_sequence<grid_dim<E0>::dimension>{};
			}
		}

		template<typename E0>
		static constexpr auto _infer_dimension(OpOperator<E0>)
		{
			return std::index_sequence<0>{};
		}

		template<typename E0, typename E1>
		static constexpr auto _infer_dimension(OpOperatorChain<E0, E1>)
		{
			return std::index_sequence<grid_dim<E0>::dimension, grid_dim<E1>::dimension>{};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(any_vector_t<T, D>)
		{
			return {};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(Grid<T, D>)
		{
			return {};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(BoundaryGrid<T, D>)
		{
			return {};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(RegionalGrid<T, D>)
		{
			return {};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(MultiBlock<D, T>)
		{
			return {};
		}

		template<typename T, size_t D>
		static constexpr std::index_sequence<D> _infer_dimension(GridData<T, D>)
		{
			return {};
		}

		template<typename T0, typename... Ts>
		static constexpr auto _infer_dimension(Substitution<T0, Ts...>)
		{
			return std::index_sequence<grid_dim<T0>::value>{};
		}

		template<typename T>
		static constexpr auto _infer_dimension(SymbolicData<T>)
		{
			return std::index_sequence<grid_dim<T>::value>{};
		}

		template<Axis ax, typename G>
		static constexpr auto _infer_dimension(VectorComponent<ax, G>)
		{
			return std::index_sequence<grid_dim<G>::value>{};
		}

		static constexpr std::index_sequence<0> _infer_dimension(...)
		{
			return {};
		}

		static constexpr size_t infer_dimension()
		{
			return symphas::lib::seq_index_value<-1, symphas::lib::sorted_seq<decltype(_infer_dimension(std::declval<E>()))>>::value;
		}

	public:

		static const size_t dimension = infer_dimension();
		static const size_t value = dimension;
	};

	//! Constructs the list of all data types.
	/*!
	 * Constructs a types list where the template arguments are (uniquely
	 * listed) template arguments of the OpTerms from the given
	 * expression.
	 */
	template<typename E>
	struct op_types
	{
		using type = symphas::lib::types_list<>;
	};

	template<typename E>
	using op_types_t = typename op_types<E>::type;

	//! Determines the underlying data in the expression.
	/*!
	 * Type trait returns as its type a list of types of data that is found in the
	 * expression.
	 */
	template<typename E>
	struct term_types;
	template<typename E>
	using term_types_t = typename term_types<E>::type;

	//! Determines the type that the given expression evaluates to.
	/*!
	 * Determines the type that the given expression evaluates to.
	 */
	template<typename E>
	struct eval_type;
	template<typename E>
	using eval_type_t = typename eval_type<E>::type;


	//! Determines the grid which would enclose the data of the expression.
	/*!
	 * Determines the evaluation type of the expression, and then builds the
	 * type based on the type of the grid that first appears in the given 
	 * expression type. This is distinct from expr::eval_type, in that
	 * the resulting type here is the full grid.
	 */
	template<typename E>
	struct storage_type;
	template<typename E>
	using storage_type_t = typename storage_type<E>::type;

	template<typename G>
	struct parent_storage_type;
	template<typename E>
	using parent_storage_t = typename parent_storage_type<E>::type;

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
		template<typename T, typename = int>
		struct HasAx : std::false_type {};

		template<typename T>
		struct HasAx<T, decltype((void)T::ax, 0)> : std::true_type { };

		template<typename T, typename = int>
		struct GetAx {};

		template<typename T>
		struct GetAx<T, decltype((void)T::ax, 0)> 
		{ 
			static const Axis ax = T::ax;
		};

	public:

		static const bool value = HasAx<E>::value;
	};

	//! Provides the value of the derivative order of the given operator or expression.
	template<typename E>
	struct derivative_order
	{
		static const size_t value = 0;
	};


	template<size_t O, typename V, typename Sp>
	struct derivative_order<OpOperatorDerivative<O, V, Sp>>
	{
		const static size_t value = O;
	};

	template<Axis ax, size_t O, typename V, typename Sp>
	struct derivative_order<OpOperatorDirectionalDerivative<ax, O, V, Sp>>
	{
		const static size_t value = O;
	};

	template<typename V, typename Sp, size_t... Os>
	struct derivative_order<OpOperatorMixedDerivative<V, Sp, Os...>>
	{
		const static size_t value = (Os + ...);
	};

	template<typename Dd, typename V, typename Sp, typename E>
	struct derivative_order<OpDerivative<Dd, V, Sp, E>>
	{
		const static size_t value = OpDerivative<Dd, V, Sp, E>::order;
	};

	template<typename E1, typename E2>
	struct derivative_order<OpOperatorChain<E1, E2>>
	{
		const static size_t value = derivative_order<E1>::value + derivative_order<E2>::value;
	};

	template<typename E1, typename E2>
	struct derivative_order<OpOperatorCombination<E1, E2>>
	{
		const static size_t value = fixed_min<derivative_order<E1>::value, derivative_order<E2>::value>;
	};

	template<typename E1, typename E2, typename E>
	struct derivative_order<OpChain<E1, E2, E>>
	{
		const static size_t value = derivative_order<OpOperatorChain<E1, E2>>::value;
	};

	template<typename E1, typename E2, typename E>
	struct derivative_order<OpCombination<E1, E2, E>>
	{
		const static size_t value = derivative_order<OpOperatorCombination<E1, E2>>::value;
	};

	template<typename A, typename B>
	struct derivative_order<OpBinaryMul<A, B>>
	{
		const static size_t value = derivative_order<A>::value + derivative_order<B>::value;
	};

	template<typename E>
	struct derivative_order<OpOptimized<E>>
	{
		const static size_t value = derivative_order<E>::order;
	};

	template<typename... Es>
	struct derivative_order<OpAdd<Es...>>
	{
		const static size_t value = fixed_max<derivative_order<Es>::value...>;
	};



	//! Combines type traits is_linear and is_combination.
	template<typename E>
	struct expression_predicate;


	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs>
	const auto& get(OpTerms<V, Term<Gs, Xs>...> const& e);
	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs>
	auto& get(OpTerms<V, Term<Gs, Xs>...>& e);
	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs>
	const auto& get(OpTermsList<V, Term<Gs, Xs>...> const& e);
	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs>
	auto& get(OpTermsList<V, Term<Gs, Xs>...>& e);

	template<size_t N, typename... Es>
	const auto& get(OpAdd<Es...> const& e);
	template<size_t N, typename... Es>
	auto& get(OpAdd<Es...>& e);
	template<size_t N, typename... Es>
	const auto& get(OpAddList<Es...> const& e);
	template<size_t N, typename... Es>
	auto& get(OpAddList<Es...>& e);

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
	template<typename... Es, size_t... Is>
	auto reevaluate(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return (expr::get<Is>(e) + ...);
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename... Es>
	auto reevaluate(OpAdd<Es...> const& e)
	{
		return reevaluate(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinaryMul<E1, E2> const& e)
	{
		return e.a * e.b;
	}

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<expr::exp_key_t X, typename V, typename E>
	auto reevaluate(OpPow<X, V, E> const& e);

	//! Specialization based on reevaluate(OpExpression<E> const&).
	template<typename E1, typename E2>
	auto reevaluate(OpBinaryDiv<E1, E2> const& e)
	{
		return e.a / e.b;
	}



}




namespace expr
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



template<typename T>
struct expr::clear_named_data
{
	using type = T;
	static const bool value = false;
};

template<typename T>
struct expr::clear_named_data<Term<T>>
{
	using type = Term<typename clear_named_data<T>::type>;
	static const bool value = clear_named_data<T>::value;
};

template<typename T>
struct expr::clear_named_data<DynamicVariable<T>>
{
	using type = DynamicVariable<typename clear_named_data<T>::type>;
	static const bool value = clear_named_data<T>::value;
};

template<size_t Z, typename T>
struct expr::clear_named_data<Variable<Z, T>>
{
	using type = Variable<Z, typename clear_named_data<T>::type>;
	static const bool value = clear_named_data<T>::value;
};

template<Axis ax, typename T>
struct expr::clear_named_data<VectorComponent<ax, T>>
{
	using type = VectorComponent<ax, typename clear_named_data<T>::type>;
	static const bool value = clear_named_data<T>::value;
};

template<typename T>
struct expr::clear_named_data<NamedData<T>>
{
	using type = typename clear_named_data<T>::type;
	static const bool value = true;
};


//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<const E>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<E&>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<OpExpression<E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<>
struct expr::grid_dim<OpVoid>
{
	static const size_t dimension = 0;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A>
struct expr::grid_dim<symphas::ref<A>>
{
	static const size_t dimension = expr::grid_dim<A>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename T, size_t D>
struct expr::grid_dim<GridData<T, D>>
{
	static const size_t dimension = D;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G>
struct expr::grid_dim<NamedData<G*>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G>
struct expr::grid_dim<NamedData<G>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<size_t Z, typename G>
struct expr::grid_dim<Variable<Z, G>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<size_t Z, typename G>
struct expr::grid_dim<Variable<Z, G*>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G>
struct expr::grid_dim<DynamicVariable<G>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Get the expression that the OpConvolution applies to.
template<typename T>
struct expr::grid_dim<SymbolicData<T>>
{
	static const size_t dimension = expr::grid_dim<T>::dimension;
	static const size_t value = dimension;
};

//! Get the expression that the OpConvolution applies to.
template<typename T>
struct expr::grid_dim<SymbolicDataArray<T>>
{
	static const size_t dimension = expr::grid_dim<T>::dimension;
	static const size_t value = dimension;
};

template<Axis ax, typename G>
struct expr::grid_dim<VectorComponent<ax, G>>
{
	static const size_t dimension = expr::grid_dim<G>::dimension;
	static const size_t value = dimension;
};

//! Get the expression that the OpConvolution applies to.
// template<typename... Ts>
// struct expr::grid_dim<SymbolicDataArray<std::tuple<Term<Ts>...>>>
// {
// 	static const size_t dimension = fixed_max<expr::grid_dim<Ts>::dimension...>;
// 	static const size_t value = dimension;
// };

//! Specialization based on expr::grid_dim.
template<typename T>
struct expr::grid_dim<OpLiteral<T>>
{
	static const size_t dimension = expr::grid_dim<T>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename T, typename I>
struct expr::grid_dim<OpCoeff<T, I>>
{
	static const size_t dimension = expr::grid_dim<OpLiteral<T>>::dimension;
	static const size_t value = expr::grid_dim<OpLiteral<T>>::value;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct expr::grid_dim<OpTerms<V, Term<Gs, Xs>...>>
{
	static const size_t dimension = fixed_max<expr::grid_dim<Gs>::dimension...>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename... Es>
struct expr::grid_dim<OpAdd<Es...>>
{
public:
	static const size_t dimension = fixed_max<expr::grid_dim<Es>::dimension...>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinaryMul<A, B>>
{
protected:
	static const size_t Ad = expr::grid_dim<A>::dimension;
	static const size_t Bd = expr::grid_dim<B>::dimension;

public:
	static const size_t dimension = fixed_max<Ad, Bd>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A, typename B>
struct expr::grid_dim<OpBinaryDiv<A, B>>
{
protected:
	static const size_t Ad = expr::grid_dim<A>::dimension;
	static const size_t Bd = expr::grid_dim<B>::dimension;

public:
	static const size_t dimension = fixed_max<Ad, Bd>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E1, typename E2>
struct expr::grid_dim<OpConvolution<V, E1, E2>>
{
protected:
	static const size_t Ad = expr::grid_dim<E1>::dimension;
	static const size_t Bd = expr::grid_dim<E2>::dimension;

public:
	static const size_t dimension = fixed_max<Ad, Bd>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<size_t D>
struct expr::grid_dim<GaussianSmoothing<D>>
{
	static const size_t dimension = D;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename G, typename V, typename E>
struct expr::grid_dim<OpMap<G, V, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct expr::grid_dim<OpFunction<V, E, F, Arg0, Args...>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<auto f, typename V, typename E>
struct expr::grid_dim<OpFunctionApply<f, V, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::grid_dim<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	static const size_t dimension = fixed_max<expr::grid_dim<sub_t>::dimension, expr::grid_dim<E>::dimension>;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::grid_dim<OpDerivative<Dd, V, E, Sp>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E, typename T>
struct expr::grid_dim<OpIntegral<V, E, T>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A1, typename A2, typename E>
struct expr::grid_dim<OpCombination<A1, A2, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename A1, typename A2, typename E>
struct expr::grid_dim<OpChain<A1, A2, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E>
struct expr::grid_dim<OpExponential<V, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<expr::exp_key_t X, typename V, typename E>
struct expr::grid_dim<OpPow<X, V, E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

//! Specialization based on expr::grid_dim.
template<typename E>
struct expr::grid_dim<OpOptimized<E>>
{
	static const size_t dimension = expr::grid_dim<E>::dimension;
	static const size_t value = dimension;
};

template<expr::NoiseType nt, typename T, size_t D>
struct expr::grid_dim<NoiseData<nt, T, D>>
{
	static const size_t dimension = D;
	static const size_t value = dimension;
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
struct expr::op_types<E&&>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename E>
struct expr::op_types<OpExpression<E>>
{
	using type = typename expr::op_types<E>::type;
};


//! Get the expression that the OpConvolution applies to.
template<typename T>
struct expr::op_types<SymbolicData<T>>
{
	using type = typename symphas::lib::types_list<T>;
};


//! Get the expression that the OpConvolution applies to.
template<typename... Ts>
struct expr::op_types<SymbolicDataArray<std::tuple<Ts...>>>
{
	using type = typename symphas::lib::combine_types_unique<Ts...>::type;
};

//! Get the expression that the OpConvolution applies to.
template<typename T>
struct expr::op_types<SymbolicDataArray<T>>
{
	using type = T;
};

//! Get the expression that the OpConvolution applies to.
template<typename... Ts>
struct expr::op_types<Substitution<Ts...>>
{
	using type = typename symphas::lib::combine_types_unique<typename op_types<Ts>::type...>::type;
};

//! Get the expression that the OpConvolution applies to.
template<expr::NoiseType nt, typename T, size_t D>
struct expr::op_types<NoiseData<nt, T, D>>
{
	using type = NoiseData<nt, T, D>;
};

//! Specialization based on expr::op_types.
template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct expr::op_types<OpTerms<V, Term<Gs, Xs>...>>
{
	using type = typename symphas::lib::combine_types_unique<Gs...>::type;
};

//! Specialization based on expr::op_types.
template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
struct expr::op_types<OpTerms<Term<G0, X0>, Term<Gs, Xs>...>>
{
	using type = typename symphas::lib::combine_types_unique<G0, Gs...>;
};

//! Specialization based on expr::op_types.
template<typename... Es>
struct expr::op_types<OpAdd<Es...>>
{
	using type = typename symphas::lib::combine_types_unique<typename expr::op_types<Es>::type...>::type;
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
template<typename V, typename E1, typename E2>
struct expr::op_types<OpConvolution<V, E1, E2>>
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
	using type = decltype(GaussianSmoothing<D>::data);
};

//! Specialization based on expr::op_types.
template<typename G, typename V, typename E>
struct expr::op_types<OpMap<G, V, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct expr::op_types<OpFunction<V, E, F, Arg0, Args...>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::grid_dim.
template<auto f, typename V, typename E>
struct expr::op_types<OpFunctionApply<f, V, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Get the expression that the OpConvolution applies to.
template<typename E, typename... Ts>
struct expr::op_types<SymbolicFunction<E, Ts...>>
{
	using type = typename symphas::lib::combine_types_unique<
		typename expr::op_types<E>::type, typename expr::op_types<Ts>::type...>::type;
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::op_types<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
protected:
	using At = typename expr::op_types<sub_t>::type;
	using Bt = typename expr::op_types<SymbolicFunction<E, Ts...>>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename E, typename... Ts, typename Seq, typename A, typename B, typename E0, typename... T0s>
struct expr::op_types<OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, Seq, A, B, SymbolicFunction<E0, T0s...>>>
{
protected:
	using At = typename expr::op_types<E>::type;
	using Bt = typename symphas::lib::expand_types_list<typename expr::op_types<SymbolicDataArray<Ts>>::type...>;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::op_types<OpDerivative<Dd, V, E, Sp>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename V, typename E, typename T>
struct expr::op_types<OpIntegral<V, E, T>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2>
struct expr::op_types<OpOperatorCombination<A1, A2>>
{
protected:
	using At = typename expr::op_types<A1>::type;
	using Bt = typename expr::op_types<A2>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2>
struct expr::op_types<OpOperatorChain<A1, A2>>
{
protected:
	using At = typename expr::op_types<A1>::type;
	using Bt = typename expr::op_types<A2>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2, typename E>
struct expr::op_types<OpCombination<A1, A2, E>>
{
protected:
	using At = typename expr::op_types<OpOperatorCombination<A1, A2>>::type;
	using Bt = typename expr::op_types<E>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename A1, typename A2, typename E>
struct expr::op_types<OpChain<A1, A2, E>>
{
protected:
	using At = typename expr::op_types<OpOperatorChain<A1, A2>>::type;
	using Bt = typename expr::op_types<E>::type;

public:
	using type = typename symphas::lib::combine_types_unique<At, Bt>::type;
};

//! Specialization based on expr::op_types.
template<typename V, typename E>
struct expr::op_types<OpExponential<V, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<expr::exp_key_t X, typename V, typename E>
struct expr::op_types<OpPow<X, V, E>>
{
	using type = typename expr::op_types<E>::type;
};

//! Specialization based on expr::op_types.
template<typename E>
struct expr::op_types<OpOptimized<E>>
{
	using type = typename expr::op_types<E>::type;
};


// *******************************************************************************


template<typename... Gs>
struct expr::term_types<symphas::lib::types_list<Gs...>>
{
	using type = symphas::lib::types_list<typename expr::original_data_type<Gs>::type...>;
};

template<typename E>
struct expr::term_types
{
	using type = typename term_types<typename expr::op_types<E>::type>::type;
};



// *******************************************************************************


namespace symphas::internal
{

	//! Get the first type of types that are in a type list.
	/*!
	 * If the given type is not a list of types, then type trait always maps
	 * to void.
	 */
	template<typename G>
	struct get_first_storage_type
	{
		using type = void;
	};

	//! Specialization based on get_first_type.
	template<template<typename, size_t> typename G, typename T, size_t D, typename... Gs>
	struct get_first_storage_type<symphas::lib::types_list<G<T, D>, Gs...>>
	{
		using type = G<T, D>;
	};

	//! Specialization based on get_first_type.
	template<typename G0, typename... Gs>
	struct get_first_storage_type<symphas::lib::types_list<G0, Gs...>>
	{
		using type = typename get_first_storage_type<symphas::lib::types_list<Gs...>>::type;
	};
}

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

}


//! Return the type of the grid in the expression as a complete grid.
template<typename E>
struct expr::storage_type
{
	using data_t = symphas::lib::type_at_index<0, symphas::lib::unroll_types_list<expr::term_types_t<E>>>;
	using storage_type_t = typename symphas::internal::get_first_storage_type<expr::term_types_t<E>>::type;
	using check_t = std::conditional_t<
		std::is_same<storage_type_t, void>::value, 
			std::conditional_t<std::is_same<data_t, void>::value, expr::symbols::Symbol, data_t>, 
			storage_type_t>;

	template<template<typename, size_t> typename enc_type>
	struct grid_class_wrap 
	{
		using type = expr::base_data_t<enc_type<
			expr::eval_type_t<E>,
			expr::grid_dim<E>::value>>;
	};

	struct block_class_wrap
	{
		using type = Block<expr::eval_type_t<E>>;
	};

	template<typename T>
	struct type_wrap
	{
		using type = T;
	};

	template<template<typename, size_t> typename G, typename T, size_t D>
	static grid_class_wrap<G> _pack_grid(G<T, D>)
	{
		return {};
	}

	template<typename T, size_t D>
	static grid_class_wrap<Grid> _pack_grid(GridData<T, D> const& g)
	{
		return {};
	}

	template<typename T>
	static auto _pack_grid(SymbolicData<T> const& g)
	{
		return _pack_grid(*g.data);
	}

	template<Axis ax, size_t N, typename T>
	static block_class_wrap _pack_grid(VectorComponent<ax, MultiBlock<N, T>> const& g)
	{
		return {};
	}

	template<typename T>
	static block_class_wrap _pack_grid(T*)
	{
		return {};
	}

	static type_wrap<expr::symbols::Symbol> _pack_grid(expr::symbols::Symbol)
	{
		return {};
	}

	static type_wrap<expr::eval_type_t<E>> _pack_grid(...)
	{
		return {};
	}


	static auto pack_grid(check_t g)
	{
		return _pack_grid(g);
	}


public:
	using type = typename std::invoke_result_t<decltype(&expr::storage_type<E>::pack_grid), check_t>::type;
};

template<typename E, typename... Ts>
struct expr::storage_type<SymbolicFunction<E, Ts...>>
{
	using type = typename storage_type<E>::type;
};


template<typename G>
struct expr::parent_storage_type
{
protected:

	template<typename T>
	static Block<T> _cast_grid(Block<T>)
	{
		return {};
	}

	template<typename T, size_t D>
	static Grid<T, D> _cast_grid(Grid<T, D>)
	{
		return {};
	}

	static auto cast_grid() -> decltype(_cast_grid(std::declval<G>()))
	{
		return {};
	}

public:

	using type = std::invoke_result_t<decltype(&expr::parent_storage_type<G>::cast_grid)>;


};

//! Specialization based on expr::has_state.
template<typename V, typename E1, typename E2>
struct expr::has_state<OpConvolution<V, E1, E2>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::has_state<OpDerivative<Dd, V, E, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename V, typename E, typename T>
struct expr::has_state<OpIntegral<V, E, T>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename Dd, typename V, typename G, typename Sp>
struct expr::has_state<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>>
{
	static const bool value = false;
};

template<size_t O, typename V, typename E, typename G0>
struct expr::has_state<OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>>>
{
	static const bool value = false;
};

//! Specialization based on expr::has_state.
template<typename G, typename V, typename E>
struct expr::has_state<OpMap<G, V, E>>
{
	static const bool value = true;
};

//! Specialization based on expr::grid_dim.
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct expr::has_state<OpFunction<V, E, F, Arg0, Args...>>
{
	static const bool value = expr::has_state<E>::value;
};

//! Specialization based on expr::grid_dim.
template<auto f, typename V, typename E>
struct expr::has_state<OpFunctionApply<f, V, E>>
{
	static const bool value = expr::has_state<E>::value;
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::has_state<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	static const bool value = expr::has_state<E>::value;
};

//! Specialization based on expr::has_state.
template<typename E>
struct expr::has_state<OpOptimized<E>>
{
	static const bool value = true;
};

//! Specialization based on expr::has_state.
template<typename V, typename E>
struct expr::has_state<OpExponential<V, E>>
{
	static const bool value = expr::has_state<E>::value;
};

template<expr::exp_key_t X, typename V, typename E>
struct expr::has_state<OpPow<X, V, E>>
{
	static const bool value = expr::has_state<E>::value;
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
template<typename... Es>
struct expr::has_state<OpAdd<Es...>>
{
	static const bool value = (expr::has_state<Es>::value || ...);
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
template<typename... Rs>
struct expr::has_state<symphas::lib::types_list<Rs...>>
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
template<typename T, typename I>
struct expr::is_operator_like<OpCoeff<T, I>>
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
template<size_t N, size_t D>
struct expr::is_operator_like<OpFractionLiteral<N, D>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<size_t N, size_t D>
struct expr::is_operator_like<OpNegFractionLiteral<N, D>>
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

template<Axis ax, size_t O, typename V, typename Sp>
struct expr::is_operator_like<OpOperatorDirectionalDerivative<ax, O, V, Sp>>
{
	static const bool value = true;
};

template<typename V, typename Sp, size_t... Os>
struct expr::is_operator_like<OpOperatorMixedDerivative<V, Sp, Os...>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_operator_like<OpDerivative<Dd, V, E, Sp>>
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
template<typename... Es>
struct expr::is_operator_like<OpAdd<Es...>>
{
	static const bool value = (expr::is_operator_like<Es>::value && ...);
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::is_operator_like<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	static const bool value = expr::is_operator_like<E>::value;
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
template<typename T, typename I>
struct expr::is_linear<OpCoeff<T, I>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<size_t N, size_t D>
struct expr::is_linear<OpFractionLiteral<N, D>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_operator_like.
template<size_t N, size_t D>
struct expr::is_linear<OpNegFractionLiteral<N, D>>
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
struct expr::is_linear<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_linear<OpDerivative<Dd, V, E, Sp>>
{
	static const bool value = expr::is_linear<E>::value;
};

//! Specialization based on expr::is_linear.
template<size_t O, typename V, typename Sp>
struct expr::is_linear<OpOperatorDerivative<O, V, Sp>>
{
	static const bool value = true;
};

template<Axis ax, size_t O, typename V, typename Sp>
struct expr::is_linear<OpOperatorDirectionalDerivative<ax, O, V, Sp>>
{
	static const bool value = true;
};

template<typename V, typename Sp, size_t... Os>
struct expr::is_linear<OpOperatorMixedDerivative<V, Sp, Os...>>
{
	static const bool value = true;
};

//! Specialization based on expr::is_linear.
template<typename V, typename E1, typename E2>
struct expr::is_linear<OpConvolution<V, E1, E2>>
{
	static const bool value = expr::is_linear<E1>::value && expr::is_linear<E2>::value;
};

//! Specialization based on expr::is_linear.
template<typename T, typename G>
struct expr::is_linear<OpTerm<T, G>>
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
template<typename... Es>
struct expr::is_linear<OpAdd<Es...>>
{
	static const bool value = (expr::is_linear<Es>::value && ...);
};

//! Get the expression that the OpConvolution applies to.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::is_linear<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	static const bool value = expr::is_operator_like<E>::value;
};



// *******************************************************************************


//! Specialization based on expr::is_nonlinear.
template<typename E>
struct expr::is_nonlinear<OpExpression<E>>
{
	static const bool value = expr::is_nonlinear<E>::value;
};

//! Specialization based on expr::is_nonlinear.
template<typename... Es>
struct expr::is_nonlinear<OpAdd<Es...>>
{
	static const bool value = (expr::is_nonlinear<Es>::value && ...);
};





//! Specialization based on expr::is_combination.
template<typename... Es>
struct expr::is_combination<OpAdd<Es...>>
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

template<typename V, typename E, typename T, typename Seq, typename A, typename B, typename C>
struct expr::is_combination<OpSum<V, E, T, Seq, A, B, C>>
{
	static const bool value = true;
};

template<typename V, typename E>
struct expr::is_combination<OpMap<void, V, E>>
{
	static const bool value = true;
};

template<typename V, typename E>
struct expr::is_combination<OpMap<symphas::internal::STHC, V, E>>
{
	static const bool value = true;
};

template<typename V, typename E>
struct expr::is_combination<OpMap<symphas::internal::HCTS, V, E>>
{
	static const bool value = true;
};

template<typename S, typename T, size_t D, typename V, typename E>
struct expr::is_combination<OpMap<MapGridFourier<S, T, D>, V, E>>
{
	static const bool value = true;
};

template<typename S, typename T, size_t D, typename V, typename E>
struct expr::is_combination<OpMap<MapGridInverseFourier<S, T, D>, V, E>>
{
	static const bool value = true;
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

template<typename V, typename Sp, size_t... Os>
struct expr::is_directional_derivative<OpOperatorMixedDerivative<V, Sp, Os...>>
{
	static const bool value = true;
};

template<typename Dd, typename V, typename E, typename Sp>
struct expr::is_directional_derivative<OpDerivative<Dd, V, E, Sp>>
{
	static const bool value = expr::is_directional_derivative<Dd>::value;
};


namespace expr
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


//! Specialization based on expr::vars.
template<size_t Z, typename G>
struct expr::vars<Variable<Z, G>>
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

//! Specialization based on expr::vars.
template<typename G>
struct expr::vars<DynamicVariable<G>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<G>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<G>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<G>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<G>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<Axis ax, typename G>
struct expr::vars<VectorComponent<ax, G>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<G>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<G>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<G>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<G>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename T>
struct expr::vars<SymbolicData<T>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<T>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<T>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<T>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<T>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename T>
struct expr::vars<SymbolicDataArray<T>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<T>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<T>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<T>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<T>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename... Ts>
struct expr::vars<SymbolicDataArray<std::tuple<Term<Ts>...>>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<Ts>::get_ids()...));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return ((expr::vars<Ts>::template has_id<Y>() || ...));
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return ((expr::vars<Ts>::template has_id<Y>() && ...));
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return ((expr::vars<Ts>::template each_id<Y>() && ...));
	}
};

//! Specialization based on expr::vars.
template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct expr::vars<OpTerms<V, Term<Gs, Xs>...>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<Gs>::get_ids()...));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return ((expr::vars<Gs>::template has_id<Y>() || ...));
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return ((expr::vars<Gs>::template has_id<Y>() && ...));
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return has_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<auto f, typename V, typename E>
struct expr::vars<OpFunctionApply<f, V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::vars.
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct expr::vars<OpFunction<V, E, F, Arg0, Args...>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename V, typename sub_t, typename E, typename... Ts>
struct expr::vars<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename V, typename E, typename... Ts, typename Seq, typename A, typename B, typename C>
struct expr::vars<OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, Seq, A, B, C>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(
			expr::vars<E>::get_ids(), expr::vars<SymbolicDataArray<Ts>>::get_ids()...));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return (expr::vars<SymbolicDataArray<Ts>>::template has_id<Y>() || ... || expr::vars<E>::template has_id<Y>());
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return (expr::vars<SymbolicDataArray<Ts>>::template only_id<Y>() && ... && expr::vars<E>::template only_id<Y>());
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return (expr::vars<SymbolicDataArray<Ts>>::template each_id<Y>() && ... && expr::vars<E>::template each_id<Y>());
	}
};


//! Specialization based on expr::vars.
template<typename G, typename V, typename E>
struct expr::vars<OpMap<G, V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename Dd, typename V, typename E, typename Sp>
struct expr::vars<OpDerivative<Dd, V, E, Sp>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename V, typename E, typename T>
struct expr::vars<OpIntegral<V, E, T>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename A1, typename A2, typename E>
struct expr::vars<OpChain<A1, A2, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::vars.
template<typename A1, typename A2, typename E>
struct expr::vars<OpCombination<A1, A2, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};


//! Specialization based on expr::vars.
template<typename V, typename E1, typename E2>
struct expr::vars<OpConvolution<V, E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<E1>::get_ids(), expr::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E1>::template has_id<Y>() || expr::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E1>::template only_id<Y>() && expr::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E1>::template each_id<Y>() && expr::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename V, size_t D, typename E>
struct expr::vars<OpConvolution<V, GaussianSmoothing<D>, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename V, typename E>
struct expr::vars<OpExponential<V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<expr::exp_key_t X, typename V, typename E>
struct expr::vars<OpPow<X, V, E>>
{
	static constexpr auto get_ids()
	{
		return expr::vars<E>::get_ids();
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename... Es>
struct expr::vars<OpAdd<Es...>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<Es>::get_ids()...));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return (expr::vars<Es>::template has_id<Y>() || ...);
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return (expr::vars<Es>::template only_id<Y>() && ...);
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return (expr::vars<Es>::template each_id<Y>() && ...);
	}
};

//! Specialization based on expr::vars.
template<typename E1, typename E2>
struct expr::vars<OpBinaryMul<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<E1>::get_ids(), expr::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E1>::template has_id<Y>() || expr::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E1>::template only_id<Y>() && expr::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E1>::template each_id<Y>() || expr::vars<E2>::template each_id<Y>();
	}
};

//! Specialization based on expr::vars.
template<typename E1, typename E2>
struct expr::vars<OpBinaryDiv<E1, E2>>
{
	static constexpr auto get_ids()
	{
		return symphas::lib::fold_unique_ids(symphas::lib::seq_join(expr::vars<E1>::get_ids(), expr::vars<E2>::get_ids()));
	}

	template<size_t Y>
	static constexpr bool has_id()
	{
		return expr::vars<E1>::template has_id<Y>() || expr::vars<E1>::template has_id<Y>();
	}

	template<size_t Y>
	static constexpr auto only_id()
	{
		return expr::vars<E1>::template only_id<Y>() && expr::vars<E2>::template only_id<Y>();
	}

	template<size_t Y>
	static constexpr auto each_id()
	{
		return expr::vars<E1>::template each_id<Y>() || expr::vars<E2>::template each_id<Y>();
	}
};





//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2>
struct expr::derivative_index<L, OpOperatorCombination<A1, A2>>
{
protected:
	static const size_t raw_value = 
		(expr::derivative_index<L, A1>::value 
			| expr::derivative_index<L, A2>::value);

public:
	static const size_t value = raw_value;
};

//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2>
struct expr::derivative_index<L, OpOperatorChain<A1, A2>>
{
protected:

	static const size_t raw_value_1 = expr::derivative_index<L, A1>::value;
	static const size_t raw_value_2 = expr::derivative_index<L, A2>::value;

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


//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2, typename E>
struct expr::derivative_index<L, OpCombination<A1, A2, E>>
{
public:
	static const size_t value = expr::derivative_index<L, OpOperatorCombination<A1, A2>>::value;
};

//! Specialization based on expr::derivative_index.
template<size_t L, typename A1, typename A2, typename E>
struct expr::derivative_index<L, OpChain<A1, A2, E>>
{
public:
	static const size_t value = expr::derivative_index<L, OpOperatorChain<A1, A2>>::value;
};



//! Specialization based on expr::derivative_index.
template<size_t L, size_t O, typename V, typename Sp>
struct expr::derivative_index<L, OpOperatorDerivative<O, V, Sp>>
{
	static const size_t value = expr::derivative_index_raw<O>::value;
};

//! Specialization based on expr::derivative_index.
template<size_t L, Axis ax, size_t O, typename V, typename Sp>
struct expr::derivative_index<L, OpOperatorDirectionalDerivative<ax, O, V, Sp>>
{
	static const size_t value = expr::derivative_index_raw<O>::value;
};

template<size_t L, typename V, typename Sp, size_t... Os>
struct expr::derivative_index<L, OpOperatorMixedDerivative<V, Sp, Os...>>
{
	static const size_t value = expr::derivative_index_raw<(Os + ...)>::value;
};


//! Specialization based on expr::derivative_index.
template<size_t L, typename Dd, typename V, typename E, typename Sp>
struct expr::derivative_index<L, OpDerivative<Dd, V, E, Sp>>
{
protected:
	static const size_t raw_value = expr::derivative_index_raw<OpDerivative<Dd, V, E, Sp>::order>::value;

public:
	static const size_t value = raw_value;
};

//! Specialization based on expr::derivative_index.
template<size_t L, typename... Es>
struct expr::derivative_index<L, OpAdd<Es...>>
{
protected:
	static const size_t value_1 = fixed_min<expr::derivative_index<L, Es>::value...>;
	static const size_t value_2 = fixed_max<expr::derivative_index<L, Es>::value...>;

public:
	static const size_t value = (value_1 < L) ? value_2 : value_1;
};

//! Specialization based on expr::derivative_index.
template<size_t L, typename A, typename B>
struct expr::derivative_index<L, OpBinaryMul<A, B>>
{
	static const size_t value = expr::derivative_index<L, OpOperatorChain<A, B>>::value;
};

// ******************************************************************************************
// Identity expressions.
// ******************************************************************************************



//! Additive identity.
/*!
 * Expression representing the additive identity, a special type of operator
 * used primarily in the simplification of expressions.
 */
struct OpVoid : OpExpression<OpVoid>
{
	constexpr scalar_t eval(iter_type = 0) const
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

	constexpr operator scalar_t() const
	{
		return eval();
	}

	constexpr auto operator--(int) const;
};


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

	size_t print(char* out) const
	{
		return sprintf(out, "1");
	}

	size_t print_length() const
	{
		return 1;
	}

	constexpr auto operator--(int) const;

#endif

	constexpr auto operator-() const;

	operator int() const
	{
		return int(eval());
	}

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

	constexpr auto operator-() const;

	constexpr auto operator--(int) const;

	operator int() const
	{
		return int(eval());
	}

};

inline constexpr auto OpIdentity::operator-() const
{
	return OpNegIdentity{};
}

inline constexpr auto OpNegIdentity::operator-() const
{
	return OpIdentity{};
}



namespace expr
{
	namespace symbols
	{
		//! The value one.
		constexpr OpIdentity one{};

		//! The value zero.
		constexpr OpVoid zero{};
	}
}


// ******************************************************************************************
// Extract the coefficient from an expression.
// ******************************************************************************************



namespace expr
{

	template<typename E,
		typename std::enable_if_t<(!has_coeff<E> && !(expr::is_fraction<E> || expr::is_identity<E>)), int> = 0>
	auto coeff(E const& e)
	{
		return OpIdentity{};
	}

	template<typename T, typename I>
	OpCoeff<T, I> const& coeff(OpCoeff<T, I> const& e)
	{
		return e;
	}

	template<typename T, typename I>
	OpCoeff<T, I>& coeff(OpCoeff<T, I>& e)
	{
		return e;
	}

	inline auto coeff(OpVoid)
	{
		return OpVoid{};
	}

	template<typename T, size_t... Ns>
	decltype(auto) coeff(OpTensor<T, Ns...> const& tensor);
	template<typename... Es, size_t... Is>
	decltype(auto) coeff(OpAdd<Es...> const& e, std::index_sequence<Is...>);

	template<typename E, typename std::enable_if_t<has_coeff<E>, int> = 0>
	decltype(auto) coeff(OpExpression<E> const& e)
	{
		return expr::make_literal((*static_cast<E const*>(&e)).value);
	}

	template<typename E, typename std::enable_if_t<has_coeff<E>, int> = 0>
	decltype(auto) coeff(OpOperator<E> const& e)
	{
		return expr::make_literal((*static_cast<E const*>(&e)).value);
	}

	template<typename E, typename std::enable_if_t<!has_coeff<E>, int> = 0>
	decltype(auto) coeff(OpOperator<E> const& e)
	{
		return OpIdentity{};
	}

	template<typename E,
		typename std::enable_if_t<(!has_coeff<E> && expr::is_coeff<E>), int> = 0>
	constexpr decltype(auto) coeff(OpExpression<E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

    template<typename G0, expr::exp_key_t X0>
    decltype(auto) coeff(Term<G0, X0> const& term)
    {
        return OpIdentity{};
    }

	template<typename... Vs>
	decltype(auto) coeff(OpTerms<Vs...> const& e)
	{
		return expr::coeff(e.term);
	}

	template<typename A, typename B>
	constexpr decltype(auto) coeff(OpBinaryMul<A, B> const& e)
	{
		return coeff(e.a);
	}

	//template<typename... Es, size_t... Is>
	//auto coeff(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	//{
	//	return (coeff(std::get<Is>(e)).value + ...);
	//}

	template<typename... Es, typename std::enable_if_t<(is_coeff<Es> && ...), int> = 0>
	decltype(auto) coeff(OpAdd<Es...> const& e)
	{
		return e;
	}

	template<typename E, typename Enable>
	struct coeff_type;

	template<typename E>
	struct coeff_type<E, std::enable_if_t<(!has_coeff<E> && !(is_identity<E> || is_fraction<E> || is_tensor<E>)), int>>
	{
		using type = OpIdentity;
	};

	template<typename E>
	struct coeff_type<E, std::enable_if_t<has_coeff<E>, int>>
	{
		using type = decltype(expr::coeff(std::declval<E>()));
	};

	template<typename coeff_t>
	struct coeff_type<coeff_t, std::enable_if_t<(is_identity<coeff_t> || is_fraction<coeff_t> || is_tensor<coeff_t>), int>>
	{
		using type = coeff_t;
	};

	//template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	//struct coeff_type<OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>,
	//	std::enable_if_t<has_coeff<OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>>, int>>
	//{
	//	using type = V;
	//};

	//template<typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	//struct coeff_type<OpTerms<Term<G0, X0>, Term<Gs, Xs>...>, int>;

	template<typename E>
	using coeff_t = typename coeff_type<E, int>::type;

}


// *******************************************************************************

namespace symphas::internal
{

	template<typename E0>
	struct test_eval
	{
		using type = decltype(std::declval<OpIdentity>() * std::declval<E0>());
	};

    template<>
    struct test_eval<void>
    {
        using type = void;
    };
}


template<typename E>
struct expr::eval_type
{
	template<typename T, size_t D>
	static constexpr std::index_sequence<1> _get_rank_1(any_vector_t<T, D> const&)
	{
		return {};
	}

	template<typename T, size_t D>
	static constexpr std::index_sequence<D> _get_rank_1(any_row_vector_t<T, D> const&)
	{
		return {};
	}

	template<typename T, size_t N, size_t M>
	static constexpr std::index_sequence<M> _get_rank_1(any_matrix_t<T, N, M> const&)
	{
		return {};
	}

	static constexpr std::index_sequence<0> _get_rank_1(...)
	{
		return {};
	}

	template<typename T, size_t D>
	static constexpr std::index_sequence<D> _get_rank(any_vector_t<T, D> const&)
	{
		return {};
	}

	template<typename T, size_t D>
	static constexpr std::index_sequence<1> _get_rank(any_row_vector_t<T, D> const&)
	{
		return {};
	}

	template<typename T, size_t N, size_t M>
	static constexpr std::index_sequence<N> _get_rank(any_matrix_t<T, N, M> const&)
	{
		return {};
	}

	static constexpr std::index_sequence<0> _get_rank(...)
	{
		return {};
	}

	template<typename E0>
	static constexpr std::index_sequence<expr::eval_type<E0>::rank> _get_rank(OpExpression<E0> const& e)
	{
		return {};
	}

	template<typename E0>
	static constexpr std::index_sequence<expr::eval_type<E0>::template rank_<1>> _get_rank_1(OpExpression<E0> const& e)
	{
		return {};
	}

	template<typename E0>
	static constexpr auto get_rank()
	{
		if constexpr (std::is_same<E0, void>::value)
		{
			return std::index_sequence<0>{};
		}
		else
		{
			return decltype(_get_rank(std::declval<E0>())){};
		}
	}

	template<typename E0>
	static constexpr auto get_rank_1()
	{
		if constexpr (std::is_same<E0, void>::value)
		{
			return std::index_sequence<0>{};
		}
		else
		{
			return decltype(_get_rank_1(std::declval<E0>())){};
		}
	}

	template<typename T0>
	static auto _get_eval(SymbolicDataArray<T0> const& e)
	{
		return typename eval_type<OpTerm<OpIdentity, DynamicVariable<T0>>>::type{};
	}

	template<typename T0>
	static auto _get_eval(SymbolicDataArray<NamedData<T0*>> const& e)
	{
		return typename eval_type<OpTerm<OpIdentity, T0>>::type{};
	}

	template<typename T0, typename... Ts>
	static auto _get_eval(SymbolicDataArray<std::tuple<T0, Ts...>> const& e)
	{
		return typename eval_type<OpTerm<OpIdentity, T0>>::type{};
	}

	template<typename T0, typename... Ts>
	static auto _get_eval(Substitution<T0, Ts...> const& e) -> typename eval_type<T0>::type
	{
		return {};
	}

	template<typename E0>
	static auto _get_eval(OpExpression<E0> const& e)
	{
		return static_cast<E0 const*>(&e)->eval(0);
	}

	template<typename G, expr::exp_key_t X>
	static auto _get_eval(Term<G, X> const& e)
	{
		return e.eval(0);
	}

	template<typename E0>
	static auto _get_eval(OpOperator<E0> const& e)
	{
		return static_cast<E0 const*>(&e)->eval(0);
	}

	static auto _get_eval(...)
	{
		return expr::symbols::Symbol{};
	}

	template<typename E0>
	static auto get_eval(E0 const& e0)
	{
		return _get_eval(e0);
	}

	template<typename E0>
	using eval_t = std::invoke_result_t<decltype(&eval_type<E>::template get_eval<E0>), E0>;
        
public:
	
	using type = typename symphas::internal::test_eval<eval_t<E>>::type;
	static constexpr size_t rank = symphas::lib::seq_index_value<0, std::invoke_result_t<decltype(&eval_type<E>::get_rank<type>)>>::value;
	
protected:

    static constexpr size_t rank_1 = symphas::lib::seq_index_value<0, std::invoke_result_t<decltype(&expr::eval_type<E>::get_rank_1<type>)>>::value;

public:

	template<size_t D>
	static constexpr size_t rank_ = (D == 0) ? rank : (D == 1) ? rank_1 : 0;
};





////////////////////////////////////////////////////////////////////////////////////////////////
//
//  library functions for symbolic expressions
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////


namespace expr
{

	namespace symbols
	{

		namespace
		{
			template<typename T>
			struct trait_arg_type
			{
				using type = std::conditional_t<is_simple_data<T>, OpLiteral<T>, T>;
			};
		}

		//! Substitutable variable for use in templates.
		/*!
		 * Object used as an enumerated argument for templates.
		 *
		 * Templates can be created without specifying an argument list, and
		 * the created template object will count the given expr::symbols::arg types used in the
		 * expression to understand how many arguments can be passed. When expressions are
		 * substituted in the template in order to expand it, their order in the template
		 * substitution list (the paramters to the `operator()` method of the template object)
		 * will be matched in order of increasing index of the argument objects used to construct
		 * the template expression.
		 */
		template<size_t N, typename T = expr::symbols::Symbol>
		using arg_t = OpTerm<OpIdentity, Variable<N, T>>;


	}

	namespace
	{
		template<size_t... Ns>
		struct SymbolicTemplateDef;
		template<typename... Gs>
		struct SymbolicTemplateDefSwap;
		template<typename... Gs>
		struct SymbolicFunctionDef;

		template<size_t... Ns>
		constexpr bool all_ne = true;
		template<size_t N0, size_t N1>
		constexpr bool all_ne<N0, N1> = N0 != N1;
		template<size_t N0, size_t N1, size_t N2, size_t... Ns>
		constexpr bool all_ne<N0, N1, N2, Ns...> = (N0 != N1) ? all_ne<N0, N2, Ns...> && all_ne<N1, N2, Ns...> : false;

		template<typename... Ts>
		constexpr bool all_different = true;
		template<typename T0, typename T1>
		constexpr bool all_different<T0, T1> = !std::is_same_v<T0, T1>;
		template<typename T0, typename T1, typename T2, typename... Ts>
		constexpr bool all_different<T0, T1, T2, Ts...> = (!std::is_same_v<T0, T1>) ? all_different<T0, T2, Ts...> && all_different<T1, T2, Ts...> : false;

	}

	//! Create a template that can create an expression by substituting parameters to placeholders.
	/*!
	 * Expression templates are essentially patterns that other expressions are substituted into, in
	 * order to build another expression.
	 *
	 * These templates can be created without specifying an argument list, and
	 * the created function object will count the given expr::symbols::arg_t types used in the
	 * expression to understand how many arguments can be passed. When expressions are
	 * substituted in the function in order to expand it, their order in the function
	 * substitution list (the paramters to the `operator()` method of the function object)
	 * will be matched in order of increasing index of the argument objects used to construct
	 * the function expression.
	 */
	template<typename... Ts>
	decltype(auto) template_of(Ts&&... ts);

	//! Define a symbolic function using the given argument types.
	/*!
	 * Functions are objects that substitute given data into their arguments so that the
	 * enclosed expression can be evaluated with any data that is swapped on the fly. See
	 * SymbolicFunction.
	 *
	 * Functions can be created without specifying an argument list, and
	 * the created function object will not take any arguments as placeholders. When expressions 
	 * are substituted in the function in order to expand it, their order in the function
	 * substitution list (the paramters to the `operator()` method of the function object)
	 * will be matched in order of increasing index of the argument objects used to construct
	 * the function expression.
	 */
	template<typename... Ts>
	decltype(auto) function_of(Ts&&... ts);


}



namespace expr
{
	//! Get the expression that the OpMap applies to.
	template<typename G, typename V, typename E>
	auto const& get_enclosed_expression(OpMap<G, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto const& get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
	{
		return e.e;
	}
	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto const& get_enclosed_expression(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename V, typename E>
	auto const& get_enclosed_expression(OpMap<void, V, E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename G, typename Sp>
	decltype(auto) get_enclosed_expression(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename G, typename G0>
	auto const& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G>, SymbolicDerivative<G0>> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename E, typename G0>
	auto const& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename E, typename Sp>
	auto const& get_enclosed_expression(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename V, typename E, typename T>
	auto const& get_enclosed_expression(OpIntegral<V, E, T> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpChain applies to.
	template<typename A1, typename A2, typename E>
	auto const& get_enclosed_expression(OpChain<A1, A2, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpCombination applies to.
	template<typename A1, typename A2, typename E>
	auto const& get_enclosed_expression(OpCombination<A1, A2, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpFunction applies to.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto const& get_enclosed_expression(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<auto f, typename V, typename E>
	auto const& get_enclosed_expression(OpFunctionApply<f, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename sub_t, typename E, typename... Ts>
	const auto& get_enclosed_expression(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		return e.f.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename G>
	auto get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}



	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename E>
	auto const& get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return e.e;
	}


	//! Get the expression that the OpExponential applies to.
	template<typename V, typename E>
	auto const& get_enclosed_expression(OpExponential<V, E> const& e)
	{
		return e.e;
	}

	template<expr::exp_key_t X, typename V, typename E>
	auto const& get_enclosed_expression(OpPow<X, V, E> const& e)
	{
		return e.e;
	}

	//! Get the expression that the OpOptimized applies to.
	template<typename E>
	decltype(auto) get_enclosed_expression(OpOptimized<E> const& e)
	{
		return e.get_expression();
	}

	//! Get the expression that the OpMap applies to.
	template<typename G, typename V, typename E>
	auto& get_enclosed_expression(OpMap<G, V, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto& get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename E>
	auto& get_enclosed_expression(OpMap<symphas::internal::STHC, OpIdentity, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpMap applies to.
	template<typename V, typename E>
	auto& get_enclosed_expression(OpMap<void, V, E>& e)
	{
		return *static_cast<E*>(&e);
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename G, typename Sp>
	decltype(auto) get_enclosed_expression(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename G, typename G0>
	auto& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G>, SymbolicDerivative<G0>>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<typename Dd, typename V, typename E, typename Sp>
	auto& get_enclosed_expression(OpDerivative<Dd, V, E, Sp>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpDerivative applies to.
	template<size_t O, typename V, typename E, typename G0>
	auto& get_enclosed_expression(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpIntegral applies to.
	template<typename V, typename E, typename T>
	auto& get_enclosed_expression(OpIntegral<V, E, T>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpChain applies to.
	template<typename A1, typename A2, typename E>
	auto& get_enclosed_expression(OpChain<A1, A2, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpCombination applies to.
	template<typename A1, typename A2, typename E>
	auto& get_enclosed_expression(OpCombination<A1, A2, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto& get_enclosed_expression(OpFunction<V, E, F, Arg0, Args...>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<auto f, typename V, typename E>
	auto& get_enclosed_expression(OpFunctionApply<f, V, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, typename sub_t, typename E, typename... Ts>
	auto& get_enclosed_expression(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>& e)
	{
		return e.f.e;
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename G>
	auto get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e)
	{
		return OpTerm<OpIdentity, G>(OpIdentity{}, e.data);
	}

	//! Get the expression that the OpConvolution applies to.
	template<typename V, size_t D, typename E>
	auto& get_enclosed_expression(OpConvolution<V, GaussianSmoothing<D>, E>& e)
	{
		return e.e;
	}

	template<expr::exp_key_t X, typename V, typename E>
	auto& get_enclosed_expression(OpPow<X, V, E> &e)
	{
		return e.e;
	}

	//! Get the expression that the OpExponential applies to.
	template<typename V, typename E>
	auto& get_enclosed_expression(OpExponential<V, E>& e)
	{
		return e.e;
	}

	//! Get the expression that the OpOptimized applies to.
	template<typename E>
	decltype(auto) get_enclosed_expression(OpOptimized<E>& e)
	{
		return e.get_expression();
	}


	template<typename E>
	auto& get_result_data(OpExpression<E>& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename G, typename Sp>
	auto& get_result_data(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename E, typename Sp>
	auto& get_result_data(OpDerivative<Dd, V, E, Sp>& e)
	{
		return e.grid;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<size_t O, typename V, typename E, typename G0>
	auto& get_result_data(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>>&) = delete;

	//! Get the grid storing the underlying data of the OpIntegral.
	template<typename V, typename E, typename T>
	auto& get_result_data(OpIntegral<V, E, T>& e)
	{
		return e.data;
	}

	//! Get the grid storing the underlying data of the OpIntegral.
	template<typename V, typename E, typename T>
	auto& get_result_data(OpIntegral<V, E, SymbolicDerivative<T>>& e) = delete;

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, typename E1, typename E2>
	auto& get_result_data(OpConvolution<V, E1, E2>& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename G>
	auto& get_result_data(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename E>
	auto& get_result_data(OpConvolution<V, GaussianSmoothing<D>, E>& e)
	{
		return e.g0;
	}

	//! Get the expression that the OpOptimized applies to.
	template<typename E>
	decltype(auto) get_result_data(OpOptimized<E>& e)
	{
		return e.grid;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename E, typename Sp>
	auto const& get_result_data(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return e.grid;
	}

	//! Get the grid storing the underlying data of the OpDerivative.
	template<typename Dd, typename V, typename G, typename Sp>
	auto const& get_result_data(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e) = delete;

	//! Get the grid storing the underlying data of the OpDerivative.
	template<size_t O, typename V, typename E, typename G0>
	auto const& get_result_data(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G0>> const&) = delete;

	//! Get the grid storing the underlying data of the OpIntegral.
	template<typename V, typename E, typename T>
	auto const& get_result_data(OpIntegral<V, E, T> const& e)
	{
		return e.data;
	}

	//! Get the grid storing the underlying data of the OpIntegral.
	template<typename V, typename E, typename T>
	auto const& get_result_data(OpIntegral<V, E, SymbolicDerivative<T>> const& e) = delete;

	//! Get the grid storing the underlying data of the OpExponential.
	template<typename V, typename E>
	auto const& get_result_data(OpExponential<V, E> const& e)
	{
		return e.data;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, typename E1, typename E2>
	auto const& get_result_data(OpConvolution<V, E1, E2> const& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename G>
	auto const& get_result_data(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e)
	{
		return e.g0;
	}

	//! Get the grid storing the underlying data of the OpConvolution.
	template<typename V, size_t D, typename E>
	auto const& get_result_data(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return e.g0;
	}

	//! Get the expression that the OpOptimized applies to.
	template<typename E>
	decltype(auto) get_result_data(OpOptimized<E> const& e)
	{
		return e.grid;
	}



	template<typename E>
	auto get_solver(E const& e);


	template<typename G0, exp_key_t X0>
	const auto& terms_after_first(OpTerms<Term<G0, X0>> const& e);
	template<typename G0, exp_key_t X0>
	auto& terms_after_first(OpTerms<Term<G0, X0>>& e);
	template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_first(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...> const& e);
	template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_first(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>& e);
	template<typename V>
	decltype(auto) terms_after_first(OpTerms<V> const& e);
	template<typename V>
	decltype(auto) terms_after_first(OpTerms<V>& e);
    template<typename... Ts>
    decltype(auto) terms_after_first(OpTermsList<Ts...> const& e);
    template<typename... Ts>
    decltype(auto) terms_after_first(OpTermsList<Ts...>& e);
	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTerms<V, Term<Gs, Xs>...> const& e);
	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTerms<V, Term<Gs, Xs>...>& e);
    template<typename... Ts>
    decltype(auto) terms_after_n(OpTermsList<Ts...> const& e);
    template<typename... Ts>
    decltype(auto) terms_after_n(OpTermsList<Ts...>& e);

	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...> const& e);
	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...>& e);

	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...> const& e);
	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...>& e);
    template<typename... Ts>
    decltype(auto) terms_after_n(OpAddList<Ts...> const& e);
    template<typename... Ts>
    decltype(auto) terms_after_n(OpAddList<Ts...>& e);
	template<typename E0, typename E1, typename E2, typename... Es>
	const auto& terms_after_first(OpAdd<E0, E1, E2, Es...> const& e);
	template<typename E0, typename E1>
	const auto& terms_after_first(OpAdd<E0, E1> const& e);
	template<typename E0, typename E1, typename E2, typename... Es>
	auto& terms_after_first(OpAdd<E0, E1, E2, Es...>& e);
	template<typename E0, typename E1>
	auto& terms_after_first(OpAdd<E0, E1>& e);
    template<typename... Ts>
    decltype(auto) terms_after_first(OpAddList<Ts...> const& e);
    template<typename... Ts>
    decltype(auto) terms_after_first(OpAddList<Ts...>& e);


    template<typename T>
    auto eval(T const& value);

}
