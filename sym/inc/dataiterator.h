
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
 * PURPOSE: Defines an iterator that is used for iterating over data so that
 * expressions can be evaluated into the data.
 *
 * ***************************************************************************
 */

#pragma once

#include "grid.h"
#include "expressionlib.h"



namespace symphas::internal
{
	template<typename T>
	struct data_value_type
	{
		using type = T;
		using ref = type&;

		ref operator()(T* data, iter_type n)
		{
			return data[n];
		}
	};

	template<template<typename, size_t> typename G, typename T, size_t D>
	struct data_value_type<G<T, D>>
	{
		using type = T;
		using ref = type&;

		ref operator()(G<T, D> *data, iter_type n)
		{
			return (*data)[n];
		}
	};

	template<typename T>
	struct data_value_type<Block<T>>
	{
		using type = T;
		using ref = type&;

		ref operator()(Block<T> *data, iter_type n)
		{
			return (*data)[n];
		}
	};

	template<size_t N, typename T>
	struct data_value_type<MultiBlock<N, T>>
	{
		using type = any_vector_t<T, N>;
		using ref = multi_value<N, T>;

		ref operator()(MultiBlock<N, T>* data, iter_type n)
		{
			return (*data)[n];
		}
	};

	template<template<typename, size_t> typename G, typename T, size_t D>
	struct data_value_type<G<any_vector_t<T, D>, D>> : data_value_type<MultiBlock<D, T>>
	{
		using parent_type = data_value_type<MultiBlock<D, T>>;
		using typename parent_type::type;
		using typename parent_type::ref;
	};

	template<typename G>
	struct data_value_type<symphas::ref<G>> : data_value_type<G>
	{
		using parent_type = data_value_type<G>;
		using typename parent_type::type;
		using typename parent_type::ref;

		decltype(auto) operator()(symphas::ref<G>* data, iter_type n)
		{
			return parent_type::operator()(&data->get(), n);
		}
	};

	template<typename G>
	struct data_value_type<Term<G, 1>> : data_value_type<G>
	{
		using parent_type = data_value_type<G>;
		using typename parent_type::type;
		using typename parent_type::ref;

		decltype(auto) operator()(Term<G, 1>* data, iter_type n)
		{
			return parent_type::operator()(static_cast<G*>(data), n);
		}
	};


	template<>
	struct data_value_type<expr::symbols::Symbol>
	{
		using type = expr::symbols::Symbol;
		using ref = expr::symbols::Symbol;

		decltype(auto) operator()(expr::symbols::Symbol* data, iter_type n)
		{
			return expr::symbols::Symbol{};
		}
	};


	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam G The expression type which is evaluated.
	 */
	template<typename G>
	class data_iterator
	{
	public:

		using iterator_category = std::random_access_iterator_tag;
		using value_type = typename data_value_type<G>::type;
		using reference_type = typename data_value_type<G>::ref;
		using difference_type = int;
		using pointer = value_type*;
		using reference = int;

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position.
		 *
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		//data_iterator(difference_type pos = 0)
		//	: e_ptr{ nullptr }, pos{ pos } {
		//	printf("N%d\n", pos);
		//	printf("N%p\n", e_ptr);
		//}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator(G& data, difference_type pos = 0)
			: e_ptr{ &data }, pos{pos} {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator(G* data, difference_type pos = 0)
			: e_ptr{ data }, pos{ pos } {}

		data_iterator(data_iterator<G> const& other) :
			data_iterator(other.e_ptr, other.pos) {}
		data_iterator(data_iterator<G>&& other) :
			data_iterator(other.e_ptr, other.pos) {}
		data_iterator<G>& operator=(data_iterator<G> other)
		{
			e_ptr = other.e_ptr;
			pos = other.pos;
			return *this;
		}


		//! Prefix increment, returns itself.
		data_iterator<G>& operator++()
		{
			++pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		data_iterator<G> operator++(difference_type)
		{
			data_iterator<G> it = *this;
			++pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		data_iterator<G>& operator--()
		{
			--pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		data_iterator<G> operator--(difference_type)
		{
			data_iterator<G> it = *this;
			--pos;
			return it;
		}


		data_iterator<G>& operator+=(difference_type offset)
		{
			pos += offset;
			return *this;
		}

		data_iterator<G>& operator-=(difference_type offset)
		{
			pos -= offset;
			return *this;
		}



		//! Dereference the iterator.
		inline decltype(auto) operator*()
		{
			return data_value_type<G>{}(e_ptr, pos);
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos)
		{
			return data_value_type<G>{}(e_ptr, pos + given_pos);
		}

		//! Member access of the iterated expression.
		inline G* operator->()
		{
			return e_ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(data_iterator<G> const& other) const
		{
			return pos == other.pos
				&& e_ptr == other.e_ptr;
		}

		//! Inequality comparison with another iterator.
		/*!
		 * Inequality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator!=(data_iterator<G> const& other) const
		{
			return !(*this == other);
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(data_iterator<G> const& other) const
		{
			return pos > other.pos
				&& e_ptr == other.e_ptr;
		}

		//! Comparison with another iterator.
		/*!
		 * Less than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<(data_iterator<G> const& other) const
		{
			return other > *this;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>=(data_iterator<G> const& other) const
		{
			return !(*this < other);
		}

		//! Comparison with another iterator.
		/*!
		 * Less than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<=(data_iterator<G> const& other) const
		{
			return !(*this > other);
		}


		//! Convertible to the difference type of two iterators.
		operator difference_type() const
		{
			return pos;
		}

		//! Add two iterators.
		difference_type operator+(data_iterator<G> const& rhs)
		{
			return pos + rhs;
		}

		//! Subtract two iterators.
		difference_type operator-(data_iterator<G> const& rhs)
		{
			return pos - rhs;
		}

		//! Add an offset from the iterator.
		data_iterator<G> operator+(difference_type offset)
		{
			data_iterator<G> it = *this;
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		data_iterator<G> operator-(difference_type offset)
		{
			data_iterator<G> it = *this;
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend difference_type operator+(difference_type offset, data_iterator<G> rhs)
		{
			return offset + rhs.pos;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend difference_type operator-(difference_type offset, data_iterator<G> rhs)
		{
			return offset - rhs.pos;
		}


		G* e_ptr;				//!< Pointer to the expression that is iterated.
		difference_type pos;	//!< Current index of iteration.
	};

	template<typename G>
	data_iterator(G*)->data_iterator<G>;

	template<typename G>
	data_iterator(G*, int)->data_iterator<G>;
}

