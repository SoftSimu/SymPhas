
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
 * PURPOSE: Defines an iterator that is used for evaluating expressions.
 * The primary purpose is to parallelize the evaluation of an expression.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionprototypes.h"
#include "dataiterator.h"


namespace symphas::internal
{

	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 * 
	 * \tparam E The expression type which is evaluated.
	 */
	template<typename E>
	class expression_iterator
	{

	public:

		using iterator_category = std::random_access_iterator_tag;
		using value_type = std::invoke_result_t<decltype(&E::eval), E, iter_type>;
		using difference_type = iterator_difference_type<const E>;
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
		expression_iterator(difference_type ptr = {})
			: ptr{ ptr } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 * 
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator(OpExpression<E> const& e, iter_type pos = 0)
			: ptr{ static_cast<E const*>(&e), pos } {}


		expression_iterator(expression_iterator<E> const& other) :
			expression_iterator(other.ptr) {}
		expression_iterator(expression_iterator<E>&& other) :
			expression_iterator(other.ptr) {}
		expression_iterator<E>& operator=(expression_iterator<E> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}


		//! Prefix increment, returns itself.
		expression_iterator<E>& operator++()
		{
			++ptr.pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		expression_iterator<E> operator++(int)
		{
			expression_iterator<E> it = *this; 
			++ptr.pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		expression_iterator<E>& operator--()
		{
			--ptr.pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		expression_iterator<E> operator--(int)
		{
			expression_iterator<E> it = *this;
			--ptr.pos;
			return it;
		}


		template<typename E0>
		expression_iterator<E>& operator+=(iterator_difference_type<E0> offset)
		{
			ptr.pos += offset.pos;
			return *this;
		}

		template<typename E0>
		expression_iterator<E>& operator-=(iterator_difference_type<E0> offset)
		{
			ptr.pos -= offset.pos;
			return *this;
		}



		//! Dereference the iterator.
		inline value_type operator*() const
		{
			return (ptr.ptr)->eval(ptr.pos);
		};

		//! Dereference past the iterator.
		inline value_type operator[](difference_type given_pos) const
		{
			return (ptr.ptr)->eval(ptr.pos + ptr.given_pos);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(expression_iterator<E> const& other) const
		{
			return ptr == other.ptr;
		}

		//! Inequality comparison with another iterator.
		/*!
		 * Inequality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator!=(expression_iterator<E> const& other) const
		{
			return !(*this == other);
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(expression_iterator<E> const& other) const
		{
			return ptr > other.ptr;
		}

		//! Comparison with another iterator.
		/*!
		 * Less than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<(expression_iterator<E> const& other) const
		{
			return other > *this;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>=(expression_iterator<E> const& other) const
		{
			return !(*this < other);
		}

		//! Comparison with another iterator.
		/*!
		 * Less than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<=(expression_iterator<E> const& other) const
		{
			return !(*this > other);
		}



		//! Convertible to the difference type of two iterators.
		operator difference_type() const
		{
			return ptr;
		}

		//! Add two iterators.
		difference_type operator+(expression_iterator<E> const& rhs) const
		{
			return ptr + rhs;
		}

		//! Subtract two iterators.
		difference_type operator-(expression_iterator<E> const& rhs) const
		{
			return ptr - rhs;
		}

		//! Add an offset from the iterator.
		template<typename E0>
		expression_iterator<E> operator+(iterator_difference_type<E0> offset) const
		{
			expression_iterator<E> it = *this;
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		template<typename E0>
		expression_iterator<E> operator-(iterator_difference_type<E0> offset) const
		{
			expression_iterator<E> it = *this;
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend difference_type operator+(difference_type offset, expression_iterator<E> rhs)
		{
			return offset + rhs.ptr;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend difference_type operator-(difference_type offset, expression_iterator<E> rhs)
		{
			return offset - rhs.ptr;
		}


		//E const* e_ptr;		//!< Pointer to the expression that is iterated.
		difference_type ptr;	//!< Current index of iteration.
	};


	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 * 
	 * \tparam E The expression type which is evaluated.
	 */
	template<typename E>
	class expression_iterator_selection
	{

	public:

		using iterator_category = std::random_access_iterator_tag;
		using value_type = std::invoke_result_t<decltype(&E::eval), E, iter_type>;
		using difference_type = iterator_selection_difference_type<const E>;
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
		expression_iterator_selection(difference_type ptr = {})
			: ptr{ ptr } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 * 
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator_selection(OpExpression<E> const& e, iter_type* iters, iter_type pos = 0)
			: ptr{ static_cast<E const*>(&e), iters, pos } {}


		expression_iterator_selection(expression_iterator_selection<E> const& other) :
			expression_iterator_selection(other.ptr) {}
		expression_iterator_selection(expression_iterator_selection<E>&& other) :
			expression_iterator_selection(other.ptr) {}
		expression_iterator_selection<E>& operator=(expression_iterator_selection<E> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}


		//! Prefix increment, returns itself.
		expression_iterator_selection<E>& operator++()
		{
			++ptr.pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		expression_iterator_selection<E> operator++(int)
		{
			expression_iterator_selection<E> it = *this; 
			++ptr.pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		expression_iterator_selection<E>& operator--()
		{
			--ptr.pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		expression_iterator_selection<E> operator--(int)
		{
			expression_iterator_selection<E> it = *this;
			--ptr.pos;
			return it;
		}


		template<typename E0>
		expression_iterator_selection<E>& operator+=(iterator_difference_type<E0> offset)
		{
			ptr.pos += offset.pos;
			return *this;
		}

		template<typename E0>
		expression_iterator_selection<E>& operator-=(iterator_difference_type<E0> offset)
		{
			ptr.pos -= offset.pos;
			return *this;
		}



		//! Dereference the iterator.
		inline value_type operator*() const
		{
			return (ptr.ptr)->eval(ptr.iters[ptr.pos]);
		};

		//! Dereference past the iterator.
		inline value_type operator[](difference_type given_pos) const
		{
			return (ptr.ptr)->eval(ptr.iters[ptr.pos + ptr.given_pos]);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(expression_iterator_selection<E> const& other) const
		{
			return ptr == other.ptr;
		}

		//! Inequality comparison with another iterator.
		/*!
		 * Inequality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator!=(expression_iterator_selection<E> const& other) const
		{
			return !(*this == other);
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(expression_iterator_selection<E> const& other) const
		{
			return ptr > other.ptr;
		}

		//! Comparison with another iterator.
		/*!
		 * Less than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<(expression_iterator_selection<E> const& other) const
		{
			return other > *this;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>=(expression_iterator_selection<E> const& other) const
		{
			return !(*this < other);
		}

		//! Comparison with another iterator.
		/*!
		 * Less than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<=(expression_iterator_selection<E> const& other) const
		{
			return !(*this > other);
		}



		//! Convertible to the difference type of two iterators.
		operator difference_type() const
		{
			return ptr;
		}

		//! Add two iterators.
		difference_type operator+(expression_iterator_selection<E> const& rhs) const
		{
			return ptr + rhs;
		}

		//! Subtract two iterators.
		difference_type operator-(expression_iterator_selection<E> const& rhs) const
		{
			return ptr - rhs;
		}

		//! Add an offset from the iterator.
		template<typename E0>
		expression_iterator_selection<E> operator+(iterator_difference_type<E0> offset) const
		{
			expression_iterator_selection<E> it = *this;
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		template<typename E0>
		expression_iterator_selection<E> operator-(iterator_difference_type<E0> offset) const
		{
			expression_iterator_selection<E> it = *this;
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend difference_type operator+(difference_type offset, expression_iterator_selection<E> rhs)
		{
			return offset + rhs.ptr;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend difference_type operator-(difference_type offset, expression_iterator_selection<E> rhs)
		{
			return offset - rhs.ptr;
		}


		//E const* e_ptr;		//!< Pointer to the expression that is iterated.
		difference_type ptr;	//!< Current index of iteration.
		iter_type* iters;
	};

}


