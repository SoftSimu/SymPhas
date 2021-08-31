
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

#include "expressionlib.h"


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
		using value_type = typename expr::eval_type<E>::type;
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
		expression_iterator(difference_type pos = 0)
			: e_ptr{ nullptr }, pos{ pos }  {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 * 
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator(OpExpression<E> const& e, difference_type pos = 0)
			: e_ptr{ static_cast<E const*>(&e) }, pos{ pos } {}


		expression_iterator(expression_iterator<E> const& other) :
			expression_iterator(*other.e_ptr, other.pos) {}
		expression_iterator(expression_iterator<E>&& other) :
			expression_iterator(*other.e_ptr, other.pos) {}
		expression_iterator<E>& operator=(expression_iterator<E> other)
		{
			e_ptr = other.e_ptr;
			pos = other.pos;
			return *this;
		}


		//! Prefix increment, returns itself.
		expression_iterator<E>& operator++()
		{
			++pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		expression_iterator<E> operator++(difference_type)
		{
			expression_iterator<E> it = *this; 
			++pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		expression_iterator<E>& operator--()
		{
			--pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		expression_iterator<E> operator--(difference_type)
		{
			expression_iterator<E> it = *this;
			--pos;
			return it;
		}


		expression_iterator<E>& operator+=(difference_type offset)
		{
			pos += offset;
			return *this;
		}

		expression_iterator<E>& operator-=(difference_type offset)
		{
			pos -= offset;
			return *this;
		}



		//! Dereference the iterator.
		inline value_type operator*() const
		{
			return e_ptr->eval(pos);
		};

		//! Dereference past the iterator.
		inline value_type operator[](difference_type given_pos) const
		{
			return e_ptr->eval(pos + given_pos);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return e_ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(expression_iterator<E> const& other) const
		{
			return pos == other.pos 
				&& e_ptr == other.e_ptr;
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
			return pos > other.pos
				&& e_ptr == other.e_ptr;
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
			return pos;
		}

		//! Add two iterators.
		difference_type operator+(expression_iterator<E> const& rhs)
		{
			return pos + rhs;
		}

		//! Subtract two iterators.
		difference_type operator-(expression_iterator<E> const& rhs)
		{
			return pos - rhs;
		}

		//! Add an offset from the iterator.
		expression_iterator<E> operator+(difference_type offset)
		{
			expression_iterator<E> it = *this;
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		expression_iterator<E> operator-(difference_type offset)
		{
			expression_iterator<E> it = *this;
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend difference_type operator+(difference_type offset, expression_iterator<E> rhs)
		{
			return offset + rhs.pos;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend difference_type operator-(difference_type offset, expression_iterator<E> rhs)
		{
			return offset - rhs.pos;
		}


		E const* e_ptr;			//!< Pointer to the expression that is iterated.
		difference_type pos;	//!< Current index of iteration.
	};

}

