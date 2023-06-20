
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
#include "expressionprototypes.h"



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

	template<typename T, size_t D>
	struct data_value_type<any_vector_t<T, D>> 
	{
		using type = any_vector_t<T, D>;
		using ref = type&;

		ref operator()(any_vector_t<T, D>* data, iter_type n)
		{
			return data[n];
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


	template<typename specialized_iterator>
	struct iterator_difference_type_impl
	{
		iterator_difference_type_impl(std::ptrdiff_t pos) : pos{ pos } {}

		bool operator<(specialized_iterator const& other) const
		{
			return pos < other.pos;
		}

		bool operator<(size_t value) const
		{
			return pos < value;
		}

		//! Dereference the iterator.
		inline decltype(auto) operator+(specialized_iterator const& other) const
		{
			specialized_iterator diff(cast());
			diff.pos += other.pos;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator-(specialized_iterator const& other) const
		{
			specialized_iterator diff(cast());
			diff.pos -= other.pos;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator+(std::ptrdiff_t offset) const
		{
			specialized_iterator diff(cast());
			diff.pos += offset;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator-(std::ptrdiff_t offset) const
		{
			specialized_iterator diff(cast());
			diff.pos += offset;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator/(size_t div) const
		{
			specialized_iterator diff(cast());
			diff.pos /= div;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator/(specialized_iterator const& div) const
		{
			specialized_iterator diff(cast());
			diff.pos /= div.pos;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator%(size_t div) const
		{
			specialized_iterator diff(cast());
			diff.pos %= div;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator%(specialized_iterator const& div) const
		{
			specialized_iterator diff(cast());
			diff.pos %= div.pos;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator*(size_t div) const
		{
			specialized_iterator diff(cast());
			diff.pos *= div;
			return diff;
		};

		//! Dereference the iterator.
		inline decltype(auto) operator*(specialized_iterator const& div) const
		{
			specialized_iterator diff(cast());
			diff.pos *= div.pos;
			return diff;
		};

		//! Dereference the iterator.
		specialized_iterator& operator+=(specialized_iterator const& other)
		{
			pos += other.pos;
			return cast();
		};

		//! Dereference the iterator.
		specialized_iterator& operator-=(specialized_iterator const& other)
		{
			pos -= other.pos;
			return cast();
		};

		//! Dereference the iterator.
		specialized_iterator& operator+=(std::ptrdiff_t offset)
		{
			pos += offset;
			return cast();
		};

		//! Dereference the iterator.
		specialized_iterator& operator-=(std::ptrdiff_t offset)
		{
			pos += offset;
			return cast();
		};

		//! Prefix decrement, returns itself.
		specialized_iterator& operator--()
		{
			--pos;
			return cast();
		}

		//! Prefix decrement, returns itself.
		specialized_iterator& operator++()
		{
			++pos;
			return cast();
		}

        specialized_iterator operator-() const
        {
            specialized_iterator result(cast());
            result.pos = -result.pos;
            return result;
        }

		explicit operator size_t() const
		{
			return size_t(pos);
		}

		explicit operator int() const
		{
			return pos;
		}

		explicit operator std::ptrdiff_t() const
		{
			return std::ptrdiff_t(pos);
		}

		specialized_iterator& cast()
		{
			return *static_cast<specialized_iterator*>(this);
		}

		const specialized_iterator& cast() const
		{
			return *static_cast<specialized_iterator const*>(this);
		}

        std::ptrdiff_t pos;

	};

	template<typename G>
	struct iterator_difference_type : iterator_difference_type_impl<iterator_difference_type<G>>
	{
		using parent_type = iterator_difference_type_impl<iterator_difference_type<G>>;
		using parent_type::pos;

		iterator_difference_type(std::ptrdiff_t pos) : iterator_difference_type(nullptr, pos) {}
		iterator_difference_type(G* ptr, std::ptrdiff_t pos) : parent_type(pos), ptr{ ptr } {}
		iterator_difference_type() : iterator_difference_type(nullptr, 0) {}

		bool operator==(iterator_difference_type<G> const& other) const
		{
			return ptr == other.ptr && pos == other.pos;
		}

		bool operator==(size_t value) const
		{
			return pos == value;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(iterator_difference_type<G> const& other) const
		{
			return pos > other.pos && ptr == other.ptr;
		}

		bool operator>(size_t value) const
		{
			return pos > value;
		}

		G* ptr;

	};


	template<typename G>
	struct iterator_selection_difference_type : iterator_difference_type_impl<iterator_selection_difference_type<G>>
	{
		using parent_type = iterator_difference_type_impl<iterator_selection_difference_type<G>>;
		using parent_type::pos;

		iterator_selection_difference_type(std::ptrdiff_t pos) : iterator_selection_difference_type(nullptr, nullptr, pos) {}
		iterator_selection_difference_type(G* ptr, iter_type* iters, std::ptrdiff_t pos) : parent_type(pos), ptr{ ptr }, iters{ iters } {}
		iterator_selection_difference_type() : iterator_selection_difference_type(nullptr, nullptr, 0) {}


		bool operator==(iterator_selection_difference_type<G> const& other) const
		{
			return ptr == other.ptr && pos == other.pos;
		}

		bool operator==(size_t value) const
		{
			return pos == value;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(iterator_selection_difference_type<G> const& other) const
		{
			return pos > other.pos && ptr == other.ptr;
		}

		bool operator>(size_t value) const
		{
			return pos > value;
		}

		G* ptr;
		iter_type* iters;
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
		using difference_type = std::ptrdiff_t;//iterator_difference_type<G>;
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
        template<typename specialized_difference = iterator_difference_type<G>>
		data_iterator(iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}
		//data_iterator(difference_type ptr = {})
		//	: ptr{ ptr } {}

		data_iterator(difference_type pos = {})
			: ptr{ pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator(G& data, difference_type pos = {})
			: ptr{ &data, pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator(G* data, difference_type pos = {})
			: ptr{ data, pos } {}

		data_iterator(data_iterator<G> const& other) :
			data_iterator(other.ptr) {}
		data_iterator(data_iterator<G>&& other) :
			data_iterator(other.ptr) {}
		data_iterator<G>& operator=(data_iterator<G> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}


		//! Prefix increment, returns itself.
		data_iterator<G>& operator++()
		{
			++ptr.pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		data_iterator<G> operator++(int)
		{
			data_iterator<G> it(*this);
			++ptr.pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		data_iterator<G>& operator--()
		{
			--ptr.pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		data_iterator<G> operator--(int)
		{
			data_iterator<G> it(*this);
			--ptr.pos;
			return it;
		}

		data_iterator<G>& operator+=(difference_type offset)
		{
			ptr.pos += offset;
			return *this;
		}

		data_iterator<G>& operator-=(difference_type offset)
		{
			ptr.pos -= offset;
			return *this;
		}



		//! Dereference the iterator.
		inline decltype(auto) operator*()
		{
			return data_value_type<G>{}(ptr.ptr, ptr.pos);
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos)
		{
			return data_value_type<G>{}(ptr.ptr, iter_type(ptr.pos + given_pos));
		}

		//! Member access of the iterated expression.
		inline G* operator->()
		{
			return ptr.ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(data_iterator<G> const& other) const
		{
			return ptr == other.ptr;
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
			return ptr > other.ptr;
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
			return ptr;
		}

		//! Add two iterators.
		difference_type operator+(data_iterator<G> const& rhs) const
		{
			return difference_type(ptr + rhs.ptr);
		}

		//! Subtract two iterators.
		difference_type operator-(data_iterator<G> const& rhs) const
		{
			return difference_type(ptr - rhs.ptr);
		}

		//! Add an offset from the iterator.
		data_iterator<G> operator+(difference_type offset) const
		{
			data_iterator<G> it(*this);
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		data_iterator<G> operator-(difference_type offset) const
		{
			data_iterator<G> it(*this);
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend data_iterator<G> operator+(difference_type offset, data_iterator<G> rhs)
		{
            auto result(rhs);
            rhs.ptr = rhs.ptr + offset;
            return result;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend data_iterator<G> operator-(difference_type offset, data_iterator<G> rhs)
		{
            auto result(rhs);
            rhs.ptr = -rhs.ptr + offset;
            return result;
		}

        iterator_difference_type<G> ptr;
		//difference_type ptr;
		//G* ptr;				//!< Pointer to the expression that is iterated.
		//difference_type pos;	//!< Current index of iteration.
	};

	template<typename G>
	data_iterator(G*)->data_iterator<G>;

	template<typename G>
	data_iterator(G*, int)->data_iterator<G>;



	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam G The expression type which is evaluated.
	 */
	template<typename G>
	class data_iterator_select
	{

	public:

		using iterator_category = std::random_access_iterator_tag;
		using value_type = typename data_value_type<G>::type;
		using reference_type = typename data_value_type<G>::ref;
		using difference_type = std::ptrdiff_t;//iterator_selection_difference_type<G>;
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
		template<typename specialized_difference = iterator_selection_difference_type<G>>
		data_iterator_select(iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}
		//data_iterator_select(difference_type ptr = {})
		//	: ptr{ ptr } {}

		data_iterator_select(difference_type pos = {})
			: ptr{ pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator_select(G& data, iter_type* iters, difference_type pos = 0)
			: ptr{ &data, iters, pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit data_iterator_select(G* data, iter_type* iters, difference_type pos = 0)
			: ptr{ data, iters, pos } {}

		data_iterator_select(data_iterator_select<G> const& other) :
			data_iterator_select(other.ptr) {}
		data_iterator_select(data_iterator_select<G>&& other) :
			data_iterator_select(other.ptr) {}
		data_iterator_select<G>& operator=(data_iterator_select<G> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}


		//! Prefix increment, returns itself.
		data_iterator_select<G>& operator++()
		{
			++ptr.pos;
			return *this;
		}

		//! Postfix increment, return a copy before the increment.
		data_iterator_select<G> operator++(int)
		{
			data_iterator_select<G> it = *this;
			++ptr.pos;
			return it;
		}


		//! Prefix decrement, returns itself.
		data_iterator_select<G>& operator--()
		{
			--ptr.pos;
			return *this;
		}

		//! Postfix decrement, return a copy before the increment.
		data_iterator_select<G> operator--(int)
		{
			data_iterator_select<G> it = *this;
			--ptr.pos;
			return it;
		}

		data_iterator_select<G>& operator+=(difference_type offset)
		{
			ptr.pos += offset;
			return *this;
		}

		data_iterator_select<G>& operator-=(difference_type offset)
		{
			ptr.pos -= offset;
			return *this;
		}



		//! Dereference the iterator.
		inline decltype(auto) operator*()
		{
			return data_value_type<G>{}(ptr.ptr, ptr.iters[ptr.pos]);
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos)
		{
			return data_value_type<G>{}(ptr.ptr, ptr.iters[ptr.pos + given_pos]);
		}

		//! Member access of the iterated expression.
		inline G* operator->()
		{
			return ptr.ptr;
		};


		//! Equality comparison with another iterator.
		/*!
		 * Equality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator==(data_iterator_select<G> const& other) const
		{
			return ptr == other.ptr;
		}

		//! Inequality comparison with another iterator.
		/*!
		 * Inequality comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator!=(data_iterator_select<G> const& other) const
		{
			return !(*this == other);
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>(data_iterator_select<G> const& other) const
		{
			return ptr > other.ptr;
		}

		//! Comparison with another iterator.
		/*!
		 * Less than comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<(data_iterator_select<G> const& other) const
		{
			return other > *this;
		}

		//! Comparison with another iterator.
		/*!
		 * Greater than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator>=(data_iterator_select<G> const& other) const
		{
			return !(*this < other);
		}

		//! Comparison with another iterator.
		/*!
		 * Less than or equal to comparison with another iterator.
		 * Compares the current position.
		 */
		bool operator<=(data_iterator_select<G> const& other) const
		{
			return !(*this > other);
		}


		//! Convertible to the difference type of two iterators.
		operator difference_type() const
		{
			return ptr;
		}

		//! Add two iterators.
		difference_type operator+(data_iterator_select<G> const& rhs) const
		{
			return difference_type(ptr + rhs.ptr);
		}

		//! Subtract two iterators.
		difference_type operator-(data_iterator_select<G> const& rhs) const
		{
			return difference_type(ptr - rhs.ptr);
		}

		//! Add an offset from the iterator.
		//template<typename G0>
		data_iterator_select<G> operator+(difference_type offset) const
		{
			data_iterator_select<G> it = *this;
			return it += offset;
		}

		//! Subtract an offset from the iterator.
		//template<typename G0>
		data_iterator_select<G> operator-(difference_type offset) const
		{
			data_iterator_select<G> it = *this;
			return it -= offset;
		}

		//! Add an offset from the left hand side to an iterator.
		friend data_iterator_select<G> operator+(difference_type offset, data_iterator_select<G> rhs)
		{
            auto result(rhs);
            rhs.ptr = rhs.ptr + offset;
            return result;
			//return offset + rhs.ptr;
		}

		//! Subtract an offset from the left hand side to an iterator.
		friend data_iterator_select<G> operator-(difference_type offset, data_iterator_select<G> rhs)
		{
            auto result(rhs);
            rhs.ptr = -rhs.ptr + offset;
            return result;
			//return offset - rhs.ptr;
		}

        iterator_selection_difference_type<G> ptr;
		//difference_type ptr;
		//G* ptr;				//!< Pointer to the expression that is iterated.
		//difference_type pos;	//!< Current index of iteration.
	};

}




