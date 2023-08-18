
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

}

namespace expr
{


	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam E The expression type which is evaluated.
	 */
	template<typename E>
	struct expression_iterator : symphas::iterator_type_impl<expression_iterator<E>,
		std::random_access_iterator_tag,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::ptrdiff_t,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>*,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>
	>
	{

		using difference_type = std::ptrdiff_t;

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position.
		 *
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		template<typename specialized_difference = symphas::internal::iterator_difference_type<const E>>
		expression_iterator(symphas::internal::iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}

		expression_iterator(difference_type pos = {})
			: ptr{ pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator(OpExpression<E> const& e, difference_type pos = {})
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



		//! Dereference the iterator.
		inline decltype(auto) operator*() const
		{
			return (ptr.ptr)->eval(iter_type(ptr.pos));
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos) const
		{
			return (ptr.ptr)->eval(iter_type(ptr.pos + given_pos));
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return ptr;
		};

		symphas::internal::iterator_difference_type<const E> ptr;  //!< Current index of iteration.
	};


	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam E The expression type which is evaluated.
	 */
	template<typename E>
	struct expression_iterator_selection : symphas::iterator_type_impl<expression_iterator_selection<E>,
		std::random_access_iterator_tag,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::ptrdiff_t,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>*,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>
	>
	{

		using difference_type = std::ptrdiff_t;

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position.
		 *
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		template<typename specialized_difference = symphas::internal::iterator_selection_difference_type<const E>>
		expression_iterator_selection(symphas::internal::iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}

		expression_iterator_selection(difference_type pos = {})
			: ptr{ pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator_selection(OpExpression<E> const& e, const iter_type* iters, difference_type pos = {})
			: ptr{ static_cast<E const*>(&e), iters, pos } {}


		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param data The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		template<size_t D>
		explicit expression_iterator_selection(OpExpression<E> const& e, grid::region_index_list<D> const& list, difference_type pos = 0)
			: expression_iterator_selection(*static_cast<E const*>(&e), list.iters, pos) {}


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

		//! Dereference the iterator.
		inline decltype(auto) operator*() const
		{
			return (ptr.ptr)->eval(ptr.iters[ptr.pos]);
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos) const
		{
			return (ptr.ptr)->eval(ptr.iters[ptr.pos + given_pos]);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return ptr;
		};

		symphas::internal::iterator_selection_difference_type<const E> ptr;	//!< Current index of iteration.
	};




	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam E The expression type which is evaluated.
	 */
	template<typename E, size_t D>
	struct expression_iterator_region : symphas::iterator_type_impl<expression_iterator_region<E, D>,
		std::random_access_iterator_tag,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::ptrdiff_t,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>*,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>
	>
	{

		using difference_type = std::ptrdiff_t;

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position.
		 *
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		template<typename specialized_difference = symphas::internal::iterator_region_difference_type<const E, D>>
		expression_iterator_region(symphas::internal::iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}

		expression_iterator_region(difference_type pos = {})
			: ptr{ pos } {}

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit expression_iterator_region(OpExpression<E> const& e, grid::region_interval<D> const& interval, difference_type pos = {})
			: ptr{ static_cast<E const*>(&e), interval, pos } {}


		expression_iterator_region(expression_iterator_region<E, D> const& other) :
			expression_iterator_region(other.ptr) {}
		expression_iterator_region(expression_iterator_region<E, D>&& other) :
			expression_iterator_region(other.ptr) {}
		expression_iterator_region<E, D>& operator=(expression_iterator_region<E, D> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}



		//! Dereference the iterator.
		inline decltype(auto) operator*() const
		{
			return (ptr.ptr)->eval(*ptr);
		};

		//! Dereference past the iterator.
		inline decltype(auto) operator[](difference_type given_pos) const
		{
			return (ptr.ptr)->eval(ptr[given_pos]);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return ptr;
		};


		//E const* e_ptr;		//!< Pointer to the expression that is iterated.
		symphas::internal::iterator_region_difference_type<const E, D> ptr;//difference_type ptr;	//!< Current index of iteration.
	};


	//! An iterator used to evaluate an expression on its underlying data.
	/*!
	 * Implements the forward iterator functionality, in order
	 * to support parallelization using the standard library.
	 *
	 * \tparam G The expression type which is evaluated.
	 */
	template<typename E, size_t D>
	struct expression_iterator_group : symphas::iterator_type_impl<expression_iterator_group<E, D>,
		std::random_access_iterator_tag,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>,
		std::ptrdiff_t,
		std::invoke_result_t<decltype(&E::eval), E, iter_type>*,
		symphas::internal::iterator_group_expression<const E, D>
	>
	{

		using difference_type = std::ptrdiff_t;//iterator_regionion_difference_type<G>;


		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position.
		 *
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		template<typename specialized_difference = symphas::internal::iterator_region_difference_type<E, D>>
		expression_iterator_group(symphas::internal::iterator_difference_type_impl<specialized_difference> const& ptr)
			: ptr{ *static_cast<specialized_difference const*>(&ptr) } {}
		//expression_iterator_group(difference_type ptr = {})
		//	: ptr{ ptr } {}

		expression_iterator_group(difference_type pos = {})
			: ptr{ pos } {}
		explicit expression_iterator_group(OpExpression<E> const& e, grid::region_interval<D> const& interval, difference_type pos = {})
			: ptr{ static_cast<E const*>(&e), interval, pos } {}


		expression_iterator_group(expression_iterator_group<E, D> const& other) :
			expression_iterator_group(other.ptr) {}
		expression_iterator_group(expression_iterator_group<E, D>&& other) :
			expression_iterator_group(other.ptr) {}
		expression_iterator_group<E, D>& operator=(expression_iterator_group<E, D> other)
		{
			using std::swap;
			swap(ptr, other.ptr);
			return *this;
		}

		//! Dereference the iterator.
		inline auto operator*() const
		{
			return symphas::internal::iterator_group_expression(ptr);// symphas::internal::data_value_type<E>{}(ptr.ptr, ptr.get_index());
		};

		//! Dereference past the iterator.
		inline auto operator[](difference_type given_pos) const
		{
			return symphas::internal::iterator_group_expression(ptr, given_pos);//symphas::internal::data_value_type<E>{}(ptr.ptr, ptr.get_index(given_pos));
		}

		//! Member access of the iterated expression.
		inline E* operator->()
		{
			return ptr.ptr;
		};


		symphas::internal::iterator_group_difference_type<const E, D> ptr;
	};

}


namespace symphas::internal
{
	template<typename operation_t, typename A, typename B>
	struct iterator_group_operation;

	template<typename E, size_t D>
	struct iterator_group_expression
	{
		using eval_type = std::invoke_result_t<decltype(&E::eval), E, iter_type>;

		iterator_group_expression(symphas::internal::iterator_group_difference_type<E, D> const& ptr, iter_type offset = 0) : ptr{ ptr }, offset{ offset } {}
		iterator_group_expression() : ptr{}, offset{} {}

		inline len_type len(grid::region_group<D> const& region) const
		{
			return region.dims[0];
		}

		inline len_type index(grid::region_group<D> const& region, iter_type pos0) const
		{
			return region[pos0 * len(region)];
		}

		template<typename E0>
		auto operator+(iterator_group_expression<E0, D> const& other) const
		{
			eval_type reduce{};
			auto n0 = index(ptr.region, ptr.pos);
			auto n1 = index(other.ptr.region, other.ptr.pos);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				reduce += (ptr.ptr)->eval(n0 + n) + (other.ptr.ptr)->eval(n1 + n);
			}
			return reduce;
		}

		void accumulate(eval_type& value) const
		{
			auto n0 = index(ptr.region, ptr.pos);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				value += (ptr.ptr)->eval(n0 + n);
			}
		}

		operator eval_type() const
		{
			eval_type value{};
			accumulate(value);
			return value;
		}

		symphas::internal::iterator_group_difference_type<E, D> ptr;
		iter_type offset;
	};

	template<typename E, size_t D>
	auto operator+(iterator_group_expression<E, D> const& group, typename iterator_group_expression<E, D>::eval_type const& value)
	{
		auto accumulate(value);
		group.accumulate(accumulate);
		return accumulate;
	}

	template<typename E, size_t D>
	auto operator+(typename iterator_group_expression<E, D>::eval_type const& value, iterator_group_expression<E, D> const& group)
	{
		auto accumulate(value);
		group.accumulate(accumulate);
		return accumulate;
	}

	template<typename E, size_t D>
	typename iterator_group_expression<E, D>::eval_type& operator+=(typename iterator_group_expression<E, D>::eval_type& value, iterator_group_expression<E, D> const& iterator)
	{
		value = value + iterator;
		return value;
	}

	template<typename G, size_t D>
	struct iterator_group
	{
		using eval_type = typename symphas::internal::data_value_type<G>::type;

		iterator_group(symphas::internal::iterator_group_difference_type<G, D> const& ptr, iter_type offset = 0) : ptr{ ptr }, offset{ offset } {}

		inline len_type len(grid::region_group<D> const& region) const
		{
			return region.dims[0];
		}

		inline len_type index(grid::region_group<D> const& region) const
		{
			return region[ptr.pos * len(region)];
		}

		template<typename G0>
		auto operator=(iterator_group<G0, D> const& other) const
		{
			auto n0 = index(ptr.region);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				symphas::internal::data_value_type<G>{}(ptr.ptr, n0 + n) = symphas::internal::data_value_type<G0>{}(other.ptr.ptr, n0 + n);
			}
		}

		template<typename E>
		auto operator=(iterator_group_expression<E, D> const& other) const
		{
			auto n0 = index(ptr.region);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				symphas::internal::data_value_type<G>{}(ptr.ptr, n0 + n) = (other.ptr.ptr)->eval(n0 + n);
			}
		}

		template<typename operation_t, typename A, typename B>
		auto operator=(iterator_group_operation<operation_t, A, B> const& other) const
		{
			auto n0 = index(ptr.region);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				symphas::internal::data_value_type<G>{}(ptr.ptr, n0 + n) = other[n0 + n];
			}
		}

		auto operator+(iterator_group<G, D> const& other) const
		{
			eval_type reduce{};
			auto n0 = index(ptr.region, ptr.pos);
			auto n1 = index(other.ptr.region, other.ptr.pos);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				reduce += symphas::internal::data_value_type<G>{}(ptr.ptr, n0 + n) + symphas::internal::data_value_type<G>{}(other.ptr.ptr, n1 + n);
			}
			return reduce;
		}

		void accumulate(eval_type& value) const
		{
			auto n0 = index(ptr.region, ptr.pos);
			for (iter_type n = 0; n < len(ptr.region); ++n)
			{
				value += (ptr.ptr)->eval(n0 + n);
			}
		}


		symphas::internal::iterator_group_difference_type<G, D> ptr;
		iter_type offset;
	};

	template<typename G, size_t D>
	auto operator+(iterator_group<G, D> const& group, typename iterator_group<G, D>::eval_type const& value)
	{
		auto accumulate(value);
		group.accumulate(accumulate);
		return accumulate;
	}

	template<typename G, size_t D>
	auto operator+(typename iterator_group<G, D>::eval_type const& value, iterator_group<G, D> const& group)
	{
		auto accumulate(value);
		group.accumulate(accumulate);
		return accumulate;
	}

	template<typename G, size_t D>
	auto operator+=(typename iterator_group<G, D>::eval_type& value, iterator_group<G, D> const& iterator)
	{
		value = value + iterator;
		return value;
	}



	template<typename G, size_t D>
	iterator_group(symphas::internal::iterator_region_difference_type<G, D>, iter_type = 0) -> iterator_group<G, D>;
	template<typename G, size_t D>
	iterator_group(symphas::internal::iterator_group_difference_type<G, D>, iter_type = 0) -> iterator_group<G, D>;

	template<typename E, size_t D>
	iterator_group_expression(symphas::internal::iterator_region_difference_type<E, D>, iter_type = 0) -> iterator_group_expression<E, D>;




	template<typename operation_t, typename G0, typename G1, size_t D>
	struct iterator_group_operation<operation_t, iterator_group<G0, D>, iterator_group<G1, D>>
	{
		iterator_group_operation(iterator_group<G0, D> const& a, iterator_group<G1, D> const& b) : a{ a }, b{ b } {}

		auto operator[](iter_type n) const
		{
			return operation_t{}(symphas::internal::data_value_type<G0>{}(a.ptr.ptr, n), symphas::internal::data_value_type<G1>{}(b.ptr.ptr, n));
		}

		iterator_group<G0, D> a;
		iterator_group<G1, D> b;
	};

	template<typename operation_t, typename G, typename E, size_t D>
	struct iterator_group_operation<operation_t, iterator_group<G, D>, iterator_group_expression<E, D>>
	{
		iterator_group_operation(iterator_group<G, D> const& a, iterator_group_expression<E, D> const& b) : a{ a }, b{ b } {}

		auto operator[](iter_type n) const
		{
			return operation_t{}(symphas::internal::data_value_type<G>{}(a.ptr.ptr, n), (b.ptr.ptr)->eval(n));
		}

		iterator_group<G, D> a;
		iterator_group_expression<E, D> b;
	};

	template<typename operation_t, typename G, typename E, size_t D>
	struct iterator_group_operation<operation_t, iterator_group_expression<E, D>, iterator_group<G, D>>
	{
		iterator_group_operation(iterator_group_expression<E, D> const& a, iterator_group<G, D> const& b) : a{ a }, b{ b } {}

		auto operator[](iter_type n) const
		{
			return operation_t{}((a.ptr.ptr)->eval(n), symphas::internal::data_value_type<G>{}(b.ptr.ptr, n));
		}

		iterator_group_expression<E, D> a;
		iterator_group<G, D> b;
	};

	template<typename operation_t, typename E0, typename E1, size_t D>
	struct iterator_group_operation<operation_t, iterator_group_expression<E0, D>, iterator_group_expression<E1, D>>
	{
		iterator_group_operation(iterator_group_expression<E0, D> const& a, iterator_group_expression<E1, D> const& b) : a{ a }, b{ b } {}

		auto operator[](iter_type n) const
		{
			return operation_t{}((a.ptr.ptr)->eval(n), (b.ptr.ptr)->eval(n));
		}

		iterator_group_expression<E0, D> a;
		iterator_group_expression<E1, D> b;
	};


	template<typename G0, typename G1, size_t D>
	auto operator+(iterator_group<G0, D> const& a, iterator_group<G1, D> const& b)
	{
		return iterator_group_operation<std::plus<>, iterator_group<G0, D>, iterator_group_expression<G1, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator+(iterator_group<G, D> const& a, iterator_group_expression<E, D> const& b)
	{
		return iterator_group_operation<std::plus<>, iterator_group<G, D>, iterator_group_expression<E, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator+(iterator_group_expression<E, D> const& a, iterator_group<G, D> const& b)
	{
		return iterator_group_operation<std::plus<>, iterator_group_expression<E, D>, iterator_group<G, D>>(a, b);
	}

	template<typename E0, typename E1, size_t D>
	auto operator+(iterator_group_expression<E0, D> const& a, iterator_group_expression<E1, D> const& b)
	{
		return iterator_group_operation<std::plus<>, iterator_group_expression<E0, D>, iterator_group_expression<E1, D>>(a, b);
	}

	template<typename G0, typename G1, size_t D>
	auto operator-(iterator_group<G0, D> const& a, iterator_group<G1, D> const& b)
	{
		return iterator_group_operation<std::minus<>, iterator_group<G0, D>, iterator_group_expression<G1, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator-(iterator_group<G, D> const& a, iterator_group_expression<E, D> const& b)
	{
		return iterator_group_operation<std::minus<>, iterator_group<E, D>, iterator_group_expression<E, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator-(iterator_group_expression<E, D> const& a, iterator_group<G, D> const& b)
	{
		return iterator_group_operation<std::minus<>, iterator_group_expression<E, D>, iterator_group<E, D>>(a, b);
	}

	template<typename E0, typename E1, size_t D>
	auto operator-(iterator_group_expression<E0, D> const& a, iterator_group_expression<E1, D> const& b)
	{
		return iterator_group_operation<std::minus<>, iterator_group_expression<E0, D>, iterator_group_expression<E1, D>>(a, b);
	}

	template<typename G0, typename G1, size_t D>
	auto operator*(iterator_group<G0, D> const& a, iterator_group<G1, D> const& b)
	{
		return iterator_group_operation<std::multiplies<>, iterator_group<G0, D>, iterator_group_expression<G1, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator*(iterator_group<G, D> const& a, iterator_group_expression<E, D> const& b)
	{
		return iterator_group_operation<std::multiplies<>, iterator_group<E, D>, iterator_group_expression<E, D>>(a, b);
	}

	template<typename G, typename E, size_t D>
	auto operator*(iterator_group_expression<E, D> const& a, iterator_group<G, D> const& b)
	{
		return iterator_group_operation<std::multiplies<>, iterator_group_expression<E, D>, iterator_group<E, D>>(a, b);
	}

	template<typename E0, typename E1, size_t D>
	auto operator*(iterator_group_expression<E0, D> const& a, iterator_group_expression<E1, D> const& b)
	{
		return iterator_group_operation<std::multiplies<>, iterator_group_expression<E0, D>, iterator_group_expression<E1, D>>(a, b);
	}
}




