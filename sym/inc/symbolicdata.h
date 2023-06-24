
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
 * PURPOSE: Implements data used in symbolic constructs.
 *
 * ***************************************************************************
 */

#pragma once

#include "symbolicprototypes.h"
#include "expressions.h"


/*!
 * \defgroup substitutables Symbolic Constructs and Expression Building
 * Symbolic constructs are used to simplify expression definitions and easily substitute data.
 * 
 * The symbolic algebra library defines two types of symbolic algebra groups: Classes with a name
 * beginning in **Op**, which indicates that they are a symbolic algebra expression which can
 * be _evaluated_; and, Classes which begin with **Symbolic**, which are constructs that extend
 * the functionality/usability of symbolic expressions and cannot be evaluated.
 * 
 * The **Symbolic** constructs include an evaluable function (SymbolicFunction), a pattern 
 * or template type that allows substitution of expressions (SymbolicTemplate), and a series 
 * object (SymbolicSeries). 
 * @{
 */

namespace expr
{

	namespace symbols
	{



		template<int N, int P = 0>
		using i_op_type = OpTerm<OpIdentity, i_<N, P>>;

		template<typename I>
		constexpr bool is_index = false;
		template<int N, int P>
		constexpr bool is_index<i_<N, P>> = true;
		template<int N, int P>
		constexpr bool is_index<i_op_type<N, P>> = true;



		template<int N, int P>
		struct i_ : Symbol
		{
			explicit constexpr operator int() const
			{
				return 0;
			}

			constexpr operator i_op_type<N, P>() const
			{
				return {};
			}

			auto operator-() const
			{
				return -i_op_type<N, P>{};
			}

			template<int M, int Q>
			auto operator=(i_<M, Q>) const
			{
				return index_eq<i_<N, P>, i_<M, Q>>{};
			}

			auto operator=(OpVoid) const
			{
				return index_eq_N<i_<N, P>, 0>{};
			}

			auto operator=(OpIdentity) const
			{
				return index_eq_N<i_<N, P>, 1>{};
			}

			auto operator=(OpNegIdentity) const
			{
				return index_eq_N<i_<N, P>, -1>{};
			}

			template<size_t N0>
			auto operator=(OpFractionLiteral<N0, 1>) const
			{
				return index_eq_N<i_<N, P>, N0>{};
			}

			template<size_t N0>
			auto operator=(OpNegFractionLiteral<N0, 1>) const
			{
				return index_eq_N<i_<N, P>, -int(N0)>{};
			}

			auto operator!=(OpVoid) const
			{
				return index_neq_N<i_<N, P>, 0>{};
			}

			auto operator!=(OpIdentity) const
			{
				return index_neq_N<i_<N, P>, 1>{};
			}

			auto operator!=(OpNegIdentity) const
			{
				return index_neq_N<i_<N, P>, -1>{};
			}

			template<size_t N0>
			auto operator!=(OpFractionLiteral<N0, 1>) const
			{
				return index_neq_N<i_<N, P>, N0>{};
			}

			template<size_t N0>
			auto operator!=(OpNegFractionLiteral<N0, 1>) const
			{
				return index_neq_N<i_<N, P>, -int(N0)>{};
			}

			template<int M, int Q>
			auto operator!=(i_<M, Q>) const
			{
				return index_neq<i_<N, P>, i_<M, Q>>{};
			}

			template<int N0, int P0, int M>
			auto operator=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpFractionLiteral<M, 1>>) const
			{
				return index_eq<i_<N, P>, expr::symbols::i_<N0, P0 + M>>{};
			}

			template<int N0, int P0, int M>
			auto operator=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpNegFractionLiteral<M, 1>>) const
			{
				return index_eq<i_<N, P>, expr::symbols::i_<N0, P0 - M>>{};
			}

			template<int N0, int P0>
			auto operator=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpIdentity>) const
			{
				return index_eq<i_<N, P>, expr::symbols::i_<N0, P0 + 1>>{};
			}

			template<int N0, int P0>
			auto operator=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpNegIdentity>) const
			{
				return index_eq<i_<N, P>, expr::symbols::i_<N0, P0 - 1>>{};
			}

			template<typename... E0s>
			auto operator=(OpAdd<E0s...>) const
			{
				return index_eq<i_<N, P>, OpAdd<E0s...>>{};
			}

			template<int N0, int P0, int M>
			auto operator!=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpFractionLiteral<M, 1>>) const
			{
				return index_neq<i_<N, P>, expr::symbols::i_<N0, P0 + M>>{};
			}

			template<int N0, int P0, int M>
			auto operator!=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpNegFractionLiteral<M, 1>>) const
			{
				return index_neq<i_<N, P>, expr::symbols::i_<N0, P0 - M>>{};
			}

			template<int N0, int P0>
			auto operator!=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpIdentity>) const
			{
				return index_neq<i_<N, P>, expr::symbols::i_<N0, P0 + 1>>{};
			}

			template<int N0, int P0>
			auto operator!=(OpAdd<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpNegIdentity>) const
			{
				return index_neq<i_<N, P>, expr::symbols::i_<N0, P0 - 1>>{};
			}

			template<typename... E0s>
			auto operator!=(OpAdd<E0s...>) const
			{
				return index_neq<i_<N, P>, OpAdd<E0s...>>{};
			}

			template<typename E>
			auto operator=(E) const
			{
				return index_eq<i_<N, P>, E>{};
			}

			template<typename E>
			auto operator!=(E) const
			{
				return index_neq<i_<N, P>, E>{};
			}



			template<typename E>
			auto operator>(E) const;
			template<typename E>
			auto operator>=(E) const;
			template<typename E>
			auto operator<(E) const;
			template<typename E>
			auto operator<=(E) const;

		};


		template<int N0, int P0, int... Ns, int... Ps>
		struct v_id_type<i_<N0, P0>, i_<Ns, Ps>...> : Symbol {};



		template<int P0, int N, int P>
		auto make_series_variable(i_<N, P>)
		{
			return v_<i_<N, P0>>{};
		}

	}

}


template<int N, int P>
template<typename E>
auto expr::symbols::i_<N, P>::operator>=(E) const
{
	return this->operator=(E{});
}

template<int N, int P>
template<typename E>
auto expr::symbols::i_<N, P>::operator>(E) const
{
	return this->operator>(E{} + OpIdentity{});
}

template<int N, int P>
template<typename E>
auto expr::symbols::i_<N, P>::operator<=(E) const
{
	return this->operator=(OpIdentity{}) && this->operator=(E{});
}

template<int N, int P>
template<typename E>
auto expr::symbols::i_<N, P>::operator<(E) const
{
	return this->operator<=(E{} - OpIdentity{});
}


template<typename I0, typename T0, typename T1>
auto operator&&(expr::symbols::index_eq<I0, T0>, expr::symbols::index_eq<I0, T1>)
{
	using namespace expr::symbols;
	return symphas::lib::types_list<index_eq<I0, T0>, index_eq<I0, T1>>{};
}

template<typename I0, typename T0, int N1>
auto operator&&(expr::symbols::index_eq<I0, T0>, expr::symbols::index_eq_N<I0, N1>)
{
	using namespace expr::symbols;
	return symphas::lib::types_list<index_eq<I0, T0>, index_eq_N<I0, N1>>{};
}

template<typename I0, int N0, typename T1>
auto operator&&(expr::symbols::index_eq_N<I0, N0>, expr::symbols::index_eq<I0, T1>)
{
	using namespace expr::symbols;
	return symphas::lib::types_list<index_eq_N<I0, N0>, index_eq<I0, T1>>{};
}

template<typename I0, int N0, int N1>
auto operator&&(expr::symbols::index_eq_N<I0, N0>, expr::symbols::index_eq_N<I0, N1>)
{
	using namespace expr::symbols;
	return symphas::lib::types_list<index_eq_N<I0, N0>, index_eq_N<I0, N1>>{};
}


DEFINE_SYMBOL_ID((size_t N), (expr::symbols::arg_t<N>), static char* name = expr::print_with_subscript<N>("arg").new_str(); return name;)
ALLOW_COMBINATION((size_t N), (expr::symbols::arg_t<N>))

DEFINE_SYMBOL_ID((int N, int P), (expr::symbols::i_<N, P>), static char* name = expr::print_with_subscript<N>("i").new_str(); return name;)
ALLOW_COMBINATION((int N, int P), (expr::symbols::i_<N, P>))
DEFINE_BASE_DATA((int N, int P), (expr::symbols::i_<N, P>), expr::symbols::Symbol{}, expr::symbols::Symbol{})

//DEFINE_SYMBOL_ID((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>), static char* name = expr::print_with_subscripts<Ns...>("v").new_str(); return name;)
DEFINE_SYMBOL_ID((int... Ns, int... Ps),
	(expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>),
	static char* name = expr::print_with_subscript("v", print_list(expr::symbols::i_<Ns, Ps>{}...)).new_str(); return name;)
ALLOW_COMBINATION((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>))
DEFINE_BASE_DATA((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>), expr::symbols::Symbol{}, expr::symbols::Symbol{})

//! Specialization based on expr::base_data_type.
template<int N0, int P0>
struct expr::base_data_type<expr::symbols::i_<N0, P0>>
{
	using type = expr::symbols::i_<N0, P0>;
};

namespace expr::symbols
{

	template<int N, int P, typename E>
	auto operator+(expr::symbols::i_<N, P>, E const& e)
	{
		return expr::symbols::i_op_type<N, P>{} + e;
	}

	template<int N, int P, typename E>
	auto operator+(E const& e, expr::symbols::i_<N, P>)
	{
		return e + expr::symbols::i_op_type<N, P>{};
	}

	template<int N0, int P0, int N1, int P1>
	auto operator+(expr::symbols::i_<N0, P0>, expr::symbols::i_<N1, P1>)
	{
		return expr::symbols::i_op_type<N0, P0>{} + expr::symbols::i_op_type<N1, P1>{};
	}

	template<int N, int P, typename E>
	auto operator-(expr::symbols::i_<N, P>, E const& e)
	{
		return expr::symbols::i_op_type<N, P>{} - e;
	}

	template<int N, int P, typename E>
	auto operator-(E const& e, expr::symbols::i_<N, P>)
	{
		return e - expr::symbols::i_op_type<N, P>{};
	}

	template<int N0, int P0, int N1, int P1>
	auto operator-(expr::symbols::i_<N0, P0>, expr::symbols::i_<N1, P1>)
	{
		return expr::symbols::i_op_type<N0, P0>{} - expr::symbols::i_op_type<N1, P1>{};
	}

	template<int N, int P, typename E>
	auto operator*(expr::symbols::i_<N, P>, E const& e)
	{
		return expr::symbols::i_op_type<N, P>{} *e;
	}

	template<int N, int P, typename E>
	auto operator*(E const& e, expr::symbols::i_<N, P>)
	{
		return e * expr::symbols::i_op_type<N, P>{};
	}

	template<int N0, int P0, int N1, int P1>
	auto operator*(expr::symbols::i_<N0, P0>, expr::symbols::i_<N1, P1>)
	{
		return expr::symbols::i_op_type<N0, P0>{} *expr::symbols::i_op_type<N1, P1>{};
	}

	template<int N, int P, typename E>
	auto operator/(expr::symbols::i_<N, P>, E const& e)
	{
		return expr::symbols::i_op_type<N, P>{} / e;
	}

	template<int N, int P, typename E>
	auto operator/(E const& e, expr::symbols::i_<N, P>)
	{
		return e / expr::symbols::i_op_type<N, P>{};
	}

	template<int N0, int P0, int N1, int P1>
	auto operator/(expr::symbols::i_<N0, P0>, expr::symbols::i_<N1, P1>)
	{
		return expr::symbols::i_op_type<N0, P0>{} / expr::symbols::i_op_type<N1, P1>{};
	}

	template<int N, int P>
	auto operator*(expr::symbols::i_<N, P>, OpIdentity)
	{
		return expr::symbols::i_<N, P>{};
	}

	template<int N, int P>
	auto operator*(OpIdentity, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_<N, P>{};
	}


	template<int N, int P>
	auto operator+(expr::symbols::i_<N, P>, OpVoid)
	{
		return expr::symbols::i_op_type<N, P>{};
	}

	template<int N, int P>
	auto operator+(OpVoid, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P>{};
	}

	template<int N, int P>
	auto operator-(expr::symbols::i_<N, P>, OpVoid)
	{
		return expr::symbols::i_op_type<N, P>{};
	}

	template<int N, int P>
	auto operator-(OpVoid, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P>{};
	}

	template<int N, int P>
	auto operator*(expr::symbols::i_<N, P>, OpVoid)
	{
		return OpVoid{};
	}

	template<int N, int P>
	auto operator*(OpVoid, expr::symbols::i_<N, P>)
	{
		return OpVoid{};
	}

	template<int N, int P>
	auto operator/(OpVoid, expr::symbols::i_<N, P>)
	{
		return OpVoid{};
	}


	template<int N, int P>
	auto operator+(expr::symbols::i_<N, P>, OpIdentity)
	{
		return expr::symbols::i_op_type<N, P + 1>{};
	}

	template<int N, int P>
	auto operator+(expr::symbols::i_<N, P>, OpNegIdentity)
	{
		return expr::symbols::i_op_type<N, P - 1>{};
	}

	template<int N, int P, size_t N0>
	auto operator+(expr::symbols::i_<N, P>, OpFractionLiteral<N0, 1>)
	{
		return expr::symbols::i_op_type<N, P + N0>{};
	}

	template<int N, int P, size_t N0>
	auto operator+(expr::symbols::i_<N, P>, OpNegFractionLiteral<N0, 1>)
	{
		return expr::symbols::i_op_type<N, P - N0>{};
	}

	template<int N, int P>
	auto operator-(expr::symbols::i_<N, P>, OpIdentity)
	{
		return expr::symbols::i_op_type<N, P>{} + OpNegIdentity{};
	}

	template<int N, int P>
	auto operator-(expr::symbols::i_<N, P>, OpNegIdentity)
	{
		return expr::symbols::i_op_type<N, P>{} + OpIdentity{};
	}

	template<int N, int P, size_t N0>
	auto operator-(expr::symbols::i_<N, P>, OpFractionLiteral<N0, 1>)
	{
		return expr::symbols::i_op_type<N, P>{} + OpNegFractionLiteral<N0, 1>{};
	}

	template<int N, int P, size_t N0>
	auto operator-(expr::symbols::i_<N, P>, OpNegFractionLiteral<N0, 1>)
	{
		return expr::symbols::i_op_type<N, P>{} + OpFractionLiteral<N0, 1>{};
	}

	template<int N, int P>
	auto operator+(OpIdentity, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P + 1>{};
	}

	template<int N, int P>
	auto operator+(OpNegIdentity, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P - 1>{};
	}

	template<int N, int P, size_t N0>
	auto operator+(OpFractionLiteral<N0, 1>, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P + N0>{};
	}

	template<int N, int P, size_t N0>
	auto operator+(OpNegFractionLiteral<N0, 1>, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, P - N0>{};
	}

	template<int N, int P>
	auto operator-(OpIdentity, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, -P + 1>{};
	}

	template<int N, int P>
	auto operator-(OpNegIdentity, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, -P - 1>{};
	}

	template<int N, int P, size_t N0>
	auto operator-(OpFractionLiteral<N0, 1>, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, -P + N0>{};
	}

	template<int N, int P, size_t N0>
	auto operator-(OpNegFractionLiteral<N0, 1>, expr::symbols::i_<N, P>)
	{
		return expr::symbols::i_op_type<N, -P - N0>{};
	}
}





//! Stores a data that is used symbolically.
/*!
 * Data can be stored symbolically in an expression so that it can be conveniently updated
 * outside of the expression. In other cases, a set value be passed, which will be
 * copied into the SymbolicData object and managed locally.
 *
 * In order to hold a copy, a reference to a data must be passed using std::ref().
 */
template<typename T>
struct SymbolicData
{
	constexpr SymbolicData() : data{ nullptr }, is_local{ false } {}

	//! Takes a copy of the data and manages it locally within the SymbolicData.
	SymbolicData(T&& data) : data{ new T(data) }, is_local{ true } {}

	//! Takes a reference to a data and stores the pointer to it.
	SymbolicData(T& data) : data{ &data }, is_local{ false } {}

	SymbolicData(T* data, bool is_local) : data{ data }, is_local{ is_local } {}

	SymbolicData(SymbolicData<T> const& other) :
		data{ (other.is_local) ? new T(*other.data) : other.data }, is_local{ other.is_local } {}

	SymbolicData(SymbolicData<T>&& other) : SymbolicData()
	{
		swap(*this, other);
	}

	SymbolicData<T> operator=(SymbolicData<T> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(SymbolicData<T>& first, SymbolicData<T>& second)
	{
		using std::swap;
		swap(first.data, second.data);
		swap(first.is_local, second.is_local);
	}

	auto set_data(T data)
	{
		if (is_local)
		{
			delete this->data;
		}

		this->data = new T(data);
		is_local = true;
	}

	auto set_data(symphas::ref<T> data)
	{
		if (is_local)
		{
			delete this->data;
		}

		this->data = &data.get();
		is_local = false;
	}

	auto set_data(symphas::ref<const T> data)
	{
		if (is_local)
		{
			delete this->data;
		}

		this->data = const_cast<T*>(&data.get());
		is_local = false;
	}

	auto set_data(SymbolicData<T> other)
	{
		if (is_local)
		{
			delete data;
		}

		is_local = other.is_local;
		if (is_local)
		{
			data = new T(*other.data);
		}
		else
		{
			data = other.data;
			is_local = false;
		}
	}

	template<typename T0>
	operator SymbolicData<T0>() const
	{
		return SymbolicData<T0>(T0(*data));
	}

	~SymbolicData()
	{
		if (is_local)
		{
			delete data;
		}
	}


	T* data;		//!< The data of this object
	bool is_local;	//!< Indicates whether the data has been allocated within the object
};


//! Stores a list data that is used symbolically.
/*!
 * A list of data can be stored symbolically in an expression so that the expression can operate
 * over multiple instances of a particular data, which can also be updated outside the expression.
 * In other cases, set values can be passed, which will be copied into the SymbolicDataArray
 * object and managed locally.
 *
 * In order to hold a copy, a reference to a data must be passed using std::ref().
 */
template<typename T>
struct SymbolicDataArray : SymbolicData<T>
{
	using parent_type = SymbolicData<T>;
	using parent_type::data;
	using parent_type::is_local;

	constexpr SymbolicDataArray() : parent_type(), len{ 0 } {}

	//! Takes a copy of the data and manages it locally within the SymbolicData.
	SymbolicDataArray(len_type len) : parent_type(get_initialized_array(len), false), len{ len } {}

	//! Takes a copy of the data and manages it locally within the SymbolicData.
	SymbolicDataArray(T* data, len_type len, bool is_local = true) :
		parent_type((is_local) ? get_initialized_array(len) : data, is_local), len{ len }
	{
		if (is_local)
		{
			std::copy(data, data + len, parent_type::data);
		}
	}

	SymbolicDataArray(SymbolicDataArray<T> const& other) :
		parent_type((other.is_local) ? get_initialized_array(other.len) : other.data, other.is_local), len{ other.len }
	{
		if (is_local)
		{
			std::copy(other.data, other.data + len, parent_type::data);
		}
	}

	SymbolicDataArray(SymbolicDataArray<T>&& other) : SymbolicDataArray()
	{
		swap(*this, other);
	}

	SymbolicDataArray<T> operator=(SymbolicDataArray<T> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(SymbolicDataArray<T>& first, SymbolicDataArray<T>& second)
	{
		using std::swap;
		swap(static_cast<SymbolicData<T>&>(first), static_cast<SymbolicData<T>&>(second));
		swap(first.len, second.len);
	}

	auto length() const
	{
		return len;
	}

protected:

	T* get_initialized_array(len_type len)
	{
		T* arr = (T*)std::calloc(len, sizeof(T));
		//for (iter_type i = 0; i < len; ++i)
		//{
		//	arr[i] = T{ 0 };
		//}
		return arr;
	}

public:

	auto set_data(T* data, len_type len = -1)
	{
		if (is_local)
		{
			delete[] data;
		}

		parent_type::data = new T(data);
		std::copy(data, data + len, parent_type::data);
		is_local = true;

		if (len > 0)
		{
			this->len = len;
		}
	}

	auto set_data(symphas::ref<T*> data, len_type len = -1)
	{
		if (is_local)
		{
			delete[] data;
		}

		parent_type::data = &data.get();
		is_local = false;

		if (len > 0)
		{
			this->len = len;
		}
	}

	auto set_data(symphas::ref<const T*> data, len_type len = -1)
	{
		set_data(std::ref(const_cast<T*>(data.get())), len);
	}

	auto set_data(SymbolicDataArray<T> other)
	{
		if (is_local)
		{
			delete[] data;
		}

		is_local = other.is_local;
		len = other.len;
		if (is_local)
		{
			parent_type::data = new T[len];
			std::copy(other.data, other.data + len, parent_type::data);
		}
		else
		{
			data = other.data;
			is_local = false;
		}
	}

	auto operator-() const
	{
		SymbolicDataArray<T> result(len);
		for (iter_type i = 0; i < data.len; ++i)
		{
			result.data[i] = -data[i];
		}
		return result;
	}


	template<typename E, typename T0>
	auto result(SymbolicDataArray<T0> const& other, SymbolicData<T> lhs, SymbolicData<T0> rhs, OpExpression<E> const& e) const
	{
		auto len0 = std::min(other.len, len);
		SymbolicDataArray<expr::eval_type_t<E>> result(len0);

		for (iter_type i = 0; i < len0; ++i)
		{
			lhs.set_data(std::ref(data[i]));
			rhs.set_data(std::ref(other.data[i]));
			expr::result(*static_cast<E const*>(&e), result.data[i]);
		}
		return result;
	}

	~SymbolicDataArray()
	{
		if (is_local)
		{
			free(data);
		}
		is_local = false;
	}

	len_type len;
};

namespace symphas::internal
{

	template<typename T>
	struct to_term_base_impl
	{
		using type = T;
	};

	template<typename T>
	struct to_term_base_impl<Term<T>>
	{
		using type = T;
	};

	template<typename T>
	using to_term_base = typename to_term_base_impl<T>::type;

}

template<typename... Ts>
struct SymbolicDataArray<std::tuple<Ts...>>
	: SymbolicDataArray<SymbolicData<expr::storage_type_t<OpAdd<OpTerm<OpIdentity, Ts>...>>>>
{
	using tuple_type = std::tuple<Ts...>;
	using storage_type = expr::storage_type_t<OpAdd<OpTerm<OpIdentity, Ts>...>>;
	using parent_type = SymbolicDataArray<SymbolicData<storage_type>>;

	using parent_type::parent_type;
	using parent_type::is_local;
	using parent_type::data;
	using parent_type::len;

	SymbolicDataArray(len_type len = sizeof...(Ts)) : parent_type(len) {}

	template<typename... T0s>
	SymbolicDataArray(std::tuple<T0s...> const& list, len_type len = sizeof...(Ts)) : parent_type(len)
	{
		set_data(list, std::make_index_sequence<sizeof...(Ts)>{});
	}

	template<size_t N, size_t Z>
	auto set_data_1(Term<Variable<Z>> const& term)
	{
		data[N] = SymbolicData<storage_type>{};
	}

	template<size_t N, typename T0>
	auto set_data_1(T0 const& term)
	{
		using type_n = symphas::lib::type_at_index<N, Ts...>;
		data[N] = SymbolicData<storage_type>(&expr::BaseData<T0>::get(static_cast<type_n&>(const_cast<T0&>(term))), false);
	}

	template<size_t N, typename T0>
	auto set_data_1(Term<T0> const& term)
	{
		using type_n = symphas::lib::type_at_index<N, Ts...>;
		data[N] = SymbolicData<storage_type>(&expr::BaseData<T0>::get(static_cast<type_n&>(const_cast<T0&>(term.data()))), false);
	}

	template<typename... T0s, size_t... Is>
	auto set_data(std::tuple<T0s...> const& list, std::index_sequence<Is...>)
	{
		(set_data_1<Is>(std::get<Is>(list)), ...);
	}

	template<typename... T0s>
	auto set_data(std::tuple<T0s...> const& list)
	{
		set_data(list, std::make_index_sequence<sizeof...(Ts)>{});
	}

	auto operator-() const
	{
		SymbolicData<storage_type> rhs0;
		auto expr = -expr::make_term(rhs0);

		SymbolicDataArray<tuple_type> result;

		for (iter_type i = 0; i < len; ++i)
		{
			rhs0.set_data(data[i]);
			expr::result(expr, result.data[i]);
		}

		return result;
	}

	template<typename E, typename T0>
	auto result(SymbolicDataArray<T0> const& other, SymbolicData<storage_type> lhs,
		SymbolicData<T0> rhs, OpExpression<E> const& e) const
	{
		auto len0 = std::min(other.len, len);
		SymbolicDataArray<expr::eval_type_t<E>> result(len0);

		for (iter_type i = 0; i < len0; ++i)
		{
			lhs.set_data(data[i]);
			rhs.set_data(other.data[i]);
			expr::result(*static_cast<E const*>(&e), result.data[i]);
		}
		return result;
	}

	template<size_t... Is>
	auto get_data_tuple(std::index_sequence<Is...>) const
	{
		return std::make_tuple(Ts(*data[Is].data)...);
	}

	auto get_data_tuple() const
	{
		return get_data_tuple(std::make_index_sequence<sizeof...(Ts)>{});
	}
};


template<typename... Ts>
struct SymbolicDataArray<std::tuple<NamedData<Ts>...>>
	: SymbolicDataArray<NamedData<SymbolicData<expr::storage_type_t<OpAdd<OpTerm<OpIdentity, Ts>...>>>>>
{
	using tuple_type = std::tuple<Ts...>;
	using storage_type = expr::storage_type_t<OpAdd<OpTerm<OpIdentity, Ts>...>>;
	using parent_type = SymbolicDataArray<NamedData<SymbolicData<storage_type>>>;

	using parent_type::parent_type;
	using parent_type::is_local;
	using parent_type::data;
	using parent_type::len;

	SymbolicDataArray(len_type len = sizeof...(Ts)) : parent_type(len)//, names{} 
	{}

	template<typename... T0s>
	SymbolicDataArray(std::tuple<T0s...> const& list, len_type len = sizeof...(Ts)) : parent_type(len)//, names{}
	{
		//set_names(list, std::make_index_sequence<sizeof...(T0s)>{});
		set_data(list, std::make_index_sequence<sizeof...(Ts)>{});
	}


	SymbolicDataArray(SymbolicDataArray<std::tuple<NamedData<Ts>...>> const& other) :
		parent_type(other.len)//, names{} 
	{
		for (iter_type i = 0; i < len; ++i)
		{
			//names[i] = new char[std::strlen(other.names[i]) + 1];
			//std::strcpy(names[i], other.names[i]);
			data[i] = other.data[i];
		}


	}

	SymbolicDataArray(SymbolicDataArray<std::tuple<NamedData<Ts>...>>&& other) noexcept : SymbolicDataArray()
	{
		swap(*this, other);
	}

	SymbolicDataArray<std::tuple<NamedData<Ts>...>> operator=(SymbolicDataArray<std::tuple<NamedData<Ts>...>> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(SymbolicDataArray<std::tuple<NamedData<Ts>...>>& first,
		SymbolicDataArray<std::tuple<NamedData<Ts>...>>& second)
	{
		using std::swap;
		swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
		//swap(first.names, second.names);
	}

	template<size_t N, typename T0>
	auto set_data_1(T0 const& term)
	{
#ifdef PRINTABLE_EQUATIONS
		if constexpr (expr::clear_named_data<T0>::value)
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term)), false), term.name);
		}
		else
#endif
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term)), false), "");
		}
	}

	template<size_t N, typename T0>
	auto set_data_1(Term<T0> const& term)
	{
#ifdef PRINTABLE_EQUATIONS
		if constexpr (expr::clear_named_data<T0>::value)
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term.data())), false), term.name);
		}
		else
#endif
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term.data())), false), "");
		}
	}

	template<typename... T0s, size_t... Is>
	auto set_data(std::tuple<T0s...> const& list, std::index_sequence<Is...>)
	{
		(set_data_1<Is>(std::get<Is>(list)), ...);
	}

	template<typename... T0s>
	auto set_data(std::tuple<T0s...> const& list)
	{
		set_data(list, std::make_index_sequence<sizeof...(Ts)>{});
	}

	auto operator-() const
	{
		SymbolicData<storage_type> rhs0;
		auto expr = -expr::make_term(rhs0);

		SymbolicDataArray<tuple_type> result;

		for (iter_type i = 0; i < len; ++i)
		{
			rhs0.set_data(data[i]);
			expr::result(expr, result.data[i]);
		}

		return result;
	}

	template<typename E, typename T0>
	auto result(SymbolicDataArray<T0> const& other, SymbolicData<storage_type> lhs,
		SymbolicData<T0> rhs, OpExpression<E> const& e) const
	{
		auto len0 = std::min(other.len, len);
		SymbolicDataArray<expr::eval_type_t<E>> result(len0);

		for (iter_type i = 0; i < len0; ++i)
		{
			lhs.set_data(data[i]);
			rhs.set_data(other.data[i]);
			expr::result(*static_cast<E const*>(&e), result.data[i]);
		}
		return result;
	}

	template<size_t... Is>
	auto get_data_tuple(std::index_sequence<Is...>) const
	{
		return std::make_tuple(Ts(*data[Is].data)...);
	}

	auto get_data_tuple() const
	{
		return get_data_tuple(std::make_index_sequence<sizeof...(Ts)>{});
	}

};

template<typename... T0s>
SymbolicDataArray(std::tuple<T0s...>) -> SymbolicDataArray<
	std::conditional_t<
	(!expr::clear_named_data<T0s>::value && ...),
	std::tuple<symphas::internal::to_term_base<T0s>...>,
	std::tuple<NamedData<symphas::internal::to_term_base<typename expr::clear_named_data<T0s>::type>>...>>>;
//template<size_t... zs, typename... ts>
//symbolicdataarray(std::tuple<term<variable<zs, nameddata<ts>>>...>) 
//	-> symbolicdataarray<std::tuple<term<variable<zs, nameddata<ts>>>...>>;

template<>
struct SymbolicDataArray<expr::symbols::Symbol>
{
	auto operator-() const
	{
		return *this;
	}

	template<typename T0>
	auto operator+(T0 const& other) const
	{
		return *this;
	}

	template<typename T0>
	auto operator-(T0 const& other) const
	{
		return *this;
	}

	template<typename T0>
	auto operator*(T0 const& other) const
	{
		return *this;
	}

	template<typename T0>
	auto operator/(T0 const& other) const
	{
		return *this;
	}
};


template<typename T, typename T0>
auto operator+(SymbolicDataArray<T> const& lhs, SymbolicDataArray<T0> const& rhs)
{
	SymbolicData<T> lhs0;
	SymbolicData<T0> rhs0;
	auto expr = expr::make_term(lhs0) + expr::make_term(rhs0);
	return lhs.result(rhs, lhs0, rhs0, expr);
}

template<typename T, typename T0>
auto operator-(SymbolicDataArray<T> const& lhs, SymbolicDataArray<T0> const& rhs)
{
	SymbolicData<T> lhs0;
	SymbolicData<T0> rhs0;
	auto expr = expr::make_term(lhs0) - expr::make_term(rhs0);
	return lhs.result(rhs, lhs0, rhs0, expr);
}

template<typename T, typename T0>
auto operator*(SymbolicDataArray<T> const& lhs, SymbolicDataArray<T0> const& rhs)
{
	SymbolicData<T> lhs0;
	SymbolicData<T0> rhs0;
	auto expr = expr::make_term(lhs0) * expr::make_term(rhs0);
	return lhs.result(rhs, lhs0, rhs0, expr);
}

template<typename T, typename T0>
auto operator/(SymbolicDataArray<T> const& lhs, SymbolicDataArray<T0> const& rhs)
{
	SymbolicData<T> lhs0;
	SymbolicData<T0> rhs0;
	auto expr = expr::make_term(lhs0) / expr::make_term(rhs0);
	return lhs.result(rhs, lhs0, rhs0, expr);
}


template<typename T, typename E>
auto operator+(SymbolicDataArray<T> const& data, E const& e)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = expr::make_term(el) + op;

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator+(E const& e, SymbolicDataArray<T> const& data)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = op + expr::make_term(el);

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator-(SymbolicDataArray<T> const& data, E const& e)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = expr::make_term(el) - op;

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator-(E const& e, SymbolicDataArray<T> const& data)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = op - expr::make_term(el);

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator*(SymbolicDataArray<T> const& data, E const& e)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = expr::make_term(el) * op;

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator*(E const& e, SymbolicDataArray<T> const& data)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = op * expr::make_term(el);

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator/(SymbolicDataArray<T> const& data, E const& e)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = expr::make_term(el) / op;

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}

template<typename T, typename E>
auto operator/(E const& e, SymbolicDataArray<T> const& data)
{
	SymbolicData<T> el;
	auto op = expr::make_term(e);
	auto expr = op / expr::make_term(el);

	SymbolicDataArray<expr::eval_type_t<decltype(expr)>> result(data.len);

	for (iter_type i = 0; i < data.len; ++i)
	{
		el.set_data(std::ref(data.data[i]));
		expr::result(expr, result.data[i]);
	}
	return result;
}



template<typename T>
auto operator+(SymbolicDataArray<T> const& data, expr::symbols::Symbol)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator+(expr::symbols::Symbol, SymbolicDataArray<T> const& data)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator-(SymbolicDataArray<T> const& data, expr::symbols::Symbol)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator-(expr::symbols::Symbol, SymbolicDataArray<T> const& data)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator*(SymbolicDataArray<T> const& data, expr::symbols::Symbol)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator*(expr::symbols::Symbol, SymbolicDataArray<T> const& data)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator/(SymbolicDataArray<T> const& data, expr::symbols::Symbol)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator/(expr::symbols::Symbol, SymbolicDataArray<T> const& data)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}


template<typename T>
auto operator+(T const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator-(T const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator*(T const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator/(T const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}


template<typename T>
auto operator+(SymbolicDataArray<T> const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator-(SymbolicDataArray<T> const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator*(SymbolicDataArray<T> const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}

template<typename T>
auto operator/(SymbolicDataArray<T> const&, SymbolicDataArray<expr::symbols::Symbol> const&)
{
	return SymbolicDataArray<expr::symbols::Symbol>{};
}


namespace expr
{

	template<typename T>
	auto eval(SymbolicData<T> const& value)
	{
		return expr::eval(*value.data);
	}

	template<typename T>
	auto eval(SymbolicDataArray<T> const& value)
	{
		return expr::eval(*value.data);
	}
}



 //! @}
 
namespace symphas::internal
{


}



