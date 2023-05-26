
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

#include "expressionlogic.h"


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
			explicit constexpr operator int() 
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

		template<typename I0, typename... Is>
		using v_ = OpTerm<OpIdentity, v_id_type<I0, Is...>>;



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
	static char* name = expr::print_with_subscript("v", expr::print_list(expr::symbols::i_<Ns, Ps>{}...)).new_str(); return name;)
ALLOW_COMBINATION((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>))
DEFINE_BASE_DATA((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>), expr::symbols::Symbol{}, expr::symbols::Symbol{})

//! Specialization based on expr::base_data_type.
template<int N0, int P0>
struct expr::base_data_type<expr::symbols::i_<N0, P0>>
{
	using type = expr::symbols::i_<N0, P0>;
};

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
	return expr::symbols::i_op_type<N0, P0>{} * expr::symbols::i_op_type<N1, P1>{};
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




template<typename E0, typename... T0s>
struct SymbolicFunctionArray
{
	using f_type = SymbolicFunction<E0, T0s...>;
	using eval_type = std::invoke_result_t<decltype(&f_type::eval), f_type, iter_type>;

	SymbolicFunctionArray(len_type n = 0, len_type len = 0) :
		data{ (len > 0) ? new f_type * [len] : nullptr },
		offsets{ (len > 0) ? new iter_type * [len] {} : nullptr },
		len{ len }, n{ n }
	{
		for (iter_type i = 0; i < len; ++i)
		{
			offsets[i] = new iter_type[n]{};
			data[i] = nullptr;
		}
	}

	void init(f_type const& f)
	{
		if (data != nullptr)
		{
			if (data[0] == nullptr)
			{
				for (iter_type i = 0; i < len; ++i)
				{
					data[i] = new f_type(f);
				}
			}
		}
	}

	SymbolicFunctionArray(SymbolicFunctionArray const& other) :
		data{ (other.len > 0) ? new f_type * [other.len] : nullptr },
		offsets{ (other.len > 0) ? new iter_type * [other.len] : nullptr },
		len{ other.len }, n{ other.n }
	{
		for (iter_type i = 0; i < len; ++i)
		{
			data[i] = (other.data[i] != nullptr) ? new f_type(*other.data[i]) : nullptr;
			offsets[i] = new iter_type[n];
			std::copy(other.offsets[i], other.offsets[i] + n, offsets[i]);
		}
	}

	SymbolicFunctionArray(SymbolicFunctionArray&& other) : SymbolicFunctionArray()
	{
		swap(*this, other);
	}

	SymbolicFunctionArray<E0, T0s...> operator=(SymbolicFunctionArray<E0, T0s...> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(SymbolicFunctionArray<E0, T0s...>& first, SymbolicFunctionArray<E0, T0s...>& second)
	{
		std::swap(first.data, second.data);
		std::swap(first.offsets, second.offsets);
		std::swap(first.len, second.len);
		std::swap(first.n, second.n);
	}


	const auto& operator[](iter_type i) const
	{
		return *data[i];
	}

	template<size_t N>
	const auto& operator()(DynamicIndex(&index)[N], iter_type i) const
	{
		for (iter_type j = 0; j < n; ++j)
		{
			index[j] = offsets[i][j];
		}
		return *data[i];
	}

	~SymbolicFunctionArray()
	{
		if (data != nullptr)
		{
			for (iter_type i = 0; i < len; ++i)
			{
				delete data[i];
			}
			delete[] data;
		}
	}

	f_type** data;
	iter_type** offsets;
	len_type len;
	len_type n;


};



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
	: SymbolicDataArray<SymbolicData<expr::storage_t<OpAdd<OpTerm<OpIdentity, Ts>...>>>>
{
	using tuple_type = std::tuple<Ts...>;
	using storage_type = expr::storage_t<OpAdd<OpTerm<OpIdentity, Ts>...>>;
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
	: SymbolicDataArray<NamedData<SymbolicData<expr::storage_t<OpAdd<OpTerm<OpIdentity, Ts>...>>>>>
{
	using tuple_type = std::tuple<Ts...>;
	using storage_type = expr::storage_t<OpAdd<OpTerm<OpIdentity, Ts>...>>;
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
		if constexpr (expr::clear_named_data<T0>::value)
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term)), false), term.name);
		}
		else
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term)), false), "");
		}
	}

	template<size_t N, typename T0>
	auto set_data_1(Term<T0> const& term)
	{
		if constexpr (expr::clear_named_data<T0>::value)
		{
			data[N] = NamedData<SymbolicData<storage_type>>(
				SymbolicData<storage_type>(&expr::BaseData<T0>::get(const_cast<T0&>(term.data())), false), term.name);
		}
		else
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

	template<typename I, typename T0, typename List>
	struct select_v_i_impl;


	template<int N0, int P0, typename... Vs>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int N0, int P0, typename... Vs, int... Ns, int... Ps, typename... Rest>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>
	{
		static const bool flag = symphas::lib::index_of_value<int, N0, Ns...> >= 0;

		using type = typename select_v_i_impl<
			expr::symbols::i_<N0, P0>,
			std::conditional_t<flag,
				symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>,
				symphas::lib::types_list<Vs...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int N0, int P0, typename... Vs, int... Ns, int... Ps, size_t D, typename... Rest>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>, Rest...>>
	{
		using type = typename select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>::type;
	};

	template<int N0, int P0, typename... Vs, typename T, typename... Rest>
	struct select_v_i_impl<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_i_impl<
			expr::symbols::i_<N0, P0>,
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

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

	template<typename Is, typename T0, typename List>
	struct select_v_nested_impl;

	template<int... N0s, int... P0s, typename... Vs>
	struct select_v_nested_impl<symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int... N0s, int... P0s, typename... Vs, int... Ns, int... Ps, typename... Rest>
	struct select_v_nested_impl<symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>
	{
		static const bool flag = symphas::lib::filter_seq_t<std::integer_sequence<int, Ns...>, std::integer_sequence<int, N0s...>>::size() == 0;

		using type = typename select_v_nested_impl<
			symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>,
			std::conditional_t<flag,
				symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>,
				symphas::lib::types_list<Vs...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int... N0s, int... P0s, typename... Vs, int... Ns, int... Ps, size_t D, typename... Rest>
	struct select_v_nested_impl<symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, D>, Rest...>>
	{
		using type = typename select_v_nested_impl<symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>::type;
	};

	template<int... N0s, int... P0s, typename... Vs, typename T, typename... Rest>
	struct select_v_nested_impl<symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_nested_impl<
			symphas::lib::types_list<expr::symbols::i_<N0s, P0s>...>,
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename Is, typename List>
	using select_v_nested_ = typename select_v_nested_impl<Is, symphas::lib::types_list<>, List>::type;


	template<typename... Ts>
	struct select_v_impl;
	
	template<typename... Vs, typename I0, typename... Is, typename... Rest>
	struct select_v_impl<symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<I0, Is...>, Rest...>>
	{
		using type = typename select_v_impl<
			symphas::lib::types_list<Vs..., expr::symbols::v_id_type<I0, Is...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Vs, typename I0, typename... Is, size_t D, typename... Rest>
	struct select_v_impl<symphas::lib::types_list<Vs...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<I0, Is...>, D>, Rest...>>
	{
		using type = typename select_v_impl<symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<I0, Is...>, Rest...>>::type;
	};

	template<typename... Vs, typename T, typename... Rest>
	struct select_v_impl<symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_impl<
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Vs>
	struct select_v_impl<symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<typename List>
	using select_v_ = typename select_v_impl<symphas::lib::types_list<>, List>::type;



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



