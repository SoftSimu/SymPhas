
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
 * PURPOSE: Implements evaluation routines for symbolic constructs so they
 * are compatible in general symbolic expressions.
 *
 * ***************************************************************************
 */

#pragma once

#include "symbolicseries.h"

namespace symphas::internal
{

    template<typename K, typename E>
	auto swapped_for_dynamic(OpExpression<E> const& e, DynamicIndex const& index) 
	{
		return expr::transform::swap_grid<K>(*static_cast<E const*>(&e), index);
	}

}

template<typename E>
struct SymbolicListIndex<E, void>
{
	E e;

	SymbolicListIndex() = default;
	SymbolicListIndex(E const& e) : e{ e } {}
	SymbolicListIndex(SymbolicListIndex<E> const& e) : SymbolicListIndex(e.e) {}

	iter_type start() const
	{
		return symphas::internal::limit_dimension_start(std::integer_sequence<int>{}, std::make_tuple(), expr::series_limits(e, 0));
	}

	iter_type end() const
	{
		return symphas::internal::limit_dimension_end(std::integer_sequence<int>{}, std::make_tuple(), expr::series_limits(0, e));
	}

	iter_type length() const
	{
		return end() - start() + 1;
	}

};

template<typename E, typename K>
struct SymbolicListIndex : SymbolicListIndex<E, void>
{
	SymbolicListIndex() = default;
	SymbolicListIndex(E const& e) : SymbolicListIndex<E, void>(e) {}
	SymbolicListIndex(E const& e, K const& k) : SymbolicListIndex(e) {}
	SymbolicListIndex(SymbolicListIndex<E> const& e, K const& k) : SymbolicListIndex(e.e) {}

};

template<typename K>
struct SymbolicListIndex<void, K>
{
	template<typename E>
	auto operator()(OpExpression<E> const& e)
	{
		return expr::make_list(K{}, *static_cast<E const*>(&e));
	}
};

template<typename E>
SymbolicListIndex(E) -> SymbolicListIndex<E, void>;
template<typename E, typename K>
SymbolicListIndex(E, K) -> SymbolicListIndex<E, K>;
template<typename E, typename K>
SymbolicListIndex(SymbolicListIndex<E>, K) -> SymbolicListIndex<E, K>;


namespace symphas::internal
{
	template<int N, int P>
	auto indexed_array(expr::symbols::i_<N, P>)
	{
		return SymbolicListIndex<void, expr::symbols::i_<N, P>>{};
	}
}


template<typename V, typename... Ss, typename E, size_t... Ns, typename... Ts>
struct OpSymbolicEval<V, SymbolicSeries<Ss...>, SymbolicFunction<E, Variable<Ns, Ts>...>> :
	OpExpression<OpSymbolicEval<V, SymbolicSeries<Ss...>, SymbolicFunction<E, Variable<Ns, Ts>...>>>
{
    using sub_t = SymbolicSeries<Ss...>;
	using eval_t = SymbolicFunction<E, Variable<Ns, Ts>...>;
	using this_t = OpSymbolicEval<V, sub_t, SymbolicFunction<E, Variable<Ns, Ts>...>>;

	OpSymbolicEval() = default;

	OpSymbolicEval(V value, sub_t const& data, eval_t const& f)
		: value{ value }, f{ f }, data{ data } {}

	V value;
	eval_t f;
	sub_t data;

public:

	auto operator()() const
	{
		auto result = data();
		expr::result(expr::make_term(value, result), result);
		return result;
	}

	auto eval(iter_type n = 0) const
	{
		return expr::eval(value) * data.eval(n);
	}

	auto operator-() const
	{
		return symphas::internal::make_symbolic_eval(-value, data, f);
	}

	void update()
	{
		data.update(f);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::symbolic_eval_print<sub_t>{}(out, data, value * (f.e));
	}

	size_t print(char* out) const
	{
		return expr::symbolic_eval_print<sub_t>{}(out, data, value * (f.e));
	}

	size_t print_length() const
	{
		return expr::symbolic_eval_print<sub_t>{}(data, value * (f.e));
	}

#endif

};

template<typename V, typename sub_t, typename E, size_t... Ns, typename... Ts>
struct OpSymbolicEval<V, sub_t, SymbolicFunction<E, symphas::lib::types_list<Variable<Ns, Ts>...>>> :
	OpSymbolicEval<V, sub_t, SymbolicFunction<E, Variable<Ns, Ts>...>>
{
	using parent_type = OpSymbolicEval<V, sub_t, SymbolicFunction<E, Variable<Ns, Ts>...>>;
	using parent_type::parent_type;
};

template<typename V, int N, int P, typename E, size_t... Ns, typename... Ts>
struct OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Variable<Ns, Ts>...>> :
	OpExpression<OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Variable<Ns, Ts>...>>>
{
	using sub_t = expr::symbols::i_<N, P>;
	using eval_t = SymbolicFunction<E, Variable<Ns, Ts>...>;
	using this_t = OpSymbolicEval<V, sub_t, eval_t>;

	OpSymbolicEval() = default;

	OpSymbolicEval(V value, sub_t const& data, eval_t const& f)
		: value{ value }, f{ f }, data{ data } {}

	V value;
	SymbolicFunction<E, Variable<Ns, Ts>...> f;
	sub_t data;

	auto eval(iter_type n) const
	{
		return expr::eval(value) * f[n];
	}

	auto operator-() const
	{
		return symphas::internal::make_symbolic_eval(-value, data, f);
	}

	void update() {}

	template<typename E0>
	auto operator[](OpExpression<E0> const& index) const
	{
		return value * expr::make_list<expr::symbols::i_<N, P>>(*static_cast<E0 const*>(&index), f);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::symbolic_eval_print<SymbolicListIndex<sub_t>>{}(out, data, value * (f.e));
	}

	size_t print(char* out) const
	{
		return expr::symbolic_eval_print<SymbolicListIndex<sub_t>>{}(out, data, value * (f.e));
	}

	size_t print_length() const
	{
		return expr::symbolic_eval_print<SymbolicListIndex<sub_t>>{}(data, value * (f.e));
	}

#endif

};

template<typename V, typename E0, typename K, typename E, size_t... Ns, typename... Ts>
struct OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Variable<Ns, Ts>...>> :
	OpExpression<OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Variable<Ns, Ts>...>>>
{
	using sub_t = SymbolicListIndex<E0, K>;
	using eval_t = SymbolicFunction<E, Variable<Ns, Ts>...>;
	using this_t = OpSymbolicEval<V, sub_t, eval_t>;

	OpSymbolicEval() : OpSymbolicEval(V{}, sub_t{}, eval_t{}) {}

	OpSymbolicEval(V value, sub_t const& data, eval_t const& f) :
		value{ value }, f{ f }, data{ data },
		index{ data.start(), data.start(), data.end() },
		list{ (data.length() > 0) ? new list_t * [data.length()] {} : nullptr}
	{
		auto sw = symphas::internal::swapped_for_dynamic<K>(f.e, index);
		for (iter_type i = 0; i < data.length(); ++i)
		{
			list[i] = new list_t(sw);
			list[i]->set_data_tuple(f.data);
		}
	}

	OpSymbolicEval(V value, SymbolicListIndex<E0> const& data, eval_t const& f) :
		OpSymbolicEval(value, sub_t(data.e), f) {}

	OpSymbolicEval(OpSymbolicEval<V, sub_t, eval_t> const& other) :
		OpSymbolicEval(other.value, other.data, other.f) {}

	OpSymbolicEval(OpSymbolicEval<V, sub_t, eval_t>&& other) : OpSymbolicEval()
	{
		swap(*this, other);
	}

	OpSymbolicEval<V, sub_t, eval_t>& operator=(OpSymbolicEval<V, sub_t, eval_t> other)
	{
		using std::swap;
		swap(*this, other);
		return *this;
	}

	friend void swap(
		OpSymbolicEval<V, sub_t, eval_t>& first,
		OpSymbolicEval<V, sub_t, eval_t>& second)
	{
		using std::swap;
		swap(first.value, second.value);
		swap(first.f, second.f);
		swap(first.data, second.data);
		swap(first.list, second.list);
		swap(first.index, second.index);
	}

	V value;
	eval_t f;
	sub_t data;
	DynamicIndex index;


	using list_expr_t = std::invoke_result_t<decltype(&symphas::internal::swapped_for_dynamic<K, E>), E, DynamicIndex>;
	using list_t = SymbolicFunction<list_expr_t, Variable<Ns, Ts>...>;

	list_t** list;

	auto eval(iter_type n) const
	{
		int i = (int)data.e.eval(n) - data.start();
		return expr::eval(value) * list[i]->operator[](n);
	}

	auto operator-() const
	{
		return symphas::internal::make_symbolic_eval(-value, data, f);
	}

	void update()
	{
		index = data.e.eval();
		expr::prune::update(list[index.index() - data.start()]->e);
		//for (iter_type i = 0; i < data.length(); ++i)
		//{
		//	index = i + data.start();
		//	expr::prune::update(list[i]->e);
		//}
	}

	//template<int N, int P>
	//auto operator[](expr::symbols::i_<N, P>) const
	//{
	//	return symphas::internal::make_symbolic_eval(value, expr::symbols::i_<N, P>{}, f);
	//}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::symbolic_eval_print<sub_t>{}(out, data, value * (f.e));
	}

	size_t print(char* out) const
	{
		return expr::symbolic_eval_print<sub_t>{}(out, data, value * (f.e));
	}

	size_t print_length() const
	{
		return expr::symbolic_eval_print<sub_t>{}(data, value * (f.e));
	}

#endif

	~OpSymbolicEval()
	{
		if (list != nullptr)
		{
			for (iter_type i = 0; i < data.length(); ++i)
			{
				delete list[i];
			}
			delete[] list;
		}
	}

};

template<typename coeff_t, typename sub_t, typename V, typename... Ts,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpSymbolicEval<V, sub_t, Ts...> const& e)
{
	return symphas::internal::make_symbolic_eval(value * e.value, e.data, e.f);
}

template<typename coeff_t, typename sub_t, typename tensor_t, typename... Ts,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpSymbolicEval<tensor_t, sub_t, Ts...> const& e)
{
	return (value * e.value) * symphas::internal::make_symbolic_eval(OpIdentity{}, e.data, e.f);
}

namespace expr
{
	template<int... I0s, int... P0s, typename E>
	constexpr bool has_selected_index<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, E> =
		(symphas::lib::types_list_size<symphas::internal::select_only_i_<
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, expr::op_types_t<E>>>::value > 0);

	template<int I0, int P0, typename E>
	constexpr bool has_selected_index<expr::symbols::i_<I0, P0>, E> =
		(symphas::lib::types_list_size<symphas::internal::select_only_i_<
			symphas::lib::types_list<expr::symbols::i_<I0, P0>>, expr::op_types_t<E>>>::value > 0);
}


namespace symphas::internal
{
	template<typename E>
	constexpr bool has_index_impl<symphas::lib::types_list<E>> = 
		(symphas::lib::types_list_size<symphas::internal::select_unique_i_<expr::op_types_t<E>>>::value > 0);


	template<size_t N, typename T1, typename T2, typename... T1s, typename... T2s>
	auto update_limits(
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		return std::tuple_cat(
			symphas::lib::get_tuple_lt<N>(limits),
			std::make_tuple(limit0),
			symphas::lib::get_tuple_ge<N + 1>(limits));
	}

	////////// helpers.

	template<int I0, typename V, int I1, int P1>
	auto separate_index(OpTerm<V, expr::symbols::i_<I1, P1>> const& e)
	{
		if constexpr (I0 == I1)
		{
			return std::make_pair(-expr::coeff(e), expr::coeff(e) * val<P1>);
		}
		else
		{
			return std::make_pair(OpVoid{}, e);
		}
	}

	template<int I0, typename T/*, typename = std::enable_if_t<!has_index<T>, int>*/>
	auto separate_index(T const& e)
	{
		return std::make_pair(OpVoid{}, e);
	}

	template<typename... Ls, typename... Rs>
	auto combine_adds(std::pair<Ls, Rs> const& ...pairs)
	{
		return std::make_pair((std::get<0>(pairs) + ...), (std::get<1>(pairs) + ...));
	}

	template<int I0, typename... Es, size_t... Is>
	auto separate_index(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return combine_adds(separate_index<I0>(expr::get<Is>(e))...);
	}

	template<int I0, typename... Es>
	auto separate_index(OpAdd<Es...> const& e)
	{
		return separate_index<I0>(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<int I0, int P0, typename E>
	auto separate_index(expr::symbols::i_<I0, P0>, E const& e)
	{
		return separate_index<I0>(e);
	}


	template<typename T, typename T0>
	struct filter_placeholder_impl;
	

	template<size_t... Ns>
	struct filter_placeholder_impl<
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>...>,
		symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>...>;
	};

	template<size_t... Ns, size_t N0, typename... Rest>
	struct filter_placeholder_impl<
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>...>,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>, Rest...>>
	{
		using type = typename filter_placeholder_impl<
			symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>..., 
				expr::symbols::placeholder_N_symbol_<N0>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<size_t... Ns, typename T, typename... Rest>
	struct filter_placeholder_impl<
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>...>,
		symphas::lib::types_list<T, Rest...>>
	{
		using type = typename filter_placeholder_impl<
			symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<Ns>...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename T>
	using filter_placeholder_t = typename filter_placeholder_impl<symphas::lib::types_list<>, T>::type;

	template<size_t N0, size_t N1, typename V>
	auto separate_placeholder(OpTerm<V, expr::symbols::placeholder_N_symbol_<N1>> const& e)
	{
		if constexpr (N0 == N1)
		{
			return std::make_pair(-expr::coeff(e), OpVoid{});
		}
		else
		{
			return std::make_pair(OpVoid{}, e);
		}
	}

	template<size_t N0>
	auto separate_placeholder(DynamicIndex const& index)
	{
		return std::make_pair(OpVoid{}, index);
	}

	template<size_t N0>
	auto separate_placeholder(OpVoid)
	{
		return std::make_pair(OpVoid{}, OpVoid{});
	}

	template<size_t N0, typename T, typename = std::enable_if_t<(expr::is_coeff<T>), int>>
	auto separate_placeholder(T const& e)
	{
		return std::make_pair(OpVoid{}, e);
	}

	template<size_t N0, typename... Es, size_t... Is>
	auto separate_placeholder(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return combine_adds(separate_placeholder<N0>(expr::get<Is>(e))...);
	}

	template<size_t N0, typename... Es>
	auto separate_placeholder(OpAdd<Es...> const& e)
	{
		return separate_placeholder<N0>(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<size_t N0, typename E>
	auto separate_placeholder(expr::symbols::placeholder_N_symbol_<N0>, E const& e)
	{
		return separate_placeholder<N0>(e);
	}

	template<typename E>
	auto separate_placeholder(symphas::lib::types_list<>, E const& e)
	{
		return std::make_pair(OpVoid{}, e);
	}

	template<size_t N0, typename E>
	auto separate_placeholder(symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>, E const& e)
	{
		return separate_placeholder<N0>(e);
	}

	template<typename E>
	auto separate_placeholder(E const& e)
	{
		return separate_placeholder(filter_placeholder_t<expr::op_types_t<E>>{}, e);
	}

	template<typename T, typename R>
	struct index_id_match;

	template<int... I0s, int... P0s, int... I00s, int... P00s>
	struct index_id_match<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, symphas::lib::types_list<expr::symbols::i_<I00s, P00s>...>>
	{
		using mask_t = std::integer_sequence<bool, symphas::lib::is_value_in_seq<int, I0s, std::integer_sequence<int, I00s...>>::value...>;
		using type = symphas::lib::filter_types<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, mask_t>;
	};

	template<typename E, typename Is>
	struct minor_index;

	template<int... I1s, int... P1s>
	struct minor_index<symphas::lib::types_list<>,
		symphas::lib::types_list<expr::symbols::i_<I1s, P1s>...>>
	{
		using type = OpVoid;
		static const int value = -1;
	};

	template<int I0, int P0, int... I0s, int... P0s, int... I1s, int... P1s>
	struct minor_index<symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::symbols::i_<I1s, P1s>...>>
	{
		using type = symphas::lib::type_at_index<sizeof...(I0s), expr::symbols::i_<I0, P0>, expr::symbols::i_<I0s, P0s>...>;
		static const int value = symphas::lib::index_of_type<type, expr::symbols::i_<I1s, P1s>...>;
	};

	template<typename E, int... I0s, int... P0s>
	struct minor_index<void, symphas::lib::types_list<E, expr::symbols::i_<I0s, P0s>...>>
	{
	protected:

		using index_list = symphas::internal::select_all_i_<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, expr::op_types_t<E>>;
		using in_list_t = typename index_id_match<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, index_list>::type;

	public:

		using type = typename minor_index<in_list_t, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::type;
		static const int value = minor_index<in_list_t, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
	};

	template<int I0, int P0, int... I0s, int... P0s>
	struct minor_index<void, symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I0s, P0s>...>>
	{

	public:

		using type = expr::symbols::i_<I0, P0>;
		static const int value = symphas::lib::index_of_value<int, I0, I0s...>;
	};

	template<typename E, typename Is>
	using minor_index_t = typename minor_index<E, Is>::type;


	////////// update upper limits by propagating a lower limit from the initial index
	// propagating a lower limit updates the upper limits for the more
	// major indices.


	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T2, typename E,
			typename... T1s, typename... T2s, typename = std::enable_if_t<
		(!has_index<E> && !has_selected_index<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, T2>), int>>
		auto propagate_lower_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			expr::symbols::i_<I1, P1>, E const& lower,
			expr::series_limits<T1, T2> const& limit0,
			std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s,
        typename = std::enable_if_t<
			(has_index<E> || has_selected_index<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, T2>), int>>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& lower,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<typename offset_t, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		offset_t const& offset,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T2, typename E,
		typename... T1s, typename... T2s, typename = std::enable_if_t<
			(!has_index<E> && !has_selected_index<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, T1>), int>>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s,
		typename = std::enable_if_t<
			(has_index<E> || has_selected_index<symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, T1>), int>>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<typename offset_t, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		offset_t const& offset,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, std::tuple<expr::series_limits<T1s, T2s>...> const& limits);



	template<bool flag_left, typename coeff_t, typename La_t, typename Ra_t>
	auto compute_limit_side_b(
		coeff_t const& coeff, La_t const& la, Ra_t const& ra)
	{
		bool side = (coeff.eval() > 0);
		if constexpr (flag_left)
		{
			return (side) ? expr::make_literal(la.eval()) : expr::make_literal(ra.eval());
		}
		else
		{
			return (side) ? expr::make_literal(ra.eval()) : expr::make_literal(la.eval());
		}
	}

	template<bool flag_left, typename V>
	auto compute_limit_side_a(
		OpVoid, V const& la, V const& ra)
	{
		if constexpr (flag_left)
		{
			return la;
		}
		else
		{
			return ra;
		}
	}

	template<bool flag_left, typename coeff_t, typename V>
	auto compute_limit_side_a(
		coeff_t const& coeff, OpVoid const& la, V const& ra)
	{
		if constexpr (flag_left)
		{
			return ra;
		}
		else
		{
			return la;
		}
	}

	template<bool flag_left, typename coeff_t, typename V>
	auto compute_limit_side_a(
		coeff_t const& coeff, V const& la, OpVoid const& ra)
	{
		if constexpr (flag_left)
		{
			return la;
		}
		else
		{
			return ra;
		}
	}

	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>);
	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>,
		symphas::lib::types_list<>);
	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>);
	template<bool flag_left, typename L, typename R>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<>);


	template<bool flag_left, typename L, typename R, typename Other0, typename Other1>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other0>,
		symphas::lib::types_list<Other1>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other1>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<Other1>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other0>,
		symphas::lib::types_list<>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other00, typename Other01,
		typename... Other0s, typename Other10, typename Other11, typename... Other1s>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other00, Other01, Other0s...>,
		symphas::lib::types_list<Other10, Other11, Other1s...>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other00, 
		typename Other10, typename Other11, typename... Other1s>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other00>,
		symphas::lib::types_list<Other10, Other11, Other1s...>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other00, typename Other01,
		typename... Other0s, typename Other10>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other00, Other01, Other0s...>,
		symphas::lib::types_list<Other10>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, 
		typename Other10, typename Other11, typename... Other1s>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<Other10, Other11, Other1s...>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, typename Other00, typename Other01,
		typename... Other0s>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<Other00, Other01, Other0s...>,
		symphas::lib::types_list<>)
	{
		return std::make_pair(left, right);
	}

	

	template<bool flag_left, typename L, typename R>
	auto compute_limit_side(L const& left, R const& right);
	template<bool flag_left, typename R>
	auto compute_limit_side(int left, R const& right);
	template<bool flag_left, typename L>
	auto compute_limit_side(L const& left, int right);
	template<bool flag_left, typename A, typename B, typename R>
	auto compute_limit_side(std::pair<A, B> const& left, R const& right);
	template<bool flag_left, typename L, typename A, typename B>
	auto compute_limit_side(L const& left, std::pair<A, B> const& right);




	template<bool flag_left, typename L, typename A, typename B>
	auto _compute_limit_side(L const& left, std::pair<A, B> const& right)
	{
		return right;
	}

	template<bool flag_left, typename A, typename B, typename R>
	auto _compute_limit_side(std::pair<A, B> const& left, R const& right)
	{
		return left;
	}

	template<bool flag_left, typename L, typename R>
	auto _compute_limit_side(L const& left, R const& right)
	{
		return compute_limit_side<flag_left>(left, right, expr::op_types_t<L>{}, expr::op_types_t<R>{});
	}

	template<bool flag_left, typename L, typename R>
	auto compute_limit_side(L const& left, R const& right)
	{
		return compute_limit_side<flag_left>(left, right, expr::op_types_t<L>{}, expr::op_types_t<R>{});
	}

	template<bool flag_left, typename R>
	auto compute_limit_side(int left, R const& right)
	{
		if constexpr (std::is_same_v<expr::symbols::Symbol, expr::eval_type_t<R>>)
		{
			return std::make_pair(left, right);
		}
		else
		{
			return right;
		}
	}

	template<bool flag_left, typename L>
	auto compute_limit_side(L const& left, int right)
	{
		if constexpr (std::is_same_v<expr::symbols::Symbol, expr::eval_type_t<L>>)
		{
			return std::make_pair(left, right);
		}
		else
		{
			return left;
		}
	}

	template<bool flag_left>
	auto compute_limit_side(int left, int right)
	{
		if constexpr (flag_left)
		{
			return std::max(left, right);
		}
		else
		{
			return std::min(left, right);
		}
	}

	template<bool flag_left, typename A, typename B, typename R>
	auto compute_limit_side(std::pair<A, B> const& left, R const& right)
	{
		auto [a, b] = left;
		auto lefta = compute_limit_side<flag_left>(a, right);
		auto leftb = compute_limit_side<flag_left>(b, right);
		return _compute_limit_side<flag_left>(lefta, leftb);
	}

	template<bool flag_left, typename L, typename A, typename B>
	auto compute_limit_side(L const& left, std::pair<A, B> const& right)
	{
		auto [a, b] = right;
		auto righta = compute_limit_side<flag_left>(left, a);
		auto rightb = compute_limit_side<flag_left>(left, b);
		return _compute_limit_side<flag_left>(righta, rightb);
	}

	template<bool flag_left, typename A, typename B>
	auto compute_limit_side(std::pair<A, B> const& left, int right)
	{
		auto [a, b] = left;
		auto lefta = compute_limit_side<flag_left>(a, right);
		auto leftb = compute_limit_side<flag_left>(b, right);
		return _compute_limit_side<flag_left>(lefta, leftb);
	}

	template<bool flag_left, typename A, typename B>
	auto compute_limit_side(int left, std::pair<A, B> const& right)
	{
		auto [a, b] = right;
		auto righta = compute_limit_side<flag_left>(left, a);
		auto rightb = compute_limit_side<flag_left>(left, b);
		return _compute_limit_side<flag_left>(righta, rightb);
	}

	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right, 
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>, 
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>)
	{
		auto [la, lb] = separate_placeholder(expr::symbols::placeholder_N_symbol_<N0>{}, left);
		auto [ra, rb] = separate_placeholder(expr::symbols::placeholder_N_symbol_<N0>{}, right);

		auto da = la - ra;
		auto db = lb - rb;

		auto a = compute_limit_side_a<flag_left>(-da, -la, -ra);
		auto b = compute_limit_side_b<flag_left>(db, lb, rb);
		return a * expr::symbols::placeholder_N_<N0>{} + b;
	}

	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>,
		symphas::lib::types_list<>)
	{
		return std::make_pair(left, right);
	}

	template<bool flag_left, typename L, typename R, size_t N0>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<expr::symbols::placeholder_N_symbol_<N0>>)
	{
		return std::make_pair(right, left);
	}

	template<bool flag_left, typename L, typename R>
	auto compute_limit_side(L const& left, R const& right,
		symphas::lib::types_list<>,
		symphas::lib::types_list<>)
	{
		if constexpr (flag_left)
		{
			return expr::make_literal(std::max((int)expr::eval(left), (int)expr::eval(right)));
		}
		else
		{
			return expr::make_literal(std::min((int)expr::eval(left), (int)expr::eval(right)));
		}
	}

	template<typename L, typename V>
	auto make_limit(L const& min, std::pair<OpLiteral<V>, int> const& max)
	{
		return expr::series_limits(min, std::min(int(max.first), max.second));
	}

	template<typename LA, typename LB, typename V>
	auto make_limit(std::pair<LA, LB> const& min, std::pair<OpLiteral<V>, int> const& max)
	{
		auto [la, lb] = min;
		return expr::series_limits(
			compute_limit_side<true>(la, lb),
			std::min(int(max.first), max.second));
	}

	template<typename L, typename A, typename B>
	auto make_limit(L const& min, std::pair<A, B> const& max)
	{
		auto [a, b] = max;
		return expr::series_limits(min, 
			compute_limit_side<false>(a, b));
	}

	template<typename V, typename R>
	auto make_limit(std::pair<OpLiteral<V>, int> const& min, R const& max)
	{
		return expr::series_limits(std::max(int(min.first), min.second), max);
	}

	template<typename V, typename RA, typename RB>
	auto make_limit(std::pair<OpLiteral<V>, int> const& min, std::pair<RA, RB> const& max)
	{
		auto [ra, rb] = max;
		return expr::series_limits(
			std::max(int(min.first), min.second),
			compute_limit_side<true>(ra, rb));
	}

	template<typename A, typename B, typename R>
	auto make_limit(std::pair<A, B> const& min, R const& max)
	{
		auto [a, b] = min;
		return expr::series_limits(compute_limit_side<true>(a, b), max);
	}

	template<typename LA, typename LB, typename RA, typename RB>
	auto make_limit(std::pair<LA, LB> const& min, std::pair<RA, RB> const& max)
	{
		auto [la, lb] = min;
		auto [ra, rb] = max;
		return expr::series_limits(
			compute_limit_side<true>(la, lb),
			compute_limit_side<true>(ra, rb));
	}

	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T2, typename E,
		typename... T1s, typename... T2s, typename>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& lower,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		return update_limits<size_t(N)>(
			make_limit(expr::limit_0(limit0), std::make_pair(lower, expr::limit_1(limit0))),
			limits);
	}

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s, typename>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& lower,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = minor_index<void, symphas::lib::types_list<E, expr::symbols::i_<I0s, P0s>...>>::value;
		constexpr int L = symphas::lib::index_of_value<int, I1, I0s...>;
		
		// The most minor index of the expression of the updated lower bound of an upstream index
		//	contains the current index. The current index is isolated and the propagation is resumed.
		if constexpr (N == L)
		{
			// we have that lhs <= rhs, so the upper limit is rhs
			auto [lhs, rhs] = separate_index<I1>(lower - expr::symbols::i_<I1, P1>{});
			return propagate_lower_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				expr::symbols::i_<I1, P1>{}, rhs / lhs,
				limit0, limits);
		}
		// Triggered when the given lower limit includes an index. For an updated lower limit that
		//	carries an index, the current index (position `L`) upper bound is updated and set immediately
		else if constexpr (N >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(N), expr::symbols::i_<I0s, P0s>...>;

			//auto limit1 = std::get<size_t(L)>(limits);
			auto updated_limits = update_limits<size_t(L)>(
				make_limit(expr::limit_0(limit0), std::make_pair(lower, expr::limit_1(limit0))), limits);

			return propagate_lower_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, lower - expr::limit_0(limit0) + index_t{},
				std::get<N>(updated_limits), updated_limits);

		}
		else
		{
			// Look for an index in the upper value of the current limit.
			// The upper value of the current limit is the one being updated after the lower limit 
			//	of a more minor index was fixed. The current upper limit defining another index
			//	means that the updated lower must be propogated to it.
			constexpr int L0 = minor_index<void, symphas::lib::types_list<T2, expr::symbols::i_<I0s, P0s>...>>::value;
			if constexpr (L0 >= 0)
			{
				using index_t = symphas::lib::type_at_index<size_t(L0), expr::symbols::i_<I0s, P0s>...>;
				auto [lhs, rhs] = separate_index(index_t{}, lower - expr::limit_1(limit0) + index_t{});

				return propagate_lower_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					index_t{}, rhs / lhs,
					std::get<L0>(limits), limits);
			}
			else
			{
				return propagate_lower_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					expr::symbols::i_<I1, P1>{}, lower,
					limit0, limits);
			}
		}
	}


	// entry point
	template<typename offset_t, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		offset_t const& offset,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, int start, int end,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		auto limit0 = std::get<size_t(N)>(limits);

		auto offset_start = symphas::internal::compute_offset(
			std::integer_sequence<int, I0s...>{}, limits, limit0);

		auto offset_delta = offset - offset_start;
		auto index_value = expr::make_literal(start) + offset_delta + offset_start;

		// update the lower index for this starting index
		auto updated_lower = propagate_upper_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			expr::symbols::i_<I1, P1>{}, index_value,
			//limit0,
			make_limit(std::make_pair(index_value, expr::limit_0(limit0)), expr::limit_1(limit0)),
			limits);

		constexpr int L = minor_index<void, symphas::lib::types_list<decltype(expr::limit_0(limit0)), expr::symbols::i_<I0s, P0s>...>>::value;
		if constexpr (L >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;
			auto [lhs, rhs] = separate_index(index_t{}, expr::limit_0(limit0) - index_value);

			return propagate_lower_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, rhs / lhs,
				std::get<L>(updated_lower), updated_lower);
		}
		else
		{
			return updated_lower;
		}
	}

	inline auto to_dynamic_index(DynamicIndex const& index)
	{
		return index;
	}

	template<typename V, size_t N>
	auto to_dynamic_index(OpTerm<V, expr::symbols::placeholder_N_symbol_<N>> const& value)
	{
		return value;
	}

	template<typename E>
	auto to_dynamic_index(E const& value)
	{
		return DynamicIndex(expr::eval(value), 0);
	}

	template<typename... Es, size_t... Is>
	auto to_dynamic_index(OpAdd<Es...> const& add, std::index_sequence<Is...>)
	{
		return (to_dynamic_index(expr::get<Is>(add)) + ...);
	}

	template<typename... Es>
	auto to_dynamic_index(OpAdd<Es...> const& add)
	{
		return to_dynamic_index(add, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T2, typename E,
		typename... T1s, typename... T2s, typename>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;

		int start = symphas::internal::limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<N>(limits));
		//auto offset_start = symphas::internal::compute_offset(std::integer_sequence<int, I0s...>{}, limits, std::get<N>(limits));
		auto offset_delta = upper - start;// (start + offset_start);

		return update_limits<size_t(N)>(
			make_limit(std::make_pair(start + to_dynamic_index(offset_delta), expr::limit_0(limit0)), expr::limit_1(limit0)),
			limits);
	}

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s, typename>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = minor_index<void, symphas::lib::types_list<E, expr::symbols::i_<I0s, P0s>...>>::value;
		constexpr int L = symphas::lib::index_of_value<int, I1, I0s...>;

		if constexpr (N == L)
		{
			// we have that lhs <= rhs, so the upper limit is rhs
			auto [lhs, rhs] = separate_index<I1>(upper - expr::symbols::i_<I1, P1>{});
			return propagate_upper_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				expr::symbols::i_<I1, P1>{}, rhs / lhs,
				limit0, limits);
		}
		else if constexpr (N >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(N), expr::symbols::i_<I0s, P0s>...>;

			int start = symphas::internal::limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<N>(limits));
			//auto offset_start = symphas::internal::compute_offset(std::integer_sequence<int, I0s...>{}, limits, std::get<N>(limits));
			auto offset_delta = upper - start;// (start + offset_start);

			//auto limit1 = std::get<size_t(L)>(limits);
			auto updated_limits = update_limits<size_t(L)>(
				make_limit(std::make_pair(start + to_dynamic_index(offset_delta), expr::limit_0(limit0)), expr::limit_1(limit0)), limits);


			return propagate_upper_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, upper - expr::limit_1(limit0) + index_t{},
				std::get<N>(updated_limits), updated_limits);

		}
		else
		{
			constexpr int L0 = minor_index<void, symphas::lib::types_list<T1, expr::symbols::i_<I0s, P0s>...>>::value;
			if constexpr (L0 >= 0)
			{
				using index_t = symphas::lib::type_at_index<size_t(L0), expr::symbols::i_<I0s, P0s>...>;
				auto [lhs, rhs] = separate_index(index_t{}, upper - expr::limit_0(limit0) + index_t{});

				return propagate_upper_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					index_t{}, rhs / lhs,
					std::get<L0>(limits), limits);
			}
			else
			{
				return propagate_upper_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					expr::symbols::i_<I1, P1>{}, upper,
					limit0, limits);
			}
		}
	}

	// entry point
	template<typename offset_t, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		offset_t const& offset,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, int start, int end,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		auto limit0 = std::get<size_t(N)>(limits);

		//auto offset_start = symphas::internal::compute_offset(std::integer_sequence<int, I0s...>{}, limits, limit0);
		//auto offset_delta = offset - offset_start;
		//auto index_value = expr::make_literal(start) + offset_delta + offset_start;
		auto index_value = start + offset;

		// update the upper index for this starting index
		auto updated_upper = propagate_lower_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			expr::symbols::i_<I1, P1>{}, index_value,
			//limit0,
			make_limit(expr::limit_0(limit0), std::make_pair(index_value, expr::limit_1(limit0))),
			limits);

		constexpr int L = minor_index<void, symphas::lib::types_list<decltype(expr::limit_1(limit0)), expr::symbols::i_<I0s, P0s>...>>::value;
		if constexpr (L >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;

			auto [lhs, rhs] = separate_index(index_t{}, expr::limit_1(limit0) - index_value /*-expr::limit_1(limit0)*/);

			//auto [lhs_limit, rhs_limit] = separate_index(index_t{}, -expr::limit_1(limit0));
			//auto [lhs_plhlr, rhs_plhlr] = separate_placeholder(rhs_limit);

			//auto rhs_corrected = rhs_limit - rhs_plhlr
			//	+ DynamicIndex(int(rhs_plhlr)) 
			//	+ (start) - offset_start + offset;

			return propagate_upper_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, rhs / lhs,
				std::get<L>(updated_upper), updated_upper);
		}
		else
		{
			return updated_upper;
		}
	}


	////////// fix a specific offset and find new limits that satisfy it


	template<size_t N, typename offset_t, int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		offset_t const& offset,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		using current_ind_t = symphas::lib::type_at_index<N, expr::symbols::i_<I0s, P0s>...>;

		len_type offsets[sizeof...(I0s)]{};
		auto limit0 = std::get<size_t(N)>(limits);
		int start = symphas::internal::limit_start(std::integer_sequence<int, I0s...>{}, limits, limit0/*, offsets*/);
		int end = symphas::internal::limit_end(std::integer_sequence<int, I0s...>{}, limits, limit0/*, offsets*/);

		auto lower_limits = propagate_lower_limit(offset,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			current_ind_t{}, start, end, limits);
		auto upper_limits = propagate_upper_limit(offset,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			current_ind_t{}, start, end, lower_limits);

		return upper_limits;
	}



	template<int I0, int P0, int offset, int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		expr::symbols::index_eq_N<expr::symbols::i_<I0, P0>, offset>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		constexpr int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N >= 0)
		{
			return limits_at_offset<size_t(N)>(
				expr::val<offset> - val<P0>, limits,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
		}
		else
		{
			return limits;
		}
	}

	template<int I0, int P0, int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		expr::symbols::index_eq<expr::symbols::i_<I0, P0>, expr::symbols::placeholder_N>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		constexpr int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N >= 0)
		{
			return limits_at_offset<size_t(N)>(
				expr::symbols::placeholder_N{} - val<P0>, limits,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
		}
		else
		{
			return limits;
		}
	}

	template<int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		symphas::lib::types_list<>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		return limits;
	}

	template<int I00, int P00, int offset, int... I00s, int... P00s, int... offsets,
		int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<I00, P00>, offset>,
			expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		auto fixed_limits = limits_at_offset(
			expr::symbols::index_eq_N<expr::symbols::i_<I00, P00>, offset>{}, 
			limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
		return limits_at_offset(
			symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...>{},
			fixed_limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
	}
	

	template<int I00, int P00, int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		symphas::lib::types_list<expr::symbols::index_eq<expr::symbols::i_<I00, P00>, expr::symbols::placeholder_N>>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		return limits_at_offset(
			expr::symbols::index_eq<expr::symbols::i_<I00, P00>, expr::symbols::placeholder_N>{},
			limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
		//return limits_at_offset(
		//	symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...>{},
		//	fixed_limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
	}

}

namespace expr
{
	namespace
	{

		template<typename... Ts>
		struct filter_grid_symbol_impl;

		template<typename... Ts>
		struct filter_grid_symbol_impl<symphas::lib::types_list<Ts...>, symphas::lib::types_list<>>
		{
			using type = symphas::lib::types_list<Ts...>;
		};

		template<typename... Ts, typename T, typename... Rest>
		struct filter_grid_symbol_impl<symphas::lib::types_list<Ts...>, symphas::lib::types_list<T, Rest...>>
		{
			using type = typename filter_grid_symbol_impl<symphas::lib::types_list<Ts...>, symphas::lib::types_list<Rest...>>::type;
		};

		template<typename... Ts, int N0, int P0, size_t D, typename... Rest>
		struct filter_grid_symbol_impl<symphas::lib::types_list<Ts...>, symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<N0, P0>>, D>, Rest...>>
		{
			using type = typename filter_grid_symbol_impl<symphas::lib::types_list<Ts..., GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<N0, P0>>, D>>, symphas::lib::types_list<Rest...>>::type;
		};

		template<typename T>
		using filter_grid_symbol_t = typename filter_grid_symbol_impl<symphas::lib::types_list<>, T>::type;
		
		template<typename V0, size_t D0, typename... Vs, size_t... Ds, typename E0>
		auto with_vs(
			OpExpression<E0> const& e0,
			symphas::lib::types_list<GridSymbol<V0, D0>, GridSymbol<Vs, Ds>...>)
		{
			return expr::transform::swap_grid<GridSymbol<V0, D0>, GridSymbol<Vs, Ds>...>(*static_cast<E0 const*>(&e0), V0{}, Vs{}...);
		}

		template<typename E0>
		auto with_vs(
			OpExpression<E0> const& e0,
			symphas::lib::types_list<>)
		{
			return *static_cast<E0 const*>(&e0);
		}

		template<typename E0, typename... T1s, typename... T2s, typename TT0, typename... TTs,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			OpExpression<E0> const& e0,
			std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<E,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
				symphas::lib::types_list<limit_ts...>,
				symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series,
			TT0&& data0, TTs&&... datas)
		{
			return expr::series<Op, expr::symbols::i_<I0s, P0s>...>(
				with_vs(*static_cast<E0 const*>(&e0), filter_grid_symbol_t<op_types_t<E0>>{})
			)(limits)(std::forward<TT0>(data0), std::forward<TTs>(datas)...);
		}

		template<typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<E,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
				symphas::lib::types_list<limit_ts...>,
				symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(*static_cast<E0 const*>(&e0), series.limits, series, series.substitution);
		}

		template<int... I00s, int... P00s, int... offsets, typename E0, typename Op, typename... Ts, typename E,
			int... I0s, int... P0s, typename... limit_ts, typename... Is, typename TT0, typename... TTs>
		auto recreate_series_impl(
			std::tuple<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...> const& conditions,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<E,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
				symphas::lib::types_list<limit_ts...>,
				symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series,
			TT0&& data0, TTs&&... datas)
		{
			constexpr size_t num_ind = (
				symphas::lib::types_list_size<symphas::internal::select_all_i_<expr::symbols::i_<I0s, P0s>, op_types_t<E0>>>::value + ...);
			if constexpr (sizeof...(I0s) == sizeof...(I00s)
				&& ((symphas::lib::index_of_value<int, I00s, I0s...> >= 0) && ...)
				&& num_ind == 0)
			{
				return *static_cast<E0 const*>(&e0);
			}
			else
			{
				using condition_ts = symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...>;
				auto fixed_limits = symphas::internal::limits_at_offset(
					condition_ts{}, series.limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});

				return expr::series<Op, expr::symbols::i_<I0s, P0s>...>(
					with_vs(*static_cast<E0 const*>(&e0), filter_grid_symbol_t<op_types_t<E0>>{})
				)(fixed_limits)(std::forward<TT0>(data0), std::forward<TTs>(datas)...);
			}
		}

		template<int... I00s, int... P00s, int... offsets, typename E0, typename Op, typename... Ts, typename E,
			int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			std::tuple<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...> const& conditions,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(conditions, *static_cast<E0 const*>(&e0), series, series.substitution);
		}

		template<int... I00s, int... P00s, int... offsets, typename Op, typename... Ts, typename E, 
			int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			std::tuple<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...> const& conditions,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<E,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
				symphas::lib::types_list<limit_ts...>,
				symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(conditions, series.e, series, series.substitution);
		}
		
		template<int I0, int P0, int offset, typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			expr::symbols::index_eq_N<expr::symbols::i_<I0, P0>, offset> const& condition,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(std::make_tuple(condition), *static_cast<E0 const*>(&e0), series);
		}

		template<int I0, int P0, int offset, typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			expr::symbols::index_eq_N<expr::symbols::i_<I0, P0>, offset> const& condition,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(condition, series.e, series);
		}
		
		template<int I0, int P0, int offset, typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			expr::symbols::index_neq_N<expr::symbols::i_<I0, P0>, offset> const& condition,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(std::make_tuple(condition), *static_cast<E0 const*>(&e0), series);
		}

		template<int I0, int P0, int offset, typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			expr::symbols::index_neq_N<expr::symbols::i_<I0, P0>, offset> const& condition,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(condition, series.e, series);
		}

		template<int I0, int P0, typename E0, typename Op, typename... Ts, typename E,
			int... I0s, int... P0s, typename... limit_ts, typename... Is, typename TT0, typename... TTs>
		auto recreate_series_impl(
			expr::symbols::index_eq<expr::symbols::i_<I0, P0>, expr::symbols::placeholder_N> const& condition,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
				symphas::lib::types_list<limit_ts...>,
				symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series,
			TT0&& data0, TTs&&... datas)
		{
			using condition_ts = symphas::lib::types_list<expr::symbols::index_eq<expr::symbols::i_<I0, P0>, expr::symbols::placeholder_N>>;
			auto fixed_limits = symphas::internal::limits_at_offset(
				condition_ts{}, series.limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});

			auto e00 = expr::transform::swap_grid<OpCoeffSwap<expr::symbols::i_<I0, P0>>>(*static_cast<E0 const*>(&e0), expr::symbols::placeholder_N_symbol{});
			
			constexpr size_t num_ind = (
				symphas::lib::types_list_size<symphas::internal::select_all_i_<expr::symbols::i_<I0s, P0s>, op_types_t<E0>>>::value + ...);
			if constexpr (sizeof...(I0s) == 1
				&& symphas::lib::index_of_value<int, I0, I0s...> == 0
				&& num_ind == 0)
			{
				return e00;
			}
			else
			{
				return expr::series<Op, expr::symbols::i_<I0s, P0s>...>(
					with_vs(e00, filter_grid_symbol_t<op_types_t<E0>>{})
				)(fixed_limits)(std::forward<TT0>(data0), std::forward<TTs>(datas)...);
			}
		}
		
		template<int I0, int P0, typename E0,
			typename Op, typename... Ts, typename E, int... I0s, int... P0s, typename... limit_ts, typename... Is>
		auto recreate_series_impl(
			expr::symbols::index_eq<expr::symbols::i_<I0, P0>, expr::symbols::placeholder_N> const& condition,
			OpExpression<E0> const& e0,
			SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<E,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
					symphas::lib::types_list<limit_ts...>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>> const& series)
		{
			return recreate_series_impl(condition, *static_cast<E0 const*>(&e0), series, series.substitution);
		}

	}

	template<typename... Ts>
	auto recreate_series(Ts&&... args)
	{
		return recreate_series_impl(std::forward<Ts>(args)...);
	}

	template<int N, int P, typename E>
	auto make_list(expr::symbols::i_<N, P>, OpExpression<E> const& e)
	{
		return symphas::internal::make_symbolic_eval(OpIdentity{}, expr::symbols::i_<N, P>{}, *static_cast<E const*>(&e));
	}

	template<typename K, typename E0, typename E>
	auto make_list(OpExpression<E0> const& index, OpExpression<E> const& e)
	{
		return symphas::internal::make_symbolic_eval(
			OpIdentity{}, SymbolicListIndex<E0, K>{ *static_cast<E0 const*>(&index) }, * static_cast<E const*>(&e));
	}
}

namespace symphas::internal
{
	template<typename V, typename Op, typename E, typename Inds, typename... Ts>
	auto make_symbolic_eval_impl(V const& value, SymbolicSeries<Op, E, Inds> const& data, Ts const&... ts)
	{
		return OpSymbolicEval<V, SymbolicSeries<Op, E, Inds>, Ts...>(value, data, ts...);
	}

	template<typename V, typename E0, typename... T0s, typename... Ts>
	auto make_symbolic_eval_impl(V const& value, SymbolicFunction<E0, T0s...> const& data, Ts const&... ts)
	{
		return OpSymbolicEval<V, SymbolicFunction<E0, T0s...>, Ts...>(value, data, ts...);
	}

	//template<typename V, typename E, typename T0, typename... Ts>
	//auto make_symbolic_eval_impl(V const& value, OpExpression<E> const& e, T0 const& t0, Ts const&... ts)
	//{
	//	return value * (*static_cast<E const*>(&e));
	//}

	template<typename V, expr::NoiseType nt, typename T, size_t D, typename E, size_t... Ns, typename... Ts>
	auto make_symbolic_eval_impl(V const& value, NoiseData<nt, T, D> const& noise, SymbolicFunction<E, Variable<Ns, Ts>...> const& f)
	{
		return OpSymbolicEval(value, noise, f);
	}

	template<typename V, expr::NoiseType nt, typename T, size_t D, typename E>
	auto make_symbolic_eval_impl(V const& value, NoiseData<nt, T, D> const& noise, OpExpression<E> const& e)
	{
		return make_symbolic_eval_impl(value, noise, expr::function_of() = *static_cast<E const*>(&e));
	}

	template<typename V, int N, int P, typename E, size_t... Ns, typename... Ts>
	auto make_symbolic_eval_impl(V const& value, expr::symbols::i_<N, P>, SymbolicFunction<E, Variable<Ns, Ts>...> const& f)
	{
		return OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Variable<Ns, Ts>...>>(value, expr::symbols::i_<N, P>{}, f);
	}

	template<typename V, int N, int P, typename E>
	auto make_symbolic_eval_impl(V const& value, expr::symbols::i_<N, P>, OpExpression<E> const& e)
	{
		return make_symbolic_eval_impl(value, expr::symbols::i_<N, P>{}, expr::function_of() = *static_cast<E const*>(&e));
	}

	template<typename V, typename E0, typename K, typename E, size_t... Ns, typename... Ts, std::enable_if_t<!std::is_same<K, void>::value, int> = 0>
	auto make_symbolic_eval_impl(V const& value, SymbolicListIndex<E0, K> const& data, SymbolicFunction<E, Variable<Ns, Ts>...> const& f)
	{
		return OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Variable<Ns, Ts>...>>(value, data, f);
	}

	template<typename V, typename E0, typename K, typename E, std::enable_if_t<!std::is_same<K, void>::value, int> = 0>
	auto make_symbolic_eval_impl(V const& value, SymbolicListIndex<E0, K> const& data, OpExpression<E> const& e)
	{
		return make_symbolic_eval_impl(value, data, expr::function_of() = *static_cast<E const*>(&e));
	}

	template<typename V, int N, int P, typename E, size_t... Ns, typename... Ts>
	auto make_symbolic_eval_impl(V const& value, SymbolicListIndex<expr::symbols::i_<N, P>> const& data, SymbolicFunction<E, Variable<Ns, Ts>...> const& f)
	{
		return make_symbolic_eval_impl(value, data.e, f);
	}

	template<typename V, int N, int P, typename E>
	auto make_symbolic_eval_impl(V const& value, SymbolicListIndex<expr::symbols::i_<N, P>> const& data, OpExpression<E> const& e)
	{
		return make_symbolic_eval_impl(value, data.e, *static_cast<E const*>(&e));
	}

	//template<typename V, typename E0, typename E, typename... Ts>
	//auto make_symbolic_eval_impl(V const& value, OpExpression<E0> const& data, SymbolicFunction<E, Ts...> const& f)
	//{
	//	return make_symbolic_eval(value, SymbolicListIndex{ *static_cast<E0 const*>(&data) }, f);
	//}

	//template<typename V, typename E0, typename E>
	//auto make_symbolic_eval_impl(V const& value, OpExpression<E0> const& data, OpExpression<E> const& e)
	//{
	//	return make_symbolic_eval(value, SymbolicListIndex{ *static_cast<E0 const*>(&data) }, *static_cast<E const*>(&e));
	//}


	template<typename V, typename sub_t, typename... Ts>
	auto make_symbolic_eval(V const& value, sub_t const& data, Ts&&... ts)
	{
		return make_symbolic_eval_impl(value, data, std::forward<Ts>(ts)...);
	}
}

template<typename E, size_t... ArgNs, typename... Ts>
template<typename... T0s>
auto SymbolicFunction<E, Variable<ArgNs, Ts>...>::substitute(Substitution<T0s...> const& s) const
{
	auto tmpl = (expr::template_of() = *this);
	return s.substitute(tmpl);
}

template<typename E, size_t... ArgNs>
template<typename... T0s>
auto SymbolicTemplate<E, ArgNs...>::substitute(Substitution<T0s...> const& s) const
{
	return s.substitute(*this);
}




