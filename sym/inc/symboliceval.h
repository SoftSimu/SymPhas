
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


template<typename V, typename sub_t, typename E, typename... Ts>
struct OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> : OpExpression<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>>
{
	using eval_t = SymbolicFunction<E, Ts...>;
	using this_t = OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>;

	OpSymbolicEval(V value, sub_t const& data, eval_t const& f)
		: value{ value }, f{ f }, data{ data } {}

	V value;
	eval_t f;
	sub_t data;

public:

	auto eval(iter_type n) const
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
		return f.e.print(out);
	}

	size_t print(char* out) const
	{
		return f.e.print(out);
	}

	size_t print_length() const
	{
		return f.e.print_length();
	}

#endif

};

template<typename V, typename sub_t, typename E, typename... Ts>
struct OpSymbolicEval<V, sub_t, SymbolicFunction<E, symphas::lib::types_list<Ts...>>> : OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>
{
	using parent_type = OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>;
	using parent_type::parent_type;
};


template<typename coeff_t, typename sub_t, typename V, typename... Ts,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpSymbolicEval<V, sub_t, Ts...> const& e)
{
	return symphas::internal::make_symbolic_eval(value * e.value, e.data, e.f);
}

template<typename coeff_t, typename sub_t, typename tensor_t, typename... Ts,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpSymbolicEval<tensor_t, sub_t, Ts...>& e)
{
	return (value * e.value) * symphas::internal::make_symbolic_eval(OpIdentity{}, e.data, e.f);
}




namespace symphas::internal
{


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
			return std::make_pair(-expr::coeff(e), OpVoid{});
		}
		else
		{
			return std::make_pair(OpVoid{}, e);
		}
	}

	template<int I0, typename T, typename = std::enable_if_t<expr::is_coeff<T>, int>>
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
	struct minor_index<E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>
	{

	protected:

		using index_list = expr::op_types_t<E>;
		using mask_t = std::integer_sequence<bool, (index_of_type<expr::symbols::i_<I0s, P0s>, index_list> >= 0)...>;
		using in_list_t = symphas::lib::filter_types<types_list<expr::symbols::i_<I0s, P0s>...>, mask_t>;

	public:

		using type = typename minor_index<in_list_t, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::type;
		static const int value = minor_index<in_list_t, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
	};

	template<int I0, int P0, int... I0s, int... P0s>
	struct minor_index<expr::symbols::i_<I0, P0>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>
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


	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T,
		typename... T1s, typename... T2s, typename = std::enable_if_t<expr::is_coeff<T>, int>>
		auto propagate_lower_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			expr::symbols::i_<I1, P1>, T const& lower,
			expr::series_limits<T1, int> const& limit0,
			std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& lower,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int P0, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename T2, typename T,
		typename... T1s, typename... T2s, typename = std::enable_if_t<expr::is_coeff<T>, int>>
		auto propagate_upper_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			expr::symbols::i_<I1, P1>, T const& upper,
			expr::series_limits<int, T2> const& limit0,
			std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits);

	template<int P0, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, std::tuple<expr::series_limits<T1s, T2s>...> const& limits);




	template<int... I0s, int... P0s, int I1, int P1, typename T1, typename T,
		typename... T1s, typename... T2s, typename>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, T const& lower,
		expr::series_limits<T1, int> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		return update_limits<size_t(N)>(
			expr::series_limits(limit0._0, std::min(int(lower.eval()), limit0._1)),
			limits);
	}

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& lower,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = minor_index<E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
		constexpr int L = symphas::lib::index_of_value<int, I1, I0s...>;

		if constexpr (N == L)
		{
			// we have that lhs <= rhs, so the upper limit is rhs
			auto [lhs, rhs] = separate_index<I1>(lower - expr::symbols::i_<I1, P1>{});
			return propagate_lower_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				expr::symbols::i_<I1, P1>{}, rhs / lhs,
				limit0, limits);
		}
		else if constexpr (N >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(N), expr::symbols::i_<I0s, P0s>...>;

			auto limit1 = std::get<size_t(L)>(limits);
			auto updated_limits = update_limits<size_t(L)>(
				expr::series_limits(limit1._0, lower), limits);

			return propagate_lower_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, lower - limit1._0 + index_t{},
				std::get<N>(updated_limits), updated_limits);

		}
		else
		{

			constexpr int L = minor_index<decltype(limit0._1), symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
			if constexpr (L >= 0)
			{
				using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;
				auto [lhs, rhs] = separate_index(index_t{}, lower - limit0._1 + index_t{});

				return propagate_lower_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					index_t{}, rhs / lhs,
					std::get<L>(limits), limits);
			}
			else
			{
				return propagate_lower_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					expr::symbols::i_<I1, P1>{}, expr::make_literal(lower),
					limit0, limits);
			}
		}
	}


	// entry point
	template<int P0, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_lower_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, int start, int end,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		auto limit0 = std::get<size_t(N)>(limits);

		// update the lower index for this starting index
		auto updated_lower = propagate_upper_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			expr::symbols::i_<I1, P1>{}, expr::make_literal(start + P0),
			expr::series_limits(std::min(start + P0, end + 1), limit0._1), limits);

		constexpr int L = minor_index<decltype(limit0._0), symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
		if constexpr (L >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;
			auto [lhs, rhs] = separate_index(index_t{}, limit0._0 - expr::make_literal(start + P0));

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



	template<int... I0s, int... P0s, int I1, int P1, typename T2, typename T,
		typename... T1s, typename... T2s, typename>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, T const& upper,
		expr::series_limits<int, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		return update_limits<size_t(N)>(
			expr::series_limits(std::max(int(upper.eval()), limit0._0), limit0._1),
			limits);
	}

	template<int... I0s, int... P0s, int I1, int P1, typename E,
		typename T1, typename T2, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, E const& upper,
		expr::series_limits<T1, T2> const& limit0,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = minor_index<E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
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

			auto limit1 = std::get<size_t(L)>(limits);
			auto updated_limits = update_limits<size_t(L)>(
				expr::series_limits(upper, std::get<size_t(L)>(limits)._1), limits);

			return propagate_upper_limit(
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
				index_t{}, upper - limit1._1 + index_t{},
				std::get<N>(updated_limits), updated_limits);

		}
		else
		{
			constexpr int L = minor_index<decltype(limit0._0), symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
			if constexpr (L >= 0)
			{
				using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;
				auto [lhs, rhs] = separate_index(index_t{}, upper - limit0._0 + index_t{});

				return propagate_upper_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					index_t{}, rhs / lhs,
					std::get<L>(limits), limits);
			}
			else
			{
				return propagate_upper_limit(
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
					expr::symbols::i_<I1, P1>{}, expr::make_literal(upper),
					limit0, limits);
			}
		}
	}


	// entry point
	template<int P0, int... I0s, int... P0s, int I1, int P1, typename... T1s, typename... T2s>
	auto propagate_upper_limit(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		expr::symbols::i_<I1, P1>, int start, int end,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		constexpr int N = symphas::lib::index_of_value<int, I1, I0s...>;
		auto limit0 = std::get<size_t(N)>(limits);

		// update the upper index for this starting index
		auto updated_upper = propagate_lower_limit(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			expr::symbols::i_<I1, P1>{}, expr::make_literal(start + P0),
			expr::series_limits(limit0._0, std::min(start + P0, end)), limits);

		constexpr int L = minor_index<decltype(limit0._1), symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>::value;
		if constexpr (L >= 0)
		{
			using index_t = symphas::lib::type_at_index<size_t(L), expr::symbols::i_<I0s, P0s>...>;
			auto [lhs, rhs] = separate_index(index_t{}, limit0._1 - expr::make_literal(P0 + start));

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


	template<size_t N, int P0, int... I0s, int... P0s, typename... T1s, typename... T2s>
	auto limits_at_offset(
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		using current_ind_t = symphas::lib::type_at_index<N, expr::symbols::i_<I0s, P0s>...>;
		using seq_i = std::integer_sequence<int, I0s...>;

		auto limit0 = std::get<size_t(N)>(limits);
		int start = symphas::internal::limit_start(seq_i{}, limits, limit0);
		int end = symphas::internal::limit_end(seq_i{}, limits, limit0);

		auto lower_limits = propagate_lower_limit<P0>(
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{},
			current_ind_t{}, start, end, limits);
		auto upper_limits = propagate_upper_limit<P0>(
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
			return limits_at_offset<size_t(N), offset>(limits,
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

}

namespace expr
{
	namespace
	{
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
			return expr::series<Op, expr::symbols::i_<I0s, P0s>...>(*static_cast<E0 const*>(&e0))(limits)
				(std::forward<TT0>(data0), std::forward<TTs>(datas)...);
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
			return recreate_series_impl(*static_cast<E0 const*>(&e0), series.limits, series);
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
			using condition_ts = symphas::lib::types_list<expr::symbols::index_eq_N<expr::symbols::i_<I00s, P00s>, offsets>...>;
			auto fixed_limits = symphas::internal::limits_at_offset(
				condition_ts{}, series.limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
			
			return expr::series<Op, expr::symbols::i_<I0s, P0s>...>(*static_cast<E0 const*>(&e0))(fixed_limits)
				(std::forward<TT0>(data0), std::forward<TTs>(datas)...);
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

	}

	template<typename... Ts>
	auto recreate_series(Ts&&... args)
	{
		return recreate_series_impl(std::forward<Ts>(args)...);
	}

}

template<typename T, size_t... Zs>
auto test_array_sub(SymbolicVariableArray<T, Zs...> const& sub)
{
	return sub;
}

template<int I0>
inline auto test_array_sub(SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<I0, 0>>> const& sub)
{
	return sub;
}

template<typename... T1s, typename... T2s, int... I0s, int... P0s>
inline auto test_array_limits(std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
	symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
{
	using ii = symphas::lib::type_at_index<1, expr::symbols::i_<I0s, P0s>...>;
	return symphas::internal::limits_at_offset(
		ii{} = expr::val<2>,
		limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
}

template<typename V, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename C>
auto test_sum(OpSum<V, E, 
	Substitution<SymbolicDataArray<Ts>...>, 
	symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, C> const& sum)
{
	auto subs = sum.data.substitution;
	auto e = sum.data.e;
	return test_array_limits(sum.data.limits, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
	//return test_array_sub(std::get<1>(subs));
}




template<typename V, typename sub_t, typename... Ts>
auto symphas::internal::make_symbolic_eval(V const& value, sub_t const& data, Ts const&... ts)
{
	return OpSymbolicEval<V, sub_t, Ts...>(value, data, ts...);
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




