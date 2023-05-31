
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
 * PURPOSE: Defines the foundational symbolic algebra elements, such as the
 * base expression object and the binary operations.
 *
 * ***************************************************************************
 */

#pragma once

#include <iostream>


#include "expressiontransforms.h"


ADD_EXPR_TYPE_SYMBOL(x)
ADD_EXPR_TYPE_SYMBOL(y)
ADD_EXPR_TYPE_SYMBOL(z)
ADD_EXPR_TYPE_SYMBOL(h)


namespace expr::symbols
{
	namespace internal
	{
		template<int N>
		struct S1_symbol : Symbol {};

		template<int N0, int N1>
		struct S2_symbol : Symbol {};

		template<int N0, int N1, int N2>
		struct S3_symbol : Symbol {};

		template<int N>
		using S1 = OpTerm<OpIdentity, S1_symbol<N>>;

		template<int N0, int N1>
		using S2 = OpTerm<OpIdentity, S2_symbol<N0, N1>>;

		template<int N0, int N1, int N2>
		using S3 = OpTerm<OpIdentity, S3_symbol<N0, N1, N2>>;
	}
}

ALLOW_COMBINATION((int N), (expr::symbols::internal::S1_symbol<N>))
ALLOW_COMBINATION((int N0, int N1), (expr::symbols::internal::S2_symbol<N0, N1>))
ALLOW_COMBINATION((int N0, int N1, int N2), (expr::symbols::internal::S3_symbol<N0, N1, N2>))


DEFINE_SYMBOL_ID((int N), (expr::symbols::internal::S1_symbol<N>), static char* name = expr::print_with_subscript<N>("S").new_str(); return name;)
DEFINE_SYMBOL_ID((int N0, int N1), (expr::symbols::internal::S2_symbol<N0, N1>),
	static auto name1 = expr::print_with_subscript<N0>("S");
	static char* name2 = expr::print_with_subscript<N1>(name1.value).new_str(); return name2;)
DEFINE_SYMBOL_ID((int N0, int N1, int N2), (expr::symbols::internal::S3_symbol<N0, N1, N2>),
	static auto name1 = expr::print_with_subscript<N0>("S");
	static auto name2 = expr::print_with_subscript<N1>(name1.value);
	static char* name3 = expr::print_with_subscript<N2>(name2.value).new_str(); return name3;)

namespace symphas::internal
{

	using namespace symphas::lib;

	template<int... Ns>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, Ns...>)
	{
		return types_list<expr::symbols::internal::S1_symbol<Ns>...>{};
	}

	template<int... N0s, int N1>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, N0s...>, std::integer_sequence<int, N1>)
	{
		return types_list<expr::symbols::internal::S2_symbol<N0s, N1>...>{};
	}

	template<int... N0s, int N11, int N12, int... N1s>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, N0s...>, std::integer_sequence<int, N11, N12, N1s...>)
	{
		return expand_types_list<
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N11>{})),
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N12, N1s...>{}))>{};
	}


	template<int... N0s, int N1, int N2>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, N0s...>, std::integer_sequence<int, N1>, std::integer_sequence<int, N2>)
	{
		return types_list<expr::symbols::internal::S3_symbol<N0s, N1, N2>...>{};
	}

	template<int... N0s, int N11, int N12, int... N1s, int N21>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, N0s...>,
		std::integer_sequence<int, N11, N12, N1s...>, std::integer_sequence<int, N21>)
	{
		return expand_types_list<
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N11>{}, std::integer_sequence<int, N21>{})),
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N12, N1s...>{}, std::integer_sequence<int, N21>{}))>{};
	}


	template<int... N0s, int N11, int N12, int... N1s, int N21, int N22, int... N2s>
	constexpr auto get_stencil_symbols(std::integer_sequence<int, N0s...>,
		std::integer_sequence<int, N11, N12, N1s...>, std::integer_sequence<int, N21, N22, N2s...>)
	{
		return expand_types_list<
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N11, N12, N1s...>{}, std::integer_sequence<int, N21>{})),
			decltype(get_stencil_symbols(std::integer_sequence<int, N0s...>{}, std::integer_sequence<int, N11, N12, N1s...>{}, std::integer_sequence<int, N22, N2s...>{}))>{};
	}


	template<size_t D, int R>
	constexpr auto get_stencil_symbols()
	{
		constexpr auto seq = seq_add(std::make_index_sequence<2 * R + 1>{}, seq_repeating_value_t<2 * R + 1, int, -R>{});
		if constexpr (D == 1)
		{
			return get_stencil_symbols(seq);
		}
		else if constexpr (D == 2)
		{
			return get_stencil_symbols(seq, seq);
		}
		else
		{
			return get_stencil_symbols(seq, seq, seq);
		}
	}


	template<int I0, int I1>
	constexpr auto compare_symbols(expr::symbols::internal::S1_symbol<I0>, expr::symbols::internal::S1_symbol<I1>)
	{
		return I0 < I1;
	}

	template<int I0, int J0, int I1, int J1>
	constexpr auto compare_symbols(expr::symbols::internal::S2_symbol<I0, J0>, expr::symbols::internal::S2_symbol<I1, J1>)
	{
		return (J0 == J1) ? (I0 < I1) : (J0 < J1);
	}

	template<int I0, int J0, int K0, int I1, int J1, int K1>
	constexpr auto compare_symbols(expr::symbols::internal::S3_symbol<I0, J0, K0>, expr::symbols::internal::S3_symbol<I1, J1, K1>)
	{
		return (K0 == K1) ? ((J0 == J1) ? (I0 < I1) : (J0 < J1)) : (K0 < K1);
	}



	//template<int... Is, int... Js, typename... Es>
	//constexpr auto get_stencil_coeffs(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict)
	//{
	//	return typename to_tuples<decltype(get_rows(dict))>::type{};
	//}

	//template<int... Is, typename... Es>
	//constexpr auto get_stencil_coeffs(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...> const& dict)
	//{
	//	return std::tuple<Es...>{};
	//}

	//template<int... Is>
	//constexpr auto get_stencil_coeffs(types_list<expr::symbols::internal::S1_symbol<Is>...>)
	//{
	//	return std::tuple<OpTerm<OpIdentity, expr::symbols::internal::S1_symbol<Is>>...>{};
	//}

	//template<size_t D, int R>
	//constexpr auto get_stencil_coeffs()
	//{
	//	return get_stencil_coeffs(get_stencil_symbols<D, R>());
	//}


	template<typename T>
	auto neg_one_element(T const&)
	{
		return std::invoke_result_t<decltype(&T::operator-), T>{};
	}

	template<typename... Ts>
	auto neg_elements(types_list<Ts...> const&)
	{
		return types_list<std::invoke_result_t<decltype(&Ts::operator-), Ts>...>{};
	}

	template<typename List, typename... Seqs>
	struct pack_dictionary;

	template<>
	struct pack_dictionary<types_list<>, std::integer_sequence<int>>
	{
		using type = types_list<>;
	};

	template<>
	struct pack_dictionary<types_list<>, std::integer_sequence<int>, std::integer_sequence<int>>
	{
		using type = types_list<>;
	};

	template<>
	struct pack_dictionary<types_list<>, std::integer_sequence<int>, std::integer_sequence<int>, std::integer_sequence<int>>
	{
		using type = types_list<>;
	};

	template<typename E0, typename... Es, int I0, int... Is>
	struct pack_dictionary<types_list<E0, Es...>, std::integer_sequence<int, I0, Is...>>
	{
		using type = expand_types_list<
			std::pair<expr::symbols::internal::S1_symbol<I0>, E0>,
			typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>>::type>;
	};

	template<typename... Es, int I0, int... Is>
	struct pack_dictionary<types_list<OpVoid, Es...>, std::integer_sequence<int, I0, Is...>>
	{
		using type = typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>>::type;
	};


	template<typename E0, typename... Es, int I0, int... Is, int J0, int... Js>
	struct pack_dictionary<types_list<E0, Es...>, std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>>
	{
		using type = expand_types_list<
			std::pair<expr::symbols::internal::S2_symbol<I0, J0>, E0>,
			typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type>;
	};

	template<typename... Es, int I0, int... Is, int J0, int... Js>
	struct pack_dictionary<types_list<OpVoid, Es...>, std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>>
	{
		using type = typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type;
	};


	template<typename E0, typename... Es, int I0, int... Is, int J0, int... Js, int K0, int... Ks>
	struct pack_dictionary<types_list<E0, Es...>, std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>, std::integer_sequence<int, K0, Ks...>>
	{
		using type = expand_types_list<
			std::pair<expr::symbols::internal::S3_symbol<I0, J0, K0>, E0>,
			typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type>;
	};

	template<typename... Es, int I0, int... Is, int J0, int... Js, int K0, int... Ks>
	struct pack_dictionary<types_list<OpVoid, Es...>, std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>, std::integer_sequence<int, K0, Ks...>>
	{
		using type = typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
	};


	//template<int I0, int J0, typename E0>
	//constexpr auto get_rows(types_list<std::pair<expr::symbols::internal::S2_symbol<I0, J0>, E0>>)
	//{
	//	return types_list<types_list<E0>>{};
	//}

	//template<int I0, int I1, int... Is, int J0, int J1, int... Js, typename E0, typename E1, typename... Es>
	//constexpr auto get_rows(types_list<
	//	std::pair<expr::symbols::internal::S2_symbol<I0, J0>, E0>,
	//	std::pair<expr::symbols::internal::S2_symbol<I1, J1>, E1>,
	//	std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>)
	//{
	//	if constexpr (I1 <= I0)
	//	{
	//		using next_list_t = decltype(get_rows(types_list<
	//			std::pair<expr::symbols::internal::S2_symbol<I1, J1>, E1>,
	//			std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>{}));

	//		return expand_types_list<types_list<types_list<E0>>, next_list_t>{};
	//	}
	//	else
	//	{
	//		using next_list_t = decltype(get_rows(types_list<
	//			std::pair<expr::symbols::internal::S2_symbol<I1, J1>, E1>,
	//			std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>{}));

	//		using first_t = type_at_index<0, unroll_types_list<next_list_t>>;
	//		using rest_t = types_after_at_index<1, unroll_types_list<next_list_t>>;

	//		return expand_types_list<
	//			types_list<expand_types_list<E0, first_t>>,
	//			rest_t>{};
	//	}
	//}


	//template<size_t>
	//using zero_t = OpVoid;
	//template<typename>
	//using zero_tt = OpVoid;

	//template<size_t O, typename... E0s, typename... Es, size_t... Is, size_t N2 = sizeof...(Es) / 2>
	//auto set_zeros_2d_row(types_list<E0s...>, types_list<Es...>, std::index_sequence<Is...>)
	//{
	//	if constexpr (O % 2 == 0)
	//	{
	//		using first_half_t = types_between_index<sizeof...(Is) + 1, N2 + 1, E0s...>;
	//		using center_t = types_list<type_at_index<N2, Es...>>;

	//		using new_row_t = expand_types_list<
	//			types_list<zero_t<Is>...>,
	//			first_half_t,
	//			center_t,
	//			reverse_types_list<first_half_t>,
	//			types_list<zero_t<Is>...>>;

	//		return new_row_t{};
	//	}
	//	else
	//	{
	//		using first_half_t = types_between_index<sizeof...(Is), N2, Es...>;

	//		using new_row_t = expand_types_list <
	//			types_list<zero_t<Is>...>,
	//			first_half_t,
	//			types_list<OpVoid>,
	//			reverse_types_list<decltype(neg_elements(first_half_t{})) > ,
	//			types_list<zero_t<Is>... >> ;

	//		return new_row_t{};
	//	}
	//}


	//template<size_t O, typename... E0s, typename... Es, size_t... Is, size_t N2 = sizeof...(Es) / 2>
	//auto set_zeros_2d_row(types_list<Es...>, std::index_sequence<Is...>)
	//{
	//	if constexpr (O % 2 == 0)
	//	{
	//		using first_half_t = types_between_index<sizeof...(Is), N2, Es...>;
	//		using center_t = types_list<type_at_index<N2, Es...>>;

	//		using new_row_t = expand_types_list<
	//			types_list<zero_t<Is>...>,
	//			first_half_t,
	//			center_t,
	//			reverse_types_list<first_half_t>,
	//			types_list<zero_t<Is>...>>;

	//		return new_row_t{};
	//	}
	//	else
	//	{
	//		using first_half_t = types_between_index<sizeof...(Is), N2, Es...>;

	//		using new_row_t = expand_types_list<
	//			types_list<zero_t<Is>...>,
	//			first_half_t,
	//			types_list<OpVoid>,
	//			reverse_types_list<first_half_t>,
	//			types_list<zero_t<Is>...>>;

	//		return new_row_t{};
	//	}
	//}


	//template<typename T>
	//struct combine_rows
	//{
	//	using type = types_list<>;
	//};

	//template<typename... Es, typename... Ts>
	//struct combine_rows<types_list<types_list<Es...>, Ts...>>
	//{
	//	using type = expand_types_list<
	//		types_list<Es...>,
	//		typename combine_rows<types_list<Ts...>>::type>;
	//};

	//template<size_t O, size_t N, typename... Es>
	//auto set_zeros_2d(types_list<Es...> last_row, types_list<> const& rows, std::index_sequence<>)
	//{
	//	return types_list<>{};
	//}

	//template<size_t O, size_t N, typename... Es, typename T0, typename... Ts>
	//auto set_zeros_2d(types_list<Es...> last_row, types_list<T0, Ts...> const& rows, std::index_sequence<>)
	//{
	//	using row_t = T0;
	//	using rest_t = types_after_at_index<1, T0, Ts...>;

	//	auto new_row = set_zeros_2d_row<O>(last_row, row_t{}, std::index_sequence<>{});

	//	return expand_types_list<
	//		types_list<decltype(new_row)>,
	//		decltype(set_zeros_2d<O, N>(new_row, rest_t{}, std::index_sequence<>{}))>{};
	//}

	//template<size_t O, size_t N, typename... Es, typename... Ts, size_t I0, size_t... Is>
	//auto set_zeros_2d(types_list<Es...> last_row, types_list<Ts...> const& rows, std::index_sequence<I0, Is...>)
	//{
	//	using row_t = type_at_index<0, Ts...>;
	//	using rest_t = types_after_at_index<1, Ts...>;

	//	constexpr size_t N1 = (N < I0) ? 0 : N - I0;
	//	auto new_row = set_zeros_2d_row<O>(last_row, row_t{}, std::make_index_sequence<N1>{});

	//	return expand_types_list<
	//		types_list<decltype(new_row)>,
	//		decltype(set_zeros_2d<O, N>(new_row, rest_t{}, std::index_sequence<Is...>{}))>{};
	//}

	//template<size_t O, size_t N, typename... Ts, size_t I0, size_t... Is>
	//auto set_zeros_2d(types_list<Ts...> const& rows, std::index_sequence<I0, Is...>)
	//{
	//	using row_t = type_at_index<0, Ts...>;
	//	using rest_t = types_after_at_index<1, Ts...>;

	//	constexpr size_t N1 = (N < I0) ? 0 : N - I0;
	//	auto new_row = set_zeros_2d_row<O>(row_t{}, std::make_index_sequence<N1>{});

	//	return expand_types_list<
	//		types_list<decltype(new_row)>,
	//		decltype(set_zeros_2d<O, N>(new_row, rest_t{}, std::index_sequence<Is...>{}))>{};
	//}

	//

	template<int R, int... As>
	constexpr bool test_in_radius = (fixed_abs<As> + ...) <= R;

	template<int I, int J, typename Ts, typename SeqI, typename SeqJ>
	struct at_stencil_index_2d;

	template<int I, int J, typename T0, typename... Ts, int I0, int... Is, int J0, int... Js>
	struct at_stencil_index_2d<I, J, types_list<T0, Ts...>,
		std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>>
	{
		using type = typename at_stencil_index_2d<I, J, types_list<Ts...>,
			std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type;
	};

	template<int I, int J, typename T0, typename... Ts, int... Is, int... Js>
	struct at_stencil_index_2d<I, J, types_list<T0, Ts...>,
		std::integer_sequence<int, I, Is...>, std::integer_sequence<int, J, Js...>>
	{
		using type = T0;
	};
	

// This number is subtracted from the whole radius to get
// the radius of values populated in a 2d central stencil 
// (not including axis)
#define NEG_RADIUS_CENTRAL_2D 0
#define USE_RADIUS_CENTRAL_2D true

// this number is subtracted from the whole radius to get
// the radius of values populated in a 2d central stencil 
// (not including axis)
// A value of -1 indicates there is no radius and only
// axes are used.
#define NEG_RADIUS_CENTRAL_3D 0
#define USE_RADIUS_CENTRAL_3D true

	template<int... Is, int... Js>
	auto simplify_axes_2d_odd(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>>...>)
	{
		using Ts = types_list<OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>...>;
		using nTs = types_list<OpTerm<OpNegIdentity, expr::symbols::internal::S2_symbol<Is, Js>>...>;

		return types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>,
			typename at_stencil_index_2d<
				-fixed_abs<Is>,
				-fixed_abs<Js>,
				std::conditional_t<(Is < 0), Ts, nTs>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type>...
			>{};
	}

	template<int... Is, int... Js>
	auto simplify_axes_2d_even(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>>...>)
	{
		using Ts = types_list<OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>...>;

		return types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>,
			typename at_stencil_index_2d<
				-fixed_max<fixed_abs<Js>, fixed_abs<Is>>,
				-fixed_min<fixed_abs<Js>, fixed_abs<Is>>,
				Ts, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type>...
			>{};
	}


	template<size_t O, int... Is, int... Js>
	auto keep_axes_2d(types_list<OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>...>, 
		std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>)
	{
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int RJ = JM - Jm + 1;

		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int RI = IM - Im + 1;

		constexpr int R = fixed_min<RI, RJ>;

		if constexpr (O == 0)
		{
			return typename pack_dictionary<
				types_list<std::conditional_t<
					Is == 0 || Js == 0 || (
						(USE_RADIUS_CENTRAL_2D) 
						? (test_in_radius<R / 2 - NEG_RADIUS_CENTRAL_2D + 1, Is, Js>)
						: false),
					OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>
			>::type{};
		}
		else if constexpr (O % 2 == 0)
		{
			return simplify_axes_2d_even(
				typename pack_dictionary<
					types_list<std::conditional_t<
						Is == 0 || Js == 0 || (
							(USE_RADIUS_CENTRAL_2D) 
							? (test_in_radius<R / 2 - NEG_RADIUS_CENTRAL_2D, Is, Js>)
							: false),
						OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>,
						OpVoid>...>,
					std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>
				>::type{});
		}
		else
		{
			return simplify_axes_2d_odd(
				typename pack_dictionary<
					types_list<std::conditional_t<
						((Js > Jm && Js < JM) && Is != 0 && test_in_radius<RJ / 2, Is, Js>),
						OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>,
						OpVoid>...>,
					std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>
				>::type{});
		}
	}


	template<int I, int J, int K, typename Ts, typename SeqI, typename SeqJ, typename SeqK>
	struct at_stencil_index_3d;

	template<int I, int J, int K, typename T0, typename... Ts, int I0, int... Is, int J0, int... Js, int K0, int... Ks>
	struct at_stencil_index_3d<I, J, K, types_list<T0, Ts...>,
		std::integer_sequence<int, I0, Is...>, std::integer_sequence<int, J0, Js...>, std::integer_sequence<int, K0, Ks...>>
	{
		using type = typename at_stencil_index_3d<I, J, K, types_list<Ts...>,
			std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
	};

	template<int I, int J, int K, typename T0, typename... Ts, int... Is, int... Js, int... Ks>
	struct at_stencil_index_3d<I, J, K, types_list<T0, Ts...>,
		std::integer_sequence<int, I, Is...>, std::integer_sequence<int, J, Js...>, std::integer_sequence<int, K, Ks...>>
	{
		using type = T0;
	};

	template<int... Is, int... Js, int... Ks>
	auto simplify_axes_3d_odd(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, 
		OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>>...>)
	{
		using Ts = types_list<OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>...>;
		using nTs = types_list<OpTerm<OpNegIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>...>;

		return types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>,
			typename at_stencil_index_3d<
				-fixed_abs<Is>,
				-fixed_max<fixed_abs<Js>, fixed_abs<Ks>>,
				-fixed_min<fixed_abs<Js>, fixed_abs<Ks>>,
				std::conditional_t<(Is < 0), Ts, nTs>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type>...
			>{};
	}

	template<int... Is, int... Js, int... Ks>
	auto simplify_axes_3d_even(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>>...>)
	{
		using Ts = types_list<OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>...>;

		return types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>,
			typename at_stencil_index_3d<
				-fixed_max<fixed_abs<Is>, fixed_abs<Js>, fixed_abs<Ks>>,
				-fixed_min<fixed_max<fixed_abs<Is>, fixed_abs<Js>>, fixed_max<fixed_abs<Is>, fixed_abs<Ks>>, fixed_max<fixed_abs<Js>, fixed_abs<Ks>>>,
				-fixed_min<fixed_abs<Is>, fixed_abs<Js>, fixed_abs<Ks>>,
				Ts,	std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type>...
			>{};
	}

	template<size_t O, int... Is, int... Js, int... Ks>
	auto keep_axes_3d(types_list<OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>...>, 
		std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>)
	{
		constexpr int Im = fixed_min<Js...>;
		constexpr int IM = fixed_max<Js...>;
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Km = fixed_min<Ks...>;
		constexpr int KM = fixed_max<Ks...>;
		constexpr int RI = IM - Im + 1;
		constexpr int RJ = JM - Jm + 1;
		constexpr int RK = KM - Km + 1;
		constexpr int R = fixed_min<RI, RJ, RK>;

		if constexpr (O % 2 == 0)
		{

			return simplify_axes_3d_even(
				typename pack_dictionary<
					types_list<std::conditional_t<
						((Is == 0 && Js == 0) || (Is == 0 && Ks == 0) || (Js == 0 && Ks == 0)) 
						|| ((USE_RADIUS_CENTRAL_3D) ? test_in_radius<R / 2 - NEG_RADIUS_CENTRAL_3D, Is, Js, Ks> : false),
						OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>, 
						OpVoid>...>,
					std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>
				>::type{});
		}
		else
		{
			return simplify_axes_3d_odd(
				typename pack_dictionary<
					types_list<std::conditional_t<
						(((Js > Jm && Js < JM) && (Ks > Km && Ks < KM) && Is != 0)
							&& test_in_radius<R / 2, Is, Js, Ks>),
						OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>, 
						OpVoid>...>,
					std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>
				>::type{});
		}
	}


	template<typename Symbol, typename E, typename IE,
		typename std::enable_if_t<(symphas::lib::index_of_type<Symbol, expr::op_types_t<E>> < 0), int> = 0>
	auto make_substitution(E const& e, OpExpression<IE> const&)
	{
		return E{};
	}

	template<typename Symbol, typename IE, typename E, 
		typename std::enable_if_t<(symphas::lib::index_of_type<Symbol, expr::op_types_t<E>> >= 0), int> = 0>
	auto make_substitution(E const& e, OpExpression<IE> const& expr)
	{
		return expr::transform::swap_grid<Symbol>(E{}, IE{});
	}

	template<typename E>
	auto make_all_substitutions(OpExpression<E> const& e, types_list<> const&)
	{
		return E{};
	}

	template<typename E, typename Symbol0, typename E0, typename... Symbols, typename... Es>
	auto make_all_substitutions(OpExpression<E> const& e, types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>)
	{
		return make_all_substitutions(make_substitution<Symbol0>(E{}, E0{}), types_list<std::pair<Symbols, Es>...>{});
	}

	template<typename F, typename T, typename M, typename L>
	struct match_select_impl;

	template<typename T, typename... Ms, typename... Ls>
	struct match_select_impl<std::integral_constant<int, -1>, T, types_list<Ms...>, types_list<Ls...>>
	{
		using type = types_list<>;
	};

	template<int N, typename T, typename... Ms, typename... Ls>
	struct match_select_impl<std::integral_constant<int, N>, T, types_list<Ms...>, types_list<Ls...>>
	{
		using type = type_at_index<size_t(N), Ls...>;
	};

	template<typename T, typename... Ms, typename... Ls>
	struct match_select_impl<void, T, types_list<Ms...>, types_list<Ls...>>
	{
		using type = typename match_select_impl<std::integral_constant<int, index_of_type<T, Ms...>>, T, types_list<Ms...>, types_list<Ls...>>::type;
	};

	template<typename T, typename M, typename L>
	using match_select = typename match_select_impl<void, T, M, L>::type;

	template<typename... Symbols, typename... Es, typename... Ts, typename E>
	auto make_all_substitutions(types_list<std::pair<Symbols, Es>...> const& dict, types_list<Ts...>, OpExpression<E> const& e)
	{
		using list_substitutions_t = expand_types_list<match_select<Ts, types_list<Symbols...>, types_list<std::pair<Symbols, Es>...>>...>;
		return make_all_substitutions(E{}, list_substitutions_t{});
	}

	template<typename... Symbols, typename... Es, typename E>
	auto make_all_substitutions(types_list<std::pair<Symbols, Es>...> const& dict, OpExpression<E> const& e)
	{
		return make_all_substitutions(dict, expr::op_types_t<E>{}, E{});
	}

	//! Perform back substitution on the dictionary to update it.
	template<typename SymbolKey, int... Is, typename... Es, typename E>
	auto back_substitution(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...> const& dict, OpExpression<E> const& expr)
	{
		return types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, decltype(make_substitution<SymbolKey>(Es{}, E{}))>...>{};
	}

	//! Perform back substitution on the dictionary to update it.
	template<typename SymbolKey, int... Is, int... Js, typename... Es, typename E>
	auto back_substitution(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict, OpExpression<E> const& expr)
	{
		return types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, decltype(make_substitution<SymbolKey>(Es{}, E{}))>...>{};
	}

	//! Perform back substitution on the dictionary to update it.
	template<typename SymbolKey, int... Is, int... Js, int... Ks, typename... Es, typename E>
	auto back_substitution(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict, OpExpression<E> const& expr)
	{
		return types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, decltype(make_substitution<SymbolKey>(Es{}, E{}))>...>{};
	}

	//! Perform back substitution on the dictionary to update it.
	template<typename... Ss, typename... S0s, typename... Es>
	auto self_complete(types_list<std::pair<Ss, Es>...> const& dict)
	{
		return types_list<std::pair<Ss, decltype(make_all_substitutions(dict, Es{}))>...>{};
	}


	////! Perform back substitution on the dictionary to update it.
	//template<typename SymbolKey, typename... Symbols, typename... Es>
	//auto back_substitution(types_list<std::pair<Symbols, Es>...> const& dict, OpVoid)
	//{
	//	return types_list<std::pair<Symbols, Es>...>{};
	//}



	namespace
	{
		template<size_t>
		using empty_symbol_expr_entry = OpVoid;

		template<typename... Es>
		struct to_list_position;

		template<typename S, typename E, size_t... Is>
		struct to_list_position<S, E, std::index_sequence<Is...>>
		{
			using elem = decltype(expr::split::factor<sizeof...(Is), S>(std::declval<E>()).second);
			using type = types_list<empty_symbol_expr_entry<Is>..., elem>;
		};


		template<typename T, typename Seq>
		struct extend_symbol_list;

		template<typename... Es, size_t... Is>
		struct extend_symbol_list<types_list<Es...>, std::index_sequence<Is...>>
		{
			using type = types_list<Es..., empty_symbol_expr_entry<Is>...>;
		};

		template<size_t I, typename... Ts>
		struct add_types_at_index
		{
			using type = add_result_t<type_at_index<I, unroll_types_list<Ts>>...>;
		};


		//! Added this type trait in order to remove && from OpAdd when generated in the
		//! merging of types list.
		template<typename T>
		struct filter_amp
		{
			using type = T;
		};

		template<typename T>
		struct filter_amp<T&&>
		{
			using type = T;
		};

		template<typename... Ts>
		struct merge_symbol_exprs;

		template<typename... Ts, size_t... Is>
		struct merge_symbol_exprs<std::index_sequence<Is...>, Ts...>
		{
			using type = reverse_types_list<types_list<typename filter_amp<typename add_types_at_index<Is, Ts...>::type>::type...>>;
		};

		template<typename... E0s, typename... Ts>
		struct merge_symbol_exprs<types_list<E0s...>, Ts...>
		{
			static const size_t N = fixed_max<sizeof...(E0s), types_list_size<Ts>::value...>;

			using type = typename merge_symbol_exprs<
				std::make_index_sequence<N>,
				typename extend_symbol_list<types_list<E0s...>, std::make_index_sequence<N - sizeof...(E0s)>>::type,
				typename extend_symbol_list<Ts, std::make_index_sequence<N - types_list_size<Ts>::value>>::type...>::type;
		};

		template<typename T, typename E>
		struct divide_exprs;

		template<typename... Es, typename E>
		struct divide_exprs<types_list<Es...>, E>
		{
			//using type = types_list<OpBinaryDiv<Es, E>...>;
			using type = types_list<decltype(std::declval<Es>() / std::declval<E>())...>;
		};
	}

	//! Generates a list of types of expressions from the given expression by sorting the powers
	//! of the symbol in each term into list of types indices, where the list index represents the 
	//! multiplicative order of the given variable `S` in that subexpression.
	template<typename S, typename E>
	struct collect_symbol_exprs
	{
		static const size_t I = expr::factor_count<S, E>::value;

		using type = typename to_list_position<S, E, std::make_index_sequence<I>>::type;
	};

	template<typename S, typename... Es>
	struct collect_symbol_exprs<S, OpAdd<Es...>>
	{
		using type = typename merge_symbol_exprs<typename collect_symbol_exprs<S, Es>::type...>::type;
	};

	template<typename S, typename E1, typename E2>
	struct collect_symbol_exprs<S, OpBinaryDiv<E1, E2>>
	{
		using type = typename divide_exprs<typename collect_symbol_exprs<S, E1>::type, E2>::type;
	};



	template<typename S, typename E>
	struct split_by
	{
		using type = typename collect_symbol_exprs<S, E>::type;
	};

	template<typename S, typename... Es>
	struct split_by<S, types_list<Es...>>
	{
		using type = types_list<typename split_by<S, Es>::type...>;
	};



	//! Finds the index of the first undefined symbol in the dictionary.
	/*!
	 * Finds the index of the first undefined symbol in the dictionary. A symbol is undefined
	 * if it is associated with a variable of itself. I.e. given a symbol `S`, it is undefined
	 * if its associated dictionary entry is `S, OpTerm<OpIdentity, S>`.
	 */
	template<typename M, typename T>
	struct first_undefined_symbol
    {
		static const int value = -1;
    };

	template<typename Symbol0, typename... Symbols, typename E0, typename... Es>
	struct first_undefined_symbol<void, types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>>
	{
	protected:

		using dict_t = types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>;

		template<int I>
		static constexpr int get_index_no_matching()
		{
			if constexpr (I >= sizeof...(Symbols) + 1)
			{
				return -1;
			}
			else
			{
				using ith_element = type_at_index<I, unroll_types_list<dict_t>>;
				using ith_symbol = type_at_index<I, Symbol0, Symbols...>;

				if constexpr (std::is_same<ith_element, std::pair<ith_symbol, OpTerm<OpIdentity, ith_symbol>>>::value)
				{
					return I;
				}
				else
				{
					return get_index_no_matching<I + 1>();
				}
			}
		}

		template<int I>
		static constexpr auto get_index()
		{
			return -1;
		}

		template<int I, typename S0, typename... Ss>
		static constexpr auto get_index()
		{
			if constexpr (I >= sizeof...(Symbols) + 1)
			{
				return -1;
			}
			else
			{
				using ith_element = type_at_index<I, unroll_types_list<dict_t>>;
				using ith_symbol = type_at_index<I, Symbol0, Symbols...>;

				using pack_t = types_list<S0, Ss...>;
				constexpr size_t N = index_of_type<ith_symbol, pack_t>;
				if constexpr (N >= 0
					&& std::is_same<ith_element, std::pair<ith_symbol, OpTerm<OpIdentity, ith_symbol>>>::value)
				{
					return I;
				}
				else
				{
					return get_index<I + 1, S0, Ss...>();
				}
			}
		}

	public:

		static const int value = get_index_no_matching<0>();
	};

	template<typename S0, typename... Ss, typename Symbol0, typename... Symbols, typename E0, typename... Es>
	struct first_undefined_symbol<types_list<S0, Ss...>, types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>>
        : first_undefined_symbol<void, types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>>
    {
        using parent_type = first_undefined_symbol<void, types_list<std::pair<Symbol0, E0>, std::pair<Symbols, Es>...>>;
		static const int value = parent_type::template get_index<0, S0, Ss...>();
    };


	//! Determines whether the dictionary is complete or not.
	/*!
	 * Determines whether the dictionary is complete or not. It is determined complete if
	 * there are no undefined symbols.
	 */
	template<typename T>
	struct dict_complete;

	template<typename... Symbols, typename... Es>
	struct dict_complete<types_list<std::pair<Symbols, Es>...>>
	{
	protected:

		using dict_type = types_list<std::pair<Symbols, Es>...>;

		//template<size_t... Is>
		//static constexpr bool get_value(std::index_sequence<Is...>)
		//{
		//	return !((
		//		std::is_same<
		//			type_at_index<Is, Es...>,
		//			OpTerm<OpIdentity, type_at_index<Is, Symbols...>>
		//		>::value || ...));
		//}

	public:

		//static const bool value = get_value(std::make_index_sequence<sizeof...(Symbols)>{});
		static const bool value = !(std::is_same<Es, OpTerm<OpIdentity, Symbols>>::value || ...);
	};


	template<size_t O, size_t N, size_t R0 = O + N - 2>
	constexpr size_t R_ = (R0 <= 1) ? 1 : (R0 % 2 == 1) ? (R0 + 1) / 2 : R0 / 2;

}

	


namespace expr
{

	template<size_t N, size_t I, typename E>
	struct stencil_vector_type {};

	// Represents the vector which is a row vector when I = 0, column when I = 1 and depth when I = 2.
	template<size_t I, typename... Es>
	struct stencil_vector_type<0, I, symphas::lib::types_list<Es...>> {};


	// Represents the 2d matrix, composed of list of row vectors.
	template<size_t I, typename... Es>
	struct stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, I, Es>...>> {};

	template<size_t I0, size_t I, typename... Es>
	auto get_stencil_vector_entry(expr::stencil_vector_type<0, I, symphas::lib::types_list<Es...>>)
	{
		return symphas::lib::type_at_index<I0, Es...>{};
	}

	template<size_t I0, size_t I1, size_t I, typename... Es>
	auto get_stencil_vector_entry(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>...>>)
	{
		return get_stencil_vector_entry<I0>(symphas::lib::type_at_index<I1, stencil_vector_type<0, 0, Es>...>{});
	}

	template<size_t N, size_t I, typename E>
	struct stencil_vector_entry_impl
	{
		using type = decltype(expr::choose<N, N - I>() * expr::pow<N - I>(E{}));
	};

	template<size_t N, size_t I, typename V, typename G0, typename... Gs>
	struct stencil_vector_entry_impl<N, I, OpTerms<V, Term<G0, 1>, Term<Gs, 1>...>>
	{
		using type = OpTerms<typename stencil_vector_entry_impl<N, I, V>::type, Term<G0, expr::Xk<N - I>>, Term<Gs, expr::Xk<N - I>>...>;
	};

	template<size_t N, typename V, typename G0, typename... Gs>
	struct stencil_vector_entry_impl<N, N, OpTerms<V, Term<G0, 1>, Term<Gs, 1>...>>
	{
		using type = typename stencil_vector_entry_impl<N, N, V>::type;
	};

	template<size_t N, size_t I, typename E>
	using stencil_vector_entry = typename stencil_vector_entry_impl<N, I, E>::type;



	template<typename E, size_t... Is>
	auto stencil_vector(expr::symbols::x, E, std::index_sequence<Is...>)
	{
		return stencil_vector_type<0, 0, symphas::lib::types_list<stencil_vector_entry<sizeof...(Is) - 1, Is, E>...>>{};
	}

	template<typename E, size_t... Is>
	auto stencil_vector(expr::symbols::y, E, std::index_sequence<Is...>)
	{
		return stencil_vector_type<0, 1, symphas::lib::types_list<stencil_vector_entry<sizeof...(Is) - 1, Is, E>...>>{};
	}

	template<typename E, size_t... Is>
	auto stencil_vector(expr::symbols::z, E, std::index_sequence<Is...>)
	{
		return stencil_vector_type<0, 2, symphas::lib::types_list<stencil_vector_entry<sizeof...(Is) - 1, Is, E>...>>{};
	}

	template<int I0, typename T, typename Seq>
	struct empty_stencil_vector_type;

	template<int I0, int I, typename... Es, size_t... Is>
	struct empty_stencil_vector_type<I0, expr::stencil_vector_type<0, I, symphas::lib::types_list<Es...>>, std::index_sequence<Is...>>
	{
		using type = expr::stencil_vector_type<0, I, symphas::lib::types_list<symphas::internal::ttype<Es, OpVoid>...>>;
	};

	template<int I0, size_t N, int I, typename... E0s, typename... Es, size_t... Is>
	struct empty_stencil_vector_type<I0, expr::stencil_vector_type<N, I, symphas::lib::types_list<expr::stencil_vector_type<N - 1, I, symphas::lib::types_list<E0s...>>, Es...>>, std::index_sequence<Is...>>
	{
		using empties_t = typename empty_stencil_vector_type<I0, expr::stencil_vector_type<N - 1, I, symphas::lib::types_list<E0s...>>, std::make_index_sequence<sizeof...(E0s)>>::type;

		using type = expr::stencil_vector_type<N, I, symphas::lib::types_list<empties_t, symphas::internal::ttype<Es, empties_t>...>>;
	};

	template<int I0, size_t N, typename T>
	using empty_stencil_vector_t = typename empty_stencil_vector_type<I0, T, std::make_index_sequence<N>>::type;

	template<size_t I, typename... Es, size_t... Is>
	auto fill_stencil_vector(expr::stencil_vector_type<0, I, symphas::lib::types_list<Es...>>, std::index_sequence<Is...>)
	{
		return expr::stencil_vector_type<0, I, symphas::lib::types_list<Es..., symphas::internal::itype<Is, OpVoid>...>>{};
	}

	//template<size_t I, typename... Es, size_t... Is, size_t... Ns>
	//auto fill_stencil_vector(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>...>>, std::index_sequence<Is...>, std::index_sequence<Ns...>)
	//{
	//	using empty_t = stencil_vector_type<0, 0, symphas::lib::types_list<symphas::internal::itype<Ns, OpVoid>...>>;
	//	return expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>..., 
	//		symphas::internal::itype<Is, empty_t>...>>{};
	//}

	//template<size_t I, typename... Es, size_t... Is>
	//auto fill_stencil_vector(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>...>>, std::index_sequence<Is...>)
	//{
	//	constexpr size_t N = symphas::lib::types_list_size<symphas::lib::type_at_index<0, Es...>>::value;
	//	return fill_stencil_vector(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>...>>{}, std::index_sequence<Is...>{}, std::make_index_sequence<N>{});
	//}

	//template<size_t I, typename... Es, size_t... Is, typename T0, typename... Ts>
	//auto fill_stencil_vector(expr::stencil_vector_type<2, I, symphas::lib::types_list<stencil_vector_type<1, I, Es>...>>, T0, Ts...)
	//{
	//	return expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>..., Ts...>>{};
	//}

	template<size_t N, size_t I, typename... Es, size_t... Is>
	auto fill_stencil_vector(expr::stencil_vector_type<N, I, symphas::lib::types_list<stencil_vector_type<N - 1, I, Es>...>>, std::index_sequence<Is...>)
	{
		using E0 = symphas::lib::type_at_index<0, stencil_vector_type<N - 1, I, Es>...>;
		return expr::stencil_vector_type<N, I, symphas::lib::types_list<stencil_vector_type<N - 1, I, Es>..., empty_stencil_vector_t<Is, sizeof...(Is), E0>...>>{};
	}


	template<typename V, expr::exp_key_t X, size_t... N0s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X>>, std::index_sequence<N0s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return stencil_vector_type<0, 0, symphas::lib::types_list<std::conditional_t<N == N0s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X, size_t... N0s, size_t... N1s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X>>, std::index_sequence<N0s...>, std::index_sequence<0, N1s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return 
			stencil_vector_type<0, 0, symphas::lib::types_list<std::conditional_t<N == N0s, V, OpVoid>...>>{}
			* stencil_vector_type<0, 1, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N1s, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X, size_t... N0s, size_t... N1s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return
			stencil_vector_type<0, 0, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N0s, OpVoid>...>>{}
			* stencil_vector_type<0, 1, symphas::lib::types_list<std::conditional_t<N == N1s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, size_t... N0s, size_t... N1s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::x_symbol, X0>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		return
			stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
			* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, size_t... N0s, size_t... N1s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::y_symbol, X1>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		return
			stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
			* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{};
	}


	template<typename V, expr::exp_key_t X, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X>>, std::index_sequence<N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<std::conditional_t<N == N0s, V, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N1s, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N2s, OpVoid>...>>{};
	}
	
	template<typename V, expr::exp_key_t X, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N0s, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<std::conditional_t<N == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N2s, OpVoid>...>>{};
	}
	
	template<typename V, expr::exp_key_t X, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::z_symbol, X>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<N2s...>)
	{
		constexpr size_t N = expr::_Xk_t<X>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N0s, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N1s, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<std::conditional_t<N == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::y_symbol, X1>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N2s, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::x_symbol, X0>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N2s, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::z_symbol, X2>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N1s, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::x_symbol, X0>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N1s, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::z_symbol, X2>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N0s, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, OpIdentity, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::y_symbol, X1>>, std::index_sequence<0, N0s...>, std::index_sequence<N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpIdentity, symphas::internal::itype<N0s, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, OpIdentity, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::z_symbol, X2>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::y_symbol, X1>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::z_symbol, X2>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::x_symbol, X0>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::x_symbol, X0>, Term<expr::symbols::y_symbol, X1>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<typename V, expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t X2, size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpTerms<V, Term<expr::symbols::z_symbol, X2>, Term<expr::symbols::y_symbol, X1>, Term<expr::symbols::x_symbol, X0>>, std::index_sequence<0, N0s...>, std::index_sequence<0, N1s...>, std::index_sequence<0, N2s...>)
	{
		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t N1 = expr::_Xk_t<X1>::N;
		constexpr size_t N2 = expr::_Xk_t<X2>::N;
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<OpVoid, std::conditional_t<N0 == N0s, OpIdentity, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<OpVoid, std::conditional_t<N1 == N1s, V, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<OpVoid, std::conditional_t<N2 == N2s, V, OpVoid>...>>{};
	}

	template<size_t... N0s>
	auto as_stencil_vector_impl(OpVoid, std::index_sequence<N0s...>)
	{
		return stencil_vector_type<0, 0, symphas::lib::types_list<symphas::internal::itype<N0s, OpVoid>...>>{};
	}

	template<size_t... N0s, size_t... N1s>
	auto as_stencil_vector_impl(OpVoid, std::index_sequence<N0s...>, std::index_sequence<N1s...>)
	{
		return
			stencil_vector_type<0, 0, symphas::lib::types_list<symphas::internal::itype<N0s, OpVoid>...>>{}
			* stencil_vector_type<0, 1, symphas::lib::types_list<symphas::internal::itype<N1s, OpVoid>...>>{};
	}

	template<size_t... N0s, size_t... N1s, size_t... N2s>
	auto as_stencil_vector_impl(OpVoid, std::index_sequence<N0s...>, std::index_sequence<N1s...>, std::index_sequence<N2s...>)
	{
		return
			(stencil_vector_type<0, 0, symphas::lib::types_list<symphas::internal::itype<N0s, OpVoid>...>>{}
				* stencil_vector_type<0, 1, symphas::lib::types_list<symphas::internal::itype<N1s, OpVoid>...>>{})
			* stencil_vector_type<0, 2, symphas::lib::types_list<symphas::internal::itype<N2s, OpVoid>...>>{};
	}

	template<typename... Es, typename... Seqs>
	auto as_stencil_vector_impl(OpAdd<Es...>, Seqs...)
	{
		return (as_stencil_vector_impl(Es{}, Seqs{}...) + ...);
	}

	template<typename E, size_t... Ns>
	auto as_stencil_vector(OpExpression<E>, std::index_sequence<Ns...>)
	{
		return as_stencil_vector_impl(E{}, std::make_index_sequence<Ns>{}...);
	}

	template<size_t I, typename... Es>
	auto unpack_stencil_vector(expr::stencil_vector_type<0, I, symphas::lib::types_list<Es...>>)
	{
		return std::make_tuple(Es{}...);
	}

	template<size_t I, typename... Es>
	auto unpack_stencil_vector(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>...>>)
	{
		return std::make_tuple(unpack_stencil_vector(stencil_vector_type<0, 0, Es>{})...);
	}

	template<size_t I, typename... Es>
	auto unpack_stencil_vector(expr::stencil_vector_type<1, I, symphas::lib::types_list<stencil_vector_type<0, 0, Es>&&...>>)
	{
		return std::make_tuple(unpack_stencil_vector(stencil_vector_type<0, 0, Es>{})...);
	}

	template<size_t I, typename... Es>
	auto unpack_stencil_vector(expr::stencil_vector_type<2, I, symphas::lib::types_list<stencil_vector_type<1, 0, Es>...>>)
	{
		return std::make_tuple(unpack_stencil_vector(stencil_vector_type<1, 0, Es>{})...);
	}

	template<typename... Es>
	auto unpack_stencil_vector(symphas::lib::types_list<Es...>)
	{
		return std::tuple_cat(unpack_stencil_vector(Es{})...);
	}
}

namespace symphas::internal
{
	using expr::stencil_vector_type;

	template<typename T>
	struct filter_zeros
	{
		using type = symphas::lib::types_list<T>;
	};

	template<>
	struct filter_zeros<OpVoid>
	{
		using type = symphas::lib::types_list<>;
	};

	template<typename... Ts>
	struct filter_zeros<symphas::lib::types_list<Ts...>>
	{
		using type = symphas::lib::expand_types_list<typename filter_zeros<Ts>::type...>;
	};

	template<size_t N, size_t I, typename... Es>
	struct filter_zeros<expr::stencil_vector_type<N, I, symphas::lib::types_list<Es...>>>
	{
		using type = symphas::lib::expand_types_list<typename filter_zeros<symphas::lib::types_list<Es...>>::type>;
	};

	template<typename T>
	using filter_zeros_t = typename filter_zeros<T>::type;

}


template<typename V, typename G, size_t I, size_t N, typename... Es>
auto operator*(OpTerm<V, G> const& term, expr::stencil_vector_type<N, I, symphas::lib::types_list<Es...>>)
{
	return expr::stencil_vector_type<N, I, symphas::lib::types_list<mul_result_t<OpTerm<V, G>, Es>...>>{};
}

template<size_t I, size_t N, typename... Es, typename E>
auto operator*(expr::stencil_vector_type<I, N, symphas::lib::types_list<Es...>>, OpExpression<E> const& term)
{
	return expr::stencil_vector_type<I, N, symphas::lib::types_list<mul_result_t<Es, E>...>>{};
}

template<size_t I, size_t N, typename... Es>
auto operator*(expr::stencil_vector_type<I, N, symphas::lib::types_list<Es...>>, OpVoid)
{
	return expr::stencil_vector_type<I, N, symphas::lib::types_list<symphas::internal::ttype<Es, OpVoid>...>>{};
}

template<size_t I, size_t N, typename... Es>
auto operator*(expr::stencil_vector_type<I, N, symphas::lib::types_list<Es...>>, OpIdentity)
{
	return expr::stencil_vector_type<I, N, symphas::lib::types_list<Es...>>{};
}

template<typename... E0s, typename... E1s>
auto operator*(expr::stencil_vector_type<0, 0, symphas::lib::types_list<E0s...>>,
	expr::stencil_vector_type<0, 1, symphas::lib::types_list<E1s...>>)
{
	return expr::stencil_vector_type<1, 0,
		symphas::lib::types_list<mul_result_t<expr::stencil_vector_type<0, 0, symphas::lib::types_list<E0s...>>, E1s>...>>{};
}

template<typename... E0s, typename... E1s>
auto operator*(expr::stencil_vector_type<1, 0, symphas::lib::types_list<E0s...>>,
	expr::stencil_vector_type<0, 2, symphas::lib::types_list<E1s...>>)
{
	return expr::stencil_vector_type<2, 0,
		symphas::lib::types_list<mul_result_t<expr::stencil_vector_type<1, 0, symphas::lib::types_list<E0s...>>, E1s>...>>{};
}

template<size_t N, size_t I, typename... E0s, typename... E1s>
auto operator+(expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>,
	expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>)
{
	if constexpr (sizeof...(E0s) == sizeof...(E1s))
	{
		//auto test0 = unpack_stencil_vector(expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>{});
		//auto test1 = unpack_stencil_vector(expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>{});
		//auto test = unpack_stencil_vector(expr::stencil_vector_type<N, I, symphas::lib::types_list<add_result_t<E0s, E1s>...>>{});
		return expr::stencil_vector_type<N, I, symphas::lib::types_list<add_result_t<E0s, E1s>...>>{};
	}
	else if constexpr (sizeof...(E0s) < sizeof...(E1s))
	{
		return expr::fill_stencil_vector(
				expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>{},
				std::make_index_sequence<sizeof...(E1s) - sizeof...(E0s)>{})
		+ expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>{};
	}
	else
	{
		return expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>{}
			+ expr::fill_stencil_vector(
				expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>{},
				std::make_index_sequence<sizeof...(E0s) - sizeof...(E1s)>{});
	}
}

template<size_t N, size_t I, typename... E0s, typename... E1s>
auto operator-(expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>,
	expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>)
{
	if constexpr (sizeof...(E0s) == sizeof...(E1s))
	{
		return expr::stencil_vector_type<N, I, symphas::lib::types_list<sub_result_t<E0s, E1s>...>>{};
	}
	else if constexpr (sizeof...(E0s) < sizeof...(E1s))
	{
		return expr::fill_stencil_vector(
			expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>{},
			std::make_index_sequence<sizeof...(E1s) - sizeof...(E0s)>{})
			- expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>{};
	}
	else
	{
		return expr::stencil_vector_type<N, I, symphas::lib::types_list<E0s...>>{}
		- expr::fill_stencil_vector(
			expr::stencil_vector_type<N, I, symphas::lib::types_list<E1s...>>{},
			std::make_index_sequence<sizeof...(E0s) - sizeof...(E1s)>{});
	}
}

namespace expr
{

	using namespace symphas::internal;

	template<size_t K, int... Is, typename... Es, typename E>
	auto setup_stencil_equation(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...>, E const& d_op)
	{
		using namespace expr::symbols;

		return ((Es{} * expr::stencil_vector(x{}, expr::make_integer<Is>() * h{}, std::make_index_sequence<K + 1>{})) + ...)
			+ expr::as_stencil_vector(-expr::apply_operators(d_op * (expr::pow<K>(x{}))), std::index_sequence<K + 1>{});
		//return ((Es{} * expr::pow<K>(x{} + expr::make_integer<Is>() * h{})) + ...)
		//	- expr::apply_operators(d_op * (expr::pow<K>(x{})));
	}

	template<size_t K, size_t L, int... Is, int... Js, typename... Es, typename E>
	auto setup_stencil_equation(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>, E const& d_op)
	{
		using namespace expr::symbols;
		return (
			(Es{}
				* expr::stencil_vector(x{}, expr::make_integer<Is>() * h{}, std::make_index_sequence<K + 1>{})
				* expr::stencil_vector(y{}, expr::make_integer<Js>() * h{}, std::make_index_sequence<L + 1>{})) + ...)
			+ expr::as_stencil_vector(-expr::apply_operators(d_op * (expr::pow<K>(x{}) * expr::pow<L>(y{}))), std::index_sequence<K + 1, L + 1>{});

		//return ((Es{} * expr::pow<K>(x{} + expr::make_integer<Is>() * h{}) * expr::pow<L>(y{} + expr::make_integer<Js>() * h{})) + ...)
		//	- expr::apply_operators(d_op * (expr::pow<K>(x{}) * expr::pow<L>(y{})));
	}


	template<size_t K, size_t L, size_t M, int... Is, int... Js, int... Ks, typename... Es, typename E>
	auto setup_stencil_equation(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...>, E const& d_op)
	{
		using namespace expr::symbols;
		return (
			(Es{}
				* (expr::stencil_vector(x{}, expr::make_integer<Is>() * h {}, std::make_index_sequence<K + 1>{})
					* expr::stencil_vector(y{}, expr::make_integer<Js>() * h {}, std::make_index_sequence<L + 1>{}))
				* expr::stencil_vector(z{}, expr::make_integer<Ks>() * h{}, std::make_index_sequence<M + 1>{})) + ...)
			+ expr::as_stencil_vector(-expr::apply_operators(d_op * (expr::pow<K>(x{}) * expr::pow<L>(y{}) * expr::pow<M>(z{}))), std::index_sequence<K + 1, L + 1, M + 1>{});

		//return ((Es{} * expr::pow<K>(x{} + expr::make_integer<Is>() * h{}) * expr::pow<L>(y{} + expr::make_integer<Js>() * h{}) * expr::pow<M>(z{} + expr::make_integer<Ks>() * h{})) + ...)
		//	- expr::apply_operators(d_op * (expr::pow<K>(x{}) * expr::pow<L>(y{}) * expr::pow<M>(z{})));
	}

	template<typename... Symbols>
	auto make_central_stencil_dictionary(types_list<Symbols...>)
	{
		return types_list<std::pair<Symbols, OpTerm<OpIdentity, Symbols>>...>{};
	}


	template<size_t D, int R>
	auto make_central_stencil_dictionary()
	{
		return make_central_stencil_dictionary(get_stencil_symbols<D, R>());
	}

	// *************************************************************************************************
	// Simplifies a stencil dictionary based on the symmetry of a central difference operator.
	// In particular, for even-ordered derivatives, applies "rotational" symmetry, and for
	// odd-ordered derivatives, makes the right side equal to the negative of the left side.
	// *************************************************************************************************

	template<size_t O, int... Is>
	auto simplify_central_stencils(
		types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, OpTerm<OpIdentity, expr::symbols::internal::S1_symbol<Is>>>...> const& dict)
	{
		using expression_first_half_t = types_between_index<0, sizeof...(Is) / 2, OpTerm<OpIdentity, expr::symbols::internal::S1_symbol<Is>>...>;
		using expression_at_half_t = type_at_index<sizeof...(Is) / 2, OpTerm<OpIdentity, expr::symbols::internal::S1_symbol<Is>>...>;

		if constexpr (O == 0)
		{
			return dict;
		}
		else if constexpr (O % 2 == 0)
		{
			return typename pack_dictionary<
				expand_types_list<
					expression_first_half_t,
					types_list<expression_at_half_t>,
					reverse_types_list<expression_first_half_t>>,
				std::integer_sequence<int, Is...>
			>::type{};
		}
		else if constexpr (O % 2 == 1)
		{
			return typename pack_dictionary<
				expand_types_list<
					expression_first_half_t,
					types_list<OpVoid>,
					decltype(neg_elements(reverse_types_list<expression_first_half_t>{}))>,
				std::integer_sequence<int, Is...>
			>::type{};

		}
	}

	template<size_t O, int... Is, int... Js>
	auto simplify_central_stencils(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>>...> const& dict)
	{
		return keep_axes_2d<O>(
			types_list<OpTerm<OpIdentity, expr::symbols::internal::S2_symbol<Is, Js>>...>{},
			std::integer_sequence<int, Is...>{}, std::integer_sequence<int, Js...>{});
	}

	template<size_t Oy, size_t N, int... Is, int... Js, typename... Es>
	auto simplify_central_mixed_stencils(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict)
	{
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int R = fixed_max<JM - Jm + 1, IM - Im + 1>;

		constexpr int RJ_2 = R_<Oy, N>;

		if constexpr (Oy == 0)
		{
			return typename pack_dictionary<types_list<
				std::conditional_t<(test_in_radius<RJ_2 - 1, Js>), Es, OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type{};
		}
		else if constexpr (Oy % 2 == 1)
		{
			using Ts = types_list<Es...>;
			using nTs = types_list<std::invoke_result_t<decltype(&neg_one_element<Es>), Es>...>;

			return typename pack_dictionary<types_list<
				std::conditional_t<
					((Js != 0 && test_in_radius<RJ_2, Is, Js>) || (!test_in_radius<R / 2, Is, Js> && test_in_radius<R / 2 + 1, Is, Js>)),
					typename at_stencil_index_2d<
						Is, -fixed_abs<Js>,
						std::conditional_t<(Js < 0), Ts, nTs>, 
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js... >> ::type{};
		}
		else
		{
			return typename pack_dictionary<types_list<
				std::conditional_t<
					((test_in_radius<RJ_2 - NEG_RADIUS_CENTRAL_2D, Is, Js> || fixed_abs<Js> == fixed_abs<Is>)),
					typename at_stencil_index_2d<
						(Is != 0 && Js != 0) ? fixed_sign<Is> * fixed_min<fixed_abs<Js>, fixed_abs<Is>> : Is,
						(Is != 0 && Js != 0) ? -fixed_max<fixed_abs<Is>, fixed_abs<Js>> : -fixed_abs<Js>,
						types_list<Es...>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type{};
		}
	}

	template<size_t Ox, size_t Oy, size_t N, int... Is, int... Js, typename... Es>
	auto simplify_central_mixed_stencils(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict)
	{
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int R = fixed_max<JM - Jm + 1, IM - Im + 1>;

		constexpr int RI_2 = R_<Ox, N>;

		if constexpr (Ox == 0)
		{
			using dictx = typename pack_dictionary<
				types_list<std::conditional_t<(test_in_radius<RI_2 - 1, Is>), Es, OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type;
			return simplify_central_mixed_stencils<Oy, N>(dictx{});
		}
		else if constexpr (Ox % 2 == 1)
		{
			using Ts = types_list<Es...>;
			using nTs = types_list<std::invoke_result_t<decltype(&neg_one_element<Es>), Es>...>;

			using dictx = typename pack_dictionary<types_list<
				std::conditional_t<
					((Is != 0 && test_in_radius<RI_2, Is, Js>) || (!test_in_radius<R / 2, Is, Js> && test_in_radius<R / 2 + 1, Is, Js>)),
					typename at_stencil_index_2d<
						-fixed_abs<Is>, Js,
						std::conditional_t<(Is < 0), Ts, nTs>, 
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type;
			return simplify_central_mixed_stencils<Oy, N>(dictx{});
		}
		else
		{
			using dictx = typename pack_dictionary<types_list<
				std::conditional_t<
					((test_in_radius<RI_2 - NEG_RADIUS_CENTRAL_2D, Is, Js> || fixed_abs<Js> == fixed_abs<Is>)),
					typename at_stencil_index_2d<
						(Is != 0 && Js != 0) ? -fixed_max<fixed_abs<Js>, fixed_abs<Is>> : -fixed_abs<Is>,
						(Is != 0 && Js != 0) ? fixed_sign<Js> * fixed_min<fixed_abs<Js>, fixed_abs<Is>> : Js,
						types_list<Es...>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type;
			return simplify_central_mixed_stencils<Oy, N>(dictx{});
		}

	}


	
	template<size_t Oz, size_t N, int... Is, int... Js, int... Ks, typename... Es>
	auto simplify_central_mixed_stencils(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict)
	{
		constexpr int Km = fixed_min<Ks...>;
		constexpr int KM = fixed_max<Ks...>;
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int R = fixed_max<JM - Jm + 1, IM - Im + 1, KM - Km + 1>;

		constexpr int RK_2 = R_<Oz, N>;

		if constexpr (Oz == 0)
		{
			return typename pack_dictionary<types_list<
				std::conditional_t<(test_in_radius<RK_2 - 1, Ks>), Es, OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type{};
		}
		else if constexpr (Oz % 2 == 1)
		{
			using Ts = types_list<Es...>;
			using nTs = types_list<std::invoke_result_t<decltype(&neg_one_element<Es>), Es>...>;

			return typename pack_dictionary<types_list<
				std::conditional_t<
					((Ks != 0 && test_in_radius<RK_2, Ks>) || !test_in_radius<R / 2, Is, Js, Ks>),
					typename at_stencil_index_3d<
						Is, Js, -fixed_abs<Ks>,
						std::conditional_t<(Ks < 0), Ts, nTs>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js... >, std::integer_sequence<int, Ks...>>::type{};
		}
		else
		{
			return typename pack_dictionary<types_list<
				std::conditional_t<
					((test_in_radius<R / 2 - 1 - NEG_RADIUS_CENTRAL_3D, Is, Js, Ks> || fixed_abs<Ks> == fixed_abs<Js> || fixed_abs<Ks> == fixed_abs<Is>)),
					typename at_stencil_index_3d<
						Is, Js, -fixed_abs<Ks>,
						types_list<Es...>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type{};
		}
	}

	template<size_t Oy, size_t Oz, size_t N, int... Is, int... Js, int... Ks, typename... Es>
	auto simplify_central_mixed_stencils(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict)
	{
		constexpr int Km = fixed_min<Ks...>;
		constexpr int KM = fixed_max<Ks...>;
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int R = fixed_max<JM - Jm + 1, IM - Im + 1, KM - Km + 1>;

		constexpr int RJ_2 = R_<Oy, N>;

		if constexpr (Oy == 0)
		{
			using dicty = typename pack_dictionary<types_list<
				std::conditional_t<(test_in_radius<RJ_2 - 1, Js>), Es, OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oz, N>(dicty{});
		}
		else if constexpr (Oy % 2 == 1)
		{
			using Ts = types_list<Es...>;
			using nTs = types_list<std::invoke_result_t<decltype(&neg_one_element<Es>), Es>...>;

			using dicty = typename pack_dictionary<types_list<
				std::conditional_t<
					((Js != 0 && test_in_radius<RJ_2, Js>) || !test_in_radius<R / 2, Is, Js, Ks>),
					typename at_stencil_index_3d<
						Is, -fixed_abs<Js>, Ks,
						std::conditional_t<(Js < 0), Ts, nTs>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oz, N>(dicty{});
		}
		else
		{
			using dicty = typename pack_dictionary<types_list<
				std::conditional_t<
					((test_in_radius<R / 2 - 1 - NEG_RADIUS_CENTRAL_3D, Is, Js, Ks> || fixed_abs<Js> == fixed_abs<Is> || fixed_abs<Js> == fixed_abs<Ks>)),
					typename at_stencil_index_3d<
						Is, -fixed_abs<Js>, Ks,
						types_list<Es...>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oz, N>(dicty{});
		}
	}

	template<size_t Ox, size_t Oy, size_t Oz, size_t N, int... Is, int... Js, int... Ks, typename... Es>
	auto simplify_central_mixed_stencils(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict)
	{
		constexpr int Jm = fixed_min<Js...>;
		constexpr int JM = fixed_max<Js...>;
		constexpr int Im = fixed_min<Is...>;
		constexpr int IM = fixed_max<Is...>;
		constexpr int R = fixed_max<JM - Jm + 1, IM - Im + 1>;

		constexpr int RI_2 = R_<Ox, N>;

		if constexpr (Ox == 0)
		{
			using dictx = typename pack_dictionary<
				types_list<std::conditional_t<(test_in_radius<RI_2 - 1, Is>), Es, OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oy, Oz, N>(dictx{});
		}
		else if constexpr (Ox % 2 == 1)
		{
			using Ts = types_list<Es...>;
			using nTs = types_list<std::invoke_result_t<decltype(&neg_one_element<Es>), Es>...>;

			using dictx = typename pack_dictionary<types_list<
				std::conditional_t<
					((Is != 0 && test_in_radius<RI_2, Is>) || !test_in_radius<R / 2, Is, Js, Ks>),
					typename at_stencil_index_3d<
						-fixed_abs<Is>, Js, Ks,
						std::conditional_t<(Is < 0), Ts, nTs>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oy, Oz, N>(dictx{});
		}
		else
		{
			using dictx = typename pack_dictionary<types_list<
				std::conditional_t<
					((test_in_radius<R / 2 - 1 - NEG_RADIUS_CENTRAL_3D, Is, Js, Ks> || fixed_abs<Is> == fixed_abs<Js> || fixed_abs<Is> == fixed_abs<Ks>)),
					typename at_stencil_index_3d<
						-fixed_abs<Is>, Js, Ks,
						types_list<Es...>,
						std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type,
					OpVoid>...>,
				std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type;
			return simplify_central_mixed_stencils<Oy, Oz, N>(dictx{});
		}

	}


	template<size_t O, int... Is, int... Js, int... Ks>
	auto simplify_central_stencils(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>>...> const& dict)
	{
		return keep_axes_3d<O>(
				types_list<OpTerm<OpIdentity, expr::symbols::internal::S3_symbol<Is, Js, Ks>>...>{},
				std::integer_sequence<int, Is...>{}, std::integer_sequence<int, Js...>{}, std::integer_sequence<int, Ks...>{});
	}

	template<typename... Ss>
	auto simplify_central_stencils(types_list<Ss...> const& dict)
	{
		return simplify_central_stencils<0>(dict);
	}

	// *************************************************************************************************
	// `update_stencil_dictionary`: updates the dictionary of stencil coefficients using equations generated
	// from the equation setup phase. `update_stencil_dictionary` is applied
	// for every equation that it is passed. On each pass, the given equation has substitutions
	// made, and then the dictionary is updated by isolating for the first unknown symbol. The 
	// updated dictionary is used on subsequent calls.
	// *************************************************************************************************

	template<typename Symbol, typename E>
	auto get_term_with_symbol(OpExpression<E>)
	{
		return OpVoid{};
	}

	template<typename Symbol, typename T, typename... Gs, expr::exp_key_t... Xs>
	auto get_term_with_symbol(OpTerms<T, Term<Gs, Xs>...>)
	{
		if constexpr (symphas::lib::index_of_type<Symbol, Gs...> >= 0)
		{
			return OpTerms<T, Term<Gs, Xs>...>{};
		}
		else
		{
			return OpVoid{};
		}
	}

	template<typename Symbol, typename... Es>
	auto get_term_with_symbol(OpAdd<Es...>);

	template<typename Symbol, typename A, typename B>
	auto get_term_with_symbol(OpBinaryDiv<A, B>)
	{
		return expr::make_div(get_term_with_symbol<Symbol>(A{}), B{});
	}

	template<typename Symbol, typename... Es>
	auto get_term_with_symbol(OpAdd<Es...>)
	{
		return (get_term_with_symbol<Symbol>(Es{}) + ...);
	}

	template<size_t N, typename Ts>
	struct symbol_match
	{
		static const size_t value = N;
	};

	template<size_t N, typename... Ts, int I>
	struct symbol_match<N, types_list<expr::symbols::internal::S1_symbol<I>, Ts...>>
	{
		static const size_t value = N;
	};

	template<size_t N, typename... Ts, int I, int J>
	struct symbol_match<N, types_list<expr::symbols::internal::S2_symbol<I, J>, Ts...>>
	{
		static const size_t value = N;
	};

	template<size_t N, typename... Ts, int I, int J, int K>
	struct symbol_match<N, types_list<expr::symbols::internal::S3_symbol<I, J, K>, Ts...>>
	{
		static const size_t value = N;
	};

	template<size_t N, typename T0, typename... Ts>
	struct symbol_match<N, types_list<T0, Ts...>>
	{
		static const size_t value = symbol_match<N + 1, types_list<Ts...>>::value;
	};

	template<typename... Symbols, typename... Es>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, OpVoid)
	{
		return dict;
	}

	template<typename... Symbols, typename... Es, typename E0>//, typename F = first_undefined_symbol<types_list<std::pair<Symbols, Es>...>>>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, OpExpression<E0> const& expr)
	{
		auto substitution = make_all_substitutions(dict, E0{});

		using ops_t = typename expr::op_types<decltype(substitution)>::type;
		constexpr size_t I = symbol_match<0, ops_t>::value;

		if constexpr (I < symphas::lib::types_list_size<ops_t>::value)
		{
			using symbol_t = type_at_index<I, unroll_types_list<ops_t>>;

			auto op = OpTerm<OpIdentity, symbol_t>{};
			auto term = get_term_with_symbol<symbol_t>(substitution);
			auto rhs = -(substitution - term) / (term / op);

			return back_substitution<symbol_t>(dict, rhs);
		}
		else
		{
			return dict;
		}
	}

	//template<int... Is, int... Js, int... Ks, typename... Es>
	//auto update_stencil_dictionary(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict, types_list<> const& exprs)
	//{
	//	return typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>, std::integer_sequence<int, Ks...>>::type{};
	//}

	//template<int... Is, int... Js, typename... Es>
	//auto update_stencil_dictionary(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict, types_list<> const& exprs)
	//{
	//	return typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>, std::integer_sequence<int, Js...>>::type{};
	//}

	//template<int... Is, typename... Es>
	//auto update_stencil_dictionary(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...> const& dict, types_list<> const& exprs)
	//{
	//	return typename pack_dictionary<types_list<Es...>, std::integer_sequence<int, Is...>>::type{};
	//}

	template<typename... Symbols, typename... Es, typename E0, typename... E0s>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, types_list<E0, E0s...> const& exprs);

	template<typename... Symbols, typename... Es>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, types_list<> const& exprs)
	{
		return self_complete(dict);
	}

	template<typename... Symbols, typename... Es, typename... E0s>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, types_list<OpVoid, E0s...> const& exprs)
	{
		return update_stencil_dictionary(dict, types_list<E0s...>{});
	}

	template<typename... Symbols, typename... Es, typename E0, typename... E0s>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict, types_list<E0, E0s...> const& exprs)
	{
		if constexpr (!dict_complete<types_list<std::pair<Symbols, Es>...>>::value)
		{
			return update_stencil_dictionary(update_stencil_dictionary(dict, E0{}), types_list<E0s...>{});
		}
		else
		{
			return dict;
		}
	}


	template<typename... Symbols, typename... Es, typename... E0s>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict,
		stencil_vector_type<0, 0, symphas::lib::types_list<E0s...>>)
	{
		return update_stencil_dictionary(dict, symphas::lib::types_list<E0s...>{});
	}

	template<typename... Symbols, typename... Es, size_t N, typename... E0s>
	auto update_stencil_dictionary(types_list<std::pair<Symbols, Es>...> const& dict,
		stencil_vector_type<N, 0, symphas::lib::types_list<stencil_vector_type<N - 1, 0, E0s>...>>)
	{
		return update_stencil_dictionary(dict, symphas::lib::types_list<stencil_vector_type<N - 1, 0, E0s>...>{});
	}

}

namespace symphas::internal
{

	template<size_t Q, size_t P, size_t L, int... Is, int... Js, int... Ks, typename... Es, typename E>
	auto setup_stencil_equation_permute(
		types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict,
		E const& d_op)
	{
		return types_list<
			decltype(expr::setup_stencil_equation<Q, P, L>(dict, d_op)),
			decltype(expr::setup_stencil_equation<L, Q, P>(dict, d_op)),
			decltype(expr::setup_stencil_equation<P, L, Q>(dict, d_op))>{};
	}


	template<size_t Q0, size_t Q, size_t P, size_t L, int... Is, int... Js, int... Ks, typename... Es, typename E>
	auto setup_stencil_equation_expand_2(
		types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict,
		E const& d_op)
	{
		if constexpr (Q == Q0 / 2 && Q0 % 2 == 0)
		{
			return types_list<>{};
		}
		else
		{
			return setup_stencil_equation_permute<Q, P, L>(dict, d_op);
		}

	}

	template<size_t Q0, size_t Q, size_t P, size_t... Ps, int... Is, int... Js, int... Ks, typename... Es, typename E>
	auto setup_stencil_equation_expand_1(
		types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict,
		E const& d_op, std::index_sequence<Ps...>)
	{
		return expand_types_list<decltype(setup_stencil_equation_expand_2<Q0, Q, P - Ps, Ps>(dict, d_op))...>{};
	}

}

namespace expr
{

	template<size_t Ox>
	auto get_mixed_derivative()
	{
		return expr::make_operator_derivative<Ox, expr::symbols::x_symbol>();
	}

	template<size_t Ox, size_t Oy>
	auto get_mixed_derivative()
	{
		return expr::make_operator_derivative<Ox, expr::symbols::x_symbol>() *
			expr::make_operator_derivative<Oy, expr::symbols::y_symbol>();
	}

	template<size_t Ox, size_t Oy, size_t Oz>
	auto get_mixed_derivative()
	{
		return expr::make_operator_derivative<Ox, expr::symbols::x_symbol>() *
			(expr::make_operator_derivative<Oy, expr::symbols::y_symbol>() *
				expr::make_operator_derivative<Oz, expr::symbols::z_symbol>());
	}

	//! Separates a list of equations based on the order of the x and y variables.
	template<size_t Q, int... Is, int... Js, typename... Es, size_t... Ls, typename E>
	auto setup_stencil_equation_list(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict,
		E const& d_op, std::index_sequence<Ls...>)
	{
		using symphas::internal::split_by;
		using symphas::internal::reverse_types_list;
		
		using expr_types = types_list<decltype(expr::setup_stencil_equation<Q - Ls, Ls>(dict, d_op))...>;
		return symphas::internal::filter_zeros_t<expr_types>{};
	}

	//! Separates a list of equations based on the order of the x and y variables.
	template<size_t Q, int R, int... Is, int... Js, int... Ks, typename... Es, size_t... Ls, typename E>
	auto setup_stencil_equation_list(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict,
		E const& d_op, std::index_sequence<Ls...>)
	{
		using symphas::internal::split_by;

		//auto exprs = std::make_tuple(expr::unpack_stencil_vector(symphas::internal::setup_stencil_equation_expand_1<Q, Q - Ls, Ls>(dict, d_op, std::make_index_sequence<Ls>{}))...);
		using expr_types = types_list<decltype(
			symphas::internal::setup_stencil_equation_expand_1<Q, Q - Ls, Ls>(dict, d_op, std::make_index_sequence<Ls>{}))...>;
		return symphas::internal::filter_zeros_t<expr_types>{};
	}


	// *************************************************************************************************
	// `get_central_space_stencil`: entry point for setting up the equations and recursively updating the 
	// dictionary based on those equations.
	// *************************************************************************************************

	template<size_t Q, size_t N, int R, int... Is,
		typename... Es, typename E, size_t O = expr::derivative_order<E>::value>
	auto construct_stencil(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...> const& dict, E const& d_op)
	{
		//expr::printe(setup_stencil_equation_<Q>(dict, d_op));
		//auto terms = expr::unpack_stencil_vector(setup_stencil_equation<Q>(dict, d_op));
		return update_stencil_dictionary(dict, setup_stencil_equation<Q>(dict, d_op));
	}

	template<size_t Q, size_t N, int R, int... Is, int... Js, typename... Es,
		typename E, size_t O = expr::derivative_order<E>::value>
	auto construct_stencil(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> const& dict, E const& d_op)
	{
		using symphas::internal::split_by;
		using symphas::internal::reverse_types_list;

		constexpr size_t QQ = (N % 2 == 0) ? Q / 2 + 1 : Q + 1;
		return update_stencil_dictionary(dict, setup_stencil_equation_list<Q>(dict, d_op, std::make_index_sequence<QQ>{}));
	}


	template<size_t Q, size_t N, int R, int... Is, int... Js, int... Ks, typename... Es,
		typename E, size_t O = expr::derivative_order<E>::value>
	auto construct_stencil(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...> const& dict, E const& d_op)
	{
		constexpr size_t QQ = (N % 2 == 0) ? Q / 2 + 1 : Q + 1;
		return update_stencil_dictionary(dict, setup_stencil_equation_list<Q, R>(dict, d_op, std::make_index_sequence<QQ>{}));
	}

	template<size_t O1, size_t O2>
	using mixed_deriv_2d_t = std::invoke_result_t<decltype(&get_mixed_derivative<O1, O2>)>;

	template<int I, int J, size_t O, typename T>
	using dict_entry_2d_t = std::pair<expr::symbols::internal::S2_symbol<I, J>, decltype(T{} / expr::pow<O>(expr::symbols::h{}))>;

	//template<size_t N, size_t D, typename E, size_t O = expr::derivative_order<E>::value>
	//auto get_central_space_stencil(E const& d_op)
	//{
	//	constexpr int R = symphas::internal::R_<O, N>;

	//	auto dict = make_central_stencil_dictionary<D, R>();

	//	constexpr size_t Q = O + N - 1;
	//	return get_central_space_stencil<Q, N, R>(dict, d_op);
	//}

	template<size_t N, size_t D, typename E, typename... Ss, size_t O = expr::derivative_order<E>::value>
	auto get_stencil(E const& d_op, types_list<Ss...> const& dict)
	{
		constexpr int R = symphas::internal::R_<O, N>;

		constexpr size_t Q = O + N - 1;
		auto dict0 = construct_stencil<Q, N, R>(dict, d_op);
		//if constexpr (!dict_complete<decltype(dict0)>::value)
		//{
		//	return construct_stencil<Q + 1, N, R>(dict0, d_op);
		//}
		//else
		{
			return dict0;
		}
	}

	template<size_t O, size_t N, size_t D>
	auto get_central_space_stencil()
	{
		constexpr int R = symphas::internal::R_<O, N>;
		auto dict = simplify_central_stencils<O>(make_central_stencil_dictionary<D, R>());
		
		if constexpr (D == 1)
		{
			return get_stencil<N, 1>(expr::make_operator_derivative<O, expr::symbols::x_symbol>(), dict);
		}
		else if constexpr (D == 2)
		{
			if constexpr (O == 1)
			{
				return get_stencil<N, 2>(expr::make_operator_derivative<1, expr::symbols::x_symbol>(), dict);
			}
			if constexpr (O % 2 == 0)
			{
				return get_stencil<N, 2>(
					expr::pow<O / 2>(
						expr::make_operator_derivative<2, expr::symbols::x_symbol>() + expr::make_operator_derivative<2, expr::symbols::y_symbol>()), dict);
			}
			else
			{
				return get_stencil<N, 2>(
					expr::make_operator_derivative<1, expr::symbols::x_symbol>()
					* expr::pow<O / 2>(expr::make_operator_derivative<2, expr::symbols::x_symbol>()
						+ expr::make_operator_derivative<2, expr::symbols::y_symbol>()), dict);
			}
		}
		else if constexpr (D == 3)
		{
			if constexpr (O == 1)
			{
				return get_stencil<N, 3>(expr::make_operator_derivative<1, expr::symbols::x>(), dict);
			}
			if constexpr (O % 2 == 0)
			{
				return get_stencil<N, 3>(
					expr::pow<O / 2>(
						expr::make_operator_derivative<2, expr::symbols::x_symbol>()
						+ expr::make_operator_derivative<2, expr::symbols::y_symbol>()
						+ expr::make_operator_derivative<2, expr::symbols::z_symbol>()), dict);
			}
			else
			{
				return get_stencil<N, 3>(
					expr::make_operator_derivative<1, expr::symbols::x_symbol>()
					* expr::pow<O / 2>(
						expr::make_operator_derivative<2, expr::symbols::x_symbol>()
						+ expr::make_operator_derivative<2, expr::symbols::y_symbol>()
						+ expr::make_operator_derivative<2, expr::symbols::z_symbol>()), dict);
			}
		}
	}


	template<size_t N, size_t O1, size_t O2>
	auto get_central_space_mixed_stencil(std::index_sequence<O1, O2>)
	{
		constexpr int R = symphas::internal::R_<fixed_max<O1, O2>, N>;
		auto dict = simplify_central_mixed_stencils<O1, O2, N>(make_central_stencil_dictionary<2, R>());
		return get_stencil<N, 2>(expr::get_mixed_derivative<O1, O2>(), dict);
	}

	template<size_t N, size_t O1, size_t O2, size_t O3>
	auto get_central_space_mixed_stencil(std::index_sequence<O1, O2, O3>)
	{
		constexpr int R = symphas::internal::R_<fixed_max<O1, O2, O3>, N>;
		auto dict = simplify_central_mixed_stencils<O1, O2, O3, N>(make_central_stencil_dictionary<3, R>());
		return get_stencil<N, 3>(expr::get_mixed_derivative<O1, O2, O3>(), dict);
	}

	template<Axis ax, size_t O, size_t N, size_t D>
	auto get_central_space_directional_stencil()
	{
		constexpr int R = symphas::internal::R_<O, N>;
		if constexpr (D == 1)
		{
			return get_central_space_stencil<O, N, 1>();
		}
		if constexpr (D == 2)
		{
			return get_central_space_mixed_stencil<N>(std::index_sequence<(ax == Axis::X) ? O : 0, (ax == Axis::Y) ? O : 0>{});
		}
		else
		{
			return get_central_space_mixed_stencil<N>(std::index_sequence<(ax == Axis::X) ? O : 0, (ax == Axis::Y) ? O : 0, (ax == Axis::Z) ? O : 0>{});
		}
	}

}



namespace symphas::internal
{



	// *************************************************************************************************
	// Applies the result of `get_*_stencil` (which is in the form of the stencil dictionary)
	// and can either print it or apply it as a stencil to a grid.
	// *************************************************************************************************

	template<typename T>
	struct StencilCoeff;


#define BOX_SPACING 26
#define ENTRY_SPACING 9

	template<typename S, typename E>
	struct StencilCoeff<std::pair<S, E>>
	{
		static const size_t h_exponent = 0;

		template<typename T, size_t D>
		constexpr auto operator()(T const* v, const len_type(&stride)[D])
		{
			return E{}.eval(0);
		}

		size_t print(char* out)
		{
			size_t n = 0;
			n += E{}.print(out);
			n += OpTerm<OpIdentity, S>{}.print(out + n);
			return n;
		}

		size_t print(FILE* out)
		{
			size_t n = 0;
			n += E{}.print(out);
			n += OpTerm<OpIdentity, S>{}.print(out);
			return n;
		}

		size_t print_length()
		{
			return E{}.print_length() + OpTerm<OpIdentity, S>{}.print_length();
		}
	};

	template<typename S>
	struct StencilCoeff<std::pair<S, OpVoid>>
	{
		static const size_t h_exponent = 0;

		template<typename T>
		constexpr auto operator()(T const* v, len_type)
		{
			return OpVoid{}.eval();
		}

		template<typename T, size_t D>
		constexpr auto operator()(T const* v, const len_type(&stride)[D])
		{
			return OpVoid{}.eval();
		}

		size_t print(char* out)
		{
			size_t n = 0;
			n += OpVoid{}.print(out);
			n += OpTerm<OpIdentity, S>{}.print(out + n);
			return n;
		}

		size_t print(FILE* out)
		{
			size_t n = 0;
			n += OpVoid{}.print(out);
			n += OpTerm<OpIdentity, S>{}.print(out);
			return n;
		}

		size_t print_length()
		{
			return OpVoid{}.print_length() + OpTerm<OpIdentity, S>{}.print_length();
		}
	};

	template<int N, typename Nt, typename Dt, expr::exp_key_t P>
	struct StencilCoeff<std::pair<
		expr::symbols::internal::S1_symbol<N>,
		OpBinaryDiv<Nt, OpTerms<Dt, Term<expr::symbols::h_symbol, P>>>>>
	{

		static const size_t h_exponent = P;

		template<typename T>
		auto operator()(T const* v, len_type stride)
		{
			return (Nt{} / Dt{}).eval() * v[N * stride];
		}

		size_t print(char* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S1<N>{}.print(out + n);
			return n;
		}

		size_t print(FILE* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S1<N>{}.print(out);
			return n;
		}

		size_t print_length()
		{
			return OpBinaryDiv<Nt, Dt>{}.print_length() + expr::symbols::internal::S1<N>{}.print_length();
		}
	};

	template<int N0, int N1, typename Nt, typename Dt, expr::exp_key_t P>
	struct StencilCoeff<std::pair<
		expr::symbols::internal::S2_symbol<N0, N1>,
		OpBinaryDiv<Nt, OpTerms<Dt, Term<expr::symbols::h_symbol, P>>>>>
	{
		static const size_t h_exponent = P;

		template<typename T>
		auto operator()(T const* v, len_type const (&stride)[2])
		{
			return (Nt{} / Dt{}).eval() * v[stride[0] * N0 + stride[1] * N1];
		}

		size_t print(char* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S2<N0, N1>{}.print(out + n);
			return n;
		}

		size_t print(FILE* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S2<N0, N1>{}.print(out);
			return n;
		}

		size_t print_length()
		{
			return OpBinaryDiv<Nt, Dt>{}.print_length() + expr::symbols::internal::S2<N0, N1>{}.print_length();
		}
	};

	template<int N0, int N1, int N2, typename Nt, typename Dt, expr::exp_key_t P>
	struct StencilCoeff<std::pair<
		expr::symbols::internal::S3_symbol<N0, N1, N2>,
		OpBinaryDiv<Nt, OpTerms<Dt, Term<expr::symbols::h_symbol, P>>>>>
	{
		static const size_t h_exponent = P;

		template<typename T>
		auto operator()(T const* v, len_type const (&stride)[3])
		{
			return (Nt{} / Dt{}).eval() * v[stride[0] * N0 + stride[1] * N1 + stride[2] * N2];
		}

		size_t print(char* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S3<N0, N1, N2>{}.print(out + n);
			return n;
		}

		size_t print(FILE* out)
		{
			size_t n = 0;
			n += OpBinaryDiv<Nt, Dt>{}.print(out);
			n += expr::symbols::internal::S3<N0, N1, N2>{}.print(out);
			return n;
		}

		size_t print_length()
		{
			return OpBinaryDiv<Nt, Dt>{}.print_length() + expr::symbols::internal::S3<N0, N1, N2>{}.print_length();
		}
	};

	template<typename T>
	struct GeneratedStencilApply;

	template<int Im, int Jm, int I>
	void print_stencil(types_list<>)
	{
		return printf("|\n");
	}

	template<int Im, int Jm, int I, int I0, int J0, typename E0, int... Is, int... Js, typename... Es>
	void print_stencil(types_list<std::pair<expr::symbols::internal::S2_symbol<I0, J0>, E0>, std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> content)
	{
		using entry_t = StencilCoeff<std::pair<expr::symbols::internal::S2_symbol<I0, J0>, E0>>;

		size_t n = 0;
		if constexpr (I0 < I && J0 != Jm)
		{
			n += printf("|\n");
			n += print_stencil<Im, Jm, Im>(content);
			return n;
		}
		else
		{
			if constexpr (I0 == I)
			{
				char* buffer = nullptr;
				size_t len = entry_t{}.print_length();
				buffer = new char[len + 1];
				entry_t{}.print(buffer);
				n += printf("| %*s ", BOX_SPACING, buffer);
				delete[] buffer;
				n += print_stencil<Im, Jm, I + 1>(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>{});
				return n;
			}
			else
			{
				n += printf("| %*s ", BOX_SPACING, "");
				n += print_stencil<Im, Jm, I + 1>(content);
				return n;
			}
		}
	}

	template<int... Is, int... Js, typename... Es>
	void print_stencil(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...> content)
	{
		return print_stencil<fixed_min<Is...>, fixed_min<Js...>, fixed_min<Is...>>(content);
	}

	template<int... Is, int... Js, typename... Es>
	void print_stencil(GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>>)
	{
		return print_stencil(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>{});
	}


	template<typename E0>
	struct coeff_type_divh { using type = E0; };
	template<>
	struct coeff_type_divh<OpVoid> { using type = OpVoid; };
	template<typename A, typename B>
	struct coeff_type_divh<OpBinaryDiv<A, B>> { using type = B; };


	template<>
	struct GeneratedStencilApply<types_list<>>
	{
		GeneratedStencilApply(...) {}

		size_t print(char* out) const
		{
			return sprintf(out, "x");
		}

		size_t print(FILE* out) const
		{
			return fprintf(out, "x");
		}
	};

	template<int... Is, typename... Es>
	struct GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...>>
	{
		static const size_t h_exponent = fixed_max<expr::factor_count<expr::symbols::h_symbol, typename coeff_type_divh<Es>::type>::value...>;

		GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...>) {}

		template<typename T>
		auto operator()(T const* v, len_type stride, double divh) const
		{
			return std::pow(divh, h_exponent) * (StencilCoeff<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>>{}(v, stride) + ...);
		}

		size_t print(char* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>>{}.print(out + n)),
				(n += sprintf(out + n, "\n"))), ...);
			return n;
		}

		size_t print(FILE* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>>{}.print(out)),
				(n += fprintf(out, "\n"))), ...);
			return n;
		}
	};

	template<int... Is, int... Js, typename... Es>
	struct GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>>
	{
		static const size_t h_exponent = fixed_max<expr::factor_count<expr::symbols::h_symbol, typename coeff_type_divh<Es>::type>::value...>;

		GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>) {}

		template<typename T>
		auto operator()(T const* v, len_type const (&stride)[2], double divh) const
		{
			return std::pow(divh, h_exponent) * (StencilCoeff<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>>{}(v, stride) + ...);
		}

		size_t print(char* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>>{}.print(out + n)),
				(n += sprintf(out + n, "\n"))), ...);
			return n;
		}

		size_t print(FILE* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>>{}.print(out)),
				(n += fprintf(out, "\n"))), ...);
			return n;
		}
	};

	template<int... Is, int... Js, int... Ks, typename... Es>
	struct GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...>>
	{
		static const size_t h_exponent = fixed_max<expr::factor_count<expr::symbols::h_symbol, typename coeff_type_divh<Es>::type>::value...>;

		GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...>) {}

		template<typename T>
		auto operator()(T const* v, len_type const (&stride)[3], double divh) const
		{
			return std::pow(divh, h_exponent) * (StencilCoeff<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>>{}(v, stride) + ...);
		}

		size_t print(char* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>>{}.print(out + n)),
				(n += sprintf(out + n, "\n"))), ...);
			return n;
		}

		size_t print(FILE* out) const
		{
			size_t n = 0;
			(((n += StencilCoeff<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>>{}.print(out)),
				(n += fprintf(out, "\n"))), ...);
			return n;
		}
	};

	template<int... Is, typename... Es>
	GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...>) ->
		GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S1_symbol<Is>, Es>...>>;
	template<int... Is, int... Js, typename... Es>
	GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>) ->
		GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S2_symbol<Is, Js>, Es>...>>;
	template<int... Is, int... Js, int... Ks, typename... Es>
	GeneratedStencilApply(types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...>) ->
		GeneratedStencilApply<types_list<std::pair<expr::symbols::internal::S3_symbol<Is, Js, Ks>, Es>...>>;
}

namespace expr
{
	template<typename Stt>
	size_t print_stencil(GeneratedStencilApply<Stt> const& stencil)
	{
		size_t n = 0;
		n += stencil.print(stdout);
		n += printf("\n");
		return n;
	}

	template<typename Stt>
	size_t print_stencil(Stt const& stencil)
	{
		size_t n = 0;
		n += symphas::internal::print_stencil(stencil);
		n += printf("\n");
		return n;
	}
}

/*
New 21-point stencil:

0S{-2}{-2}
((1/72)) / (h^2)S{-1}{-2}
(-(1/9)) / (h^2)S{0}{-2}
((1/72)) / (h^2)S{1}{-2}
0S{2}{-2}
((1/72)) / (h^2)S{-2}{-1}
(-(1/9)) / (h^2)S{-1}{-1}
((55/36)) / (h^2)S{0}{-1}
(-(1/9)) / (h^2)S{1}{-1}
((1/72)) / (h^2)S{2}{-1}
(-(1/9)) / (h^2)S{-2}{0}
((55/36)) / (h^2)S{-1}{0}
(-(16/3)) / (h^2)S{0}{0}
((55/36)) / (h^2)S{1}{0}
(-(1/9)) / (h^2)S{2}{0}
((1/72)) / (h^2)S{-2}{1}
(-(1/9)) / (h^2)S{-1}{1}
((55/36)) / (h^2)S{0}{1}
(-(1/9)) / (h^2)S{1}{1}
((1/72)) / (h^2)S{2}{1}
0S{-2}{2}
((1/72)) / (h^2)S{-1}{2}
(-(1/9)) / (h^2)S{0}{2}
((1/72)) / (h^2)S{1}{2}
0S{2}{2}




*/
