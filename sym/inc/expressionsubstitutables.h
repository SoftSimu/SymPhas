
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
 * PURPOSE: Implements substitutable functions and series, including sum
 * and product.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressiontransforms.h"


#ifdef LATEX_PLOT

#define SYEX_SUM_FMT_AA "\\sum_{"
#define SYEX_SUM_FMT_AB "}"
#define SYEX_SUM_FMT_BA "{"
#define SYEX_SUM_FMT_BB "}"

#else


#define SYEX_SUM_FMT_AA "SUM["
#define SYEX_SUM_FMT_AB "]"
#define SYEX_SUM_FMT_BA "("
#define SYEX_SUM_FMT_BB ")"

#endif


#define SYEX_SUM_FMT SYEX_SUM_FMT_AA "%s" SYEX_SUM_FMT_AB SYEX_SUM_FMT_BA "%s" SYEX_SUM_FMT_BB
#define SYEX_MUL_FMT_LEN (STR_ARR_LEN(SYEX_SUM_FMT_AA SYEX_SUM_FMT_AB SYEX_SUM_FMT_BA SYEX_SUM_FMT_BB) - 1)


namespace expr
{
	namespace symbols
	{
		//! Substitutable variable for use in functions.
		/*!
		 * Object used as an enumerated argument for functions.
		 *
		 * Functions can be created without specifying an argument list, and
		 * the created function object will count the given expr::symbols::arg types used in the
		 * expression to understand how many arguments can be passed. When expressions are
		 * substituted in the function in order to expand it, their order in the function
		 * substitution list (the paramters to the `operator()` method of the function object)
		 * will be matched in order of increasing index of the argument objects used to construct
		 * the function expression.
		 */
		template<size_t N>
		using arg = OpTerm<OpIdentity, Variable<N>>;

	}
}

DEFINE_SYMBOL_ID((size_t N), (expr::symbols::arg<N>), static char* name = expr::print_with_subscript<N>("arg").new_str(); return name;)
ALLOW_COMBINATION((size_t N), (expr::symbols::arg<N>))

//! A function into which other expressions can be substituted when called.
/*!
 * A function into which other expressions can be substituted when called.
 * The function is created with a given number parameters, which may be zero.
 *
 * \tparam E The type of the expression function.
 * \tparam ArgsNs... The indices of the independent variables, v<N>.
 */
	template<typename E, size_t... ArgNs>
struct OpFuncSubstitutable
{
	OpFuncSubstitutable(E const& e) : e{ e } {}

protected:

	template<size_t... Ss>
	inline auto swap_x()
	{
		return e;
	}

	template<size_t S, size_t... Ss, typename X, typename... Xs>
	auto swap_x(X const& x, Xs const& ...xs)
	{
		return expr::transform::swap_grid<S>(swap_x<Ss...>(xs...), x);
	}

public:

	E operator()() const
	{
		return e;
	}

	template<typename... Xs, std::enable_if_t<(sizeof...(Xs) == sizeof...(ArgNs)), int> = 0>
	auto operator()(OpExpression<Xs> const& ...xs)
	{
		return swap_x<ArgNs...>(*static_cast<const Xs*>(&xs)...);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

	E e;				//!< The substitutable function.
};



namespace expr
{

	namespace
	{

		template<typename G>
		struct tuple_list_to_ptrs
		{
			using type = G;
		};

		template<typename... Gs>
		struct tuple_list_to_ptrs<std::tuple<Gs...>>
		{
			using type = std::tuple<Gs*...>;
		};

		template<size_t N>
		constexpr std::true_type _is_subvar(Variable<N>)
		{
			return {};
		}

		constexpr std::false_type _is_subvar(...)
		{
			return {};
		}

		template<typename T>
		constexpr auto is_subvar(T t)
		{
			return _is_subvar(t);
		}

		template<size_t N>
		constexpr std::index_sequence<N> subvar_id(Variable<N>*)
		{
			return {};
		}

		constexpr auto get_subvar_list(std::tuple<>)
		{
			return std::index_sequence<>{};
		}

		template<typename G0, typename... Gs>
		constexpr auto get_subvar_list(std::tuple<G0*, Gs*...>)
		{
			if constexpr (std::invoke_result_t<decltype(is_subvar<G0>), G0>::value)
			{
				return symphas::lib::fold_unique_ids(symphas::lib::seq_join(subvar_id((G0*)nullptr), get_subvar_list(std::tuple<Gs*...>{})));
			}
			else
			{
				return get_subvar_list(std::tuple<Gs*...>{});
			}
		}
	}

	template<typename E>
	constexpr auto get_independent_variables(OpExpression<E> const& e)
	{
		return get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
	}

	template<typename E, size_t... ArgNs>
	constexpr auto get_independent_variables(OpFuncSubstitutable<E, ArgNs...> const& e)
	{
		return std::index_sequence<ArgNs...>{};
	}

	namespace
	{


		template<typename E, size_t... ArgNs>
		OpFuncSubstitutable<E, ArgNs...> make_new_function(OpExpression<E> const& e, std::index_sequence<ArgNs...>)
		{
			return { *static_cast<const E*>(&e) };
		}


		//! Create a substitutable function.
		template<typename E>
		auto make_substitutable_function(OpExpression<E> const& e)
		{
			constexpr auto all_ids = get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
			return make_new_function(*static_cast<const E*>(&e), all_ids);
		}

		//! Create a substitutable function.
		template<typename E, size_t... Ns>
		auto make_substitutable_function(OpExpression<E> const& e, std::index_sequence<Ns...> filter)
		{
			constexpr auto all_ids = get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
			return make_new_function(*static_cast<const E*>(&e), symphas::lib::filter_seq(all_ids, filter));
		}


		template<size_t... Ns>
		struct OpSymbolicFunction;

		template<size_t N0, size_t... Ns>
		struct OpSymbolicFunction<N0, Ns...>
		{

		protected:

			template<typename S>
			struct convert_vars;

			template<int... Ms>
			struct convert_vars<std::integer_sequence<int, Ms...>>
			{
				using type = std::index_sequence<size_t(Ms)...>;
			};

			template<typename S>
			struct has_enough_vars;

			template<size_t... Ms>
			struct has_enough_vars<std::index_sequence<Ms...>>
			{
				using check_seq = symphas::lib::filter_seq_t<std::index_sequence<N0, Ns...>, std::index_sequence<Ms...>>;

				static const bool value = std::is_same<std::index_sequence<>, check_seq>::value;
			};

		public:

			template<typename E, typename Seq = decltype(expr::get_independent_variables(std::declval<E>())),
				std::enable_if_t<has_enough_vars<Seq>::value, int> = 0>
			auto operator=(OpExpression<E> const& e)
			{
				return make_substitutable_function(*static_cast<const E*>(&e), symphas::lib::filter_seq(Seq{}, std::index_sequence<N0, Ns...>{}));
			}

			OpSymbolicFunction(OpSymbolicFunction<N0, Ns...> const&) = delete;
			OpSymbolicFunction(OpSymbolicFunction<N0, Ns...>&&) = delete;
		};

		template<>
		struct OpSymbolicFunction<>
		{
			template<typename E>
			auto operator=(OpExpression<E> const& e)
			{
				return make_substitutable_function(*static_cast<const E*>(&e));
			}

			OpSymbolicFunction(OpSymbolicFunction<> const&) = delete;
			OpSymbolicFunction(OpSymbolicFunction<>&&) = delete;
		};


		template<typename... Gs>
		struct OpSymbolicFunctionSwap
		{

		public:

			template<size_t N, typename E>
			inline auto swap_all_grids(OpExpression<E> const& e)
			{
				return *static_cast<E const*>(&e);
			}

			template<size_t N, typename E, typename X, typename... Xs, typename Seq = decltype(expr::get_independent_variables(std::declval<E>()))>
			auto swap_all_grids(OpExpression<E> const& e)
			{
				constexpr size_t M = N - sizeof...(Xs) - 1;
				if constexpr (!symphas::lib::is_value_in_seq<size_t, M, Seq>::value)
				{
					return expr::transform::swap_grid<X>(
						swap_all_grids<N, E, Xs...>(*static_cast<E const*>(&e)), expr::symbols::arg<M>{ OpVoid{} });
				}
				else
				{
					return swap_all_grids<N + 1, E, X, Xs...>(*static_cast<const E*>(&e));
				}
			}

		public:

			template<typename E, typename Seq = decltype(expr::get_independent_variables(std::declval<E>()))>
			auto operator=(OpExpression<E> const& e)
			{
				return make_substitutable_function(
					swap_all_grids<sizeof...(Gs), E, Gs...>(*static_cast<const E*>(&e)), Seq{});
			}

			OpSymbolicFunctionSwap(OpSymbolicFunctionSwap<Gs...> const&) = delete;
			OpSymbolicFunctionSwap(OpSymbolicFunctionSwap<Gs...>&&) = delete;
		};

	}

	template<size_t... Ns>
	OpSymbolicFunction<Ns...> fn(symbols::arg<Ns>...) { return {}; }

	//! Create a function such that parameters are interpreted from the expression.
	/*!
	 * Functions can be created without specifying an argument list, and
	 * the created function object will count the given expr::symbols::arg types used in the
	 * expression to understand how many arguments can be passed. When expressions are
	 * substituted in the function in order to expand it, their order in the function
	 * substitution list (the paramters to the `operator()` method of the function object)
	 * will be matched in order of increasing index of the argument objects used to construct
	 * the function expression.
	 */
	OpSymbolicFunction<> fn() { return {}; }

	//! Create a function such that the parameters are ::Variable types of the prescribed indices.
	/*!
	 * The linear variables (with OpIdentity coefficients) are given as arguments of the function,
	 * and they must be defined in terms of a ::Variable.
	 */
	template<size_t... Ns, typename... Gs>
	OpSymbolicFunctionSwap<Variable<Ns, Gs>...> fn(OpTerm<OpIdentity, Variable<Ns, Gs>> const&...) { return {}; }

	//! Create a function such that the parameters are ::Variable types of the prescribed indices.
	/*!
	 * Linear variables of the associated data types which have their data as one of the
	 * types defined in the function argument list will be substituted. The grid types must
	 * satisfy expr::is_combinable; it is not possible to use any data type as a substitution.
	 */
	template<typename... Gs, typename std::enable_if_t<(expr::is_combinable<Gs> && ...), int> = 0>
	OpSymbolicFunctionSwap<Gs...> fn(OpTerm<OpIdentity, Gs> const&...) { return {}; }

	//! Create a function such that the parameters are expr::symbols::Symbol types.
	/*!
	 * Linear variables of the associated function expression which have their data as one of the
	 * symbols defined in the function argument list will be substituted.
	 */
	template<size_t... Ns, typename... Symbols, typename std::enable_if_t<(std::is_convertible<Symbols, expr::symbols::Symbol>::value && ...), int> = 0>
	OpSymbolicFunctionSwap<Symbols...> fn(Symbols const&...) { return {}; }

}

namespace expr
{

	namespace symbols
	{
		template<int N, int P = 0>
		struct i_;

		template<int N, int P = 0>
		using i_op_type = OpTerm<OpIdentity, i_<N, P>>;


		template<int N, int P>
		struct i_ : Symbol
		{
			constexpr operator i_op_type<N, P>()
			{
				return {};
			}

			template<int M, int Q>
			auto operator+(i_<M, Q>)
			{
				return i_op_type<N, P>() + i_op_type<M, Q>();
			}

			template<int M, int Q>
			auto operator-(i_<M, Q>)
			{
				return i_op_type<N, P>() - i_op_type<M, Q>();
			}

			template<int M, int Q>
			auto operator*(i_<M, Q>)
			{
				return i_op_type<N, P>() * i_op_type<M, Q>();
			}

			template<int M, int Q>
			auto operator/(i_<M, Q>)
			{
				return i_op_type<N, P>() / i_op_type<M, Q>();
			}

			template<int M, int Q>
			auto operator=(i_<M, Q>)
			{
				return symphas::lib::types_list<i_<N, P>, i_<M, Q>>{};
			}

			template<int M, int Q>
			auto operator!=(i_<M, Q>)
			{
				return symphas::lib::types_list<std::false_type, i_<N, P>, i_<M, Q>>{};
			}
		};


		template<typename, typename...>
		struct v_id_type;
		template<int N0, int P0, int... Ns, int... Ps>
		struct v_id_type<i_<N0, P0>, i_<Ns, Ps>...> : Symbol {};

		template<typename I0, typename... Is>
		using v_ = OpTerm<OpIdentity, v_id_type<I0, Is...>>;


		template<int N0, int P0, int... Ns, int... Ps>
		auto make_sum_variable(i_<N0, P0>, i_<Ns, Ps>...)
		{
			return v_<i_<N0, P0>, i_<Ns, Ps>...>{};
		};

		template<int P0, int N, int P>
		auto make_sum_variable(i_<N, P>)
		{
			return v_<i_<N, P + P0>>{};
		};
	}


	template<typename I, typename E,
		typename = std::enable_if_t<std::is_convertible<I, expr::symbols::Symbol>::value, int>>
		auto sum(OpExpression<E> const& e);

	template<typename I0, typename I1, typename... Is, typename E,
		typename = std::enable_if_t<(
			std::is_convertible<I0, expr::symbols::Symbol>::value&&
			std::is_convertible<I1, expr::symbols::Symbol>::value &&
			(std::is_convertible<Is, expr::symbols::Symbol>::value && ...)), int>>
		auto sum(OpExpression<E> const& e);

	template<typename I0, typename... Is, int N>
	auto sum(expr::symbols::i_<N> const& e)
	{
		return sum<I0, Is...>(expr::symbols::i_op_type<N>{});
	}

	template<typename I, typename E,
		typename = std::enable_if_t<std::is_convertible<I, expr::symbols::Symbol>::value, int>>
		auto prod(OpExpression<E> const& e);

	template<typename I0, typename I1, typename... Is, typename E,
		typename = std::enable_if_t<(
			std::is_convertible<I0, expr::symbols::Symbol>::value&&
			std::is_convertible<I1, expr::symbols::Symbol>::value &&
			(std::is_convertible<Is, expr::symbols::Symbol>::value && ...)), int>>
		auto prod(OpExpression<E> const& e);

	template<typename I0, typename... Is, int N>
	auto prod(expr::symbols::i_<N> const& e)
	{
		return prod<I0, Is...>(expr::symbols::i_op_type<N>{});
	}

}

DEFINE_SYMBOL_ID((int N, int P), (expr::symbols::i_<N, P>), static char* name = expr::print_with_subscript<N>("i").new_str(); return name;)
ALLOW_COMBINATION((int N, int P), (expr::symbols::i_<N, P>))

DEFINE_SYMBOL_ID((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>), static char* name = expr::print_with_subscripts<Ns...>("v").new_str(); return name;)
ALLOW_COMBINATION((int... Ns, int... Ps), (expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>))

//! Specialization based on expr::base_data_type.
template<int N0, int P0>
struct expr::base_data_type<expr::symbols::i_<N0, P0>>
{
	using type = expr::symbols::i_<N0>;
};

template<int N, int P, typename E>
auto operator+(expr::symbols::i_<N, P>, OpExpression<E> const& e)
{
	return expr::symbols::i_op_type<N, P>{} + (*static_cast<E const*>(&e));
}

template<int N, int P, typename E>
auto operator+(OpExpression<E> const& e, expr::symbols::i_<N, P>)
{
	return (*static_cast<E const*>(&e)) + expr::symbols::i_op_type<N, P>{};
}

template<int N, int P, typename E>
auto operator-(expr::symbols::i_<N, P>, OpExpression<E> const& e)
{
	return expr::symbols::i_op_type<N, P>{} - (*static_cast<E const*>(&e));
}

template<int N, int P, typename E>
auto operator-(OpExpression<E> const& e, expr::symbols::i_<N, P>)
{
	return (*static_cast<E const*>(&e)) - expr::symbols::i_op_type<N, P>{};
}

template<int N, int P, typename E>
auto operator*(expr::symbols::i_<N, P>, OpExpression<E> const& e)
{
	return expr::symbols::i_op_type<N, P>{} * (*static_cast<E const*>(&e));
}

template<int N, int P, typename E>
auto operator*(OpExpression<E> const& e, expr::symbols::i_<N, P>)
{
	return (*static_cast<E const*>(&e)) * expr::symbols::i_op_type<N, P>{};
}

template<int N, int P, typename E>
auto operator/(expr::symbols::i_<N, P>, OpExpression<E> const& e)
{
	return expr::symbols::i_op_type<N, P>{} / (*static_cast<E const*>(&e));
}

template<int N, int P, typename E>
auto operator/(OpExpression<E> const& e, expr::symbols::i_<N, P>)
{
	return (*static_cast<E const*>(&e)) / expr::symbols::i_op_type<N, P>{};
}

namespace symphas::internal
{
	enum struct CompoundOp
	{
		ADD, SUB, MUL, DIV
	};

	template<CompoundOp>
	struct identity_op;

	template<>
	struct identity_op<CompoundOp::ADD>
	{
		using type = OpVoid;
	};

	template<>
	struct identity_op<CompoundOp::SUB>
	{
		using type = OpVoid;
	};

	template<>
	struct identity_op<CompoundOp::MUL>
	{
		using type = OpIdentity;
	};

	template<>
	struct identity_op<CompoundOp::DIV>
	{
		using type = OpIdentity;
	};


	template<CompoundOp>
	struct recursive_op_apply;


	template<>
	struct recursive_op_apply<CompoundOp::ADD>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (std::forward<Es>(es) + ... + OpVoid{});
		}
	};

	template<>
	struct recursive_op_apply<CompoundOp::SUB>
	{
		template<typename E0, typename... Es>
		auto operator()(E0&& e0, Es&&... es)
		{
			return std::forward<E0>(e0) - (std::forward<Es>(es) + ... + OpVoid{});
		}
	};

	template<>
	struct recursive_op_apply<CompoundOp::MUL>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (OpIdentity{} * ... * std::forward<Es>(es));
		}
	};

	template<>
	struct recursive_op_apply<CompoundOp::DIV>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (std::forward<Es>(es) / ... / OpIdentity{});
		}
	};



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

	template<typename... Ts>
	struct select_v_;

	template<int N0, int P0, typename... Vs>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int N0, int P0, typename... Vs, int... Ns, int... Ps, typename... Rest>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>, Rest...>>
	{
		static const bool flag = symphas::lib::is_value_in_seq<int, N0, std::integer_sequence<int, Ns...>>::value;

		using type = typename select_v_<
			expr::symbols::i_<N0, P0>,
			std::conditional_t<flag, 
				symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>...>>,
				symphas::lib::types_list<Vs...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int N0, int P0, typename... Vs, typename T, typename... Rest>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_<
			expr::symbols::i_<N0, P0>,
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int N0, int P0, typename... Vs>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>>
	{
		using type = typename select_v_<
			expr::symbols::i_<N0, P0>,
			symphas::lib::types_list<>,
			symphas::lib::types_list<Vs...>>::type;
	};


	template<typename... Vs, typename... Is, typename... Rest>
	struct select_v_<symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<Is...>, Rest...>>
	{
		using type = typename select_v_<
			symphas::lib::types_list<Vs..., expr::symbols::v_id_type<Is...>>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Vs, typename T, typename... Rest>
	struct select_v_<symphas::lib::types_list<Vs...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_v_<
			symphas::lib::types_list<Vs...>,
			symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... Vs>
	struct select_v_<symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<typename... Vs>
	struct select_v_<symphas::lib::types_list<Vs...>>
	{
		using type = typename select_v_<
			symphas::lib::types_list<>,
			symphas::lib::types_list<Vs...>>::type;
	};

	
	template<typename... Ts>
	struct select_all_i_;

	template<int N, typename... Is>
	struct select_all_i_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Is...>;
	};

	template<int N, typename... Is, int N0, int P0, typename... Rest>
	struct select_all_i_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N0, P0>>, Rest...>>
	{
		using type = typename select_all_i_<expr::symbols::i_<N, 0>, 
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, int P0, typename... Rest>
	struct select_all_i_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N, P0>>, Rest...>>
	{
		using type = typename select_all_i_<expr::symbols::i_<N, 0>, 
			symphas::lib::types_list<Is..., expr::symbols::i_<N, P0>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, typename... Is, typename Other, typename... Rest>
	struct select_all_i_<expr::symbols::i_<N, 0>, symphas::lib::types_list<Is...>, symphas::lib::types_list<Other, Rest...>>
	{
		using type = typename select_all_i_<expr::symbols::i_<N, 0>,
			symphas::lib::types_list<Is...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename... Vs>
	struct select_all_i_<expr::symbols::i_<N, P>, symphas::lib::types_list<Vs...>>
	{
		using type = typename select_all_i_<expr::symbols::i_<N, 0>, 
			symphas::lib::types_list<>, symphas::lib::types_list<Vs...>>::type;
	};
}

//! Defines a summation of an expression with indices and substitutable expressions.
/*!
 * A summation is performed over an index, which is provided by the user, and is of type
 * expr::symbols::i_, with an `int` template parameter indicating the "value" of the index. This
 * is to make it possible to perform multiple summations or nested summations.
 *
 * A sum is expanded by calling its member function OpSum::expand.
 * Terms or other expressions can also be substituted into sums when they are expanded.
 *
 * A sum is written like a regular expression, with the additional use of expr::symbols::i_ to
 * substitute the current index that is being summed, as well of the expr::symbols::v_, which takes
 * the sum's corresponding index, expr::symbols::i_, in order to allow terms or expressions
 * to be substituted into the expansion of the sum.
 */
template<symphas::internal::CompoundOp Op, typename I, typename E, typename... Is>
struct OpCompound;

template<symphas::internal::CompoundOp Op, typename I, typename E>
struct OpCompound<Op, I, E> : OpExpression<OpCompound<Op, I, E>>
{
protected:

	using all_var_types = typename expr::op_types<E>::type;
	using v_types = typename symphas::internal::select_v_<I, all_var_types>::type;
	using all_indices_of_id = typename symphas::internal::select_all_i_<I, all_var_types>::type;

	template<typename V>
	static constexpr int index_of_v = symphas::lib::index_of_type<V, v_types>;

	template<typename V>
	static constexpr int index_of_type = symphas::lib::index_of_type<V, all_var_types>;


public:

	OpCompound(E const& e) : e{ e } {}

	auto eval(iter_type n = 0)
	{
		return e.eval(n);
	}

	//! Expand the sum between the given indices.
	/*!
	 * Expand the sum between the given indices, where `I0` is the first index, and `I1` is the
	 * final index. The final index is included in the sum (**not** C-style, where the last index of
	 * a loop is omitted).
	 *
	 * This expand function can only be used if there are no substitutable terms
	 * corresponding to this index.
	 */
	template<int I0, int I1, typename V = expr::symbols::v_id_type<I>, 
		typename std::enable_if_t<(index_of_v<V> < 0), int> = 0>
	auto expand() const
	{
		if constexpr (I1 < I0)
		{
			return typename symphas::internal::identity_op<Op>::type{};
		}
		else
		{
			return expand_integer(
				all_indices_of_id{},
				symphas::lib::seq_add(
					std::make_integer_sequence<int, I1 - I0>{},
					symphas::lib::seq_repeating_value_t<I1 - I0, int, I0>{}));
		}
	}


	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int I0, typename E0, typename... Es, typename V = expr::symbols::v_id_type<I>,
		typename std::enable_if_t<(index_of_v<V> >= 0), int> = 0>
	auto expand(OpExpression<E0> const& sub, OpExpression<Es> const& ...subs) const
	{
		return expand_function<I0>(
			std::make_index_sequence<sizeof...(Es) + 1>{}, all_indices_of_id{},
			std::make_tuple(*static_cast<E0 const*>(&sub), *static_cast<Es const*>(&subs)...));
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int I0, typename E0, typename... Es, typename V = expr::symbols::v_id_type<I>,
		typename std::enable_if_t<(index_of_v<V> >= 0), int> = 0>
	auto expand(std::tuple<E0, Es...> const& subs) const
	{
		return expand_function<I0>(std::make_index_sequence<sizeof...(Es) + 1>{}, all_indices_of_id{}, subs);
	}

	auto operator-() const
	{
		return expr::sum<I>(-e);
	}

	//// TODO
#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

	E e;				//!< The expression which is summed.

protected:

	template<int I0, int I00, int... IJ0s>
	auto swap_ii(symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>)
	{
		return expr::transform::swap_grid<expr::symbols::i_<I00, IJ0s>...>(e, expr::make_integer<I0 + IJ0s>()...);
	}

	template<int I00, int... IJ0s, int... Is>
	auto expand_integer(symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>, std::integer_sequence<int, Is...>) const
	{
		if constexpr (((index_of_type<expr::symbols::i_<I00, IJ0s>> >= 0) || ...))
		{
			return symphas::internal::recursive_op_apply<Op>{}(swap_ii<Is>(symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>{})...);
		}
		else
		{
			auto ee = [&] (int) { return e; };
			return symphas::internal::recursive_op_apply<Op>{}(ee(Is)...);
		}
	}

	template<int N, typename... Es>
	const auto& get_substituted(std::tuple<Es...> const& subs) const
	{
		if constexpr (N < 0)
		{
			return std::get<sizeof...(Es) + N>(subs);
		}
		else if constexpr (N >= sizeof...(Es))
		{
			return std::get<N - sizeof...(Es)>(subs);
		}
		else
		{
			return std::get<N>(subs);
		}
	}

	template<size_t Is, int I0, int P0, typename Ex, typename... Es>
	auto substitute_functions(
		expr::symbols::i_<I0, P0>,
		symphas::lib::types_list<>,
		Ex const& ex,
		std::tuple<Es...> const& subs) const
	{
		return ex;
	}

	template<size_t Is, int I0, int P0, typename Ex, typename... Es, int... Ps>
	auto substitute_functions(
		expr::symbols::i_<I0, P0>,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, Ps>>...>, 
		Ex const& ex, 
		std::tuple<Es...> const& subs) const
	{
		return expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<I0, Ps>>...>(ex, get_substituted<int(Is) + P0 + Ps>(subs)...);
	}

	template<int I0, int I00, int... IJ0s, typename... Es, size_t... Is>
	auto expand_function(std::index_sequence<Is...>, 
		symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>, 
		std::tuple<Es...> const& subs) const
	{
		if constexpr (((index_of_type<expr::symbols::i_<I00, IJ0s>> >= 0) || ...))
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				substitute_functions<Is>(
					I{},
					v_types{},
					swap_ii<Is>(symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>{}),
					subs)...);
		}
		else
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				substitute_functions<Is>(I{}, v_types{}, e, subs)...);
		}
	}

	template<int I0>
	auto expand() const
	{
		return typename symphas::internal::identity_op<Op>::type{};
	}

};

//! Defines a summation of an expression with indices and substitutable expressions.
/*!
 * A summation is performed over an index, which is provided by the user, and is of type
 * expr::symbols::i_, with an `int` template parameter indicating the "value" of the index. This
 * is to make it possible to perform multiple summations or nested summations.
 *
 * A sum is expanded by calling its member function OpSum::expand.
 * Terms or other expressions can also be substituted into sums when they are expanded.
 *
 * A sum is written like a regular expression, with the additional use of expr::symbols::i_ to
 * substitute the current index that is being summed, as well of the expr::symbols::v_, which takes
 * the sum's corresponding index, expr::symbols::i_, in order to allow terms or expressions
 * to be substituted into the expansion of the sum.
 */
template<symphas::internal::CompoundOp Op, typename I0, typename E, typename I1, typename... Is>
struct OpCompound<Op, I0, E, I1, Is...> : OpExpression<OpCompound<Op, I0, E, I1, Is...>>
{

protected:

	using all_var_types = typename expr::op_types<E>::type;
	using v_types = typename symphas::internal::select_v_<I0, all_var_types>::type;
	using all_v_types = typename symphas::internal::select_v_<all_var_types>::type;

	template<typename II>
	using all_indices_of_id = typename symphas::internal::select_all_i_<II, all_var_types>::type;

	template<typename V>
	static constexpr int index_of_type = symphas::lib::index_of_type<V, all_var_types>;



public:

	OpCompound(E const& e) : e{ e } {}

	auto eval(iter_type n = 0)
	{
		return e.eval(n);
	}

	//! Expand the sum between the given indices.
	/*!
	 * Expand the sum between the given indices, where `I0` is the first index, and `I1` is the
	 * final index. The final index is included in the sum (**not** C-style, where the last index of
	 * a loop is omitted).
	 *
	 * This expand function can only be used if there are no substitutable terms
	 * corresponding to this index.
	 */
	template<int I00, int I01, typename V0 = expr::symbols::v_id_type<I0>, typename Vs = expr::symbols::v_id_type<I0, I1, Is...>,
		typename std::enable_if_t<(index_of_type<V0> < 0 && index_of_type<Vs> < 0), int> = 0>
	auto expand() const
	{
		auto expr = OpCompound<Op, I0, E>(e).template expand<I00, I01>();
		return OpCompound<Op, I1, decltype(expr), Is...>(expr);
	}


	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int I00, typename... S, typename Vs = expr::symbols::v_id_type<I0, I1, Is...>,
		typename std::enable_if_t<(index_of_type<Vs> < 0), int> = 0>
	auto expand(S&&... s) const
	{
		auto expr = OpCompound<Op, I0, E>(e).template expand<I00>(std::forward<S>(s)...);
		return OpCompound<Op, I1, decltype(expr), Is...>(expr);
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int... I0s, typename T0, typename... Ts, typename Vs = expr::symbols::v_id_type<I0, I1, Is...>,
		typename std::enable_if_t<(index_of_type<Vs> >= 0 && sizeof...(I0s) == sizeof...(Is) + 2), int> = 0>
	auto expand(std::tuple<T0, Ts...> const& sub) const
	{
		return expand_nested<0>(sub, e, std::integer_sequence<int, I0s...>{}, std::make_index_sequence<sizeof...(Ts) + 1>{});
	}

	//! Start the iterated sum expansion using a single list of functions.
	/*!
	 * Expand the sum, starting at the given index. The rest of the indices have offsets that
	 * indicate the starting position of the substitutions in the list. The offsets are not based
	 * on the primary index, can only be positive, and do not loop around the list.
	 * 
	 * The list of functions is non-nested; it is assumed that the iterated sum has expressions 
	 * that are indexed by only a single value, e.g. $\psi_i$.
	 * 
	 * \param subs The list of expressions that is substituted, only a single list
	 * since it is assumed that the iterated sum has expressions that are indexed by only
	 * a single value, e.g. $\psi_i$.
	 * 
	 * \tparam I The starting value of the first iteration index.
	 * \tparam ISs... The starting value of the other indices is set to the current running index
	 * of the specified index.
	 */
	template<int I, typename... ISs, typename T0, typename... Ts, typename Vs = expr::symbols::v_id_type<I0>,
		typename std::enable_if_t<(index_of_type<Vs> >= 0 && sizeof...(ISs) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<T0, Ts...> const& subs) const
	{
		return expand_iterated<I>(subs, symphas::lib::types_list<ISs...>{}, std::make_index_sequence<sizeof...(Ts) + 1>{});
	}

	//! Start the iterated sum expansion using a single list of functions.
	template<int I, typename... ISs, typename T0, typename... Ts, typename Vs = expr::symbols::v_id_type<I0>,
		typename std::enable_if_t<(index_of_type<Vs> >= 0 && sizeof...(ISs) == sizeof...(Is) + 1), int> = 0>
	auto expand(T0 const& sub0, Ts const&... subs) const
	{
		return expand<I, ISs...>(std::make_tuple(sub0, subs...));
	}

	auto operator-() const
	{
		return expr::sum<I0, I1, Is...>(-e);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

protected:

	template<int Ip, typename... Ts>
	auto substitute_from_list(std::tuple<Ts...> const& subs) const
	{
		if constexpr (Ip < 0 || Ip >= sizeof...(Ts))
		{
			return OpVoid{};
		}
		else
		{
			return std::get<size_t(Ip)>(subs);
		}
	}

	template<size_t N, typename... Ts, typename EE, int I00, int... I00s>
	auto expand_nested(std::tuple<> const& sub, EE const& expr0, std::integer_sequence<int, I00, I00s...>) const
	{
		return typename symphas::internal::identity_op<Op>::type{};
	}
	
	template<size_t N, typename E0, typename... Es, typename EE, int I00, size_t... Js, typename Vs = expr::symbols::v_<I0, I1, Is...>,
		typename In = std::tuple_element_t<N, std::tuple<I0, I1, Is...>>, typename V0 = expr::symbols::v_<In>>
	auto expand_nested(std::tuple<E0, Es...> const& sub, EE const& expr0, 
		std::integer_sequence<int, I00>, std::index_sequence<Js...>) const
	{
		if constexpr (index_of_type<In> >= 0)
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				expr::transform::swap_grid<Vs>(
					expr::transform::swap_grid<In>(expr0, expr::make_integer<int(Js) + I00>()),
					std::get<Js>(sub))...);
		}
		else
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				expr::transform::swap_grid<Vs>(expr0, std::get<Js>(sub))...);
		}
	}

	template<size_t N, typename... T0s, typename... Ts, typename EE, int I00, int... I00s, size_t... Js,
		typename In = std::tuple_element_t<N, std::tuple<I0, I1, Is...>>>
	auto expand_nested(std::tuple<std::tuple<T0s...>, Ts...> const& sub, EE const& expr0,
		std::integer_sequence<int, I00, I00s...>, std::index_sequence<Js...>) const
	{
		using sub_t = std::tuple<std::tuple<T0s...>, Ts...>;

		if constexpr (index_of_type<In> >= 0)
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				expand_nested<N + 1>(
					std::get<Js>(sub),
					expr::transform::swap_grid<In>(expr0, expr::make_integer<int(Js) + I00>()),
					std::integer_sequence<int, I00s...>{},
					std::make_index_sequence<std::tuple_size<std::tuple_element_t<Js, sub_t>>::value>{})...
				);
		}
		else
		{
			return symphas::internal::recursive_op_apply<Op>{}(
				expand_nested<N + 1>(
					std::get<Js>(sub),
					expr0,
					std::integer_sequence<int, I00s...>{},
					std::make_index_sequence<std::tuple_size<std::tuple_element_t<Js, sub_t>>::value>{})...
				);
		}

	}

	template<int JI, int J0, int JP, typename T0, typename... Ts, typename EE, int... Ns, int... Ps>
	auto construct_iterated(
		std::tuple<T0, Ts...> const& subs, 
		EE const& expr0, 
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>) const
	{
		return OpCompound<Op, expr::symbols::i_<J0, JP>, EE, Is...>(expr0).template expand<JI, expr::symbols::i_<Ns, Ps>...>(subs);
	}

	template<int JI, int J0, int JP, typename... Ts, typename EE, int... Ns, int... Ps>
	auto construct_iterated(
		std::tuple<> const& subs,
		EE const& expr0,
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>) const
	{
		return typename symphas::internal::identity_op<Op>::type{};
	}

	//! Processes the first index for an iterated sum.
	/*!
	 * If any index is set to the current index that is being expanded, it will be set equal to
	 * its own ID value in order to be substituted correctly. Its global offset will be
	 * set to the current value of this index, plus its own offset.
	 * 
	 * \param subs The list of functions that is substituted.
	 * 
	 * \tparam I is the value that the first index starts at.
	 * \tparam Ic is the current position in the list for the first index.
	 * \tparam N0 The ID of the index which `I1` will start at.
	 */
	template<int I, size_t Ic, int I00, int... IJ0s, int... IJ1s, int N0, int P0, int... Ns, int... Ps,
		int IP0, int J00, int JP0, int... Ls, int... Rs, typename... Ts>
	auto construct_iterated(
		std::tuple<Ts...> const& subs,
		symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>,							// all offsets of I0 in the expression
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I00, IJ1s>>...>,// all index functions to substitute, v_<i_...>>
		symphas::lib::types_list<expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>,	// All iterated index starting values
		symphas::lib::types_list<
			expr::symbols::i_<I00, IP0>,								// I0
			expr::symbols::i_<J00, JP0>,								// I1
			expr::symbols::i_<Ls, Rs>...>) const						// All indices to be substituted, Is...
	{
		constexpr int Ip = int(Ic) + IP0;								// position from which to start picking from substitution list
		constexpr int JI = (N0 == I00) ? 0 : P0;						// starting value of the next index
		constexpr int J0 = (N0 == I00) ? J00 : N0;						// the identity of the next index
		constexpr int JP = (N0 == I00) ? Ip + P0 : P0;

		using J0s = symphas::lib::types_list<
			expr::symbols::i_<
				(Ns == I00) ? Ls : Ns,		// if the index starts at I0, set it back to its own ID
				(Ns == I00) ? Ip + Ps : Ps>	// if the index starts at I0, add the offset of I0 to it
			...>;

		if constexpr (((index_of_type<expr::symbols::v_id_type<expr::symbols::i_<I00, IJ0s>>> >= 0) || ...))
		{
			if constexpr (((index_of_type<expr::symbols::i_<I00, IJ0s>> >= 0) || ...))
			{
				return construct_iterated<JI, J0, JP>(
					subs, 
					expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<I00, IJ1s>>...>(
						expr::transform::swap_grid<expr::symbols::i_<I00, IJ0s>...>(e, expr::make_integer<I + int(Ic) + IJ0s>()...),
						substitute_from_list<Ip + IJ1s>(subs)...),
					J0s{});
			}

			return construct_iterated<JI, J0, JP>(
				subs, 
				expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<I00, IJ1s>>...>(e, substitute_from_list<Ip + IJ1s>(subs)...),
				J0s{});
		}
		else 
		{
			return construct_iterated<JI, J0, JP>(subs, e, J0s{});
		}
	}

	//! Start the iterated sum expansion.
	/*!
	 * Start the iterated sum expansion.
	 * 
	 * \param subs The list of expressions that is substituted, only a single list
	 * since it is assumed that the iterated sum has expressions that are indexed by only
	 * a single value, e.g. $\psi_i$. 
	 * 
	 * \tparam Ns... The starting index ID of an index from Is...
	 * \tparam Ps... The starting index value of an index from Is...
	 */
	template<int I, typename... Ts, int... Ns, int... Ps, size_t... Ics>
	auto expand_iterated(std::tuple<Ts...> const& subs, 
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>, 
		std::index_sequence<Ics...>) const
	{
		return symphas::internal::recursive_op_apply<Op>{}(
			construct_iterated<I, Ics>(subs, all_indices_of_id<I0>{}, v_types{},
				symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>{}, // All iterated index starting values
				symphas::lib::types_list<I0, I1, Is...>{})...); // All indices to be substituted, Is
	}

public:

	E e;				//!< The expression which is summed.

};


template<typename coeff_t, typename I, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t, OpCompound<symphas::internal::CompoundOp::ADD, I, E> const& sum)
{
	return expr::sum<I>(coeff_t{} *sum.e);
}

template<typename T, typename I, typename E>
auto operator*(OpLiteral<T> const& a, OpCompound<symphas::internal::CompoundOp::ADD, I, E> const& sum)
{
	return expr::sum<I>(a * sum.e);
}

template<typename coeff_t, typename I0, typename I1, typename... Is, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t, OpCompound<symphas::internal::CompoundOp::ADD, I0, E, I1, Is...> const& sum)
{
	return expr::sum<I0, I1, Is...>(coeff_t{} *sum.e);
}

template<typename T, typename I0, typename I1, typename... Is, typename E>
auto operator*(OpLiteral<T> const& a, OpCompound<symphas::internal::CompoundOp::ADD, I0, E, I1, Is...> const& sum)
{
	return expr::sum<I0, I1, Is...>(a * sum.e);
}


namespace symphas::internal
{


	template<typename... Is>
	struct index_conditions;

	template<>
	struct index_conditions<> {};

	template<typename... Is>
	struct index_conditions
	{
		index_conditions(Is...) {}
		index_conditions() {}

		template<typename E>
		struct sum : index_conditions<E, Is...>, OpExpression<sum<E>>, index_conditions<>
		{
			auto eval(iter_type n) const
			{
				return e.eval(n);
			}

			auto operator-() const
			{
				return index_conditions<Is...>(Is{}...)(-e);
			}

			sum(E const& e) : e{ e } {}
			E e;
		};

		template<typename E>
		auto operator()(OpExpression<E> const& e)
		{
			return construct_sum(*static_cast<E const*>(&e), parse_condition(Is{})...);
		}

	protected:

		template<typename E, typename... I0s>
		auto construct_sum(OpExpression<E> const& e, I0s...)
		{
			return typename index_conditions<I0s...>::template sum<E>(*static_cast<E const*>(&e));
		}

		template<typename C>
		C parse_condition(C)
		{
			return {};
		}

		template<size_t N0, size_t I0, size_t M>
		auto parse_condition(OpAdd<expr::symbols::i_<N0, I0>, OpFractionLiteral<M, 1>>)
		{
			return expr::symbols::i_<N0, I0 + M>{};
		}

		template<size_t N0, size_t I0, size_t M>
		auto parse_condition(OpAdd<expr::symbols::i_<N0, I0>, OpNegFractionLiteral<M, 1>>)
		{
			return expr::symbols::i_<N0, I0 - M>{};
		}

		template<size_t N0, size_t I0>
		auto parse_condition(OpAdd<expr::symbols::i_<N0, I0>, OpIdentity>)
		{
			return expr::symbols::i_<N0, I0 + 1>{};
		}

		template<size_t N0, size_t I0>
		auto parse_condition(OpAdd<expr::symbols::i_<N0, I0>, OpNegIdentity>)
		{
			return expr::symbols::i_<N0, I0 - 1>{};
		}

		template<size_t N0, size_t I0, size_t N1, size_t I1, size_t M>
		auto parse_condition(types_list<expr::symbols::i_<N0, I0>, OpAdd<expr::symbols::i_<N1, I1>, OpFractionLiteral<M, 1>>>)
		{
			return types_list<expr::symbols::i_<N0, I0>, expr::symbols::i_<N1, I1 + M>>{};
		}

		template<size_t N0, size_t I0, size_t N1, size_t I1, size_t M>
		auto parse_condition(types_list<expr::symbols::i_<N0, I0>, OpAdd<expr::symbols::i_<N1, I1>, OpNegFractionLiteral<M, 1>>>)
		{
			return types_list<expr::symbols::i_<N0, I0>, expr::symbols::i_<N1, I1 - M>>{};
		}

		template<size_t N0, size_t I0, size_t N1, size_t I1>
		auto parse_condition(types_list<expr::symbols::i_<N0, I0>, OpAdd<expr::symbols::i_<N1, I1>, OpIdentity>>)
		{
			return types_list<expr::symbols::i_<N0, I0>, expr::symbols::i_<N1, I1 + 1>>{};
		}

		template<size_t N0, size_t I0, size_t N1, size_t I1>
		auto parse_condition(types_list<expr::symbols::i_<N0, I0>, OpAdd<expr::symbols::i_<N1, I1>, OpNegIdentity>>)
		{
			return types_list<expr::symbols::i_<N0, I0>, expr::symbols::i_<N1, I1 - 1>>{};
		}
	};

	template<typename E, typename... Is>
	constexpr bool test_index_conditions(index_conditions<E, Is...>*)
	{
		return true;
	}

	inline constexpr bool test_index_conditions(...)
	{
		return false;
	}


	template<typename E, typename... Is>
	using sum_conditions = typename symphas::internal::index_conditions<Is...>::template sum<E>;


	template<typename... Ts>
	using index_relation = symphas::lib::types_list<Ts...>;
}


namespace expr
{

	//! Defines a summation of an expression with indices and substitutable expressions.
	/*!
	 * A summation is performed over an index, which is provided by the user, and is of type
	 * expr::symbols::i_, with an `int` template parameter indicating the "value" of the index. This
	 * is to make it possible to perform multiple summations or nested summations.
	 *
	 * A sum is expanded by calling its member function OpSum::expand.
	 * Terms or other expressions can also be substituted into sums when they are expanded.
	 *
	 * A sum is written like a regular expression, with the additional use of expr::symbols::i_ to
	 * substitute the current index that is being summed, as well of the expr::symbols::v_, which takes
	 * the sum's corresponding index, expr::symbols::i_, in order to allow terms or expressions
	 * to be substituted into the expansion of the sum.
	 */
	template<typename I, typename E, typename>
	auto sum(OpExpression<E> const& e)
	{
		return OpCompound<symphas::internal::CompoundOp::ADD, I, E>(*static_cast<E const*>(&e));
	}

	//! Defines a summation of an expression with indices and substitutable expressions.
	/*!
	 * A summation is performed over an index, which is provided by the user, and is of type
	 * expr::symbols::i_, with an `int` template parameter indicating the "value" of the index. This
	 * is to make it possible to perform multiple summations or nested summations.
	 *
	 * A sum is expanded by calling its member function OpSum::expand.
	 * Terms or other expressions can also be substituted into sums when they are expanded.
	 *
	 * A sum is written like a regular expression, with the additional use of expr::symbols::i_ to
	 * substitute the current index that is being summed, as well of the expr::symbols::v_, which takes
	 * the sum's corresponding index, expr::symbols::i_, in order to allow terms or expressions
	 * to be substituted into the expansion of the sum.
	 *
	 * Multiple indices can be defined, which indicates multiple sums.
	 */
	template<typename I0, typename I1, typename... Is, typename E, typename>
	auto sum(OpExpression<E> const& e)
	{
		return OpCompound<symphas::internal::CompoundOp::ADD, I0, E, I1, Is...>(*static_cast<E const*>(&e));
	}


	//! Defines a product of an expression with indices and substitutable expressions.
	/*!
	 * See expr::sum.
	 */
	template<typename I, typename E, typename>
	auto prod(OpExpression<E> const& e)
	{
		return OpCompound<symphas::internal::CompoundOp::MUL, I, E>(*static_cast<E const*>(&e));
	}

	//! Defines a product of an expression with indices and substitutable expressions.
	/*!
	 * See expr::sum.
	 */
	template<typename I0, typename I1, typename... Is, typename E, typename>
	auto prod(OpExpression<E> const& e)
	{
		return OpCompound<symphas::internal::CompoundOp::MUL, I0, E, I1, Is...>(*static_cast<E const*>(&e));
	}


	using namespace expr::symbols;



	template<typename E, typename I0, typename... Is, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, I0, Is...> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, I0, Is...> const*>(&sum)).e;
		return expr::sum<I0, Is...>(e).template expand<0, Is...>(std::make_tuple(std::get<Ns>(ops)...));
	}

	template<typename E, int I, int N, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, i_<I, N>> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, i_<I, N>> const*>(&sum)).e;
		return expr::sum<i_<I, N>>(e).template expand<0>(std::get<Ns>(ops)...);
	}

	template<typename E, int I0, int N0, int I1, int N1, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, i_<I0, N0>, i_<I1, N1>> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, i_<I0, N0>, i_<I1, N1>> const*>(&sum)).e;
		return expr::sum<i_<I0, N0>, i_<I1, N1>>(e).template expand<0>(std::get<Ns>(ops)...);
	}

	template<typename E, int I0, int N0, int I1, int N1, int N2, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, i_<I0, N0>, 
		symphas::internal::index_relation<i_<I1, N1>, i_<I0, N2>>> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, i_<I0, N0>, 
			symphas::internal::index_relation<i_<I1, N1>, i_<I0, N2>>> const*>(&sum)).e;

		return expr::sum<i_<I0, N0>, i_<I1, N1>>(e).template expand<0, i_<I1, N2>>(std::make_tuple(std::get<Ns>(ops)...));
	}

	namespace
	{
		template<size_t Ns, size_t NN, int I0, int I1, int N0, int N1, int N2, typename E, typename... Us>
		auto expand_cancelled_term(OpExpression<E> const& e, i_<I0, N0>, i_<I1, N1>, i_<I0, N2>, std::tuple<Us...> const& ops)
		{
			if constexpr (N1 <= Ns + N0 + N2 && Ns + N0 + N2 < NN)
			{
				auto e_j = expr::sum<i_<I1, 0>>(e).template expand<Ns + N0 + N2 - N1>(std::get<Ns + N0 + N2>(ops));
				return expr::sum<i_<I0, 0>>(e_j).template expand<Ns + N0>(std::get<Ns + N0 + N2>(ops));;
			}
			else
			{
				return OpVoid{};
			}
		}
	}

	//! Expand the sum where the condition is: jj != ii + n.
	/*!
	 * The full sum is expanded, where the expansion of jj assumes the iteration value starts at 0.
	 * Then sums are carried out with only jj = ii + n is expanded, they are all summed, and then
	 * ii is substituted as its usual value. This is then subtracted from the full sum, in order
	 * to get an expression where jj != ii. Note that N0 refers to the starting point in the
	 * list of functions that v_ii starts substituting. N2 is means that jj is not supposed to be
	 * ii + N0 + N2.
	 *
	 * \tparam N0 The global offset of index I0
	 * \tparam N1 The global offset of index I1
	 */
	template<typename E, int I0, int N0, int I1, int N1, int N2, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, i_<I0, N0>, 
		symphas::internal::index_relation<std::false_type, i_<I1, N1>, i_<I0, N2>>> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, i_<I0, N0>, 
			symphas::internal::index_relation<std::false_type, i_<I1, N1>, i_<I0, N2>>> const*>(&sum)).e;

		auto cancelled = (expand_cancelled_term<Ns, sizeof...(Ns)>(e, i_<I0, N0>{}, i_<I1, N1>{}, i_<I0, N2>{}, ops) + ...);
		return expr::sum<i_<I0, N0>, i_<I1, N1>>(e).template expand<0, i_<I1, 0>>(std::make_tuple(std::get<Ns>(ops)...)) - cancelled;
	}

	//! Expand the sum where the first index does not have one of the values substituted.
	/*!
	 * The full sum is expanded, where the expansion of jj assumes the iteration value starts at 0.
	 * Then sums are carried out with only jj = ii + n is expanded, they are all summed, and then
	 * ii is substituted as its usual value. This is then subtracted from the full sum, in order
	 * to get an expression where jj != ii. Note that N0 refers to the starting point in the
	 * list of functions that v_ii starts substituting. N2 is means that jj is not supposed to be
	 * ii + N0 + N2.
	 *
	 * \tparam N0 The global offset of index I0
	 * \tparam N1 The global offset of index I1
	 */
	template<typename E, int I0, int N0, int I, int N, typename... Us, size_t... Ns>
	auto expand_sum(symphas::internal::index_conditions<E, symphas::internal::index_relation<std::false_type, i_<I0, N0>, i_<I, N>>> const& sum,
		std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		auto e = (*static_cast<symphas::internal::sum_conditions<E, 
			symphas::internal::index_relation<std::false_type, i_<I0, N0>, i_<I, N>>> const*>(&sum)).e;

		auto cancelled = v_<i_<I, N>>{};
		return expr::sum<i_<I0, N0>>(e).template expand<0>(std::make_tuple(std::get<Ns>(ops)...)) - cancelled;
	}

	template<typename E, typename... Us, size_t... Ns,
		typename std::enable_if_t<!std::is_convertible_v<E, symphas::internal::index_conditions<>>, int> = 0>
	auto expand_sum(OpExpression<E> const& e, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return (*static_cast<E const*>(&e));
	}

	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpBinaryMul<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);
	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpOperatorCombination<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);
	template<typename A, typename B, typename E, typename... Us, size_t... Ns>
	auto expand_sum(OpCombination<A, B, E> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);
	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpOperatorChain<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);
	template<typename A, typename B, typename E, typename... Us, size_t... Ns>
	auto expand_sum(OpChain<A, B, E> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);
	template<typename Dd, typename V, typename E, typename Sp, typename... Us, size_t... Ns>
	auto expand_sum(OpFuncDerivative<Dd, V, E, Sp> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>);

	template<typename... As, size_t... Is, typename... Us, size_t... Ns>
	auto expand_sum_adds(OpAdd<As...> const& sums, std::index_sequence<Is...>, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return (expand_sum(expr::get<Is>(sums), ops, std::index_sequence<Ns...>{}) + ...);
	}

	template<typename... As, typename... Us, size_t... Ns>
	auto expand_sum(OpAdd<As...> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return expand_sum_adds(sums, std::make_index_sequence<sizeof...(As)>{}, ops, std::index_sequence<Ns...>{});
	}

	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpOperatorCombination<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return OpOperatorCombination(expand_sum(sums.f, ops, std::index_sequence<Ns...>{}), expand_sum(sums.g, ops, std::index_sequence<Ns...>{}));
	}

	template<typename A, typename B, typename E, typename... Us, size_t... Ns>
	auto expand_sum(OpCombination<A, B, E> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return expand_sum(sums.combination, ops, std::index_sequence<Ns...>{})(expand_sum(expr::get_enclosed_expression(sums), ops, std::index_sequence<Ns...>{}));
	}

	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpOperatorChain<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return OpOperatorChain(expand_sum(sums.f, ops, std::index_sequence<Ns...>{}), expand_sum(sums.g, ops, std::index_sequence<Ns...>{}));
	}

	template<typename A, typename B, typename E, typename... Us, size_t... Ns>
	auto expand_sum(OpChain<A, B, E> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return expand_sum(sums.combination, ops, std::index_sequence<Ns...>{})(expand_sum(expr::get_enclosed_expression(sums), ops, std::index_sequence<Ns...>{}));
	}

	template<typename Dd, typename V, typename E, typename Sp, typename... Us, size_t... Ns>
	auto expand_sum(OpFuncDerivative<Dd, V, E, Sp> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return expr::make_derivative<Dd>(expr::coeff(sums), expand_sum(expr::get_enclosed_expression(sums), ops, std::index_sequence<Ns...>{}));
	}

	template<typename A, typename B, typename... Us, size_t... Ns>
	auto expand_sum(OpBinaryMul<A, B> const& sums, std::tuple<Us...> const& ops, std::index_sequence<Ns...>)
	{
		return expand_sum(sums.a, ops, std::index_sequence<Ns...>{}) * expand_sum(sums.b, ops, std::index_sequence<Ns...>{});
	}

}


/* multiplication of a derivative object by a literal
 */
template<typename coeff_t, typename... Is, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& value, symphas::internal::index_conditions<E, Is...> const& sum)
{
	return symphas::internal::index_conditions(Is{}...)(value * 
		(*static_cast<symphas::internal::sum_conditions<E, Is...> const*>(&sum)).e);
}



