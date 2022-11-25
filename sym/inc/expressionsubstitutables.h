
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
				if constexpr (!symphas::lib::is_value_in_seq<M, Seq>::value)
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
		};


		template<typename, typename...>
		struct v_id_type;
		template<int N0, int P0, int... Ns, int... Ps>
		struct v_id_type<i_<N0, P0>, i_<Ns, Ps>...> : Symbol {};

		template<typename I0, typename... Is>
		using v_ = OpTerm<OpIdentity, v_id_type<I0, Is...>>;

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

namespace expr::internal
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


	template<typename... Ts>
	struct select_v_;

	template<int N0, int P0, typename... Vs>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<Vs...>;
	};

	template<int N0, int P0, typename... Vs, int... Ps, typename... Rest>
	struct select_v_<expr::symbols::i_<N0, P0>, symphas::lib::types_list<Vs...>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N0, Ps>...>, Rest...>>
	{
		using type = typename select_v_<
			expr::symbols::i_<N0, P0>,
			symphas::lib::types_list<Vs..., expr::symbols::v_id_type<expr::symbols::i_<N0, Ps>...>>,
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
template<expr::internal::CompoundOp Op, typename I, typename E, typename... Is>
struct OpCompound;

template<expr::internal::CompoundOp Op, typename I, typename E>
struct OpCompound<Op, I, E> : OpExpression<OpCompound<Op, I, E>>
{
protected:

	using all_var_types = typename expr::op_types<E>::type;
	template<typename V>
	static constexpr int index_of_type = symphas::lib::index_of_type<V, all_var_types>;
	
	using v_types = typename expr::internal::select_v_<I, all_var_types>::type;

	template<typename V>
	static constexpr int index_of_v = symphas::lib::index_of_type<V, v_types>;

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
			return typename expr::internal::identity_op<Op>::type{};
		}
		else
		{
			return expand_integer(
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
			std::make_index_sequence<sizeof...(Es) + 1>{}, 
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
		return expand_function<I0>(std::make_index_sequence<sizeof...(Es) + 1>{}, subs);
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

//protected:

	template<int... Is>
	auto expand_integer(std::integer_sequence<int, Is...>) const
	{
		if constexpr (index_of_type<I> >= 0)
		{
			return expr::internal::recursive_op_apply<Op>{}(expr::transform::swap_grid<I>(e, expr::make_integer<Is>())...);
		}
		else
		{
			auto ee = [&] (int) { return e; };
			return expr::internal::recursive_op_apply<Op>{}(ee(Is)...);
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

	template<size_t Is, typename Ex, typename... Es>
	auto substitute_functions(
		symphas::lib::types_list<>,
		Ex const& ex,
		std::tuple<Es...> const& subs) const
	{
		return ex;
	}

	template<size_t Is, typename Ex, typename... Es, int N0, int P0, int... Ns, int... Ps>
	auto substitute_functions(
		symphas::lib::types_list<
			expr::symbols::v_id_type<expr::symbols::i_<N0, P0>>, 
			expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>>...>, 
		Ex const& ex, 
		std::tuple<Es...> const& subs) const
	{
		return substitute_functions<Is>(
			symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<Ns, Ps>>...>{},
			expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<N0, P0>>>(ex, get_substituted<int(Is) + P0>(subs)),
			subs);
	}

	template<int I0, typename... Es, size_t... Is>
	auto expand_function(std::index_sequence<Is...>, std::tuple<Es...> const& subs) const
	{
		if constexpr (index_of_type<I> >= 0)
		{
			return expr::internal::recursive_op_apply<Op>{}(
				substitute_functions<Is>(
					v_types{},
					expr::transform::swap_grid<I>(e, expr::make_integer<int(Is) + I0>()),
					subs)...);
		}
		else
		{
			return expr::internal::recursive_op_apply<Op>{}(
				substitute_functions<Is>(v_types{}, e, subs)...);
		}
	}

	template<int I0>
	auto expand() const
	{
		return typename expr::internal::identity_op<Op>::type{};
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
template<expr::internal::CompoundOp Op, typename I0, typename E, typename I1, typename... Is>
struct OpCompound<Op, I0, E, I1, Is...> : OpExpression<OpCompound<Op, I0, E, I1, Is...>>
{

protected:

	using all_var_types = typename expr::op_types<E>::type;
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

	//! Expand the sum on a single list of functions using indices offset from the first.
	/*!
	 * Expand the sum, starting at the given index. The rest of the indices are given as
	 * offsets of the first, and indicate the starting position of the substitutions in the list.
	 * The offsets can only be positive, and do not loop around the list.
	 */
	template<int I, typename... ISs, typename T0, typename... Ts, typename Vs = expr::symbols::v_id_type<I0>,
		typename std::enable_if_t<(index_of_type<Vs> >= 0 && (std::is_same<expr::base_data_t<I0>, expr::base_data_t<ISs>>::value && ...)), int> = 0>
	auto expand(std::tuple<T0, Ts...> const& subs) const
	{
		return expand_iterated<I>(subs, symphas::lib::types_list<ISs...>{}, std::make_index_sequence<sizeof...(Ts) + 1>{});
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

	template<size_t N, typename... Ts, typename EE, int I00, int... I00s>
	auto expand_nested(std::tuple<> const& sub, EE const& expr0, std::integer_sequence<int, I00, I00s...>) const
	{
		return typename expr::internal::identity_op<Op>::type{};
	}

	template<size_t N, typename E0, typename... Es, typename EE, int I00, size_t... Js, typename Vs = expr::symbols::v_<I0, I1, Is...>,
		typename In = std::tuple_element_t<N, std::tuple<I0, I1, Is...>>, typename V0 = expr::symbols::v_<In>>
	auto expand_nested(std::tuple<E0, Es...> const& sub, EE const& expr0, std::integer_sequence<int, I00>, std::index_sequence<Js...>) const
	{
		if constexpr (index_of_type<In> >= 0)
		{
			return expr::internal::recursive_op_apply<Op>{}(
				expr::transform::swap_grid<Vs>(
					expr::transform::swap_grid<In>(expr0, expr::make_integer<int(Js) + I00>()),
					std::get<Js>(sub))...);
		}
		else
		{
			return expr::internal::recursive_op_apply<Op>{}(
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
			return expr::internal::recursive_op_apply<Op>{}(
				expand_nested<N + 1>(
					std::get<Js>(sub),
					expr::transform::swap_grid<In>(expr0, expr::make_integer<int(Js) + I00>()),
					std::integer_sequence<int, I00s...>{},
					std::make_index_sequence<std::tuple_size<std::tuple_element_t<Js, sub_t>>::value>{})...
				);
		}
		else
		{
			return expr::internal::recursive_op_apply<Op>{}(
				expand_nested<N + 1>(
					std::get<Js>(sub),
					expr0,
					std::integer_sequence<int, I00s...>{},
					std::make_index_sequence<std::tuple_size<std::tuple_element_t<Js, sub_t>>::value>{})...
				);
		}

	}

	template<int I, typename T0, typename... Ts, typename EE, int... Ns, int... Ps>
	auto construct_iterated(
		std::tuple<T0, Ts...> const& subs, 
		EE const& expr0, 
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>) const
	{
		return OpCompound<Op, I1, EE, Is...>(expr0).template expand<I, expr::symbols::i_<Ns, Ps>...>(subs);
	}

	template<int I, typename... Ts, typename EE, int... Ns, int... Ps>
	auto construct_iterated(
		std::tuple<> const& subs,
		EE const& expr0,
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>) const
	{
		return typename expr::internal::identity_op<Op>::type{};
	}

	template<int I, size_t J, int M0, int Q0, int N0, int P0, int... Ns, int... Ps, typename... Ts, typename V = expr::symbols::v_id_type<I0>>
	auto construct_iterated(
		std::tuple<Ts...> const& subs,
		expr::symbols::i_<M0, Q0>, 
		symphas::lib::types_list<expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>) const
	{
		using J0s = symphas::lib::types_list<expr::symbols::i_<M0, Ps - I>...>;

		if constexpr (index_of_type<V> >= 0)
		{
			if constexpr (index_of_type<I0> >= 0)
			{
				return construct_iterated<I + int(J)>(
					symphas::lib::get_tuple_ge<J + P0>(subs),
					expr::transform::swap_grid<V>(
						expr::transform::swap_grid<I0>(e, expr::make_integer<I + (int)J>()), 
						std::get<J>(subs)), 
					J0s{});
			}

			return construct_iterated<I + int(J)>(
				symphas::lib::get_tuple_ge<J + P0>(subs),
				expr::transform::swap_grid<V>(e, std::get<J>(subs)),
				J0s{});
		}
		else
		{
			return construct_iterated(symphas::lib::get_tuple_ge<J + P0>(subs), e, J0s{});
		}
	}

	template<int I, typename... Ts, int... Ns, int... Ps, size_t... Js>
	auto expand_iterated(std::tuple<Ts...> const& subs, 
		symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>, 
		std::index_sequence<Js...>) const
	{
		return expr::internal::recursive_op_apply<Op>{}(
			construct_iterated<I, Js>(subs, I1{}, symphas::lib::types_list<expr::symbols::i_<Ns, Ps>...>{})...);
	}

public:

	E e;				//!< The expression which is summed.

};


template<typename coeff_t, typename I, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t, OpCompound<expr::internal::CompoundOp::ADD, I, E> const& sum)
{
	return expr::sum<I>(coeff_t{} *sum.e);
}

template<typename T, typename I, typename E>
auto operator*(OpLiteral<T> const& a, OpCompound<expr::internal::CompoundOp::ADD, I, E> const& sum)
{
	return expr::sum<I>(a * sum.e);
}

template<typename T, typename I, typename E>
auto operator-(OpCompound<expr::internal::CompoundOp::ADD, I, E> const& sum)
{
	return expr::sum<I>(-sum.e);
}

template<typename coeff_t, typename I0, typename I1, typename... Is, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t, OpCompound<expr::internal::CompoundOp::ADD, I0, E, I1, Is...> const& sum)
{
	return expr::sum<I0, I1, Is...>(coeff_t{} *sum.e);
}

template<typename T, typename I0, typename I1, typename... Is, typename E>
auto operator*(OpLiteral<T> const& a, OpCompound<expr::internal::CompoundOp::ADD, I0, E, I1, Is...> const& sum)
{
	return expr::sum<I0, I1, Is...>(a * sum.e);
}

template<typename T, typename I0, typename I1, typename... Is, typename E>
auto operator-(OpCompound<expr::internal::CompoundOp::ADD, I0, E, I1, Is...> const& sum)
{
	return expr::sum<I0, I1, Is...>(-sum.e);
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
		return OpCompound<expr::internal::CompoundOp::ADD, I, E>(*static_cast<E const*>(&e));
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
		return OpCompound<expr::internal::CompoundOp::ADD, I0, E, I1, Is...>(*static_cast<E const*>(&e));
	}


	//! Defines a product of an expression with indices and substitutable expressions.
	/*!
	 * See expr::sum.
	 */
	template<typename I, typename E, typename>
	auto prod(OpExpression<E> const& e)
	{
		return OpCompound<expr::internal::CompoundOp::MUL, I, E>(*static_cast<E const*>(&e));
	}

	//! Defines a product of an expression with indices and substitutable expressions.
	/*!
	 * See expr::sum.
	 */
	template<typename I0, typename I1, typename... Is, typename E, typename>
	auto prod(OpExpression<E> const& e)
	{
		return OpCompound<expr::internal::CompoundOp::MUL, I0, E, I1, Is...>(*static_cast<E const*>(&e));
	}
}




