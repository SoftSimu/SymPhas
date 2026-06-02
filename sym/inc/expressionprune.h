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
 * PURPOSE: Defines ways to transform expressions, usually from one
 * expression to the other based on some rules or the type of expression.
 *
 * ***************************************************************************
 */

#pragma once

#include <fstream>
#include <iostream>
#include <type_traits>

#include "expressionswap.h"

// ******************************************************************************************

struct EvalHandler : BaseEvalHandler<EvalHandler> {
  template <typename E, typename assign_type, typename interval_type>
  void result(OpEvaluable<E> const& e, assign_type&& data,
              interval_type&& interval) const {
    expr::result(*static_cast<E const*>(&e), std::forward<assign_type>(data),
                 std::forward<interval_type>(interval));
  }
  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const& e, assign_type&& data) const {
    expr::result(*static_cast<E const*>(&e), std::forward<assign_type>(data));
  }
  template <typename E, typename assign_type, typename interval_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                  interval_type&& interval) const {
    expr::result_sum(*static_cast<E const*>(&e),
                     std::forward<assign_type>(assign),
                     std::forward<interval_type>(interval));
  }
  template <typename E, typename assign_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign) const {
    return expr::result_sum(*static_cast<E const*>(&e),
                            std::forward<assign_type>(assign));
  }
};
#ifdef USING_CUDA

struct KernelEvalHandler : BaseEvalHandler<KernelEvalHandler> {
  template <typename E, typename assign_type, typename interval_type>
  void result(OpEvaluable<E> const& e, assign_type&& assign,
              interval_type&& interval) const {
    using grid_type = expr::storage_type_t<E>;
    expr::evaluate_expression_trait<grid_type>(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign),
        std::forward<interval_type>(interval), false);
  }
  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const& e, assign_type&& assign) const {
    using grid_type = expr::storage_type_t<E>;
    expr::evaluate_expression_trait<grid_type>(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign),
        expr::iterable_domain(*static_cast<E const*>(&e)), false);
  }
  template <typename E, typename assign_type, typename interval_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                  interval_type&& interval) const {
    using grid_type = expr::storage_type_t<E>;
    expr::result_sum_trait<grid_type>{}.result(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign),
        std::forward<interval_type>(interval), false);
  }
  template <typename E, typename assign_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign) const {
    using grid_type = expr::storage_type_t<E>;
    expr::result_sum_trait<grid_type>{}.result(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign),
        expr::iterable_domain(*static_cast<E const*>(&e)), false);
  }

  // template <typename condition_t, typename... condition_ts, typename E,
  //           typename assign_type>
  // void result_by_term(OpEvaluable<E> const& e, assign_type&& assign) const {
  //   result(*static_cast<E const*>(&e), std::forward<assign_type>(assign));
  // }

  // template <typename condition_t, typename... condition_ts, typename... Es,
  //           typename assign_type>
  // void result_by_term(OpAdd<Es...> const& e, assign_type&& assign) const {
  //   if constexpr ((expr::satisfies<Es,
  //                                  expr::or_<condition_t, condition_ts...>>
  //                                  ||
  //                  ...)) {
  //     result_only<expr::not_<expr::or_<condition_t, condition_ts...>>>(
  //         e, std::forward<assign_type>(assign));
  //     result_by_term(e, std::forward<assign_type>(assign),
  //                    symphas::lib::types_list<condition_t,
  //                    condition_ts...>{},
  //                    std::make_index_sequence<sizeof...(Es)>{});
  //   } else {
  //     result(e, std::forward<assign_type>(assign));
  //   }
  // }

  ~KernelEvalHandler() { CHECK_CUDA_ERROR(cudaDeviceSynchronize()); }
};
#endif

namespace expr {

namespace internal {
template <typename grid_type>
struct select_eval_handler_trait_on_grid {
  using type = EvalHandler;
};
#ifdef USING_CUDA

template <typename T, size_t D>
struct select_eval_handler_trait_on_grid<GridCUDA<T, D>> {
  using type = KernelEvalHandler;
};
template <typename T, size_t D>
struct select_eval_handler_trait_on_grid<BoundaryGridCUDA<T, D>> {
  using type = KernelEvalHandler;
};
template <typename T, size_t D>
struct select_eval_handler_trait_on_grid<RegionalGridCUDA<T, D>> {
  using type = KernelEvalHandler;
};

#endif
}  // namespace internal

template <typename E>
struct select_eval_handler_trait {
  using type = typename expr::internal::select_eval_handler_trait_on_grid<
      expr::storage_type_t<E>>::type;
};
template <typename E>
using eval_handler_type = typename select_eval_handler_trait<E>::type;
}  // namespace expr

namespace expr::prune {
namespace {
inline void _update(...) {}

template <typename E, typename eval_handler_type, typename... condition_ts>
inline void _update(OpExpression<E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {}

template <typename T, typename eval_handler_type, typename... condition_ts>
inline void _update(Block<T>&, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {}
template <size_t N, typename T, typename eval_handler_type,
          typename... condition_ts>
inline void _update(MultiBlock<N, T>&, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {}

template <typename V, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpExponential<V, E>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpExponential<V, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    expr::exp_key_t X, typename V, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpPow<X, V, E>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(OpPow<X, V, E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename V, typename E1, typename E2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpConvolution<V, E1, E2>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpConvolution<V, E1, E2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename V, typename E, typename F, typename... Args,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpFunction<V, E, F, Args...>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpFunction<V, E, F, Args...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <auto f, typename V, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpFunctionApply<f, V, E>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpFunctionApply<f, V, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename V, size_t D, template <typename, size_t> typename grid_type, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D, grid_type>, E>,
                                     expr::or_<condition_ts...>>,
                     int> = 0>
inline void _update(OpConvolution<V, GaussianSmoothing<D, grid_type>, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename V, size_t D, template <typename, size_t> typename grid_type, typename G, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D, grid_type>,
                                                   OpTerm<OpIdentity, G>>,
                                     expr::or_<condition_ts...>>,
                     int> = 0>
inline void _update(
    OpConvolution<V, GaussianSmoothing<D, grid_type>, OpTerm<OpIdentity, G>>& e,
    eval_handler_type const& eval_handler,
    symphas::lib::types_list<condition_ts...>);
template <size_t O, typename V, typename E, typename G,
          typename eval_handler_type, typename... condition_ts>
inline void _update(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>& e,
    eval_handler_type const& eval_handler,
    symphas::lib::types_list<condition_ts...>);
template <typename Dd, typename V, typename G, typename Sp,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>,
                              expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename Dd, typename V, typename E, typename Sp,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, E, Sp>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpDerivative<Dd, V, E, Sp>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename V, typename E, typename T, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpIntegral<V, E, T>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpIntegral<V, E, T>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename... Es, typename eval_handler_type, typename... condition_ts,
    size_t... Is,
    std::enable_if_t<
        expr::satisfies<OpAddList<Es...>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(OpAddList<Es...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>,
                    std::index_sequence<Is...>);
template <
    typename... Es, typename eval_handler_type, typename... condition_ts,
    std::enable_if_t<(expr::satisfies<Es, expr::or_<condition_ts...>> || ...),
                     int> = 0>
inline void _update(OpAdd<Es...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpBinaryMul<A1, A2>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpBinaryMul<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpBinaryDiv<A1, A2>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpBinaryDiv<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpOperatorCombination<A1, A2>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpOperatorCombination<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpCombination<A1, A2, E>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpCombination<A1, A2, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpOperatorChain<A1, A2>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpOperatorChain<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename A1, typename A2, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpChain<A1, A2, E>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpChain<A1, A2, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename G, typename V, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpMap<G, V, E>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(OpMap<G, V, E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename... Ts, typename eval_handler_type, typename... condition_ts,
          size_t... Is,
          std::enable_if_t<
              expr::satisfies<OpTermsList<Ts...>, expr::or_<condition_ts...>>,
              int> = 0>
inline void _update(OpTermsList<Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>,
                    std::index_sequence<Is...>);
template <
    typename... Ts, typename eval_handler_type, typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpTerms<Ts...>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(OpTerms<Ts...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename G, expr::exp_key_t X, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<Term<G, X>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(Term<G, X>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename V, typename sub_t, typename eval_t,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpSymbolicEval<V, sub_t, eval_t>,
                                           expr::or_<condition_ts...>>,
                           int> = 0>
inline void _update(OpSymbolicEval<V, sub_t, eval_t>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename E, typename eval_handler_type, typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpOptimized<E>, expr::or_<condition_ts...>>, int> = 0>
inline void _update(OpOptimized<E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <
    typename E, typename... Ts, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<expr::satisfies<E, expr::or_<condition_ts...>>, int> = 0>
inline void _update(SymbolicFunction<E, Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename... Ts, typename eval_handler_type, typename... condition_ts>
inline void _update(SymbolicCase<Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);

// template<typename Op, typename E, typename Inds>
// inline void _update(SymbolicSeries<Op, E, Inds>& e);
// template<typename Op, typename... Ts, typename E, int... I0s, int... P0s,
//	typename... T1s, typename... T2s, typename... Is>
// inline void _update(SymbolicSeries<Op,
// Substitution<SymbolicDataArray<Ts>...>, 	symphas::lib::types_list<E,
//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
//		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
//		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
//	>>& e);

// template<typename... Ts>
// inline void _update(SymbolicCase<Ts...>& e);
// template<size_t Z, typename G>
// inline void _update(Variable<Z, G>& e);
// template<typename G>
// inline void _update(NamedData<G>& e);
// template<Axis ax, typename G>
// inline void _update(VectorComponent<ax, G>& e);
// template<typename G>
// inline void _update(symphas::ref<G>& e);

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type,
          typename eval_handler_type, typename... condition_ts>
inline void _update(NoiseData<nt, T, D, grid_type>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);
template <typename T, typename eval_handler_type, typename... condition_ts>
inline void _update(SymbolicData<T>&, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>);

template <
    typename V, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpExponential<V, E>, expr::or_<condition_ts...>>, int>>
inline void _update(OpExponential<V, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
}

template <expr::exp_key_t X, typename V, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpPow<X, V, E>, expr::or_<condition_ts...>>, int>>
inline void _update(OpPow<X, V, E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
}

template <typename V, typename E1, typename E2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpConvolution<V, E1, E2>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpConvolution<V, E1, E2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.a, eval_handler, symphas::lib::types_list<condition_ts...>{});
  _update(e.b, eval_handler, symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename V, typename E, typename F, typename... Args,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpFunction<V, E, F, Args...>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpFunction<V, E, F, Args...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
}

template <auto f, typename V, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpFunctionApply<f, V, E>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpFunctionApply<f, V, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
}

template <
    typename V, size_t D, template <typename, size_t> typename grid_type, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D, grid_type>, E>,
                                     expr::or_<condition_ts...>>,
                     int>>
inline void _update(OpConvolution<V, GaussianSmoothing<D, grid_type>, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <
    typename V, size_t D, template <typename, size_t> typename grid_type, typename G, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D, grid_type>,
                                                   OpTerm<OpIdentity, G>>,
                                     expr::or_<condition_ts...>>,
                     int>>
inline void _update(
    OpConvolution<V, GaussianSmoothing<D, grid_type>, OpTerm<OpIdentity, G>>& e,
    eval_handler_type const& eval_handler,
    symphas::lib::types_list<condition_ts...>) {
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

/* derivative pruning
 */

template <size_t O, typename V, typename E, typename G,
          typename eval_handler_type, typename... condition_ts>
inline void _update(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>& e,
    eval_handler_type const& eval_handler,
    symphas::lib::types_list<condition_ts...>) {}

template <typename Dd, typename V, typename G, typename Sp,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>,
                              expr::or_<condition_ts...>>,
              int>>
inline void _update(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename Dd, typename V, typename E, typename Sp,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, E, Sp>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpDerivative<Dd, V, E, Sp>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <
    typename V, typename E, typename T, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpIntegral<V, E, T>, expr::or_<condition_ts...>>, int>>
inline void _update(OpIntegral<V, E, T>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

/* binary op pruning
 */

template <
    typename... Es, typename eval_handler_type, typename... condition_ts,
    size_t... Is,
    std::enable_if_t<
        expr::satisfies<OpAddList<Es...>, expr::or_<condition_ts...>>, int>>
inline void _update(OpAddList<Es...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>,
                    std::index_sequence<Is...>) {
  (_update(expr::get<Is>(e), eval_handler,
           symphas::lib::types_list<condition_ts...>{}),
   ...);
}

template <typename... Es, typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<
              (expr::satisfies<Es, expr::or_<condition_ts...>> || ...), int>>
inline void _update(OpAdd<Es...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(*static_cast<OpAddList<Es...>*>(&e), eval_handler,
          symphas::lib::types_list<condition_ts...>{},
          std::make_index_sequence<sizeof...(Es)>{});
}

template <
    typename A1, typename A2, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpBinaryMul<A1, A2>, expr::or_<condition_ts...>>, int>>
inline void _update(OpBinaryMul<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.a, eval_handler, symphas::lib::types_list<condition_ts...>{});
  _update(e.b, eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <
    typename A1, typename A2, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpBinaryDiv<A1, A2>, expr::or_<condition_ts...>>, int>>
inline void _update(OpBinaryDiv<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.a, eval_handler, symphas::lib::types_list<condition_ts...>{});
  _update(e.b, eval_handler, symphas::lib::types_list<condition_ts...>{});
}

/* operator pruning
 */

template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpOperatorCombination<A1, A2>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpOperatorCombination<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.f, eval_handler, symphas::lib::types_list<condition_ts...>{});
  _update(e.g, eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename A1, typename A2, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpCombination<A1, A2, E>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpCombination<A1, A2, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename A1, typename A2, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpOperatorChain<A1, A2>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpOperatorChain<A1, A2>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.f, eval_handler, symphas::lib::types_list<condition_ts...>{});
  _update(e.g, eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <
    typename A1, typename A2, typename E, typename eval_handler_type,
    typename... condition_ts,
    std::enable_if_t<
        expr::satisfies<OpChain<A1, A2, E>, expr::or_<condition_ts...>>, int>>
inline void _update(OpChain<A1, A2, E>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  //_update(e.combination);
  // if constexpr (expr::has_state<E>::value)
  //{
  //	_update(expr::get_enclosed_expression(e));
  //}
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename G, typename V, typename E, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpMap<G, V, E>, expr::or_<condition_ts...>>, int>>
inline void _update(OpMap<G, V, E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(expr::get_enclosed_expression(e), eval_handler,
          symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <
    typename... Ts, typename eval_handler_type, typename... condition_ts,
    size_t... Is,
    std::enable_if_t<
        expr::satisfies<OpTermsList<Ts...>, expr::or_<condition_ts...>>, int>>
inline void _update(OpTermsList<Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>,
                    std::index_sequence<Is...>) {
  (_update(expr::get<Is>(e), eval_handler,
           symphas::lib::types_list<condition_ts...>{}),
   ...);
}

template <typename... Ts, typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpTerms<Ts...>, expr::or_<condition_ts...>>, int>>
inline void _update(OpTerms<Ts...>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  //_update(*static_cast<OpTermsList<Ts...>*>(&e),
  // std::make_index_sequence<sizeof...(Ts)>{});
}

template <typename G, expr::exp_key_t X, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<Term<G, X>, expr::or_<condition_ts...>>, int>>
inline void _update(Term<G, X>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  //_update(e.data());
}

template <typename V, typename sub_t, typename eval_t,
          typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<expr::satisfies<OpSymbolicEval<V, sub_t, eval_t>,
                                           expr::or_<condition_ts...>>,
                           int>>
inline void _update(OpSymbolicEval<V, sub_t, eval_t>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.f, eval_handler, symphas::lib::types_list<condition_ts...>{});
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename E, typename eval_handler_type, typename... condition_ts,
          std::enable_if_t<
              expr::satisfies<OpOptimized<E>, expr::or_<condition_ts...>>, int>>
inline void _update(OpOptimized<E>& e, eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  e.update(eval_handler, symphas::lib::types_list<condition_ts...>{});
}

template <typename E, typename... Ts, typename eval_handler_type,
          typename... condition_ts,
          std::enable_if_t<expr::satisfies<E, expr::or_<condition_ts...>>, int>>
inline void _update(SymbolicFunction<E, Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {
  _update(e.e, eval_handler, symphas::lib::types_list<condition_ts...>{});
}

// template<typename Op, typename E, typename Inds>
// inline void _update(SymbolicSeries<Op, E, Inds>& e) {}

// template<typename Op, typename... Ts, typename E, int... I0s, int... P0s,
//	typename... T1s, typename... T2s, typename... Is>
// inline void _update(SymbolicSeries<Op,
// Substitution<SymbolicDataArray<Ts>...>, 	symphas::lib::types_list<E,
//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
//		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
//		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
//	>>& e)
//{
//	//e.update();
// }

template <typename... Ts, typename eval_handler_type, typename... condition_ts>
inline void _update(SymbolicCase<Ts...>& e,
                    eval_handler_type const& eval_handler,
                    symphas::lib::types_list<condition_ts...>) {}

// template<size_t Z, typename G>
// inline void _update(Variable<Z, G>& e)
//{
//     _update(*static_cast<G*>(&e));
// }

// template<typename G>
// inline void _update(NamedData<G>& e)
//{
//     _update(*static_cast<G*>(&e));
// }

// template<Axis ax, typename G>
// inline void _update(VectorComponent<ax, G>& e)
//{
//     _update(*static_cast<G*>(&e));
// }

// template<typename G>
// inline void _update(symphas::ref<G>& e)
//{
//     _update(e.get());
// }

// template<expr::NoiseType nt, typename T, size_t D>
// inline void _update(NoiseData<nt, T, D>& data)
//{
//     data.update();
// }

//     template<typename T>
//     inline void _update(SymbolicData<T>& data)
//     {
// if (data.data != nullptr)
//{
//	_update(*data.data);
//}
//     }

}  // namespace

//! Update underlying the given expression.
/*!
 * For expressions which store intermediate data, such as derivatives, this
 * data must be updated before the expression can be evaluated. This will
 * traverse the expression tree and perform all necessary updating.
 *
 * \param e The expression that is updated.
 */
template <typename... condition_ts, typename E, typename eval_handler_type>
inline void update(E& e,
                   BaseEvalHandler<eval_handler_type> const& eval_handler) {
  //e.allocate();
  TIME_THIS_CONTEXT_LIFETIME(expression_update)
  _update(e, *static_cast<eval_handler_type const*>(&eval_handler),
          symphas::lib::types_list<condition_ts...>{});
}

template <typename... condition_ts, typename E>
inline void update(E& e) {
  e.allocate();
  _update(e, EvalHandler{}, symphas::lib::types_list<condition_ts...>{});
}

}  // namespace expr::prune
