
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

#include "expressionconditions.h"
#include "expressionproperties.h"
#include "expressionsprint.h"

// ******************************************************************************************

namespace expr {
namespace compile_trait {
template <typename condition_t>
struct separate_by_trait {
  template <typename E>
  auto operator()(E&&);
};
}  // namespace compile_trait

template <typename L, typename R>
struct equation_ptr_list_type {
  L* left;
  R* right;
};

#define PARALLELIZATION_CUTOFF_COUNT 1000

struct forward_value {
  template <typename T>
  decltype(auto) operator()(T && value) {
    return std::forward<T>(value);
  }
};

template <typename assign_type>
struct evaluate_expression_trait {
  template <typename E, size_t D>
  evaluate_expression_trait(OpEvaluable<E> const& e, assign_type& data,
                            grid::region_interval<D> const& interval) {
    symphas::data_iterator_group it(data, interval);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      std::transform(
          std::execution::par_unseq,
          static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
          static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
          forward_value{});
    else
#endif
      std::transform(
          static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
          static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
          forward_value{});
  }
  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const& e, assign_type& data,
                            grid::region_interval<0> const& interval) {
    auto data_region = expr::iterable_domain(data);
    if (grid::length(data_region) > 1) {
      result(*static_cast<E const*>(&e), data, data_region);
    } else {
      result(*static_cast<E const*>(&e), data, 1);
    }
  }

  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const& e, assign_type& data,
                            len_type len) {
    symphas::data_iterator it(data);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      std::transform(std::execution::par_unseq,
                     static_cast<const E*>(&e)->begin(),
                     static_cast<const E*>(&e)->end(len), it, forward_value{});
    else
#endif
      std::transform(static_cast<const E*>(&e)->begin(),
                     static_cast<const E*>(&e)->end(len), it, forward_value{});
  }

  template <typename E>
  evaluate_expression_trait(equation_ptr_list_type<assign_type, E>* eq_list,
                            len_type len) {
    auto range = symphas::parallel::get_index_range(len);
    iter_type start = range.first;
    iter_type end = range.second;

    SYMPHAS_OMP_PARALLEL_DIRECTIVE
    for (iter_type i = start; i < end; ++i) {
      evaluate_expression_trait(*eq_list[i].right, *eq_list[i].left,
                                expr::iterable_domain(*eq_list[i].left));
    }
  }
};

template <typename assign_type>
struct evaluate_expression_trait<assign_type&>
    : evaluate_expression_trait<assign_type> {
  using evaluate_expression_trait<assign_type>::evaluate_expression_trait;
};

//! Evaluate the expression into the underlying data member.
/*!
 * The expression must be iterable over the entire given length.
 *
 * \param e Expression that is evaluated.
 * \param data The array containing the result of the expression.
 * \param len The number of elements in the array.
 */
template <typename E, typename assign_type>
void result(OpEvaluable<E> const& e, assign_type&& data, len_type len) {
  evaluate_expression_trait<assign_type>(*static_cast<E const*>(&e),
                                         std::forward<assign_type>(data), len);
}

template <typename E, typename assign_type, size_t D>
void result(OpEvaluable<E> const& e, assign_type&& data,
            grid::region_interval<D> const& interval) {
  evaluate_expression_trait<assign_type>(
      *static_cast<E const*>(&e), std::forward<assign_type>(data), interval);
  /*symphas::data_iterator_group it(std::forward<assign_type>(data),
interval);

  if (params::parallelization)
    std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par_unseq,
#endif
        static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
        static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
        forward_value{});
  else
    std::transform(static_cast<const E*>(&e)->begin(symphas::it_grp,
interval), static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
forward_value{});*/
}

template <typename E, typename assign_type, size_t D>
void result(OpEvaluable<E> const& e, assign_type&& data,
            grid::region_interval_multiple<D> const& regions) {
  for (grid::region_interval<D> region : regions) {
    result(*static_cast<E const*>(&e), std::forward<assign_type>(data), region);
  }
}

template <typename E, typename assign_type>
void result(OpEvaluable<E> const& e, assign_type&& data,
            grid::region_interval<0> const& interval) {
  evaluate_expression_trait<assign_type>(
      *static_cast<E const*>(&e), std::forward<assign_type>(data), interval);
  // auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
  // if (grid::length(data_region) > 1) {
  //  result(*static_cast<E const*>(&e), std::forward<assign_type>(data),
  //         data_region);
  //} else {
  //  result(*static_cast<E const*>(&e), std::forward<assign_type>(data), 1);
  //}
}

template <typename E, typename assign_type>
void result(OpEvaluable<E> const& e, assign_type&& data, grid::region_empty) {}

template <typename E, typename assign_type>
void result(OpEvaluable<E> const& e, assign_type&& data) {
  result(*static_cast<E const*>(&e), std::forward<assign_type>(data),
         expr::iterable_domain(*static_cast<E const*>(&e)));
}

template <typename L, typename R>
void result(std::pair<L, R>& evaluate) {
  auto&& [lhs, rhs] = evaluate;
  result(rhs, lhs);
}

template <typename L, typename R>
void result(equation_ptr_list_type<L, R>* evaluate, len_type len) {
  evaluate_expression_trait<L>(evaluate, len);
}

template <typename G, typename R>
void result(std::pair<OpTerm<OpIdentity, G>, R> const& evaluate) {
  auto&& [lhs, rhs] = evaluate;
  result(rhs, BaseData<G>::get(expr::get<1>(lhs).data()));
}

template <typename assign_type>
struct accumulate_expression_trait {
  //! Add the result of the expression into the underlying data member.
  /*!
   * The expression is evaluated and the result is added to the existing
   * values in the data array.
   *
   * \param e Expression that is evaluated.
   * \param data The array of data.
   * \param len The length of the array.
   */
  template <typename E>
  accumulate_expression_trait(OpEvaluable<E> const& e, assign_type& data,
                              len_type len) {
    symphas::data_iterator it(std::forward<assign_type>(data));

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      std::transform(std::execution::par_unseq,
                     static_cast<const E*>(&e)->begin(),
                     static_cast<const E*>(&e)->end(len), it, it,
                     [](auto expr_value, auto data_value) {
                       return data_value + expr_value;
                     });
    else
#endif
      std::transform(static_cast<const E*>(&e)->begin(),
                     static_cast<const E*>(&e)->end(len), it, it,
                     [](auto expr_value, auto data_value) {
                       return data_value + expr_value;
                     });
  }

  template <typename E, size_t D>
  accumulate_expression_trait(OpEvaluable<E> const& e, assign_type& data,
                              grid::region_interval<D> const& interval) {
    symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

    if (grid::length<D>(interval) <= PARALLELIZATION_CUTOFF_COUNT) {
      auto start = static_cast<const E*>(&e)->begin(symphas::it_grp, interval);
      auto end = static_cast<const E*>(&e)->end(symphas::it_grp, interval);

      for (auto eit = start; eit < end; ++eit, ++it) {
        *it = *it + *eit;
      }
    } else {
#ifdef EXECUTION_HEADER_AVAILABLE
      if (params::parallelization)
        std::transform(
            std::execution::par_unseq,
            static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
            static_cast<const E*>(&e)->end(symphas::it_grp, interval), it, it,
            [](auto expr_value, auto data_value) {
              return data_value + expr_value;
            });
      else
#endif
        std::transform(
            static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
            static_cast<const E*>(&e)->end(symphas::it_grp, interval), it, it,
            [](auto expr_value, auto data_value) {
              return data_value + expr_value;
            });
    }
  }

  template <typename V, typename... Gs, expr::exp_key_t... Xs, size_t D>
  accumulate_expression_trait(OpTerms<V, Term<Gs, Xs>...> const& e,
                              assign_type& data,
                              grid::region_interval<D> const& interval) {
    symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

    auto start = e.begin(symphas::it_grp, interval);
    auto end = e.end(symphas::it_grp, interval);

    for (auto eit = start; eit < end; ++eit, ++it) {
      *it = *it + *eit;
    }
  }
};

template <typename assign_type>
struct accumulate_expression_trait<assign_type&>
    : accumulate_expression_trait<assign_type> {
  using accumulate_expression_trait<assign_type>::accumulate_expression_trait;
};

//! Add the result of the expression into the underlying data member.
/*!
 * The expression is evaluated and the result is added to the existing
 * values in the data array.
 *
 * \param e Expression that is evaluated.
 * \param data The array of data.
 * \param len The length of the array.
 */
template <typename E, typename assign_type>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data,
                       len_type len) {
  accumulate_expression_trait<assign_type>(
      *static_cast<E const*>(&e), std::forward<assign_type>(data), len);
}

template <typename E, typename assign_type, size_t D>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data,
                       grid::region_interval<D> const& interval) {
  accumulate_expression_trait<assign_type>(
      *static_cast<E const*>(&e), std::forward<assign_type>(data), interval);
}

template <typename V, typename... Gs, expr::exp_key_t... Xs,
          typename assign_type, size_t D>
void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data,
                       grid::region_interval<D> const& interval) {
  accumulate_expression_trait<assign_type>(e, std::forward<assign_type>(data),
                                           interval);
}

template <typename assign_type>
void result_accumulate(OpVoid, assign_type&& data) {}
template <typename assign_type>
void result_accumulate(OpVoid, assign_type&& data, len_type) {}
template <typename assign_type>
void result_accumulate(OpVoid, assign_type&& data, grid::region_interval<0>) {}
template <typename assign_type, size_t D>
void result_accumulate(OpVoid, assign_type&& data, grid::region_interval<D>) {}
template <typename assign_type, size_t D>
void result_accumulate(OpVoid, assign_type&& data,
                       grid::region_interval_multiple<D>) {}

template <typename V, typename... Gs, expr::exp_key_t... Xs,
          typename assign_type, size_t D>
void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data,
                       grid::region_interval_multiple<D> const& regions) {
  for (grid::region_interval<D> region : regions) {
    result_accumulate(e, std::forward<assign_type>(data), region);
  }
}

template <typename E, typename assign_type, size_t D>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data,
                       grid::region_interval_multiple<D> const& regions) {
  for (grid::region_interval<D> region : regions) {
    result_accumulate(*static_cast<E const*>(&e),
                      std::forward<assign_type>(data), region);
  }
}

template <typename E, typename assign_type>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data,
                       grid::region_interval<0> const& interval) {
  auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
  if (grid::length(data_region) > 1) {
    result_accumulate(*static_cast<E const*>(&e),
                      std::forward<assign_type>(data), data_region);
  } else {
    result_accumulate(*static_cast<E const*>(&e),
                      std::forward<assign_type>(data), 1);
  }
}

template <typename V, typename... Gs, expr::exp_key_t... Xs,
          typename assign_type>
void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data,
                       grid::region_interval<0> const& interval) {
  auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
  if (grid::length(data_region) > 1) {
    result_accumulate(e, std::forward<assign_type>(data), data_region);
  } else {
    result_accumulate(e, std::forward<assign_type>(data), 1);
  }
}

template <typename E, typename assign_type>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data,
                       grid::region_empty) {}

template <typename E, typename assign_type>
void result_accumulate(OpEvaluable<E> const& e, assign_type&& data) {
  result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data),
                    expr::iterable_domain(*static_cast<E const*>(&e)));
}

template <typename data_type>
struct result_sum_trait {
  template <typename E, typename assign_type, size_t D>
  void result(OpEvaluable<E> const& e, assign_type&& assign,
              grid::region_interval<D> const& interval) {
    if (grid::length<D>(interval) <= PARALLELIZATION_CUTOFF_COUNT) {
      auto start = static_cast<const E*>(&e)->begin(symphas::it_grp, interval);
      auto end = static_cast<const E*>(&e)->end(symphas::it_grp, interval);

      auto reduce = expr::eval_type_t<E>{};
      for (auto it = start; it < end; ++it) {
        reduce = reduce + *it;
      }
      assign = reduce;
    } else {
      auto reduce =
#ifdef EXECUTION_HEADER_AVAILABLE
          (params::parallelization)
              ? std::reduce(
                    std::execution::par_unseq,
                    static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
                    static_cast<const E*>(&e)->end(symphas::it_grp, interval))
              :
#endif
              std::reduce(
                  static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
                  static_cast<const E*>(&e)->end(symphas::it_grp, interval));
      assign = reduce;
    }
  }

  template <typename V, typename... Gs, expr::exp_key_t... Xs,
            typename assign_type, size_t D>
  void result(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& assign,
              grid::region_interval<D> const& interval) {
    auto start = e.begin(symphas::it_grp, interval);
    auto end = e.end(symphas::it_grp, interval);

    auto reduce = expr::eval_type_t<OpTerms<V, Term<Gs, Xs>...>>{};
    for (auto it = start; it < end; ++it) {
      reduce = reduce + *it;
    }
    assign = reduce;
  }

  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const& e, assign_type&& assign, len_type len) {
#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      assign = std::reduce(std::execution::par_unseq,
                           static_cast<const E*>(&e)->begin(),
                           static_cast<const E*>(&e)->end(len));
    else
#endif
      assign = std::reduce(static_cast<const E*>(&e)->begin(),
                           static_cast<const E*>(&e)->end(len));
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const& e, assign_type&& assign, len_type len,
              std::index_sequence<Is...>) {
    auto start = symphas::reduce_seq_iterator(expr::get<Is>(e).begin()...);
    auto end = symphas::reduce_seq_iterator(expr::get<Is>(e).end(len)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      assign = std::reduce(std::execution::par_unseq, start, end);
    else
#endif
      assign = std::reduce(start, end);
  }

  template <typename... Es, typename assign_type, size_t D, size_t... Is>
  void result(OpAdd<Es...> const& e, assign_type&& assign,
              grid::region_interval<D> const& interval,
              std::index_sequence<Is...>) {
    auto start = symphas::reduce_iterator(
        expr::get<Is>(e).begin(symphas::it_grp, interval)...);
    auto end = symphas::reduce_iterator(
        expr::get<Is>(e).end(symphas::it_grp, interval)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      assign = std::reduce(std::execution::par_unseq, start, end);
    else
#endif
      assign = std::reduce(start, end);
  }
};

//! Evaluate the expression into the underlying data member.
/*!
 * The expression must be iterable over the entire given length.
 *
 * \param e Expression that is evaluated.
 * \param data The array containing the result of the expression.
 * \param len The number of elements in the array.
 */
template <typename E, typename assign_type>
void result_sum(OpEvaluable<E> const& e, assign_type&& assign, len_type len) {
  result_sum_trait<expr::storage_type_t<E>>{}.result(
      *static_cast<E const*>(&e), std::forward<assign_type>(assign), len);
}

template <typename assign_type>
void result_sum(OpVoid, assign_type&& assign) {
  std::forward<assign_type>(assign) = OpVoid{};
}

template <typename assign_type>
void result_sum(OpVoid, assign_type&& assign, len_type) {
  std::forward<assign_type>(assign) = OpVoid{};
}
template <typename assign_type>
void result_sum(OpVoid, assign_type&& assign, grid::region_interval<0>) {
  std::forward<assign_type>(assign) = OpVoid{};
}
template <typename assign_type, size_t D>
void result_sum(OpVoid, assign_type&& assign, grid::region_interval<D>) {
  std::forward<assign_type>(assign) = OpVoid{};
}
template <size_t D, typename assign_type>
void result_sum(OpVoid, assign_type&& assign,
                grid::region_interval_multiple<D>) {
  std::forward<assign_type>(assign) = OpVoid{};
}

template <typename E, typename assign_type>
void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                grid::region_interval<0> const& interval) {
  result_sum(*static_cast<E const*>(&e), std::forward<assign_type>(assign), 1);
}

template <typename V, typename... Gs, expr::exp_key_t... Xs,
          typename assign_type, size_t D>
void result_sum(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& assign,
                grid::region_interval<D> const& interval) {
  result_sum_trait<expr::storage_type_t<OpTerms<V, Term<Gs, Xs>...>>>{}.result(
      e, std::forward<assign_type>(assign), interval);
}

template <typename V, typename... Gs, expr::exp_key_t... Xs,
          typename assign_type, size_t D>
void result_sum(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& assign,
                grid::region_interval_multiple<D> const& regions) {
  expr::eval_type_t<OpTerms<V, Term<Gs, Xs>...>> sum{};
  std::forward<assign_type>(assign) = sum;
  for (grid::region_interval<D> region : regions) {
    result_sum(e, sum, region);
    std::forward<assign_type>(assign) += sum;
  }
}

template <typename E, typename assign_type, size_t D>
void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                grid::region_interval<D> const& interval) {
  result_sum_trait<expr::storage_type_t<E>>{}.result(
      *static_cast<E const*>(&e), std::forward<assign_type>(assign), interval);
}

template <typename E, typename assign_type, size_t D>
void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                grid::region_interval_multiple<D> const& regions) {
  expr::eval_type_t<E> sum{};
  std::forward<assign_type>(assign) = sum;
  for (grid::region_interval<D> region : regions) {
    result_sum(*static_cast<E const*>(&e), sum, region);
    std::forward<assign_type>(assign) += sum;
  }
}

template <typename E, typename assign_type>
void result_sum(OpEvaluable<E> const& e, assign_type&& assign) {
  TIME_THIS_EXPRESSION_LIFETIME(
      iterable_domain,
      auto r = expr::iterable_domain(*static_cast<E const*>(&e));)
  result_sum(*static_cast<E const*>(&e), std::forward<assign_type>(assign),
             expr::iterable_domain(*static_cast<E const*>(&e)));
}

struct matches_series;
struct matches_mul;
struct matches_div;
struct matches_term;
struct matches_integral;
struct matches_derivative;
struct matches_operator;

template <typename matches_t>
struct matching_in_mul;

//
// template <typename condition_t, typename... condition_ts, typename E,
//          typename assign_type>
// void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign,
//                        grid::region_interval<0> const& interval) {
//  using namespace symphas::internal;
//  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
//                 expression_satisfies_condition<E, condition_t>)) {
//    result_sum_by_term(*static_cast<E const*>(&e),
//                       std::forward<assign_type>(assign), 1);
//  } else {
//    expr::eval_type_t<E> sum{};
//    std::forward<assign_type>(assign) = sum;
//  }
//}
//
// template <typename condition_t, typename... condition_ts, typename E,
//          typename assign_type, size_t D>
// void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign,
//                        grid::region_interval<D> const& interval) {
//  using namespace symphas::internal;
//  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
//                 expression_satisfies_condition<E, condition_t>)) {
//    expr::eval_type_t<E> sum{};
//    std::forward<assign_type>(assign) = sum;
// #ifdef EXECUTION_HEADER_AVAILABLE
//    if (params::parallelization)
//      std::forward<assign_type>(assign) = std::reduce(
//          std::execution::par_unseq,
//          static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
//          static_cast<const E*>(&e)->end(symphas::it_grp, interval));
//    else
// #endif
//      std::forward<assign_type>(assign) = std::reduce(
//          static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
//          static_cast<const E*>(&e)->end(symphas::it_grp, interval));
//  } else {
//    expr::eval_type_t<E> sum{};
//    std::forward<assign_type>(assign) = sum;
//  }
//}
//
// template <typename condition_t, typename... condition_ts, typename E,
//          typename assign_type, size_t D>
// void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign,
//                        grid::region_interval_multiple<D> const& regions) {
//  expr::eval_type_t<E> sum{};
//  std::forward<assign_type>(assign) = sum;
//  for (grid::region_interval<D> region : regions) {
//    result_sum_by_term<condition_t, condition_ts...>(*static_cast<E
//    const*>(&e),
//                                                     sum, region);
//    std::forward<assign_type>(assign) += sum;
//  }
//}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_of_matching(OpEvaluable<E> const& e, assign_type&& data) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 expression_satisfies_condition<E, condition_t>)) {
    result_accumulate(*static_cast<E const*>(&e),
                      std::forward<assign_type>(data));
  }
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, size_t... Is>
void result_of_matching(OpAdd<Es...> const& e, assign_type&& data,
                        std::index_sequence<Is...>) {
  (result_of_matching<condition_t, condition_ts...>(
       expr::get<Is>(e), std::forward<assign_type>(data)),
   ...);
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_of_matching(OpAdd<Es...> const& e, assign_type&& data) {
  result_of_matching<condition_t, condition_ts...>(
      e, std::forward<assign_type>(data),
      std::make_index_sequence<sizeof...(Es)>{});
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, typename region_type>
void result_of_matching(OpEvaluable<E> const& e, assign_type&& data,
                        region_type&& region) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 expression_satisfies_condition<E, condition_t>)) {
    result_accumulate(*static_cast<E const*>(&e),
                      std::forward<assign_type>(data),
                      std::forward<region_type>(region));
  } else {
    expr::eval_type_t<E> sum{};
    std::forward<assign_type>(assign) = sum;
  }
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, typename region_type, size_t... Is>
void result_of_matching(OpAdd<Es...> const& e, assign_type&& data,
                        region_type&& region, std::index_sequence<Is...>) {
  (result_of_matching<condition_t, condition_ts...>(
       expr::get<Is>(e), std::forward<assign_type>(data),
       std::forward<region_type>(region)),
   ...);
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, typename region_type>
void result_of_matching(OpAdd<Es...> const& e, assign_type&& data,
                        region_type&& region) {
  result_of_matching<condition_t, condition_ts...>(
      e, std::forward<assign_type>(data), std::forward<region_type>(region),
      std::make_index_sequence<sizeof...(Es)>{});
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, typename region_type>
void result_sum_of_matching(OpEvaluable<E> const& e, assign_type&& assign,
                            region_type&& region) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 expression_satisfies_condition<E, condition_t>)) {
    result_sum(*static_cast<E const*>(&e), std::forward<assign_type>(assign),
               std::forward<region_type>(region));
  } else {
    std::forward<assign_type>(assign) = expr::eval_type_t<E>{};
  }
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, typename region_type>
void result_sum_of_matching(OpAdd<Es...> const& e, assign_type&& assign,
                            region_type&& region, std::index_sequence<>) {
  std::forward<assign_type>(assign) = expr::eval_type_t<OpAdd<Es...>>{};
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, typename region_type, size_t I0, size_t... Is>
void result_sum_of_matching(OpAdd<Es...> const& e, assign_type&& assign,
                            region_type&& region,
                            std::index_sequence<I0, Is...>) {
  expr::eval_type_t<OpAdd<Es...>> sum1{}, sum2{};
  result_sum_of_matching<condition_t, condition_ts...>(
      expr::get<I0>(e), sum1, std::forward<region_type>(region));
  result_sum_of_matching<condition_t, condition_ts...>(
      e, std::forward<assign_type>(assign), sum2, std::index_sequence<Is...>{});
  std::forward<assign_type>(assign) = sum1 + sum2;
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type, typename region_type>
void result_sum_of_matching(OpAdd<Es...> const& e, assign_type&& assign,
                            region_type&& region) {
  result_sum_of_matching<condition_t, condition_ts...>(
      e, std::forward<assign_type>(assign), std::forward<region_type>(region),
      std::make_index_sequence<sizeof...(Es)>{});
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename assign_type,
          typename E>
void result_sum_of_matching(OpEvaluable<E> const& e, assign_type&& assign) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 expression_satisfies_condition<E, condition_t>)) {
    result_sum(*static_cast<E const*>(&e), std::forward<assign_type>(assign),
               expr::iterable_domain(*static_cast<E const*>(&e)));
  } else {
    std::forward<assign_type>(assign) = expr::eval_type_t<E>{};
  }
}

//! Accumulates the result if the given expression matches the condition.
template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_sum_of_matching(OpAdd<Es...> const& e, assign_type&& assign) {
  result_sum_of_matching<condition_t, condition_ts...>(
      e, std::forward<assign_type>(assign), expr::iterable_domain(e));
}

template <typename condition_t>
struct result_by_term_apply {
  template <typename E, typename assign_type, typename region_type>
  void operator()(OpEvaluable<E> const& e, assign_type&& data,
                  region_type&& region) const {
    result_of_matching<condition_t>(*static_cast<E const*>(&e),
                                    std::forward<assign_type>(data),
                                    std::forward<region_type>(region));
  }

  template <typename E, typename assign_type>
  void operator()(OpEvaluable<E> const& e, assign_type&& data) const {
    result_of_matching<condition_t>(*static_cast<E const*>(&e),
                                    std::forward<assign_type>(data));
  }
};

template <>
struct result_by_term_apply<expr::matches_series> : result_by_term_apply<void> {
  using result_by_term_apply<void>::operator();

  template <typename V0, typename E, typename... Ts, int... I0s, int... P0s,
            typename A, typename B, typename... Vs, typename assign_type,
            typename region_type>
  void operator()(
      OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
            symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
            symphas::lib::types_list<Vs...>> const& sum,
      assign_type&& data, region_type&& region) const {
    result_accumulate(sum, std::forward<assign_type>(data),
                      std::forward<region_type>(region));
  }

  template <typename V0, typename E, typename... Ts, int... I0s, int... P0s,
            typename A, typename B, typename... Vs, typename assign_type>
  void operator()(
      OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
            symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
            symphas::lib::types_list<Vs...>> const& sum,
      assign_type&& data) const {
    TIME_THIS_CONTEXT_LIFETIME(apply_summation);
    auto coeff = expr::coeff(sum);
    for (iter_type i = 0; i < sum.data.persistent.len; ++i) {
      auto region = expr::iterable_domain(sum.data.persistent[i].e);
      result_accumulate(coeff * sum.data.persistent[i].e,
                        std::forward<assign_type>(data), region);
    }
  }
};

template <typename condition_t = void>
struct result_sum_by_term_apply {
  template <typename E, typename assign_type, typename region_type>
  void operator()(OpEvaluable<E> const& e, assign_type&& assign,
                  region_type&& region) const {
    result_sum_of_matching<condition_t>(*static_cast<E const*>(&e),
                                        std::forward<assign_type>(assign),
                                        std::forward<region_type>(region));
  }

  template <typename E, typename assign_type>
  void operator()(OpEvaluable<E> const& e, assign_type&& assign) const {
    result_sum_of_matching<condition_t>(*static_cast<E const*>(&e),
                                        std::forward<assign_type>(assign));
  }
};

template <>
struct result_sum_by_term_apply<expr::matches_series>
    : result_sum_by_term_apply<void> {
  using result_sum_by_term_apply<void>::operator();

  static const len_type group_size = 6;

  template <typename V0, typename E, typename... Ts, int... I0s, int... P0s,
            typename A, typename B, typename... Vs, typename assign_type,
            typename region_type>
  void operator()(
      OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
            symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
            symphas::lib::types_list<Vs...>> const& sum,
      assign_type&& assign, region_type&& region) const {
    result_sum_by_term(sum, std::forward<assign_type>(assign),
                       std::forward<region_type>(region));
  }

  template <typename V0, typename E, typename... Ts, int... I0s, int... P0s,
            typename A, typename B, typename... Vs, typename assign_type>
  void operator()(
      OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
            symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
            symphas::lib::types_list<Vs...>> const& sum,
      assign_type&& assign) const {
    if (sum.data.persistent.len > 0) {
      auto reduced = expr::eval_type_t<
          OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
                symphas::lib::types_list<Vs...>>>{};
      result_sum(expr::coeff(sum) * sum.data.persistent[0].e, reduced);
      std::forward<assign_type>(assign) = reduced;

      auto coeff = expr::coeff(sum);
      for (iter_type i = 1; i < sum.data.persistent.len; ++i) {
        result_sum(coeff * sum.data.persistent[i].e, reduced);
        std::forward<assign_type>(assign) += reduced;
      }
    } else {
      using expr_type =
          OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
                symphas::lib::types_list<Vs...>>;
      std::forward<assign_type>(assign) = expr::eval_type_t<expr_type>{};
    }
  }
};

template <>
struct result_sum_by_term_apply<expr::matching_in_mul<expr::matches_series>>
    : result_sum_by_term_apply<void> {
  using result_sum_by_term_apply<void>::operator();

  template <typename A, typename B, typename assign_type, typename region_type,
            std::enable_if_t<symphas::internal::expression_satisfies_condition<
                                 OpBinaryMul<A, B>,
                                 expr::matching_in_mul<expr::matches_series>>,
                             int> = 0>
  void operator()(OpBinaryMul<A, B> const& e, assign_type&& assign,
                  region_type&& region) const {
    result_sum_by_term(e, std::forward<assign_type>(assign),
                       std::forward<region_type>(region));
  }

  template <typename A, typename B, typename assign_type,
            std::enable_if_t<symphas::internal::expression_satisfies_condition<
                                 OpBinaryMul<A, B>,
                                 expr::matching_in_mul<expr::matches_series>>,
                             int> = 0>
  void operator()(OpBinaryMul<A, B> const& e, assign_type&& assign) const {
    apply_mul(std::forward<assign_type>(assign), e.a, e.b);
  }

  template <typename assign_type, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs,
            typename... Es>
  void apply_mul(assign_type&& assign,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>> const& sum,
                 Es&&... terms) const {
    using expr_type = mul_result_t<
        mul_result_t<Es...>,
        OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
              symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
              symphas::lib::types_list<Vs...>>>;

    if (sum.data.persistent.len > 0) {
      auto e = (expr::coeff(sum) * ... * std::forward<Es>(terms));
      auto region = expr::iterable_domain(e);
      auto reduced = expr::eval_type_t<expr_type>{};
      std::forward<assign_type>(assign) = reduced;

      TIME_THIS_CONTEXT_LIFETIME(reduce_mul_summation);
      for (iter_type i = 0; i < sum.data.persistent.len; ++i) {
        auto region0 = region / expr::iterable_domain(sum.data.persistent[i].e);
        TIME_THIS_EXPRESSION_LIFETIME(
            iterable_domain_reduce,
            auto r = expr::iterable_domain(sum.data.persistent[i].e);)

        if (!grid::is_empty(region0)) {
          result_sum(expr::make_mul(e, sum.data.persistent[i].e), reduced,
                     region0);
          std::forward<assign_type>(assign) += reduced;
        }
      }
    } else {
      std::forward<assign_type>(assign) = expr::eval_type_t<expr_type>{};
    }
  }

  template <typename assign_type, typename E0, typename V0, typename E,
            typename... Ts, int... I0s, int... P0s, typename A, typename B,
            typename... Vs, typename... Es>
  void apply_mul(assign_type&& assign, OpEvaluable<E0> const& e,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>> const& sum,
                 Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), sum,
              *static_cast<E0 const*>(&e), std::forward<Es>(terms)...);
  }

  template <typename assign_type, typename E0, typename V0, typename E,
            typename... Ts, int... I0s, int... P0s, typename A, typename B,
            typename... Vs, typename... Es>
  void apply_mul(assign_type&& assign, OpOperator<E0> const& e,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>> const& sum,
                 Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), sum,
              *static_cast<E0 const*>(&e), std::forward<Es>(terms)...);
  }

  template <typename assign_type, typename A, typename B, typename E,
            typename... Es>
  void apply_mul(assign_type&& assign, OpEvaluable<E> const& e0,
                 OpBinaryMul<A, B> const& e1, Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), e1.a, e1.b,
              *static_cast<E const*>(&e0), std::forward<Es>(terms)...);
  }

  template <typename assign_type, typename A, typename B, typename E,
            typename... Es>
  void apply_mul(assign_type&& assign, OpOperator<E> const& e0,
                 OpBinaryMul<A, B> const& e1, Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), e1.a, e1.b,
              *static_cast<E const*>(&e0), std::forward<Es>(terms)...);
  }

  template <typename assign_type, typename A, typename B, typename E,
            typename... Es>
  void apply_mul(assign_type&& assign, OpBinaryMul<A, B> const& e0,
                 OpEvaluable<E> const& e1, Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), *static_cast<E const*>(&e1),
              e0, std::forward<Es>(terms)...);
  }

  template <typename assign_type, typename A, typename B, typename C,
            typename D, typename... Es>
  void apply_mul(assign_type&& assign, OpBinaryMul<A, B> const& e0,
                 OpBinaryMul<C, D> const& e1, Es&&... terms) const {
    apply_mul(std::forward<assign_type>(assign), e0 * e1.a, e1.b,
              std::forward<Es>(terms)...);
  }
};

template <typename E, typename assign_type>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     symphas::lib::types_list<>) {
  result_accumulate(*static_cast<E const*>(&e),
                    std::forward<assign_type>(data));
}

template <typename E, typename assign_type, typename condition_t,
          typename... condition_ts>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     symphas::lib::types_list<condition_t, condition_ts...>) {
  auto&& [eval, rest] = compile_trait::separate_by_trait<condition_t>{}(
      *static_cast<E const*>(&e));
  result_accumulate(eval, std::forward<assign_type>(data));
  result_by_group(rest, std::forward<assign_type>(data),
                  symphas::lib::types_list<condition_ts...>{});
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, typename E0>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     E0 const& init) {
  result(init, std::forward<assign_type>(data));
  result_by_group(*static_cast<E const*>(&e), std::forward<assign_type>(data),
                  symphas::lib::types_list<condition_t, condition_ts...>{});
}

template <typename E, typename assign_type, typename region_type>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     region_type&& region, symphas::lib::types_list<>) {
  result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data),
                    std::forward<region_type>(region));
}

template <typename E, typename assign_type, typename region_type,
          typename condition_t, typename... condition_ts>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     region_type&& region,
                     symphas::lib::types_list<condition_t, condition_ts...>) {
  auto&& [eval, rest] = compile_trait::separate_by_trait<condition_t>{}(
      *static_cast<E const*>(&e));
  result_accumulate(eval, std::forward<assign_type>(data));
  result_by_group(rest, std::forward<assign_type>(data),
                  std::forward<region_type>(region),
                  symphas::lib::types_list<condition_ts...>{});
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, typename region_type, typename E0>
void result_by_group(OpEvaluable<E> const& e, assign_type&& data,
                     region_type&& region, E0 const& init) {
  result(init, std::forward<assign_type>(data),
         std::forward<region_type>(region));
  result_by_group(*static_cast<E const*>(&e), std::forward<assign_type>(data),
                  std::forward<region_type>(region),
                  symphas::lib::types_list<condition_t, condition_ts...>{});
}

template <typename assign_type>
struct result_only_trait {
  template <typename... Es, size_t D, size_t... Is>
  result_only_trait(OpAdd<Es...> const& e, assign_type&& data,
                    grid::region_interval<D> const& interval,
                    std::index_sequence<Is...>) {
    symphas::data_iterator_group it(std::forward<assign_type>(data), interval);
    auto start = symphas::reduce_iterator(
        expr::get<Is>(e).begin(symphas::it_grp, interval)...);
    auto end = symphas::reduce_iterator(
        expr::get<Is>(e).end(symphas::it_grp, interval)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      std::transform(std::execution::par_unseq, start, end, it,
                     forward_value{});
    else
#endif
      std::transform(start, end, it, forward_value{});
  }

  template <typename... Es, size_t... Is>
  result_only_trait(OpAdd<Es...> const& e, assign_type&& data, len_type len,
                    std::index_sequence<Is...>) {
    symphas::data_iterator it(std::forward<assign_type>(data));
    auto start = symphas::reduce_seq_iterator(expr::get<Is>(e).begin()...);
    auto end = symphas::reduce_seq_iterator(expr::get<Is>(e).end(len)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      std::transform(std::execution::par_unseq, start, end, it,
                     forward_value{});
    else
#endif
      std::transform(start, end, it, forward_value{});
  }
};

template <typename... Es, typename assign_type, size_t D, size_t... Is>
result_only_trait(OpAdd<Es...>, assign_type&&, grid::region_interval<D>,
                  std::index_sequence<Is...>) -> result_only_trait<assign_type>;
template <typename... Es, typename assign_type, size_t... Is>
result_only_trait(OpAdd<Es...>, assign_type&&, len_type,
                  std::index_sequence<Is...>) -> result_only_trait<assign_type>;

//! Evaluate only the terms of the expression matching the condition.
/*!
 * The expression must be iterable over the entire given length.
 *
 * \param e Expression that is evaluated.
 * \param data The array containing the result of the expression.
 * \param len The number of elements in the array.
 */
template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_only(OpEvaluable<E> const& e, assign_type&& data, len_type len) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result(*static_cast<E const*>(&e), std::forward<assign_type>(data), len);
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, size_t D>
void result_only(OpEvaluable<E> const& e, assign_type&& data,
                 grid::region_interval<D> const& interval) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result(*static_cast<E const*>(&e), std::forward<assign_type>(data),
           interval);
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, size_t D>
void result_only(OpEvaluable<E> const& e, assign_type&& data,
                 grid::region_interval_multiple<D> const& regions) {
  for (grid::region_interval<D> region : regions) {
    result_only<condition_t, condition_ts...>(
        *static_cast<E const*>(&e), std::forward<assign_type>(data), region);
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_only(OpEvaluable<E> const& e, assign_type&& data,
                 grid::region_interval<0> const& interval) {
  auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
  if (grid::length(data_region) > 1) {
    result_only<condition_t, condition_ts...>(*static_cast<E const*>(&e),
                                              std::forward<assign_type>(data),
                                              data_region);
  } else {
    result_only<condition_t, condition_ts...>(
        *static_cast<E const*>(&e), std::forward<assign_type>(data), 1);
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_only(OpEvaluable<E> const& e, assign_type&& data,
                 grid::region_empty) {}

template <typename... Es, typename assign_type, size_t... Is>
void result_only(OpAdd<Es...> const& e, assign_type&& data, len_type len,
                 std::index_sequence<Is...>) {
  result_only_trait(e, std::forward<assign_type>(data), len,
                    std::index_sequence<Is...>{});
}

template <typename... Es, typename assign_type, size_t... Is>
void result_only(OpAdd<Es...> const& e, assign_type&& data,
                 grid::region_interval<0> const& interval,
                 std::index_sequence<Is...>) {
  auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
  if (grid::length(data_region) > 1) {
    result_only(e, std::forward<assign_type>(data), data_region,
                std::index_sequence<Is...>{});
  } else {
    result_only(e, std::forward<assign_type>(data), 1,
                std::index_sequence<Is...>{});
  }
}

template <typename... Es, typename assign_type, size_t D, size_t... Is>
void result_only(OpAdd<Es...> const& e, assign_type&& data,
                 grid::region_interval<D> const& interval,
                 std::index_sequence<Is...>) {
  result_only_trait<assign_type>(e, std::forward<assign_type>(data), interval,
                                 std::index_sequence<Is...>{});
}

template <typename... Es, typename assign_type, size_t D, size_t... Is>
void result_only(OpAdd<Es...> const& e, assign_type&& data,
                 grid::region_interval_multiple<D> const& regions,
                 std::index_sequence<Is...>) {
  for (grid::region_interval<D> region : regions) {
    result_only(e, std::forward<assign_type>(data), region,
                std::index_sequence<Is...>{});
  }
}

template <typename... Es, typename assign_type>
void result_only(OpAdd<Es...> const& e, assign_type&& data,
                 std::index_sequence<>) {
  result(OpVoid{}, std::forward<assign_type>(data),
         expr::iterable_domain(std::forward<assign_type>(data)));
}

template <typename... Es, typename assign_type, size_t I0, size_t... Is>
void result_only(OpAdd<Es...> const& e, assign_type&& data,
                 std::index_sequence<I0, Is...>) {
  auto region = (expr::iterable_domain(expr::get<Is>(e)) + ... +
                 expr::iterable_domain(expr::get<I0>(e)));
  result_only(e, std::forward<assign_type>(data), region,
              std::index_sequence<I0, Is...>{});
}

template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_only(OpAdd<Es...> const& e, assign_type&& data) {
  using namespace symphas::internal;
  using mask = std::integer_sequence<
      bool, expression_satisfies_condition<
                Es, expr::or_<condition_t, condition_ts...>>...>;
  using seq = std::make_index_sequence<sizeof...(Es)>;
  result_only(e, std::forward<assign_type>(data),
              symphas::lib::filter_seq_t<seq, mask>{});
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_only(OpEvaluable<E> const& e, assign_type&& data) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result(*static_cast<E const*>(&e), std::forward<assign_type>(data));
  }
}

template <typename... Es, typename assign_type, size_t... Is>
void result_by_term(OpAdd<Es...> const& e, assign_type&& data,
                    symphas::lib::types_list<>, std::index_sequence<Is...>) {}

template <typename... Es, typename assign_type, typename condition_t,
          typename... condition_ts, size_t... Is>
void result_by_term(OpAdd<Es...> const& e, assign_type&& data,
                    symphas::lib::types_list<condition_t, condition_ts...>,
                    std::index_sequence<Is...>) {
  (result_by_term_apply<condition_t>{}(expr::get<Is>(e),
                                       std::forward<assign_type>(data)),
   ...);
  result_by_term(e, std::forward<assign_type>(data),
                 symphas::lib::types_list<condition_ts...>{},
                 std::index_sequence<Is...>{});
}

template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_by_term(OpAdd<Es...> const& e, assign_type&& data) {
  if constexpr ((expr::satisfies<Es, expr::or_<condition_t, condition_ts...>> ||
                 ...)) {
    result_only<expr::not_<expr::or_<condition_t, condition_ts...>>>(
        e, std::forward<assign_type>(data));
    result_by_term(e, std::forward<assign_type>(data),
                   symphas::lib::types_list<condition_t, condition_ts...>{},
                   std::make_index_sequence<sizeof...(Es)>{});
  } else {
    result(e, std::forward<assign_type>(data));
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_by_term(OpEvaluable<E> const& e, assign_type&& data) {
  result(*static_cast<E const*>(&e), std::forward<assign_type>(data));
}

template <typename assign_type>
struct result_sum_only_trait {
  template <typename... Es, size_t... Is>
  void result(OpAdd<Es...> const& e, assign_type& assign, len_type len,
              std::index_sequence<Is...>) {
    auto start = symphas::reduce_seq_iterator(expr::get<Is>(e).begin()...);
    auto end = symphas::reduce_seq_iterator(expr::get<Is>(e).end(len)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      assign = std::reduce(std::execution::par_unseq, start, end);
    else
#endif
      assign = std::reduce(start, end);
  }

  template <typename... Es, size_t D, size_t... Is>
  void result(OpAdd<Es...> const& e, assign_type& assign,
              grid::region_interval<D> const& interval,
              std::index_sequence<Is...>) {
    auto start = symphas::reduce_iterator(
        expr::get<Is>(e).begin(symphas::it_grp, interval)...);
    auto end = symphas::reduce_iterator(
        expr::get<Is>(e).end(symphas::it_grp, interval)...);

#ifdef EXECUTION_HEADER_AVAILABLE
    if (params::parallelization)
      assign = std::reduce(std::execution::par_unseq, start, end);
    else
#endif
      assign = std::reduce(start, end);
  }
};

//! Evaluate only the terms of the expression matching the condition.
/*!
 * The expression must be iterable over the entire given length.
 *
 * \param e Expression that is evaluated.
 * \param data The array containing the result of the expression.
 * \param len The number of elements in the array.
 */
template <typename condition_t, typename... condition_ts, typename assign_type,
          typename E>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign,
                     len_type len) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result_sum_trait<expr::storage_type_t<E>>{}.result(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign), len);
  } else {
    std::forward<assign_type>(assign) = expr::eval_type_t<E>{};
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, size_t D>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign,
                     grid::region_interval<D> const& interval) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result_sum_trait<expr::storage_type_t<E>>{}.result(
        *static_cast<E const*>(&e), std::forward<assign_type>(assign),
        interval);
  } else {
    std::forward<assign_type>(assign) = expr::eval_type_t<E>{};
  }
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, size_t D>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign,
                     grid::region_interval_multiple<D> const& regions) {
  expr::eval_type_t<E> sum{};
  std::forward<assign_type>(assign) = sum;
  for (grid::region_interval<D> region : regions) {
    result_sum_only<condition_t, condition_ts...>(*static_cast<E const*>(&e),
                                                  result, sum);
    std::forward<assign_type>(assign) += sum;
  }
}

template <typename condition_t, typename... condition_ts, typename assign_type,
          typename E>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign,
                     grid::region_interval<0> const& interval) {
  result_sum_only(*static_cast<E const*>(&e), std::forward<assign_type>(assign),
                  1);
}

template <typename condition_t, typename... condition_ts, typename assign_type,
          typename E>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign,
                     grid::region_empty) {
  std::forward<assign_type>(assign) = expr::eval_type_t<E>{};
}

template <typename... Es, typename assign_type, size_t... Is>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign, len_type len,
                     std::index_sequence<Is...>) {
  result_sum_only_trait<expr::storage_type_t<OpAdd<Es...>>>{}.result(
      e, std::forward<assign_type>(assign), len, std::index_sequence<Is...>{});
}

template <typename... Es, typename assign_type, size_t... Is>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign,
                     grid::region_interval<0> const& interval,
                     std::index_sequence<Is...>) {
  return result_sum_only(e, std::forward<assign_type>(assign), 1,
                         std::index_sequence<Is...>{});
}

template <typename... Es, typename assign_type, size_t D, size_t... Is>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign,
                     grid::region_interval<D> const& interval,
                     std::index_sequence<Is...>) {
  return result_sum_only_trait<expr::storage_type_t<OpAdd<Es...>>>{}.result(
      e, std::forward<assign_type>(assign), interval,
      std::index_sequence<Is...>{});
}

template <typename... Es, typename assign_type, size_t D, size_t... Is>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign,
                     grid::region_interval_multiple<D> const& regions,
                     std::index_sequence<Is...>) {
  expr::eval_type_t<OpAdd<symphas::lib::direct_type_at_index<Is, Es...>...>>
      result{};
  std::forward<assign_type>(assign) = result;
  for (grid::region_interval<D> region : regions) {
    result_sum_only(e, result, region, std::index_sequence<Is...>{});
    std::forward<assign_type>(assign) += result;
  }
}

template <typename... Es, typename assign_type>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign,
                     std::index_sequence<>) {
  std::forward<assign_type>(assign) = expr::eval_type_t<OpAdd<Es...>>{};
}

template <typename... Es, typename assign_type, size_t I0, size_t... Is>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign,
                     std::index_sequence<I0, Is...>) {
  auto region = (expr::iterable_domain(expr::get<Is>(e)) + ... +
                 expr::iterable_domain(expr::get<I0>(e)));
  result_sum_only(e, std::forward<assign_type>(assign), region,
                  std::index_sequence<I0, Is...>{});
}

template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_sum_only(OpAdd<Es...> const& e, assign_type&& assign) {
  using namespace symphas::internal;
  using mask = std::integer_sequence<
      bool, expression_satisfies_condition<
                Es, expr::or_<condition_t, condition_ts...>>...>;
  using seq = std::make_index_sequence<sizeof...(Es)>;
  return result_sum_only(e, std::forward<assign_type>(assign),
                         symphas::lib::filter_seq_t<seq, mask>{});
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_sum_only(OpEvaluable<E> const& e, assign_type&& assign) {
  if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>) {
    result_sum(*static_cast<E const*>(&e), std::forward<assign_type>(assign));
  } else {
    std::forward<assign_type>(assign) = OpVoid{};
  }
}

template <typename E, typename assign_type>
void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign,
                        symphas::lib::types_list<>) {
  expr::eval_type_t<E> sum{};
  std::forward<assign_type>(assign) = sum;
}

template <typename E, typename assign_type, typename condition_t,
          typename... condition_ts>
void result_sum_by_term(
    OpEvaluable<E> const& e, assign_type&& assign,
    symphas::lib::types_list<condition_t, condition_ts...>) {
  expr::eval_type_t<E> sum{};
  std::forward<assign_type>(assign) = sum;
  result_sum_by_term_apply<condition_t>{}(*static_cast<E const*>(&e), sum);

  result_sum_by_term(*static_cast<E const*>(&e), sum,
                     symphas::lib::types_list<condition_ts...>{});
  std::forward<assign_type>(assign) += sum;
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type>
void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign) {
  result_sum_by_term(*static_cast<E const*>(&e),
                     std::forward<assign_type>(assign),
                     symphas::lib::types_list<condition_t, condition_ts...>{});
}

template <typename condition_t, typename... condition_ts, typename E,
          typename assign_type, typename interval_type>
void result_sum_by_term(OpEvaluable<E> const& e, assign_type&& assign,
                        interval_type&& interval) {
  result_sum<condition_t, condition_ts...>(
      *static_cast<E const*>(&e), std::forward<assign_type>(assign),
      std::forward<interval_type>(interval));
}

template <typename... Es, typename assign_type, typename condition_t,
          typename... condition_ts>
void result_sum_by_term(OpAdd<Es...> const& e, assign_type&& assign,
                        symphas::lib::types_list<condition_t, condition_ts...>,
                        std::index_sequence<>) {
  expr::eval_type_t<OpAdd<Es...>> sum{};
  std::forward<assign_type>(assign) = sum;
}

template <typename... Es, typename assign_type, typename condition_t,
          typename... condition_ts, size_t I0, size_t... Is>
void result_sum_by_term(OpAdd<Es...> const& e, assign_type&& assign,
                        symphas::lib::types_list<condition_t, condition_ts...>,
                        std::index_sequence<I0, Is...>) {
  expr::eval_type_t<OpAdd<Es...>> sum{};
  result_sum_by_term(expr::get<I0>(e), sum,
                     symphas::lib::types_list<condition_t, condition_ts...>{});
  std::forward<assign_type>(assign) = sum;

  result_sum_by_term(e, sum,
                     symphas::lib::types_list<condition_t, condition_ts...>{},
                     std::index_sequence<Is...>{});
  std::forward<assign_type>(assign) += sum;
}

template <typename condition_t, typename... condition_ts, typename... Es,
          typename assign_type>
void result_sum_by_term(OpAdd<Es...> const& e, assign_type&& assign) {
  expr::eval_type_t<OpAdd<Es...>> sum{};

  result_sum_only<expr::not_<expr::or_<condition_t, condition_ts...>>>(e, sum);
  std::forward<assign_type>(assign) = sum;

  result_sum_by_term(e, sum,
                     symphas::lib::types_list<condition_t, condition_ts...>{},
                     std::make_index_sequence<sizeof...(Es)>{});
  std::forward<assign_type>(assign) += sum;
}

}  // namespace expr
