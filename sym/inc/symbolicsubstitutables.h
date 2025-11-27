
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
#include "symbolicdata.h"

/*!
 * \addtogroup substitutables Symbolic Constructs and Expression Building
 * Symbolic constructs are used to simplify expression definitions and easily
 * substitute data.
 *
 * The symbolic algebra library defines two types of symbolic algebra groups:
 * Classes with a name beginning in **Op**, which indicates that they are a
 * symbolic algebra expression which can be _evaluated_; and, Classes which
 * begin with **Symbolic**, which are constructs that extend the
 * functionality/usability of symbolic expressions and cannot be evaluated.
 *
 * The **Symbolic** constructs include an evaluable function (SymbolicFunction),
 * a pattern or template type that allows substitution of expressions
 * (SymbolicTemplate), and a series object (SymbolicSeries).
 * @{
 */

namespace symphas::internal {

template <typename T0>
auto _make_symbolic_data(T0&& arg) {
  return SymbolicData<T0>(arg);
}

template <typename T0>
auto _make_symbolic_data(T0& arg) {
  return SymbolicData<T0>(arg);
}

template <typename T0>
auto _make_symbolic_data(SymbolicData<T0>&& arg) {
  return arg;
}

template <typename T0>
auto _make_symbolic_data(SymbolicData<T0>& arg) {
  return arg;
}

template <typename T0>
auto _make_symbolic_data(SymbolicDataArray<T0>&& arg) {
  return arg;
}

template <typename T0>
auto _make_symbolic_data(SymbolicDataArray<T0>& arg) {
  return arg;
}

template <typename T0>
auto make_symbolic_data(T0&& arg) {
  return _make_symbolic_data(arg);
}

template <typename T0>
auto make_symbolic_data(T0& arg) {
  return _make_symbolic_data(arg);
}

template <typename T0>
auto make_symbolic_data(T0 const& arg) {
  return _make_symbolic_data(const_cast<T0&>(arg));
}

template <typename T0>
auto as_arg(T0&& arg) {
  if constexpr (expr::is_expression<T0>) {
    return std::forward<T0>(arg);
  } else {
    return make_symbolic_data(std::forward<T0>(arg));
  }
}

template <typename T0>
using const_arr_type = const T0*;
template <typename T0>
using arr_type = T0*;

template <typename T0>
auto _make_symbolic_data(const_arr_type<T0>&& arg, len_type len) {
  return SymbolicDataArray<T0>(const_cast<T0*>(arg), len, true);
}

template <typename T0>
auto _make_symbolic_data(const_arr_type<T0>& arg, len_type len) {
  return SymbolicDataArray<T0>(const_cast<T0*&>(arg), len, false);
}

/*template<typename T0>
auto _make_symbolic_data(const_arr_type<T0> const& arg, len_type len)
{
        return SymbolicDataArray<T0>(const_cast<T0*&>(arg), len, false);
}*/

template <typename T0>
auto _make_symbolic_data(NamedData<T0*> arg, len_type len) {
  return SymbolicDataArray<NamedData<T0*>>(arg, len, false);
}

template <typename T0>
auto _make_symbolic_data(arr_type<T0>&& arg, len_type len) {
  return SymbolicDataArray<T0>(arg, len, true);
}

template <typename T0>
auto _make_symbolic_data(arr_type<T0>& arg, len_type len) {
  return SymbolicDataArray<T0>(arg, len, false);
}

// template<typename T0>
// auto _make_symbolic_data(arr_type<T0> const& arg, len_type len)
//{
//	return SymbolicDataArray<T0>(const_cast<arr_type<T0>*&>(arg), len,
// false);
// }

template <typename... T0s>
auto _make_symbolic_data(std::tuple<T0s...>&& arg, len_type len) {
  return SymbolicDataArray(arg);
}

template <typename... T0s>
auto _make_symbolic_data(std::tuple<T0s...>& arg, len_type len) {
  return SymbolicDataArray(arg);
}

template <typename T0>
auto _make_symbolic_data(SymbolicDataArray<T0>&& arg, len_type len) {
  return arg;
}

template <typename T0>
auto _make_symbolic_data(SymbolicDataArray<T0>& arg, len_type len) {
  return arg;
}

template <typename... Vs, typename... T0s, size_t... Is>
auto _make_symbolic_data(std::tuple<OpTerm<Vs, T0s>...>& arg,
                         std::index_sequence<Is...>, len_type len) {
  return _make_symbolic_data(
      std::make_tuple(expr::get<1>(std::get<Is>(arg))...), len);
}

template <typename... Vs, typename... T0s>
auto _make_symbolic_data(std::tuple<OpTerm<Vs, T0s>...>& arg, len_type len) {
  return _make_symbolic_data(arg, std::make_index_sequence<sizeof...(Vs)>{},
                             len);
}

template <typename... Vs, typename... T0s, size_t... Is>
auto _make_symbolic_data(std::tuple<OpTerm<Vs, T0s>...>&& arg,
                         std::index_sequence<Is...>, len_type len) {
  return _make_symbolic_data(
      std::make_tuple(expr::get<1>(std::get<Is>(arg))...), len);
}

template <typename... Vs, typename... T0s>
auto _make_symbolic_data(std::tuple<OpTerm<Vs, T0s>...>&& arg, len_type len) {
  return _make_symbolic_data(arg, std::make_index_sequence<sizeof...(Vs)>{},
                             len);
}

template <typename T0>
auto make_symbolic_data(T0&& arg, len_type len) {
  return _make_symbolic_data(arg, len);
}

template <typename T0, size_t D>
auto make_symbolic_data(T0 (&arg)[D], len_type len) {
  return SymbolicDataArray<T0>(&(arg[0]), len, false);
}

template <typename T0, size_t D>
auto make_symbolic_data(const T0 (&arg)[D], len_type len) {
  return SymbolicDataArray<T0>(const_cast<T0*>(&(arg[0])), len, false);
}

template <typename G>
auto make_symbolic_data(DynamicVariable<G> const& arg, len_type len) {
  return _make_symbolic_data(arg.data, len);
}

template <typename G>
auto make_symbolic_data(OpTerm<OpIdentity, DynamicVariable<G>> const& arg,
                        len_type len) {
  const auto& [v, t] = arg;
  return make_symbolic_data(t, len);
}

template <typename T0>
auto make_symbolic_data(T0& arg, len_type len) {
  return _make_symbolic_data(arg, len);
}

template <typename T0>
auto make_symbolic_data(T0 const& arg, len_type len) {
  return make_symbolic_data(const_cast<T0&>(arg), len);
}

template <typename T0>
auto as_array_arg(T0&& arg, len_type len) {
  return make_symbolic_data(std::forward<T0>(arg), len);
}

template <typename... T0s>
auto make_substitution(T0s... ts) {
  return Substitution<T0s...>(ts...);
}

}  // namespace symphas::internal

namespace expr::symbols {

//! A variable of argument type that can be directly used to define a symbolic
//! construct.
template <size_t N, typename T = expr::symbols::Symbol>
inline arg_t<N, typename trait_arg_type<T>::type> arg;

inline arg_t<1, expr::symbols::Symbol> arg_1;
inline arg_t<2, expr::symbols::Symbol> arg_2;
inline arg_t<3, expr::symbols::Symbol> arg_3;
inline arg_t<4, expr::symbols::Symbol> arg_4;

template <typename T>
inline arg_t<1, typename trait_arg_type<T>::type> arg_1_;
template <typename T>
inline arg_t<2, typename trait_arg_type<T>::type> arg_2_;
template <typename T>
inline arg_t<3, typename trait_arg_type<T>::type> arg_3_;
template <typename T>
inline arg_t<4, typename trait_arg_type<T>::type> arg_4_;

}  // namespace expr::symbols

//! A template into which other expressions can be substituted into.
/*!
 * A symbolic template for expressions to be substituted into.
 * The template is created with a specific number of substitutable parameters.
 * If there are no substitutable parameters, then it is equivalent to an
 * expression.
 *
 * Symbolic templates are defined by using the function expr::template_of() in
 * the following way:
 *
 *     auto tmpl = (expr::template_of() = ...);
 *
 * The elipses, `...`, represent an expression containing #expr::arg terms that
 * act as placeholders, and are where the substitutions will take place. Also,
 * arguments of type Variable<N> or of OpTerm type with Variable<N> data can be
 * passed to set the specific arguments and the order they will be substituted.
 *
 * Subsequently, by calling the initialized `tmpl` variable with OpExpression
 * type arguments, these will be substituted in the pattern that was defined,
 * i.e. `tmpl(x, y)` will generate a new expression with `x` and `y` substituted
 * into the placeholders arg<0> and arg<1>.
 *
 * \tparam E The type of the expression function.
 * \tparam ArgsNs... The indices of the independent variables, v<N>.
 */
template <typename E, size_t... ArgNs>
struct SymbolicTemplate {
  SymbolicTemplate(E const& e) : e{e} {}

 public:
  E operator()() const { return e; }

  template <typename... Xs, typename std::enable_if_t<
                                (sizeof...(Xs) == sizeof...(ArgNs)), int> = 0>
  auto operator()(Xs const&... xs) const {
    return expr::transform::swap_grid<
        Variable<ArgNs, symphas::internal::exclusive_swap>...>(e, xs...);
  }

  template <typename... T0s>
  auto substitute(Substitution<T0s...> const&) const;

  //! Allows this template to be converted into an expression.
  /*!
   * Allows this template to be converted into an expression if there are no
   * substitutable
   */
  template <size_t N = sizeof...(ArgNs),
            typename = std::enable_if_t<(N == 0), int>>
  operator E() const {
    return e;
  }

  E e;  //!< The substitutable function.
};

namespace symphas::internal {

template <typename T>
auto _construct_arg(SymbolicData<T> arg) {
  return arg;
}

template <typename T>
auto _construct_arg(NamedData<SymbolicData<T>> arg) {
  return arg;
}

template <typename T>
auto _construct_arg(SymbolicDataArray<T> arg) {
  return arg;
}

template <typename T>
auto _construct_arg(T arg) {
  if constexpr (expr::is_expression<T> || expr::is_symbol<T>) {
    return arg;
  } else {
    return SymbolicData<T>(&arg, true);
  }
}

template <typename T>
auto construct_arg(T arg) {
  return _construct_arg(arg);
}

template <size_t N, typename T>
auto substitute_arg(T const& data0) {
  if constexpr (expr::is_symbol<T>) {
    return expr::make_term<N>(T{});
  } else {
    return expr::make_term<N>(std::ref(const_cast<T&>(data0)));
  }
}

template <typename E0, typename... Ts, size_t... ArgNs, size_t... Is>
auto substitute_args(OpExpression<E0> const& e0, std::tuple<Ts...> const& data,
                     std::index_sequence<ArgNs...>,
                     std::index_sequence<Is...>) {
  if constexpr (sizeof...(Ts) > 0) {
    return expr::transform::swap_grid<Variable<ArgNs>...>(
        *static_cast<E0 const*>(&e0),
        substitute_arg<ArgNs>(std::get<Is>(data))...);
  } else {
    return *static_cast<E0 const*>(&e0);
  }
}

template <typename E0, size_t... ArgNs, typename... Ts, size_t... Is>
auto substitute_args(SymbolicTemplate<E0, ArgNs...> const& tmpl,
                     std::tuple<Ts...> const& data,
                     std::index_sequence<Is...>) {
  return expr::apply_operators(
      tmpl(substitute_arg<ArgNs>(std::get<Is>(data))...));
}

template <typename E0, size_t... ArgNs, typename... Ts, size_t... Is>
auto substitute_args(SymbolicFunction<E0, Variable<ArgNs, Ts>...> const& func,
                     std::tuple<Ts...> const& data,
                     std::index_sequence<Is...>) {
  return substitute_args(func.e, data, std::index_sequence<ArgNs...>{},
                         std::index_sequence<Is...>{});
}

template <typename E0, size_t... ArgNs, typename... Ts>
auto construct_function(SymbolicTemplate<E0, ArgNs...> const& tmpl,
                        std::tuple<Ts...> const& data) {
  return substitute_args(tmpl, data,
                         std::make_index_sequence<sizeof...(ArgNs)>{});
}

template <typename T, typename R>
auto construct_function_resolve(T const& tmpl, R const& data) {
  return construct_function(tmpl, data);
}

}  // namespace symphas::internal

//! A function into which data can be substituted and be automatically
//! evaluated.
/*!
 * A symbolic function for values to be substituted and be automatically
 * evaluated. The evaluation automatically determines the evaluation type and
 * size of the evaluation grid or block.
 *
 * Symbolic templates are defined by using the function expr::template_of() in
 * the following way:
 *
 *     auto func = (expr::function_of() = ...);
 *
 * The elipses, `...`, represent an expression containing #expr::arg terms that
 * act as placeholders, and are where the substitutions of the data will take
 * place. Also, arguments of type Variable<N, T> or of OpTerm type containing
 * data can be passed to set the specific arguments and the order they will be
 * substituted. The arguments must have some type associated with them so that
 * the function can perform automatic simplifications/symbolic rules based on
 * the data type. This can easily be achieved by using the expr::arg variables
 * with the second template parameter defining the type, i.e.:
 *
 *     auto u = expr::arg<1, double>;
 *     auto func = (expr::function_of(u) = u * u);
 *
 * Subsequently, by calling the initialized `func` variable with any type of
 * data will be substituted in the pattern that was defined, i.e. `func(x, y)`
 * will generate a new expression with `x` and `y` data substituted into the
 * placeholders arg<0> and arg<1>.
 *
 * \tparam E The type of the expression function.
 * \tparam Ts... The types of the variables of the function.
 */
template <typename E, size_t... ArgNs, typename... Ts>
struct SymbolicFunction<E, Variable<ArgNs, Ts>...> {
  using this_type = SymbolicFunction<E, Variable<ArgNs, Ts>...>;

  SymbolicFunction()
      : data{},
        e{symphas::internal::substitute_args(
            E{}, data, std::index_sequence<ArgNs...>{},
            std::make_index_sequence<sizeof...(Ts)>{})} {}

  template <typename T>
  using arg_type = typename std::invoke_result_t<
      decltype(&symphas::internal::construct_arg<T>), T>;
  using data_type = std::tuple<arg_type<Ts>...>;

  //! The linked substitutable data.
  data_type data;

  using expr_type = typename std::invoke_result_t<
      decltype(&symphas::internal::construct_function_resolve<
               SymbolicTemplate<E, ArgNs...>, data_type>),
      SymbolicTemplate<E, ArgNs...>, data_type>;
  using type = SymbolicFunction<expr_type, Variable<ArgNs, Ts>...>;

  using result_type = expr::storage_type_t<E>;

  SymbolicFunction(OpExpression<E> const& e)
      : data{}, e{*static_cast<E const*>(&e)} {}

  template <typename E0>
  SymbolicFunction(SymbolicTemplate<E0, ArgNs...> const& tmpl)
      : data{},
        e{symphas::internal::substitute_args(
            tmpl, data, std::make_index_sequence<sizeof...(Ts)>{})} {}

  SymbolicFunction(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& other)
      : data{},
        e{symphas::internal::substitute_args(
            other, data, std::make_index_sequence<sizeof...(Ts)>{})} {
    set_data_tuple(other.data);
  }

  SymbolicFunction(SymbolicFunction<E, Variable<ArgNs, Ts>...>&& other)
      : data{},
        e{symphas::internal::substitute_args(
            other, data, std::make_index_sequence<sizeof...(Ts)>{})} {
    set_data_tuple(other.data);
    // swap(*this, other);
  }

  SymbolicFunction<E, Variable<ArgNs, Ts>...> operator=(
      SymbolicFunction<E, Variable<ArgNs, Ts>...> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(this_type& first, this_type& second) {
    using std::swap;
    // swap(static_cast<parent_type&>(first),
    // static_cast<parent_type&>(second));
    swap(first.data, second.data);
    auto first_expr = first.e;
    first.e = symphas::internal::substitute_args(
        second.e, first.data, std::index_sequence<ArgNs...>{},
        std::make_index_sequence<sizeof...(Ts)>{});
    second.e = symphas::internal::substitute_args(
        first_expr, second.data, std::index_sequence<ArgNs...>{},
        std::make_index_sequence<sizeof...(Ts)>{});
  }

 protected:
  template <typename T, typename T0>
  void set_data_1(SymbolicData<T0>& data0, T const& arg) {
    data0.set_data(arg);
  }

  template <typename T, typename T0>
  void set_data_1(T0& data0, T const& arg) {
    data0 = arg;
  }

  template <typename... T0s, size_t... Is>
  void set_data(std::index_sequence<Is...>, T0s&&... args) {
    (set_data_1<T0s>(std::get<Is>(data), std::forward<T0s>(args)), ...);
    // expr::prune::update(e);
  }

  template <typename... T0s, size_t... Is>
  result_type evaluate(std::tuple<T0s...> const& args,
                       std::index_sequence<Is...>) {
    set_data(std::index_sequence<Is...>{}, std::get<Is>(args)...);
    result_type result(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));

    expr::result(e, result);
    return result;
  }

  template <typename... T0s, size_t... Is>
  result_type evaluate(Substitution<SymbolicDataArray<T0s>...> const& args,
                       std::index_sequence<Is...>, iter_type i) {
    set_data(std::index_sequence<Is...>{}, std::get<Is>(args).data[i]...);
    result_type result(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));

    expr::result(e, result);
    return result;
  }

  template <typename... T0s, size_t... Is>
  SymbolicDataArray<result_type> evaluate(
      Substitution<SymbolicDataArray<T0s>...> const& args,
      std::index_sequence<Is...>) {
    SymbolicDataArray<result_type> result_arr(std::get<0>(args).len);

    for (iter_type i = 0; i < result_arr.len; ++i) {
      result_arr.data[i] = evaluate(args, std::index_sequence<Is...>{}, i);
    }
    return result_arr;
  }

  template <typename... T0s, size_t... Is>
  auto swap_vars(std::tuple<T0s...> const& args,
                 std::index_sequence<Is...>) const {
    return expr::transform::swap_grid<ArgNs...>(e, std::get<Is>(args)...);
  }

 public:
  template <typename... T0s>
  void set_data(T0s&&... args) {
    set_data(std::make_index_sequence<sizeof...(T0s)>{},
             std::forward<T0s>(args)...);
  }

  template <typename... T0s, size_t... Is>
  void set_data_tuple(std::index_sequence<Is...>,
                      std::tuple<T0s...> const& args) {
    set_data(std::index_sequence<Is...>{}, std::get<Is>(args)...);
  }

  template <typename... T0s>
  void set_data_tuple(std::tuple<T0s...> const& args) {
    set_data_tuple(std::make_index_sequence<sizeof...(T0s)>{}, args);
  }

  template <size_t N = sizeof...(Ts), std::enable_if_t<(N > 0), int> = 0>
  result_type operator()(Ts const&... args) {
    set_data(std::make_index_sequence<sizeof...(Ts)>{}, std::ref(args)...);
    result_type result(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));

    expr::result(e, result);
    return result;
  }

  template <size_t N = sizeof...(Ts), std::enable_if_t<(N > 0), int> = 0>
  result_type operator()() const {
    result_type result(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));

    expr::result(e, result);
    return result;
  }

  template <size_t N = sizeof...(Ts), std::enable_if_t<(N == 0), int> = 0>
  result_type operator()() const {
    result_type result(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));
    expr::result(e, result);
    return result;
  }

  template <typename... T0s, typename std::enable_if_t<
                                 (sizeof...(T0s) == sizeof...(Ts)), int> = 0>
  auto operator()(Substitution<T0s...> const& args) {
    return evaluate(args, std::make_index_sequence<sizeof...(T0s)>{});
  }

  template <typename... T0s, typename std::enable_if_t<
                                 (sizeof...(T0s) == sizeof...(Ts)), int> = 0>
  auto operator()(std::tuple<T0s...> const& args) {
    return evaluate(args, std::make_index_sequence<sizeof...(T0s)>{});
  }

  template <
      typename... T0s,
      typename std::enable_if_t<
          (sizeof...(T0s) == sizeof...(Ts) && sizeof...(T0s) > 0), int> = 0>
  auto operator[](std::tuple<T0s...> const& args) const {
    return swap_vars(args, std::make_index_sequence<sizeof...(T0s)>{});
  }

  template <
      typename... T0s,
      typename std::enable_if_t<
          (sizeof...(T0s) == sizeof...(Ts) && sizeof...(T0s) == 0), int> = 0>
  E const& operator[](std::tuple<T0s...> const& args) const {
    return e;
  }

  // Evaluate the function using the last used arguments.
  result_type result() const {
    result_type r(
        expr::construct_result_data<result_type>{}(expr::data_list(e)));
    expr::result(e, r);
    return r;
  }

  // Evaluate the function at a particular index using the last used arguments.
  auto operator[](iter_type n) const { return e.eval(n); }

  // Evaluate the function at a particular index using the last used arguments.
  auto eval(iter_type n) const { return this->operator[](n); }

  template <typename... T0s>
  auto substitute(Substitution<T0s...> const&) const;

  E e;  //!< The substitutable function.
};

template <typename E, typename... Ts>
struct SymbolicFunction<E, symphas::lib::types_list<Ts...>>
    : SymbolicFunction<E, Ts...> {
  using parent_type = SymbolicFunction<E, Ts...>;
  using parent_type::parent_type;
  using this_type = SymbolicFunction<E, symphas::lib::types_list<Ts...>>;

  friend void swap(this_type& first, this_type& second) {}

  SymbolicFunction() : parent_type() {}
  SymbolicFunction(parent_type const& parent) : parent_type(parent) {}
  SymbolicFunction(parent_type&& parent) : parent_type(parent) {}
};

namespace symphas::internal {
template <size_t N, size_t D, typename E, typename... Ts>
auto swapped_coeffs(OpExpression<E> const& e, DynamicIndex (&index)[D],
                    symphas::lib::types_list<>) {
  return *static_cast<E const*>(&e);
}

template <size_t N, typename E, size_t D, int I0, int P0, int... I0s,
          int... P0s>
auto swapped_coeffs(OpExpression<E> const& e, DynamicIndex (&index)[D],
                    symphas::lib::types_list<expr::symbols::i_<I0, P0>,
                                             expr::symbols::i_<I0s, P0s>...>) {
  return swapped_coeffs<N + 1>(
      expr::transform::swap_grid<OpCoeffSwap<expr::symbols::i_<I0, P0>>>(
          *static_cast<E const*>(&e), index[N]),
      index, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
}

}  // namespace symphas::internal

template <typename E0, typename... T0s>
struct SymbolicFunctionArray {
  using f_type = SymbolicFunction<E0, T0s...>;
  using eval_type =
      std::invoke_result_t<decltype(&f_type::eval), f_type, iter_type>;

  SymbolicFunctionArray(len_type n = 0, len_type len = 0)
      : data{nullptr}, offsets_data{nullptr}, len{len}, n{n} {}

  template <typename E, size_t... Ns, typename... Ts, int... I0s, int... P0s>
  void init_one(iter_type i, SymbolicFunction<E, Variable<Ns, Ts>...> const& f,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>) {
    DynamicIndex index[sizeof...(I0s)]{};
    for (iter_type j = 0; j < sizeof...(I0s); ++j) {
      index[j] = DynamicIndex(offsets_data[i * n + j]);
    }

    auto ff =
        (expr::function_of(Variable<Ns, Ts>{}...) =
             symphas::internal::swapped_coeffs<0>(
                 f.e, index,
                 symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{}));
    ff.set_data_tuple(f.data);
    swap(data[i], ff);
  }

  template <typename E, size_t... Ns, typename... Ts, int... I0s, int... P0s>
  void reallocate(SymbolicFunction<E, Variable<Ns, Ts>...> const& f,
                  len_type len,
                  symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>) {
    this->len = len;
    delete[] data;
    delete[] offsets_data;

    data = new f_type[len]{};
    offsets_data = new iter_type[n * len]{};
    for (iter_type i = 0; i < len; ++i) {
      init_one(i, f,
               symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{});
    }
  }

  SymbolicFunctionArray(SymbolicFunctionArray const& other)
      : data{(other.data != nullptr && other.len > 0) ? new f_type[other.len]{}
                                                      : nullptr},
        offsets_data{(other.offsets_data != nullptr && other.len * other.n > 0)
                         ? new iter_type[other.len * other.n]{}
                         : nullptr},
        len{other.len},
        n{other.n} {
    if (other.offsets_data != nullptr) {
      std::copy(other.offsets_data, other.offsets_data + len * n, offsets_data);
      std::copy(other.data, other.data + len, data);
      for (iter_type i = 0; i < len; ++i) {
        for (iter_type n0 = 0; n0 < n; ++n0) {
          expr::fix_index(data[i].e,
                          DynamicIndex(other.offsets_data[i * n + n0]) ==
                              (offsets_data + i * n + n0));
        }
      }
    }
  }

  SymbolicFunctionArray(SymbolicFunctionArray&& other)
      : SymbolicFunctionArray() {
    swap(*this, other);
  }

  SymbolicFunctionArray<E0, T0s...> operator=(
      SymbolicFunctionArray<E0, T0s...> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(SymbolicFunctionArray<E0, T0s...>& first,
                   SymbolicFunctionArray<E0, T0s...>& second) {
    std::swap(first.data, second.data);
    std::swap(first.offsets_data, second.offsets_data);
    std::swap(first.len, second.len);
    std::swap(first.n, second.n);
  }

  const auto& operator[](iter_type i) const { return data[i]; }

  template <size_t N>
  const auto& operator()(DynamicIndex (&index)[N], iter_type i) const {
    for (iter_type j = 0; j < n; ++j) {
      index[j] = offsets_data[i * n + j];
    }
    return data[i];
  }

  ~SymbolicFunctionArray() {
    if (data != nullptr) {
      delete[] data;
      delete[] offsets_data;
    }
  }

  f_type* data;
  iter_type* offsets_data;
  len_type len;
  len_type n;
};

//! The purpose of this persistent object is storing only offsets related to a
//! series
/*!
 * The idea is that a series boils down to being an object that represents the
 * sum of a list of expressions which only differ by the indices of summation.
 * By storing only those indicies, a.k.a offsets, the process of evaluating and
 * summing the list should be relatively quick compared to other methods.
 */
template <typename E0, typename... T0s>
struct SymbolicFunctionOffsets {
  using f_type = SymbolicFunction<E0, T0s...>;
  using eval_type =
      std::invoke_result_t<decltype(&f_type::eval), f_type, iter_type>;

  SymbolicFunctionOffsets() : /*e{}, */ offsets_data{nullptr}, len{0}, n{0} {}

  template <int... I0s, int... P0s>
  void reallocate(len_type len,
                  symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>) {
    n = len_type(sizeof...(I0s));
    this->len = len;

    delete[] offsets_data;
    offsets_data = new iter_type[len * n]{};
  }

  SymbolicFunctionOffsets(SymbolicFunctionOffsets const& other)
      : /*e{ other.e }, offsets{ (other.len > 0) ? new iter_type*[other.len] :
           nullptr },*/
        offsets_data{(other.len * other.n > 0)
                         ? new iter_type[other.len * other.n]
                         : nullptr},
        len{other.len},
        n{other.n} {
    std::copy(other.offsets_data, other.offsets_data + len * n, offsets_data);
  }

  SymbolicFunctionOffsets(SymbolicFunctionOffsets&& other)
      : SymbolicFunctionOffsets() {
    swap(*this, other);
  }

  SymbolicFunctionOffsets<E0, T0s...> operator=(
      SymbolicFunctionOffsets<E0, T0s...> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(SymbolicFunctionOffsets<E0, T0s...>& first,
                   SymbolicFunctionOffsets<E0, T0s...>& second) {
    std::swap(first.offsets_data, second.offsets_data);
    // std::swap(first.e, second.e);
    std::swap(first.len, second.len);
    std::swap(first.n, second.n);
  }

  ~SymbolicFunctionOffsets() { delete[] offsets_data; }

  iter_type* offsets_data;
  // f_type e;
  len_type len;
  len_type n;
};

template <typename... Ts>
struct Substitution : std::tuple<Ts...> {
  using parent_type = std::tuple<Ts...>;
  using parent_type::parent_type;

  template <typename E>
  auto eval(SymbolicFunction<E, Ts...> const& function) const {
    return substitute(function, std::make_index_sequence<sizeof...(Ts)>{});
  }

  template <typename E, size_t... Is>
  auto substitute(SymbolicFunction<E, Ts...> const& function) const {
    return substitute(function, std::make_index_sequence<sizeof...(Ts)>{});
  }

  template <typename E, size_t... Ns>
  auto substitute(SymbolicTemplate<E, Ns...> const& tmpl) const {
    return tmpl(std::move(std::get<Ns>(*this))...);
  }

  template <typename... T0s>
  operator Substitution<T0s...>() const {
    return cast<T0s...>(std::make_index_sequence<sizeof...(T0s)>{});
  }

  template <typename E0, typename... T0s>
  auto update(SymbolicFunction<E0, T0s...> const& f) {
    return f;
  }

  template <typename E0, typename... T0s>
  auto update(SymbolicFunction<E0, T0s...> const& f,
              SymbolicFunction<E0, T0s...>&) {
    return f;
  }

  template <typename E0, typename... T0s>
  void update(SymbolicFunction<E0, T0s...>& persistent) {}

 protected:
  template <typename... T0s, size_t... Is>
  Substitution<T0s...> cast(std::index_sequence<Is...>) const {
    return {std::get<Is>(*this)...};
  }

  template <typename E, size_t... Is>
  auto substitute(SymbolicFunction<E, Ts...> const& function,
                  std::index_sequence<Is...>) const {
    return function(std::get<Is>(*this)...);
  }
};

template <typename... Ts>
Substitution(Ts...) -> Substitution<Ts...>;
template <typename... Ts>
Substitution(std::tuple<Ts...>) -> Substitution<Ts...>;

namespace expr {

namespace {

template <typename E, size_t... ArgNs>
// auto get_substitutable_template(OpExpression<E> const& e,
// std::index_sequence<ArgNs...>)
SymbolicTemplate<E, ArgNs...> get_substitutable_template(
    OpExpression<E> const& e, std::index_sequence<ArgNs...>) {
  // return std::index_sequence<ArgNs...>{};
  return {*static_cast<const E*>(&e)};
  // auto ee = (expr::make_term<ArgNs>() + ...);
  // return SymbolicTemplate<decltype(ee), ArgNs...>(ee);
}

//! Create a substitutable function.
template <typename E>
auto to_template_expression(OpExpression<E> const& e, std::index_sequence<>,
                            std::index_sequence<>) {
  return *static_cast<const E*>(&e);
}

//! Create a substitutable function.
template <typename E, size_t N0, size_t... Ns, size_t I0, size_t... Is>
auto to_template_expression(OpExpression<E> const& e,
                            std::index_sequence<N0, Ns...>,
                            std::index_sequence<I0, Is...>) {
  return expr::transform::swap_grid<N0, Ns...>(
      *static_cast<const E*>(&e), Variable<I0>{}, Variable<Is>{}...);
}

//! Create a substitutable function.
template <typename E>
auto make_substitutable_template(OpExpression<E> const& e) {
  auto all_ids = independent_variables_of<E>;
  using seq = std::make_index_sequence<decltype(all_ids)::size()>;
  return get_substitutable_template(
      to_template_expression(*static_cast<const E*>(&e),
                             symphas::lib::sorted_seq<decltype(all_ids)>{},
                             seq{}),
      seq{});
}

//! Create a substitutable template.
/*!
 * Creates a substitutable template. The filter sequence list are the variables
 * that will not be used to generate the template. The remaining independent
 * variables are sorted and used to generate the template.
 */
template <typename E, size_t... Ns>
auto make_substitutable_template(OpExpression<E> const& e,
                                 std::index_sequence<Ns...> ids) {
  return get_substitutable_template(
      to_template_expression(*static_cast<const E*>(&e), ids, ids), ids);
}

template <size_t N0, size_t... Ns>
struct SymbolicTemplateDef<N0, Ns...> {
 protected:
  template <typename S>
  struct convert_vars;

  template <int... Ms>
  struct convert_vars<std::integer_sequence<int, Ms...>> {
    using type = std::index_sequence<size_t(Ms)...>;
  };

 public:
  auto operator=(OpVoid) { return OpVoid{}; }

  template <typename E>
  auto operator=(OpExpression<E> const& e) {
    return make_substitutable_template(*static_cast<const E*>(&e),
                                       std::index_sequence<N0, Ns...>{});
  }

  SymbolicTemplateDef(SymbolicTemplateDef<N0, Ns...> const&) = delete;
  SymbolicTemplateDef(SymbolicTemplateDef<N0, Ns...>&&) = delete;
};

template <>
struct SymbolicTemplateDef<> {
  auto operator=(OpVoid) { return OpVoid{}; }

  template <typename E>
  auto operator=(OpExpression<E> const& e) {
    return make_substitutable_template(*static_cast<const E*>(&e));
  }

  template <typename E, typename... Ts>
  auto operator=(SymbolicFunction<E, Ts...> const& e) {
    return make_substitutable_template(e.e);
  }

  SymbolicTemplateDef(SymbolicTemplateDef<> const&) = delete;
  SymbolicTemplateDef(SymbolicTemplateDef<>&&) = delete;
};

template <typename... Gs>
struct SymbolicTemplateDefSwap {
 public:
  template <typename E>
  auto swap_all_grids(OpExpression<E> const& e, std::index_sequence<>) {
    return *static_cast<E const*>(&e);
  }

  template <typename X, typename... Xs, typename E, size_t N0, size_t... Ns>
  auto swap_all_grids(OpExpression<E> const& e,
                      std::index_sequence<N0, Ns...>) {
    return expr::transform::swap_grid<X, Xs...>((*static_cast<E const*>(&e)),
                                                expr::symbols::arg_t<N0>{},
                                                expr::symbols::arg_t<Ns>{}...);
  }

 public:
  auto operator=(OpVoid) { return OpVoid{}; }

  template <
      typename E,
      typename Seq = decltype(expr::get_independent_variables(
          std::declval<E>() +
          std::declval<add_result_t<
              OpVoid, decltype(expr::make_term(std::declval<Gs&>()))...>>()))>
  auto operator=(OpExpression<E> const& e) {
    using seq_pick = std::make_index_sequence<sizeof...(Gs) * 2>;
    using seq_filtered = symphas::lib::filter_seq_t<seq_pick, Seq>;
    using seq_cut = symphas::lib::seq_lt_t<sizeof...(Gs), seq_filtered>;

    auto swapped = swap_all_grids<Gs...>(*static_cast<const E*>(&e), seq_cut{});
    return make_substitutable_template(swapped, seq_cut{});
  }

  SymbolicTemplateDefSwap(SymbolicTemplateDefSwap<Gs...> const&) = delete;
  SymbolicTemplateDefSwap(SymbolicTemplateDefSwap<Gs...>&&) = delete;
};

//! Create a template that can create an expression by substituting parameters
//! to placeholders.
/*!
 * Defines a template using the given placeholders.
 */
template <size_t... Ns, typename std::enable_if_t<all_ne<Ns...>, int> = 0>
SymbolicTemplateDef<Ns...> template_of_apply(symbols::arg_t<Ns>...) {
  return {};
}

//! Create a symbolic template where parameters are ::Variable types of inferred
//! indices.
/*!
 * The placeholders of the template are inferred from the existing variables in
 * the expression.
 */
SymbolicTemplateDef<> template_of_apply() { return {}; }

//! Create a symbolic template where parameters are ::Variable types of given
//! indices.
/*!
 * The linear variables (with OpIdentity coefficients) are given as arguments of
 * the function, and they must be defined in terms of a ::Variable.
 */
template <size_t... Ns, typename... Gs,
          typename std::enable_if_t<(all_ne<Ns...> && (is_symbol<Gs> && ...)),
                                    int> = 0>
SymbolicTemplateDef<Ns...> template_of_apply(
    OpTerm<OpIdentity, Variable<Ns, Gs>> const&...) {
  return {};
}

//! Create a symbolic template where parameters are ::Variable types of given
//! indices.
/*!
 * The linear variables (with OpIdentity coefficients) are given as arguments of
 * the function, and they must be defined in terms of a ::Variable.
 */
template <size_t... Ns, typename... Gs,
          typename std::enable_if_t<(all_ne<Ns...> && (is_symbol<Gs> && ...)),
                                    int> = 0>
SymbolicTemplateDef<Ns...> template_of_apply(Variable<Ns, Gs> const&...) {
  return {};
}

//! Create a function such that the parameters are expr::symbols::Symbol types.
/*!
 * Linear variables of the associated function expression which have their data
 * as one of the symbols defined in the function argument list will be
 * substituted.
 */
template <size_t... Ns, typename... symbol_ts,
          typename std::enable_if_t<
              (all_different<symbol_ts...> &&
               ((is_symbol<symbol_ts> && !is_id_variable<symbol_ts>) && ...)),
              int> = 0>
SymbolicTemplateDefSwap<symbol_ts...> template_of_apply(symbol_ts const&...) {
  return {};
}

//! Create a function such that the parameters are expr::symbols::Symbol types.
/*!
 * Linear variables of the associated function expression which have their data
 * as one of the symbols defined in the function argument list will be
 * substituted.
 */
template <size_t... Ns, typename... symbol_ts,
          typename std::enable_if_t<(is_symbol<symbol_ts> && ...), int> = 0>
SymbolicTemplateDefSwap<symbol_ts...> template_of_apply(
    OpTerms<OpIdentity, Term<symbol_ts, 1>> const&...) {
  return {};
}

template <size_t N0, typename G0, size_t... Ns, typename... Gs>
struct SymbolicFunctionDef<Variable<N0, G0>, Variable<Ns, Gs>...> {
  SymbolicFunctionDef(G0 const& arg0, Gs const&... args)
      : data{arg0, args...} {}
  SymbolicFunctionDef() : data{} {}

  template <typename E>
  auto get_function(SymbolicTemplate<E, N0, Ns...> const& tmpl) {
    using function_type =
        SymbolicFunction<E, Variable<N0, G0>, Variable<Ns, Gs>...>;
    auto f = typename function_type::type(tmpl);
    f.set_data_tuple(data);
    return f;
  }

  template <typename E>
  auto get_function(OpExpression<E> const& e) {
    using function_type =
        SymbolicFunction<E, Variable<N0, G0>, Variable<Ns, Gs>...>;
    auto f = typename function_type::type(*static_cast<E const*>(&e));
    f.set_data_tuple(data);
    return f;
  }

  template <typename E>
  auto operator=(OpExpression<E> const& e) {
    auto tmpl = (template_of(Variable<N0>{}, Variable<Ns>{}...) =
                     *static_cast<E const*>(&e));
    return get_function(tmpl);
  }

  template <typename E, size_t M0, size_t... Ms>
  auto operator=(SymbolicTemplate<E, M0, Ms...> const& e) {
    return (*this = e(Variable<N0>{}, Variable<Ns>{}...));
  }

  SymbolicFunctionDef(SymbolicFunctionDef<Variable<N0, G0>,
                                          Variable<Ns, Gs>...> const&) = delete;
  SymbolicFunctionDef(
      SymbolicFunctionDef<Variable<N0, G0>, Variable<Ns, Gs>...>&&) = delete;

  std::tuple<G0, Gs...> data;
};

template <typename G0, typename... Gs>
struct SymbolicFunctionDef<G0, Gs...> {
  using var_type = OpLiteral<double>;

  template <typename E, size_t N0, size_t... Ns>
  auto get_function(SymbolicTemplate<E, N0, Ns...> const& tmpl) {
    using function_type =
        SymbolicFunction<E, Variable<N0, var_type>, Variable<Ns, var_type>...>;
    return typename function_type::type(tmpl);
  }

  template <typename E, size_t I0, size_t... Is>
  auto get_function(OpExpression<E> const& e, std::index_sequence<I0, Is...>) {
    auto tmpl = (template_of(G0{}, Gs{}...) = *static_cast<E const*>(&e));
    return get_function(tmpl);
  }

  template <typename E>
  auto operator=(OpExpression<E> const& e) {
    return get_function(*static_cast<E const*>(&e),
                        std::make_index_sequence<sizeof...(Gs) + 1>{});
  }

  SymbolicFunctionDef(SymbolicFunctionDef<G0, Gs...> const&) = delete;
  SymbolicFunctionDef(SymbolicFunctionDef<G0, Gs...>&&) = delete;
};

template <>
struct SymbolicFunctionDef<> {
  template <typename E>
  auto operator=(OpExpression<E> const& e) {
    return SymbolicFunction<E>(*static_cast<E const*>(&e));
  }

  SymbolicFunctionDef(SymbolicFunctionDef<> const&) = delete;
  SymbolicFunctionDef(SymbolicFunctionDef<>&&) = delete;
};

//! Create a function such that parameters are interpreted from the expression.
/*!
 * There are no arguments passed to the function, and it is just a way to call
 * the expression evaluation.
 */
inline SymbolicFunctionDef<> function_of_apply() { return {}; }

//! Create a symbolic template where parameters are ::Variable types of given
//! indices.
/*!
 * The linear variables (with OpIdentity coefficients) are given as arguments of
 * the function, and they must be defined in terms of a ::Variable.
 */
template <size_t N0, typename G0, size_t... Ns, typename... Gs,
          typename std::enable_if_t<(all_ne<N0, Ns...> &&
                                     !(is_symbol<G0> || ... || is_symbol<Gs>)),
                                    int> = 0>
SymbolicFunctionDef<Variable<N0, G0>, Variable<Ns, Gs>...> function_of_apply(
    OpTerm<OpIdentity, Variable<N0, G0>> const& arg0,
    OpTerm<OpIdentity, Variable<Ns, Gs>> const&... args) {
  return {expr::get<1>(arg0).data(), expr::get<1>(args).data()...};
}

//! Create a function such that the parameters are expr::symbols::Symbol types.
/*!
 * Linear variables of the associated function expression which have their data
 * as one of the symbols defined in the function argument list will be
 * substituted.
 */
template <size_t N0, typename G0, size_t... Ns, typename... Gs,
          typename std::enable_if_t<
              (all_ne<N0, Ns...> /* && !(is_symbol<Gs> || ...)*/), int> = 0>
SymbolicFunctionDef<Variable<N0, G0>, Variable<Ns, Gs>...> function_of_apply(
    Variable<N0, G0> const arg0, Variable<Ns, Gs> const&... args) {
  return {arg0, args...};
}

inline SymbolicFunctionDef<> function_of_apply(symphas::lib::types_list<>) {
  return function_of_apply();
}

template <size_t N0, typename G0, size_t... Ns, typename... Gs,
          typename std::enable_if_t<(all_ne<N0, Ns...> &&
                                     !(is_symbol<G0> || ... || is_symbol<Gs>)),
                                    int> = 0>
SymbolicFunctionDef<Variable<N0, G0>, Variable<Ns, Gs>...> function_of_apply(
    symphas::lib::types_list<Variable<N0, G0>, Variable<Ns, Gs>...>) {
  return {};
}

//! Create a function such that the parameters are expr::symbols::Symbol types.
/*!
 * Linear variables of the associated function expression which have their data
 * as one of the symbols defined in the function argument list will be
 * substituted.
 */
template <typename symbol_t, typename... symbol_ts,
          typename std::enable_if_t<(all_different<symbol_t, symbol_ts...> &&
                                     (is_symbol<symbol_t> && ... &&
                                      is_symbol<symbol_ts>)),
                                    int> = 0>
SymbolicFunctionDef<symbol_ts...> function_of_apply(symbol_t const,
                                                    symbol_ts const&...) {
  return {};
}

//! Create a function such that the parameters are expr::symbols::Symbol types.
/*!
 * Linear variables of the associated function expression which have their data
 * as one of the symbols defined in the function argument list will be
 * substituted.
 */
template <typename symbol_t, typename... symbol_ts,
          typename std::enable_if_t<(all_different<symbol_t, symbol_ts...> &&
                                     (is_symbol<symbol_t> && ... &&
                                      is_symbol<symbol_ts>)),
                                    int> = 0>
SymbolicFunctionDef<symbol_ts...> function_of_apply(
    OpTerms<OpIdentity, Term<symbol_t, 1>> const,
    OpTerms<OpIdentity, Term<symbol_ts, 1>> const&...) {
  return {};
}
}  // namespace

template <typename... Ts>
decltype(auto) template_of(Ts&&... ts) {
  return template_of_apply(std::forward<Ts>(ts)...);
}

template <typename... Ts>
decltype(auto) function_of(Ts&&... ts) {
  return function_of_apply(std::forward<Ts>(ts)...);
}

}  // namespace expr

namespace expr {

template <typename... T0s>
auto arg_list(T0s&&... args) {
  return symphas::internal::make_substitution(
      symphas::internal::as_arg(std::forward<T0s>(args))...);
}

template <typename... T0s>
auto arg_array_list(len_type len, T0s&&... args) {
  return symphas::internal::make_substitution(
      symphas::internal::as_array_arg(std::forward<T0s>(args), len)...);
}

template <typename T0>
auto array_arg(len_type len, T0&& arg) {
  return symphas::internal::as_array_arg(std::forward<T0>(arg), len);
}

template <typename T>
std::reference_wrapper<T> bind(T const& arg) {
  return std::ref(const_cast<T&>(arg));
}

template <typename T>
T unbind(T const& arg) {
  return arg;
}
}  // namespace expr

namespace expr::symbols {

template <typename... Ts>
auto fn(Ts&&... ts) {
  return expr::function_of(std::forward<Ts>(ts)...);
}

template <typename E, size_t Z0, typename G0>
auto D(SymbolicFunction<E, Variable<Z0, G0>> const& f) {
  return expr::apply_operators(expr::make_operator_derivative<1, Z0>()(f.e));
}

template <typename E, size_t Z0, typename G0, size_t O>
auto D(SymbolicFunction<E, Variable<Z0, G0>> const& f,
       OpFractionLiteral<O, 1>) {
  return expr::apply_operators(expr::make_operator_derivative<O, Z0>()(f.e));
}

template <typename E, size_t Z0, typename G0>
auto D(SymbolicFunction<E, Variable<Z0, G0>> const& f, OpIdentity) {
  return expr::apply_operators(expr::make_operator_derivative<1, Z0>()(f.e));
}

template <size_t O, typename E, size_t... Zs, size_t... Is,
          size_t N = sizeof...(Is)>
auto D(OpExpression<E> const& e, std::index_sequence<Zs...>,
       std::index_sequence<Is...>) {
  return ((expr::make_column_vector<Is, N>() *
           expr::apply_operators(expr::make_operator_derivative<O, Zs>()(
               *static_cast<E const*>(&e)))) +
          ...);
}

template <typename E, size_t Z0, typename G0, size_t Z1, typename G1,
          size_t... Zs, typename... Gs>
auto D(SymbolicFunction<E, Variable<Z0, G0>, Variable<Z1, G1>,
                        Variable<Zs, Gs>...> const& f) {
  return D<1>(f.e, std::index_sequence<Z0, Z1, Zs...>{},
              std::make_index_sequence<sizeof...(Zs) + 2>{});
}

template <typename E, size_t Z0, typename G0, size_t Z1, typename G1,
          size_t... Zs, typename... Gs, size_t O>
auto D(SymbolicFunction<E, Variable<Z0, G0>, Variable<Z1, G1>,
                        Variable<Zs, Gs>...> const& f,
       OpFractionLiteral<O, 1>) {
  return D<O>(f.e, std::index_sequence<Z0, Z1, Zs...>{},
              std::make_index_sequence<sizeof...(Zs) + 2>{});
}

template <typename E, size_t Z0, typename G0, size_t Z1, typename G1,
          size_t... Zs, typename... Gs>
auto D(SymbolicFunction<E, Variable<Z0, G0>, Variable<Z1, G1>,
                        Variable<Zs, Gs>...> const& f,
       OpIdentity) {
  return D<1>(f.e, std::index_sequence<Z0, Z1, Zs...>{},
              std::make_index_sequence<sizeof...(Zs) + 2>{});
}
}  // namespace expr::symbols

template <typename V, typename F, typename... Args>
struct OpCallable : OpExpression<OpCallable<V, F, Args...>> {
 protected:
  template <size_t... Is>
  auto _eval(iter_type n, std::index_sequence<Is...>) const {
    return expr::eval(value) * f(std::get<Is>(args)...);
  }

 public:
  OpCallable(V const& value, F const& f, Args const&... args)
      : value(value), f(f), args(args...) {}

  auto eval(iter_type n) const {
    return _eval(n, std::make_index_sequence<sizeof...(Args)>{});
  }

  auto operator-() const;

  V value;
  F f;
  std::tuple<Args...> args;
};

template <typename V, typename F, typename... Args>
struct OpCallable<V, F*, Args...> : OpExpression<OpCallable<V, F*, Args...>> {
 protected:
  template <size_t... Is>
  auto _eval(iter_type n, std::index_sequence<Is...>) const {
    return expr::eval(value) * f->operator()(std::get<Is>(args)...);
  }

 public:
  OpCallable(V const& value, F* f, Args const&... args)
      : value(value), f(f), args(args...) {}

  auto eval(iter_type n) const {
    return _eval(n, std::make_index_sequence<sizeof...(Args)>{});
  }

  auto operator-() const;

  V value;
  F* f;
  std::tuple<Args...> args;
};

template <typename V, typename F, typename... Args>
OpCallable(V, F, Args...) -> OpCallable<V, F, Args...>;

template <typename V, typename F, typename... Args>
OpCallable(V, F*, Args...) -> OpCallable<V, F*, Args...>;

template <typename F>
struct SymbolicCallable {
  SymbolicCallable(F const& f) : f(f) {}
  template <typename... Args>
  auto operator()(Args&&... args) const {
    return OpCallable(OpIdentity{}, f, std::forward<Args>(args)...);
  }
  F f;
};

template <typename F>
struct SymbolicCallable<F*> {
  SymbolicCallable(F* f) : f(f) {}
  template <typename... Args>
  auto operator()(Args&&... args) const {
    return OpCallable(OpIdentity{}, f, std::forward<Args>(args)...);
  }
  F* f;
};

/* multiplication of a derivative object by a literal
 */

namespace expr {
template <typename F>
auto callable_of(F const& f) {
  return SymbolicCallable<F>(f);
}

template <typename F>
auto callable_of(F* f) {
  return SymbolicCallable<F*>(f);
}

template <typename V, typename F, typename... Args, size_t... Is>
auto make_callable_of(V&& coeff, F&& f, std::tuple<Args...> const& args,
                      std::index_sequence<Is...>) {
  return OpCallable(std::forward<V>(coeff), std::forward<F>(f),
                    std::get<Is>(args)...);
}
template <typename V, typename F, typename... Args>
auto make_callable_of(V&& coeff, F&& f, std::tuple<Args...> const& args) {
  return make_callable_of(std::forward<V>(coeff), std::forward<F>(f), args,
                          std::make_index_sequence<sizeof...(Args)>{});
}
}  // namespace expr

template <typename V, typename F, typename... Args>
auto OpCallable<V, F, Args...>::operator-() const {
  return expr::make_callable_of(-value, f, args);
}

template <typename coeff_t, typename V, typename F, typename... Args,
          typename std::enable_if_t<
              (expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpCallable<V, F, Args...> const& b) {
  return expr::make_callable_of(value * b.value, b.f, b.args);
}

template <typename coeff_t, typename tensor_t, typename F, typename... Args,
          typename std::enable_if_t<
              (expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value,
               OpCallable<tensor_t, F, Args...> const& b) {
  return (value * b.value) * expr::make_callable_of(OpIdentity{}, b.f, b.args);
}

//! @}
