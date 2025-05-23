
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
 * MODULE:  sol
 * PURPOSE: Defines the model used for phase field crystal problems.
 *
 * ***************************************************************************
 */

#pragma once

#include <array>
#include <set>
#include <utility>

#include "solver.h"

namespace symphas::internal {
template <typename T>
struct field_array_t;

namespace parameterized {

struct ExpressionStats {
 protected:
  template <typename T, size_t D>
  static auto max(any_vector_t<T, D> const &v1, any_vector_t<T, D> const &v2) {
    using std::max;
    any_vector_t<T, D> vm;
    for (iter_type i = 0; i < D; ++i) {
      vm[i] = max(v1[i], v2[i]);
    }
    return vm;
  }

  template <typename T, size_t D>
  static auto min(any_vector_t<T, D> const &v1, any_vector_t<T, D> const &v2) {
    using std::min;
    any_vector_t<T, D> vm;
    for (iter_type i = 0; i < D; ++i) {
      vm[i] = min(v1[i], v2[i]);
    }
    return vm;
  }

 public:
  template <typename E>
  auto sum(OpExpression<E> const &e) const {
    expr::eval_type_t<E> sum;
    expr::result(OpVoid{}, sum);

    auto sumop = expr::make_term(sum);

    len_type len = expr::data_length(*static_cast<E const *>(&e));
    for (iter_type i = 0; i < len; ++i) {
      sum += (*static_cast<E const *>(&e)).eval(i);
    }

    return expr::make_literal(sum);
  }

  template <typename E>
  auto mean(OpExpression<E> const &e) const {
    len_type len = expr::data_length(*static_cast<E const *>(&e));
    return sum(*static_cast<E const *>(&e)) * (1. / len);
  }

  template <typename E>
  auto max(OpExpression<E> const &e) const {
    using std::max;

    expr::eval_type_t<E> result;
    expr::result((*static_cast<E const *>(&e)).eval(0), result);

    len_type len = expr::data_length(*static_cast<E const *>(&e));
    for (iter_type i = 1; i < len; ++i) {
      result = max((*static_cast<E const *>(&e)).eval(i), result);
    }

    return expr::make_literal(result);
  }

  template <typename E>
  auto min(OpExpression<E> const &e) const {
    using std::max;

    expr::eval_type_t<E> result;
    expr::result((*static_cast<E const *>(&e)).eval(0), result);

    len_type len = expr::data_length(*static_cast<E const *>(&e));
    for (iter_type i = 1; i < len; ++i) {
      result = min((*static_cast<E const *>(&e)).eval(i), result);
    }

    return expr::make_literal(result);
  }
};

//! Real valued order parameter.
/*!
 * The order parameter will be of type ::scalar_t.
 */
using SCALAR = scalar_t;

//! Complex value order parameter.
/*!
 * The order parameter will be of type ::complex_t.
 */
using COMPLEX = complex_t;

//! Vector valued order parameter of real elements.
/*!
 * The order parameter will a vector type, specified by ::vector_t.
 */
template <size_t Dm>
using VECTOR_D = vector_t<Dm>;

//! Vector valued order parameter of real elements.
/*!
 * The order parameter will a vector type, specified by ::vector_t, this type
 * in particular is used to infer the vector.
 */
using VECTOR = vector_t<0>;

//! Indicates that the field number is selected from configuration.
/*!
 * The number of fields of the type this option is passed to is selected using
 * the configuration. The configuration will select the total number of fields.
 * Only a single field can have this option applied.
 */
constexpr int CONFIGURATION = -1;

//! The imaginary number for an expression.
/*!
 * An imaginary number literal that can be used in an expression.
 */
inline OpLiteral<complex_t> Ii = expr::make_literal(complex_t(0, 1));

//! Access statistics about the system, such as the mean.
/*!
 * Included statistics are: mean, max, min and sum.
 */
inline ExpressionStats STATS;

constexpr expr::symbols::i_<0, 0> ii;
constexpr expr::symbols::i_<1, 0> jj;
constexpr expr::symbols::i_<2, 0> kk;

//! Take the imaginary part of a complex number resulting from an expression.
/*!
 * Take the imaginary part of the complex number that is computed by
 * evaluating the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
template <typename E>
decltype(auto) Im(E &&e) {
  return expr::imag(std::forward<E>(e));
}

//! Take the imaginary part of a complex number resulting from an expression.
/*!
 * Take the real part of the complex number that is computed by
 * evaluating the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
template <typename E>
decltype(auto) Re(E &&e) {
  return expr::real(std::forward<E>(e));
}

//! Use a number in an expression.
/*!
 * Create a number with the given value to use in an expression.
 * This is a constant.
 *
 * \param V The value of the number.
 */
template <typename E>
decltype(auto) lit(E &&e) {
  return expr::make_literal(std::forward<E>(e));
}

//! Define an array that will be iterated using the given index.
/*!
 * This object is only applicable in the context of a index-defined EOM or
 * free energy. The array is defined using an expression which can use the same
 * index, in which case the index of the array will take over. The array is
 * generally used when multiple distinct instances have to be generated for each
 * field, such as a noise term.
 */
template <int N, int P>
decltype(auto) ARRAY(expr::symbols::i_<N, P>) {
  return symphas::internal::indexed_array(expr::symbols::i_<N, P>{});
}

//! Define a domain integral of the expression.
/*!
 * Writes the integral of the given expression.
 */
template <typename E>
decltype(auto) INT(E &&e) {
  return expr::make_domain_integral(std::forward<E>(e));
}

template <typename E>
decltype(auto) T(E &&e) {
  return expr::transpose(std::forward<E>(e));
}

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
decltype(auto) NOISE(symphas::grid_info const &info, const double *dt,
                     enclosing_class_wrap<grid_type>) {
  return NoiseData<nt, T, D, grid_type>(info.get_dims(), info.get_widths(), dt);
}

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
decltype(auto) NOISE(symphas::grid_info const &info,
                     grid::region_interval<D> const &region, const double *dt,
                     enclosing_class_wrap<grid_type>) {
  len_type dims[D]{};
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = region[i][1] - region[i][0];
  }
  return NoiseData<nt, T, D, grid_type>(dims, info.get_widths(), dt);
}

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
decltype(auto) NOISE(symphas::grid_info const &info,
                     grid::region_size const &region, const double *dt,
                     enclosing_class_wrap<grid_type>) {
  return NoiseData<nt, T, D, grid_type>(info.get_dims(), info.get_widths(), dt);
}

template <size_t D, typename T>
struct dimensionalized {
  using type = T;
};

template <size_t D, typename T>
struct dimensionalized<D, any_vector_t<T, 0>> {
  using type = any_vector_t<T, D>;
};

template <size_t D, typename T>
using dimensionalized_t = typename dimensionalized<D, T>::type;

}  // namespace parameterized

template <typename T, typename S = void>
struct parameterized_type;

template <typename M, typename T>
struct model_field_parameters;

template <>
struct parameterized_type<void, void> {
  parameterized_type() : parameters{nullptr}, len{0} {}

  template <typename... Ts>
  parameterized_type(double parameter0, Ts... parameters)
      : parameters{new double[sizeof...(Ts) + 1]{}}, len{sizeof...(Ts) + 1} {
    populate_parameters(std::make_index_sequence<sizeof...(Ts) + 1>{},
                        parameter0, parameters...);
  }

  parameterized_type(parameterized_type<void, void> const &other)
      : parameters{(other.len > 0) ? new double[other.len]{} : nullptr},
        len{other.len} {
    std::copy(other.parameters, other.parameters + other.len, parameters);
  }

  parameterized_type(parameterized_type<void, void> &&other)
      : parameterized_type() {
    swap(*this, other);
  }

  template <typename M, typename T>
  parameterized_type(parameterized_type<M, T> const &other)
      : parameters{(other.len > 0) ? new double[other.len]{} : nullptr},
        len{other.len} {
    std::copy(other.parameters, other.parameters + other.len, parameters);
  }

  template <typename M, typename T>
  parameterized_type(parameterized_type<M, T> &&other) : parameterized_type() {
    swap(*this, *static_cast<parameterized_type<void, void> *>(&other));
  }

  parameterized_type<void, void> &operator=(
      parameterized_type<void, void> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(parameterized_type<void, void> &first,
                   parameterized_type<void, void> &second) {
    using std::swap;
    swap(first.parameters, second.parameters);
    swap(first.len, second.len);
  }

  const parameterized_type<void, void> &operator()() const { return *this; }

  parameterized_type<void, void> &operator()() { return *this; }

  double operator[](iter_type i) const {
    if (i >= len) {
      return DEFAULT_COEFF_VALUE;
    } else {
      return parameters[i];
    }
  }

  double *parameters;
  len_type len;

  ~parameterized_type() { delete[] parameters; }

 protected:
  auto as_parameter(double parameter) { return parameter; }

  template <size_t I0, size_t... Is, typename T0, typename... Ts>
  void populate_parameters(std::index_sequence<I0, Is...>, T0 &&parameter0,
                           Ts &&...parameters) {
    this->parameters[I0] = as_parameter(std::forward<T0>(parameter0));
    ((this->parameters[Is] = as_parameter(std::forward<Ts>(parameters))), ...);
  }

  inline void populate_parameters(std::index_sequence<>) {}
};

template <typename T, typename S>
struct parameterized_type<parameterized_type<T, S>, void>
    : parameterized_type<T, S> {
  using parent_type = parameterized_type<T, S>;
  using parent_type::parent_type;
};

template <typename T>
struct parameterized_type<field_array_t<T>, void>
    : parameterized_type<T, void> {
  using parent_type = parameterized_type<T, void>;
  using parent_type::parent_type;
};

template <typename S>
struct parameterized_type<S, void> : parameterized_type<void, void> {
  template <typename M>
  parameterized_type(M const &model)
      : parameterized_type<void>(model_field_parameters<M, S>(model)) {}

  parameterized_type() : parameterized_type<void>() {}

  const parameterized_type<void, void> &operator()() const { return *this; }

  parameterized_type<void, void> &operator()() { return *this; }
};

template <typename M, typename T>
struct model_field_parameters : model_field_parameters<M, void>,
                                parameterized_type<void, void> {
  using parent_type = model_field_parameters<M, void>;
  model_field_parameters(M const &m)
      : parent_type{m}, parameterized_type<void, void>{m.param(0)} {}
};

template <typename M>
struct model_field_parameters<M, void> {
  model_field_parameters(M const &m) : model{&m} {}

  template <size_t I = 0>
  auto op() const {
    return model->template op<I>();
  }

  template <size_t I = 0>
  auto dop() const {
    return model->template dop<I>();
  }

  template <int N, int P>
  auto op_i(expr::symbols::i_<N, P>) const {
    return model->op_i(expr::symbols::i_<N, P>{});
  }

  auto param(size_t I) const { return model->param(I); }

  template <size_t N>
  auto param_matrix() const {
    return model->template param_matrix<N>();
  }

  auto c(size_t I) const { return param(I - 1); }

  template <expr::NoiseType nt, typename T, typename... T0s>
  auto make_noise(T0s &&...args) const {
    return model->template make_noise<nt, T>(std::forward<T0s>(args)...);
  }

  template <expr::NoiseType nt, typename G, typename... T0s>
  auto make_noise(OpTerm<OpIdentity, G> const &term, T0s &&...args) const {
    return model->template make_noise<nt>(term, std::forward<T0s>(args)...);
  }
  const M *model;
};

template <size_t D, typename T>
struct non_parameterized_type_impl {
  using type = T;
};

template <size_t D, typename T, typename S>
struct non_parameterized_type_impl<D, parameterized_type<T, S>> {
  using type = S;
};

template <size_t D, typename T>
struct non_parameterized_type_impl<
    D, parameterized_type<T, parameterized::VECTOR>> {
  using type = parameterized::VECTOR_D<D>;
};

template <size_t D, typename S>
struct non_parameterized_type_impl<D, symphas::internal::field_array_t<S>> {
  using type = typename non_parameterized_type_impl<D, S>::type;
};

template <size_t D, typename T>
using non_parameterized_type = typename non_parameterized_type_impl<D, T>::type;

}  // namespace symphas::internal

template <template <typename> typename M>
struct model_field_name_format {
#ifdef LATEX_PLOT
  static constexpr const char *value = "\\phi_{%d}";
#else
  static constexpr const char *value = "phi_%d";
#endif
};

template <template <typename> typename M>
struct model_field_name_builder {
  const char *operator()(int index) const {
    static auto names = build_names(model_field_name_format<M>::value);
    return names[index];
  }

  decltype(auto) operator()() const {
    return model_field_name_format<M>::value;
  }

  char **build_names(const char *value) const {
    constexpr size_t len = 100;
    char **names = new char *[len];
    char buffer[100];
    for (size_t i = 1; i <= len; ++i) {
      sprintf(buffer, value, int(i));
      names[i - 1] = new char[std::strlen(buffer) + 1];
      std::strcpy(names[i - 1], buffer);
    }
    return names;
  }

  template <size_t N>
  decltype(auto) build_names(const char *(&list)[N]) const {
    return list;
  }
};

//! Used to obtain names of phase fields for a given model.
/*!
 * Used for getting the name of the phase field at the given index, typically
 * for printing to output or for using the name in expressions.
 *
 * Specialization of this class by model type is used to associate different
 * phase field names than the default ones. All models use default phase field
 * names unless specialized.
 */
template <template <typename> typename M>
struct model_field_name_default : model_field_name_builder<M> {
  using model_field_name_builder<M>::operator();

  //! Returns the name of the \f$i\f$-th phase field.
  /*!
   * Get the name of the phase field at the given index for the prescribed
   * model.
   *
   * \param i The index of the phase field.
   */
  const char *operator()(int i) {
    constexpr size_t MAX_NAME_COUNT =
        sizeof(ORDER_PARAMETER_NAMES) / sizeof(*ORDER_PARAMETER_NAMES);

    if (i < MAX_NAME_COUNT) {
      return ORDER_PARAMETER_NAMES[i];
    } else {
      static size_t EXTRA_NAME_COUNT = 0;
      static std::vector<char *> EXTRA_NAMES;

      if (i - MAX_NAME_COUNT < EXTRA_NAME_COUNT) {
        return EXTRA_NAMES[i - MAX_NAME_COUNT];
      } else {
        char *name = new char[BUFFER_LENGTH];
        snprintf(name, BUFFER_LENGTH, ORDER_PARAMETER_NAME_EXTRA_FMT, i);
        EXTRA_NAMES.push_back(name);
        return EXTRA_NAMES.back();
      }
    }
  }
};

template <template <typename> typename M>
struct model_field_name : model_field_name_default<M> {};

// *****************************************************************************************

#ifdef VTK_ON
#include <thread>

#include "colormap.h"
#endif

template <size_t D, typename Sp, typename... Ts>
using ArrayModel = Model<D, Sp, symphas::internal::field_array_t<void>, Ts...>;

template <size_t D, typename Sp, typename... S>
struct Model;

namespace symphas::internal {

template <typename S>
constexpr bool is_field_array_type = false;
template <typename S>
constexpr bool is_field_array_type<field_array_t<S>> = true;
template <typename T, typename S>
constexpr bool is_field_array_type<parameterized_type<T, S>> =
    is_field_array_type<S>;

template <typename S>
struct modifier_save {
  virtual void update(const S *_s, len_type len) = 0;
  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir) const = 0;
  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir, const char *name) const = 0;
  virtual const S *get_first(const S *_s) const = 0;
  virtual S *get_first(S *_s) = 0;
  virtual ~modifier_save() {}
};

template <typename S, ModelModifiers modifier = ModelModifiers::PLOT_DEFAULT>
struct modifier_save_apply : modifier_save<S> {
  modifier_save_apply(const S *_s, len_type len,
                      const symphas::interval_data_type *vdata, size_t N)
      : N{N} {}

  virtual void update(const S *_s, len_type len) override {}

  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir) const override {
    if (len > 0) {
      for (iter_type i = 0; i < len; ++i) {
        _s[i].save(dir, index);
      }
    }
  }

  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir, const char *name) const override {
    if (len > 0) {
      char **names = new char *[len];
      for (iter_type i = 0; i < len; ++i) {
        names[i] =
            new char[std::strlen(name) + symphas::lib::num_digits(i) + 1];
        snprintf(names[i], BUFFER_LENGTH, "#%zd-%d-%s", N, i, name);
      }

      for (iter_type i = 0; i < len; ++i) {
        _s[i].save(dir, names[i], index);
      }

      for (iter_type i = 0; i < len; ++i) {
        delete[] names[i];
      }
      delete[] names;
    }
  }

  virtual const S *get_first(const S *_s) const override { return _s; }

  virtual S *get_first(S *_s) override { return _s; }

  size_t N;
};

#ifdef USING_CUDA
template <typename T, size_t D, typename grid_type>
void update_max_grid(const SolverSystemFDwSDCUDA<T, D> *_s,
                     grid_type const &s_max, len_type len);

template <typename T, size_t D, typename grid_type>
void update_max_grid(const SolverSystemFDCUDA<T, D> *_s, grid_type const &s_max,
                     len_type len);
#endif

template <typename S>
struct modifier_save_apply<S, ModelModifiers::PLOT_MAX> : modifier_save<S> {
  modifier_save_apply(const S *_s, len_type len,
                      const symphas::interval_data_type *vdata, size_t N)
      : s_max{symphas::init_data_type({{Axis::NONE, Inside::CONSTANT}}), *vdata,
              symphas::b_data_type(), (len > 0) ? _s[0].id : 0},
        N{N} {}

  virtual void update(const S *_s, len_type len) override { _update(_s, len); }

  template <template <typename, size_t> typename solver_system_type, typename T,
            size_t D>
  void _update(const solver_system_type<T, D> *_s, len_type len) {
    auto s_max_it_start = symphas::data_iterator_region(
        s_max.as_grid(), grid::get_iterable_domain(s_max));
    auto s_max_it_end = s_max_it_start + grid::length_interior(s_max);

    for (auto it(s_max_it_start); it < s_max_it_end; ++it) {
      *it = OpVoid{};
    }
    for (iter_type i = 0; i < len; ++i) {
      auto it0 = symphas::data_iterator_region(
          _s[i].as_grid(), grid::get_iterable_domain(s_max));
      for (auto it(s_max_it_start); it < s_max_it_end; ++it, ++it0) {
        *it = std::max(*it, *it0);
      }
    }
  }

#ifdef USING_CUDA

  template <typename T, size_t D>
  void _update(const SolverSystemFDwSDCUDA<T, D> *_s, len_type len) {
    update_max_grid(_s, s_max, len);
  }
  template <typename T, size_t D>
  void _update(const SolverSystemFDCUDA<T, D> *_s, len_type len) {
    update_max_grid(_s, s_max, len);
  }

#endif

  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir) const override {
    if (len > 0) {
      s_max.save(dir, index);
    }
  }

  virtual void operator()(const S *_s, len_type len, iter_type index,
                          const char *dir, const char *name) const override {
    if (len > 0) {
      char *save_name =
          new char[std::strlen(name) + symphas::lib::num_digits(N) + 1];
      snprintf(save_name, BUFFER_LENGTH, "#%zd-max-%s", N, name);

      s_max.save(dir, save_name, index);
      delete[] save_name;
    }
  }

  virtual const S *get_first(const S *_s) const override { return &s_max; }

  virtual S *get_first(S *_s) override { return &s_max; }

  S s_max;
  size_t N;
};

template <size_t D, typename... S>
struct list_array_types;

template <size_t D, typename... array_ts, typename T, typename... Rest>
struct list_array_types<D, symphas::lib::types_list<array_ts...>, T, Rest...> {
  using type = typename list_array_types<
      D,
      symphas::lib::expand_types_list<
          array_ts..., std::conditional_t<
                           symphas::internal::is_field_array_type<T>,
                           symphas::lib::types_list<
                               symphas::internal::non_parameterized_type<D, T>>,
                           symphas::lib::types_list<>>>,
      Rest...>::type;
};

template <size_t D, typename... array_ts>
struct list_array_types<D, symphas::lib::types_list<array_ts...>> {
  using type = symphas::lib::types_list<array_ts...>;
};

template <size_t D, typename... S>
using list_array_t =
    typename list_array_types<D, symphas::lib::types_list<>, S...>::type;

template <typename Sp, size_t D, typename... S>
struct modifier_save_types;

template <typename Sp, size_t D, typename... array_ts>
struct modifier_save_types<Sp, D, symphas::lib::types_list<array_ts...>> {
  using type = std::tuple<modifier_save<typename symphas::solver_system_type<
      Sp>::template type<array_ts, D>> *...>;
};

template <typename Sp, size_t D, typename... S>
struct modifier_save_types {
  using type = typename modifier_save_types<Sp, D, list_array_t<D, S...>>::type;
};

template <typename Sp, size_t D, typename... S>
using modifier_save_t = typename modifier_save_types<Sp, D, S...>::type;

template <size_t D, typename Sp, typename... Ts>
auto num_fields_as_literal(ArrayModel<D, Sp, Ts...> const &model) {
  return expr::make_literal(model.len);
}

template <size_t D, typename Sp, typename... S>
auto num_fields_as_literal(Model<D, Sp, S...> const &model) {
  return expr::val<sizeof...(S)>;
}

}  // namespace symphas::internal

namespace symphas {

template <size_t D, typename Sp, typename... Ts>
auto model_num_fields(ArrayModel<D, Sp, Ts...> const &model) {
  size_t num_types = sizeof...(Ts);
  return (model.len > num_types) ? model.len : num_types;
}

template <size_t D, typename Sp, typename... S>
auto model_num_fields(Model<D, Sp, S...> const &model) {
  return sizeof...(S);
}
}  // namespace symphas

//! A representation of the phase field crystal model.
/*!
 * A representation of the phase field crystal model. Implements the phase
 * field crystal model equation, and allows
 *
 * \tparam PFC The phase field crystal parameters specialized object.
 * \tparam D The dimension of the phase field crystal problem.
 * \tparam Sp The solver type for numerically solving the phase field crystal
 * problem.
 * \tparam Ts... The types of the order parameters.
 */
template <size_t D, typename Sp, typename... Ts>
struct Model<D, Sp, symphas::internal::field_array_t<void>, Ts...> {
  //! Get the type of the phase-fields, takes an argument for comptability with
  //! other Model.
  /*!
   * Get the type of the `N`th system, but since they are all the same type,
   * will return the type directly regardless of the value of `N`.
   *
   * \tparam N The phase field index to get the type.
   */
  template <size_t N>
  using type_of_S = symphas::lib::type_at_index<
      N, symphas::internal::non_parameterized_type<D, Ts>...>;

  using first_type_of_S = symphas::lib::type_at_index<
      0, symphas::internal::non_parameterized_type<D, Ts>...>;

  using all_field_types = symphas::lib::types_list<
      symphas::internal::non_parameterized_type<D, Ts>...>;

  template <typename Type>
  static constexpr bool model_has_type =
      symphas::lib::index_of_type<
          Type, symphas::lib::expand_types_list<all_field_types>> >= 0;

  template <typename Type, size_t I = 0>
  static constexpr int index_of_type = symphas::lib::nth_index_of_type<
      Type, I, symphas::lib::expand_types_list<all_field_types>>;

  static constexpr size_t num_array_types =
      (size_t(symphas::internal::is_field_array_type<Ts>) + ...);

  //! The type of the system storing the phase fields, used by the solver.
  template <size_t N = 0>
  using SolverSystemApplied =
      typename symphas::solver_system_type<Sp>::template type<type_of_S<N>, D>;
  using save_method_t = symphas::internal::modifier_save_t<Sp, D, Ts...>;

 protected:
  template <size_t I>
  void construct_save_type(const symphas::interval_data_type *vdata) {
    size_t n = num_field_offset<I>();
    len_type len = num_type_fields[I];

    using namespace symphas::internal;
    switch (plot_type[I]) {
      case symphas::ModelModifiers::PLOT_DEFAULT:
        (std::get<I>(save_method) =
             new modifier_save_apply<SolverSystemApplied<I>,
                                     symphas::ModelModifiers::PLOT_DEFAULT>(
                 _s + n, len, vdata, I));
        break;
      case symphas::ModelModifiers::PLOT_MAX:
        (std::get<I>(save_method) =
             new modifier_save_apply<SolverSystemApplied<I>,
                                     symphas::ModelModifiers::PLOT_MAX>(
                 _s + n, len, vdata, I));
        break;
      case symphas::ModelModifiers::PLOT_MIN:
        (std::get<I>(save_method) =
             new modifier_save_apply<SolverSystemApplied<I>,
                                     symphas::ModelModifiers::PLOT_MIN>(
                 _s + n, len, vdata, I));
        break;
      case symphas::ModelModifiers::PLOT_CONTOURS:
        (std::get<I>(save_method) =
             new modifier_save_apply<SolverSystemApplied<I>,
                                     symphas::ModelModifiers::PLOT_CONTOURS>(
                 _s + n, len, vdata, I));
        break;
      default:
        (std::get<I>(save_method) =
             new modifier_save_apply<SolverSystemApplied<I>,
                                     symphas::ModelModifiers::PLOT_DEFAULT>(
                 _s + n, len, vdata, I));
    }
  }

  template <size_t I>
  void apply_save_type(iter_type index, const char *dir) const {
    size_t n = num_field_offset<I>();
    len_type len = num_type_fields[I];

    std::get<I>(save_method)->operator()(_s + n, len, index, dir);
  }

  template <size_t... Is>
  void apply_save_types(iter_type index, const char *dir,
                        std::index_sequence<Is...>) const {
    (apply_save_type<Is>(index, dir), ...);
  }

  void apply_save_types(iter_type index, const char *dir) const {
    if (symphas::parallel::is_host_node()) {
      apply_save_types(index, dir, std::make_index_sequence<num_array_types>{});
    }
  }

  template <size_t I>
  void apply_save_type(iter_type index, const char *dir,
                       const char *name) const {
    size_t n = num_field_offset<I>();
    len_type len = num_type_fields[I];

    std::get<I>(save_method)->operator()(_s + n, len, index, dir, name);
  }

  template <size_t... Is>
  void apply_save_types(iter_type index, const char *dir, const char *name,
                        std::index_sequence<Is...>) const {
    (apply_save_type<Is>(index, dir, name), ...);
  }

  void apply_save_types(iter_type index, const char *dir,
                        const char *name) const {
    if (symphas::parallel::is_host_node()) {
      apply_save_types(index, dir, name,
                       std::make_index_sequence<num_array_types>{});
    }
  }

  template <size_t I>
  void update_save_type() const {
    size_t n = num_field_offset<I>();
    len_type len = num_type_fields[I];

    std::get<I>(const_cast<save_method_t &>(save_method))->update(_s + n, len);
  }

  template <size_t... Is>
  void update_save_types(std::index_sequence<Is...>) const {
    (update_save_type<Is>(), ...);
  }

  void update_save_types() const {
    SolverSystemApplied<>::synchronize(
        (void *)NULL,
        const_cast<Model<D, Sp, symphas::internal::field_array_t<void>, Ts...>
                       *>(this));
    if (symphas::parallel::is_host_node()) {
      update_save_types(std::make_index_sequence<num_array_types>{});
    }
  }

  template <size_t... Is>
  void construct_save_types(const symphas::interval_data_type *vdata,
                            std::index_sequence<Is...>) {
    (construct_save_type<Is>(vdata), ...);
  }

  void construct_save_types(const symphas::interval_data_type *vdata) {
    construct_save_types(vdata, std::make_index_sequence<num_array_types>{});
    update_save_types();
  }

  template <size_t I>
  void delete_save_type() {
    delete std::get<I>(save_method);
  }

  template <size_t... Is>
  void delete_save_types(std::index_sequence<Is...>) {
    (delete_save_type<Is>(), ...);
  }

  void delete_save_types() {
    delete_save_types(std::make_index_sequence<num_array_types>{});
  }

  Model()
      : len{0},
        index{params::start_index},
        _s{construct_systems({}, {}, {}, 0, 0)},
        solver{Sp::make_solver()},
        coeff{nullptr},
        num_coeff{0},
        time{0},
        num_type_fields{},
        plot_type{},
        save_method{}
#ifdef VTK_ON
        ,
        viz_update{nullptr}
#endif
  {
    for (iter_type i = 0; i < num_array_types; ++i) {
      plot_type[i] = symphas::ModelModifiers::PLOT_DEFAULT;
      num_type_fields[i] = 0;
    }
    construct_save_types({});
  }

 public:
  //! Create a new model.
  /*!
   * A new model is created using the problem parameters. The model
   * coefficients for the dynamical equations are also stored.
   *
   * \param coeff The list of coefficients for this model.
   * \param num_coeff The number of coefficients provided.
   * \param parameters The problem parameters for the phase field problem
   * represented by the model.
   */
  Model(double const *coeff, size_t num_coeff,
        symphas::problem_parameters_type const &parameters)
      : len{static_cast<len_type>(parameters.length())},
        _s{construct_systems(parameters.get_initial_data(),
                             parameters.get_interval_data(),
                             parameters.get_boundary_data(),
                             parameters.length(), parameters.length())},
        index{parameters.index},
        solver{Sp::make_solver(get_updated_parameters(parameters))},
        coeff{(num_coeff > 0) ? new double[num_coeff] : nullptr},
        num_coeff{num_coeff},
        time{parameters.time},
        num_type_fields{},
        plot_type{},
        save_method{}
#ifdef VTK_ON
        ,
        viz_update{nullptr}
#endif
  {
    iter_type i0;
    i0 = static_cast<iter_type>(
        std::min(num_array_types, parameters.get_num_fields_len()));
    for (iter_type i = 0; i < i0; ++i) {
      plot_type[i] = parameters.get_modifiers()[i];
      num_type_fields[i] = parameters.get_num_fields()[i];
    }
    i0 = static_cast<iter_type>(
        std::min(num_array_types, parameters.get_num_fields_len()));
    for (iter_type i = i0; i < num_array_types; ++i) {
      plot_type[i] = symphas::ModelModifiers::PLOT_DEFAULT;
      num_type_fields[i] = 0;
    }

    construct_save_types(parameters.get_interval_data());
    std::copy(coeff, coeff + num_coeff, this->coeff);
    visualize();
  }

  //! Create a new model without coefficients.
  /*!
   * The coefficient list is assumed to be empty, so no coefficients are
   * saved.
   *
   * \param parameters The problem parameters for the phase field problem
   * represented by the model.
   */
  Model(symphas::problem_parameters_type const &parameters)
      : Model(nullptr, 0, parameters) {}

  Model(ArrayModel<D, Sp, Ts...> const &other)
      : len{other.len},
        _s{construct_systems(other._s, other.len)},
        solver{other.solver},
        coeff{(other.num_coeff > 0) ? new double[other.num_coeff] : nullptr},
        num_coeff{other.num_coeff},
        index{other.index},
        time{other.time},
        num_type_fields{},
        plot_type{},
        save_method{other.save_method}
#ifdef VTK_ON
        ,
        viz_update{nullptr}
#endif
  {
    std::copy(other._s, other._s + other.len, _s);
    std::copy(other.coeff, other.coeff + other.num_coeff, coeff);
    std::copy(other.plot_type, other.plot_type + num_array_types, plot_type);
    std::copy(other.num_type_fields, other.num_type_fields + num_array_types,
              num_type_fields);
    visualize();
  }

  Model(ArrayModel<D, Sp, Ts...> &&other) noexcept : Model() {
    swap(*this, other);
  }

  ArrayModel<D, Sp, Ts...> &operator=(ArrayModel<D, Sp, Ts...> other) {
    swap(*this, other);
    return *this;
  }

  ~Model() {
    for (iter_type i = len - 1; i >= 0; --i) {
      std::destroy_at(&_s[i]);
    }
    operator delete[](_s);
    delete_save_types();

    devisualize();
    delete[] coeff;
  }

  friend void swap(ArrayModel<D, Sp, Ts...> &first,
                   ArrayModel<D, Sp, Ts...> &second) {
    using std::swap;

    swap(first.solver, second.solver);
    swap(first._s, second._s);
    swap(first.coeff, second.coeff);
    swap(first.len, second.len);
    swap(first.num_coeff, second.num_coeff);
    swap(first.index, second.index);
    swap(first.time, second.time);
    swap(first.num_type_fields, second.num_type_fields);
    swap(first.plot_type, second.plot_type);
    swap(first.save_method, second.save_method);

#ifdef VTK_ON
    swap(first.viz_thread, second.viz_thread);
    swap(first.viz_update, second.viz_update);
#endif
  }

  template <typename Sp0,
            typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
  friend void swap(ArrayModel<D, Sp, Ts...> &first,
                   ArrayModel<D, Sp0, Ts...> &second) {
    using std::swap;

    swap(first._s, second._s);
    swap(first.coeff, second.coeff);
    swap(first.len, second.len);
    swap(first.num_coeff, second.num_coeff);
    swap(first.index, second.index);
    swap(first.time, second.time);
    swap(first.num_type_fields, second.num_type_fields);
    swap(first.plot_type, second.plot_type);
    swap(first.save_method, second.save_method);

#ifdef VTK_ON
    swap(first.viz_thread, second.viz_thread);
    swap(first.viz_update, second.viz_update);
#endif
  }

  //! Returns the list of coefficients.
  /*!
   * Returns a pointer to the array of coefficients used in the dynamic
   * equations.
   */
  const double *get_coeff() const { return coeff; }

  //! Get the current solution index.
  /*!
   * At each complete iteration of the solver, the index is incremented to
   * track the total number of solver iterations performed. Returns this
   * value.
   */
  iter_type get_index() const { return index; }

  //! Return the number of coefficients saved by the model.
  /*!
   * Returns the number of coefficients that are used by the model in
   * the dynamic equations.
   */
  auto get_num_coeff() const { return num_coeff; }

  //! Returns the current simulation time of the model.
  /*!
   * Returns the current simulation time of the model.
   */
  auto get_time() const { return time; }

  void set_index(iter_type index) { this->index = index; }

  void set_time(double time) { this->time = time; }

  //! Returns a symbolic variable for the time value of the model.
  /*!
   * Returns a symbolic variable for the time value of the model. The
   * variable maintains the time value as a pointer to the time value managed
   * by the model.
   */
  auto get_time_var() const { return expr::make_term(TimeValue{&time}); }

  //! Get the number of fields of the given types.
  /*!
   * Return the number of fields which correspond to the given types for this
   * model.
   *
   * \tparam S0 The type which is counted from the list of all types.
   */
  template <typename S0>
  size_t num_fields() const {
    return num_fields<S0>(std::make_index_sequence<sizeof...(Ts)>{});
  }

  //! Get the number of fields of the given types.
  /*!
   * Return the number of fields which correspond to the given types for this
   * model. The number returned is a sum of the number of times each type
   * appears in the list.
   *
   * \tparam S0 The type which is counted from the list of all types.
   * \tparam S1 The second type which is counted.
   */
  template <typename S0, typename S1, typename... Ss>
  size_t num_fields() const {
    return num_fields<S0>() + num_fields<S1, Ss...>();
  }

  //! Get the number of fields for the ith type.
  /*!
   * Return the number of fields which are used by type i.
   */
  template <size_t N>
  size_t num_fields() const {
    return num_field<N>();
  }

  //! Get the number of fields for the ith type.
  /*!
   * Return the number of fields which are used by type i.
   */
  size_t num_fields(iter_type i) const { return num_field(i); }

  //! Get the index of the `I`-th given type.
  /*!
   * Get the index of the selected type from the list of types for this model.
   * The index is the position in the list, starting from 0. Determines the
   * index of the type after appearing `I` times.
   *
   * \tparam Type The type to find in the list and return its position index.
   * \tparam I Chose the `I`-th appearance of the chosen type.
   */
  // template<typename Type, size_t I = 0>
  // int index_of_type = (std::is_same<Type, S>::value) ? I : -1;

  //! Execute a function for all fields of the given type.
  /*!
   * Execute a function for all fields of the given type.
   * The function expects as parameters the axes values, the data values and
   * the index of the corresponding field, followed by desired arguments.
   *
   * \param f The function to execute on the field of the given type.
   * \param args The arguments to pass to the function.
   *
   * \tparam Type The field type that is selected to be passed to the
   * function.
   */
  template <typename Type, typename F, typename... Args,
            typename std::enable_if_t<model_has_type<Type>, int> = 0>
  void do_for_field_type(F f, Args &&...args) const {
    do_for_field_type<Type, 0>(f, std::forward<Args>(args)...);
  }

  //! Execute a function for the `N`-th field.
  /*!
   * Execute a function for all fields of the given type.
   * The function expects as parameters the axes values, the data values and
   * the index of the corresponding field, followed by desired arguments.
   *
   * \param f The function to execute on the field of the given type.
   * \param args The arguments to pass to the function.
   *
   * \tparam N the index of the field that is passed to the function.
   */
  template <size_t N, typename F, typename... Args>
  void do_for_field(F f, Args &&...args) const {
    iter_type i0 = num_field_offset<N>();
    for (iter_type i = 0; i < num_field<N>(); ++i) {
      f(_s[i + i0], std::forward<Args>(args)...);
    }
  }

  //! Copy the values of the `N`-th field into the given array.
  /*!
   * Copies the values from the field into the given array, ensuring that
   * only the true phase field values are copied. I.e. this means that for
   * systems which use boundaries, the boundary domain is not included.
   *
   * \param into The array into which to copy the values.
   *
   * \tparam N The index of the phase field which is copied from.
   */
  template <size_t N>
  void copy_field_values(first_type_of_S *into) const {
    _s[N].persist(into);
  }

  //! Copy the values of the `I`-th field of the given type.
  /*!
   * Copy values from the field corresponding to the given type. The `I`-th
   * field of that type is chosen from the list. Depending on the system type,
   * the boundary domain will not be copied.
   *
   * \param[in] into The array into which the field is copied into.
   *
   * \tparam Type is the type to match and copy into.
   */
  template <typename Type, size_t I,
            typename std::enable_if_t<std::is_same<Type, type_of_S<I>>::value,
                                      int> = 0>
  void copy_field_type_values(Type *into) const {
    copy_field_values<I>(into);
  }

  //! Copy the values of the `I`-th field of the given type.
  /*!
   * Copy values from the field corresponding to the given type. The `I`-th
   * field of that type is chosen from the list. Depending on the system type,
   * the boundary domain will not be copied.
   *
   * \param[in] into The array into which the field is copied into.
   *
   * \tparam Type is the type to match and copy into.
   */
  template <typename Type, size_t I,
            typename std::enable_if_t<!std::is_same<Type, type_of_S<I>>::value,
                                      int> = 0>
  void copy_field_type_values(Type *) const {}

  //! Return the values of the `N`-th field as a Grid.
  /*!
   * Copies the values from the field and returns a new Grid, ensuring that
   * only the true phase field values are copied. I.e. this means that for
   * systems which use boundaries, the boundary domain is not included.
   *
   * \tparam N The index of the phase field which is copied from.
   */
  template <size_t N>
  auto get_field() const {
    Grid<first_type_of_S, D> out(_s[N].get_info().get_dims().get());
    _s[N].persist(out.values);
    return out;
  }

  //! Persists all the systems managed by the model.
  /*!
   * Persists the phase field systems to disk, saving the entire array.
   * This function has no effect without the **io** module.
   *
   * \param dir The directory into which to store the resulting data file.
   */
  void save_systems(const char *dir = nullptr) const {
    update_save_types();
    apply_save_types(index, dir);

    iter_type n = num_field_offset<num_array_types>();
    for (iter_type i = 0; i < len - n; ++i) {
      _s[i + n].save(dir, index);
    }
  }

  //! Persists all the systems managed by the model.
  /*!
   * Persists the phase field systems to disk, saving the entire array.
   * This function has no effect without the **io** module.
   *
   * \param dir The directory into which to store the resulting data file.
   * \param name The suffix name, to which the system id number is appended.
   */
  void save_systems(const char *dir, const char *name) const {
    update_save_types();
    apply_save_types(index, dir, name);

    iter_type n = num_field_offset<num_array_types>();

    char *names[sizeof...(Ts) - num_array_types];
    for (iter_type i = 0; i < sizeof...(Ts) - num_array_types; ++i) {
      names[i] =
          new char[std::strlen(name) + symphas::lib::num_digits(n + i) + 1];
      snprintf(names[i], BUFFER_LENGTH, "%d%s", n + i, name);
    }

    for (iter_type i = 0; i < sizeof...(Ts) - num_array_types; ++i) {
      _s[i + n].save(dir, names[i], index);
    }

    for (iter_type i = 0; i < sizeof...(Ts) - num_array_types; ++i) {
      delete[] names[i];
    }
  }

  //! Updates each of the phase field systems.
  /*!
   * Updates the systems by directly calling their `update` function.
   *
   * \param time The current time of the solution.
   */
  void update_systems(double time) {
    this->time = time;

    SYMPHAS_OMP_PARALLEL_DIRECTIVE
    for (iter_type i = 0; i < len; ++i) {
      _s[i].update(index, time);
    }
  }

  //! Typically implemented for parallelization routines.
  /*!
   * Synchronizes the state of the systems using the list of solver objects, and
   * the algorithm is selected based on the implemented solver and the
   * implemented system.
   *
   * \param es The list of equation objects generated and used by the solver.
   */
  template <typename E>
  void synchronize_systems(E *es) {
    SolverSystemApplied<>::synchronize(es, this);
  }

  //! Advances to the next solution iteration.
  /*!
   * Calls the `step` function of the solver, stepping forward the solution
   * iteration. This also increments Moodel::index.
   *
   * \param dt The time increment.
   */
  void step(double dt) {
    ++index;
    solver.dt = dt;

    auto range = symphas::parallel::get_index_range(len);
    iter_type start = range.first;
    iter_type end = range.second;

    SYMPHAS_OMP_PARALLEL_DIRECTIVE
    for (iter_type i = start; i < end; ++i) {
      solver.step(system(i));
    }

#ifdef VTK_ON
    if (viz_update != nullptr && index % params::viz_interval == 0) {
      get_viz_field();
      viz_update->update();
    }
#endif
  }

  //! Returns the system at the given index.
  /*!
   * Returns the reference to the system at the given index from the list of
   * phase field systems of the model.
   *
   * \tparam I The index of the phase field to return.
   */
  template <size_t I>
  const auto &system() const {
    return _s[I];
  }

  //! Returns the system at the given index.
  /*!
   * Returns the reference to the system at the given index from the list of
   * phase field systems of the model.
   *
   * \tparam I The index of the phase field to return.
   */
  template <size_t I>
  auto &system() {
    return const_cast<SolverSystemApplied<> &>(
        static_cast<const ArrayModel<D, Sp, Ts...> &>(*this)
            .template system<I>());
  }

  //! Returns the underlying grid of the system at the given index.
  /*!
   * Returns a reference to the grid of the system at the given index.
   *
   * \tparam I The index of the phase field to return.
   */
  template <size_t I>
  auto &grid() {
    return system<I>().as_grid();
  }

  //! Returns the underlying grid of the system at the given index.
  /*!
   * Returns a reference to the grid of the system at the given index.
   *
   * \tparam I The index of the phase field to return.
   */
  template <size_t I>
  const auto &grid() const {
    return system<I>().as_grid();
  }

  //! Returns the system at the given index.
  /*!
   * Returns the reference to the system at the given index from the list of
   * phase field systems of the model.
   *
   * \tparam I The index of the phase field to return.
   */
  const auto &system(iter_type i) const { return _s[i]; }

  //! Returns the system at the given index.
  /*!
   * Returns the reference to the system at the given index from the list of
   * phase field systems of the model.
   *
   * \tparam I The index of the phase field to return.
   */
  auto &system(iter_type i) { return _s[i]; }

  //! Returns the underlying grid of the system at the given index.
  /*!
   * Returns a reference to the grid of the system at the given index.
   *
   * \tparam I The index of the phase field to return.
   */
  auto &grid(iter_type i) { return system(i).as_grid(); }

  //! Returns the underlying grid of the system at the given index.
  /*!
   * Returns a reference to the grid of the system at the given index.
   *
   * \tparam I The index of the phase field to return.
   */
  const auto &grid(iter_type i) const { return system(i).as_grid(); }

  //! Returns the tuple of all model systems.
  /*!
   * Directly returns the list of all systems of this model.
   */
  decltype(auto) systems_tuple() const { return std::make_pair(_s, len); }

  //! Returns the tuple of all model systems.
  /*!
   * Directly returns the list of all systems of this model.
   */
  decltype(auto) systems_tuple() { return std::make_pair(_s, len); }

  //! Returns the tuple of all model systems.
  /*!
   * Directly returns the list of all systems of this model.
   */
  const auto &systems() const { return _s; }

  //! Returns the tuple of all model systems.
  /*!
   * Directly returns the list of all systems of this model.
   */
  auto &systems() { return _s; }

  //! Print information about this model to the given output.
  /*!
   * The information about this model is written to the given output.
   * Data about this model includes the number of coefficients and the number
   * of fields.
   *
   * Since the name is not directly associated with the model, the name is
   * not provided in this information.
   *
   * \param out The output file or stream to which to print the
   * information.
   */
  void print_info(FILE *out) const {
    fprintf(out, OUTPUT_BANNER);
    fprintf(out, "Phase field problem of %d system%s:\n", len,
            (len > 1) ? "s" : "");
    print_all_dimensions(out);
    print_all_coeff(out);
    fprintf(out, OUTPUT_BANNER);
  }

  auto generate_parameters() const {
    symphas::problem_parameters_type parameters(len);
    populate_parameters(parameters);

    parameters.time = get_time();
    ;
    parameters.index = index;

    return parameters;
  }

  len_type len;  //!< Number of fields.

 protected:
  //! The latest index of simulation.
  /*
   * Internally tracks the current index in the simulation. It is
   * incremented for every call to Model::step(), and can be accessed
   * (but not modified) through the function Model::index().
   */
  iter_type index;

  SolverSystemApplied<>
      *_s;  //!< Container managing pointers to phase field data.

  Sp solver;  //! Solver for determining the phase field solution.

  //! Coefficients in the equations of motion.
  /*!
   * List of coefficients used in the equations of motion. If
   * the number of coefficients is insufficient, then the default coefficient
   * value #DEFAULT_COEFF_VALUE should be used.
   */
  double *coeff;

  //! Number of supplied coefficients.
  /*!
   * The length of the coefficient list, coeff.
   */
  size_t num_coeff;

  //! The current simulation time.
  double time;

  len_type num_type_fields[num_array_types];  //!< Number of fields of each
                                              //!< array type.
  symphas::ModelModifiers
      plot_type[num_array_types];  //!< How array fields are plotted.
  save_method_t save_method;

#ifdef VTK_ON

  std::thread viz_thread;
  ColourPlotUpdater *viz_update;

  SolverSystemApplied<> *get_viz_field(std::index_sequence<>) {
    update_save_type<num_array_types - 1>();
    return std::get<num_array_types - 1>(save_method)
        ->get_first(_s + num_field_offset<num_array_types - 1>());
  }

  template <size_t I0, size_t... Is>
  SolverSystemApplied<> *get_viz_field(std::index_sequence<I0, Is...>) {
    if (I0 == params::viz_index) {
      update_save_type<I0>();
      return std::get<I0>(save_method)->get_first(_s + num_field_offset<I0>());
    } else {
      return get_viz_field(std::index_sequence<Is...>{});
    }
  }

  SolverSystemApplied<> *get_viz_field() {
    return get_viz_field(std::make_index_sequence<num_array_types>{});
  }

  void visualize() {
    if constexpr (std::is_same<type_of_S<>, scalar_t>::value) {
      SolverSystemApplied<> *ss = get_viz_field();
      if (params::viz_enabled && params::viz_interval > 0) {
        viz_thread = std::thread(
            [](auto *viz_grid, int *index, auto **viz_update) {
              ColourPlot2d viz{};
              viz.init(viz_grid->values, grid::get_data_domain(*viz_grid),
                       *index, *viz_update);
            },
            &(ss->as_grid()), &index, &viz_update);
      }
    }
  }

  void devisualize() {
    if constexpr (std::is_same<type_of_S<>, scalar_t>::value) {
      if (params::viz_enabled && params::viz_interval > 0) {
        viz_thread.join();
      }
    }
  }

#else
  void visualize() {}
  void devisualize() {}
#endif

  auto get_updated_parameters(
      symphas::problem_parameters_type parameters) const {
    parameters.extend(len);
    for (iter_type i = 0; i < len; ++i) {
      parameters.set_interval_data(system(i).get_info().intervals, i);
    }
    return parameters;
  }

  template <size_t... Is>
  auto populate_parameters(symphas::problem_parameters_type &parameters) const {
    for (iter_type i = 0; i < len; ++i) {
      fill_interval_data(system(i), parameters, i);
      fill_boundary_data(system(i), parameters, i);
      fill_initial_data(system(i), parameters, i);
    }
  }

  void fill_interval_data(
      PhaseFieldSystem<BoundaryGrid, first_type_of_S, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const {
    auto intervals = system.get_info().intervals;
    for (iter_type i = 0; i < D; ++i) {
      Axis ax = symphas::index_to_axis(i);
      auto &interval = intervals.at(ax);

      interval.set_count(interval.left(), interval.right(), system.dims[i]);
    }

    parameters.set_interval_data(intervals, n);
  }

  void fill_interval_data(
      PhaseFieldSystem<RegionalGrid, first_type_of_S, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const {
    auto intervals = system.get_info().intervals;
    for (iter_type i = 0; i < D; ++i) {
      Axis ax = symphas::index_to_axis(i);
      auto &interval = intervals.at(ax);

      interval.set_count(interval.left(), interval.right(), system.dims[i]);
    }

    parameters.set_interval_data(intervals, n);
  }

#ifdef USING_CUDA
  void fill_interval_data(
      PhaseFieldSystem<BoundaryGridCUDA, first_type_of_S, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const;

  void fill_interval_data(
      PhaseFieldSystem<RegionalGridCUDA, first_type_of_S, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const;
#endif

  template <template <typename, size_t> typename G>
  void fill_interval_data(PhaseFieldSystem<G, first_type_of_S, D> const &system,
                          symphas::problem_parameters_type &parameters,
                          iter_type n) const {
    parameters.set_interval_data(system.get_info().intervals, n);
  }

  template <template <typename, size_t> typename G, typename T0>
  void fill_boundary_data(PhaseFieldSystem<G, T0, D> const &system,
                          symphas::problem_parameters_type &parameters,
                          iter_type n) const {
    symphas::b_data_type bdata;

    for (iter_type i = 0; i < D * 2; ++i) {
      Side side = symphas::index_to_side(i);
      bdata[side] = symphas::b_element_type(BoundaryType::PERIODIC);
    }

    parameters.set_boundary_data(bdata, n);
  }

  template <typename T0>
  void fill_boundary_data(PhaseFieldSystem<BoundaryGrid, T0, D> const &system,
                          symphas::problem_parameters_type &parameters,
                          iter_type n) const {
    symphas::b_data_type bdata;

    for (iter_type i = 0; i < D * 2; ++i) {
      Side side = symphas::index_to_side(i);
      BoundaryType type = system.types[i];
      bdata[side] = system.boundaries[i]->get_parameters();
    }

    parameters.set_boundary_data(bdata, n);
  }

  template <typename T0>
  void fill_boundary_data(PhaseFieldSystem<RegionalGrid, T0, D> const &system,
                          symphas::problem_parameters_type &parameters,
                          iter_type n) const {
    symphas::b_data_type bdata;

    for (iter_type i = 0; i < D * 2; ++i) {
      Side side = symphas::index_to_side(i);
      BoundaryType type = system.types[i];
      bdata[side] = system.boundaries[i]->get_parameters();
    }

    parameters.set_boundary_data(bdata, n);
  }

#ifdef USING_CUDA

  template <typename T0>
  void fill_boundary_data(
      PhaseFieldSystem<BoundaryGridCUDA, T0, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const;

  template <typename T0>
  void fill_boundary_data(
      PhaseFieldSystem<RegionalGridCUDA, T0, D> const &system,
      symphas::problem_parameters_type &parameters, iter_type n) const;
#endif

  template <template <typename, size_t> typename G, typename T0>
  void fill_initial_data(PhaseFieldSystem<G, T0, D> const &system,
                         symphas::problem_parameters_type &parameters,
                         iter_type n) const {
    parameters.set_initial_data(init_data_functor(system.as_grid()), n);
  }

  auto construct_systems(const symphas::init_data_type *tdata,
                         const symphas::interval_data_type *vdata,
                         const symphas::b_data_type *bdata, size_t data_len,
                         size_t len) const {
    void *systems = operator new[](len * sizeof(SolverSystemApplied<>));
    SolverSystemApplied<> *ptr = static_cast<SolverSystemApplied<> *>(systems);

    auto [start, end] = symphas::parallel::get_index_range(len);

    iter_type end0 =
        std::min(static_cast<iter_type>(data_len), static_cast<iter_type>(len));
    for (iter_type i = 0; i < end0; ++i) {
      auto tdata_i = tdata[i];
      if (i < start || i >= end) {
        for (auto &[axis, init] : tdata_i) {
          init = symphas::init_entry_type();
        }
      }

      new (&ptr[i]) SolverSystemApplied<>(tdata_i, vdata[i], bdata[i], i);
      ptr[i].update(index, 0);
    }
    for (iter_type i = end0; i < len; ++i) {
      auto tdata_i = tdata[0];
      if (i < start || i >= end) {
        for (auto &[axis, init] : tdata_i) {
          init = symphas::init_entry_type();
        }
      }

      new (&ptr[i]) SolverSystemApplied<>(tdata_i, vdata[0], bdata[0], i);
      ptr[i].update(index, 0);
    }

    return ptr;
  }

  template <
      typename other_sys_type,
      typename std::enable_if_t<
          !std::is_same<other_sys_type, SolverSystemApplied<>>::value, int> = 0>
  SolverSystemApplied<> construct_system(other_sys_type const &system,
                                         iter_type n) const {
    symphas::b_data_type boundaries{};
    for (iter_type i = 0; i < D * 2; ++i) {
      boundaries[symphas::index_to_side(i)] = BoundaryType::PERIODIC;
    }

    bool extend_boundary = params::extend_boundary;
    params::extend_boundary = true;

    auto new_system =
        SolverSystemApplied<>(symphas::init_data_type{}, system.info.intervals,
                              boundaries, system.id);

    params::extend_boundary = extend_boundary;

    first_type_of_S *swap_arr = new first_type_of_S[system.length()];
    system.persist(swap_arr);
    new_system.fill(swap_arr);

    delete[] swap_arr;

    return new_system;
  }

  template <typename S0>
  SolverSystemApplied<> *construct_systems(S0 *const &systems,
                                           len_type len) const {
    if (len > 0) {
      void *new_systems = operator new[](len * sizeof(SolverSystemApplied<>));
      SolverSystemApplied<> *ptr =
          static_cast<SolverSystemApplied<> *>(new_systems);

      for (iter_type i = 0; i < len; ++i) {
        new (&ptr[i]) SolverSystemApplied<>(systems[i]);
        ptr[i].update(index, 0);
      }
      return ptr;
    } else {
      return nullptr;
    }
  }

  //! Returns the number of fields corresponding to the field type at index N.
  template <size_t N>
  size_t num_field() const {
    using type_n = symphas::lib::type_at_index<N, Ts...>;
    if constexpr (symphas::internal::is_field_array_type<type_n>) {
      return num_type_fields[N];
    } else {
      return 1;
    }
  }

  template <size_t... Is>
  size_t num_field(std::index_sequence<Is...>, iter_type i) const {
    return (((i == Is) ? ((symphas::internal::is_field_array_type<
                              symphas::lib::type_at_index<Is, Ts...>>)
                              ? num_type_fields[i]
                              : 1)
                       : 0) +
            ...);
  }

  //! Returns the number of fields corresponding to the field type at index i.
  size_t num_field(iter_type i) const {
    return num_field(std::make_index_sequence<sizeof...(Ts)>{}, i);
  }

  template <size_t... Is>
  size_t num_field_offset(std::index_sequence<Is...>) const {
    return (num_field<Is>() + ... + 0);
  }

  //! Returns the number of fields corresponding to the field type at index N.
  template <size_t N>
  size_t num_field_offset() const {
    return num_field_offset(std::make_index_sequence<N>{});
  }

  template <typename S0, size_t... Is>
  size_t num_fields(std::index_sequence<Is...>) const {
    return (((std::is_same<S0, type_of_S<Is>>::value) ? num_field<Is>() : 0) +
            ...);
  }

  template <typename Type, size_t... Is, typename F, typename... Args>
  void do_for_field_type(std::index_sequence<Is...>, F &&f,
                         Args &&...args) const {
    (do_for_field<index_of_type<Type, Is>>(f, std::forward<Args>(args)...),
     ...);
  }

  void print_all_dimensions(FILE *out) const {
    fprintf(out, "\tsystem<*> ... ");
    for (iter_type d = 0; d < D; ++d) {
      fprintf(out, "%d", grid(0).dims[d]);
      if (d < D - 1) {
        fprintf(out, " x ");
      }
    }
    fprintf(out, "\n");
  }

  void print_all_coeff(FILE *out) const {
    if (num_coeff > 0) {
      fprintf(out, "Number of coefficients provided: %zd.\n", num_coeff);
      fprintf(out, "List of those coefficients:\n");
      for (iter_type i = 0; i < num_coeff; ++i) {
        fprintf(out, "\tc%-8d ... % .4E\n", i + 1, coeff[i]);
      }
    } else {
      fprintf(out, "No coefficients provided.\n");
    }
    fprintf(out, "Default coefficient value: % .2lf.\n", DEFAULT_COEFF_VALUE);
  }
};
