
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
 * MODULE:  lib
 * PURPOSE: Defines random noise, which represents, for instance, thermal
 * fluctuations in a field.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"
#include "expressiontypek.h"
#include "symboliceval.h"

//! \cond

#define NOISE_MEAN 0
#define NOISE_STD 1

//! \endcond

// **************************************************************************************

namespace expr {

template <typename T, size_t D, template <typename, size_t> typename grid_type,
          typename E>
auto make_poisson_event(
    NoiseData<expr::NoiseType::POISSON, T, D, grid_type> const& noise,
    OpExpression<E> const& e);

template <size_t D, template <typename, size_t> typename grid_type, typename E>
auto make_poisson_event(
    expr::noise_data<expr::NoiseType::POISSON, D, grid_type> const& noise,
    OpExpression<E> const& e) {
  return make_poisson_event(
      NoiseData<expr::NoiseType::POISSON, scalar_t, D, grid_type>(noise),
      *static_cast<E const*>(&e));
}

template <typename T, size_t D, template <typename, size_t> typename grid_type,
          typename E, typename... Ts>
auto make_poisson_event(
    NoiseData<expr::NoiseType::POISSON, T, D, grid_type> const& noise,
    SymbolicFunction<E, Ts...> const& f);

template <size_t D, template <typename, size_t> typename grid_type, typename E,
          typename... Ts>
auto make_poisson_event(
    expr::noise_data<expr::NoiseType::POISSON, D, grid_type> const& noise,
    SymbolicFunction<E, Ts...> const& f) {
  return make_poisson_event(
      NoiseData<expr::NoiseType::POISSON, scalar_t, D, grid_type>(noise), f);
}

template <NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename... Ts>
auto make_noise(NoiseData<nt, T, D, grid_type> const& noise,
                SymbolicFunction<E, Ts...> const& f);

template <NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E>
auto make_noise(NoiseData<nt, T, D, grid_type> const& noise,
                OpExpression<E> const& e);

template <NoiseType nt, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename... Ts>
auto make_noise(expr::noise_data<nt, D, grid_type> const& noise,
                SymbolicFunction<E, Ts...> const& f) {
  return make_noise(NoiseData<nt, scalar_t, D, grid_type>(noise), f);
}

template <NoiseType nt, typename T, size_t D>
auto make_noise(const len_type* dimensions, const double* h, const double* dt);

}  // namespace expr

namespace symphas::internal {

template <size_t... Is, typename E, typename... Ts>
auto build_function_for_noise(std::index_sequence<Is...>,
                              OpExpression<E> const& e, Ts const&... args) {
  auto f = expr::function_of(Variable<Is, Ts>{}...) =
      *static_cast<E const*>(&e);
  f.set_data(args...);
  return f;
}

template <typename E, typename... Ts>
auto build_function_for_noise(OpExpression<E> const& e, Ts&&... args) {
  return build_function_for_noise(std::make_index_sequence<sizeof...(Ts)>{},
                                  *static_cast<E const*>(&e),
                                  std::forward<Ts>(args)...);
}
}  // namespace symphas::internal

namespace expr {

template <typename G>
struct random_state {
  void allocate() {}
  random_state(len_type = 0) {}
  void* states;
};

template <NoiseType nt>
struct random_seed {
  random_seed(int index = -1) : index{index} {}

  auto get_seed() const {
    if (index < 0) {
      return std::mt19937{std::random_device{}()};
    } else {
      static std::map<int, std::mt19937> gen_map;
      if (gen_map.count(index) > 0) {
        return gen_map[index];
      } else {
        auto [it, _] =
            gen_map.emplace(index, std::mt19937{std::random_device{}()});
        auto [k, v] = *it;
        return v;
      }
    }
  }

  void set_seed_index(int index) { this->index = index; }

  int index;
};

template <typename F, size_t D, template <typename, size_t> typename grid_type>
struct noise_data_with_function : grid_type<complex_t, D>,
                                  random_state<grid_type<complex_t, D>> {
  using parent_type = grid_type<complex_t, D>;
  using state_type = random_state<grid_type<complex_t, D>>;
  using parent_type::dims;
  using parent_type::len;
  using parent_type::values;

  void allocate() {
    if (parent_type::len == 0) {
      parent_type grid(parent_type::dims);
      using std::swap;
      swap(*static_cast<parent_type*>(this), grid);
    }

    state_type::allocate();
  }

  noise_data_with_function(F f, const len_type* dims, const double* h,
                           const double* dt)
      : parent_type(nullptr),
        state_type(grid::length<D>(dims)),
        dt{dt},
        h{0},
        f{f} {
    std::copy(h, h + D, this->h);
    std::copy(dims, dims + D, parent_type::dims);
  }

  scalar_t operator[](iter_type n) const { return values[n].real(); }

 protected:
  template <NoiseType nt, typename eval_handler_type>
  void update(random_seed<nt> const& seed, double intensity, bool fourier_space,
              eval_handler_type const& eval_handler) {
    using std::exp;
    using std::sqrt;

    auto gen = seed.get_seed();
    std::normal_distribution<> dis(NOISE_MEAN, NOISE_STD);

    k_field<D>::template fill<Axis::NONE, 3>(*static_cast<parent_type*>(this),
                                             dims, h);

    double H = 1;
    for (iter_type i = 0; i < D; ++i) H *= h[i];
    scalar_t R = sqrt(*dt) / H;

    for (iter_type n = 0; n < len; ++n) {
      double E_n = f(intensity, values[n].real());
      complex_t rIi{dis(gen) * R, dis(gen) * R};
      values[n] = rIi * E_n;
    }

    if (!fourier_space) {
#ifdef USING_FFTW
      auto plan = symphas::dft::new_fftw_plan<D, complex_t, complex_t>{}(
          values, values, dims, false, true);
      symphas::dft::fftw_execute(plan);
      symphas::dft::fftw_destroy_plan(plan);
#endif
    } else {
      // auto plan = symphas::dft::new_fftw_plan<D, complex_t,
      // complex_t>{}(values, values, dims, false);
      // symphas::dft::fftw_execute(plan);
      // symphas::dft::fftw_destroy_plan(plan);
    }
  }

  const double* dt;
  double h[D];
  F f;
};

inline auto eigen_exponential(scalar_t intensity, scalar_t lambda_n) {
  return std::exp(-intensity * lambda_n / 2);
}

inline auto eigen_polynomial(scalar_t intensity, scalar_t lambda_n) {
  return std::pow(lambda_n, -intensity);
}

inline constexpr auto eigen_one(scalar_t intensity, scalar_t lambda_n) {
  return intensity * lambda_n;
}

template <size_t D, template <typename, size_t> typename grid_type>
struct noise_data<NoiseType::WHITE, D, grid_type>
    : grid_type<scalar_t, D>,
      random_seed<NoiseType::WHITE>,
      random_state<grid_type<scalar_t, D>> {
  using parent_type = grid_type<scalar_t, D>;
  using parent_type::as_grid;
  using parent_type::dims;
  using parent_type::len;
  using parent_type::operator[];
  using seed_type = random_seed<NoiseType::WHITE>;
  using state_type = random_state<grid_type<scalar_t, D>>;
  using state_type::states;

  void allocate() {
    if (parent_type::len == 0) {
      parent_type grid(parent_type::dims);
      using std::swap;
      swap(*static_cast<parent_type*>(this), grid);
    }
    state_type::allocate();
  }

  noise_data(const len_type* dims, const double* h, const double* dt)
      : parent_type(nullptr),
        seed_type(),
        state_type(grid::length<D>(dims)),
        H{1},
        dt{dt},
        gen{get_seed()} {
    for (iter_type i = 0; i < D; ++i) H *= h[i];
    std::copy(dims, dims + D, parent_type::dims);
  }

  noise_data()
      : parent_type(nullptr),
        seed_type(),
        state_type(grid::length<D>(dims)),
        H{1},
        dt{nullptr},
        gen{get_seed()} {}

  auto operator()(double mean, double std_dev) {
    std::normal_distribution<> dis(mean, std_dev);
    return dis(gen);
  }

  template <typename eval_handler_type>
  void update(double intensity, double mean, double std_dev,
              eval_handler_type const& eval_handler) {
    scalar_t sq_corr = intensity * std::sqrt(*dt) / H;
    auto ex =
        expr::make_literal(sq_corr) * expr::callable_of(this)(mean, std_dev);
    eval_handler.result(ex, this->as_grid(), len);
  }

  template <typename R, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<> const& args,
              eval_handler_type const& eval_handler) {
    update(static_cast<R const*>(&intensity)->eval(), NOISE_MEAN, NOISE_STD,
           eval_handler);
  }

  template <typename R, typename T0, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<T0> const& args,
              eval_handler_type const& eval_handler) {
    update(static_cast<R const*>(&intensity)->eval(), NOISE_STD,
           std::get<0>(args), eval_handler);
  }

  template <typename R, typename E, typename T0, typename T1, typename... Ts,
            typename eval_handler_type>
  void update(OpExpression<R> const& intensity,
              std::tuple<T0, T1, Ts...> const& args,
              eval_handler_type const& eval_handler) {
    update(static_cast<R const*>(&intensity)->eval(), std::get<0>(args),
           std::get<1>(args), eval_handler);
  }

  double H;
  const double* dt;
  std::mt19937 gen;
};

template <size_t D, template <typename, size_t> typename grid_type>
struct noise_data<NoiseType::DECAY_EXP, D, grid_type>
    : noise_data_with_function<decltype(&eigen_exponential), D, grid_type>,
      random_seed<NoiseType::DECAY_EXP> {
  using parent_type =
      noise_data_with_function<decltype(&eigen_exponential), D, grid_type>;
  using seed_type = random_seed<NoiseType::DECAY_EXP>;
  using parent_type::operator[];

  void allocate() {}

  noise_data(const len_type* dims, const double* h, const double* dt)
      : parent_type(&eigen_exponential, dims, h, dt), seed_type() {}

  noise_data()
      : parent_type(&eigen_exponential, nullptr, nullptr, nullptr),
        seed_type() {}

  template <typename R, typename eval_handler_type>
  void update(OpExpression<R> const& intensity,
              eval_handler_type const& eval_handler) {
    parent_type::update(*this, static_cast<R const*>(&intensity)->eval(), false,
                        eval_handler);
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    update(*static_cast<R const*>(&intensity), eval_handler);
  }
};

template <size_t D, template <typename, size_t> typename grid_type>
struct noise_data<NoiseType::DECAY_POLY, D, grid_type>
    : noise_data_with_function<decltype(&eigen_polynomial), D, grid_type>,
      random_seed<NoiseType::DECAY_POLY> {
  using parent_type =
      noise_data_with_function<decltype(&eigen_polynomial), D, grid_type>;
  using seed_type = random_seed<NoiseType::DECAY_POLY>;
  using parent_type::operator[];

  void allocate() {}

  noise_data(const len_type* dims, const double* h, const double* dt)
      : parent_type(&eigen_polynomial, dims, h, dt), seed_type() {}

  noise_data()
      : parent_type(&eigen_polynomial, nullptr, nullptr, nullptr, 1.0),
        seed_type() {}

  template <typename R, typename eval_handler_type>
  void update(OpExpression<R> const& intensity,
              eval_handler_type const& eval_handler) {
    parent_type::update(*this, static_cast<R const*>(&intensity)->eval(), false,
                        eval_handler);
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    update(*static_cast<R const*>(&intensity), eval_handler);
  }
};

template <size_t D, template <typename, size_t> typename grid_type>
struct noise_data<NoiseType::NONE, D, grid_type>
    : noise_data_with_function<decltype(&eigen_one), D, grid_type>,
      random_seed<NoiseType::NONE> {
  using parent_type =
      noise_data_with_function<decltype(&eigen_one), D, grid_type>;
  using seed_type = random_seed<NoiseType::NONE>;
  using parent_type::operator[];

  void allocate() {}

  noise_data(const len_type* dims, const double* h, const double* dt)
      : parent_type(&eigen_one, dims, h, dt), seed_type() {}

  noise_data()
      : parent_type(&eigen_one, nullptr, nullptr, nullptr, 1.0), seed_type() {}

  template <typename R, typename eval_handler_type>
  void update(OpExpression<R> const& intensity,
              eval_handler_type const& eval_handler) {
    parent_type::update(*this, static_cast<R const*>(&intensity)->eval(), false,
                        eval_handler);
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    update(*static_cast<R const*>(&intensity), eval_handler);
  }
};

template <size_t D, template <typename, size_t> typename grid_type>
struct noise_data<NoiseType::POISSON, D, grid_type>
    : random_seed<NoiseType::POISSON>, random_state<scalar_t> {
  using seed_type = random_seed<NoiseType::POISSON>;
  using state_type = random_state<scalar_t>;

  void allocate() { state_type::allocate(); }

  noise_data(const len_type* dims, const double* h, const double* dt)
      : seed_type(), state_type(grid::length<D>(dims)), next{0}, value{0} {}

  noise_data(const len_type* dims) : noise_data(dims, nullptr, nullptr) {}

  noise_data() : seed_type(), state_type(), next{0}, value{0} {}

  auto get_index(expr::symbols::Symbol) { return -1; }

  auto get_index(int index) { return index; }

  template <typename R, typename T0, typename T1, typename... Ts,
            typename eval_handler_type>
  void update(OpExpression<R> const& intensity,
              std::tuple<T0, T1, Ts...> const& args,
              eval_handler_type const& eval_handler) {
    auto current = expr::eval(std::get<0>(args));
    if (current >= next) {
      index = get_index(expr::eval(std::get<1>(args)));

      auto gen = get_seed();
      std::exponential_distribution<> dis(
          static_cast<R const*>(&intensity)->eval());
      next += dis(gen);

      std::uniform_real_distribution<> dis0(0.0, 1.0);
      value = dis0(gen);
    }
  }

  auto operator[](iter_type n) const { return next; }

 public:
  double next;
  double value;
};

template <size_t D0, NoiseType nt, size_t D,
          template <typename, size_t> typename grid_type>
struct noise_data_axis;

template <NoiseType nt, size_t D,
          template <typename, size_t> typename grid_type>
struct noise_data_axis<3, nt, D, grid_type>
    : noise_data_axis<2, nt, D, grid_type> {
  using parent_type = noise_data_axis<2, nt, D, grid_type>;
  using parent_type::parent_type;

  void allocate() {
    parent_type::allocate();
    data.allocate();
  }

  noise_data_axis(const len_type* dims, const double* h, const double* dt)
      : parent_type(dims, h, dt), data(dims, h, dt) {}
  noise_data_axis() : parent_type(), data() {}

  vector_t<3> operator[](iter_type n) const {
    return {parent_type::operator[](n)[0], parent_type::operator[](n)[1],
            data[n]};
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    data.update(*static_cast<R const*>(&intensity), args, eval_handler);
    parent_type::update(*static_cast<R const*>(&intensity), args, eval_handler);
  }

  noise_data<nt, D, grid_type> data;
};

template <NoiseType nt, size_t D,
          template <typename, size_t> typename grid_type>
struct noise_data_axis<2, nt, D, grid_type>
    : noise_data_axis<1, nt, D, grid_type> {
  using parent_type = noise_data_axis<1, nt, D, grid_type>;

  void allocate() {
    parent_type::allocate();
    data.allocate();
  }

  noise_data_axis(const len_type* dims, const double* h, const double* dt)
      : parent_type(dims, h, dt), data(dims, h, dt) {}
  noise_data_axis() : parent_type(), data() {}

  vector_t<2> operator[](iter_type n) const {
    return {parent_type::operator[](n)[0], data[n]};
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    data.update(*static_cast<R const*>(&intensity), args, eval_handler);
    parent_type::update(*static_cast<R const*>(&intensity), args, eval_handler);
  }

  noise_data<nt, D, grid_type> data;
};

template <NoiseType nt, size_t D,
          template <typename, size_t> typename grid_type>
struct noise_data_axis<1, nt, D, grid_type> : noise_data<nt, D, grid_type> {
  using parent_type = noise_data<nt, D, grid_type>;

  void allocate() { parent_type::allocate(); }

  noise_data_axis(const len_type* dims, const double* h, const double* dt)
      : parent_type(dims, h, dt) {}
  noise_data_axis() : parent_type() {}

  vector_t<1> operator[](iter_type n) const {
    return {parent_type::operator[](n)};
  }

  template <typename R, typename... Ts, typename eval_handler_type>
  void update(OpExpression<R> const& intensity, std::tuple<Ts...> const& args,
              eval_handler_type const& eval_handler) {
    parent_type::update(*static_cast<R const*>(&intensity), args, eval_handler);
  }
};
}  // namespace expr

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
struct NoiseData : expr::noise_data<nt, D, grid_type> {
  using parent_type = expr::noise_data<nt, D, grid_type>;
  using parent_type::parent_type;
  using parent_type::operator[];
  using parent_type::update;

  NoiseData(expr::noise_data<nt, D, grid_type> const& noise)
      : parent_type(noise) {}
  NoiseData(expr::noise_data<nt, D, grid_type>&& noise) : parent_type(noise) {}
  NoiseData() : parent_type() {}

  template <typename E, typename... Ts>
  auto operator()(OpExpression<E> const& e, Ts&&... args) const {
    auto f = symphas::internal::build_function_for_noise(
        *static_cast<E const*>(&e), std::forward<Ts>(args)...);
    return expr::make_noise(*this, f);
  }

  template <typename T0, typename... Ts>
  auto operator()(T0 const& arg0, Ts&&... args) const {
    auto f = symphas::internal::build_function_for_noise(
        expr::make_literal(arg0), std::forward<Ts>(args)...);
    return expr::make_noise(*this, f);
  }

  auto operator()() const { return expr::make_noise(*this, OpIdentity{}); }
};

template <expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
struct NoiseData<nt, any_vector_t<T, D>, D, grid_type>
    : expr::noise_data_axis<D, nt, D, grid_type> {
  using parent_type = expr::noise_data_axis<D, nt, D, grid_type>;
  using parent_type::allocate;
  using parent_type::parent_type;
  using parent_type::update;
  using parent_type::operator[];

  NoiseData(expr::noise_data_axis<D, nt, D, grid_type> const& noise)
      : parent_type(noise) {}
  NoiseData(expr::noise_data_axis<D, nt, D, grid_type>&& noise)
      : parent_type(noise) {}
  NoiseData() : parent_type() {}

  template <typename E, typename... Ts>
  auto operator()(OpExpression<E> const& e, Ts&&... args) const {
    auto f = symphas::internal::build_function_for_noise(
        *static_cast<E const*>(&e), std::forward<Ts>(args)...);
    return expr::make_noise(*this, f);
  }

  template <typename T0, typename... Ts>
  auto operator()(T0 const& arg0, Ts&&... args) const {
    auto f = symphas::internal::build_function_for_noise(
        expr::make_literal(arg0), std::forward<Ts>(args)...);
    return expr::make_noise(*this, f);
  }
};

//
// DEFINE_BASE_DATA_ARRAY((expr::NoiseType nt, typename T, size_t D),
// (NoiseData<nt, T, D>)) DEFINE_SYMBOL_ID((expr::NoiseType nt, typename T,
// size_t D), (NoiseData<nt, T, D>), { return SYEX_NOISE_OP_STR; })
// DEFINE_SYMBOL_ID((typename T, size_t D), (NoiseData<expr::NoiseType::NONE, T,
// D>), { return SYEX_NOISE_NONE_OP_STR; })

template <typename V, expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          size_t... Ns, typename... Ts>
struct OpSymbolicEval<V, NoiseData<nt, T, D, grid_type>,
                      SymbolicFunction<E, Variable<Ns, Ts>...>>
    : OpExpression<OpSymbolicEval<V, NoiseData<nt, T, D, grid_type>,
                                  SymbolicFunction<E, Variable<Ns, Ts>...>>> {
  using sub_t = NoiseData<nt, T, D, grid_type>;
  using eval_t = SymbolicFunction<E, Variable<Ns, Ts>...>;
  using this_t = OpSymbolicEval<V, sub_t, eval_t>;

  void allocate() { data.allocate(); }

  OpSymbolicEval() = default;

  OpSymbolicEval(V value, sub_t const& data, eval_t const& f)
      : value{value}, f{f}, data{data} {}

  V value;
  eval_t f;
  sub_t data;

 public:
  auto eval(iter_type n) const { return expr::eval(value) * data[n]; }

  auto operator-() const {
    return symphas::internal::make_symbolic_eval(-value, data, f);
  }

  template <typename eval_handler_type, typename... condition_ts>
  void update(eval_handler_type const& eval_handler,
              symphas::lib::types_list<condition_ts...>) {
    data.update(f.e, f.data, eval_handler);
  }

  template <typename eval_handler_type>
  void update(eval_handler_type const& eval_handler) {
    update(eval_handler, symphas::lib::types_list<>{});
  }

#ifdef PRINTABLE_EQUATIONS

  size_t print(FILE* out) const {
    size_t n = 0;
    n += expr::print_with_coeff(out, value);
    n += expr::symbolic_eval_print<sub_t>{}(out, data, (f.e));
    return n;
  }

  size_t print(char* out) const {
    size_t n = 0;
    n += expr::print_with_coeff(out + n, value);
    n += expr::symbolic_eval_print<sub_t>{}(out + n, data, (f.e));
    return n;
  }

  size_t print_length() const {
    return expr::coeff_print_length(value) +
           expr::symbolic_eval_print<sub_t>{}(data, (f.e));
  }

#endif
};

template <typename V, expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          size_t... Ns, typename... Ts>
OpSymbolicEval(V, NoiseData<nt, T, D, grid_type>,
               SymbolicFunction<E, Variable<Ns, Ts>...>)
    -> OpSymbolicEval<V, NoiseData<nt, T, D, grid_type>,
                      SymbolicFunction<E, Variable<Ns, Ts>...>>;

namespace expr {

template <typename T, size_t D, template <typename, size_t> typename grid_type,
          typename E, typename... Ts>
auto make_poisson_event(
    NoiseData<expr::NoiseType::POISSON, T, D, grid_type> const& noise,
    SymbolicFunction<E, Ts...> const& f) {
  return symphas::internal::make_symbolic_eval(OpIdentity{}, noise, f);
}

template <typename T, size_t D, template <typename, size_t> typename grid_type,
          typename E, typename L>
auto make_poisson_event(
    NoiseData<expr::NoiseType::POISSON, T, D, grid_type> const& noise,
    OpExpression<E> const& e, OpExpression<L> const& lambda) {
  return make_poisson_event(noise, function_of(*static_cast<L>(&lambda)) =
                                       *static_cast<E const*>(&e));
}

template <NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename... Ts>
auto make_noise(NoiseData<nt, T, D, grid_type> const& noise,
                SymbolicFunction<E, Ts...> const& f) {
  return symphas::internal::make_symbolic_eval(OpIdentity{}, noise, f);
}

template <NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename E>
auto make_noise(NoiseData<nt, T, D, grid_type> const& noise,
                OpExpression<E> const& e) {
  return make_noise(noise, function_of() = *static_cast<E const*>(&e));
}

template <NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
auto make_noise(const len_type* dimensions, const double* h, const double* dt) {
  return make_noise(
      NoiseData<nt, T, D, grid_type>(dimensions, h, const_cast<double*>(dt)),
      OpIdentity{});
}

template <typename T, size_t D, template <typename, size_t> typename grid_type>
auto make_white_noise(const len_type* dimensions, const double* h,
                      const double* dt, double intensity = 1.0) {
  return make_noise(NoiseData<NoiseType::WHITE, T, D, grid_type>(dimensions, h,
                                                                 dt, intensity),
                    OpVoid{});
}

template <Axis ax, NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
auto resolve_axis_component(
    NoiseData<nt, any_vector_t<T, D>, D, grid_type> const& data) {
  if constexpr (ax == Axis::X) {
    return as_component_data<ax, D>(
        static_cast<noise_data<nt, D, grid_type> const*>(&data)->values);
  } else if constexpr (ax == Axis::Y) {
    return as_component_data<ax, D>(
        static_cast<noise_data_axis<2, nt, D, grid_type> const*>(&data)
            ->data.values);
  } else if constexpr (ax == Axis::Z) {
    return as_component_data<ax, D>(
        static_cast<noise_data_axis<3, nt, D, grid_type> const*>(&data)
            ->data.values);
  } else {
    return as_component_data<ax, D>((double*)nullptr);
  }
}

template <Axis ax, NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type>
auto resolve_axis_component(
    NoiseData<nt, any_vector_t<T, D>, D, grid_type>& data) {
  if constexpr (ax == Axis::X) {
    return as_component_data<ax, D>(
        static_cast<noise_data<nt, D, grid_type>*>(&data)->values);
  } else if constexpr (ax == Axis::Y) {
    return as_component_data<ax, D>(
        static_cast<noise_data_axis<2, nt, D, grid_type>*>(&data)->data.values);
  } else if constexpr (ax == Axis::Z) {
    return as_component_data<ax, D>(
        static_cast<noise_data_axis<3, nt, D, grid_type>*>(&data)->data.values);
  } else {
    return as_component_data<ax, D>((double*)nullptr);
  }
}

}  // namespace expr

namespace expr::transform {

template <size_t D, NoiseType nt, typename T,
          template <typename, size_t> typename grid_type>
auto to_ft(NoiseData<nt, T, D, grid_type> const& e, double const* h,
           const len_type* dims) {
  return e;
}

template <size_t D, typename V, typename T,
          template <typename, size_t> typename grid_type, typename E,
          size_t... Ns, typename... Ts>
auto to_ft(OpSymbolicEval<V, NoiseData<NoiseType::WHITE, T, D, grid_type>,
                          SymbolicFunction<E, Variable<Ns, Ts>...>> const& e,
           double const* h, const len_type* dims) {
  auto noise = to_ft<D>(e.data, h, dims);
  return expr::coeff(e) * expr::make_noise(noise, e.f);
}
}  // namespace expr::transform
