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
 * PURPOSE: Implements the convolution operation. Represents taking the
 * convolution between two expressions, with specializations based on
 * when one of the terms is a variable or smoothing kernel.
 *
 * ***************************************************************************
 */

#pragma once

#include "convolutionlib.h"

// Forward declarations for CUDA implementations
#ifdef USING_CUDA
namespace symphas::internal {

template <typename V, typename E1, typename E2, size_t D0, typename G_T>
void update_impl_cuda(expr::ConvolutionDataPairCUDA<D0>& compute_data, 
                      Grid<G_T, D0>& result_grid);

template <typename V, size_t D, template <typename, size_t> typename grid_type, 
          typename E, size_t D0, typename G_T>
void update_impl_cuda(expr::ConvolutionDataCUDA<D0>& compute_data, 
                      grid_type<G_T, D>& result_grid,
                      const GaussianSmoothing<D, grid_type>& smoother);

} // namespace symphas::internal
#endif

//! \cond

#ifdef LATEX_PLOT
#define SYEX_CONVOLUTION_FMT_A "\\left("
#define SYEX_CONVOLUTION_FMT_B "\\right)"
#define SYEX_CONVOLUTION_FMT_SEP "*"
#else
#define SYEX_CONVOLUTION_FMT_A "("
#define SYEX_CONVOLUTION_FMT_B ")"
#define SYEX_CONVOLUTION_FMT_SEP " x "
#endif

#define SYEX_CONVOLUTION_FMT                           \
  SYEX_CONVOLUTION_FMT_A "%s" SYEX_CONVOLUTION_FMT_SEP \
                         "%s" SYEX_CONVOLUTION_FMT_B
#define SYEX_CONVOLUTION_FMT_LEN                             \
  (STR_ARR_LEN(SYEX_CONVOLUTION_FMT_A SYEX_CONVOLUTION_FMT_B \
                   SYEX_CONVOLUTION_FMT_SEP) -               \
   1)

//! \endcond

namespace symphas::internal {
//! Implementation of convolution expression construction.
/*!
 * Implementation of functions which generalize and simplify the way to
 * construct convolution expressions. Wraps the template deduction necessary
 * to initialize a convolution expression.
 */
struct make_convolution {
  //! Constructs the convolution with the identity coefficient.
  template <typename A, typename B>
  static auto get(A&&, B&&);

  //! Constructs the convolution applied to the given expression.
  template <typename V, typename E2>
  static auto get(V const&, OpVoid, OpExpression<E2> const&) {
    return OpVoid{};
  }

  //! Constructs the convolution applied to the given expression.
  template <typename V, typename E1>
  static auto get(V const&, OpExpression<E1> const&, OpVoid) {
    return OpVoid{};
  }

  //! Constructs the convolution applied to the given expression.
  template <typename V, typename E1, typename E2>
  static auto get(V const&, OpExpression<E1> const&, OpExpression<E2> const&);

  //! Constructs the convolution applied to the given expression.
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type, typename E>
  static auto get(V const&, GaussianSmoothing<D, grid_type> const&,
                  OpExpression<E> const&);

  //! Constructs the convolution applied to the given expression.
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type>
  static auto get(V const&, GaussianSmoothing<D, grid_type> const&, OpVoid);

  //! Constructs the convolution applied to the given expression.
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type, typename E>
  static auto get(V const& value, OpExpression<E> const& e,
                  GaussianSmoothing<D, grid_type> const& smoother) {
    return get(value, smoother, *static_cast<E const*>(&e));
  }

  //! Constructs the convolution between a variable and Gaussian kernel.
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type, typename S,
            typename G>
  static auto get(V const&, GaussianSmoothing<D, grid_type> const&,
                  OpTerm<S, G> const&);

  //! Constructs the convolution between a variable and Gaussian kernel.
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type, typename S,
            typename G>
  static auto get(V const& value, OpTerm<S, G> const& e,
                  GaussianSmoothing<D, grid_type> const& smoother) {
    return get(value, smoother, e);
  }

  //! Constructs the convolution using a grid instead of an expression.
  /*!
   * Used for the convolution specialization for the OpTerm.
   */
  template <typename V, size_t D,
            template <typename, size_t> typename grid_type, typename G>
  static auto get_g(V const&, GaussianSmoothing<D, grid_type> const&, G g);
};
}  // namespace symphas::internal

namespace expr {
//! Create a convolution expression.
/*!
 * Create a convolution expression.
 *
 * \param a The left hand side of the convolution operation.
 * \param b The right hand side of the convolution operation.
 */
template <typename A, typename B>
auto make_convolution(A&& a, B&& b) {
  return symphas::internal::make_convolution::get(std::forward<A>(a),
                                                  std::forward<B>(b));
}

//! Create a convolution expression.
/*!
 * Create a convolution expression.
 *
 * \param v The coefficient of the convolution expression.
 * \param a The left hand side of the convolution operation.
 * \param b The right hand side of the convolution operation.
 */
template <typename V, typename A, typename B>
auto make_convolution(V&& v, A&& a, B&& b) {
  return symphas::internal::make_convolution::get(
      std::forward<V>(v), std::forward<A>(a), std::forward<B>(b));
}

}  // namespace expr

namespace expr::internal {
template <typename G>
struct select_convolution_data_type {
 protected:
  template <typename T>
  struct wrap_type {
    using type = T;
  };

  template <typename T, size_t D>
  auto _get(Grid<T, D>) {
    return wrap_type<ConvolutionData<D>>{};
  }

#ifdef USING_CUDA
  template <typename T, size_t D>
  auto _get(GridCUDA<T, D>) {
    return wrap_type<ConvolutionDataCUDA<D>>{};
  }
#endif

  auto get(G const& g) { return _get(g); }

 public:
  using type = typename std::invoke_result_t<
      decltype(&select_convolution_data_type::get),
      select_convolution_data_type, G>::type;
};

template <typename G>
using convolution_data_type = typename select_convolution_data_type<G>::type;

template <typename G>
struct select_convolution_data_pair_type {
 protected:
  template <typename T>
  struct wrap_type {
    using type = T;
  };

  template <typename T, size_t D>
  auto _get(Grid<T, D>) {
    return wrap_type<ConvolutionDataPair<D>>{};
  }

#ifdef USING_CUDA
  template <typename T, size_t D>
  auto _get(GridCUDA<T, D>) {
    return wrap_type<ConvolutionDataPairCUDA<D>>{};
  }
#endif

  auto get(G const& g) { return _get(g); }

 public:
  using type = typename std::invoke_result_t<
      decltype(&select_convolution_data_pair_type::get),
      select_convolution_data_pair_type, G>::type;
};

template <typename G>
using convolution_data_pair_type =
    typename select_convolution_data_pair_type<G>::type;

template <typename G>
struct select_block_type;

template <typename T, size_t D>
struct select_block_type<Grid<T, D>> {
  using type = Block<T>;
};

#ifdef USING_CUDA
template <typename T, size_t D>
struct select_block_type<GridCUDA<T, D>> {
  using type = BlockCUDA<T>;
};
#endif

template <typename G>
using block_type = typename select_block_type<G>::type;
}  // namespace expr::internal

//! Convolution of two arbitrary expressions.
/*!
 * An implementation of the convolution operator; the convolution takes two
 * functions and upon evaluating, computes the convolution.
 * In order to evaluate the convolution, the update function must be called in
 * the pruning step.
 *
 * In order to compute the convolution, the convolution operator uses the
 * convolution theorem, evaluating the terms in Fourier space before performing
 * an inverse Fourier transform back into real space.
 *
 * \tparam V Type of the coefficient.
 * \tparam E1 Left hand expression type.
 * \tparam E2 Right hand expression type.
 */
template <typename V, typename E1, typename E2>
struct OpConvolution : OpExpression<OpConvolution<V, E1, E2>> {
  static const size_t e1_U = expr::grid_dim<E1>::dimension;
  static const size_t e2_U = expr::grid_dim<E2>::dimension;

  static const size_t D = fixed_max<e1_U, e1_U>;

  // identifying the type that comes out as a result of the derivative
  // get the type of the system holding the intermediate results; based on
  // the multiplication of the 'real' space terms
  using e1_T = typename expr::eval_type<E1>::type;
  using e2_T = typename expr::eval_type<E2>::type;
  using G_T = mul_result_t<e1_T, e2_T>;

  using storage_type_1 = expr::storage_type_t<E1>;
  using storage_type_2 = expr::storage_type_t<E2>;

  template <typename T, size_t D0>
  using grid_type = typename expr::enclosing_parent_storage_t<
      storage_type_1>::template type<T, D0>;
  using result_type = grid_type<G_T, D>;
  using convolution_type =
      typename expr::internal::convolution_data_pair_type<result_type>;

  using data_a_t = storage_type_1;
  using data_b_t = storage_type_2;

  void allocate() {
    b.allocate();
    a.allocate();

    if (g0.len == 0) {
      g0 = result_type(expr::data_dimensions(a, b));
      data_a = data_a_t(expr::data_dimensions(a, b));
      data_b = data_b_t(expr::data_dimensions(a, b));
      compute = convolution_type(data_a.values, data_b.values, g0);
    }
  }

  OpConvolution()
      : data_a{0}, data_b{0}, g0{0}, value{V{}}, a{}, b{}, compute{} {}

  //! Generate the convolution expression.
  /*!
   * The convolution expression is generated using the given value as the
   * coefficient and the two terms, with the order representing the
   * side they appear on relative to the convolution operator.
   *
   * \param value The coefficient value of the convolution term.
   * \param a The left hand side of the convolution operator.
   * \param b The right hand side of the convolution operator.
   */
  OpConvolution(V value, E1 const& a, E2 const& b)
      : data_a{0},
        data_b{0},
        g0{0},
        value{value},
        a{a},
        b{b},
        compute{} { /*update();*/ }

  OpConvolution(OpConvolution<V, E1, E2> const& other)
      : data_a{0},
        data_b{0},
        g0{0},
        value{other.value},
        a{other.a},
        b{other.b},
        compute{} { /*update();*/
    ;
  }

  OpConvolution(OpConvolution<V, E1, E2>&& other) noexcept : OpConvolution() {
    swap(*this, other);
  }

  friend void swap(OpConvolution<V, E1, E2>& first,
                   OpConvolution<V, E1, E2>& second) {
    using std::swap;
    swap(first.data_a, second.data_a);
    swap(first.data_b, second.data_b);
    swap(first.g0, second.g0);
    swap(first.value, second.value);
    swap(first.a, second.a);
    swap(first.b, second.b);
    swap(first.compute, second.compute);
  }

  inline auto eval(iter_type n) const { return expr::eval(value) * g0[n]; }

  auto operator-() const {
    return symphas::internal::make_convolution::get(-value, a, b);
  }

  // template<typename V1,
  //	typename std::enable_if_t<expr::is_combinable<E1>, int> = 0>
  // auto operator+(OpConvolution<V1, E1, E2> const& other) const
  //{
  //	return symphas::internal::make_convolution::get(value + other.value, a,
  // b);
  // }

  // template<typename V1>
  // auto operator+(OpConvolution<V1, E2, E1> const& other) const
  //{
  //	return symphas::internal::make_convolution::get(value + other.value, a,
  // b);
  // }

  // template<typename V1>
  // auto operator-(OpConvolution<V1, E1, E2> const& other) const
  //{
  //	return symphas::internal::make_convolution::get(value - other.value, a,
  // b);
  // }

  // template<typename V1>
  // auto operator-(OpConvolution<V1, E2, E1> const& other) const
  //{
  //	return symphas::internal::make_convolution::get(value - other.value, a,
  // b);
  // }

#ifdef PRINTABLE_EQUATIONS

  size_t print(FILE* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
    n += a.print(out);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
    n += b.print(out);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print(char* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
    n += a.print(out + n);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
    n += b.print(out + n);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print_length() const {
    return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN +
           a.print_length() + b.print_length();
  }

#endif

 protected:
  data_a_t data_a;  //!< Data of the result of the first expression.
  data_b_t data_b;  //!< Data of the result of the second expression.
  result_type g0;   //!< Grid storing the final result of the convolution.

 public:
  V value;  //!< Coefficient of this convolution expression.
  E1 a;     //!< First expression in the convolution.
  E2 b;     //!< Second expression in the convolution.

  template <typename V0, typename E10, typename E20>
  friend auto const& expr::get_result_data(OpConvolution<V0, E10, E20> const&);
  template <typename V0, typename E10, typename E20>
  friend auto& expr::get_result_data(OpConvolution<V0, E10, E20>&);

  //! Update the convolution by computing the result into the stored grid.
  /*!
   * Update the convolution by computing the result into the stored grid.
   */
  template <typename eval_handler_type, typename... condition_ts>
  void update(eval_handler_type const& eval_handler,
              symphas::lib::types_list<condition_ts...>) {
    // compute the result of the expressions and update the grids
    eval_handler.result(a, data_a.values, data_a.len);
    eval_handler.result(b, data_b.values, data_b.len);

    compute.transform_in_out_0(&expr::BaseData<data_a_t>::get(data_a)[0],
                               g0.dims);
    compute.transform_in_out_1(&expr::BaseData<data_a_t>::get(data_b)[0],
                               g0.dims);

    // Dispatch to CPU or CUDA implementation based on convolution type
    update_impl(compute, g0);
  }

 private:
  //! CPU implementation of convolution update
  template <size_t D0>
  void update_impl(expr::ConvolutionDataPair<D0>& compute_data, result_type& result_grid) {
    len_type len = symphas::dft::length<G_T, D>(result_grid.dims);
    if constexpr (std::is_same<G_T, complex_t>::value) {
      grid::scale(compute_data.out_0, len);
      grid::scale(compute_data.out_1, len);
    }

#pragma omp parallel for
    for (iter_type i = 0; i < len; ++i) {
      symphas::math::real(compute_data.in_2[i]) =
          symphas::math::real(compute_data.out_0[i]) *
              symphas::math::real(compute_data.out_1[i]) -
          symphas::math::imag(compute_data.out_0[i]) *
              symphas::math::imag(compute_data.out_1[i]);
      symphas::math::imag(compute_data.in_2[i]) =
          symphas::math::real(compute_data.out_0[i]) *
              symphas::math::imag(compute_data.out_1[i]) +
          symphas::math::imag(compute_data.out_0[i]) *
              symphas::math::real(compute_data.out_1[i]);
    }

    compute_data.transform_out_in(result_grid.values, result_grid.dims);

    if constexpr (std::is_same<G_T, scalar_t>::value) {
      grid::scale(result_grid);
    }
  }

#ifdef USING_CUDA
  //! CUDA implementation of convolution update
  template <size_t D0>
  void update_impl(expr::ConvolutionDataPairCUDA<D0>& compute_data, result_type& result_grid) {
    symphas::internal::update_impl_cuda<V, E1, E2, D0, G_T>(compute_data, result_grid);
  }
#endif

 public:

  template <typename eval_handler_type>
  void update(eval_handler_type const& eval_handler) {
    update(eval_handler, symphas::lib::types_list<>{});
  }

 protected:
  convolution_type compute;
};

template <typename coeff_t, typename V2, typename E1, typename E2,
          typename std::enable_if_t<(expr::is_coeff<coeff_t> &&
                                     !expr::is_tensor<coeff_t> &&
                                     !expr::is_tensor<V2>),
                                    int> = 0>
auto operator*(coeff_t const& value, OpConvolution<V2, E1, E2> const& b) {
  return symphas::internal::make_convolution::get(value * b.value, b.a, b.b);
}

template <typename coeff_t, typename tensor_t, typename E1, typename E2,
          typename std::enable_if_t<
              (expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpConvolution<tensor_t, E1, E2> const& b) {
  return (value * b.value) *
         symphas::internal::make_convolution::get(OpIdentity{}, b.a, b.b);
}

template <
    typename tensor_t, typename V, typename E1, typename E2,
    typename std::enable_if_t<
        (expr::is_tensor<tensor_t> && !expr::is_tensor<V> &&
         expr::eval_type<E1>::rank == 0 && expr::eval_type<E2>::rank == 0),
        int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b) {
  return symphas::internal::make_convolution::get(tensor * b.value, b.a, b.b);
}

template <typename tensor_t, typename V, typename E1, typename E2,
          typename std::enable_if_t<
              (expr::is_tensor<tensor_t> && !expr::is_tensor<V> &&
               expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank == 0),
              int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b) {
  return symphas::internal::make_convolution::get(b.value, tensor * b.a, b.b);
}

template <typename tensor_t, typename V, typename E1, typename E2,
          typename std::enable_if_t<
              (expr::is_tensor<tensor_t> && !expr::is_tensor<V> &&
               expr::eval_type<E1>::rank == 0 && expr::eval_type<E2>::rank > 0),
              int> = 0>
auto operator*(tensor_t const& tensor, OpConvolution<V, E1, E2> const& b) {
  return symphas::internal::make_convolution::get(b.value, b.a, tensor * b.b);
}

// ******************************************************************************************************************

//! Convolution of an expression with a Gaussian smoothing kernel.
/*!
 * Specialization of the convolution of two arbitrary expressions.
 *
 * \tparam V Type of the coefficient.
 * \tparam D the dimension of the Gaussian smoothing kernel.
 * \tparam E The type of the expression convoluted with the Gaussian smoother.
 */
template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename E>
struct OpConvolution<V, GaussianSmoothing<D, grid_type>, E>
    : OpExpression<OpConvolution<V, GaussianSmoothing<D, grid_type>, E>> {
  using G_T = typename expr::eval_type<E>::type;
  using result_type = grid_type<G_T, D>;
  using convolution_type = expr::internal::convolution_data_type<result_type>;
  using data_type = result_type;

  void allocate() {
    e.allocate();

    if (g0.len == 0) {
      smoother.allocate();
      g0 = result_type(expr::data_dimensions(smoother));
      data = data_type(expr::data_dimensions(smoother));

      using std::swap;
      convolution_type allocated_compute(data.values, g0);
      swap(compute, allocated_compute);
    }
  }

  OpConvolution()
      : g0{0},
        data{0},
        value{V{}},
        smoother{GaussianSmoothing<D, grid_type>()},
        compute{} {}

  //! Generate the convolution expression.
  /*!
   * The convolution expression is generated using the given value as the
   * coefficient and the two terms, with the order representing the
   * side they appear on relative to the convolution operator.
   *
   * \param value The coefficient value of the convolution term.
   * \param smoother The Gaussian smoothing kernel.
   * \param e The expression that is smoothed with this expression.
   */
  OpConvolution(V value, GaussianSmoothing<D, grid_type> const& smoother,
                E const& e)
      : g0{0},
        data{0},
        value{value},
        e{e},
        smoother{smoother},
        compute{} { /*update();*/ }

  OpConvolution(
      OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& other)
      : g0{0},
        data{0},
        value{other.value},
        e{other.e},
        smoother{other.smoother},
        compute{} { /*update();*/ }

  OpConvolution(
      OpConvolution<V, GaussianSmoothing<D, grid_type>, E>&& other) noexcept
      : OpConvolution() {
    swap(*this, other);
  }

  friend void swap(
      OpConvolution<V, GaussianSmoothing<D, grid_type>, E>& first,
      OpConvolution<V, GaussianSmoothing<D, grid_type>, E>& second) {
    using std::swap;
    swap(first.g0, second.g0);
    swap(first.data, second.data);
    swap(first.value, second.value);
    swap(first.e, second.e);
    swap(first.smoother, second.smoother);
    swap(first.compute, second.compute);
  }

  inline auto eval(iter_type n) const { return expr::eval(value) * g0[n]; }

  auto operator-() const {
    return symphas::internal::make_convolution::get(-value, smoother, e);
  }

  template <typename V1>
  auto operator+(OpConvolution<V1, GaussianSmoothing<D, grid_type>, E> const&
                     other) const {
    return symphas::internal::make_convolution::get(value + other.value,
                                                    smoother, e);
  }

  template <typename V1>
  auto operator-(OpConvolution<V1, GaussianSmoothing<D, grid_type>, E> const&
                     other) const {
    return symphas::internal::make_convolution::get(value - other.value,
                                                    smoother, e);
  }

#ifdef PRINTABLE_EQUATIONS

  size_t print(FILE* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
    n += smoother.print(out);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
    n += e.print(out);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print(char* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
    n += smoother.print(out + n);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
    n += e.print(out + n);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print_length() const {
    return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN +
           e.print_length() + smoother.print_length();
  }

#endif

 protected:
  result_type g0;  //!< The grid storing the result of this convolution.
  data_type data;  //!< Grid storing the result of the second expression.

 public:
  V value;  //!< Value multiplying the result of this convolution.
  E e;      //!< The expression to populate the first grid.
  GaussianSmoothing<D, grid_type> smoother;  //!< The smoothing kernel.

  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename E0>
  friend auto const& expr::get_enclosed_expression(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>, E0> const&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename E0>
  friend auto& expr::get_enclosed_expression(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>, E0>&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename E0>
  friend auto const& expr::get_result_data(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>, E0> const&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename E0>
  friend auto& expr::get_result_data(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>, E0>&);

  //! Update the convolution by computing the result into the stored grid.
  /*!
   * Update the convolution by computing the result into the stored grid.
   */
  template <typename eval_handler_type, typename... condition_ts>
  void update(eval_handler_type const& eval_handler,
              symphas::lib::types_list<condition_ts...>) {
    eval_handler.result(e, data, data.len);
    compute.transform_in_out(&expr::BaseData<data_type>::get(data)[0],
                             g0.dims);

    // Dispatch to CPU or CUDA implementation
    update_impl(compute, g0);
  }

 private:
  //! CPU implementation of Gaussian smoothing convolution
  template <size_t D0>
  void update_impl(expr::ConvolutionData<D0>& compute_data, result_type& result_grid) {
    auto f = [&](iter_type i, iter_type ft_i) {
      symphas::math::real(compute_data.in_1[ft_i]) =
          smoother.eval(i) * symphas::math::real(compute_data.out_0[ft_i]);
      symphas::math::imag(compute_data.in_1[ft_i]) =
          smoother.eval(i) * symphas::math::imag(compute_data.out_0[ft_i]);
    };

    symphas::dft::iterate_rc<G_T, D>(f, result_grid.dims);
    compute_data.transform_out_in(result_grid.values, result_grid.dims);
    grid::scale(result_grid);
  }

#ifdef USING_CUDA
  //! CUDA implementation of Gaussian smoothing convolution
  template <size_t D0>
  void update_impl(expr::ConvolutionDataCUDA<D0>& compute_data, result_type& result_grid) {
    symphas::internal::update_impl_cuda<V, D, grid_type, E, D0, G_T>(compute_data, result_grid, smoother);
  }
#endif

 public:

  template <typename eval_handler_type>
  void update(eval_handler_type const& eval_handler) {
    update(eval_handler, symphas::lib::types_list<>{});
  }

 protected:
  convolution_type compute;
  // expr::ConvolutionData<D> compute;
};

template <typename coeff_t, typename V, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<(expr::is_coeff<coeff_t> &&
                                     !expr::is_tensor<coeff_t> &&
                                     !expr::is_tensor<V>),
                                    int> = 0>
auto operator*(coeff_t const& value,
               OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& b) {
  return symphas::internal::make_convolution::get(
      value * b.value, b.smoother, expr::get_enclosed_expression(b));
}

template <typename coeff_t, typename tensor_t, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<
              (expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(
    coeff_t const& value,
    OpConvolution<tensor_t, GaussianSmoothing<D, grid_type>, E> const& b) {
  return (value * b.value) *
         symphas::internal::make_convolution::get(
             OpIdentity{}, b.smoother, expr::get_enclosed_expression(b));
}

template <typename tensor_t, typename V, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<(expr::is_tensor<tensor_t> &&
                                     !expr::is_tensor<V> &&
                                     expr::eval_type<E>::rank == 0),
                                    int> = 0>
auto operator*(tensor_t const& tensor,
               OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& b) {
  return symphas::internal::make_convolution::get(
      tensor * b.value, b.smoother, expr::get_enclosed_expression(b));
}

template <typename tensor_t, typename V, size_t D,
          template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<(expr::is_tensor<tensor_t> &&
                                     !expr::is_tensor<V> &&
                                     expr::eval_type<E>::rank > 0),
                                    int> = 0>
auto operator*(tensor_t const& tensor,
               OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& b) {
  return symphas::internal::make_convolution::get(
      b.value, b.smoother, tensor * expr::get_enclosed_expression(b));
}

// ******************************************************************************************************************

//! Convolution of an OpTerm with a Gaussian smoothing kernel.
/*!
 * Specialization of the convolution of two arbitrary expressions.
 *
 * \tparam V Type of the coefficient.
 * \tparam D the dimension of the Gaussian smoothing kernel.
 * \tparam grid_type The grid type template used by the Gaussian kernel.
 * \tparam G The type of variable data.
 */
template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename G>
struct OpConvolution<V, GaussianSmoothing<D, grid_type>, OpTerm<OpIdentity, G>>
    : OpExpression<OpConvolution<V, GaussianSmoothing<D, grid_type>,
                                 OpTerm<OpIdentity, G>>> {
  using E = OpTerm<V, G>;
  using G_T = expr::eval_type_t<E>;
  using result_type = grid_type<G_T, D>;
  using convolution_type = typename expr::internal::convolution_data_type<result_type>;

  void allocate() {
    if (g0.len == 0) {
      smoother.allocate();
      g0 = result_type(expr::data_dimensions(smoother));

      using std::swap;
      convolution_type allocated_compute(expr::BaseData<G>::get(data), g0);
      swap(compute, allocated_compute);
    }
  }

  OpConvolution() : g0{0}, data{0}, value{V{}}, smoother{} {}

  //! Generate the convolution expression.
  /*!
   * The convolution expression is generated using the given value as the
   * coefficient and the two terms, with the order representing the
   * side they appear on relative to the convolution operator.
   *
   * \param value The coefficient value of the convolution term.
   * \param smoother The Gaussian smoothing kernel.
   * \param a The variable which is smoothed.
   */
  template <typename V0, typename V1,
            typename std::enable_if_t<
                std::is_convertible<mul_result_t<V0, V1>, V>::value, int> = 0>
  OpConvolution(V0 value, GaussianSmoothing<D, grid_type> const& smoother,
                OpTerm<V1, G> const& a)
      : g0{0},
        data{expr::data(a)},
        value{value * expr::coeff(a)},
        smoother{smoother},
        compute{} { /*update();*/
    ;
  }
  OpConvolution(V value, GaussianSmoothing<D, grid_type> const& smoother,
                G grid)
      : g0{0},
        data{grid},
        value{value},
        smoother{smoother},
        compute{} { /*update();*/ }

  OpConvolution(OpConvolution<V, GaussianSmoothing<D, grid_type>,
                              OpTerm<OpIdentity, G>> const& other)
      : g0{0},
        data{other.data},
        value{other.value},
        smoother{other.smoother},
        compute{} { /*update();*/ }

  OpConvolution(OpConvolution<V, GaussianSmoothing<D, grid_type>,
                              OpTerm<OpIdentity, G>>&& other) noexcept
      : OpConvolution() {
    swap(*this, other);
  }

  friend void swap(OpConvolution<V, GaussianSmoothing<D, grid_type>,
                                 OpTerm<OpIdentity, G>>& first,
                   OpConvolution<V, GaussianSmoothing<D, grid_type>,
                                 OpTerm<OpIdentity, G>>& second) {
    using std::swap;
    swap(first.g0, second.g0);
    swap(first.data, second.data);
    swap(first.value, second.value);
    swap(first.smoother, second.smoother);
    swap(first.compute, second.compute);
  }

  inline auto eval(iter_type n) const { return expr::eval(value) * g0[n]; }

  auto operator-() const {
    return symphas::internal::make_convolution::get_g(-value, smoother, data);
  }

  template <typename V1>
  auto operator+(OpConvolution<V1, GaussianSmoothing<D, grid_type>, E> const&
                     other) const {
    return symphas::internal::make_convolution::get_g(value + other.value,
                                                      smoother, data);
  }

  template <typename V1>
  auto operator-(OpConvolution<V1, GaussianSmoothing<D, grid_type>, E> const&
                     other) const {
    return symphas::internal::make_convolution::get_g(value - other.value,
                                                      smoother, data);
  }

#ifdef PRINTABLE_EQUATIONS

  size_t print(FILE* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_A);
    n += smoother.print(out);
    n += fprintf(out, SYEX_CONVOLUTION_FMT_SEP);
    n += fprintf(out, "%s", expr::get_op_name(data));
    n += fprintf(out, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print(char* out) const {
    size_t n = expr::print_with_coeff(out, value);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_A);
    n += smoother.print(out + n);
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_SEP);
    n += sprintf(out + n, "%s", expr::get_op_name(data));
    n += sprintf(out + n, SYEX_CONVOLUTION_FMT_B);
    return n;
  }

  size_t print_length() const {
    return expr::coeff_print_length(value) + SYEX_CONVOLUTION_FMT_LEN +
           smoother.print_length() + std::strlen(expr::get_op_name(data));
  }

#endif

 protected:
  result_type g0;  //!< Grid storing the result of the convolution.
  G data;  //!< The data from the OpTerm.

 public:
  V value;  //!< Value multiplying the result of this convolution.
  GaussianSmoothing<D, grid_type> smoother;  //!< The smoothing kernel.

  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename G0>
  friend auto expr::get_enclosed_expression(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>,
                    OpTerm<OpIdentity, G0>> const&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename G0>
  friend auto expr::get_enclosed_expression(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>,
                    OpTerm<OpIdentity, G0>>&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename G0>
  friend auto const& expr::get_result_data(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>,
                    OpTerm<OpIdentity, G0>> const&);
  template <typename V0, size_t D0,
            template <typename, size_t> typename grid_type0, typename G0>
  friend auto& expr::get_result_data(
      OpConvolution<V0, GaussianSmoothing<D0, grid_type0>,
                    OpTerm<OpIdentity, G0>>&);

  template <typename eval_handler_type, typename... condition_ts>
  void update(eval_handler_type const& eval_handler,
              symphas::lib::types_list<condition_ts...>) {
    compute.transform_in_out(&expr::BaseData<G>::get(data)[0], g0.dims);

    // Dispatch to CPU or CUDA implementation
    update_impl(compute, g0);
  }

 private:
  //! CPU implementation of OpTerm Gaussian smoothing convolution
  template <size_t D0>
  void update_impl(expr::ConvolutionData<D0>& compute_data, result_type& result_grid) {
    using G_T = typename expr::eval_type<OpTerm<OpIdentity, G>>::type;
    auto f = [&](iter_type ft_i, iter_type i) {
      symphas::math::real(compute_data.in_1[i]) =
          smoother.eval(ft_i) * symphas::math::real(compute_data.out_0[i]);
      symphas::math::imag(compute_data.in_1[i]) =
          smoother.eval(ft_i) * symphas::math::imag(compute_data.out_0[i]);
    };

    symphas::dft::iterate_rc<G_T, D>(f, result_grid.dims);
    compute_data.transform_out_in(result_grid.values, result_grid.dims);
    grid::scale(result_grid);
  }

#ifdef USING_CUDA
  //! CUDA implementation of OpTerm Gaussian smoothing convolution
  template <size_t D0>
  void update_impl(expr::ConvolutionDataCUDA<D0>& compute_data, result_type& result_grid) {
    symphas::internal::update_impl_cuda<V, D, grid_type, G, D0>(compute_data, result_grid, smoother);
  }
#endif

 public:

  template <typename eval_handler_type>
  void update(eval_handler_type const& eval_handler) {
    update(eval_handler, symphas::lib::types_list<>{});
  }

 protected:
  convolution_type compute;
};

template <typename V0, size_t D, template <typename, size_t> typename grid_type,
          typename T, typename G>
OpConvolution(V0, GaussianSmoothing<D, grid_type>, OpTerm<T, G>)
    -> OpConvolution<mul_result_t<V0, T>, GaussianSmoothing<D, grid_type>,
                     OpTerm<OpIdentity, G>>;
template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename E>
OpConvolution(V, GaussianSmoothing<D, grid_type>, OpExpression<E>)
    -> OpConvolution<V, GaussianSmoothing<D, grid_type>, E>;
template <typename V, typename E1, typename E2>
OpConvolution(V, OpExpression<E1>, OpExpression<E2>)
    -> OpConvolution<V, E1, E2>;

template <typename coeff_t, typename V2, size_t D,
          template <typename, size_t> typename grid_type, typename G2,
          typename std::enable_if_t<(expr::is_coeff<coeff_t> &&
                                     !expr::is_tensor<coeff_t> &&
                                     !expr::is_tensor<V2>),
                                    int> = 0>
auto operator*(coeff_t const& value,
               OpConvolution<V2, GaussianSmoothing<D, grid_type>,
                             OpTerm<OpIdentity, G2>> const& b) {
  return symphas::internal::make_convolution::get_g(value * b.value, b.smoother,
                                                    b.op_g);
}

template <typename coeff_t, typename tensor_t, size_t D,
          template <typename, size_t> typename grid_type, typename G2,
          typename std::enable_if_t<
              (expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value,
               OpConvolution<tensor_t, GaussianSmoothing<D, grid_type>,
                             OpTerm<OpIdentity, G2>> const& b) {
  return (value * b.value) * symphas::internal::make_convolution::get_g(
                                 OpIdentity{}, b.smoother, b.op_g);
}

template <typename tensor_t, typename V, size_t D,
          template <typename, size_t> typename grid_type, typename G2,
          typename std::enable_if_t<
              (expr::is_tensor<tensor_t> && !expr::is_tensor<V> &&
               expr::eval_type<OpTerm<OpIdentity, G2>>::rank == 0),
              int> = 0>
auto operator*(tensor_t const& tensor,
               OpConvolution<V, GaussianSmoothing<D, grid_type>,
                             OpTerm<OpIdentity, G2>> const& b) {
  return symphas::internal::make_convolution::get(tensor * b.value, b.smoother,
                                                  b.op_g);
}

template <typename tensor_t, typename V, size_t D,
          template <typename, size_t> typename grid_type, typename G2,
          typename std::enable_if_t<
              (expr::is_tensor<tensor_t> && !expr::is_tensor<V> &&
               expr::eval_type<OpTerm<OpIdentity, G2>>::rank > 0),
              int> = 0>
auto operator*(tensor_t const& tensor,
               OpConvolution<V, GaussianSmoothing<D, grid_type>,
                             OpTerm<OpIdentity, G2>> const& b) {
  return symphas::internal::make_convolution::get(b.value, b.smoother,
                                                  tensor * b.op_g);
}

// ************************************************************************************************

namespace symphas::internal {

template <typename V, typename E1, typename E2>
auto construct_convolution(V const& value, OpExpression<E1> const& e1,
                           OpExpression<E2> const& e2) {
  return OpConvolution<V, E1, E2>(value, *static_cast<const E1*>(&e1),
                                  *static_cast<const E2*>(&e2));
}

template <typename V, typename E2>
auto construct_convolution(V const& value, OpVoid, OpExpression<E2> const&) {
  return OpVoid{};
}

template <typename V, typename E1>
auto construct_convolution(V const& value, OpExpression<E1> const&, OpVoid) {
  return OpVoid{};
}

template <typename V, typename E1, typename E2, size_t... Rs,
          size_t R = sizeof...(Rs),
          typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 &&
                                     expr::eval_type<E2>::rank == 0),
                                    int> = 0>
auto make_convolution_tensor(V const& value, OpExpression<E1> const& e1,
                             OpExpression<E2> const& e2,
                             std::index_sequence<Rs...>) {
  return (construct_convolution(
              expr::make_column_vector<Rs, R>() * value,
              expr::make_row_vector<Rs, R>() * (*static_cast<const E1*>(&e1)),
              *static_cast<const E2*>(&e2)) +
          ...);
}

template <typename V, typename E1, typename E2, size_t... Rs,
          size_t R = sizeof...(Rs),
          typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 &&
                                     expr::eval_type<E2>::rank > 0),
                                    int> = 0>
auto make_convolution_tensor(V const& value, OpExpression<E1> const& e1,
                             OpExpression<E2> const& e2,
                             std::index_sequence<Rs...>) {
  return (construct_convolution(
              expr::make_column_vector<Rs, R>() * value,
              *static_cast<const E1*>(&e1),
              expr::make_row_vector<Rs, R>() * (*static_cast<const E2*>(&e2))) +
          ...);
}

template <typename A, typename B>
auto make_convolution::get(A&& a, B&& b) {
  return get(OpIdentity{}, std::forward<A>(a), std::forward<B>(b));
}

template <typename V, typename E1, typename E2>
auto make_convolution::get(V const& value, OpExpression<E1> const& e1,
                           OpExpression<E2> const& e2) {
  constexpr size_t R =
      fixed_max<expr::eval_type<E1>::rank, expr::eval_type<E2>::rank>;
  if constexpr (R > 0) {
    return make_convolution_tensor(value, *static_cast<E1 const*>(&e1),
                                   *static_cast<E2 const*>(&e2),
                                   std::make_index_sequence<R>{});
  } else {
    return OpConvolution<V, E1, E2>(value, *static_cast<const E1*>(&e1),
                                    *static_cast<const E2*>(&e2));
  }
}

template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename E>
auto make_convolution::get(V const& value,
                           GaussianSmoothing<D, grid_type> const& smoother,
                           OpExpression<E> const& e) {
  constexpr size_t R = expr::eval_type<E>::rank;
  if constexpr (R > 0) {
    return make_convolution_tensor(value, smoother, *static_cast<E const*>(&e),
                                   std::make_index_sequence<R>{});
  } else {
    return OpConvolution<V, GaussianSmoothing<D, grid_type>, E>(
        value, smoother, *static_cast<E const*>(&e));
  }
}

template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename T, typename G>
auto make_convolution::get(V const& value,
                           GaussianSmoothing<D, grid_type> const& smoother,
                           OpTerm<T, G> const& e) {
  constexpr size_t R = expr::eval_type<OpTerm<T, G>>::rank;
  if constexpr (R > 0) {
    return make_convolution_tensor(value, smoother, e,
                                   std::make_index_sequence<R>{});
  } else {
    return OpConvolution<V, GaussianSmoothing<D, grid_type>,
                         OpTerm<OpIdentity, G>>(value, smoother, e);
  }
}

template <typename V, size_t D, template <typename, size_t> typename grid_type>
auto make_convolution::get(V const& value,
                           GaussianSmoothing<D, grid_type> const& smoother,
                           OpVoid) {
  return OpVoid{};
}

template <typename V, size_t D, template <typename, size_t> typename grid_type,
          typename G>
auto make_convolution::get_g(V const& value,
                             GaussianSmoothing<D, grid_type> const& smoother,
                             G grid) {
  return OpConvolution<V, GaussianSmoothing<D, grid_type>,
                       OpTerm<OpIdentity, G>>(value, smoother, grid);
}
}  // namespace symphas::internal
