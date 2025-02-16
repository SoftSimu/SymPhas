
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
 * PURPOSE: Forward declaration of all expression objects.
 *
 * ***************************************************************************
 */

#pragma once

#include "spslib.h"

template <typename Sp, size_t N = 0>
struct Solver;

//! Contains all elements constituting the symbolic algebra functionality.
/*!
 * Defines elements which support the symbolic algebra functionality, such
 * as creating new expression terms, manipulating and transforming expressions,
 * and printing expressions.
 */
namespace expr {
using exp_key_t = unsigned int;
template <typename G>
struct variational_t {};

}  // namespace expr

namespace symphas::internal {
//! Convert values of the exponent into an exponent key.
/*!
 * Convert values of the exponent into an exponent key. If \p b is
 * given as zero, then it will consider a to be a key, and return the
 * real value of the exponent from the key.
 *
 * \tparam a The numerator of the exponent.
 * \tparam b The denominator of the exponent.
 * \tparam sign If true, then it means the exponent is negative.
 */
template <expr::exp_key_t a, expr::exp_key_t b = 1, bool sign = false>
struct exp_compute_key;

using expr::exp_key_t;

inline constexpr unsigned int xS =
    sizeof(exp_key_t) * 8;  //<! Size of the exponent key.
inline constexpr unsigned int xsm =
    (1u << (xS - 1u));  //<! Mask of the sign bit.
inline constexpr unsigned int xam =
    (~0u >> (xS >> 1u));  //<! Mask of the numerator value.
inline constexpr unsigned int xbm =
    ~(xam | xsm);  //<! Mask of the denominator value.

template <exp_key_t a, exp_key_t b, bool sign>
struct exp_compute_key {
 protected:
  static const exp_key_t GCD = GCD_of<a, b>;

 public:
  static const exp_key_t N = a / GCD;
  static const exp_key_t D = b / GCD;

  static const exp_key_t value =
      ((sign) ? xsm : 0ul) | (N & xam) | (((D - 1) << (xS >> 1ul)) & xbm);
};

template <exp_key_t v>
struct exp_compute_key<v, 0, false> {
  static const exp_key_t N = (v & xam);
  static const exp_key_t D = ((v & xbm) >> (xS >> 1ul)) + 1;
  static const bool sign = (v & xsm) >> (xS - 1);

  static constexpr double value = (double(N) / D) * ((sign) ? -1 : 1);
  using type = exp_compute_key<N, D, sign>;
};

struct tensor_fold;
struct tensor_cast;

//! A mapping to map half complex plane (FFTW format) to sequential complex
//! values.
struct HCTS {};
//! A sequential complex values to half complex plane (FFTW format).
struct STHC {};
}  // namespace symphas::internal

namespace expr {
template <exp_key_t a, exp_key_t b = 1, bool sign = false>
using Xk_t = symphas::internal::exp_compute_key<a, b, sign>;

template <exp_key_t a>
using _Xk_t = symphas::internal::exp_compute_key<a, 0, false>;

//! Exponent key helper to convert fraction into an exponent key.
template <exp_key_t a, exp_key_t b = 1, bool sign = false>
constexpr auto Xk = Xk_t<a, b, sign>::value;

//! Exponent key helper to convert an exponent key into an exponent value.
template <exp_key_t value>
constexpr auto _Xk = Xk<value, 0, false>;

template <bool sign1, bool sign2, exp_key_t X1, exp_key_t X2>
constexpr int compute_XXk_N = int(expr::_Xk_t<X1>::N * expr::_Xk_t<X2>::D +
                                  expr::_Xk_t<X2>::N * expr::_Xk_t<X1>::D);
template <exp_key_t X1, exp_key_t X2>
constexpr int compute_XXk_N<true, true, X1, X2> =
    -compute_XXk_N<false, false, X1, X2>;
template <exp_key_t X1, exp_key_t X2>
constexpr int compute_XXk_N<false, true, X1, X2> =
    int(expr::_Xk_t<X1>::N * expr::_Xk_t<X2>::D) -
    int(expr::_Xk_t<X2>::N * expr::_Xk_t<X1>::D);
template <exp_key_t X1, exp_key_t X2>
constexpr int compute_XXk_N<true, false, X1, X2> =
    int(expr::_Xk_t<X2>::N * expr::_Xk_t<X1>::D) -
    int(expr::_Xk_t<X1>::N * expr::_Xk_t<X2>::D);

template <exp_key_t X1, exp_key_t X2>
constexpr int switch_XXk_N =
    compute_XXk_N<_Xk_t<X1>::sign, _Xk_t<X2>::sign, X1, X2>;

template <bool flag, exp_key_t X1, exp_key_t X2>
constexpr exp_key_t abs_XXk_N = exp_key_t(switch_XXk_N<X1, X2>);
template <exp_key_t X1, exp_key_t X2>
constexpr exp_key_t abs_XXk_N<false, X1, X2> = exp_key_t(-switch_XXk_N<X1, X2>);

template <exp_key_t X1, exp_key_t X2>
constexpr exp_key_t _XXk_N = abs_XXk_N<(switch_XXk_N<X1, X2> >= 0), X1, X2>;

template <exp_key_t X1, exp_key_t X2>
constexpr exp_key_t _XXk_D = expr::_Xk_t<X1>::D * expr::_Xk_t<X2>::D;

template <exp_key_t X1, exp_key_t X2>
using XXk_t = Xk_t<_XXk_N<X1, X2>, _XXk_D<X1, X2>, (switch_XXk_N<X1, X2> < 0)>;

namespace symbols {
struct Symbol {
  auto eval(iter_type = 0) const { return *this; }

  inline auto operator+(expr::symbols::Symbol) const {
    return expr::symbols::Symbol{};
  }

  inline auto operator-(expr::symbols::Symbol) const {
    return expr::symbols::Symbol{};
  }

  inline auto operator*(expr::symbols::Symbol) const {
    return expr::symbols::Symbol{};
  }

  inline auto operator/(expr::symbols::Symbol) const {
    return expr::symbols::Symbol{};
  }

  template <typename E>
  auto operator=(E) {
    return expr::symbols::Symbol{};
  }

  template <typename E>
  auto operator+(E) const {
    return expr::symbols::Symbol{};
  }

  template <typename E>
  auto operator-(E) const {
    return expr::symbols::Symbol{};
  }

  template <typename E>
  auto operator*(E) const {
    return expr::symbols::Symbol{};
  }

  template <typename E>
  auto operator/(E) const {
    return expr::symbols::Symbol{};
  }

  auto& operator+=(Symbol) { return *this; }

  auto& operator-=(Symbol) { return *this; }
};

template <typename T>
struct SymbolType {};
}  // namespace symbols

template <size_t N>
auto pow(symbols::Symbol const&) {
  return symbols::Symbol{};
}

}  // namespace expr

template <typename E>
auto operator+(E, expr::symbols::Symbol);
template <typename E>
auto operator-(E, expr::symbols::Symbol);
template <typename E>
auto operator*(E, expr::symbols::Symbol);
template <typename E>
auto operator/(E, expr::symbols::Symbol);
namespace expr::symbols {
using ::operator+;
using ::operator-;
using ::operator*;
using ::operator/;
}  // namespace expr::symbols

/*!
 * \defgroup Op Symbolic Algebra Objects
 * @{
 */

template <typename T, size_t... Ns>
struct OpTensor;
struct OpIdentity;
struct OpNegIdentity;
struct OpVoid;
template <typename T>
struct OpLiteral;
template <typename T, typename I = void>
struct OpCoeff;
template <size_t N, size_t D>
struct OpFractionLiteral;
template <size_t N, size_t D>
struct OpNegFractionLiteral;

template <typename... Es>
struct OpAdd;
template <typename... Es>
struct OpAddList;
template <typename E1, typename E2>
struct OpBinaryMul;
template <typename E1, typename E2>
struct OpBinaryDiv;

template <typename E>
struct OpEvaluable;
template <typename E>
struct OpExpression;
template <typename E>
struct OpOperator;

template <typename T, size_t D>
struct GridData;
template <size_t Z, typename G = expr::symbols::Symbol>
struct Variable;
template <typename G = expr::symbols::Symbol>
struct DynamicVariable;
template <typename G>
struct NamedData;
struct DynamicIndex;
struct DynamicIndexSet;

#ifdef USING_CUDA

template <typename T, size_t D>
struct GridDataCUDA;
#endif

template <typename T, size_t D>
using GridSymbol = GridData<expr::symbols::SymbolType<T>, D>;

template <typename G, expr::exp_key_t X = expr::Xk<1>>
struct Term;
template <typename... Ts>
struct OpTermsList;
template <typename... Ts>
struct OpTerms;
template <typename V, typename G>
using OpTerm = OpTerms<V, Term<G, 1>>;

template <typename G>
struct SymbolicDerivative;
template <typename Dd, typename V, typename E, typename Sp>
struct OpDerivative;
template <typename V, typename E, typename T>
struct OpIntegral;
template <size_t O, typename V, typename T>
struct OpOperatorDerivative;
template <Axis ax, size_t O, typename V, typename Sp>
struct OpOperatorDirectionalDerivative;
template <typename V, typename Sp, size_t... Os>
struct OpOperatorMixedDerivative;

template <typename V, typename E1, typename E2>
struct OpConvolution;
template <size_t D>
struct GaussianSmoothing;

template <typename A1, typename A2, typename E>
struct OpCombination;
template <typename A1, typename A2>
struct OpOperatorCombination;
template <typename A1, typename A2, typename E>
struct OpChain;
template <typename A1, typename A2>
struct OpOperatorChain;
template <typename E>
struct OpOperator;
template <typename V, typename E>
struct OpExponential;
template <expr::exp_key_t X, typename V, typename E>
struct OpPow;
template <typename G, typename V, typename E>
struct OpMap;

template <typename V, typename E, typename F, typename Arg, typename... Args>
struct OpFunction;
template <auto f, typename V, typename E>
struct OpFunctionApply;
template <typename V, typename sub_t, typename... Ts>
struct OpSymbolicEval;

template <typename E>
struct OpOptimized;

template <typename T>
struct SymbolicData;

template <size_t D, Axis ax>
struct GridAxis;

template <typename I>
using OpCoeffSwap = OpCoeff<void, I>;

template <typename G>
using SymbolicFunctionalDerivative = SymbolicDerivative<expr::variational_t<G>>;

//! @}

template <typename specialized_eval_handler>
struct BaseEvalHandler {
  template <typename E, typename assign_type, typename interval_type>
  void result(OpEvaluable<E> const& e, assign_type&& assign,
              interval_type&& interval) const = delete;
  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const& e, assign_type&& assign) const = delete;
  template <typename E, typename assign_type, typename interval_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign,
                  interval_type&& interval) const = delete;
  template <typename E, typename assign_type>
  void result_sum(OpEvaluable<E> const& e, assign_type&& assign) const = delete;
};

namespace expr::prune {

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
                   BaseEvalHandler<eval_handler_type> const& eval_handler);

template <typename... condition_ts, typename E>
inline void update(E& e);
}  // namespace expr::prune

namespace expr {

template <Axis ax, typename V, typename E>
auto to_axis(V const& value, OpExpression<E> const& e);
template <Axis ax, typename E>
auto to_axis(OpExpression<E> const& e);

enum class NoiseType { WHITE, NONE, DECAY_EXP, DECAY_POLY, POISSON };

}  // namespace expr

template <expr::NoiseType nt, typename T, size_t D>
struct NoiseData;

#ifdef USING_CUDA

template <typename E>
struct CuEvaluable {};

#endif

template <typename T>
struct CUDADataType {};

namespace expr {
//! Defines elements which transform an expression into another one.
/*!
 * Implements the functions that transform an expression of a given
 * form into another form. The functions typically constructing another
 * expression based on some rules.
 */
namespace transform {}
//! Prune an expression to update the state of all nested terms.
/*!
 * Specifies algorithms for different expression types that will update
 * the underlying data before the expression can be evaluated.
 */
namespace prune {}

//! Defines ways to split functions based on criteria.
/*!
 * Implements all functions that will split an expression based
 * on certain criteria.
 */
namespace split {}
}  // namespace expr
