
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
#include <utility>

#include "modelpfctraits.h"

#define DEFAULT_ALPHA_COEFF -0.3

#define SN sizeof...(S)

//! A representation of the phase field crystal model.
/*!
 * A representation of the phase field crystal model. Implements the phase
 * field crystal model equation, and allows
 *
 * \tparam PFC The phase field crystal parameters specialized object.
 * \tparam D The dimension of the phase field crystal problem.
 * \tparam Sp The solver type for numerically solving the phase field crystal
 * problem.
 * \tparam S... The types of the order parameters.
 */
template <template <size_t, typename> typename PFC, size_t D, typename Sp,
          typename... S>
struct ModelPFCEquation : Model<D, Sp, S...>,
                          PFCTraitDynamic<PFC<D, Sp>>,
                          PFCTraitMode<PFC<D, Sp>> {
  using trait_dynamic = PFCTraitDynamic<PFC<D, Sp>>;
  using trait_mode = PFCTraitMode<PFC<D, Sp>>;

  using parent_type = Model<D, Sp, S...>;
  using parent_type::_s;
  using parent_type::coeff;
  using parent_type::num_coeff;
  using parent_type::solver;

  ModelPFCEquation(double const* coeff, size_t num_coeff,
                   symphas::problem_parameters_type const& parameters)
      : Model<D, Sp, S...>(coeff, num_coeff, parameters) {
    set_parameters();
  }

  ModelPFCEquation(symphas::problem_parameters_type const& parameters)
      : ModelPFCEquation(nullptr, 0, parameters) {}

  ModelPFCEquation(ModelPFCEquation<PFC, D, Sp, S...> const& other)
      : Model<D, Sp, S...>(*static_cast<Model<D, Sp, S...> const*>(&other)) {
    set_parameters();
  }

  ModelPFCEquation(ModelPFCEquation<PFC, D, Sp, S...>&& other)
      : Model<D, Sp, S...>(static_cast<Model<D, Sp, S...>&&>(other)) {
    set_parameters();
  }

  void set_parameters() {
    size_t num_interactions = (SN - 1) * ((SN - 1) + 1) / 2;
    size_t base_num_coeff = SN * 5;
    size_t total_num_coeff = base_num_coeff + num_interactions * 6;

    if (num_coeff < total_num_coeff) {
      symphas::lib::expand_append_array(DEFAULT_COEFF_VALUE, coeff,
                                        total_num_coeff - num_coeff, num_coeff);

      for (size_t i = ((num_coeff + 5) - base_num_coeff) / 6;
           i < num_interactions; ++i) {
        size_t offset_interaction = base_num_coeff + i * 6;
        coeff[offset_interaction] = DEFAULT_ALPHA_COEFF;
      }

      num_coeff = total_num_coeff;
    }

    for (iter_type i = 0; i < SN; ++i) {
      alpha[i][i] = coeff[i];
      beta[i][i] = coeff[i + SN];
      nu[i][i] = coeff[i + 2 * SN];
      gamma[i][i] = coeff[i + 3 * SN];
      delta[i][i] = coeff[i + 4 * SN];
      eps[i][i] = 1.0;

      for (iter_type j = i + 1; j < SN; ++j) {
        size_t offset_group = (SN * (SN - 1)) / 2;  // offset from groups
        size_t offset_last = offset_group - ((SN - i) * (SN - (i + 1))) / 2;
        size_t offset_noninteraction = SN * 5;
        size_t offset_interaction =
            offset_noninteraction  // offset from the noninteracting
                                   // coefficients
            +
            (j -
             (i +
              1))  // offset from the previous coupling index for the current i
            + offset_last;  // offset from the coupling index of the previous
                            // fields

        alpha[i][j] = alpha[j][i] = coeff[offset_interaction];
        beta[i][j] = beta[j][i] = coeff[offset_interaction + offset_group];
        nu[i][j] = nu[j][i] = coeff[offset_interaction + 2 * offset_group];
        gamma[i][j] = gamma[j][i] =
            coeff[offset_interaction + 3 * offset_group];
        delta[i][j] = delta[j][i] =
            coeff[offset_interaction + 4 * offset_group];
        eps[i][j] = eps[j][i] = coeff[offset_interaction + 5 * offset_group];
      }
    }
  }

  // phase field crystal coefficients

  double alpha[SN][SN];  //!< Coefficients representing \f$\alpha_{ij}\f$
  double beta[SN][SN];   //!< Coefficients representing \f$\beta_{ij}\f$
  double nu[SN][SN];     //!< Coefficients representing \f$\nu_{ij}\f$
  double gamma[SN][SN];  //!< Coefficients representing \f$\gamma_{ij}\f$
  double delta[SN][SN];  //!< Coefficients representing \f$\delta_{ij}\f$
  double eps[SN][SN];    //!< Coefficients representing \f$\epsilon_{ij}\f$

  template <size_t N>
  decltype(auto) get_field_op() {
    std::ostringstream ss;
    ss << "n_" << N;
    return expr::make_term<N>(
        NamedData(parent_type::template grid<N>(), ss.str()));
  }

  // build the bulk phase equation
  // getting the grid this is associated with
  template <size_t N>
  auto bulk_pfc_dynamic_N() {
    auto n = get_field_op<N>();
    return (expr::make_literal(gamma[N][N]) +
            expr::make_literal(delta[N][N]) * n) *
               n * n +
           get_mode_N<N>() * n;
  }

  // build the bulk free energy equation
  // getting the grid this is associated with
  template <size_t N>
  auto bulk_pfc_fe_N() {
    auto n = get_field_op<N>();
    return (expr::make_literal(gamma[N][N] / 3) +
            expr::make_literal(delta[N][N] / 4) * n) *
               n * n * n +
           n * get_mode_N<N>() * (n * expr::make_literal(0.5));
  }

  template <size_t N>
  auto coupled_pfc_dynamic_N() {
    return coupled_pfc_dynamic_N<N>(std::make_index_sequence<SN - 1>{});
  }

  // the coupled EVOLUTION construct a convolution term using the existing
  // convolution term eta, and hijack that grid to save time and space
  template <size_t N, size_t M>
  auto coupled_pfc_dynamic_NM() {
    auto n = get_field_op<N>();
    auto m = get_field_op<M>();
    auto dims = parent_type::template system<M>().dims;
    auto widths = parent_type::template system<M>().get_info().get_widths();

    auto G = GaussianSmoothing<D>{dims, widths, 1.0};
    return expr::make_literal(alpha[N][M]) * m +
           expr::make_literal(beta[N][M]) * get_mode_NM<N, M>() * m +
           expr::make_literal(0.5 * gamma[N][M]) *
               (expr::make_literal(2.) * n * m + m * m) +
           expr::make_literal(eps[N][M]) *
               expr::make_convolution(G, eta_N<M>());
  }

  template <size_t N, size_t M>
  auto coupled_pfc_fe_NM() {
    auto n = get_field_op<N>();
    auto m = get_field_op<M>();

    return expr::make_literal(alpha[N][M]) * m +
           expr::make_literal(beta[N][M]) * get_mode_NM<N, M>() * m +
           expr::make_literal(0.5 * gamma[N][M]) * n * (n * m + m * m) +
           expr::make_literal(eps[N][M]) * eta_N<N>() * eta_N<M>();
  }

  template <size_t N>
  auto eta_N() {
    auto dims = parent_type::template system<N>().dims;
    auto widths = parent_type::template system<N>().get_info().get_widths();

    auto G = GaussianSmoothing<D>{dims, widths, 1.0};
    return expr::make_convolution(G, get_field_op<N>());
  }

  // get the dynamic equation dictating evolution of field N
  template <size_t N>
  auto get_dynamic_N() {
    return trait_dynamic::template construct_dynamic<
        N, PFC<D, Sp>::template dynamic_val<N>()>();
  }

  // get the mode operator for the field N
  template <size_t N, size_t M>
  auto get_mode_NM() {
    return trait_mode::template construct_mode<
        N, M, PFC<D, Sp>::template mode_val<N>()>();
  }

  // get the mode operator for the field N
  template <size_t N>
  auto get_mode_N() {
    return get_mode_NM<N, N>();
  }

  auto make_equations() {
    return _make_equations(std::make_index_sequence<sizeof...(S)>{});
  }

 protected:
  // constructs an array of SN - 1 terms, skipping N
  // equivalent to the indices covered by the sum in the coupling of the
  // dynamical equation
  template <size_t N>
  static constexpr auto make_arr_N() {
    std::array<int, SN - 1> result = {};
    for (int i = 0; i < SN - 1; ++i) {
      result[i] = (i < N) ? i : (i + 1);
    }
    return result;
  }

  // implicit function for generating the coupled terms
  template <size_t N, size_t... Is>
  auto coupled_pfc_dynamic_N(std::index_sequence<Is...>) {
    return (
        (OpVoid{} + ... + coupled_pfc_dynamic_NM<N, make_arr_N<N>()[Is]>()));
  }

  template <size_t... Is>
  auto _make_equations(std::index_sequence<Is...>) {
    ((..., expr::printe(std::get<1>(get_dynamic_N<Is>()), "given equation")));
    return solver.template form_expr_all<SN>(_s, get_dynamic_N<Is>()...);
  }
};

template <template <size_t, typename> typename PFC, size_t D, typename Sp,
          typename... S>
using GeneralizedPFCModel = typename symphas::internal::expand_types_to_model_t<
    D, Sp, S...>::template Specialized<symphas::internal::MakeEquation,
                                       ModelPFCEquation<PFC, D, Sp, S...>>;

#undef SN

template <template <size_t, typename> typename PFC, size_t D, typename Sp,
          typename... S>
struct ModelPFCEquation<PFC, D, Sp, symphas::lib::types_list<S...>>
    : ModelPFCEquation<PFC, D, Sp, S...> {
  using parent_type = ModelPFCEquation<PFC, D, Sp, S...>;
  using parent_type::parent_type;
};
