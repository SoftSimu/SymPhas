
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
 * PURPOSE: Defines properties about a phase field crystal model.
 *
 * ***************************************************************************
 */

#pragma once

#include "modelspecialized.h"

/*
 * modifiers that affect the behaviour of the PFC model, in terms of different
 * qualities like conserved/non conserved dynamics
 * the recipe is, the pfc model inherits from these classes through the CRTP
 * construct, and these classes expose a function that the pfc model must use to
 * apply various functionalities
 */

/* the generic modifiers class that contains default parameters
 * it uses expression templates to find the correct parameter for the given
 * field
 */
template <typename Parameters>
struct PFCParametersDefault {
  template <size_t N>
  static constexpr symphas::internal::DynamicType dynamic_val_apply() {
    return Parameters::DEFAULT_DYNAMIC;
  }
  template <size_t N>
  static constexpr symphas::internal::DynamicType dynamic_val() {
    return Parameters::template dynamic_val_apply<N>();
  }

  template <size_t N>
  static constexpr size_t mode_val_apply() {
    return Parameters::DEFAULT_MODE_N;
  }
  template <size_t N>
  static constexpr size_t mode_val() {
    return Parameters::template mode_val_apply<N>();
  }

 protected:
  /* default parameter values
   */

  static const symphas::internal::DynamicType DEFAULT_DYNAMIC =
      symphas::internal::DynamicType::NONCONSERVED;
  static const size_t DEFAULT_MODE_N = 1;
};

/*
 * specializations of modifiers based on the modifiers parameters
 */

/* specialization of the dynamics for the field
 */

template <typename PFC>
struct PFCTraitDynamic {
  template <size_t N, symphas::internal::DynamicType dd>
  auto construct_dynamic() {
    return symphas::internal::apply_dynamics<dd>{}(
        OpLHS(expr::as_variable<N>(cast().template system<N>())),
        cast().template bulk_pfc_dynamic_N<N>() +
            cast().template coupled_pfc_dynamic_N<N>(),
        cast().solver);

    // return std::make_pair(OpLHS(expr::as_variable<N>(cast().template
    // system<N>())),
    //	/*OpIdentity{});*/ ((cast().template bulk_pfc_dynamic_N<N>())));
  }

 protected:
  PFC& cast() { return *static_cast<PFC*>(this); }
};

/* specialization of the mode for the field
 */

template <typename PFC>
struct PFCTraitMode {
  template <size_t N, size_t M, size_t M_N>
  auto construct_mode() {
    return construct_mode_apply<N, M>(std::make_index_sequence<M_N>{});
  }

 protected:
  template <size_t N, size_t M, size_t... Is>
  auto construct_mode_apply(std::index_sequence<Is...>) {
    return ((... * get_mode<N, M, Is>()));
  }

  template <size_t N, size_t M, size_t I>
  auto get_mode() {
    // auto qd = (1. + expr::make_operator_derivative<2>(cast().solver));
    // return qd * qd;
    auto qd = (expr::make_literal(cast().nu[N][M] * cast().nu[N][M]) +
               expr::make_operator_derivative<2>(cast().solver));
    return expr::make_literal(cast().alpha[N][M]) +
           expr::make_literal(cast().beta[N][N]) * qd * qd;
  }

  PFC& cast() { return *static_cast<PFC*>(this); }
};
