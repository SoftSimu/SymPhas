
/* ***************************************************************************
 * This file is part of the SymPhas package, containing a framework for
 * implementing solvers for phase-field problems with compile-time symbolic
 * algebra.
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
 * This file is part of the distributed solvers.
 *
 * FORWARD-TIME CENTRAL-SPACE METHOD (FORWARD EULER)
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"

#ifndef SYMPHAS_MPI_PROFILE_SCOPE
#define SYMPHAS_MPI_PROFILE_SCOPE(name) ((void)0)
#endif
#include "solver.h"

//! Finite difference solver.
/*!
 * The boundaries are updated after the provisional variables are computed
 * from the corresponding equations. The boundary data is typically the same
 * as the first field that is given, but in practise it does not matter
 * for the numerical results unless finite difference approximations are
 * applied extensively to the provisional variables. This is not recommended.
 */
START_NEW_SOLVER_WITH_STENCIL(SolverFT)

/*
 * forward euler method for actually updating the grid
 */
template <typename S>
void step(S& sys) const {
  SYMPHAS_MPI_PROFILE_SCOPE("step_euler");
  expr::result(expr::make_term(sys.as_grid()) + expr::make_term(dt, sys.dframe),
               sys.as_grid(), expr::iterable_domain(sys.as_grid()));
}

/*
 * the parameter is a tuple with the first parameter being a Variable object
 * with ref base type, and the ref is a reference to the system; like this, the
 * variable index is packaged and the equality operator is deferred to the
 * oplvariable
 */
template <typename S, typename E>
inline void equation(std::pair<S, E>* r, len_type len) const {
  TIME_THIS_CONTEXT_LIFETIME(solverft_equation_list);

  using L = decltype(r[0].first.get().dframe);
  using R = E;

  {
    expr::eval_handler_type<E> handler;

    SYMPHAS_OMP_PARALLEL_DIRECTIVE
    for (iter_type i = 0; i < len; ++i) {
      expr::prune::update<expr::not_<expr::matches_series>>(r[i].second,
                                                            handler);
    }
  }
  {
    expr::eval_handler_type<E> handler;
    for (iter_type i = 0; i < len; ++i) {
      handler.result(r[i].second, r[i].first.get().dframe);
    }
  }

  // expr::result_by_term<expr::matches_series>(r.second, r.first.get().dframe);
}

/*
 * the parameter is a tuple with the first parameter being a Variable object
 * with ref base type, and the ref is a reference to the system; like this, the
 * variable index is packaged and the equality operator is deferred to the
 * oplvariable
 */
template <typename S, typename E>
inline void equation(std::pair<S, E>& r) const {
  TIME_THIS_CONTEXT_LIFETIME(solverft_equation);
  SYMPHAS_MPI_PROFILE_SCOPE("equation_total");

  {
    SYMPHAS_MPI_PROFILE_SCOPE("equation_prune_update");
    expr::eval_handler_type<E> handler;
    expr::prune::update<expr::not_<expr::matches_series>>(r.second, handler);
  }

  {
    SYMPHAS_MPI_PROFILE_SCOPE("equation_result");
    expr::result(r.second, r.first.get().dframe);
  }
}

/*
 * some solvers require the equation that is provided
 * to be in a slightly different form
 * the forward euler solver has no such requirement and so
 * it returns its argument directly
 */

template <size_t En, typename SS, typename S, typename E>
auto form_expr_one(SS&&, std::pair<S, E> const& e) const {
  auto [sys, equation] = e;
  auto eq_ft = expr::transform::optimize(expr::apply_operators(equation));
  expr::prune::update(eq_ft);
  expr::printe(eq_ft, "scheme");
  return std::make_pair(sys, eq_ft);
}

static auto make_solver(symphas::problem_parameters_type const& parameters) {
  if (parameters.length()) {
    double h = parameters.get_interval_data()[0].at(Axis::X).width();
    size_t dim = parameters.get_dimension();

    // integrity check: all grid widths must be the same
    for (iter_type i = 0; i < parameters.length(); ++i) {
      for (iter_type n = 0; n < dim; ++n) {
        if (h != parameters.get_interval_data()[i]
                     .at(symphas::index_to_axis(n))
                     .width()) {
          char axis = (symphas::index_to_axis(n) == Axis::X)   ? 'x'
                      : (symphas::index_to_axis(n) == Axis::Y) ? 'y'
                      : (symphas::index_to_axis(n) == Axis::Z) ? 'z'
                                                               : '?';
          fprintf(
              SYMPHAS_WARN,
              "the grid spacing of system %d for axis '%c' is "
              "not consistent, results will not reflect the given problem!\n",
              i, axis);
        }
      }
    }

    /* the dimensions of the problem are taken from the first system
     * since the solver assumes they are homogeneous
     */
    len_type* dims = new len_type[dim];
    for (iter_type i = 0; i < dim; ++i) {
      Axis side = symphas::index_to_axis(i);
      dims[i] = parameters.get_interval_data()[0].at(side).get_count() +
                2 * BOUNDARY_DEPTH;
    }

    auto s = this_type{dims, h};
    delete[] dims;
    return s;
  } else {
    return this_type{grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3)};
  }
}

// given the grid/equation pair, evaluate the equation into the grid
// element by element
template <typename G, typename E>
void evaluate_one(std::pair<G, E>& r) const {
  auto& [grid, equation] = r;
  if constexpr (expr::has_state<E>::value) {
    expr::prune::update(equation);
  }

  expr::result(equation, expr::BaseData<G>::get(grid));
}

template <typename G, typename E>
void evaluate_one(std::pair<G, E>&& r) const {
  auto& [grid, equation] = r;
  if constexpr (expr::has_state<E>::value) {
    expr::prune::update(equation);
  }

  expr::result(equation, expr::BaseData<G>::get(grid));
}
END_SOLVER

template <size_t D>
bool check_overlapping_domain(grid::region_interval<D> const& region0,
                              const iter_type (&dims)[D],
                              const iter_type (&origin)[D]) {
  grid::region_interval<D> region1(origin, dims);
  return grid::is_overlapping(region0, region1);
}

template <size_t D>
bool check_overlapping_domain(grid::region_interval_multiple<D> region,
                              const iter_type (&dims)[D],
                              const iter_type (&origin)[D]) {
  grid::region_interval<D> region0(origin, dims);
  region /= region0;
  if (region.regions.empty()) {
    return false;
  } else {
    return true;
  }
}

#ifdef USING_MPI
#include "solverftmpisync.h"

template <typename T, size_t D>
template <typename... Ts>
void PhaseFieldSystem<RegionalGridMPI, T, D>::synchronize(Ts&&... args) {
  //::synchronize(std::forward<Ts>(args)...);
}

#endif

// Build-time guard: under an MPI build, RegionalGrid variants are NOT
// registered as selectable solver systems. RegionalGrid assumes a single
// process owns the full sparse domain, which is incompatible with MPI's
// spatial decomposition. Valid `solver_variation` values:
//   serial build   : 0 = SolverSystemFD (Grid), 1 = SolverSystemFDwSD (RegionalGrid)
//   MPI build      : 0 = SolverSystemFD (Grid) only
// The `SolverSystemFDwSDMPI` class body itself carries a static_assert
// (see sol/inc/solversystem.h) that fires if anything tries to instantiate it.
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFD)
#ifndef USING_MPI
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDwSD)
#endif

#ifdef USING_CUDA
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDCUDA)
ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDwSDCUDA)
#endif

// NOTE: SolverSystemFDwSDMPI (RegionalGrid + MPI) is intentionally NOT
// registered. See the static_assert guard in sol/inc/solversystem.h.

ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverFT, ProvisionalSystemFD)

#ifdef USING_CUDA
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE_FOR_SPECIFIC(SolverFT, SolverSystemFDCUDA,
                                               ProvisionalSystemFDCUDA)
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE_FOR_SPECIFIC(SolverFT, SolverSystemFDwSDCUDA,
                                               ProvisionalSystemFDCUDA)
#endif

SYMPHAS_SOLVER_ALL_SUPPORTED(SolverFT)
