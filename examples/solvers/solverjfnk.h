
// Uncomment to enable verbose JFNK debugging output
// #define JFNK_VERBOSE
// #define JFNK_VERBOSE
#define JFNK_VERBOSE

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
 * JACOBIAN-FREE NEWTON-KRYLOV (JFNK) SOLVER
 *
 * Implements a fully implicit solver using Newton-GMRES with Jacobian-free
 * matrix-vector products computed via finite differences:
 *
 *   J*v ≈ (G(u + ε*v) - G(u)) / ε
 *
 * The Newton iteration solves the nonlinear system:
 *
 *   G(u^{n+1}) = u^{n+1} - u^n - dt * F(u^{n+1}) = 0
 *
 * where F(u) is the RHS of the evolution equations (du/dt = F(u)).
 *
 * The inner linear solve uses restarted GMRES(m) with m=20 by default.
 *
 * KEY DESIGN DECISION: Uses SolverSystemFD (real-space BoundaryGrid dframe),
 * NOT SolverSystemSpectral. The PoissonSolver operator embedded in
 * expressions (e.g., MagneticPFC2013) handles its own FFTs internally and
 * requires BoundaryGrid<T,D> as input. expr::result(expr, sys.dframe)
 * writes F(u) as real-valued data in sys.dframe.values.
 *
 * Design within SymPhas architecture:
 * - form_expr_one: stores the optimized RHS expression (like SolverFT)
 * - equation: On first field, collects all field info. On last field,
 *   triggers the full Newton-GMRES solve. All subsequent calls are no-ops.
 * - step: no-op (solution already written into grid data by Newton solve)
 *
 * VECTOR FIELD LAYOUT:
 * Grid<vector_t<D>,D> inherits MultiBlock<D,scalar_t>, which stores D
 * SEPARATE scalar_t* arrays (NOT interleaved). sys.values is really
 * scalar_t* values[D]. The JFNK solver registers each axis component
 * as a separate DOF entry. For a 2D vector field M, this means 2 DOFs
 * (Mx, My), each of length sys.len (= number of grid points).
 *
 * BOUNDARY REFRESH:
 * Between Newton iterations and during Jacobian-free matvec perturbations,
 * sys.update() is called to refresh boundary ghost cells. This is critical
 * because the main model loop only calls update() once at the start of each
 * timestep, and JFNK modifies grid interior data multiple times.
 *
 * Reference: Validated against CUDA JFNK implementation for MagneticPFC2013
 * (Faghihi et al., PRE 88, 032407, 2013) at 512×512 producing correct
 * multi-domain magnetic structures.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"
#include "solver.h"

#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

// =============================================================================
// JFNK Solver Parameters
// =============================================================================

namespace jfnk {

//! Default GMRES restart dimension.
constexpr int GMRES_M = 20;

//! Maximum Newton iterations per timestep.
constexpr int MAX_NEWTON_ITERS = 50;

//! Absolute tolerance for Newton convergence (RMS of ||G||).
constexpr double NEWTON_ATOL = 1e-6;

//! Relative tolerance for Newton convergence.
constexpr double NEWTON_RTOL = 1e-4;

//! Base for finite-difference perturbation (≈ sqrt(machine epsilon)).
constexpr double EPSILON_BASE = 1.49011611938476e-08;

//! Maximum GMRES iterations (with restarts).
constexpr int MAX_GMRES_ITERS = 5000;

//! GMRES convergence tolerance (relative to Newton residual norm).
constexpr double GMRES_RTOL = 1e-3;

// =============================================================================
// Type Traits for Vector Field Detection
// =============================================================================

//! Primary template: T is scalar.
template <typename T>
struct field_traits {
    static constexpr int num_components = 1;
};

//! Specialization for any_vector_t<T,N> — stored in MultiBlock<N,T>
//! with N separate scalar_t* arrays (values[0]..values[N-1]).
template <typename T, size_t N>
struct field_traits<any_vector_t<T, N>> {
    static constexpr int num_components = static_cast<int>(N);
};

// =============================================================================
// Linear Algebra Helpers
// =============================================================================

inline double dot_product(const double* a, const double* b, len_type n) {
    double sum = 0.0;
    for (len_type i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

inline double norm2(const double* v, len_type n) {
    return std::sqrt(dot_product(v, v, n));
}

inline void axpy(double alpha, const double* x, double* y, len_type n) {
    for (len_type i = 0; i < n; ++i) {
        y[i] += alpha * x[i];
    }
}

inline void scale_vec(double alpha, double* v, len_type n) {
    for (len_type i = 0; i < n; ++i) {
        v[i] *= alpha;
    }
}

inline void copy_vec(const double* src, double* dst, len_type n) {
    std::memcpy(dst, src, sizeof(double) * static_cast<size_t>(n));
}

inline void zero_vec(double* v, len_type n) {
    std::memset(v, 0, sizeof(double) * static_cast<size_t>(n));
}

inline void apply_givens(double& h_i, double& h_ip1, double cs_i, double sn_i) {
    double temp = cs_i * h_i + sn_i * h_ip1;
    h_ip1 = -sn_i * h_i + cs_i * h_ip1;
    h_i = temp;
}

inline void compute_givens(double a, double b, double& cs, double& sn) {
    if (std::abs(b) < 1e-30) {
        cs = 1.0; sn = 0.0;
    } else if (std::abs(b) > std::abs(a)) {
        double tau = a / b;
        sn = 1.0 / std::sqrt(1.0 + tau * tau);
        cs = tau * sn;
    } else {
        double tau = b / a;
        cs = 1.0 / std::sqrt(1.0 + tau * tau);
        sn = tau * cs;
    }
}

}  // namespace jfnk


// =============================================================================
// SolverJFNK Definition
// =============================================================================
//
// Uses SolverSystemFD (BoundarySystem with real-space BoundaryGrid dframe).
// Expression evaluation via expr::result(expr, sys.dframe) produces F(u)
// in real space, which is exactly what JFNK needs.
//

START_NEW_SOLVER_WITH_STENCIL(SolverJFNK)

    // =========================================================================
    // Per-component DOF info (type-erased for the Newton solver)
    // =========================================================================

    //! Type-erased DOF data per scalar component.
    //!
    //! For scalar fields: one FieldDOF per equation.
    //! For vector_t<D> fields: D FieldDOFs per equation (one per axis
    //! component). All share the same eval_rhs callback (which evaluates
    //! the full vector F(u) into sys.dframe), but grid_values and
    //! dframe_values point to individual component arrays from the
    //! MultiBlock<D,scalar_t> storage.
    //!
    //! The update_boundaries callback refreshes ghost cells after writing
    //! to grid_values (needed for stencil correctness in JFNK iterations).
    //!
    //! INTERIOR INDEXING: Grid arrays (grid_values, dframe_values) include
    //! boundary ghost cells. The JFNK DOF vector operates only on interior
    //! cells. `interior_map[j]` gives the grid flat index for DOF j.
    //! All workspace arrays (u_n, F_u, residual, saved) are sized to
    //! `interior_len` and indexed by dense DOF indices.
    struct FieldDOF {
        double* grid_values;    //!< Component's scalar_t* from sys (or MultiBlock).
        double* dframe_values;  //!< Component's scalar_t* from sys.dframe.
        len_type grid_len;      //!< Total grid array length (including boundary).
        len_type interior_len;  //!< Number of interior cells (DOF size for this comp).

        //! Map from dense DOF index j to grid flat index.
        //! interior_map[j] ∈ [0, grid_len) for j ∈ [0, interior_len).
        std::vector<len_type> interior_map;

        //! Callback: evaluate F(u) into sys.dframe (shared for all
        //! components of the same vector field — call once per field).
        std::function<void()> eval_rhs;

        //! Callback: refresh boundary ghost cells after grid modification.
        std::function<void()> update_boundaries;

        //! Per-DOF group tag. All components from the same field share a tag.
        //! eval_rhs is called once per group during evaluate_all_rhs().
        int field_group;

        // Newton workspace — sized to interior_len (allocated in allocate()).
        double* u_n;        //!< Saved solution at start of timestep.
        double* F_u;        //!< F(u) at current Newton iterate.
        double* residual;   //!< Newton residual G(u) = u - u_n - dt*F(u).

        // Pre-allocated save buffer for matvec.
        double* saved;      //!< Saved grid values during matvec perturbation.

        FieldDOF() : grid_values{nullptr}, dframe_values{nullptr},
                     grid_len{0}, interior_len{0},
                     field_group{-1},
                     u_n{nullptr}, F_u{nullptr}, residual{nullptr},
                     saved{nullptr} {}

        void allocate(len_type n_interior) {
            interior_len = n_interior;
            u_n = new double[n_interior]();
            F_u = new double[n_interior]();
            residual = new double[n_interior]();
            saved = new double[n_interior]();
        }

        void deallocate() {
            delete[] u_n;       u_n = nullptr;
            delete[] F_u;       F_u = nullptr;
            delete[] residual;  residual = nullptr;
            delete[] saved;     saved = nullptr;
        }

        //! Gather interior cells from a grid array into a dense buffer.
        void gather(const double* grid_arr, double* dense) const {
            for (len_type j = 0; j < interior_len; ++j) {
                dense[j] = grid_arr[interior_map[j]];
            }
        }

        //! Scatter dense buffer values INTO a grid array's interior cells.
        void scatter_add(const double* dense, double alpha,
                         double* grid_arr) const {
            for (len_type j = 0; j < interior_len; ++j) {
                grid_arr[interior_map[j]] += alpha * dense[j];
            }
        }

        //! Overwrite interior cells of grid_arr from dense buffer.
        void scatter_set(const double* dense, double* grid_arr) const {
            for (len_type j = 0; j < interior_len; ++j) {
                grid_arr[interior_map[j]] = dense[j];
            }
        }

        ~FieldDOF() { deallocate(); }

        FieldDOF(FieldDOF const& o)
            : grid_values{o.grid_values}, dframe_values{o.dframe_values},
              grid_len{o.grid_len}, interior_len{o.interior_len},
              interior_map{o.interior_map},
              eval_rhs{o.eval_rhs},
              update_boundaries{o.update_boundaries},
              field_group{o.field_group},
              u_n{nullptr}, F_u{nullptr}, residual{nullptr},
              saved{nullptr} {
            if (o.interior_len > 0) {
                allocate(o.interior_len);
                if (o.u_n) jfnk::copy_vec(o.u_n, u_n, interior_len);
                if (o.F_u) jfnk::copy_vec(o.F_u, F_u, interior_len);
                if (o.residual) jfnk::copy_vec(o.residual, residual, interior_len);
                if (o.saved) jfnk::copy_vec(o.saved, saved, interior_len);
            }
        }
        FieldDOF& operator=(FieldDOF const& o) {
            if (this != &o) {
                deallocate();
                grid_values = o.grid_values; dframe_values = o.dframe_values;
                grid_len = o.grid_len; interior_len = o.interior_len;
                interior_map = o.interior_map;
                eval_rhs = o.eval_rhs;
                update_boundaries = o.update_boundaries;
                field_group = o.field_group;
                if (o.interior_len > 0) {
                    allocate(o.interior_len);
                    if (o.u_n) jfnk::copy_vec(o.u_n, u_n, interior_len);
                    if (o.F_u) jfnk::copy_vec(o.F_u, F_u, interior_len);
                    if (o.residual) jfnk::copy_vec(o.residual, residual, interior_len);
                    if (o.saved) jfnk::copy_vec(o.saved, saved, interior_len);
                }
            }
            return *this;
        }
        FieldDOF(FieldDOF&& o) noexcept
            : grid_values{o.grid_values}, dframe_values{o.dframe_values},
              grid_len{o.grid_len}, interior_len{o.interior_len},
              interior_map{std::move(o.interior_map)},
              eval_rhs{std::move(o.eval_rhs)},
              update_boundaries{std::move(o.update_boundaries)},
              field_group{o.field_group},
              u_n{o.u_n}, F_u{o.F_u}, residual{o.residual},
              saved{o.saved} {
            o.u_n = nullptr; o.F_u = nullptr;
            o.residual = nullptr; o.saved = nullptr; o.interior_len = 0;
        }
        FieldDOF& operator=(FieldDOF&& o) noexcept {
            if (this != &o) {
                deallocate();
                grid_values = o.grid_values; dframe_values = o.dframe_values;
                grid_len = o.grid_len; interior_len = o.interior_len;
                interior_map = std::move(o.interior_map);
                eval_rhs = std::move(o.eval_rhs);
                update_boundaries = std::move(o.update_boundaries);
                field_group = o.field_group;
                u_n = o.u_n; F_u = o.F_u; residual = o.residual;
                saved = o.saved;
                o.u_n = nullptr; o.F_u = nullptr;
                o.residual = nullptr; o.saved = nullptr; o.interior_len = 0;
            }
            return *this;
        }
    };

    //! Build the interior index map for a 2D grid.
    //! Interior indices: (ix + iy*W) for ix in [B, W-B), iy in [B, H-B).
    static std::vector<len_type> build_interior_map_2d(
            len_type W, len_type H, len_type B) {
        std::vector<len_type> map;
        map.reserve((W - 2 * B) * (H - 2 * B));
        for (len_type iy = B; iy < H - B; ++iy) {
            for (len_type ix = B; ix < W - B; ++ix) {
                map.push_back(ix + iy * W);
            }
        }
        return map;
    }

    //! Build the interior index map for a 3D grid.
    //! Interior indices: (ix + iy*W + iz*W*H) for each axis in [B, dim-B).
    static std::vector<len_type> build_interior_map_3d(
            len_type W, len_type H, len_type D, len_type B) {
        std::vector<len_type> map;
        map.reserve((W - 2 * B) * (H - 2 * B) * (D - 2 * B));
        for (len_type iz = B; iz < D - B; ++iz) {
            for (len_type iy = B; iy < H - B; ++iy) {
                for (len_type ix = B; ix < W - B; ++ix) {
                    map.push_back(ix + iy * W + iz * W * H);
                }
            }
        }
        return map;
    }

    //! Build interior index map, dispatching on grid dimensionality.
    //! @param dims  Grid dimensions array (length ndim, including boundaries).
    //! @param ndim  Number of spatial dimensions (1, 2, or 3).
    //! @param B     Boundary depth (BOUNDARY_DEPTH).
    static std::vector<len_type> build_interior_map(
            const len_type* dims, size_t ndim, len_type B) {
        if (ndim == 3) {
            return build_interior_map_3d(dims[0], dims[1], dims[2], B);
        } else if (ndim == 2) {
            return build_interior_map_2d(dims[0], dims[1], B);
        } else {
            // 1D fallback.
            std::vector<len_type> map;
            map.reserve(dims[0] - 2 * B);
            for (len_type ix = B; ix < dims[0] - B; ++ix) {
                map.push_back(ix);
            }
            return map;
        }
    }

    // =========================================================================
    // Solver state
    // =========================================================================

    mutable bool solve_triggered{false};         //!< Set once Newton solve is done this step.
    mutable int num_fields_expected{0};      //!< equation() calls per timestep.
    mutable int num_fields_registered{0};    //!< equation() calls received so far.
    mutable int next_field_group{0};         //!< Counter for field group tags.
    mutable std::vector<FieldDOF> dofs;   //!< Per-component DOF data.
    mutable len_type total_len{0};           //!< Sum of all DOF lengths.

    // GMRES workspace (allocated lazily).
    mutable std::vector<double*> krylov_V;
    mutable double* krylov_w{nullptr};
    mutable std::vector<double> H;
    mutable std::vector<double> gmres_g;
    mutable std::vector<double> gmres_cs;
    mutable std::vector<double> gmres_sn;
    mutable std::vector<double> gmres_y;
    mutable double* delta_u{nullptr};
    mutable bool workspace_allocated{false};

    // =========================================================================
    // Destructor: free GMRES workspace
    // =========================================================================

    ~SolverJFNK() {
        for (auto* p : krylov_V) delete[] p;
        krylov_V.clear();
        delete[] krylov_w; krylov_w = nullptr;
        delete[] delta_u;  delta_u = nullptr;
    }

    // =========================================================================
    // step: no-op — solution was written directly during Newton solve.
    // =========================================================================
    //
    // SolverFT's step does: u += dt * dframe. We do NOT want that.
    // The Newton solve has already written the converged u^{n+1} into
    // sys.values directly.
    //

    template <typename S>
    void step(S&&) const {
        // Reset for next timestep so equation() runs again.
        solve_triggered = false;
        num_fields_registered = 0;
    }

    // =========================================================================
    // form_expr_one: optimize expression (same as SolverFT)
    // =========================================================================

    template <size_t En, typename SS, typename S, typename E>
    auto form_expr_one(SS&&, std::pair<S, E> const& e) const {
        auto [sys, equation] = e;
        auto eq_opt = expr::transform::optimize(expr::apply_operators(equation));
        expr::prune::update(eq_opt);
        expr::printe(eq_opt, "JFNK scheme");
        return std::make_pair(sys, eq_opt);
    }

    // =========================================================================
    // equation: collect fields, trigger Newton on last one
    // =========================================================================

    //! Helper: register DOFs for a SCALAR field.
    //! Block<scalar_t> has a single scalar_t* values array.
    template <size_t Z, typename S, typename E>
    void register_scalar_field(
            std::pair<Variable<Z, symphas::ref<S>>, E>& r,
            int group) const
    {
        auto& sys = r.first.get();

        FieldDOF dof;
        dof.grid_values = reinterpret_cast<double*>(sys.values);
        dof.dframe_values = reinterpret_cast<double*>(sys.dframe.values);
        dof.grid_len = sys.len;
        dof.field_group = group;

        // Build interior map based on grid dims and BOUNDARY_DEPTH.
        constexpr size_t ndim = sizeof(sys.dims) / sizeof(sys.dims[0]);
        dof.interior_map = build_interior_map(sys.dims, ndim, BOUNDARY_DEPTH);
        dof.interior_len = static_cast<len_type>(dof.interior_map.size());
        dof.allocate(dof.interior_len);

        // Callback: evaluate F(u) into sys.dframe.
        // During Newton iterations, noise must NOT be regenerated — the
        // stochastic forcing is fixed once per timestep.  We exclude
        // OpSymbolicEval (which wraps noise) from the prune::update so
        // cached noise grids are read but not overwritten.
        dof.eval_rhs = [&r]() {
            auto& sys_ref = r.first.get();
            {
                using Expr = std::remove_reference_t<decltype(r.second)>;
                expr::eval_handler_type<Expr> handler;
                expr::prune::update<expr::not_<expr::or_<
                        expr::matches_series,
                        expr::matches_symbolic_eval>>>(
                        r.second, handler);
            }
            expr::result(r.second, sys_ref.dframe);
        };

        // Callback: refresh boundaries.
        dof.update_boundaries = [&r]() {
            r.first.get().update();
        };

        dofs.push_back(std::move(dof));
    }

    //! Helper: register DOFs for a VECTOR field (any_vector_t<T,N>).
    //! MultiBlock<N,T> stores N separate T* arrays: values[0]..values[N-1].
    //! Similarly, sys.dframe (BoundaryGrid<vector_t<N>,N>) stores N arrays.
    //! We register N per-component DOFs, all sharing the same eval_rhs.
    template <size_t Z, typename S, typename E, typename T, size_t VN>
    void register_vector_field(
            std::pair<Variable<Z, symphas::ref<S>>, E>& r,
            int group,
            any_vector_t<T, VN>*) const
    {
        auto& sys = r.first.get();

        // Build the shared eval_rhs callback.
        auto shared_eval = std::make_shared<std::function<void()>>(
            [&r]() {
                auto& sys_ref = r.first.get();
                {
                    using Expr = std::remove_reference_t<decltype(r.second)>;
                    expr::eval_handler_type<Expr> handler;
                    expr::prune::update<expr::not_<expr::or_<
                            expr::matches_series,
                            expr::matches_symbolic_eval>>>(
                            r.second, handler);
                }
                expr::result(r.second, sys_ref.dframe);
            }
        );

        // Shared boundary update callback.
        auto shared_update = std::make_shared<std::function<void()>>(
            [&r]() {
                r.first.get().update();
            }
        );

        // Build interior map once (all components share same grid dims).
        constexpr size_t ndim = sizeof(sys.dims) / sizeof(sys.dims[0]);
        auto imap = build_interior_map(sys.dims, ndim, BOUNDARY_DEPTH);

        for (int comp = 0; comp < static_cast<int>(VN); ++comp) {
            FieldDOF dof;
            // MultiBlock<VN,T>: sys.values is T* values[VN].
            // Each values[comp] is a contiguous scalar_t[sys.len] array.
            dof.grid_values = reinterpret_cast<double*>(sys.values[comp]);
            dof.dframe_values = reinterpret_cast<double*>(
                    sys.dframe.values[comp]);
            dof.grid_len = sys.len;
            dof.field_group = group;
            dof.interior_map = imap;
            dof.interior_len = static_cast<len_type>(imap.size());
            dof.allocate(dof.interior_len);

            // Share the eval_rhs callback (only called once per group).
            auto se = shared_eval;
            dof.eval_rhs = [se]() { (*se)(); };

            auto su = shared_update;
            dof.update_boundaries = [su]() { (*su)(); };

            dofs.push_back(std::move(dof));
        }
    }

    //! Detect value type and dispatch to scalar or vector registration.
    template <size_t Z, typename S, typename E>
    void register_field_dofs(
            std::pair<Variable<Z, symphas::ref<S>>, E>& r,
            int group) const
    {
        auto& sys = r.first.get();
        using value_type = typename grid::value_type_of<
                std::remove_reference_t<decltype(sys)>>::type;
        constexpr int nc = jfnk::field_traits<value_type>::num_components;

        if constexpr (nc == 1) {
            register_scalar_field(r, group);
        } else {
            // Pass a typed null pointer for tag dispatch.
            register_vector_field(r, group,
                    static_cast<value_type*>(nullptr));
        }
    }

    template <size_t Z, typename S, typename E>
    inline void equation(std::pair<Variable<Z, symphas::ref<S>>, E>& r) const {
        // Skip if Newton solve already done this timestep.
        if (solve_triggered) return;

        // First equation() call: reset.
        if (num_fields_registered == 0) {
            dofs.clear();
            next_field_group = 0;
            solve_triggered = false;
        }

        auto& sys = r.first.get();

        // Ensure expression sub-data (noise grids, etc.) are allocated.
        r.second.allocate();

        // Evaluate F(u^current) into sys.dframe (real-space BoundaryGrid).
        // This uses not_<matches_series> so noise IS updated here — this
        // generates the stochastic forcing for this timestep.
        {
            expr::eval_handler_type<E> handler;
            expr::prune::update<expr::not_<expr::matches_series>>(
                    r.second, handler);
        }
        expr::result(r.second, sys.dframe);

        // Register DOFs for this field (1 for scalar, D for vector).
        register_field_dofs(r, next_field_group);
        next_field_group++;
        num_fields_registered++;

        // If all fields registered, run Newton-GMRES.
        if (num_fields_registered >= num_fields_expected) {
            perform_newton_solve();
            solve_triggered = true;
            num_fields_registered = 0;
        }
    }

    template <typename S, typename E>
    inline void equation(std::pair<S, E>* r, len_type len) const {
        for (iter_type i = 0; i < len; ++i) {
            equation(r[i]);
        }
    }

    // =========================================================================
    // evaluate_one: for provisional variables
    // =========================================================================

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

    // =========================================================================
    // make_solver
    // =========================================================================

    static auto make_solver(symphas::problem_parameters_type const& parameters) {
        if (parameters.length()) {
            double h = parameters.get_interval_data()[0].at(Axis::X).width();
            size_t dim = parameters.get_dimension();

            len_type* dims = new len_type[dim];
            for (iter_type i = 0; i < static_cast<iter_type>(dim); ++i) {
                Axis side = symphas::index_to_axis(i);
                dims[i] = parameters.get_interval_data()[0].at(side).get_count()
                        + 2 * BOUNDARY_DEPTH;
            }

            auto s = this_type{dims, h};
            s.num_fields_expected = parameters.length();
            s.solve_triggered = false;
            s.num_fields_registered = 0;
            s.next_field_group = 0;
            s.total_len = 0;
            s.krylov_w = nullptr;
            s.delta_u = nullptr;
            s.workspace_allocated = false;
            delete[] dims;
            return s;
        } else {
            auto s = this_type{grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3)};
            s.num_fields_expected = 0;
            s.solve_triggered = false;
            s.num_fields_registered = 0;
            s.next_field_group = 0;
            s.total_len = 0;
            s.krylov_w = nullptr;
            s.delta_u = nullptr;
            s.workspace_allocated = false;
            return s;
        }
    }

private:

    // =========================================================================
    // GMRES Workspace Allocation
    // =========================================================================

    void ensure_workspace(len_type new_total) const {
        if (workspace_allocated && total_len == new_total) return;

        // Free old.
        for (auto* p : krylov_V) delete[] p;
        krylov_V.clear();
        delete[] krylov_w; krylov_w = nullptr;
        delete[] delta_u; delta_u = nullptr;

        total_len = new_total;
        int m = jfnk::GMRES_M;

        krylov_V.resize(m + 1);
        for (int i = 0; i <= m; ++i) {
            krylov_V[i] = new double[total_len]();
        }
        krylov_w = new double[total_len]();
        delta_u = new double[total_len]();

        H.resize(static_cast<size_t>((m + 1) * m), 0.0);
        gmres_g.resize(m + 1, 0.0);
        gmres_cs.resize(m, 0.0);
        gmres_sn.resize(m, 0.0);
        gmres_y.resize(m, 0.0);

        workspace_allocated = true;
    }

    // =========================================================================
    // Newton-GMRES Core
    // =========================================================================

    void perform_newton_solve() const {
        if (dofs.empty()) return;

        // On the very first timestep, dt=0 (framework sets dt after step()).
        // G(u) = u - u_n - 0*F(u) = 0, so skip immediately.
        if (dt <= 0.0) return;

        // Compute total interior DOFs across all components.
        len_type tot = 0;
        for (auto const& d : dofs) tot += d.interior_len;
        ensure_workspace(tot);

        // Step 1: Gather u_n and F(u_n) from grid arrays into dense DOF buffers.
        // F(u_n) was already evaluated into sys.dframe by equation() calls.
        for (auto& d : dofs) {
            d.gather(d.grid_values, d.u_n);    // u_n[j] = grid[interior_map[j]]
            d.gather(d.dframe_values, d.F_u);  // F_u[j] = dframe[interior_map[j]]
        }

        // Step 2: Compute initial Newton residual.
        // G(u) = u - u_n - dt*F(u).  Since u = u_n here, G = -dt*F(u_n).
        compute_newton_residual();
        double res0 = compute_residual_norm();
        double res = res0;

#ifdef JFNK_VERBOSE
        fprintf(stdout, "[JFNK] initial ||G|| = %.6e  (%d interior DOFs)\n",
                res0, static_cast<int>(tot));
#endif

        // NOTE: Do NOT early-exit here based on absolute tolerance.
        // The initial residual G(u_old) = -dt*F(u_old) scales with dt,
        // so small dt will always have small initial residual — but the
        // Newton solve must still proceed to apply the update.
        // Convergence is checked AFTER at least one Newton iteration.

        // Step 3: Newton iteration.
        for (int k = 0; k < jfnk::MAX_NEWTON_ITERS; ++k) {
            // Solve J * delta_u = -G(u) via GMRES.
            jfnk::zero_vec(delta_u, total_len);
            double gmres_tol = jfnk::GMRES_RTOL * res;
            int gmres_iters = 0;
            bool gmres_ok = solve_gmres(gmres_tol, gmres_iters);

            // Armijo backtracking line search: accept the largest alpha in
            // {1, 1/2, 1/4, ...} for which ||G(u + alpha*du)|| < (1 - 1e-4*alpha)*||G(u)||.
            // Without this, full Newton on PFC (6th-order, ill-conditioned)
            // diverges immediately because the GMRES direction overshoots.
            constexpr int LS_MAX = 20;
            constexpr double LS_SIGMA = 1e-4;
            double alpha = 1.0;
            double res_trial = res;
            double res_prev = res;
            bool ls_ok = false;
            for (int ls = 0; ls < LS_MAX; ++ls) {
                scatter_update(delta_u, alpha);
                update_all_boundaries();
                evaluate_all_rhs();
                for (auto& d : dofs) {
                    d.gather(d.dframe_values, d.F_u);
                }
                compute_newton_residual();
                res_trial = compute_residual_norm();
                if (res_trial <= (1.0 - LS_SIGMA * alpha) * res_prev) {
                    ls_ok = true;
                    break;
                }
                // Reject: undo step, halve alpha, retry.
                scatter_update(delta_u, -alpha);
                alpha *= 0.5;
            }
            if (!ls_ok) {
                // No descent direction found. Keep the last (rejected) trial
                // applied so state advances; warn and break.
                scatter_update(delta_u, alpha);
                update_all_boundaries();
                evaluate_all_rhs();
                for (auto& d : dofs) {
                    d.gather(d.dframe_values, d.F_u);
                }
                compute_newton_residual();
                res_trial = compute_residual_norm();
#ifdef JFNK_VERBOSE
                fprintf(stdout, "[JFNK] newton %2d: LINE SEARCH FAILED, "
                        "accepting alpha=%.3e, ||G||=%.6e\n",
                        k + 1, alpha, res_trial);
#endif
            }
            res = res_trial;

#ifdef JFNK_VERBOSE
            fprintf(stdout, "[JFNK] newton %2d: ||G|| = %.6e, alpha = %.3e, "
                    "GMRES iters = %d (%s)\n",
                    k + 1, res, alpha, gmres_iters,
                    gmres_ok ? "conv" : "maxiter");
#endif

            if (res < jfnk::NEWTON_ATOL ||
                res < jfnk::NEWTON_RTOL * res0) {
#ifdef JFNK_VERBOSE
                fprintf(stdout, "[JFNK] converged in %d Newton iterations\n",
                        k + 1);
#endif
                return;
            }
        }

        fprintf(SYMPHAS_WARN, "[JFNK] did NOT converge after %d Newton "
                "iterations (||G|| = %.6e)\n",
                jfnk::MAX_NEWTON_ITERS, res);
    }

    // =========================================================================
    // F(u) evaluation and boundary refresh
    // =========================================================================

    //! Re-evaluate all RHS expressions. Each callback calls
    //! expr::result(expr, sys.dframe), writing F(u) into sys.dframe.
    //! For vector fields, multiple DOFs share the same eval_rhs callback
    //! via field_group — we call it only once per group.
    void evaluate_all_rhs() const {
        int last_group = -1;
        for (auto& d : dofs) {
            if (d.field_group != last_group) {
                d.eval_rhs();
                last_group = d.field_group;
            }
        }
    }

    //! Refresh boundary ghost cells for all fields.
    //! Called after modifying grid interior data during Newton iterations.
    //! Only calls update once per field group (vector components share it).
    void update_all_boundaries() const {
        int last_group = -1;
        for (auto& d : dofs) {
            if (d.field_group != last_group) {
                if (d.update_boundaries) d.update_boundaries();
                last_group = d.field_group;
            }
        }
    }

    // =========================================================================
    // Newton Residual
    // =========================================================================

    //! Compute G(u) = u - u_n - dt * F(u) per component.
    //! u is read from grid_values via interior_map, while u_n and F_u are
    //! dense DOF-indexed arrays (already gathered from grid interior).
    void compute_newton_residual() const {
        for (auto& d : dofs) {
            for (len_type j = 0; j < d.interior_len; ++j) {
                double u_curr = d.grid_values[d.interior_map[j]];
                d.residual[j] = u_curr - d.u_n[j] - dt * d.F_u[j];
            }
        }
    }

    //! RMS norm of Newton residual across all interior DOFs.
    double compute_residual_norm() const {
        double sq = 0.0;
        len_type tot = 0;
        for (auto const& d : dofs) {
            sq += jfnk::dot_product(d.residual, d.residual, d.interior_len);
            tot += d.interior_len;
        }
        return std::sqrt(sq / tot);
    }

    // =========================================================================
    // GMRES(m) with restart
    // =========================================================================

    //! Solve J * x = -G(u) via restarted GMRES(m).
    //! The solution x is accumulated into delta_u.
    bool solve_gmres(double tol, int& total_iters) const {
        int m = jfnk::GMRES_M;
        len_type Ntot = total_len;
        total_iters = 0;

        // Outer restart loop.
        for (int restart = 0;
             restart < jfnk::MAX_GMRES_ITERS / m + 1 &&
             total_iters < jfnk::MAX_GMRES_ITERS;
             ++restart) {

            // Compute r = b - A*x_current = -G(u) - J*delta_u.
            // On first pass, delta_u = 0, so r = -G(u).
            // On restarts, we recompute the residual properly.
            if (restart == 0) {
                gather_neg_residual(krylov_V[0]);
            } else {
                // r = b - A*delta_u: compute J*delta_u, then r = b - J*delta_u.
                jacobian_free_matvec(delta_u, krylov_V[0]);
                // krylov_V[0] now has J*delta_u; we want b - J*delta_u.
                // b = -G(u_current), but G was computed at the start of Newton.
                // Re-gather -G and subtract J*delta_u.
                // Actually: for GMRES restart, we need to recompute exactly.
                // Simple approach: gather -G and subtract matvec result.
                double* tmp = krylov_w;  // reuse workspace temporarily.
                gather_neg_residual(tmp);
                for (len_type i = 0; i < Ntot; ++i) {
                    krylov_V[0][i] = tmp[i] - krylov_V[0][i];
                }
            }

            double beta = jfnk::norm2(krylov_V[0], Ntot);
            if (beta < 1e-30) return true;

            jfnk::scale_vec(1.0 / beta, krylov_V[0], Ntot);

            std::fill(H.begin(), H.end(), 0.0);
            std::fill(gmres_g.begin(), gmres_g.end(), 0.0);
            gmres_g[0] = beta;
            std::fill(gmres_cs.begin(), gmres_cs.end(), 0.0);
            std::fill(gmres_sn.begin(), gmres_sn.end(), 0.0);

            int j;
            for (j = 0; j < m && total_iters < jfnk::MAX_GMRES_ITERS;
                 ++j, ++total_iters) {

                // Arnoldi: w = J * V[j].
                jacobian_free_matvec(krylov_V[j], krylov_w);

                // Modified Gram-Schmidt.
                for (int i = 0; i <= j; ++i) {
                    H[i + j * (m + 1)] =
                        jfnk::dot_product(krylov_w, krylov_V[i], Ntot);
                    jfnk::axpy(-H[i + j * (m + 1)], krylov_V[i], krylov_w, Ntot);
                }
                double h_jp1_j = jfnk::norm2(krylov_w, Ntot);
                H[(j + 1) + j * (m + 1)] = h_jp1_j;

                if (h_jp1_j < 1e-30) { j++; break; }

                jfnk::copy_vec(krylov_w, krylov_V[j + 1], Ntot);
                jfnk::scale_vec(1.0 / h_jp1_j, krylov_V[j + 1], Ntot);

                // Apply previous Givens rotations.
                for (int i = 0; i < j; ++i) {
                    jfnk::apply_givens(H[i + j * (m + 1)],
                                       H[(i + 1) + j * (m + 1)],
                                       gmres_cs[i], gmres_sn[i]);
                }

                // Compute new Givens rotation.
                jfnk::compute_givens(H[j + j * (m + 1)],
                                     H[(j + 1) + j * (m + 1)],
                                     gmres_cs[j], gmres_sn[j]);

                jfnk::apply_givens(H[j + j * (m + 1)],
                                   H[(j + 1) + j * (m + 1)],
                                   gmres_cs[j], gmres_sn[j]);
                jfnk::apply_givens(gmres_g[j], gmres_g[j + 1],
                                   gmres_cs[j], gmres_sn[j]);

                if (std::abs(gmres_g[j + 1]) < tol) { j++; break; }
            }

            // Back-substitution.
            int kk = j;
            for (int i = kk - 1; i >= 0; --i) {
                gmres_y[i] = gmres_g[i];
                for (int jj = i + 1; jj < kk; ++jj) {
                    gmres_y[i] -= H[i + jj * (m + 1)] * gmres_y[jj];
                }
                double diag = H[i + i * (m + 1)];
                if (std::abs(diag) > 1e-30) {
                    gmres_y[i] /= diag;
                }
            }

            // Accumulate delta_u += V * y.
            for (int i = 0; i < kk; ++i) {
                jfnk::axpy(gmres_y[i], krylov_V[i], delta_u, Ntot);
            }

            // Check convergence.
            if (kk <= m && std::abs(gmres_g[kk]) < tol) {
                break;  // Converged.
            }
            if (total_iters >= jfnk::MAX_GMRES_ITERS) {
                break;  // Hit iteration limit.
            }
            // Otherwise: restart with updated residual (next pass of outer loop).
        }

        return total_iters < jfnk::MAX_GMRES_ITERS;
    }

    // =========================================================================
    // Jacobian-Free Matrix-Vector Product
    // =========================================================================

    //! Compute J*v ≈ (G(u + ε*v) - G(u)) / ε.
    //!
    //! Since G(u) = u - u_n - dt*F(u), J = I - dt*(dF/du), so:
    //!   J*v = v - dt * (F(u + ε*v) - F(u)) / ε
    //!
    //! F(u) at the current iterate is stored in dof.F_u (dense, interior-indexed).
    //!
    //! Pre-allocated dof.saved buffers are used to save/restore grid interior
    //! data (no heap allocation during matvec).
    void jacobian_free_matvec(const double* v, double* Jv) const {
        len_type Ntot = total_len;

        // If v is (near-)zero, J*v = v (the identity part dominates).
        double v_norm = jfnk::norm2(v, Ntot);
        if (v_norm < 1e-30) {
            jfnk::copy_vec(v, Jv, Ntot);
            return;
        }

        // Compute ε using Walker-Pernice formula (over interior DOFs only).
        double u_norm_sq = 0.0;
        for (auto const& d : dofs) {
            // Gather interior u values for norm computation.
            for (len_type j = 0; j < d.interior_len; ++j) {
                double val = d.grid_values[d.interior_map[j]];
                u_norm_sq += val * val;
            }
        }
        double u_norm = std::sqrt(u_norm_sq);
        double eps = jfnk::EPSILON_BASE * (1.0 + u_norm) / v_norm;

        // Save current grid interior values (using pre-allocated buffers).
        for (auto& d : dofs) {
            d.gather(d.grid_values, d.saved);
        }

        // Perturb: grid_interior ← grid_interior + ε*v.
        len_type offset = 0;
        for (auto& d : dofs) {
            for (len_type j = 0; j < d.interior_len; ++j) {
                d.grid_values[d.interior_map[j]] += eps * v[offset + j];
            }
            offset += d.interior_len;
        }

        // Refresh boundaries after perturbation.
        update_all_boundaries();

        // Evaluate F(u + ε*v).
        evaluate_all_rhs();

        // Compute J*v = v - dt * (F(u+εv) - F(u)) / ε.
        // F(u+εv) is now in dframe (grid layout). F(u) is in dof.F_u (dense).
        // We gather F(u+εv) from interior and subtract F_u.
        offset = 0;
        for (auto& d : dofs) {
            for (len_type j = 0; j < d.interior_len; ++j) {
                double F_perturbed = d.dframe_values[d.interior_map[j]];
                double dF = (F_perturbed - d.F_u[j]) / eps;
                Jv[offset + j] = v[offset + j] - dt * dF;
            }
            offset += d.interior_len;
        }

        // Restore grid interior values from saved buffers.
        for (auto& d : dofs) {
            d.scatter_set(d.saved, d.grid_values);
        }

        // Restore boundaries to unperturbed state.
        update_all_boundaries();
    }

    // =========================================================================
    // Scatter/Gather
    // =========================================================================

    //! Apply delta_u to all field grids: u_interior += alpha * du.
    //! Scatters flat DOF vector entries to grid interior cells.
    void scatter_update(const double* du, double alpha) const {
        len_type offset = 0;
        for (auto& d : dofs) {
            for (len_type j = 0; j < d.interior_len; ++j) {
                d.grid_values[d.interior_map[j]] += alpha * du[offset + j];
            }
            offset += d.interior_len;
        }
    }

    //! Gather -G(u) into flat vector (dense, interior-indexed).
    void gather_neg_residual(double* flat) const {
        len_type offset = 0;
        for (auto const& d : dofs) {
            for (len_type j = 0; j < d.interior_len; ++j) {
                flat[offset + j] = -d.residual[j];
            }
            offset += d.interior_len;
        }
    }

END_SOLVER

ASSOCIATE_SOLVER_SYSTEM_TYPE(SolverJFNK, SolverSystemFD)
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverJFNK, ProvisionalSystemFD)
SYMPHAS_SOLVER_ALL_SUPPORTED(SolverJFNK)
