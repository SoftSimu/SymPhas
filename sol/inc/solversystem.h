
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
 * PURPOSE: Manages a group of provisional systems, used by the phase field
 * model in order to manage the provisional variables.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundarysystem.h"
#include "spslibfftw.h"
#include "spsmpi.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

//! The default phase field system.
/*!
 * The representation of a phase field system, storing the values of the
 * order parameter that is defined in a phase field problem.
 *
 * Unless explicitly specified, this phase field system type will be
 * used by the solver. It does not manage boundaries or other data.
 *cid
 * \tparam T The order parameter type.
 * \tparam D The order parameter dimension.
 */
template <typename T, size_t D>
using SolverSystem = System<T, D>;

template <typename T, size_t D>
struct SolverSystemFD : BoundarySystem<T, D> {
  using BoundarySystem<T, D>::dims;

  BoundaryGrid<T, D> dframe;  // the working grid for the solver

#if defined(USING_MPI) && defined(SYMPHAS_MPI_LOCAL_STORAGE)
  // SYMPHAS_MPI_LOCAL_STORAGE: each rank allocates only its sub-region
  // of the global grid (with halo). Phase 1-4 of the local-storage
  // refactor live here, on the actually-registered FD system class
  // (see ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE in solverft.h).
  symphas::parallel::domain_info<D> dinfo;

 private:
  // Compute global+halo dims from the user-facing vdata (matches what
  // the legacy global-storage parent ctor would produce).
  static grid::dim_list compute_global_dims(
      symphas::interval_data_type const& vdata) {
    auto v = vdata;
    if (params::extend_boundary) {
      for (auto& [_, iv] : v) {
        iv.set_count(iv.get_count() + 2 * BOUNDARY_DEPTH);
        iv.interval_to_domain();
      }
    }
    return symphas::grid_info(v).get_dims();
  }

  // domain_info wants a fixed-size array; convert from dim_list and tag
  // the result as local-storage (so halo + iterable_domain consumers know
  // grid.dims is local, not global).
  static symphas::parallel::domain_info<D> make_dinfo(
      grid::dim_list const& gd) {
    len_type fixed[D];
    for (size_t i = 0; i < D; ++i) fixed[i] = gd[i];
    auto di = symphas::parallel::domain_info<D>(fixed, BOUNDARY_DEPTH);
    di.local_storage = true;
    return di;
  }

  // Build a vdata describing this rank's local sub-region INCLUDING halo
  // on each side. Total per-axis count = local_n + 2*BOUNDARY_DEPTH;
  // the physical interval is positioned at the rank's slice of the user
  // domain, extended by BOUNDARY_DEPTH cells of halo.
  static symphas::interval_data_type make_local_vdata(
      symphas::interval_data_type const& vdata,
      symphas::parallel::domain_info<D> const& di) {
    if constexpr (D < 2) {
      return vdata;
    } else {
      if (di.num_ranks <= 1) return vdata;
      symphas::interval_data_type local = vdata;
      Axis xa = symphas::index_to_axis(0);
      Axis ya = symphas::index_to_axis(1);
      auto& xv = local.at(xa);
      auto& yv = local.at(ya);
      double xh = xv.width();
      double yh = yv.width();
      len_type local_total_x = di.local_nx + 2 * BOUNDARY_DEPTH;
      len_type local_total_y = di.local_ny + 2 * BOUNDARY_DEPTH;
      double x_left = xv.domain_left() +
                      (di.local_x_start - BOUNDARY_DEPTH) * xh;
      double x_right = x_left + (local_total_x - 1) * xh;
      double y_left = yv.domain_left() +
                      (di.local_y_start - BOUNDARY_DEPTH) * yh;
      double y_right = y_left + (local_total_y - 1) * yh;
      xv.set_domain(x_left, x_right, xh);
      xv.set_interval(x_left, x_right);
      yv.set_domain(y_left, y_right, yh);
      yv.set_interval(y_left, y_right);
      return local;
    }
  }

  // Replace the four periodic boundary updaters with MPI exchange
  // updaters carrying the rank-local domain_info. Cart topology has
  // periods={1,1}; every rank has neighbors on every side. Side::TOP
  // is the one that actually triggers exchange_halos at update time;
  // the other three are no-ops, gating exchange_halos to fire exactly
  // once per step.
  void rewire_mpi_boundaries() {
    if constexpr (D < 2) return;
    if (dinfo.num_ranks <= 1) return;
    int rk = dinfo.rank;
    constexpr iter_type TOP_IDX = static_cast<iter_type>(Side::TOP);
    constexpr iter_type BOT_IDX = static_cast<iter_type>(Side::BOTTOM);
    constexpr iter_type LEFT_IDX = static_cast<iter_type>(Side::LEFT);
    constexpr iter_type RIGHT_IDX = static_cast<iter_type>(Side::RIGHT);
    auto& bg = static_cast<BoundaryGroup<T, D>&>(*this);

    auto install_mpi = [&](iter_type idx, int neighbor) {
      delete bg.boundaries[idx];
      bg.types[idx] = BoundaryType::MPI;
      auto* b = new grid::BoundaryApplied<T, D - 1, BoundaryType::MPI>(
          neighbor, rk);
      // Attach the local-storage dinfo so the halo updater uses it
      // directly instead of constructing one from grid.dims.
      b->dinfo_storage = new symphas::parallel::domain_info<D>(dinfo);
      bg.boundaries[idx] = b;
    };

    install_mpi(TOP_IDX, dinfo.neighbor_above);
    install_mpi(BOT_IDX, dinfo.neighbor_below);
    if (dinfo.dims_cart[0] > 1) {
      install_mpi(LEFT_IDX, dinfo.neighbor_left);
      install_mpi(RIGHT_IDX, dinfo.neighbor_right);
    }
  }

  // Delegating ctor: builds dinfo from global dims, then constructs
  // parent with rank-local vdata so the BoundaryGrid storage is sized
  // (local_n + 2*bdepth) per axis instead of global N.
  SolverSystemFD(symphas::init_data_type const& tdata,
                 symphas::interval_data_type const& vdata,
                 symphas::b_data_type const& bdata, size_t id,
                 grid::dim_list const& global_dims)
      : BoundarySystem<T, D>(
            tdata,
            make_local_vdata(vdata, make_dinfo(global_dims)),
            bdata, id),
        dframe{dims},
        dinfo{make_dinfo(global_dims)} {
    if (dinfo.num_ranks > 1) {
      fprintf(stderr,
              "[SolverSystemFD::ctor LOCAL_STORAGE] rank=%d global_dims=%dx%d "
              "local_nx=%d local_ny=%d after-ctor dims=%dx%d\n",
              dinfo.rank, (int)global_dims[0], (int)global_dims[1],
              (int)dinfo.local_nx, (int)dinfo.local_ny,
              (int)dims[0], (int)dims[1]);
    }
    rewire_mpi_boundaries();
  }

 public:
  SolverSystemFD(symphas::init_data_type const& tdata,
                 symphas::interval_data_type const& vdata,
                 symphas::b_data_type const& bdata, size_t id = 0)
      : SolverSystemFD(tdata, vdata, bdata, id,
                       compute_global_dims(vdata)) {}
  SolverSystemFD() : BoundarySystem<T, D>(), dframe{0}, dinfo{} {}
#else
  SolverSystemFD(symphas::init_data_type const& tdata,
                 symphas::interval_data_type const& vdata,
                 symphas::b_data_type const& bdata, size_t id = 0)
      : BoundarySystem<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFD() : BoundarySystem<T, D>(), dframe{0} {}
#endif
};

template <typename T, size_t D>
struct SolverSystemFDwSD : RegionalSystem<T, D> {
  using RegionalSystem<T, D>::dims;

  RegionalGrid<T, D> dframe;  // the working grid for the solver
  SolverSystemFDwSD(symphas::init_data_type const& tdata,
                    symphas::interval_data_type const& vdata,
                    symphas::b_data_type const& bdata, size_t id = 0)
      : RegionalSystem<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFDwSD() : RegionalSystem<T, D>(), dframe{0} {}

  inline void update(iter_type index, double time) {
    RegionalSystem<T, D>::update(index, time);
    dframe.adjust(RegionalSystem<T, D>::region);
  }
};

#ifdef USING_MPI

// =====================================================================
// BUILD-TIME GUARD — SolverSystemFDwSDMPI is structurally unsound.
//
// SolverSystemFDwSDMPI combines:
//   * RegionalSystemMPI<T,D>   (MPI-decomposed sparse region storage)
//   * RegionalGrid<T,D>        (single-process sparse working frame)
//
// These assumptions conflict: RegionalGrid assumes one process owns the
// full sparse domain, while RegionalSystemMPI decomposes it across ranks.
// In practice, using this class produces access violations and
// MPI_Bcast: Invalid buffer pointer errors.
//
// The `static_assert` below is a template-instantiation guard. The class
// declaration is kept so existing headers still parse, but any attempt
// to instantiate it (e.g. through SolverFT's selectable-system machinery)
// triggers a compile-time error at the instantiation site. The registration
// macro ASSOCIATE_SELECTABLE_SOLVER_SYSTEM_TYPE(SolverFT, SolverSystemFDwSDMPI)
// in examples/solvers/solverft.h is also intentionally omitted.
// =====================================================================

namespace symphas::internal {
template <typename, size_t>
struct solver_system_fdwsdmpi_forbidden : std::false_type {};
}  // namespace symphas::internal

template <typename T, size_t D>
struct SolverSystemFDwSDMPI : RegionalSystemMPI<T, D> {
  static_assert(symphas::internal::solver_system_fdwsdmpi_forbidden<T, D>::value,
                "SolverSystemFDwSDMPI is forbidden: RegionalGrid is not "
                "compatible with MPI spatial decomposition. Use "
                "SolverSystemFD (variation 0) under MPI builds, or switch "
                "to a serial build to use SolverSystemFDwSD (variation 1).");

  using RegionalSystemMPI<T, D>::dims;
  using RegionalSystemMPI<T, D>::thr_info;

  RegionalGrid<T, D> dframe;  // the working grid for the solver
  SolverSystemFDwSDMPI(symphas::init_data_type const& tdata,
                       symphas::interval_data_type const& vdata,
                       symphas::b_data_type const& bdata, size_t id = 0)
      : RegionalSystemMPI<T, D>(tdata, vdata, bdata, id),
        dframe{(thr_info.is_in_node()) ? dims : nullptr} {}
  SolverSystemFDwSDMPI() : RegionalSystemMPI<T, D>(), dframe{0} {}

  inline void update(iter_type index, double time) {
    RegionalSystemMPI<T, D>::update(index, time);
    if (thr_info.is_in_node()) {
      dframe.adjust(RegionalSystemMPI<T, D>::region);
    }
  }
};

template <typename T, size_t D>
struct SolverSystemFDMPI : BoundarySystem<T, D> {
  using BoundarySystem<T, D>::dims;

  BoundaryGrid<T, D> dframe;
  symphas::parallel::domain_info<D> dinfo;

 private:
  // Compute global+halo dims from the user-facing vdata. This mirrors
  // what BoundarySystem's parent ctor would produce so we can construct
  // a global domain_info before delegating to the parent.
  static grid::dim_list compute_global_extended_dims(
      symphas::interval_data_type const& vdata) {
    auto v = vdata;
    if (params::extend_boundary) {
      for (auto& [_, iv] : v) {
        iv.set_count(iv.get_count() + 2 * BOUNDARY_DEPTH);
        iv.interval_to_domain();
      }
    }
    return symphas::grid_info(v).get_dims();
  }

  // Compute a vdata that describes ONLY this rank's local sub-region of
  // the global grid PLUS the boundary halo on each side. The result has
  // count = local_n + 2*BOUNDARY_DEPTH on each axis (matching what the
  // grid will allocate as `dims`), and the physical interval positioned
  // at the rank's slice of the user domain.
  //
  // Why +2*BOUNDARY_DEPTH: SymPhas's BoundaryGrid stores `dims` as the
  // total array size (interior + halo on both sides). The interval
  // `count` is interpreted by `grid::construct<BoundaryGrid>(vdata)` as
  // exactly the number of cells to allocate, so we must include halo.
  //
  // For non-MPI / single-rank fallback or when D < 2, returns vdata
  // unchanged.
  static symphas::interval_data_type make_local_vdata(
      symphas::interval_data_type const& vdata,
      symphas::parallel::domain_info<D> const& di) {
    if constexpr (D < 2) {
      return vdata;
    } else {
      if (di.num_ranks <= 1) return vdata;
      symphas::interval_data_type local = vdata;
      Axis xa = symphas::index_to_axis(0);
      Axis ya = symphas::index_to_axis(1);
      auto& xv = local.at(xa);
      auto& yv = local.at(ya);
      double xh = xv.width();
      double yh = yv.width();
      // Total local extent (interior + 2*halo) on each axis.
      len_type local_total_x = di.local_nx + 2 * BOUNDARY_DEPTH;
      len_type local_total_y = di.local_ny + 2 * BOUNDARY_DEPTH;
      // Position the slice in physical coords. domain_info splits the
      // GLOBAL INTERIOR; the interior starts at xv.left() in physical
      // coords (xv from vdata represents the user-facing interior).
      // Each cell is xh wide; we extend by BOUNDARY_DEPTH cells on each
      // side to cover halo too.
      double x_left = xv.domain_left() +
                      (di.local_x_start - BOUNDARY_DEPTH) * xh;
      double x_right = x_left + (local_total_x - 1) * xh;
      double y_left = yv.domain_left() +
                      (di.local_y_start - BOUNDARY_DEPTH) * yh;
      double y_right = y_left + (local_total_y - 1) * yh;
      xv.set_domain(x_left, x_right, xh);
      xv.set_interval(x_left, x_right);
      yv.set_domain(y_left, y_right, yh);
      yv.set_interval(y_left, y_right);
      return local;
    }
  }

#ifdef SYMPHAS_MPI_LOCAL_STORAGE
  // Delegating ctor: builds dinfo from the global vdata, then constructs
  // the parent with the rank-local vdata. Storage allocates only
  // (local_n + 2*bdepth) per axis.
  SolverSystemFDMPI(symphas::init_data_type const& tdata,
                    symphas::interval_data_type const& vdata,
                    symphas::b_data_type const& bdata, size_t id,
                    grid::dim_list const& global_dims)
      : BoundarySystem<T, D>(
            tdata,
            make_local_vdata(
                vdata,
                make_domain_info(global_dims)),
            bdata, id),
        dframe{dims},
        dinfo{make_domain_info(global_dims)} {
    fprintf(stderr,
            "[SOLVER_LOCALSTOR] rank=%d global_dims=%dx%d local_nx=%d "
            "local_ny=%d after-ctor dims=%dx%d\n",
            dinfo.rank, (int)global_dims[0], (int)global_dims[1],
            (int)dinfo.local_nx, (int)dinfo.local_ny,
            (int)dims[0], (int)dims[1]);
    rewire_mpi_boundaries();
  }

  // Helper: domain_info wants a fixed-size array. Convert from dim_list.
  static symphas::parallel::domain_info<D> make_domain_info(
      grid::dim_list const& gd) {
    len_type fixed[D];
    for (size_t i = 0; i < D; ++i) fixed[i] = gd[i];
    auto di = symphas::parallel::domain_info<D>(fixed, BOUNDARY_DEPTH);
    di.local_storage = true;
    return di;
  }

  void rewire_mpi_boundaries() {
    if constexpr (D < 2) return;
    if (dinfo.num_ranks <= 1) return;
    int rk = dinfo.rank;
    constexpr iter_type TOP_IDX = static_cast<iter_type>(Side::TOP);
    constexpr iter_type BOT_IDX = static_cast<iter_type>(Side::BOTTOM);
    constexpr iter_type LEFT_IDX = static_cast<iter_type>(Side::LEFT);
    constexpr iter_type RIGHT_IDX = static_cast<iter_type>(Side::RIGHT);
    auto& bg = static_cast<BoundaryGroup<T, D>&>(*this);

    delete bg.boundaries[TOP_IDX];
    delete bg.boundaries[BOT_IDX];
    bg.types[TOP_IDX] = BoundaryType::MPI;
    bg.types[BOT_IDX] = BoundaryType::MPI;
    bg.boundaries[TOP_IDX] =
        new grid::BoundaryApplied<T, D - 1, BoundaryType::MPI>(
            dinfo.neighbor_above, rk);
    bg.boundaries[BOT_IDX] =
        new grid::BoundaryApplied<T, D - 1, BoundaryType::MPI>(
            dinfo.neighbor_below, rk);
    if (dinfo.dims_cart[0] > 1) {
      delete bg.boundaries[LEFT_IDX];
      delete bg.boundaries[RIGHT_IDX];
      bg.types[LEFT_IDX] = BoundaryType::MPI;
      bg.types[RIGHT_IDX] = BoundaryType::MPI;
      bg.boundaries[LEFT_IDX] =
          new grid::BoundaryApplied<T, D - 1, BoundaryType::MPI>(
              dinfo.neighbor_left, rk);
      bg.boundaries[RIGHT_IDX] =
          new grid::BoundaryApplied<T, D - 1, BoundaryType::MPI>(
              dinfo.neighbor_right, rk);
    }
  }
#endif

 public:
#ifdef SYMPHAS_MPI_LOCAL_STORAGE
  SolverSystemFDMPI(symphas::init_data_type const& tdata,
                    symphas::interval_data_type const& vdata,
                    symphas::b_data_type const& bdata, size_t id = 0)
      : SolverSystemFDMPI(tdata, vdata, bdata, id,
                          compute_global_extended_dims(vdata)) {}
#else
  SolverSystemFDMPI(symphas::init_data_type const& tdata,
                    symphas::interval_data_type const& vdata,
                    symphas::b_data_type const& bdata, size_t id = 0)
      : BoundarySystem<T, D>(tdata, vdata, bdata, id),
        dframe{dims},
        dinfo{dims, BOUNDARY_DEPTH} {}
#endif

  SolverSystemFDMPI() : BoundarySystem<T, D>(), dframe{0}, dinfo{} {}

  inline void update(iter_type index, double time) {
    // BoundarySystem::update -> BoundaryGroup::update_boundaries iterates
    // the four sides and dispatches per-side updaters. Under MPI the TOP
    // side updater calls exchange_halos (Y+X in a single 2-step exchange).
    // No additional exchange is needed here.
    //
    // Phase 1: removed redundant exchange_halos(dframe.values, dinfo).
    // Phase 3 prep: removed redundant exchange_halos(values, dinfo) — that
    // exchange was already performed by the MPI boundary updater chain
    // above, so calling it again was a 2x halo-traffic wastage.
    BoundarySystem<T, D>::update(index, time);
  }
};

#endif

#ifdef USING_FFTW

//! The phase field system used by the spectral solver.
/*!
 * The spectral solver requires the Fourier transforms of all variables
 * when computing the solution, so in addition to storing the real-space
 * provisional system, the Fourier space transform is also stored and updated.
 *
 * The phase field does not use or define boundary conditions, although
 * they are still supplied to the constructor in order to use the
 * workflow.
 *
 * \tparam T The provisional system type, in real space.
 * \tparam D The provisional system dimension.
 */
template <typename T, size_t D>
struct SolverSystemSpectral : System<T, D> {
  using System<T, D>::System;
  using Grid<T, D>::values;

  //! Create the phase field data for the spectral solver implementation.
  SolverSystemSpectral(symphas::init_data_type const& tdata,
                       symphas::interval_data_type const& vdata,
                       symphas::b_data_type const& bdata, size_t id = 0)
      : System<T, D>(tdata, vdata, id) {}
  SolverSystemSpectral() : System<T, D>() {}
};

//! The phase field system used by the spectral solver.
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a real-valued order parameter. It defines the Fourier
 * transform of the order parameter with dimensions that do not include
 * the duplicated half.
 *
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template <size_t D>
struct SolverSystemSpectral<scalar_t, D> : System<scalar_t, D> {
  using System<scalar_t, D>::System;
  using Grid<scalar_t, D>::dims;
  using Grid<scalar_t, D>::values;

  len_type transformed_len;  //!< Length of the transformed array.
  complex_t* frame_t;        //!< The values of the transformed grid.
  complex_t* dframe;  //!< Accumulates the values of the solver computation.

  fftw_plan p;  //!< Fourier transform plan of the order parameter.

  //! Create the order parameter data used by the spectral solver.
  /*!
   * Create the order parameter data used by the spectral solver. The
   * boundaries are provided but not used in the implementation, as they are
   * imposed by the solver instead.
   *
   * \param tdata The initial conditions data of the system.
   * \param vdata The interval data of the system.
   * \param id The ID value of the system.
   */
  SolverSystemSpectral(symphas::init_data_type const& tdata,
                       symphas::interval_data_type const& vdata,
                       symphas::b_data_type const& bdata, size_t id);
  SolverSystemSpectral(SolverSystemSpectral<scalar_t, D> const& other);
  SolverSystemSpectral(SolverSystemSpectral<scalar_t, D>&& other) noexcept
      : SolverSystemSpectral() {
    swap(*this, other);
  }

  SolverSystemSpectral<scalar_t, D>& operator=(
      SolverSystemSpectral<scalar_t, D> other) {
    swap(*this, other);
    return *this;
  }

  //! Compute the Fourier transform and update the system.
  /*!
   * Compute the Fourier transform and update the system. On some iterations,
   * the complex transform is recomputed in order to eliminate numerical
   * error from accumulating. The values computed
   * by the solver are copied to the Fourier transformed data, and then
   * the Fourier transform is inverted to recover the real phase field values.
   *
   * \param index The index of the solution.
   */
  void update(iter_type index, double) {
    if (index % 100 == 0) {
      symphas::dft::arrange_fftw_stip<D>(
          values, reinterpret_cast<scalar_t*>(dframe), dims);
      symphas::dft::fftw_execute(p_to_t);
    }
    std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par,
#endif
        dframe, dframe + transformed_len, frame_t);

    symphas::dft::fftw_execute(p);
    symphas::dft::arrange_fftw_ipts<D>(reinterpret_cast<scalar_t*>(dframe),
                                       values, dims);
    grid::scale(System<scalar_t, D>::as_grid());
  }

  friend void swap(SolverSystemSpectral<scalar_t, D>& first,
                   SolverSystemSpectral<scalar_t, D>& second) {
    using std::swap;

    swap(static_cast<System<scalar_t, D>&>(first),
         static_cast<System<scalar_t, D>&>(second));
    swap(first.transformed_len, second.transformed_len);
    swap(first.frame_t, second.frame_t);
    swap(first.dframe, second.dframe);
    swap(first.p, second.p);
    swap(first.p_to_t, second.p_to_t);
  }

  SolverSystemSpectral()
      : System<scalar_t, D>(),
        transformed_len{0},
        frame_t{nullptr},
        dframe{nullptr},
        p{0},
        p_to_t{0} {}

  ~SolverSystemSpectral();

  fftw_plan p_to_t;  // Made public for SolverSP access.
};

//! The phase field system used by the spectral solver.
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a complex-valued order parameter.
 *
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template <size_t D>
struct SolverSystemSpectral<complex_t, D> : System<complex_t, D> {
  using System<complex_t, D>::System;
  using Grid<complex_t, D>::dims;
  using Grid<complex_t, D>::values;

  len_type transformed_len;  //!< Length of the transformed array.
  complex_t* frame_t;        //!< The values of the transformed grid.
  complex_t* dframe;  //!< Accumulates the values of the solver computation.
  fftw_plan p;        //!< Fourier transform plan of the order parameter.

  //! Create the order parameter data used by the spectral solver.
  /*!
   * Create the order parameter data used by the spectral solver. The
   * boundaries are provided but not used in the implementation, as they are
   * imposed by the solver instead.
   *
   * \param tdata The initial conditions data of the system.
   * \param vdata The interval data of the system.
   * \param id The ID value of the system.
   */
  SolverSystemSpectral(symphas::init_data_type const& tdata,
                       symphas::interval_data_type const& vdata,
                       symphas::b_data_type const& bdata, size_t id);
  SolverSystemSpectral(SolverSystemSpectral<complex_t, D> const& other);
  SolverSystemSpectral(SolverSystemSpectral<complex_t, D>&& other) noexcept
      : SolverSystemSpectral() {
    swap(*this, other);
  }
  SolverSystemSpectral<complex_t, D>& operator=(
      SolverSystemSpectral<complex_t, D> other) {
    swap(*this, other);
    return *this;
  }

  //! Compute the Fourier transform and update the system.
  /*!
   * Compute the Fourier transform and update the system. The values computed
   * by the solver are copied to the Fourier transformed data, and then
   * the Fourier transform is inverted to recover the real phase field values.
   *
   * \param index The index of the solution.
   */
  inline void update(iter_type, double) {
    std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par,
#endif
        dframe, dframe + transformed_len, frame_t);

    symphas::dft::fftw_execute(p);
    std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par,
#endif
        dframe, dframe + transformed_len, values);
    grid::scale(*this);
  }

  friend void swap(SolverSystemSpectral<complex_t, D>& first,
                   SolverSystemSpectral<complex_t, D>& second) {
    using std::swap;

    swap(static_cast<System<complex_t, D>&>(first),
         static_cast<System<complex_t, D>&>(second));
    swap(first.transformed_len, second.transformed_len);
    swap(first.frame_t, second.frame_t);
    swap(first.dframe, second.dframe);
    swap(first.p, second.p);
    swap(first.p_to_t, second.p_to_t);
  }

  SolverSystemSpectral()
      : System<complex_t, D>(),
        transformed_len{0},
        frame_t{nullptr},
        dframe{nullptr},
        p{0},
        p_to_t{0} {}

  ~SolverSystemSpectral();

 protected:
  fftw_plan p_to_t;
};

//! The phase field system used by the spectral solver.
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a vector-valued order parameter. It defines the Fourier
 * transform of the order parameter with dimensions that do not include
 * the duplicated half.
 *
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template <size_t D>
struct SolverSystemSpectral<vector_t<D>, D> : System<vector_t<D>, D> {
  using base_type = vector_t<D>;

  using System<base_type, D>::System;
  using Grid<base_type, D>::dims;
  using Grid<base_type, D>::axis;

  len_type transformed_len;          //!< Length of the transformed array.
  MultiBlock<D, complex_t> frame_t;  //!< The values of the transformed grid.
  MultiBlock<D, complex_t>
      dframe;  //!< Accumulates the values of the solver computation.

  fftw_plan p[D];  //!< Fourier transform plan of the order parameter.

  //! Create the order parameter data used by the spectral solver.
  /*!
   * Create the order parameter data used by the spectral solver. The
   * boundaries are provided but not used in the implementation, as they are
   * imposed by the solver instead.
   *
   * \param tdata The initial conditions data of the system.
   * \param vdata The interval data of the system.
   * \param id The ID value of the system.
   */
  SolverSystemSpectral(symphas::init_data_type const& tdata,
                       symphas::interval_data_type const& vdata,
                       symphas::b_data_type const& bdata, size_t id);
  SolverSystemSpectral(SolverSystemSpectral<base_type, D> const& other);
  SolverSystemSpectral(SolverSystemSpectral<base_type, D>&& other) noexcept
      : SolverSystemSpectral() {
    swap(*this, other);
  }
  SolverSystemSpectral<base_type, D>& operator=(
      SolverSystemSpectral<base_type, D> other) {
    swap(*this, other);
    return *this;
  }

  //! Compute the Fourier transform and update the system.
  /*!
   * Compute the Fourier transform and update the system. On some iterations,
   * the complex transform is recomputed in order to eliminate numerical
   * error from accumulating. The values computed
   * by the solver are copied to the Fourier transformed data, and then
   * the Fourier transform is inverted to recover the real phase field values.
   *
   * \param index The index of the solution.
   */
  void update(iter_type index, double) {
    if (index % 100 == 0) {
      for (iter_type i = 0; i < D; ++i) {
        Axis ax = symphas::index_to_axis(i);
        symphas::dft::arrange_fftw_stip<D>(
            axis(ax), reinterpret_cast<scalar_t*>(dframe(ax)), dims);
        symphas::dft::fftw_execute(p_to_t[i]);
      }
    }
    for (iter_type i = 0; i < D; ++i) {
      Axis ax = symphas::index_to_axis(i);

      std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
          std::execution::par,
#endif
          dframe(ax), dframe(ax) + transformed_len, frame_t(ax));

      symphas::dft::fftw_execute(p[i]);
      symphas::dft::arrange_fftw_ipts<D>(
          reinterpret_cast<scalar_t*>(dframe(ax)), axis(ax), dims);
    }
    grid::scale(System<base_type, D>::as_grid());
  }

  friend void swap(SolverSystemSpectral<base_type, D>& first,
                   SolverSystemSpectral<base_type, D>& second) {
    using std::swap;

    swap(static_cast<System<base_type, D>&>(first),
         static_cast<System<base_type, D>&>(second));
    swap(first.transformed_len, second.transformed_len);
    swap(first.frame_t, second.frame_t);
    swap(first.dframe, second.dframe);
    swap(first.p, second.p);
    swap(first.p_to_t, second.p_to_t);
  }

  SolverSystemSpectral()
      : System<base_type, D>(),
        transformed_len{0},
        frame_t{0},
        dframe{0},
        p{0},
        p_to_t{0} {}

  ~SolverSystemSpectral();

 protected:
  fftw_plan p_to_t[D];
};

using symphas::dft::new_fftw_plan;

template <>
inline SolverSystemSpectral<scalar_t, 1>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<scalar_t, 1>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 1>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<1, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<1, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<1>(Grid<scalar_t, 1>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<scalar_t, 2>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<scalar_t, 2>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 2>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<2, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<2, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<2>(Grid<scalar_t, 2>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<scalar_t, 3>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<scalar_t, 3>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 3>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<3, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<3, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<3>(Grid<scalar_t, 3>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 1>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<complex_t, 1>(tdata, vdata, id),
      transformed_len{symphas::dft::length<complex_t, 1>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 1>::values,
            Grid<complex_t, 1>::values + Grid<complex_t, 1>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 2>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<complex_t, 2>(tdata, vdata, id),
      transformed_len{symphas::dft::length<complex_t, 2>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 2>::values,
            Grid<complex_t, 2>::values + Grid<complex_t, 2>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 3>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<complex_t, 3>(tdata, vdata, id),
      transformed_len{symphas::dft::length<complex_t, 3>(dims)},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 3>::values,
            Grid<complex_t, 3>::values + Grid<complex_t, 3>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<any_vector_t<scalar_t, 1>, 1>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 1>(dims)},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<1, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims)},
      p_to_t{new_fftw_plan<1, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims)} {
  symphas::dft::arrange_fftw_stip<1>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<any_vector_t<scalar_t, 2>, 2>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 2>(dims)},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims),
        new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::Y),
                                                dframe(Axis::Y), dims)},
      p_to_t{new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims),
             new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::Y),
                                                     dframe(Axis::Y), dims)} {
  symphas::dft::arrange_fftw_stip<2>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::arrange_fftw_stip<2>(
      axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
  symphas::dft::fftw_execute(p_to_t[1]);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3>::SolverSystemSpectral(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata, symphas::b_data_type const&,
    size_t id)
    : System<any_vector_t<scalar_t, 3>, 3>(tdata, vdata, id),
      transformed_len{symphas::dft::length<scalar_t, 3>(dims)},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims),
        new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Y),
                                                dframe(Axis::Y), dims),
        new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Z),
                                                dframe(Axis::Z), dims)},
      p_to_t{new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims),
             new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Y),
                                                     dframe(Axis::Y), dims),
             new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Z),
                                                     dframe(Axis::Z), dims)} {
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::Z), reinterpret_cast<scalar_t*>(dframe(Axis::Z)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
  symphas::dft::fftw_execute(p_to_t[1]);
  symphas::dft::fftw_execute(p_to_t[2]);
}

template <>
inline SolverSystemSpectral<scalar_t, 1>::SolverSystemSpectral(
    SolverSystemSpectral<scalar_t, 1> const& other)
    : System<scalar_t, 1>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<1, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<1, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<1>(Grid<scalar_t, 1>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<scalar_t, 2>::SolverSystemSpectral(
    SolverSystemSpectral<scalar_t, 2> const& other)
    : System<scalar_t, 2>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<2, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<2, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<2>(Grid<scalar_t, 2>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<scalar_t, 3>::SolverSystemSpectral(
    SolverSystemSpectral<scalar_t, 3> const& other)
    : System<scalar_t, 3>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<3, complex_t, scalar_t>{}(dframe, dframe, dims)},
      p_to_t{new_fftw_plan<3, scalar_t, complex_t>{}(dframe, dframe, dims)} {
  symphas::dft::arrange_fftw_stip<3>(Grid<scalar_t, 3>::values,
                                     reinterpret_cast<scalar_t*>(dframe), dims);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 1>::SolverSystemSpectral(
    SolverSystemSpectral<complex_t, 1> const& other)
    : System<complex_t, 1>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 1>::values,
            Grid<complex_t, 1>::values + Grid<complex_t, 1>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 2>::SolverSystemSpectral(
    SolverSystemSpectral<complex_t, 2> const& other)
    : System<complex_t, 2>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 2>::values,
            Grid<complex_t, 2>::values + Grid<complex_t, 2>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<complex_t, 3>::SolverSystemSpectral(
    SolverSystemSpectral<complex_t, 3> const& other)
    : System<complex_t, 3>(other),
      transformed_len{other.transformed_len},
      frame_t{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      dframe{reinterpret_cast<complex_t*>(
          symphas::dft::fftw_alloc_complex(transformed_len))},
      p{new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false,
                                                 true)},
      p_to_t{new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims,
                                                      false, false)} {
  std::copy(Grid<complex_t, 3>::values,
            Grid<complex_t, 3>::values + Grid<complex_t, 3>::len, dframe);
  symphas::dft::fftw_execute(p_to_t);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1>::SolverSystemSpectral(
    SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1> const& other)
    : System<any_vector_t<scalar_t, 1>, 1>(other),
      transformed_len{other.transformed_len},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<1, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims)},
      p_to_t{new_fftw_plan<1, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims)} {
  symphas::dft::arrange_fftw_stip<1>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2>::SolverSystemSpectral(
    SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2> const& other)
    : System<any_vector_t<scalar_t, 2>, 2>(other),
      transformed_len{other.transformed_len},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims),
        new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::Y),
                                                dframe(Axis::Y), dims)},
      p_to_t{new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims),
             new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::Y),
                                                     dframe(Axis::Y), dims)} {
  symphas::dft::arrange_fftw_stip<2>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::arrange_fftw_stip<2>(
      axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
  symphas::dft::fftw_execute(p_to_t[1]);
}

template <>
inline SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3>::SolverSystemSpectral(
    SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3> const& other)
    : System<any_vector_t<scalar_t, 3>, 3>(other),
      transformed_len{other.transformed_len},
      frame_t{transformed_len},
      dframe{transformed_len},
      p{new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::X),
                                                dframe(Axis::X), dims),
        new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Y),
                                                dframe(Axis::Y), dims),
        new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Z),
                                                dframe(Axis::Z), dims)},
      p_to_t{new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::X),
                                                     dframe(Axis::X), dims),
             new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Y),
                                                     dframe(Axis::Y), dims),
             new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Z),
                                                     dframe(Axis::Z), dims)} {
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
  symphas::dft::arrange_fftw_stip<3>(
      axis(Axis::Z), reinterpret_cast<scalar_t*>(dframe(Axis::Z)), dims);
  symphas::dft::fftw_execute(p_to_t[0]);
  symphas::dft::fftw_execute(p_to_t[1]);
  symphas::dft::fftw_execute(p_to_t[2]);
}

template <size_t D>
inline SolverSystemSpectral<scalar_t, D>::~SolverSystemSpectral() {
  symphas::dft::fftw_destroy_plan(p);
  symphas::dft::fftw_destroy_plan(p_to_t);
  symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(dframe));
  symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(frame_t));
}

#if defined(USING_FFTW) && defined(USING_MPI) && defined(USING_FFTW_MPI)

//! Distributed spectral system for MPI.
/*!
 * Uses fftw_mpi distributed r2c/c2r transforms. Data is distributed
 * along the y-axis (leftmost dimension in FFTW row-major convention).
 * Each rank holds local_n0 rows of k-space data.
 *
 * The real-space grid (values[]) holds the FULL global grid, but only the
 * local slab [local_0_start .. local_0_start+local_n0) is kept up-to-date
 * during the solve loop. No MPI communication occurs per timestep.
 * Call sync_full_grid() before I/O to reconstruct the full grid via Allgatherv.
 */
template <size_t D>
struct SolverSystemSpectralMPI : System<scalar_t, D> {
  using System<scalar_t, D>::System;
  using Grid<scalar_t, D>::dims;
  using Grid<scalar_t, D>::values;

  static_assert(D == 2 || D == 3,
                "distributed spectral MPI system supports D=2 and D=3");

  len_type transformed_len;   //!< Local k-space length.
  complex_t* frame_t;         //!< Local k-space snapshot (previous step).
  complex_t* dframe;          //!< Local k-space accumulator (solver writes here).
  scalar_t* real_work;        //!< Padded real workspace for MPI FFT.

  fftw_plan p;       //!< Inverse FFT plan (c2r).
  fftw_plan p_to_t;  //!< Forward FFT plan (r2c).
  bool owns_plans;   //!< Whether this instance owns (and should destroy) the plans.

  ptrdiff_t local_n0;       //!< Number of slab rows/planes on this rank.
  ptrdiff_t local_0_start;  //!< Global slab offset of this rank.
  ptrdiff_t alloc_local;    //!< Allocation size returned by fftw_mpi.

  SolverSystemSpectralMPI(symphas::init_data_type const& tdata,
                          symphas::interval_data_type const& vdata,
                          symphas::b_data_type const& bdata, size_t id);

  SolverSystemSpectralMPI()
      : System<scalar_t, D>(),
        transformed_len{0}, frame_t{nullptr}, dframe{nullptr},
        real_work{nullptr}, p{0}, p_to_t{0}, owns_plans{false},
        local_n0{0}, local_0_start{0}, alloc_local{0} {}

  SolverSystemSpectralMPI(SolverSystemSpectralMPI&& other) noexcept
      : SolverSystemSpectralMPI() {
    swap(*this, other);
  }

  SolverSystemSpectralMPI(SolverSystemSpectralMPI const& other);

  SolverSystemSpectralMPI& operator=(SolverSystemSpectralMPI other) {
    swap(*this, other);
    return *this;
  }

  //! Number of real cells per slab index (a y-row in 2D, a z-plane in 3D).
  len_type plane_size() const {
    if constexpr (D == 2) {
      return dims[0];
    } else {
      return dims[0] * dims[1];
    }
  }

  //! Local slab offset and length in the values[] array.
  len_type local_real_start() const {
    return static_cast<len_type>(local_0_start) * plane_size();
  }
  len_type local_real_len() const {
    return static_cast<len_type>(local_n0) * plane_size();
  }

  //! Copy local real-space slab into padded MPI workspace.
  void scatter_real_to_work() {
    len_type Nx = dims[0];
    len_type row_pad = 2 * (Nx / 2 + 1);
    if constexpr (D == 2) {
      for (ptrdiff_t j = 0; j < local_n0; ++j) {
        ptrdiff_t global_j = local_0_start + j;
        for (len_type i = 0; i < Nx; ++i) {
          real_work[j * row_pad + i] = values[global_j * Nx + i];
        }
      }
    } else {
      len_type Ny = dims[1];
      for (ptrdiff_t k = 0; k < local_n0; ++k) {
        ptrdiff_t global_k = local_0_start + k;
        for (len_type j = 0; j < Ny; ++j) {
          for (len_type i = 0; i < Nx; ++i) {
            real_work[(k * Ny + j) * row_pad + i] =
                values[(global_k * Ny + j) * Nx + i];
          }
        }
      }
    }
  }

  //! Unpad inverse FFT result into the local slab of values[].
  void unpad_local_slab() {
    len_type Nx = dims[0];
    len_type row_pad = 2 * (Nx / 2 + 1);
    if constexpr (D == 2) {
      for (ptrdiff_t j = 0; j < local_n0; ++j) {
        ptrdiff_t global_j = local_0_start + j;
        for (len_type i = 0; i < Nx; ++i) {
          values[global_j * Nx + i] = real_work[j * row_pad + i];
        }
      }
    } else {
      len_type Ny = dims[1];
      for (ptrdiff_t k = 0; k < local_n0; ++k) {
        ptrdiff_t global_k = local_0_start + k;
        for (len_type j = 0; j < Ny; ++j) {
          for (len_type i = 0; i < Nx; ++i) {
            values[(global_k * Ny + j) * Nx + i] =
                real_work[(k * Ny + j) * row_pad + i];
          }
        }
      }
    }
  }

  //! Allgatherv to reconstruct full grid on all ranks. Call before I/O only.
  void sync_full_grid() {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    std::vector<int> recvcounts(nprocs), displs(nprocs);
    int local_count = static_cast<int>(local_n0 * plane_size());
    MPI_Allgather(&local_count, 1, MPI_INT,
                  recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    displs[0] = 0;
    for (int r = 1; r < nprocs; ++r)
      displs[r] = displs[r-1] + recvcounts[r-1];
    MPI_Allgatherv(MPI_IN_PLACE, local_count, MPI_DOUBLE,
                   values, recvcounts.data(), displs.data(),
                   MPI_DOUBLE, MPI_COMM_WORLD);
  }

  void update(iter_type index, double) {
    SYMPHAS_MPI_PROFILE_SCOPE("sp_update");
    // Periodically refresh k-space from local real-space to prevent drift.
    if (index % 100 == 0) {
      SYMPHAS_MPI_PROFILE_SCOPE("sp_update_refresh");
      scatter_real_to_work();
      symphas::dft::fftw_execute(p_to_t);
    }
    // dframe has the updated k-space solution; copy to frame_t.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp_update_copy_dframe");
      std::copy(dframe, dframe + transformed_len, frame_t);
    }
    // Inverse FFT (c2r): dframe → real_work.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp_update_inv_fft");
      symphas::dft::fftw_execute(p);
    }
    // Unpad only local slab — no MPI communication.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp_update_unpad");
      unpad_local_slab();
    }
    // Scale only the local slab by 1/(total grid points).
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp_update_scale");
      len_type start = local_real_start();
      len_type len = local_real_len();
      // The inverse FFT is unnormalized; divide by the GLOBAL total real-grid
      // point count (product of every global dimension). The local `len` above
      // is only this rank's slab and must not be used here. Looping over all D
      // dimensions is correct for any dimensionality (including D == 1).
      double total = 1.0;
      for (size_t d = 0; d < D; ++d) total *= static_cast<double>(dims[d]);
      double scale = 1.0 / total;
      for (len_type i = start; i < start + len; ++i)
        values[i] *= scale;
    }
  }

  friend void swap(SolverSystemSpectralMPI& a, SolverSystemSpectralMPI& b) {
    using std::swap;
    swap(static_cast<System<scalar_t, D>&>(a),
         static_cast<System<scalar_t, D>&>(b));
    swap(a.transformed_len, b.transformed_len);
    swap(a.frame_t, b.frame_t);
    swap(a.dframe, b.dframe);
    swap(a.real_work, b.real_work);
    swap(a.p, b.p);
    swap(a.p_to_t, b.p_to_t);
    swap(a.owns_plans, b.owns_plans);
    swap(a.local_n0, b.local_n0);
    swap(a.local_0_start, b.local_0_start);
    swap(a.alloc_local, b.alloc_local);
  }

  ~SolverSystemSpectralMPI() {
    if (owns_plans) {
      if (p) symphas::dft::fftw_destroy_plan(p);
      if (p_to_t) symphas::dft::fftw_destroy_plan(p_to_t);
    }
    if (dframe) symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(dframe));
    if (frame_t) symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(frame_t));
    if (real_work) symphas::dft::fftw_free(real_work);
  }
};

//! Detects whether a solver-system type is a distributed spectral MPI system.
/*!
 * Matches any \c SolverSystemSpectralMPI<D> regardless of dimension so that
 * dispatch sites (solver equation branch, I/O gather) can recognize the
 * distributed spectral path uniformly across 2D and 3D.
 */
template <typename S>
struct is_spectral_mpi_system : std::false_type {};
template <size_t D>
struct is_spectral_mpi_system<SolverSystemSpectralMPI<D>> : std::true_type {};
template <typename S>
inline constexpr bool is_spectral_mpi_system_v =
    is_spectral_mpi_system<std::decay_t<S>>::value;

template <size_t D>
inline SolverSystemSpectralMPI<D>::SolverSystemSpectralMPI(
    symphas::init_data_type const& tdata,
    symphas::interval_data_type const& vdata,
    symphas::b_data_type const&, size_t id)
    : System<scalar_t, D>(tdata, vdata, id),
      transformed_len{0}, frame_t{nullptr}, dframe{nullptr},
      real_work{nullptr}, p{0}, p_to_t{0}, owns_plans{true},
      local_n0{0}, local_0_start{0}, alloc_local{0}
{
  // FFTW MPI slab-decomposes along the first (slowest) logical axis. The grid
  // is stored x-contiguous, so the FFTW logical dimensions are (Ny,Nx) in 2D
  // and (Nz,Ny,Nx) in 3D; the contiguous half-spectrum axis is Nx.
  ptrdiff_t Nx = dims[0];
  if constexpr (D == 2) {
    ptrdiff_t Ny = dims[1];
    alloc_local = symphas::dft::fftw_mpi_local_size_2d(
        Ny, Nx, MPI_COMM_WORLD, &local_n0, &local_0_start);
    transformed_len = static_cast<len_type>(local_n0 * (Nx / 2 + 1));
  } else {
    ptrdiff_t Ny = dims[1], Nz = dims[2];
    alloc_local = symphas::dft::fftw_mpi_local_size_3d(
        Nz, Ny, Nx, MPI_COMM_WORLD, &local_n0, &local_0_start);
    transformed_len = static_cast<len_type>(local_n0 * Ny * (Nx / 2 + 1));
  }

  // Allocate k-space arrays. Must use alloc_local (not transformed_len)
  // because FFTW MPI may need extra space for the internal transpose.
  frame_t = reinterpret_cast<complex_t*>(
      symphas::dft::fftw_alloc_complex(alloc_local));
  dframe = reinterpret_cast<complex_t*>(
      symphas::dft::fftw_alloc_complex(alloc_local));

  // Allocate padded real workspace for MPI r2c/c2r.
  // alloc_local is in units of complex, so 2*alloc_local doubles.
  real_work = symphas::dft::fftw_alloc_real(2 * alloc_local);

  // Create MPI plans (in-place on real_work / cast to dframe).
  if constexpr (D == 2) {
    ptrdiff_t Ny = dims[1];
    p_to_t = symphas::dft::fftw_mpi_plan_r2c_2d(
        Ny, Nx, real_work, reinterpret_cast<fftw_complex*>(dframe),
        MPI_COMM_WORLD);
    p = symphas::dft::fftw_mpi_plan_c2r_2d(
        Ny, Nx, reinterpret_cast<fftw_complex*>(dframe), real_work,
        MPI_COMM_WORLD);
  } else {
    ptrdiff_t Ny = dims[1], Nz = dims[2];
    p_to_t = symphas::dft::fftw_mpi_plan_r2c_3d(
        Nz, Ny, Nx, real_work, reinterpret_cast<fftw_complex*>(dframe),
        MPI_COMM_WORLD);
    p = symphas::dft::fftw_mpi_plan_c2r_3d(
        Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(dframe), real_work,
        MPI_COMM_WORLD);
  }

  // Initial forward FFT.
  scatter_real_to_work();
  symphas::dft::fftw_execute(p_to_t);
}

template <size_t D>
inline SolverSystemSpectralMPI<D>::SolverSystemSpectralMPI(
    SolverSystemSpectralMPI const& other)
    : System<scalar_t, D>(other),
      transformed_len{other.transformed_len}, frame_t{nullptr}, dframe{nullptr},
      real_work{nullptr}, p{0}, p_to_t{0}, owns_plans{true},
      local_n0{other.local_n0}, local_0_start{other.local_0_start},
      alloc_local{other.alloc_local}
{
  // Allocate own buffers. Must use alloc_local for k-space arrays
  // (FFTW MPI may need extra space beyond transformed_len for transpose).
  frame_t = reinterpret_cast<complex_t*>(
      symphas::dft::fftw_alloc_complex(alloc_local));
  dframe = reinterpret_cast<complex_t*>(
      symphas::dft::fftw_alloc_complex(alloc_local));
  real_work = symphas::dft::fftw_alloc_real(2 * alloc_local);

  // Create own MPI plans bound to our buffers.
  // fftw_mpi_plan_* are collective — all ranks must call simultaneously.
  ptrdiff_t Nx = dims[0];
  if constexpr (D == 2) {
    ptrdiff_t Ny = dims[1];
    p_to_t = symphas::dft::fftw_mpi_plan_r2c_2d(
        Ny, Nx, real_work, reinterpret_cast<fftw_complex*>(dframe),
        MPI_COMM_WORLD);
    p = symphas::dft::fftw_mpi_plan_c2r_2d(
        Ny, Nx, reinterpret_cast<fftw_complex*>(dframe), real_work,
        MPI_COMM_WORLD);
  } else {
    ptrdiff_t Ny = dims[1], Nz = dims[2];
    p_to_t = symphas::dft::fftw_mpi_plan_r2c_3d(
        Nz, Ny, Nx, real_work, reinterpret_cast<fftw_complex*>(dframe),
        MPI_COMM_WORLD);
    p = symphas::dft::fftw_mpi_plan_c2r_3d(
        Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(dframe), real_work,
        MPI_COMM_WORLD);
  }

  // Copy k-space state from the original.
  std::copy(other.frame_t, other.frame_t + transformed_len, frame_t);
  std::copy(other.dframe, other.dframe + transformed_len, dframe);
}

#endif  // USING_FFTW && USING_MPI && USING_FFTW_MPI

template <size_t D>
inline SolverSystemSpectral<complex_t, D>::~SolverSystemSpectral() {
  symphas::dft::fftw_destroy_plan(p);
  symphas::dft::fftw_destroy_plan(p_to_t);
  symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(dframe));
  symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(frame_t));
}

template <size_t D>
inline SolverSystemSpectral<vector_t<D>, D>::~SolverSystemSpectral() {
  for (iter_type i = 0; i < D; ++i) {
    symphas::dft::fftw_destroy_plan(p[i]);
    symphas::dft::fftw_destroy_plan(p_to_t[i]);
  }
}

#endif

#ifdef USING_MPI
DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFDwSDMPI<T, D>),
                           (RegionalGridMPI<T, D>))
#endif

DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFDwSD<T, D>),
                           (RegionalGrid<T, D>))
DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFD<T, D>),
                           (BoundaryGrid<T, D>))

#ifdef USING_FFTW
DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemSpectral<T, D>),
                           (Grid<T, D>))
#endif

#ifdef USING_CUDA
template <typename T, size_t D>
struct SolverSystemFDwSDCUDA;
template <typename T, size_t D>
struct SolverSystemFDCUDA;

// Forward declaration of the GPU spectral system (defined in solversystem.cuh,
// which carries device syntax and is only included by nvcc translation units).
// Declaring it here — in the host-parseable header — lets the SP solver's
// system-type association and the is_spectral_cuda_system trait below be named
// by ordinary (host-compiled) translation units without pulling in cuFFT.
template <size_t D>
struct SolverSystemSpectralCUDA;

//! Detects a GPU spectral solver system (any dimension) for solver dispatch.
//! Visible to both host and device translation units (no CUDA syntax here).
template <typename S>
struct is_spectral_cuda_system : std::false_type {};
template <size_t D>
struct is_spectral_cuda_system<SolverSystemSpectralCUDA<D>> : std::true_type {};
template <typename S>
inline constexpr bool is_spectral_cuda_system_v =
    is_spectral_cuda_system<std::decay_t<S>>::value;
#endif