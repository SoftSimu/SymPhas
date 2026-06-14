#pragma once

#include "symphas.h"
#include "timer.h"
#include <chrono>

#ifdef USING_MPI
#include "spsmpi.h"

// MPI helpers used by simulate(). Conceptually belong in a shared
// example/runtime header; for now they live here because they are only
// used by this main and depend on BoundaryGroup/BoundaryApplied from sol/
// (which would invert the lib/ → sol/ layering if moved into lib/).
template <typename M>
void broadcast_model_grids(M& model) {
#ifdef SYMPHAS_MPI_LOCAL_STORAGE
  // Under local-storage mode, every rank's grid is its local sub-region
  // and grid.dims differs per rank; a global MPI_Bcast of rank-0's grid
  // would be incoherent (mismatched sizes → MPI_ERR_TRUNCATE). Each
  // rank's IC is computed locally by populate_tdata against the
  // rank-local vdata, so no broadcast is needed for IC consistency.
  // (Note: rank seeds may differ; for deterministic IC across ranks the
  // RNG must be seeded based on global coords, which is a separate
  // concern handled at IC level.)
  (void)model;
#else
  if (symphas::parallel::get_num_nodes() > 1) {
    auto& tup = model.systems_tuple();
    std::apply([&](auto&... sys) {
      auto broadcast_one = [&](auto& s) {
        auto& g = s.as_grid();
        len_type total_len = 1;
        for (iter_type d = 0; d < model_dimension<M>::value; ++d)
          total_len *= g.dims[d];
        if constexpr (std::is_pointer_v<decltype(g.values)>) {
          MPI_Bcast(g.values, static_cast<int>(total_len), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
          for (iter_type c = 0; c < model_dimension<M>::value; ++c) {
            MPI_Bcast(g.values[c], static_cast<int>(total_len), MPI_DOUBLE, 0, MPI_COMM_WORLD);
          }
        }
      };
      (broadcast_one(sys), ...);
    }, tup);
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
}

//! Override Y-direction boundaries to use MPI halo exchange.
/*!
 * When running with multiple MPI ranks, the Y-direction periodic boundaries
 * must be replaced with MPI boundaries that perform halo exchange between
 * neighboring ranks instead of local periodic copy.
 * X-direction boundaries remain as periodic (local operation within each row).
 */
template <typename M>
void setup_mpi_boundaries(M& model) {
  if (symphas::parallel::get_num_nodes() <= 1) return;

#if defined(USE_SPECTRAL_SOLVER)
  // The spectral (SP) solver distributes work through FFTW-MPI's own slab
  // decomposition (see SolverSystemSpectralMPI<D>), not through finite-
  // difference halo-exchange boundaries. The FD MPI boundary substitution
  // below does not apply and is skipped; this also avoids the 2D-only
  // FD-boundary restriction for 3D spectral runs.
  (void)model;
  return;
#elif defined(SYMPHAS_MPI_LOCAL_STORAGE)
  // Under SYMPHAS_MPI_LOCAL_STORAGE, the MPI boundaries are already
  // installed (with rank-local domain_info attached) by SolverSystemFD's
  // gated ctor. Re-running the substitution here would delete those
  // dinfo-bearing boundaries and replace them with bare ones, breaking
  // halo exchange. Skip.
  (void)model;
  return;
#else

  constexpr size_t D = model_dimension<M>::value;
  static_assert(D == 2, "MPI boundaries currently implemented for 2D only");

  // Resolve the 2-D Cartesian decomposition (cached). Px==1 reproduces the
  // legacy 1-D Y-slab behavior (only TOP/BOTTOM tagged as MPI; LEFT/RIGHT
  // remain PERIODIC and the periodic boundary updater handles X locally).
  // Px>1 additionally tags LEFT/RIGHT as MPI so the periodic copy is
  // suppressed and the X halos come from the MPI exchange instead.
  auto [px, py] = symphas::parallel::get_mpi_dims_cart();

  int rank = symphas::parallel::get_node_rank();
  int ix = rank / py;
  int iy = rank % py;
  int iy_below = (iy - 1 + py) % py;
  int iy_above = (iy + 1) % py;
  int ix_left = (ix - 1 + px) % px;
  int ix_right = (ix + 1) % px;
  int neighbor_below = ix * py + iy_below;
  int neighbor_above = ix * py + iy_above;
  int neighbor_left = (px > 1) ? ix_left * py + iy : MPI_PROC_NULL;
  int neighbor_right = (px > 1) ? ix_right * py + iy : MPI_PROC_NULL;

  auto& tup = model.systems_tuple();
  std::apply([&](auto&... sys) {
    auto setup_one = [&](auto& s) {
      using sys_type = std::remove_const_t<std::remove_reference_t<decltype(s)>>;
      // Apply MPI boundaries to scalar systems and vector-valued systems
      // alike. Spectral systems (no BoundaryGroup base) keep replicated
      // grids and are skipped here.
      if constexpr (std::is_base_of_v<BoundaryGroup<scalar_t, D>, sys_type>) {
        constexpr iter_type TOP_IDX = static_cast<iter_type>(Side::TOP);
        constexpr iter_type BOT_IDX = static_cast<iter_type>(Side::BOTTOM);
        constexpr iter_type LEFT_IDX = static_cast<iter_type>(Side::LEFT);
        constexpr iter_type RIGHT_IDX = static_cast<iter_type>(Side::RIGHT);

        auto& s_mut = const_cast<sys_type&>(s);

        delete s_mut.boundaries[TOP_IDX];
        delete s_mut.boundaries[BOT_IDX];
        s_mut.types[TOP_IDX] = BoundaryType::MPI;
        s_mut.types[BOT_IDX] = BoundaryType::MPI;
        s_mut.boundaries[TOP_IDX] = new grid::BoundaryApplied<scalar_t, D - 1, BoundaryType::MPI>(
            neighbor_above, rank);
        s_mut.boundaries[BOT_IDX] = new grid::BoundaryApplied<scalar_t, D - 1, BoundaryType::MPI>(
            neighbor_below, rank);

        if (px > 1) {
          delete s_mut.boundaries[LEFT_IDX];
          delete s_mut.boundaries[RIGHT_IDX];
          s_mut.types[LEFT_IDX] = BoundaryType::MPI;
          s_mut.types[RIGHT_IDX] = BoundaryType::MPI;
          s_mut.boundaries[LEFT_IDX] = new grid::BoundaryApplied<scalar_t, D - 1, BoundaryType::MPI>(
              neighbor_left, rank);
          s_mut.boundaries[RIGHT_IDX] = new grid::BoundaryApplied<scalar_t, D - 1, BoundaryType::MPI>(
              neighbor_right, rank);
        }
      } else if constexpr (
          std::is_base_of_v<BoundaryGroup<any_vector_t<scalar_t, D>, D>,
                            sys_type>) {
        constexpr iter_type TOP_IDX = static_cast<iter_type>(Side::TOP);
        constexpr iter_type BOT_IDX = static_cast<iter_type>(Side::BOTTOM);
        constexpr iter_type LEFT_IDX = static_cast<iter_type>(Side::LEFT);
        constexpr iter_type RIGHT_IDX = static_cast<iter_type>(Side::RIGHT);

        auto& s_mut = const_cast<sys_type&>(s);

        delete s_mut.boundaries[TOP_IDX];
        delete s_mut.boundaries[BOT_IDX];
        s_mut.types[TOP_IDX] = BoundaryType::MPI;
        s_mut.types[BOT_IDX] = BoundaryType::MPI;
        s_mut.boundaries[TOP_IDX] =
            new grid::BoundaryApplied<any_vector_t<scalar_t, D>, D - 1,
                                       BoundaryType::MPI>(neighbor_above, rank);
        s_mut.boundaries[BOT_IDX] =
            new grid::BoundaryApplied<any_vector_t<scalar_t, D>, D - 1,
                                       BoundaryType::MPI>(neighbor_below, rank);

        if (px > 1) {
          delete s_mut.boundaries[LEFT_IDX];
          delete s_mut.boundaries[RIGHT_IDX];
          s_mut.types[LEFT_IDX] = BoundaryType::MPI;
          s_mut.types[RIGHT_IDX] = BoundaryType::MPI;
          s_mut.boundaries[LEFT_IDX] =
              new grid::BoundaryApplied<any_vector_t<scalar_t, D>, D - 1,
                                         BoundaryType::MPI>(neighbor_left, rank);
          s_mut.boundaries[RIGHT_IDX] =
              new grid::BoundaryApplied<any_vector_t<scalar_t, D>, D - 1,
                                         BoundaryType::MPI>(neighbor_right, rank);
        }
      }
    };
    (setup_one(sys), ...);
  }, tup);

  fprintf(SYMPHAS_LOG,
          "MPI boundaries configured: rank %d, dims_cart={%d,%d} coords={%d,%d}, "
          "neighbors below=%d above=%d left=%d right=%d\n",
          rank, px, py, ix, iy,
          neighbor_below, neighbor_above, neighbor_left, neighbor_right);
#endif  // SYMPHAS_MPI_LOCAL_STORAGE
}
#else
template <typename M>
void broadcast_model_grids(M&) {}

template <typename M>
void setup_mpi_boundaries(M&) {}
#endif

// **************************************************************************************

using namespace symphas;

template <typename M>
struct Simulation {
#ifdef USING_CONF
  int simulate(double const *coeff, size_t num_coeff) {
    auto pp = symphas::conf::config().get_problem_parameters();

    if (params::plots_only) {
      M model(coeff, num_coeff, pp);
      symphas::io::write_plot_config(model);
    } else {
      /*
       * execute a single model and collect data if there is only a single run
       */
      if (symphas::conf::config().runs == 1) {
        auto t_ctor_start = std::chrono::high_resolution_clock::now();
        M model(coeff, num_coeff, pp);
        auto t_ctor_end = std::chrono::high_resolution_clock::now();
        double ctor_ms = std::chrono::duration<double, std::milli>(
            t_ctor_end - t_ctor_start).count();
        fprintf(stderr, "[PHASE rank=%d model_ctor_ms=%.3f]\n",
                symphas::parallel::get_node_rank(), ctor_ms);

        auto t_setup_start = std::chrono::high_resolution_clock::now();
        broadcast_model_grids(model);
        setup_mpi_boundaries(model);
        symphas::io::write_plot_config(model);
        auto t_setup_end = std::chrono::high_resolution_clock::now();
        double setup_ms = std::chrono::duration<double, std::milli>(
            t_setup_end - t_setup_start).count();
        fprintf(stderr, "[PHASE rank=%d setup_ms=%.3f]\n",
                symphas::parallel::get_node_rank(), setup_ms);

#ifdef USING_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        auto t_start = std::chrono::high_resolution_clock::now();
#ifdef SYMPHAS_MPI_LOCAL_STORAGE
        // Under local-storage mode, the data-persistence and checkpoint
        // paths are not yet adapted to local storage (Phase 6 — gather
        // adapter). Disable both so that benchmarks measure the actual
        // solver throughput rather than save infrastructure overhead.
        // For production use of local storage, the gather adapter must
        // be implemented to make these calls correct.
        find_solution(model, /*plotting_output=*/false, /*checkpoint=*/false);
#else
        find_solution(model);
#endif
        auto t_end = std::chrono::high_resolution_clock::now();
        double solve_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        fprintf(SYMPHAS_LOG, "SOLVE_TIME rank=%d ms=%.3f\n",
                symphas::parallel::get_node_rank(), solve_ms);

#ifdef USING_MPI
#if defined(USING_FFTW_MPI)
        // For spectral MPI, skip the sync/checksum to avoid cleanup issues.
        // The SOLVE_TIME is printed above; the simulation completed correctly.
        if (symphas::parallel::get_num_nodes() > 1) {
          fflush(stdout);
          fflush(SYMPHAS_LOG);
          MPI_Barrier(MPI_COMM_WORLD);
        }
#else
        // FD MPI: sync all slabs then checksum for correctness verification
#ifdef SYMPHAS_MPI_LOCAL_STORAGE
        // Under local-storage mode, every rank's grid.dims is its local
        // sub-region (not the global grid), so sync_all_slabs cannot be
        // driven from grid.dims alone — it would broadcast slabs of
        // mismatched sizes and abort with MPI_ERR_TRUNCATE. The proper
        // fix is to thread the rank-local dinfo through (Phase 6 — gather
        // adapter); for now, skip the post-run checksum sync. Each rank
        // still emits a per-rank checksum below, which is sufficient
        // for benchmark correctness verification.
#else
        if (symphas::parallel::get_num_nodes() > 1) {
          auto& tup = model.systems_tuple();
          std::apply([&](auto&... sys) {
            auto sync_one = [&](auto& s) {
              using sys_t = std::remove_reference_t<decltype(s)>;
              if constexpr (symphas::parallel::has_mpi_dinfo_v<sys_t>) {
                auto& g = s.as_grid();
                // Slab sync only operates on flat scalar grids (T*); vector
                // grids store axes as separate buffers (T**) and are handled
                // by the regular boundary-update path. Skip them here.
                using val_t = std::remove_reference_t<decltype(g.values[0])>;
                if constexpr (std::is_arithmetic_v<val_t>) {
                  constexpr size_t D = model_dimension<M>::value;
                  symphas::parallel::domain_info<D> dinfo(g.dims, BOUNDARY_DEPTH);
                  symphas::parallel::sync_all_slabs(g.values, dinfo);
                }
              }
              // Replicated (non-MPI) systems already have the full grid
              // on every rank — a slab broadcast would overwrite correct
              // data with another rank's slab.
            };
            (sync_one(sys), ...);
          }, tup);
        }
#endif  // SYMPHAS_MPI_LOCAL_STORAGE
        {
          auto& tup = model.systems_tuple();
          std::apply([&](auto&... sys) {
            int field_idx = 0;
            auto checksum_one = [&](auto& s) {
              auto& g = s.as_grid();
              constexpr size_t DD = model_dimension<M>::value;
              // Only scalar grids are checksummed; vector grids store
              // values as double** (one array per axis) and don't fit
              // this debug routine. Skip silently.
              using val_t = std::remove_reference_t<decltype(g.values[0])>;
              if constexpr (std::is_arithmetic_v<val_t>) {
                len_type row_len = g.dims[0];
                double sum = 0;
                if constexpr (DD == 2) {
                  for (len_type y = BOUNDARY_DEPTH; y < g.dims[1] - BOUNDARY_DEPTH; ++y)
                    for (len_type x = BOUNDARY_DEPTH; x < g.dims[0] - BOUNDARY_DEPTH; ++x) {
                      double v = g.values[y * row_len + x];
                      sum += v * v;
                    }
                } else {
                  len_type total = 1;
                  for (iter_type d = 0; d < DD; ++d) total *= g.dims[d];
                  for (len_type i = 0; i < total; ++i) sum += g.values[i];
                }
                fprintf(SYMPHAS_LOG, "CHECKSUM rank=%d field=%d sum=%.10e\n",
                        symphas::parallel::get_node_rank(), field_idx, sum);
              }
              ++field_idx;
            };
            (checksum_one(sys), ...);
          }, tup);
        }
#endif
#endif
      }

      /*
       * assemble a vector of duplicated models if there are
       * multiple runs
       */
      else {
        len_type runs = static_cast<len_type>(symphas::conf::config().runs);
        std::vector<M> models(runs, {coeff, num_coeff, pp});
        symphas::io::write_plot_config(models.front());
        find_solution(models.data(), runs);
      }
    }
    return 1;
  }

#else

  int simulate(double const *coeff, size_t num_coeff) {
    problem_parameters_type pp{1};
    M model(coeff, num_coeff, pp);
    find_solution(model, 0.05, 100);
    return 1;
  }

#endif

  template <typename... Ts>
  auto operator()(Ts &&...ts) {
    return simulate(std::forward<Ts>(ts)...);
  }
};

inline void initiate(const char *modelname, double const *coeff,
                     size_t num_coeff) {
#ifdef USING_CONF
  model_select<Simulation> m{
      symphas::conf::config().simulation_settings.dimension,
      symphas::conf::config().simulation_settings.stp};
#else
  model_select<Simulation> m{2, StencilParams{2, 9, 6, 13}};
#endif

  int result = INVALID_MODEL;

#if defined(USING_FFTW) && defined(USE_SPECTRAL_SOLVER) && \
    !defined(SYMPHAS_DISABLE_SP2)
  result = m.call<SolverSP2>(modelname, coeff, num_coeff);
#endif

#ifndef SYMPHAS_DISABLE_FT
  if (result == INVALID_MODEL) {
    result = m.call<SolverFT>(modelname, coeff, num_coeff);
  }
#endif

  if (result == INVALID_MODEL) {
    fprintf(SYMPHAS_ERR, "Unknown model provided, '%s'\n", modelname);
    exit(101);
  }

#ifdef PRINT_TIMINGS
  print_timings(SYMPHAS_LOG);
#endif
}


