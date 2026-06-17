
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
 * PURPOSE: Defines functions for MPI functionality.
 *
 * ***************************************************************************
 */

#pragma once
#include "definitions.h"

#define SYMPHAS_MPI_HOST_RANK 0
#ifdef USING_MPI

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef SYMPHAS_MPI_PROFILE
#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <string>
#endif

#ifdef _MSC_VER
#include <process.h>
#else
#include <unistd.h>
#endif

namespace symphas {
namespace parallel {}
}  // namespace symphas

namespace symphas::parallel {

// ===========================================================================
// Phase 3.5 — per-step phase timers (CMake-gated by SYMPHAS_MPI_PROFILE).
// ===========================================================================
//
// Lightweight RAII timer that accumulates wall time into a named bucket.
// When SYMPHAS_MPI_PROFILE is not defined the timer expands to a no-op so
// the production build sees zero overhead. Buckets are dumped to stderr
// at MPI_Finalize (or via mpi_profile_dump()).
//
// Buckets are explicit string keys to keep the implementation trivial and
// thread-safe-by-construction (lookup happens at start; the per-bucket
// accumulator is a single atomic). Suggested keys for the per-step path:
//   "halo_y"          — Y-direction Isend/Irecv + Waitall
//   "halo_x"          — X-direction Isend/Irecv + Waitall (no-op for Px=1)
//   "halo_total"      — sum of halo_y + halo_x (wraps both)
//   "io_sync"         — sync_all_slabs (per I/O event)
//   "io_gather"       — gather_field_to_host (per I/O event)

#ifdef SYMPHAS_MPI_PROFILE
struct mpi_phase_bucket {
  std::atomic<long long> nanoseconds{0};
  std::atomic<long long> count{0};
};

inline std::map<std::string, mpi_phase_bucket>& mpi_profile_buckets() {
  static std::map<std::string, mpi_phase_bucket> buckets;
  return buckets;
}

inline std::mutex& mpi_profile_mutex() {
  static std::mutex m;
  return m;
}

inline mpi_phase_bucket& mpi_profile_get(const char* name) {
  std::lock_guard<std::mutex> lk(mpi_profile_mutex());
  return mpi_profile_buckets()[name];
}

struct mpi_profile_scope {
  mpi_phase_bucket& b;
  std::chrono::high_resolution_clock::time_point t0;
  mpi_profile_scope(const char* name)
      : b(mpi_profile_get(name)),
        t0(std::chrono::high_resolution_clock::now()) {}
  ~mpi_profile_scope() {
    auto dt = std::chrono::high_resolution_clock::now() - t0;
    b.nanoseconds.fetch_add(
        std::chrono::duration_cast<std::chrono::nanoseconds>(dt).count(),
        std::memory_order_relaxed);
    b.count.fetch_add(1, std::memory_order_relaxed);
  }
};

inline void mpi_profile_dump(FILE* fp = stderr) {
  std::lock_guard<std::mutex> lk(mpi_profile_mutex());
  int rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(fp, "\n[MPI_PROFILE rank=%d]  %-22s %12s %12s %14s\n",
          rank, "phase", "calls", "total_ms", "mean_us");
  for (auto& kv : mpi_profile_buckets()) {
    long long ns = kv.second.nanoseconds.load();
    long long n = kv.second.count.load();
    double total_ms = ns / 1.0e6;
    double mean_us = (n > 0) ? (ns / 1.0e3) / n : 0.0;
    fprintf(fp, "[MPI_PROFILE rank=%d]  %-22s %12lld %12.3f %14.3f\n",
            rank, kv.first.c_str(), n, total_ms, mean_us);
  }
  fflush(fp);
}

#define SYMPHAS_MPI_PROFILE_SCOPE_CONCAT(a, b) a##b
#define SYMPHAS_MPI_PROFILE_SCOPE_NAME(line) SYMPHAS_MPI_PROFILE_SCOPE_CONCAT(_mpi_prof_scope_, line)
#define SYMPHAS_MPI_PROFILE_SCOPE(name) \
    ::symphas::parallel::mpi_profile_scope SYMPHAS_MPI_PROFILE_SCOPE_NAME(__LINE__)(name)
#define SYMPHAS_MPI_PROFILE_DUMP() ::symphas::parallel::mpi_profile_dump()
#else
#define SYMPHAS_MPI_PROFILE_SCOPE(name) ((void)0)
#define SYMPHAS_MPI_PROFILE_DUMP()      ((void)0)
inline void mpi_profile_dump(FILE* /*fp*/ = nullptr) {}
#endif

// ===========================================================================
inline int get_node_pid() {
#ifdef _MSC_VER
  return _getpid();
#else
  return getpid();
#endif
}

inline int get_node_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

inline bool is_host_node() { return get_node_rank() == SYMPHAS_MPI_HOST_RANK; }

inline bool is_host_node(int rank) { return rank == SYMPHAS_MPI_HOST_RANK; }

inline int get_num_nodes() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

inline std::pair<iter_type, iter_type> get_index_range(len_type num_fields,
                                                       iter_type rank) {
  int N = get_num_nodes();
  int size = num_fields / N;
  int R = num_fields % N;

  int lower = size * rank + std::min(R, rank);
  int upper = size * (rank + 1) + std::min(R, rank + 1);

  return {lower, upper};
}

inline std::pair<iter_type, iter_type> get_index_range(len_type num_fields) {
  int rank = get_node_rank();
  return get_index_range(num_fields, rank);
}

inline bool index_in_node_range(int index, int num_fields) {
  auto [lower, upper] = get_index_range(num_fields);
  return (index >= lower && index < upper);
}

struct info_type {
  int rank;
  int index;
  int num_fields;
  int owning_node;

  info_type(int index, int num_fields, int owning_node = -1)
      : rank{get_node_rank()},
        index{index},
        num_fields{num_fields},
        owning_node{(owning_node >= 0) ? owning_node
                                       : num_fields / get_num_nodes()} {}

  info_type() : rank{0}, index{0}, num_fields{1} {}

  inline bool is_in_node() const {
    if (owning_node >= 0)
      return get_node_rank() == owning_node;
    else
      return false;
  }
};

//! Decode the 2-D Cartesian decomposition (Px, Py) for the current MPI run.
//! Resolved once per process and cached. Px is the X-axis (column) factor,
//! Py is the Y-axis (row) factor; Px*Py == num_ranks.
//!
//! Resolution order:
//!   1. Env var SYMPHAS_MPI_DIMS=PxxPy (e.g. "2x4"). Must satisfy Px*Py==N.
//!   2. Env var SYMPHAS_MPI_DECOMP=1d   -> {1, N}    (legacy slab decomp)
//!   3. Default                          -> MPI_Dims_create(N, 2, ...)
//!      (Phase 2.4: auto 2-D balanced decomposition).
//!
//! Phase 2.4: default is now auto-balanced via MPI_Dims_create. The
//! larger factor is placed along Y (Py >= Px) so the X strip stays short
//! when the factorization is non-square; this keeps the non-contiguous
//! X halo communication minimal. Set SYMPHAS_MPI_DECOMP=1d to recover
//! the legacy 1-D Y-slab path (e.g. for benchmark comparisons).
inline std::pair<int, int> get_mpi_dims_cart() {
  static int cached_px = -1;
  static int cached_py = -1;
  if (cached_px > 0) return {cached_px, cached_py};

  int N = get_num_nodes();
  int px = 1, py = N;
  bool resolved = false;

  const char* env_dims = std::getenv("SYMPHAS_MPI_DIMS");
  if (env_dims && env_dims[0]) {
    int pa = 0, pb = 0;
    if (sscanf(env_dims, "%dx%d", &pa, &pb) == 2 && pa > 0 && pb > 0 && pa * pb == N) {
      px = pa;
      py = pb;
      resolved = true;
    } else {
      fprintf(stderr,
              "SYMPHAS_MPI_DIMS='%s' invalid for N=%d (expected PxxPy with Px*Py==N); falling back to auto\n",
              env_dims, N);
    }
  }

  if (!resolved) {
    const char* env_decomp = std::getenv("SYMPHAS_MPI_DECOMP");
    if (env_decomp && (strcmp(env_decomp, "1d") == 0)) {
      px = 1;
      py = N;
    } else {
      // Phase 2.4 default: auto-balanced 2-D decomposition.
      int dims2[2] = {0, 0};
      MPI_Dims_create(N, 2, dims2);
      px = std::min(dims2[0], dims2[1]);
      py = std::max(dims2[0], dims2[1]);
    }
  }

  cached_px = px;
  cached_py = py;
  return {px, py};
}

//! Domain decomposition information for MPI slab / 2-D Cartesian
//! decomposition. The legacy 1-D Y-slab decomposition corresponds to
//! dims_cart == {1, num_ranks}. The 2-D decomposition with Px > 1 splits
//! both X and Y axes simultaneously (Px*Py == num_ranks).
//!
//! Memory layout for 2-D blocks: the host grid storage covers the full
//! global domain on every rank, but only the rank's local block interior
//! is updated each step. Block extents (interior coords):
//!   x: [local_x_start, local_x_end)    nx = local_x_end - local_x_start
//!   y: [local_y_start, local_y_end)    ny = local_y_end - local_y_start
//! Halo strips of width boundary_depth surround the block on the four
//! sides; corners are filled by ordering the exchange Y-first then X
//! (the X-strip extends across the freshly-updated Y halos).
template <size_t D>
struct domain_info {
  int rank;
  int num_ranks;
  len_type global_dims[D];
  len_type interior_dims[D];

  // 2-D Cartesian topology metadata (unused for D==1; sized for D==2).
  // dims_cart[0] = Px (X-axis factor), dims_cart[1] = Py (Y-axis factor).
  // coords[0]    = ix in [0, Px),       coords[1]    = iy in [0, Py).
  // Cart rank arithmetic: rank == ix * Py + iy (row-major over Py first).
  int dims_cart[2];
  int coords[2];

  // Local block extents (interior coords) within the 2-D block.
  // For Px=1 (legacy slab), local_x_start=0 and local_nx=interior_dims[0].
  len_type local_x_start;
  len_type local_x_end;
  len_type local_nx;
  len_type local_y_start;
  len_type local_y_end;
  len_type local_ny;

  // Local block extents along Z (3-D Z-slab decomposition only). For D<3 these
  // are unused. The 3-D path decomposes ONLY the slowest (Z) axis across all
  // ranks, keeping the full X and Y extent on every rank — the direct analog of
  // the bit-exact 1-D Y-slab (Px==1) used in 2-D.
  len_type local_z_start;
  len_type local_z_end;
  len_type local_nz;

  // Periodic neighbor ranks for Y halos (always valid for periodic BCs).
  int neighbor_below;     // -Y
  int neighbor_above;     // +Y
  // Periodic neighbor ranks for X halos. MPI_PROC_NULL when Px==1; the
  // X-direction halo is then provided by the local periodic boundary
  // updater on every rank (no MPI traffic). MPI calls with PROC_NULL
  // complete immediately as no-ops, so exchange_halos is correct in
  // either case.
  int neighbor_left;      // -X
  int neighbor_right;     // +X

  // Periodic neighbor ranks for Z halos (3-D Z-slab only). For the Z-slab the
  // X and Y halos are handled locally by the periodic boundary updater (every
  // rank owns the full X,Y extent), so only the Z halos cross the wire.
  int neighbor_front;     // -Z
  int neighbor_back;      // +Z

  len_type boundary_depth;

  // When true, each rank allocates only (local_n+2b) × (local_n+2b)
  // storage and the rank's interior begins at LOCAL coords
  // (bdepth, bdepth) within its own buffer. When false (legacy), every
  // rank holds the full N² grid and the interior block begins at GLOBAL
  // coords (bdepth + local_x_start, bdepth + local_y_start). The flag
  // is consulted by exchange_halos / gather / sync to compute the
  // correct row stride and starting offsets.
  bool local_storage;

  domain_info()
      : rank{}, num_ranks{}, global_dims{}, interior_dims{},
        dims_cart{1, 1}, coords{0, 0},
        local_x_start{}, local_x_end{}, local_nx{},
        local_y_start{}, local_y_end{}, local_ny{},
        local_z_start{}, local_z_end{}, local_nz{},
        neighbor_below{MPI_PROC_NULL}, neighbor_above{MPI_PROC_NULL},
        neighbor_left{MPI_PROC_NULL}, neighbor_right{MPI_PROC_NULL},
        neighbor_front{MPI_PROC_NULL}, neighbor_back{MPI_PROC_NULL},
        boundary_depth{}, local_storage{false} {}

  domain_info(const len_type (&dims)[D], len_type bdepth)
      : rank{get_node_rank()}, num_ranks{get_num_nodes()},
        global_dims{}, interior_dims{},
        dims_cart{1, get_num_nodes()},
        coords{0, get_node_rank()},
        local_x_start{}, local_x_end{}, local_nx{},
        local_y_start{}, local_y_end{}, local_ny{},
        local_z_start{}, local_z_end{}, local_nz{},
        neighbor_below{MPI_PROC_NULL}, neighbor_above{MPI_PROC_NULL},
        neighbor_left{MPI_PROC_NULL}, neighbor_right{MPI_PROC_NULL},
        neighbor_front{MPI_PROC_NULL}, neighbor_back{MPI_PROC_NULL},
        boundary_depth{bdepth}, local_storage{false} {
    for (iter_type i = 0; i < D; ++i) {
      global_dims[i] = dims[i];
      interior_dims[i] = dims[i] - 2 * bdepth;
    }

    if constexpr (D == 3) {
      // 3-D Z-slab: decompose ONLY the slowest (Z) axis across all ranks;
      // every rank keeps the full X and Y extent. This is the direct analog of
      // the bit-exact 1-D Y-slab (Px==1) 2-D path: the rank's block is a
      // contiguous range of Z-planes, so halo exchange and I/O sync are plain
      // contiguous sends (no MPI_Type_vector). X and Y periodicity is handled
      // locally by the periodic boundary updater on every rank.
      local_x_start = 0;
      local_x_end = interior_dims[0];
      local_nx = interior_dims[0];
      local_y_start = 0;
      local_y_end = interior_dims[1];
      local_ny = interior_dims[1];

      len_type nz_total = interior_dims[2];
      len_type z_base = nz_total / num_ranks;
      len_type z_rem = nz_total % num_ranks;
      local_z_start = z_base * rank + std::min((len_type)rank, z_rem);
      local_z_end = z_base * (rank + 1) + std::min((len_type)(rank + 1), z_rem);
      local_nz = local_z_end - local_z_start;

      // Periodic Z neighbors (front = -Z, back = +Z).
      neighbor_front = (rank - 1 + num_ranks) % num_ranks;
      neighbor_back = (rank + 1) % num_ranks;
    } else if constexpr (D == 2) {
      auto [px, py] = get_mpi_dims_cart();
      dims_cart[0] = px;
      dims_cart[1] = py;
      // Row-major mapping: rank = ix * Py + iy.
      coords[0] = rank / py;     // ix
      coords[1] = rank % py;     // iy

      // X-axis split.
      len_type nx_total = interior_dims[0];
      len_type x_base = nx_total / px;
      len_type x_rem = nx_total % px;
      local_x_start = x_base * coords[0] + std::min((len_type)coords[0], x_rem);
      local_x_end = x_base * (coords[0] + 1) + std::min((len_type)(coords[0] + 1), x_rem);
      local_nx = local_x_end - local_x_start;

      // Y-axis split.
      len_type ny_total = interior_dims[1];
      len_type y_base = ny_total / py;
      len_type y_rem = ny_total % py;
      local_y_start = y_base * coords[1] + std::min((len_type)coords[1], y_rem);
      local_y_end = y_base * (coords[1] + 1) + std::min((len_type)(coords[1] + 1), y_rem);
      local_ny = local_y_end - local_y_start;

      // Periodic Cartesian neighbors.
      int ix_left = (coords[0] - 1 + px) % px;
      int ix_right = (coords[0] + 1) % px;
      int iy_below = (coords[1] - 1 + py) % py;
      int iy_above = (coords[1] + 1) % py;
      // Y neighbors always set (periodic).
      neighbor_below = coords[0] * py + iy_below;
      neighbor_above = coords[0] * py + iy_above;
      // X neighbors set only when Px > 1; otherwise leave as PROC_NULL so
      // the X exchange is skipped and the periodic boundary updater on
      // every rank handles X periodicity locally.
      if (px > 1) {
        neighbor_left = ix_left * py + coords[1];
        neighbor_right = ix_right * py + coords[1];
      }
    } else {
      // D < 2: keep historical 1-D behavior (rare/untested path).
      len_type total = interior_dims[0];
      len_type base = total / num_ranks;
      len_type rem = total % num_ranks;
      local_y_start = base * rank + std::min((len_type)rank, rem);
      local_y_end = base * (rank + 1) + std::min((len_type)(rank + 1), rem);
      local_ny = local_y_end - local_y_start;
      neighbor_below = (rank > 0) ? rank - 1 : num_ranks - 1;
      neighbor_above = (rank < num_ranks - 1) ? rank + 1 : 0;
    }
  }

  len_type halo_strip_len() const {
    return boundary_depth * global_dims[0];
  }
};

template <typename T>
inline MPI_Datatype mpi_type() {
  if constexpr (std::is_same_v<T, double>) {
    return MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, float>) {
    return MPI_FLOAT;
  } else if constexpr (std::is_same_v<T, int>) {
    return MPI_INT;
  } else {
    static MPI_Datatype dt = MPI_DATATYPE_NULL;
    if (dt == MPI_DATATYPE_NULL) {
      MPI_Type_contiguous(sizeof(T), MPI_BYTE, &dt);
      MPI_Type_commit(&dt);
    }
    return dt;
  }
}

//! Cache an MPI_Type_vector for a strided X-column halo strip, keyed by
//! (T, strip_count, bdepth, row_len). Built once per unique geometry per
//! process and reused across every exchange. This removes the per-step
//! Type_create / Type_commit / Type_free overhead that otherwise dominates
//! the X-halo cost (~5–10 µs per step). Phase 4 stretch (H2/H5 in
//! docs/MPI_HOTSPOTS.md). The cached datatype is intentionally never
//! MPI_Type_free'd — handles persist for the lifetime of the process.
template <typename T>
inline MPI_Datatype get_cached_col_strip(int strip_count, int bdepth, int row_len) {
  struct key_t {
    int strip_count;
    int bdepth;
    int row_len;
    bool operator==(const key_t& o) const {
      return strip_count == o.strip_count && bdepth == o.bdepth && row_len == o.row_len;
    }
  };
  // Small linear cache (typical run has ≤2 distinct geometries per T).
  static thread_local std::vector<std::pair<key_t, MPI_Datatype>> cache;
  key_t k{strip_count, bdepth, row_len};
  for (auto const& e : cache) {
    if (e.first == k) return e.second;
  }
  MPI_Datatype dt;
  MPI_Type_vector(strip_count, bdepth, row_len, mpi_type<T>(), &dt);
  MPI_Type_commit(&dt);
  cache.emplace_back(k, dt);
  return dt;
}

//! Cache an MPI_Type_vector for the Y-direction halo strip when Px > 1,
//! keyed by (T, bdepth, local_nx, row_len). Without this the Y exchange
//! sends bdepth × global_dims[0] doubles per step (full row width) even
//! though only bdepth × local_nx of those columns are actually owned by
//! the rank — a factor-Px bandwidth waste. The X-halo columns at the
//! Y-halo rows are subsequently populated by the X exchange (whose strip
//! count already covers local_ny + 2*bdepth rows for corner correctness).
//! Cache entries are intentionally never MPI_Type_free'd.
template <typename T>
inline MPI_Datatype get_cached_y_strip(int bdepth, int local_nx, int row_len) {
  struct key_t {
    int bdepth;
    int local_nx;
    int row_len;
    bool operator==(const key_t& o) const {
      return bdepth == o.bdepth && local_nx == o.local_nx && row_len == o.row_len;
    }
  };
  static thread_local std::vector<std::pair<key_t, MPI_Datatype>> cache;
  key_t k{bdepth, local_nx, row_len};
  for (auto const& e : cache) {
    if (e.first == k) return e.second;
  }
  MPI_Datatype dt;
  // count = bdepth rows, blocklength = local_nx, stride = row_len.
  MPI_Type_vector(bdepth, local_nx, row_len, mpi_type<T>(), &dt);
  MPI_Type_commit(&dt);
  cache.emplace_back(k, dt);
  return dt;
}

//! Exchange halos for the 2-D slab/block decomposition.
//!
//! Implements two-step exchange (Y first, then X) so that the X-direction
//! halo strip extends across the freshly-updated Y halos and corner cells
//! are correctly populated. The X-direction exchange is automatically a
//! no-op when Px == 1 (neighbor_left/right == MPI_PROC_NULL); MPI calls
//! with PROC_NULL complete immediately and produce no traffic, so this
//! is correct (and bit-identical) for the legacy 1-D Y-slab case.
template <typename T, size_t D>
inline void exchange_halos(T* values, const domain_info<D>& dinfo) {
  static_assert(D == 2 || D == 3,
                "Halo exchange implemented for 2D and 3D");
  SYMPHAS_MPI_PROFILE_SCOPE("halo_total");

  if constexpr (D == 3) {
    // 3-D Z-slab: contiguous front/back (−Z/+Z) halo exchange. The rank owns a
    // contiguous range of Z-planes spanning the full X×Y extent, so a Z-halo of
    // width bdepth is exactly bdepth whole XY-planes of contiguous memory — the
    // direct analog of the 2-D Px==1 contiguous full-row send. X and Y
    // periodicity is filled locally by the periodic boundary updater on every
    // rank (each rank owns the full X,Y extent), so no X/Y MPI traffic.
    len_type bdepth = dinfo.boundary_depth;
    len_type plane = dinfo.global_dims[0] * dinfo.global_dims[1];
    len_type strip = bdepth * plane;
    len_type my_first = bdepth + dinfo.local_z_start;
    len_type my_last = bdepth + dinfo.local_z_start + dinfo.local_nz;
    T* send_front = values + my_first * plane;
    T* send_back = values + (my_last - bdepth) * plane;
    T* recv_front = values + (my_first - bdepth) * plane;
    T* recv_back = values + my_last * plane;
    MPI_Request requests[4];
    MPI_Isend(send_front, static_cast<int>(strip), mpi_type<T>(),
              dinfo.neighbor_front, 0, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_front, static_cast<int>(strip), mpi_type<T>(),
              dinfo.neighbor_front, 1, MPI_COMM_WORLD, &requests[1]);
    MPI_Isend(send_back, static_cast<int>(strip), mpi_type<T>(),
              dinfo.neighbor_back, 1, MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(recv_back, static_cast<int>(strip), mpi_type<T>(),
              dinfo.neighbor_back, 0, MPI_COMM_WORLD, &requests[3]);
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    return;
  } else {
  len_type bdepth = dinfo.boundary_depth;
  // Row stride and per-rank origin depend on storage layout:
  //   - Legacy: every rank holds full N² → row_len = global_dims[0],
  //     rank's interior block at GLOBAL coords (bdepth + local_x_start,
  //     bdepth + local_y_start).
  //   - Local: each rank's buffer is (local_nx + 2b) × (local_ny + 2b) →
  //     row_len = local_nx + 2*bdepth, rank's interior begins at LOCAL
  //     coords (bdepth, bdepth) inside its own storage.
  len_type row_len = dinfo.local_storage
                         ? (dinfo.local_nx + 2 * bdepth)
                         : dinfo.global_dims[0];
  len_type rank_x0 = dinfo.local_storage ? 0 : dinfo.local_x_start;
  len_type rank_y0 = dinfo.local_storage ? 0 : dinfo.local_y_start;

  // -------- Y-direction exchange --------
  // For Px==1 each rank owns the full global row width (the legacy 1-D
  // Y-slab decomposition), so the contiguous bdepth × global_dims[0] send
  // is correct and tight. For Px>1 only the rank's local_nx columns need
  // to cross the wire; the X-halo columns at the Y-halo rows are then
  // filled by the X exchange below (which already includes the Y halo
  // rows in its strip for corner correctness).
  {
    SYMPHAS_MPI_PROFILE_SCOPE("halo_y");
    len_type my_first_row = bdepth + rank_y0;
    len_type my_last_row = bdepth + rank_y0 + dinfo.local_ny;
    MPI_Request requests[4];

    if (dinfo.dims_cart[0] == 1) {
      // Px == 1: contiguous full-row send (bit-identical legacy path).
      len_type strip_len = bdepth * row_len;
      T* send_down = values + my_first_row * row_len;
      T* send_up = values + (my_last_row - bdepth) * row_len;
      T* recv_below = values + (my_first_row - bdepth) * row_len;
      T* recv_above = values + my_last_row * row_len;
      MPI_Isend(send_down, strip_len, mpi_type<T>(),
                dinfo.neighbor_below, 0, MPI_COMM_WORLD, &requests[0]);
      MPI_Irecv(recv_below, strip_len, mpi_type<T>(),
                dinfo.neighbor_below, 1, MPI_COMM_WORLD, &requests[1]);
      MPI_Isend(send_up, strip_len, mpi_type<T>(),
                dinfo.neighbor_above, 1, MPI_COMM_WORLD, &requests[2]);
      MPI_Irecv(recv_above, strip_len, mpi_type<T>(),
                dinfo.neighbor_above, 0, MPI_COMM_WORLD, &requests[3]);
    } else {
      // Px > 1: strided send carrying only local_nx columns at the rank's
      // X extent. Reduces Y-halo bandwidth from O(global_dims[0]) to
      // O(local_nx) per direction (factor Px improvement).
      MPI_Datatype y_strip = get_cached_y_strip<T>(
          static_cast<int>(bdepth),
          static_cast<int>(dinfo.local_nx),
          static_cast<int>(row_len));
      len_type x0 = bdepth + rank_x0;
      T* send_down = values + my_first_row * row_len + x0;
      T* send_up = values + (my_last_row - bdepth) * row_len + x0;
      T* recv_below = values + (my_first_row - bdepth) * row_len + x0;
      T* recv_above = values + my_last_row * row_len + x0;
      MPI_Isend(send_down, 1, y_strip,
                dinfo.neighbor_below, 0, MPI_COMM_WORLD, &requests[0]);
      MPI_Irecv(recv_below, 1, y_strip,
                dinfo.neighbor_below, 1, MPI_COMM_WORLD, &requests[1]);
      MPI_Isend(send_up, 1, y_strip,
                dinfo.neighbor_above, 1, MPI_COMM_WORLD, &requests[2]);
      MPI_Irecv(recv_above, 1, y_strip,
                dinfo.neighbor_above, 0, MPI_COMM_WORLD, &requests[3]);
    }
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
  }

  // -------- X-direction exchange (column strip; no-op when Px==1) --------
  if (dinfo.dims_cart[0] > 1) {
    SYMPHAS_MPI_PROFILE_SCOPE("halo_x");
    // Column strip extent in Y: include both Y halo bands so corners are
    // populated by the X exchange. count = local_ny + 2*bdepth rows;
    // each row contributes bdepth contiguous elements; stride = row_len.
    len_type strip_count = dinfo.local_ny + 2 * bdepth;
    // Phase 4: cached MPI_Type_vector (built once per geometry, never freed).
    MPI_Datatype col_strip = get_cached_col_strip<T>(
        static_cast<int>(strip_count),
        static_cast<int>(bdepth),
        static_cast<int>(row_len));

    // Column origins: y at (rank_y0) i.e. starting at the bottom Y
    // halo of the block, x at the four locations below.
    len_type y0 = rank_y0;  // includes bottom Y halo because
                            // first row index = bdepth + y0 - bdepth = y0
    len_type x_send_left = bdepth + rank_x0;                        // left interior
    len_type x_send_right = bdepth + rank_x0 + dinfo.local_nx
                            - bdepth;                               // right interior
    len_type x_recv_left = rank_x0;                                 // left halo
    len_type x_recv_right = bdepth + rank_x0 + dinfo.local_nx;      // right halo

    T* send_left = values + y0 * row_len + x_send_left;
    T* send_right = values + y0 * row_len + x_send_right;
    T* recv_left = values + y0 * row_len + x_recv_left;
    T* recv_right = values + y0 * row_len + x_recv_right;

    MPI_Request requests[4];
    MPI_Isend(send_left, 1, col_strip,
              dinfo.neighbor_left, 2, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_left, 1, col_strip,
              dinfo.neighbor_left, 3, MPI_COMM_WORLD, &requests[1]);
    MPI_Isend(send_right, 1, col_strip,
              dinfo.neighbor_right, 3, MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(recv_right, 1, col_strip,
              dinfo.neighbor_right, 2, MPI_COMM_WORLD, &requests[3]);
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    // No MPI_Type_free: datatype is cached for the lifetime of the run.
  }
  }  // end else (D == 2)
}

//! Recompute another rank's local (x_start, x_end, y_start, y_end) given
//! the current dinfo. Used by sync_all_slabs and gather_field_to_host to
//! address remote ranks' blocks without rebuilding a full domain_info.
template <size_t D>
inline void rank_block_extents(const domain_info<D>& dinfo, int r,
                               len_type& rx0, len_type& rx1,
                               len_type& ry0, len_type& ry1) {
  int px = dinfo.dims_cart[0];
  int py = dinfo.dims_cart[1];
  int rix = r / py;
  int riy = r % py;
  len_type nx_total = dinfo.interior_dims[0];
  len_type x_base = nx_total / px;
  len_type x_rem = nx_total % px;
  rx0 = x_base * rix + std::min((len_type)rix, x_rem);
  rx1 = x_base * (rix + 1) + std::min((len_type)(rix + 1), x_rem);
  len_type ny_total = dinfo.interior_dims[1];
  len_type y_base = ny_total / py;
  len_type y_rem = ny_total % py;
  ry0 = y_base * riy + std::min((len_type)riy, y_rem);
  ry1 = y_base * (riy + 1) + std::min((len_type)(riy + 1), y_rem);
}

template <typename T, size_t D>
inline void sync_all_slabs(T* values, const domain_info<D>& dinfo) {
  static_assert(D == 2 || D == 3, "Slab sync implemented for 2D and 3D");
  SYMPHAS_MPI_PROFILE_SCOPE("io_sync");

  if constexpr (D == 3) {
    // 3-D Z-slab: each rank's block is a contiguous range of whole XY-planes,
    // so the sync is a plain contiguous broadcast per rank (mirrors the 2-D
    // Px==1 contiguous path). X and Y periodic ghosts are already filled in
    // each plane by the per-rank periodic boundary updater; only the Z (front/
    // back) ghost planes need a wrap refresh afterwards.
    len_type bdepth = dinfo.boundary_depth;
    len_type plane = dinfo.global_dims[0] * dinfo.global_dims[1];
    len_type nz_total = dinfo.interior_dims[2];
    for (int r = 0; r < dinfo.num_ranks; ++r) {
      len_type z_base = nz_total / dinfo.num_ranks;
      len_type z_rem = nz_total % dinfo.num_ranks;
      len_type rz0 = z_base * r + std::min((len_type)r, z_rem);
      len_type rz1 = z_base * (r + 1) + std::min((len_type)(r + 1), z_rem);
      len_type r_nz = rz1 - rz0;
      T* block_ptr = values + (bdepth + rz0) * plane;
      MPI_Bcast(block_ptr, static_cast<int>(r_nz * plane), mpi_type<T>(), r,
                MPI_COMM_WORLD);
    }
    // Refresh Z periodic ghost planes (front = -Z, back = +Z).
    for (len_type b = 0; b < bdepth; ++b) {
      T* front_ghost = values + b * plane;
      T* front_src = values + (bdepth + nz_total - bdepth + b) * plane;
      std::copy(front_src, front_src + plane, front_ghost);
      T* back_ghost = values + (bdepth + nz_total + b) * plane;
      T* back_src = values + (bdepth + b) * plane;
      std::copy(back_src, back_src + plane, back_ghost);
    }
    return;
  }

  len_type row_len = dinfo.global_dims[0];
  len_type bdepth = dinfo.boundary_depth;

  // Each rank broadcasts its local block (interior, no halos) to all others.
  // For Px==1 the block is contiguous (full row * local_ny), preserving the
  // historical fast path. For Px>1 the block is a 2-D rectangle with stride
  // row_len; we use MPI_Type_vector to describe it.
  for (int r = 0; r < dinfo.num_ranks; ++r) {
    len_type rx0, rx1, ry0, ry1;
    rank_block_extents(dinfo, r, rx0, rx1, ry0, ry1);
    len_type r_nx = rx1 - rx0;
    len_type r_ny = ry1 - ry0;
    T* block_ptr = values + (bdepth + ry0) * row_len + (bdepth + rx0);
    if (dinfo.dims_cart[0] == 1) {
      // Contiguous: full row * r_ny rows.
      MPI_Bcast(block_ptr, static_cast<int>(r_ny * row_len),
                mpi_type<T>(), r, MPI_COMM_WORLD);
    } else {
      MPI_Datatype block_t;
      MPI_Type_vector(static_cast<int>(r_ny), static_cast<int>(r_nx),
                      static_cast<int>(row_len), mpi_type<T>(), &block_t);
      MPI_Type_commit(&block_t);
      MPI_Bcast(block_ptr, 1, block_t, r, MPI_COMM_WORLD);
      MPI_Type_free(&block_t);
    }
  }

  // Refresh periodic ghost rows in Y (still needed for I/O completeness).
  // For the 1-D Y-slab path we kept this unchanged. For 2-D blocks the
  // X-axis ghost columns still come from the periodic boundary updater
  // applied by every rank to its local block; after the broadcasts above
  // the full interior is consistent, so we just refresh Y ghost rows.
  len_type ny_interior = dinfo.interior_dims[1];
  for (len_type b = 0; b < bdepth; ++b) {
    T* dst = values + b * row_len;
    T* src = values + (bdepth + ny_interior - bdepth + b) * row_len;
    std::copy(src, src + row_len, dst);
  }
  for (len_type b = 0; b < bdepth; ++b) {
    T* dst = values + (bdepth + ny_interior + b) * row_len;
    T* src = values + (bdepth + b) * row_len;
    std::copy(src, src + row_len, dst);
  }
  // Refresh periodic ghost columns in X for full I/O completeness.
  len_type nx_interior = dinfo.interior_dims[0];
  for (len_type r2 = 0; r2 < dinfo.global_dims[1]; ++r2) {
    T* row_ptr = values + r2 * row_len;
    for (len_type b = 0; b < bdepth; ++b) {
      row_ptr[b] = row_ptr[bdepth + nx_interior - bdepth + b];
      row_ptr[bdepth + nx_interior + b] = row_ptr[bdepth + b];
    }
  }
}

template <typename T, size_t D>
inline void gather_field_to_host(T* values, const domain_info<D>& dinfo) {
  static_assert(D == 2 || D == 3, "Gather implemented for 2D and 3D");
  SYMPHAS_MPI_PROFILE_SCOPE("io_gather");

  if constexpr (D == 3) {
    // 3-D Z-slab gather: host receives each rank's contiguous Z-slab of whole
    // XY-planes. Mirrors the 2-D Px==1 contiguous path.
    len_type bdepth = dinfo.boundary_depth;
    len_type plane = dinfo.global_dims[0] * dinfo.global_dims[1];
    len_type nz_total = dinfo.interior_dims[2];
    if (is_host_node()) {
      for (int r = 1; r < dinfo.num_ranks; ++r) {
        len_type z_base = nz_total / dinfo.num_ranks;
        len_type z_rem = nz_total % dinfo.num_ranks;
        len_type rz0 = z_base * r + std::min((len_type)r, z_rem);
        len_type rz1 = z_base * (r + 1) + std::min((len_type)(r + 1), z_rem);
        len_type r_nz = rz1 - rz0;
        T* block_ptr = values + (bdepth + rz0) * plane;
        MPI_Recv(block_ptr, static_cast<int>(r_nz * plane), mpi_type<T>(), r,
                 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    } else {
      T* my_block = values + (bdepth + dinfo.local_z_start) * plane;
      MPI_Send(my_block, static_cast<int>(dinfo.local_nz * plane), mpi_type<T>(),
               SYMPHAS_MPI_HOST_RANK, 100, MPI_COMM_WORLD);
    }
    return;
  }

  len_type row_len = dinfo.global_dims[0];
  len_type bdepth = dinfo.boundary_depth;
  if (is_host_node()) {
    for (int r = 1; r < dinfo.num_ranks; ++r) {
      len_type rx0, rx1, ry0, ry1;
      rank_block_extents(dinfo, r, rx0, rx1, ry0, ry1);
      len_type r_nx = rx1 - rx0;
      len_type r_ny = ry1 - ry0;
      T* block_ptr = values + (bdepth + ry0) * row_len + (bdepth + rx0);
      if (dinfo.dims_cart[0] == 1) {
        MPI_Recv(block_ptr, static_cast<int>(r_ny * row_len), mpi_type<T>(),
                 r, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Datatype block_t;
        MPI_Type_vector(static_cast<int>(r_ny), static_cast<int>(r_nx),
                        static_cast<int>(row_len), mpi_type<T>(), &block_t);
        MPI_Type_commit(&block_t);
        MPI_Recv(block_ptr, 1, block_t,
                 r, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Type_free(&block_t);
      }
    }
  } else {
    len_type my_x0 = dinfo.local_x_start;
    len_type my_y0 = dinfo.local_y_start;
    len_type my_nx = dinfo.local_nx;
    len_type my_ny = dinfo.local_ny;
    T* my_block = values + (bdepth + my_y0) * row_len + (bdepth + my_x0);
    if (dinfo.dims_cart[0] == 1) {
      MPI_Send(my_block, static_cast<int>(my_ny * row_len), mpi_type<T>(),
               SYMPHAS_MPI_HOST_RANK, 100, MPI_COMM_WORLD);
    } else {
      MPI_Datatype block_t;
      MPI_Type_vector(static_cast<int>(my_ny), static_cast<int>(my_nx),
                      static_cast<int>(row_len), mpi_type<T>(), &block_t);
      MPI_Type_commit(&block_t);
      MPI_Send(my_block, 1, block_t,
               SYMPHAS_MPI_HOST_RANK, 100, MPI_COMM_WORLD);
      MPI_Type_free(&block_t);
    }
  }
}

//! Trait: does a solver system type carry MPI domain-decomposition metadata?
//!
//! SolverSystemFDMPI has a `dinfo` member of type domain_info<D>; the
//! non-MPI replicated SolverSystemFD does not. This distinguishes
//! MPI-decomposed systems (whose per-rank slab must be gathered before
//! I/O) from replicated systems (whose full grid already exists on
//! every rank, so no sync is needed or correct).
template <typename, typename = void>
struct has_mpi_dinfo : std::false_type {};

template <typename S>
struct has_mpi_dinfo<S, std::void_t<decltype(std::declval<S&>().dinfo)>>
    : std::true_type {};

template <typename S>
inline constexpr bool has_mpi_dinfo_v = has_mpi_dinfo<S>::value;

}  // namespace symphas::parallel

#else

namespace symphas::parallel {

inline int get_node_rank() { return SYMPHAS_MPI_HOST_RANK; }

inline bool is_host_node() { return true; }

inline bool is_host_node(int rank) { return true; }

inline int get_num_nodes() { return 1; }

inline std::pair<int, int> get_index_range(size_t num_fields, int rank) {
  return {0, static_cast<int>(num_fields)};
}

inline std::pair<int, int> get_index_range(size_t num_fields) {
  return get_index_range(static_cast<int>(num_fields), get_node_rank());
}

using info_type = size_t;
}  // namespace symphas::parallel

#endif

#ifndef SYMPHAS_MPI_PROFILE_SCOPE
#define SYMPHAS_MPI_PROFILE_SCOPE(name) ((void)0)
#endif
#ifndef SYMPHAS_MPI_PROFILE_DUMP
#define SYMPHAS_MPI_PROFILE_DUMP()      ((void)0)
#endif

namespace symphas {
using multi_thr_info_type = symphas::parallel::info_type;
}