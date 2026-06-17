
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

#include "solversystem.h"

#ifdef USING_CUDA

#include <vector>

#include "boundarysystem.cuh"
#include "systemlib.cuh"

//! The default phase field system.
/*!
 * The representation of a phase field system, storing the values of the
 * order parameter that is defined in a phase field problem.
 *
 * Unless explicitly specified, this phase field system type will be
 * used by the solver. It does not manage boundaries or other data.
 *
 * \tparam T The order parameter type.
 * \tparam D The order parameter dimension.
 */
template <typename T, size_t D>
using SolverSystemCUDA = SystemCUDA<T, D>;

template <typename T, size_t D>
struct SolverSystemFDCUDA : BoundarySystemCUDA<T, D> {
  using BoundarySystemCUDA<T, D>::dims;

  BoundaryGridCUDA<T, D> dframe;  // the working grid for the solver
  SolverSystemFDCUDA(symphas::init_data_type const& tdata,
                     symphas::interval_data_type const& vdata,
                     symphas::b_data_type const& bdata, size_t id = 0)
      : BoundarySystemCUDA<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFDCUDA() : BoundarySystemCUDA<T, D>(), dframe{0} {}
};

template <typename T, size_t D>
struct SolverSystemFDwSDCUDA : RegionalSystemCUDA<T, D> {
  using RegionalSystemCUDA<T, D>::dims;

  RegionalGridCUDA<T, D> dframe;  // the working grid for the solver
  SolverSystemFDwSDCUDA(symphas::init_data_type const& tdata,
                        symphas::interval_data_type const& vdata,
                        symphas::b_data_type const& bdata, size_t id = 0)
      : RegionalSystemCUDA<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFDwSDCUDA() : RegionalSystemCUDA<T, D>(), dframe{0} {}

  inline void update(iter_type index, double time) {
    RegionalSystemCUDA<T, D>::update(index, time);
    dframe.adjust(RegionalSystemCUDA<T, D>::region);
  }
};

DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFDCUDA<T, D>),
                           (BoundaryGridCUDA<T, D>))

DEFINE_BASE_DATA_INHERITED((typename T, size_t D),
                           (SolverSystemFDwSDCUDA<T, D>),
                           (RegionalGridCUDA<T, D>))

// ===========================================================================
// GPU spectral (semi-implicit Fourier) solver system.
//
// Device-resident counterpart of SolverSystemSpectral<scalar_t, D>. The order
// parameter lives on the GPU (GridCUDA values), and the per-step transforms
// run with cuFFT (real-to-complex D2Z forward, complex-to-real Z2D inverse).
// The real axis dims[0] (the contiguous x-axis in SymPhas storage) is the
// half-spectrum axis, exactly matching the CPU layout
// (symphas::dft::length<scalar_t,D> = (dims[0]/2+1)*dims[1][*dims[2]]), so the
// A(k)/B(k) operator arrays built by the SP solver transfer to the device
// unchanged.
// ===========================================================================

#include <cufft.h>

namespace symphas::internal {

//! Scale a real device array in place (normalizes the cuFFT inverse).
__global__ inline void sp_cuda_scale_kernel(scalar_t* v, len_type n,
                                            double scale) {
  len_type i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) v[i] *= scale;
}

//! Spectral linear update: dframe = A * frame_t + B * dframe (per k-mode).
__global__ inline void sp_cuda_spectral_mul_kernel(
    cufftDoubleComplex* dframe, const cufftDoubleComplex* frame_t,
    const cufftDoubleComplex* A, const cufftDoubleComplex* B, len_type n) {
  len_type i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  // complex: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
  cufftDoubleComplex a = A[i], ft = frame_t[i], b = B[i], df = dframe[i];
  double re = (a.x * ft.x - a.y * ft.y) + (b.x * df.x - b.y * df.y);
  double im = (a.x * ft.y + a.y * ft.x) + (b.x * df.y + b.y * df.x);
  dframe[i].x = re;
  dframe[i].y = im;
}

}  // namespace symphas::internal

//! Distributed-memory-free GPU spectral system (single device).
template <size_t D>
struct SolverSystemSpectralCUDA : SystemCUDA<scalar_t, D> {
  using parent_type = SystemCUDA<scalar_t, D>;
  using parent_type::dims;
  using parent_type::values;  // device real field
  using parent_type::len;

  static_assert(D == 1 || D == 2 || D == 3,
                "GPU spectral system supports D=1,2,3");

  len_type transformed_len;        //!< Half-spectrum complex length.
  cufftDoubleComplex* frame_t;     //!< k-space field (device).
  cufftDoubleComplex* dframe;      //!< k-space accumulator (device).
  scalar_t* nl_work;               //!< device scratch holding the field copy.
  cufftHandle p;                   //!< inverse (Z2D) plan.
  cufftHandle p_to_t;              //!< forward (D2Z) plan.
  bool owns;

  //! cuFFT logical dimensions: SymPhas stores x-contiguous and halves dims[0],
  //! so the cuFFT contiguous (last) axis must be dims[0]. Pass axes reversed.
  static void plan_dims(const len_type* d, int* out) {
    if constexpr (D == 1) {
      out[0] = static_cast<int>(d[0]);
    } else if constexpr (D == 2) {
      out[0] = static_cast<int>(d[1]);
      out[1] = static_cast<int>(d[0]);
    } else {
      out[0] = static_cast<int>(d[2]);
      out[1] = static_cast<int>(d[1]);
      out[2] = static_cast<int>(d[0]);
    }
  }

  void make_plans() {
    int n[D];
    plan_dims(dims, n);
    if constexpr (D == 1) {
      CHECK_CUFFT_ERROR(cufftPlan1d(&p_to_t, n[0], CUFFT_D2Z, 1));
      CHECK_CUFFT_ERROR(cufftPlan1d(&p, n[0], CUFFT_Z2D, 1));
    } else if constexpr (D == 2) {
      CHECK_CUFFT_ERROR(cufftPlan2d(&p_to_t, n[0], n[1], CUFFT_D2Z));
      CHECK_CUFFT_ERROR(cufftPlan2d(&p, n[0], n[1], CUFFT_Z2D));
    } else {
      CHECK_CUFFT_ERROR(cufftPlan3d(&p_to_t, n[0], n[1], n[2], CUFFT_D2Z));
      CHECK_CUFFT_ERROR(cufftPlan3d(&p, n[0], n[1], n[2], CUFFT_Z2D));
    }
  }

  SolverSystemSpectralCUDA(symphas::init_data_type const& tdata,
                           symphas::interval_data_type const& vdata,
                           symphas::b_data_type const&, size_t id = 0)
      : parent_type(tdata, vdata, id),
        transformed_len{symphas::dft::length<scalar_t, D>(dims)},
        frame_t{nullptr}, dframe{nullptr}, nl_work{nullptr},
        p{0}, p_to_t{0}, owns{true} {
    CHECK_CUDA_ERROR(cudaMalloc(&frame_t,
                                transformed_len * sizeof(cufftDoubleComplex)));
    CHECK_CUDA_ERROR(cudaMalloc(&dframe,
                                transformed_len * sizeof(cufftDoubleComplex)));
    CHECK_CUDA_ERROR(cudaMalloc(&nl_work, len * sizeof(scalar_t)));
    make_plans();
    // Seed BOTH the k-space snapshot (frame_t) and accumulator (dframe) with
    // the forward transform of the initial field. The time loop calls
    // update() BEFORE the first equation() (see model_iteration); update()
    // copies dframe -> frame_t and inverse-transforms dframe -> values, so
    // dframe must already hold FFT(IC) on entry or the field would be zeroed.
    CHECK_CUFFT_ERROR(cufftExecD2Z(p_to_t, values, frame_t));
    CHECK_CUDA_ERROR(cudaMemcpy(dframe, frame_t,
                                transformed_len * sizeof(cufftDoubleComplex),
                                cudaMemcpyDeviceToDevice));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  SolverSystemSpectralCUDA()
      : parent_type(), transformed_len{0}, frame_t{nullptr}, dframe{nullptr},
        nl_work{nullptr}, p{0}, p_to_t{0}, owns{false} {}

  SolverSystemSpectralCUDA(SolverSystemSpectralCUDA&& o) noexcept
      : SolverSystemSpectralCUDA() {
    swap(*this, o);
  }
  SolverSystemSpectralCUDA& operator=(SolverSystemSpectralCUDA o) {
    swap(*this, o);
    return *this;
  }

  //! Deep copy: clone the device buffers and rebuild the cuFFT plans (plans
  //! are not copyable). Required because Model<...> holds systems by value and
  //! has a copy constructor; the parent_type copy clones the device field.
  SolverSystemSpectralCUDA(SolverSystemSpectralCUDA const& o)
      : parent_type(o),
        transformed_len{o.transformed_len}, frame_t{nullptr}, dframe{nullptr},
        nl_work{nullptr}, p{0}, p_to_t{0}, owns{true} {
    CHECK_CUDA_ERROR(cudaMalloc(&frame_t,
                                transformed_len * sizeof(cufftDoubleComplex)));
    CHECK_CUDA_ERROR(cudaMalloc(&dframe,
                                transformed_len * sizeof(cufftDoubleComplex)));
    CHECK_CUDA_ERROR(cudaMalloc(&nl_work, len * sizeof(scalar_t)));
    CHECK_CUDA_ERROR(cudaMemcpy(frame_t, o.frame_t,
                                transformed_len * sizeof(cufftDoubleComplex),
                                cudaMemcpyDeviceToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(dframe, o.dframe,
                                transformed_len * sizeof(cufftDoubleComplex),
                                cudaMemcpyDeviceToDevice));
    make_plans();
  }

  friend void swap(SolverSystemSpectralCUDA& a, SolverSystemSpectralCUDA& b) {
    using std::swap;
    swap(static_cast<parent_type&>(a), static_cast<parent_type&>(b));
    swap(a.transformed_len, b.transformed_len);
    swap(a.frame_t, b.frame_t);
    swap(a.dframe, b.dframe);
    swap(a.nl_work, b.nl_work);
    swap(a.p, b.p);
    swap(a.p_to_t, b.p_to_t);
    swap(a.owns, b.owns);
  }

  //! One-time host->device upload of the constant spectral operators A(k),
  //! B(k) into device cufftDoubleComplex arrays. No-op after first upload.
  //! host_A/host_B are host arrays of length transformed_len in the same
  //! (dims[0]/2+1)-major layout as the field FFT. The device buffers are held
  //! in shared_ptr<void> with a cudaFree deleter so they release when the last
  //! SpectralDataSP copy is destroyed.
  void upload_operators(std::shared_ptr<void>& A_dev,
                        std::shared_ptr<void>& B_dev,
                        const complex_t* host_A, const complex_t* host_B) {
    if (A_dev) return;
    size_t bytes =
        static_cast<size_t>(transformed_len) * sizeof(cufftDoubleComplex);
    void* a = nullptr;
    void* b = nullptr;
    CHECK_CUDA_ERROR(cudaMalloc(&a, bytes));
    CHECK_CUDA_ERROR(cudaMalloc(&b, bytes));
    CHECK_CUDA_ERROR(cudaMemcpy(a, host_A, bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(b, host_B, bytes, cudaMemcpyHostToDevice));
    A_dev = std::shared_ptr<void>(a, [](void* p) { if (p) cudaFree(p); });
    B_dev = std::shared_ptr<void>(b, [](void* p) { if (p) cudaFree(p); });
  }

  //! Forward-transform the (already NL-evaluated) device field into dframe and
  //! apply the spectral linear step dframe = A*frame_t + B*dframe on device.
  //! Keeps all cuFFT calls and kernel launches inside this .cuh so the shared
  //! solver header (compiled by both nvcc and the host compiler) stays clean.
  void solve_spectral_step(void* A_dev, void* B_dev) {
    CHECK_CUFFT_ERROR(cufftExecD2Z(p_to_t, values, dframe));
    int threads = 256;
    int blocks =
        static_cast<int>((transformed_len + threads - 1) / threads);
    symphas::internal::sp_cuda_spectral_mul_kernel<<<blocks, threads>>>(
        dframe, frame_t, reinterpret_cast<cufftDoubleComplex*>(A_dev),
        reinterpret_cast<cufftDoubleComplex*>(B_dev), transformed_len);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  //! Save/restore the device field into nl_work across NL evaluation.
  void save_field() {
    CHECK_CUDA_ERROR(cudaMemcpy(nl_work, values, len * sizeof(scalar_t),
                                cudaMemcpyDeviceToDevice));
  }
  void restore_field() {
    CHECK_CUDA_ERROR(cudaMemcpy(values, nl_work, len * sizeof(scalar_t),
                                cudaMemcpyDeviceToDevice));
  }

  //! Apply cross-field spectral contributions on the GPU path. The GPU spectral
  //! step (solve_spectral_step) only applies the self-field A/B operators; for
  //! multi-field models with LINEAR cross-coupling (e.g. Model D, KKS) the
  //! serial/MPI paths additionally accumulate  dframe += sum_c K_c(k)*frame_t_j(k).
  //! Without this the GPU silently drops those terms and diverges from serial.
  //!
  //! K_host[c] are host arrays (length transformed_len) and other_frame_t_dev[c]
  //! are DEVICE pointers to the coupled field's k-space snapshot (each coupled
  //! field is itself a SolverSystemSpectralCUDA, so its frame_t lives on device).
  //! We copy this field's dframe and each coupled frame_t to the host, apply the
  //! multiply-add there, and copy dframe back. num_cross is tiny (<= a few) and
  //! this runs once per field per step, so the transfer cost is acceptable and
  //! the result is bit-consistent with the host solver's complex arithmetic.
  void apply_cross_fields_host(size_t num_cross,
                               complex_t* const* K_host,
                               void* const* other_frame_t_dev) {
    if (num_cross == 0) return;
    size_t n = static_cast<size_t>(transformed_len);
    std::vector<complex_t> df(n), ft(n);
    CHECK_CUDA_ERROR(cudaMemcpy(df.data(), dframe,
                                n * sizeof(cufftDoubleComplex),
                                cudaMemcpyDeviceToHost));
    for (size_t c = 0; c < num_cross; ++c) {
      CHECK_CUDA_ERROR(cudaMemcpy(ft.data(), other_frame_t_dev[c],
                                  n * sizeof(cufftDoubleComplex),
                                  cudaMemcpyDeviceToHost));
      const complex_t* K = K_host[c];
      for (size_t i = 0; i < n; ++i) df[i] += K[i] * ft[i];
    }
    CHECK_CUDA_ERROR(cudaMemcpy(dframe, df.data(),
                                n * sizeof(cufftDoubleComplex),
                                cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  //! Device pointer to this system's k-space snapshot (for cross-field coupling
  //! from another field's GPU spectral step).
  void* frame_t_device() { return reinterpret_cast<void*>(frame_t); }


  //! Recover the real field: inverse FFT (dframe -> values) and normalize.
  void update(iter_type index, double) {
    // Periodically rebuild the k-space accumulator from the real field to
    // prevent drift AND to make update() idempotent: cufftExecZ2D below
    // DESTROYS dframe (cuFFT c2r overwrites its input), so if update() is
    // called again before equation() regenerates dframe (e.g. an initial
    // checkpoint save's update followed by the first time-step update at the
    // same index), the second inverse transform would read garbage. Rebuilding
    // dframe = FFT(values) here from the still-valid field guards against that.
    // Matches the CPU SolverSystemSpectral::update re-sync (index % 100 == 0).
    if (index % 100 == 0) {
      CHECK_CUFFT_ERROR(cufftExecD2Z(p_to_t, values, dframe));
    }
    // Refresh k-space field snapshot from the accumulated solution.
    CHECK_CUDA_ERROR(cudaMemcpy(frame_t, dframe,
                                transformed_len * sizeof(cufftDoubleComplex),
                                cudaMemcpyDeviceToDevice));
    CHECK_CUFFT_ERROR(cufftExecZ2D(p, dframe, values));
    // Inverse FFT (cuFFT Z2D) is unnormalized; divide by the total real-grid
    // point count. `len` (= product of all dims, the real field allocation
    // length) is the correct factor for every dimension, including D == 1.
    int threads = 256;
    int blocks = static_cast<int>((len + threads - 1) / threads);
    symphas::internal::sp_cuda_scale_kernel<<<blocks, threads>>>(
        values, len, 1.0 / static_cast<double>(len));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  ~SolverSystemSpectralCUDA() {
    if (owns) {
      if (p) cufftDestroy(p);
      if (p_to_t) cufftDestroy(p_to_t);
    }
    if (frame_t) cudaFree(frame_t);
    if (dframe) cudaFree(dframe);
    if (nl_work) cudaFree(nl_work);
  }
};

//! Teach grid::value_type_of about the GPU spectral system. The generic
//! implementation casts to Block<T>/MultiBlock<N,T>, but the device system's
//! storage is GridCUDA (BlockCUDA), which is not a host Block. The order
//! parameter is real-valued, so the underlying value type is scalar_t.
namespace grid {
template <size_t D>
struct value_type_of<SolverSystemSpectralCUDA<D>> {
  using type = scalar_t;
};

//! Likewise teach grid::dimension_of (the generic version casts to Grid<T,D>,
//! which the device system is not).
template <size_t D>
struct dimension_of<SolverSystemSpectralCUDA<D>> {
  static const size_t value = D;
};
}  // namespace grid

#endif



