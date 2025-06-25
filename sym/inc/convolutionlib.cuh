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
 * PURPOSE: Implements CUDA version of the Gaussian smoothing kernel and 
 * elements used in computing the convolution on GPU using cuFFT.
 *
 * ***************************************************************************
 */

#pragma once

#include "convolutionlib.h"
#include "gridfunctions.cuh"
#include "dft.cuh"
#include <cufft.h>
#include <cuda_runtime.h>
#include <iostream>


namespace symphas::internal {

template <size_t D, size_t D0 = 0>
__global__ void complex_gaussian_kernel_fill_cuda(double* kernel, len_type const* dims,
                                                  len_type const* stride, iter_type* index,
                                                  const double* h, double sigma) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int total_elements = 1;
  for (int i = 0; i < D; ++i) {
    total_elements *= dims[i];
  }
  
  if (idx < total_elements) {
    // Convert linear index to multi-dimensional index
    iter_type temp_idx = idx;
    iter_type pos[D];
    for (int i = D - 1; i >= 0; --i) {
      pos[i] = temp_idx % dims[i];
      temp_idx /= dims[i];
    }

    double k[D]{};
    for (size_t i = 0; i < D; ++i) {
      double factor = 2.0 * symphas::PI * h[i] / dims[i];
      k[i] = ((pos[i] < dims[i] / 2) ? pos[i] : pos[i] - dims[i]) * factor;
    }
    double kk = 0.0;
    for (size_t i = 0; i < D; ++i) {
      kk += k[i] * k[i];
    }
    kernel[idx] = exp(-kk / (2.0 * sigma * sigma));
  }
}

template <size_t D, size_t D0 = 0>
__global__ void real_gaussian_kernel_fill_cuda(double* kernel, len_type const* dims,
                                               len_type const* stride, iter_type* index,
                                               const double* h, double sigma) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int total_elements = 1;
  for (int i = 0; i < D; ++i) {
    total_elements *= dims[i];
  }
  
  if (idx < total_elements) {
    // Convert linear index to multi-dimensional index
    iter_type temp_idx = idx;
    iter_type pos[D];
    for (int i = D - 1; i >= 0; --i) {
      pos[i] = temp_idx % dims[i];
      temp_idx /= dims[i];
    }

    double r[D]{};
    for (size_t i = 0; i < D; ++i) {
      r[i] = pos[i] - dims[i] / 2;
    }
    double rr = 0.0;
    for (size_t i = 0; i < D; ++i) {
      rr += r[i] * r[i];
    }
    double value = (1. / (2.0 * symphas::PI * sigma * sigma)) *
                   exp(-rr / (2.0 * sigma * sigma));
    kernel[idx] = value;
  }
}

template <size_t D>
void complex_gaussian_kernel_cuda(double* kernel, len_type const* dims,
                                 const double* h, double sigma) {
  // Allocate device memory for dimensions and h
  len_type* d_dims;
  double* d_h;
  len_type* d_stride;
  iter_type* d_index;
  
  CHECK_CUDA_ERROR(cudaMalloc(&d_dims, D * sizeof(len_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_h, D * sizeof(double)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_stride, D * sizeof(len_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_index, D * sizeof(iter_type)));
  
  CHECK_CUDA_ERROR(cudaMemcpy(d_dims, dims, D * sizeof(len_type), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(d_h, h, D * sizeof(double), cudaMemcpyHostToDevice));
  
  len_type total_elements = 1;
  for (size_t i = 0; i < D; ++i) {
    total_elements *= dims[i];
  }
  
  int blockSize = BLOCK_SIZE;
  int numBlocks = (total_elements + blockSize - 1) / blockSize;
  
  complex_gaussian_kernel_fill_cuda<D><<<numBlocks, blockSize>>>(
      kernel, d_dims, d_stride, d_index, d_h, sigma);
  
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(d_dims));
  CHECK_CUDA_ERROR(cudaFree(d_h));
  CHECK_CUDA_ERROR(cudaFree(d_stride));
  CHECK_CUDA_ERROR(cudaFree(d_index));
}

template <size_t D>
void real_gaussian_kernel_cuda(double* kernel, len_type const* dims,
                              const double* h, double sigma) {
  // Allocate device memory for dimensions and h
  len_type* d_dims;
  double* d_h;
  len_type* d_stride;
  iter_type* d_index;
  
  CHECK_CUDA_ERROR(cudaMalloc(&d_dims, D * sizeof(len_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_h, D * sizeof(double)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_stride, D * sizeof(len_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_index, D * sizeof(iter_type)));
  
  CHECK_CUDA_ERROR(cudaMemcpy(d_dims, dims, D * sizeof(len_type), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(d_h, h, D * sizeof(double), cudaMemcpyHostToDevice));
  
  len_type total_elements = 1;
  for (size_t i = 0; i < D; ++i) {
    total_elements *= dims[i];
  }
  
  int blockSize = BLOCK_SIZE;
  int numBlocks = (total_elements + blockSize - 1) / blockSize;
  
  real_gaussian_kernel_fill_cuda<D><<<numBlocks, blockSize>>>(
      kernel, d_dims, d_stride, d_index, d_h, sigma);
  
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(d_dims));
  CHECK_CUDA_ERROR(cudaFree(d_h));
  CHECK_CUDA_ERROR(cudaFree(d_stride));
  CHECK_CUDA_ERROR(cudaFree(d_index));
}

// Overloaded functions for CUDA Grid types
template <typename T, size_t D>
void complex_gaussian_kernel(GridCUDA<T, D>& grid,
                             grid::region_interval<D> const& domain,
                           const double* h, double sigma) {
  complex_gaussian_kernel_cuda<D>(grid.values, grid.dims, h, sigma);
}

template <typename T, size_t D>
void real_gaussian_kernel(GridCUDA<T, D>& grid,
                          grid::region_interval<D> const& domain,
                         const double* h, double sigma) {
  real_gaussian_kernel_cuda<D>(grid.values, grid.dims, h, sigma);
}

} // namespace symphas::internal

//! CUDA version of Gaussian smoothing kernel.
/*!
 * CUDA expression encapsulating a grid that contains values of a Gaussian smoothing
 * kernel stored on GPU device memory. The values are precomputed upon construction
 * of the object using CUDA kernels.
 */
template <size_t D>
using GaussianSmoothingCUDA = GaussianSmoothing<D, GridCUDA>;

// ******************************************************************************************************************

namespace expr {

//! CUDA delegation method for computing convolution of different types.
/*!
 * A CUDA computation object for the convolution using cuFFT routines 
 * to compute transforms on GPU device memory.
 */
template <size_t D>
struct ConvolutionDataCUDA {
  ConvolutionDataCUDA() : out_0{nullptr}, in_1{nullptr}, plan_forward{0}, plan_backward{0} {}

  template <typename T, typename S>
  ConvolutionDataCUDA(T* in_0, S* out_1, len_type* dims, len_type len);
  template <typename T, typename S>
  ConvolutionDataCUDA(GridCUDA<T, D> const& in_0, S* out_1)
      : ConvolutionDataCUDA(in_0.values, out_1, in_0.dims, in_0.len) {}
  template <typename T, typename S>
  ConvolutionDataCUDA(T* in_0, GridCUDA<S, D>& out_1)
      : ConvolutionDataCUDA(in_0, out_1.values, out_1.dims, out_1.len) {}
  template <typename T, typename S>
  ConvolutionDataCUDA(GridCUDA<T, D> const& in_0, GridCUDA<S, D>& out_1)
      : ConvolutionDataCUDA(in_0.values, out_1.values, out_1.dims, out_1.len) {}

  ConvolutionDataCUDA(ConvolutionDataCUDA<D> const& other) = delete;
  ConvolutionDataCUDA(ConvolutionDataCUDA<D>&& other) { swap(*this, other); }
  ConvolutionDataCUDA<D>& operator=(ConvolutionDataCUDA<D> const&) = delete;

  friend void swap(ConvolutionDataCUDA<D>& first, ConvolutionDataCUDA<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.in_1, second.in_1);
    swap(first.plan_forward, second.plan_forward);
    swap(first.plan_backward, second.plan_backward);
  }

  void transform_in_out(...) const { 
    CHECK_CUFFT_ERROR(cufftExecZ2Z(plan_forward, out_0, out_0, CUFFT_FORWARD));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  void transform_out_in(...) const { 
    CHECK_CUFFT_ERROR(cufftExecZ2Z(plan_backward, in_1, in_1, CUFFT_INVERSE));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  cufftDoubleComplex* out_0;  //!< Fourier transform of the given data.
  cufftDoubleComplex* in_1;   //!< Input of Fourier data.
  cufftHandle plan_forward;   //!< Forward cuFFT plan.
  cufftHandle plan_backward;  //!< Backward cuFFT plan.

  ~ConvolutionDataCUDA() {
    if (plan_forward) cufftDestroy(plan_forward);
    if (plan_backward) cufftDestroy(plan_backward);
    if (out_0) CHECK_CUDA_ERROR(cudaFree(out_0));
    if (in_1) CHECK_CUDA_ERROR(cudaFree(in_1));
  }
};

// Specializations for different dimensionalities
template <>
template <>
inline ConvolutionDataCUDA<1>::ConvolutionDataCUDA(scalar_t* in_0, scalar_t* out_1,
                                                   len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_1, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan1d(&plan_forward, dims[0], CUFFT_Z2Z, 1));
  CHECK_CUFFT_ERROR(cufftPlan1d(&plan_backward, dims[0], CUFFT_Z2Z, 1));
}

template <>
template <>
inline ConvolutionDataCUDA<2>::ConvolutionDataCUDA(scalar_t* in_0, scalar_t* out_1,
                                                   len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_1, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan2d(&plan_forward, dims[0], dims[1], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan2d(&plan_backward, dims[0], dims[1], CUFFT_Z2Z));
}

template <>
template <>
inline ConvolutionDataCUDA<3>::ConvolutionDataCUDA(scalar_t* in_0, scalar_t* out_1,
                                                   len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_1, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan3d(&plan_forward, dims[0], dims[1], dims[2], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan3d(&plan_backward, dims[0], dims[1], dims[2], CUFFT_Z2Z));
}

//! CUDA delegation method for computing convolution with two input fields.
/*!
 * A CUDA computation object for the convolution using cuFFT routines when there
 * are two input fields and one output field.
 */
template <size_t D>
struct ConvolutionDataPairCUDA {
  ConvolutionDataPairCUDA()
      : out_0{nullptr},
        out_1{nullptr},
        in_2{nullptr},
        plan_forward_0{0},
        plan_forward_1{0},
        plan_backward{0} {}

  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPairCUDA(T_0* in_0, T_1* in_1, R* out_2, len_type* dims,
                         len_type len);
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPairCUDA(GridCUDA<T_0, D> const& in_0, GridCUDA<T_1, D> const& in_1,
                         R* out_2)
      : ConvolutionDataPairCUDA(in_0.values, in_1.values, out_2, in_0.dims,
                               in_0.len) {}
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPairCUDA(T_0* in_0, T_1* in_1, GridCUDA<R, D>& out_2)
      : ConvolutionDataPairCUDA(in_0, in_1, out_2.values, out_2.dims, out_2.len) {}

  ConvolutionDataPairCUDA(ConvolutionDataPairCUDA<D> const&) = delete;
  ConvolutionDataPairCUDA(ConvolutionDataPairCUDA<D>&& other) : ConvolutionDataPairCUDA() {
    swap(*this, other);
  }
  ConvolutionDataPairCUDA<D>& operator=(ConvolutionDataPairCUDA<D> const&) = delete;

  friend void swap(ConvolutionDataPairCUDA<D>& first,
                   ConvolutionDataPairCUDA<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.out_1, second.out_1);
    swap(first.in_2, second.in_2);
    swap(first.plan_forward_0, second.plan_forward_0);
    swap(first.plan_forward_1, second.plan_forward_1);
    swap(first.plan_backward, second.plan_backward);
  }

  void transform_in_out_0(...) const { 
    CHECK_CUFFT_ERROR(cufftExecZ2Z(plan_forward_0, out_0, out_0, CUFFT_FORWARD));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  void transform_in_out_1(...) const { 
    CHECK_CUFFT_ERROR(cufftExecZ2Z(plan_forward_1, out_1, out_1, CUFFT_FORWARD));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  void transform_out_in(...) const { 
    CHECK_CUFFT_ERROR(cufftExecZ2Z(plan_backward, in_2, in_2, CUFFT_INVERSE));
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  ~ConvolutionDataPairCUDA() {
    if (plan_backward) cufftDestroy(plan_backward);
    if (plan_forward_0) cufftDestroy(plan_forward_0);
    if (plan_forward_1) cufftDestroy(plan_forward_1);
    if (out_0) CHECK_CUDA_ERROR(cudaFree(out_0));
    if (out_1) CHECK_CUDA_ERROR(cudaFree(out_1));
    if (in_2) CHECK_CUDA_ERROR(cudaFree(in_2));
  }

  cufftDoubleComplex* out_0;   //!< Fourier transform of the first given data.
  cufftDoubleComplex* out_1;   //!< Fourier transform of the second given data.
  cufftDoubleComplex* in_2;    //!< Input of Fourier data.
  cufftHandle plan_forward_0;  //!< Forward cuFFT plan for first input.
  cufftHandle plan_forward_1;  //!< Forward cuFFT plan for second input.
  cufftHandle plan_backward;   //!< Backward cuFFT plan for output.
};

template <>
template <>
inline ConvolutionDataPairCUDA<1>::ConvolutionDataPairCUDA(scalar_t* in_0,
                                                           scalar_t* in_1,
                                                           scalar_t* out_2,
                                                           len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&out_1, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_2, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan1d(&plan_forward_0, dims[0], CUFFT_Z2Z, 1));
  CHECK_CUFFT_ERROR(cufftPlan1d(&plan_forward_1, dims[0], CUFFT_Z2Z, 1));
  CHECK_CUFFT_ERROR(cufftPlan1d(&plan_backward, dims[0], CUFFT_Z2Z, 1));
}

template <>
template <>
inline ConvolutionDataPairCUDA<2>::ConvolutionDataPairCUDA(scalar_t* in_0,
                                                           scalar_t* in_1,
                                                           scalar_t* out_2,
                                                           len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&out_1, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_2, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan2d(&plan_forward_0, dims[0], dims[1], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan2d(&plan_forward_1, dims[0], dims[1], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan2d(&plan_backward, dims[0], dims[1], CUFFT_Z2Z));
}

template <>
template <>
inline ConvolutionDataPairCUDA<3>::ConvolutionDataPairCUDA(scalar_t* in_0,
                                                           scalar_t* in_1,
                                                           scalar_t* out_2,
                                                           len_type* dims, len_type len) {
  // Allocate device memory for complex arrays
  CHECK_CUDA_ERROR(cudaMalloc(&out_0, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&out_1, len * sizeof(cufftDoubleComplex)));
  CHECK_CUDA_ERROR(cudaMalloc(&in_2, len * sizeof(cufftDoubleComplex)));
  
  // Create cuFFT plans
  CHECK_CUFFT_ERROR(cufftPlan3d(&plan_forward_0, dims[0], dims[1], dims[2], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan3d(&plan_forward_1, dims[0], dims[1], dims[2], CUFFT_Z2Z));
  CHECK_CUFFT_ERROR(cufftPlan3d(&plan_backward, dims[0], dims[1], dims[2], CUFFT_Z2Z));
}

} // namespace expr



template <typename T1, size_t D, typename T2>
auto expr::make_convolution_data(const GridCUDA<T1, D>& in_0, T2* out_1) {
  return ConvolutionDataCUDA<D>(in_0.values, out_1, in_0.dims, in_0.len);
}

template <typename T1, typename T2, size_t D, typename T3>
auto expr::make_convolution_data(const GridCUDA<T1, D>& in_0,
                           const GridCUDA<T2, D>& in_1,
                           T3* out_1) {
  return ConvolutionDataPairCUDA<D>(in_0, in_1, out_1);
}

