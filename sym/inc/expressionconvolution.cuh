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
 * PURPOSE: CUDA kernels for convolution operations in expressionconvolution.h
 *
 * ***************************************************************************
 */

#pragma once

#ifdef USING_CUDA

#include <cuda_runtime.h>
#include <cufft.h>
#include "definitions.h"
#include "convolutionlib.cuh"

namespace symphas::internal::cuda {

//! CUDA kernel for complex multiplication in convolution computation
/*!
 * Performs element-wise complex multiplication of two arrays and stores 
 * the result in a third array. This implements the convolution theorem
 * multiplication step in Fourier space.
 *
 * \param out_0 First complex array (Fourier transform of first expression)
 * \param out_1 Second complex array (Fourier transform of second expression) 
 * \param in_2 Output array for the complex multiplication result
 * \param len Number of elements to process
 */
__global__ void complex_multiply_kernel(cufftDoubleComplex* out_0,
                                      cufftDoubleComplex* out_1, 
                                      cufftDoubleComplex* in_2,
                                      len_type len) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx < len) {
    // Complex multiplication: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
    double real_0 = out_0[idx].x;
    double imag_0 = out_0[idx].y;
    double real_1 = out_1[idx].x;
    double imag_1 = out_1[idx].y;
    
    in_2[idx].x = real_0 * real_1 - imag_0 * imag_1;
    in_2[idx].y = real_0 * imag_1 + imag_0 * real_1;
  }
}

//! CUDA kernel for complex scaling in convolution computation
/*!
 * Scales a complex array by a real scalar value.
 *
 * \param data Complex array to scale
 * \param len Number of elements to process
 * \param scale Scaling factor
 */
__global__ void complex_scale_kernel(cufftDoubleComplex* data,
                                   len_type len,
                                   double scale) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx < len) {
    data[idx].x *= scale;
    data[idx].y *= scale;
  }
}

//! CUDA kernel for Gaussian smoothing convolution
/*!
 * Applies Gaussian smoothing in Fourier space by multiplying the 
 * Fourier-transformed data with the Gaussian kernel values.
 *
 * \param out_0 Fourier-transformed input data
 * \param in_1 Output array for smoothed data
 * \param smoother_values Gaussian kernel values
 * \param len Number of elements to process
 */
__global__ void gaussian_smooth_kernel(cufftDoubleComplex* out_0,
                                     cufftDoubleComplex* in_1,
                                     double* smoother_values,
                                     len_type len) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx < len) {
    double smooth_val = smoother_values[idx];
    in_1[idx].x = smooth_val * out_0[idx].x;
    in_1[idx].y = smooth_val * out_0[idx].y;
  }
}

//! Host function to launch complex multiplication kernel
inline void launch_complex_multiply(cufftDoubleComplex* out_0,
                                   cufftDoubleComplex* out_1,
                                   cufftDoubleComplex* in_2,
                                   len_type len) {
  constexpr int blockSize = 256;
  int numBlocks = (len + blockSize - 1) / blockSize;
  
  complex_multiply_kernel<<<numBlocks, blockSize>>>(out_0, out_1, in_2, len);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

//! Host function to launch complex scaling kernel
inline void launch_complex_scale(cufftDoubleComplex* data,
                                len_type len,
                                double scale) {
  constexpr int blockSize = 256;
  int numBlocks = (len + blockSize - 1) / blockSize;
  
  complex_scale_kernel<<<numBlocks, blockSize>>>(data, len, scale);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

//! Host function to launch Gaussian smoothing kernel  
inline void launch_gaussian_smooth(cufftDoubleComplex* out_0,
                                  cufftDoubleComplex* in_1,
                                  double* smoother_values,
                                  len_type len) {
  constexpr int blockSize = 256;
  int numBlocks = (len + blockSize - 1) / blockSize;
  
  gaussian_smooth_kernel<<<numBlocks, blockSize>>>(out_0, in_1, smoother_values, len);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

} // namespace symphas::internal::cuda

// Forward declarations for OpConvolution classes
template <typename V, typename E1, typename E2>
struct OpConvolution;

template <typename V, size_t D, template <typename, size_t> typename grid_type, typename E>
struct OpConvolution<V, GaussianSmoothing<D, grid_type>, E>;

template <typename V, size_t D, template <typename, size_t> typename grid_type, typename G>
struct OpConvolution<V, GaussianSmoothing<D, grid_type>, OpTerm<OpIdentity, G>>;

namespace symphas::internal {

//! CUDA implementation of convolution update for general two-expression convolution
template <typename V, typename E1, typename E2, size_t D0, typename G_T>
void update_impl_cuda(expr::ConvolutionDataPairCUDA<D0>& compute_data, 
                      Grid<G_T, D0>& result_grid) {
  len_type len = symphas::dft::length<G_T, D0>(result_grid.dims);
  if constexpr (std::is_same<G_T, complex_t>::value) {
    symphas::internal::cuda::launch_complex_scale(compute_data.out_0, len, 1.0 / len);
    symphas::internal::cuda::launch_complex_scale(compute_data.out_1, len, 1.0 / len);
  }

  // Perform complex multiplication in CUDA
  symphas::internal::cuda::launch_complex_multiply(compute_data.out_0, compute_data.out_1, 
                                                  compute_data.in_2, len);

  compute_data.transform_out_in(result_grid.values, result_grid.dims);

  if constexpr (std::is_same<G_T, scalar_t>::value) {
    grid::scale(result_grid);
  }
}

//! CUDA implementation of Gaussian smoothing convolution for expressions
template <typename V, size_t D, template <typename, size_t> typename grid_type, 
          typename E, size_t D0, typename G_T>
void update_impl_cuda(expr::ConvolutionDataCUDA<D0>& compute_data, 
                      grid_type<G_T, D>& result_grid,
                      const GaussianSmoothing<D, grid_type>& smoother) {
  len_type len = symphas::dft::length<G_T, D>(result_grid.dims);
  
  // Copy smoother values to device if needed (this assumes smoother.data.values is on device)
  symphas::internal::cuda::launch_gaussian_smooth(compute_data.out_0, compute_data.in_1,
                                                  smoother.data.values, len);
  
  compute_data.transform_out_in(result_grid.values, result_grid.dims);
  grid::scale(result_grid);
}

} // namespace symphas::internal

#endif // USING_CUDA
