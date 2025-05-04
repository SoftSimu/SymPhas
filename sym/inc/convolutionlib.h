
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
 * PURPOSE: Implements the Gaussian smoothing kernel and elements used in
 * computing the convolution.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressiontypek.h"
#include "gridfunctions.h"
#include "spslibfftw.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

#ifdef LATEX_PLOT
#define SYEX_GAUSSIAN_KERNEL "G"
#else
#define SYEX_GAUSSIAN_KERNEL "G"
#endif

namespace symphas::internal {

template <size_t D, size_t D0 = 0>
void complex_gaussian_kernel_fill(double* kernel, len_type const (&dims)[D],
                                  len_type (&stride)[D], iter_type (&index)[D],
                                  const double* h, double sigma) {
  constexpr size_t Dn = D - D0 - 1;
  if constexpr (D0 == D) {
    double k[D]{};
    for (size_t i = 0; i < D; ++i) {
      double factor = 2.0 * symphas::PI * h[i] / dims[i];
      k[i] =
          ((index[i] < dims[i] / 2) ? index[i] : index[i] - dims[i]) * factor;
    }
    double kk = std::inner_product(k, k + D, k, 0.0);
    iter_type n = std::inner_product(index, index + D, stride, 0.0);
    kernel[n] = std::exp(-kk / (2.0 * sigma * sigma));
  } else {
    for (int i = 0; i < dims[D0]; ++i) {
      index[Dn] = i;
      complex_gaussian_kernel_fill<D, D0 + 1>(kernel, dims, stride, index, h,
                                              sigma);
    }
  }
}

template <size_t D>
void complex_gaussian_kernel(double* kernel, len_type const (&dims)[D],
                             const double* h, double sigma) {
  len_type stride[D]{};
  iter_type index[D]{};
  grid::get_stride(stride, dims);
  complex_gaussian_kernel_fill<D>(kernel, dims, stride, index, h, sigma);
}

template <size_t D, size_t D0 = 0>
void real_gaussian_kernel_fill(double* kernel, len_type const (&dims)[D],
                               len_type (&stride)[D], iter_type (&index)[D],
                               const double* h, double sigma) {
  constexpr size_t Dn = D - D0 - 1;
  if constexpr (D0 == D) {
    double r[D]{};
    for (size_t i = 0; i < D; ++i) {
      r[i] = index[i] - dims[i] / 2;
    }
    double rr = std::inner_product(r, r + D, r, 0.0);
    iter_type n = std::inner_product(index, index + D, stride, 0.0);
    double value = (1. / (2.0 * symphas::PI * sigma * sigma)) *
                   std::exp(-rr / (2.0 * sigma * sigma));
    kernel[n] = value;
  } else {
    for (int i = 0; i < dims[Dn]; ++i) {
      index[Dn] = i;
      real_gaussian_kernel_fill<D, D0 + 1>(kernel, dims, stride, index, h,
                                           sigma);
    }
  }
}

template <size_t D>
void real_gaussian_kernel(double* kernel, len_type const (&dims)[D],
                          const double* h, double sigma) {
  len_type stride[D]{};
  iter_type index[D]{};
  grid::get_stride(stride, dims);
  real_gaussian_kernel_fill<D>(kernel, dims, stride, index, h, sigma);

  // double sum = std::accumulate(kernel, kernel + grid::length<D>(dims), 0.0);
  // std::transform(kernel, kernel + grid::length<D>(dims), kernel,
  //                [sum](double val) { return val / sum; });
}
}  // namespace symphas::internal

//! Gaussian smoothing kernel.
/*!
 * Expression encapsulating a grid that contains values of a Gaussian smoothing
 * kernel. The values are precomputed upon construction of the object using the
 * class k_field, which computes the wavenumber for each point in the grid
 * according to its position.
 */
template <size_t D>
struct GaussianSmoothing : OpExpression<GaussianSmoothing<D>> {
  void allocate() {
    if (data.len == 0) {
      data = Grid<scalar_t, D>(dims);
      initialize_values(h, sigma);
    }
  }

  GaussianSmoothing() : data{0} {}

  GaussianSmoothing(symphas::interval_data_type const& intervals,
                    double sigma = 1.0, bool fourier_space = false)
      : data{0}, fourier_space{fourier_space}, sigma{sigma}, h{}, dims{} {
    for (iter_type i = 0; i < D; ++i) {
      h[i] = intervals.at(symphas::index_to_axis(i)).width();
      dims[i] = intervals.at(symphas::index_to_axis(i)).get_interval_count();
    }
  }

  GaussianSmoothing(const len_type* dims, const double* h, double sigma = 1.0,
                    bool fourier_space = false)
      : data{0}, fourier_space{fourier_space}, sigma{sigma}, h{}, dims{} {
    std::copy(h, h + D, this->h);
    std::copy(dims, dims + D, this->dims);
  }

  inline double eval(iter_type n) const { return data.values[n]; }

#ifdef PRINTABLE_EQUATIONS
  /* print statements
   */

  size_t print(FILE* out) const {
    if (fourier_space) {
      return fprintf(
          out, SYEX_FT_OF_OP_FMT_A SYEX_GAUSSIAN_KERNEL SYEX_FT_OF_OP_FMT_B);
    } else {
      return fprintf(out, SYEX_GAUSSIAN_KERNEL);
    }
  }

  size_t print(char* out) const {
    if (fourier_space) {
      return sprintf(
          out, SYEX_FT_OF_OP_FMT_A SYEX_GAUSSIAN_KERNEL SYEX_FT_OF_OP_FMT_B);
    } else {
      return sprintf(out, SYEX_GAUSSIAN_KERNEL);
    }
  }

  size_t print_length() const {
    if (fourier_space) {
      return STR_ARR_LEN(
          SYEX_FT_OF_OP_FMT_A SYEX_GAUSSIAN_KERNEL SYEX_FT_OF_OP_FMT_B);
    } else {
      return STR_ARR_LEN(SYEX_GAUSSIAN_KERNEL);
    }
  }
#endif

  auto operator*(GaussianSmoothing<D> const& other) const {
    GaussianSmoothing<D> g2(*this);
    std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par,
#endif
        data.values, data.values + data.len, g2.data.values, g2.data.values,
        [](auto a, auto b) { return a * b; });
    return g2;
  }

  auto operator-() const {
    GaussianSmoothing g{*this};
    std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par,
#endif
        g.data.values, g.data.values + g.data.len, g.data.values,
        [](auto v) { return -v; });
    return g;
  }

  template <typename coeff_t, typename T, size_t U1,
            std::enable_if_t<expr::is_coeff<coeff_t>, int>>
  friend auto operator*(coeff_t const& a, GaussianSmoothing<U1> const& b);

  Grid<scalar_t, D> data;
  bool fourier_space;
  double sigma;
  double h[D];
  len_type dims[D];

 protected:
  void scale(scalar_t value) {
    for (iter_type i = 0; i < data.len; ++i) {
      data[i] *= value;
    }
  }

  void initialize_values(const double* h, double sigma) {
    if (fourier_space) {
      symphas::internal::complex_gaussian_kernel<D>(data, data.dims, h, sigma);
    } else {
      symphas::internal::real_gaussian_kernel<D>(data, data.dims, h, sigma);
    }
  }
};

template <typename coeff_t, typename T, size_t U1,
          std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& a, GaussianSmoothing<U1> const& b) {
  GaussianSmoothing<U1> scaled(b);
  scaled.scale(expr::eval(a));
  return scaled;
}

// ******************************************************************************************************************

#ifdef USING_FFTW

namespace expr {

//! Delegation method for computing convolution of different types.
/*!
 * A computation object for the convolution
 * so that the FFTW routines can be used on the data to compute
 * from real space to Fourier space and back.
 */
template <size_t D>
struct ConvolutionData {
  ConvolutionData() : out_0{0}, in_1{0}, p_in_out{0}, p_out_in{0} {}

  template <typename T, typename S>
  ConvolutionData(T* in_0, S* out_1, len_type* dims, len_type len);
  template <typename T, typename S>
  ConvolutionData(Grid<T, D> const& in_0, S* out_1)
      : ConvolutionData(in_0.values, out_1, in_0.dims, in_0.len) {}
  template <typename T, typename S>
  ConvolutionData(T* in_0, Grid<S, D>& out_1)
      : ConvolutionData(in_0, out_1.values, out_1.dims, out_1.len) {}
  template <typename T, typename S>
  ConvolutionData(Grid<T, D> const& in_0, Grid<S, D>& out_1)
      : ConvolutionData(in_0.values, out_1.values, out_1.dims, out_1.len) {}

  ConvolutionData(ConvolutionData<D> const& other) = delete;
  ConvolutionData(ConvolutionData<D>&& other) { swap(*this, other); }
  ConvolutionData<D>& operator=(ConvolutionData<D> const&) = delete;

  friend void swap(ConvolutionData<D>& first, ConvolutionData<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.in_1, second.in_1);
    swap(first.p_in_out, second.p_in_out);
    swap(first.p_out_in, second.p_out_in);
  }

  void transform_in_out(...) const { symphas::dft::fftw_execute(p_in_out); }

  void transform_out_in(...) const { symphas::dft::fftw_execute(p_out_in); }

  fftw_complex* out_0;  //!< Fourier transform of the given data.
  fftw_complex* in_1;   //!< Input of Fourier data.
  fftw_plan p_in_out;   //!< Transforms the given input data.
  fftw_plan
      p_out_in;  //!< Transforms member Fourier data into given output data.

  ~ConvolutionData() {
    symphas::dft::fftw_destroy_plan(p_out_in);
    symphas::dft::fftw_destroy_plan(p_in_out);
    symphas::dft::fftw_free(out_0);
    symphas::dft::fftw_free(in_1);
  }
};

using symphas::dft::new_fftw_plan;

template <>
template <>
inline ConvolutionData<1>::ConvolutionData(scalar_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 1>(dims))},
      in_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 1>(dims))},
      p_in_out{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<1, complex_t, scalar_t>{}(in_1, out_1, dims)} {}

template <>
template <>
inline ConvolutionData<2>::ConvolutionData(scalar_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 2>(dims))},
      in_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 2>(dims))},
      p_in_out{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<2, complex_t, scalar_t>{}(in_1, out_1, dims)} {}

template <>
template <>
inline ConvolutionData<3>::ConvolutionData(scalar_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 3>(dims))},
      in_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 3>(dims))},
      p_in_out{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<3, complex_t, scalar_t>{}(in_1, out_1, dims)} {}

template <>
template <>
inline ConvolutionData<1>::ConvolutionData(complex_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<1, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{new_fftw_plan<1, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<2>::ConvolutionData(complex_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<2, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{new_fftw_plan<2, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<3>::ConvolutionData(complex_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<3, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{new_fftw_plan<3, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<1>::ConvolutionData(scalar_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<1, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<2>::ConvolutionData(scalar_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<2, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<3>::ConvolutionData(scalar_t* in_0, complex_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_out_in{new_fftw_plan<3, complex_t, complex_t>{}(in_1, out_1, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionData<1>::ConvolutionData(complex_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<1, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{
          new_fftw_plan<1, complex_t, scalar_t>{}(in_1, out_1, dims, false)} {}

template <>
template <>
inline ConvolutionData<2>::ConvolutionData(complex_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<2, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{
          new_fftw_plan<2, complex_t, scalar_t>{}(in_1, out_1, dims, false)} {}

template <>
template <>
inline ConvolutionData<3>::ConvolutionData(complex_t* in_0, scalar_t* out_1,
                                           len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      in_1{symphas::dft::fftw_alloc_complex(len)},
      p_in_out{new_fftw_plan<3, complex_t, complex_t>{}(in_0, out_0, dims, true,
                                                        false)},
      p_out_in{
          new_fftw_plan<3, complex_t, scalar_t>{}(in_1, out_1, dims, false)} {}

//! Delegation method for computing convolution of different types.
/*!
 * A computation object for the convolution
 * so that the FFTW routines can be used on the data to compute
 * from real space to Fourier space and back. Performs convolution when there
 * are two input fields and one output field.
 */
template <size_t D>
struct ConvolutionDataPair {
  ConvolutionDataPair()
      : out_0{0},
        out_1{0},
        in_2{0},
        p_in_out_0{0},
        p_in_out_1{0},
        p_out_in{0} {}

  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(T_0* in_0, T_1* in_1, R* out_2, len_type* dims,
                      len_type len);
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(Grid<T_0, D> const& in_0, Grid<T_1, D> const& in_1,
                      R* out_2)
      : ConvolutionDataPair(in_0.values, in_1.values, out_2, in_0.dims,
                            in_0.len) {}
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(T_0* in_0, T_1* in_1, Grid<R, D>& out_2)
      : ConvolutionDataPair(in_0, in_1, out_2.values, out_2.dims, out_2.len) {}

  ConvolutionDataPair(ConvolutionDataPair<D> const&) = delete;
  ConvolutionDataPair(ConvolutionDataPair<D>&& other) : ConvolutionDataPair() {
    swap(*this, other);
  }
  ConvolutionDataPair<D>& operator=(ConvolutionDataPair<D> const&) = delete;

  friend void swap(ConvolutionDataPair<D>& first,
                   ConvolutionDataPair<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.out_1, second.out_1);
    swap(first.in_2, second.in_2);
    swap(first.p_in_out_0, second.p_in_out_0);
    swap(first.p_in_out_1, second.p_in_out_1);
    swap(first.p_out_in, second.p_out_in);
  }

  void transform_in_out_0(...) const { symphas::dft::fftw_execute(p_in_out_0); }

  void transform_in_out_1(...) const { symphas::dft::fftw_execute(p_in_out_1); }

  void transform_out_in(...) const { symphas::dft::fftw_execute(p_out_in); }

  ~ConvolutionDataPair() {
    symphas::dft::fftw_destroy_plan(p_out_in);
    symphas::dft::fftw_destroy_plan(p_in_out_0);
    symphas::dft::fftw_destroy_plan(p_in_out_1);
    symphas::dft::fftw_free(out_0);
    symphas::dft::fftw_free(out_1);
    symphas::dft::fftw_free(in_2);
  }

  /* for finding the convolution from real space, back into real space, we need
   * to first compute it in fourier space
   */
  fftw_complex* out_0;   //!< Fourier transform of the first given data.
  fftw_complex* out_1;   //!< Fourier transform of the second given data.
  fftw_complex* in_2;    //!< Input of Fourier data.
  fftw_plan p_in_out_0;  //!< Transforms the first of the given input data.
  fftw_plan p_in_out_1;  //!< Transforms the second of the given input data.
  fftw_plan
      p_out_in;  //!< Transforms member Fourier data into given output data.
};

template <>
template <>
inline ConvolutionDataPair<1>::ConvolutionDataPair(scalar_t* in_0,
                                                   scalar_t* in_1,
                                                   scalar_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 1>(dims))},
      out_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 1>(dims))},
      in_2{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 1>(dims))},
      p_in_out_0{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{
          new_fftw_plan<1, complex_t, scalar_t>{}(in_2, out_2, dims, false)} {}

template <>
template <>
inline ConvolutionDataPair<2>::ConvolutionDataPair(scalar_t* in_0,
                                                   scalar_t* in_1,
                                                   scalar_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 2>(dims))},
      out_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 2>(dims))},
      in_2{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 2>(dims))},
      p_in_out_0{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{
          new_fftw_plan<2, complex_t, scalar_t>{}(in_2, out_2, dims, false)} {}

template <>
template <>
inline ConvolutionDataPair<3>::ConvolutionDataPair(scalar_t* in_0,
                                                   scalar_t* in_1,
                                                   scalar_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 3>(dims))},
      out_1{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 3>(dims))},
      in_2{symphas::dft::fftw_alloc_complex(
          symphas::dft::length<scalar_t, 3>(dims))},
      p_in_out_0{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{
          new_fftw_plan<3, complex_t, scalar_t>{}(in_2, out_2, dims, false)} {}

template <>
template <>
inline ConvolutionDataPair<1>::ConvolutionDataPair(scalar_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{new_fftw_plan<1, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<1, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<2>::ConvolutionDataPair(scalar_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{new_fftw_plan<2, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<2, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<3>::ConvolutionDataPair(scalar_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_0, out_0, dims, true)},
      p_in_out_1{new_fftw_plan<3, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<3, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<1>::ConvolutionDataPair(complex_t* in_0,
                                                   scalar_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<1, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{
          new_fftw_plan<1, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{new_fftw_plan<1, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<2>::ConvolutionDataPair(complex_t* in_0,
                                                   scalar_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<2, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{
          new_fftw_plan<2, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{new_fftw_plan<2, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<3>::ConvolutionDataPair(complex_t* in_0,
                                                   scalar_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<3, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{
          new_fftw_plan<3, scalar_t, complex_t>{}(in_1, out_1, dims, true)},
      p_out_in{new_fftw_plan<3, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<1>::ConvolutionDataPair(complex_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<1, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{new_fftw_plan<1, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<1, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<2>::ConvolutionDataPair(complex_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<2, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{new_fftw_plan<2, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<2, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

template <>
template <>
inline ConvolutionDataPair<3>::ConvolutionDataPair(complex_t* in_0,
                                                   complex_t* in_1,
                                                   complex_t* out_2,
                                                   len_type* dims, len_type len)
    : out_0{symphas::dft::fftw_alloc_complex(len)},
      out_1{symphas::dft::fftw_alloc_complex(len)},
      in_2{symphas::dft::fftw_alloc_complex(len)},
      p_in_out_0{new_fftw_plan<3, complex_t, complex_t>{}(in_0, out_0, dims,
                                                          true, false)},
      p_in_out_1{new_fftw_plan<3, complex_t, complex_t>{}(in_1, out_1, dims,
                                                          true, false)},
      p_out_in{new_fftw_plan<3, complex_t, complex_t>{}(in_2, out_2, dims,
                                                        false, true)} {}

}  // namespace expr

#else

namespace expr {

//! Delegation method for computing convolution of different types.
/*!
 * A computation object for the convolution
 * so that the FFTW routines can be used on the data to compute
 * from real space to Fourier space and back.
 */
template <size_t D>
struct ConvolutionData {
  ConvolutionData() : out_0{0}, in_1{0} {}

  template <typename T, typename S>
  ConvolutionData(T* in_0, S* out_1, len_type* dims, len_type len)
      : out_0{new complex_t[len]{}}, in_1{new complex_t[len]{}} {}
  template <typename T, typename S>
  ConvolutionData(Grid<T, D> const& in_0, S* out_1)
      : ConvolutionData(in_0.values, out_1, in_0.dims, in_0.len) {}
  template <typename T, typename S>
  ConvolutionData(T* in_0, Grid<S, D>& out_1)
      : ConvolutionData(in_0, out_1.values, out_1.dims, out_1.len) {}
  template <typename T, typename S>
  ConvolutionData(Grid<T, D> const& in_0, Grid<S, D>& out_1)
      : ConvolutionData(in_0.values, out_1.values, out_1.dims, out_1.len) {}

  ConvolutionData(ConvolutionData<D> const& other) = delete;
  ConvolutionData(ConvolutionData<D>&& other) { swap(*this, other); }
  ConvolutionData<D>& operator=(ConvolutionData<D> const&) = delete;

  friend void swap(ConvolutionData<D>& first, ConvolutionData<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.in_1, second.in_1);
  }

  template <typename T_0>
  void transform_in_out(T_0* in_0, const len_type* dims) const {
    symphas::dft::dft<D>(in_0, out_0, dims, true);
  }

  template <typename T_0>
  void transform_out_in(T_0* out_1, const len_type* dims) const {
    symphas::dft::dft<D>(in_1, out_1, dims, true);
  }

  complex_t* out_0;  //!< Fourier transform of the given data.
  complex_t* in_1;   //!< Input of Fourier data.

  ~ConvolutionData() {
    delete[] out_0;
    delete[] in_1;
  }
};

//! Delegation method for computing convolution of different types.
/*!
 * A computation object for the convolution
 * so that the FFTW routines can be used on the data to compute
 * from real space to Fourier space and back. Performs convolution when there
 * are two input fields and one output field.
 */
template <size_t D>
struct ConvolutionDataPair {
  ConvolutionDataPair() : out_0{0}, out_1{0}, in_2{0} {}

  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(T_0* in_0, T_1* in_1, R* out_2, len_type* dims,
                      len_type len)
      : out_0{new complex_t[len]{}},
        out_1{new complex_t[len]{}},
        in_2{new complex_t[len]{}} {}
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(Grid<T_0, D> const& in_0, Grid<T_1, D> const& in_1,
                      R* out_2)
      : ConvolutionDataPair(in_0.values, in_1.values, out_2, in_0.dims,
                            in_0.len) {}
  template <typename T_0, typename T_1, typename R>
  ConvolutionDataPair(T_0* in_0, T_1* in_1, Grid<R, D>& out_2)
      : ConvolutionDataPair(in_0, in_1, out_2.values, out_2.dims, out_2.len) {}

  ConvolutionDataPair(ConvolutionDataPair<D> const&) = delete;
  ConvolutionDataPair(ConvolutionDataPair<D>&& other) : ConvolutionDataPair() {
    swap(*this, other);
  }
  ConvolutionDataPair<D>& operator=(ConvolutionDataPair<D> const&) = delete;

  friend void swap(ConvolutionDataPair<D>& first,
                   ConvolutionDataPair<D>& second) {
    using std::swap;
    swap(first.out_0, second.out_0);
    swap(first.out_1, second.out_1);
    swap(first.in_2, second.in_2);
  }

  template <typename T_0>
  void transform_in_out_0(T_0* in_0, const len_type* dims) const {
    symphas::dft::dft<D>(in_0, out_0, dims, true);
  }

  template <typename T_0>
  void transform_in_out_1(T_0* in_1, const len_type* dims) const {
    symphas::dft::dft<D>(in_1, out_1, dims, true);
  }

  template <typename T_0>
  void transform_out_in(T_0* out_2, const len_type* dims) const {
    symphas::dft::dft<D>(in_2, out_2, dims, true);
  }

  ~ConvolutionDataPair() {
    delete[] out_0;
    delete[] out_1;
    delete[] in_2;
  }

  /* for finding the convolution from real space, back into real space, we need
   * to first compute it in fourier space
   */
  complex_t* out_0;  //!< Fourier transform of the first given data.
  complex_t* out_1;  //!< Fourier transform of the second given data.
  complex_t* in_2;   //!< Input of Fourier data.
};

}  // namespace expr

#endif
