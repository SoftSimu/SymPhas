
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
 * PURPOSE: Defines Fourier transform algorithms that can be used in solvers.
 *
 * ***************************************************************************
 */

#pragma once

#include "dft.h"

namespace symphas::dft {

//! Arrange FFTW data from real dft output to in place format.
/*!
 * From a real to real type of arrangement, copy the data into the
 * target array as an in-place type of arrangement.
 *
 * This does not support in-place.
 *
 * \param src The array which is arranged sequentially, and will be copied
 * from into target.
 * \param target The data here is copied from src according to an
 * arrangement that allows FFTW to perform an in place transform.
 * \param dims The underlying dimensions of the complex data.
 */
template <size_t D>
void arrange_fftw_rtip(const scalar_t* src, scalar_t* target,
                       const len_type* dims);

//! Arrange FFTW data from in-place to real dft input format.
/*!
 * From an real to real type of arrangement, copy the data into the
 * target array as an in-place type of arrangement.
 *
 * This does not support in-place.
 *
 * \param src The array which is arranged sequentially, and will be copied
 * from into target.
 * \param target The data here is copied from src according to an
 * arrangement that allows FFTW to perform an in place transform.
 * \param dims The underlying dimensions of the complex data.
 */
template <size_t D>
void arrange_fftw_iptr(const scalar_t* src, scalar_t* target,
                       const len_type* dims);

//! Arrange FFTW data from sequential to in place format.
/*!
 * From source, copy the data to the target array. This arranges the data
 * in target so that the FFTW transform for an in place r2c transform
 * in target can produce the correct results.
 *
 * This algorithm supports the case when src and target are the same array.
 *
 * See arrange_fftw_ipts, which performs the opposite transposition.
 *
 * \param src The array which is arranged sequentially, and will be copied
 * from into target.
 * \param target The data here is copied from src according to an
 * arrangement that allows FFTW to perform an in place transform.
 * \param dims The underlying dimensions of the complex data.
 */
template <size_t D>
void arrange_fftw_stip(const scalar_t* src, scalar_t* target,
                       const len_type* dims);

//! Arrange FFTW data from in place format to sequential format.
/*!
 * From source, copy the data to the target array. This arranges the data
 * in target so that it is sequential, as it assumes that the source
 * array is in the format given by an in place FFTW transform.
 *
 * This algorithm supports the case when src and target are the same array.
 *
 * See arrange_fftw_stip, which performs the opposite transposition.
 *
 * \param src The array which is arranged according to the FFTW format for
 * in place Fourier transforms.
 * \param target The data here is copied from src into a sequential
 * arrangement suitable for regular array use cases.
 * \param dims The underlying dimensions of the complex data.
 */
template <size_t D>
void arrange_fftw_ipts(scalar_t* src, scalar_t* target, const len_type* dims);

//! Arranges FFTW half complex output into a (full) sequential array.
/*!
 * When an FFTW r2c plan has been executed, the complex array storing the
 * result will only store the non-redundant half. In some cases, the whole
 * complex array should be represented, and this algorithm rearranges the
 * values of the array to fit from the result into the desired usable shape.
 *
 * The source and target array may be the same
 *
 * \param src The source array which is in half complex output.
 * \param target The array that is rearranged so that it is expanded into
 * the full dimensions.
 * \param dims The underlying dimensions of the data.
 */
template <size_t D>
void arrange_fftw_hcts(complex_t* src, complex_t* target, const len_type* dims);

//! Arranges a (full) sequential array into FFTW half complex output.
/*!
 * When an FFTW r2c plan has been executed, the complex array storing the
 * result will only store the non-redundant half. In some cases, the whole
 * complex array should be represented, and this algorithm rearranges the
 * values of the array to fit from the result into the desired usable shape.
 *
 * The source and target array may be the same
 *
 * \param src The source array which is the source complex output.
 * \param target The array that is rearranged so that it is expanded into
 * the full dimensions.
 * \param dims The underlying dimensions of the data.
 */
template <size_t D>
void arrange_fftw_sthc(complex_t* src, complex_t* target, const len_type* dims);

namespace {

//! Executes a function with consideration of an FFTW arrangement.
/*!
 * When using FFTW transformations between complex- and real-valued
 * arrays, the data of the real-valued grid corresponds to only the
 * non-duplicated half of the complex (Fourier space) array. Thus, an
 * iteration must iterate to skip over the duplicated values when
 * iterating over the complex array.
 *
 * The template parameter specializes for either complex or real types,
 * and when the real type is provided, the iteration will skip the
 * duplicated part, and otherwise there will be no skipping.
 *
 * The function given to the operator member is expected to have two
 * indices, the first is the index in the complex array which will be
 * skipped on duplicated rows when the given type is ::scalar_t.
 * The second index is the index is always sequential.
 *
 * If the iteration is for a complex type, it will be entirely
 * sequential (both indices will always match).
 *
 * \tparam T The type of the grid values.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct iterate_complex_dup {
  template <typename F>
  void operator()(F f, const len_type* dim) {
    len_type len = grid::length<D>(dim);
#pragma omp parallel for
    for (iter_type i = 0; i < len; ++i) {
      f(i, i);
    }
  }

  template <typename F>
  void operator()(T* values, F f, const len_type* dim) {
    len_type len = grid::length<D>(dim);
#pragma omp parallel for
    for (iter_type i = 0; i < len; ++i) {
      values[i] = f(i);
    }
  }
};

}  // namespace

//! Applies the binary function over matching indices.
/*!
 * See symphas::dft::iterate_complex_dup. Iterates through by skipping over
 * duplicated data in each row, which is chosen by the source type.
 */
template <typename T, size_t D, typename F>
void iterate_rc(F&& f, len_type const* dims) {
  iterate_complex_dup<T, D>{}(std::forward<F>(f), dims);
}

//! Applies the binary function over matching indices.
/*!
 * See symphas::dft::iterate_complex_dup. Iterates through by skipping over
 * duplicated data in each row, which is chosen by the source type.
 */
template <typename T, size_t D, typename F>
void iterate_rc(T* values, F&& f, len_type const* dims) {
  iterate_complex_dup<T, D>{}(values, std::forward<F>(f), dims);
}
}  // namespace symphas::dft
