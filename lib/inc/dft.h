
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
 * PURPOSE: Adds additional necessary features for applying Fourier transforms.
 *
 * ***************************************************************************
 */

#pragma once

#include <vector>

#include "spslib.h"

namespace symphas {
//! Contains elements related to discrete fourier transforms.
/*!
 * Defines functions for computing the discrete fourier transform, as well
 * as functions related to the manipulation of discrete fourier transform data.
 */
namespace dft {}
}  // namespace symphas

namespace symphas::dft {

iter_type next_block_i(iter_type i, len_type full_len, len_type block_count);

void long_dft(scalar_t* data, complex_t* out, iter_type L,
              bool backward = false);
void long_dft(scalar_t* data, complex_t* out, iter_type L, iter_type M,
              bool backward = false);
void long_dft(scalar_t* data, complex_t* out, iter_type L, iter_type M,
              iter_type N, bool backward = false);
void long_dft(complex_t* data, complex_t* out, iter_type L,
              bool backward = false);
void long_dft(complex_t* data, complex_t* out, iter_type L, iter_type M,
              bool backward = false);
void long_dft(complex_t* data, complex_t* out, iter_type L, iter_type M,
              iter_type N, bool backward = false);

template <size_t D>
struct long_dft_call;
//
//	template<>
//	struct long_dft_call<1>
//	{
//		template<typename T>
//		void operator()(const axis_nd_t<1>*, const T* data_y, complex_t*
// out, len_type len, bool backward = false)
//		{
//			long_dft(data_y, out, len, backward);
//		}
//	};
//
//	template<>
//	struct long_dft_call<2>
//	{
//		template<typename T>
//		void operator()(const axis_nd_t<2>* data_x, const T* data_y,
// complex_t* out, len_type len, bool backward = false)
//		{
//			auto [L, M] = symphas::lib::get_dimensions<2>(data_x,
// len)._2(); 			long_dft(data_y, out, L, M, backward);
//		}
//	};
//
//	template<>
//	struct long_dft_call<3>
//	{
//		template<typename T>
//		void operator()(const axis_nd_t<3>* data_x, const T* data_y,
// complex_t* out, len_type len, bool backward = false)
//		{
//			auto [L, M, N] = symphas::lib::get_dimensions<3>(data_x,
// len)._3(); 			long_dft(data_y, out, L, M, N, backward);
//		}
//	};
//
template <size_t D, typename T>
std::vector<std::pair<axis_nd_t<D>, complex_t>> long_dft(
    std::vector<std::pair<axis_nd_t<D>, T>> const& data,
    bool backward = false) {
  complex_t* out_y = new complex_t[data.size()];
  T* in_y = new T[data.size()];
  axis_nd_t<D>* out_x = new axis_nd_t<D>[data.size()];

  fill_x_axis(out_x, out_x, data.size());
  long_dft_call<D>{}(out_x, in_y, out_y, data.size(), backward);

  for (iter_type i = 0; i < data.size(); ++i) {
    in_y = data[i].first;
    out_x = data[i].second;
  }

  std::vector<std::pair<typename axis_nd<D>::type, complex_t>> out(data.size());
  for (iter_type i = 0; i < data.size(); ++i) {
    out[i] = {out_x[i], out_y[i]};
  }
  return out;
}
}  // namespace symphas::dft

// ****************************************************************************************

namespace grid {

/*
 * functions to scale the data after fourier transforms
 */

template <typename T>
void scale(T* data_y, len_type len, double r) {
#pragma omp parallel for
  for (iter_type i = 0; i < len; ++i) {
    data_y[i] *= r;
  }
}

template <typename T>
void scale(T (*data_y)[2], len_type len, double r) {
#pragma omp parallel for
  for (iter_type i = 0; i < len; ++i) {
    data_y[i][0] *= r;
    data_y[i][1] *= r;
  }
}

template <typename T>
void scale(T* data_y, len_type len) {
  scale(data_y, len, 1.0 / len);
}

template <typename T>
void scale(T* data_y, len_type L, len_type M) {
  scale(data_y, L * M);
}

template <typename T>
void scale(T* data_y, len_type L, len_type M, len_type N) {
  scale(data_y, L * M * N);
}

template <size_t D, typename T>
void scale(T* data_y, const len_type* dims) {
  len_type len = 1;
  for (iter_type i = 0; i < D; ++i) {
    len *= dims[i];
  }

  scale(data_y, len);
}

}  // namespace grid
