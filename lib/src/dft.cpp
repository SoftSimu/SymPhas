
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
 */

#include "dft.h"

iter_type symphas::dft::next_block_i(iter_type i, len_type full_len,
                                     len_type block_count) {
  iter_type count = full_len / block_count,
            rem = full_len - (count * block_count);

  return count * i + std::min(i, rem);
}

template <typename T>
inline void long_dft(T* data, complex_t* out, len_type len,
                     bool backward = false) {
  // int const dL = static_cast<int>(len) / 2;
  for (iter_type ii = 0; ii < len; ++ii) {
    double ti = 0, p_ti = (2.0 * symphas::PI / len) * (ii);
    for (iter_type i = 0; i < len; ++i, ti += p_ti) {
      out[ii] += data[i] * complex_t{cos(ti), -sin(ti)};
    }
  }
}

template <typename T>
inline void long_dft(T* data, complex_t* out, iter_type L, iter_type M,
                     bool backward = false) {
  int sign = (backward) ? 1 : -1;

  len_type len = L * M;
  // int const
  //	dL = L / 2,
  //	dM = M / 2;

  auto long_dft = [&](iter_type thr_i, len_type thr_n) {
    iter_type start_i = symphas::dft::next_block_i(thr_i, len, thr_n),
              end_i = symphas::dft::next_block_i(thr_i + 1, len, thr_n);
    int i0 = start_i % L, j0 = start_i / L;

    iter_type jj = j0, ii = i0;
    iter_type index = start_i;
    for (; jj < M; ++jj) {
      for (; ii < L; ++ii) {
        complex_t yy;
        double tj = 0, p_tj = (2.0 * symphas::PI / M) * (jj);
        for (iter_type j = 0; j < M; ++j, tj += p_tj) {
          complex_t xx;
          double ti = 0, p_ti = (2.0 * symphas::PI / L) * (ii);
          for (iter_type i = 0; i < L; ++i, ti += p_ti) {
            xx += data[i + j * L] * complex_t{cos(ti), sign * sin(ti)};
          }
          yy += xx * complex_t{cos(tj), sign * sin(tj)};
        }
        out[ii + jj * L] = yy;
        if (++index == end_i) return;
      }
      ii = 0;
    }
  };

#pragma omp parallel for
  for (iter_type i = 0; i < L; ++i) {
    long_dft(i, L);
  }
}

template <typename T>
inline void long_dft(T* data, complex_t* out, iter_type L, iter_type M,
                     iter_type N, bool backward = false) {
  int sign = (backward) ? 1 : -1;

  len_type len = L * M * N;

  auto long_dft = [&](iter_type thr_i, len_type thr_n) {
    iter_type start_i = symphas::dft::next_block_i(thr_i, len, thr_n),
              end_i = symphas::dft::next_block_i(thr_i + 1, len, thr_n);
    int i0 = (start_i % L), j0 = ((start_i / L) % M), k0 = (start_i / (L * M));

    iter_type kk = k0, jj = j0, ii = i0;
    iter_type index = start_i;
    for (; kk < N; ++kk) {
      for (; jj < M; ++jj) {
        for (; ii < L; ++ii) {
          complex_t zz;
          double p_tk = (2.0 * symphas::PI / N) * (kk), tk = 0;
          for (iter_type k = 0; k < N; ++k, tk += p_tk) {
            complex_t yy;
            double p_tj = (2.0 * symphas::PI / M) * (jj), tj = 0;

            for (iter_type j = 0; j < M; ++j, tj += p_tj) {
              complex_t xx;
              double p_ti = (2.0 * symphas::PI / L) * (ii), ti = 0;

              for (iter_type i = 0; i < L; ++i, ti += p_ti) {
                xx += data[i + j * L + k * L * M] *
                      complex_t(cos(ti), sign * sin(ti));
              }
              yy += xx * complex_t(cos(tj), sign * sin(tj));
            }
            zz += yy * complex_t(cos(tk), sign * sin(tk));
          }
          out[ii + jj * L + kk * L * M] = zz;
          if (++index == end_i) return;
        }
        ii = 0;
      }
      jj = 0;
    }
  };

#pragma omp parallel for
  for (iter_type i = 0; i < L; ++i) {
    long_dft(i, L);
  }
}

void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L,
                            bool backward) {
  ::long_dft(data, out, L, backward);
}
void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L,
                            iter_type M, bool backward) {
  ::long_dft(data, out, L, M, backward);
}
void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L,
                            iter_type M, iter_type N, bool backward) {
  ::long_dft(data, out, L, M, N, backward);
}

void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L,
                            bool backward) {
  ::long_dft(data, out, L, backward);
}
void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L,
                            iter_type M, bool backward) {
  ::long_dft(data, out, L, M, backward);
}
void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L,
                            iter_type M, iter_type N, bool backward) {
  ::long_dft(data, out, L, M, N, backward);
}
