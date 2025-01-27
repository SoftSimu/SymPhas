
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
 * PURPOSE: Defines a custom FFT implementation (unused).
 *
 * ***************************************************************************
 */

#pragma once

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include <complex>
#include <vector>

#include "dft.h"
#include "symphasthread.h"

#ifdef __NO_DEF__

//! Macro for calling custom FFT implementation.
#define FFT_2D_BASE(IN, OUT, L_2_POW, M_2_POW)    \
  FFT_ID(L_2_POW, M_2_POW, 0) : {                 \
    symphas::dft::fft<L_2_POW, M_2_POW>(IN, OUT); \
    break;                                        \
  }

//! Macro for calling custom FFT implementation.
#define FFT_3D_BASE(IN, OUT, L_2_POW, M_2_POW, N_2_POW)    \
  FFT_ID(L_2_POW, M_2_POW, N_2_POW) : {                    \
    symphas::dft::fft<L_2_POW, M_2_POW, N_2_POW>(IN, OUT); \
    break;                                                 \
  }

#define _MM_TRANSPOSE4_PD(row0, row1, row2, row3)      \
  {                                                    \
    __m256d tmp3, tmp2, tmp1, tmp0;                    \
                                                       \
    tmp0 = _mm256_shuffle_pd((row0), (row1), 0x0);     \
    tmp2 = _mm256_shuffle_pd((row0), (row1), 0xF);     \
    tmp1 = _mm256_shuffle_pd((row2), (row3), 0x0);     \
    tmp3 = _mm256_shuffle_pd((row2), (row3), 0xF);     \
                                                       \
    (row0) = _mm256_permute2f128_pd(tmp0, tmp1, 0x20); \
    (row1) = _mm256_permute2f128_pd(tmp2, tmp3, 0x20); \
    (row2) = _mm256_permute2f128_pd(tmp0, tmp1, 0x31); \
    (row3) = _mm256_permute2f128_pd(tmp2, tmp3, 0x31); \
  }

#define FFT_ID(L_DIM, M_DIM, N_DIM)                   \
  (symphas::internal::fft_depth(N_DIM) << 8ull) +     \
      (symphas::internal::fft_depth(M_DIM) << 4ull) + \
      (symphas::internal::fft_depth(L_DIM))

#define FFT_ID_RUN_3D(L_DIM, M_DIM, N_DIM)                   \
  ((static_cast<size_t>(log2(N_DIM)) - log2(N_DIM) == 0)     \
       ? (static_cast<size_t>(log2(N_DIM)) << 8)             \
       : 0) +                                                \
      ((static_cast<size_t>(log2(M_DIM)) - log2(M_DIM) == 0) \
           ? (static_cast<size_t>(log2(M_DIM)) << 4)         \
           : 0) +                                            \
      ((static_cast<size_t>(log2(L_DIM)) - log2(L_DIM) == 0) \
           ? (static_cast<size_t>(log2(L_DIM)))              \
           : 0)

#define FFT_ID_RUN_2D(L_DIM, M_DIM)                          \
  ((static_cast<size_t>(log2(M_DIM)) - log2(M_DIM) == 0)     \
       ? (static_cast<size_t>(log2(M_DIM)) << 4)             \
       : 0) +                                                \
      ((static_cast<size_t>(log2(L_DIM)) - log2(L_DIM) == 0) \
           ? (static_cast<size_t>(log2(L_DIM)))              \
           : 0)

namespace symphas::dft {

void fft(double* data_real, double* data_imag, std::complex<double>* out,
         size_t LL, size_t MM);
void fft(double* data_real, double* data_imag, std::complex<double>* out,
         size_t LL, size_t MM, size_t NN);

template <size_t LL, size_t MM, typename T>
void fft(T* data, complex_t* out);
template <size_t LL, size_t MM, size_t NN, typename T>
void fft(T* data, complex_t* out);
}  // namespace symphas::dft

/*
 *
 *
 *
 *
 *
 *
 *
 *
 * implementation of the symphas fast fourier transform
 */

namespace symphas::internal {

/*!
 * Determines the depth of the recursion required for
 * the fft algorithm at runtime.
 */
inline constexpr size_t fft_depth(size_t N) {
  if (N > 1) {
    return 1 + fft_depth(N >> 1);
  } else {
    return 0;
  }
}

/*!
 * Inplace calculation of the discrete fourier transform using
 * the Cooley�Tukey FFT algorithm.
 */
inline void fft_alg(double* const(&r), double* const(&c), size_t N) {
  size_t M = fft_depth(N);
  for (size_t l = 0; l < M; ++l) {
    size_t const le = 1ull << (M - static_cast<size_t>(l)), le2 = le >> 1;

    std::complex<double> u{1.0, 0.0}, s{cos(symphas::PI / le2), -sin(PI / le2)};

    for (size_t j = 0; j < le2; ++j) {
      for (size_t i = j; i < N; i += le) {
        size_t ip = i + le2;

        double tr = r[i] + r[ip];
        double tc = c[i] + c[ip];

        auto w = u * std::complex<double>{r[i] - r[ip], c[i] - c[ip]};
        r[ip] = w.real();
        c[ip] = w.imag();

        r[i] = tr;
        c[i] = tc;
      }
      u *= s;
    }
  }

  size_t nd2 = N >> 1ull, nm1 = N - 1;

  size_t j = 0;

  for (size_t i = 0; i < nm1; ++i) {
    if (i < j) {
      double tr = r[j];
      r[j] = r[i];
      r[i] = tr;

      double tc = c[j];
      c[j] = c[i];
      c[i] = tc;
    }

    size_t k = nd2;
    while (k < j + 1) {
      j -= k;
      k = k >> 1;
    }

    j += k;
  }
}

/*!
 * Uses the SSE library to quickly compute the transpose
 * of a 4x4 square of doubles.
 */
inline void transpose_4x4_SSE(double* A, double* B, size_t N, size_t M) {
  __m256d row1 = _mm256_load_pd(&A[0 * N]);
  __m256d row2 = _mm256_load_pd(&A[1 * N]);
  __m256d row3 = _mm256_load_pd(&A[2 * N]);
  __m256d row4 = _mm256_load_pd(&A[3 * N]);
  _MM_TRANSPOSE4_PD(row1, row2, row3, row4);
  _mm256_store_pd(&B[0 * M], row1);
  _mm256_store_pd(&B[1 * M], row2);
  _mm256_store_pd(&B[2 * M], row3);
  _mm256_store_pd(&B[3 * M], row4);
}

/*!
 * Determine the tranpose of first matrix and put it in
 * the second matrix.
 */
inline void tr_matrix(double* A, double* B, size_t N, size_t M,
                      const int block_size) {
  for (size_t i = 0; i < N; i += block_size) {
    for (size_t j = 0; j < M; j += block_size) {
      size_t max_i2 = (i + block_size < N) ? i + block_size : N;
      size_t max_j2 = (j + block_size < M) ? j + block_size : M;
      for (size_t i2 = i; i2 < max_i2; i2 += 4) {
        for (size_t j2 = j; j2 < max_j2; j2 += 4) {
          transpose_4x4_SSE(&A[i2 * N + j2], &B[j2 * M + i2], N, M);
        }
      }
    }
  }
}

/*
 * runs the fft algorithm on the offset array, a helper for
 * threads
 */
inline void fft_thr_run(iter_type i, double* const(&r), double* const(&c),
                        size_t B) {
  fft_alg(r + i * B, c + i * B, B);
}

/*
 * unrolls all the required fft loops, a helper for threads
 */
inline void fft_thr_loop(double* const(&r), double* const(&c), size_t A,
                         size_t B, iter_type i) {
  for (iter_type n = ThreadPool::thr_start_i(A, i);
       n < ThreadPool::thr_start_i(A, i + 1); ++n) {
    fft_thr_run(n, r, c, B);
  }
}

/*
 * the fft job done by one of the threads
 * it performs the row-wise fft for each slice of the 3d system to recover the
 * full fourier transform
 */
inline void fft_thr_job_3(double*((&r)[2]), double*((&c)[2]), ThreadPool& thr,
                          std::complex<double>* out, size_t LL, size_t MM,
                          size_t NN, iter_type i) {
  for (iter_type a = 0; a < NN; ++a) {
    fft_thr_loop(r[0] + MM * LL * a, c[0] + MM * LL * a, MM, LL, a);
  }
  thr.sync(i, 1);
  for (iter_type a = 0; a < NN; ++a) {
    fft_thr_loop(r[1] + MM * LL * a, c[1] + MM * LL * a, LL, MM, a);
  }
  thr.sync(i, 2);
  for (iter_type a = 0; a < LL; ++a) {
    fft_thr_loop(r[0] + MM * NN * a, c[0] + MM * NN * a, NN, MM, a);
  }
  thr.sync(i, 3);

  size_t n = ThreadPool::thr_start_i(LL * MM, i),
         e = ThreadPool::thr_start_i(LL * MM, i + 1),
         offset_z = ((LL * MM * NN) >> 1ull), offset_y = ((LL * MM) >> 1ull),
         offset_x = LL >> 1ull;

  auto *it = out + n + offset_z + offset_y + offset_x,
       *last = out + LL * MM * NN;

  // the indexing is changed here in order to align the values with the correct
  // axes in the native list, the dft of the center is put at the corners

  iter_type nn = static_cast<iter_type>(n),
            ii = static_cast<iter_type>(n) % static_cast<iter_type>(LL),
            jj = (static_cast<iter_type>(n) / static_cast<iter_type>(LL)) %
                 static_cast<iter_type>(MM);
  for (; nn < e && it < last + offset_x + offset_y; ++nn) {
    if (jj >= offset_y) {
      if (ii++ >= offset_x) {
        *(it - LL * MM - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      } else {
        *(it - LL * MM) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      }
    } else {
      if (ii++ >= offset_x) {
        *(it - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      } else {
        *(it) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      }
    }

    ++it;
    ii %= LL;
    if (ii == 0) ++jj;
    jj %= MM;
  }

  if (it >= last) {
    it = out + (it - last);
  }
  for (; nn < e; ++nn) {
    if (jj >= offset_y) {
      if (ii++ >= offset_x) {
        *(it - LL * MM - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      } else {
        *(it - LL * MM) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      }
    } else {
      if (ii++ >= offset_x) {
        *(it - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      } else {
        *(it) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
      }
    }

    ++it;
    ii %= LL;
    if (ii == 0) ++jj;
    jj %= MM;
  }
  thr.sync(i, 4);
}

/*
 * the fft job done by one of the threads thread
 * it takes the row-wise slice, transposes it and then takes the row-wise
 * (originally columns) slice
 */
inline void fft_thr_job_2(double*((&r)[2]), double*((&c)[2]), ThreadPool& thr_p,
                          std::complex<double>* out, size_t LL, size_t MM,
                          iter_type i) {
  fft_thr_loop(r[0], c[0], MM, LL, i);
  thr_p.sync(i, 1);
  fft_thr_loop(r[1], c[1], LL, MM, i);
  thr_p.sync(i, 2);

  size_t n = ThreadPool::thr_start_i(LL * MM, i),
         e = ThreadPool::thr_start_i(LL * MM, i + 1),
         offset_y = ((LL * MM) >> 1ull), offset_x = LL >> 1ull;

  auto *it = out + n + offset_y + offset_x, *last = out + LL * MM;

  // the indexing is changed here in order to align the values with the correct
  // axes in the native list, the dft of the center is put at the corners

  iter_type nn = static_cast<iter_type>(n),
            ii = static_cast<iter_type>(n) % static_cast<iter_type>(LL);
  for (; nn < e && it < last + offset_x; ++nn) {
    if (ii++ >= offset_x) {
      *(it - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
    } else {
      *it = std::complex<double>{(r[1])[nn], (c[1])[nn]};
    }

    ++it;
    ii %= LL;
  }

  if (it >= last) {
    it = out + (it - last);
  }
  for (; nn < e; ++nn) {
    if (ii++ >= offset_x) {
      *(it - LL) = std::complex<double>{(r[1])[nn], (c[1])[nn]};
    } else {
      *it = std::complex<double>{(r[1])[nn], (c[1])[nn]};
    }

    ++it;
    ii %= LL;
  }
  thr_p.sync(i, 3);
}

template <size_t LEN>
inline void copy_in(const scalar_t* data, scalar_t real_out[LEN],
                    scalar_t imag_out[LEN]) {
  size_t i = 0;
  for (const scalar_t* it = data; it < data + LEN; ++it, ++i) {
    real_out[i] = *it;
    imag_out[i] = 0;
  }
}

template <size_t LEN>
inline void copy_in(const complex_t* data, scalar_t real_out[LEN],
                    scalar_t imag_out[LEN]) {
  size_t i = 0;
  for (const complex_t* it = data; it < data + LEN; ++it, ++i) {
    real_out[i] = it->real();
    imag_out[i] = it->imag();
  }
}
}  // namespace symphas::internal

inline void symphas::dft::fft(double* data_real, double* data_imag,
                              std::complex<double>* out, size_t LL, size_t MM) {
  symphas::Time t("time for FFT (2d)");
  ThreadPool thr_p;

  double* dftr_tr = new double[LL * MM];
  double* dftc_tr = new double[LL * MM];

  using arr = double* [2];
  arr pr{data_real, dftr_tr};
  arr pc{data_imag, dftc_tr};

  thr_p.start_thr(symphas::internal::fft_thr_job_2, std::ref(pr), std::ref(pc),
                  std::ref(thr_p), out, LL, MM);
  for (iter_type i = 0; i < 3; ++i) {
    thr_p.main_thr_sync([&]() {
      if (i == 0) {
        symphas::internal::tr_matrix(data_real, dftr_tr, LL, MM, 16);
        symphas::internal::tr_matrix(data_imag, dftc_tr, LL, MM, 16);
      }
    });
  }
  thr_p.close_thr();
  delete[] dftr_tr;
  delete[] dftc_tr;
}

inline void symphas::dft::fft(double* data_real, double* data_imag,
                              std::complex<double>* out, size_t LL, size_t MM,
                              size_t NN) {
  symphas::Time t("time for FFT (3d)");
  ThreadPool thr_p;

  double* dftr_tr = new double[LL * MM * NN];
  double* dftc_tr = new double[LL * MM * NN];

  using arr = double* [2];
  arr pr{data_real, dftr_tr};
  arr pc{data_imag, dftc_tr};

  thr_p.start_thr(symphas::internal::fft_thr_job_3, std::ref(pr), std::ref(pc),
                  std::ref(thr_p), out, LL, MM, NN);
  for (iter_type i = 0; i < 4; ++i) {
    thr_p.main_thr_sync([&]() {
      if (i == 0) {
        for (iter_type n = 0; n < NN; ++n) {
          symphas::internal::tr_matrix(data_real + MM * LL * n,
                                       dftr_tr + MM * LL * n, LL, MM, 16);
          symphas::internal::tr_matrix(data_imag + MM * LL * n,
                                       dftc_tr + MM * LL * n, LL, MM, 16);
        }
      } else if (i == 1) {
        for (iter_type x = 0; x < NN; ++x) {
          for (iter_type y = 0; y < MM; ++y) {
            for (iter_type z = 0; z < NN; ++z) {
              (data_real + MM * NN * x)[NN * y + z] =
                  dftr_tr[x + LL * MM * z + LL * y];
              (data_imag + MM * NN * x)[NN * y + z] =
                  dftc_tr[x + LL * MM * z + LL * y];
            }
          }
        }
      }
    });
  }
  thr_p.close_thr();
  delete[] dftr_tr;
  delete[] dftc_tr;
}

template <size_t LL, size_t MM, typename T>
void symphas::dft::fft(T* data, complex_t* out) {
  scalar_t* data_real = new scalar_t[LL * MM];
  scalar_t* data_imag = new scalar_t[LL * MM];

  symphas::internal::copy_in<LL * MM>(data, data_real, data_imag);
  symphas::dft::fft(data_real, data_imag, out, LL, MM);

  delete[] data_real;
  delete[] data_imag;
}

template <size_t LL, size_t MM, size_t NN, typename T>
void symphas::dft::fft(T* data, complex_t* out) {
  scalar_t* data_real = new scalar_t[LL * MM * NN];
  scalar_t* data_imag = new scalar_t[LL * MM * NN];

  symphas::internal::copy_in<LL * MM * NN>(data, data_real, data_imag);
  symphas::dft::fft(data_real, data_imag, out, LL, MM, NN);

  delete[] data_real;
  delete[] data_imag;
}

#endif
