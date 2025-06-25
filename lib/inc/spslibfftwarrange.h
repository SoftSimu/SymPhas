
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


//! Specialization of iterate_complex_dup.
template <>
struct symphas::dft::iterate_complex_dup<scalar_t, 2> {
  template <typename F>
  void operator()(F&& f, const len_type* dim) {
    len_type Lh = dim[0] / 2 + 1;

#pragma omp parallel for
    for (iter_type j = 0; j < dim[1]; ++j) {
      iter_type skip_index =
          j * dim[0];                // Iterated across the whole dimension.
      iter_type seq_index = j * Lh;  // Iterated sequentially.
      for (iter_type i = 0; i < Lh; ++i) {
        f(skip_index + i, seq_index + i);
      }
    }
  }

  template <typename F>
  void operator()(scalar_t* values, F&& f, const len_type* dim) {
    const len_type Lh = dim[0] / 2 + 1;

#pragma omp parallel for shared(f)
    for (iter_type j = 0; j < dim[1]; ++j) {
      iter_type skip_index = j * dim[0];
      iter_type seq_index = j * Lh;
      for (iter_type i = 0; i < Lh; ++i) {
        values[seq_index + i] = f(skip_index + i);
      }
    }
  }
};

//! Specialization of iterate_complex_dup.
template <>
struct symphas::dft::iterate_complex_dup<scalar_t, 3> {
  template <typename F>
  void operator()(F&& f, const len_type* dim) {
    len_type Lh = dim[0] / 2 + 1;

#pragma omp parallel for
    for (iter_type kj = 0; kj < dim[1] * dim[2]; ++kj) {
      iter_type skip_index = kj * dim[0];
      iter_type seq_index = kj * Lh;
      for (iter_type i = 0; i < Lh; ++i) {
        f(skip_index + i, seq_index + i);
      }
    }
  }

  template <typename F>
  void operator()(scalar_t* values, F&& f, const len_type* dim) {
    len_type Lh = dim[0] / 2 + 1;

#pragma omp parallel for
    for (iter_type kj = 0; kj < dim[1] * dim[2]; ++kj) {
      iter_type skip_index = kj * dim[0];
      iter_type seq_index = kj * Lh;
      for (iter_type i = 0; i < Lh; ++i) {
        values[seq_index + i] = f(skip_index + i);
      }
    }
  }
};

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_rtip<1>(const scalar_t* src, scalar_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_rtip<2>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  int hsize = dims[0] / 2 + 1;
  for (int j = 0; j < dims[1]; ++j) {
    for (int i = 0; i < hsize; ++i) {
      if (j <= hsize) {
        target[i + j * hsize] = src[i + j * dims[0]];  // Real part
      } else {
        target[i + j * hsize] = src[i + (dims[1] - j) * dims[0]];  // Real part
      }
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_rtip<3>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  int hsize = dims[0] / 2 + 1;
  for (int k = 0; k < dims[2]; ++k) {
    for (int j = 0; j < dims[1]; ++j) {
      for (int i = 0; i < hsize; ++i) {
        if (k <= hsize) {
          if (j <= hsize) {
            target[i + j * hsize + k * dims[1] * hsize] =
                src[i + j * dims[0] + k * dims[0] * dims[1]];  // Real part
          } else {
            target[i + j * hsize + k * dims[1] * hsize] =
                src[i + (dims[1] - j) * dims[0] +
                    k * dims[0] * dims[1]];  // Real part
          }
        } else {
          if (j <= hsize) {
            target[i + j * hsize + k * dims[1] * hsize] =
                src[i + j * dims[0] +
                    (dims[2] - k) * dims[0] * dims[1]];  // Real part
          } else {
            target[i + j * hsize + k * dims[1] * hsize] =
                src[i + (dims[1] - j) * dims[0] +
                    (dims[2] - k) * dims[0] * dims[1]];  // Real part
          }
        }
      }
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_iptr<1>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_iptr<2>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  int hsize = dims[0] / 2 + 1;
  for (int j = 0; j < dims[1]; ++j) {
    for (int i = 0; i < hsize; ++i) {
      if (j <= hsize) {
        target[i + j * hsize] = src[i + j * dims[0]];  // Real part
      } else {
        target[i + j * hsize] = src[i + (dims[1] - j) * dims[0]];  // Real part
      }
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_iptr<3>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  int hsize = dims[0] / 2 + 1;
  for (int k = 0; k < dims[2]; ++k) {
    for (int j = 0; j < dims[0]; ++j) {
      for (int i = 0; i < hsize; ++i) {
        if (k <= hsize) {
          if (j <= hsize) {
            target[i + j * dims[0] + k * dims[0] * dims[1]] =
                src[i + j * hsize + k * dims[1] * hsize];  // Real part
          } else {
            target[i + (dims[1] - j) * dims[0] + k * dims[0] * dims[1]] =
                src[i + j * hsize + k * dims[0] * hsize];  // Real part
          }
        } else {
          if (j <= hsize) {
            target[i + j * dims[0] + (dims[2] - k) * dims[1] * dims[0]] =
                src[i + j * hsize + k * dims[1] * hsize];  // Real part
          } else {
            target[i + (dims[10] - j) * dims[0] +
                   (dims[2] - k) * dims[1] * dims[0]] =
                src[i + j * hsize + k * dims[1] * hsize];  // Real part
          }
        }
      }
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_sthc.
template <>
inline void symphas::dft::arrange_fftw_sthc<1>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_sthc.
template <>
inline void symphas::dft::arrange_fftw_sthc<2>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = dn;

  for (iter_type j = 1; j < dims[1]; ++j) {
    iter_type n = dims[0] * j;
    for (iter_type i = 0; i < dn; ++i) {
      target[ii++] = src[n++];
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_sthc.
template <>
inline void symphas::dft::arrange_fftw_sthc<3>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = dn;

  for (iter_type j = 1; j < dims[1]; ++j) {
    iter_type n = dims[0] * j;
    for (iter_type i = 0; i < dn; ++i) {
      target[ii++] = src[n++];
    }
  }
  for (iter_type k = 1; k < dims[2]; ++k) {
    for (iter_type j = 0; j < dims[1]; ++j) {
      iter_type n = dims[0] * j + dims[0] * dims[1] * k;
      for (iter_type i = 0; i < dn; ++i) {
        target[ii++] = src[n++];
      }
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_hcts.
template <>
inline void symphas::dft::arrange_fftw_hcts<1>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;

  if (src != target) {
#pragma omp parallel for
    for (iter_type i = 0; i < dn; ++i) {
      target[i] = src[i];
    }
  }

#pragma omp parallel for
  for (iter_type i = dn; i < dims[0]; ++i) {
    target[i] = symphas::math::conj(src[2 * dn - 2 - i]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_hcts.
template <>
inline void symphas::dft::arrange_fftw_hcts<2>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = grid::length<2>(dims) - dims[0] + dn;
  iter_type n = dn * dims[1];

  /* Iterate over the source half.
   */
  for (iter_type j = 0; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dn; ++i) {
      target[--ii] = symphas::math::conj(src[--n]);
    }
    ii -= dims[0] - dn;
  }

  /* Iterate over the hermitian half.
   */
  n = dims[0] * (dims[1] - 1) + dn - 2;
  ii = dims[0] + dn;
  for (iter_type j = 1; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = symphas::math::conj(target[n--]);
    }
    n -= dn;
    ii += dn;
  }
  ii = dn;
  for (iter_type i = 0; i < dims[0] - dn; ++i) {
    target[ii++] = symphas::math::conj(target[n--]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_hcts.
template <>
inline void symphas::dft::arrange_fftw_hcts<3>(complex_t* src,
                                               complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = grid::length<3>(dims) - dims[0] + dn;
  iter_type n = dn * dims[1] * dims[2];

  for (iter_type k = 0; k < dims[2]; ++k) {
    for (iter_type j = 0; j < dims[1]; ++j) {
      for (iter_type i = 0; i < dn; ++i) {
        target[--ii] = symphas::math::conj(src[--n]);
      }
      ii -= dims[0] - dn;
    }
  }

  iter_type layer_offset = 0, layer_length = dims[0] * dims[1],
            last_layer_offset = layer_length * (dims[2] - 1);
  for (iter_type k = 1; k < dims[2]; ++k) {
    n = (last_layer_offset - layer_offset) + dims[0] * (dims[1] - 1) + dn - 2;
    ii = layer_length + layer_offset + dims[0] + dn;
    for (iter_type j = 1; j < dims[1]; ++j) {
      for (iter_type i = 0; i < dims[0] - dn; ++i) {
        target[ii++] = symphas::math::conj(target[n--]);
      }
      n -= dn;
      ii += dn;
    }
    ii = layer_length + layer_offset + dn;
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = symphas::math::conj(target[n--]);
    }
    layer_offset += layer_length;
  }
  n = dims[0] * (dims[1] - 1) + dn - 2;
  ii = dims[0] + dn;
  for (iter_type j = 1; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = symphas::math::conj(target[n--]);
    }
    n -= dn;
    ii += dn;
  }
  ii = dn;
  for (iter_type i = 0; i < dims[0] - dn; ++i) {
    target[ii++] = symphas::math::conj(target[n--]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_ipts.
template <>
inline void symphas::dft::arrange_fftw_ipts<1>(scalar_t* src, scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = 0;
  iter_type n = 0;

  for (iter_type i = 0; i < dims[0]; ++i) {
    target[ii++] = src[n++];
  }
}

//! Specialization based on symphas::dft::arrange_fftw_ipts.
template <>
inline void symphas::dft::arrange_fftw_ipts<2>(scalar_t* src, scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = 0;
  iter_type n = 0;
  iter_type c = 2 - dims[0] % 2;

  for (iter_type j = 0; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dims[0]; ++i) {
      target[ii++] = src[n++];
    }
    n += c;
  }
}

//! Specialization based on symphas::dft::arrange_fftw_ipts.
template <>
inline void symphas::dft::arrange_fftw_ipts<3>(scalar_t* src, scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = 0;
  iter_type n = 0;
  iter_type c = 2 - dims[0] % 2;

  for (iter_type k = 0; k < dims[2]; ++k) {
    for (iter_type j = 0; j < dims[1]; ++j) {
      for (iter_type i = 0; i < dims[0]; ++i) {
        target[ii++] = src[n++];
      }
      n += c;
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_stip<1>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = grid::length<1>(dims) - 1;  //!< The index in the source array.
  iter_type c = 2 - dims[0] % 2;  //!< The jump to the next row as according to
                                  //!< the FFTW arrangement.
  iter_type n = (dims[0] + c) - 1;  //!< Index in the transformed array.

  n -= c;

  for (iter_type i = 0; i < dims[0]; ++i) {
    target[n--] = src[ii--];
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_stip<2>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = grid::length<2>(dims) - 1;
  iter_type c = 2 - dims[0] % 2;
  iter_type n = (dims[0] + c) * dims[1] - 1;

  for (iter_type j = 0; j < dims[1]; ++j) {
    n -= c;
    for (iter_type i = 0; i < dims[0]; ++i) {
      target[n--] = src[ii--];
    }
  }
}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
inline void symphas::dft::arrange_fftw_stip<3>(const scalar_t* src,
                                               scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = grid::length<3>(dims) - 1;
  iter_type c = 2 - dims[0] % 2;
  iter_type n = (dims[0] + c) * dims[1] * dims[2] - 1;

  for (iter_type k = 0; k < dims[2]; ++k) {
    for (iter_type j = 0; j < dims[1]; ++j) {
      n -= c;
      for (iter_type i = 0; i < dims[0]; ++i) {
        target[n--] = src[ii--];
      }
    }
  }
}