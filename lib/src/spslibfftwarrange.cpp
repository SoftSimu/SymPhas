
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

#include "spslibfftwarrange.h"

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
void symphas::dft::arrange_fftw_rtip<1>(const scalar_t* src, scalar_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
void symphas::dft::arrange_fftw_rtip<2>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_rtip<3>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_iptr<1>(const scalar_t* src, scalar_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_stip.
template <>
void symphas::dft::arrange_fftw_iptr<2>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_iptr<3>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_sthc<1>(complex_t* src, complex_t* target,
                                        const len_type* dims) {}

//! Specialization based on symphas::dft::arrange_fftw_sthc.
template <>
void symphas::dft::arrange_fftw_sthc<2>(complex_t* src, complex_t* target,
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
void symphas::dft::arrange_fftw_sthc<3>(complex_t* src, complex_t* target,
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
void symphas::dft::arrange_fftw_hcts<1>(complex_t* src, complex_t* target,
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
    target[i] = std::conj(src[2 * dn - 2 - i]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_hcts.
template <>
void symphas::dft::arrange_fftw_hcts<2>(complex_t* src, complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = grid::length<2>(dims) - dims[0] + dn;
  iter_type n = dn * dims[1];

  /* Iterate over the source half.
   */
  for (iter_type j = 0; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dn; ++i) {
      target[--ii] = std::conj(src[--n]);
    }
    ii -= dims[0] - dn;
  }

  /* Iterate over the hermitian half.
   */
  n = dims[0] * (dims[1] - 1) + dn - 2;
  ii = dims[0] + dn;
  for (iter_type j = 1; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = std::conj(target[n--]);
    }
    n -= dn;
    ii += dn;
  }
  ii = dn;
  for (iter_type i = 0; i < dims[0] - dn; ++i) {
    target[ii++] = std::conj(target[n--]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_hcts.
template <>
void symphas::dft::arrange_fftw_hcts<3>(complex_t* src, complex_t* target,
                                        const len_type* dims) {
  iter_type dn = dims[0] / 2 + 1;
  iter_type ii = grid::length<3>(dims) - dims[0] + dn;
  iter_type n = dn * dims[1] * dims[2];

  for (iter_type k = 0; k < dims[2]; ++k) {
    for (iter_type j = 0; j < dims[1]; ++j) {
      for (iter_type i = 0; i < dn; ++i) {
        target[--ii] = std::conj(src[--n]);
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
        target[ii++] = std::conj(target[n--]);
      }
      n -= dn;
      ii += dn;
    }
    ii = layer_length + layer_offset + dn;
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = std::conj(target[n--]);
    }
    layer_offset += layer_length;
  }
  n = dims[0] * (dims[1] - 1) + dn - 2;
  ii = dims[0] + dn;
  for (iter_type j = 1; j < dims[1]; ++j) {
    for (iter_type i = 0; i < dims[0] - dn; ++i) {
      target[ii++] = std::conj(target[n--]);
    }
    n -= dn;
    ii += dn;
  }
  ii = dn;
  for (iter_type i = 0; i < dims[0] - dn; ++i) {
    target[ii++] = std::conj(target[n--]);
  }
}

//! Specialization based on symphas::dft::arrange_fftw_ipts.
template <>
void symphas::dft::arrange_fftw_ipts<1>(scalar_t* src, scalar_t* target,
                                        const len_type* dims) {
  iter_type ii = 0;
  iter_type n = 0;

  for (iter_type i = 0; i < dims[0]; ++i) {
    target[ii++] = src[n++];
  }
}

//! Specialization based on symphas::dft::arrange_fftw_ipts.
template <>
void symphas::dft::arrange_fftw_ipts<2>(scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_ipts<3>(scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_stip<1>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_stip<2>(const scalar_t* src, scalar_t* target,
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
void symphas::dft::arrange_fftw_stip<3>(const scalar_t* src, scalar_t* target,
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
//
////! Specialization of iterate_complex_dup.
// template struct symphas::dft::iterate_complex_dup<scalar_t, 2>;
//
////! Specialization of iterate_complex_dup.
// template struct symphas::dft::iterate_complex_dup<scalar_t, 3>;
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_rtip<1>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_rtip<2>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_rtip<3>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_iptr<1>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_iptr<2>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_iptr<3>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_sthc.
// template void symphas::dft::arrange_fftw_sthc<1>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_sthc.
// template void symphas::dft::arrange_fftw_sthc<2>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
////! Specialization based on symphas::dft::arrange_fftw_sthc.
// template void symphas::dft::arrange_fftw_sthc<3>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_hcts.
// template void symphas::dft::arrange_fftw_hcts<1>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_hcts.
// template void symphas::dft::arrange_fftw_hcts<2>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_hcts.
// template void symphas::dft::arrange_fftw_hcts<3>(complex_t* src,
//                                                  complex_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_ipts.
// template void symphas::dft::arrange_fftw_ipts<1>(scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
////! Specialization based on symphas::dft::arrange_fftw_ipts.
// template void symphas::dft::arrange_fftw_ipts<2>(scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_ipts.
// template void symphas::dft::arrange_fftw_ipts<3>(scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_stip<1>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
//
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_stip<2>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
////! Specialization based on symphas::dft::arrange_fftw_stip.
// template void symphas::dft::arrange_fftw_stip<3>(const scalar_t* src,
//                                                  scalar_t* target,
//                                                  const len_type* dims);
