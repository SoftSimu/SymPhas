
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
