
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

#include "spslibfftw.h"

#ifdef USING_FFTW

#include <fftw3.h>

fftw_plan symphas::dft::new_fftw_plan<1, scalar_t, scalar_t>::operator()(
    scalar_t* src, scalar_t* target, const len_type* dims, bool estimate,
    bool backward) {
  return fftw_plan_r2r_1d(dims[0], src, target,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<2, scalar_t, scalar_t>::operator()(
    scalar_t* src, scalar_t* target, const len_type* dims, bool estimate,
    bool backward) {
  return fftw_plan_r2r_2d(dims[1], dims[0], src, target,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<3, scalar_t, scalar_t>::operator()(
    scalar_t* src, scalar_t* target, const len_type* dims, bool estimate,
    bool backward) {
  return fftw_plan_r2r_3d(dims[2], dims[1], dims[0], src, target,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (!backward) ? FFTW_R2HC : FFTW_HC2R,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}
fftw_plan symphas::dft::new_fftw_plan<1, complex_t, complex_t>::operator()(
    fftw_complex* src, fftw_complex* target, const len_type* dims,
    bool estimate, bool backward) {
  return fftw_plan_dft_1d(dims[0], reinterpret_cast<fftw_complex*>(src),
                          reinterpret_cast<fftw_complex*>(target),
                          (!backward) ? FFTW_FORWARD : FFTW_BACKWARD,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<2, complex_t, complex_t>::operator()(
    fftw_complex* src, fftw_complex* target, const len_type* dims,
    bool estimate, bool backward) {
  return fftw_plan_dft_2d(dims[1], dims[0],
                          reinterpret_cast<fftw_complex*>(src),
                          reinterpret_cast<fftw_complex*>(target),
                          (!backward) ? FFTW_FORWARD : FFTW_BACKWARD,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<3, complex_t, complex_t>::operator()(
    fftw_complex* src, fftw_complex* target, const len_type* dims,
    bool estimate, bool backward) {
  return fftw_plan_dft_3d(dims[2], dims[1], dims[0],
                          reinterpret_cast<fftw_complex*>(src),
                          reinterpret_cast<fftw_complex*>(target),
                          (!backward) ? FFTW_FORWARD : FFTW_BACKWARD,
                          (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<1, scalar_t, complex_t>::operator()(
    scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_r2c_1d(dims[0], reinterpret_cast<scalar_t*>(src),
                              reinterpret_cast<fftw_complex*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<2, scalar_t, complex_t>::operator()(
    scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_r2c_2d(dims[1], dims[0],
                              reinterpret_cast<scalar_t*>(src),
                              reinterpret_cast<fftw_complex*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<3, scalar_t, complex_t>::operator()(
    scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_r2c_3d(dims[2], dims[1], dims[0],
                              reinterpret_cast<scalar_t*>(src),
                              reinterpret_cast<fftw_complex*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<1, complex_t, scalar_t>::operator()(
    fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_c2r_1d(dims[0], reinterpret_cast<fftw_complex*>(src),
                              reinterpret_cast<scalar_t*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<2, complex_t, scalar_t>::operator()(
    fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_c2r_2d(dims[1], dims[0],
                              reinterpret_cast<fftw_complex*>(src),
                              reinterpret_cast<scalar_t*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

fftw_plan symphas::dft::new_fftw_plan<3, complex_t, scalar_t>::operator()(
    fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate) {
  return fftw_plan_dft_c2r_3d(dims[2], dims[1], dims[0],
                              reinterpret_cast<fftw_complex*>(src),
                              reinterpret_cast<scalar_t*>(target),
                              (estimate) ? FFTW_ESTIMATE : FFTW_MEASURE);
}

void symphas::dft::fftw_execute(fftw_plan p) { ::fftw_execute(p); }

void symphas::dft::fftw_destroy_plan(fftw_plan p) {
  if (p != NULL) {
    ::fftw_destroy_plan(p);
  }
}

fftw_complex* symphas::dft::fftw_alloc_complex(size_t n) {
  return ::fftw_alloc_complex(n);
}

scalar_t* symphas::dft::fftw_alloc_real(size_t n) {
  return (scalar_t*)fftw_malloc(sizeof(scalar_t) * n);
}

void symphas::dft::fftw_free(fftw_complex*& arr) { ::fftw_free(arr); }

#endif
