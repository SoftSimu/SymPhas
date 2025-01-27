
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
#include "spslibfftwarrange.h"

struct fftw_plan_s;
typedef double fftw_complex[2];
typedef fftw_plan_s* fftw_plan;

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

namespace symphas::dft {

#ifdef USING_FFTW

/* Subset of FFTW functions used in SymPhas. FFTW is not exposed to
 * the whole program, but a limited number of functions are specified
 * and the rest are adapted to the specific workings of SymPhas.
 */

//! Executes the FFTW plan.
void fftw_execute(fftw_plan p);
//! Destroys the FFTW plan.
void fftw_destroy_plan(fftw_plan p);
//! Allocates memory aligned complex array of the prescribed size.
fftw_complex* fftw_alloc_complex(size_t);
//! Allocates memory aligned real-valued array of the prescribed size.
scalar_t* fftw_alloc_real(size_t);
//! Frees the memory associated with an aligned FFTW array.
void fftw_free(fftw_complex*&);

//! Creates a new FFTW plan from the given types and dimension.
/*!
 * Using the given dimension, construct a new FFTW plan based on
 * transforming a source array of type `T_src` to a transform target
 * of type `T_target`. The template parameters are used to specialize the
 * plan type which is selected.
 *
 * \tparam D The dimension of the system that is transformed.
 * \tparam T_src The source system value type.
 * \tparam T_target The type of the Fourier transform output.
 */
template <size_t D, typename T_src, typename T_target>
struct new_fftw_plan;

/* specialized plans, in particular implementing real to complex and vice versa
 * note that these plans (when only complex values are involved) are purely
 * forward, hence if a "backward" fourier transform needs to be performed, the
 * output grid needs to be passed as the second argument
 */

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<1, scalar_t, scalar_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(scalar_t* src, scalar_t* target, const len_type* dims,
                       bool estimate = false, bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<1, scalar_t, scalar_t>{}(
        reinterpret_cast<scalar_t*>(src), reinterpret_cast<scalar_t*>(target),
        dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<2, scalar_t, scalar_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(scalar_t* src, scalar_t* target, const len_type* dims,
                       bool estimate = false, bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<2, scalar_t, scalar_t>{}(
        reinterpret_cast<scalar_t*>(src), reinterpret_cast<scalar_t*>(target),
        dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<3, scalar_t, scalar_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(scalar_t* src, scalar_t* target, const len_type* dims,
                       bool estimate = false, bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<3, scalar_t, scalar_t>{}(
        reinterpret_cast<scalar_t*>(src), reinterpret_cast<scalar_t*>(target),
        dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<1, complex_t, complex_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(fftw_complex* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false,
                       bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<1, complex_t, complex_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<2, complex_t, complex_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(fftw_complex* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false,
                       bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<2, complex_t, complex_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<3, complex_t, complex_t> {
  //! Creates a new complex to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to complex plan is considered to be
   * forward by default, unless \p backward is provided equal to true.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  fftw_plan operator()(fftw_complex* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false,
                       bool backward = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   * \param backward Whether or not the plan refers to taking a forward
   * Fourier transform, or backward.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false, bool backward = false) {
    return new_fftw_plan<3, complex_t, complex_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<1, scalar_t, complex_t> {
  //! Creates a new real to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the real to complex plan is always forward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(scalar_t* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<1, scalar_t, complex_t>{}(
        reinterpret_cast<scalar_t*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<2, scalar_t, complex_t> {
  //! Creates a new real to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the real to complex plan is always forward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(scalar_t* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<2, scalar_t, complex_t>{}(
        reinterpret_cast<scalar_t*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<3, scalar_t, complex_t> {
  //! Creates a new real to complex FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the real to complex plan is always forward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(scalar_t* src, fftw_complex* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<3, scalar_t, complex_t>{}(
        reinterpret_cast<scalar_t*>(src),
        reinterpret_cast<fftw_complex*>(target), dims, estimate);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<1, complex_t, scalar_t> {
  //! Creates a new complex to real FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to real plan is always backward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(fftw_complex* src, scalar_t* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<1, complex_t, scalar_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<scalar_t*>(target), dims, estimate);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<2, complex_t, scalar_t> {
  //! Creates a new complex to real FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to real plan is always backward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(fftw_complex* src, scalar_t* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<2, complex_t, scalar_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<scalar_t*>(target), dims, estimate);
  }
};

//! Specialization based on symphas::dft::new_fftw_plan.
template <>
struct new_fftw_plan<3, complex_t, scalar_t> {
  //! Creates a new complex to real FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   * The direction of the complex to real plan is always backward.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  fftw_plan operator()(fftw_complex* src, scalar_t* target,
                       const len_type* dims, bool estimate = false);
  //! Forwards arbitrary arguments to create a specific FFTW plan.
  /*!
   * Creates a plan based on transforming \p src to \p target.
   *
   * \param src The source array.
   * \param target The target array.
   * \param dims The dimension of the source array.
   * \param estimate Whether or not to estimate the plan. If the plan
   * is not estimated, the source and target arrays will be modified.
   */
  template <typename T1, typename T2>
  fftw_plan operator()(T1* src, T2* target, const len_type* dims,
                       bool estimate = false) {
    return new_fftw_plan<3, complex_t, scalar_t>{}(
        reinterpret_cast<fftw_complex*>(src),
        reinterpret_cast<scalar_t*>(target), dims, estimate);
  }
};

#endif

// **************************************************************************************

//! Computes the number of array elements for the given type.
/*!
 * Using the FFTW algorithm means that the number of elements in the
 * complex array may be shortened since the transform will be duplicated
 * by way of Hermitian conjugate.
 *
 * When the given type is complex, the true length of the data array
 * represented by the given dimensions is returned. When the given type
 * is real, then only the length of the non-duplicated half is returned.
 *
 * \param dims The dimensions of the grid which is transformed.
 *
 * \tparam T The type which the length is adjusted to. A real type will
 * return a length half of the total dimensions.
 * \tparam D the dimension of the corresponding dimension.
 */
template <typename T, size_t D>
len_type length(len_type const* dims) {
  return 1;
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<scalar_t, 1>(len_type const* dims) {
  return (dims[0] / 2 + 1);
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<scalar_t, 2>(len_type const* dims) {
  return (dims[0] / 2 + 1) * dims[1];
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<scalar_t, 3>(len_type const* dims) {
  return (dims[0] / 2 + 1) * dims[1] * dims[2];
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<complex_t, 1>(len_type const* dims) {
  return grid::length<1>(dims);
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<complex_t, 2>(len_type const* dims) {
  return grid::length<2>(dims);
}

//! Specialization of symphas::dft::length().
template <>
inline len_type length<complex_t, 3>(len_type const* dims) {
  return grid::length<3>(dims);
}

#ifdef USING_FFTW

//! Executes the given plan to compute the Fourier transform.
/*!
 * An overload that takes an instance of the FFTW plan type and executes it
 * to compute the Fourier transform.
 *
 * \param p An instance of the FFTW plan type.
 */
inline void dft(fftw_plan p) { symphas::dft::fftw_execute(p); }

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given real-valued system and store
 * the result in the given complex array. The given data is a sequential
 * array representing a `D`-dimensional grid with dimensions \p dims. The
 * input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. Since the input is real-valued, the result then has to be
 * rearranged to represent the whole system that will include duplicated
 * values. The plan will always do a forward type transformation. The plan
 * is then destroyed.
 *
 * \param[in] in The input to the Fourier transform, as real-valued data.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const scalar_t* in, complex_t* out, const len_type* dims, bool) {
  fftw_plan p = new_fftw_plan<D, scalar_t, complex_t>{}(
      const_cast<scalar_t*>(in), out, dims, true);
  symphas::dft::fftw_execute(p);
  arrange_fftw_hcts<D>(out, out, dims);
  symphas::dft::fftw_destroy_plan(p);
}

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given complex-valued system and
 * store the result in the given complex array. The given data is a
 * sequential array representing a `D`-dimensional grid with dimensions
 * \p dims. The input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. The plan is then destroyed.
 *
 * \param[in] in The input to the Fourier transform.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const complex_t* in, complex_t* out, const len_type* dims,
         bool backward = false) {
  fftw_plan p = new_fftw_plan<D, complex_t, complex_t>{}(
      const_cast<complex_t*>(in), out, dims, true, backward);
  symphas::dft::fftw_execute(p);
  symphas::dft::fftw_destroy_plan(p);
}

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given complex-valued system and
 * store the result in the given real-valued array. The given data is a
 * sequential array representing a `D`-dimensional grid with dimensions
 * \p dims. The input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. Since the output is real-valued, the input is first rearranged
 * because this type of transformation assumes that the input array has
 * values which repeat in the two halves. Once the plan is executed, the
 * input is restored to the original arrangement. Therefore, this plan
 * always assumes a backward type transformation. That is, it is the
 * opposite of symphas::dft::dft(const scalar_t*, complex_t*, const len_type*).
 * The plan is subsequently destroyed.
 *
 * **NOTE**: If the input array does not have repeated values between each
 * half, then this routine will destroy the original data.

 * \param[in] in The input to the Fourier transform, as complex-valued data.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const complex_t* in, scalar_t* out, const len_type* dims, bool) {
  complex_t* hf = const_cast<complex_t*>(in);
  arrange_fftw_sthc<D>(hf, hf, dims);
  fftw_plan p = new_fftw_plan<D, complex_t, scalar_t>{}(hf, out, dims, true);
  symphas::dft::fftw_execute(p);
  symphas::dft::fftw_destroy_plan(p);
  arrange_fftw_hcts<D>(hf, hf, dims);
}

#else

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given real-valued system and store
 * the result in the given complex array. The given data is a sequential
 * array representing a `D`-dimensional grid with dimensions \p dims. The
 * input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. Since the input is real-valued, the result then has to be
 * rearranged to represent the whole system that will include duplicated
 * values. The plan will always do a forward type transformation. The plan
 * is then destroyed.
 *
 * \param[in] in The input to the Fourier transform, as real-valued data.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const scalar_t* in, complex_t* out, const len_type* dims, bool) {
  if constexpr (D == 1) {
    long_dft(const_cast<scalar_t*>(in), out, dims[0]);
  } else if constexpr (D == 2) {
    long_dft(const_cast<scalar_t*>(in), out, dims[0], dims[1]);
  } else if constexpr (D == 3) {
    long_dft(const_cast<scalar_t*>(in), out, dims[0], dims[1], dims[2]);
  }
}

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given complex-valued system and
 * store the result in the given complex array. The given data is a
 * sequential array representing a `D`-dimensional grid with dimensions
 * \p dims. The input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. The plan is then destroyed.
 *
 * \param[in] in The input to the Fourier transform.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const complex_t* in, complex_t* out, const len_type* dims,
         bool backward = false) {
  if constexpr (D == 1) {
    long_dft(const_cast<complex_t*>(in), out, dims[0], backward);
  } else if constexpr (D == 2) {
    long_dft(const_cast<complex_t*>(in), out, dims[0], dims[1], backward);
  } else if constexpr (D == 3) {
    long_dft(const_cast<complex_t*>(in), out, dims[0], dims[1], dims[2],
             backward);
  }
}

//! Compute the Fouier transform of a `D`-dimensional system.
/*!
 * Compute the Fourier transform of the given complex-valued system and
 * store the result in the given real-valued array. The given data is a
 * sequential array representing a `D`-dimensional grid with dimensions
 * \p dims. The input and output both have the same dimensions.
 *
 * The FFTW routine is used by creating an estimated plan, and then
 * executed. Since the output is real-valued, the input is first rearranged
 * because this type of transformation assumes that the input array has
 * values which repeat in the two halves. Once the plan is executed, the
 * input is restored to the original arrangement. Therefore, this plan
 * always assumes a backward type transformation. That is, it is the
 * opposite of symphas::dft::dft(const scalar_t*, complex_t*, const len_type*).
 * The plan is subsequently destroyed.
 *
 * **NOTE**: If the input array does not have repeated values between each
 * half, then this routine will destroy the original data.

 * \param[in] in The input to the Fourier transform, as complex-valued data.
 * \param[out] out The output array of the Fourier transform, where the
 * result is written to. This has dimensions equal to the input.
 * \param dims The dimensions of the input array of length `D`.
 *
 * \tparam D The dimension of the input system.
 */
template <size_t D>
void dft(const complex_t* in, scalar_t* out, const len_type* dims, bool) {
  if constexpr (D == 1) {
    long_dft(out, const_cast<complex_t*>(in), dims[0], true);
  } else if constexpr (D == 2) {
    long_dft(out, const_cast<complex_t*>(in), dims[0], dims[1], true);
  } else if constexpr (D == 3) {
    long_dft(out, const_cast<complex_t*>(in), dims[0], dims[1], dims[2], true);
  }
}

#endif

//! An arbitrary Fourier transform chosen based on the argument types.
/*!
 * The type of transform is chosen based on the argument types, which
 * must be compatible with the ::complex_t and ::scalar_t types.
 *
 * \param[in] in The data which is transformed.
 * \param[out] out The result of the Fourier transform of \p in.
 * \param dims The dimensions of the input array.
 *
 * \tparam D The dimension of the system.
 */
template <size_t D, typename T, typename C>
void dft(T&& in, C&& out, std::initializer_list<len_type> dims,
         bool backward = false) {
  dft<D>(std::forward<T>(in), std::forward<C>(out), dims.begin(), backward);
}

//! An arbitrary Fourier transform chosen based on the argument types.
/*!
 * The type of transform is chosen based on the argument types, which
 * must be compatible with the ::complex_t and ::scalar_t types. This
 * transform takes the input system to be a 1-dimensional system.
 *
 * \param[in] in The data which is transformed.
 * \param[out] out The result of the Fourier transform of \p in.
 * \param L The length of the input array.
 */
template <typename T, typename C>
void dft(T&& in, C&& out, len_type L, bool backward = false) {
  dft<1>(std::forward<T>(in), std::forward<C>(out), &L, backward);
}

//! An arbitrary Fourier transform chosen based on the argument types.
/*!
 * The type of transform is chosen based on the argument types, which
 * must be compatible with the ::complex_t and ::scalar_t types. This
 * transform takes the input system to be a 2-dimensional system.
 *
 * \param[in] in The data which is transformed.
 * \param[out] out The result of the Fourier transform of \p in.
 * \param L The length of the \f$x\f$-axis.
 * \param M The length of the \f$y\f$-axis.
 */
template <typename T, typename C>
void dft(T&& in, C&& out, len_type L, len_type M, bool backward = false) {
  len_type dims[2]{L, M};
  dft<2>(std::forward<T>(in), std::forward<C>(out), dims, backward);
}

//! An arbitrary Fourier transform chosen based on the argument types.
/*!
 * The type of transform is chosen based on the argument types, which
 * must be compatible with the ::complex_t and ::scalar_t types. This
 * transform takes the input system to be a 3-dimensional system.
 *
 * \param[in] in The data which is transformed.
 * \param[out] out The result of the Fourier transform of \p in.
 * \param L The length of the \f$x\f$-axis.
 * \param M The length of the \f$y\f$-axis.
 * \param N The length of the \f$z\f$-axis.
 */
template <typename T, typename C>
void dft(T&& in, C&& out, len_type L, len_type M, len_type N,
         bool backward = false) {
  len_type dims[3]{L, M, N};
  dft<3>(std::forward<T>(in), std::forward<C>(out), dims, backward);
}
}  // namespace symphas::dft

namespace symphas::math {

#ifdef USING_FFTW

inline scalar_t real(fftw_complex const& value) { return value[0]; }

inline scalar_t imag(fftw_complex const& value) { return value[1]; }

inline scalar_t& real(fftw_complex& value) { return value[0]; }

inline scalar_t& imag(fftw_complex& value) { return value[1]; }
#endif

inline scalar_t real(complex_t const& value) { return value.real(); }

inline scalar_t imag(complex_t const& value) { return value.imag(); }

inline scalar_t& real(complex_t& value) {
  return reinterpret_cast<scalar_t(&)[2]>(value)[0];
}

inline scalar_t& imag(complex_t& value) {
  return reinterpret_cast<scalar_t(&)[2]>(value)[1];
}
}  // namespace symphas::math
