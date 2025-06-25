
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
 * PURPOSE: Defines the vector class used as the phase field order parameter
 * type.
 *
 * ***************************************************************************
 */

#pragma once

#include <cmath>
#include <complex>

#include "configured-defs.h"

#ifdef USING_CUDA

#define CHECK_CUDA_ERROR(call)                                              \
  {                                                                         \
    cudaError_t err = call;                                                 \
    if (err != cudaSuccess) {                                               \
      printf("CUDA Error: %s @ %s:%d\n", cudaGetErrorString(err), __FILE__, \
             __LINE__);                                                     \
      exit(err);                                                            \
    }                                                                       \
  }

#endif

#ifdef USING_CUDA
#include <cuda_runtime.h>
#else
#undef __host__
#undef __device__
#define __host__
#define __device__
#endif

//! A D-dimensional vector value supporting various operations.
/*!
 * This is a basic datatype in SymPhas, used primarily to represent values of
 * an order parameter with multiple components, such as for a convective flow.
 */
template <typename T, size_t D>
struct VectorValue {
  //! Data member.
  /*!
   * An array of D elements which is the raw data of VectorValue.
   */
  T v[D]{};

  //! Data access operator.
  /*!
   * Access operator for the data memeber `v`. The position in the array is
   * passed. The member access only returns a modifiable reference to the member
   * data array element. In the context of position information on a grid, the
   * first index is typically used for the horizontal component, \f$x\f$.
   *
   * \param i The position in the member data array to access.
   */
  __host__ __device__ T& operator[](int i) { return v[i]; }

  //! Data access operator.
  /*!
   * Const access operator for the data member `v`. See
   * VectorValue::operator[](int). The access is an unmodifiable reference to
   * the member data array element.
   *
   * \param i The position in the data array to access.
   */
  __host__ __device__ const T& operator[](int i) const { return v[i]; }

  //! Compound assignment operator.
  /*!
   * Compound assignment operator for addition.
   *
   * Operator is called by the left hand instance, taking data from the right
   * hand VectorValue to perform point-wise addition between left and right
   * sides and assigning the result to the left hand side.
   *
   * \param rhs The VectorValue instance from where data is added point-wise
   * to the vector of the VectorValue on the left of the operator.
   */
  __host__ __device__ VectorValue<T, D>& operator+=(
      const VectorValue<T, D>& rhs) {
    for (int i = 0; i < D; ++i) {
      v[i] += rhs.v[i];
    }
    return *this;  // return the result by reference
  }

  //! Addition operator.
  /*!
   * Addition operator, which returns a new VectorValue instance with its data
   * equal to the result of point-wise addition of the data in the given
   * operands.
   *
   * \param lhs The `VectorValue` instance on the left hand side of the addition
   * operator. \param rhs The `VectorValue` instance on the right hand side of
   * the addition operator.
   */
  __host__ __device__ friend VectorValue<T, D> operator+(
      VectorValue<T, D> lhs, const VectorValue<T, D>& rhs) {
    lhs += rhs;
    return lhs;
  }

  //! Negative compound assignment operator.
  /*!
   * Compound assignment operator for subtraction.
   *
   * Assignment operator to subtract the values from the calling instance.
   * Operator is called by the left hand instance, taking data from the right
   * hand VectorValue to perform point-wise subtraction between left and right
   * sides, assigning the result to the left hand side.
   *
   * \param rhs The VectorValue instance where data is subtracted point-wise
   * from the vector of the VectorValue on the left of the operator.
   */
  __host__ __device__ VectorValue<T, D>& operator-=(
      const VectorValue<T, D>& rhs) {
    for (int i = 0; i < D; ++i) {
      v[i] -= rhs.v[i];
    }
    return *this;
  }

  //! Subtraction operator.
  /*!
   * Subtraction operator, which returns a new VectorValue instance with its
   * data equal to the result of point-wise addition of the data in the given
   * operands.
   *
   * \param lhs The `VectorValue` instance on the left hand side of the
   * subtraction operator. \param rhs The `VectorValue` instance on the right
   * hand side of the subtraction operator.
   */
  __host__ __device__ friend VectorValue<T, D> operator-(
      VectorValue<T, D> lhs, const VectorValue<T, D>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  //! Unary negative operator.
  /*!
   * Operator applied to the left side of an instance that returns a new
   * VectorValue instance where all components of the data are determined by
   * applying the unary negative operator to the components of the calling
   * instance.
   */
  __host__ __device__ VectorValue<T, D> operator-() const {
    VectorValue<T, D> neg = *this;
    for (int i = 0; i < D; ++i) {
      neg.v[i] = -neg.v[i];
    }
    return neg;
  }

  //! Dot product operator.
  /*!
   * Operator for computing the dot product of two VectorValue instances. A
   * vector multiplied with another vector of the same dimensions is assumed to
   * be the dot product. This would be equivalent of multiplying a row vector by
   * a column vector.
   *
   * Assignment operator to multiply all the components of the data from the
   * calling instance with the scalar quantity on the left hand side.
   *
   * \param lhs First VectorValue instance used to compute the dot product.
   * \param rhs Second VectorValue instance used to compute the dot product.
   */
  template <typename S>
  __host__ __device__ friend auto dot(VectorValue<T, D> const& lhs,
                                      VectorValue<S, D> const& rhs) {
    auto sum = lhs[0] * rhs[0];
    for (int i = 1; i < D; ++i) {
      sum += lhs[i] * rhs[i];
    }
    return sum;
  }
};

template <typename T>
VectorValue(T) -> VectorValue<T, 1>;
template <typename T>
VectorValue(T, T) -> VectorValue<T, 2>;
template <typename T>
VectorValue(T, T, T) -> VectorValue<T, 3>;

//
// extern template struct VectorValue<double, 1>;
// extern template struct VectorValue<double, 2>;
// extern template struct VectorValue<double, 3>;
// extern template struct VectorValue<std::complex<double>, 1>;
// extern template struct VectorValue<std::complex<double>, 2>;
// extern template struct VectorValue<std::complex<double>, 3>;
// extern template struct VectorValue<float, 1>;
// extern template struct VectorValue<float, 2>;
// extern template struct VectorValue<float, 3>;
// extern template struct VectorValue<std::complex<float>, 1>;
// extern template struct VectorValue<std::complex<float>, 2>;
// extern template struct VectorValue<std::complex<float>, 3>;

//! A matrix with `N` rows and `M` columns.
template <typename T, size_t N, size_t M>
using MatrixValue = VectorValue<VectorValue<T, M>, N>;
template <typename T, size_t D>
using RowVectorValue = VectorValue<VectorValue<T, D>, 1>;

/* default implementations of distance and abs for VectorValue
 */

template <typename T, size_t D>
__host__ __device__ T abs(VectorValue<T, D> const& c) {
  double sum = 0;
  for (int i = 0; i < D; ++i) {
    sum += c.v[i] * c.v[i];
  }
  using std::sqrt;
  return sqrt(sum);
}

template <typename T, size_t D>
__host__ __device__ T distance(VectorValue<T, D> const& c1,
                               VectorValue<T, D> const& c2) {
  using std::abs;
  return abs(c2 - c1);
}

//
// template<>
// inline double distance(VectorValue<double, 1> const& c1, VectorValue<double,
// 1> const& c2); template<> inline double distance(VectorValue<double, 2>
// const& c1, VectorValue<double, 2> const& c2); template<> inline double
// distance(VectorValue<double, 3> const& c1, VectorValue<double, 3> const& c2);
// template<>
// inline double abs(VectorValue<double, 1> const& c);
// template<>
// inline double abs(VectorValue<double, 2> const& c);
// template<>
// inline double abs(VectorValue<double, 3> const& c);
//
//
// template<>
// inline float distance(VectorValue<float, 1> const& c1, VectorValue<float, 1>
// const& c2); template<> inline float distance(VectorValue<float, 2> const& c1,
// VectorValue<float, 2> const& c2); template<> inline float
// distance(VectorValue<float, 3> const& c1, VectorValue<float, 3> const& c2);
// template<>
// inline float abs(VectorValue<float, 1> const& c);
// template<>
// inline float abs(VectorValue<float, 2> const& c);
// template<>
// inline float abs(VectorValue<float, 3> const& c);

//! Scalar multiplication operator.
/*!
 * Operator for multiplication with a scalar quantity.
 *
 * Assignment operator to multiply all the components of the data from the
 calling
 * instance with the scalar quantity on the right hand side.
 *
 * \param lhs The VectorValue instance from which data is multiplied point-wise
 by the scalar
 * quantity on the right of the operator.
 * \param rhs The scalar value used which multiplies the data in the VectorValue
 instance.

 * \tparam V The template parameter for the scalar type. SFINAE strategy is used
 so that
 * the operator will only be compiled for types `V` which satisfy the condition
 that the
 * resultant type of multiplication between instances of types `T` and `V` is
 implicitly
 * convertible to `T`.
 */
// template<typename T, size_t D>
// auto operator*(VectorValue<T, D> const& lhs, double rhs)
//{
//	VectorValue<mul_result_t<T, double>, D> lhs0{ lhs };
//	for (int i = 0; i < D; ++i)
//	{
//		lhs0[i] *= rhs;
//	}
//	lhs *= rhs;
//	return lhs;
// }

//! Scalar multiplication operator.
/*!
 * Operator for multiplication with a scalar quantity with scalar quantity on
 left hand side.
 *
 * Assignment operator to multiply all the components of the data from the
 calling
 * instance with the scalar quantity on the left hand side.
 *
 * \param lhs The scalar value used which multiplies the data in the VectorValue
 instance.
 * \param rhs The VectorValue instance from which data is multiplied point-wise
 by the scalar
 * quantity on the left of the operator.

 * \tparam V The template parameter for the scalar type. SFINAE strategy is used
 so that
 * the operator will only be compiled for types `V` which satisfy the condition
 that the
 * resultant type of multiplication between instances of types `T` and `V` is
 implicitly
 * convertible to `T`.
 */
// friend auto operator*(T const& lhs, VectorValue<T, D> const& rhs)
//{
//	return rhs * lhs;
// }
