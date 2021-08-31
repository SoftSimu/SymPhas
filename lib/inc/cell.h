
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


//! Flag for implicit conversion of product result.
/*!
 * Flag indicating if multiplying objects of two types, `T` and `V` can be 
 * implicitly cast to the first type `T`.
 */
template<typename T, typename V>
constexpr bool is_mult_convertible = std::is_convertible<
	T, decltype(std::declval<T>() * std::declval<V>())>::value;


//! Used to conditionally compile functions satisfying scalar multiplication.
/*!
 * Operators which involve a type `V` that can be multiplied by the underlying
 * type of the ValueVector `T` and then be implicitly convertible to `T` will
 * be compiled for that type `V`. Otherwise the template function will not
 * compile. This strategy is based on SFINAE. This class evaluates whether
 * this condition holds true, and stores the result as member `value`.
 *
 * \tparam T The underlying type of VectorValue, checked in multiplication with
 * `V`.
 * \tparam V The type to multiply with `T` so it can be checked if the result 
 * is a type that is implicitly convertible back to `T`.
 * \tparam else_type The type which should be returned by the function if the 
 * condition is satisfied.
 */
//! This is `true` if the multiplication converts to `T`, else `false`.
template<typename T, typename V, typename else_type>
constexpr bool is_scalar_mult =
	(std::is_same<V, typename std::remove_reference<else_type>::type>::value) ? false : is_mult_convertible<T, V>;

//! Type alias for enabling multiplication on scalar.
/*!
 * Aliases the `enable_if` construct to use SFINAE on member operators so that
 * the given type `V` can be such that the product of a type `T` and type `V`
 * is implicitly convertible back to type `T`.
 */
template<typename T, typename V, typename else_type>
using mult_enable_if = std::enable_if<is_scalar_mult<T, V, else_type>, else_type>;




//! A D-dimensional vector value supporting various operations.
/*!
 * This is a basic datatype in SymPhas, used primarily to represent values of
 * an order parameter with multiple components, such as for a convective flow.
 */
template<typename T, size_t D>
struct VectorValue
{
	//! Data member.
	/*!
	 * An array of D elements which is the raw data of VectorValue.
	 */
 	T v[D]{ 0 };

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
	VectorValue<T, D>& operator+=(const VectorValue<T, D>& rhs)
	{
		for (int i = 0; i < D; ++i)
		{
			v[i] += rhs.v[i];
		}
		return *this; // return the result by reference
	}

	//! Addition operator.
	/*!
	 * Addition operator, which returns a new VectorValue instance with its data equal
	 * to the result of point-wise addition of the data in the given operands.
	 *
	 * \param lhs The `VectorValue` instance on the left hand side of the addition operator. 
	 * \param rhs The `VectorValue` instance on the right hand side of the addition operator.
	 */
	friend VectorValue<T, D> operator+(VectorValue<T, D> lhs, const VectorValue<T, D>& rhs)
	{
		lhs += rhs;
		return lhs;
	}

	//! Negative compound assignment operator.
	/*!
	 * Compound assignment operator for subtraction.
	 * 
	 * Assignment operator to subtract the values from the calling instance. Operator 
	 * is called by the left hand instance, taking data from the right hand VectorValue 
	 * to perform point-wise subtraction between left and right sides, assigning the result
	 * to the left hand side.
	 *
	 * \param rhs The VectorValue instance where data is subtracted point-wise from the vector
	 * of the VectorValue on the left of the operator.
	 */
	VectorValue<T, D>& operator-=(const VectorValue<T, D>& rhs)
	{
		for (int i = 0; i < D; ++i)
		{
			v[i] -= rhs.v[i];
		}
		return *this;
	}

	//! Subtraction operator.
	/*!
	 * Subtraction operator, which returns a new VectorValue instance with its data equal
	 * to the result of point-wise addition of the data in the given operands.
	 *
	 * \param lhs The `VectorValue` instance on the left hand side of the subtraction operator. 
	 * \param rhs The `VectorValue` instance on the right hand side of the subtraction operator.
	 */
	friend VectorValue<T, D> operator-(VectorValue<T, D> lhs, const VectorValue<T, D>& rhs)
	{
		lhs -= rhs;
		return lhs;
	}

	//! Unary negative operator.
	/*!
	 * Operator applied to the left side of an instance that returns a new VectorValue 
	 * instance where all components of the data are determined by applying the unary 
	 * negative operator to the components of the calling instance.
	 */
	VectorValue<T, D> operator-() const
	{
		VectorValue<T, D> neg = *this;
		for (int i = 0; i < D; ++i)
		{
			neg.v[i] = -neg.v[i];
		}
		return neg;
	}

	//! Scalar multiplication compound assignment operator.
	/*!
	 * Compound assignment operator for multiplication with a scalar quantity.
	 *
	 * Assignment operator to multiply all the components of the data from the calling 
	 * instance with the scalar quantity on the right hand side.
	 *
	 * \param rhs The VectorValue instance where data is multiplied point-wise by the scalar quantity
	 * on the left of the operator.

	 * \tparam V The template parameter for the scalar type. SFINAE strategy is used so that
	 * the operator will only be compiled for types `V` which satisfy the condition that the
	 * resultant type of multiplication between instances of types `T` and `V` is implicitly
	 * convertible to `T`.
	 */
	template<typename V>
	auto operator*=(V const& rhs) -> typename mult_enable_if<T, V, VectorValue<T, D>&>::type
	{
		for (int i = 0; i < D; ++i)
		{
			v[i] *= rhs;
		}
		return *this;
	}

	//! Scalar multiplication operator.
	/*!
	 * Operator for multiplication with a scalar quantity.
	 *
	 * Assignment operator to multiply all the components of the data from the calling
	 * instance with the scalar quantity on the right hand side.
	 *
	 * \param lhs The VectorValue instance from which data is multiplied point-wise by the scalar 
	 * quantity on the right of the operator.
	 * \param rhs The scalar value used which multiplies the data in the VectorValue instance.

	 * \tparam V The template parameter for the scalar type. SFINAE strategy is used so that
	 * the operator will only be compiled for types `V` which satisfy the condition that the
	 * resultant type of multiplication between instances of types `T` and `V` is implicitly
	 * convertible to `T`.
	 */
	template<typename V>
	friend auto operator*(VectorValue<T, D> lhs, V const& rhs) -> 
		typename mult_enable_if<T, V, VectorValue<T, D>>::type
	{
		lhs *= rhs;
		return lhs;
	}

	//! Scalar multiplication operator.
	/*!
	 * Operator for multiplication with a scalar quantity with scalar quantity on left hand side.
	 *
	 * Assignment operator to multiply all the components of the data from the calling
	 * instance with the scalar quantity on the left hand side.
	 *
	 * \param lhs The scalar value used which multiplies the data in the VectorValue instance.
	 * \param rhs The VectorValue instance from which data is multiplied point-wise by the scalar 
	 * quantity on the left of the operator.

	 * \tparam V The template parameter for the scalar type. SFINAE strategy is used so that
	 * the operator will only be compiled for types `V` which satisfy the condition that the
	 * resultant type of multiplication between instances of types `T` and `V` is implicitly
	 * convertible to `T`.
	 */
	template<typename V>
	friend auto operator*(V const& lhs, VectorValue<T, D> const& rhs) -> 
		typename mult_enable_if<T, V, VectorValue<T, D>>::type
	{
		return rhs * lhs;
	}

	//! Dot product operator.
	/*!
	 * Operator for computing the dot product of two VectorValue instances.
	 *
	 * Assignment operator to multiply all the components of the data from the calling
	 * instance with the scalar quantity on the left hand side.
	 *
	 * \param lhs First VectorValue instance used to compute the dot product.
	 * \param rhs Second VectorValue instance used to compute the dot product.
	 */
	friend T operator*(VectorValue<T, D> const& lhs, VectorValue<T, D> const& rhs)
	{
		T sum = 0;
		for (int i = 0; i < D; ++i)
		{
			sum += lhs.v[i] * rhs.v[i];
		}
		return sum;
	}

	//! Scalar division compound assignment operator.
	/*!
	 * Operator for dividing by a scalar quantity, which is simply multiplying with the inverse
	 * of the scalar quantity (multiplicative inverse).
	 *
	 * Assignment operator to divide all the components of the data from the calling
	 * instance with the scalar quantity on the right hand side.
	 *
	 * \param rhs The scalar quantity of which the inverse will be naively determined by dividing
	 * the real valued number 1.0 by the scalar quantity.
	 * 
	 * \tparam V The template parameter for the scalar type. SFINAE strategy is used so that 
	 * the operator will only be compiled for types `V` which satisfy the condition that the 
	 * resultant type of multiplication between instances of types `T` and `V` is implicitly 
	 * convertible to `T`.
	 */
	template<typename V>
	auto operator/=(V const& rhs) -> 
		typename mult_enable_if<T, V, VectorValue<T, D>&>::type
	{
		*this *= (1.0 / rhs);
		return *this;
	}

	//! Scalar division operator.
	/*!
	 * Operator for dividing by a scalar quantity, which is equivalent multiplying with the
	 * inverse of the scalar quantity (multiplicative inverse). In particular, the inverse 
	 * is naively determined by dividing the real valued number 1.0 by the scalar quantity.
	 *
	 * This operator is only applicable for scalar quantities which are on the right hand side,
	 * because division is not commutative.
	 *
	 * \param lhs The VectorValue instance from where the data will be divided by the inverse
	 * of the scalar quantity.
	 * \param rhs The scalar quantity that divides the data of the VectorValue instance.
	 *
	 * \tparam V The template parameter for the scalar type. SFINAE strategy is used so that
	 * the operator will only be compiled for types `V` which satisfy the condition that the
	 * resultant type of multiplication between instances of types `T` and `V` is implicitly
	 * convertible to `T`.
	 */
	template<typename V>
	friend auto operator/(VectorValue<T, D> lhs, V const& rhs) -> 
		typename mult_enable_if<T, V, VectorValue<T, D>>::type
	{
		lhs /= rhs;
		return lhs;
	}

	//! Cross product compound assignment operator.
	/*!
	 * Compound assignment operator to calculate the conventional cross product and assign it
	 * to the calling instance.
	 * 
	 * Assignment operator to compute the cross product, which is only applicable to the
	 * 3-dimensional VectorValue template specialization. Operator is called by the left hand 
	 * instance, taking data from the right hand VectorValue to compute the cross product and
	 * assigning the result to the left hand side.
	 *
	 * \param rhs The VectorValue instance where data is taken from to compute the cross product
	 * the data from the left hand instance. 
	 */
	VectorValue<T, D>& operator%=(const VectorValue<T, D>& rhs)
	{
		if constexpr (D == 3)
		{
			T v0 = v[0], v1 = v[1], v2 = v[2];
			v[0] = (v1 * rhs.v[2] - v2 * rhs.v[1]);
			v[1] = (v2 * rhs.v[0] - v0 * rhs.v[2]);
			v[2] = (v0 * rhs.v[1] - v1 * rhs.v[0]);
			
			return *this;
		}
		else
		{
			// send an error message to indicate the dot product is being returned
			throw "attempting to use the cross product in a non 3-dimensional model";
		}
	}

	//! Cross product operator.
	/*!
	 * Assignment operator to calculate the conventional cross product.
	 *
	 * Assignment operator to compute the cross product, which is only applicable to the
	 * 3-dimensional VectorValue template specialization.
	 *
	 * \param lhs The VectorValue instance on the left hand side of the cross product
	 * operator.
	 * \param rhs The VectorValue instance on the right hand side of the cross product
	 * operator.
	 */
	friend VectorValue<T, D> operator%(VectorValue<T, D> lhs, const VectorValue<T, D>& rhs)
	{
		lhs %= rhs;
		return lhs;
	}

	//! Data access operator.
	/*!
	 * Access operator for the data memeber `v`. The position in the array is passed.
	 * The member access only returns a modifiable reference to the member data array 
	 * element. In the context of position information on a grid, the first index is 
	 * typically used for the horizontal component, \f$x\f$.
	 * 
	 * \param i The position in the member data array to access.
	 */
	T& operator[](int i)
	{
		return v[i];
	}

	//! Data access operator.
	/*!
	 * Const access operator for the data member `v`. See VectorValue::operator[](int). 
	 * The access is an unmodifiable reference to the member data array element.
	 *
	 * \param i The position in the data array to access.
	 */
	const T& operator[](int i) const
	{
		return v[i];
	}


};
template<typename T>
VectorValue(T)->VectorValue<T, 1>;
template<typename T>
VectorValue(T, T)->VectorValue<T, 2>;
template<typename T>
VectorValue(T, T, T)->VectorValue<T, 3>;


extern template struct VectorValue<double, 1>;
extern template struct VectorValue<double, 2>;
extern template struct VectorValue<double, 3>;
extern template struct VectorValue<std::complex<double>, 1>;
extern template struct VectorValue<std::complex<double>, 2>;
extern template struct VectorValue<std::complex<double>, 3>;
extern template struct VectorValue<float, 1>;
extern template struct VectorValue<float, 2>;
extern template struct VectorValue<float, 3>;
extern template struct VectorValue<std::complex<float>, 1>;
extern template struct VectorValue<std::complex<float>, 2>;
extern template struct VectorValue<std::complex<float>, 3>;



/* default implementations of distance and abs for VectorValue
 */

template<typename T, size_t D>
T abs(VectorValue<T, D> const& c) 
{ 
	double sum = 0;
	for (int i = 0; i < D; ++i)
	{
		sum += c.v[i];
	}
	using std::sqrt;
	return sqrt(sum);
}

template<typename T, size_t D>
T distance(VectorValue<T, D> const& c1, VectorValue<T, D> const& c2)
{
	using std::abs;
	return abs(c2 - c1);
}



template<>
inline double distance(VectorValue<double, 1> const& c1, VectorValue<double, 1> const& c2);
template<>
inline double distance(VectorValue<double, 2> const& c1, VectorValue<double, 2> const& c2);
template<>
inline double distance(VectorValue<double, 3> const& c1, VectorValue<double, 3> const& c2);
template<>
inline double abs(VectorValue<double, 1> const& c);
template<>
inline double abs(VectorValue<double, 2> const& c);
template<>
inline double abs(VectorValue<double, 3> const& c);


template<>
inline float distance(VectorValue<float, 1> const& c1, VectorValue<float, 1> const& c2);
template<>
inline float distance(VectorValue<float, 2> const& c1, VectorValue<float, 2> const& c2);
template<>
inline float distance(VectorValue<float, 3> const& c1, VectorValue<float, 3> const& c2);
template<>
inline float abs(VectorValue<float, 1> const& c);
template<>
inline float abs(VectorValue<float, 2> const& c);
template<>
inline float abs(VectorValue<float, 3> const& c);














