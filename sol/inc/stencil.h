
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 * 
 * Copyright (c) 2018-2020 by Steven A. Silber and Mikko Karttunen 
 * 
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published 
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details
 *
 * ***************************************************************************
 *
 * MODULE:  sol
 * PURPOSE:
 *
 * The stencil object defines the three finite difference functions, and
 * uses templates to specialize each one by dimension and order, for the
 * gradient, laplacian, bilaplacian and gradlaplacian.
 *
 * A specialized template for the stencil object in turn specializes each 
 * function to take advantage of generic programming.
 * 
 * A generic stencil can also be specialized for any stencils which are not
 * included in the specialization of first to fourth order derivatives, and
 * allows higher order derivatives to be programmatically or manually computed
 * by specializing the dimension, accuracy and derivative order.
 * 
 * ***************************************************************************
 */

#pragma once

#include "definitions.h"

#define DEFAULT_STENCIL_ACCURACY 2



/*!
 * \defgroup stencil Finite Difference Stencils
 * @{
 */


// ******************************************************************************************

//!
/*!
 * Basic stencil object on which concrete stencils are inherited. This stencil
 * is used for uniform grids where the spatial discretization is identical across
 * dimensions.
 */
template<typename Sp>
struct Stencil
{
	template <typename T>
	auto laplacian(T * const v) const
	{
		return (*static_cast<const Sp*>(this)).laplacian(v);
	}
	template <typename T>
	auto bilaplacian(T * const v) const
	{
		return (*static_cast<const Sp*>(this)).bilaplacian(v);
	}
	template <typename T>
	auto gradlaplacian(T * const v) const
	{
		return (*static_cast<const Sp*>(this)).gradlaplacian(v);
	}

	template<typename T>
	auto gradient(T * const v) const
	{
		return (*static_cast<const Sp*>(this)).gradient(v);
	}
};

//! A generalized derivative stencil.
/*!
 * A generalized derivative object that will produce stencils for any order of
 * derivative with a maximum boundary expansion of 3.
 * 
 * \tparam DD Dimension the stencil is applied in. 
 * \tparam OA The order of accuracy for the stencil.
 */
template<size_t DD, size_t OA = DEFAULT_STENCIL_ACCURACY>
struct GeneralizedStencil
{
	len_type dims[DD];
	double divh;

	GeneralizedStencil(const len_type* dims, double h) : divh{ 1. / h }
	{
		std::copy(dims, dims + DD, this->dims);
	}

	//! Determine the finite difference approximation to the derivative.
	/*!
	 * Determine the finite difference approximation to the derivative of
	 * the given order.
	 * 
	 * \tparam OD Order of the derivative.
	 */
	template<size_t OD, typename T>
	auto apply(T* const v) const
	{
		/* future implementation of a generalized stencil
		 * this function would use precomputed constants to approximate the derivative, given
		 * the dimension of the system and order of accuracy which the containing object is
		 * parameterized on
		 *
		 * currently this function returns the same value which is passed to it
		 */

		return *v;
	}

};

//! A specialization of the generalized derivative which returns the original value.
/*!
 * This specialization is meant to match the return value type for the 
 * derivative of the specified order. It does not actually compute the derivative
 * and is primarily used for testing by substituting it where a functional
 * stencil is intended. 
 */
template<>
struct GeneralizedStencil<0, DEFAULT_STENCIL_ACCURACY>
{
	GeneralizedStencil(...) {}

	template<size_t OD, typename T>
	auto apply(T* const v) const
	{
		using namespace std;
		using namespace symphas::math;
		if constexpr (OD % 2 == 1)
		{
			return *v;
		}
		else
		{
			return abs(*v);
		}
	}

};



// ******************************************************************************************


//! @}

// \cond

/* translations
 */

#define __DZ lenX * lenY
#define __DY lenX
#define __DX 1

#define __2DZ (2 * lenX * lenY)
#define __2DY (lenX + lenX)
#define __2DX 2

#define __3DZ (3 * lenX * lenY)
#define __3DY (lenX + lenX + lenX)
#define __3DX 3


// primary index
#define v0 v[0]


/* single combinations
 */

#define vx v[__DX]
#define vx_ v[-__DX]
#define vy v[__DY]
#define vy_ v[-__DY]
#define vz v[__DZ]
#define vz_ v[-__DZ]

#define vx2 v[__2DX]
#define vx2_ v[-__2DX]
#define vy2 v[__2DY]
#define vy2_ v[-__2DY]
#define vz2 v[__2DZ]
#define vz2_ v[-__2DZ]

#define vx3 v[__3DX]
#define vx3_ v[-__3DX]
#define vy3 v[__3DY]
#define vy3_ v[-__3DY]
#define vz3 v[__3DZ]
#define vz3_ v[-__3DZ]


/*
 * all xy combinations
 */

#define vxy v[__DX + __DY]
#define vxy_ v[__DX - __DY]
#define vx_y v[-__DX + __DY]
#define vx_y_ v[-__DX - __DY]

#define vx2y v[__2DX + __DY]
#define vx2y_ v[__2DX - __DY]
#define vx2_y v[-__2DX + __DY]
#define vx2_y_ v[-__2DX - __DY]

#define vxy2 v[__DX + __2DY]
#define vxy2_ v[__DX - __2DY]
#define vx_y2 v[-__DX + __2DY]
#define vx_y2_ v[-__DX - __2DY]

#define vx3y v[__3DX + __DY]
#define vx3y_ v[__3DX - __DY]
#define vx3_y v[-__3DX + __DY]
#define vx3_y_ v[-__3DX - __DY]

#define vxy3 v[__DX + __3DY]
#define vxy3_ v[__DX - __3DY]
#define vx_y3 v[-__DX + __3DY]
#define vx_y3_ v[-__DX - __3DY]

#define vx2y2 v[__2DX + __2DY]
#define vx2y2_ v[__2DX - __2DY]
#define vx2_y2 v[-__2DX + __2DY]
#define vx2_y2_ v[-__2DX - __2DY]

#define vx3y3 v[__3DX + __3DY]
#define vx3y3_ v[__3DX - __3DY]
#define vx3_y3 v[-__3DX + __3DY]
#define vx3_y3_ v[-__3DX - __3DY]

/*
* all xz combinations
*/

#define vxz v[__DX + __DZ]
#define vxz_ v[__DX - __DZ]
#define vx_z v[-__DX + __DZ]
#define vx_z_ v[-__DX - __DZ]

#define vx2z v[__2DX + __DZ]
#define vx2z_ v[__2DX - __DZ]
#define vx2_z v[-__2DX + __DZ]
#define vx2_z_ v[-__2DX - __DZ]

#define vxz2 v[__DX + __2DZ]
#define vxz2_ v[__DX - __2DZ]
#define vx_z2 v[-__DX + __2DZ]
#define vx_z2_ v[-__DX - __2DZ]

#define vx2z2 v[__2DX + __2DZ]
#define vx2z2_ v[__2DX - __2DZ]
#define vx2_z2 v[-__2DX + __2DZ]
#define vx2_z2_ v[-__2DX - __2DZ]


/*
* all yz combinations
*/

#define vyz v[__DY + __DZ]
#define vyz_ v[__DY - __DZ]
#define vy_z v[-__DY + __DZ]
#define vy_z_ v[-__DY - __DZ]

#define vyz2 v[__DY + __2DZ]
#define vyz2_ v[__DY - __2DZ]
#define vy_z2 v[-__DY + __2DZ]
#define vy_z2_ v[-__DY - __2DZ]

#define vy2z v[__2DY + __DZ]
#define vy2z_ v[__2DY - __DZ]
#define vy2_z v[-__2DY + __DZ]
#define vy2_z_ v[-__2DY - __DZ]

#define vy2z2 v[__2DY + __2DZ]
#define vy2z2_ v[__2DY - __2DZ]
#define vy2_z2 v[-__2DY + __2DZ]
#define vy2_z2_ v[-__2DY - __2DZ]


/*
* xyz combinations
*/

#define vxyz v[__DX + __DY + __DZ]
#define vxyz_ v[__DX + __DY - __DZ]
#define vxy_z v[__DX - __DY + __DZ]
#define vxy_z_ v[__DX - __DY - __DZ]
#define vx_yz v[-__DX + __DY + __DZ]
#define vx_yz_ v[-__DX + __DY - __DZ]
#define vx_y_z v[-__DX - __DY + __DZ]
#define vx_y_z_ v[-__DX - __DY - __DZ]

#define vx2yz v[__2DX + __DY + __DZ]
#define vx2yz_ v[__2DX + __DY - __DZ]
#define vx2y_z v[__2DX - __DY + __DZ]
#define vx2y_z_ v[__2DX - __DY - __DZ]
#define vx2_yz v[-__2DX + __DY + __DZ]
#define vx2_yz_ v[-__2DX + __DY - __DZ]
#define vx2_y_z v[-__2DX - __DY + __DZ]
#define vx2_y_z_ v[-__2DX - __DY - __DZ]

#define vxy2z v[__DX + __2DY + __DZ]
#define vxy2z_ v[__DX + __2DY - __DZ]
#define vxy2_z v[__DX - __2DY + __DZ]
#define vxy2_z_ v[__DX - __2DY - __DZ]
#define vx_y2z v[-__DX + __2DY + __DZ]
#define vx_y2z_ v[-__DX + __2DY - __DZ]
#define vx_y2_z v[-__DX - __2DY + __DZ]
#define vx_y2_z_ v[-__DX - __2DY - __DZ]

#define vxyz2 v[__DX + __DY + __2DZ]
#define vxyz2_ v[__DX + __DY - __2DZ]
#define vxy_z2 v[__DX - __DY + __2DZ]
#define vxy_z2_ v[__DX - __DY - __2DZ]
#define vx_yz2 v[-__DX + __DY + __2DZ]
#define vx_yz2_ v[-__DX + __DY - __2DZ]
#define vx_y_z2 v[-__DX - __DY + __2DZ]
#define vx_y_z2_ v[-__DX - __DY - __2DZ]

#define vx2y2z2 v[__2DX + __2DY + __2DZ]
#define vx2y2z2_ v[__2DX + __2DY - __2DZ]
#define vx2y2_z2 v[__2DX - __2DY + __2DZ]
#define vx2y2_z2_ v[__2DX - __2DY - __2DZ]
#define vx2_y2z2 v[-__2DX + __2DY + __2DZ]
#define vx2_y2z2_ v[-__2DX + __2DY - __2DZ]
#define vx2_y2_z2 v[-__2DX - __2DY + __2DZ]
#define vx2_y2_z2_ v[-__2DX - __2DY - __2DZ]


// \endcond

