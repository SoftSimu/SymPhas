
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
 * MODULE:  sol
 * PURPOSE: Stencil of 6th order accuracy, 2 dimensions.
 * Isotropic stencils for the Laplacian, bilaplacian, and
 * gradlaplacian on a square lattice.
 *
 * ***************************************************************************
 */

#pragma once

#include "stencil.h"

template<size_t L>
struct apply_laplacian_2d6h;
template<size_t B>
struct apply_bilaplacian_2d6h;
template<size_t G>
struct apply_gradlaplacian_2d6h;

// =====================================================================
// 6th-order isotropic 2D Laplacian
// =====================================================================

//! 33-point 6th-order isotropic Laplacian.
/*!
 * Weights (multiply by 1/(2520 h^2)):
 *   center:     -11308
 *   (+-1,0):      2232   (x4)
 *   (+-1,+-1):     972   (x4)
 *   (+-2,+-1):    -216   (x8)
 *   (+-2,+-2):      27   (x4)
 *   (+-3,0):        -8   (x4)
 *   (+-3,+-1):      18   (x8)
 *
 * Note: (+-2,0) has weight 0 and is excluded.
 */
template<>
struct apply_laplacian_2d6h<33>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh2, const len_type(&stride)[2])
	{
		return (
			-8. * (vx3 + vx3_ + vy3 + vy3_)
			+ 18. * (vx3y + vx3_y + vx3y_ + vx3_y_ + vxy3 + vx_y3 + vxy3_ + vx_y3_)
			+ 27. * (vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			- 216. * (vx2y + vx2_y + vx2y_ + vx2_y_ + vxy2 + vx_y2 + vxy2_ + vx_y2_)
			+ 972. * (vxy + vx_y + vxy_ + vx_y_)
			+ 2232. * (vx + vx_ + vy + vy_)
			- 11308. * v0
			) * divh2 * (1.0 / 2520);
	}
};

// =====================================================================
// 6th-order isotropic 2D Bilaplacian
// =====================================================================

//! 37-point 6th-order isotropic bilaplacian.
/*!
 * Weights (multiply by 1/(180 h^4)):
 *   center:       3116
 *   (+-1,0):      -764   (x4)
 *   (+-1,+-1):    -374   (x4)
 *   (+-2,0):        42   (x4)
 *   (+-2,+-1):     188   (x8)
 *   (+-2,+-2):     -29   (x4)
 *   (+-3,0):         4   (x4)
 *   (+-3,+-1):     -17   (x8)
 */
template<>
struct apply_bilaplacian_2d6h<37>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh4, const len_type(&stride)[2])
	{
		return (
			-17. * (vx3y + vx3_y + vx3y_ + vx3_y_ + vxy3 + vx_y3 + vxy3_ + vx_y3_)
			+ 4. * (vx3 + vx3_ + vy3 + vy3_)
			- 29. * (vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			+ 188. * (vx2y + vx2_y + vx2y_ + vx2_y_ + vxy2 + vx_y2 + vxy2_ + vx_y2_)
			+ 42. * (vx2 + vx2_ + vy2 + vy2_)
			- 374. * (vxy + vx_y + vxy_ + vx_y_)
			- 764. * (vx + vx_ + vy + vy_)
			+ 3116. * v0
			) * divh4 * (1.0 / 180);
	}
};

// =====================================================================
// 6th-order isotropic 2D Gradlaplacian (x-component)
// =====================================================================

//! 30-point 6th-order isotropic gradlaplacian.
/*!
 * Weights for x-component (multiply by 1/(720 h^3)):
 *   (+-1,0):      -382   (x2 antisym)
 *   (+-1,+-1):    -374   (x4 antisym)
 *   (+-2,0):        42   (x2 antisym)
 *   (+-2,+-1):     376   (x4 antisym)
 *   (+-1,+-2):     188   (x4 antisym)
 *   (+-2,+-2):     -58   (x4 antisym)
 *   (+-3,0):         6   (x2 antisym)
 *   (+-3,+-1):     -51   (x4 antisym)
 *   (+-1,+-3):     -17   (x4 antisym)
 *
 * y-component obtained by 90-degree rotation.
 */
template<>
struct apply_gradlaplacian_2d6h<30>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh3, const len_type(&stride)[2])
	{
		auto x = (
			6. * (vx3 - vx3_)
			- 51. * (vx3y - vx3_y + vx3y_ - vx3_y_)
			- 17. * (vxy3 - vx_y3 + vxy3_ - vx_y3_)
			- 58. * (vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_)
			+ 376. * (vx2y - vx2_y + vx2y_ - vx2_y_)
			+ 188. * (vxy2 - vx_y2 + vxy2_ - vx_y2_)
			+ 42. * (vx2 - vx2_)
			- 374. * (vxy - vx_y + vxy_ - vx_y_)
			- 382. * (vx - vx_)
			) * divh3 * (1.0 / 720);

		auto y = (
			6. * (vy3 - vy3_)
			- 51. * (vxy3 - vxy3_ + vx_y3 - vx_y3_)
			- 17. * (vx3y - vx3y_ + vx3_y - vx3_y_)
			- 58. * (vx2y2 - vx2y2_ + vx2_y2 - vx2_y2_)
			+ 376. * (vxy2 - vxy2_ + vx_y2 - vx_y2_)
			+ 188. * (vx2y - vx2y_ + vx2_y - vx2_y_)
			+ 42. * (vy2 - vy2_)
			- 374. * (vxy - vxy_ + vx_y - vx_y_)
			- 382. * (vy - vy_)
			) * divh3 * (1.0 / 720);

		return VectorValue<T, 2>{ x, y };
	}
};
