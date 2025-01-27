
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
 * SymPhas is distributed with the faith that it will be divhelpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 *
 * MODULE:  sol
 * PURPOSE: Stencil of 2nd order accuracy, both 2 and 3 dimensions.
 * The stencil is defined by returning values corresponding
 * to the number of points used in the finite difference approximation.
 *
 * ***************************************************************************
 */

#pragma once


#include "stencil.h"



template<size_t L>
struct apply_laplacian_2d4h;
template<size_t B>
struct apply_bilaplacian_2d4h;
template<size_t G>
struct apply_gradlaplacian_2d4h;


template<>
struct apply_laplacian_2d4h<9>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh2, const len_type(&stride)[2])
	{
		return (-vx2_ - vx2 - vy2 - vy2_
			+ 16. * (vx_ + vx + vy + vy_) - 60. * v0
			) * divh2 * (1.0 / 12);
	}
};

template<>
struct apply_laplacian_2d4h<17>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh2, const len_type(&stride)[2])
	{
		return (-(vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			- 8. * (vx2 + vx2_ + vy2 + vy2_) +
			16. * (vxy + vx_y + vxy_ + vx_y_) +
			128. * (vx + vx_ + vy + vy_)
			- 540. * v0
			) * divh2 * (1.0 / 120);
	}
};

template<>
struct apply_laplacian_2d4h<21>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh2, const len_type(&stride)[2])
	{
		return (-2. * (vx2y + vx2_y + vx2y_ + vx2_y_ +
			vxy2 + vx_y2 + vxy2_ + vx_y2_)
			- (vx2 + vx2_ + vy2 + vy2_) +
			16. * (vxy + vx_y + vxy_ + vx_y_) +
			52. * (vx + vx_ + vy + vy_)
			- 252. * v0
			) * divh2 * (1.0 / 60);
	}
};


template<>
struct apply_bilaplacian_2d4h<21>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh4, const len_type(&stride)[2])
	{
		return (-4. * (vx3 + vx3_ + vy3 + vy3_)
			- (vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			+ 50. * (vx2 + vx2_ + vy2 + vy2_)
			+ 64. * (vxy + vx_y + vxy_ + vx_y_)
			- 284. * (vx + vx_ + vy + vy_)
			+ 700. * v0
			) * divh4 * (1.0 / 24);
	}
};

template<>
struct apply_bilaplacian_2d4h<25>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh4, const len_type(&stride)[2])
	{
		return (-(vx3 + vx3_ + vy3 + vy3_)
			- (vx2y + vx2_y + vx2y_ + vx2_y_ + vxy2 + vx_y2 + vxy2_ + vx_y2_) +
			14. * (vx2 + vx2_ + vy2 + vy2_) +
			20. * (vxy + vx_y + vxy_ + vx_y_)
			- 77. * (vx + vx_ + vy + vy_) +
			184. * v0
			) * divh4 * (1.0 / 6);
	}
};

template<>
struct apply_bilaplacian_2d4h<33>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh4, const len_type(&stride)[2])
	{
		return (-17. * (vx3y3 + vx3_y3 + vx3y3_ + vx3_y3_)
			- 236. * (vx3 + vx3_ + vy3 + vy3_)
			+ 351. * (vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			- 756. * (vx2y + vx2_y + vx2y_ + vx2_y_ + vxy2 + vx_y2 + vxy2_ + vx_y2_)
			+ 4050. * (vx2 + vx2_ + vy2 + vy2_)
			+ 5049. * (vxy + vx_y + vxy_ + vx_y_)
			- 19116. * (vx + vx_ + vy + vy_)
			+ 45724. * v0
			) * divh4 * (1.0 / 1620);
	}
};

template<>
struct apply_bilaplacian_2d4h<37>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh4, const len_type(&stride)[2])
	{
		return (-17. * (vx3y + vx3_y + vx3y_ + vx3_y_ + vxy3 + vx_y3 + vxy3_ + vx_y3_)
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

template<>
struct apply_gradlaplacian_2d4h<14>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh3, const len_type(&stride)[2])
	{
		auto x = (-6. * (vx3 - vx3_)
			- (vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_)
			+ 50. * (vx2 - vx2_)
			+ 32. * (vxy - vx_y + vxy_ - vx_y_)
			- 142. * (vx - vx_)
			) * divh3 * (1.0 / 48);
		return x;
	}
};

template<>
struct apply_gradlaplacian_2d4h<18>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh3, const len_type(&stride)[2])
	{
		auto x = (-3. * (vx3)
			-2. * (vx2y - vx2_y + vx2y_ - vx2_y_)
			+ 28. * (vx2)
			-(vxy2 - vx_y2 + vxy2_ - vx_y2_)
			+ 20. * (vxy - vx_y + vxy_ - vx_y_)
			- 77. * (vx)
			) * divh3 * (1.0 / 24);
		return x;
	}
};

template<>
struct apply_gradlaplacian_2d4h<26>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh3, const len_type(&stride)[2])
	{
		auto x = (-17. * (vx3y3 - vx3_y3 + vx3y3_ - vx3_y3_)
			- 236. * (vx3)
			+234. * (vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_)
			- 504. * (vx2y - vx2_y + vx2y_ - vx2_y_)
			+ 2700. * (vx2)
			-252. * (vxy2 - vx_y2 + vxy2_ - vx_y2_)
			+ 1683. * (vxy - vx_y + vxy_ - vx_y_)
			- 6372. * (vx)
			) / (2160. * divh3);
		return x;
	}
};

template<>
struct apply_gradlaplacian_2d4h<30>
{
	template<typename T>
	__device__ __host__ auto operator()(T* const v, double divh3, const len_type(&stride)[2])
	{
		auto x = (-51. * (vx3y - vx3_y + vx3y_ - vx3_y_)
			+ 12. * (vx3)
			-58. * (vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_)
			+ 376. * (vx2y - vx2_y + vx2y_ - vx2_y_)
			+ 84. * (vx2)
			-17. * (vxy3 - vx_y3 + vxy3_ - vx_y3_)
			+ 188. * (vxy2 - vx_y2 + vxy2_ - vx_y2_)
			- 374. * (vxy - vx_y + vxy_ - vx_y_)
			- 764. * (vx)
			) * divh3 * (1.0 / 720);
		return x;
	}
};


