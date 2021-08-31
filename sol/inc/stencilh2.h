
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
 * PURPOSE: Stencil of 2nd order accuracy, both 2 and 3 dimensions.
 * The stencil is defined by returning values corresponding
 * to the number of points used in the finite difference approximation.
 *
 * ***************************************************************************
 */

#pragma once

/*
 * author:	Steven Silber
 * date:	June, 2018
 *
 * 
 */

#include "stencil.h"


template<size_t L>
struct apply_laplacian_1d2h;
template<size_t B>
struct apply_bilaplacian_1d2h;
template<size_t G>
struct apply_gradlaplacian_1d2h;


template<>
struct apply_laplacian_1d2h<3>
{
	template<typename T>
	auto operator()(T* const v, double h2)
	{
		return (vx_ + vx - 2. * v0) * h2;
	}
};

template<>
struct apply_bilaplacian_1d2h<5>
{
	template<typename T>
	auto operator()(T* const v, double h4)
	{
		return (vx2 - 4 * vx + 6 * v0 - 4 * vx_ + vx2_) * h4;
	}
};

template<>
struct apply_gradlaplacian_1d2h<4>
{
	template<typename T>
	auto operator()(T* const v, double h3)
	{
		return VectorValue<T, 1>{ (vx2 - 2 * vx + 2 * vx_ - vx2_) * h3 };
	}
};



template<size_t L>
struct apply_laplacian_2d2h;
template<size_t B>
struct apply_bilaplacian_2d2h;
template<size_t G>
struct apply_gradlaplacian_2d2h;

template<>
struct apply_laplacian_2d2h<5>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX)
	{
		return (vx_ + vx + vy + vy_ - 4. * v0) * h2;
	}
};

template<>
struct apply_laplacian_2d2h<9>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX)
	{
		return (vx_y + vxy_ + vxy + vx_y_
			+ 4. * vx + 4. * vy + 4. * vy_ + 4. * vx_ - 20. * v0
			) * h2 * (1.0 / 6);
	}
};




template<>
struct apply_bilaplacian_2d2h<13>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX)
	{
		return (
			vy2_ + vy2 + vx2_ + vx2 +
			2. * (vx_y_ + vxy_ + vx_y + vxy)
			- 8. * (vx_ + vx + vy_ + vy)
			+ 20. * v0
			) * h4;
	}
};

template<>
struct apply_bilaplacian_2d2h<17>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX)
	{
		return (
			(vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_)
			+ 10. * (vx2 + vx2_ + vy2 + vy2_)
			+ 8. * (vxy + vx_y + vxy_ + vx_y_)
			- 64. * (vx + vx_ + vy + vy_)
			+ 180. * v0
			) * h4 * (1.0 / 12);
	}
};

template<>
struct apply_bilaplacian_2d2h<21>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX)
	{
		return (
			(vx2y + vx2_y + vx2y_ + vx2_y_ + vxy2 + vx_y2 + vxy2_ + vx_y2_)
			+ (vx2 + vx2_ + vy2 + vy2_)
			- 2. * (vxy + vx_y + vxy_ + vx_y_)
			- 10. * (vx + vx_ + vy + vy_)
			+ 36. * v0
			) * h4 * (1.0 / 12);
	}
};



template<>
struct apply_gradlaplacian_2d2h<6>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX)
	{
		auto x = ((vx2y - vx2_y + vx2y_ - vx2_y_)
			- 4. * (vx - vx_)
			) * h3 * (1.0 / 4);
		auto y = ((vxy2 - vxy2_ + vx_y2 - vx_y2_)
			- 4. * (vy - vy_)
			) * h3 * (1.0 / 4);
		return VectorValue<T, 2>{ x, y };
	}
};

template<>
struct apply_gradlaplacian_2d2h<8>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX)
	{
		auto x = ((vx2 - vx2_)
			+ (vxy - vx_y + vxy_ - vx_y_)
			- 4. * (vx - vx_)
			) * h3 * (1.0 / 2);
		auto y = ((vy2 - vy2_)
			+ (vxy - vxy_ + vx_y - vx_y_)
			- 4. * (vy - vy_)
			) * h3 * (1.0 / 2);
		return VectorValue<T, 2>{ x, y };
	}
};

template<>
struct apply_gradlaplacian_2d2h<12>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX)
	{
		auto x = ((vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_)
			+ 10. * (vx2 - vx2_)
			+ 4. * (vxy - vx_y + vxy_ - vx_y_)
			- 32. * (vx - vx_)
			) * h3 * (1.0 / 24);
		auto y = ((vx2y2 - vx2y2_ + vx2_y2 - vx2_y2_)
			+ 10. * (vy2 - vy2_)
			+ 4. * (vxy - vxy_ + vx_y - vx_y_)
			- 32. * (vy - vy_)
			) * h3 * (1.0 / 24);
		return VectorValue<T, 2>{ x, y };
	}

};

template<>
struct apply_gradlaplacian_2d2h<16>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX)
	{
		auto x = (2. * (vx2y - vx2_y + vx2y_ - vx2_y_)
			+ 2. * (vx2 - vx2_)
			+ (vxy2 - vx_y2 + vxy2_ - vx_y2_)
			- 2. * (vxy - vx_y + vxy_ - vx_y_)
			- 10. * (vx - vx_)
			) * h3 * (1.0 / 12);
		auto y = (2. * (vxy2 - vxy2_ + vx_y2 - vx_y2_)
			+ 2. * (vy2 - vy2_)
			+ (vx2y - vx2y_ + vx2_y - vx2_y_)
			- 2. * (vxy - vxy_ + vx_y - vx_y_)
			- 10. * (vy - vy_)
			) * h3 * (1.0 / 12);
		return VectorValue<T, 2>{ x, y };
	}
};






template<size_t L>
struct apply_laplacian_3d2h;
template<size_t B>
struct apply_bilaplacian_3d2h;
template<size_t G>
struct apply_gradlaplacian_3d2h;

template<>
struct apply_laplacian_3d2h<7>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX, len_type lenY)
	{
		return (
			vx_ + vx + vy + vy_ +
			vz + vz_ - 6. * v0
			) * h2;
	}
};

template<>
struct apply_laplacian_3d2h<15>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX, len_type lenY)
	{
		return (
			vx_y_z_ + vxy_z_ + vx_yz_ +
			vxyz_ + vx_y_z + vxy_z + vx_yz + vxyz +
			8. * (vx_ + vx + vy + vy_ + vz + vz_ - 7. * v0)
			) * h2 * (1.0 / 12);
	}
};

template<>
struct apply_laplacian_3d2h<19>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX, len_type lenY)
	{
		return (
			vx_y_ + vxy_ + vx_y + vxy +
			vx_z_ + vxz_ + vx_z + vxz +
			vy_z_ + vyz_ + vy_z + vyz +
			2. * (vx_ + vx + vy + vy_ +
				vz + vz) - 24. * v0
			) * h2 * (1.0 / 6);
	}
};

template<>
struct apply_laplacian_3d2h<21>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX, len_type lenY)
	{
		return (
			-(vx_y_z_ + vxy_z_ + vx_yz_ +
				vxyz_ + vx_y_z + vxy_z + vx_yz + vxyz) +
			4. * (vx_y_ + vxy_ + vx_y + vxy +
				vx_z_ + vxz_ + vx_z + vxz +
				vy_z_ + vyz_ + vy_z + vyz) +
			-40. * v0
			) * h2 * (1.0 / 12);
	}
};

template<>
struct apply_laplacian_3d2h<27>
{
	template<typename T>
	auto operator()(T* const v, double h2, len_type lenX, len_type lenY)
	{
		return (
			vx_y_z_ + vxy_z_ + vx_yz_ +
			vxyz_ + vx_y_z + vxy_z + vx_yz + vxyz +
			3. * (vx_y_ + vxy_ + vx_y + vxy +
				vx_z_ + vxz_ + vx_z + vxz +
				vy_z_ + vyz_ + vy_z + vyz) +
			14. * (vx_ + vx + vy + vy_ + vz + vz_) - 128. * v0
			) * h2 * (1.0 / 30);
	}
};

template<>
struct apply_bilaplacian_3d2h<21>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX, len_type lenY)
	{
		return (
			(vx2y2 + vx2_y2 + vx2y2_ + vx2_y2_ +
				vx2z2_ + vx2_z2_ + vx2z2 + vx2_z2 +
				vy2z2_ + vy2z2 + vy2_z2_ + vy2_z2)
			- 4. * (vxyz + vx_yz + vxy_z + vx_y_z +
				vxyz_ + vx_yz_ + vxy_z_ + vx_y_z_)
			+ 20. * v0
			) * h4 * (1.0 / 4);
	}
};

template<>
struct apply_bilaplacian_3d2h<25>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX, len_type lenY)
	{
		return (
			2. * (vxy + vx_y + vxy_ + vx_y_ +
				vxz_ + vx_z_ + vxz + vx_z +
				vyz_ + vyz + vy_z_ + vy_z)
			+ (vx2 + vx2_ + vy2 + vy2_ + vz2 + vz2_)
			- 12. * (vx + vx_ + vy + vy_ + vz + vz_)
			+ 42. * v0
			) * h4;
	}
};

template<>
struct apply_bilaplacian_3d2h<41>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX, len_type lenY)
	{
		return (
			(vx2y2z2 + vx2_y2z2 + vx2y2_z2 + vx2_y2_z2 +
				vx2y2z2_ + vx2_y2z2_ + vx2y2_z2_ + vx2_y2_z2_)
			- 40. * (vxyz + vx_yz + vxy_z + vx_y_z +
				vxyz_ + vx_yz_ + vxy_z_ + vx_y_z_)
			+ 96. * (vxy + vx_y + vxy_ + vx_y_ +
				vxz_ + vx_z_ + vxz + vx_z +
				vyz_ + vyz + vy_z_ + vy_z)
			+ 20. * (vx2 + vx2_ + vy2 + vy2_ + vz2 + vz2_)
			- 320. * (vx + vx_ + vy + vy_ + vz + vz_)
			+ 960. * v0
			) * h4 * (1.0 / 24);
	}
};

template<>
struct apply_bilaplacian_3d2h<52>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX, len_type lenY)
	{
		return (
			-(vx2y2z2 + vx2_y2z2 + vx2y2_z2 + vx2_y2_z2 +
				vx2y2z2_ + vx2_y2z2_ + vx2y2_z2_ + vx2_y2_z2_)
			+ 10. * (vxyz2 + vx_yz2 + vxy_z2 + vx_y_z2 +
				vxyz2_ + vx_yz2_ + vxy_z2_ + vx_y_z2_ +
				vxy2z + vx_y2z + vxy2_z + vx_y2_z +
				vxy2z_ + vx_y2z_ + vxy2_z_ + vx_y2_z_ +
				vx2yz + vx2_yz + vx2y_z + vx2_y_z +
				vx2yz_ + vx2_yz_ + vx2y_z_ + vx2_y_z_)
			- 20. * (vxyz + vx_yz + vxy_z + vx_y_z +
				vxyz_ + vx_yz_ + vxy_z_ + vx_y_z_)
			- 36. * (vxy + vx_y + vxy_ + vx_y_ +
				vxz_ + vx_z_ + vxz + vx_z +
				vyz_ + vyz + vy_z_ + vy_z)
			+ 360. * v0
			) * h4 * (1.0 / 36);
	}
};

template<>
struct apply_bilaplacian_3d2h<57>
{
	template<typename T>
	auto operator()(T* const v, double h4, len_type lenX, len_type lenY)
	{
		return (
			(vyz2 + vy_z2 + vyz2_ + vy_z2_ +
				vy2z + vy2_z + vy2z_ + vy2_z_ +
				vx2z + vx2_z + vx2z_ + vx2_z_ +
				vxz2 + vx_z2 + vxz2_ + vx_z2_ +
				vxy2 + vx_y2 + vxy2_ + vx_y2_ +
				vx2y + vx2_y + vx2y_ + vx2_y_)
			+ 3. * (vxyz + vx_yz + vxy_z + vx_y_z +
				vxyz_ + vx_yz_ + vxy_z_ + vx_y_z_)
			- 8. * (vxy + vx_y + vxy_ + vx_y_ +
				vxz_ + vx_z_ + vxz + vx_z +
				vyz_ + vyz + vy_z_ + vy_z)
			- (vx2 + vx2_ + vy2 + vy2_ + vz2 + vz2_)
			+ 4. * (vx + vx_ + vy + vy_ + vz + vz_)
			+ 30. * v0
			) * h4 * (1.0 / 3);
	}
};

template<>
struct apply_gradlaplacian_3d2h<10>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX, len_type lenY)
	{
		auto x = ((vx2yz - vx2_yz + vx2y_z - vx2_y_z +
			vx2yz_ - vx2_yz_ + vx2y_z_ - vx2_y_z_)
			- 8. * (vx - vx_)
			) * h3 * (1.0 / 8);
		auto y = ((vxy2z - vxy2_z + vx_y2z - vx_y2_z +
			vxy2z_ - vxy2_z_ + vx_y2z_ - vx_y2_z_)
			- 8. * (vy - vy_)
			) * h3 * (1.0 / 8);
		auto z = ((vxyz2 - vxyz2_ + vxy_z2 - vxy_z2_ +
			vx_yz2 - vx_yz2_ + vx_y_z2 - vx_y_z2_)
			- 8. * (vz - vz_)
			) * h3 * (1.0 / 8);
		return VectorValue<T, 3>{ x, y, z };
	}
};

template<>
struct apply_gradlaplacian_3d2h<12>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX, len_type lenY)
	{
		auto x = ((vx2 - vx2_)
			+ (vxy - vx_y + vxy_ - vx_y_ +
				vxz_ - vx_z_ + vxz - vx_z)
			- 6. * (vx - vx_)
			) * h3 * (1.0 / 2);
		auto y = ((vy2 - vy2_)
			+ (vxy - vxy_ + vx_y - vx_y_ +
				vyz_ - vy_z_ + vyz - vy_z)
			- 6. * (vy - vy_)
			) * h3 * (1.0 / 2);
		auto z = ((vz2 - vz2_)
			+ (vyz - vyz_ + vy_z - vy_z_ +
				vx_z - vx_z_ + vxz - vxz_)
			- 6. * (vz - vz_)
			) * h3 * (1.0 / 2);
		return VectorValue<T, 3>{ x, y, z };
	}
};

template<>
struct apply_gradlaplacian_3d2h<28>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX, len_type lenY)
	{
		auto x = ((vx2y2 - vx2_y2 + vx2y2_ - vx2_y2_ +
			vx2z2_ - vx2_z2_ + vx2z2 - vx2_z2)
			+ 8. * (vx2 - vx2_)
			+ 6. * (vxyz - vx_yz + vxy_z - vx_y_z +
				vxyz_ - vx_yz_ + vxy_z_ - vx_y_z_)
			- 8. * (vxy - vx_y + vxy_ - vx_y_ +
				vxz_ - vx_z_ + vxz - vx_z)
			- 16. * (vx - vx_)
			) * h3 * (1.0 / 24);
		auto y = ((vx2y2 - vx2y2_ + vx2_y2 - vx2_y2_ +
			vy2z2_ - vy2_z2_ + vy2z2 - vy2_z2)
			+ 8. * (vy2 - vy2_)
			+ 6. * (vxyz - vxy_z + vx_yz - vx_y_z +
				vxyz_ - vxy_z_ + vx_yz_ - vx_y_z_)
			- 8. * (vxy - vxy_ + vx_y - vx_y_ +
				vyz_ - vy_z_ + vyz - vy_z)
			- 16. * (vy - vy_)
			) * h3 * (1.0 / 24);
		auto z = ((vx2y2 - vx2y2_ + vx2_y2 - vx2_y2_ +
			vx2_z2 - vx2_z2_ + vx2z2 - vx2z2_)
			+ 8. * (vz2 - vz2_)
			+ 6. * (vxyz - vxyz_ + vxy_z - vxy_z_ +
				vx_yz - vx_yz_ + vx_y_z - vx_y_z_)
			- 8. * (vxy - vxy_ + vx_y - vx_y_ +
				vx_z - vx_z_ + vxz - vxz_)
			- 16. * (vz - vz_)
			) * h3 * (1.0 / 24);
		return VectorValue<T, 3>{ x, y, z };
	}
};

template<>
struct apply_gradlaplacian_3d2h<36>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX, len_type lenY)
	{
		auto x = (2. * (vx2z - vx2_z + vx2z_ - vx2_z_ +
			vx2y - vx2_y + vx2y_ - vx2_y_)
			- 2. * (vx2 - vx2_)
			+ (vxz2 - vx_z2 + vxz2_ - vx_z2_ +
				vxy2 - vx_y2 + vxy2_ - vx_y2_)
			+ 3. * (vxyz - vx_yz + vxy_z - vx_y_z +
				vxyz_ - vx_yz_ + vxy_z_ - vx_y_z_)
			- 8. * (vxy - vx_y + vxy_ - vx_y_ +
				vxz_ - vx_z_ + vxz - vx_z)
			+ 4. * (vx - vx_)
			) * h3 * (1.0 / 12);
		auto y = (2. * (vy2z - vy2_z + vy2z_ - vy2_z_ +
			vxy2 - vxy2_ + vx_y2 - vx_y2_)
			- 2. * (vy2 - vy2_)
			+ (vyz2 - vy_z2 + vyz2_ - vy_z2_ +
				vx2y - vx2y_ + vx2_y - vx2_y_)
			+ 3. * (vxyz - vxy_z + vx_yz - vx_y_z +
				vxyz_ - vxy_z_ + vx_yz_ - vx_y_z_)
			- 8. * (vxy - vxy_ + vx_y - vx_y_ +
				vyz_ - vy_z_ + vyz - vy_z)
			+ 4. * (vy - vy_)
			) * h3 * (1.0 / 12);
		auto z = (2. * (vxz2 - vxz2_ + vx_z2 - vx_z2_ +
			vxy2 - vxy2_ + vx_y2 - vx_y2_)
			- 2. * (vz2 - vz2_)
			+ (vx2z - vx2z_ + vx2_z - vx2_z_ +
				vx2y - vx2y_ + vx2_y - vx2_y_)
			+ 3. * (vxyz - vxyz_ + vxy_z - vxy_z_ +
				vx_yz - vx_yz_ + vx_y_z - vx_y_z_)
			- 8. * (vxy - vxy_ + vx_y - vx_y_ +
				vx_z - vx_z_ + vxz - vxz_)
			+ 4. * (vz - vz_)
			) * h3 * (1.0 / 12);
		return VectorValue<T, 3>{ x, y, z };
	}
};

template<>
struct apply_gradlaplacian_3d2h<40>
{
	template<typename T>
	auto operator()(T* const v, double h3, len_type lenX, len_type lenY)
	{
		auto x = (16. * (vx2z - vx2_z + vx2z_ - vx2_z_ +
			vx2y - vx2_y + vx2y_ - vx2_y_)
			- 16. * (vx2 - vx2_)
			+ 3. * (vxy2z - vx_y2z + vxy2_z - vx_y2_z +
				vxy2z_ - vx_y2z_ + vxy2_z_ - vx_y2_z_ +
				vxyz2 - vx_yz2 + vxy_z2 - vx_y_z2 +
				vxyz2_ - vx_yz2_ + vxy_z2_ - vx_y_z2_)
			+ 2. * (vxz2 - vx_z2 + vxz2_ - vx_z2_ +
				vxy2 - vx_y2 + vxy2_ - vx_y2_)
			- 22. * (vxy - vx_y + vxy_ - vx_y_ +
				vxz_ - vx_z_ + vxz - vx_z)
			- 40. * (vx - vx_)
			) * h3 * (1.0 / 96);
		auto y = (16. * (vy2z - vy2_z + vy2z_ - vy2_z_ +
			vxy2 - vxy2_ + vx_y2 - vx_y2_)
			- 16. * (vy2 - vy2_)
			+ 3. * (vx2yz - vx2y_z + vx2_yz - vx2_y_z +
				vx2yz_ - vx2y_z_ + vx2_yz_ - vx2_y_z_ +
				vxyz2 - vxy_z2 + vx_yz2 - vx_y_z2 +
				vxyz2_ - vxy_z2_ + vx_yz2_ - vx_y_z2_)
			+ 2. * (vyz2 - vy_z2 + vyz2_ - vy_z2_ +
				vx2y - vx2y_ + vx2_y - vx2_y_)
			- 22. * (vxy - vxy_ + vx_y - vx_y_ +
				vyz_ - vy_z_ + vyz - vy_z)
			- 40. * (vy - vy_)
			) * h3 * (1.0 / 96);
		auto z = (16. * (vxz2 - vxz2_ + vx_z2 - vx_z2_ +
			vxy2 - vxy2_ + vx_y2 - vx_y2_)
			- 16. * (vz2 - vz2_)
			+ 3. * (vxy2z - vxy2z_ + vxy2_z - vxy2_z_ +
				vx_y2z - vx_y2z_ + vx_y2_z - vx_y2_z_ +
				vx2yz - vx2yz_ + vx2y_z - vx2y_z_ +
				vx2_yz - vx2_yz_ + vx2_y_z - vx2_y_z_)
			+ 2. * (vx2z - vx2z_ + vx2_z - vx2_z_ +
				vx2y - vx2y_ + vx2_y - vx2_y_)
			- 22. * (vxy - vxy_ + vx_y - vx_y_ +
				vx_z - vx_z_ + vxz - vxz_)
			- 40. * (vz - vz_)
			) * h3 * (1.0 / 96);
		return VectorValue<T, 3>{ x, y, z };
	}
};





