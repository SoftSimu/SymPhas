
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
 * PURPOSE: Defines all the finite difference stencils which will be available 
 * to the solvers to approximate the derivatives. These stencils assume
 * that the grid is uniform with equal discretization along all axes.
 *
 * ***************************************************************************
 */



#pragma once


#include "stencilh2.h"
#include "stencilh4.h"


/*!
 * \addtogroup stencil
 * @{
 */

//! Characterizes finite difference stencils.
/*!
 * Defines the stencil sizes of the Laplacian, bilaplacian
 * and gradlaplacian stencils, as well as their accuracy. Used
 * to characterize the stencil configuration used for the solution.
 */
struct StencilParams
{
	StencilParams() : ord{ 2 }, ptl{ 5 }, ptb{ 13 }, ptg{ 6 } {}

	size_t
		ord,	//!< Order of accuracy for the stencils.
		ptl,	//!< Number of cells involved in the Laplacian stencil.
		ptb,	//!< Number of cells involved in the bilaplacian stencil.
		ptg;	//!< Number of cells involved in the gradlaplacian stencil.
};

namespace symphas::internal
{

	//! Implements the gradient for the 1-dimensional stencil.
	/*!
	 * 1-dimensional implementation for gradient with 2nd order of accuracy.
	 */
	struct StencilBase1d2h : GeneralizedStencil<1, 2>
	{
		using parent_type = GeneralizedStencil<1, 2>;
		using parent_type::dims;

		//! Construct a new stencil from the system dimensions and \f$h\f$.
		StencilBase1d2h(const len_type* dims, double h) : parent_type(dims, h),
			divh2{ divh * divh }, divh3{ divh * divh2 }, divh4{ divh2 * divh2 } {}

		double divh2, divh3, divh4;

		//! Gradient of the field.
		template<typename T>
		auto gradient(T* const) const;
	};


	//! Implements the gradient for the 2-dimensional stencil.
	/*!
	 * 2-dimensional implementation for gradient with 2nd order of accuracy.
	 */
	struct StencilBase2d2h : GeneralizedStencil<2, 2>
	{
		using parent_type = GeneralizedStencil<2, 2>;
		using parent_type::dims;

		//! Construct a new stencil from the system dimensions and \f$h\f$.
		StencilBase2d2h(const len_type* dims, double h) : parent_type(dims, h),
			divh2{ divh * divh }, divh3{ divh * divh2 }, divh4{ divh2 * divh2 } {}

		double divh2, divh3, divh4;

		//! Gradient of the field.
		template<typename T>
		auto gradient(T* const) const;
	};


	//! Implements the gradient for the 3-dimensional stencil.
	/*!
	 * 3-dimensional implementation for gradient with 2nd order of accuracy.
	 */
	struct StencilBase3d2h : GeneralizedStencil<3, 2>
	{
		using parent_type = GeneralizedStencil<3, 2>;
		using parent_type::dims;

		//! Construct a new stencil from the system dimensions and \f$h\f$.
		StencilBase3d2h(const len_type* dims, double h) : parent_type(dims, h),
			divh2{ divh * divh }, divh3{ divh * divh2 }, divh4{ divh2 * divh2 } {}

		double divh2, divh3, divh4;

		//! Gradient of the field.
		template<typename T>
		auto gradient(T* const) const;
	};



	template<>
	inline auto StencilBase1d2h::gradient<scalar_t>(scalar_t* const v) const
	{
		return VectorValue<scalar_t, 1>{ (vx_ + vx - 2. * v0) * divh };
	}

	template<>
	inline auto StencilBase2d2h::gradient<scalar_t>(scalar_t* const v) const
	{
		const len_type& lenX = dims[0];
		return VectorValue<scalar_t, 2>{
			(vx_ + vx - 2. * v0) * divh,
				(vy_ + vy - 2. * v0) * divh };
	}

	template<>
	inline auto StencilBase3d2h::gradient<scalar_t>(scalar_t* const v) const
	{
		const len_type& lenX = dims[0];
		const len_type& lenY = dims[1];
		return VectorValue<scalar_t, 3>{
			(vx_ + vx - 2. * v0) * divh,
				(vy_ + vy - 2. * v0) * divh,
				(vz_ + vz - 2. * v0) * divh };
	}


	template<typename T>
	auto StencilBase1d2h::gradient(T* const) const
	{
		fprintf(SYMPHAS_ERR, "using gradient on incompatible type.\n");
		exit(2);
	}

	template<typename T>
	auto StencilBase2d2h::gradient(T* const) const
	{
		fprintf(SYMPHAS_ERR, "using gradient on incompatible type.\n");
		exit(2);
	}

	template<typename T>
	auto StencilBase3d2h::gradient(T* const) const
	{
		fprintf(SYMPHAS_ERR, "using gradient on incompatible type.\n");
		exit(2);
	}



	//! Implements the gradient for the 2-dimensional stencil.
	/*!
	 * 2-dimensional implementation for gradient with 4th order of accuracy.
	 */
	struct StencilBase2d4h : GeneralizedStencil<2, 4>
	{
		using parent_type = GeneralizedStencil<2, 4>;
		using parent_type::dims;

		//! Construct a new stencil from the system dimensions and \f$h\f$.
		StencilBase2d4h(const len_type* dims, double h) : parent_type(dims, h),
			divh2{ divh * divh }, divh3{ divh * divh2 }, divh4{ divh2 * divh2 } {}

		double divh2, divh3, divh4;

		//! Gradient of the field.
		template<typename T>
		auto gradient(T* const) const;
	};

	template<>
	inline auto StencilBase2d4h::gradient<scalar_t>(scalar_t* const v) const
	{
		const len_type& lenX = dims[0];
		return VectorValue<scalar_t, 2>{
			-(vx_ + vx - 2. * v0) * divh,
				-(vy_ + vy - 2. * v0) * divh };
	}

	template<typename T>
	auto StencilBase2d4h::gradient(T* const) const
	{
		fprintf(SYMPHAS_ERR, "using gradient on incompatible type.\n");
		exit(2);
	}
}

//! 1-dimensional stencil with 2nd order of accuracy.
/*!
 * Implements the 1-dimensional stencil of 2nd order accuracy for all derivative
 * orders up to fourth order. For higher orders, the generalized stencil is
 * applied.
 *
 * \tparam L The number of points for the laplacian.
 * \tparam B The number of points for the bilaplacian.
 * \tparam G The number of points for the gradlaplacian.
 */
template <size_t L, size_t B, size_t G>
struct Stencil1d2h : symphas::internal::StencilBase1d2h, Stencil<Stencil1d2h<L, B, G>>
{
	using parent_type = Stencil<Stencil1d2h<L, B, G>>;
	using base_type = symphas::internal::StencilBase1d2h;
	
	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	auto laplacian(T* const v) const
	{
		return apply_laplacian_1d2h<L>{}(v, divh2);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	auto bilaplacian(T* const v) const
	{
		return apply_bilaplacian_1d2h<B>{}(v, divh4);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	auto gradlaplacian(T* const v) const
	{
		return apply_gradlaplacian_1d2h<G>{}(v, divh3);
	}

};


//! 2-dimensional stencil with 2nd order of accuracy.
/*!
 * Implements the 2-dimensional stencil of 2nd order accuracy for all derivative
 * orders up to fourth order. For higher orders, the generalized stencil is
 * applied.
 *
 * \tparam L The number of points for the laplacian.
 * \tparam B The number of points for the bilaplacian.
 * \tparam G The number of points for the gradlaplacian.
 */
template <size_t L, size_t B, size_t G>
struct Stencil2d2h : symphas::internal::StencilBase2d2h, Stencil<Stencil2d2h<L, B, G>>
{
	using base_type = symphas::internal::StencilBase2d2h;
	using parent_type = Stencil<Stencil2d2h<L, B, G>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v) const
	{
		return apply_laplacian_2d2h<L>{}(v, divh2, dims[0]);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v) const
	{
		return apply_bilaplacian_2d2h<B>{}(v, divh4, dims[0]);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v) const
	{
		return apply_gradlaplacian_2d2h<G>{}(v, divh3, dims[0]);
	}

};

//! 3-dimensional stencil with 2nd order of accuracy.
/*!
 * Implements the 3-dimensional stencil of 2nd order accuracy for all derivative
 * orders up to fourth order. For higher orders, the generalized stencil is
 * applied.
 *
 * \tparam L The number of points for the laplacian.
 * \tparam B The number of points for the bilaplacian.
 * \tparam G The number of points for the gradlaplacian.
 */
template <size_t L, size_t B, size_t G>
struct Stencil3d2h : symphas::internal::StencilBase3d2h, Stencil<Stencil3d2h<L, B, G>>
{
	using base_type = symphas::internal::StencilBase3d2h;
	using parent_type = Stencil<Stencil3d2h<L, B, G>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v) const
	{
		return apply_laplacian_3d2h<L>{}(v, divh2, dims[0], dims[1]);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v) const
	{
		return apply_bilaplacian_3d2h<B>{}(v, divh4, dims[0], dims[1]);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v) const
	{
		return apply_gradlaplacian_3d2h<G>{}(v, divh3, dims[0], dims[1]);
	}
};




//! 2-dimensional stencil with 4th order of accuracy.
/*!
 * Implements the 2-dimensional stencil of 4th order accuracy for all derivative
 * orders up to fourth order. For higher orders, the generalized stencil is
 * applied.
 *
 * \tparam L The number of points for the laplacian.
 * \tparam B The number of points for the bilaplacian.
 * \tparam G The number of points for the gradlaplacian.
 */
template <size_t L, size_t B, size_t G>
struct Stencil2d4h : symphas::internal::StencilBase2d4h, Stencil<Stencil2d4h<L, B, G>>
{
	using base_type = symphas::internal::StencilBase2d4h;
	using parent_type = Stencil<Stencil2d4h<L, B, G>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;


	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v) const
	{
		return apply_laplacian_2d4h<L>{}(v, divh3, dims[0]);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v) const
	{
		return apply_bilaplacian_2d4h<B>{}(v, divh3, dims[0]);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v) const
	{
		return apply_gradlaplacian_2d4h<G>{}(v, divh3, dims[0]);
	}
};

//! @}


#undef v0
#undef vx
#undef vx_
#undef vy
#undef vy_
#undef vz
#undef vz_
#undef vx2
#undef vx2_
#undef vy2
#undef vy2_
#undef vz2
#undef vz2_
#undef vx3
#undef vx3_
#undef vy3
#undef vy3_
#undef vz3
#undef vz3_
#undef vxy
#undef vxy_
#undef vx_y
#undef vx_y_
#undef vx2y
#undef vx2y_
#undef vx2_y
#undef vx2_y_
#undef vxy2
#undef vxy2_
#undef vx_y2
#undef vx_y2_
#undef vx3y
#undef vx3y_
#undef vx3_y
#undef vx3_y_
#undef vxy3
#undef vxy3_
#undef vx_y3
#undef vx_y3_
#undef vx2y2
#undef vx2y2_
#undef vx2_y2
#undef vx2_y2_
#undef vx3y3
#undef vx3y3_
#undef vx3_y3
#undef vx3_y3_
#undef vxz
#undef vxz_
#undef vx_z
#undef vx_z_
#undef vx2z
#undef vx2z_
#undef vx2_z
#undef vx2_z_
#undef vxz2
#undef vxz2_
#undef vx_z2
#undef vx_z2_
#undef vx2z2
#undef vx2z2_
#undef vx2_z2
#undef vx2_z2_
#undef vyz
#undef vyz_
#undef vy_z
#undef vy_z_
#undef vyz2
#undef vyz2_
#undef vy_z2
#undef vy_z2_
#undef vy2z
#undef vy2z_
#undef vy2_z
#undef vy2_z_
#undef vy2z2
#undef vy2z2_
#undef vy2_z2
#undef vy2_z2_
#undef vxyz
#undef vxyz_
#undef vxy_z
#undef vxy_z_
#undef vx_yz
#undef vx_yz_
#undef vx_y_z
#undef vx_y_z_
#undef vx2yz
#undef vx2yz_
#undef vx2y_z
#undef vx2y_z_
#undef vx2_yz
#undef vx2_yz_
#undef vx2_y_z
#undef vx2_y_z_
#undef vxy2z
#undef vxy2z_
#undef vxy2_z
#undef vxy2_z_
#undef vx_y2z
#undef vx_y2z_
#undef vx_y2_z
#undef vx_y2_z_
#undef vxyz2
#undef vxyz2_
#undef vxy_z2
#undef vxy_z2_
#undef vx_yz2
#undef vx_yz2_
#undef vx_y_z2
#undef vx_y_z2_
#undef vx2y2z2
#undef vx2y2z2_
#undef vx2y2_z2
#undef vx2y2_z2_
#undef vx2_y2z2
#undef vx2_y2z2_
#undef vx2_y2_z2
#undef vx2_y2_z2_



