
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

#include "gridinfo.h"
#include "expressionstencils.h"

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

	template <typename T, size_t D>
	auto laplacian(T* const v, const len_type(&stride)[D]) const
	{
		return cast().laplacian(v, stride);
	}
	template <typename T, size_t D>
	auto bilaplacian(T* const v, const len_type(&stride)[D]) const
	{
		return cast().bilaplacian(v, stride);
	}
	template <typename T, size_t D>
	auto gradlaplacian(T* const v, const len_type(&stride)[D]) const
	{
		return cast().gradlaplacian(v, stride);
	}

	template<typename T, size_t D>
	auto gradient(T* const v, const len_type (&stride)[D]) const
	{
		return cast().gradient(v, stride);
	}

	template<Axis ax, size_t O, template<typename, size_t> typename G, typename T, size_t D, 
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_generalized_directional_derivative(G<T, D> const& grid, iter_type n) const
	{
		return cast().template apply_directional<ax, O>(grid.values + n);
	}

	template<size_t... Os, template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_generalized_mixed_derivative(G<T, D> const& grid, iter_type n) const
	{
		return cast().template apply_mixed<Os...>(grid.values + n);
	}

	template<Axis ax, size_t O, template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_generalized_derivative(G<T, D> const& grid, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return cast().template apply<O>(grid.values + n, stride);
	}

	template<template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_laplacian(G<T, D> const& grid, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<Axis::X>(stride, cast().dims);
		return laplacian(grid.values + n, stride);
	}

	template<template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_bilaplacian(G<T, D> const& grid, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<Axis::X>(stride, cast().dims);
		return bilaplacian(grid.values + n, stride);
	}

	template<Axis ax, template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_gradlaplacian(G<T, D> const& grid, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return gradlaplacian(grid.values + n, stride);
	}

	template<Axis ax, template<typename, size_t> typename G, typename T, size_t D,
		typename = std::enable_if_t<std::is_convertible<G<T, D>, Block<T>>::value, int>>
	auto applied_gradient(G<T, D> const& grid, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return gradient(grid.values + n, stride);
	}


	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_directional_derivative(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			return cast().template apply_directional<ax, O>(v.value);
		}
		else
		{
			return *v.value;
		}
	}

	template<size_t... Os, typename T, size_t D>
	auto applied_generalized_mixed_derivative(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			return cast().template apply_mixed<Os...>(v.value);
		}
		else
		{
			return *v.value;
		}
	}

	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_derivative(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			len_type stride[D];
			grid::get_stride<ax>(stride, grid.region.dims);
			return cast().template apply<O>(v.value, stride);
		}
		else
		{
			return *v.value;
		}
	}

	template<typename T, size_t D>
	auto applied_laplacian(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			return laplacian(v.value, grid.region.stride);
		}
		else
		{
			return *v.value;
		}
	}

	template<typename T, size_t D>
	auto applied_bilaplacian(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			return bilaplacian(v.value, grid.region.stride);
		}
		else
		{
			return *v.value;
		}
	}

	template<Axis ax, typename T, size_t D>
	auto applied_gradlaplacian(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			len_type stride[D];
			grid::get_stride<ax>(stride, grid.region.dims);
			return gradlaplacian(v.value, stride);
		}
		else
		{
			return *v.value;
		}
	}


	template<Axis ax, typename T, size_t D>
	auto applied_gradient(RegionalGrid<T, D> const& grid, iter_type n) const
	{
		auto v = grid[n];
		if (v.is_valid())
		{
			len_type stride[D];
			grid::get_stride<ax>(stride, grid.region.dims);
			return gradient(v.value, stride);
		}
		else
		{
			return *v.value;
		}
	}


	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_directional_derivative(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		return cast().template apply_directional<ax, O>(data.values + n);
	}

	template<size_t... Os, Axis ax, typename T, size_t D>
	auto applied_generalized_mixed_derivative(VectorComponentData<ax, T*, D> const& grid, iter_type n) const
	{
		return cast().template apply_mixed<Os...>(grid.values + n);
	}

	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_derivative(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return cast().template apply<O>(data.values + n, stride);
	}

	template<Axis ax, typename T, size_t D>
	auto applied_laplacian(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<Axis::X>(stride, cast().dims);
		return laplacian(data.values + n, stride);
	}

	template<Axis ax, typename T, size_t D>
	auto applied_bilaplacian(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<Axis::X>(stride, cast().dims);
		return bilaplacian(data.values + n, stride);
	}

	template<Axis ax, typename T, size_t D>
	auto applied_gradlaplacian(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return gradlaplacian(data.values + n, stride);
	}

	template<Axis ax, typename T, size_t D>
	auto applied_gradient(VectorComponentData<ax, T*, D> const& data, iter_type n) const
	{
		len_type stride[D];
		grid::get_stride<ax>(stride, cast().dims);
		return gradient(data.values + n, stride);
	}


	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_directional_derivative(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<size_t... Os, typename T, size_t D>
	auto applied_generalized_mixed_derivative(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<Axis ax, size_t O, typename T, size_t D>
	auto applied_generalized_derivative(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<typename T, size_t D>
	auto applied_laplacian(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<typename T, size_t D>
	auto applied_bilaplacian(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<Axis ax, typename T, size_t D>
	auto applied_gradlaplacian(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}

	template<Axis ax, typename T, size_t D>
	auto applied_gradient(GridSymbol<T, D> const& grid, iter_type n) const
	{
		return T{};
	}


	template<Axis ax, size_t O>
	auto applied_generalized_directional_derivative(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	template<size_t... Os>
	auto applied_generalized_mixed_derivative(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	template<Axis ax, size_t O>
	auto applied_generalized_derivative(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	auto applied_laplacian(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	auto applied_bilaplacian(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	template<Axis ax>
	auto applied_gradlaplacian(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	template<Axis ax>
	auto applied_gradient(...) const
	{
		return symphas::lib::get_identity<scalar_t>();
	}

	const Sp& cast() const
	{
		return (*static_cast<const Sp*>(this));
	}

	Sp& cast()
	{
		return (*static_cast<Sp*>(this));
	}
};

namespace symphas::internal
{

#ifdef GENERATE_UNDEFINED_STENCILS_ON

	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type (&stride)[1])
	{
		static symphas::internal::GeneratedStencilApply stencil{ expr::get_central_space_stencil<OD, OA, 1>() };
		return stencil(v, stride, divh);
	}

	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type(&stride)[2])
	{
		static symphas::internal::GeneratedStencilApply stencil{ expr::get_central_space_stencil<OD, OA, 2>() };
		return stencil(v, stride, divh);
	}

	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type(&stride)[3])
	{
		static symphas::internal::GeneratedStencilApply stencil{ expr::get_central_space_stencil<OD, OA, 3>() };
		return stencil(v, stride, divh);
	}

#else

	struct no_derivative_message_printed
	{
		no_derivative_message_printed(size_t OD, size_t OA, size_t D)
		{
			fprintf(SYMPHAS_ERR, "no derivative of order %zd accuracy %zd available in dimension %zd\n", OD, OA, D);
		}
	};

	template<size_t OD, size_t OA, size_t D>
	void print_no_derivative_message()
	{
		static no_derivative_message_printed message{ OD, OA, D };
	}


	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type(&stride)[1])
	{
		print_no_derivative_message<OD, OA, 1>();
		return T{};
	}

	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type(&stride)[2])
	{
		print_no_derivative_message<OD, OA, 2>();
		return T{};
	}

	template<size_t OD, size_t OA, typename T>
	auto apply_generalized_derivative(T* const v, double divh, const len_type(&stride)[3])
	{
		print_no_derivative_message<OD, OA, 3>();
		return T{};
	}

#endif
}


namespace symphas::internal
{
}

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

    GeneralizedStencil() : dims{}, divh{ 0 } {}

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
	auto apply(T* const v, const len_type(&stride)[DD]) const
	{
		return symphas::internal::apply_generalized_derivative<OD, OA>(v, divh, stride);
	}

	//! Determine the finite difference approximation to the directional derivative.
	/*!
	 * Determine the finite difference approximation to the directional derivative of
	 * the given order.
	 *
	 * \tparam OD Order of the derivative.
	 * \tparam ax Axis of the directional derivative.
	 */
	template<Axis ax, size_t OD, typename T/*, typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_stencil<OD, OA, 1>)>,
		typename = std::enable_if_t<((ax == Axis::Z) ? (DD == 3) : (ax == Axis::Y) ? (DD >= 2) : (ax == Axis::X) ? (DD >= 1) : false), int>*/>
	auto apply_directional(T* const v) const
	{
		len_type stride = (ax == Axis::Z) ? (dims[0] * dims[1]) : (ax == Axis::Y) ? dims[0] : 1;
		static symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_stencil<OD, OA, 1>());

#ifdef DEBUG
		static size_t n = expr::print_stencil(stencil);
#endif

		return stencil(v, stride, divh);

	}

	//! Determine the finite difference approximation to the mixed derivative.
	/*!
	 * Determine the finite difference approximation to the mixed derivative of
	 * the given order.
	 *
	 * \tparam OD Order of the derivative.
	 * \tparam ax Axis of the directional derivative.
	 */
	template<size_t OD1, size_t OD2, size_t OD3, typename T/*,
		typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_mixed_stencil<OA, OD1, OD2, OD3>), std::index_sequence<OD1, OD2, OD3>>*/>
	auto apply_mixed(T* const v) const
	{
		len_type stride[DD];
		grid::get_stride<Axis::X>(stride, dims);
		static symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_mixed_stencil<OA>(std::index_sequence<OD1, OD2, OD3>{}));

#ifdef DEBUG
		static size_t n = expr::print_stencil(stencil);
#endif

		return stencil(v, stride, divh);
		//return symphas::internal::GeneratedStencilApply<Stt>{ stride, divh }(v);
	}

	//! Determine the finite difference approximation to the mixed derivative.
	/*!
	 * Determine the finite difference approximation to the mixed derivative of
	 * the given order.
	 *
	 * \tparam OD Order of the derivative.
	 * \tparam ax Axis of the directional derivative.
	 */
	template<size_t OD1, size_t OD2, typename T/*,
		typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_mixed_stencil<OA, OD1, OD2>), std::index_sequence<OD1, OD2>>*/>
	auto apply_mixed(T* const v) const
	{
		len_type stride[DD];
		grid::get_stride<Axis::X>(stride, dims);
		static symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_mixed_stencil<OA>(std::index_sequence<OD1, OD2>{}));

#ifdef DEBUG
		static size_t n = expr::print_stencil(stencil);
#endif

		return stencil(v, stride, divh);
		//return symphas::internal::GeneratedStencilApply<Stt>{ stride, divh }(v);
	}

	//! Determine the finite difference approximation to the mixed derivative.
	/*!
	 * Determine the finite difference approximation to the mixed derivative of
	 * the given order.
	 *
	 * \tparam OD Order of the derivative.
	 * \tparam ax Axis of the directional derivative.
	 */
	template<size_t OD1, typename T/*, typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_stencil<OD1, OA, 1>)>*/>
	auto apply_mixed(T* const v) const
	{
		len_type stride[DD];
		grid::get_stride<Axis::X>(stride, dims);
		static symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_stencil<OD1, OA, 1>());

#ifdef DEBUG
		static size_t n = expr::print_stencil(stencil);
#endif

		return stencil(v, stride, divh);
		//return symphas::internal::GeneratedStencilApply<Stt>{ stride, divh }(v);
	}


	template<typename T>
	auto gradient(T* const v, const len_type(&stride)[DD]) const
	{
		return apply<1>(v, stride);
	}

	template<typename T>
	auto gradient(T* const v) const
	{
		len_type stride[DD]{};
		grid::get_stride<Axis::X>(stride, dims);
		return apply<1>(v, stride);
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

#define __DZ stride[2]
#define __DY stride[1]
#define __DX stride[0]

#define __2DZ (2 * __DZ)
#define __2DY (2 * __DY)
#define __2DX (2 * __DX)

#define __3DZ (3 * __DZ)
#define __3DY (3 * __DY)
#define __3DX (3 * __DX)


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

