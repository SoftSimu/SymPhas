
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
	constexpr StencilParams(unsigned short ord, unsigned short ptl, unsigned short ptg, unsigned short ptb)
		: ord{ ord }, ptl{ ptl }, ptg{ ptg }, ptb{ ptb } {}
	StencilParams() : ord{ 0 }, ptl{ 0 }, ptg{ 0 }, ptb{ 0 } {}

	unsigned short
		ord,	//!< Order of accuracy for the stencils.
		ptl,	//!< Number of cells involved in the Laplacian stencil.
		ptg,	//!< Number of cells involved in the gradlaplacian stencil.
		ptb;	//!< Number of cells involved in the bilaplacian stencil.

	constexpr operator size_t() const
	{
		return (ord & 7ull)
			| ((ptl & 63ull) << 3ull)
			| ((ptg & 63ull) << 9ull)
			| ((ptb & 63ull) << 15ull);
	}
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
		auto gradient(T* const, const len_type(&)[1]) const = delete;
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
		auto gradient(T* const, const len_type(&)[2]) const = delete;
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
		auto gradient(T* const, const len_type(&)[3]) const = delete;
	};



	template<>
	inline auto StencilBase1d2h::gradient<scalar_t>(scalar_t* const v, const len_type (&stride)[1]) const
	{
		return -0.5 * (vx_ - vx) * divh;
	}

	template<>
	inline auto StencilBase2d2h::gradient<scalar_t>(scalar_t* const v, const len_type (&stride)[2]) const
	{
		return -0.5 * (vx_ - vx) * divh;
	}

	template<>
	inline auto StencilBase3d2h::gradient<scalar_t>(scalar_t* const v, const len_type (&stride)[3]) const
	{
		return -0.5 * (vx_ - vx) * divh;
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
		auto gradient(T* const, const len_type(&)[2]) const;
	};

	template<>
	inline auto StencilBase2d4h::gradient<scalar_t>(scalar_t* const v, const len_type (&stride)[2]) const
	{
		return (1 / 12.) * ((vx2_ - vx2) - 8 * (vx_ - vx)) * divh;
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
template <size_t L, size_t G, size_t B>
struct Stencil1d2h : symphas::internal::StencilBase1d2h, Stencil<Stencil1d2h<L, G, B>>
{
	using parent_type = Stencil<Stencil1d2h<L, G, B>>;
	using base_type = symphas::internal::StencilBase1d2h;
	
	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	auto laplacian(T* const v, const len_type (&stride)[1]) const
	{
		return apply_laplacian_1d2h<L>{}(v, divh2, stride);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	auto bilaplacian(T* const v, const len_type (&stride)[1]) const
	{
		return apply_bilaplacian_1d2h<B>{}(v, divh4, stride);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	auto gradlaplacian(T* const v, const len_type (&stride)[1]) const
	{
		return apply_gradlaplacian_1d2h<G>{}(v, divh3, stride);
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
template <size_t L, size_t G, size_t B>
struct Stencil2d2h : symphas::internal::StencilBase2d2h, Stencil<Stencil2d2h<L, G, B>>
{
	using base_type = symphas::internal::StencilBase2d2h;
	using parent_type = Stencil<Stencil2d2h<L, G, B>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v, const len_type(&stride)[2]) const
	{
		return apply_laplacian_2d2h<L>{}(v, divh2, stride);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v, const len_type(&stride)[2]) const
	{
		return apply_bilaplacian_2d2h<B>{}(v, divh4, stride);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v, const len_type (&stride)[2]) const
	{
		return apply_gradlaplacian_2d2h<G>{}(v, divh3, stride);
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
template <size_t L, size_t G, size_t B>
struct Stencil3d2h : symphas::internal::StencilBase3d2h, Stencil<Stencil3d2h<L, G, B>>
{
	using base_type = symphas::internal::StencilBase3d2h;
	using parent_type = Stencil<Stencil3d2h<L, G, B>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;

	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v, const len_type (&stride)[3]) const
	{
		return apply_laplacian_3d2h<L>{}(v, divh2, stride);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v, const len_type (&stride)[3]) const
	{
		return apply_bilaplacian_3d2h<B>{}(v, divh4, stride);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v, const len_type (&stride)[3]) const
	{
		return apply_gradlaplacian_3d2h<G>{}(v, divh3, stride);
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
template <size_t L, size_t G, size_t B>
struct Stencil2d4h : symphas::internal::StencilBase2d4h, Stencil<Stencil2d4h<L, G, B>>
{
	using base_type = symphas::internal::StencilBase2d4h;
	using parent_type = Stencil<Stencil2d4h<L, G, B>>;

	using base_type::base_type;
	using base_type::gradient;
	using base_type::dims;


	//! Laplacian (2nd order derivative) of the field.
	template<typename T>
	inline auto laplacian(T* const v, const len_type (&stride)[2]) const
	{
		return apply_laplacian_2d4h<L>{}(v, divh3, stride);
	}

	//! Bilaplacian (4th order derivative) of the field.
	template<typename T>
	inline auto bilaplacian(T* const v, const len_type (&stride)[2]) const
	{
		return apply_bilaplacian_2d4h<B>{}(v, divh3, stride);
	}

	//! Gradlaplacian (gradient of the laplacian) of the field.
	template<typename T>
	inline auto gradlaplacian(T* const v, const len_type (&stride)[2]) const
	{
		return apply_gradlaplacian_2d4h<G>{}(v, divh3, stride);
	}
};

template<size_t DD, size_t OA = 2>
struct SelfSelectingStencil : GeneralizedStencil<DD, OA>, Stencil<SelfSelectingStencil<DD, OA>>
{
	using base_type = GeneralizedStencil<DD, OA>;
	using base_type::base_type;
};



template<>
struct SelfSelectingStencil<1, 2> : GeneralizedStencil<1, 2>, Stencil<SelfSelectingStencil<1, 2>>
{
	using base_type = GeneralizedStencil<1, 2>;
	using base_type::base_type;
	
	template<size_t... Ps>
	using Points = Stencil1d2h<Ps...>;
};

template<>
struct SelfSelectingStencil<2, 2> : GeneralizedStencil<2, 2>, Stencil<SelfSelectingStencil<2, 2>>
{
	using base_type = GeneralizedStencil<2, 2>;
	using base_type::base_type;

	template<size_t... Ps>
	using Points = Stencil2d2h<Ps...>;
};

template<>
struct SelfSelectingStencil<2, 4> : GeneralizedStencil<2, 4>, Stencil<SelfSelectingStencil<2, 4>>
{
	using base_type = GeneralizedStencil<2, 4>;
	using base_type::base_type;

	template<size_t... Ps>
	using Points = Stencil2d4h<Ps...>;
};


template<>
struct SelfSelectingStencil<3, 2> : GeneralizedStencil<3, 2>, Stencil<SelfSelectingStencil<3, 2>>
{
	using base_type = GeneralizedStencil<3, 2>;
	using base_type::base_type;

	template<size_t... Ps>
	using Points = Stencil3d2h<Ps...>;
};


//! @}

namespace symphas::internal
{

	/*!
	 * Defines the "point number" for the given finite difference approximation,
	 * based on the index of the derivative that is approximated, as well as its
	 * dimension and order of accuracy.
	 *
	 * \tparam N The order of the derivative for the point list.
	 * \tparam D Dimension of derivative.
	 * \tparam O The order of accuracy of the derivative.
	 */
	template<size_t N, size_t D, size_t O>
	struct StencilPointList;

	template<size_t D>
	struct OrderList
	{
		using type = std::index_sequence<2>;
	};
}


// **************************************************************************************


#define MAKE_STENCIL_POINT_LIST(N, DIM, ORD, PS) \
template<> struct symphas::internal::StencilPointList<N, DIM, ORD> { using type = std::index_sequence<SINGLE_ARG PS>; };

#define MAKE_AVAILABLE_ORDER_LIST(DIM, ORDS) \
template<> struct symphas::internal::OrderList<DIM> { using type = std::index_sequence<SINGLE_ARG ORDS>; };

#if !defined(DEBUG) && defined(ALL_STENCILS)

MAKE_STENCIL_POINT_LIST(2, 1, 2, (3))
MAKE_STENCIL_POINT_LIST(4, 1, 2, (5))
MAKE_STENCIL_POINT_LIST(3, 1, 2, (4))

MAKE_STENCIL_POINT_LIST(2, 2, 2, (5, 9))
MAKE_STENCIL_POINT_LIST(4, 2, 2, (13, 17, 21))
MAKE_STENCIL_POINT_LIST(3, 2, 2, (6, 8, 12, 16))
MAKE_STENCIL_POINT_LIST(2, 2, 4, (9, 17, 21))
MAKE_STENCIL_POINT_LIST(4, 2, 4, (21, 25, 33, 37))
MAKE_STENCIL_POINT_LIST(3, 2, 4, (14, 18, 26, 30))

MAKE_STENCIL_POINT_LIST(2, 3, 2, (7, 15, 19, 21, 27))
MAKE_STENCIL_POINT_LIST(4, 3, 2, (21, 25, 41, 52, 57))
MAKE_STENCIL_POINT_LIST(3, 3, 2, (10, 12, 28, 36, 40))

MAKE_AVAILABLE_ORDER_LIST(1, (2))
MAKE_AVAILABLE_ORDER_LIST(2, (2, 4))
MAKE_AVAILABLE_ORDER_LIST(3, (2))

#else

MAKE_STENCIL_POINT_LIST(2, 1, 2, (3))
MAKE_STENCIL_POINT_LIST(4, 1, 2, (5))
MAKE_STENCIL_POINT_LIST(3, 1, 2, (4))

MAKE_STENCIL_POINT_LIST(2, 2, 2, (9))
MAKE_STENCIL_POINT_LIST(4, 2, 2, (13))
MAKE_STENCIL_POINT_LIST(3, 2, 2, (6))

MAKE_STENCIL_POINT_LIST(2, 3, 2, (7))
MAKE_STENCIL_POINT_LIST(4, 3, 2, (21))
MAKE_STENCIL_POINT_LIST(3, 3, 2, (10))

#endif

#define AVAILABLE_DIMENSIONS 2


namespace symphas::lib::internal
{
	template<size_t I, typename Seq0, typename... Seqs>
	constexpr size_t get_search_offset(CrossProductList<Seq0, Seqs...>)
	{
		using cl_type = CrossProductList<Seq0, Seqs...>;
		if constexpr (I + 1 >= cl_type::rank)
		{
			return 1;
		}
		else
		{
			return cl_type::template size<I + 1> *get_search_offset<I + 1>(CrossProductList<Seq0, Seqs...>{});
		}
	}

	template<typename>
	struct SearchCrossListReturn
	{
		bool operator()()
		{
			return false;
		}
	};

	template<size_t... As>
	struct SearchCrossListReturn<std::index_sequence<As...>>
	{
		bool operator()()
		{
			return true;
		}
	};

	template<size_t I0, size_t I, size_t L>
	auto matches(const size_t(&)[L], std::index_sequence<>)
	{
		return false;
	}

	template<size_t I0, size_t I, size_t L, size_t V, size_t... Vs>
	auto matches(const size_t(&parameters)[L], std::index_sequence<V, Vs...>)
	{
		if constexpr (I0 == 0)
		{
			return (parameters[I] == V);
		}
		else
		{
			return matches<I0 - 1, I>(parameters, std::index_sequence<Vs...>{});
		}
	}

	template<size_t I, size_t Pos, template<typename> typename F, size_t L,
		typename Seq0, typename... Seqs, typename... Ts,
		typename = std::enable_if_t<(L == CrossProductList<Seq0, Seqs...>::rank), int>>
		auto search(const size_t(&parameters)[L], CrossProductList<Seq0, Seqs...>, Ts&& ...args)
	{
		using cl_type = CrossProductList<Seq0, Seqs...>;

		if constexpr (Pos >= cl_type::count)
		{
			return F<void>{}(std::forward<Ts>(args)...);
		}
		else
		{
			using row_type = typename cl_type::template row<Pos>;

			if constexpr (I >= L)
			{
				return F<row_type>{}(std::forward<Ts>(args)...);
			}
			else
			{
				if (matches<I, I>(parameters, row_type{}))
				{
					return search<I + 1, Pos, F>(parameters, cl_type{}, std::forward<Ts>(args)...);
				}
				else
				{
					constexpr size_t offset = get_search_offset<I>(cl_type{});
					return search<I, Pos + offset, F>(parameters, cl_type{}, std::forward<Ts>(args)...);
				}
			}
		}
	}

	template<template<typename> typename F, size_t L, typename Seq0, typename... Seqs, typename... Ts>
	auto search(const size_t(&parameters)[L], CrossProductList<Seq0, Seqs...>, Ts&& ...args)
	{
		return search<0, 0, F>(parameters, CrossProductList<Seq0, Seqs...>{}, std::forward<Ts>(args)...);
	}

	template<size_t L, typename Seq0, typename... Seqs, typename... Ts>
	auto search(const size_t(&parameters)[L], CrossProductList<Seq0, Seqs...>, Ts&& ...args)
	{
		return search<0, 0, SearchCrossListReturn>(parameters, CrossProductList<Seq0, Seqs...>{}, std::forward<Ts>(args)...);
	}

}

namespace symphas::internal
{

	template<size_t... Ds>
	using dim_ord_list_t = symphas::lib::types_list<
		symphas::lib::types_list<std::index_sequence<Ds>, typename symphas::internal::OrderList<Ds>::type>...>;

	template<size_t D, size_t O>
	using cross_list_t = symphas::lib::CrossProductList<std::index_sequence<D>, std::index_sequence<O>,
		typename symphas::internal::StencilPointList<2, D, O>::type,
		typename symphas::internal::StencilPointList<3, D, O>::type,
		typename symphas::internal::StencilPointList<4, D, O>::type>;

}


/*!
 * \addtogroup stencil
 * @{
 */



struct DefaultStencil
{
protected:

	template<size_t D>
	auto search_ord(std::index_sequence<>) const
	{
		return StencilParams{};
	}

	template<size_t D, size_t O, size_t... Os>
	auto search_ord(std::index_sequence<O, Os...>) const
	{
		if (parameters[1] == O)
		{
			using pl = typename symphas::internal::cross_list_t<D, O>::template row<0>;
			return StencilParams{ O,
				static_cast<unsigned short>(symphas::lib::internal::seq_value<2>(pl{})),
				static_cast<unsigned short>(symphas::lib::internal::seq_value<3>(pl{})),
				static_cast<unsigned short>(symphas::lib::internal::seq_value<4>(pl{})) };
		}
		else
		{
			return search_ord<D>(std::index_sequence<Os...>{});
		}
	}

	template<typename... Ts>
	auto search_dim(symphas::lib::types_list<>) const
	{
		return StencilParams{};
	}

	template<size_t D, size_t... Ns, typename... Seqs>
	auto search_dim(symphas::lib::types_list<symphas::lib::types_list<std::index_sequence<D>, std::index_sequence<Ns...>>, Seqs...>) const
	{
		if (parameters[0] == D)
		{
			return search_ord<D>(std::index_sequence<Ns...>{});
		}
		else
		{
			return search_dim(symphas::lib::types_list<Seqs...>{});
		}
	}

	size_t parameters[2];

public:

	DefaultStencil(size_t dimension, size_t order)
		: parameters{ dimension, order } {}

	auto operator()() const
	{
		return search_dim(symphas::internal::dim_ord_list_t<AVAILABLE_DIMENSIONS>{});
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



