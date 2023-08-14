
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
 * MODULE:  datatypes
 * PURPOSE: Adds macros for iterating over grids, including grids with
 * boundaries.
 *
 * ***************************************************************************
 */

#pragma once

#include "macros.h"
#include "gridinfo.h"



namespace grid
{

	//! Compute the stride for each dimension when getting the index in a flattened grid using position.
	template<size_t D>
	inline void init_region_stride(len_type(&stride)[D], const len_type(&dims)[D], len_type boundary_size = 0)
	{
		stride[0] = 1;
		for (iter_type i = 1; i < D; ++i)
		{
			stride[i] = stride[i - 1] * (dims[i - 1] + boundary_size * 2);
		}
	}

	template<size_t D, size_t... Is>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D], const len_type(&offset)[D], len_type boundary_size, std::index_sequence<Is...>)
	{
		return (((pos[Is] + offset[Is] + boundary_size) * stride[Is]) + ...);
	}

	//! Compute the stride for each dimension when getting the index in a flattened grid using position.
	template<size_t D>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D], const len_type(&offset)[D], len_type boundary_size = 0)
	{
		return index_from_position(pos, stride, offset, boundary_size, std::make_index_sequence<D>{});
	}

	template<size_t D, size_t... Is>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D], len_type boundary_size, std::index_sequence<Is...>)
	{
		return (((pos[Is] + boundary_size) * stride[Is]) + ...);
	}

	//! Compute the stride for each dimension when getting the index in a flattened grid using position.
	template<size_t D>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D], len_type boundary_size)
	{
		return index_from_position(pos, stride, boundary_size, std::make_index_sequence<D>{});
	}

	template<size_t D, size_t... Is>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D], std::index_sequence<Is...>)
	{
		return ((pos[Is] * stride[Is]) + ...);
	}

	//! Compute the stride for each dimension when getting the index in a flattened grid using position.
	template<size_t D>
	inline iter_type index_from_position(const len_type(&pos)[D], const len_type(&stride)[D])
	{
		return index_from_position(pos, stride, std::make_index_sequence<D>{});
	}

	template<size_t D, size_t... Is>
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[D], const len_type(&origin)[D], const len_type(&stride)[D], std::index_sequence<Is...>)
	{
		return ((((n / stride[Is]) % dims[Is] + origin[Is]) * stride[Is]) + ...);
	}

	//! Compute a new index by shifting the current index by origin inside the grid of given dimensions.
	template<size_t D>
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[D], const len_type(&origin)[D], const len_type(&stride)[D])
	{
		return index_from_position(n, dims, origin, stride);
	}

	//! Compute a new index by shifting the current index by origin inside the grid of given dimensions.
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[2], const len_type(&origin)[2], const len_type(&stride)[2])
	{
		return (n % dims[0] + origin[0]) * stride[0] + (n / dims[0] + origin[1]) * stride[1];
	}

	//! Compute the flattened position with an offset.
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[2], const len_type(&stride)[2], const iter_type(&offset)[2], iter_type delta)
	{
		return n * stride[0] - (n / dims[0]) * delta + offset[0] + offset[1];
	}


	//! Used for optimized calculation of the index from the flattened position index.
	/*!
	 * \param dims The dimensions of the region to iterate over.
	 * \param stride The stride of the global region.
	 */
	inline iter_type compute_delta(const len_type(&dims)[2], const len_type(&stride)[2])
	{
		return dims[0] * stride[0] - stride[1];
	}

	//! Compute the flattened position with an offset.
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[3], const len_type(&stride)[3], const iter_type(&offset)[3], iter_type delta)
	{
		auto s2 = dims[0] * dims[1];
		auto z = (n / s2) ;
		auto nxy = n - z * s2;
		auto mxy = nxy * stride[0] - (nxy / dims[0]) * delta + offset[0] + offset[1] + offset[2];
		return mxy + z * stride[2];
	}

	//! Used for optimized calculation of the index from the flattened position index.
	/*!
	 * \param dims The dimensions of the region to iterate over.
	 * \param stride The stride of the global region.
	 */
	inline iter_type compute_delta(const len_type(&dims)[3], const len_type(&stride)[3])
	{
		return dims[0] * stride[0] - stride[1];
	}

	//! Compute the flattened position with an offset.
	inline iter_type index_from_position(iter_type n, const len_type(&dims)[1], const len_type(&stride)[1], const iter_type(&offset)[1], iter_type delta)
	{
		return n * stride[0] + offset[0];
	}

	//! Used for optimized calculation of the index from the flattened position index.
	/*!
	 * \param dims The dimensions of the region to iterate over.
	 * \param stride The stride of the global region.
	 */
	inline iter_type compute_delta(const len_type(&dims)[1], const len_type(&stride)[1])
	{
		return 1;
	}
}


/*
 * the following defines dictate how to iterate over the one dimensional array owned
 * by a grid object to correctly visit each cell corresponding to the given boundary
 * ..it does not actually visit a grid object, only provides the loop structure for
 * getting there
 *
 * the defines are called with the desired operation, and two variables are provided
 * required is the dimensions of the grid
 *
 * there are separate iteration routines for edges/faces and corners
 * 
 * the indexing starts with TOP (so that the top left corner is the very first index)
 */

 /*
  * ENTRY keeps track of the actual loop iteration (i.e. always increments by 1)
  * INDEX is the grid array index corresponding to the ENTRY-th value on the chosen side
  */

  /* THREE dimensional iterations
   */

// faces, for when the whole system is a periodic system

//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ2 + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT(OP, __L, __M) \
for (int iter_k = BOUNDARY_DEPTH - 1, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + BOUNDARY_DEPTH, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + __L * __M * (__N - BOUNDARY_DEPTH) + BOUNDARY_DEPTH, ENTRY = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }


// full faces, for when there is a periodic side only on the opposite side

//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = 0, ENTRY = 0; iter_i < __N; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ2, ENTRY = 0; iter_i < __N; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_ALL(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = 0, ENTRY = 0; iter_j < __N; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2, ENTRY = 0; iter_j < __N; ++iter_j, INDEX += __SZ2) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_ALL(OP, __L, __M) \
for (int iter_k = BOUNDARY_DEPTH - 1, INDEX = 0, ENTRY = 0; iter_k >= 0; --iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_ALL(OP, __L, __M, __N) \
for (int iter_k = 0, INDEX = __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_k < BOUNDARY_DEPTH; ++iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }



//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_3A(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_3A(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_3A(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_3A(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ2 + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_3A(OP, __L, __M) \
for (int iter_k = BOUNDARY_DEPTH - 1, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + BOUNDARY_DEPTH, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_3A(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + __L * __M * (__N - BOUNDARY_DEPTH) + BOUNDARY_DEPTH, ENTRY = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }


//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_3AA(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_3AA(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_3AA(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_3AA(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ2 + BOUNDARY_DEPTH, ENTRY = 0; iter_j < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = BOUNDARY_DEPTH - 1; iter_k >= 0; --iter_k, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_3AA(OP, __L, __M) \
for (int iter_k = BOUNDARY_DEPTH - 1, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + BOUNDARY_DEPTH, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_3AA(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * BOUNDARY_DEPTH, INDEX = __SZ1 + __L * __M * (__N - BOUNDARY_DEPTH) + BOUNDARY_DEPTH, ENTRY = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }





// edges all periodic

//! Iterates over the edge at the top and front faces of a 3d grid.
#define ITER_GRID3_FRONT_TOP(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = BOUNDARY_DEPTH, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_k = 0; iter_k < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and front faces of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2 + BOUNDARY_DEPTH, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_k = 0; iter_k < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and top faces of a 3d grid.
#define ITER_GRID3_LEFT_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and bottom faces of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ2, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the back and top faces of a 3d grid.
#define ITER_GRID3_BACK_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * (__N - BOUNDARY_DEPTH) + BOUNDARY_DEPTH, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_k = 0; iter_k < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and back faces of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * (__N - BOUNDARY_DEPTH) + __SZ2 + BOUNDARY_DEPTH, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_k = 0; iter_k < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and top faces of a 3d grid.
#define ITER_GRID3_RIGHT_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and right faces of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * BOUNDARY_DEPTH + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and left faces of a 3d grid.
#define ITER_GRID3_FRONT_LEFT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ1, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and right faces of a 3d grid.
#define ITER_GRID3_FRONT_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and back faces of a 3d grid.
#define ITER_GRID3_BACK_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ1 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and back faces of a 3d grid.
#define ITER_GRID3_BACK_RIGHT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * BOUNDARY_DEPTH, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ1 + __SZ2 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }



// edges, full edges for the 3A and 3AA type updates



//! Iterates over the edge at the top and front faces of a 3d grid.
#define ITER_GRID3_FRONT_TOP_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = 0, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and front faces of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and top faces of a 3d grid.
#define ITER_GRID3_LEFT_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = 0, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and bottom faces of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the back and top faces of a 3d grid.
#define ITER_GRID3_BACK_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and back faces of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * (__N - BOUNDARY_DEPTH) + __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and top faces of a 3d grid.
#define ITER_GRID3_RIGHT_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and right faces of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and left faces of a 3d grid.
#define ITER_GRID3_FRONT_LEFT_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = 0, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and right faces of a 3d grid.
#define ITER_GRID3_FRONT_RIGHT_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and back faces of a 3d grid.
#define ITER_GRID3_BACK_LEFT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and back faces of a 3d grid.
#define ITER_GRID3_BACK_RIGHT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L - BOUNDARY_DEPTH, INDEX = __SZ2 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }



// corners all periodic

//! Iterates over the corner at the front top right intersection of a 3d grid.
#define ITER_GRID3_FRONT_TOP_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the front bottom right of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the right top back intersection of a 3d grid.
#define ITER_GRID3_RIGHT_TOP_BACK(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the right bottom back intersection of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM_BACK(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ1 + __SZ2 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the back top left intersection of a 3d grid.
#define ITER_GRID3_BACK_TOP_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the back bottom left intersection of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2 + __L * __M * (__N - BOUNDARY_DEPTH), ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the left top front intersection of a 3d grid.
#define ITER_GRID3_LEFT_TOP_FRONT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = 0, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the left bottom front intersection of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM_FRONT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - BOUNDARY_DEPTH, __SZ2 = __L * (__M - BOUNDARY_DEPTH), INDEX = __SZ2, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < BOUNDARY_DEPTH; ++iter_k, ++INDEX, ++ENTRY) { OP; }



/* TWO dimensional iterations
 */

// edges

//! Iterates over the left edge of a 2d grid.
#define ITER_GRID2_LEFT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * BOUNDARY_DEPTH, ENTRY = 0; iter_i < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = BOUNDARY_DEPTH - 1; iter_j >= 0; --iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right edge of a 2d grid.
#define ITER_GRID2_RIGHT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * BOUNDARY_DEPTH + __L - BOUNDARY_DEPTH, ENTRY = 0; iter_i < __M - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top edge of a 2d grid.
#define ITER_GRID2_TOP(OP, __L) \
for (int iter_j = 0, INDEX = BOUNDARY_DEPTH, ENTRY = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom edge of a 2d grid.
#define ITER_GRID2_BOTTOM(OP, __L, __M) \
for (int iter_j = BOUNDARY_DEPTH - 1, INDEX = __L * (__M - BOUNDARY_DEPTH) + BOUNDARY_DEPTH, ENTRY = 0; iter_j >= 0; --iter_j, INDEX += BOUNDARY_DEPTH + BOUNDARY_DEPTH) \
for (int iter_i = 0; iter_i < __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

// full edge, for mirroring only the direct opposite side

//! Iterates over the left edge of a 2d grid.
#define ITER_GRID2_LEFT_ALL(OP, __L, __M) \
for (int iter_i = 0, INDEX = 0, ENTRY = 0; iter_i < __M; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = BOUNDARY_DEPTH - 1; iter_j >= 0; --iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right edge of a 2d grid.
#define ITER_GRID2_RIGHT_ALL(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L - BOUNDARY_DEPTH, ENTRY = 0; iter_i < __M; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top edge of a 2d grid.
#define ITER_GRID2_TOP_ALL(OP, __L) \
for (int iter_j = 0, INDEX = 0, ENTRY = 0; iter_j < BOUNDARY_DEPTH; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom edge of a 2d grid.
#define ITER_GRID2_BOTTOM_ALL(OP, __L, __M) \
for (int iter_j = BOUNDARY_DEPTH - 1, INDEX = __L * (__M - BOUNDARY_DEPTH), ENTRY = 0; iter_j >= 0; --iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }



// corners across periodic

//! Iterates over the top left corner of a 2d grid.
#define ITER_GRID2_LEFT_TOP(OP, __L) \
for (int iter_i = 0, INDEX = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX) { OP; }

//! Iterates over the bottom right corner of a 2d grid.
#define ITER_GRID2_RIGHT_BOTTOM(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * (__M - BOUNDARY_DEPTH) + __L - BOUNDARY_DEPTH; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX) { OP; }

//! Iterates over the top right corner of a 2d grid.
#define ITER_GRID2_TOP_RIGHT(OP, __L) \
for (int iter_i = 0, INDEX = __L - BOUNDARY_DEPTH; iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX) { OP; }

//! Iterates over the bottom left corner of a 2d grid.
#define ITER_GRID2_BOTTOM_LEFT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * (__M - BOUNDARY_DEPTH); iter_i < BOUNDARY_DEPTH; ++iter_i, INDEX += __L - BOUNDARY_DEPTH) \
for (int iter_j = 0; iter_j < BOUNDARY_DEPTH; ++iter_j, ++INDEX) { OP; }


 /* ONE dimensional iterations
  */

  //! Iterates over the left side of a 1d grid.
#define ITER_GRID1_LEFT(OP) \
for (int iter_i = 0, INDEX = 0, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

  //! Iterates over the right side of a 1d grid.
#define ITER_GRID1_RIGHT(OP, __L) \
for (int iter_i = 0, INDEX = __L - BOUNDARY_DEPTH - BOUNDARY_DEPTH, ENTRY = 0; iter_i < BOUNDARY_DEPTH; ++iter_i, ++INDEX, ++ENTRY) { OP; }

 // **************************************************************************************



/*
 * executes the given function by passing it the index of each cell inside the grid
 * iterates from the start index to the end index in the grid
 * supports passing arguments to the function
 */

#define PARALLEL_ITER_GRID _Pragma("omp parallel for")


//! Iterates over the interior values of a 3d grid.
#define ITER_GRID3(OP, __L, __M, __N) { \
PARALLEL_ITER_GRID \
for (int iter_k = BOUNDARY_DEPTH; iter_k < __N - BOUNDARY_DEPTH; ++iter_k) { \
for (int iter_j = BOUNDARY_DEPTH; iter_j < __M - BOUNDARY_DEPTH; ++iter_j) { \
for (int iter_i = BOUNDARY_DEPTH; iter_i < __L - BOUNDARY_DEPTH; ++iter_i) { \
int INDEX = iter_i + iter_j * __L + iter_k * __L * __M; \
int ENTRY = (iter_i - BOUNDARY_DEPTH) + (iter_j - BOUNDARY_DEPTH) * (__L - BOUNDARY_DEPTH * 2) + (iter_k - BOUNDARY_DEPTH) * (__L - BOUNDARY_DEPTH * 2) * (__M - BOUNDARY_DEPTH * 2); \
{ OP; } } } } }

//! Iterates over the interior values of a 2d grid.
#define ITER_GRID2(OP, __L, __M) { \
PARALLEL_ITER_GRID \
for (int iter_j = BOUNDARY_DEPTH; iter_j < __M - BOUNDARY_DEPTH; ++iter_j) { \
for (int iter_i = BOUNDARY_DEPTH; iter_i < __L - BOUNDARY_DEPTH; ++iter_i) { \
int INDEX = iter_i + iter_j * __L; \
int ENTRY = (iter_i - BOUNDARY_DEPTH) + (iter_j - BOUNDARY_DEPTH) * (__L - BOUNDARY_DEPTH * 2); \
{ OP; } } } }

//! Iterates over the interior values of a 1d grid.
#define ITER_GRID1(OP, __L) { \
PARALLEL_ITER_GRID \
for (int iter_i = BOUNDARY_DEPTH; iter_i < __L - BOUNDARY_DEPTH; ++iter_i) { \
int INDEX = iter_i; int ENTRY = (iter_i - BOUNDARY_DEPTH); \
{ OP; } } }


//! Iterates over all the values of a 3d grid.
#define ITER_GRID3_ENTIRE(OP, __L, __M, __N) { \
PARALLEL_ITER_GRID \
for (int iter_k = 0; iter_k < __N; ++iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i) { \
int INDEX = iter_i + iter_j * __L + iter_k * __L * __M; \
{ OP; } } }

//! Iterates over all the values of a 2d grid.
#define ITER_GRID2_ENTIRE(OP, __L, __M) { \
PARALLEL_ITER_GRID \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i)  { \
int INDEX = iter_i + iter_j * __L; \
{ OP; } } }

//! Iterates over all the values of a 1d grid.
#define ITER_GRID1_ENTIRE(OP, __L) { \
PARALLEL_ITER_GRID \
for (int iter_i = 0; iter_i < __L; ++iter_i) { \
int INDEX = iter_i; \
 { OP; } } }





