
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
for (int iter_i = 0, __SZ1 = __L * THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __L - THICKNESS) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ2 + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT(OP, __L, __M) \
for (int iter_k = THICKNESS - 1, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + THICKNESS, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + __L * __M * (__N - THICKNESS) + THICKNESS, ENTRY = 0; iter_k < THICKNESS; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }


// full faces, for when there is a periodic side only on the opposite side

//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, INDEX = 0, ENTRY = 0; iter_i < __N; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __L - THICKNESS) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __SZ2, ENTRY = 0; iter_i < __N; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_ALL(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = 0, ENTRY = 0; iter_j < __N; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2, ENTRY = 0; iter_j < __N; ++iter_j, INDEX += __SZ2) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_ALL(OP, __L, __M) \
for (int iter_k = THICKNESS - 1, INDEX = 0, ENTRY = 0; iter_k >= 0; --iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_ALL(OP, __L, __M, __N) \
for (int iter_k = 0, INDEX = __L * __M * (__N - THICKNESS), ENTRY = 0; iter_k < THICKNESS; ++iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }



//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_3A(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __L - THICKNESS) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_3A(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_3A(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_3A(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ2 + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_3A(OP, __L, __M) \
for (int iter_k = THICKNESS - 1, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + THICKNESS, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_3A(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + __L * __M * (__N - THICKNESS) + THICKNESS, ENTRY = 0; iter_k < THICKNESS; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }


//! Iterates over the left face of a 3d grid.
#define ITER_GRID3_LEFT_3AA(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __L - THICKNESS) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right face of a 3d grid.
#define ITER_GRID3_RIGHT_3AA(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __L * __M * THICKNESS + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top face of a 3d grid.
#define ITER_GRID3_TOP_3AA(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom face of a 3d grid.
#define ITER_GRID3_BOTTOM_3AA(OP, __L, __M, __N) \
for (int iter_j = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ2 + THICKNESS, ENTRY = 0; iter_j < __N - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = THICKNESS - 1; iter_k >= 0; --iter_k, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the front face of a 3d grid.
#define ITER_GRID3_FRONT_3AA(OP, __L, __M) \
for (int iter_k = THICKNESS - 1, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + THICKNESS, ENTRY = 0; iter_k >= 0; --iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the back face of a 3d grid.
#define ITER_GRID3_BACK_3AA(OP, __L, __M, __N) \
for (int iter_k = 0, __SZ1 = __L * THICKNESS, INDEX = __SZ1 + __L * __M * (__N - THICKNESS) + THICKNESS, ENTRY = 0; iter_k < THICKNESS; ++iter_k, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }





// edges all periodic

//! Iterates over the edge at the top and front faces of a 3d grid.
#define ITER_GRID3_FRONT_TOP(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = THICKNESS, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_k = 0; iter_k < __L - THICKNESS - THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and front faces of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2 + THICKNESS, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_k = 0; iter_k < __L - THICKNESS - THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and top faces of a 3d grid.
#define ITER_GRID3_LEFT_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and bottom faces of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ2, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the back and top faces of a 3d grid.
#define ITER_GRID3_BACK_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * (__N - THICKNESS) + THICKNESS, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_k = 0; iter_k < __L - THICKNESS - THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and back faces of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * (__N - THICKNESS) + __SZ2 + THICKNESS, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_k = 0; iter_k < __L - THICKNESS - THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and top faces of a 3d grid.
#define ITER_GRID3_RIGHT_TOP(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ1, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and right faces of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * THICKNESS + __SZ1 + __SZ2, ENTRY = 0; iter_i < __N - THICKNESS - THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and left faces of a 3d grid.
#define ITER_GRID3_FRONT_LEFT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __SZ1, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and right faces of a 3d grid.
#define ITER_GRID3_FRONT_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and back faces of a 3d grid.
#define ITER_GRID3_BACK_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __SZ1 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and back faces of a 3d grid.
#define ITER_GRID3_BACK_RIGHT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L * THICKNESS, __SZ2 = __L - THICKNESS, INDEX = __SZ1 + __SZ2 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ1 + __SZ1) \
for (int iter_j = 0; iter_j < __M - THICKNESS - THICKNESS; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }



// edges, full edges for the 3A and 3AA type updates



//! Iterates over the edge at the top and front faces of a 3d grid.
#define ITER_GRID3_FRONT_TOP_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = 0, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and front faces of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and top faces of a 3d grid.
#define ITER_GRID3_LEFT_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = 0, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and bottom faces of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the back and top faces of a 3d grid.
#define ITER_GRID3_BACK_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and back faces of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * (__N - THICKNESS) + __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j) \
for (int iter_k = 0; iter_k < __L; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and top faces of a 3d grid.
#define ITER_GRID3_RIGHT_TOP_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the bottom and right faces of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < __N; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and left faces of a 3d grid.
#define ITER_GRID3_FRONT_LEFT_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L - THICKNESS, INDEX = 0, ENTRY = 0; iter_i < THICKNESS; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the front and right faces of a 3d grid.
#define ITER_GRID3_FRONT_RIGHT_ALL(OP, __L, __M) \
for (int iter_i = 0, __SZ2 = __L - THICKNESS, INDEX = __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the left and back faces of a 3d grid.
#define ITER_GRID3_BACK_LEFT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L - THICKNESS, INDEX = __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the edge at the right and back faces of a 3d grid.
#define ITER_GRID3_BACK_RIGHT_ALL(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ2 = __L - THICKNESS, INDEX = __SZ2 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i) \
for (int iter_j = 0; iter_j < __M; ++iter_j, INDEX += __SZ2) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }



// corners all periodic

//! Iterates over the corner at the front top right intersection of a 3d grid.
#define ITER_GRID3_FRONT_TOP_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the front bottom right of a 3d grid.
#define ITER_GRID3_FRONT_BOTTOM_RIGHT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1 + __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the right top back intersection of a 3d grid.
#define ITER_GRID3_RIGHT_TOP_BACK(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the right bottom back intersection of a 3d grid.
#define ITER_GRID3_RIGHT_BOTTOM_BACK(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ1 + __SZ2 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the back top left intersection of a 3d grid.
#define ITER_GRID3_BACK_TOP_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the back bottom left intersection of a 3d grid.
#define ITER_GRID3_BACK_BOTTOM_LEFT(OP, __L, __M, __N) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2 + __L * __M * (__N - THICKNESS), ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the left top front intersection of a 3d grid.
#define ITER_GRID3_LEFT_TOP_FRONT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = 0, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the corner at the left bottom front intersection of a 3d grid.
#define ITER_GRID3_LEFT_BOTTOM_FRONT(OP, __L, __M) \
for (int iter_i = 0, __SZ1 = __L - THICKNESS, __SZ2 = __L * (__M - THICKNESS), INDEX = __SZ2, ENTRY = 0; iter_i < THICKNESS; ++iter_i, INDEX += __SZ2) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, INDEX += __SZ1) \
for (int iter_k = 0; iter_k < THICKNESS; ++iter_k, ++INDEX, ++ENTRY) { OP; }



/* TWO dimensional iterations
 */

// edges

//! Iterates over the left edge of a 2d grid.
#define ITER_GRID2_LEFT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * THICKNESS, ENTRY = 0; iter_i < __M - THICKNESS - THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = THICKNESS - 1; iter_j >= 0; --iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right edge of a 2d grid.
#define ITER_GRID2_RIGHT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * THICKNESS + __L - THICKNESS, ENTRY = 0; iter_i < __M - THICKNESS - THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top edge of a 2d grid.
#define ITER_GRID2_TOP(OP, __L) \
for (int iter_j = 0, INDEX = THICKNESS, ENTRY = 0; iter_j < THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom edge of a 2d grid.
#define ITER_GRID2_BOTTOM(OP, __L, __M) \
for (int iter_j = THICKNESS - 1, INDEX = __L * (__M - THICKNESS) + THICKNESS, ENTRY = 0; iter_j >= 0; --iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = 0; iter_i < __L - THICKNESS - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

// full edge, for mirroring only the direct opposite side

//! Iterates over the left edge of a 2d grid.
#define ITER_GRID2_LEFT_ALL(OP, __L, __M) \
for (int iter_i = 0, INDEX = 0, ENTRY = 0; iter_i < __M; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = THICKNESS - 1; iter_j >= 0; --iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the right edge of a 2d grid.
#define ITER_GRID2_RIGHT_ALL(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L - THICKNESS, ENTRY = 0; iter_i < __M; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the top edge of a 2d grid.
#define ITER_GRID2_TOP_ALL(OP, __L) \
for (int iter_j = 0, INDEX = 0, ENTRY = 0; iter_j < THICKNESS; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the bottom edge of a 2d grid.
#define ITER_GRID2_BOTTOM_ALL(OP, __L, __M) \
for (int iter_j = THICKNESS - 1, INDEX = __L * (__M - THICKNESS), ENTRY = 0; iter_j >= 0; --iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX, ++ENTRY) { OP; }



// corners across periodic

//! Iterates over the top left corner of a 2d grid.
#define ITER_GRID2_LEFT_TOP(OP, __L) \
for (int iter_i = 0, INDEX = 0; iter_i < THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX) { OP; }

//! Iterates over the bottom right corner of a 2d grid.
#define ITER_GRID2_RIGHT_BOTTOM(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * (__M - THICKNESS) + __L - THICKNESS; iter_i < THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX) { OP; }

//! Iterates over the top right corner of a 2d grid.
#define ITER_GRID2_TOP_RIGHT(OP, __L) \
for (int iter_i = 0, INDEX = __L - THICKNESS; iter_i < THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX) { OP; }

//! Iterates over the bottom left corner of a 2d grid.
#define ITER_GRID2_BOTTOM_LEFT(OP, __L, __M) \
for (int iter_i = 0, INDEX = __L * (__M - THICKNESS); iter_i < THICKNESS; ++iter_i, INDEX += __L - THICKNESS) \
for (int iter_j = 0; iter_j < THICKNESS; ++iter_j, ++INDEX) { OP; }


 /* ONE dimensional iterations
  */

  //! Iterates over the left side of a 1d grid.
#define ITER_GRID1_LEFT(OP) \
for (int iter_i = 0, INDEX = 0, ENTRY = 0; iter_i < THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

  //! Iterates over the right side of a 1d grid.
#define ITER_GRID1_RIGHT(OP, __L) \
for (int iter_i = 0, INDEX = __L - THICKNESS - THICKNESS, ENTRY = 0; iter_i < THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

 // **************************************************************************************



/*
 * executes the given function by passing it the index of each cell inside the grid
 * iterates from the start index to the end index in the grid
 * supports passing arguments to the function
 */


//! Iterates over the interior values of a 3d grid.
#define ITER_GRID3(OP, __L, __M, __N) \
for (int iter_k = THICKNESS, INDEX = THICKNESS + __L * THICKNESS + __L * __M * THICKNESS, ENTRY = 0, __SZ1 = (THICKNESS + THICKNESS) * __L; iter_k < __N - THICKNESS; ++iter_k, INDEX += __SZ1) \
for (int iter_j = THICKNESS; iter_j < __M - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = THICKNESS; iter_i < __L - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the interior values of a 2d grid.
#define ITER_GRID2(OP, __L, __M) \
for (int iter_j = THICKNESS, INDEX = THICKNESS + __L * THICKNESS, ENTRY = 0; iter_j < __M - THICKNESS; ++iter_j, INDEX += THICKNESS + THICKNESS) \
for (int iter_i = THICKNESS; iter_i < __L - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }

//! Iterates over the interior values of a 1d grid.
#define ITER_GRID1(OP, __L) \
for (int iter_i = THICKNESS, INDEX = THICKNESS, ENTRY = 0; iter_i < __L - THICKNESS; ++iter_i, ++INDEX, ++ENTRY) { OP; }


//! Iterates over all the values of a 3d grid.
#define ITER_GRID3_ENTIRE(OP, __L, __M, __N) \
for (int iter_k = 0, INDEX = 0; iter_k < __N; ++iter_k) \
for (int iter_j = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX) { OP; }

//! Iterates over all the values of a 2d grid.
#define ITER_GRID2_ENTIRE(OP, __L, __M) \
for (int iter_j = 0, INDEX = 0; iter_j < __M; ++iter_j) \
for (int iter_i = 0; iter_i < __L; ++iter_i, ++INDEX) { OP; }

//! Iterates over all the values of a 1d grid.
#define ITER_GRID1_ENTIRE(OP, __L) \
for (int iter_i = 0, INDEX = 0; iter_i < __L; ++iter_i, ++INDEX) { OP; }





