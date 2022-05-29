
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
 * PURPOSE: Defines the grid types used in the finite difference and other
 * semi-discrete numerical solvers.
 *
 * ***************************************************************************
 */

#pragma once

#include <iostream>
#include <cassert>

#include "griditer.h"
#include "gridinfo.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif


//! \cond
#define INNERS_ARR_START 1
#define INNERS_ARR_INCREMENT 1
//! \endcond


/*!
 * \defgroup grid Grid Objects
 * @{
 */


namespace grid
{
	//! Compute the number of interior points.
	/*!
	 * For grids that include boundaries, only the interior values are typically
	 * considered relevant data. Returns the number of points that are interior.
	 * 
	 * \param dimensions The dimensions of the grid.
	 * 
	 * \tparam D The dimension of the grid.
	 */
	template<size_t D>
	len_type length_interior(len_type const* dimensions);

	//! Specialization of grid::length_interior(len_type const*).
	template<>
	inline len_type length_interior<1>(len_type const* dimensions)
	{
		return (dimensions != nullptr) ? std::max(0, dimensions[0] - THICKNESS - THICKNESS) : 0;
	}

	template<size_t D>
	len_type length_interior(len_type const* dimensions)
	{
		return ((dimensions != nullptr) ? std::max(0, (dimensions[D - 1] - THICKNESS - THICKNESS)) : 0) * length_interior<D - 1>(dimensions);
	}


	//! Return the total number of points in a grid.
	/*!
	 * Using the interval information about a grid, return the total number
	 * of points.
	 * 
	 * \param intervals The intervals that the grid spans.
	 * 
	 * \tparam D The dimension of the grid.
	 */
	template<size_t D>
	len_type length(symphas::interval_data_type intervals)
	{
		len_type dimensions[D];
		for (iter_type i = 0; i < D; ++i)
		{
			dimensions[i] = intervals.at(symphas::index_to_axis(i)).count();
		}

		return length<D>(dimensions);
	}

}




namespace grid
{

	namespace
	{


		/*
		 * helper functions to create an array for iterating over the interior values of a grid
		 *
		 * works by creating function with variables of static duration; calling the get function will
		 * create the array of interior index values the first time, then each subsequent call is made
		 * to immediately return the inner indices that was last used (based on the product of the
		 * dimensions)
		 *
		 * hence, each grid creation will immediately have access to the same list of inner values, which
		 * increase runtime and more importantly, decreases memory cost
		 */

		template<size_t D>
		using dim_arr_type = len_type[D];


		template<size_t I>
		void fill_inner_arr(iter_type* (&inner_i), len_type const* dimensions, iter_type& index);

		template<>
		inline void fill_inner_arr<0>(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			for (iter_type i = 0; i < dimensions[0] - THICKNESS - THICKNESS; ++i)
			{
				*(inner_i++) = index++;
			}
		}

		template<>
		inline void fill_inner_arr<1>(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			iter_type inc = THICKNESS + THICKNESS;
			for (iter_type i = 0; i < dimensions[1] - THICKNESS - THICKNESS; ++i, index += inc)
			{
				fill_inner_arr<0>(inner_i, dimensions, index);
			}
		}
		
		template<size_t I>
		void fill_inner_arr(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			iter_type inc = (THICKNESS + THICKNESS) * (dimensions[I - 2]);
			for (iter_type i = 0; i < dimensions[I] - THICKNESS - THICKNESS; ++i, index += inc)
			{
				fill_inner_arr<I - 1>(inner_i, dimensions, index);
			}
		}



		template<size_t D>
		iter_type get_init_index(len_type const* dimensions);

		template<>
		inline iter_type get_init_index<1>(len_type const*)
		{
			return THICKNESS;
		}

		template<size_t D>
		iter_type get_init_index(len_type const* dimensions)
		{
			return grid::length<D - 1>(dimensions) * THICKNESS + get_init_index<D - 1>(dimensions);
		}



		template<size_t D>
		iter_type* get_inner_arr(len_type const* dimensions)
		{
			iter_type* inner_i = new iter_type[grid::length_interior<D>(dimensions)];
			iter_type* iter = inner_i;
			iter_type index = get_init_index<D>(dimensions);
			fill_inner_arr<D - 1>(iter, dimensions, index);
			return inner_i;
		}

		/* append a new element to the array at the given index
		 */

		template<size_t D>
		void insert_dimensions_arr(dim_arr_type<D>* (&dimensions_arr), len_type index, len_type const* dimensions)
		{
			for (iter_type i = 0; i < D; ++i)
			{
				dimensions_arr[index][i] = dimensions[i];
			}
		}

		template<size_t D>
		void insert_inner_arr(iter_type** (&inner_arr), len_type index, len_type const* dimensions)
		{
			inner_arr[index] = get_inner_arr<D>(dimensions);
		}


		template<size_t D>
		dim_arr_type<D>* new_dimensions_arr(len_type const* dimensions)
		{
			dim_arr_type<D>* dimensions_arr = new dim_arr_type<D>[INNERS_ARR_START];
			insert_dimensions_arr<D>(dimensions_arr, 0, dimensions);
			return dimensions_arr;
		}

		template<size_t D>
		iter_type** new_inner_arr(len_type const* dimensions)
		{
			iter_type** inner_arr = new iter_type*[INNERS_ARR_START];
			insert_inner_arr<D>(inner_arr, 0, dimensions);
			return inner_arr;
		}


		template<size_t D>
		void extend_dimensions_arr(dim_arr_type<D>* (&dimensions_arr), len_type len)
		{
			dim_arr_type<D>* dimensions_arr_2 = new dim_arr_type<D>[len + INNERS_ARR_INCREMENT];
			for (iter_type n = 0; n < len; ++n)
			{
				insert_dimensions_arr(dimensions_arr_2, n, dimensions_arr[n]);
			}
			delete[] dimensions_arr;
			dimensions_arr = dimensions_arr_2;
		}

		template<size_t D>
		void extend_inner_arr(iter_type** (&inner_arr), len_type len)
		{
			iter_type** inner_arr_2 = new iter_type*[len + INNERS_ARR_INCREMENT];
			for (iter_type n = 0; n < len; ++n)
			{
				inner_arr_2[n] = inner_arr[n];
			}
			delete[] inner_arr;
			inner_arr = inner_arr_2;
		}

	}


	template<size_t D>
	iter_type* interior_indices_list(len_type const* dimensions)
	{
		/* on static variable creation, one list of indexes of inner values is created
		 */

		static len_type len = 1;							// the current count of static arrays
		static len_type capacity = INNERS_ARR_START;		// the total capacity of the arrays
		static iter_type index = 0;							// the last index accessed

		static dim_arr_type<D> *dimensions_arr = new_dimensions_arr<D>(dimensions);
		static iter_type** inner_arr = new_inner_arr<D>(dimensions);
		

		if (std::equal(dimensions, dimensions + D, dimensions_arr[index]))
		{
			return inner_arr[index];
		}
		else
		{
			for (iter_type i = (index + 1) % len; i != index; i = (i + 1) % len)
			{
				if (std::equal(dimensions, dimensions + D, dimensions_arr[index]))
				{
					index = i;
					return inner_arr[index];
				}
			}
		}

		if (len == capacity)
		{
			extend_dimensions_arr<D>(dimensions_arr, len);
			extend_inner_arr<D>(inner_arr, len);
			capacity += INNERS_ARR_INCREMENT;
		}

		len += 1;
		index = len - 1;
		insert_dimensions_arr<D>(dimensions_arr, index, dimensions);
		insert_inner_arr<D>(inner_arr, index, dimensions);

		return inner_arr[index];
	}
}


// ***********************************************************************************************

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 * 
 * \tparam T The value type of the underlying array.
 */
template<typename T>
struct Block
{

	T* values;		//!< The list of values managed by this object.
	len_type len;	//!< The number of values in the list.


	//! Create this object with \p len values.
	/*!
	 * Create this object with \p len values. The length can be 0, in which case
	 * no memory will be allocated, no values can be accessed, but the
	 * object can still be constructed.
	 * 
	 * \param len The number of values to create.
	 */
	Block(len_type len) : values{ (len > 0) ? new T[len] : nullptr }, len { len }
	{
		std::fill(values, values + len, T{});
	}

	explicit Block(const len_type* len, size_t dimensions = 1) : Block((len != nullptr) ? grid::length(len, dimensions) : 0) {}

	Block(Block<T> const& other) : Block(other.len)
	{
		std::copy(other.values, other.values + other.len, values);
	}

	Block(Block<T>&& other) noexcept : Block()
	{
		swap(*this, other);
	}

	Block& operator=(Block<T> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(Block<T> &first, Block<T> &second)
	{
		using std::swap;
		swap(first.len, second.len);
		swap(first.values, second.values);
	}

	//! Return the value at the \p i index in the list.
	/*!
	 * Return the value from the list from index \p i.
	 * 
	 * \param i The index of the value to return.
	 */
	T& operator[](iter_type i) const
	{
		return values[i];
	}

	~Block()
	{
		delete[] values;
	}

protected:

	Block() : values{ nullptr }, len{ 0 } {}


};

// ***********************************************************************************************

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type. The grid
 * forms the basis of a finite difference grid. Only manages its own 
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 * 
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template<typename T, size_t D>
struct Grid : Block<T>
{
	
protected:

	void set_dimensions(const len_type* dimensions)
	{
		if (dimensions == nullptr)
		{
			std::fill(dims, dims + D, 0);
		}
		else
		{
			for (iter_type i = 0; i < D; ++i)
			{
				dims[i] = dimensions[i];
			}
		}
	}

public:

	len_type dims[D];		//!< Dimensions of the grid, arranged beginning from the horizontal coordinate.

	//! Create a grid of the prescribed dimensions.
	/*!
	 * Creates a new grid using the given dimensions. The number of values in
	 * the grid directly correspond to the dimensions, equal to the product
	 * of all dimensions.
	 * 
	 * \param dimensions The dimensions of the grid.
	 */
	Grid(std::initializer_list<len_type> dimensions) : Block<T>{ grid::length<D>(dimensions.begin()) }, dims{ 0 }
	{
		set_dimensions(dimensions.begin());
	}

	//! Create a grid of the prescribed dimensions.
	/*!
	 * Creates a new grid using the given dimensions. The number of values in
	 * the grid directly correspond to the dimensions, equal to the product
	 * of all dimensions.
	 *
	 * \param dimensions The dimensions of the grid.
	 */
	Grid(const len_type* dimensions) : Block<T>{ grid::length<D>(dimensions) }, dims{ 0 }
	{
		set_dimensions(dimensions);
	}

	const Grid<T, D>& as_grid() const
	{
		return *this;
	}

	Grid<T, D>& as_grid()
	{
		return *this;
	}

protected:

	Grid() : Grid(nullptr) {}

};


// ***********************************************************************************************

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type. 
 * This grid implementation is an extension of the base ::Grid, but it
 * contains virtual boundary conditions.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template<typename T, size_t D>
struct BoundaryGrid : Grid<T, D>
{
	using Block<T>::values;
	using Grid<T, D>::dims;

#ifdef MULTITHREAD

	iter_type* inners;
	len_type len_inner;

public:
	BoundaryGrid(std::initializer_list<len_type> dimensions) : BoundaryGrid(dimensions.begin()) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, D>{ dimensions }, len_inner{ grid::length_interior<D>(dims) }, inners{ grid::interior_indices_list<D>(dims) } {}
	BoundaryGrid(Grid<T, D> const& other) : Grid<T, D>{ other }, len_inner{ grid::length_interior<D>(dims) }, inners{ grid::interior_indices_list<D>(dims) } {}
	BoundaryGrid(Grid<T, D>&& other) noexcept : Grid<T, D>{ other }, len_inner{ grid::length_interior<D>(dims) }, inners{ grid::interior_indices_list<D>(dims) } {}


	T& invalue(iter_type i)
	{
		return values[inners[i]];
	}

#else

	BoundaryGrid(std::initializer_list<len_type> dimensions) : Grid<T, D>{ dimensions } {}
	BoundaryGrid(Grid<T, D> const& other) : Grid<T, D>{ other } {}
	BoundaryGrid(Grid<T, D>&& other) noexcept : Grid<T, D>{ other } {}

	T& invalue(iter_type i)
	{
		iter_type index = 0;
		for (iter_type n = 0; n < D; ++n)
		{
			len_type sz = 1;
			len_type md = (n > 0) ? dims[D - n - 1] : 1;

			for (iter_type m = 0; m < D - n - 1; ++m)
			{
				sz *= dims[m];
			}
			
			iter_type ci = (i / sz) % md;
			index += ci + THICKNESS;
		}

		return values[index];
	}

#endif


	const BoundaryGrid<T, D>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, D>& as_grid()
	{
		return *this;
	}

protected:

#ifdef MULTITHREAD
	BoundaryGrid() : Grid<T, D>(), inners{ nullptr }, len_inner{ 0 } {}
#else
	BoundaryGrid() : Grid<T, D>() {}
#endif

};


template <typename T>
struct BoundaryGrid<T, 3> : Grid<T, 3>
{
	using Block<T>::values;
	using Grid<T, 3>::dims;

#ifdef MULTITHREAD

	iter_type* inners;
	len_type len_inner;

public:
	BoundaryGrid(std::initializer_list<len_type> dimensions) : BoundaryGrid(dimensions.begin()) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 3>(dimensions), len_inner{ grid::length_interior<3>(dims) }, inners{ grid::interior_indices_list<3>(dims) } {}
	BoundaryGrid(Grid<T, 3> const& other) : Grid<T, 3>(other), len_inner{ grid::length_interior<3>(dims) }, inners{ grid::interior_indices_list<3>(dims) } {}
	BoundaryGrid(Grid<T, 3>&& other) noexcept : Grid<T, 3>(other), len_inner{ grid::length_interior<3>(dims) }, inners{ grid::interior_indices_list<3>(dims) } {}


	T& invalue(iter_type i)
	{
		return values[inners[i]];
	}

#else

	BoundaryGrid(std::initializer_list<len_type> dimensions) : Grid<T, 3>(dimensions) {}
	BoundaryGrid(Grid<T, 3> const& other) : Grid<T, 3>(other) {}
	BoundaryGrid(Grid<T, 3>&& other) noexcept : Grid<T, 3>(other) {}


	T& invalue(iter_type i)
	{
		iter_type xi = i % dims[0];
		iter_type yi = (i / dims[0]) % dims[1];
		iter_type zi = i / (dims[0] * dims[1]);
		return values[(zi + THICKNESS) * dims[0] * dims[1] + (yi + THICKNESS) * dims[0] + (xi + THICKNESS)];
	}

#endif




	const BoundaryGrid<T, 3>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 3>& as_grid()
	{
		return *this;
	}

	void copy_ptr_left(T** into);
	void copy_ptr_right(T** into);
	void copy_ptr_top(T** into);
	void copy_ptr_bottom(T** into);
	void copy_ptr_front(T** into);
	void copy_ptr_back(T** into);


protected:

#ifdef MULTITHREAD
	BoundaryGrid() : Grid<T, 3>(), inners{ nullptr }, len_inner{ 0 } {}
#else
	BoundaryGrid() : Grid<T, 3>() {}
#endif

};


template <typename T>
struct BoundaryGrid<T, 2> : Grid<T, 2>
{
	using Block<T>::values;
	using Grid<T, 2>::dims;

#ifdef MULTITHREAD

	iter_type* inners;
	len_type len_inner;

public:
	BoundaryGrid(std::initializer_list<len_type> dimensions) : BoundaryGrid(dimensions.begin()) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 2>(dimensions), len_inner{ grid::length_interior<2>(dims) }, inners{ grid::interior_indices_list<2>(dims) } {}
	BoundaryGrid(Grid<T, 2> const& other) : Grid<T, 2>(other), len_inner{ grid::length_interior<2>(dims) }, inners{ grid::interior_indices_list<2>(dims) } {}
	BoundaryGrid(Grid<T, 2>&& other) noexcept : Grid<T, 2>(other), len_inner{ grid::length_interior<2>(dims) }, inners{ grid::interior_indices_list<2>(dims) } {}

	T& invalue(iter_type i)
	{
		return values[inners[i]];
	}

#else

	BoundaryGrid(std::initializer_list<len_type> dimensions) : Grid<T, 2>(dimensions) {}
	BoundaryGrid(Grid<T, 2> const& other) : Grid<T, 2>(other) {}
	BoundaryGrid(Grid<T, 2>&& other) noexcept : Grid<T, 2>(other) {}

	T& invalue(iter_type i)
	{
		iter_type xi = i % dims[0];
		iter_type yi = i / dims[0];
		return values[(yi + THICKNESS) * dims[0] + (xi + THICKNESS)];
	}

#endif

	void copy_ptr_left(T** into);
	void copy_ptr_right(T** into);
	void copy_ptr_top(T** into);
	void copy_ptr_bottom(T** into);


	const BoundaryGrid<T, 2>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 2>& as_grid()
	{
		return *this;
	}

protected:

#ifdef MULTITHREAD
	BoundaryGrid() : Grid<T, 2>(), inners{ nullptr }, len_inner{ 0 } {}
#else
	BoundaryGrid() : Grid<T, 2>() {}
#endif

};


template <typename T>
struct BoundaryGrid<T, 1> : Grid<T, 1>
{
	using Block<T>::values;
	using Grid<T, 1>::dims;

#ifdef MULTITHREAD

	iter_type* inners;
	len_type len_inner;

public:
	BoundaryGrid(std::initializer_list<len_type> dimensions) : BoundaryGrid(dimensions.begin()) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 1>(dimensions), len_inner{ grid::length_interior<1>(dims) }, inners{ grid::interior_indices_list<1>(dims) } {}
	BoundaryGrid(Grid<T, 1> const& other) : Grid<T, 1>(other), len_inner{ grid::length_interior<1>(dims) }, inners{ grid::interior_indices_list<1>(dims) } {}
	BoundaryGrid(Grid<T, 1>&& other) noexcept : Grid<T, 1>(other), len_inner{ grid::length_interior<1>(dims) }, inners{ grid::interior_indices_list<1>(dims) } {}

	T& invalue(iter_type i)
	{
		return values[inners[i]];
	}

#else

	BoundaryGrid(std::initializer_list<len_type> dimensions) : Grid<T, 1>(dimensions) {}
	BoundaryGrid(Grid<T, 1> const& other) : Grid<T, 1>(other) {}
	BoundaryGrid(Grid<T, 1>&& other) noexcept : Grid<T, 1>(other) {}

	T& invalue(iter_type i)
	{
		return values[i + THICKNESS];
	}

#endif


	const BoundaryGrid<T, 1>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 1>& as_grid()
	{
		return *this;
	}

	void copy_ptr_left(T** into);
	void copy_ptr_right(T** into);


protected:


#ifdef MULTITHREAD
	BoundaryGrid() : Grid<T, 1>(), inners{ nullptr }, len_inner{ 0 } {}
#else
	BoundaryGrid() : Grid<T, 1>() {}
#endif

};




// **************************************************************************************



template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_left(T** into)
{
	ITER_GRID3_LEFT(into[ENTRY] = values + INDEX, dims[0], dims[1], dims[2]);
}

template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_right(T** into)
{
	ITER_GRID3_RIGHT(into[ENTRY] = values + INDEX, dims[0], dims[1], dims[2]);
}

template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_top(T** into)
{
	ITER_GRID3_TOP(into[ENTRY] = values + INDEX, dims[0], dims[1], dims[2]);
}

template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_bottom(T** into)
{
	ITER_GRID3_BOTTOM(into[ENTRY] = values + INDEX, dims[0], dims[1], dims[2]);
}

template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_front(T** into)
{
	ITER_GRID3_FRONT(into[ENTRY] = values + INDEX, dims[0], dims[1]);
}

template<typename T>
void BoundaryGrid<T, 3>::copy_ptr_back(T** into)
{
	ITER_GRID3_BACK(into[ENTRY] = values + INDEX, dims[0], dims[1], dims[2]);
}

template<typename T>
void BoundaryGrid<T, 2>::copy_ptr_left(T** into)
{

	ITER_GRID2_LEFT(into[ENTRY] = values + INDEX, dims[0], dims[1]);
}

template<typename T>
void BoundaryGrid<T, 2>::copy_ptr_right(T** into)
{
	ITER_GRID2_RIGHT(into[ENTRY] = values + INDEX, dims[0], dims[1]);
}

template<typename T>
void BoundaryGrid<T, 2>::copy_ptr_top(T** into)
{
	ITER_GRID2_TOP(into[ENTRY] = values + INDEX, dims[0]);
}

template<typename T>
void BoundaryGrid<T, 2>::copy_ptr_bottom(T** into)
{
	ITER_GRID2_BOTTOM(into[ENTRY] = values + INDEX, dims[0], dims[1]);
}

template<typename T>
void BoundaryGrid<T, 1>::copy_ptr_left(T** into)
{
	ITER_GRID1_LEFT(into[ENTRY] = values + INDEX);
}

template<typename T>
void BoundaryGrid<T, 1>::copy_ptr_right(T** into)
{
	ITER_GRID1_RIGHT(into[ENTRY] = values + INDEX, dims[0]);
}




// **************************************************************************************

namespace grid
{
	//! Copy the interior values of the grid into an array.
	/*!
	 * The interior values of the given grid are copied into an array.
	 * It is assumed that the array has enough space to store all the interior
	 * values. For computing the number of interior points, see 
	 * grid::length_interior(len_type const*). The grid is 1-dimensional.
	 * 
	 * \param g The grid from which the interior values are copied.
	 * \param t The array into which the values are copied.
	 */
	template<typename T>
	void copy_interior(Grid<T, 1> const& g, T* t)
	{
		ITER_GRID1(t[ENTRY] = g.values[INDEX], g.dims[0])
	}


	//! Copy the interior values of the grid into an array.
	/*!
	 * Implementation of copying interior values for a 2-dimensional Grid, see
	 * grid::copy_interior(Grid<T, 1> const&, T*).
	 *
	 * \param g The grid from which the interior values are copied.
	 * \param t The array into which the values are copied.
	 */
	template<typename T>
	void copy_interior(Grid<T, 2> const& g, T* t)
	{
		ITER_GRID2(t[ENTRY] = g.values[INDEX], g.dims[0], g.dims[1]);
	}

	//! Copy the interior values of the grid into an array.
	/*!
	 * Implementation of copying interior values for a 3-dimensional Grid, see
	 * grid::copy_interior(Grid<T, 1> const&, T*).
	 *
	 * \param g The grid from which the interior values are copied.
	 * \param t The array into which the values are copied.
	 */
	template<typename T>
	void copy_interior(Grid<T, 3> const& g, T* t)
	{
		ITER_GRID3(t[ENTRY] = g.values[INDEX], g.dims[0], g.dims[1], g.dims[2]);
	}

	namespace
	{

		//! Create a new grid of the given primary type.
		/*!
		 * This is used to generate a new grid of the prescribed dimensions
		 * and of the prescribed source type, for the given primary type. The
		 * primary type is not a complete type. It defines the class of
		 * this structure which generates the grid.
		 * 
		 * \tparam G The primary grid type, templated on the value type and 
		 * the dimension, is provided.
		 */
		template<template<typename, size_t> typename G>
		struct make
		{
			//! Obtain a newly initialized grid using the interval data.
			/*
			 * This is used to generate a new grid of the prescribed dimensions
			 * and of the prescribed source type. The dimensions are assumed
			 * to all be zero, thereby creating an empty grid.
			 */
			template<typename T, size_t D>
			static G<T, D> apply()
			{
				len_type dims[D]{ 0 };
				return { dims };
			}

			//! Obtain a newly initialized grid using the interval data.
			/*
			 * This is used to generate a new grid of the prescribed dimensions
			 * and of the prescribed source type. The given interval data is
			 * used to obtain the dimensions.
			 *
			 * \param vdata The interval data of the grid to generate.
			 */
			template<typename T, size_t D>
			static G<T, D> apply(symphas::interval_data_type const& vdata)
			{
				if (vdata.size() < D)
				{
					return apply<T, D>();
				}

				len_type dims[D]{ 0 };
				for (iter_type i = 0; i < D; ++i)
				{
					Axis side = symphas::index_to_axis(i);
					dims[i] = vdata.at(side).count();
				}
				return { dims };
			}
		};
	}

	//! Create a new grid with of the given base type and dimension.
	/*!
	 * This is used to generate a new grid of the prescribed dimensions
	 * and of the prescribed source type.
	 * 
	 * \param vdata The intervals of the grid, from which the dimensions
	 * are taken.
	 */
	template<template<typename, size_t> typename G, typename T, size_t D>
	G<T, D> construct(symphas::interval_data_type const& vdata)
	{
		return make<G>::template apply<T, D>(vdata);
	}

	//! Create a new grid with of the given base type and dimension.
	/*!
	 * This is used to generate a new grid of the prescribed dimensions
	 * and of the prescribed source type.
	 *
	 * \param vdata The intervals of the grid, from which the dimensions
	 * are taken.
	 */
	template<template<typename, size_t> typename G, typename T, size_t D>
	G<T, D> construct()
	{
		return make<G>::template apply<T, D>();
	}


	//! From the given parameter of base type Block, infer the underlying type.
	/*!
	 * Type traits technique is used to infer the underlying type of the given
	 * template parameter, which must be a child type of Block. If the given type
	 * cannot be implicitly cast to Block, then a compile error will be generated.
	 *
	 * \tparam G The grid type to infer the underlying type from.
	 */
	template<typename G>
	struct value_type_of
	{
	protected:

		template<typename T>
		struct wrap_type
		{
			using type = T;
		};

		template<typename T>
		static constexpr wrap_type<T> cast(Block<T>)
		{
			return {};
		}


		static constexpr auto call_wrap(G g)
		{
			return cast(g);
		}

		using wrapped_type = typename std::invoke_result_t<decltype(&grid::value_type_of<G>::call_wrap), G>;

	public:

		using type = typename wrapped_type::type;

	};

	//! Specialization of grid::value_type_of.
	/*!
	 * Specialization based on the reference of the given type.
	 */
	template<typename G>
	struct value_type_of<G&>
	{
		using type = typename grid::value_type_of<G>::type;
	};




	//! From the given parameter of base type Grid, infer the underlying dimension.
	/*!
	 * Type traits technique is used to infer the underlying dimension of the given
	 * template parameter, which must be a child type of Grid. If the given type
	 * cannot be implicitly cast to Grid, then a compile error will be generated.
	 *
	 * \tparam G The grid type to infer the underlying dimension from.
	 */
	template<typename G>
	struct dimension_of
	{
	protected:

		template<size_t D>
		struct wrap_type
		{
			static const size_t value = D;
		};

		template<typename T, size_t D>
		static constexpr wrap_type<D> cast(Grid<T, D>)
		{
			return {};
		}


		static constexpr auto call_wrap(G g)
		{
			return cast(g);
		}

		using wrapped_type = typename std::invoke_result_t<decltype(&grid::dimension_of<G>::call_wrap), G>;

	public:

		static const size_t value = wrapped_type::value;

	};

}


//! @}





