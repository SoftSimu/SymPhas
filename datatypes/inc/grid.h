
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
#include <random>
#include <omp.h>

#include "spsmpi.h"
#include "griditer.h"

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
	inline bool is_overlapping_compare(iter_type left0, iter_type left1, iter_type right0, iter_type right1)
	{
		return (right0 > left1 && right1 > left0);
		//return ((std::max(region0[Is][0], region1[Is][0]) <= std::min(region0[Is][1], region1[Is][1])) && ...);
	}

	template<size_t D, size_t... Is>
	bool is_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2], std::index_sequence<Is...>)
	{
		return (is_overlapping_compare(region0[Is][0], region1[Is][0], region0[Is][1], region1[Is][1]) && ...);
	}

	template<size_t D>
	bool is_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2])
	{
		return is_overlapping(region0, region1, std::make_index_sequence<D>{});
	}

	template<size_t D, size_t... Is>
	bool is_fully_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2], std::index_sequence<Is...>)
	{
		return ((region0[Is][0] <= region1[Is][0] && region0[Is][1] >= region1[Is][1]) && ...);
	}

	template<size_t D>
	bool is_fully_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2])
	{
		return is_fully_overlapping(region0, region1, std::make_index_sequence<D>{});
	}


	inline bool is_contact_overlapping_compare(iter_type left0, iter_type left1, iter_type right0, iter_type right1)
	{
		return (right0 >= left1 && right1 >= left0);
		//return ((std::max(region0[Is][0], region1[Is][0]) <= std::min(region0[Is][1], region1[Is][1])) && ...);
	}

	template<size_t D, size_t... Is>
	bool is_contact_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2], std::index_sequence<Is...>)
	{
		return (is_contact_overlapping_compare(region0[Is][0], region1[Is][0], region0[Is][1], region1[Is][1]) && ...);
	}

	template<size_t D>
	bool is_contact_overlapping(const iter_type(&region0)[D][2], const iter_type(&region1)[D][2])
	{
		return is_contact_overlapping(region0, region1, std::make_index_sequence<D>{});
	}

	template<size_t D>
	void get_region_union(iter_type(&result)[D][2], const iter_type(&region0)[D][2], const iter_type(&region1)[D][2])
	{
		for (iter_type i = 0; i < D; ++i)
		{
			result[i][0] = std::min(region0[i][0], region1[i][0]);
			result[i][1] = std::max(region0[i][1], region1[i][1]);
		}
	}

	template<size_t D>
	void get_region_intersection(iter_type(&result)[D][2], const iter_type(&region0)[D][2], const iter_type(&region1)[D][2])
	{
		for (iter_type i = 0; i < D; ++i)
		{
			result[i][0] = std::max(region0[i][0], region1[i][0]);
			result[i][1] = std::min(region0[i][1], region1[i][1]);
		}
	}



	// Explicit forms to allow implicit casting.

	inline bool is_overlapping(const iter_type(&region0)[1][2], const iter_type(&region1)[1][2])
	{
		return is_overlapping<1>(region0, region1);
	}
	inline bool is_overlapping(const iter_type(&region0)[2][2], const iter_type(&region1)[2][2])
	{
		return is_overlapping<2>(region0, region1);
	}
	inline bool is_overlapping(const iter_type(&region0)[3][2], const iter_type(&region1)[3][2])
	{
		return is_overlapping<3>(region0, region1);
	}

	inline bool is_fully_overlapping(const iter_type(&region0)[1][2], const iter_type(&region1)[1][2])
	{
		return is_fully_overlapping<1>(region0, region1);
	}
	inline bool is_fully_overlapping(const iter_type(&region0)[2][2], const iter_type(&region1)[2][2])
	{
		return is_fully_overlapping<2>(region0, region1);
	}
	inline bool is_fully_overlapping(const iter_type(&region0)[3][2], const iter_type(&region1)[3][2])
	{
		return is_fully_overlapping<3>(region0, region1);
	}

	inline bool is_contact_overlapping(const iter_type(&region0)[1][2], const iter_type(&region1)[1][2])
	{
		return is_contact_overlapping<1>(region0, region1);
	}
	inline bool is_contact_overlapping(const iter_type(&region0)[2][2], const iter_type(&region1)[2][2])
	{
		return is_contact_overlapping<2>(region0, region1);
	}
	inline bool is_contact_overlapping(const iter_type(&region0)[3][2], const iter_type(&region1)[3][2])
	{
		return is_contact_overlapping<3>(region0, region1);
	}

	inline void get_region_union(iter_type(&result)[1][2], const iter_type(&region0)[1][2], const iter_type(&region1)[1][2])
	{
		get_region_union<1>(result, region0, region1);
	}
	inline void get_region_union(iter_type(&result)[2][2], const iter_type(&region0)[2][2], const iter_type(&region1)[2][2])
	{
		get_region_union<2>(result, region0, region1);
	}
	inline void get_region_union(iter_type(&result)[3][2], const iter_type(&region0)[3][2], const iter_type(&region1)[3][2])
	{
		get_region_union<3>(result, region0, region1);
	}

	inline void get_region_intersection(iter_type(&result)[1][2], const iter_type(&region0)[1][2], const iter_type(&region1)[1][2])
	{
		get_region_intersection<1>(result, region0, region1);
	}
	inline void get_region_intersection(iter_type(&result)[2][2], const iter_type(&region0)[2][2], const iter_type(&region1)[2][2])
	{
		get_region_intersection<2>(result, region0, region1);
	}
	inline void get_region_intersection(iter_type(&result)[3][2], const iter_type(&region0)[3][2], const iter_type(&region1)[3][2])
	{
		get_region_intersection<3>(result, region0, region1);
	}



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
	len_type length_interior(len_type const* dimensions, len_type boundary_size = BOUNDARY_DEPTH);

	//! Specialization of grid::length_interior(len_type const*).
	template<>
	inline len_type length_interior<1>(len_type const* dimensions, len_type boundary_size)
	{
		return (dimensions != nullptr) ? std::max(0, dimensions[0] - boundary_size - boundary_size) : 0;
	}

	template<size_t D>
	len_type length_interior(len_type const* dimensions, len_type boundary_size)
	{
		return ((dimensions != nullptr) ? std::max(0, (dimensions[D - 1] - boundary_size - boundary_size)) : 0) * length_interior<D - 1>(dimensions, boundary_size);
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
			dimensions[i] = intervals.at(symphas::index_to_axis(i)).get_count();
		}

		return length<D>(dimensions);
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
	len_type length_interior(symphas::interval_data_type intervals)
	{
		len_type dimensions[D];
		for (iter_type i = 0; i < D; ++i)
		{
			dimensions[i] = intervals.at(symphas::index_to_axis(i)).get_interval_count();
		}

		return length<D>(dimensions);
	}

}


template<size_t N, typename T>
struct multi_value;


template<typename T>
struct carry_value
{
	carry_value() : value{ &fallback }, fallback{ 0 } {}
	carry_value(T fallback) : value{ &this->fallback }, fallback{ fallback } {}
	carry_value(T* value, T fallback) : value{ value }, fallback{ fallback } {}
	carry_value(T* value, T fallback, bool valid) : value{ this->value = (valid) ? value : &this->fallback }, fallback{ fallback } {}
	carry_value(carry_value<T> const& other) : carry_value{ other.value, other.fallback, other.is_valid() } {}
	carry_value(carry_value<T>&& other) : carry_value() 
	{
		swap(*this, other);
	}

	friend void swap(carry_value<T>& first, carry_value<T>& second)
	{
		using std::swap;
		T* tmp = second.value;
		if (first.is_valid())
		{
			second.value = first.value;
		}
		else
		{
			second.value = &second.fallback;
		}

		if (second.is_valid())
		{
			first.value = tmp;
		}
		else
		{
			first.value = &first.fallback;
		}
		swap(first.fallback, second.fallback);
	}

	carry_value& operator=(T const& other)
	{
		*value = (is_valid()) ? other : *value;
		return *this;
	}

	carry_value& operator=(carry_value<T> const& other)
	{
		*value = (is_valid()) ? *other.value : *value;
		return *this;
	}

	T* value;
	T fallback;

	operator T () const
	{
		return *value;
	}

	inline bool is_valid() const
	{
		return value != &fallback;
	}
};


namespace grid
{
	template<size_t D, size_t... Is>
	inline bool is_same(const len_type(&pos0)[D], const len_type(&pos1)[D], std::index_sequence<Is...>)
	{
		return ((pos0[Is] == pos1[Is]) && ...);
	}

	template<size_t D>
	inline bool is_same(const len_type(&pos0)[D], const len_type(&pos1)[D])
	{
		return is_same(pos0, pos1, std::make_index_sequence<D>{});
	}

	inline iter_type position(iter_type n, len_type dim, const iter_type stride)
	{
		return (n / stride) % dim;
	}

	template<size_t D>
	inline void get_grid_position(iter_type(&pos)[D], const len_type(&dims)[D], const len_type(&stride)[D], iter_type n)
	{
		for (iter_type i = 0; i < D; ++i)
		{
			pos[i] = (n / stride[i]) % dims[i];
		}
	}

	inline void get_grid_position(iter_type(&pos)[1], const len_type(&dims)[1], iter_type n)
	{
		pos[0] = n;
	}

	inline void get_grid_position(iter_type(&pos)[2], const len_type(&dims)[2], iter_type n)
	{
		pos[1] = n / dims[0];
		pos[0] = n - dims[0] * pos[1];
	}

	inline void get_grid_position(iter_type(&pos)[3], const len_type(&dims)[3], iter_type n)
	{
		pos[0] = n % dims[0];
		pos[1] = (n / dims[0]) % dims[1];
		pos[2] = (n / (dims[0] * dims[1]));
	}
	
	template<size_t D>
	inline void get_grid_position(iter_type(&pos)[D], const len_type* dims, iter_type n)
	{
		len_type stride[D]{};
		get_stride(stride, dims);
		for (iter_type i = 0; i < D; ++i)
		{
			pos[i] = (n / stride[i]) % dims[i];
		}
	}

	inline void get_grid_position_offset(iter_type(&pos)[1], const len_type(&dims)[1], const len_type(&offset)[1], iter_type n)
	{
		pos[0] = n + offset[0];
	}

	inline void get_grid_position_offset(iter_type(&pos)[3], const len_type(&dims)[3], const len_type(&offset)[3], iter_type n)
	{
		pos[0] = n % dims[0] + offset[0];
		pos[1] = (n / dims[0]) % dims[1] + offset[1];
		pos[2] = (n / (dims[0] * dims[1])) + offset[2];
	}

	template<size_t D>
	void get_grid_position_offset(iter_type(&pos)[D], const len_type(&dims)[D], const len_type(&offset)[D], iter_type n)
	{
		get_grid_position(pos, dims, n);
		for (iter_type i = 0; i < D; ++i)
		{
			pos[i] += offset[i];
		}
	}

	template<size_t D>
	void get_grid_position_offset(iter_type(&pos)[D], const len_type(&dims)[D], const len_type(&stride)[D], const len_type(&offset)[D], iter_type n)
	{
		get_grid_position(pos, dims, stride, n);
		for (iter_type i = 0; i < D; ++i)
		{
			pos[i] += offset[i];
		}
	}

	inline void get_grid_position(pos_list& pos, dim_list const& dims, dim_list const& stride, iter_type n)
	{
		for (iter_type i = 0; i < pos.n; ++i)
		{
			pos[i] = (n / stride[i]) % dims[i];
		}
	}

	inline void get_grid_position_offset(pos_list& pos, dim_list const& dims, dim_list const& stride, pos_list const& offset, iter_type n)
	{
		get_grid_position(pos, dims, stride, n);
		for (iter_type i = 0; i < pos.n; ++i)
		{
			pos[i] += offset[i];
		}
	}

	template<Axis ax = Axis::X>
	void get_stride(len_type(&stride)[1], len_type const* dims)
	{
		stride[0] = 1;
	}

	template<Axis ax = Axis::X>
	void get_stride(len_type(&stride)[2], len_type const* dims)
	{
		if constexpr (ax == Axis::X)
		{
			stride[0] = 1;
			stride[1] = dims[0];
		}
		else
		{
			stride[0] = dims[0];
			stride[1] = -1;
		}
	}

	template<Axis ax = Axis::X>
	void get_stride(len_type(&stride)[3], len_type const* dims)
	{
		if constexpr (ax == Axis::X)
		{
			stride[0] = 1;
			stride[1] = dims[0];
			stride[2] = dims[0] * dims[1];
		}
		else if constexpr (ax == Axis::Y)
		{
			stride[0] = dims[0];
			stride[1] = -1;
			stride[2] = dims[0] * dims[1];
		}
		else
		{
			stride[0] = -dims[0] * dims[1];
			stride[1] = dims[0];
			stride[2] = -1;
		}
	}

	template<typename T, size_t D>
	void adjust_region_to_from(
		T* new_values, const iter_type(&new_origin)[D], const len_type(&new_dims)[D],
		const T* old_values, const iter_type(&old_origin)[D], const len_type(&old_dims)[D],
		const len_type(&global_dims)[D], T empty, len_type boundary_size = 0);

	template<typename T, size_t D>
	void adjust_origin_to_from(
		T* (&values), const iter_type(&new_origin)[D], const iter_type(&old_origin)[D],
		const len_type(&dims)[D], const len_type(&global_dims)[D], T empty, len_type boundary_size = 0);


	namespace
	{

		template<size_t D>
		using dim_arr_type = len_type[D];
		
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
		template<size_t I>
		void fill_inner_arr(iter_type* (&inner_i), len_type const* dimensions, iter_type& index);

		template<>
		inline void fill_inner_arr<0>(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			for (iter_type i = 0; i < dimensions[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++i)
			{
				*(inner_i++) = index++;
			}
		}

		template<>
		inline void fill_inner_arr<1>(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			iter_type inc = BOUNDARY_DEPTH + BOUNDARY_DEPTH;
			for (iter_type i = 0; i < dimensions[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++i, index += inc)
			{
				fill_inner_arr<0>(inner_i, dimensions, index);
			}
		}
		
		template<size_t I>
		void fill_inner_arr(iter_type* (&inner_i), len_type const* dimensions, iter_type& index)
		{
			iter_type inc = (BOUNDARY_DEPTH + BOUNDARY_DEPTH) * (dimensions[I - 2]);
			for (iter_type i = 0; i < dimensions[I] - BOUNDARY_DEPTH - BOUNDARY_DEPTH; ++i, index += inc)
			{
				fill_inner_arr<I - 1>(inner_i, dimensions, index);
			}
		}



		template<size_t D>
		iter_type get_init_index(len_type const* dimensions);

		template<>
		inline iter_type get_init_index<1>(len_type const*)
		{
			return BOUNDARY_DEPTH;
		}

		template<size_t D>
		iter_type get_init_index(len_type const* dimensions)
		{
			return grid::length<D - 1>(dimensions) * BOUNDARY_DEPTH + get_init_index<D - 1>(dimensions);
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


	template<size_t D>
	struct select_grid_index;

	template<>
	struct select_grid_index<1>
	{
		select_grid_index(const len_type(&dims)[1]) {}

		template<typename T>
		T& operator()(T* values, iter_type n)
		{
			return values[n];
		}

		template<typename T>
		const T& operator()(T const* values, iter_type n)
		{
			return values[n];
		}

		template<typename T>
		multi_value<1, T> operator()(T* (&values)[1], iter_type n)
		{
			multi_value<1, T> value;
			value[0] = values[0][n];
			return value;
		}

		template<typename T>
		multi_value<1, T> operator()(const T* (&values)[1], iter_type n)
		{
			multi_value<1, T> value;
			value[0] = values[0][n];
			return value;
		}
	};

	template<>
	struct select_grid_index<2>
	{
		select_grid_index(const len_type(&dims)[2]) : stride{ dims[0] } {}
		len_type stride;

		template<typename T>
		T& operator()(T* values, iter_type x, iter_type y) const
		{
			return values[y * stride + x];
		}

		template<typename T>
		const T& operator()(T const* values, iter_type x, iter_type y) const
		{
			return values[y * stride + x];
		}

		template<typename T>
		multi_value<2, T> operator()(T* (&values)[2], iter_type x, iter_type y)
		{
			multi_value<2, T> value;
			value[0] = operator()(values[0], x, y);
			value[1] = operator()(values[1], x, y);
			return value;
		}

		template<typename T>
		multi_value<2, T> operator()(const T* (&values)[2], iter_type x, iter_type y)
		{
			multi_value<2, T> value;
			value[0] = operator()(values[0], x, y);
			value[1] = operator()(values[1], x, y);
			return value;
		}
	};

	template<>
	struct select_grid_index<3>
	{
		select_grid_index(const len_type(&dims)[3]) : stride{ dims[0], dims[0] * dims[1] } {}
		len_type stride[2];

		template<typename T>
		T& operator()(T* values, iter_type x, iter_type y, iter_type z) const
		{
			return values[z * stride[1] + y * stride[0] + x];
		}

		template<typename T>
		const T& operator()(T const* values, iter_type x, iter_type y, iter_type z) const
		{
			return values[z * stride[1] + y * stride[0] + x];
		}

		template<typename T>
		multi_value<3, T> operator()(T* (&values)[3], iter_type x, iter_type y, iter_type z)
		{
			multi_value<3, T> value;
			value[0] = operator()(values[0], x, y, z);
			value[1] = operator()(values[1], x, y, z);
			value[2] = operator()(values[2], x, y, z);
			return value;
		}

		template<typename T>
		multi_value<3, T> operator()(const T* (&values)[3], iter_type x, iter_type y, iter_type z)
		{
			multi_value<3, T> value;
			value[0] = operator()(values[0], x, y, z);
			value[1] = operator()(values[1], x, y, z);
			value[2] = operator()(values[2], x, y, z);
			return value;
		}
	};

	template<size_t D>
	select_grid_index(const len_type(&)[D]) -> select_grid_index<D>;


	inline double distance(iter_type pos0, iter_type pos1)
	{
		double delta = pos0 - pos1;
		return delta * delta;
	}

	template<size_t D, size_t... Is>
	double distance(const iter_type(&pos0)[D], const iter_type(&pos1)[D], std::index_sequence<Is...>)
	{
		return std::sqrt((distance(pos0[Is], pos1[Is]) + ...));
	}

	template<size_t D>
	double distance(const iter_type(&pos0)[D], const iter_type(&pos1)[D])
	{
		return distance(pos0, pos1, std::make_index_sequence<D>{});
	}

	inline bool is_same_point(iter_type pos0, iter_type pos1)
	{
		return pos0 == pos1;
	}

	template<size_t D, size_t... Is>
	bool is_same_point(const iter_type(&pos0)[D], const iter_type(&pos1)[D], std::index_sequence<Is...>)
	{
		return (is_same_point(pos0[Is], pos1[Is]) && ...);
	}

	template<size_t D>
	bool is_same_point(const iter_type(&pos0)[D], const iter_type(&pos1)[D])
	{
		return is_same_point(pos0, pos1, std::make_index_sequence<D>{});
	}

	inline bool is_in_region(iter_type pos, const iter_type dims)
	{
		return (pos >= 0 && pos < dims);
	}

	template<size_t D, size_t... Is>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], std::index_sequence<Is...>)
	{
		return (is_in_region(pos[Is], dims[Is]) && ...);
	}

	template<size_t D>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D])
	{
		return is_in_region(pos, dims, std::make_index_sequence<D>{});
	}

	inline bool is_in_region_boundary(iter_type pos, iter_type dims, len_type boundary_size)
	{
		return (pos >= boundary_size && pos < dims - boundary_size);
	}

	template<size_t D, size_t... Is>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], len_type boundary_size, std::index_sequence<Is...>)
	{
		return (is_in_region_boundary(pos[Is], dims[Is], boundary_size) && ...);
	}

	template<size_t D>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], len_type boundary_size)
	{
		return is_in_region(pos, dims, boundary_size, std::make_index_sequence<D>{});
	}

	inline bool is_in_region_boundary(iter_type pos, iter_type dims, iter_type offset, len_type boundary_size)
	{
		return (pos >= offset + boundary_size && pos < offset + dims - boundary_size);
	}

	template<size_t D, size_t... Is>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], const iter_type(&offset)[D], len_type boundary_size, std::index_sequence<Is...>)
	{
		return (is_in_region_boundary(pos[Is], dims[Is], offset[Is], boundary_size) && ...);
	}

	template<size_t D>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], const iter_type(&offset)[D], len_type boundary_size)
	{
		return is_in_region(pos, dims, boundary_size, offset, std::make_index_sequence<D>{});
	}

	inline bool is_in_region(iter_type pos, iter_type left, iter_type right)
	{
		return (pos >= left && pos < right);
	}

	template<size_t D, size_t... Is>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&intervals)[D][2], std::index_sequence<Is...>)
	{
		return (is_in_region(pos[Is], intervals[Is][0], intervals[Is][1]) && ...);
	}

	template<size_t D>
	bool is_in_region(const iter_type(&pos)[D], const iter_type(&intervals)[D][2])
	{
		return is_in_region(pos, intervals, std::make_index_sequence<D>{});
	}

	template<size_t D, size_t... Is>
	bool is_in_region(const iter_type(&intervals0)[D][2], const iter_type(&intervals1)[D][2], std::index_sequence<Is...>)
	{
		return is_in_region((iter_type[D]) { intervals0[Is][0]... }, intervals1, std::make_index_sequence<D>{})
			&& is_in_region((iter_type[D]) { intervals0[Is][1]... }, intervals1, std::make_index_sequence<D>{});
	}

	template<size_t D>
	bool is_in_region(const iter_type(&intervals0)[D][2], const iter_type(&intervals1)[D][2])
	{
		return is_in_region(intervals0, intervals1, std::make_index_sequence<D>{});
	}

	inline bool not_in_region(const iter_type pos, const iter_type dims)
	{
		return (pos < 0 || pos >= dims);
	}

	template<size_t D, size_t... Is>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], std::index_sequence<Is...>)
	{
		return (not_in_region(pos[Is], dims[Is]) && ...);
	}

	template<size_t D>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D])
	{
		return not_in_region(pos, dims, std::make_index_sequence<D>{});
	}

	inline bool not_in_region_boundary(const iter_type pos, const iter_type dims, len_type boundary_size)
	{
		return (pos < boundary_size || pos >= dims - boundary_size);
	}

	template<size_t D, size_t... Is>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], len_type boundary_size, std::index_sequence<Is...>)
	{
		return (not_in_region_boundary(pos[Is], dims[Is], boundary_size) || ...);
	}

	template<size_t D>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&dims)[D], len_type boundary_size)
	{
		return not_in_region(pos, dims, boundary_size, std::make_index_sequence<D>{});
	}

	inline bool not_in_region(iter_type pos, iter_type left, iter_type right)
	{
		return (pos < left || pos >= right);
	}

	template<size_t D, size_t... Is>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&intervals)[D][2], std::index_sequence<Is...>)
	{
		return (not_in_region(pos[Is], intervals[Is][0], intervals[Is][1]) || ...);
	}

	template<size_t D>
	bool not_in_region(const iter_type(&pos)[D], const iter_type(&intervals)[D][2])
	{
		return not_in_region(pos, intervals, std::make_index_sequence<D>{});
	}

	template<size_t D, bool... Fs>
	auto compute_periodic_stride(len_type(&stride)[D], const len_type(&dims)[D], std::integer_sequence<bool, Fs...>)
	{
		iter_type extend[]{ ((Fs) ? 2 : 1)... };

		stride[0] = 1;
		for (iter_type i = 1; i < D; ++i)
		{
			stride[i] = stride[i - 1] * dims[i - 1] * extend[i - 1];
		}
	}

	template<size_t D, size_t N, size_t... Is>
	auto compute_periodic_strides(len_type(&strides)[N][D], const len_type(&dims)[D], std::index_sequence<Is...>)
	{
		(compute_periodic_stride(strides[Is], dims, symphas::lib::nth_periodic_shift_t<Is, D>{}), ...);
	}

	template<size_t D, size_t N>
	auto compute_periodic_strides(len_type(&strides)[N][D], const len_type(&dims)[D])
	{
		compute_periodic_strides(strides, dims, std::make_index_sequence<N>{});
	}

	//template<size_t I, size_t D, bool... Fs>
	//auto compute_periodic_offset0(len_type(&offset), const len_type(&dims)[D], std::integer_sequence<bool, Fs...>)
	//{
	//	iter_type extend[]{ ((Fs) ? 1 : 0)... };
	//	iter_type stride[D];
	//	iter_type stride0[D];

	//	compute_periodic_stride(stride, dims, std::integer_sequence<bool, Fs...>{});
	//	get_stride(stride0, dims);

	//	for (iter_type i = 0; i < D; ++i)
	//	{
	//		offset += dims[i] * stride[i]; //dims[i] * extend[i] * stride[i];
	//	}
	//}

	//template<size_t D, size_t N, size_t... Is>
	//auto compute_periodic_offsets(len_type(&offsets)[N][2], const len_type(&dims)[D], std::index_sequence<Is...>)
	//{
	//	(compute_periodic_offset0(offsets[Is][0], dims, get_periodic_shift<Is, D>()), ...);
	//	(compute_periodic_offset1(offsets[Is][1], dims, get_periodic_shift<Is, D>()), ...);
	//}

	//template<size_t D, size_t N>
	//auto compute_periodic_offsets(len_type(&offsets)[N][2], const len_type(&dims)[D])
	//{
	//	compute_periodic_offsets(offsets, dims, std::make_index_sequence<N>{});
	//}

	template<size_t D>
	struct select_region
	{
		select_region() : offset{}, stride{}, origin{}, dims{}, len{}, boundary_size{} {}

		select_region(const len_type(&origin)[D], const len_type(&dims)[D], len_type boundary_size = 0) :
			offset{}, stride{}, origin{}, dims{}, len{}, boundary_size{ boundary_size }
		{
			update(origin, dims);
		}

		select_region(const len_type(&dims)[D], len_type boundary_size = 0) :
			stride{}, origin{}, len{}, boundary_size{ boundary_size }
		{
			iter_type origin[D]{};
			update(origin, dims);
		}
		
		auto update(const iter_type(&origin)[D])
		{
			for (iter_type i = 0; i < D; ++i)
			{
				this->origin[i] = origin[i];
			}
			offset = grid::index_from_position(origin, stride);
		}

		auto update(const iter_type(&origin)[D], const iter_type(&dims)[D])
		{
			this->stride[0] = 1;
			this->origin[0] = origin[0];
			this->dims[0] = dims[0];
			for (iter_type i = 1; i < D; ++i)
			{
				this->stride[i] = stride[i - 1] * dims[i - 1];
				this->origin[i] = origin[i];
				this->dims[i] = dims[i];
			}
			len = grid::length<D>(dims);
			offset = grid::index_from_position(origin, stride);
		}


		template<typename T>
		inline carry_value<T> operator()(iter_type(&pos)[D], T* values, const iter_type(&domain_dims)[D], const T& empty) const
		{
			for (iter_type i = 0; i < D; ++i)
			{
				pos[i] = (pos[i] >= origin[i] + boundary_size) ? pos[i] - origin[i] : pos[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
			}
			if (is_in_region(pos, dims, boundary_size))
			{
				return { &values[grid::index_from_position(pos, stride)], empty };
			}
			else
			{
				return { empty };
			}
		}

		template<typename T>
		inline carry_value<const T> operator()(iter_type(&pos)[D], const T* values, const iter_type(&domain_dims)[D], const T& empty) const
		{
			for (iter_type i = 0; i < D; ++i)
			{
				pos[i] = (pos[i] >= origin[i] + boundary_size) ? pos[i] - origin[i] : pos[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
			}
			if (is_in_region(pos, dims, boundary_size))
			{
				return { &values[grid::index_from_position(pos, stride)], empty };
			}
			else
			{
				return { empty };
			}
		}

		template<typename T>
		inline T& operator()(iter_type n, T* values, const iter_type(&domain_dims)[D], const T& empty) const
		{
			return values[(n - offset) % grid::length<D>(domain_dims)];
		}

		template<typename T>
		inline const T& operator()(iter_type n, const T* values, const iter_type(&domain_dims)[D], const T& empty) const
		{
			return values[(n - offset) % grid::length<D>(domain_dims)];
		}

		template<typename T>
		inline decltype(auto) operator()(T* values, const iter_type(&domain_dims)[D], const T& empty, const iter_type(&pos)[D]) const
		{
			iter_type pos0[D]{};
			for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
			return operator()(pos0, values, domain_dims, empty);
		}

		template<typename T>
		inline decltype(auto) operator()(const T* values, const iter_type(&domain_dims)[D], const T& empty, const iter_type(&pos)[D]) const
		{
			iter_type pos0[D]{};
			for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
			return operator()(pos0, values, domain_dims, empty);
		}

		//template<typename T>
		//carry_value<const T> operator()(T const* values, const iter_type(&origin)[D], const iter_type(&dims)[D], const T& empty, iter_type (&pos)[D]) const
		//{
		//	for (iter_type i = 0; i < D; ++i)
		//	{
		//		pos[i] -= origin[i];
		//	}
		//
		//	for (iter_type i = 0; i < D; ++i)
		//	{
		//		pos[i] = (pos[i] >= 0) ? pos[i] : ((pos[i] - boundary_size + domain_dims[i]) % domain_dims[i] + boundary_size);
		//	}
		//
		//	if (is_in_region(pos, dims, boundary_size))
		//	{
		//		return &values[grid::index_from_position(pos, stride)];
		//	}
		//	else
		//	{
		//		return empty;
		//	}
		//}

		template<typename T>
		inline multi_value<D, T> operator()(iter_type n, T* (&values)[D], const iter_type(&domain_dims)[D], T(&empty)[D]) const
		{
			multi_value<D, T> value;
			for (iter_type i = 0; i < D; ++i)
			{
				value[i] = operator()(n, values[i], domain_dims, empty[i]);
			}
			return value;
		}

		template<typename T>
		inline multi_value<D, T> operator()(iter_type n, T* const (&values)[D], const iter_type(&domain_dims)[D], const T(&empty)[D]) const
		{
			multi_value<D, T> value;
			for (iter_type i = 0; i < D; ++i)
			{
				value[i] = operator()(n, values[i], domain_dims, empty[i]);
			}
			return value;
		}
		

		template<typename T>
		inline multi_value<D, T> operator()(iter_type(&pos)[D], T* (&values)[D], const iter_type(&domain_dims)[D], const T(&empty)[D]) const
		{
			multi_value<D, T> value;
			for (iter_type i = 0; i < D; ++i)
			{
				pos[i] = (pos[i] >= origin[i] + boundary_size) ? pos[i] - origin[i] : pos[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
			}
			if (is_in_region(pos, dims, boundary_size))
			{
				return multi_value<D, T>(values, grid::index_from_position(pos, stride));
			}
			else
			{
				return { empty };
			}
		}

		template<typename T>
		inline multi_value<D, T> operator()(iter_type(&pos)[D], T* const (&values)[D], const iter_type(&domain_dims)[D], const T(&empty)[D]) const
		{
			for (iter_type i = 0; i < D; ++i)
			{
				pos[i] = (pos[i] >= origin[i] + boundary_size) ? pos[i] - origin[i] : pos[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
			}
			if (is_in_region(pos, dims, boundary_size))
			{
				return multi_value<D, T>(values, grid::index_from_position(pos, stride));
			}
			else
			{
				return { empty };
			}
		}

		template<typename T>
		inline multi_value<D, T> operator()(T* (&values)[D], const iter_type(&domain_dims)[D], const T(&empty)[D], const iter_type(&pos)[D]) const
		{
			iter_type pos0[D]{};
			for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
			return operator()(pos0, values, domain_dims, empty);
		}

		template<typename T>
		inline multi_value<D, T> operator()(const T* (&values)[D], const iter_type(&domain_dims)[D], const T(&empty)[D], const iter_type(&pos)[D]) const
		{
			iter_type pos0[D]{};
			for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
			return operator()(pos0, values, domain_dims, empty);
		}


		template<typename T, typename T0, typename... Ts, std::enable_if_t<(std::is_same<Ts, iter_type>::value && ...), int> = 0>
		decltype(auto) operator()(T&& values, const iter_type(&domain_dims)[D], T0 empty, iter_type coord0, Ts&&... coords) const
		{
			iter_type pos[D]{ coord0, std::forward<Ts>(coords)... };
			return operator()(pos, std::forward<T>(values), domain_dims, empty);
		}

		len_type offset;
		len_type stride[D];
		iter_type origin[D];
		len_type dims[D];
		len_type len;
		len_type boundary_size;
	};

	template<size_t D>
	select_region(const len_type(&)[D], const len_type(&)[D], len_type boundary_size = 0) -> select_region<D>;

	template<size_t D>
	len_type length_interior(select_region<D> const& region)
	{
		return length_interior<D>(region.dims, region.boundary_size);
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
	explicit Block(grid::dim_list const& dims) : Block(dims.data, dims.n) {}

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
	const T& operator[](iter_type i) const
	{
		return values[i];
	}

	//! Return the value at the \p i index in the list.
	/*!
	 * Return the value from the list from index \p i.
	 *
	 * \param i The index of the value to return.
	 */
	T& operator[](iter_type i)
	{
		return values[i];
	}

	operator const T* () const
	{
		return values;
	}

	operator T* ()
	{
		return values;
	}

	~Block()
	{
		delete[] values;
	}

protected:

	Block() : values{ nullptr }, len{ 0 } {}


};

template<size_t N, typename T>
struct multi_value
{
	T* value[N];

	multi_value(T* (&value)[N], iter_type offset = 0) : multi_value()
	{
		for (iter_type i = 0; i < N; ++i)
		{
			this->value[i] = value[i] + offset;
		}
	}

	multi_value(T* const (&value)[N], iter_type offset = 0) : multi_value()
	{
		for (iter_type i = 0; i < N; ++i)
		{
			this->value[i] = value[i] + offset;
		}
	}

	multi_value(T* const value) : multi_value()
	{
		for (iter_type i = 0; i < N; ++i)
		{
			this->value[i] = &value[i];
		}
	}

	multi_value(T(&value)[N]) : multi_value()
	{
		for (iter_type i = 0; i < N; ++i)
		{
			this->value[i] = &value[i];
		}
	}

	multi_value(const T(&value)[N]) : multi_value()
	{
		for (iter_type i = 0; i < N; ++i)
		{
			this->value[i] = &const_cast<T&>(value[i]);
		}
	}

	const T& operator[](iter_type i) const
	{
		return *(value[i]);
	}

	T& operator[](iter_type i)
	{
		return *(value[i]);
	}

	bool is_valid()
	{
		return true;
	}

	//! Return a multi_value as a vector for compatibility.
	operator any_vector_t<T, N>() const;

	//! Set the values of the multi_value from a vector.
	/*!
	 * Assigning the value from a vector can update the values
	 * in the MultiBlock list, using the pointers stored by
	 * the multi_value object.
	 */
	multi_value& operator=(any_vector_t<T, N> const& vector);
	multi_value& operator=(any_vector_t<T, N>&& vector);
	multi_value& operator=(multi_value<N, T> const& other);
	multi_value& operator=(multi_value<N, T>&& other);
	multi_value& operator=(T const& other);

	multi_value(multi_value<N, T> const& other) : multi_value(other.value) {}
	multi_value(multi_value<N, T>&& other) : multi_value(other.value) {}

	//! Set the values of the multi_value from a vector.
	/*!
	 * Assigning the value from a vector can update the values
	 * in the MultiBlock list, using the pointers stored by
	 * the multi_value object.
	 */
	 //multi_value& operator=(multi_value other);

	multi_value() : value{} {}


};

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 *
 * \tparam T The value type of the underlying array.
 */
template<size_t N, typename T>
struct MultiBlock
{

	T* values[N];	//!< The list of values managed by this object.
	len_type len;	//!< The number of values in the list.


	//! Create this object with \p len values.
	/*!
	 * Create this object with \p len values. The length can be 0, in which case
	 * no memory will be allocated, no values can be accessed, but the
	 * object can still be constructed.
	 *
	 * \param len The number of values to create.
	 */
	MultiBlock(len_type len) : values{ 0 }, len{ len }
	{
		for (iter_type n = 0; n < N; ++n)
		{
			values[n] = (len > 0) ? new T[len] : nullptr;
			std::fill(values[n], values[n] + len, T{});
		}
	}

	explicit MultiBlock(const len_type* len, size_t dimensions = 1) : MultiBlock((len != nullptr) ? grid::length(len, dimensions) : 0) {}

	MultiBlock(MultiBlock<N, T> const& other) : MultiBlock(other.len)
	{
		for (iter_type n = 0; n < N; ++n)
		{
			std::copy(other.values[n], other.values[n] + other.len, values[n]);
		}
	}

	MultiBlock(MultiBlock<N, T>&& other) noexcept : MultiBlock()
	{
		swap(*this, other);
	}

	MultiBlock& operator=(MultiBlock<N, T> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(MultiBlock<N, T>& first, MultiBlock<N, T>& second)
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
	multi_value<N, T> operator[](iter_type i) const
	{
		multi_value<N, T> value;
		for (iter_type n = 0; n < N; ++n)
		{
			value.value[n] = values[n] + i;
		}
		return value;
	}

	T* operator()(Axis ax) const
	{
		return values[symphas::axis_to_index(ax)];
	}

	~MultiBlock()
	{
		for (iter_type n = 0; n < N; ++n)
		{
			delete[] values[n];
		}
	}


protected:

	MultiBlock() : values{ 0 }, len{ 0 } {}

};


template<size_t N, typename T>
multi_value<N, T>::operator any_vector_t<T, N>() const
{
	any_vector_t<T, N> vector;

	for (iter_type i = 0; i < N; ++i)
	{
		vector[i] = *(value[i]);
	}

	return vector;
}


template<size_t N, typename T>
multi_value<N, T>& multi_value<N, T>::operator=(any_vector_t<T, N> const& vector)
{
	for (iter_type i = 0; i < N; ++i)
	{
		*(value[i]) = vector[i];
	}

	return *this;
}

template<size_t N, typename T>
multi_value<N, T>& multi_value<N, T>::operator=(any_vector_t<T, N>&& vector)
{
	for (iter_type i = 0; i < N; ++i)
	{
		*(value[i]) = vector[i];
	}

	return *this;
}

template<size_t N, typename T>
multi_value<N, T>& multi_value<N, T>::operator=(multi_value<N, T> const& other)
{
	any_vector_t<T, N> vector = other;
	*this = vector;
	return *this;
}

template<size_t N, typename T>
multi_value<N, T>& multi_value<N, T>::operator=(multi_value<N, T>&& other)
{
	any_vector_t<T, N> vector = other;
	*this = vector;
	return *this;
}

template<size_t N, typename T>
multi_value<N, T>& multi_value<N, T>::operator=(T const& other)
{
	any_vector_t<T, N> vector;;
	for (iter_type i = 0; i < N; ++i)
	{
		vector[i] = other;
	}
	*this = vector;
	return *this;
}


template<size_t N, typename T, typename V>
auto operator*(multi_value<N, T> const& a, V&& b)
{
	return any_vector_t<T, N>(a) * std::forward<V>(b);
}

template<size_t N, typename T, typename V>
auto operator*(V&& a, multi_value<N, T> const& b)
{
	return std::forward<V>(a) * any_vector_t<T, N>(b);
}

template<size_t N, typename T, typename V>
auto operator+(multi_value<N, T> const& a, multi_value<N, V> const& b)
{
	return any_vector_t<T, N>(a) + any_vector_t<V, N>(a);
}

template<size_t N, typename T, typename V>
auto operator-(multi_value<N, T> const& a, multi_value<N, V> const& b)
{
	return any_vector_t<T, N>(a) -any_vector_t<V, N>(a);
}

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
private:

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

	void set_dimensions(grid::dim_list dimensions)
	{
		std::fill(dims, dims + D, 0);
		for (iter_type i = 0; i < dimensions.n; ++i)
		{
			dims[i] = dimensions[i];
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
	Grid(grid::dim_list dimensions) : Block<T>{ grid::length<D>(dimensions) }, dims{ 0 }
	{
		set_dimensions(dimensions);
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

	template<typename... Ts, std::enable_if_t<(sizeof...(Ts) == D), int> = 0>
	decltype(auto) operator()(Ts&&... indices) const
	{
		return grid::select_grid_index(dims)(Block<T>::values, std::forward<Ts>(indices)...);
	}

protected:

	constexpr Grid() : Grid(nullptr) {}

};

template<typename T, size_t D>
struct Grid<any_vector_t<T, D>, D> : MultiBlock<D, T>
{
	using parent_type = MultiBlock<D, T>;
	using parent_type::parent_type;

private:

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
	Grid(grid::dim_list dimensions) :
		parent_type{ grid::length<D>(dimensions) }, dims{ 0 }
	{
		set_dimensions(dimensions);
	}

	//! Create a grid of the prescribed dimensions.
	/*!
	 * Creates a new grid using the given dimensions. The number of values in
	 * the grid directly correspond to the dimensions, equal to the product
	 * of all dimensions.
	 *
	 * \param dimensions The dimensions of the grid.
	 */
	Grid(const len_type* dimensions) :
		parent_type{ grid::length<D>(dimensions) }, dims{ 0 }
	{
		set_dimensions(dimensions);
	}

	const Grid<any_vector_t<T, D>, D>& as_grid() const
	{
		return *this;
	}

	Grid<any_vector_t<T, D>, D>& as_grid()
	{
		return *this;
	}

	const T* axis(Axis ax) const
	{
		return parent_type::values[symphas::axis_to_index(ax)];
	}

	T* axis(Axis ax)
	{
		return parent_type::values[symphas::axis_to_index(ax)];
	}

	template<typename... Ts, std::enable_if_t<(sizeof...(Ts) == D), int> = 0>
	decltype(auto) operator()(Ts&&... rest) const
	{
		return grid::select_grid_index(dims)(MultiBlock<D, T>::values, std::forward<Ts>(rest)...);
	}

protected:

	constexpr Grid() : Grid(nullptr) {}

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
	using parent_type = Grid<T, D>;
	using parent_type::dims;

	BoundaryGrid(grid::dim_list dimensions) : BoundaryGrid(static_cast<const len_type*>(dimensions)) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, D>{ dimensions } {}
	BoundaryGrid(Grid<T, D> const& other) : Grid<T, D>{ other } {}
	BoundaryGrid(Grid<T, D>&& other) noexcept : Grid<T, D>{ other } {}

	const BoundaryGrid<T, D>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, D>& as_grid()
	{
		return *this;
	}

	template<typename... Ts, std::enable_if_t<(sizeof...(Ts) == D), int> = 0>
	decltype(auto) operator()(Ts&&... rest) const
	{
		return grid::select_grid_index(dims)(Block<T>::values, (std::forward<Ts>(rest) + BOUNDARY_DEPTH)...);
	}

	constexpr BoundaryGrid() : Grid<T, D>() {}
};


template <typename T>
struct BoundaryGrid<T, 3> : Grid<T, 3>
{
	using parent_type = Grid<T, 3>;
	using parent_type::dims;

	BoundaryGrid(grid::dim_list dimensions) : BoundaryGrid(static_cast<const len_type*>(dimensions)) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 3>(dimensions) {}
	BoundaryGrid(Grid<T, 3> const& other) : Grid<T, 3>(other) {}
	BoundaryGrid(Grid<T, 3>&& other) noexcept : Grid<T, 3>(other) {}

	const BoundaryGrid<T, 3>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 3>& as_grid()
	{
		return *this;
	}

	decltype(auto) operator()(iter_type x, iter_type y, iter_type z) const
	{
		return grid::select_grid_index(dims)(parent_type::values, x + BOUNDARY_DEPTH, y + BOUNDARY_DEPTH, z + BOUNDARY_DEPTH);
	}

	constexpr BoundaryGrid() : Grid<T, 3>() {}

};


template <typename T>
struct BoundaryGrid<T, 2> : Grid<T, 2>
{
	using parent_type = Grid<T, 2>;
	using parent_type::dims;

	BoundaryGrid(grid::dim_list dimensions) : Grid<T, 2>(dimensions) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 2>(dimensions) {}
	BoundaryGrid(Grid<T, 2> const& other) : Grid<T, 2>(other) {}
	BoundaryGrid(Grid<T, 2>&& other) noexcept : Grid<T, 2>(other) {}


	const BoundaryGrid<T, 2>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 2>& as_grid()
	{
		return *this;
	}

	decltype(auto) operator()(iter_type x, iter_type y) const
	{
		return grid::select_grid_index(dims)(parent_type::values, x + BOUNDARY_DEPTH, y + BOUNDARY_DEPTH);
	}

	constexpr BoundaryGrid() : Grid<T, 2>() {}
};


template <typename T>
struct BoundaryGrid<T, 1> : Grid<T, 1>
{
	using parent_type = Grid<T, 1>;
	using parent_type::dims;

	BoundaryGrid(grid::dim_list dimensions) : BoundaryGrid(static_cast<const len_type*>(dimensions)) {}
	BoundaryGrid(const len_type* dimensions) : Grid<T, 1>(dimensions) {}
	BoundaryGrid(Grid<T, 1> const& other) : Grid<T, 1>(other) {}
	BoundaryGrid(Grid<T, 1>&& other) noexcept : Grid<T, 1>(other) {}

	decltype(auto) invalue(iter_type i)
	{
		return parent_type::operator[](i + BOUNDARY_DEPTH);
	}

	const BoundaryGrid<T, 1>& as_grid() const
	{
		return *this;
	}

	BoundaryGrid<T, 1>& as_grid()
	{
		return *this;
	}


	decltype(auto) operator()(iter_type x) const
	{
		return grid::select_grid_index(dims)(parent_type::values, x + BOUNDARY_DEPTH);
	}


	constexpr BoundaryGrid() : Grid<T, 1>() {}
};



// ***********************************************************************************************

#define REGIONAL_GRID_OUTER_VALUE 0

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type.
 * This grid implementation is an extension of the base ::Grid, but it
 * is meant to sub-domain a larger domain when there are many fields.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template<typename T, size_t D>
struct RegionalGrid : Grid<T, D>
{
	using parent_type = Grid<T, D>;

	grid::select_region<D> region;
	T empty;

public:
	RegionalGrid(grid::dim_list dimensions, T empty = REGIONAL_GRID_OUTER_VALUE, len_type boundary_size = BOUNDARY_DEPTH) :
		RegionalGrid(static_cast<const len_type*>(dimensions), empty, boundary_size) {}
	RegionalGrid(const len_type* dimensions, T empty = REGIONAL_GRID_OUTER_VALUE, len_type boundary_size = BOUNDARY_DEPTH) :
		parent_type{ dimensions }, region{ parent_type::dims, boundary_size }, empty{ empty } {}

	RegionalGrid(Grid<T, D> const& other) :
		parent_type{ other }, region{ parent_type::dims, BOUNDARY_DEPTH }, empty{ REGIONAL_GRID_OUTER_VALUE } {}
	RegionalGrid(Grid<T, D>&& other) noexcept :
		parent_type{ other }, region{ parent_type::dims, BOUNDARY_DEPTH }, empty{ REGIONAL_GRID_OUTER_VALUE } {}

	RegionalGrid(RegionalGrid<T, D> const& other) : RegionalGrid(nullptr, other.empty, other.region.boundary_size)
	{
		using std::swap;

		parent_type::len = other.len;
		std::copy(other.dims, other.dims + D, parent_type::dims);
		region = other.region;


		Block<T> tmp(other.region.len);
		std::copy(other.values, other.values + other.region.len, tmp.values);
		swap(tmp.values, parent_type::values);
	}

	RegionalGrid(RegionalGrid<T, D>&& other) noexcept : RegionalGrid()
	{
		swap(*this, other);
	}
	RegionalGrid& operator=(RegionalGrid<T, D> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(RegionalGrid<T, D>& first, RegionalGrid<T, D>& second)
	{
		using std::swap;
		swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
		swap(first.empty, second.empty);
		swap(first.region, second.region);
	}


	const RegionalGrid<T, D>& as_grid() const
	{
		return *this;
	}

	RegionalGrid<T, D>& as_grid()
	{
		return *this;
	}

	template<typename... Ts, std::enable_if_t<(sizeof...(Ts) == D), int> = 0>
	decltype(auto) operator()(Ts&&... rest) const
	{
		return region(Block<T>::values, parent_type::dims, empty, (std::forward<Ts>(rest) + region.boundary_size)... );
	}

	decltype(auto) operator[](iter_type n) const
	{
		iter_type pos[D]{};
		grid::get_grid_position(pos, parent_type::dims, n);
		return region(pos, Block<T>::values, parent_type::dims, empty);
	}

	T const& get_unsafe(iter_type n) const
	{
		return region(n, Block<T>::values, parent_type::dims, empty);
	}

	T& get_unsafe(iter_type n)
	{
		return region(n, Block<T>::values, parent_type::dims, empty);
	}

	void adjust(const iter_type(&new_origin)[D])
	{
		grid::adjust_origin_to_from(parent_type::values, new_origin, region.origin, region.dims, parent_type::dims, empty, region.boundary_size);
		region.update(new_origin);
	}

	void adjust(const iter_type(&new_origin)[D], const len_type(&new_dims)[D])
	{
		T* new_values = new T[grid::length<D>(new_dims)]{};
		std::fill(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			new_values, new_values + grid::length<D>(new_dims), empty);
		grid::adjust_region_to_from(new_values, new_origin, new_dims, parent_type::values, region.origin, region.dims, parent_type::dims, empty, region.boundary_size);

		using std::swap;
		swap(parent_type::values, new_values);
		delete[] new_values;
		
		region.update(new_origin, new_dims);
	}

	void adjust(grid::select_region<D> const& other)
	{
		if (!grid::is_same(region.dims, other.dims)
			|| !grid::is_same(region.origin, other.origin))
		{
			if (grid::is_same(region.dims, other.dims))
			{
				adjust(other.origin);
			}
			else
			{
				adjust(other.origin, other.dims);
			}
		}
	}

protected:

	constexpr RegionalGrid() : parent_type(), region{}, empty{} {}
};

template<typename T, size_t D>
struct RegionalGrid<any_vector_t<T, D>, D> : Grid<any_vector_t<T, D>, D>
{
	using parent_type = Grid<any_vector_t<T, D>, D>;

	grid::select_region<D> region;
	T empty[D];

public:

	template<size_t... Is>
	RegionalGrid(const len_type* dimensions, const T(&empty)[D], len_type boundary_size, std::index_sequence<Is...>) :
		parent_type{ dimensions }, region{ parent_type::dims, boundary_size }, empty{ empty[Is]... } {}
	RegionalGrid(const len_type* dimensions, const T(&empty)[D], len_type boundary_size = BOUNDARY_DEPTH) :
		RegionalGrid(static_cast<const len_type*>(dimensions), empty, boundary_size, std::make_index_sequence<D>{}) {}
	RegionalGrid(grid::dim_list dimensions, const T(&empty)[D], len_type boundary_size = BOUNDARY_DEPTH) :
		RegionalGrid(static_cast<const len_type*>(dimensions), empty, boundary_size) {}

	template<size_t... Is>
	RegionalGrid(const len_type* dimensions, T empty, len_type boundary_size, std::index_sequence<Is...>) :
		RegionalGrid(dimensions, (T[D]) { symphas::internal::repeat_value<Is>(empty)... }, boundary_size, std::index_sequence<Is...>{}) {}
	RegionalGrid(const len_type* dimensions, T empty = REGIONAL_GRID_OUTER_VALUE, len_type boundary_size = BOUNDARY_DEPTH) :
		RegionalGrid(dimensions, empty, boundary_size, std::make_index_sequence<D>{}) {}
	RegionalGrid(grid::dim_list dimensions, T empty = REGIONAL_GRID_OUTER_VALUE, len_type boundary_size = BOUNDARY_DEPTH) :
		RegionalGrid(static_cast<const len_type*>(dimensions), empty, boundary_size) {}

	template<size_t... Is>
	RegionalGrid(Grid<T, D> const& other, std::index_sequence<Is...>) :
		parent_type{ other }, region{ parent_type::dims, BOUNDARY_DEPTH }, empty{ symphas::internal::repeat_value<Is>(REGIONAL_GRID_OUTER_VALUE)... } {}
	template<size_t... Is>
	RegionalGrid(Grid<T, D>&& other, std::index_sequence<Is...>) noexcept :
		parent_type{ other }, region{ parent_type::dims, BOUNDARY_DEPTH }, empty{ symphas::internal::repeat_value<Is>(REGIONAL_GRID_OUTER_VALUE)... } {}

	RegionalGrid(Grid<T, D> const& other) : RegionalGrid(other, std::make_index_sequence<D>{}) {}
	RegionalGrid(Grid<T, D>&& other) noexcept : RegionalGrid(other, std::make_index_sequence<D>{}) {}

	RegionalGrid(RegionalGrid<any_vector_t<T, D>, D> const& other) : RegionalGrid(nullptr, other.empty, other.region.boundary_size)
	{
		using std::swap;

		parent_type::len = other.len;
		std::copy(other.dims, other.dims + D, parent_type::dims);
		region = other.region;


		MultiBlock<D, T> tmp(other.region.len);
		for (iter_type i = 0; i < D; ++i)
		{
			std::copy(other.values[i], other.values[i] + other.region.len, tmp.values[i]);
		}
		swap(tmp.values, parent_type::values);
	}


	const RegionalGrid<any_vector_t<T, D>, D>& as_grid() const
	{
		return *this;
	}

	RegionalGrid<any_vector_t<T, D>, D>& as_grid()
	{
		return *this;
	}

	const T* axis(Axis ax) const
	{
		return parent_type::values[symphas::axis_to_index(ax)];
	}

	T* axis(Axis ax)
	{
		return parent_type::values[symphas::axis_to_index(ax)];
	}

	template<typename... Ts, std::enable_if_t<(sizeof...(Ts) == D), int> = 0>
	decltype(auto) operator()(Ts&&... rest) const
	{
		return region(MultiBlock<D, T>::values, parent_type::dims, empty, (std::forward<Ts>(rest) + region.boundary_size)...);
	}

	decltype(auto) operator[](iter_type n) const
	{
		iter_type pos[D]{};
		grid::get_grid_position(pos, parent_type::dims, n);
		return region(pos, MultiBlock<D, T>::values, parent_type::dims, empty);
	}

	T const& get_unsafe(iter_type n) const
	{
		return region(n, MultiBlock<D, T>::values, parent_type::dims, empty);
	}

	T& get_unsafe(iter_type n)
	{
		return region(n, MultiBlock<D, T>::values, parent_type::dims, empty);
	}

	void adjust(const iter_type(&new_origin)[D])
	{
		for (iter_type i = 0; i < D; ++i)
		{
			grid::adjust_origin_to_from(MultiBlock<D, T>::values[i], new_origin, region.origin, region.dims, parent_type::dims, empty[i], region.boundary_size);
		}
		region.update(new_origin);
	}

	void adjust(const iter_type(&new_origin)[D], const len_type(&new_dims)[D])
	{
		for (iter_type i = 0; i < D; ++i)
		{
			T* new_values = new T[grid::length<D>(new_dims)]{};

			std::fill(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par_unseq,
#endif
				new_values, new_values + grid::length<D>(new_dims), empty[i]);


			grid::adjust_region_to_from(new_values, new_origin, new_dims,
				MultiBlock<D, T>::values[i], region.origin, region.dims, parent_type::dims, empty[i], region.boundary_size);

			using std::swap;
			swap(MultiBlock<D, T>::values[i], new_values);
			delete[] new_values;
		}

		region.update(new_origin, new_dims);
	}

	void adjust(grid::select_region<D> const& other)
	{
		if (!grid::is_same(region.dims, other.dims)
			|| !grid::is_same(region.origin, other.origin))
		{
			if (grid::is_same(region.dims, other.dims))
			{
				adjust(other.origin);
			}
			else
			{
				adjust(other.origin, other.dims);
			}
		}
	}

	constexpr RegionalGrid() : parent_type(), region{}, empty{} {}
};

#ifdef USING_MPI

template<typename T, size_t D>
struct RegionalGridMPI : RegionalGrid<T, D>
{
	using parent_type = RegionalGrid<T, D>;
	using parent_type::parent_type;
	using parent_type::dims;
	using parent_type::len;
	using parent_type::region;
	
	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions) :
		parent_type(mpi_dims(info, dimensions))
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, T empty) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, T empty, len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions) :
		parent_type(mpi_dims(info, dimensions))
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, T empty) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, T empty, len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(RegionalGridMPI<T, D> const& other) : RegionalGridMPI(nullptr, other.empty, other.region.boundary_size)
	{
		using std::swap;

		parent_type::len = other.len;
		std::copy(other.dims, other.dims + D, parent_type::dims);
		region = other.region;

		if (other.values != nullptr)
		{
			Block<T> tmp(other.region.len);
			std::copy(other.values, other.values + other.region.len, tmp.values);
			swap(tmp.values, parent_type::values);
		}
	}

	auto operator[](iter_type n) const -> std::invoke_result_t<decltype(&parent_type::operator[]), parent_type, iter_type>
	{
		if (parent_type::values != nullptr) return parent_type::operator[](n); else return parent_type::empty;
	}

	const RegionalGridMPI<T, D>& as_grid() const
	{
		return *this;
	}

	RegionalGridMPI<T, D>& as_grid()
	{
		return *this;
	}
	
protected:

	void update_dimensions(const len_type* dimensions)
	{
		if (dimensions != nullptr)
		{
			std::copy(dimensions, dimensions + D, parent_type::dims);
			parent_type::len = grid::length<D>(dimensions);
			parent_type::region = { parent_type::dims, parent_type::region.boundary_size };
		}
	}

	auto mpi_dims(symphas::multi_thr_info_type const& info, const len_type* dimensions)
	{
		return (info.index_in_node()) ? dimensions : nullptr;
	}
};



template<typename T, size_t D>
struct RegionalGridMPI<any_vector_t<T, D>, D> : RegionalGrid<any_vector_t<T, D>, D>
{
	using parent_type = RegionalGrid<any_vector_t<T, D>, D>;
	using parent_type::parent_type;
	using parent_type::dims;
	using parent_type::len;
	using parent_type::region;


	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions) :
		parent_type(mpi_dims(info, dimensions))
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, T empty) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, const T(&empty)[D]) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, T empty, len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, const len_type* dimensions, const T(&empty)[D], len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions) :
		parent_type(mpi_dims(info, dimensions))
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, T empty) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, const T(&empty)[D]) :
		parent_type(mpi_dims(info, dimensions), empty)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, T empty, len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(symphas::multi_thr_info_type const& info, grid::dim_list dimensions, const T(&empty)[D], len_type boundary_size) :
		parent_type(mpi_dims(info, dimensions), empty, boundary_size)
	{
		update_dimensions(dimensions);
	}

	RegionalGridMPI(RegionalGridMPI<any_vector_t<T, D>, D> const& other) : RegionalGridMPI(nullptr, other.empty, other.region.boundary_size)
	{
		using std::swap;

		parent_type::len = other.len;
		std::copy(other.dims, other.dims + D, parent_type::dims);
		region = other.region;

		if (other.values != nullptr)
		{
			MultiBlock<D, T> tmp(other.region.len);
			for (iter_type i = 0; i < D; ++i)
			{
				std::copy(other.values[i], other.values[i] + other.region.len, tmp.values[i]);
			}
			swap(tmp.values, parent_type::values);
		}
	}

	auto operator[](iter_type n) const->std::invoke_result_t<decltype(&parent_type::operator[]), parent_type, iter_type>
	{
		if (parent_type::values != nullptr) return parent_type::operator[](n); else return parent_type::empty;
	}

	const RegionalGridMPI<any_vector_t<T, D>, D>& as_grid() const
	{
		return *this;
	}

	RegionalGridMPI<any_vector_t<T, D>, D>& as_grid()
	{
		return *this;
	}

protected:

	void update_dimensions(const len_type* dimensions)
	{
		if (dimensions != nullptr)
		{
			std::copy(dimensions, dimensions + D, parent_type::dims);
			parent_type::len = grid::length<D>(dimensions);
			parent_type::region = { parent_type::dims, parent_type::region.boundary_size };
		}
	}

	auto mpi_dims(symphas::multi_thr_info_type const& info, const len_type* dimensions)
	{
		return (info.index_in_node()) ? dimensions : nullptr;
	}

};

#endif

// **************************************************************************************



namespace grid
{
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
				return { symphas::grid_info(vdata).get_dims() };
			}
		};

#ifdef USING_MPI

		template<>
		struct make<RegionalGridMPI>
		{
			//! Obtain a newly initialized grid using the interval data.
			/*
			 * This is used to generate a new grid of the prescribed dimensions
			 * and of the prescribed source type. The dimensions are assumed
			 * to all be zero, thereby creating an empty grid.
			 */
			template<typename T, size_t D>
			static RegionalGridMPI<T, D> apply()
			{
				len_type dims[D]{ 0 };
				return { {}, dims };
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
			static RegionalGridMPI<T, D> apply(symphas::multi_thr_info_type const& thr_info, symphas::interval_data_type const& vdata)
			{
				if (vdata.size() < D)
				{
					return apply<T, D>();
				}
				return { thr_info, symphas::grid_info(vdata).get_dims()};
			}
		};

#endif 
	}


	//! Create a new grid with of the given base type and dimension.
	/*!
	 * This is used to generate a new grid of the prescribed dimensions
	 * and of the prescribed source type.
	 * 
	 * \param vdata The intervals of the grid, from which the dimensions
	 * are taken.
	 */
	template<template<typename, size_t> typename G, typename T, size_t D, typename... Ts>
	G<T, D> construct(Ts&&... args)
	{
		return make<G>::template apply<T, D>(std::forward<Ts>(args)...);
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

		template<size_t N, typename T>
		static constexpr wrap_type<any_vector_t<T, N>> cast(MultiBlock<N, T>)
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





