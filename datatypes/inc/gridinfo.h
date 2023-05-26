
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
 * PURPOSE: Adds information about a grid, such as representation of sides
 * of a grid and axes of a grid.
 *
 * ***************************************************************************
 */

#pragma once

#include <map>

#include "definitions.h"
#include "params.h"

//! The thickness of a boundary on a BoundaryGrid.
/*! 
 * The number of layers on the outside of a BoundaryGrid which is considered
 * not to be part of the data of the grid itself, but the which contains and is
 * outside the data. The boundary is the same thickness on all sides of a grid.
 *
 * One layer is one cell deep into the logical arrangement of values of a grid.
 * The arrangement of elements constituting the boundary is equal to the exclusion
 * of the interior grid, which is a grid of the elements with dimensions equal to 
 * that of the containing grid less twice the #THICKNESS, from the containing grid
 */
#define THICKNESS 3


//! Values for labeling the sides of a grid.
/*!
 * Global names that are used to refer to the sides of a grid, in use cases such
 * as labeling the boundaries of a grid. The order of the enumerated values
 * 
 * The sides listed in this enumeration apply for grids up to 3 dimensions.
 */
enum class Side
{
	LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
};

namespace symphas
{


	//! Get the iteration index for the given ::Side.
	/*!
	 * Returns the index corresponding to the given ::Side, used
	 * in the context of obtaining an iterable variable for arranging
	 * a list according to ::Side values. Each side has a
	 * unique index between 0 and the total number of sides minus 1.
	 */
	inline constexpr iter_type side_to_index(Side side)
	{
		return static_cast<iter_type>(side);
	}

	//! Get the iteration index for the given ::Axis.
	/*!
	 * Returns the index corresponding to the given ::Axis, used
	 * in the context of obtaining an iterable variable for arranging
	 * a list according to ::Axis values. Each axis has a
	 * unique index between 0 and the total number of axes minus 1.
	 */
	inline constexpr iter_type axis_to_index(Axis axis)
	{
		return static_cast<iter_type>(axis);
	}

	//! Get the ::Side value corresponding to the given index.
	/*!
	 * Each side is associated an index, see symphas::side_to_index(Side). Given
	 * an iterated index between 0 and the number of sides minus 1, the
	 * corresponding side will be returned.
	 */
	inline constexpr Side index_to_side(iter_type i)
	{
		return static_cast<Side>(i);
	}

	//! Get the ::Axis value corresponding to the given index.
	/*!
	 * Each axis is associated an index, see symphas::axis_to_index(Side). Given
	 * an iterated index between 0 and the number of axes minus 1, the
	 * corresponding axis will be returned.
	 */
	inline constexpr Axis index_to_axis(iter_type i)
	{
		return static_cast<Axis>(i);
	}

	//! Representation of an interval along a grid edge.
	/*!
	 * Contains an array representing the left and right endpoints of an interval,
	 * used to characterize the axis positioning of the edges of a grid. Intervals
	 * are considered to all start from the same point; i.e. \f$x_0, y_0, z_0\f$
	 * for three dimensional systems. An interval corresponds to one axis.
	 *
	 * The interval will also return the number of discrete elements used to
	 * represent it, i.e. the size of the grid dimension corresponding to the
	 * interval. 
	 */
	struct interval_element_type
	{
	protected:
		double data[2]{ 0 };	//!< Left and right endpoints of the interval.
		double h = 1.0;			//!< Grid spacing of this interval.

		len_type direct_count() const
		{
			return std::lround(std::abs(data[1] - data[0]) / h) + 1;
		}

	public:




		//! Return number of discrete elements in the interval.
		/*!
		 * Computes and returns the number of discrete elements in this interval
		 * in a standardized way. The computation assumes that an interval always
		 * contains its endpoints; for example, an interval between 0 and 2 with
		 * a spatial width of 1 means there are 3 discrete elements in the
		 * interval (elements 0, 1 and 2).
		 */
		len_type get_count() const
		{
				return direct_count();
		}

		//! Returns the left endpoint of the interval.
		/*!
		 * The endpoint that is returned is always the physical end point with which
		 * the interval is initialized with.
		 */
		axis_coord_t left() const
		{
			return data[0];
		}

		//! Returns the right endpoint of the interval.
		/*!
		 * The endpoint that is returned is always the physical end point with which
		 * the interval is initialized with.
		 */
		axis_coord_t right() const
		{
			return data[1];
		}

		//! Returns the grid spacing.
		/*!
		 * The grid spacing is always computed with respect
		 * to the number of elements in the domain.
		 * That is, whether the boundaries are extended.
		 * Specifying the width directly always returns the correct
		 * width.
		 */
		double width() const
		{
			return h;
		}

		//! Return the length of the interval on the real axis.
		double length() const
		{
			return right() - left() + width();
		}

		//! Sets the beginning and end values of the interval. 
		/*
		 * Given the start and the end of the interval, which are the left and
		 * right endpoints, respectively, the interval data is initialized.
		 *
		 * The design of this initialization is such that intervals can always be
		 * used to initialize each other using the values returned by left() and
		 * right() member functions. 
		 * 
		 * \param left The left endpoint of the interval.
		 * \param right The right endpoint of the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_interval(double left, double right, double width = 0)
		{
			if (width > 0)
			{
				h = width;
				data[0] = left;
				data[1] = right;
			}
			else
			{
				set_count(left, right, get_count());
			}
		}

		//! Sets the beginning and end values of the interval. 
		/*
		 * See set_interval(double, double, double).
		 *
		 * This overload is given the number of points and will determine the
		 * spatial width of elements from the given number of elements
		 * constituting the interval.
		 * 
		 * The number of given elements is the same value that is
		 * produced by the function count.
		 *
		 * \param count The number of points in the interval.
		 * \param left The left endpoint of the interval.
		 * \param right The right endpoint of the interval.
		 */
		void set_count(double left, double right, len_type count)
		{
			double width = std::abs(right - left) / (count - 1);
			set_interval(left, right, width);
		}

		//! Resize the interval with the given number of elements, with constant width.
		/*!
		 * See set_interval(double, double, double).
		 *
		 * This overload is given the number of points (without the width),
		 * and will resize the interval such that a point is fixed
		 * along the interval that corresponds to the value
		 * of the given ratio (between 0 and 1).
		 *
		 * The number of given elements is the same value that is
		 * produced by the function count().
		 *
		 * \param count The number of points in the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_count_from_r(double width, len_type count, double ratio)
		{
			double length0 = ((double)count / get_count()) * (width / h) * length();
			double pos = ratio * length();

			double left = pos - length0 * ratio;
			double right = left + (count - 1) * width;
			set_interval(left, right, width);
		}

		//! Resize the interval with the given number of elements, with constant width.
		/*!
		 * See set_interval(double, double, double).
		 *
		 * This overload is given the number of points (without the width),
		 * and will resize the interval such that a point is fixed 
		 * along the interval that corresponds to the value
		 * of the given ratio (between 0 and 1).
		 *
		 * The number of given elements is the same value that is
		 * produced by the function count().
		 *
		 * \param count The number of points in the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_count_from_r(len_type count, double ratio)
		{
			set_count_from_r(h, count, ratio);
		}

		//! Resize the interval with the given number of elements, with constant width.
		/*!
		 * See set_interval(double, double, double).
		 *
		 * This overload is given the number of points (without the width),
		 * and will fix the interval at the center when resizing
		 *
		 * The number of given elements is the same value that is
		 * produced by the function count().
		 *
		 * \param count The number of points in the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_count(len_type count)
		{
			set_count_from_r(count, 0.5);
		}

		//! Resize the interval with the given number of elements, with constant width.
		/*!
		 * See set_interval(double, double, double).
		 *
		 * This overload is given the number of points (without the width),
		 * and will fix the interval at the center when resizing
		 *
		 * The number of given elements is the same value that is
		 * produced by the function count().
		 *
		 * \param count The number of points in the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_count(double width, len_type count)
		{
			set_count_from_r(width, count, 0.5);
		}

	};

	using interval_data_type = std::map<Axis, interval_element_type>;

	struct grid_info;
}

void swap(symphas::grid_info& first, symphas::grid_info& second);


namespace grid
{

	template<typename T>
	struct box_list
	{
		T* data;
		size_t n;

		box_list(T dim0, T dim1, T dim2) :
			box_list(dim0, dim1, dim2, (dim1 > 0 && dim2 > 0) ? 3 : (dim1 > 0 || dim2 > 0) ? 2 : 1) {}
		box_list(T dim0, T dim1) :
			box_list(dim0, dim1, 0, (dim1 == 0) ? 1 : 2) {}
		box_list(T dim0) : box_list(dim0, 0, 0, 1) {}
		box_list() : box_list(0, 0, 0, 0) {}

		box_list(symphas::interval_data_type const& intervals) :
			box_list(
				(intervals.find(Axis::X) != intervals.end()) ? intervals.at(Axis::X).get_count() : 0,
				(intervals.find(Axis::Y) != intervals.end()) ? intervals.at(Axis::Y).get_count() : 0,
				(intervals.find(Axis::Z) != intervals.end()) ? intervals.at(Axis::Z).get_count() : 0) {}

		box_list(box_list<T> const& other) noexcept : box_list(other.data, other.n) {}

		box_list(const T* data, size_t n) : data{ (n > 0) ? new T[n]{ 0 } : nullptr }, n{ n }
		{
			if (data)
			{
				for (iter_type i = 0; i < n; ++i)
				{
					this->data[i] = data[i];
				}
			}
		}

		//dim_list(len_type* data, size_t D) : dim_list(
		//	(data && D > 0) ? data[0] : 0,
		//	(data && D > 1) ? data[1] : 0,
		//	(data && D > 2) ? data[2] : 0,
		//	D) {}

		box_list(box_list<T>&& other) noexcept : box_list()
		{
			swap(*this, other);
		}

		~box_list()
		{
			delete[] data;
		}

		operator T* () const
		{
			return data;
		}

		operator T () const
		{
			return data[0];
		}

		std::tuple<T, T, T> _3() const
		{
			return { data[0], data[1], data[2] };
		}

		std::tuple<T, T> _2() const
		{
			return { data[0], data[1] };
		}

		std::tuple<T> _1() const
		{
			return { data[0] };
		}

		T& operator[](iter_type i)
		{
			return data[i];
		}

		const T& operator[](iter_type i) const
		{
			return data[i];
		}

		friend void swap(box_list<T>& first, box_list<T>& second)
		{
			using std::swap;
			swap(first.data, second.data);
			swap(first.n, second.n);
		}

		const T* begin() const
		{
			return data;
		}

		const T* end() const
		{
			return data + n;
		}

		T* begin()
		{
			return data;
		}

		T* end()
		{
			return data + n;
		}

	protected:


		box_list(T dim0, T dim1, T dim2, size_t n) : data{ (n > 0) ? new T[n] : nullptr }, n{ n }
		{
			if (n > 0) data[0] = dim0;
			if (n > 1) data[1] = dim1;
			if (n > 2) data[2] = dim2;
		}

	};

	using dim_list = box_list<len_type>;
	using h_list = box_list<double>;
}



//! Represents parameters of a grid.
/*! 
 * Stores the information about a grid, namely the data which fully defines the 
 * characteristics of a grid.
 */
struct symphas::grid_info
{
protected:

	grid_info() : intervals{} {}

	interval_data_type make_intervals(const len_type* dims, size_t dimension)
	{
		if (!dims)
		{
			return interval_data_type{};
		}
		else
		{
			interval_data_type v;
			for (iter_type i = 0; i < dimension; ++i)
			{
				symphas::interval_element_type interval;
				interval.set_count(1.0, dims[i]);
				v[symphas::index_to_axis(i)] = interval;
			}

			return v;
		}
	}

public:

	symphas::interval_data_type intervals;		//!< Extent of the grid in the spatial axes.

	//! Create the grid information using the intervals of the system.
	/*!
	 * Create the grid information using the intervals of the system.
	 * The dimension is equivalent to the
	 * number of axes represented by the interval data.
	 * 
	 * \param intervals The intervals which define the system.
	 */
	grid_info(interval_data_type const& intervals) : intervals{ intervals } {}

	//! Create the grid information using the intervals of the system.
	/*!
	 * Create the grid information using the dimensions of the system.
	 *
	 * \param dims The side lengths (number of elements in each direction)
	 * which define the system.
	 * \param dimension The dimension of the system
	 */
	grid_info(const len_type* dims, size_t dimension) : grid_info(make_intervals(dims, dimension)) {}

	//! Create the grid information using the intervals of the system.
	/*!
	 * Create the grid information using the dimensions of the system.
	 * Converts the dimension type to be the acceptable type.
	 *
	 * \param dims The side lengths (number of elements in each direction)
	 * which define the system.
	 * \param dimension The dimension of the system
	 */
	grid_info(const len_type* dims, int dimension) : grid_info(dims, static_cast<size_t>(dimension)) {}


	grid_info(grid_info const& other) : grid_info{ other.intervals } {}
	grid_info(grid_info&& other) noexcept : grid_info()
	{
		::swap(*this, other);
	}
	grid_info& operator=(grid_info other)
	{
		::swap(*this, other);
		return *this;
	}

	//! Return the interval on the given axis.
	/*!
	 * Return the interval on the given axis. This provides the end point
	 * values of the system on that axis.
	 * 
	 * \param ax The axis to get the interval from.
	 */
	auto& at(Axis ax)
	{
		return intervals.at(ax);
	}

	//! Return the interval on the given axis.
	/*!
	 * Return the interval on the given axis. This provides the end point
	 * values of the system on that axis.
	 *
	 * \param ax The axis to get the interval from.
	 */
	const auto& at(Axis ax) const
	{
		return intervals.at(ax);
	}

	//! Return the interval on the given axis.
	/*!
	 * Return the interval on the given axis. This provides the end point
	 * values of the system on that axis.
	 *
	 * \param ax The axis to get the interval from.
	 */
	const auto& operator[](Axis ax) const
	{
		return at(ax);
	}

	//! Return the interval on the given axis.
	/*!
	 * Return the interval on the given axis. This provides the end point
	 * values of the system on that axis.
	 *
	 * \param ax The axis to get the interval from.
	 */
	auto& operator[](Axis ax)
	{
		return at(ax);
	}

	//! Return the dimension of the system.
	/*!
	 * Return the dimension of the system, equivalent to the
	 * number of axes represented by the interval data.
	 */
	len_type dimension() const
	{
		return static_cast<len_type>(intervals.size());
	}

	//! Gives the number of cells in the entire grid in discrete space.
	/*!
	 * Gives the number of cells in the entire grid in discrete space.
	 * This is computed using the lengths of all the intervals.
	 */
	len_type num_points()
	{
		len_type length = 1;
		for (iter_type i = 0; i < dimension(); ++i)
		{
			length *= at(symphas::index_to_axis(i)).get_count();
		}
		return length;
	}


	//! Gives the spacing of the intervals.
	/*!
	 * Returns the spatial discretization of each of the intervals. This is the
	 * list of widths between discretized elements in each of the dimensions.
	 */
	grid::h_list get_widths() const
	{
		grid::h_list widths(nullptr, dimension());
		for (iter_type i = 0; i < dimension(); ++i)
		{
			widths[i] = at(symphas::index_to_axis(i)).width();
		}
		return widths;
	}


	//! Gives the dimensions of the grid.
	/*!
	 * Gives the number of cells along each dimension.
	 */
	grid::dim_list get_dims() const
	{
		grid::dim_list dims(nullptr, dimension());
		for (iter_type i = 0; i < dimension(); ++i)
		{
			dims[i] = at(symphas::index_to_axis(i)).get_count();
		}
		return dims;
	}

	//! Gives the dimensions of the grid.
	/*!
	 * Gives the number of cells along each dimension.
	 */
	grid::box_list<double> get_lengths() const
	{
		grid::box_list<double> lengths(nullptr, dimension());
		for (iter_type i = 0; i < dimension(); ++i)
		{
			lengths[i] = at(symphas::index_to_axis(i)).length();
		}
		return lengths;
	}

	len_type area() const
	{
		double area = 1;
		for (auto length : get_lengths())
		{
			area *= length;
		}
		return area;
	}

	operator symphas::interval_data_type() const
	{
		return intervals;
	}

	auto begin() const
	{
		return intervals.begin();
	}

	auto end() const
	{
		return intervals.end();
	}

	auto begin()
	{
		return intervals.begin();
	}

	auto end()
	{
		return intervals.end();
	}

};

inline void swap(symphas::grid_info& first, symphas::grid_info& second)
{
	using std::swap;
	swap(first.intervals, second.intervals);
}


//! \cond

/* definitions for accessing the given interval from the map named
 * intervals which has key type Axis and value type double[2]
 * the implementation of the value type is given below
 * this is used in the configuration and used as well in the writing utility
 */
#define INTERVAL_X0_AT(i) intervals[i].at(Axis::X).left()
#define INTERVAL_Xn_AT(i) intervals[i].at(Axis::X).right()
#define INTERVAL_Y0_AT(i) intervals[i].at(Axis::Y).left()
#define INTERVAL_Yn_AT(i) intervals[i].at(Axis::Y).right()
#define INTERVAL_Z0_AT(i) intervals[i].at(Axis::Z).left()
#define INTERVAL_Zn_AT(i) intervals[i].at(Axis::Z).right()
#define INTERVAL_X0 intervals.at(Axis::X).left()
#define INTERVAL_Xn intervals.at(Axis::X).right()
#define INTERVAL_Y0 intervals.at(Axis::Y).left()
#define INTERVAL_Yn intervals.at(Axis::Y).right()
#define INTERVAL_Z0 intervals.at(Axis::Z).left()
#define INTERVAL_Zn intervals.at(Axis::Z).right()

#define INTERVAL_Xh_AT(i) intervals[i].at(Axis::X).width()
#define INTERVAL_Yh_AT(i) intervals[i].at(Axis::Y).width()
#define INTERVAL_Zh_AT(i) intervals[i].at(Axis::Z).width()
#define INTERVAL_Xh intervals.at(Axis::X).width()
#define INTERVAL_Yh intervals.at(Axis::Y).width()
#define INTERVAL_Zh intervals.at(Axis::Z).width()

#define INTERVAL_Xc_AT(i) intervals[i].at(Axis::X).get_count()
#define INTERVAL_Yc_AT(i) intervals[i].at(Axis::Y).get_count()
#define INTERVAL_Zc_AT(i) intervals[i].at(Axis::Z).get_count()
#define INTERVAL_Xc intervals.at(Axis::X).get_count()
#define INTERVAL_Yc intervals.at(Axis::Y).get_count()
#define INTERVAL_Zc intervals.at(Axis::Z).get_count()


//! \endcond

