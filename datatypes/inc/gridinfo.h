
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
 * that of the containing grid less twice the #BOUNDARY_DEPTH, from the containing grid
 */
#define BOUNDARY_DEPTH 3

namespace symphas::internal
{

	inline len_type range_count(const double(&interval)[2], double h)
	{
		return std::lround(std::abs(interval[1] - interval[0]) / h) + 1;
	}

	//! Return the length of the interval on the real axis.
	inline double range_length(const double(&interval)[2], double h)
	{
		return interval[1] - interval[0] + h;
	}

	inline void set_interval(double(&interval)[2], double& h, double left, double right, double width = 0)
	{
		h = (width > 0) ? width : (h > 0) ? h : 1;
		interval[0] = left;
		interval[1] = right;
	}

	inline void set_interval(double(&interval)[2], double left, double right)
	{
		interval[0] = left;
		interval[1] = right;
	}

	inline void set_count(double(&interval)[2], double& h, double left, double right, len_type count)
	{
		double width = std::abs(right - left) / (count - 1);
		set_interval(interval, h, left, right, width);
	}

	inline void set_count_from_r(double(&interval)[2], double& h, len_type count, double width, double ratio)
	{
		double length0 = ((double)count / range_count(interval, h)) * (width / h) * range_length(interval, h);
		double pos = interval[0] + ratio * range_length(interval, h);

		double left = pos - length0 * ratio;
		double right = left + (count - 1) * width;
		set_interval(interval, h, left, right, width);
	}

	inline void set_count_from_r(double(&interval)[2], double& h, len_type count, double ratio)
	{
		double length0 = (double(count) / range_count(interval, h)) * range_length(interval, h);
		double pos = interval[0] + ratio * range_length(interval, h);

		double left = pos - length0 * ratio;
		double right = left + (count - 1) * h;
		set_interval(interval, left, right);
	}
}

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

	//! Get the iteration index for the given ::Axis.
	/*!
	 * Returns the index corresponding to the given ::Axis, used
	 * in the context of obtaining an iterable variable for arranging
	 * a list according to ::Axis values. Each axis has a
	 * unique index between 0 and the total number of axes minus 1.
	 */
	inline constexpr Axis axis_of_side(Side side)
	{
		return static_cast<Axis>(static_cast<iter_type>(side) / 2);
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
		double domain[2];		//!< Left and right endpoints of the domain.
		double data[2];			//!< Interval of the axis element grid in the domain.
		double h;				//!< Grid spacing of this interval.

	public:

		interval_element_type(const double(&data)[2], double h = 1.) : interval_element_type(data, data, h) {}
		interval_element_type(const double(&domain)[2], const double(&data)[2], double h) : 
			domain{ domain[0], domain[1] }, data{ data[0], data[1] }, h{ h } {}
		interval_element_type(double left, double right, double h = 1.) :
			domain{ left, right }, data{ left, right }, h{ h } {}
		interval_element_type(len_type length) : interval_element_type(0, length - 1) {}
		interval_element_type() : interval_element_type(0, 1) {}


		//! Return number of discrete elements in the interval.
		/*!
		 * Computes and returns the number of discrete elements in this interval
		 * in a standardized way. The computation assumes that an interval always
		 * contains its endpoints; for example, an interval between 0 and 2 with
		 * a spatial width of 1 means there are 3 discrete elements in the
		 * interval (elements 0, 1 and 2).
		 */
		len_type get_interval_count() const
		{
			return symphas::internal::range_count(data, width());
		}

		//! Return number of discrete elements in the domain interval.
		/*!
		 * Computes and returns the number of discrete elements in this domain
		 * in a standardized way. The computation assumes that an interval always
		 * contains its endpoints; for example, an interval between 0 and 2 with
		 * a spatial width of 1 means there are 3 discrete elements in the
		 * interval (elements 0, 1 and 2).
		 */
		len_type get_domain_count() const
		{
			return symphas::internal::range_count(domain, width());
		}

		//! By default, this will return the number of elements in the domain.
		len_type get_count() const
		{
			return get_domain_count();
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

		//! Returns the left endpoint of the domain.
		/*!
		 * The endpoint that is returned is always the physical end point with which
		 * the interval is initialized with.
		 */
		axis_coord_t domain_left() const
		{
			return domain[0];
		}

		//! Returns the right endpoint of the domain.
		/*!
		 * The endpoint that is returned is always the physical end point with which
		 * the interval is initialized with.
		 */
		axis_coord_t domain_right() const
		{
			return domain[1];
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
		double interval_length() const
		{
			return symphas::internal::range_length(data, width());
		}

		//! Return the length of the domain on the real axis.
		double domain_length() const
		{
			return symphas::internal::range_length(domain, width());
		}

		//! Return the length of the interval on the real axis.
		double length() const
		{
			return interval_length();
		}

		//! Sets the beginning and end values of the domain. 
		/*
		 * Chooses new values for the start and endpoints of the domain, keeping
		 * the existing interval within the domain.
		 *
		 * \param left The left endpoint of the interval.
		 * \param right The right endpoint of the interval.
		 * \param width The spatial distance between points in the interval.
		 */
		void set_domain(double left, double right, double width = 0)
		{
			symphas::internal::set_interval(domain, h, left, right, width);
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
		void set_interval(double left, double right)
		{
			symphas::internal::set_interval(data, left, right);
		}

		void set_interval_fraction(double left, double right)
		{
			double left0 = domain[0] + left * get_domain_count() * h;
			double right0 = left0 + (right - left) * get_domain_count() * h;
			symphas::internal::set_interval(data, left0, right0);
		}

		void set_interval_count(double left, len_type count)
		{
			symphas::internal::set_interval(data, left, (left + count * width()) / domain_length());
		}

		void set_interval_count(len_type count)
		{
			symphas::internal::set_count_from_r(data, h, count, 0.5);
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
		void set_domain_count(double left, double right, len_type count)
		{
			symphas::internal::set_count(domain, h, left, right, count);
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
			set_domain_count(left, right, count);
			set_interval_count(left, count);
		}

		//! Resize the domain with the given number of elements, with constant width.
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
		 * \param The point (from 0 to 1) along which the new interval is centered 
		 * around the old interval.
		 */
		void set_domain_count_from_r(len_type count, double width, double ratio)
		{
			symphas::internal::set_count_from_r(domain, h, count, width, ratio);
		}

		//! Resize the domain with the given number of elements, with constant width.
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
		 * \param The point (from 0 to 1) along which the new interval is centered 
		 * around the old interval.
		 */
		void set_domain_count_from_r(len_type count, double ratio)
		{
			symphas::internal::set_count_from_r(domain, h, count, ratio);
		}

		//! Resize the domain with the given number of elements, with constant width.
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
		 * \param The point (from 0 to 1) along which the new interval is centered 
		 * around the old interval.
		 */
		void set_count_from_r(len_type count, double width, double ratio)
		{
			set_domain_count_from_r(count, width, ratio);
		}

		//! Resize the domain with the given number of elements, with constant width.
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
			set_domain_count_from_r(count, ratio);
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
		void set_count(len_type count, double width)
		{
			set_count_from_r(width, count, 0.5);
		}

		void domain_to_interval()
		{
			domain[0] = data[0];
			domain[1] = data[1];
		}

		void interval_to_domain()
		{
			data[0] = domain[0];
			data[1] = domain[1];
		}
	};

	struct interval_data_type : std::map<Axis, interval_element_type>
	{
		using parent_type = std::map<Axis, interval_element_type>;
		using parent_type::parent_type;
		using parent_type::operator[];
		interval_data_type(size_t dim, interval_element_type const& interval) : parent_type()
		{
			for (iter_type i = 0; i < dim; ++i)
			{
				auto axis = symphas::index_to_axis(i);
				this->operator[](axis) = interval;
			}
		}
	};

	struct grid_info;
}

void swap(symphas::grid_info& first, symphas::grid_info& second);


namespace grid
{

	template<typename T>
	struct box_list : symphas::lib::array_container<T>
	{
		using parent_type = symphas::lib::array_container<T>;
		using parent_type::parent_type;
		using parent_type::data;
		using parent_type::n;
		using parent_type::operator T*;
		using parent_type::operator[];

		box_list(T dim0, T dim1, T dim2) :
			box_list(dim0, dim1, dim2, (dim1 > 0 && dim2 > 0) ? 3 : (dim1 > 0 || dim2 > 0) ? 2 : 1) {}
		box_list(T dim0, T dim1) :
			box_list(dim0, dim1, 0, (dim1 == 0) ? 1 : 2) {}
		box_list(T dim0) : box_list(dim0, 0, 0, 1) {}
		box_list() : box_list(0, 0, 0, 0) {}

		template<typename TT = T, std::enable_if_t<std::is_same<TT, len_type>::value, int> = 0>
		box_list(symphas::interval_data_type const& intervals) :
			box_list(
				(intervals.find(Axis::X) != intervals.end()) ? intervals.at(Axis::X).get_count() : 0,
				(intervals.find(Axis::Y) != intervals.end()) ? intervals.at(Axis::Y).get_count() : 0,
				(intervals.find(Axis::Z) != intervals.end()) ? intervals.at(Axis::Z).get_count() : 0) {}

		template<typename TT = T, std::enable_if_t<std::is_same<TT, double>::value, int> = 0>
		box_list(symphas::interval_data_type const& intervals) :
			box_list(
				(intervals.find(Axis::X) != intervals.end()) ? intervals.at(Axis::X).width() : 0,
				(intervals.find(Axis::Y) != intervals.end()) ? intervals.at(Axis::Y).width() : 0,
				(intervals.find(Axis::Z) != intervals.end()) ? intervals.at(Axis::Z).width() : 0) {}

		box_list(const T* data, size_t n) : parent_type(data, n) {}

		operator T() const
		{
			if (n > 0)
			{
				return data[0];
			}
			else
			{
				return 0;
			}
		}

		auto operator[] (iter_type i) const
		{
			if (i < n)
			{
				return data[i];
			}
			else
			{
				return 0;
			}
		}

		std::tuple<T, T, T> _3() const
		{
			return { operator[](0), operator[](1), operator[](2) };
		}

		std::tuple<T, T> _2() const
		{
			return { operator[](0), operator[](1) };
		}

		std::tuple<T> _1() const
		{
			return { operator[](0) };
		}


	protected:


		box_list(T dim0, T dim1, T dim2, size_t n) : parent_type(nullptr, n)
		{
			if (n > 0) data[0] = dim0;
			if (n > 1) data[1] = dim1;
			if (n > 2) data[2] = dim2;
		}

	};

	using dim_list = box_list<len_type>;
	using pos_list = symphas::lib::array_container<iter_type>;
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
				v[symphas::index_to_axis(i)] = symphas::interval_element_type(dims[i]);
			}

			return v;
		}
	}

public:
	
	symphas::interval_data_type intervals;				//!< Extent of the grid in the spatial axes.


	template<size_t... Is>
	grid_info(const double(&intervals)[sizeof...(Is)][2], double width, std::index_sequence<Is...>) :
		intervals{ { symphas::index_to_axis(Is), interval_element_type(intervals[Is][0], intervals[Is][1], width) }... } {}

	template<size_t D>
	grid_info(const double(&intervals)[D][2], double width = 1.) : grid_info(intervals, width, std::make_index_sequence<D>{}) {}

	template<size_t... Is>
	grid_info(const len_type(&intervals)[sizeof...(Is)][2], std::index_sequence<Is...>) :
		intervals{ { symphas::index_to_axis(Is), interval_element_type(intervals[Is][0], intervals[Is][1]) }... } {}

	template<size_t D>
	grid_info(const len_type(&intervals)[D][2]) : grid_info(intervals, std::make_index_sequence<D>{}) {}

	template<size_t... Is>
	grid_info(const double(&domain)[sizeof...(Is)][2], const double(&intervals)[sizeof...(Is)][2], double width, std::index_sequence<Is...>) :
		intervals{ { symphas::index_to_axis(Is), interval_element_type(domain[Is], intervals[Is], width)}...} {}

	template<size_t D>
	grid_info(const double(&domain)[D][2], const double(&intervals)[D][2], double width = 1.) : grid_info(domain, intervals, width, std::make_index_sequence<D>{}) {}

	template<size_t... Is>
	grid_info(const len_type(&dims)[sizeof...(Is)], const len_type(&intervals)[sizeof...(Is)][2], std::index_sequence<Is...>) :
		intervals{ { symphas::index_to_axis(Is), interval_element_type((double[2]) { 0, double(dims[Is] - 1) }, (double[2]) { double(intervals[Is][0]), double(intervals[Is][1] - 1) }, 1.) }... } {}

	template<size_t D>
	grid_info(const len_type(&dims)[D], const len_type(&intervals)[D][2]) : grid_info(dims, intervals, std::make_index_sequence<D>{}) {}

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

	grid_info(grid::dim_list const& dims) : grid_info(dims, dims.n) {}

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
		return intervals[ax];
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
	len_type num_points() const
	{
		len_type length = 1;
		for (iter_type i = 0; i < dimension(); ++i)
		{
			length *= at(symphas::index_to_axis(i)).get_count();
		}
		return length;
	}

	//! Gives the number of cells in the entire grid in discrete space.
	/*!
	 * Gives the number of cells in the entire grid in discrete space.
	 * This is computed using the lengths of all the intervals.
	 */
	len_type num_interval_points()
	{
		len_type length = 1;
		for (iter_type i = 0; i < dimension(); ++i)
		{
			length *= at(symphas::index_to_axis(i)).get_interval_count();
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
	grid::dim_list get_interval_dims() const
	{
		grid::dim_list dims(nullptr, dimension());
		for (iter_type i = 0; i < dimension(); ++i)
		{
			dims[i] = at(symphas::index_to_axis(i)).get_interval_count();
		}
		return dims;
	}

	//! Gives the get_stride of the grid.
	/*!
	 * Gives the number of cells along each dimension.
	 */
	grid::dim_list get_stride() const
	{
		auto dims = get_dims();
		grid::dim_list stride(nullptr, dimension());
		stride[0] = 1;
		for (iter_type i = 1; i < dimension(); ++i)
		{
			stride[i] = dims[i - 1] * stride[i - 1];
		}
		return stride;
	}

	//! Gives the get_stride of the grid.
	/*!
	 * Gives the number of cells along each dimension.
	 */
	grid::dim_list get_interval_stride() const
	{
		auto dims = get_interval_dims();
		grid::dim_list stride(nullptr, dimension());
		stride[0] = 1;
		for (iter_type i = 1; i < dimension(); ++i)
		{
			stride[i] = dims[i - 1] * stride[i - 1];
		}
		return stride;
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

	double area() const
	{
		double area = 1;
		for (auto length : get_lengths())
		{
			area *= length;
		}
		return area;
	}

	double element_area() const
	{
		double r = 1;
		for (auto width : get_widths())
		{
			r *= width;
		}
		return r;
	}

	//! Gives the first position of the interval.
	/*!
	 * Gives the left element count for each interval as a position list.
	 */
	grid::pos_list left() const
	{
		grid::pos_list pos(dimension());
		for (iter_type i = 1; i < dimension(); ++i)
		{
			auto interval = at(symphas::index_to_axis(i));
			pos[i] = (interval.left() - interval.domain_left()) / interval.width();
		}
		return pos;
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

namespace symphas
{
	inline bool is_valid(grid_info const& ginfo)
	{
		return (ginfo.dimension() > 0);
	}
}


//! \cond

/* definitions for accessing the given interval from the map named
 * intervals which has key type Axis and value type double[2]
 * the implementation of the value type is given below
 * this is used in the configuration and used as well in the writing utility
 */

#define INTERVAL_Xh_AT(i) intervals[i].at(Axis::X).width()
#define INTERVAL_Yh_AT(i) intervals[i].at(Axis::Y).width()
#define INTERVAL_Zh_AT(i) intervals[i].at(Axis::Z).width()
#define INTERVAL_Xh intervals.at(Axis::X).width()
#define INTERVAL_Yh intervals.at(Axis::Y).width()
#define INTERVAL_Zh intervals.at(Axis::Z).width()

#define DOMAIN_X0_AT(i) intervals[i].at(Axis::X).domain_left()
#define DOMAIN_Xn_AT(i) intervals[i].at(Axis::X).domain_right()
#define DOMAIN_Y0_AT(i) intervals[i].at(Axis::Y).domain_left()
#define DOMAIN_Yn_AT(i) intervals[i].at(Axis::Y).domain_right()
#define DOMAIN_Z0_AT(i) intervals[i].at(Axis::Z).domain_left()
#define DOMAIN_Zn_AT(i) intervals[i].at(Axis::Z).domain_right()
#define DOMAIN_X0 intervals.at(Axis::X).domain_left()
#define DOMAIN_Xn intervals.at(Axis::X).domain_right()
#define DOMAIN_Y0 intervals.at(Axis::Y).domain_left()
#define DOMAIN_Yn intervals.at(Axis::Y).domain_right()
#define DOMAIN_Z0 intervals.at(Axis::Z).domain_left()
#define DOMAIN_Zn intervals.at(Axis::Z).domain_right()

#define DOMAIN_Xc_AT(i) intervals[i].at(Axis::X).get_count()
#define DOMAIN_Yc_AT(i) intervals[i].at(Axis::Y).get_count()
#define DOMAIN_Zc_AT(i) intervals[i].at(Axis::Z).get_count()
#define DOMAIN_Xc intervals.at(Axis::X).get_count()
#define DOMAIN_Yc intervals.at(Axis::Y).get_count()
#define DOMAIN_Zc intervals.at(Axis::Z).get_count()

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

#define INTERVAL_Xc_AT(i) intervals[i].at(Axis::X).get_interval_count()
#define INTERVAL_Yc_AT(i) intervals[i].at(Axis::Y).get_interval_count()
#define INTERVAL_Zc_AT(i) intervals[i].at(Axis::Z).get_interval_count()
#define INTERVAL_Xc intervals.at(Axis::X).get_interval_count()
#define INTERVAL_Yc intervals.at(Axis::Y).get_interval_count()
#define INTERVAL_Zc intervals.at(Axis::Z).get_interval_count()


//! \endcond

