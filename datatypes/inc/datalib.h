
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
 * PURPOSE: Defines objects for managing data representing independent
 * and dependent data. This also applies to managing two arbitrary lists of
 * data, where one is independent and the other is dependent, which can be
 * considered x and y data on a Cartesian grid. Extensions are also added
 * where these objects also allow managing up to 3-dimensional phase field
 * data. Also adds functions for computing properties about the data, such
 * as the radial average.
 *
 * ***************************************************************************
 */


#pragma once

#include <memory>

#include "grid.h"
#include "spslibfftw.h"

namespace symphas
{
	namespace lib
	{
		//! Returns a list of each axis position using the given intervals.
		/*!
		 * Returns a list of each axis position using the given intervals.
		 *
		 * \param intervals Intervals used to construct the list of grid
		 * axis positions.
		 */
		template<size_t D>
		auto new_system_axis_list(symphas::interval_data_type const& intervals);

		//! Returns a list of each axis position using the given dimensions.
		/*!
		 * Returns a list of each axis position using the given dimensions.
		 *
		 * \param dims Dimensions of the system for which the axis list
		 * is constructed.
		 */
		template<size_t D>
		auto new_system_axis_list(const len_type* dims);

		template<size_t D>
		grid::dim_list get_dimensions(const axis_nd_t<D>* axis_values, len_type len);

		template<size_t D, typename T>
		grid::dim_list get_dimensions(std::vector<std::pair<axis_nd_t<D>, T>> const& data);
	}



	namespace dft
	{
		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 3-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<3>* dft_x, iter_type L, iter_type M, iter_type N);

		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 2-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<2>* dft_x, iter_type L, iter_type M);

		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 1-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<1>* dft_x, iter_type L);

		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 3-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<3>* dft_x, const axis_nd_t<3>* data_x, len_type len);

		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 2-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<2>* dft_x, const axis_nd_t<2>* data_x, len_type len);

		//! Fills the array with axis values of Fourier space. 
		/*!
		 * Compute \f$x\f$-axis values of a 1-dimensional Fourier transform, and
		 * store them in the provided array.
		 */
		void fill_x_axis(axis_nd_t<1>* dft_x, const axis_nd_t<1>* data_x, len_type len);
	}



	//! Representation of a data series with dependent and independent axes.
	/*!
	 * A container for data in order to allow manipulation of that data.
	 * The basic field contains a pair of the give _X_ and _Y_ data. The data, once
	 * passed to the constructor, is owned by this object. This does not apply
	 * to data which are arrays.
	 * 
	 * \tparam X The type of the independent data.
	 * \tparam Y The type of the dependent data.
	 */
	template<typename X, typename Y>
	struct Field
	{

		//! Create a new field with the data points.
		/*!
		 * The first data point is an independent data point, and the second
		 * is dependent data. These are copied into the object.
		 * 
		 * \param data_x The independent data.
		 * \param data_y The dependent data.
		 * \param len The length of the underlying values of the given data.
		 */
		Field(X data_x, Y data_y, len_type len = 1) : data{ data_x, data_y }, len{ len } {}
		Field(len_type len = 1) : data{ X{}, Y{} }, len{ len } {}


		//! Get the reference to the independent data.
		const X& data_x() const
		{
			return data.first;
		}

		//! Get the reference to the dependent data.
		const Y& data_y() const
		{
			return data.second;
		}

		//! Get the reference to the independent data.
		X& data_x()
		{
			return data.first;
		}

		//! Get the reference to the dependent data.
		Y& data_y()
		{
			return data.second;
		}

		//! Get the length of the dependent data array.
		len_type length() const
		{
			return len;
		}

	protected:

		std::pair<X, Y> data;
		len_type len;
	};



	//! Representation of data arrays with dependent and independent axes.
	/*!
	 * The given dependent data is an array of values. 
	 * 
	 * Specialization of symphas::Field to take into account pointer memory management. 
	 * The symphas::Field will move the given independent data into the shared pointer, then the
	 * field can be manipulated at will with no special consideration for memory
	 * management.
	 * 
	 * This means that if a raw pointer is given, it must 
	 * not be deleted.
	 * 
	 * \tparam X The type of the independent data.
	 * \tparam Y The type of the dependent data.
	 */
	template<typename X, typename Y>
	struct Field<X, Y*>
	{
		//! Create a new field out of the data. 
		/*!
		 * The first list is the independent data point, and the second list is
		 * the dependent data of length \p len. Deleting the pointer which
		 * is given will corrupt the data in this symphas::Field as well. 
		 * The original array becomes owned by the shared
		 * pointer, preventing any memory leak.
		 *
		 * \param data_x The independent data.
		 * \param data_y The dependent data.
		 * \param len The length of the data.
		 */
		Field(X data_x, Y* data_y, len_type len) : 
			data{ data_x, data_y }, len{ len }, y{ data.second.get() } {}

		//! Create a new field out of the data. 
		/*!
		 * The first list is the independent data point, and the second list is 
		 * the dependent data of length \p len.
		 *
		 * \param data_x The independent data.
		 * \param data_y The dependent data.
		 * \param len The length of the data.
		 */
		Field(X data_x, std::shared_ptr<Y[]> data_y, len_type len) : 
			data{ data_x, data_y }, len{ len }, y{ data.second.get() } {}

		//! Create a new field out of the block data. 
		/*!
		 * The first list is the independent data, and the second list is the
		 * dependent data. **The data is copied from the block**, so there is no
		 * longer a data dependency.
		 *
		 * \param data_x The independent data.
		 * \param block The dependent data.
		 */
		Field(X data_x, Block<Y> const& block) : 
			data{ data_x, new Y[block.len] }, len{ block.len }, y{ data.second.get() }
		{
			for (iter_type i = 0; i < block.len; ++i)
			{
				data.second[i] = block.values[i];
			}
		}

		Field(len_type len = 1) : data{ X{}, new Y[len]{} }, len{ len } {}


		Field(Field<X, Y*> const& other) : Field(other.data.first, other.data.second, other.len) {}
		Field(Field<X, Y*>&& other) : Field(other.data.first, other.data.second, other.len) {}
		Field<X, Y*>& operator=(Field<X, Y*> other)
		{
			std::swap(data, other.data);
			std::swap(len, other.len);
			std::swap(y, other.y);
			return *this;
		}


		//! Get the reference to the independent data.
		X& data_x()
		{
			return data.first;
		}

		//! Get the reference to the independent data.
		const X& data_x() const
		{
			return data.first;
		}

		//! Get the reference to the dependent data array.
		std::shared_ptr<Y[]>& data_y()
		{
			return data.second;
		}

		//! Get the reference to the dependent data array.
		std::shared_ptr<Y[]> const& data_y() const
		{
			return data.second;
		}

		//! Get the reference to the dependent data array.
		Y* field_y()
		{
			return data.second.get();
		}

		//! Get the reference to the dependent data array.
		const Y* field_y() const
		{
			return data.second.get();
		}

		//! Get the length of the dependent data array.
		len_type length() const
		{
			return len;
		}

		//! Access an element in the dependent data array.
		Y& operator[](iter_type i)
		{
			return data_y()[i];
		}

		//! Access an element in the dependent data array.
		const Y& operator[](iter_type i) const
		{
			return data_y()[i];
		}

	protected:

		std::pair<X, std::shared_ptr<Y[]>> data;
		len_type len;

	public:

		Y* y;

	};

	template<size_t D, typename Y>
	struct FieldAxis;

	//! Representation of data arrays with dependent and independent axes.
	/*!
	 * A list of arrays representing independent and dependent data. The data
	 * is stored in a shared pointer to avoid data duplication, so a copy
	 * of the field will manage the same data.
	 * 
	 * Specialization of symphas::Field to take into account pointer memory management. 
	 * The symphas::Field will move the given arrays into the shared pointers, then the
	 * field can be manipulated at will with no special consideration for memory
	 * management. 
	 * 
	 * This means that if a raw pointer is given, it must 
	 * not be deleted.
	 * 
	 * Also, important point: If copying only one data array from one
	 * field into another one, then you have to use the function data_x() or
	 * data_y(), because this will take the shared pointer. There are members
	 * to simplify memory access to the underlying pointers, but these 
	 * should never be used in initializing another Field.
	 *
	 * In order to create a different symphas::Field for the same data without
	 * data dependence, new arrays will have to be created.
	 * 
	 * \tparam X The type of the independent data.
	 * \tparam Y The type of the dependent data.
	 */
	template<typename X, typename Y>
	struct Field<X*, Y*>
	{
		//! Create a new field out of the lists of data. 
		/*!
		 * The first list is the independent data, and the second list is the
		 * dependent data. They should be the same length \p len. Since the
		 * arrays which are given are raw pointers, deallocating them
		 * would corrupt the data. The original arrays become owned by the shared
		 * pointers, preventing any memory leak.
		 *
		 * \param data_x The independent data.
		 * \param data_y The dependent data.
		 * \param len The length of the data.
		 */
		Field(X* data_x, Y* data_y, len_type len) :
			data{ std::shared_ptr<X[]>(data_x), std::shared_ptr<Y[]>(data_y) }, len{ len },
			x{ data.first.get() }, y{ data.second.get() } {}

		//! Create a new field out of the lists of data. 
		/*!
		 * The first list is the independent data, and the second list is the
		 * dependent data. They should be the same length \p len.
		 *
		 * \param data_x The independent data.
		 * \param data_y The dependent data.
		 * \param len The length of the data.
		 */
		Field(std::shared_ptr<X[]> data_x, std::shared_ptr<Y[]> data_y, len_type len) :
			data{ data_x, data_y }, len{ len },
			x{ data.first.get() }, y{ data.second.get() } {}

		Field(X* data_x, std::shared_ptr<Y[]> data_y, len_type len) :
			data{ data_x, data_y }, len{ len },
			x{ data.first.get() }, y{ data.second.get() } {}
		Field(std::shared_ptr<X[]> data_x, Y* data_y, len_type len) :
			data{ data_x, data_y }, len{ len },
			x{ data.first.get() }, y{ data.second.get() } {}

		Field(len_type len = 1) : Field(new X[len]{}, new Y[len]{}, len) {}

		//! Create a new field out of the block data. 
		/*!
		 * The first list is the independent data, and the second list is the
		 * dependent data. **The data is copied from the block**, so there is no
		 * longer a data dependency.
		 *
		 * \param data_x The independent data.
		 * \param block The dependent data.
		 */
		Field(X* data_x, Block<Y> const& block) : Field(block.len)
		{
			for (iter_type i = 0; i < block.len; ++i)
			{
				data.second[i] = block.values[i];
			}
		}

		//! Create a new field out of the block data. 
		/*!
		 * The first list is the independent data, and the second list is the
		 * dependent data. **The data is copied from the block**, so there is no
		 * longer a data dependency.
		 *
		 * \param data_x The independent data.
		 * \param block The dependent data.
		 */
		template<size_t D>
		Field(Grid<Y, D> const& grid) : 
			Field(symphas::lib::new_system_axis_list<D>(grid.dims), new Y[grid.len], grid.len)
		{
			for (iter_type i = 0; i < grid.len; ++i)
			{
				data.second[i] = grid.values[i];
			}
		}


		Field(Field<X*, Y*> const& other) : Field(other.data.first, other.data.second, other.len) {}
		Field(Field<X*, Y*>&& other) : Field(other.data.first, other.data.second, other.len) {}
		Field<X*, Y*>& operator=(Field<X*, Y*> other)
		{
			std::swap(data, other.data);
			std::swap(len, other.len);
			std::swap(x, other.x);
			std::swap(y, other.y);
			return *this;
		}

		//! Get the reference to the independent data array.
		std::shared_ptr<X[]>& data_x()
		{
			return data.first;
		}

		//! Get the reference to the dependent data array.
		std::shared_ptr<Y[]>& data_y()
		{
			return data.second;
		}

		//! Get the reference to the dependent data array.
		std::shared_ptr<X[]> const& data_x() const
		{
			return data.first;
		}

		//! Get the reference to the independent data array.
		std::shared_ptr<Y[]> const& data_y() const
		{
			return data.second;
		}

		//! Get the reference to the dependent data array.
		X* field_x()
		{
			return data.first.get();
		}

		//! Get the reference to the dependent data array.
		const X* field_x() const
		{
			return data.first.get();
		}

		//! Get the reference to the dependent data array.
		Y* field_y()
		{
			return data.second.get();
		}

		//! Get the reference to the dependent data array.
		const Y* field_y() const
		{
			return data.second.get();
		}

		//! Get the length of the dependent data array.
		len_type length() const
		{
			return len;
		}

		//! Access an element in the dependent data array.
		Y& operator[](iter_type i)
		{
			return data_y()[i];
		}

		//! Access an element in the dependent data array.
		const Y& operator[](iter_type i) const
		{
			return data_y()[i];
		}

	protected:

		std::pair<std::shared_ptr<X[]>, std::shared_ptr<Y[]>> data;
		len_type len;

	public:

		X* x;
		Y* y;
	};



	template<typename X, typename Y>
	Field(X, Y*, len_type)->Field<X, Y*>;
	template<typename X, typename Y>
	Field(X, std::shared_ptr<Y[]>, len_type)->Field<X, Y*>;

	template<typename Y, size_t D>
	Field(Grid<Y, D>)->Field<axis_nd_t<D>*, Y*>;

	template<typename X, typename Y>
	Field(X*, Y*, len_type)->Field<X*, Y*>;
	template<typename X, typename Y>
	Field(std::shared_ptr<X[]>, std::shared_ptr<Y[]>, len_type)->Field<X*, Y*>;


	//! Representation of data arrays with standardized dependent axis.
	/*!
	 * A derived field in which the \f$x\f$ values are axis types.
	 * Therefore this object is parameterized on the dimension.
	 */
	template<size_t D, typename Y>
	struct FieldAxis : Field<axis_nd_t<D>*, Y>
	{
		using parent_type = Field<axis_nd_t<D>*, Y>;
		using parent_type::parent_type;
		using parent_type::x;
		using parent_type::len;

		auto& sort(const len_type*);
		auto& sort();
	};


	namespace internal
	{
		template<typename Y>
		auto& sort_data(FieldAxis<1, Y*> &data, const len_type* dims)
		{
			symphas::lib::sort_data(
				data.x, data.y, dims[0], 
				symphas::lib::regular_sort);
			return data;
		}

		template<typename Y>
		auto& sort_data(FieldAxis<2, Y*> &data, const len_type* dims)
		{
			symphas::lib::sort_data(
				data.x, data.y, dims[0], dims[1], 
				symphas::lib::regular_sort);
			return data;
		}

		template<typename Y>
		auto& sort_data(FieldAxis<3, Y*> &data, const len_type* dims)
		{
			symphas::lib::sort_data(
				data.x, data.y, dims[0], dims[1], dims[1], 
				symphas::lib::regular_sort);
			return data;
		}
	}


	template<size_t D, typename Y>
	auto& FieldAxis<D, Y>::sort(const len_type *dims)
	{
		internal::sort_data(*this, dims);
		return *this;
	}

	template<size_t D, typename Y>
	auto& FieldAxis<D, Y>::sort()
	{
		len_type dims[D]{ 0 };
		if constexpr (D == 1)
		{
			auto L = symphas::lib::get_dimensions<1>(x, len);
			dims[0] = L;
		}
		else if constexpr (D == 2)
		{
			auto [L, M] = symphas::lib::get_dimensions<2>(x, len)._2();
			dims[0] = L;
			dims[1] = M;
		}
		else if constexpr (D == 3)
		{
			auto [L, M, N] = symphas::lib::get_dimensions<3>(x, len)._3();
			dims[0] = L;
			dims[1] = M;
			dims[2] = N;
		}
		internal::sort_data(*this, dims);
		return *this;
	}


	//! Compatibility for 1-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_1d_type*, Y*, len_type)->FieldAxis<1, Y*>;

	//! Compatibility for 2-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_2d_type*, Y*, len_type)->FieldAxis<2, Y*>;

	//! Compatibility for 3-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_3d_type*, Y*, len_type)->FieldAxis<3, Y*>;

	//! Compatibility for 1-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_1d_type*, Block<Y> const&)->FieldAxis<1, Y*>;

	//! Compatibility for 2-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_2d_type*, Block<Y> const&)->FieldAxis<2, Y*>;

	//! Compatibility for 3-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_3d_type*, Block<Y> const&)->FieldAxis<3, Y*>;

	//! Compatibility for 1-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_1d_type*, Grid<Y, 1> const&)->FieldAxis<1, Y*>;

	//! Compatibility for 2-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_2d_type*, Grid<Y, 2> const&)->FieldAxis<2, Y*>;

	//! Compatibility for 3-dimensional phase field data type.
	template<typename Y>
	FieldAxis(axis_3d_type*, Grid<Y, 3> const&)->FieldAxis<3, Y*>;

	template<typename Y>
	FieldAxis(Field<axis_1d_type*, Y*>)->FieldAxis<1, Y*>;
	template<typename Y>
	FieldAxis(Field<axis_2d_type*, Y*>)->FieldAxis<2, Y*>;
	template<typename Y>
	FieldAxis(Field<axis_3d_type*, Y*>)->FieldAxis<3, Y*>;


	template<typename Y, size_t D>
	FieldAxis(Grid<Y, D>)->FieldAxis<D, Y*>;

	template<typename F>
	struct field_data_type;
	

	template<typename X, typename Y>
	struct field_data_type<Field<X, Y>>
	{
		using x = X;
		using y = Y;
	};

	template<size_t D, typename Y>
	struct field_data_type<FieldAxis<D, Y>>
	{
		using x = typename field_data_type<typename FieldAxis<D, Y>::parent_type>::x;
		using y = typename field_data_type<typename FieldAxis<D, Y>::parent_type>::y;
	};

	//! Obtain the data type of the x-axis of the Field.
	/*!
	 * Obtains the data type of the x-axis of the given Field type. This is always the
	 * direct type, which in the case of a Field of an array in `X` (e.g. `Field<X*, _>`)
	 * would be a pointer type. This is likewise for FieldAxis, which inherits from
	 * `Field<axis_nd_t<D>, _>`.
	 */
	template<typename F>
	using field_x_t = typename field_data_type<F>::x;

	//! Obtain the data type of the y-axis of the Field.
	/*!
	 * Obtains the data type of the y-axis of the given Field type. This is always the
	 * direct type, which in the case of a Field of an array in `Y` (e.g. `Field<_, Y*>`)
	 * would be a pointer type. 
	 */
	template<typename F>
	using field_y_t = typename field_data_type<F>::y;

}

// ****************************************************************************************


namespace symphas::lib
{

	//! Get the Fourier transform of the given 1-dimensional data.
	/*!
	 * The primary FT computation routine, combining the algorithms and
	 * using the FFTW routine. Places the newly computed arrays into 
	 * a field.
	 * 
	 * \param data_x The \f$x\f$-axis of the data to be transformed.
	 * \param data_y The data which is transformed.
	 * \param len The length of the data series that is transformed.
	 */
	template<typename Y,
		typename = decltype(symphas::dft::dft(std::declval<Y*>(), std::declval<complex_t*>(), std::declval<len_type>(), std::declval<bool>()))>
	auto fourier_transform(const axis_1d_type* data_x, const Y* data_y, const len_type len, bool backward = false)
	{
		axis_1d_type* dft_x = new axis_1d_type[len];
		complex_t* dft_y = new complex_t[len];

		symphas::dft::dft<1>(data_y, dft_y, { len }, backward);
		symphas::dft::fill_x_axis(dft_x, data_x, len);

		return Field<axis_1d_type*, complex_t*>{ dft_x, dft_y, len };

	}



	//! Get the Fourier transform of the given 2-dimensional data.
	/*!
	 * The primary FT computation routine, combining the algorithms and
	 * using the FFTW routine. Places the newly computed arrays into 
	 * a field.
	 *
	 * \param data_x The \f$x\f$-axis of the data to be transformed.
	 * \param data_y The data which is transformed.
	 * \param len The length of the data series that is transformed.
	 */
	template<typename Y,
		typename = decltype(symphas::dft::dft(std::declval<Y*>(), std::declval<complex_t*>(), std::declval<len_type>(), std::declval<len_type>(), std::declval<bool>()))>
	auto fourier_transform(const axis_2d_type* data_x, const Y* data_y, const len_type len, bool backward = false)
	{
		axis_2d_type* dft_x = new axis_2d_type[len];
		complex_t* dft_y = new complex_t[len];

		auto [L, M] = symphas::lib::get_dimensions<2>(data_x, len)._2();
		symphas::dft::dft(data_y, dft_y, L, M, backward);
		symphas::dft::fill_x_axis(dft_x, data_x, len);

		return FieldAxis<2, complex_t*>{ dft_x, dft_y, len };
	}

	//! Get the Fourier transform of the given 3-dimensional data.
	/*!
	 * The primary FT computation routine, combining the algorithms and
	 * using the FFTW routine. Places the newly computed arrays into 
	 * a field.
	 *
	 * \param data_x The \f$x\f$-axis of the data to be transformed.
	 * \param data_y The data which is transformed.
	 * \param len The length of the data series that is transformed.
	 */
	template<typename Y,
		typename = decltype(symphas::dft::dft(std::declval<Y*>(), std::declval<complex_t*>(), std::declval<len_type>(), std::declval<len_type>(), std::declval<len_type>(), std::declval<bool>()))>
	auto fourier_transform(const axis_3d_type* data_x, const Y* data_y, const len_type len, bool backward = false)
	{
		axis_3d_type* dft_x = new axis_3d_type[len];
		complex_t* dft_y = new complex_t[len];

		auto [L, M, N] = symphas::lib::get_dimensions<3>(data_x, len)._3();
		symphas::dft::dft(data_y, dft_y, L, M, N, backward);
		symphas::dft::fill_x_axis(dft_x, data_x, len);

		return FieldAxis<3, complex_t*>{ dft_x, dft_y, len };
	}


	//! Get the Fourier transform of the given 3-dimensional data.
	/*!
	 * The primary FT computation routine, combining the algorithms and
	 * using the FFTW routine. Places the newly computed arrays into
	 * a field.
	 *
	 * \param data_x The \f$x\f$-axis of the data to be transformed.
	 * \param data_y The data which is transformed.
	 * \param len The length of the data series that is transformed.
	 */
	template<typename Y,
		typename = decltype(symphas::dft::dft(std::declval<Y*>(), std::declval<complex_t*>(), std::declval<len_type>(), std::declval<bool>()))>
	auto fourier_transform(const iter_type* data_x, const Y* data_y, const len_type len, bool backward = false)
	{
		complex_t* dft_y = new complex_t[len];
		iter_type* dft_x = new iter_type[len];

		symphas::dft::dft(data_y, dft_y, len, backward);
		std::copy(data_x, data_x + len, dft_x);

		return Field<iter_type*, complex_t*>{ dft_x, dft_y, len };
	}

	//! Get the Fourier transform of the given 3-dimensional data.
	/*!
	 * The primary FT computation routine, combining the algorithms and
	 * using the FFTW routine. Places the newly computed arrays into
	 * a field.
	 *
	 * \param data_x The \f$x\f$-axis of the data to be transformed.
	 * \param data_y The data which is transformed.
	 * \param len The length of the data series that is transformed.
	 */
	template<typename X, typename Y, typename std::enable_if<!std::is_pointer<X>::value, int>::type = 0, 
		typename = decltype(symphas::dft::dft(std::declval<Y*>(), std::declval<complex_t*>(), std::declval<len_type>(), std::declval<bool>()))>
	auto fourier_transform(const X data_x, const Y* data_y, const len_type len, bool backward = false) 
	{
		complex_t* dft_y = new complex_t[len];
		symphas::dft::dft(data_y, dft_y, len, backward);
		return Field<X, complex_t*>{ data_x, dft_y, len };
	}

	// ****************************************************************************************

	//! Computes a new Fourier transformed field.
	/*!
	 * The given data field is Fourier transformed and a new field is returned.
	 * This is primarily an overload of symphas::lib::fourier_transform(const X, const Y*, const len_type).
	 * 
	 * \param data The data which is transformed.
	 */
	template<typename X, typename Y>
	auto fourier_transform(Field<X, Y*> const& data, bool backward = false)
	{
		return fourier_transform(data.data_x(), data.field_y(), data.length(), backward);
	}

	//! Computes a new Fourier transformed field.
	/*!
	 * The given data field is Fourier transformed and a new field is returned.
	 * This is primarily an overload of symphas::lib::fourier_transform(const X, const Y*, const len_type).
	 *
	 * \param data The data which is transformed.
	 */
	template<typename X, typename Y>
	auto fourier_transform(Field<X*, Y*> const& data, bool backward = false)
	{
		return fourier_transform(data.field_x(), data.field_y(), data.length(), backward);
	}
}



//! Compound addition assignment operator.
/*!
 * Compound assignment operator for addition, performing point-wise
 * addition of all elements returned by symphas::Field<X, Y>::data_y().
 *
 * \param rhs The symphas::Field which is added pointwise to the symphas::Field on the
 * left hand side.
 */
template<typename X, typename Y1, typename Y2>
symphas::Field<X, Y1*>& operator+=(symphas::Field<X, Y1*>& lhs, symphas::Field<X, Y2*>& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());
	for (int i = 0; i < len; ++i)
	{
		lhs[i] += rhs[i];
	}
	return lhs;
}

//! Addition operator.
/*!
 * Addition operator, which returns a new symphas::Field by adding the left
 * hand side symphas::Field to the right hand side symphas::Field.
 *
 * \param lhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 * \param rhs The symphas::Field<X, Y> instance on the right hand side of the
 * addition operator.
 */
template<typename X, typename Y1, typename Y2>
auto operator+(symphas::Field<X*, Y1*> const& lhs, symphas::Field<X*, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());

	using Y = decltype(std::declval<Y1>() + std::declval<Y2>());
	symphas::Field<X*, Y*> res{ len };
	for (int i = 0; i < len; ++i)
	{
		res.x[i] = lhs.x[i];
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}

//! Addition operator.
/*!
 * Addition operator, which returns a new symphas::Field by adding the left
 * hand side symphas::Field to the right hand side symphas::Field.
 *
 * \param lhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 * \param rhs The symphas::Field<X, Y> instance on the right hand side of the
 * addition operator.
 */
template<typename X, typename Y1, typename Y2>
auto operator+(symphas::Field<X, Y1*> const& lhs, symphas::Field<X, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());

	using Y = decltype(std::declval<Y1>() + std::declval<Y2>());
	symphas::Field<X, Y*> res{ len };
	for (int i = 0; i < len; ++i)
	{
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}


//! Compound subtraction assignment operator.
/*!
 * Compound assignment operator for addition, performing point-wise
 * subtraction of all elements returned by symphas::Field<X, Y>::data_y().
 *
 * \param rhs The symphas::Field which subtracts point-wise from the to the symphas::Field
 * on the left hand side.
 */
template<typename X, typename Y1, typename Y2>
symphas::Field<X, Y1*>& operator-=(symphas::Field<X, Y1*>& lhs, symphas::Field<X, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());
	for (int i = 0; i < len; ++i)
	{
		lhs[i] -= rhs[i];
	}
	return lhs;
}

//! Subtraction operator.
/*!
 * Subtraction operator, which returns a new symphas::Field by subtracting the
 * right hand side symphas::Field from the left hand side symphas::Field.
 *
 * \param lhs The symphas::Field instance on the left hand side of the
 * addition operator.
 * \param rhs The symphas::Field instance on the right hand side of the
 * addition operator.
 */
template<typename X, typename Y1, typename Y2>
auto operator-(symphas::Field<X, Y1*> const& lhs, symphas::Field<X, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());

	using Y = decltype(std::declval<Y1>() - std::declval<Y2>());
	symphas::Field<X, Y*> res{ len };
	for (int i = 0; i < len; ++i)
	{
		res.x[i] = lhs.x[i];
		res[i] = lhs[i] - rhs[i];
	}
	return res;
}

//! Unary negative operator.
/*!
 * Operator applied to the a symphas::Field, returning a new symphas::Field with all
 * values with the unary negative operator applied to them.
 */
template<typename X, typename Y>
auto operator-(symphas::Field<X, Y*> const& lhs)
{
	symphas::Field<X, decltype(-std::declval<Y>())*> neg{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		neg.x[i] = lhs.x[i];
		neg[i] = -lhs[i];
	}
	return neg;
}

//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param lhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 * \param scalar The scalar to multiply all the values.
 */
template<typename X, typename Y>
auto operator*(symphas::Field<X, Y*> const& lhs, Y const& scalar)
{
	symphas::Field<X, Y*> scaled{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		scaled.x[i] = lhs.x[i];
		scaled[i] *= scalar;
	}
	return scaled;
}

//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param lhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 * \param scalar The scalar to multiply all the values.
 */
template<typename X, typename Y>
auto operator*(symphas::Field<X, Y*> const& lhs, scalar_t scalar)
{
	symphas::Field<X, mul_result_t<scalar_t, Y>*> scaled{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		scaled.x[i] = lhs.x[i];
		scaled[i] = scalar * lhs[i];
	}
	return scaled;
}


//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param scalar The scalar to multiply all the values.
 * \param rhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 */
template<typename T, typename X, typename Y>
auto operator*(T const& scalar, symphas::Field<X, Y*> const& rhs)
{
	return rhs * scalar;
}

//! Scalar multiplication compound assignment for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param scalar The scalar to multiply all the values.
 * \param rhs The symphas::Field<X, Y> instance on the left hand side of the
 * addition operator.
 */
template<typename X, typename Y, typename T>
auto operator*=(symphas::Field<X, Y*>& lhs, T scalar)
{
	for (int i = 0; i < lhs.length(); ++i)
	{
		lhs[i] *= scalar;
	}
	return lhs;
}




//! Compound addition assignment operator.
/*!
 * Compound assignment operator for addition, performing point-wise
 * addition of all elements returned by symphas::FieldAxis<D, Y>::data_y().
 *
 * \param rhs The symphas::Field which is added pointwise to the symphas::Field on the
 * left hand side.
 */
template<size_t D, typename Y1, typename Y2>
symphas::FieldAxis<D, Y1*>& operator+=(symphas::FieldAxis<D, Y1*>& lhs, symphas::FieldAxis<D, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());
	for (int i = 0; i < len; ++i)
	{
		lhs[i] += rhs[i];
	}
	return lhs;
}

//! Addition operator.
/*!
 * Addition operator, which returns a new symphas::Field by adding the left
 * hand side symphas::Field to the right hand side symphas::Field.
 *
 * \param lhs The symphas::FieldAxis<D, Y> instance on the left hand side of the
 * addition operator.
 * \param rhs The symphas::FieldAxis<D, Y> instance on the right hand side of the
 * addition operator.
 */
template<size_t D, typename Y1, typename Y2>
auto operator+(symphas::FieldAxis<D, Y1*> const& lhs, symphas::FieldAxis<D, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());

	using Y = decltype(std::declval<Y1>() + std::declval<Y2>());
	symphas::FieldAxis<D, Y*> res{ len };
	for (int i = 0; i < len; ++i)
	{
		res.x[i] = lhs.x[i];
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}

//! Compound subtraction assignment operator.
/*!
 * Compound assignment operator for addition, performing point-wise
 * subtraction of all elements returned by symphas::FieldAxis<D, Y>::data_y().
 *
 * \param rhs The symphas::Field which subtracts point-wise from the to the symphas::Field
 * on the left hand side.
 */
template<size_t D, typename Y1, typename Y2>
symphas::FieldAxis<D, Y1*>& operator-=(symphas::FieldAxis<D, Y1*>& lhs, const symphas::FieldAxis<D, Y2*>& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());
	for (int i = 0; i < len; ++i)
	{
		lhs[i] -= rhs[i];
	}
	return lhs;
}

//! Subtraction operator.
/*!
 * Subtraction operator, which returns a new symphas::Field by subtracting the
 * right hand side symphas::Field from the left hand side symphas::Field.
 *
 * \param lhs The symphas::Field instance on the left hand side of the
 * addition operator.
 * \param rhs The symphas::Field instance on the right hand side of the
 * addition operator.
 */
template<size_t D, typename Y1, typename Y2>
auto operator-(symphas::FieldAxis<D, Y1*> const& lhs, symphas::FieldAxis<D, Y2*> const& rhs)
{
	len_type len = std::min(lhs.length(), rhs.length());

	using Y = decltype(std::declval<Y1>() - std::declval<Y2>());
	symphas::FieldAxis<D, Y*> res{ len };
	for (int i = 0; i < len; ++i)
	{
		res.x[i] = lhs.x[i];
		res[i] = lhs[i] - rhs[i];
	}
	return res;
}

//! Unary negative operator.
/*!
 * Operator applied to the a symphas::Field, returning a new symphas::Field with all
 * values with the unary negative operator applied to them.
 */
template<size_t D, typename Y>
auto operator-(symphas::FieldAxis<D, Y*> const& lhs)
{
	symphas::FieldAxis<D, decltype(-std::declval<Y>())*> neg{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		neg.x[i] = lhs.x[i];
		neg[i] = -lhs[i];
	}
	return neg;
}

//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param lhs The symphas::FieldAxis<D, Y> instance on the left hand side of the
 * addition operator.
 * \param scalar The scalar to multiply all the values.
 */
template<size_t D, typename Y>
auto operator*(symphas::FieldAxis<D, Y*> const& lhs, Y const& scalar)
{
	symphas::FieldAxis<D, Y*> scaled{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		scaled.x[i] = lhs.x[i];
		scaled[i] *= scalar;
	}
	return scaled;
}

//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param lhs The symphas::FieldAxis<D, Y> instance on the left hand side of the
 * addition operator.
 * \param scalar The scalar to multiply all the values.
 */
template<size_t D, typename Y>
auto operator*(symphas::FieldAxis<D, Y*> const& lhs, scalar_t scalar)
{
	symphas::FieldAxis<D, mul_result_t<scalar_t, Y>*> scaled{ lhs.length() };
	for (int i = 0; i < lhs.length(); ++i)
	{
		scaled.x[i] = lhs.x[i];
		scaled[i] = scalar * lhs[i];
	}
	return scaled;
}


//! Scalar multiplication for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param scalar The scalar to multiply all the values.
 * \param rhs The symphas::FieldAxis<D, Y> instance on the left hand side of the
 * addition operator.
 */
template<typename T, size_t D, typename Y>
auto operator*(T const& scalar, symphas::FieldAxis<D, Y*> const& rhs)
{
	return rhs * scalar;
}

//! Scalar multiplication compound assignment for a symphas::Field.
/*!
 * Scalar multiplication, which multiplies all values in the field by
 * a given value.
 *
 * \param scalar The scalar to multiply all the values.
 * \param rhs The symphas::FieldAxis<D, Y> instance on the left hand side of the
 * addition operator.
 */
template<size_t D, typename Y, typename T>
auto operator*=(symphas::FieldAxis<D, Y*>& lhs, T scalar)
{
	for (int i = 0; i < lhs.length(); ++i)
	{
		lhs[i] *= scalar;
	}
	return lhs;
}


namespace symphas::internal
{
	template<typename T>
	auto new_abs_array(const T* data_y, const len_type len)
	{
		using abs_type = decltype(std::abs(std::declval<T>()));
		abs_type* data_abs = new abs_type[len];

		using std::abs;
		for (iter_type i = 0; i < len; ++i)
		{
			data_abs[i] = abs(data_y[i]);
		}
		return data_abs;
	
	
	}


	template<typename T>
	auto new_radial_avg_field(const T* data_y, len_type L, double dx)
	{
		using Y = std::remove_const_t<T>;

		iter_type cL = (L % 2 == 0) ? 1 : 0;
		iter_type
			L_left = L / 2,
			L_right = L / 2 - cL;

		len_type dims[] = { L / 2 + 1 };
		len_type len = symphas::lib::radial_length(dims);

		iter_type* dmap = symphas::lib::make_radial_index_map(dims);
		auto [data_x, cmap] = symphas::lib::make_radial_arrays(dims, dx);
		delete[] cmap;
		
		
		Y* avgs = new Y[len]{ 0 };
		size_t* cs = new size_t[len]{ 0 };


		iter_type center_index = L_left;
		avgs[dmap[0]] = data_y[center_index];

		iter_type n = 1;
		for (iter_type i = 0; i < L_right; ++i, ++n)
		{
			avgs[dmap[n]] = (data_y[center_index + n] + data_y[center_index - n]);
			cs[dmap[n]] = 2;
		}

		if (cL > 0)
		{
			avgs[dmap[n]] = data_y[center_index - n];
			cs[dmap[n]] = 1;
		}


		return Field{ data_x, avgs, len };
	}

	template<typename T>
	auto new_radial_avg_field(const T* data_y, len_type L, len_type M, double dx, double dy)
	{
		using Y = std::remove_const_t<T>;

		iter_type cL = (L % 2 == 0) ? 1 : 0;
		iter_type cM = (M % 2 == 0) ? 1 : 0;
		iter_type
			L_left = L / 2,
			L_right = L / 2 - cL,
			M_left = M / 2,
			M_right = M / 2 - cM;


		len_type dims[] = { L / 2 + 1, M / 2 + 1 };
		len_type len = symphas::lib::radial_length(dims);

		iter_type* dmap = symphas::lib::make_radial_index_map(dims);
		auto [data_x, cmap] = symphas::lib::make_radial_arrays(dims, dx, dy);
		delete[] cmap;

		Y* avgs = new Y[len]{ 0 };
		size_t* cs = new size_t[len]{ 0 };


		iter_type center_index = L * M_left + L_left;
		avgs[dmap[0]] = data_y[center_index];
		cs[dmap[0]] += 1;

		iter_type n = 1;
		for (iter_type i = 0; i < L_right; ++i, ++n)
		{
			avgs[dmap[n]] += (data_y[center_index + n] + data_y[center_index - n]);
			cs[dmap[n]] += 2;
		}

		if (cL > 0)
		{
			avgs[dmap[n]] += data_y[center_index - n];
			cs[dmap[n]] += 1;
		}


		iter_type
			bl = L * (M_left - 1) + (L_left - 1),
			tl = L * (M_left + 1) + (L_left - 1),
			br = L * (M_left - 1) + (L_left + 1),
			tr = L * (M_left + 1) + (L_left + 1);



		// Add the data which is offset from the axes.
		for (iter_type j = 0; j < M_right; ++j)
		{
			iter_type m = (j + 1) * dims[0];
			avgs[dmap[m]] += (data_y[bl + 1] + data_y[tl + 1]);
			cs[dmap[m]] += 2;

			n = m + 1;
			for (iter_type i = 0; i < L_right; ++i, ++n)
			{
				Y value = 
					data_y[bl--] + data_y[tl--]
					+ data_y[br++] + data_y[tr++];

				avgs[dmap[n]] += value;
				cs[dmap[n]] += 4;
			}

			if (cL > 0)
			{
				avgs[dmap[n]] += (data_y[bl--] + data_y[tl--]);
				cs[dmap[n]] += 2;
			}

			bl -= L_right + 1;
			tl += L + L_left;
			br -= L + L_right;
			tr += L_left + 1;
		}

		if (cM > 0)
		{
			iter_type m = dims[0] * (dims[1] - 1);
			avgs[dmap[m]] += data_y[bl + 1];
			cs[dmap[m]] += 1;

			n = m + 1;
			for (iter_type i = 0; i < L_right; ++i, ++n)
			{
				avgs[dmap[n]] += (data_y[bl--] + data_y[br++]);
				cs[dmap[n]] += 2;
			}
			if (cL > 0)
			{
				avgs[dmap[n]] += data_y[bl--];
				cs[dmap[n]] += 1;
			}
		}

		iter_type duplicates = 0;

		for (iter_type i = 0; i < len - duplicates; ++i)
		{
			cs[i] = cs[i];
			iter_type j = i + 1;
			while (std::abs(data_x[i] - data_x[j]) < 1E-3 
				&& j < len - duplicates)
			{
				avgs[i] += avgs[j];
				cs[i] += cs[j];
				++j;
			}

			if (j > i + 1)
			{
				iter_type dist = j - i - 1;
				for (iter_type k = j; k < len - duplicates; ++k)
				{
					avgs[k - dist] = avgs[k];
					cs[k - dist] = cs[k];
					data_x[k - dist] = data_x[k];
				}
				duplicates += dist;
			}
			avgs[i] = avgs[i] * (1.0 / cs[i]);
		}
		delete[] dmap;
		delete[] cs;

		return Field{ data_x, avgs, len - duplicates };
	}



	template<typename T>
	auto new_radial_avg_field(const T* data_y, len_type L, len_type M, len_type N, double dx, double dy, double dz)
	{
		using Y = std::remove_const_t<T>;

		iter_type cL = (L % 2 == 0) ? 1 : 0;
		iter_type cM = (M % 2 == 0) ? 1 : 0;
		iter_type cN = (N % 2 == 0) ? 1 : 0;
		iter_type
			L_left = L / 2,
			L_right = L / 2 - cL,
			M_left = M / 2,
			M_right = M / 2 - cM,
			N_left = N / 2,
			N_right = N / 2 - cN;


		len_type dims[] = { L / 2 + 1, M / 2 + 1, N / 2 + 1 };
		len_type len = symphas::lib::radial_length(dims);

		iter_type* dmap = symphas::lib::make_radial_index_map(dims);
		auto [data_x, cmap] = symphas::lib::make_radial_arrays(dims, dx, dy, dz);
		delete[] cmap;

		Y* avgs = new Y[len]{ 0 };
		size_t* cs = new size_t[len]{ 0 };


		/* Iteration of the data along the x axis.
		 */
		iter_type center_index = L * M * N_left + L * M_left + L_left;
		avgs[dmap[0]] += data_y[center_index];
		cs[dmap[0]] += 1;

		iter_type n = 1;
		for (iter_type i = 0; i < L_right; ++i, ++n)
		{
			avgs[dmap[n]] += (data_y[center_index + n] + data_y[center_index - n]);
			cs[dmap[n]] += 2;
		}

		if (cL > 0)
		{
			avgs[dmap[n]] += data_y[center_index - n];
			cs[dmap[n]] += 1;
			n += 1;
		}



		/* Iteration of the data along the x-y plane, excluding the x axis itself.
		 */
		iter_type center_index_u = center_index;
		iter_type center_index_d = center_index;
		for (iter_type j = 0; j < M_right; ++j)
		{
			center_index_u += L;
			center_index_d -= L;

			avgs[dmap[n]] += (data_y[center_index_u] + data_y[center_index_d]);
			cs[dmap[n]] += 2;
			n += 1;

			for (iter_type i = 0; i < L_right; ++i, ++n)
			{
				avgs[dmap[n]] +=
					(data_y[center_index_u + i + 1] + data_y[center_index_u - i - 1]
						+ data_y[center_index_d + i + 1] + data_y[center_index_d - i - 1]);
				cs[dmap[n]] += 4;
			}

			if (cL > 0)
			{
				avgs[dmap[n]] += (data_y[center_index_u - L_right - 1] 
					+ data_y[center_index_d - L_right - 1]);
				cs[dmap[n]] += 2;
				n += 1;
			}
		}

		if (cM > 0)
		{
			center_index_d -= L;
			avgs[dmap[n]] += data_y[center_index_d];
			cs[dmap[n]] += 1;
			n += 1;

			for (iter_type i = 0; i < L_right; ++i, ++n)
			{
				avgs[dmap[n]] += (data_y[center_index_d + i + 1] + data_y[center_index_d - i - 1]);
				cs[dmap[n]] += 2;
			}

			if (cL > 0)
			{
				avgs[dmap[n]] += data_y[center_index_d - L_right - 1];
				cs[dmap[n]] += 1;
			}
		}




		iter_type
			fbl = L * M * (N_left - 1) + L * (M_left - 1) + (L_left - 1),
			ftl = L * M * (N_left - 1) + L * (M_left + 1) + (L_left - 1),
			fbr = L * M * (N_left - 1) + L * (M_left - 1) + (L_left + 1),
			ftr = L * M * (N_left - 1) + L * (M_left + 1) + (L_left + 1),
			bbl = L * M * (N_left + 1) + L * (M_left - 1) + (L_left - 1),
			btl = L * M * (N_left + 1) + L * (M_left + 1) + (L_left - 1),
			bbr = L * M * (N_left + 1) + L * (M_left - 1) + (L_left + 1),
			btr = L * M * (N_left + 1) + L * (M_left + 1) + (L_left + 1);


		for (iter_type k = 0; k < N_right; ++k)
		{

			/* Iteration of dmap starts at first row in the next depth layer.
			 * This targets the cells which are in z-x plane, first on the z
			 * axis and then on the rest of the plane.
			 */

			iter_type m = (k + 1) * dims[1] * dims[0];

			avgs[dmap[m]] += (data_y[fbl + L + 1] + data_y[bbl + L + 1]);
			cs[dmap[m]] += 2;
			m += 1;

			for (iter_type i = 0; i < L_right; ++i, ++m)
			{
				avgs[dmap[m]] +=
					(data_y[fbl + L - i] + data_y[fbr + L + i]
						+ data_y[bbl + L - i] + data_y[bbr + L + i]);
				cs[dmap[m]] += 4;
			}

			if (cL > 0)
			{
				avgs[dmap[m]] += (data_y[fbl + L - L_right] + data_y[bbl + L - L_right]);
				cs[dmap[m]] += 2;
				m += 1;
			}

			for (iter_type j = 0; j < M_right; ++j)
			{
				/*
				 * Continue iteration on the y-z plane.
				 */

				avgs[dmap[m]] +=
					(data_y[fbl + 1] + data_y[ftl + 1]
						+ data_y[bbl + 1] + data_y[btl + 1]);
				cs[dmap[m]] += 4;

				n = m + 1;
				for (iter_type i = 0; i < L_right; ++i, ++n)
				{
					Y value =
						data_y[fbl--] + data_y[ftl--]
						+ data_y[fbr++] + data_y[ftr++]
						+ data_y[bbl--] + data_y[btl--]
						+ data_y[bbr++] + data_y[btr++];

					avgs[dmap[n]] += value;
					cs[dmap[n]] += 8;
				}

				if (cL > 0)
				{
					avgs[dmap[n]] +=
						(data_y[fbl--] + data_y[ftl--]
							+ data_y[bbl--] + data_y[btl--]);
					cs[dmap[m]] += 4;
					n += 1;
				}
				m = n;

				bbl -= L_right + 1;
				btl += L + L_left;
				bbr -= L + L_right;
				btr += L_left + 1;

				fbl -= L_right + 1;
				ftl += L + L_left;
				fbr -= L + L_right;
				ftr += L_left + 1;
			}


			if (cM > 0)
			{
				avgs[dmap[m]] += (data_y[fbl + 1] + data_y[bbl + 1]);
				cs[dmap[m]] += 2;

				n = m + 1;
				for (iter_type i = 0; i < L_right; ++i, ++n)
				{
					Y value =
						data_y[fbl--] + data_y[bbl--]
						+ data_y[fbr++] + data_y[bbr++];

					avgs[dmap[n]] += value;
					cs[dmap[n]] += 4;
				}

				if (cL > 0)
				{
					avgs[dmap[n]] += (data_y[bbl--] + data_y[fbl--]);
					cs[dmap[n]] += 2;

					n += 1;
				}
				m = n;

				bbl -= L_right + 1;
				bbr -= L + L_right;
				fbl -= L_right + 1;
				fbr -= L + L_right;
			}


			bbl += L * M + L * M_left;
			bbr += L * M + L * M_left;
			btl += L * (M_left + 1);
			btr += L * (M_left + 1);

			fbl -= L * (M_right + 1);
			fbr -= L * (M_right + 1);
			ftl -= L * M + L * M_right;
			ftr -= L * M + L * M_right;
		}

		if (cN > 0)
		{
			iter_type m = (N_right + 1) * dims[1] * dims[0];

			avgs[dmap[m]] += data_y[fbl + L + 1];
			cs[dmap[m]] += 1;
			m += 1;

			for (iter_type i = 0; i < L_right; ++i, ++m)
			{
				avgs[dmap[m]] += (data_y[fbl + L - i] + data_y[fbr + L + i]);
				cs[dmap[m]] += 2;
			}

			if (cL > 0)
			{
				avgs[dmap[m]] += data_y[fbl + L - L_right];
				cs[dmap[m]] += 1;
				m += 1;
			}

			for (iter_type j = 0; j < M_right; ++j)
			{
				/*
				 * Continue iteration on the y-z plane.
				 */

				avgs[dmap[m]] += (data_y[fbl + 1] + data_y[ftl + 1]);
				cs[dmap[m]] += 2;

				n = m + 1;
				for (iter_type i = 0; i < L_right; ++i, ++n)
				{
					Y value =
						data_y[fbl--] + data_y[ftl--]
						+ data_y[fbr++] + data_y[ftr++];

					avgs[dmap[n]] += value;
					cs[dmap[n]] += 4;
				}

				if (cL > 0)
				{
					avgs[dmap[n]] += (data_y[fbl--] + data_y[ftl--]);
					cs[dmap[n]] += 2;

					n += 1;
				}
				m = n;

				fbl -= L_right + 1;
				ftl += L + L_left;
				fbr -= L + L_right;
				ftr += L_left + 1;
			}

			if (cM > 0)
			{
				avgs[dmap[m]] += data_y[fbl + 1];
				cs[dmap[m]] += 1;

				n = m + 1;
				for (iter_type i = 0; i < L_right; ++i, ++n)
				{
					Y value = data_y[fbl--] + data_y[fbr++];
					avgs[dmap[n]] += value;
					cs[dmap[n]] += 2;
				}

				if (cL > 0)
				{
					avgs[dmap[n]] += data_y[fbl--];
					cs[dmap[n]] += 1;
					n += 1;
				}
				m = n;

				fbl -= L_right + 1;
				fbr -= L + L_right;
			}
		}

		iter_type duplicates = 0;
		for (iter_type i = 0; i < len - duplicates; ++i)
		{
			cs[i] = cs[i];
			iter_type j = i + 1;
			while (std::abs(data_x[i] - data_x[j]) < 1E-3
				&& j < len - duplicates)
			{
				avgs[i] += avgs[j];
				cs[i] += cs[j];
				++j;
			}

			if (j > i + 1)
			{
				iter_type dist = j - i - 1;
				for (iter_type k = j; k < len - duplicates; ++k)
				{
					avgs[k - dist] = avgs[k];
					cs[k - dist] = cs[k];
					data_x[k - dist] = data_x[k];
				}
				duplicates += dist;
			}
			avgs[i] = avgs[i] * (1.0 / cs[i]);

		}
		delete[] dmap;
		delete[] cs;

		return Field{ data_x, avgs, len - duplicates };
	}






	template<typename T>
	auto linear_interpolate(axis_coord_t ax, T ay, axis_coord_t bx, T by, axis_coord_t point)
	{
		T slope = (by - ay) * (1.0 / (bx - ax));
		T intercept = ay - slope * ax;

		return slope * point + intercept;
	}
}



namespace symphas::lib
{

	//! Returns a new dataset by trimming part of the given data.
	/*!
	 * For a normal data set taken as a function of x-axis only, trim the
	 * sides according to the given ratio. for example, 1/3 means that a
	 * sixth of the total data on either side is cut, so that 2/3 will remain.
	 * 
	 * \param data The data that is trimmed.
	 * \param ratio The total fraction of data that is removed.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> trim_sides(std::vector<std::pair<axis_1d_type, T>> const& data, double ratio)
	{
		size_t
			len = data.size(),
			chop = static_cast<size_t>((ratio / 2.0) * len);
		auto
			first = data.begin() + chop,
			last = data.end() - chop;

		return std::vector<std::pair<double, T>>(first, last);
	}

	//! Returns a new dataset by trimming part of the given data.
	/*!
	 * For a normal data set taken as a function of x-axis only, trim the
	 * right side according to the given ratio. for example, 1/3 means that a
	 * third of the total data is cut on the right side, so 2/3 will remain.
	 * 
	 * \param data The data that is trimmed.
	 * \param ratio The total fraction of data that is removed.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> trim_right(std::vector<std::pair<axis_1d_type, T>> const& data, double ratio)
	{
		size_t
			len = data.size(),
			chop = static_cast<size_t>(ratio * len);
		auto
			first = data.begin(),
			last = data.end() - chop;

		return std::vector<std::pair<double, T>>(first, last);
	}

	//! Returns a new dataset by trimming part of the given data.
	/*!
	 * For a normal data set taken as a function of x-axis only, trim the
	 * left side according to the given ratio. for example, 1/3 means that a
	 * third of the total data is cut on the right side, so 2/3 will remain.
	 * 
	 * \param data The data that is trimmed.
	 * \param ratio The total fraction of data that is removed.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> trim_left(std::vector<std::pair<axis_1d_type, T>> const& data, double ratio)
	{
		size_t
			len = data.size(),
			chop = static_cast<size_t>(ratio * len);
		auto
			first = data.begin() + chop,
			last = data.end();

		return std::vector<std::pair<double, T>>(first, last);
	}
	


	//! Smooth the data by combining datapoints.
	/*!
	 * Reduces the number of points in the data and produces a smoothed result 
	 * by averaging over a number of x-values as specified by the function 
	 * argument a ratio of 1/2 means half the points are kept.
	 * 
	 * \param data The data which is averaged.
	 * \param ratio The proportion of data which is removed.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> smooth_data(std::vector<std::pair<axis_1d_type, T>> const& data, double ratio = 0.5)
	{
		ratio = std::min(ratio, 1.0);

		len_type len = static_cast<len_type>(data.size());
		len_type data_len = static_cast<len_type>(len * ratio);

		auto split = [&](iter_type i)
		{
			iter_type
				count = len / data_len,
				rem = len - (data_len * count);
			return count * i + std::min(i, rem);
		};

		std::vector<std::pair<axis_1d_type, T>> data_avg;
		data_avg.reserve(data_len);

		for (iter_type i = 0; i < data_len; ++i)
		{
			T value{ 0 };
			for (iter_type ii = split(i); ii < split(i + 1); ++ii)
			{
				value += data[ii].second;
			}
			value *= 1.0 / (split(i + 1) - split(i));
			data_avg.emplace_back((data[split(i)].first + data[split(i + 1) - 1].first) / 2.0, value);
		}

		return data_avg;
	}

	template<typename T>
	std::vector<std::pair<double, T>> smooth_data(std::vector<std::pair<axis_1d_type, T>> const& data, len_type points)
	{
		return smooth_data(data, static_cast<double>(points) / data.size());
	}

	//! Create a data series where the dependent axis has equal separation. 
	/*!
	 * For datasets where the separation of the \f$x\f$ axis elements is
	 * uneven, use the given separation and apply it evenly to the entire 
	 * dataset. Values between axis are interpolated linearly.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> balance_axis(std::vector<std::pair<axis_1d_type, T>> const& data, axis_coord_t dx)
	{
		std::vector<std::pair<axis_1d_type, T>> data_bal;
		data_bal.reserve(data.size());
		data_bal.emplace_back(data[0]);

		axis_coord_t cursor_x = data[0].first;
		axis_coord_t final_x = data.back().first;
		axis_coord_t next_x = data_bal.back().first;

		iter_type i = 1;
		while (cursor_x < final_x)
		{
			axis_coord_t last_x = next_x;
			next_x = last_x + dx;


			T last_y = data_bal.back().second;
			T cursor_y = data[i].second;

			// Integrate over the interval until the cursor is larger than
			// the one we need to get to. Then choose a next_y to match the
			// value of the integration.
			if (cursor_x < next_x)
			{
				T sum = cursor_y * (cursor_x - last_x);
				while (cursor_x < next_x && i + 1 < data.size())
				{
					axis_coord_t cursor_last = cursor_x;
					cursor_x = data[++i].first;
					sum += data[i].second * (cursor_x - cursor_last);
				}
				sum -= data[i].second * (cursor_x - next_x);
				
				axis_coord_t next_y = (2.0 / dx) * (sum - last_y * dx + 0.5 * dx * last_y);
				data_bal.emplace_back(next_x, next_y);
			}
			// Otherwise, interpolate the value linearly.
			else
			{
				axis_coord_t next_y = 
					symphas::internal::linear_interpolate(
						last_x, last_y,
						cursor_x, cursor_y,
						next_x);
				data_bal.emplace_back(next_x, next_y);
			}
		}

		return data_bal;
	}

	//! Create a data series where the axis is given the largest separation.
	/*!
	 * For datasets where the separation of the \f$x\f$ axis elements is
	 * uneven, choose the largest separation and apply it evenly to the entire
	 * dataset.
	 *
	 * \param data The original dataset which is averaged.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> balance_axis_ls(std::vector<std::pair<axis_1d_type, T>> const& data)
	{
		axis_coord_t x = data.front().first;
		axis_coord_t dx = std::numeric_limits<axis_coord_t>::min();

		for (iter_type i = 1; i < data.size(); ++i)
		{
			axis_coord_t dx_current = std::abs(x - data[i].first);
			if (dx_current > dx)
			{
				dx = dx_current;
			}
			x = data[i].first;
		}
		return balance_axis(data, dx);
	}

	//! Create a data series where the dependent axis has equal separation. 
	/*!
	 * For datasets where the separation of the \f$x\f$ axis elements is
	 * uneven, choose the average separation and apply it evenly to the entire
	 * dataset.
	 * 
	 * \param data The original dataset which is averaged.
	 */
	template<typename T>
	std::vector<std::pair<double, T>> balance_axis_as(std::vector<std::pair<axis_1d_type, T>> const& data)
	{
		axis_coord_t dx = (data.front().first + data.back().first) / data.size();
		return balance_axis(data, dx);
	}
	

	template<>
	inline auto new_system_axis_list<1>(symphas::interval_data_type const& intervals)
	{
		auto& ix = intervals.at(Axis::X);
		len_type len = ix.get_count();
		auto dx = ix.length() / len;

		axis_1d_type* data_x = new axis_1d_type[len];

		double x = ix.left();
		for (auto* it = data_x; it < data_x + len; ++it)
		{
			*it = axis_1d_type{ x };
			x += dx;
		}
		return data_x;
	}

	template<>
	inline auto new_system_axis_list<2>(symphas::interval_data_type const& intervals)
	{
		auto& ix = intervals.at(Axis::X);
		auto& iy = intervals.at(Axis::Y);

		len_type count_x = ix.get_interval_count();
		len_type count_y = iy.get_interval_count();

		len_type len_x = ix.interval_length();
		len_type len_y = iy.interval_length();

		auto dx = len_x / count_x;
		auto dy = len_y / count_y;

		axis_2d_type* data_x = new axis_2d_type[count_x * count_y];

		iter_type n = 0;
		double y = iy.left();
		for (iter_type j = 0; j < count_y; ++j)
		{
			double x = ix.left();
			for (iter_type i = 0; i < count_x; ++i, ++n)
			{
				data_x[n][0] = x;
				data_x[n][1] = y;
				x += dx;
				x = (x > DOMAIN_Xn) ? x - len_x : (x < DOMAIN_X0) ? x + len_x : x;
			}
			y += dy;
			y = (y > DOMAIN_Yn) ? y - len_y : (y < DOMAIN_Y0) ? y + len_y : y;
		}
		return data_x;
	}


	template<>
	inline auto new_system_axis_list<3>(symphas::interval_data_type const& intervals)
	{
		auto& ix = intervals.at(Axis::X);
		auto& iy = intervals.at(Axis::Y);
		auto& iz = intervals.at(Axis::Z);

		len_type lenx = ix.get_count();
		len_type leny = iy.get_count();
		len_type lenz = iz.get_count();

		auto dx = ix.length() / lenx;
		auto dy = iy.length() / leny;
		auto dz = iz.length() / lenz;

		axis_3d_type* data_x = new axis_3d_type[lenx * leny * lenz];


		iter_type n = 0;
		double z = iz.left();
		for (iter_type k = 0; k < lenz; ++k)
		{
			double y = iy.left();
			for (iter_type j = 0; j < leny; ++j)
			{
				double x = ix.left();
				for (iter_type i = 0; i < lenx; ++i, ++n)
				{
					data_x[n][0] = x;
					data_x[n][1] = y;
					data_x[n][2] = z;
					x += dx;
				}
				y += dy;
			}
			z += dz;
		}
		return data_x;
	}


	template<size_t D>
	auto new_system_axis_list(const len_type* dims)
	{
		symphas::interval_data_type intervals;
		for (iter_type i = 0; i < D; ++i)
		{
			symphas::interval_element_type interval(0, dims[i] - 1, 1);
			intervals[symphas::index_to_axis(i)] = interval;
		}

		return new_system_axis_list<D>(intervals);
	}





	//! For an list of data, return the dimensions.
	/*!
	 * For an list of data, return the dimensions. The data is assumed to be
	 * sorted.
	 *
	 * \param sorted_data The sorted list of data from which the dimensions
	 * are inferred.
	 */
	template<typename T>
	inline len_type get_sorted_dimensions(std::vector<std::pair<axis_1d_type, T>> sorted_data)
	{
		return static_cast<iter_type>(sorted_data.size());
	}


	template<typename T>
	inline grid::dim_list get_sorted_dimensions(std::vector<std::pair<axis_2d_type, T>> sorted_data)
	{
		len_type L = 1, M = 0;
		double y0 = sorted_data.front().first[1];
		for (auto it = sorted_data.begin() + 1; it < sorted_data.end() && y0 == it->first[1]; ++it)
		{
			++L;
		}
		M = static_cast<len_type>(sorted_data.size()) / L;

		return { L, M };

	}



	inline grid::dim_list get_sorted_dimensions(const axis_3d_type* sorted_data, len_type len)
	{
		len_type L = 0, M = 0, N = 0;
		double
			z0 = (*sorted_data)[2],
			y0 = (*sorted_data)[1];
		for (auto* it = sorted_data + 1; it < sorted_data + len; ++it)
		{
			if (y0 != (*it)[1])
			{
				L = static_cast<iter_type>(it - sorted_data);
				break;
			}
		}
		for (auto* it = sorted_data + L; it < sorted_data + len; it += L)
		{
			if (z0 != (*it)[2])
			{
				M = static_cast<iter_type>(it - sorted_data) / L;
				break;
			}
		}
		N = len / (M * L);

		return { L, M, N };

	}

	inline grid::dim_list get_sorted_dimensions(const axis_2d_type* sorted_data, len_type len)
	{
		len_type L = 0, M = 0;
		double y0 = (*sorted_data)[1];
		while (y0 == (*sorted_data++)[1])
		{
			++L;
		}
		M = len / L;

		return { L, M };

	}

	inline len_type get_sorted_dimensions(const axis_1d_type*, len_type len)
	{
		return len;
	}


	// ****************************************************************************************


	//! Construct a new field with the absolute value of the data.
	/*!
	 * Compute the absolute value of each of the dependent values and
	 * return a new field. There is no dependency with the \f$x\f$-axis values.
	 * 
	 * \param data The data which is transformed by the absolute value.
	 */
	template<typename X, typename Y>
	auto abs(Field<X, Y*> const& data)
	{
		X data_x{ data.data_x() };

		auto* abs_data = symphas::internal::new_abs_array(data.y, data.length());
		using Ya = decltype(abs_data);
		return Field<X, Ya>(data_x, abs_data, data.length());
	}

	//! Construct a new field with the absolute value of the data.
	/*!
	 * Compute the absolute value of each of the dependent values and
	 * return a new field. There is no dependency with the \f$x\f$-axis values.
	 *
	 * \param data The data which is transformed by the absolute value.
	 */
	template<typename X, typename Y>
	auto abs(Field<X*, Y*> const& data)
	{
		X* data_x = new X[data.length()];
		std::copy(data.x, data.x + data.length(), data_x);

		auto* abs_data = symphas::internal::new_abs_array(data.y, data.length());
		using Ya = decltype(abs_data);
		return Field<X*, Ya>(data_x, abs_data, data.length());
	}

	//! Construct a new field with the absolute value of the data.
	/*!
	 * Compute the absolute value of each of the dependent values and
	 * return a new field. There is no dependency with the \f$x\f$-axis values.
	 *
	 * \param data The data which is transformed by the absolute value.
	 */
	template<size_t D, typename Y>
	auto abs(FieldAxis<D, Y> const& data)
	{
		axis_nd_t<D>* data_x = new axis_nd_t<D>[data.length()];
		std::copy(data.x, data.x + data.length(), data_x);

		auto* abs_data = symphas::internal::new_abs_array(data.y, data.length());
		using Ya = decltype(abs_data);
		return FieldAxis<D, Ya>(data_x, abs_data, data.length());
	}

	//! Construct data of the radial average of the given system.
	/*!
	 * The radial average as measured from the center point (not the same
	 * as the origin), of the given 1-dimensional system is computed.
	 *
	 * \param data The 1-dimensional system of which the radial average is
	 * computed.
	 */
	template<typename Y>
	auto radial_avg(FieldAxis<1, Y> const& data)
	{
		auto L = symphas::lib::get_sorted_dimensions(data.x, data.length());

		double const
			dx = ((*(data.x + 1)) - (*data.x));

		return symphas::internal::new_radial_avg_field(data.y, L, dx);

	}

	//! Construct data of the radial average of the given system.
	/*!
	 * The radial average as measured from the center point (not the same
	 * as the origin), of the given 2-dimensional system is computed.
	 *
	 * \param data The 2-dimensional system of which the radial average is
	 * computed.
	 */
	template<typename Y>
	auto radial_avg(FieldAxis<2, Y> const& data)
	{
		auto [L, M] = symphas::lib::get_sorted_dimensions(data.x, data.length())._2();

		double const
			dx = ((*(data.x + 1))[0] - (*data.x)[0]),
			dy = ((*(data.x + L))[1] - (*data.x)[1]);

		return symphas::internal::new_radial_avg_field(data.y, L, M, dx, dy);

	}

	//! Construct data of the radial average of the given system.
	/*!
	 * The radial average as measured from the center point (not the same
	 * as the origin), of the given 3-dimensional system is computed.
	 * 
	 * \param data The 3-dimensional system of which the radial average is 
	 * computed.
	 */
	template<typename Y>
	auto radial_avg(FieldAxis<3, Y> const& data)
	{
		auto [L, M, N] = symphas::lib::get_sorted_dimensions(data.x, data.length())._3();

		double const
			dx = ((*(data.x + 1))[0] - (*data.x)[0]),
			dy = ((*(data.x + L))[1] - (*data.x)[1]),
			dz = ((*(data.x + L * M))[2] - (*data.x)[2]);

		return symphas::internal::new_radial_avg_field(data.y, L, M, N, dx, dy, dz);

	}

	template<typename T>
	grid::dim_list get_sorted_dimensions(std::vector<std::pair<axis_3d_type, T>> sorted_data)
	{
		len_type L = 0, M = 0, N = 0;
		double
			z0 = sorted_data.front().first[2],
			y0 = sorted_data.front().first[1];
		for (auto it = sorted_data.begin() + 1; it < sorted_data.end(); ++it)
		{
			if (y0 != it->first[1])
			{
				L = static_cast<iter_type>(it - sorted_data.begin());
				break;
			}
		}
		for (auto it = sorted_data.begin() + L; it < sorted_data.end(); it += L)
		{
			if (z0 != it->first[2])
			{
				M = static_cast<iter_type>(it - sorted_data.begin()) / L;
				break;
			}
		}
		N = static_cast<len_type>(sorted_data.size()) / (M * L);

		return { L, M, N };

	}

	//! Return the min, max and separation distance of the given axis.
	template<size_t D>
	auto get_axis_mms(const axis_nd_t<D>* axis_values, len_type len)
	{
		std::array<axis_coord_t, D> min{}, max{}, sep{};
		std::fill(min.begin(), min.end(), std::numeric_limits<axis_coord_t>::max());
		std::fill(max.begin(), max.end(), std::numeric_limits<axis_coord_t>::min());
		std::fill(sep.begin(), sep.end(), std::numeric_limits<axis_coord_t>::max());

		for (iter_type i = 0; i < len; ++i)
		{
			for (iter_type n = 0; n < D; ++n)
			{
				double dd = std::abs(min[n] - axis_values[i][n]);
				sep[n] = (dd > 0 && dd < sep[n]) ? dd : sep[n];

				min[n] = (axis_values[i][n] < min[n]) ? axis_values[i][n] : min[n];
				max[n] = (axis_values[i][n] > max[n]) ? axis_values[i][n] : max[n];
			}
		}

		return std::make_tuple(min, max, sep);
	}

	template<>
	inline auto get_axis_mms<1>(const axis_nd_t<1>* axis_values, len_type len)
	{
		std::array<axis_coord_t, 1> min = { std::numeric_limits<axis_coord_t>::max() };
		std::array<axis_coord_t, 1> max = { std::numeric_limits<axis_coord_t>::min() };
		std::array<axis_coord_t, 1> sep = { std::numeric_limits<axis_coord_t>::max() };

		for (iter_type i = 0; i < len; ++i)
		{
			double dd = std::abs(min[0] - axis_values[i]);
			sep[0] = (dd > 0 && dd < sep[0]) ? dd : sep[0];

			min[0] = (axis_values[i] < min[0]) ? axis_values[i] : min[0];
			max[0] = (axis_values[i] > max[0]) ? axis_values[i] : max[0];
		}

		return std::make_tuple(min, max, sep);
	}



	template<size_t D>
	grid::dim_list get_dimensions(const axis_nd_t<D>* axis_values, len_type len)
	{
		auto [min, max, sep] = get_axis_mms<D>(axis_values, len);
		len_type dim[D];
		for (iter_type n = 0; n < D; ++n)
		{
			dim[n] = static_cast<len_type>(std::round((max[n] - min[n]) / sep[n])) + 1;
		}
		return grid::dim_list(dim, D);
	}


	template<size_t D, typename T>
	grid::dim_list get_dimensions(std::vector<std::pair<axis_nd_t<D>, T>> const& data)
	{
		axis_nd_t<D>* axis_values = new axis_nd_t<D>[data.size()] {};
		for (iter_type i = 0; i < data.size(); ++i)
		{
			axis_values[i] = data[i].first;
		}

		auto&& dims = get_dimensions<D>(axis_values, static_cast<len_type>(data.size()));
		delete[] axis_values;
		return dims;
	}



	//! Return information on the interval of the given axis data.
	/*!
	 * The given data is expected to be sorted.
	 */
	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_3d_type, T>> const& data, axis_coord_t(&ranges)[6])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first[0];
			ranges[1] = data.back().first[0];
			ranges[2] = data.front().first[1];
			ranges[3] = data.back().first[1];
			ranges[4] = data.front().first[2];
			ranges[5] = data.back().first[2];
		}
	}

	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_2d_type, T>> const& data, axis_coord_t(&ranges)[4])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first[0];
			ranges[1] = data.back().first[0];
			ranges[2] = data.front().first[1];
			ranges[3] = data.back().first[1];
		}
	}

	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_1d_type, T>> const& data, axis_coord_t(&ranges)[2])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first;
			ranges[1] = data.back().first;
		}
	}

	template<size_t D, typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_nd_t<D / 2>, T>> const& data, axis_coord_t(&ranges)[D], T(&extrema)[2])
	{
		if (data.size() > 0)
		{
			fill_sorted_ranges(data, ranges);
			extrema[0] = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return a.second > b.second; })->second;
			extrema[1] = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return a.second < b.second; })->second;
		}
	}




	inline void fill_sorted_ranges(const axis_3d_type* data_x, len_type len, axis_coord_t(&ranges)[6])
	{
		if (len > 0)
		{
			auto [L, M, N] = get_sorted_dimensions(data_x, len)._3();
			ranges[0] = (data_x[0])[0];
			ranges[1] = (data_x[L - 1])[0];
			ranges[2] = (data_x[0])[1];
			ranges[3] = (data_x[L * (M - 1)])[1];
			ranges[4] = (data_x[0])[2];
			ranges[5] = (data_x[L * M * (N - 1)])[2];
		}
	}

	inline void fill_sorted_ranges(const axis_2d_type* data_x, len_type len, axis_coord_t(&ranges)[4])
	{
		if (len > 0)
		{
			auto [L, M] = get_sorted_dimensions(data_x, len)._2();
			ranges[0] = (data_x[0])[0];
			ranges[1] = (data_x[L - 1])[0];
			ranges[2] = (data_x[0])[1];
			ranges[3] = (data_x[L * (M - 1)])[1];
		}
	}

	inline void fill_sorted_ranges(const axis_1d_type* data_x, len_type len, axis_coord_t(&ranges)[2])
	{
		if (len > 0)
		{
			ranges[0] = *(data_x);
			ranges[1] = *(data_x + len - 1);
		}
	}

	template<typename T>
	void fill_sorted_ranges(const T* data, len_type len, T(&ranges)[2])
	{
		if (len > 0)
		{
			auto max = std::max_element(data, data + len, [&] (auto a, auto b) { return a < b; });
			auto min = std::max_element(data, data + len, [&] (auto a, auto b) { return a > b; });
			ranges[0] = min;
			ranges[1] = max;
		}
	}




	template<typename T>
	void fill_ranges(std::vector<std::tuple<iter_type, axis_1d_type, T>> const& data, T(&ranges)[4])
	{
		if (data.size() > 0)
		{
			auto max_it_y = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return std::get<2>(a) < std::get<2>(b); });
			auto min_it_y = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });
			auto max_it_x = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return std::get<1>(a) < std::get<1>(b); });
			auto min_it_x = std::max_element(data.begin(), data.end(), [&] (auto a, auto b) { return std::get<1>(a) > std::get<1>(b); });

			auto max_y = std::get<2>(*max_it_y);
			auto min_y = std::get<2>(*min_it_y);
			auto max_x = std::get<1>(*max_it_x);
			auto min_x = std::get<1>(*min_it_x);

			ranges[0] = static_cast<T>(max_x);
			ranges[1] = static_cast<T>(min_x);
			ranges[2] = max_y;
			ranges[3] = min_y;
		}
	}

	inline void fill_ranges(axis_3d_type const* data_x, len_type len, axis_coord_t(&ranges)[6])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[2] = +INFINITY;
		ranges[4] = +INFINITY;

		ranges[1] = -INFINITY;
		ranges[3] = -INFINITY;
		ranges[5] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < ranges[0])
			{
				ranges[0] = data_x[i][0];
			}
			if (data_x[i][0] > ranges[1])
			{
				ranges[1] = data_x[i][0];
			}
			if (data_x[i][1] < ranges[2])
			{
				ranges[2] = data_x[i][1];
			}
			if (data_x[i][1] > ranges[3])
			{
				ranges[3] = data_x[i][1];
			}
			if (data_x[i][2] < ranges[4])
			{
				ranges[4] = data_x[i][2];
			}
			if (data_x[i][2] > ranges[5])
			{
				ranges[5] = data_x[i][2];
			}
		}
	}


	inline void fill_ranges(axis_2d_type const* data_x, len_type len, axis_coord_t(&ranges)[4])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[2] = +INFINITY;

		ranges[1] = -INFINITY;
		ranges[3] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < ranges[0])
			{
				ranges[0] = data_x[i][0];
			}
			if (data_x[i][0] > ranges[1])
			{
				ranges[1] = data_x[i][0];
			}
			if (data_x[i][1] < ranges[2])
			{
				ranges[2] = data_x[i][1];
			}
			if (data_x[i][1] > ranges[3])
			{
				ranges[3] = data_x[i][1];
			}
		}
	}



	inline void fill_ranges(axis_1d_type const* data_x, len_type len, axis_coord_t(&ranges)[2])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[1] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i] < ranges[0])
			{
				ranges[0] = data_x[i];
			}
			if (data_x[i] > ranges[1])
			{
				ranges[1] = data_x[i];
			}
		}
	}


}












