
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

#include "spslibfftw.h"

namespace symphas::internal
{

	template<size_t D, typename S, typename T>
	struct make_forward_plan
	{
		auto operator()(S* from, T* to, const len_type* dims)
		{
			return symphas::dft::new_fftw_plan<D, S, T>{}(from, to, dims, false);
		}
	};

	template<size_t D>
	struct make_forward_plan<D, complex_t, complex_t>
	{
		auto operator()(complex_t* from, complex_t* to, const len_type* dims)
		{
			return symphas::dft::new_fftw_plan<D, complex_t, complex_t>{}(from, to, dims, false, false);
		}
	};

	template<size_t D, typename S, typename T>
	struct make_backward_plan
	{
		auto operator()(S* from, T* to, const len_type* dims)
		{
			return symphas::dft::new_fftw_plan<D, S, T>{}(from, to, dims, false);
		}
	};

	template<size_t D>
	struct make_backward_plan<D, complex_t, complex_t>
	{
		auto operator()(complex_t* from, complex_t* to, const len_type* dims)
		{
			return symphas::dft::new_fftw_plan<D, complex_t, complex_t>{}(from, to, dims, false, true);
		}
	};
}


//! When updated, performs an in place Fourier Transform of the current data.
/*!
 * Given the specified source type 'S', FourierGrid computes an in place 
 * Fourier transform on its owned data. The underlying data array is also  
 * appropriately sized to accommodate the r2c and c2r FFTW algorithms.
 *
 * In the case that source type is complex and the target type is scalar, 
 * complex data is transformed into real space; this type of transform uses 
 * the c2r routine. NOTE that it will assume that the complex data has an 
 * arrangement consistent with what FFTW expects, i.e. that the complex array 
 * has only half values.
 * 
 * \tparam S The source type that is transformed from.
 * \tparam T The type to transform into.
 * \tparam D The dimension of the grid.
 */
template<typename S, typename T, size_t D>
struct FourierGrid : Grid<T, D>
{
	using parent_type = Grid<T, D>;
	using parent_type::values;
	using parent_type::len;
	using parent_type::dims;

	//! Create a FourierGrid of the prescribed dimensions.
	/*!
	 * Create a FourierGrid of the prescribed dimensions.
	 * 
	 * \param dimensions The dimensions of the grid to create.
	 */
	FourierGrid(const len_type* dimensions) : 
		parent_type{ dimensions },
		p_in_out{ symphas::dft::new_fftw_plan<D, S, T>{}(values, values, dimensions) } {}
	FourierGrid(FourierGrid<S, T, D> const& other) : 
		parent_type{ other },
		p_in_out{ symphas::dft::new_fftw_plan<D, S, T>{}(values, values, other.dims) } {}

	FourierGrid(FourierGrid<S, T, D>&& other) noexcept : FourierGrid()
	{
		swap(*this, other);
	}

	FourierGrid<S, T, D>& operator=(FourierGrid<S, T, D> other) const
	{
		swap(*this, other);
		return *this;
	}


	friend void swap(FourierGrid<S, T, D>& first, FourierGrid<S, T, D>& second)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(&first), *static_cast<parent_type*>(&second));
		swap(first.p_in_out, second.p_in_out);
	}

	//! Return a pointer to the source values.
	/*!
	 * Return a pointer to the array representing the source array of the
	 * Fourier transform.
	 */
	S* values_src_cast()
	{
		return reinterpret_cast<S*>(values);
	}

	const S* values_src_cast() const
	{
		return reinterpret_cast<S*>(values);
	}

	//! Update the grid data.
	/*!
	 * Updates this grid, which computes the in place transform.
	 */
	void update();

	~FourierGrid()
	{
		symphas::dft::fftw_destroy_plan(p_in_out);
	}


protected:

	FourierGrid() : Grid<T, D>{}, p_in_out{ NULL } {}

	fftw_plan p_in_out;
};

//! When updated, performs an in place Fourier Transform of the current data.
/*!
 * Given the specified source type 'S', FourierGrid computes an in place
 * Fourier transform on its owned data. The underlying data array is also
 * appropriately sized to accommodate the r2c and c2r FFTW algorithms.
 *
 * In the case that source type is complex and the target type is scalar,
 * complex data is transformed into real space; this type of transform uses
 * the c2r routine. NOTE that it will assume that the complex data has an
 * arrangement consistent with what FFTW expects, i.e. that the complex array
 * has only half values.
 *
 * \tparam S The source type that is transformed from.
 * \tparam T The type to transform into.
 * \tparam D The dimension of the grid.
 */
template<typename S, typename T, size_t D>
struct FourierGrid<S, any_vector_t<T, D>, D> : Grid<any_vector_t<T, D>, D>
{
	using parent_type = Grid<any_vector_t<T, D>, D>;
	using parent_type::values;
	using parent_type::len;
    using parent_type::dims;

	//! Create a FourierGrid of the prescribed dimensions.
	/*!
	 * Create a FourierGrid of the prescribed dimensions.
	 *
	 * \param dimensions The dimensions of the grid to create.
	 */
	FourierGrid(const len_type* dimensions) :
		parent_type{ dimensions },
		p_in_out{ 0 }
	{
		for (iter_type i = 0; i < D; ++i)
		{
			p_in_out[i] = symphas::dft::new_fftw_plan<D, S, T>{}(values[i], values[i], dims);
		}
	}
	FourierGrid(FourierGrid<S, T, D> const& other) :
		parent_type{ other },
		p_in_out{ 0 }
	{
		for (iter_type i = 0; i < D; ++i)
		{
			p_in_out[i] = symphas::dft::new_fftw_plan<D, S, T>{}(values[i], values[i], dims);
		}
	}

	FourierGrid(FourierGrid<S, T, D>&& other) noexcept : FourierGrid()
	{
		swap(*this, other);
	}

	FourierGrid<S, T, D>& operator=(FourierGrid<S, T, D> other) const
	{
		swap(*this, other);
		return *this;
	}


	friend void swap(FourierGrid<S, T, D>& first, FourierGrid<S, T, D>& second)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(&first), *static_cast<parent_type*>(&second));
		swap(first.p_in_out, second.p_in_out);
	}

	//! Return a pointer to the source values.
	/*!
	 * Return a pointer to the array representing the source array of the
	 * Fourier transform.
	 */
	S* values_src_cast()
	{
		return reinterpret_cast<S*(&)[D]>(values);
	}

	const S* values_src_cast() const
	{
		return reinterpret_cast<S*(&)[D]>(values);
	}

	//! Update the grid data.
	/*!
	 * Updates this grid, which computes the in place transform.
	 */
	void update();

	~FourierGrid()
	{
		for (iter_type i = 0; i < D; ++i)
		{
			symphas::dft::fftw_destroy_plan(p_in_out[i]);
		}
	}


protected:

	FourierGrid() : Grid<T, D>{}, p_in_out{ NULL } {}

	fftw_plan p_in_out[D];
};


//! Base type of grids that transform given data and store the result.
/*!
 * Encompassing type for grids which are paired together in some way under a
 * transformation. These are still fundamentally grids; they inherit from the
 * Grid type.
 * 
 * MapGrid consists of source and target types. The source type is the type of 
 * the original type which will be modified or otherwise transformed by the 
 * MapGrid instance. 
 * 
 * \tparam S The type of the source data.
 * \tparam T The type of the transformed data.
 * \tparam D The dimension of the grid.
 */
template<typename S, typename T, size_t D>
struct MapGrid : Grid<T, D>
{
	using parent_type = Grid<T, D>;
	using parent_type::values;
	using parent_type::len;
	using parent_type::dims;

	MapGrid(S src, len_type const* dimensions) : parent_type{ dimensions }, src{ src } {}
	MapGrid(S src, std::initializer_list<len_type> dimensions) : parent_type{ dimensions }, src{ src } {}
	MapGrid(parent_type const& src_grid) : parent_type{ src_grid.dims }, src{ src_grid.values } {}

	friend void swap(MapGrid<S, T, D>& first, MapGrid<S, T, D>& second)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(&first), *static_cast<parent_type*>(&second));
		swap(first.src, second.src);
	}

	//! Performs an update of the data.
	virtual void update() = 0;


protected:

	MapGrid() : parent_type{}, src{ nullptr } {}

	S src;
};


template<typename T, size_t D>
struct FourierGridTarget : Grid<T, D>
{
	FourierGridTarget(len_type const* dimensions) : Grid<T, D>(dimensions) {};
	FourierGridTarget() : Grid<T, D>() {}
};

//! Computes and stores a Fourier transform of a given data set.
/*!
 * Inherently a Grid type of the prescribed dimension and underlying type.
 * Given source data, computes the Fourier transform and stores the
 * result into its own data. In the cases the source data is real, the FFTW
 * performs the Fourier transform according to the r2c algorithm, and the data
 * has to be subsequently rearranged to account for the result storing only
 * the non-redundant half.
 *
 * An important point is that copying the grid pair will retain a pointer to the
 * source data, i.e. the same data will be operated on by the copy.
 */
template<typename S, typename T, size_t D, 
	template<size_t, typename, typename> typename make_ft_helper = symphas::internal::make_forward_plan>
struct MapGridFourier : FourierGridTarget<T, D>
{
	using parent_type = FourierGridTarget<T, D>;
	using parent_type::parent_type;
	using parent_type::values;
	using parent_type::dims;
	using parent_type::len;
	using parent_type::operator[];

	MapGridFourier() : parent_type{}, src{ 0 }, p_in_out { NULL } {}

	MapGridFourier(len_type const* dimensions) :
		parent_type{ dimensions }, src{ symphas::dft::length<T, D>(dimensions) }, 
		p_in_out{ make_ft_helper<D, S, T>{}(src.values, values, dimensions) }
	{
		len = src.len;
	}

	MapGridFourier(MapGridFourier const& other) :
		parent_type{ other.dims }, src{ other.src.len },
		p_in_out{ make_ft_helper<D, S, T>{}(src.values, values, other.dims) }
	{
		len = src.len;
	}

	MapGridFourier(MapGridFourier<S, T, D>&& other) noexcept : MapGridFourier()
	{
		swap(*this, other);
	}
	MapGridFourier<S, T, D>& operator=(MapGridFourier<S, T, D> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(MapGridFourier<S, T, D, make_ft_helper>& first, MapGridFourier<S, T, D, make_ft_helper>& second)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(&first), *static_cast<parent_type*>(&second));
		swap(first.src, second.src);
		swap(first.p_in_out, second.p_in_out);
	}

	//! Performs the Fourier transform.
	/*!
	 * Performs the Fourier transformation of the given data, and the result
	 * is stored in the member data.
	 */
	void update();

	~MapGridFourier()
	{
		symphas::dft::fftw_destroy_plan(p_in_out);
	}

	Block<S> src;

protected:

	fftw_plan p_in_out;

};

//! Computes and stores a Fourier transform of a given data set.
/*!
 * Inherently a Grid type of the prescribed dimension and underlying type.
 * Given source data, computes the Fourier transform and stores the
 * result into its own data. In the cases the source data is real, the FFTW
 * performs the Fourier transform according to the r2c algorithm, and the data
 * has to be subsequently rearranged to account for the result storing only
 * the non-redundant half.
 *
 * An important point is that copying the grid pair will retain a pointer to the
 * source data, i.e. the same data will be operated on by the copy.
 */
template<typename S, typename T, size_t D, 
	template<size_t, typename, typename> typename make_ft_helper>
struct MapGridFourier<any_vector_t<S, D>, any_vector_t<T, D>, D, make_ft_helper> :
	FourierGridTarget<any_vector_t<T, D>, D>
{
	using parent_type = FourierGridTarget<any_vector_t<T, D>, D>;
	using parent_type::parent_type;
	using parent_type::values;
	using parent_type::dims;
	using parent_type::len;
	using parent_type::operator[];

	MapGridFourier() : parent_type{}, src{ 0 }, p_in_out{ NULL } {}

	MapGridFourier(len_type const* dimensions) :
		parent_type{ dimensions }, src{ symphas::dft::length<T, D>(dimensions) }, p_in_out{ 0 }
	{
		len = src.len;
		for (iter_type i = 0; i < D; ++i)
		{
			p_in_out[i] = make_ft_helper<D, S, T>{}(src.values[i], values[i], dims);
		}
	}

	MapGridFourier(MapGridFourier const& other) :
		parent_type{ other.dims }, src{ other.src.len }, p_in_out{ 0 }
	{
		len = src.len;
		for (iter_type i = 0; i < D; ++i)
		{
			p_in_out[i] = make_ft_helper<D, S, T>{}(src.values[i], values[i], dims);
		}
	}

	MapGridFourier(MapGridFourier<S, T, D>&& other) noexcept : MapGridFourier()
	{
		swap(*this, other);
	}

	MapGridFourier<S, T, D>& operator=(MapGridFourier<S, T, D> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(MapGridFourier<S, T, D>& first, MapGridFourier<S, T, D>& second)
	{
		using std::swap;

		//swap(*static_cast<source_parent_type*>(&first), *static_cast<source_parent_type*>(&second));
		swap(*static_cast<parent_type*>(&first), *static_cast<parent_type*>(&second));
		swap(first.src, second.src);
		swap(first.p_in_out, second.p_in_out);
	}

	//! Performs the Fourier transform.
	/*!
	 * Performs the Fourier transformation of the given data, and the result
	 * is stored in the member data.
	 */
	void update();

	~MapGridFourier()
	{
		for (iter_type i = 0; i < D; ++i)
		{
			symphas::dft::fftw_destroy_plan(p_in_out[i]);
		}
	}

	MultiBlock<D, S> src;

protected:

	fftw_plan p_in_out[D];

};


//! Computes and stores a Fourier transform of a given data set.
/*!
 * Inherently a Grid type of the prescribed dimension and underlying type.
 * Given source data, computes the Fourier transform and stores the
 * result into its own data. In the cases the source data is real, the FFTW
 * performs the Fourier transform according to the r2c algorithm, and the data
 * has to be subsequently rearranged to account for the result storing only
 * the non-redundant half.
 *
 * An important point is that copying the grid pair will retain a pointer to the
 * source data, i.e. the same data will be operated on by the copy.
 */
template<typename S, typename T, size_t D>
struct MapGridInverseFourier : MapGridFourier<S, T, D, symphas::internal::make_backward_plan>
{
	using parent_type = MapGridFourier<S, T, D, symphas::internal::make_backward_plan>;
	using parent_type::parent_type;
	using parent_type::values;
	using parent_type::dims;
	using parent_type::src;
	using parent_type::p_in_out;

	/*!
	 * Performs the Fourier transformation of the given data, and the result
	 * is stored in the member data.
	 */
	void update();

};

//! Computes and stores a Fourier transform of a given data set.
/*!
 * Inherently a Grid type of the prescribed dimension and underlying type.
 * Given source data, computes the Fourier transform and stores the
 * result into its own data. In the cases the source data is real, the FFTW
 * performs the Fourier transform according to the r2c algorithm, and the data
 * has to be subsequently rearranged to account for the result storing only
 * the non-redundant half.
 *
 * An important point is that copying the grid pair will retain a pointer to the
 * source data, i.e. the same data will be operated on by the copy.
 */
template<typename S, typename T, size_t D>
struct MapGridInverseFourier<any_vector_t<S, D>, any_vector_t<T, D>, D>
	: MapGridFourier<any_vector_t<S, D>, any_vector_t<T, D>, D, symphas::internal::make_backward_plan>
{
	using parent_type = MapGridFourier<any_vector_t<S, D>, any_vector_t<T, D>, D, symphas::internal::make_backward_plan>;
	using parent_type::parent_type;
	using parent_type::values;
	using parent_type::dims;
	using parent_type::src;
	using parent_type::p_in_out;

	//! Performs the Fourier transform.
	/*!
	 * Performs the Fourier transformation of the given data, and the result
	 * is stored in the member data.
	 */
	void update();

};

// ****************************************************************************


namespace symphas::internal
{
	template<size_t D>
	inline void transform_fourier_map(fftw_plan p_in_out, complex_t* in, complex_t* out, const len_type* dims)
	{
		symphas::dft::fftw_execute(p_in_out);
	}

	template<size_t D>
	inline void transform_fourier_map(fftw_plan p_in_out, scalar_t* in, complex_t* out, const len_type* dims)
	{
		symphas::dft::fftw_execute(p_in_out);
		//symphas::dft::arrange_fftw_hcts<D>(out, out, dims);
	}

	template<size_t D>
	inline void inverse_transform_fourier_map(fftw_plan p_in_out, complex_t* in, complex_t* out, const len_type* dims)
	{
		symphas::dft::fftw_execute(p_in_out);
		symphas::dft::scale<D>(out, dims);
	}

	template<size_t D>
	inline void inverse_transform_fourier_map(fftw_plan p_in_out, complex_t* in, scalar_t* out, const len_type* dims)
	{
		//symphas::dft::arrange_fftw_sthc<D>(in, in, dims);
		symphas::dft::fftw_execute(p_in_out);
		symphas::dft::scale<D>(out, dims);
	}

	template<size_t D>
	inline void transform_fourier_grid(fftw_plan p_in_out, complex_t* in, complex_t* out, const len_type* dims)
	{
		symphas::dft::fftw_execute(p_in_out);
	}

	template<size_t D>
	inline void transform_fourier_grid(fftw_plan p_in_out, scalar_t* in, complex_t* out, const len_type* dims)
	{
		//symphas::dft::arrange_fftw_stip<D>(in, in, dims);
		symphas::dft::fftw_execute(p_in_out);
	}

}


template<typename S, typename T, size_t D, template<size_t, typename, typename> typename make_ft_helper>
void MapGridFourier<S, T, D, make_ft_helper>::update()
{
	symphas::internal::transform_fourier_map<D>(p_in_out, src.values, values, dims);
}

template<typename S, typename T, size_t D, template<size_t, typename, typename> typename make_ft_helper>
void MapGridFourier<any_vector_t<S, D>, any_vector_t<T, D>, D, make_ft_helper>::update()
{
	for (iter_type i = 0; i < D; ++i)
	{
		symphas::internal::transform_fourier_map<D>(p_in_out[i], src.values[i], values[i], dims);
	}
}


template<typename S, typename T, size_t D>
void MapGridInverseFourier<S, T, D>::update()
{
	symphas::internal::inverse_transform_fourier_map<D>(p_in_out, src.values, values, dims);
}

template<typename S, typename T, size_t D>
void MapGridInverseFourier<any_vector_t<S, D>, any_vector_t<T, D>, D>::update()
{
	for (iter_type i = 0; i < D; ++i)
	{
		symphas::internal::inverse_transform_fourier_map<D>(p_in_out[i], src.values[i], values[i], dims);
	}
}

//! Updating executes the FFTW routine, potentially also rearranging data.
/*!
 * These are explicit specializations of the case when real data is transformed
 * into Fourier space; then a rearrangement of the data needs to take place.
 *
 * The type of transform scalar_t -> complex_t uses the r2c routine, from which
 * this methodology is based off.
 */
template<typename S, typename T, size_t D>
void FourierGrid<S, T, D>::update()
{
	symphas::internal::transform_fourier_grid<D>(p_in_out, values_src_cast(), values, dims);
}


//! Updating executes the FFTW routine, potentially also rearranging data.
/*!
 * These are explicit specializations of the case when real data is transformed
 * into Fourier space; then a rearrangement of the data needs to take place.
 *
 * The type of transform scalar_t -> complex_t uses the r2c routine, from which
 * this methodology is based off.
 */
template<typename S, typename T, size_t D>
void FourierGrid<S, any_vector_t<T, D>, D>::update()
{
	for (iter_type i = 0; i < D; ++i)
	{
		symphas::internal::transform_fourier_grid<D>(p_in_out[i], values_src_cast()[i], values[i], dims);
	}
}




