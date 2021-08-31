
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
struct FourierGrid : Block<T>
{
	using Block<T>::values;
	using Block<T>::len;

	//! Create a FourierGrid of the prescribed dimensions.
	/*!
	 * Create a FourierGrid of the prescribed dimensions.
	 * 
	 * \param dimensions The dimensions of the grid to create.
	 */
	FourierGrid(const len_type* dimensions) : 
		Block<T>{ symphas::dft::len<S, D>(dimensions) }, dims{ 0 },
		p_in_out{ symphas::dft::new_fftw_plan<D, S, T>{}(values, values, dimensions) }
	{
		std::copy(dimensions, dimensions + D, dims);
	}
	FourierGrid(FourierGrid<S, T, D> const& other) : 
		Block<T>{ symphas::dft::len<S, D>(other.dims) }, dims{ 0 },
		p_in_out{ symphas::dft::new_fftw_plan<D, S, T>{}(values, values, other.dims) }
	{
		std::copy(other.dims, other.dims + D, dims);
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

		swap(*static_cast<Block<T>*>(&first), *static_cast<Block<T>*>(&second));
		swap(first.dims, second.dims);
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

	FourierGrid() : Block<T>{}, p_in_out{ NULL }, dims{ 0 } {}

	len_type dims[D];
	fftw_plan p_in_out;
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
	using Grid<T, D>::values;
	using Grid<T, D>::len;
	using Grid<T, D>::dims;

	MapGrid(S* src, len_type const* dimensions) : Grid<T, D>{ dimensions }, src{ src } {}
	MapGrid(S* src, std::initializer_list<len_type> dimensions) : Grid<T, D>{ dimensions }, src{ src } {}
	MapGrid(Grid<S, D> const& src_grid) : Grid<T, D>{ src_grid.dims }, src{ src_grid.values } {}

	friend void swap(MapGrid<S, T, D>& first, MapGrid<S, T, D>& second)
	{
		using std::swap;

		swap(*static_cast<Grid<T, D>*>(&first), *static_cast<Grid<T, D>*>(&second));
		swap(first.src, second.src);
	}

	//! Performs an update of the data.
	virtual void update() = 0;

	S* source_data() const
	{
		return src;
	}

	T* transformed_data() const
	{
		return values;
	}


protected:

	MapGrid() : Grid<T, D>{}, src{ nullptr } {}

	S* src;
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
struct MapGridFourier : MapGrid<S, T, D>
{
	using MapGrid<S, T, D>::MapGrid;
	using MapGrid<S, T, D>::src;
	using Grid<T, D>::values;
	using Grid<T, D>::dims;


	MapGridFourier(S const* src, len_type const* dimensions) :
		MapGrid<S, T, D>{ src, dimensions }, p_in_out{ get_plan() } {}
	MapGridFourier(S const* src, std::initializer_list<len_type> dimensions) :
		MapGrid<S, T, D>{ src, dimensions }, p_in_out{ get_plan() } {}
	MapGridFourier(Grid<S, D> const& src_grid) :
		MapGrid<S, T, D>{ src_grid.values, src_grid.dims }, p_in_out{ get_plan() } {}

	MapGridFourier(MapGridFourier const& other) : MapGrid<S, T, D>{ other.src, other.dims }, p_in_out{ get_plan() } {}
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

		swap(*static_cast<MapGrid<S, T, D>*>(&first), *static_cast<MapGrid<S, T, D>*>(&second));
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


protected:

	fftw_plan get_plan()
	{
		return symphas::dft::new_fftw_plan<D, S, T>{}(src, values, dims);
	}

	MapGridFourier() : MapGrid<S, T, D>{}, p_in_out{ NULL } {}
	fftw_plan p_in_out;
};


// ****************************************************************************


template<typename S, typename T, size_t D>
void MapGridFourier<S, T, D>::update()
{
	symphas::dft::fftw_execute(p_in_out);
}

template<>
inline void MapGridFourier<scalar_t, complex_t, 1>::update()
{
	symphas::dft::fftw_execute(p_in_out);
	symphas::dft::arrange_fftw_hcts<1>(values, values, dims);
}

template<>
inline void MapGridFourier<scalar_t, complex_t, 2>::update()
{
	symphas::dft::fftw_execute(p_in_out);
	symphas::dft::arrange_fftw_hcts<2>(values, values, dims);
}

template<>
inline void MapGridFourier<scalar_t, complex_t, 3>::update()
{
	symphas::dft::fftw_execute(p_in_out);
	symphas::dft::arrange_fftw_hcts<3>(values, values, dims);
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
	symphas::dft::fftw_execute(p_in_out);
}

//! Specialization of FourierGrid<S, T, D>::update().
template<>
inline void FourierGrid<scalar_t, complex_t, 1>::update()
{
	symphas::dft::arrange_fftw_stip<1>(values_src_cast(), values_src_cast(), dims);
	symphas::dft::fftw_execute(p_in_out);
}

//! Specialization of FourierGrid<S, T, D>::update().
template<>
inline void FourierGrid<scalar_t, complex_t, 2>::update()
{
	symphas::dft::arrange_fftw_stip<2>(values_src_cast(), values_src_cast(), dims);
	symphas::dft::fftw_execute(p_in_out);
}

//! Specialization of FourierGrid<S, T, D>::update().
template<>
inline void FourierGrid<scalar_t, complex_t, 3>::update()
{
	symphas::dft::arrange_fftw_stip<3>(values_src_cast(), values_src_cast(), dims);
	symphas::dft::fftw_execute(p_in_out);
}





