
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
 * PURPOSE: Defines a provisional system. A provisional system stores values
 * of the provisional variables which can constitute part of the model
 * definition. Provisional variables can act as virtual variables; they
 * are computed at each step and can be used in the equations of motion, but
 * are not physical variables like the order parameters.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundarysystem.h"
#include "spslibfftw.h"

//! The default provisional system.
/*!
 * The representation of a provisional system, storing the values of the
 * provisional variable that is defined in a phase field problem.
 * 
 * Unless explicitly specified, this provisional system type will be
 * used by the solver. It does not manage boundaries or other data.
 * 
 * \tparam T The provisional system type.
 * \tparam D The provisional system dimension.
 */
template<typename T, size_t D>
using ProvisionalSystem = System<T, D>;

//! The provisional system used by the spectral solver. 
/*!
 * The spectral solver requires the Fourier transforms of all variables
 * when computing the solution, so in addition to storing the real-space
 * provisional system, the Fourier space transform is also stored and updated.
 * 
 * The provisional variable does not use boundary conditions, although
 * they are still supplied to the constructor in order to use the 
 * workflow. 
 * 
 * \tparam T The provisional system type, in real space.
 * \tparam D The provisional system dimension.
 */
template<typename T, size_t D>
struct ProvisionalSystemSpectral : Grid<T, D>
{
	using Grid<T, D>::dims;

	fftw_complex* frame_t;
	fftw_plan p_src_out;

	//! Construct the grid data of a provisional variable.
	/*!
	 * Construct the grid data of a provisional variable. The Fourier space
	 * transform of the variable is also initialized.
	 * 
	 * \param vdata Interval data about the grid system.
	 */
	ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const&);

	ProvisionalSystemSpectral(ProvisionalSystemSpectral<T, D> const& other);
	ProvisionalSystemSpectral(ProvisionalSystemSpectral<T, D>&& other) noexcept;
	ProvisionalSystemSpectral<T, D>& operator=(ProvisionalSystemSpectral<T, D> other);

	friend void swap(ProvisionalSystemSpectral<T, D>& first, ProvisionalSystemSpectral<T, D>& second)
	{
		using std::swap;

		swap(static_cast<Grid<T, D>&>(first), static_cast<Grid<T, D>&>(second));
		swap(first.frame_t, second.frame_t);
		swap(first.p_src_out, second.p_src_out);
	}

	~ProvisionalSystemSpectral();

	//! Executes the FFTW plan to transform the provisional variable.
	/*!
	 * Executes the FFTW plan to transform the provisional variable.
	 * No parameters are used in the update.
	 */
	void update(iter_type, double)
	{
		symphas::dft::fftw_execute(p_src_out);
	}
	

protected:

	ProvisionalSystemSpectral() : Grid<T, D>(), frame_t{ nullptr }, p_src_out{ nullptr } {}

};

//! The provisional system used by the spectral solver. 
/*!
 * Specialization when the provisional variable is a vector type.
 *
 * The spectral solver requires the Fourier transforms of all variables
 * when computing the solution, so in addition to storing the real-space
 * provisional system, the Fourier space transform is also stored and updated.
 *
 * \tparam D The provisional system dimension.
 */
template<size_t D>
struct ProvisionalSystemSpectral<vector_t<D>, D> : Grid<vector_t<D>, D>
{
	//! Construct the data of a vector-valued provisional variable.
	ProvisionalSystemSpectral(
		const symphas::interval_data_type vdata,
		const symphas::b_data_type) : 
		Grid<vector_t<D>, D>{ grid::construct<Grid, vector_t<D>, D>(vdata) } {}

protected:

	ProvisionalSystemSpectral() : Grid<vector_t<D>, D>() {}

};




//! The initial conditions of a provisional system.
/*!
 * Initial conditions for provisional variables are required with solvers
 * that use finite difference grids with boundaries. In that case, the
 * definition of the provisional object means that the initial conditions need
 * to be passed to the constructor of its parent object, even though
 * initial conditions have no meaning in the context of provisional variables.
 * 
 * The initial conditions of a provisional system are undefined and set to
 * `NONE` value for both parameters, thereby skipping any initial generation.
 */
inline symphas::init_data_type provisional_init{
	Inside::NONE,
	symphas::build_intag(InsideTag::NONE),
	symphas::init_data_parameters{} };


//! The representation of a provisional system.
/*!
 * The representation of a provisional system, storing the values of the
 * provisional variable that is defined in a phase field problem.
 * 
 * This provisional system implementation is used when the solver implements
 * a finite difference stencil to approximate all derivatives, thereby
 * requiring boundaries.
 * The boundaries are updated after the provisional variables are computed
 * from the corresponding equations. The boundary data is typically the same
 * as the first field that is given, but in practise it does not matter
 * for the numerical results unless finite difference approximations are
 * applied extensively to the provisional variables.
 *
 * \tparam T The provisional system type.
 * \tparam D The dimension of the provisional system, corresponding to
 * the dimension of the phase field problem.
 */
template<typename T, size_t D>
struct ProvisionalSystemFD : BoundarySystem<T, D>
{
	//! Create a provisional system.
	/*!
	 * Create a provisional system using the given intervals and boundary data.
	 */
	ProvisionalSystemFD(
		const symphas::interval_data_type vdata, 
		const symphas::b_data_type bdata) : 
		BoundarySystem<T, D>(provisional_init, vdata, bdata) {}
};




using symphas::dft::new_fftw_plan;

template<>
inline ProvisionalSystemSpectral<scalar_t, 1>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<scalar_t, 1>(grid::construct<::Grid, scalar_t, 1>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 1>(dims)) },
	p_src_out{ new_fftw_plan<1, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<scalar_t, 2>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<scalar_t, 2>(grid::construct<::Grid, scalar_t, 2>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 2>(dims)) },
	p_src_out{ new_fftw_plan<2, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<scalar_t, 3>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<scalar_t, 3>(grid::construct<::Grid, scalar_t, 3>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 3>(dims)) },
	p_src_out{ new_fftw_plan<3, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 1>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<complex_t, 1>(grid::construct<::Grid, complex_t, 1>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 1>(dims)) },
	p_src_out{ new_fftw_plan<1, complex_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 2>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<complex_t, 2>(grid::construct<::Grid, complex_t, 2>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 2>(dims)) },
	p_src_out{ new_fftw_plan<2, complex_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 3>::ProvisionalSystemSpectral(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata) :
	Grid<complex_t, 3>(grid::construct<::Grid, complex_t, 3>(vdata)),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 3>(dims)) },
	p_src_out{ new_fftw_plan<3, complex_t, complex_t>{}(values, frame_t, dims) }
{}





template<>
inline ProvisionalSystemSpectral<scalar_t, 1>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<scalar_t, 1> const& other) :
	Grid<scalar_t, 1>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 1>(dims)) },
	p_src_out{ new_fftw_plan<1, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<scalar_t, 2>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<scalar_t, 2> const& other) :
	Grid<scalar_t, 2>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 2>(dims)) },
	p_src_out{ new_fftw_plan<2, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<scalar_t, 3>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<scalar_t, 3> const& other) :
	Grid<scalar_t, 3>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<scalar_t, 3>(dims)) },
	p_src_out{ new_fftw_plan<3, scalar_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 1>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<complex_t, 1> const& other) :
	Grid<complex_t, 1>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 1>(dims)) },
	p_src_out{ new_fftw_plan<1, complex_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 2>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<complex_t, 2> const& other) :
	Grid<complex_t, 2>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 2>(dims)) },
	p_src_out{ new_fftw_plan<2, complex_t, complex_t>{}(values, frame_t, dims) }
{}

template<>
inline ProvisionalSystemSpectral<complex_t, 3>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<complex_t, 3> const& other) :
	Grid<complex_t, 3>(other),
	frame_t{ symphas::dft::fftw_alloc_complex(symphas::dft::length<complex_t, 3>(dims)) },
	p_src_out{ new_fftw_plan<3, complex_t, complex_t>{}(values, frame_t, dims) }
{}



template<typename T, size_t D>
ProvisionalSystemSpectral<T, D>::ProvisionalSystemSpectral(ProvisionalSystemSpectral<T, D>&& other) noexcept : ProvisionalSystemSpectral()
{
	swap(*this, other);
}

template<typename T, size_t D>
ProvisionalSystemSpectral<T, D>& ProvisionalSystemSpectral<T, D>::operator=(ProvisionalSystemSpectral<T, D> other)
{
	swap(*this, other);
	return *this;
}



template<typename T, size_t D>
inline ProvisionalSystemSpectral<T, D>::~ProvisionalSystemSpectral()
{
	symphas::dft::fftw_free(frame_t);
	symphas::dft::fftw_destroy_plan(p_src_out);
}