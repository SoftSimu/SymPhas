
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
 * PURPOSE: Manages a group of provisional systems, used by the phase field
 * model in order to manage the provisional variables.
 *
 * ***************************************************************************
 */

#pragma once

#include "spslibfftw.h"
#include "boundarysystem.h"


#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

 //! The default phase field system.
 /*!
  * The representation of a phase field system, storing the values of the
  * order parameter that is defined in a phase field problem.
  *
  * Unless explicitly specified, this phase field system type will be
  * used by the solver. It does not manage boundaries or other data.
  *
  * \tparam T The order parameter type.
  * \tparam D The order parameter dimension.
  */
template<typename T, size_t D>
using SolverSystem = System<T, D>;


template<typename T, size_t D>
struct SolverSystemFD : BoundarySystem<T, D>
{
	using BoundarySystem<T, D>::dims;

	BoundaryGrid<T, D> dframe;		// the working grid for the solver
	SolverSystemFD(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id = 0) :
		BoundarySystem<T, D>(tdata, vdata, bdata, id), dframe{ dims } {}
	SolverSystemFD() : BoundarySystem<T, D>(), dframe{ 0 } {}
};



//! The phase field system used by the spectral solver. 
/*!
 * The spectral solver requires the Fourier transforms of all variables
 * when computing the solution, so in addition to storing the real-space
 * provisional system, the Fourier space transform is also stored and updated.
 *
 * The phase field does not use or define boundary conditions, although
 * they are still supplied to the constructor in order to use the
 * workflow. 
 *
 * \tparam T The provisional system type, in real space.
 * \tparam D The provisional system dimension.
 */
template<typename T, size_t D>
struct SolverSystemSpectral : System<T, D> 
{
	using System<T, D>::System;
	using Grid<T, D>::values;

	//! Create the phase field data for the spectral solver implementation.
	SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id = 0) : System<T, D>(tdata, vdata, id) {}
	SolverSystemSpectral() : System<T, D>() {}
};

//! The phase field system used by the spectral solver. 
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a real-valued order parameter. It defines the Fourier
 * transform of the order parameter with dimensions that do not include
 * the duplicated half. 
 * 
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template<size_t D>
struct SolverSystemSpectral<scalar_t, D> : System<scalar_t, D>
{
	using System<scalar_t, D>::System;
	using Grid<scalar_t, D>::dims;
	using Grid<scalar_t, D>::values;
	
	len_type transformed_len;	//!< Length of the transformed array.
	complex_t* frame_t;			//!< The values of the transformed grid.
	complex_t* dframe;			//!< Accumulates the values of the solver computation.

	fftw_plan p;				//!< Fourier transform plan of the order parameter.


	//! Create the order parameter data used by the spectral solver. 
	/*!
	 * Create the order parameter data used by the spectral solver. The
	 * boundaries are provided but not used in the implementation, as they are
	 * imposed by the solver instead.
	 * 
	 * \param tdata The initial conditions data of the system.
	 * \param vdata The interval data of the system.
	 * \param id The ID value of the system.
	 */
	SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id);
	SolverSystemSpectral(SolverSystemSpectral<scalar_t, D> const& other);
	SolverSystemSpectral(SolverSystemSpectral<scalar_t, D>&& other) noexcept : SolverSystemSpectral()
	{
		swap(*this, other);
	}

	SolverSystemSpectral<scalar_t, D>& operator=(SolverSystemSpectral<scalar_t, D> other)
	{
		swap(*this, other);
		return *this;
	}

	//! Compute the Fourier transform and update the system.
	/*!
	 * Compute the Fourier transform and update the system. On some iterations,
	 * the complex transform is recomputed in order to eliminate numerical
	 * error from accumulating. The values computed
	 * by the solver are copied to the Fourier transformed data, and then 
	 * the Fourier transform is inverted to recover the real phase field values.
	 * 
	 * \param index The index of the solution.
	 */
	void update(iter_type index, double)
	{
		if (index % 100 == 0)
		{
			symphas::dft::arrange_fftw_stip<D>(values, reinterpret_cast<scalar_t*>(dframe), dims);
			symphas::dft::fftw_execute(p_to_t);
		}
		std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par,
#endif
			dframe, dframe + transformed_len, frame_t);

		symphas::dft::fftw_execute(p);
		symphas::dft::arrange_fftw_ipts<D>(reinterpret_cast<scalar_t*>(dframe), values, dims);
		grid::scale(System<scalar_t, D>::as_grid());
	}

	friend void swap(SolverSystemSpectral<scalar_t, D>& first, SolverSystemSpectral<scalar_t, D>& second)
	{
		using std::swap;

		swap(static_cast<System<scalar_t, D>&>(first), static_cast<System<scalar_t, D>&>(second));
		swap(first.transformed_len, second.transformed_len);
		swap(first.frame_t, second.frame_t);
		swap(first.dframe, second.dframe);
		swap(first.p, second.p);
		swap(first.p_to_t, second.p_to_t);
	}


	SolverSystemSpectral() : System<scalar_t, D>(), transformed_len{ 0 }, frame_t{ nullptr }, dframe{ nullptr }, p{ 0 }, p_to_t{ 0 } {}

	~SolverSystemSpectral();

protected:

	fftw_plan p_to_t;


};


//! The phase field system used by the spectral solver. 
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a complex-valued order parameter. 
 * 
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template<size_t D>
struct SolverSystemSpectral<complex_t, D> : System<complex_t, D>
{
	using System<complex_t, D>::System;
	using Grid<complex_t, D>::dims;
	using Grid<complex_t, D>::values;

	len_type transformed_len;	//!< Length of the transformed array.
	complex_t* frame_t;			//!< The values of the transformed grid.
	complex_t* dframe;			//!< Accumulates the values of the solver computation.
	fftw_plan p;				//!< Fourier transform plan of the order parameter.


	//! Create the order parameter data used by the spectral solver. 
	/*!
	 * Create the order parameter data used by the spectral solver. The
	 * boundaries are provided but not used in the implementation, as they are
	 * imposed by the solver instead.
	 *
	 * \param tdata The initial conditions data of the system.
	 * \param vdata The interval data of the system.
	 * \param id The ID value of the system.
	 */
	SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id);
	SolverSystemSpectral(SolverSystemSpectral<complex_t, D> const& other);
	SolverSystemSpectral(SolverSystemSpectral<complex_t, D>&& other) noexcept : SolverSystemSpectral()
	{
		swap(*this, other);
	}
	SolverSystemSpectral<complex_t, D>& operator=(SolverSystemSpectral<complex_t, D> other)
	{
		swap(*this, other);
		return *this;
	}

	//! Compute the Fourier transform and update the system.
	/*!
	 * Compute the Fourier transform and update the system. The values computed
	 * by the solver are copied to the Fourier transformed data, and then 
	 * the Fourier transform is inverted to recover the real phase field values.
	 *
	 * \param index The index of the solution.
	 */
	inline void update(iter_type, double)
	{

		std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par,
#endif
			dframe, dframe + transformed_len, frame_t);

		symphas::dft::fftw_execute(p);
		std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par,
#endif
			dframe, dframe + transformed_len, values);
		grid::scale(*this);
	}


	friend void swap(SolverSystemSpectral<complex_t, D>& first, SolverSystemSpectral<complex_t, D>& second)
	{
		using std::swap;

		swap(static_cast<System<complex_t, D>&>(first), static_cast<System<complex_t, D>&>(second));
		swap(first.transformed_len, second.transformed_len);
		swap(first.frame_t, second.frame_t);
		swap(first.dframe, second.dframe);
		swap(first.p, second.p);
		swap(first.p_to_t, second.p_to_t);
	}

	SolverSystemSpectral() : System<complex_t, D>(), transformed_len{ 0 }, frame_t{ nullptr }, dframe{ nullptr }, p{ 0 }, p_to_t{ 0 } {}

	~SolverSystemSpectral();

protected:

	fftw_plan p_to_t;


};

//! The phase field system used by the spectral solver. 
/*!
 * Implementation of the phase field system used for the spectral solver
 * based on a vector-valued order parameter. It defines the Fourier
 * transform of the order parameter with dimensions that do not include
 * the duplicated half.
 *
 * See SolverSystemSpectral<T, D>.
 *
 * \tparam D The provisional system dimension.
 */
template<size_t D>
struct SolverSystemSpectral<vector_t<D>, D> : System<vector_t<D>, D>
{
	using base_type = vector_t<D>;

	using System<base_type, D>::System;
	using Grid<base_type, D>::dims;
	using Grid<base_type, D>::axis;

	len_type transformed_len;				//!< Length of the transformed array.
	MultiBlock<D, complex_t> frame_t;		//!< The values of the transformed grid.
	MultiBlock<D, complex_t> dframe;		//!< Accumulates the values of the solver computation.

	fftw_plan p[D];							//!< Fourier transform plan of the order parameter.


	//! Create the order parameter data used by the spectral solver. 
	/*!
	 * Create the order parameter data used by the spectral solver. The
	 * boundaries are provided but not used in the implementation, as they are
	 * imposed by the solver instead.
	 *
	 * \param tdata The initial conditions data of the system.
	 * \param vdata The interval data of the system.
	 * \param id The ID value of the system.
	 */
	SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id);
	SolverSystemSpectral(SolverSystemSpectral<base_type, D> const& other);
	SolverSystemSpectral(SolverSystemSpectral<base_type, D>&& other) noexcept : SolverSystemSpectral()
	{
		swap(*this, other);
	}
	SolverSystemSpectral<base_type, D>& operator=(SolverSystemSpectral<base_type, D> other)
	{
		swap(*this, other);
		return *this;
	}


	//! Compute the Fourier transform and update the system.
	/*!
	 * Compute the Fourier transform and update the system. On some iterations,
	 * the complex transform is recomputed in order to eliminate numerical
	 * error from accumulating. The values computed
	 * by the solver are copied to the Fourier transformed data, and then
	 * the Fourier transform is inverted to recover the real phase field values.
	 *
	 * \param index The index of the solution.
	 */
	void update(iter_type index, double)
	{
		if (index % 100 == 0)
		{
			for (iter_type i = 0; i < D; ++i)
			{
				Axis ax = symphas::index_to_axis(i);
				symphas::dft::arrange_fftw_stip<D>(axis(ax), reinterpret_cast<scalar_t*>(dframe(ax)), dims);
				symphas::dft::fftw_execute(p_to_t[i]);
			}
		}
		for (iter_type i = 0; i < D; ++i)
		{
			Axis ax = symphas::index_to_axis(i);

			std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par,
#endif
				dframe(ax), dframe(ax) + transformed_len, frame_t(ax));

			symphas::dft::fftw_execute(p[i]);
			symphas::dft::arrange_fftw_ipts<D>(reinterpret_cast<scalar_t*>(dframe(ax)), axis(ax), dims);
		}
		grid::scale(System<base_type, D>::as_grid());
	}

	friend void swap(SolverSystemSpectral<base_type, D>& first, SolverSystemSpectral<base_type, D>& second)
	{
		using std::swap;

		swap(static_cast<System<base_type, D>&>(first), static_cast<System<base_type, D>&>(second));
		swap(first.transformed_len, second.transformed_len);
		swap(first.frame_t, second.frame_t);
		swap(first.dframe, second.dframe);
		swap(first.p, second.p);
		swap(first.p_to_t, second.p_to_t);
	}

	SolverSystemSpectral() : System<base_type, D>(), transformed_len{ 0 }, frame_t{ 0 }, dframe{ 0 }, p{ 0 }, p_to_t{ 0 } {}

	~SolverSystemSpectral();

protected:

	fftw_plan p_to_t[D];


};





using symphas::dft::new_fftw_plan;

template<>
inline SolverSystemSpectral<scalar_t, 1>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<scalar_t, 1>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 1>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<1, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<1, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<1>(Grid<scalar_t, 1>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<scalar_t, 2>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<scalar_t, 2>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 2>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<2, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<2, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<2>(Grid<scalar_t, 2>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<scalar_t, 3>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<scalar_t, 3>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 3>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<3, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<3, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<3>(Grid<scalar_t, 3>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 1>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<complex_t, 1>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<complex_t, 1>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 1>::values, Grid<complex_t, 1>::values + Grid<complex_t, 1>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 2>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<complex_t, 2>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<complex_t, 2>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 2>::values, Grid<complex_t, 2>::values + Grid<complex_t, 2>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 3>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<complex_t, 3>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<complex_t, 3>(dims) },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 3>::values, Grid<complex_t, 3>::values + Grid<complex_t, 3>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<any_vector_t<scalar_t, 1>, 1>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 1>(dims) },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<1, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims) },
	p_to_t{ new_fftw_plan<1, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims) }
{
	symphas::dft::arrange_fftw_stip<1>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
}

template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<any_vector_t<scalar_t, 2>, 2>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 2>(dims) },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims) },
	p_to_t{ new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims) }
{
	symphas::dft::arrange_fftw_stip<2>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::arrange_fftw_stip<2>(axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
	symphas::dft::fftw_execute(p_to_t[1]);
}

template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3>::SolverSystemSpectral(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id) :
	System<any_vector_t<scalar_t, 3>, 3>(tdata, vdata, id),
	transformed_len{ symphas::dft::length<scalar_t, 3>(dims) },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims), new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Z), dframe(Axis::Z), dims) },
	p_to_t{ new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims), new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Z), dframe(Axis::Z), dims) }
{
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::Z), reinterpret_cast<scalar_t*>(dframe(Axis::Z)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
	symphas::dft::fftw_execute(p_to_t[1]);
	symphas::dft::fftw_execute(p_to_t[2]);
}





template<>
inline SolverSystemSpectral<scalar_t, 1>::SolverSystemSpectral(SolverSystemSpectral<scalar_t, 1> const& other) :
	System<scalar_t, 1>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<1, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<1, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<1>(Grid<scalar_t, 1>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<scalar_t, 2>::SolverSystemSpectral(SolverSystemSpectral<scalar_t, 2> const& other) :
	System<scalar_t, 2>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<2, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<2, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<2>(Grid<scalar_t, 2>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<scalar_t, 3>::SolverSystemSpectral(SolverSystemSpectral<scalar_t, 3> const& other) :
	System<scalar_t, 3>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<3, complex_t, scalar_t>{}(dframe, dframe, dims) },
	p_to_t{ new_fftw_plan<3, scalar_t, complex_t>{}(dframe, dframe, dims) }
{
	symphas::dft::arrange_fftw_stip<3>(Grid<scalar_t, 3>::values, reinterpret_cast<scalar_t*>(dframe), dims);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 1>::SolverSystemSpectral(SolverSystemSpectral<complex_t, 1> const& other) :
	System<complex_t, 1>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<1, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 1>::values, Grid<complex_t, 1>::values + Grid<complex_t, 1>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 2>::SolverSystemSpectral(SolverSystemSpectral<complex_t, 2> const& other) :
	System<complex_t, 2>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<2, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 2>::values, Grid<complex_t, 2>::values + Grid<complex_t, 2>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}

template<>
inline SolverSystemSpectral<complex_t, 3>::SolverSystemSpectral(SolverSystemSpectral<complex_t, 3> const& other) :
	System<complex_t, 3>(other),
	transformed_len{ other.transformed_len },
	frame_t{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	dframe{ reinterpret_cast<complex_t*>(symphas::dft::fftw_alloc_complex(transformed_len)) },
	p{ new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false, true) },
	p_to_t{ new_fftw_plan<3, complex_t, complex_t>{}(dframe, dframe, dims, false, false) }
{
	std::copy(Grid<complex_t, 3>::values, Grid<complex_t, 3>::values + Grid<complex_t, 3>::len, dframe);
	symphas::dft::fftw_execute(p_to_t);
}


template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1>::SolverSystemSpectral(SolverSystemSpectral<any_vector_t<scalar_t, 1>, 1> const& other) :
	System<any_vector_t<scalar_t, 1>, 1>(other),
	transformed_len{ other.transformed_len },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<1, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims) },
	p_to_t{ new_fftw_plan<1, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims)}
{
	symphas::dft::arrange_fftw_stip<1>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
}

template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2>::SolverSystemSpectral(SolverSystemSpectral<any_vector_t<scalar_t, 2>, 2> const& other) :
	System<any_vector_t<scalar_t, 2>, 2>(other),
	transformed_len{ other.transformed_len },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<2, complex_t, scalar_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims) },
	p_to_t{ new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<2, scalar_t, complex_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims) }
{
	symphas::dft::arrange_fftw_stip<2>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::arrange_fftw_stip<2>(axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
	symphas::dft::fftw_execute(p_to_t[1]);
}

template<>
inline SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3>::SolverSystemSpectral(SolverSystemSpectral<any_vector_t<scalar_t, 3>, 3> const& other) :
	System<any_vector_t<scalar_t, 3>, 3>(other),
	transformed_len{ other.transformed_len },
	frame_t{ transformed_len },
	dframe{ transformed_len },
	p{ new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims), new_fftw_plan<3, complex_t, scalar_t>{}(dframe(Axis::Z), dframe(Axis::Z), dims) },
	p_to_t{ new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::X), dframe(Axis::X), dims), new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Y), dframe(Axis::Y), dims), new_fftw_plan<3, scalar_t, complex_t>{}(dframe(Axis::Z), dframe(Axis::Z), dims) }
{
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::X), reinterpret_cast<scalar_t*>(dframe(Axis::X)), dims);
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::Y), reinterpret_cast<scalar_t*>(dframe(Axis::Y)), dims);
	symphas::dft::arrange_fftw_stip<3>(axis(Axis::Z), reinterpret_cast<scalar_t*>(dframe(Axis::Z)), dims);
	symphas::dft::fftw_execute(p_to_t[0]);
	symphas::dft::fftw_execute(p_to_t[1]);
	symphas::dft::fftw_execute(p_to_t[2]);
}





template<size_t D>
inline SolverSystemSpectral<scalar_t, D>::~SolverSystemSpectral()
{
	symphas::dft::fftw_destroy_plan(p);
	symphas::dft::fftw_destroy_plan(p_to_t);
	symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(dframe));
	symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(frame_t));
}

template<size_t D>
inline SolverSystemSpectral<complex_t, D>::~SolverSystemSpectral()
{
	symphas::dft::fftw_destroy_plan(p);
	symphas::dft::fftw_destroy_plan(p_to_t);
	symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(dframe));
	symphas::dft::fftw_free(reinterpret_cast<fftw_complex*&>(frame_t));
}

template<size_t D>
inline SolverSystemSpectral<vector_t<D>, D>::~SolverSystemSpectral()
{
	for (iter_type i = 0; i < D; ++i)
	{
		symphas::dft::fftw_destroy_plan(p[i]);
		symphas::dft::fftw_destroy_plan(p_to_t[i]);
	}
}


DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFD<T, D>), (BoundarySystem<T, D>))
DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemSpectral<T, D>), (System<T, D>))




