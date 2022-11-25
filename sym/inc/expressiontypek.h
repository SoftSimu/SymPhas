
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
 * MODULE:  lib
 * PURPOSE: Defines the wavevector object, used for example, to construct
 * a Gaussian smoothing kernel.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"


//! \cond

#ifdef LATEX_PLOT

#define SYEX_K_COMPONENT_OP_STR "k"
#define SYEX_K_AXIS_OP_STR "\\vec{" SYEX_K_COMPONENT_OP_STR ""}"
#define SYEX_K_OP_STR "|" SYEX_K_AXIS_OP_STR "|"
#define SUBSCRIPT_AXIS "{_%c}"

#else

#define SYEX_K_COMPONENT_OP_STR "k"
#define SYEX_K_AXIS_OP_STR SYEX_K_COMPONENT_OP_STR
#define SYEX_K_COMPONENT_OP_STR "k"
#define SYEX_K_OP_STR "|" SYEX_K_AXIS_OP_STR "|"
#define SUBSCRIPT_AXIS "%c"

#endif

//! The display format for printing the wavenumber.
#define SYEX_K_EVEN_FMT SYEX_K_OP_STR SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B
#define SYEX_K_ODD_FMT SYEX_K_AXIS_OP_STR SYEX_K_OP_STR SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B
#define SYEX_K_AXIS_FMT SYEX_K_COMPONENT_OP_STR SUBSCRIPT_AXIS SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B
#define SYEX_K_AXIS_0_FMT SYEX_K_COMPONENT_OP_STR SUBSCRIPT_AXIS
#define SYEX_K_COMPONENT_FMT SYEX_K_COMPONENT_OP_STR SUBSCRIPT_AXIS SYEX_K_OP_STR SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B
//! The latex display format for printing the wavenumber.
#define SYEX_K_EVEN_FMT_LEN (STR_ARR_LEN(SYEX_K_OP_STR SYEX_POW_SEP_A SYEX_POW_SEP_B) - 1)
#define SYEX_K_ODD_FMT_LEN (STR_ARR_LEN(SYEX_K_AXIS_OP_STR SYEX_K_OP_STR SYEX_POW_SEP_A SYEX_POW_SEP_B) - 1)
#define SYEX_K_AXIS_FMT_LEN (STR_ARR_LEN(SYEX_K_AXIS_OP_STR SYEX_POW_SEP_A SYEX_POW_SEP_B))
#define SYEX_K_AXIS_0_FMT_LEN (STR_ARR_LEN(SYEX_K_AXIS_OP_STR))
#define SYEX_K_COMPONENT_FMT_LEN (STR_ARR_LEN(SYEX_K_AXIS_OP_STR SYEX_K_OP_STR SYEX_POW_SEP_A SYEX_POW_SEP_B))


//! \endcond

// **************************************************************************************

namespace expr
{



	template<size_t O, size_t D>
	using wave_vector_grid = std::conditional_t<(O % 2 == 0), Grid<scalar_t, D>, Grid<cvector_t<D>, D>>;
	template<size_t O, size_t D>
	using wave_vector_axis_grid = std::conditional_t<(O % 2 == 0), Grid<scalar_t, D>, Grid<complex_t, D>>;

	//! Manages the filling of values of a wavenumber field.
	/*!
	 * Computes the values of the wavevector field, equivalent to the Fourier
	 * space, and stores it in an array.
	 * 
	 * \tparam D The dimension of the wavevector field.
	 */
	template<size_t D>
	struct k_field;

	//! Specialization based on k_field.
	template<>
	struct k_field<1>
	{
		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<size_t O>
		static void fill(wave_vector_grid<O, 1> &into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			scalar_t dk_i = 2 * symphas::PI / (*h * L);

			std::complex Ii(0.0, 1.0);

			for (iter_type i = 0; i < L; ++i)
			{
				scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
				if (kx == 0)
				{
					kx = std::pow(symphas::EPS, 1.0 / (double)(O));
				}

				cvector_t<1> k{ Ii * kx };
				scalar_t k2 = (k * k).real();

				if constexpr (O % 2 == 0)
				{
					into[i] = std::pow(k2, O / 2);
				}
				else
				{
					into[i] = k * std::pow(k2, O / 2);
				}
			}

		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill(wave_vector_axis_grid<O, 1>& into, const len_type* dims, const double* h)
		{
			if constexpr (ax == Axis::X)
			{
				std::complex Ii(0.0, 1.0);
				len_type L = dims[0];

				scalar_t dk_i = 2 * symphas::PI / (h[0] * L);

				iter_type ii = 0;
				for (iter_type i = 0; i < L; ++i)
				{
					scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
					if (kx == 0)
					{
						kx = std::pow(symphas::EPS, 1.0 / (double)(O));
					}

					complex_t kc = Ii * ((ax == Axis::X) ? kx : 0);
					cvector_t<1> k{ Ii * kx };
					scalar_t k2 = (k * k).real();

					if constexpr (O % 2 == 0)
					{
						into[i] = std::pow(k2, O / 2);
					}
					else
					{
						into[i] = kc * std::pow(k2, O / 2);
					}
				}
			}
		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill_axis(wave_vector_axis_grid<O, 1>& into, const len_type* dims, const double* h)
		{
			if constexpr (ax == Axis::X)
			{
				std::complex Ii(0.0, 1.0);
				len_type L = dims[0];

				scalar_t dk_i = 2 * symphas::PI / (h[0] * L);

				for (iter_type i = 0; i < L; ++i)
				{
					scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
					if (kx == 0)
					{
						kx = std::pow(symphas::EPS, 1.0 / (double)(O));
					}

					complex_t k = Ii * ((ax == Axis::X) ? kx : 0);
					scalar_t k2 = (k * k).real();

					if constexpr (O % 2 == 0)
					{
						into[i] = std::pow(k2, O / 2);
					}
					else
					{
						into[i] = k * std::pow(k2, O / 2);
					}
				}
			}
		}
	};

	//! Specialization based on k_field.
	template<>
	struct k_field<2>
	{
		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<size_t O>
		static void fill(wave_vector_grid<O, 2>& into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			len_type M = dims[1];

			scalar_t dk_i = 2 * symphas::PI / (h[0] * L);
			scalar_t dk_j = 2 * symphas::PI / (h[1] * M);

			std::complex Ii(0.0, 1.0);

			iter_type ii = 0;
			for (iter_type j = 0; j < M; ++j)
			{
				scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
				if (ky == 0)
				{
					ky = std::pow(symphas::EPS, 1.0 / (double)(O));
				}
				for (iter_type i = 0; i < L; ++i)
				{
					scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
					if (kx == 0)
					{
						kx = std::pow(symphas::EPS, 1.0 / (double)(O));
					}
					cvector_t<2> k{ Ii * kx, Ii * ky };
					scalar_t k2 = (k * k).real();

					if constexpr (O % 2 == 0)
					{
						into[ii++] = std::pow(k2, O / 2);
					}
					else
					{
						into[ii++] = k * std::pow(k2, O / 2);
					}
				}
			}
		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill(wave_vector_axis_grid<O, 2> &into, const len_type* dims, const double* h)
		{
			if constexpr (ax == Axis::X || ax == Axis::Y)
			{
				std::complex Ii(0.0, 1.0);

				len_type L = dims[0];
				len_type M = dims[1];

				scalar_t dk_i = 2 * symphas::PI / (h[0] * L);
				scalar_t dk_j = 2 * symphas::PI / (h[1] * M);

				iter_type ii = 0;
				for (iter_type j = 0; j < M; ++j)
				{
					scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					if (ky == 0)
					{
						ky = std::pow(symphas::EPS, 1.0 / (double)(O));
					}

					for (iter_type i = 0; i < L; ++i)
					{
						scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						if (kx == 0)
						{
							kx = std::pow(symphas::EPS, 1.0 / (double)(O));
						}

						complex_t kc = Ii * ((ax == Axis::X) ? kx : (ax == Axis::Y) ? ky : 0);
						cvector_t<2> k{ Ii * kx, Ii * ky };
						scalar_t k2 = (k * k).real();

						if constexpr (O % 2 == 0)
						{
							into[ii++] = std::pow(k2, O / 2);
						}
						else
						{
							into[ii++] = kc * std::pow(k2, O / 2);
						}
					}
				}
			}
		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill_axis(wave_vector_axis_grid<O, 2>& into, const len_type* dims, const double* h)
		{
			if constexpr (ax == Axis::X || ax == Axis::Y)
			{
				std::complex Ii(0.0, 1.0);

				len_type L = dims[0];
				len_type M = dims[1];

				scalar_t dk_i = 2 * symphas::PI / (h[0] * L);
				scalar_t dk_j = 2 * symphas::PI / (h[1] * M);

				iter_type ii = 0;
				for (iter_type j = 0; j < M; ++j)
				{
					scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					if (ky == 0)
					{
						ky = std::pow(symphas::EPS, 1.0 / (double)(O));
					}

					for (iter_type i = 0; i < L; ++i)
					{
						scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						if (kx == 0)
						{
							kx = std::pow(symphas::EPS, 1.0 / (double)(O));
						}

						complex_t k = Ii * ((ax == Axis::X) ? kx : (ax == Axis::Y) ? ky : 0);
						scalar_t k2 = (k * k).real();

						if constexpr (O % 2 == 0)
						{
							into[ii++] = std::pow(k2, O / 2);
						}
						else
						{
							into[ii++] = k * std::pow(k2, O / 2);
						}
					}
				}
			}
		}

	};

	//! Specialization based on k_field.
	template<>
	struct k_field<3>
	{
		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<size_t O>
		static void fill(wave_vector_grid<O, 3>& into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			len_type M = dims[1];
			len_type N = dims[2];

			scalar_t
				dk_i = 2 * symphas::PI / (h[0] * L),
				dk_j = 2 * symphas::PI / (h[1] * M),
				dk_k = 2 * symphas::PI / (h[2] * N);

			std::complex Ii(0.0, 1.0);

			iter_type ii = 0;
			for (iter_type k = 0; k < N; ++k)
			{
				scalar_t kz = (k < N / 2) ? k * dk_k : (k - N) * dk_k;
				if (kz == 0)
				{
					kz = std::pow(symphas::EPS, 1.0 / (double)(O));
				}

				for (iter_type j = 0; j < M; ++j)
				{
					scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					if (ky == 0)
					{
						ky = std::pow(symphas::EPS, 1.0 / (double)(O));
					}

					for (iter_type i = 0; i < L; ++i)
					{
						scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						if (kx == 0)
						{
							kx = std::pow(symphas::EPS, 1.0 / (double)(O));
						}
						
						cvector_t<3> k{ Ii * kx, Ii * ky, Ii * kz };
						scalar_t k2 = (k * k).real();

						if constexpr (O % 2 == 0)
						{
							into[ii++] = std::pow(k2, O / 2);
						}
						else
						{
							into[ii++] = k * std::pow(k2, O / 2);
						}
					}
				}
			}
		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill(wave_vector_axis_grid<O, 3> &into, const len_type* dims, const double* h)
		{
			std::complex Ii(0.0, 1.0);

			len_type L = dims[0];
			len_type M = dims[1];
			len_type N = dims[2];

			scalar_t
				dk_i = 2 * symphas::PI / (h[0] * L),
				dk_j = 2 * symphas::PI / (h[1] * M),
				dk_k = 2 * symphas::PI / (h[2] * N);

			iter_type ii = 0;
			for (iter_type k = 0; k < N; ++k)
			{
				scalar_t kz = (k < N / 2) ? k * dk_k : (k - N) * dk_k;
				if (kz == 0)
				{
					kz = std::pow(symphas::EPS, 1.0 / (double)(O));
				}
				for (iter_type j = 0; j < M; ++j)
				{
					scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					if (ky == 0)
					{
						ky = std::pow(symphas::EPS, 1.0 / (double)(O));
					}
					for (iter_type i = 0; i < L; ++i)
					{
						scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						if (kx == 0)
						{
							kx = std::pow(symphas::EPS, 1.0 / (double)(O));
						}

						complex_t kc = Ii * ((ax == Axis::X) ? kx : (ax == Axis::Y) ? ky : (ax == Axis::Z) ? kz : 0);
						cvector_t<3> k{ Ii * kx, Ii * ky, Ii * kz };
						scalar_t k2 = (k * k).real();

						if constexpr (O % 2 == 0)
						{
							into[ii++] = std::pow(k2, O / 2);
						}
						else
						{
							into[ii++] = kc * std::pow(k2, O / 2);
						}
					}
				}
			}
		}

		//! Fills the values of the wavenumber field of the prescribed dimension.
		/*!
		 * \param into The data array of the wavenumber.
		 * \param dims The dimensions of the wavenumber.
		 * \param h The uniform grid spacing of the wavenumber.
		 *
		 * \tparam O The exponential order (power) of the wavenumber.
		 */
		template<Axis ax, size_t O>
		static void fill_axis(wave_vector_axis_grid<O, 3>& into, const len_type* dims, const double* h)
		{
			std::complex Ii(0.0, 1.0);

			len_type L = dims[0];
			len_type M = dims[1];
			len_type N = dims[2];

			scalar_t
				dk_i = 2 * symphas::PI / (h[0] * L),
				dk_j = 2 * symphas::PI / (h[1] * M),
				dk_k = 2 * symphas::PI / (h[2] * N);

			iter_type ii = 0;
			for (iter_type k = 0; k < N; ++k)
			{
				scalar_t kz = (k < N / 2) ? k * dk_k : (k - N) * dk_k;
				if (kz == 0)
				{
					kz = std::pow(symphas::EPS, 1.0 / (double)(O));
				}
				for (iter_type j = 0; j < M; ++j)
				{
					scalar_t ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					if (ky == 0)
					{
						ky = std::pow(symphas::EPS, 1.0 / (double)(O));
					}
					for (iter_type i = 0; i < L; ++i)
					{
						scalar_t kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						if (kx == 0)
						{
							kx = std::pow(symphas::EPS, 1.0 / (double)(O));
						}

						complex_t k = Ii * ((ax == Axis::X) ? kx : (ax == Axis::Y) ? ky : (ax == Axis::Z) ? kz : 0);
						scalar_t k2 = (k * k).real();

						if constexpr (O % 2 == 0)
						{
							into[ii++] = std::pow(k2, O / 2);
						}
						else
						{
							into[ii++] = k * std::pow(k2, O / 2);
						}
					}
				}
			}
		}
	};
}

namespace symphas::internal
{
	//! Predicate indicating whether the given type is a wavenumber field.
	template<typename E>
	struct is_K_type;

}

//! The grid representing the wavenumber data.
/*!
 * \tparam T The type of the wavenumber values, typically ::scalar_t.
 * \tparam D The dimension of the wavenumber field.
 */
template<size_t O, size_t D>
struct WaveVectorData;

//! The grid representing the wavenumber data.
/*!
 * \tparam T The type of the wavenumber values, typically ::scalar_t.
 * \tparam D The dimension of the wavenumber field.
 */
template<Axis ax, size_t O, size_t D>
struct WaveVectorDataComponent;

//! The grid representing the wavenumber data for one axis.
/*!
 * \tparam T The type of the wavenumber values, typically ::scalar_t.
 * \tparam D The dimension of the wavenumber field.
 */
template<Axis ax, size_t O, size_t D>
struct WaveVectorDataAxis;



//! Enclosing template of the wavenumber field.
/*!
 * \tparam O The exponential order (power) of the wavenumber.
 */
template<size_t D>
struct K
{
	K() : h{ 0 } {}


	K(double const* h) : h{ 0 }
	{
		if (h != nullptr)
		{
			std::copy(h, h + D, this->h);
		}
	}

	const double* widths() const
	{
		return h;
	}

	double width(int i) const
	{
		return h[i];
	}

	double h[D];

};

template<size_t O, size_t D>
struct WaveVectorData : expr::wave_vector_grid<O, D>, K<D>
{
	using parent_type = expr::wave_vector_grid<O, D>;

	using parent_type::dims;
	using parent_type::len;
	using K<D>::h;

	WaveVectorData(const len_type* dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill<O>(cast(), dims, h);
	}

	WaveVectorData(std::initializer_list<len_type> dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill<O>(cast(), dims, h);
	}

	template<size_t O2>
	WaveVectorData(WaveVectorData<O2, D> const& other) : parent_type(other.dims), K<D>(other.h)
	{
		expr::k_field<D>::template fill<O>(cast(), dims, h);
	}

	template<size_t O2>
	WaveVectorData(WaveVectorData<O2, D>&& other) : WaveVectorData()
	{
		swap(*this, other);
	}

	WaveVectorData<O, D>& operator=(WaveVectorData<O, D> other)
	{
		swap(*this, other);
		return *this;
	}

	parent_type& cast()
	{
		return *static_cast<parent_type*>(this);
	}

	parent_type const& cast() const
	{
		return *static_cast<parent_type const*>(this);
	}

	friend void swap(WaveVectorData<O, D>& first, WaveVectorData<O, D>& second)
	{
		using std::swap;
		swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
		swap(static_cast<K<D>&>(first), static_cast<K<D>&>(second));
	}


protected:

	WaveVectorData() : parent_type(), K<D>() {}
};

template<Axis ax, size_t O, size_t D>
struct WaveVectorDataComponent : expr::wave_vector_axis_grid<O, D>, K<D>
{
	using parent_type = expr::wave_vector_axis_grid<O, D>;

	using parent_type::dims;
	using parent_type::len;
	using K<D>::h;

	WaveVectorDataComponent(const len_type* dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill<ax, O>(cast(), dims, h);
	}

	WaveVectorDataComponent(std::initializer_list<len_type> dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill<ax, O>(cast(), dims, h);
	}

	template<Axis ax2, size_t O2>
	WaveVectorDataComponent(WaveVectorDataComponent<ax2, O2, D> const& other) : parent_type(other.dims), K<D>(other.h)
	{
		expr::k_field<D>::template fill<ax, O>(cast(), dims, h);
	}

	template<Axis ax2, size_t O2>
	WaveVectorDataComponent(WaveVectorDataComponent<ax2, O2, D>&& other) : WaveVectorDataComponent()
	{
		swap(*this, other);
	}

	WaveVectorDataComponent<ax, O, D>& operator=(WaveVectorDataComponent<ax, O, D> other)
	{
		swap(*this, other);
		return *this;
	}

	parent_type& cast()
	{
		return *static_cast<parent_type*>(this);
	}

	parent_type const& cast() const
	{
		return *static_cast<parent_type const*>(this);
	}

	friend void swap(WaveVectorDataComponent<ax, O, D>& first, WaveVectorDataComponent<ax, O, D>& second)
	{
		using std::swap;
		swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
		swap(static_cast<K<D>&>(first), static_cast<K<D>&>(second));
	}

protected:

	WaveVectorDataComponent() : parent_type(), K<D>{} {}
};

template<Axis ax, size_t O, size_t D>
struct WaveVectorDataAxis : expr::wave_vector_axis_grid<O, D>, K<D>
{
	using parent_type = expr::wave_vector_axis_grid<O, D>;

	//using Grid<T, D>::values;
	using parent_type::dims;
	using parent_type::len;
	using K<D>::h;

	WaveVectorDataAxis(const len_type* dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill_axis<ax, O>(cast(), dims, h);
	}

	WaveVectorDataAxis(std::initializer_list<len_type> dimensions, double const* h) : parent_type(dimensions), K<D>(h)
	{
		expr::k_field<D>::template fill_axis<ax, O>(cast(), dims, h);
	}

	template<Axis ax2, size_t O2>
	WaveVectorDataAxis(WaveVectorDataAxis<ax2, O2, D> const& other) : parent_type(other.dims), K<D>(other.h)
	{
		expr::k_field<D>::template fill_axis<ax, O>(cast(), dims, h);
	}

	template<Axis ax2, size_t O2>
	WaveVectorDataAxis(WaveVectorDataAxis<ax2, O2, D>&& other) : WaveVectorDataAxis()
	{
		swap(*this, other);
	}

	WaveVectorDataAxis<ax, O, D>& operator=(WaveVectorDataAxis<ax, O, D> other)
	{
		swap(*this, other);
		return *this;
	}

	parent_type& cast()
	{
		return *static_cast<parent_type*>(this);
	}

	parent_type const& cast() const
	{
		return *static_cast<parent_type const*>(this);
	}

	friend void swap(WaveVectorDataAxis<ax, O, D>& first, WaveVectorDataAxis<ax, O, D>& second)
	{
		using std::swap;
		swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
		swap(static_cast<K<D>&>(first), static_cast<K<D>&>(second));
	}


protected:

	WaveVectorDataAxis() : parent_type(), K<D>{} {}
};

template<size_t O, size_t D>
using k_grid_type = WaveVectorData<O, D>;

template<Axis ax, size_t O, size_t D>
using k_grid_component_type = WaveVectorDataComponent<ax, O, D>;

template<Axis ax, size_t O, size_t D>
using k_grid_axis_type = WaveVectorDataAxis<ax, O, D>;



namespace symphas::internal
{

	template<Axis ax, size_t O, size_t D>
	auto to_vector_component(WaveVectorData<O, D> const& data)
	{
		return WaveVectorDataComponent<ax, O, D>(data.dims, data.h);
	}

	template<Axis ax, size_t O, size_t D>
	auto to_vector_component(ref<WaveVectorData<O, D>> const& data)
	{
		return WaveVectorDataComponent<ax, O, D>(data.dims, data.h);
	}

	template<typename E>
	struct order_K_type
	{
		static const size_t value = 0;
	};

	template<size_t O, size_t D>
	struct order_K_type<k_grid_type<O, D>>
	{
		static const size_t value = O;
	};

	template<Axis ax, size_t O, size_t D>
	struct order_K_type<k_grid_component_type<ax, O, D>>
	{
		static const size_t value = O;
	};

	template<Axis ax, size_t O, size_t D>
	struct order_K_type<k_grid_axis_type<ax, O, D>>
	{
		static const size_t value = O;
	};



}


template<typename G>
struct order_K_type
{
	static const size_t value = symphas::internal::order_K_type<expr::original_data_t<G>>::value;
};

template<typename E>
struct is_K_type
{
	static const bool value = order_K_type<E>::value > 0;
};

/* type traits supporting the use of querying the enclosing type around the wave number variable
 */


//! Define the base class for the wavenumber.
/*! 
 * Define the base type for this new object to distinguish it from others.
 * The base type is used primarily in the identities.
 */
//DEFINE_BASE_TYPE((size_t O, size_t D), (WaveVectorData<O, D>), K<O>)
//DEFINE_BASE_TYPE((Axis ax, size_t O, size_t D), (WaveVectorDataComponent<ax, O, D>), K<O>)
//DEFINE_BASE_TYPE((Axis ax, size_t O, size_t D), (WaveVectorDataAxis<ax, O, D>), K<O>)

DEFINE_BASE_DATA_ARRAY((size_t O, size_t D), (WaveVectorData<O, D>))
DEFINE_BASE_DATA_ARRAY((Axis ax, size_t O, size_t D), (WaveVectorDataComponent<ax, O, D>))
DEFINE_BASE_DATA_ARRAY((Axis ax, size_t O, size_t D), (WaveVectorDataAxis<ax, O, D>))

// **************************************************************************************

namespace symphas::internal
{

	template<size_t O_, size_t O, size_t D>
	auto change_k_order(WaveVectorData<O, D> const& k)
	{
		return WaveVectorData<O_, D>(k);
	}

	template<size_t O_, Axis ax, size_t O, size_t D>
	auto change_k_order(WaveVectorDataComponent<ax, O, D> const& k)
	{
		return WaveVectorDataComponent<ax, O_, D>(k);
	}

	template<size_t O_, Axis ax, size_t O, size_t D>
	auto change_k_order(WaveVectorDataAxis<ax, O, D> const& k)
	{
		return WaveVectorDataAxis<ax, O_, D>(k);
	}

	//! Make a new string with the name representing the wave vector grid.
	/*!
	 * A new string is initialized with a name representing a wave vector grid of
	 * the prescribed order.
	 *
	 * \param degree The order of the exponent applied to the wavenumber.
	 */
	inline char* new_k_name(size_t degree)
	{
		if (degree % 2 == 0)
		{
			char* name = new char[SYEX_K_EVEN_FMT_LEN + symphas::lib::num_digits(degree) + 1];
			sprintf(name, SYEX_K_EVEN_FMT, degree);
			return name;
		}
		else
		{
			char* name = new char[SYEX_K_ODD_FMT_LEN + symphas::lib::num_digits(degree - 1) + 1];
			sprintf(name, SYEX_K_ODD_FMT, degree - 1);
			return name;
		}
	}

	inline char* new_k_component_name(size_t degree, Axis axis)
	{
		if (degree == 1)
		{
			char* name = new char[SYEX_K_AXIS_0_FMT_LEN + 1];
			sprintf(name, SYEX_K_AXIS_0_FMT, (axis == Axis::X) ? 'x' : (axis == Axis::Y) ? 'y' : (axis == Axis::Z) ? 'z' : '?');
			return name;
		}
		else
		{
			char* name = new char[SYEX_K_COMPONENT_FMT_LEN + symphas::lib::num_digits(degree - 1) + 1];
			sprintf(name, SYEX_K_COMPONENT_FMT, (axis == Axis::X) ? 'x' : (axis == Axis::Y) ? 'y' : (axis == Axis::Z) ? 'z' : '?', degree - 1);
			return name;
		}
	}

	inline char* new_k_axis_name(size_t degree, Axis axis)
	{
		if (degree == 1)
		{
			char* name = new char[SYEX_K_AXIS_0_FMT_LEN + 1];
			sprintf(name, SYEX_K_AXIS_0_FMT, (axis == Axis::X) ? 'x' : (axis == Axis::Y) ? 'y' : (axis == Axis::Z) ? 'z' : '?');
			return name;
		}
		else
		{
			char* name = new char[SYEX_K_AXIS_FMT_LEN + symphas::lib::num_digits(degree) + 1];
			sprintf(name, SYEX_K_AXIS_FMT, (axis == Axis::X) ? 'x' : (axis == Axis::Y) ? 'y' : (axis == Axis::Z) ? 'z' : '?', degree);
			return name;
		}
	}
}

DEFINE_SYMBOL_ID((size_t O, size_t D), (WaveVectorData<O, D>), 
	{ static char* name = symphas::internal::new_k_name(O); return name; })
DEFINE_SYMBOL_ID((Axis ax, size_t O, size_t D), (WaveVectorDataComponent<ax, O, D>), 
	{ static char* name = symphas::internal::new_k_component_name(O, ax); return name; })
DEFINE_SYMBOL_ID((Axis ax, size_t O, size_t D), (WaveVectorDataAxis<ax, O, D>), 
	{ static char* name = symphas::internal::new_k_axis_name(O, ax); return name; })

// **************************************************************************************

/* expression logic
*/

ALLOW_COMBINATION((size_t O, size_t D), (WaveVectorData<O, D>))
ALLOW_COMBINATION((Axis ax, size_t O, size_t D), (WaveVectorDataComponent<ax, O, D>))
ALLOW_COMBINATION((Axis ax, size_t O, size_t D), (WaveVectorDataAxis<ax, O, D>))

DEFINE_FACTOR_COUNT((size_t O1, size_t O2, size_t D), (WaveVectorData<O1, D>), (WaveVectorData<O2, D>), O2 / O1)
DEFINE_FACTOR_COUNT((Axis ax, size_t O1, size_t O2, size_t D), (WaveVectorDataComponent<ax, O1, D>), (WaveVectorDataComponent<ax, O2, D>), O2 / O1)
DEFINE_FACTOR_COUNT((Axis ax, size_t O1, size_t O2, size_t D), (WaveVectorDataAxis<ax, O1, D>), (WaveVectorDataAxis<ax, O2, D>), O2 / O1)

// **************************************************************************************



//! Overload of the identity for the multiplication between two wave vectors.
/*! 
 * Operations supporting the variables representing the wave number
 * currently only implemented when the variable contains the reference to 
 * K<O>::Grid; but this is always how the K<O>::Grid is originally constructed.
 */
template<typename K1, expr::exp_key_t X1, typename K2, expr::exp_key_t X2, 
	size_t O1 = order_K_type<K1>::value * expr::_Xk_t<X1>::N * expr::_Xk_t<X2>::D,
	size_t O2 = order_K_type<K2>::value * expr::_Xk_t<X2>::N * expr::_Xk_t<X1>::D,
	typename std::enable_if_t<(order_K_type<K1>::value > 0 && order_K_type<K2>::value > 0), int> = 0>
auto operator*(Term<K1, X1> const& a, Term<K2, X2> const& b)
{
	constexpr size_t dimension = expr::grid_dim<K1>::dimension;
	constexpr size_t D = expr::_Xk_t<X1>::D * expr::_Xk_t<X2>::D;

	if constexpr (expr::_Xk_t<X1>::sign == expr::_Xk_t<X2>::sign)
	{
		constexpr size_t O = O1 + O2;
		constexpr size_t O_reduced = O / GCD_of<O, D>;

		return symphas::internal::to_term_element<expr::Xk<1, D / GCD_of<O, D>, expr::_Xk_t<X1>::sign>>(symphas::internal::change_k_order<O_reduced>(a.data()));
	}
	else
	{
		constexpr size_t O = (O1 > O2) ? (O1 - O2) : (O2 - O1);
		constexpr size_t O_reduced = O / GCD_of<O, D>;

		constexpr bool flag = (expr::_Xk_t<X1>::sign) ? (O1 > O2) : (O2 > O1);

		if constexpr (O == 0)
		{
			return symphas::internal::to_term_element<expr::Xk<0, 1>>(symphas::internal::change_k_order<2>(a.data()));
		}
		else
		{
			return symphas::internal::to_term_element<expr::Xk<1, D / GCD_of<O, D>, flag>>(symphas::internal::change_k_order<O_reduced>(a.data()));
		}
	}

}



#undef SYEX_k_grid_type_FMT
#undef SYEX_k_grid_type_FMT_LATEX
#undef SYEX_k_grid_type_FMT_LEN
#undef SYEX_k_grid_type_FMT_LATEX_LEN



