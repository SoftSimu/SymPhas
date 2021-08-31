
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


#include "expressionlogic.h"
#include "expressionproperties.h"
#include "expressionaggregates.h"

//! \cond

#define SYEX_K_OP_STR "|k|"

//! The display format for printing the wavenumber.
#define SYEX_K_FMT SYEX_K_OP_STR SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B
//! The latex display format for printing the wavenumber.
#define SYEX_K_FMT_LEN (STR_ARR_LEN(SYEX_K_OP_STR SYEX_POW_SEP_A SYEX_POW_SEP_B) - 1)

//! \endcond

// **************************************************************************************

namespace expr
{


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
		static void fill(scalar_t* into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			double dk_i = 2 * symphas::PI / (*h * L);

			iter_type ii = 0;
			for (iter_type i = 0; i < L; ++i)
			{
				double kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
				double kk = kx * kx;
				into[ii++] = std::pow(kk, O / 2);
			}
			into[0] = symphas::EPS * std::pow(1.0, -1 * (O / 2));
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
		static void fill(scalar_t* into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			len_type M = dims[1];

			double dk_i = 2 * symphas::PI / (h[0] * L);
			double dk_j = 2 * symphas::PI / (h[1] * M);

			iter_type ii = 0;
			for (iter_type j = 0; j < M; ++j)
			{
				double ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
				for (iter_type i = 0; i < L; ++i)
				{
					double kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
					double kk = kx * kx + ky * ky;
					into[ii++] = std::pow(kk, O / 2);
				}
			}
			into[0] = symphas::EPS * std::pow(1.0, -1 * (O / 2));
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
		static void fill(scalar_t* into, const len_type* dims, const double* h)
		{
			len_type L = dims[0];
			len_type M = dims[1];
			len_type N = dims[2];

			double
				dk_i = 2 * symphas::PI / (h[0] * L),
				dk_j = 2 * symphas::PI / (h[1] * M),
				dk_k = 2 * symphas::PI / (h[2] * N);

			iter_type ii = 0;
			for (iter_type k = 0; k < N; ++k)
			{
				double kz = (k < N / 2) ? k * dk_k : (k - N) * dk_k;
				for (iter_type j = 0; j < M; ++j)
				{
					double ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
					for (iter_type i = 0; i < L; ++i)
					{
						double kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
						double kk = kx * kx + ky * ky + kz * kz;
						into[ii++] = std::pow(kk, O / 2);
					}
				}
			}
			into[0] = symphas::EPS * std::pow(1.0, -1 * (O / 2));
		}

	};


	// **************************************************************************************


	//! Predicate indicating whether the given type is a wavenumber field.
	template<typename E>
	struct is_K_type;

}


//! Enclosing template of the wavenumber field.
/*!
 * \tparam O The exponential order (power) of the wavenumber.
 */
template<size_t O>
struct K
{
	//! The grid representing the wavenumber data.
	/*!
	 * \tparam T The type of the wavenumber values, typically ::scalar_t.
	 * \tparam D The dimension of the wavenumber field.
	 */
	template<typename T, size_t D>
	struct WaveVectorData;
};

template<size_t O>
template<typename T, size_t D>
struct K<O>::WaveVectorData : Grid<T, D>, K<O>
{
	using Grid<T, D>::values;
	using Grid<T, D>::dims;
	using Grid<T, D>::len;

	WaveVectorData(const len_type* dimensions, double const* h) : Grid<T, D>{ dimensions }, h{ 0 }
	{
		if (h != nullptr)
		{
			std::copy(h, h + D, this->h);
			expr::k_field<D>::template fill<O>(values, dims, h);
		}
	}

	WaveVectorData(std::initializer_list<len_type> dimensions, double const* h) : Grid<T, D>{ dimensions }, h{ 0 }
	{
		std::copy(h, h + D, this->h);
		expr::k_field<D>::template fill<O>(values, dims, h);
	}

	template<typename K, typename std::enable_if_t<expr::is_K_type<K>::value, int> = 0>
	WaveVectorData(K const& other) : Grid<T, D>{ other }, h{ 0 }
	{
		std::copy(other.widths(), other.widths() + D, h);
	}

	WaveVectorData(K<O>::WaveVectorData<T, D> const& other) : Grid<T, D>{ other }, h{ 0 }
	{
		std::copy(other.h, other.h + D, h);
	}

	WaveVectorData(K<O>::WaveVectorData<T, D>&& other) : WaveVectorData()
	{
		swap(*this, other);
	}

	K<O>::WaveVectorData<T, D>& operator=(K<O>::WaveVectorData<T, D> other)
	{
		swap(*this, other);
		return *this;
	}



	friend void swap(K<O>::WaveVectorData<T, D>& first, K<O>::WaveVectorData<T, D>& second)
	{
		using std::swap;
		swap(static_cast<Grid<T, D>&>(first), static_cast<Grid<T, D>&>(second));
		swap(first.h, second.h);
	}

	
	const double* widths() const
	{
		return h;
	}

	double width(int i) const
	{
		return h[i];
	}

protected:

	double h[D];
	WaveVectorData() : Grid<T, D>{}, h{ 0 } {}
};

template<size_t O, size_t D>
using K_Grid = typename K<O>::template WaveVectorData<scalar_t, D>;


/* type traits supporting the use of querying the enclosing type around the wave number variable
 */

template<typename E>
struct order_K_type
{
protected:

	template<size_t O>
	struct K_return_type
	{
		static const size_t value = O;
	};

	template<size_t O>
	static constexpr auto _check_cast(K<O>*)
	{
		return K_return_type<O>{};
	}

	template<typename A>
	static constexpr auto cast(A* a) -> decltype(_check_cast(a))
	{
		return _check_cast(a);
	}

	template<typename A>
	static constexpr auto cast(...)
	{
		return K_return_type<0>{};
	}

	static constexpr auto call_wrap(E* a)
	{
		return cast<E>(a);
	}

	using wrapped_type = typename std::invoke_result_t<decltype(&order_K_type<E>::call_wrap), E*>;


public:

	static const size_t value = wrapped_type::value;
};


template<typename E>
struct expr::is_K_type
{
	static const bool value = order_K_type<E>::value > 0;
};


//! Define the base class for the wavenumber.
/*! 
 * Define the base type for this new object to distinguish it from others.
 * The base type is used primarily in the identities.
 */
DEFINE_BASE_TYPE_CONDITION((typename K_t), (K_t), expr::is_K_type<K_t>::value, K<order_K_type<K_t>::value>)

// **************************************************************************************


//! Make a new string with the name representing the wave vector grid.
/*!
 * A new string is initialized with a name representing a wave vector grid of
 * the prescribed order.
 * 
 * \param degree The order of the exponent applied to the wavenumber.
 */
inline char* new_k_name(size_t degree)
{
	char* name = new char[SYEX_K_FMT_LEN + symphas::lib::num_digits(degree) + 1];
	sprintf(name, SYEX_K_FMT, degree);
	return name;
}

DEFINE_SYMBOL_ID_CONDITION((typename K_t), (K_t), (expr::is_K_type<K_t>::value),
	{
	static char* name = new_k_name(order_K_type<K_t>::value);
	return name;
	});

// **************************************************************************************

/* expression logic
*/

ALLOW_COMBINATION_CONDITION((typename K_t), (K_t), expr::is_K_type<K_t>::value);
RESTRICT_MULTIPLICATION_CONDITION((typename K_t), (K_t), expr::is_K_type<K_t>::value);
ALLOW_DIVISION_CONDITION((size_t O1, size_t O2), (K<O1>, K<O2>), O2 > O1);

// **************************************************************************************


//! Overload of the identity for the multiplication between two wave vectors.
/*! 
 * Operations supporting the variables representing the wave number
 * currently only implemented when the variable contains the reference to 
 * K<O>::Grid; but this is always how the K<O>::Grid is originally constructed.
 */
template<size_t O1, size_t O2>
struct expr::identity_rule<expr::IdentityType::MUL, K<O1>, K<O2>>
{
	static bool get_value() { return (O1 > 0) && (O2 > 0); }
	static const bool enable = (O1 > 0) && (O2 > 0);

	template<typename E1, typename E2>
	static auto apply(E1 const& a, E2 const& b)
	{
		constexpr size_t D = expr::grid_dim<E1>::dimension;
		double const* h = a.data.widths();
		len_type const* dims = expr::property::data_dimensions(a, b);
		return expr::make_literal(a.value * b.value) * expr::make_op(K_Grid<O1 + O2, D>(dims, h));
	}
};


namespace symphas::internal
{

	template<size_t O1, size_t O2, typename E1, typename E2, typename std::enable_if_t<(O1 == O2), int> = 0>
	auto k_wave_identity(E1 const& a, E2 const& b)
	{
		constexpr size_t D = expr::grid_dim<E1>::dimension;
		return expr::make_literal(a.value / b.value);
	}

	template<size_t O1, size_t O2, typename E1, typename E2, typename std::enable_if_t<(O1 > O2), int> = 0>
	auto k_wave_identity(E1 const& a, E2 const& b)
	{
		constexpr size_t D = expr::grid_dim<E1>::dimension;
		double const* h = a.data.widths();
		len_type const* dims = expr::property::data_dimensions(a, b);
		return expr::make_literal(a.value / b.value) * expr::make_op(K_Grid<O1 - O2, D>(dims, h));
	}

	template<size_t O1, size_t O2, typename E1, typename E2, typename std::enable_if_t<(O1 < O2), int> = 0>
	auto k_wave_identity(E1 const& a, E2 const& b)
	{
		constexpr size_t D = expr::grid_dim<E1>::dimension;
		double const* h = a.data.widths();
		len_type const* dims = expr::property::data_dimensions(a, b);
		return expr::make_literal(a.value / b.value) / expr::make_op(K_Grid<O2 - O1, D>(dims, h));
	}

}

//! Overload of the identity for the division between two wave vectors.
template<size_t O1, size_t O2>
struct expr::identity_rule<expr::IdentityType::DIV, K<O1>, K<O2>>
{
	static const bool enable = (O1 > 0) && (O2 > 0);

	template<typename E1, typename E2>
	static auto apply(E1 const& a, E2 const& b)
	{
		return symphas::internal::k_wave_identity<O1, O2>(a, b);
	}

};


#undef SYEX_K_GRID_FMT
#undef SYEX_K_GRID_FMT_LATEX
#undef SYEX_K_GRID_FMT_LEN
#undef SYEX_K_GRID_FMT_LATEX_LEN



