
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
 * PURPOSE: Defines random noise, which represents, for instance, thermal
 * fluctuations in a field. 
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"
#include "expressiontypek.h"


//! \cond

#define NOISE_MEAN 0
#define NOISE_STD 1

#ifdef LATEX_PLOT

#define SYEX_NOISE_OP_STR "\\eta"
#define SYEX_NOISE_NONE_OP_STR "\\hat{" SYEX_NOISE_OP_STR "}"

#else

#define SYEX_NOISE_OP_STR "n"
#define SYEX_NOISE_NONE_OP_STR SYEX_NOISE_OP_STR "'"

#endif


//! \endcond

// **************************************************************************************

enum class NoiseType
{
	WHITE,
	NONE,
	DECAY_EXP,
	DECAY_POLY,
	POISSON
};

namespace expr
{

	template<typename F, size_t D>
	struct noise_data_with_function : Grid<complex_t, D>
	{
		using parent_type = Grid<complex_t, D>;
		using parent_type::values;
		using parent_type::len;
		using parent_type::dims;

		noise_data_with_function(F f, const len_type* dims, const double* h, double *dt, double intensity = 1.0)
			: parent_type(dims), dt{ dt }, intensity{ intensity }, h{ 0 }, f{ f }
		{
			std::copy(h, h + D, this->h);
			update();
		}

		scalar_t operator[](iter_type n) const
		{
			return values[n].real();
		}

	protected:

		void update(bool fourier_space = true)
		{
			using std::exp;
			using std::sqrt;

			std::random_device rd;
			std::mt19937 gen(rd());
			std::normal_distribution<> dis(NOISE_MEAN, NOISE_STD);

			k_field<D>::template fill<Axis::NONE, 3>(*static_cast<parent_type*>(this), dims, h);

			double H = 1;
			for (iter_type i = 0; i < D; ++i) H *= h[i];
			scalar_t R = sqrt(*dt) / H;

			for (iter_type n = 0; n < len; ++n)
			{
				double E_n = f(intensity, values[n].real());
				complex_t rIi{ dis(gen) * R, dis(gen) * R };
				values[n] = rIi * E_n;
			}

			if (!fourier_space)
			{
				auto plan = symphas::dft::new_fftw_plan<D, complex_t, complex_t>{}(values, values, dims, false, true);
				symphas::dft::fftw_execute(plan);
				symphas::dft::fftw_destroy_plan(plan);
			}
			else
			{
				//auto plan = symphas::dft::new_fftw_plan<D, complex_t, complex_t>{}(values, values, dims, false);
				//symphas::dft::fftw_execute(plan);
				//symphas::dft::fftw_destroy_plan(plan);

			}
		}

		double *dt;
		double intensity;
		double h[D];
		F f;
	};

	inline auto eigen_exponential(scalar_t intensity, scalar_t lambda_n)
	{
		return std::exp(-intensity * lambda_n / 2);
	}

	inline auto eigen_polynomial(scalar_t intensity, scalar_t lambda_n)
	{
		return std::pow(lambda_n, -intensity);
	}

	inline constexpr auto eigen_one(scalar_t intensity, scalar_t lambda_n)
	{
		return intensity * lambda_n;
	}


	template<NoiseType nt, size_t D>
	struct noise_data;

	template<size_t D>
	struct noise_data<NoiseType::WHITE, D> : Grid<scalar_t, D>
	{
		using parent_type = Grid<scalar_t, D>;
		using parent_type::values;
		using parent_type::len;
		using parent_type::dims;
		using parent_type::operator[];

		noise_data(const len_type* dims, const double* h, double *dt, double intensity = 1.0) :
			parent_type(dims), H{ 1 }, dt{ dt }, intensity{ intensity }
		{
			for (iter_type i = 0; i < D; ++i) H *= h[i];
			update();
		}

		void update() 
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::normal_distribution<> dis(NOISE_MEAN, NOISE_STD);

			scalar_t sq_corr = intensity * std::sqrt(*dt) / H;
			for (iter_type i = 0; i < len; ++i)
			{
				values[i] = sq_corr * dis(gen);
			}
		}

		double H;
		double *dt;
		scalar_t intensity;
	};

	template<size_t D>
	struct noise_data<NoiseType::POISSON, D> : Grid<scalar_t, D>
	{
		using parent_type = Grid<scalar_t, D>;
		using parent_type::values;
		using parent_type::len;
		using parent_type::dims;
		using parent_type::operator[];

		noise_data(const len_type* dims, const double* h, double* dt, double intensity = 1.0) :
			parent_type(dims), H{ 1 }, dt{ dt }, intensity{ intensity }
		{
			for (iter_type i = 0; i < D; ++i) H *= h[i];
			update();
		}

		void update()
		{
			time += *dt;
			if (time >= next_time)
			{
				std::random_device rd;
				std::mt19937 gen(rd());
				std::normal_distribution<> dis(NOISE_MEAN, NOISE_STD);

			    scalar_t sq_corr = intensity * std::sqrt(*dt) / H;
				for (iter_type i = 0; i < len; ++i)
				{
					values[i] = sq_corr * dis(gen);
				}


				// find new nexttime
			}
		}

		double H;
		double* dt;
		double time;
		double next_time;
		scalar_t intensity;
	};

	template<size_t D>
	struct noise_data<NoiseType::DECAY_EXP, D> : noise_data_with_function<decltype(&eigen_exponential), D>
	{
		using parent_type = noise_data_with_function<decltype(&eigen_exponential), D>;
		using parent_type::operator[];

		noise_data(const len_type* dims, const double* h, double *dt, double intensity = 1.0)
			: parent_type(&eigen_exponential, dims, h, dt, intensity) {}

		void update()
		{
			parent_type::update(false);
		}
	};

	template<size_t D>
	struct noise_data<NoiseType::DECAY_POLY, D> : noise_data_with_function<decltype(&eigen_polynomial), D>
	{
		using parent_type = noise_data_with_function<decltype(&eigen_polynomial), D>;
		using parent_type::operator[];

		noise_data(const len_type* dims, const double* h, double *dt, double intensity = 1.0)
			: parent_type(&eigen_polynomial, dims, h, dt, intensity) {}

		void update()
		{
			parent_type::update(false);
		}
	};


	template<size_t D>
	struct noise_data<NoiseType::NONE, D> : noise_data_with_function<decltype(&eigen_one), D>
	{
		using parent_type = noise_data_with_function<decltype(&eigen_one), D>;
		using parent_type::operator[];

		noise_data(const len_type* dims, const double* h, double* dt, double intensity = 1.0)
			: parent_type(&eigen_one, dims, h, dt, intensity) {}

		void update()
		{
			parent_type::update(true);
		}
	};

	template<size_t D0, NoiseType nt, size_t D>
	struct noise_data_axis;

	template<NoiseType nt, size_t D>
	struct noise_data_axis<3, nt, D> : noise_data_axis<2, nt, D>
	{
		using parent_type = noise_data_axis<2, nt, D>;
		using parent_type::parent_type;

		noise_data_axis(const len_type* dims, const double* h, double* dt, double intensity = 1.0)
			: parent_type(dims, h, dt, intensity), data(dims, h, dt, intensity) {}

		vector_t<3> operator[](iter_type n) const
		{
			return { parent_type::operator[](n)[0], parent_type::operator[](n)[1], data[n]};
		}

		void update()
		{
			data.update();
			parent_type::update();
		}

		noise_data<nt, D> data;
	};

	template<NoiseType nt, size_t D>
	struct noise_data_axis<2, nt, D> : noise_data_axis<1, nt, D>
	{
		using parent_type = noise_data_axis<1, nt, D>;

		noise_data_axis(const len_type* dims, const double* h, double* dt, double intensity = 1.0)
			: parent_type(dims, h, dt, intensity), data(dims, h, dt, intensity) {}

		vector_t<2> operator[](iter_type n) const
		{
			return { parent_type::operator[](n)[0], data[n]};
		}

		void update()
		{
			data.update();
			parent_type::update();
		}

		noise_data<nt, D> data;
	};

	template<NoiseType nt, size_t D>
	struct noise_data_axis<1, nt, D> : noise_data<nt, D>
	{
		using parent_type = noise_data<nt, D>;

		noise_data_axis(const len_type* dims, const double* h, double* dt, double intensity = 1.0)
			: parent_type(dims, h, dt, intensity) {}

		vector_t<1> operator[](iter_type n) const
		{
			return { parent_type::operator[](n) };
		}

		void update()
		{
			parent_type::update();
		}

	};
}

template<NoiseType nt, typename T, size_t D>
struct NoiseData : expr::noise_data<nt, D>
{
	using parent_type = expr::noise_data<nt, D>;
	using parent_type::parent_type;
	using parent_type::operator[];
	using parent_type::update;

protected:

	NoiseData() : parent_type() {}
};

template<NoiseType nt, typename T, size_t D>
struct NoiseData<nt, any_vector_t<T, D>, D> : expr::noise_data_axis<D, nt, D>
{
	using parent_type = expr::noise_data_axis<D, nt, D>;
	using parent_type::parent_type;
	using parent_type::update;
	using parent_type::operator[];


protected:

	NoiseData() : parent_type() {}
};

template<typename T, size_t D>
struct NoiseData<NoiseType::POISSON, any_vector_t<T, D>, D> : expr::noise_data_axis<D, NoiseType::POISSON, D>
{
	using parent_type = expr::noise_data_axis<D, NoiseType::POISSON, D>;
	using parent_type::update;
	using parent_type::operator[];

	NoiseData(const len_type* dims, const double* h, double* dt, double intensity = 1.0)
		: parent_type(dims, h, dt, intensity) {}



protected:

	NoiseData() : parent_type() {}
};


namespace expr
{

	template<NoiseType nt, typename T, size_t D>
	auto make_noise(const len_type* dimensions, const double* h, const double *dt, double intensity = 1.0)
	{
		return make_term(NoiseData<nt, T, D>(dimensions, h, const_cast<double*>(dt), intensity));
	}

	template<typename T, size_t D>
	auto make_white_noise(const len_type* dimensions, const double* h, const double *dt, double intensity = 1.0)
	{
		return make_term(NoiseData<NoiseType::WHITE, T, D>(dimensions, h, dt, intensity));
	}

	template<Axis ax, NoiseType nt, typename T, size_t D>
	auto resolve_axis_component(NoiseData<nt, any_vector_t<T, D>, D> const& data)
	{
		if constexpr (ax == Axis::X)
		{
			return as_component_data<ax, D>(static_cast<noise_data<nt, D> const*>(&data)->values);
		}
		else if constexpr (ax == Axis::Y)
		{
			return as_component_data<ax, D>(static_cast<noise_data_axis<2, nt, D> const*>(&data)->data.values);
		}
		else if constexpr (ax == Axis::Z)
		{
			return as_component_data<ax, D>(static_cast<noise_data_axis<3, nt, D> const*>(&data)->data.values);
		}
		else
		{
			return as_component_data<ax, D>((double*)nullptr);
		}
	}

	template<Axis ax, NoiseType nt, typename T, size_t D>
	auto resolve_axis_component(NoiseData<nt, any_vector_t<T, D>, D>& data)
	{
		if constexpr (ax == Axis::X)
		{
			return as_component_data<ax, D>(static_cast<noise_data<nt, D>*>(&data)->values);
		}
		else if constexpr (ax == Axis::Y)
		{
			return as_component_data<ax, D>(static_cast<noise_data_axis<2, nt, D>*>(&data)->data.values);
		}
		else if constexpr (ax == Axis::Z)
		{
			return as_component_data<ax, D>(static_cast<noise_data_axis<3, nt, D>*>(&data)->data.values);
		}
		else
		{
			return as_component_data<ax, D>((double*)nullptr);
		}
	}

}


DEFINE_BASE_DATA_ARRAY((NoiseType nt, typename T, size_t D), (NoiseData<nt, T, D>))
DEFINE_SYMBOL_ID((NoiseType nt, typename T, size_t D), (NoiseData<nt, T, D>), { return SYEX_NOISE_OP_STR; })
DEFINE_SYMBOL_ID((typename T, size_t D), (NoiseData<NoiseType::NONE, T, D>), { return SYEX_NOISE_NONE_OP_STR; })



