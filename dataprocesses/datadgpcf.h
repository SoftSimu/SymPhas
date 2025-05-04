#pragma once

#include <cmath>

#include "datapcf.h"
#include "processdynamic.h"

#define RUNNING_AVG 6
#define CROSSOVER_TOL 0.5


DEFINE_ALGORITHM_POINT(DGPCF)
DEFINE_ALGORITHM_DYNAMIC(DGPCF)
DEFINE_DATA(DGPCF, POINT, ALG_POINT, ALG_DYNAMIC)

inline scalar_t get_crossover(double cp, axis_coord_t *data_x, scalar_t *data_y, len_type length)
{
	for (iter_type n = 0; n < length - RUNNING_AVG; ++n)
	{
		scalar_t y1 = data_y[n + RUNNING_AVG - 1];
		if (y1 < cp)
		{
			scalar_t x0 = data_x[n];
			scalar_t x1 = data_x[n + RUNNING_AVG - 1];

			scalar_t
				avgx = 0,
				avgy = 0,
				mnum = 0,
				mden = 0;
			for (iter_type i = n; i < n + RUNNING_AVG; ++i)
			{
				avgx += data_x[i];
				avgy += data_y[i];
			}
			avgx /= RUNNING_AVG;
			avgy /= RUNNING_AVG;

			for (iter_type i = n; i < n + RUNNING_AVG; ++i)
			{
				scalar_t dx = (data_x[i] - avgx);
				mnum += dx * (data_y[i] - avgy);
				mden += dx * dx;
			}
			scalar_t m = mnum / mden;
			scalar_t b = avgy - m * avgx;

			scalar_t crossover_x = (cp - b) / m;

			if (crossover_x <= x1 && crossover_x >= x0)
			{
				return crossover_x;
			}
			else if (m * x0 + b <= cp && m < 0)
			{
				return x0;
			}
		}
	}
	return 0;
}

ALGORITHM_POINT_1D_DECLARATION(DGPCF)
{
	double const
		dx = ((*(data_x + 1)) - (*data_x));

	auto pcf_data = RUN_ALGORITHM_VECTOR_1D(PCF)(ALGORITHM_ARGS_1D);
	double r = (1.0 / len);

	T phase0 = std::reduce(std::execution::par_unseq, data_y, data_y + len);
	T p = (r * phase0);
	T p2 = p * p;

	auto pr = pcf_data[L / 2] - p2;

#	pragma omp parallel for
	for (iter_type n = 0; n < len; ++n)
	{
		pcf_data[n] -= p2;
		pcf_data[n] /= pr;
	}


	auto rd = symphas::lib::radial_avg(pcf_data);
	scalar_t r0 = get_crossover(CROSSOVER_TOL, rd.x, rd.y, rd.length());

	using pt_t = std::tuple<scalar_t, scalar_t, decltype(rd)>;
	return point_data<pt_t>{ 0, std::make_tuple(r0, r0, rd) };
}

ALGORITHM_POINT_2D_DECLARATION(DGPCF)
{

	double const
		dx = ((*(data_x + L))[0] - (*data_x)[0]),
		dy = ((*(data_x + 1))[1] - (*data_x)[1]);

	auto pcf_data = RUN_ALGORITHM_VECTOR_2D(PCF)(ALGORITHM_ARGS_2D);
	auto phase0 = std::reduce(std::execution::par, data_y, data_y + len);
	auto p = (std::pow(len, -1) * phase0);
	auto p2 = p * p;

	
	auto pr = pcf_data[(M / 2) * L + (L / 2)] - p2;
	auto rd = symphas::lib::radial_avg(pcf_data);
	scalar_t r01 = get_crossover(CROSSOVER_TOL, rd.x, rd.y, rd.length());


#	pragma omp parallel for
	for (iter_type n = 0; n < len; ++n)
	{
		pcf_data[n] -= p2;
		pcf_data[n] /= pr;
	}


	rd = symphas::lib::radial_avg(pcf_data);
	scalar_t r02 = get_crossover(CROSSOVER_TOL, rd.x, rd.y, rd.length());

	using pt_t = std::tuple<scalar_t, scalar_t, decltype(rd)>;
	return point_data<pt_t>{ 0, std::make_tuple(r01, r02, rd) };
}

ALGORITHM_POINT_3D_DECLARATION(DGPCF)
{
	double const
		dx = ((*(data_x + L * M))[0] - (*data_x)[0]),
		dy = ((*(data_x + L))[1] - (*data_x)[1]),
		dz = ((*(data_x + 1))[2] - (*data_x)[2]);

	auto pcf_data = RUN_ALGORITHM_VECTOR_3D(PCF)(ALGORITHM_ARGS_3D);
	double r = (1.0 / len);

	auto phase0 = std::reduce(std::execution::par_unseq, data_y, data_y + len);
	auto p = (r * phase0);
	auto p2 = p * p;

	auto pr = pcf_data[(N / 2) * L * M + (M / 2) * L + L / 2] - p2;

#	pragma omp parallel for
	for (iter_type n = 0; n < len; ++n)
	{
		pcf_data[n] -= p2;
		pcf_data[n] /= pr;
	}


	auto rd = symphas::lib::radial_avg(pcf_data);
	scalar_t r0 = get_crossover(CROSSOVER_TOL, rd.x, rd.y, rd.length());

	using pt_t = std::tuple<scalar_t, scalar_t, decltype(rd)>;
	return point_data<pt_t>{ 0, std::make_tuple(r0, r0, rd) };
}


ALGORITHM_DYNAMIC_DECLARATION(DGPCF)
{
	auto* ys = new std::pair<scalar_t[2], std::vector<std::pair<scalar_t, scalar_t>>>[len];
	scalar_t* xs = new scalar_t[len];

	for (iter_type i = 0; i < len; ++i)
	{
		xs[i] = static_cast<scalar_t>(data_x[i]);
		
		auto&& [r01, r02, rd] = data_y[i]->data_y();
		ys[i].first[0] = r01;
		ys[i].first[1] = r02;


		len_type lena = rd.length();
		ys[i].second.resize(lena);

		for (iter_type n = 0; n < lena; ++n)
		{
			ys[i].second[n] = std::make_pair(
				rd.x[n], 
				rd.y[n]);
		}
	}

	return scalar_data<std::pair<scalar_t[2], std::vector<std::pair<scalar_t, scalar_t>>>>{ xs, ys, len };
}



