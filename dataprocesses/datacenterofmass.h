#pragma once

#include "processdynamic.h"

DEFINE_ALGORITHM_POINT(CenterOfMass)
DEFINE_ALGORITHM_DYNAMIC(CenterOfMass)
DEFINE_DATA(CenterOfMass, POINT, ALG_POINT, ALG_DYNAMIC)



ALGORITHM_POINT_1D_DECLARATION(CenterOfMass)
{
	axis_1d_type position{};
	Y sum{};
	for (iter_type i = 0; i < len; ++i)
	{
		position += data_y[i] * i;
		sum += data_y[i];
	}

	double const dx = ((*(data_x + 1)) - (*data_x));
	position *= dx / sum;

	return symphas::lib::to_field_data(0, position);
}

ALGORITHM_POINT_2D_DECLARATION(CenterOfMass)
{
	axis_2d_type position{};
	Y sum{};
	for (iter_type i = 0; i < L; ++i)
	{
		for (iter_type j = 0; j < M; ++j)
		{
			iter_type n = i + j * L;
			position[0] += data_y[n] * i;
			position[1] += data_y[n] * j;
			sum += data_y[n];
		}
	}

	double const dx = ((*(data_x + 1))[0] - (*data_x)[0]);
	double const dy = ((*(data_x + L))[1] - (*data_x)[1]);

	position[0] *= dx / sum;
	position[1] *= dy / sum;

	return symphas::lib::to_field_data(0, position);
}

ALGORITHM_POINT_3D_DECLARATION(CenterOfMass)
{
	axis_3d_type position{};
	Y sum{};
	for (iter_type i = 0; i < L; ++i)
	{
		for (iter_type j = 0; j < M; ++j)
		{
			for (iter_type k = 0; k < N; ++k)
			{
				iter_type n = i + j * L + k * L * M;
				position[0] += data_y[n] * i;
				position[1] += data_y[n] * j;
				position[2] += data_y[n] * k;
				sum += data_y[n];
			}
		}
	}

	double const dx = ((*(data_x + 1))[0] - (*data_x)[0]);
	double const dy = ((*(data_x + L))[1] - (*data_x)[1]);
	double const dz = ((*(data_x + L * M))[2] - (*data_x)[2]);
	position[0] *= dx / sum;
	position[1] *= dy / sum;
	position[2] *= dz / sum;

	return symphas::lib::to_field_data(0, position);
}


ALGORITHM_DYNAMIC_DECLARATION(CenterOfMass)
{
	using result_type = std::tuple<axis_nd_t<D>, int>;
	result_type* ys = new result_type[len];
	scalar_t* xs = new scalar_t[len];

	iter_type id = 0;
	iter_type x_last = -1;
	for (iter_type i = 0; i < len; ++i)
	{
		id = (x_last == data_x[i]) ? id + 1 : 0;
		x_last = data_x[i];

		xs[i] = static_cast<scalar_t>(data_x[i]);
		std::get<0>(ys[i]) = data_y[i]->data_y();
		std::get<1>(ys[i]) = id;
	}

	return symphas::lib::to_field_data(xs, ys, len);
}



