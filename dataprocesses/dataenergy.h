#pragma once

#include "processdynamic.h"

DEFINE_ALGORITHM_POINT(Energy)
DEFINE_ALGORITHM_DYNAMIC(Energy)
DEFINE_DATA(Energy, POINT, ALG_POINT, ALG_DYNAMIC)



ALGORITHM_POINT_1D_DECLARATION(Energy)
{
	Y energy = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		energy += data_y[i];
	}
	return point_data<Y>{ 0, energy* (1. / len) };
}

ALGORITHM_POINT_2D_DECLARATION(Energy)
{
	Y energy = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		energy += data_y[i];
	}
	return point_data<Y>{ 0, energy* (1. / len) };
}

ALGORITHM_POINT_3D_DECLARATION(Energy)
{
	Y energy = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		energy += data_y[i];
	}
	return point_data<Y>{ 0, energy* (1. / len) };
}


ALGORITHM_DYNAMIC_DECLARATION(Energy)
{
	Y* ys = new Y[len];
	scalar_t* xs = new scalar_t[len];

	for (iter_type i = 0; i < len; ++i)
	{
		xs[i] = static_cast<scalar_t>(data_x[i]);
		ys[i] = data_y[i]->data_y();
	}

	return scalar_data<Y>{ xs, ys, len };
}



