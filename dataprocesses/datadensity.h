#pragma once

#include "processdynamic.h"

DEFINE_ALGORITHM_POINT(Density)
DEFINE_ALGORITHM_DYNAMIC(Density)
DEFINE_DATA(Density, POINT, ALG_POINT, ALG_DYNAMIC)



ALGORITHM_POINT_1D_DECLARATION(Density)
{
	Y density = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		density += data_y[i];
	}
	return point_data<Y>{ 0, density };
}


ALGORITHM_POINT_2D_DECLARATION(Density)
{
	Y density = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		density += data_y[i];
	}
	return point_data<Y>{ 0, density };
}

ALGORITHM_POINT_3D_DECLARATION(Density)
{
	Y density = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		density += data_y[i];
	}
	return point_data<Y>{ 0, density };
}


ALGORITHM_DYNAMIC_DECLARATION(Density)
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






