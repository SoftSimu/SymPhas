#pragma once

#include <algorithm>


#include "processdynamic.h"

DEFINE_ALGORITHM_POINT(DomainGrowth)
DEFINE_ALGORITHM_DYNAMIC(DomainGrowth)
DEFINE_DATA(DomainGrowth, POINT, ALG_POINT, ALG_DYNAMIC)
DEFINE_DATA_ALIAS_NAME(DomainGrowth, DG)

#define NUM_POINTS std::pow(len, 0.65)
#define NUM_RAYS 200
#define RAY_LENGTH_2 std::sqrt(L * L + M * M);
#define RAY_LENGTH_3 std::sqrt(L * L + M * M + N * N);

ALGORITHM_POINT_1D_DECLARATION(DomainGrowth)
{
	scalar_t phase = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		phase += (data_y[i] > 0) ? 1 : 0;
	}
	return point_data<Y>{ 0, phase * (1. / len) };
}

ALGORITHM_POINT_2D_DECLARATION(DomainGrowth)
{
	double Pr = std::sqrt(static_cast<double>(len) / NUM_POINTS);
	len_type PnL = static_cast<iter_type>(L / Pr);
	len_type PnM = static_cast<iter_type>(M / Pr);

	double ray_len = RAY_LENGTH_2;
	std::vector<scalar_t> lengths;
	lengths.reserve(4 * NUM_RAYS);

	std::vector<scalar_t> all_lengths;
	for (iter_type i = 0; i < PnL; ++i)
	{
		for (iter_type j = 0; j < PnM; ++j)
		{
			iter_type x0 = static_cast<iter_type>(i * Pr);
			iter_type y0 = static_cast<iter_type>(j * Pr);

			double dtheta = NUM_RAYS / (2 * symphas::PI);
			lengths.clear();

			for (iter_type t = 0; t < NUM_RAYS; ++t)
			{
				double theta = t * dtheta;
				double 
					dx = std::cos(theta),
					dy = std::sin(theta);

				double xs = x0;
				double ys = y0;

				double startx = xs;
				double starty = ys;

				iter_type n0 = x0 * L + y0;
				bool plus = (data_y[n0] > 0);
				bool first_flipped = false;

				for (iter_type d = 0; d < ray_len; ++d)
				{
					xs += dx;
					ys += dy;
					xs = (xs < 0) ? xs + L : (xs > (L - 0.5)) ? xs - L + 0.5 : xs;
					ys = (ys < 0) ? ys + M : (ys > (M - 0.5)) ? ys - M + 0.5 : ys;

					iter_type x1 = std::lround(xs);
					iter_type y1 = std::lround(ys);

					iter_type n = x1 * L + y1;

					// When the sign flips, we need to start measuring or
					//	end measuring.
					bool plus0 = (data_y[n] > 0);
					if (plus0 != plus)
					{
						if (first_flipped)
						{
							double dx0 = startx - xs;
							double dy0 = starty - ys;
							lengths.push_back(std::sqrt(dx0 * dx0 + dy0 * dy0));
						}
						else
						{
							first_flipped = true;
						}

						plus = plus0;
						startx = xs;
						starty = ys;
					}
				}
			}

			std::sort(lengths.begin(), lengths.end());
			auto value = std::reduce(lengths.begin(), lengths.end() - lengths.size() / 2);
			all_lengths.push_back(value * (1.0 / (lengths.size() / 2)));
		}
	}

	scalar_t domain_length = 0;
	for (auto&& length : all_lengths)
	{
		domain_length += length;
	}

	return point_data<scalar_t>{ 0, domain_length * (1. / all_lengths.size()) };
}

ALGORITHM_POINT_3D_DECLARATION(DomainGrowth)
{
	scalar_t phase = 0;
	for (iter_type i = 0; i < len; ++i)
	{
		phase += (data_y[i] > 0) ? 1 : 0;
	}
	return point_data<Y>{ 0, phase * (1. / len) };
}


ALGORITHM_DYNAMIC_DECLARATION(DomainGrowth)
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



