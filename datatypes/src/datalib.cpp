
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
 */


#include "datalib.h"


namespace symphas::dft
{
	void fill_x_axis(axis_nd_t<3>* dft_x, iter_type L, iter_type M, iter_type N)
	{
		iter_type
			dL = static_cast<iter_type>(L / 2),
			dM = static_cast<iter_type>(M / 2),
			dN = static_cast<iter_type>(N / 2);
		iter_type n = 0;
		for (iter_type k = 0; k < N; ++k)
		{
			for (iter_type j = 0; j < M; ++j)
			{
				for (iter_type i = 0; i < L; ++i)
				{
					dft_x[n++] =
						axis_nd_t<3>{
							(2.0 * symphas::PI / L) * ((i <= dL) ? ((dL - i)) : ((L - (i - dL)))),
							(2.0 * symphas::PI / M) * ((j <= dM) ? ((dM - j)) : ((M - (j - dM)))),
							(2.0 * symphas::PI / N) * ((k <= dN) ? ((dN - k)) : ((N - (k - dN)))) };
				}
			}
		}
	}

	void fill_x_axis(axis_nd_t<2>* dft_x, iter_type L, iter_type M)
	{
		iter_type
			dL = static_cast<iter_type>(L / 2),
			dM = static_cast<iter_type>(M / 2);
		iter_type n = 0;
		for (iter_type j = 0; j < M; ++j)
		{
			for (iter_type i = 0; i < L; ++i)
			{
				dft_x[n++] =
					axis_nd_t<2>{
						(2.0 * symphas::PI / L) * ((i <= dL) ? ((dL - i)) : ((L - (i - dL)))),
						(2.0 * symphas::PI / M) * ((j <= dM) ? ((dM - j)) : ((M - (j - dM)))) };
			}
		}
	}

	void fill_x_axis(axis_nd_t<1>* dft_x, iter_type L)
	{
		iter_type
			dL = static_cast<iter_type>(L / 2);
		iter_type n = 0;
		for (iter_type i = 0; i < L; ++i)
		{
			dft_x[n++] =
				axis_nd_t<1>{ (2.0 * symphas::PI / L) * ((i <= dL) ? ((dL - i)) : ((L - (i - dL)))) };
		}
	}

	void fill_x_axis(axis_nd_t<3>* dft_x, const axis_nd_t<3>* data_x, len_type len)
	{
		auto [L, M, N] = symphas::lib::get_dimensions<3>(data_x, len)._3();
		axis_coord_t extrema[6]{ +INFINITY, -INFINITY, +INFINITY, -INFINITY, +INFINITY, -INFINITY };
		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < extrema[0])
			{
				extrema[0] = data_x[i][0];
			}
			if (data_x[i][0] > extrema[1])
			{
				extrema[1] = data_x[i][0];
			}
			if (data_x[i][1] < extrema[2])
			{
				extrema[2] = data_x[i][1];
			}
			if (data_x[i][1] > extrema[3])
			{
				extrema[3] = data_x[i][1];
			}
			if (data_x[i][2] < extrema[4])
			{
				extrema[4] = data_x[i][2];
			}
			if (data_x[i][2] > extrema[5])
			{
				extrema[5] = data_x[i][2];
			}
		}


		axis_coord_t
			dx = ((*(data_x + 1))[0] - (*data_x)[0]),
			dy = ((*(data_x + L))[1] - (*data_x)[1]),
			dz = ((*(data_x + L * M))[2] - (*data_x)[2]);

		extrema[1] += dx;
		extrema[3] += dy;
		extrema[5] += dz;

		axis_coord_t avg[]{
			(extrema[1] + extrema[0]) / 2,
			(extrema[3] + extrema[2]) / 2,
			(extrema[5] + extrema[4]) / 2 };

		double rx = (2.0 * symphas::PI / ((extrema[1] - extrema[0]) * dx));
		double ry = (2.0 * symphas::PI / ((extrema[3] - extrema[2]) * dy));
		double rz = (2.0 * symphas::PI / ((extrema[5] - extrema[4]) * dz));
		for (iter_type i = 0; i < len; ++i)
		{
			axis_coord_t x = ((data_x[i][0] < avg[0])
				? avg[0] - data_x[i][0] - dx
				: (extrema[1] - (data_x[i][0] + dx - avg[0])));
			axis_coord_t y = ((data_x[i][1] < avg[1])
				? avg[1] - data_x[i][1] - dy
				: (extrema[3] - (data_x[i][1] + dy - avg[1])));
			axis_coord_t z = ((data_x[i][2] < avg[2])
				? avg[2] - data_x[i][2] - dz
				: (extrema[5] - (data_x[i][2] + dz - avg[2])));

			dft_x[i] =
				axis_nd_t<3>{
					rx * x,
					ry * y,
					rz * z };
		}
	}

	void fill_x_axis(axis_nd_t<2>* dft_x, const axis_nd_t<2>* data_x, len_type len)
	{
		auto [L, M] = symphas::lib::get_dimensions<2>(data_x, len)._2();
		axis_coord_t extrema[4]{ +INFINITY, -INFINITY, +INFINITY, -INFINITY };
		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < extrema[0])
			{
				extrema[0] = data_x[i][0];
			}
			if (data_x[i][0] > extrema[1])
			{
				extrema[1] = data_x[i][0];
			}
			if (data_x[i][1] < extrema[2])
			{
				extrema[2] = data_x[i][1];
			}
			if (data_x[i][1] > extrema[3])
			{
				extrema[3] = data_x[i][1];
			}
		}

		axis_coord_t
			dx = (extrema[1] - extrema[0]) / (L - 1),
			dy = (extrema[3] - extrema[2]) / (M - 1);

		extrema[1] += dx;
		extrema[3] += dy;

		axis_coord_t avg[]{
			(extrema[1] + extrema[0]) / 2,
			(extrema[3] + extrema[2]) / 2 };

		double rx = (2.0 * symphas::PI / ((extrema[1] - extrema[0]) * dx));
		double ry = (2.0 * symphas::PI / ((extrema[3] - extrema[2]) * dy));
		for (iter_type i = 0; i < len; ++i)
		{
			axis_coord_t x = ((data_x[i][0] < avg[0])
				? avg[0] - data_x[i][0] - dx
				: (extrema[1] - (data_x[i][0] + dx - avg[0])));
			axis_coord_t y = ((data_x[i][1] < avg[1])
				? avg[1] - data_x[i][1] - dy
				: (extrema[3] - (data_x[i][1] + dy - avg[1])));

			dft_x[i] =
				axis_nd_t<2>{
					rx * x,
					ry * y };
		}
	}

	void fill_x_axis(axis_nd_t<1>* dft_x, const axis_nd_t<1>* data_x, len_type len)
	{
		axis_coord_t extrema[2];
		symphas::lib::fill_sorted_ranges(data_x, len, extrema);
		axis_coord_t dx = (*(data_x + 1)) - (*data_x);
		extrema[1] += dx;

		axis_coord_t avg = (extrema[1] + extrema[0]) / 2;

		double rx = (2.0 * symphas::PI / (extrema[1] - extrema[0]));
		for (iter_type i = 0; i < len; ++i)
		{
			axis_coord_t x = ((data_x[i] < avg)
				? avg - data_x[i] - dx
				: (extrema[1] - (data_x[i] + dx - avg)));

			dft_x[i] =
				axis_nd_t<1>{ rx * x };
		}
	}


	void fill_x_axis(axis_nd_t<3>* dft_x, const axis_coord_t(&x)[2], len_type L,
		const axis_coord_t(&y)[2], len_type M, const axis_coord_t(&z)[2], len_type N)
	{
		axis_coord_t
			dx = (x[1] - x[0]) / (L - 1),
			dy = (y[1] - y[0]) / (M - 1),
			dz = (z[1] - z[0]) / (N - 1);

		axis_coord_t avg[]{
			(x[1] + dx + x[0]) / 2,
			(y[1] + dy + y[0]) / 2,
			(z[1] + dz + z[0]) / 2 };

		double rx = (2.0 * symphas::PI / ((x[1] + dx - x[0]) * dx));
		double ry = (2.0 * symphas::PI / ((y[1] + dy - y[0]) * dy));
		double rz = (2.0 * symphas::PI / ((z[1] + dz - z[0]) * dz));

		for (iter_type i = 0; i < L; ++i)
		{
			for (iter_type j = 0; j < M; ++j)
			{
				for (iter_type k = 0; k < N; ++k)
				{
					axis_nd_t<3> data{ i * dx + x[0], j * dy + y[0], k * dz + z[0] };

					axis_coord_t x0 = ((data[0] < avg[0])
						? avg[0] - data[0] - dx
						: (x[1] + dx - (data[0] + dx - avg[0])));
					axis_coord_t y0 = ((data[1] < avg[1])
						? avg[1] - data[1] - dy
						: (y[1] + dy - (data[1] + dy - avg[1])));
					axis_coord_t z0 = ((data[2] < avg[2])
						? avg[2] - data[2] - dz
						: (z[1] + dz - (data[2] + dz - avg[2])));

					dft_x[i + j * L + k * L * M] =
						axis_nd_t<3>{
							rx * x0,
							ry * y0,
							rz * z0 };
				}
			}
		}
	}

	void fill_x_axis(axis_nd_t<2>* dft_x, const axis_coord_t(&x)[2], len_type L, const axis_coord_t(&y)[2], len_type M)
	{
		axis_coord_t
			dx = (x[1] - x[0]) / (L - 1),
			dy = (y[1] - y[0]) / (M - 1);

		axis_coord_t avg[]{
			((x[1] + dx) + x[0]) / 2,
			((y[1] + dy) + y[0]) / 2 };

		double rx = (2.0 * symphas::PI / ((x[1] + dx - x[0]) * dx));
		double ry = (2.0 * symphas::PI / ((y[1] + dy - y[0]) * dy));

		for (iter_type i = 0; i < L; ++i)
		{
			for (iter_type j = 0; j < M; ++j)
			{
				axis_nd_t<2> data{ i * dx + x[0], j * dy + y[0] };

				axis_coord_t x0 = ((data[0] < avg[0])
					? avg[0] - data[0] - dx
					: (x[1] + dx - (data[0] + dx - avg[0])));
				axis_coord_t y0 = ((data[1] < avg[1])
					? avg[1] - data[1] - dy
					: (y[1] + dy - (data[1] + dy - avg[1])));

				dft_x[i + j * L] =
					axis_nd_t<2>{
						rx * x0,
						ry * y0 };
			}
		}
	}

	void fill_x_axis(axis_nd_t<1>* dft_x, const axis_coord_t(&x)[2], len_type len)
	{
		axis_coord_t dx = (x[1] - x[0]) / (len - 1);
		axis_coord_t avg = ((x[1] + dx) + x[0]) / 2;
		//axis_coord_t adj[]{ avg - x[0], (x[1] + dx) - avg };

		double rx = (2.0 * symphas::PI / ((x[1] + dx) - x[0]));
		for (iter_type i = 0; i < len; ++i)
		{
			scalar_t data_x = i * dx + x[0];
			axis_coord_t x0 = ((data_x < avg)
				? avg - data_x - dx
				: (x[1] + dx - (data_x + dx - avg)));

			dft_x[i] =
				axis_nd_t<1>{ rx * x0 };
		}
	}


}
