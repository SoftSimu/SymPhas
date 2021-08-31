
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


#include "dft.h"


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
		auto [L, M, N] = symphas::lib::get_dimensions<3>(data_x, len);
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
				: extrema[1] - (data_x[i][0] + dx - avg[0]));
			axis_coord_t y = ((data_x[i][1] < avg[1])
				? avg[1] - data_x[i][1] - dy
				: extrema[3] - (data_x[i][1] + dy - avg[1]));
			axis_coord_t z = ((data_x[i][2] < avg[2])
				? avg[2] - data_x[i][2] - dz
				: extrema[5] - (data_x[i][2] + dz - avg[2]));

			dft_x[i] =
				axis_nd_t<3>{
					rx * x,
					ry * y,
					rz * z };
		}
	}

	void fill_x_axis(axis_nd_t<2>* dft_x, const axis_nd_t<2>* data_x, len_type len)
	{
		auto [L, M] = symphas::lib::get_dimensions<2>(data_x, len);
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
				: extrema[1] - (data_x[i][0] + dx - avg[0]));
			axis_coord_t y = ((data_x[i][1] < avg[1])
				? avg[1] - data_x[i][1] - dy
				: extrema[3] - (data_x[i][1] + dy - avg[1]));

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
		axis_coord_t adj[]{ avg - extrema[0], extrema[1] - avg };

		double rx = (2.0 * symphas::PI / (extrema[1] - extrema[0]));
		for (iter_type i = 0; i < len; ++i)
		{
			axis_coord_t x = ((data_x[i] < avg)
				? avg - data_x[i] - dx
				: extrema[1] - (data_x[i] + dx - avg));

			dft_x[i] =
				axis_nd_t<1>{ rx * x };
		}
	}


}


template<typename T>
inline void long_dft(T* data, complex_t* out, len_type len, bool backward = false)
{
	int const dL = static_cast<int>(len) / 2;
	for (iter_type ii = 0; ii < len; ++ii)
	{
		std::complex xx(0.0);
		double ti = 0, p_ti = (2.0 * symphas::PI / len) * (ii);
		for (iter_type i = 0; i < len; ++i, ti += p_ti)
		{
			out[ii] += data[i] * std::complex{ cos(ti), -sin(ti) };
		}
	}
}

template<typename T>
inline void long_dft(T* data, complex_t* out, iter_type L, iter_type M, bool backward = false)
{
	int sign = (backward) ? 1 : -1;

	len_type len = L * M;
	int const
		dL = L / 2,
		dM = M / 2;

	auto long_dft = [&](iter_type thr_i, len_type thr_n)
	{
		iter_type
			start_i = symphas::lib::next_block_i(thr_i, len, thr_n),
			end_i = symphas::lib::next_block_i(thr_i + 1, len, thr_n);
		int
			i0 = start_i % L,
			j0 = start_i / L;

		iter_type jj = j0, ii = i0;
		iter_type index = start_i;
		for (; jj < M; ++jj)
		{
			for (; ii < L; ++ii)
			{
				std::complex yy(0.0);
				double tj = 0, p_tj = (2.0 * symphas::PI / M) * (jj);
				for (iter_type j = 0; j < M; ++j, tj += p_tj)
				{
					std::complex xx(0.0);
					double ti = 0, p_ti = (2.0 * symphas::PI / L) * (ii);
					for (iter_type i = 0; i < L; ++i, ti += p_ti)
					{
						xx += data[i + j * L] * std::complex{ cos(ti), sign * sin(ti) };
					}
					yy += xx * std::complex{ cos(tj), sign * sin(tj) };
				}
				out[ii + jj * L] = yy;
				if (++index == end_i) return;
			}
			ii = 0;
		}
	};

#	pragma omp parallel for
	for (iter_type i = 0; i < L; ++i)
	{
		long_dft(i, L);
	}

}

template<typename T>
inline void long_dft(T* data, complex_t* out, iter_type L, iter_type M, iter_type N, bool backward = false)
{
	int sign = (backward) ? 1 : -1;

	len_type len = L * M * N;
	int const
		dL = L / 2 + ((L % 2 == 1) ? 1 : 0),
		dM = M / 2 + ((M % 2 == 1) ? 1 : 0),
		dN = N / 2 + ((N % 2 == 1) ? 1 : 0);

	auto long_dft = [&](iter_type thr_i, len_type thr_n)
	{
		iter_type
			start_i = symphas::lib::next_block_i(thr_i, len, thr_n),
			end_i = symphas::lib::next_block_i(thr_i + 1, len, thr_n);
		int
			i0 = (start_i % L),
			j0 = ((start_i / L) % M),
			k0 = (start_i / (L * M));


		iter_type kk = k0, jj = j0, ii = i0;
		iter_type index = start_i;
		for (; kk < N; ++kk)
		{
			for (; jj < M; ++jj)
			{
				for (; ii < L; ++ii)
				{
					std::complex zz(0.0);
					double p_tk = (2.0 * symphas::PI / N) * (kk), tk = 0;
					for (iter_type k = 0; k < N; ++k, tk += p_tk)
					{
						std::complex yy(0.0);
						double p_tj = (2.0 * symphas::PI / M) * (jj), tj = 0;

						for (iter_type j = 0; j < M; ++j, tj += p_tj)
						{
							std::complex xx(0.0);
							double p_ti = (2.0 * symphas::PI / L) * (ii), ti = 0;

							for (iter_type i = 0; i < L; ++i, ti += p_ti)
							{
								xx += data[i + j * L + k * L * M] * std::complex(cos(ti), sign * sin(ti));
							}
							yy += xx * std::complex(cos(tj), sign * sin(tj));
						}
						zz += yy * std::complex(cos(tk), sign * sin(tk));
					}
					out[ii + jj * L + kk * L * M] = zz;
					if (++index == end_i) return;
				}
				ii = 0;
			}
			jj = 0;
		}
	};

#	pragma omp parallel for
	for (iter_type i = 0; i < L; ++i)
	{
		long_dft(i, L);
	}
}




void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L, bool backward) { ::long_dft(data, out, L, backward); }
void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L, iter_type M, bool backward) { ::long_dft(data, out, L, M, backward); }
void symphas::dft::long_dft(scalar_t* data, complex_t* out, iter_type L, iter_type M, iter_type N, bool backward) { ::long_dft(data, out, L, M, N, backward); }

void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L, bool backward) { ::long_dft(data, out, L, backward); }
void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L, iter_type M, bool backward) { ::long_dft(data, out, L, M, backward); }
void symphas::dft::long_dft(complex_t* data, complex_t* out, iter_type L, iter_type M, iter_type N, bool backward) { ::long_dft(data, out, L, M, N, backward); }

