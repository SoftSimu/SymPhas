#pragma once

#include "symphasthread.h"
#include "processstatic.h"

DEFINE_ALGORITHM_SCALAR(PCF)
DEFINE_ALGORITHM_VECTOR(PCF)
DEFINE_DATA(PCF, SCALAR, ALG_SCALAR, ALG_VECTOR)

/*
 * runs the function 'add' over all points
 */
template<typename T>
auto pcf_bins_1d(const T* data, len_type LL, len_type thr_i, len_type thr_n)
{
	iter_type
		n = symphas::lib::next_block_i(thr_i, LL, thr_n),
		e = symphas::lib::next_block_i(thr_i + 1, LL, thr_n);
	iter_type
		L = LL,
		L2 = L / 2,
		x = n;

	len_type dim_half[] = { LL / 2 };

	len_type len = symphas::lib::radial_length(dim_half);
	iter_type* dmap = symphas::lib::make_radial_index_map(dim_half);
	

	T* pcf_values = new T[len]{ 0 };
	for (; x < L && n < e; ++x, ++n)
	{
		int nn, ii;
		iter_type end_x;

		// nn is the index in the full array
		nn = n;

		// ii is the number of indices in the half square 
		// (that the point correlation is being measured inside)
		// required to get from n to nn
		ii = 0;

		// we measure the point correlation L2 units away, the x length
		// of the half square

		// we first go up the boundary
		end_x = std::min(x + L2, L);
		for (iter_type xx = x; xx < end_x; ++xx)
		{
			pcf_values[dmap[ii++]] += data[n] * data[nn++];
		}

		// and continue counting past the boundary
		// end_x could be either at the boundary or if the boundary wasn't reached,
		// the result will be 0
		end_x = (x + L2) - end_x;
		for (iter_type xx = 0; xx < end_x; ++xx)
		{
			pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
		}

	}

	delete[] dmap;
	return pcf_values;

}

template<typename T>
auto pcf_bins_2d(const T* data, len_type LL, len_type MM, len_type thr_i, len_type thr_n)
{
	iter_type
		n = symphas::lib::next_block_i(thr_i, LL, thr_n),
		e = symphas::lib::next_block_i(thr_i + 1, LL, thr_n);
	iter_type
		L = LL,
		M = MM,
		L2 = L / 2,
		M2 = M / 2,
		x = n % L,
		y = n / L;

	iter_type 
		cx = (L % 2 == 0) ? 0 : 1;

	len_type dim_half[] = { LL / 2, MM / 2 };

	len_type len = symphas::lib::radial_length(dim_half);
	iter_type* dmap = symphas::lib::make_radial_index_map(dim_half);

	T* pcf_values = new T[len]{ 0 };
	for (; y < M && n < e; ++y)
	{
		for (; x < L && n < e; ++x, ++n)
		{
			int nn, ii;
			iter_type end_x, end_y;

			// nn is the index in the full array
			nn = n;

			// ii is the number of indices in the half square 
			// (that the point correlation is being measured inside)
			// required to get from n to nn
			ii = 0;

			// we measure the point correlation M2 units away in the y direction,
			// the y length of the half square
			end_y = std::min(y + M2, M);
			for (iter_type yy = y; yy < end_y; ++yy)
			{
				// we measure the point correlation L2 units away, the x length
				// of the half square

				// we first go up the boundary
				end_x = std::min(x + L2, L);
				for (iter_type xx = x; xx < end_x; ++xx)
				{
					pcf_values[dmap[ii++]] += data[n] * data[nn++];
				}

				// and continue counting past the boundary
				// end_x could be either at the boundary or if the boundary wasn't reached,
				// the result will be 0
				end_x = (x + L2) - end_x;
				for (iter_type xx = 0; xx < end_x; ++xx)
				{
					pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
				}

				// move nn to the next row up
				nn += L2 + cx;
			}


			end_y = (y + M2) - end_y;
			nn -= L * M;
			for (iter_type yy = 0; yy < end_y; ++yy)
			{
				// we measure the point correlation L2 units away, the x length
				// of the half square

				// we first go up the boundary
				end_x = std::min(x + L2, L);
				for (iter_type xx = x; xx < end_x; ++xx)
				{
					pcf_values[dmap[ii++]] += data[n] * data[nn++];
				}

				// and continue counting past the boundary
				// end_x could be either at the boundary or if the boundary wasn't reached,
				// the result will be 0
				end_x = (x + L2) - end_x;
				for (iter_type xx = 0; xx < end_x; ++xx)
				{
					pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
				}

				// move nn to the next row up
				nn += L2 + cx;
			}
		}
		x = 0;
	}

	delete[] dmap;
	return pcf_values;

}

template<typename T>
auto pcf_bins_3d(const T* data, len_type LL, len_type MM, len_type NN, iter_type thr_i, len_type thr_n)
{
	iter_type
		n = symphas::lib::next_block_i(thr_i, LL, thr_n),
		e = symphas::lib::next_block_i(thr_i + 1, LL, thr_n);
	int
		L = static_cast<int>(LL),
		M = static_cast<int>(MM),
		N = static_cast<int>(NN),
		L2 = L / 2,
		M2 = M / 2,
		N2 = N / 2,
		x = static_cast<int>(n) % L,
		y = static_cast<int>(n / L) % M,
		z = static_cast<int>(n) / (L * M);

	len_type dim_half[] = { LL / 2, MM / 2, NN / 2 };

	len_type len = symphas::lib::radial_length(dim_half);
	iter_type* dmap = symphas::lib::make_radial_index_map(dim_half);

	T* pcf_values = new T[len]{ 0 };
	for (; z < N && n < e; ++z)
	{
		for (; y < M && n < e; ++y)
		{
			for (; x < L && n < e; ++x, ++n)
			{
				int nn, ii;
				iter_type end_x, end_y, end_z;

				// nn is the index in the full array
				nn = n;

				// ii is the number of indices in the half square 
				// (that the point correlation is being measured inside)
				// required to get from n to nn
				ii = 0;

				end_z = std::min(z + N2, N);
				for (iter_type zz = z; zz < end_z; ++zz)
				{
					// we measure the point correlation M2 units away in the y direction,
					// the y length of the half square
					end_y = std::min(y + M2, M);
					for (iter_type yy = y; yy < end_y; ++yy)
					{
						// we measure the point correlation L2 units away, the x length
						// of the half square

						// we first go up the boundary
						end_x = std::min(x + L2, L);
						for (iter_type xx = x; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++];
						}

						// and continue counting past the boundary
						// end_x could be either at the boundary or if the boundary wasn't reached,
						// the result will be 0
						end_x = (x + L2) - end_x;
						for (iter_type xx = 0; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
						}

						// move nn to the next row up
						nn += L2;
					}


					end_y = (y + M2) - end_y;
					nn -= L * M;
					for (iter_type yy = 0; yy < end_y; ++yy)
					{
						// we measure the point correlation L2 units away, the x length
						// of the half square

						// we first go up the boundary
						end_x = std::min(x + L2, L);
						for (iter_type xx = x; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++];
						}

						// and continue counting past the boundary
						// end_x could be either at the boundary or if the boundary wasn't reached,
						// the result will be 0
						end_x = (x + L2) - end_x;
						for (iter_type xx = 0; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
						}

						// move nn to the next row up
						nn += L2;
					}
					nn += L * M2;
				}

				end_z = (z + N2) - end_z;
				nn -= L * M * N2;
				for (iter_type zz = 0; zz < end_z; ++zz)
				{
					// we measure the point correlation M2 units away in the y direction,
					// the y length of the half square
					end_y = std::min(y + M2, M);
					for (iter_type yy = y; yy < end_y; ++yy)
					{
						// we measure the point correlation L2 units away, the x length
						// of the half square

						// we first go up the boundary
						end_x = std::min(x + L2, L);
						for (iter_type xx = x; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++];
						}

						// and continue counting past the boundary
						// end_x could be either at the boundary or if the boundary wasn't reached,
						// the result will be 0
						end_x = (x + L2) - end_x;
						for (iter_type xx = 0; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
						}

						// move nn to the next row up
						nn += L2;
					}


					end_y = (y + M2) - end_y;
					nn -= L * M;
					for (iter_type yy = 0; yy < end_y; ++yy)
					{
						// we measure the point correlation L2 units away, the x length
						// of the half square

						// we first go up the boundary
						end_x = std::min(x + L2, L);
						for (iter_type xx = x; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++];
						}

						// and continue counting past the boundary
						// end_x could be either at the boundary or if the boundary wasn't reached,
						// the result will be 0
						end_x = (x + L2) - end_x;
						for (iter_type xx = 0; xx < end_x; ++xx)
						{
							pcf_values[dmap[ii++]] += data[n] * data[nn++ - L];
						}

						// move nn to the next row up
						nn += L2;
					}
					nn += L * M2;
				}
			}
			x = 0;
		}
		y = 0;
	}

	delete[] dmap;
	return pcf_values;

}


#define OMP_BLOCK_COUNT 32

//! Returns the 2-point correlation function.
/*!
 * Computes the radial distribution function histogram
 */
ALGORITHM_SCALAR_1D_DECLARATION(PCF)
{
	len_type dim_half[] = { L / 2 };
	len_type out_len = symphas::lib::radial_length(dim_half);

	double const
		dx = ((*(data_x + 1)) - (*data_x));

	auto [pcf_x, pcf_c] = symphas::lib::make_radial_arrays(dim_half, dx);
	Y* pcf_y = new Y[out_len]{ 0 };

	len_type block_count = std::min(L, OMP_BLOCK_COUNT);
	Y** pcf_aggr = new Y*[block_count];
#	pragma omp parallel for
	for (iter_type n = 0; n < block_count; ++n)
	{
		// run the threads
		pcf_aggr[n] = pcf_bins_1d<Y>(data_y, L, n, block_count);
	}
	for (iter_type n = 0; n < block_count; ++n)
	{
		for (iter_type i = 0; i < out_len; ++i)
		{
			pcf_y[i] += pcf_aggr[n][i];
		}
		delete[] pcf_aggr[n];
	}
	delete[] pcf_aggr;

	for (iter_type i = 0; i < out_len; ++i)
	{
		pcf_y[i] = pcf_y[i] * (1.0 / (pcf_c[i] * len));
	}

	return scalar_data<Y>{ pcf_x, pcf_y, out_len };
}

ALGORITHM_SCALAR_2D_DECLARATION(PCF)
{
	len_type dim_half[] = { L / 2, M / 2 };
	len_type out_len = symphas::lib::radial_length(dim_half);

	double const
		dx = ((*(data_x + 1))[0] - (*data_x)[0]),
		dy = ((*(data_x + L))[1] - (*data_x)[1]);

	auto [pcf_x, pcf_c] = symphas::lib::make_radial_arrays(dim_half, dx, dy);
	Y* pcf_y = new Y[out_len]{ 0 };

	
	len_type block_count = std::min(L, OMP_BLOCK_COUNT);
	Y** pcf_aggr = new Y*[block_count];
#	pragma omp parallel for
	for (iter_type n = 0; n < block_count; ++n)
	{
		// run the threads
		pcf_aggr[n] = pcf_bins_2d<Y>(data_y, L, M, n, block_count);
	}
	for (iter_type n = 0; n < block_count; ++n)
	{
		for (iter_type i = 0; i < out_len; ++i)
		{
			pcf_y[i] += pcf_aggr[n][i];
		}
		delete[] pcf_aggr[n];
	}
	delete[] pcf_aggr;

	double r = (1.0 / out_len);
	for (iter_type i = 0; i < out_len; ++i)
	{
		pcf_y[i] = pcf_y[i] * (r / pcf_c[i]);
	}

	return scalar_data<Y>{ pcf_x, pcf_y, out_len };
	//return symphas::lib::radial_avg(RUN_ALGORITHM_VECTOR_2D(PCF));
}

ALGORITHM_SCALAR_3D_DECLARATION(PCF)
{
	len_type dim_half[] = { L / 2, M / 2, N / 2 };
	len_type out_len = symphas::lib::radial_length(dim_half);

	double const
		dx = ((*(data_x + 1))[0] - (*data_x)[0]),
		dy = ((*(data_x + L))[1] - (*data_x)[1]),
		dz = ((*(data_x + L * M))[2] - (*data_x)[2]);

	auto [pcf_x, pcf_c] = symphas::lib::make_radial_arrays(dim_half, dx, dy, dz);
	Y* pcf_y = new Y[out_len]{ 0 };

	len_type block_count = std::min(L, OMP_BLOCK_COUNT);
	Y** pcf_aggr = new Y * [block_count];
#	pragma omp parallel for
	for (iter_type n = 0; n < block_count; ++n)
	{
		// run the threads
		pcf_aggr[n] = pcf_bins_3d<Y>(data_y, L, M, N, n, block_count);
	}
	for (iter_type n = 0; n < block_count; ++n)
	{
		for (iter_type i = 0; i < out_len; ++i)
		{
			pcf_y[i] += pcf_aggr[n][i];
		}
		delete[] pcf_aggr[n];
	}
	delete[] pcf_aggr;

	for (iter_type i = 0; i < out_len; ++i)
	{
		pcf_y[i] = pcf_y[i] * (1.0 / (pcf_c[i] * len));
	}

	return scalar_data<Y>{ pcf_x, pcf_y, out_len };
}

/*
 * computes the point correlation function data
 */

ALGORITHM_VECTOR_1D_DECLARATION(PCF)
{
	Y* put_main_y = new Y[len];
	axis_1d_type* put_main_x = new axis_1d_type[len];

	/*
	 * here we store the length of each of the computed arrays containing
	 * the 2-point correlation of a subset of vectors. This is used to
	 * fill the vector with all values.
	 */
	double const dx = (data_x[1] - data_x[0]);
	iter_type const LL = static_cast<iter_type>(std::ceil(L / 2.0));

	auto corr_m = [&]()
	{
		//Time t("main thread 2-point correlation");
		auto put_y = put_main_y;
		auto put_x = put_main_x;

		put_y = new Y[len];
		put_x = new axis_1d_type[len];
		iter_type index = 0;

		Y sum;
		iter_type nn;

		/*
		 * count the interaction of all points with themselves
		 */
		sum = 0;
		nn = 0;
		for (iter_type xx = 0; xx < L; ++xx, ++nn)
		{
			sum += data_y[nn] * data_y[nn];
		}
		put_y[index] = (1.0 / nn) * sum;
		put_x[index] = axis_1d_type{ 0 };
		++index;


		/*
		 * count all the interactions on the x-axis
		 */
		for (iter_type x = 1; x < LL + 1; ++x)
		{
			sum = 0;
			nn = 0;
			for (iter_type xx = 0; xx < L - x; ++xx, ++nn)
			{
				sum += data_y[nn] * data_y[nn + x];
			}
			for (iter_type xx = L - x; xx < L; ++xx, ++nn)
			{
				sum += data_y[nn] * data_y[nn + x - L];
			}
			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_1d_type{ dx * x };
			++index;
			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_1d_type{ -dx * x };
			++index;
		}
	};

	corr_m();
	len_type out_len = L;

	symphas::lib::sort_data(put_main_x, put_main_y, out_len);
	return vector_data<Y, 1>{ put_main_x, put_main_y, out_len };

}

///*!
// * \param put_x \f$x\f$-axis array.
// * \param put_y \f$y\f$-axis array.
// * \param d The list of spatial widths.
// * \param s The index separations list.
// * \param index The current index into the lists to put.
// * \param
// */
//template<size_t U, typename Y>
//void insert_sum(axis_2d_type &put_x, Y &put_y, double* d, iter_type *s, double y_value)
//{
//
//}

ALGORITHM_VECTOR_2D_DECLARATION(PCF)
{
	Y* put_main_y = new Y[len];
	axis_2d_type* put_main_x = new axis_2d_type[len];

	len_type corr_0_len = L;

	//Time t("point correlation (2d)");

	double const
		dx = ((*(data_x + 1))[0] - (*data_x)[0]),
		dy = ((*(data_x + L))[1] - (*data_x)[1]);
	iter_type const
		LL = static_cast<iter_type>(std::floor(L / 2.0)),
		MM = static_cast<iter_type>(std::floor(M / 2.0));

	auto corr_m = [&]()
	{
		auto* put_y = put_main_y;
		auto* put_x = put_main_x;
		iter_type index = 0;

		Y sum;
		iter_type nn;

		/*
		 * count the interaction of all points with themselves
		 */
		sum = 0;
		nn = 0;
		for (iter_type yy = 0; yy < M; ++yy)
		{
			for (iter_type xx = 0; xx < L; ++xx, ++nn)
			{
				sum += data_y[nn] * data_y[nn];
			}
		}
		put_y[index] = (1.0 / nn) * sum;
		put_x[index] = axis_2d_type{ 0, 0 };
		++index;


		/*
		 * count all the interactions on the x-axis
		 */
		for (iter_type x = 1; x < LL + 1; ++x)
		{
			sum = 0;
			nn = 0;
			for (iter_type yy = 0; yy < M; ++yy)
			{
				for (iter_type xx = 0; xx < L - x; ++xx, ++nn)
				{
					sum += data_y[nn] * data_y[nn + x];
				}
				for (iter_type xx = L - x; xx < L; ++xx, ++nn)
				{
					sum += data_y[nn] * data_y[nn + x - L];
				}
			}
			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_2d_type{ -dx * x, 0 };
			++index;

			if (x < LL)
			{
				put_y[index] = (1.0 / nn) * sum;
				put_x[index] = axis_2d_type{ dx * x, 0 };
				++index;
			}
		}
	};

	auto corr = [&](iter_type thr_i, len_type thr_n)
	{
		iter_type
			i = symphas::lib::next_block_i(thr_i, MM, thr_n) + 1,
			e = symphas::lib::next_block_i(thr_i + 1, MM, thr_n) + 1;

		len_type offset_ptr = (i - 1) * (2 + 2 * LL + 2 * (LL - 1));
		auto* put_y = put_main_y + corr_0_len + offset_ptr;
		auto* put_x = put_main_x + corr_0_len + offset_ptr;

		Y sum;
		iter_type index = 0;
		iter_type nn, mm;

		for (iter_type y = i; y < e; ++y)
		{
			iter_type n = y * L;

			/*
			 * count all the interactions on the y-axis
			 */
			sum = 0;
			nn = 0;
			for (iter_type yy = 0; yy < M - y; ++yy)
			{
				for (iter_type xx = 0; xx < L; ++xx, ++nn)
				{
					sum += data_y[nn] * data_y[nn + n];
				}
			}
			for (iter_type yy = M - y; yy < M; ++yy)
			{
				for (iter_type xx = 0; xx < L; ++xx, ++nn)
				{
					sum += data_y[nn] * data_y[nn + n - len];
				}
			}
			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_2d_type{ 0, -dy * y };
			++index;

			if (y < MM)
			{
				put_y[index] = (1.0 / nn) * sum;
				put_x[index] = axis_2d_type{ 0, dy * y };
				++index;
			}


			iter_type m;
			for (iter_type x = 1; x < LL + 1; ++x)
			{
				Y sum_a = 0, sum_b = 0;

				n = y * L + x;
				m = y * L - x;
				nn = 0;
				mm = x;

				for (iter_type yy = 0; yy < M - y; ++yy)
				{
					for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
					{
						sum_a += data_y[nn] * data_y[nn + n];
						sum_b += data_y[mm] * data_y[mm + m];
					}
					for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
					{
						sum_a += data_y[nn] * data_y[nn + n - L];
						sum_b += data_y[mm - L] * data_y[mm + m];
					}
				}

				for (iter_type yy = M - y; yy < M; ++yy)
				{
					for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
					{
						sum_a += data_y[nn] * data_y[nn + n - len];
						sum_b += data_y[mm] * data_y[mm + m - len];
					}
					for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
					{
						sum_a += data_y[nn] * data_y[nn + n - L - len];
						sum_b += data_y[mm - L] * data_y[mm + m - len];
					}
				}
				put_y[index] = (1.0 / nn) * sum_a;
				put_x[index] = axis_2d_type{ -dx * x, -dy * y };
				++index;

				if (y < MM)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_2d_type{ -dx * x, dy * y };
					++index;
				}

				if (x < LL)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_2d_type{ dx * x, -dy * y };
					++index;

					if (y < MM)
					{
						put_y[index] = (1.0 / nn) * sum_a;
						put_x[index] = axis_2d_type{ dx * x, dy * y };
						++index;
					}
				}
			}
		}
	};

	len_type block_count = std::min(MM / 2, OMP_BLOCK_COUNT);
#	pragma omp parallel for
	for (iter_type n = 0; n < block_count; ++n)
	{
		if (n == 0)
		{
			corr_m();
		}
		corr(n, block_count);
	}


	symphas::lib::sort_data(put_main_x, put_main_y, L, M);
	return vector_data<Y, 2>{ put_main_x, put_main_y, len };

}

ALGORITHM_VECTOR_3D_DECLARATION(PCF)
{

	len_type corr_0_len = L * M;

	Y* put_main_y = new Y[len];
	axis_3d_type* put_main_x = new axis_3d_type[len];

	double
		dx = ((*(data_x + M * N))[0] - (*data_x)[0]),
		dy = ((*(data_x + N))[1] - (*data_x)[1]),
		dz = ((*(data_x + 1))[2] - (*data_x)[2]);

	iter_type LL = L / 2, MM = M / 2, NN = N / 2;

	auto corr_m = [&]()
	{
		auto* put_y = put_main_y;
		auto* put_x = put_main_x;
		iter_type index = 0;

		Y sum;
		iter_type nn, mm;

		/*
		 * count the interaction of all points with themselves
		 */
		sum = 0;
		nn = 0;
		for (iter_type zz = 0; zz < N; ++zz)
		{
			for (iter_type yy = 0; yy < M; ++yy)
			{
				for (iter_type xx = 0; xx < L; ++xx, ++nn)
				{
					sum += data_y[nn] * data_y[nn];
				}
			}
		}

		put_y[index] = (1.0 / nn) * sum;
		put_x[index] = axis_3d_type{ 0, 0, 0 };
		++index;

		/*
		 * count all the interactions on the x-axis
		 */
		for (iter_type x = 1; x < LL + 1; ++x)
		{
			sum = 0;
			nn = 0;
			for (iter_type zz = 0; zz < N; ++zz)
			{
				for (iter_type yy = 0; yy < M; ++yy)
				{
					for (iter_type xx = 0; xx < L - x; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + x];
					}
					for (iter_type xx = L - x; xx < L; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + x - L];
					}
				}
			}

			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_3d_type{ -dx * x, 0, 0 };
			++index;
			
			if (x < LL)
			{
				put_y[index] = (1.0 / nn) * sum;
				put_x[index] = axis_3d_type{ dx * x, 0, 0 };
				++index;
			}
		}

		/*
		 * count all the interactions on the y-axis
		 */
		for (iter_type y = 1, n = L; y < MM + 1; ++y, n = L * y)
		{
			sum = 0;
			nn = 0;
			for (iter_type zz = 0; zz < N; ++zz)
			{
				for (iter_type yy = 0; yy < M - y; ++yy)
				{
					for (iter_type xx = 0; xx < L; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + n];
					}
				}
				for (iter_type yy = M - y; yy < M; ++yy)
				{
					for (iter_type xx = 0; xx < L; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + n - L * M];
					}
				}
			}

			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_3d_type{ 0, -dy * y, 0 };
			++index;
			if (y < MM)
			{
				put_y[index] = (1.0 / nn) * sum;
				put_x[index] = axis_3d_type{ 0, dy * y, 0 };
				++index;
			}

			iter_type m;
			for (iter_type x = 1; x < LL + 1; ++x, ++n)
			{
				Y sum_a = 0, sum_b = 0;

				n = y * L + x;
				m = y * L - x;
				nn = 0;
				mm = x;
				for (iter_type zz = 0; zz < N; ++zz)
				{
					for (iter_type yy = 0; yy < M - y; ++yy)
					{
						for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n];
							sum_b += data_y[mm] * data_y[mm + m];
						}
						for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L];
							sum_b += data_y[mm - L] * data_y[mm + m];
						}
					}
					for (iter_type yy = M - y; yy < M; ++yy)
					{
						for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L * M];
							sum_b += data_y[mm] * data_y[mm + m - L * M];
						}
						for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L - L * M];
							sum_b += data_y[mm - L] * data_y[mm + m - L * M];
						}
					}
				}

				put_y[index] = (1.0 / nn) * sum_a;
				put_x[index] = axis_3d_type{ -dx * x, -dy * y, 0 };
				++index;
				if (y < MM)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ -dx * x, dy * y, 0 };
					++index;
				}
				if (x < LL)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ dx * x, -dy * y, 0 };
					++index;
					if (y < MM)
					{
						put_y[index] = (1.0 / nn) * sum_a;
						put_x[index] = axis_3d_type{ dx * x, dy * y, 0 };
						++index;
					}
				}
			}
		}
	};

	auto corr = [&](iter_type thr_i, len_type thr_n)
	{
		iter_type
			i = symphas::lib::next_block_i(thr_i, NN, thr_n) + 1,
			e = symphas::lib::next_block_i(thr_i + 1, NN, thr_n) + 1;

		len_type offset_ptr = 
			(i - 1) * (2 + 
				2 * LL + 2 * (LL - 1) + 
				2 * MM + 2 * (MM - 1) +
				2 * MM * LL + 2 * (MM - 1) * LL + 2 * MM * (LL - 1) + 2 * (MM - 1) * (LL - 1));
		auto* put_y = put_main_y + corr_0_len + offset_ptr;
		auto* put_x = put_main_x + corr_0_len + offset_ptr;

		Y sum;
		iter_type index = 0;
		iter_type nn, mm;

		for (iter_type z = i; z < e; ++z)
		{
			iter_type n = z * L * M;
			/*
			 * count all the interactions on the z axis
			 */
			sum = 0;
			nn = 0;
			for (iter_type zz = 0; zz < N - z; ++zz)
			{
				for (iter_type yy = 0; yy < M; ++yy)
				{
					for (iter_type xx = 0; xx < L; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + n];
					}
				}
			}
			for (iter_type zz = N - z; zz < N; ++zz)
			{
				for (iter_type yy = 0; yy < M; ++yy)
				{
					for (iter_type xx = 0; xx < L; ++xx, ++nn)
					{
						sum += data_y[nn] * data_y[nn + n - len];
					}
				}
			}

			put_y[index] = (1.0 / nn) * sum;
			put_x[index] = axis_3d_type{ 0, 0, -dz * z };
			++index;
			if (z < NN)
			{
				put_y[index] = (1.0 / nn) * sum;
				put_x[index] = axis_3d_type{ 0, 0, dz * z };
				++index;
			}

			/*
			 * count all the interactions on the x-axis
			 */
			iter_type m;
			for (iter_type x = 1; x < LL + 1; ++x)
			{
				Y sum_a = 0, sum_b = 0;

				n = z * L * M + x;
				m = z * L * M - x;
				nn = 0;
				mm = x;

				for (iter_type zz = 0; zz < N - z; ++zz)
				{
					for (iter_type yy = 0; yy < M; ++yy)
					{
						for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n];
							sum_b += data_y[mm] * data_y[mm + m];
						}
						for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L];
							sum_b += data_y[mm - L] * data_y[mm + m];
						}
					}
				}
				for (iter_type zz = N - z; zz < N; ++zz)
				{
					for (iter_type yy = 0; yy < M; ++yy)
					{
						for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - len];
							sum_b += data_y[mm] * data_y[mm + m - len];
						}
						for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L - len];
							sum_b += data_y[mm - L] * data_y[mm + m - len];
						}
					}
				}

				put_y[index] = (1.0 / nn) * sum_a;
				put_x[index] = axis_3d_type{ -dx * x, 0, -dz * z };
				++index;
				if (z < NN)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ -dx * x, 0, dz * z };
					++index;
				}
				if (x < LL)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ dx * x, 0, -dz * z };
					++index;
					if (z < NN)
					{
						put_y[index] = (1.0 / nn) * sum_a;
						put_x[index] = axis_3d_type{ dx * x, 0, dz * z };
						++index;
					}
				}
			}

			for (iter_type y = 1; y < MM + 1; ++y)
			{
				/*
				 * count all the interactions on the y-axis
				 */
				Y sum_a = 0, sum_b = 0;

				n = z * L * M + y * L;
				m = z * L * M - y * L;
				nn = 0;
				mm = y * L;

				for (iter_type zz = 0; zz < N - z; ++zz)
				{
					for (iter_type yy = 0; yy < M - y; ++yy)
					{
						for (iter_type xx = 0; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n];
							sum_b += data_y[mm] * data_y[mm + m];
						}
					}
					for (iter_type yy = M - y; yy < M; ++yy)
					{
						for (iter_type xx = 0; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L * M];
							sum_b += data_y[mm - L * M] * data_y[mm + m];
						}
					}
				}
				for (iter_type zz = N - z; zz < N; ++zz)
				{
					for (iter_type yy = 0; yy < M - y; ++yy)
					{
						for (iter_type xx = 0; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - len];
							sum_b += data_y[mm] * data_y[mm + m - len];
						}
					}
					for (iter_type yy = M - y; yy < M; ++yy)
					{
						for (iter_type xx = 0; xx < L; ++xx, ++nn, ++mm)
						{
							sum_a += data_y[nn] * data_y[nn + n - L * M - len];
							sum_b += data_y[mm - L * M] * data_y[mm + m - len];
						}
					}
				}

				put_y[index] = (1.0 / nn) * sum_a;
				put_x[index] = axis_3d_type{ 0, -dy * y, -dz * z };
				++index;

				if (z < NN)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ 0, -dy * y, dz * z };
					++index;
				}

				if (y < MM)
				{
					put_y[index] = (1.0 / nn) * sum_b;
					put_x[index] = axis_3d_type{ 0, dy * y, -dz * z };
					++index;
					if (z < NN)
					{
						put_y[index] = (1.0 / nn) * sum_a;
						put_x[index] = axis_3d_type{ 0, dy * y, dz * z };
						++index;
					}
				}

				for (iter_type x = 1; x < LL + 1; ++x)
				{


					iter_type ss, s, tt, t;
					Y sum_c = 0, sum_d = 0;
					sum_a = 0, sum_b = 0;


					n = z * L * M + y * L + x;
					m = z * L * M - y * L + x;
					s = z * L * M + y * L - x;
					t = z * L * M - y * L - x;

					nn = 0;
					mm = y * L;
					ss = x;
					tt = y * L + x;

					for (iter_type zz = 0; zz < N - z; ++zz)
					{
						for (iter_type yy = 0; yy < M - y; ++yy)
						{
							for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n];
								sum_b += data_y[mm] * data_y[mm + m];

								sum_c += data_y[ss] * data_y[ss + s];
								sum_d += data_y[tt] * data_y[tt + t];
							}
							for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L];
								sum_b += data_y[mm] * data_y[mm + m - L];

								sum_c += data_y[ss - L] * data_y[ss + s];
								sum_d += data_y[tt - L] * data_y[tt + t];
							}
						}
						for (iter_type yy = M - y; yy < M; ++yy)
						{
							for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L * M];
								sum_b += data_y[mm - L * M] * data_y[mm + m];

								sum_c += data_y[ss] * data_y[ss + s - L * M];
								sum_d += data_y[tt - L * M] * data_y[tt + t];
							}
							for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L * M - L];
								sum_b += data_y[mm - L * M] * data_y[mm + m - L];

								sum_c += data_y[ss - L] * data_y[ss + s - L * M];
								sum_d += data_y[tt - L * M - L] * data_y[tt + t];
							}
						}
					}
					for (iter_type zz = N - z; zz < N; ++zz)
					{
						for (iter_type yy = 0; yy < M - y; ++yy)
						{
							for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - len];
								sum_b += data_y[mm] * data_y[mm + m - len];

								sum_c += data_y[ss] * data_y[ss + s - len];
								sum_d += data_y[tt] * data_y[tt + t - len];
							}
							for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L - len];
								sum_b += data_y[mm] * data_y[mm + m - L - len];

								sum_c += data_y[ss - L] * data_y[ss + s - len];
								sum_d += data_y[tt - L] * data_y[tt + t - len];
							}
						}
						for (iter_type yy = M - y; yy < M; ++yy)
						{
							for (iter_type xx = 0; xx < L - x; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L * M - len];
								sum_b += data_y[mm - L * M] * data_y[mm + m - len];

								sum_c += data_y[ss] * data_y[ss + s - L * M - len];
								sum_d += data_y[tt - L * M] * data_y[tt + t - len];
							}
							for (iter_type xx = L - x; xx < L; ++xx, ++nn, ++mm, ++ss, ++tt)
							{
								sum_a += data_y[nn] * data_y[nn + n - L - L * M - len];
								sum_b += data_y[mm - L * M] * data_y[mm + m - L - len];

								sum_c += data_y[ss - L] * data_y[ss + s - L * M - len];
								sum_d += data_y[tt - L * M - L] * data_y[tt + t - len];
							}
						}
					}
					
					put_y[index] = (1.0 / nn) * sum_a;
					put_x[index] = axis_3d_type{ -dx * x, -dy * y, -dz * z };
					++index;
					if (z < NN)
					{
						put_y[index] = (1.0 / nn) * sum_d;
						put_x[index] = axis_3d_type{ -dx * x, -dy * y, dz * z };
						++index;
					}
					if (y < MM)
					{
						put_y[index] = (1.0 / nn) * sum_b;
						put_x[index] = axis_3d_type{ -dx * x, dy * y, -dz * z };
						++index;
						if (z < NN)
						{
							put_y[index] = (1.0 / nn) * sum_c;
							put_x[index] = axis_3d_type{ -dx * x, dy * y, dz * z };
							++index;
						}
					}
					if (x < LL)
					{
						put_y[index] = (1.0 / nn) * sum_c;
						put_x[index] = axis_3d_type{ dx * x, -dy * y, -dz * z };
						++index;
						if (z < NN)
						{
							put_y[index] = (1.0 / nn) * sum_b;
							put_x[index] = axis_3d_type{ dx * x, -dy * y, dz * z };
							++index;
						}
						if (y < MM)
						{
							put_y[index] = (1.0 / nn) * sum_d;
							put_x[index] = axis_3d_type{ dx * x, dy * y, -dz * z };
							++index;
							if (z < NN)
							{
								put_y[index] = (1.0 / nn) * sum_a;
								put_x[index] = axis_3d_type{ dx * x, dy * y, dz * z };
								++index;
							}
						}
					}
				}
			}
		}
	};


	len_type block_count = std::min(NN / 2, OMP_BLOCK_COUNT);
#	pragma omp parallel for
	for (iter_type n = 0; n < block_count; ++n)
	{
		if (n == 0)
		{
			corr_m();
		}
		corr(n, block_count);
	}
	

	symphas::lib::sort_data(put_main_x, put_main_y, L, M, N);
	return vector_data<Y, 3>{ put_main_x, put_main_y, len };
}







