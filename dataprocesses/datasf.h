#pragma once

#include <execution>

#include "datapcf.h"
#include "datadgpcf.h"

DEFINE_ALGORITHM_SCALAR(SF)
DEFINE_ALGORITHM_VECTOR(SF)
DEFINE_DATA(SF, SCALAR, ALG_SCALAR, ALG_VECTOR)

template<typename T>
auto get_p_hat(T* data_y,
	double dx,
	double qx,
	len_type L)
{
	std::complex<scalar_t> sum;

	double
		ddqi = qx * dx;

	iter_type n = 0;
	double ddi = 0;
	for (iter_type i = 0; i < L; ++i, ++n)
	{
		using std::sin;
		using std::cos;

		auto p = ddi;
		auto c = data_y[n] * cos(p);
		auto s = data_y[n] * sin(p);
		sum += std::complex{ c, s };
		ddi += ddqi;
	}

	return std::norm(sum);
}


template<typename T>
auto get_p_hat(T* data_y, 
	double dx, double dy, 
	double qx, double qy, 
	len_type L, len_type M)
{
	std::complex<scalar_t> sum;


	double
		ddqi = qx * dx,
		ddqj = qy * dy;

	double ddj = 0;
	iter_type n = 0;
	for (iter_type j = 0; j < M; ++j)
	{
		double ddi = 0;
		for (iter_type i = 0; i < L; ++i, ++n)
		{
			using std::sin;
			using std::cos;

			if (data_y[n] > 0)
			{
				// position of particle = (i*dx, j*dy)
				auto p = ddi + ddj;
				sum += std::complex{ cos(p), -sin(p) };
			}
			ddi += ddqi;
		}
		ddj += ddqj;
	}

	return std::norm(sum);
}

template<typename T>
auto get_p_hat(T* data_y,
	double dx, double dy, double dz,
	double qx, double qy, double qz,
	len_type L, len_type M, len_type N)
{
	std::complex<scalar_t> sum;

	double
		ddqi = qx * dx,
		ddqj = qy * dy,
		ddqk = qz * dz;

	double ddk = 0;
	iter_type n = 0;
	for (iter_type k = 0; k < N; ++k)
	{
		double ddj = 0;
		for (iter_type j = 0; j < M; ++j)
		{
			double ddi = 0;
			for (iter_type i = 0; i < L; ++i, ++n)
			{
				using std::sin;
				using std::cos;

				if (data_y[n] > 0)
				{
					auto p = ddi + ddj + ddk;
					sum += std::complex{ cos(p), -sin(p) };
					ddi += ddqi;
				}
			}
			ddj += ddqj;
		}
		ddk += ddqk;
	}

	return std::norm(sum);
}


ALGORITHM_VECTOR_1D_DECLARATION(SF)
{

	axis_1d_type* data_sfx = new axis_1d_type[len];
	symphas::dft::fill_x_axis(data_sfx, data_x, len);
	Y* data_sfy = new Y[len];

	double const
		dx = ((*(data_x + 1)) - (*data_x));

	double r = 1.0 / len;
	iter_type n = 0;
	for (auto* it = data_sfx; it < data_sfx + len; ++it)
	{
		auto& qx = *it;
		data_sfy[n++] = r * get_p_hat(data_y, dx, qx, L);
	}

	return symphas::FieldAxis{ data_sfx, data_sfy, len };
}

ALGORITHM_VECTOR_2D_DECLARATION(SF)
{
	axis_2d_type* data_sfx = new axis_2d_type[len];
	symphas::dft::fill_x_axis(data_sfx, data_x, len);
	Y* data_sfy = new Y[len];

	double r = 1.0 / len;

#	pragma omp parallel for
	for (iter_type i = 0; i < len; ++i)
	{
		data_sfy[i] = (data_y[i] > 0) ? 1 : 0;
	}

	// This implementation is verified to work, including correct placement.
	complex_t* data_sft = new complex_t[len];
	symphas::dft::dft(data_sfy, data_sft, L, M);

#	pragma omp parallel for
	for (iter_type j = 0; j < M; ++j)
	{
		iter_type b = (j <= M / 2) ? (M / 2 - j) : (M - (j - M / 2));
		b *= L;
		for (iter_type i = 0; i < L; ++i)
		{
			iter_type a = (i <= L / 2) ? ((L / 2 - i)) : ((L - (i - L / 2)));

			data_sfy[a + b] = r * std::norm(data_sft[j * L + i]);
		}
	}
	delete[] data_sft;


	// This implementation is verified to work, including correct placement.
	//complex_t* data_sft = new complex_t[len];
	//symphas::dft::dft(data_y, data_sft, L, M);
#	//pragma omp parallel for
	//for (iter_type i = 0; i < L; ++i)
	//{
	//	iter_type a = (i <= L / 2) ? ((L / 2 - i)) : ((L - (i - L / 2)));
	//	a *= L;
	//	for (iter_type j = 0; j < M; ++j)
	//	{
	//		iter_type b = (j <= M / 2) ? (M / 2 - j) : (M - (j - M / 2));
	//		data_sfy[a + b] = r * std::norm(data_sft[j * M + i]);
	//	}
	//}
	//delete[] data_sft;

	// This implementation is verified to work.
	//std::transform(std::execution::par_unseq, data_sfx, data_sfx + len, data_sfy, [&](auto& pt)
	//	{
	//		auto& [qx, qy] = pt;
	//		return r * get_p_hat(data_y, dx, dy, qx, qy, L, M);
	//	});


//	// This implementation is verified to work. OFFICIALLY.
//#	pragma omp parallel for
//	for (iter_type n = 0; n < len; ++n)
//	{
//		auto& [qx, qy] = data_sfx[n];
//		data_sfy[n] = r * get_p_hat(data_y, dx, dy, qx, qy, L, M);
//	}

	return symphas::FieldAxis{ data_sfx, data_sfy, len };

}

ALGORITHM_VECTOR_3D_DECLARATION(SF)
{
	axis_3d_type* data_sfx = new axis_3d_type[len];
	symphas::dft::fill_x_axis(data_sfx, data_x, len);
	Y* data_sfy = new Y[len];

	double r = 1.0 / len;
	

#	pragma omp parallel for
	for (iter_type i = 0; i < len; ++i)
	{
		data_sfy[i] = (data_y[i] > 0) ? 1 : 0;
	}

	// This implementation is verified to work, including correct placement.
	complex_t* data_sft = new complex_t[len];
	symphas::dft::dft(data_sfy, data_sft, L, M, N);


#	pragma omp parallel for
	for (iter_type k = 0; k < N; ++k)
	{
		iter_type c = (k <= N / 2) ? (N / 2 - k) : (N - (k - N / 2));
		c *= L * M;
		for (iter_type j = 0; j < M; ++j)
		{
			iter_type b = (j <= M / 2) ? (M / 2 - j) : (M - (j - M / 2));
			b *= L;
			for (iter_type i = 0; i < L; ++i)
			{
				iter_type a = (i <= L / 2) ? ((L / 2 - i)) : ((L - (i - L / 2)));
				data_sfy[a + b + c] = r * std::norm(data_sft[k * L * M + j * L + i]);
			}
		}
	}
	delete[] data_sft;

//
//	// This implementation is verified to work. OFFICIALLY.
//#	pragma omp parallel for
//	for (iter_type n = 0; n < len; ++n)
//	{
//		auto& [qx, qy, qz] = data_sfx[n];
//		data_sfy[n] = r * get_p_hat(data_y, dx, dy, dz, qx, qy, qz, L, M);
//	}

	//for (iter_type n = 0; n < len; ++n)
	//{
	//	auto& [qx, qy, qz] = data_sfx[n];
	//	data_sfy[n] = r * get_p_hat(data_y, dx, dy, dz, qx, qy, qz, L, M, N);
	//}

	return symphas::FieldAxis{ data_sfx, data_sfy, len };

}

ALGORITHM_SCALAR_1D_DECLARATION(SF)
{
	auto sfv_data = RUN_ALGORITHM_VECTOR_1D(SF)(ALGORITHM_ARGS_1D);
	return symphas::lib::radial_avg(sfv_data);
}

ALGORITHM_SCALAR_2D_DECLARATION(SF)
{
	auto sfv_data = RUN_ALGORITHM_VECTOR_2D(SF)(ALGORITHM_ARGS_2D);
	return symphas::lib::radial_avg(sfv_data);
}

ALGORITHM_SCALAR_3D_DECLARATION(SF)
{
	auto sfv_data = RUN_ALGORITHM_VECTOR_3D(SF)(ALGORITHM_ARGS_3D);
	return symphas::lib::radial_avg(sfv_data);
}



