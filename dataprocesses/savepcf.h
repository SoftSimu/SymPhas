#pragma once

#include "savedefines.h"

#define PCF_VEC_NAME "pv"
#define PCF_ABS_NAME "pa"

struct SavePCF {};

DLLPROC extern char const* apcfset[];
DLLPROC extern char const* pcfset[];

/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_abs<SavePCF>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char *dir, int index, size_t id)
{
	FILE *pcf = symphas::io::open_data_file(dir, PCF_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(pcf, "# spherical average of 2-point correlation data in frame %d\n", index);

	for (auto const &[x, y] : data)
	{
		fprintf(pcf, "%.6f %.6f\n", x, y);
	}
	fprintf(pcf, "\n\n");
	fclose(pcf);
	
	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);
	scalar_t dy = std::max(1e-2, std::max((extrema[1] - extrema[0]) * 0.2, std::abs(extrema[0] * 0.1)));

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { apcfset[0], apcfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, PCF_ABS_NAME, index, id, rangestr);
}


/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_abs<SavePCF>(std::vector<std::pair<axis_1d_type, complex_t>> const& data, const char* dir, int index, size_t id)
{
	FILE* pcf = symphas::io::open_data_file(dir, PCF_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(pcf, "# spherical average of 2-point correlation data in frame %d\n", index);

	for (auto const& [x, y] : data)
	{
		fprintf(pcf, "%.6f %.6lf\n", x, std::abs(y));
	}
	fprintf(pcf, "\n\n");
	fclose(pcf);

	const char* sets[] = { apcfset[0], apcfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, PCF_ABS_NAME, index, id);
}


/*
 * writes the point correlation function data to a file
 */
template<>
inline void symphas::io::save_vec<SavePCF>(std::vector<std::pair<axis_3d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *pcf = symphas::io::open_data_file(dir, PCF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(pcf, "# spherical average of 2-point correlation data in frame %d\n", index);

	auto max = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second < b.second;	});
	auto min = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second > b.second; });
	fprintf(pcf, "# limits:%lf,%lf\n", max->second, min->second);

	auto[L, M, N] = symphas::lib::get_sorted_dimensions(data)._3();
	for (iter_type k = 0; k < N; k++)
	{
		fprintf(pcf, "%d ", k);
		for (iter_type i = 0; i < L; i++)
		{
			axis_coord_t x = data[i + k * L * M].first[0];
			fprintf(pcf, "%.6f ", x);
		}
		fprintf(pcf, "\n");

		for (iter_type j = 0; j < M; j++)
		{
			iter_type jj = j * L + k * L * M;
			axis_coord_t y = data[jj].first[1];
			fprintf(pcf, "%.6f ", y);
			for (iter_type i = 0; i < L; i++)
			{
				auto v = data[jj + i].second;
				fprintf(pcf, "% 10E ", v);
			}
			fprintf(pcf, "\n");
		}
		fprintf(pcf, "\n");
	}
	fprintf(pcf, "\n");
	fclose(pcf);

	double ranges[] = {
		data[0].first[0], data[data.size() - 1].first[0],
		data[0].first[1], data[data.size() - 1].first[1] };


	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { pcfset[0], pcfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, PCF_VEC_NAME, index, id, rangestr);
}

/*
 * writes the point correlation function data to a file
 */
template<>
inline void symphas::io::save_vec<SavePCF>(std::vector<std::pair<axis_2d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *pcf = symphas::io::open_data_file(dir, PCF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(pcf, "# 2-point correlation data in frame %d\n", index);

	auto[L, M] = symphas::lib::get_sorted_dimensions(data)._2();

	fprintf(pcf, "%d ", 0);
	for (iter_type i = 0; i < L; i++)
	{
		axis_coord_t x = data[i].first[0];
		fprintf(pcf, "%.6f ", x);
	}
	fprintf(pcf, "\n");

	for (iter_type j = 0; j < M; j++)
	{
		iter_type jj = j * L;
		axis_coord_t y = data[jj].first[1];
		fprintf(pcf, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", y);
		for (iter_type i = 0; i < L; i++)
		{
			auto v = data[jj + i].second;
			fprintf(pcf, "% 10E ", v);
		}
		fprintf(pcf, "\n");
	}
	fprintf(pcf, "\n");
	fclose(pcf);

	axis_coord_t ranges[4];
	symphas::lib::fill_sorted_ranges(data, ranges);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { pcfset[0], pcfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, PCF_VEC_NAME, index, id);//, rangestr);
}


/*
 * writes the point correlation function data to a file
 */
template<>
inline void symphas::io::save_vec<SavePCF>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE* pcf = symphas::io::open_data_file(dir, PCF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(pcf, "# 2-point correlation data in frame %d\n", index);

	auto max = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second < b.second;	});
	auto min = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second > b.second; });
	fprintf(pcf, "# limits:%lf,%lf\n", max->second, min->second);

	auto L = symphas::lib::get_sorted_dimensions(data);

	for (iter_type i = 0; i < L; i++)
	{
		fprintf(pcf, "% 10.6f % 10E ", data[i].first, data[i].second);
	}
	fprintf(pcf, "\n\n");
	fclose(pcf);

	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0], extrema[1]);
	const char* sets[] = { pcfset[0], pcfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, PCF_VEC_NAME, index, id, rangestr);
}

