#pragma once

#include "savedefines.h"
#include "datalib.h"

#define SF_VEC_NAME "sv"
#define SF_ABS_NAME "sa"

struct SaveSF {};

/*
 * writes the structure factor data to a file
 */

DLLPROC extern char const* sfset[];
DLLPROC extern char const* asfset[];

template<>
inline void symphas::io::save_vec<SaveSF>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE* sf = symphas::io::open_data_file(dir, SF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(sf, "# structure data in frame %d\n", index);

	for (auto const& [q, v] : data)
	{
		double z = std::abs(v);

		fprintf(sf, "%.6f %.6f\n", q, z);
	}
	fprintf(sf, "\n\n");
	fclose(sf);

	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0], extrema[1]);
	const char* sets[] = { sfset[0], sfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, SF_VEC_NAME, index, id, rangestr);
}



template<>
inline void symphas::io::save_vec<SaveSF>(std::vector<std::pair<axis_2d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *sf = symphas::io::open_data_file(dir, SF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(sf, "# structure data in frame %d\n", index);

	auto [L, M] = symphas::lib::get_sorted_dimensions(data)._2();

	fprintf(sf, "%d ", 0);
	for (iter_type i = 0; i < L; i++)
	{
		axis_coord_t x = data[i].first[0];
		fprintf(sf, "%.6f ", x);
	}
	fprintf(sf, "\n");

	for (iter_type j = 0; j < M; j++)
	{
		iter_type jj = j * L;
		axis_coord_t y = data[jj].first[1];
		fprintf(sf, "%.6f ", y);
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(sf, "% 10E ", data[jj + i].second);
		}
		fprintf(sf, "\n");
	}

	fprintf(sf, "\n");
	fclose(sf);

	axis_coord_t ranges[4];
	symphas::lib::fill_sorted_ranges(data, ranges);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { sfset[0], sfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, SF_VEC_NAME, index, id, rangestr);
}

/*
 * writes the structure factor data to a file
 */
template<>
inline void symphas::io::save_vec<SaveSF>(std::vector<std::pair<axis_3d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *sf = symphas::io::open_data_file(dir, SF_VEC_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(sf, "# structure data in frame %d\n", index);
	

	auto [L, M, N] = symphas::lib::get_sorted_dimensions(data)._3();
	for (iter_type k = 0; k < N; k++)
	{
		fprintf(sf, "%d ", k);
		for (iter_type i = 0; i < L; i++)
		{
			iter_type ii = i + k * L * M;
			axis_coord_t x = data[ii].first[0];
			fprintf(sf, "%.6f ", x);
		}
		fprintf(sf, "\n");

		for (iter_type j = 0; j < M; j++)
		{
			iter_type jj = j * L + k * L * M;
			axis_coord_t y = data[jj].first[1];
			fprintf(sf, "%.6f ", y);
			for (iter_type i = 0; i < L; i++)
			{
				fprintf(sf, "% 10E ", data[jj * i].second);
			}
			fprintf(sf, "\n");
		}
		fprintf(sf, "\n");
	}
	fprintf(sf, "\n");
	fclose(sf);
	/*
	for (auto const &[q, v] : data)
	{
		double
			qx = q[0],
			qy = q[1],
			qz = q[2];

		fprintf(sf, "%.6f %.6f %.6f %.6f\n", qx, qy, qz, v);
	}
	fprintf(sf, "\n\n");
	fclose(sf);*/

	axis_coord_t ranges[6];
	symphas::lib::fill_sorted_ranges(data, ranges);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { sfset[0], sfset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, SF_VEC_NAME, index, id, rangestr);
}

/*
 * writes the absolute structure factor data to a file
 */
template<>
inline void symphas::io::save_abs<SaveSF>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *sf = symphas::io::open_data_file(dir, SF_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(sf, "# spherical average of structure data in frame %d\n", index);
	
	//auto trimmed_data = symphas::lib::trim_right(data, 0);
	//auto balanced_data = symphas::lib::balance_axis_as(trimmed_data);
	//auto averaged_data = symphas::lib::smooth_data(data, 0.4);



	fprintf(sf, "\n\n");
	fclose(sf);

	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);
	extrema[0] = std::max(std::abs(extrema[0]), data[1].first / 2);


	for (auto const& [x, v] : data)
	{
		if (x <= extrema[0])
		{
			fprintf(sf, "%.6f %.6E\n", extrema[0], v);
		}
		fprintf(sf, "%.6f %.6E\n", x, v);
	}


	char name[BUFFER_LENGTH];
	copy_data_file_name(DATA_DIR_RELATIVE_PLOT, SF_ABS_NAME, index, id, DataFileType::POSTPROC_DATA, name);


	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, data[1].first, ranges[1], extrema[0], extrema[1] * 1.15);
	const char* sets[] = { asfset[0], asfset[1] };

	char addstr[BUFFER_LENGTH_L2];
	snprintf(addstr, BUFFER_LENGTH_L2, asfset[2],
		name,			// filename
		2,				// dimensions
		rangestr);		// ranges

	symphas::io::write_postproc_plot_file(sets, 2, dir, SF_ABS_NAME, index, id, addstr);
}
