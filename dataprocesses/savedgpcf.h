#pragma once

#include "savedefines.h"

#define DP_VEC_NAME "dpv"
#define DP_ABS_NAME "dpa"
#define DP_PT_NAME "dpp"

struct SaveDGPCF {};

DLLPROC extern char const* dpset[];
DLLPROC extern char const* adpset[];
DLLPROC extern char const* pdpset[];
DLLPROC extern const char* fit_title_str_dp;

/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_pts<SaveDGPCF>(std::vector<std::tuple<iter_type, axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *en = symphas::io::open_data_file(dir, DP_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(en, "# domain growth curve over the time evolution");

	for (auto const &[i, _, y] : data)
	{
		fprintf(en, "%.6d %.6f\n", i, y);
	}
	fclose(en);

	double ranges[4];
	symphas::lib::fill_ranges(data, ranges);

	// set the range based on index instead
	ranges[0] = std::get<0>(data.front());
	ranges[1] = std::get<0>(data.back());

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { rangestr, dpset[0], dpset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, DP_ABS_NAME, index, id, rangestr);
}



template<>
inline void symphas::io::save_pts<SaveDGPCF>(std::vector<std::tuple<iter_type, axis_1d_type, complex_t>> const&, const char*, int, size_t)
{
}


template<>
inline void symphas::io::save_abs<SaveDGPCF>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE* def = symphas::io::open_data_file(dir, DP_PT_NAME, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6f %.6f\n", x, y);
	}
	fclose(def);


	axis_coord_t ranges[2];
	scalar_t extrema[2];
	auto data_1 = std::vector(data.begin() + ((data.size() > 0) ? 1 : 0), data.end());
	symphas::lib::fill_sorted_ranges(data_1, ranges, extrema);
	scalar_t dy = std::max((extrema[1] - extrema[0]) * 0.2, extrema[0] * 0.1);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { pdpset[0], pdpset[1] };

	char name[BUFFER_LENGTH];
	copy_data_file_name(DATA_DIR_RELATIVE_PLOT, DP_PT_NAME, index, id, DataFileType::POSTPROC_DATA, name);
	char addstr[BUFFER_LENGTH_L2];
	snprintf(addstr, BUFFER_LENGTH_L2, pdpset[2],
		rangestr,		// ranges
		name,			// filename
		fit_title_str_dp);

	symphas::io::write_postproc_plot_file(sets, 2, dir, DP_PT_NAME, index, id, addstr);
}


template<>
inline void symphas::io::save_abs<SaveDGPCF>(std::vector<std::pair<axis_1d_type, scalar_t[2]>> const& data, const char* dir, int index, size_t id)
{
	FILE* def = symphas::io::open_data_file(dir, DP_PT_NAME, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		auto x0 = (x == 0) ? x + 1 : x;
		auto y01 = (y[0] == 0) ? y[0] + 1e-3 : y[0];
		auto y02 = (y[1] == 0) ? y[1] + 1e-3 : y[1];
		fprintf(def, "%.6f %.6f %.6f\n", x0, y01, y02);
	}
	fclose(def);

	const char* sets[] = { pdpset[0], pdpset[1] };

	char name[BUFFER_LENGTH];
	copy_data_file_name(DATA_DIR_RELATIVE_PLOT, DP_PT_NAME, index, id, DataFileType::POSTPROC_DATA, name);
	char addstr[BUFFER_LENGTH_L2];


	auto fitpt = (data.back().first - data.front().first) / 4 + data.front().first;
	snprintf(addstr, BUFFER_LENGTH_L2, pdpset[2],
		symphas::io::gp::ranges_auto,		// ranges
		name,								// filename
		fitpt,								// from where fitting starts
		fit_title_str_dp);

	symphas::io::write_postproc_plot_file(sets, 2, dir, DP_PT_NAME, index, id, addstr);
}

template<>
inline void symphas::io::save_abs<SaveDGPCF>(std::vector<std::pair<axis_1d_type, std::vector<std::pair<scalar_t, scalar_t>>>> const& data, const char* dir, int index, size_t id)
{
	FILE* def = symphas::io::open_data_file(dir, DP_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);

	axis_coord_t ranges[2]{ +INFINITY, -INFINITY };
	scalar_t extrema[2]{ +INFINITY, -INFINITY };

	for (auto const& [x, y] : data)
	{
		axis_coord_t ranges0[2];
		scalar_t extrema0[2];
		symphas::lib::fill_sorted_ranges(y, ranges0, extrema0);
		
		if (ranges0[0] < ranges[0])
		{
			ranges[0] = ranges0[0];
		}
		if (ranges0[1] > ranges[1])
		{
			ranges[1] = ranges0[1];
		}
		if (extrema0[0] < extrema[0])
		{
			extrema[0] = extrema0[0];
		}
		if (extrema0[1] > extrema[1])
		{
			extrema[1] = extrema0[1];
		}

		for (auto const& [x0, y0] : y)
		{
			fprintf(def, "%.6f %.6f\n", x0, y0);
		}
		fprintf(def, "\n\n");
	}
	fclose(def);

	scalar_t dy = std::max((extrema[1] - extrema[0]) * 0.2, extrema[0] * 0.1);
	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { adpset[0], adpset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, DP_ABS_NAME, index, id, rangestr);
}

