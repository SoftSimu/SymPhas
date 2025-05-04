#pragma once

#include "savedefines.h"

#define DE_VEC_NAME "dv"
#define DE_ABS_NAME "da"

struct SaveDensity {};

DLLPROC extern char const* deset[];
DLLPROC extern char const *adeset[];

/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_pts<SaveDensity>(std::vector<std::tuple<iter_type, axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *de = symphas::io::open_data_file(dir, DE_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(de, "# energy curve over the time evolution\n");

	for (auto const &[i, _, y] : data)
	{
		fprintf(de, "%.6d %.6f\n", i, y);
	}
	fclose(de);

	double ranges[4];
	symphas::lib::fill_ranges(data, ranges);

	// set the range based on index instead
	ranges[0] = std::get<0>(data.front());
	ranges[1] = std::get<0>(data.back());

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { deset[0], deset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, DE_ABS_NAME, index, id, rangestr);
}


template<>
inline void symphas::io::save_pts<SaveDensity>(std::vector<std::tuple<iter_type, axis_1d_type, complex_t>> const&, const char*, int, size_t) {}



template<>
inline void symphas::io::save_abs<SaveDensity>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE* def = symphas::io::open_data_file(dir, DE_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6f %.6f\n", x, y);
	}
	fclose(def);


	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);
	scalar_t dy = std::max(1e-2, std::max((extrema[1] - extrema[0]) * 0.2, std::abs(extrema[0] * 0.1)));

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { adeset[0], adeset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, DE_ABS_NAME, index, id, rangestr);
}


