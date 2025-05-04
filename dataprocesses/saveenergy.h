#pragma once

#include "savedefines.h"

#define EN_VEC_NAME "ev"
#define EN_ABS_NAME "ea"

struct SaveEnergy {};

DLLPROC extern char const* enset[];

/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_pts<SaveEnergy>(std::vector<std::tuple<iter_type, axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *en = symphas::io::open_data_file(dir, EN_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
	fprintf(en, "# energy curve over the time evolution");

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
	const char* sets[] = { rangestr, enset[0], enset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, EN_ABS_NAME, index, id, rangestr);
}



template<>
inline void symphas::io::save_pts<SaveEnergy>(std::vector<std::tuple<iter_type, axis_1d_type, complex_t>> const&, const char*, int, size_t)
{
}



