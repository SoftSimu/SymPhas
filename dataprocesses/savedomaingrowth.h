#pragma once

#include "savedefines.h"

#define DG_VEC_NAME "dgv"
#define DG_ABS_NAME "dga"

struct SaveDomainGrowth {};

DLLPROC extern char const* dgset[];
DLLPROC extern char const* adgset[];
DLLPROC extern const char* fit_title_str_dg;

/*
 * writes the absolute point correlation function data to a file
 */
template<>
inline void symphas::io::save_pts<SaveDomainGrowth>(std::vector<std::tuple<iter_type, axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	FILE *en = symphas::io::open_data_file(dir, DG_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);
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
	const char* sets[] = { rangestr, dgset[0], dgset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, DG_ABS_NAME, index, id, rangestr);
}



template<>
inline void symphas::io::save_pts<SaveDomainGrowth>(std::vector<std::tuple<iter_type, axis_1d_type, complex_t>> const&, const char*, int, size_t)
{
}

template<>
inline void symphas::io::save_abs<SaveDomainGrowth>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{

	FILE* def = symphas::io::open_data_file(dir, DG_ABS_NAME, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6f %.6f\n", x, y);
	}
	fclose(def);


	axis_coord_t ranges[2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);
	scalar_t dy = std::max((extrema[1] - extrema[0]) * 0.2, extrema[0] * 0.1);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { adgset[0], adgset[1] };

	char name[BUFFER_LENGTH];
	copy_data_file_name(DATA_DIR_RELATIVE_PLOT, DG_ABS_NAME, index, id, DataFileType::POSTPROC_DATA, name);
	char addstr[BUFFER_LENGTH_L2];
	snprintf(addstr, BUFFER_LENGTH_L2, adgset[2],
		rangestr,		// ranges
		name,			// filename
		fit_title_str_dg);

	symphas::io::write_postproc_plot_file(sets, 2, dir, DG_ABS_NAME, index, id, addstr);
}




