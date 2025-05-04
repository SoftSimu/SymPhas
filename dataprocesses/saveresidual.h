#pragma once

#include "savedefines.h"
#include "datalib.h"

#define RS_VEC_NAME "rs"

struct SaveResidual {};

/*
 * writes the structure factor data to a file
 */

//DLLPROC extern char const* sfset[];
//DLLPROC extern char const* asfset[];

template<>
inline void symphas::io::save_vec<SaveResidual>(std::vector<std::pair<axis_1d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	scalar_t* data_y = new scalar_t[data.size()];
	for (iter_type i = 0; i < data.size(); ++i)
	{
		data_y[i] = data[i].second;
	}
	auto dim = symphas::lib::get_dimensions<1, scalar_t>(data);
	len_type dims[]{ dim };

	char buffer[1024];
	sprintf(buffer, "%s/%s", dir, RS_VEC_NAME);
	symphas::io::save_grid_plotting(data_y, symphas::io::write_info(buffer, index, id, DataFileType::POSTPROC_DATA), symphas::grid_info(dims, 1));
}



template<>
inline void symphas::io::save_vec<SaveResidual>(std::vector<std::pair<axis_2d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	scalar_t* data_y = new scalar_t[data.size()];
	for (iter_type i = 0; i < data.size(); ++i)
	{
		data_y[i] = data[i].second;
	}
	auto [L, M] = symphas::lib::get_dimensions<2, scalar_t>(data)._2();
	len_type dims[]{ L, M };

	char buffer[1024];
	sprintf(buffer, "%s/%s", dir, RS_VEC_NAME);
	symphas::io::save_grid_plotting(data_y, symphas::io::write_info(buffer, index, index, DataFileType::POSTPROC_DATA), symphas::grid_info(dims, 2));
}

/*
 * writes the structure factor data to a file
 */
template<>
inline void symphas::io::save_vec<SaveResidual>(std::vector<std::pair<axis_3d_type, scalar_t>> const& data, const char* dir, int index, size_t id)
{
	scalar_t* data_y = new scalar_t[data.size()];
	for (iter_type i = 0; i < data.size(); ++i)
	{
		data_y[i] = data[i].second;
	}
	auto [L, M, N] = symphas::lib::get_dimensions<3, scalar_t>(data)._3();
	len_type dims[]{ L, M, N };

	char buffer[1024];
	sprintf(buffer, "%s/%s", dir, RS_VEC_NAME);
	symphas::io::save_grid_plotting(data_y, symphas::io::write_info(buffer, index, id, DataFileType::POSTPROC_DATA), symphas::grid_info(dims, 3));
}


