
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


#include "writegp.h"

namespace symphas::io::gp
{
	template<typename helper_specialized>
	struct gp_plotting_helper_impl : gp_plotting_helper
	{
		virtual size_t get_dim() const override
		{
			return cast_const().get_dim();
		}

		virtual len_type get_len(iter_type n) const override
		{
			return cast_const().get_len(n);
		}

		virtual len_type get_index(const len_type(&coords)[3]) const override
		{
			return cast_const().get_index(&coords[0]);
		}

		virtual len_type get_index(const len_type(&coords)[2]) const override
		{
			return cast_const().get_index(&coords[0]);
		}

		virtual len_type get_index(const len_type(&coords)[1]) const override
		{
			return cast_const().get_index(&coords[0]);
		}

		virtual len_type get_index(std::initializer_list<len_type> const& list) const override
		{
			return cast_const().get_index(list.begin());
		}

		template<size_t D>
		len_type get_offset(const len_type(&coords)[D]) const
		{
			return cast_const().get_offset(&coords[0]);
		}

		virtual double get_position(Axis ax, const len_type(&coords)[1]) const override
		{
			return cast_const().get_position(ax, &coords[0]);
		}

		virtual double get_position(Axis ax, const len_type(&coords)[2]) const override
		{
			return cast_const().get_position(ax, &coords[0]);
		}

		virtual double get_position(Axis ax, const len_type(&coords)[3]) const override
		{
			return cast_const().get_position(ax, &coords[0]);
		}

		virtual double get_position(Axis ax, std::initializer_list<len_type> const& list) const override
		{
			return cast_const().get_position(ax, list.begin());
		}

		template<size_t D>
		double get_position(Axis ax, const len_type(&coords)[D]) const
		{
			return cast_const().get_position(ax, &coords[0]);
		}

		virtual double get_position(Axis ax, len_type coord) const override
		{
			return cast_const().get_position(ax, coord);
		}

		inline helper_specialized const& cast_const() const
		{
			return *static_cast<helper_specialized const*>(this);
		}

		inline helper_specialized& cast()
		{
			return *static_cast<helper_specialized*>(this);
		}

		virtual ~gp_plotting_helper_impl() override {}
	};
}


template<>
struct symphas::io::gp::gp_plotting_helper_specialized<3> : symphas::io::gp::gp_plotting_helper_impl<gp_plotting_helper_specialized<3>>
{
	using parent_type = gp_plotting_helper_impl<gp_plotting_helper_specialized<3>>;
	using parent_type::get_index;
	using parent_type::get_position;

	gp_plotting_helper_specialized(symphas::io::write_info const& winfo, symphas::grid_info const& ginfo) :
		h{}, pos0{}, len{}, stride{}, dims{}, offset{}
	{
		h[2] = ginfo.INTERVAL_Zh;
		h[1] = ginfo.INTERVAL_Yh;
		h[0] = ginfo.INTERVAL_Xh;

		len[2] = ginfo.INTERVAL_Zc;
		len[1] = ginfo.INTERVAL_Yc;
		len[0] = ginfo.INTERVAL_Xc;

		dims[1] = ginfo.DOMAIN_Yc;
		dims[1] = ginfo.DOMAIN_Zc;
		dims[0] = ginfo.DOMAIN_Xc;

		stride[0] = 1;
		stride[1] = dims[0];
		stride[2] = dims[0] * dims[1];

		len_type pos[3]{};
		offset[2] = ginfo.INTERVAL_Z0 - ginfo.DOMAIN_Z0;
		offset[1] = ginfo.INTERVAL_Y0 - ginfo.DOMAIN_Y0;
		offset[0] = ginfo.INTERVAL_X0 - ginfo.DOMAIN_X0;

		pos0[2] = ginfo.INTERVAL_Z0;
		pos0[1] = ginfo.INTERVAL_Y0;
		pos0[0] = ginfo.INTERVAL_X0;
	}

	len_type get_index(const len_type* coords) const
	{
		len_type coords_corrected[3]{};
		for (iter_type i = 0; i < 3; ++i)
		{
			coords_corrected[i] = (coords[i] + offset[i]) % dims[i];
		}
		return grid::index_from_position((len_type[3]){ coords_corrected[0], coords_corrected[1], coords_corrected[2] }, stride);
	}

	double get_position(Axis ax, const len_type* coords) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coords[n] * h[n];
	}

	double get_position(Axis ax, len_type coord) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coord * h[n];
	}

	size_t get_dim() const
	{
		return 3;
	}

	len_type get_len(iter_type n) const
	{
		return len[n];
	}

	double h[3];
	double pos0[3];
	len_type len[3];
	len_type stride[3];
	len_type dims[3];
	len_type offset[3];

};

template<>
struct symphas::io::gp::gp_plotting_helper_specialized<2> : symphas::io::gp::gp_plotting_helper_impl<gp_plotting_helper_specialized<2>>
{
	using parent_type = gp_plotting_helper_impl<gp_plotting_helper_specialized<2>>;
	using parent_type::get_index;
	using parent_type::get_position;

	gp_plotting_helper_specialized(symphas::io::write_info const& winfo, symphas::grid_info const& ginfo) :
		h{}, pos0{}, len{}, stride{}, dims{}, offset{}
	{
		h[1] = ginfo.INTERVAL_Yh;
		h[0] = ginfo.INTERVAL_Xh;

		len[1] = ginfo.INTERVAL_Yc;
		len[0] = ginfo.INTERVAL_Xc;

		dims[1] = ginfo.DOMAIN_Yc;
		dims[0] = ginfo.DOMAIN_Xc;

		stride[0] = 1;
		stride[1] = dims[0];

		offset[1] = ginfo.INTERVAL_Y0 - ginfo.DOMAIN_Y0;
		offset[0] = ginfo.INTERVAL_X0 - ginfo.DOMAIN_X0;

		pos0[1] = ginfo.INTERVAL_Y0;
		pos0[0] = ginfo.INTERVAL_X0;
	}

	len_type get_index(const len_type* coords) const
	{
		len_type coords_corrected[2]{};
		for (iter_type i = 0; i < 2; ++i)
		{
			coords_corrected[i] = (coords[i] + offset[i]) % dims[i];
		}
		return grid::index_from_position((len_type[2]) { coords_corrected[0], coords_corrected[1] }, stride);
	}

	double get_position(Axis ax, const len_type* coords) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coords[n] * h[n];
	}

	double get_position(Axis ax, len_type coord) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coord * h[n];
	}

	size_t get_dim() const
	{
		return 2;
	}

	len_type get_len(iter_type n) const
	{
		return len[n];
	}

	double h[2];
	double pos0[2];
	len_type len[2];
	len_type stride[2];
	len_type dims[2];
	len_type offset[2];
};

template<>
struct symphas::io::gp::gp_plotting_helper_specialized<1> : symphas::io::gp::gp_plotting_helper_impl<gp_plotting_helper_specialized<1>>
{
	using parent_type = gp_plotting_helper_impl<gp_plotting_helper_specialized<1>>;
	using parent_type::get_index;
	using parent_type::get_position;

	gp_plotting_helper_specialized(symphas::io::write_info const& winfo, symphas::grid_info const& ginfo) :
		h{}, pos0{}, len{}, dims{}, offset{}
	{
		h[0] = ginfo.INTERVAL_Xh;
		len[0] = ginfo.INTERVAL_Xc;
		offset[0] = ginfo.INTERVAL_X0 - ginfo.DOMAIN_X0;
		dims[0] = ginfo.DOMAIN_Xc;

		pos0[0] = ginfo.INTERVAL_X0;
	}

	len_type get_index(const len_type* coords) const
	{
		return ((coords[0] + offset[0]) % dims[0]);
	}

	double get_position(Axis ax, const len_type* coords) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coords[n] * h[n];
	}

	double get_position(Axis ax, len_type coord) const
	{
		iter_type n = symphas::axis_to_index(ax);
		return pos0[n] + coord * h[n];
	}

	size_t get_dim() const
	{
		return 1;
	}

	len_type get_len(iter_type n) const
	{
		return len[n];
	}


	double h[1];
	double pos0[1];
	len_type len[1];
	len_type dims[1];
	len_type offset[1];
};

template<>
struct symphas::io::gp::gp_plotting_helper_specialized<0> : symphas::io::gp::gp_plotting_helper_impl<gp_plotting_helper_specialized<0>>
{
	using parent_type = gp_plotting_helper_impl<gp_plotting_helper_specialized<0>>;
	using parent_type::get_index;
	using parent_type::get_position;

	gp_plotting_helper_specialized(symphas::io::write_info const& winfo, symphas::grid_info const& ginfo) :
		h{}, pos0{}, len{}, offset{}
	{
	}

	len_type get_index(const len_type* coords) const
	{
		return 0;
	}

	double get_position(Axis ax, const len_type* coords) const
	{
		return 0;
	}

	double get_position(Axis ax, len_type coord) const
	{
		return 0;
	}

	size_t get_dim() const
	{
		return 0;
	}

	len_type get_len(iter_type n) const
	{
		return 0;
	}

	double h[1];
	double pos0[1];
	len_type len[1];
	len_type offset[1];
};

symphas::io::gp::gp_plotting_helper* symphas::io::gp::new_helper(symphas::io::write_info const& winfo, symphas::grid_info const& ginfo)
{
	switch (ginfo.dimension())
	{
	case 3: return new gp_plotting_helper_specialized<3>(winfo, ginfo);
	case 2: return new gp_plotting_helper_specialized<2>(winfo, ginfo);
	case 1: return new gp_plotting_helper_specialized<3>(winfo, ginfo);
	default: return new gp_plotting_helper_specialized<0>(winfo, ginfo);
	}
}

symphas::io::gp::gp_plotting_helper* symphas::io::gp::new_helper(symphas::grid_info ginfo)
{
	//for (auto& [axis, interval] : ginfo)
	//{
	//	interval.domain_to_interval();
	//}
	return new_helper({}, ginfo);
}



DLLIO double symphas::io::gp::alignment::display_width = MULTI_OUTPUT_WIDTH_DEFAULT;
DLLIO double symphas::io::gp::alignment::mmr = MULTI_MARGIN_RATIO_DEFAULT;
DLLIO double symphas::io::gp::alignment::mas = MULTI_ALIGN_SEPARATION_DEFAULT;
DLLIO double symphas::io::gp::alignment::chr = COLORBOX_HEIGHT_RATIO_DEFAULT;
DLLIO double symphas::io::gp::alignment::cbp = COLORBOX_BOTTOM_PADDING_DEFAULT;

DLLIO std::pair<const char*, double*> symphas::io::gp::alignment::key_par_pairs[] = {
	{ STR(MULTI_OUTPUT_WIDTH_KEY), &symphas::io::gp::alignment::display_width },
	{ STR(MULTI_MARGIN_RATIO_KEY), &symphas::io::gp::alignment::mmr },
	{ STR(MULTI_ALIGN_SEPARATION_KEY), &symphas::io::gp::alignment::mas },
	{ STR(COLORBOX_HEIGHT_RATIO_KEY), &symphas::io::gp::alignment::chr },
	{ STR(COLORBOX_BOTTOM_PADDING_KEY), &symphas::io::gp::alignment::cbp } };



void symphas::io::gp::print_gp_header(int index, size_t id, symphas::grid_info const& ginfo, symphas::io::write_info const& winfo, FILE* f)
{
	static std::vector<std::tuple<std::string, size_t>> idlist;
	if (!params::single_output_file || (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir_str_ptr, id)) == idlist.end()))
	{
		fprintf(f, "%d ", ginfo.dimension());

		/* the output of the intervals always goes x, y, z
		 * keep checking size of dimension and print corresponding interval
		 */

		for (iter_type i = 0; i < ginfo.dimension(); ++i)
		{
			fprintf(f, "%d ", ginfo.at(symphas::index_to_axis(i)).get_count());
		}

		for (iter_type i = 0; i < ginfo.dimension(); ++i)
		{
			auto& interval = ginfo.intervals.at(symphas::index_to_axis(i));
			fprintf(f, "%lf %lf ", interval.domain_left(), interval.domain_right());
		}

		if (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir_str_ptr, id)) == idlist.end())
		{
			idlist.emplace_back(winfo.dir_str_ptr, id);
		}
	}
	fprintf(f, "%d", index);

	if (grid::has_subdomain(ginfo))
	{
		fprintf(f, " %d ", ginfo.dimension());
		for (auto const& [axis, interval] : ginfo)
		{
			fprintf(f, "%lf %lf ", interval.left(), interval.right());
		}
	}
	else
	{
		fprintf(f, " %c", CONFIG_OPTION_PREFIX_C);
	}

	fprintf(f, "\n");
}



/*
 * functions to write data from the grid to a file
 */
void symphas::io::gp::save_grid_plotting(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	if (ginfo.dimension() > 1)
	{
		for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
		{
			fprintf(f, "% 13d ", k);
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSX(i));
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < GP_HELPER_LENY; j++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSY(j));
				for (iter_type i = 0; i < GP_HELPER_LENX; i++)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					fprintf(f, "% 10lE ", grid[ii]);
				}
				fprintf(f, "\n");
			}
			if (k < GP_HELPER_LENZ - 1)
			{
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f % 10lE\n", GP_HELPER_POSX(i), grid[i]);
		}
		fprintf(f, "\n");
	}

	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	if (ginfo.dimension() > 1)
	{
		double
			dX = ginfo.INTERVAL_Xh,
			dY = ginfo.INTERVAL_Yh;

		for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
		{
			fprintf(f, "% 13d ", k);
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSX(i));
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < GP_HELPER_LENY; j++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSY(j));
				for (iter_type i = 0; i < GP_HELPER_LENX; i++)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					fprintf(f, "% 10lE ", std::abs(grid[ii]));
				}
				fprintf(f, "\n");
			}
			if (k < GP_HELPER_LENZ - 1)
			{
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f % 10lE\n", GP_HELPER_POSX(i), std::abs(grid[i]));
		}
		fprintf(f, "\n");
	}

	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const double_arr2* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		fprintf(f, "% 13d ", k);
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSX(i));
		}
		fprintf(f, "\n");

		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POSY(j));
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				fprintf(f, "% 10lE ", std::sqrt(grid[ii][0] * grid[ii][0] + grid[ii][1] * grid[ii][1]));
			}
			fprintf(f, "\n");
		}
		if (k < GP_HELPER_LENZ - 1)
		{
			fprintf(f, "\n");
		}
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		double z = GP_HELPER_POSZ(k);
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			double y = GP_HELPER_POSY(j);
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double x = GP_HELPER_POSX(i);
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				vector_t<3> v = grid[ii];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					x, y, z, dx, dy, dz, m);
			}
		}
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type j = 0; j < GP_HELPER_LENY; j++)
	{
		double y = GP_HELPER_POSY(j);
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			double x = GP_HELPER_POSX(i);
			iter_type ii = GP_HELPER_INDEX({ i, j });
			vector_t<2> v = grid[ii];

			double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
				dx = v.v[0] / m,
				dy = v.v[1] / m;

			fprintf(f,
				"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f\n",
				x, y, dx, dy, m);
		}
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


void symphas::io::gp::save_grid_plotting(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type i = 0; i < GP_HELPER_LENX; i++)
	{
		double x = GP_HELPER_POSX(i);
		iter_type ii = GP_HELPER_INDEX({ i });
		vector_t<1> v = grid[ii];

		fprintf(f,
			"%." AXIS_OUTPUT_ACCURACY_STR "f "
			"%" DATA_OUTPUT_ACCURACY_STR "f ",
			x, v.v[0]);
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		double z = GP_HELPER_POSZ(k);
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			double y = GP_HELPER_POSY(j);
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double x = GP_HELPER_POSX(i);
				iter_type ii = GP_HELPER_INDEX({ i, j });

				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii] + grid[2][ii] * grid[2][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m,
					dz = grid[2][ii] / m;

				fprintf(f,
					"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					x, y, z, dx, dy, dz, m);
			}
		}
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}

void symphas::io::gp::save_grid_plotting(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type j = 0; j < GP_HELPER_LENY; j++)
	{
		double y = GP_HELPER_POSY(j);
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			double x = GP_HELPER_POSX(i);
			iter_type ii = GP_HELPER_INDEX({ i, j });

			double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii]),
				dx = grid[0][ii] / m,
				dy = grid[1][ii] / m;

			fprintf(f,
				"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f "
				"%" DATA_OUTPUT_ACCURACY_STR "f\n",
				x, y, dx, dy, m);
		}
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


void symphas::io::gp::save_grid_plotting(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::gp::new_helper(winfo, ginfo);

	for (iter_type i = 0; i < GP_HELPER_LENX; i++)
	{
		double x = GP_HELPER_POSX(i);
		iter_type ii = GP_HELPER_INDEX({ i });

		fprintf(f,
			"%." AXIS_OUTPUT_ACCURACY_STR "f "
			"%" DATA_OUTPUT_ACCURACY_STR "f ",
			x, grid[0][ii]);
	}
	symphas::io::gp::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}








/*
 * functions to checkpoint the data to a file, separate from recording the data for plotting
 */



void symphas::io::gp::save_grid(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_interval_count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_interval_count();
	len_type L = ginfo.at(Axis::X).get_interval_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f,
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " ", grid[ii]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void symphas::io::gp::save_grid(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_interval_count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_interval_count();
	len_type L = ginfo.at(Axis::X).get_interval_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f,
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " ",
					grid[ii].real(), grid[ii].imag());
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void symphas::io::gp::save_grid(const double_arr2 *grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);


	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_interval_count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_interval_count();
	len_type L = ginfo.at(Axis::X).get_interval_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f, 
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " ",
					grid[ii][0], grid[ii][1]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void symphas::io::gp::save_grid(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = ginfo.at(Axis::Z).get_interval_count();
	len_type M = ginfo.at(Axis::Y).get_interval_count();
	len_type L = ginfo.at(Axis::X).get_interval_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + (j * L) + (k * L * M);
				vector_t<3> v = grid[ii];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n",
					dx, dy, dz, m);
			}
		}
	}

	fclose(f);
}

void symphas::io::gp::save_grid(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type j = 0; j < ginfo.at(Axis::Y).get_interval_count(); j++)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_interval_count(); i++)
		{
			iter_type ii = i + j * ginfo.at(Axis::X).get_interval_count();
			vector_t<2> v = grid[ii];

			double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
				dx = v.v[0] / m,
				dy = v.v[1] / m;

			fprintf(f,
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n",
				dx, dy, m);
		}
	}
	fclose(f);
}


void symphas::io::gp::save_grid(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, DataFileType::CHECKPOINT_DATA);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type i = 0; i < ginfo.at(Axis::X).get_interval_count(); i++)
	{
		fprintf(f, "%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n", grid[i].v[0]);
	}
	fclose(f);
}



void symphas::io::gp::save_grid(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = ginfo.at(Axis::Z).get_interval_count();
	len_type M = ginfo.at(Axis::Y).get_interval_count();
	len_type L = ginfo.at(Axis::X).get_interval_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + (j * L) + (k * L * M);

				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii] + grid[2][ii] * grid[2][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m,
					dz = grid[2][ii] / m;

				fprintf(f,
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
					"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n",
					dx, dy, dz, m);
			}
		}
	}

	fclose(f);
}

void symphas::io::gp::save_grid(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type j = 0; j < ginfo.at(Axis::Y).get_interval_count(); j++)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_interval_count(); i++)
		{
			iter_type ii = i + j * ginfo.at(Axis::X).get_interval_count();

			double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii]),
				dx = grid[0][ii] / m,
				dy = grid[1][ii] / m;

			fprintf(f,
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR " "
				"%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n",
				dx, dy, m);
		}
	}
	fclose(f);
}


void symphas::io::gp::save_grid(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, DataFileType::CHECKPOINT_DATA);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type i = 0; i < ginfo.at(Axis::X).get_interval_count(); i++)
	{
		fprintf(f, "%" CHECKPOINT_OUTPUT_ACCURACY_STR CHECKPOINT_OUTPUT_FORMAT_STR "\n", grid[0][i]);
	}
	fclose(f);
}

