
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
	static std::vector<std::pair<std::string, size_t>> idlist;
	if (!params::single_output_file || (std::find(idlist.begin(), idlist.end(), id) == idlist.end()))
	{
		fprintf(f, "%d ", ginfo.dimension());

		/* the output of the intervals always goes x, y, z
		 * keep checking size of dimension and print corresponding interval
		 */

		for (iter_type i = 0; i < ginfo.dimension(); ++i)
		{
			fprintf(f, "%d ", ginfo.at(symphas::index_to_axis(i)).count());
		}

		for (iter_type i = 0; i < ginfo.dimension(); ++i)
		{
			auto& interval = ginfo.intervals.at(symphas::index_to_axis(i));
			fprintf(f, "%lf %lf ", interval.left(), interval.right());
		}

		if (std::find(idlist.begin(), idlist.end(), std::make_pair(winfo.dir_str_ptr, id)) == idlist.end())
		{
			idlist.emplace_back(winfo.dir_str_ptr, id);
		}
	}
	fprintf(f, "%d", index);
	fprintf(f, "\n");
}




/*
 * functions to write data from the grid to a file
 */
template<>
void symphas::io::gp::save_grid_plotting<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	if (ginfo.dimension() > 1)
	{
		double
			dX = ginfo.INTERVAL_Xh,
			dY = ginfo.INTERVAL_Yh;

		len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
		len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
		len_type L = ginfo.at(Axis::X).count();

		for (iter_type k = 0; k < N; k++)
		{
			fprintf(f, "% 13d ", k);
			for (iter_type i = 0; i < L; i++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_X0 + dX * i);
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < M; j++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_Y0 + dY * j);
				for (iter_type i = 0; i < L; i++)
				{
					iter_type ii = i + j * L + k * L * M;
					fprintf(f, "% 10lE ", grid[ii]);
				}
				fprintf(f, "\n");
			}
			if (k < N - 1)
			{
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		double dX = ginfo.INTERVAL_Xh;
		len_type L = ginfo.at(Axis::X).count();
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f % 10lE\n", ginfo.INTERVAL_X0 + dX * i, grid[i]);
		}
		fprintf(f, "\n");
	}

	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gp::save_grid_plotting<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);



	if (ginfo.dimension() > 1)
	{
		double
			dX = ginfo.INTERVAL_Xh,
			dY = ginfo.INTERVAL_Yh;

		len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
		len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
		len_type L = ginfo.at(Axis::X).count();

		for (iter_type k = 0; k < N; k++)
		{
			fprintf(f, "% 13d ", k);
			for (iter_type i = 0; i < L; i++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_X0 + dX * i);
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < M; j++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_Y0 + dY * j);
				for (iter_type i = 0; i < L; i++)
				{
					iter_type ii = i + j * L + k * L * M;
					fprintf(f, "% 10lE ", std::abs(grid[ii]));
				}
				fprintf(f, "\n");
			}
			if (k < N - 1)
			{
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		double dX = ginfo.INTERVAL_Xh;
		len_type L = ginfo.at(Axis::X).count();
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f % 10lE\n", ginfo.INTERVAL_X0 + dX * i, std::abs(grid[i]));
		}
		fprintf(f, "\n");
	}

	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gp::save_grid_plotting<double[2]>(const double (*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double
		dX = ginfo.INTERVAL_Xh,
		dY = ginfo.INTERVAL_Yh;

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		fprintf(f, "% 13d ", k);
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_X0 + dX * i);
		}
		fprintf(f, "\n");

		for (iter_type j = 0; j < M; j++)
		{
			fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", ginfo.INTERVAL_Y0 + dY * j);
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f, "% 10lE ", std::sqrt(grid[ii][0] * grid[ii][0] + grid[ii][1] * grid[ii][1]));
			}
			fprintf(f, "\n");
		}
		if (k < N - 1)
		{
			fprintf(f, "\n");
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gp::save_grid_plotting<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double
		dX = ginfo.INTERVAL_Xh,
		dY = ginfo.INTERVAL_Yh,
		dZ = ginfo.INTERVAL_Zh;

	len_type N = ginfo.at(Axis::Z).count();
	len_type M = ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		double z = ginfo.INTERVAL_Z0 + dZ * k;
		for (iter_type j = 0; j < M; j++)
		{
			double y = ginfo.INTERVAL_Y0 + dY * j;
			for (iter_type i = 0; i < L; i++)
			{
				double x = ginfo.INTERVAL_X0 + dX * i;
				iter_type ii = i + (j * L) + (k * L * M);
				vector_t<3> v = grid[ii];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f\n",
					x, y, z, dx, dy, dz, m);
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gp::save_grid_plotting<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double
		dX = ginfo.INTERVAL_Xh,
		dY = ginfo.INTERVAL_Yh;

	for (iter_type j = 0; j < ginfo.at(Axis::Y).count(); j++)
	{
		double y = ginfo.INTERVAL_Y0 + dY * j;
		for (iter_type i = 0; i < ginfo.at(Axis::X).count(); i++)
		{
			double x = ginfo.INTERVAL_X0 + dX * i;
			iter_type ii = i + j * ginfo.at(Axis::X).count();
			vector_t<2> v = grid[ii];

			double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
				dx = v.v[0] / m,
				dy = v.v[1] / m;

			fprintf(f,
				"%." AXIS_OUTPUT_ACCURACY_STR "f %." AXIS_OUTPUT_ACCURACY_STR "f "
				"%." DATA_OUTPUT_ACCURACY_STR "f "
				"%." DATA_OUTPUT_ACCURACY_STR "f "
				"%." DATA_OUTPUT_ACCURACY_STR "f\n",
				x, y, dx, dy, m);
		}
	}
	fclose(f);
}


template<>
void symphas::io::gp::save_grid_plotting<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dX = ginfo.INTERVAL_Xh;

	for (iter_type i = 0; i < ginfo.at(Axis::X).count(); i++)
	{
		double x = ginfo.INTERVAL_X0 + dX * i;
		iter_type ii = i;
		vector_t<1> v = grid[ii];

		fprintf(f,
			"%." AXIS_OUTPUT_ACCURACY_STR "f "
			"%." DATA_OUTPUT_ACCURACY_STR "f ",
			x, v.v[0]);
	}
	fprintf(f, "\n");
	fclose(f);
}






/*
 * functions to checkpoint the data to a file, separate from recording the data for plotting
 */



template<>
void symphas::io::gp::save_grid<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f,
					"%." DATA_OUTPUT_ACCURACY_STR "f ", grid[ii]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

template<>
void symphas::io::gp::save_grid<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f,
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f ",
					grid[ii].real(), grid[ii].imag());
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

template<>
void symphas::io::gp::save_grid<double[2]>(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);


	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;
				fprintf(f, 
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f ", 
					grid[ii][0], grid[ii][1]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

template<>
void symphas::io::gp::save_grid<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	len_type N = ginfo.at(Axis::Z).count();
	len_type M = ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

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
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, dz, m);
			}
		}
	}

	fclose(f);
}

template<>
void symphas::io::gp::save_grid<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type j = 0; j < ginfo.at(Axis::Y).count(); j++)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).count(); i++)
		{
			iter_type ii = i + j * ginfo.at(Axis::X).count();
			vector_t<2> v = grid[ii];

			double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
				dx = v.v[0] / m,
				dy = v.v[1] / m;

			fprintf(f,
				"%." DATA_OUTPUT_ACCURACY_STR "f "
				"%." DATA_OUTPUT_ACCURACY_STR "f "
				"%." DATA_OUTPUT_ACCURACY_STR "f\n", 
				dx, dy, m);
		}
	}
	fclose(f);
}


template<>
void symphas::io::gp::save_grid<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, DataFileType::CHECKPOINT_DATA);
	print_gp_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type i = 0; i < ginfo.at(Axis::X).count(); i++)
	{
		fprintf(f, "%." DATA_OUTPUT_ACCURACY_STR "f\n", grid[i].v[0]);
	}
	fclose(f);
}






