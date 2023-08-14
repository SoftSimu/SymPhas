
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
 *
 * MODULE:  io
 * PURPOSE: Defines the text output functionality for plain text files.
 *
 * ***************************************************************************
 */


#pragma once


#include "writedefines.h"
#include "writegp.h"

// \cond
#define PIXEL_MULTIPLIER 4
#define MAX_RESOLUTION 2048
// \endcond

namespace symphas::io::gp
{
	//! Defines elements used in input and output for the text format.
	/*!
	 * The format of input/output is based on using it with the gnuplot
	 * utility.
	 */
	namespace mov {}
}


namespace symphas::io::gp::mov
{

	constexpr char gnuset_1[] = R"~(

reset

set x2tics scale 0
set y2tics scale 0

unset x2tics
unset y2tics
unset xtics
unset ytics
unset key
unset title
unset colorbox


set xrange [%f:%f]
set yrange [*:*]
")~";


	constexpr char gnuset[] = R"~(

reset

set x2tics scale 0
set y2tics scale 0

unset x2tics
unset y2tics
unset xtics
unset ytics
unset key
unset title
unset colorbox

set cbrange [-1:1]

set xrange [%f:%f]
set yrange [%f:%f]


)~";


	constexpr char gnu_term[] = R"~(

set terminal png size %d,%d

set tmargin at screen 0.99
set bmargin at screen 0.01
set lmargin at screen 0.01
set rmargin at screen 0.99


)~";

	constexpr char gnu_end[] = R"~(
unset output
)~";




	constexpr char gnu_1d[] = R"~(
	plot "%s" index i using 1:2 with lines title "%s")~";

	constexpr char gnu_s[] = R"~(
	plot "%s" index i matrix nonuniform with image title "%s")~";

	constexpr char gnu_2v[] = R"~(
	plot "%s" index i using 1:2:3:4:5 with vectors head size 1,20,60 filled lc palette title "%s")~";


	constexpr char gnu_3v[] = R"~(
	plot "%s" index i using 1:2:3:4:5:6:7 with vectors head size 1,20,60 filled lc palette title "%s")~";




	template<size_t D, typename S0>
	inline void plot_fmt(const char*&) {}

	template<>
	inline void plot_fmt<1, scalar_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_1d;
	}

	template<>
	inline void plot_fmt<2, scalar_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_s;
	}

	template<>
	inline void plot_fmt<3, scalar_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_s;
	}

	template<>
	inline void plot_fmt<1, complex_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_1d;
	}

	template<>
	inline void plot_fmt<2, complex_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_s;
	}

	template<>
	inline void plot_fmt<3, complex_t>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_s;
	}

	template<>
	inline void plot_fmt<2, vector_t<2>>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_2v;
	}

	template<>
	inline void plot_fmt<3, vector_t<3>>(const char*& gnu_set)
	{
		gnu_set = symphas::io::gp::mov::gnu_3v;
	}



	//! Get the string corresponding to the gnuplot input format for the data.
	template<size_t D, typename... S, size_t... Is>
	void get_plot_fmt(const char* (*gnu_set), std::index_sequence<Is...>)
	{
		((..., plot_fmt<D, S>(gnu_set[Is])));
	}
	//! Get the string corresponding to the gnuplot input format for the data.
	template<size_t D, typename... S>
	void get_plot_fmt(const char* (*gnu_set))
	{
		get_plot_fmt<D, S...>(gnu_set, std::make_index_sequence<sizeof...(S)>{});
	}


	template<size_t D>
	void print_plot_ranges(FILE* gnu, symphas::interval_data_type const& intervals)
	{
		fprintf(gnu, gnuset,
			DOMAIN_X0, DOMAIN_Xn - INTERVAL_Xh,
			DOMAIN_Y0, DOMAIN_Yn - INTERVAL_Yh
		);
	}


	template<>
	inline void print_plot_ranges<1>(FILE* gnu, symphas::interval_data_type const& intervals)
	{
		fprintf(gnu, gnuset_1,
			DOMAIN_X0, DOMAIN_Xn - INTERVAL_Xh
		);
	}


	//! Writes a data array to a file for plotting.
	/*!
	 * The given data is written to a file set by \p winfo. The information
	 * about the grid and the grid type is used to select a writing method. The
	 * output format and amount of data written is chosen to allow the resulting
	 * file to be ingested by a plotting program.
	 *
	 * \param grid The data which is written to the file.
	 * \param winfo Information about the file that is written.
	 * \param ginfo Information about the grid.
	 */
	template<typename value_type>
	void save_grid_plotting(value_type&& grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
	{
		symphas::io::gp::save_grid_plotting(std::forward<value_type>(grid), winfo, ginfo);
	}


	//! Writes a data array to a file.
	/*!
	 * The given data is written to a file set by \p winfo. The information
	 * about the grid and the grid type is used to select a writing method.
	 * 
	 * \param grid The data which is written to the file.
	 * \param winfo Information about the file that is written.
	 * \param ginfo Information about the grid.
	 */
	template<typename value_type>
	void save_grid(value_type &&grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
	{
		symphas::io::gp::save_grid(std::forward<value_type>(grid), winfo, ginfo);
	}


	template<size_t D, typename Sp, typename... S>
	void write_plot_config(Model<D, Sp, S...> const& model, const char* directory, const char* const* names, SaveParams const& save)
	{
		FILE* gnu;
		char plot_name[BUFFER_LENGTH_L2],
			data_loc[BUFFER_LENGTH];		// the data file location
		snprintf(plot_name, BUFFER_LENGTH_L2, PHASEFIELD_PLOT_LOC_FMT, directory, PHASEFIELD_DATA_NAME);

		if ((gnu = fopen(plot_name, "w")) == 0)
		{
			symphas::lib::make_directory_for_file(plot_name);
			if ((gnu = fopen(plot_name, "w")) == 0)
			{
				fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, plot_name);
				exit(ERR_CODE_FILE_OPEN);
			}
		}
		const auto intervals = model.template system<0>().get_info().intervals;
		print_plot_ranges<D>(gnu, intervals);

		fprintf(gnu, gnu_term, 
			std::min(MAX_RESOLUTION, PIXEL_MULTIPLIER * model.template system<0>().get_info().DOMAIN_Xc),
			std::min(MAX_RESOLUTION, PIXEL_MULTIPLIER * model.template system<0>().get_info().DOMAIN_Yc));

		// copy over the names of the systems for given system number
		for (size_t id = 0; id < sizeof...(S); ++id)
		{
			size_t saves = save.num_saves();
			const char* sets[sizeof...(S)];

			if (saves > 1)
			{
				if (params::single_output_file)
				{
					symphas::io::copy_data_file_name(DATA_DIR_RELATIVE_PLOT, save.get_stop(), id, DataFileType::SOLUTION_DATA, data_loc);
					get_plot_fmt<D, S...>(sets);

					fprintf(gnu, "do for [i=0:%zd] { \n", saves - 1);
					fprintf(gnu, "\tset output '" PHASEFIELD_DATA_NAME POSTFIX_ID_FMT "_'.i.'.png'", id);
					fprintf(gnu, sets[id], data_loc, names[id]);
					fprintf(gnu, "\n\tunset output\n}\n");
				}
				else
				{
					symphas::io::gp::get_plot_fmt<D, S...>(sets);

					for (iter_type index = save.get_start(), i = 0; i < saves; index = save.next_save(index), ++i)
					{
						symphas::io::copy_data_file_name(DATA_DIR_RELATIVE_PLOT, index, id, DataFileType::SOLUTION_DATA, data_loc);
						fprintf(gnu, "\tset output '" PHASEFIELD_DATA_NAME POSTFIX_ID_FMT "_%d.png'", id, i);
						fprintf(gnu, sets[id], data_loc, 0, names[id]);
						fprintf(gnu, "\n\tunset output\n");
					}
				}
			}
			else
			{
				/* in the case of a single output file, the correct index from the
				 * list of saved data needs to be loaded
				 */
				iter_type index = (params::single_output_file)
					? std::max(0, save.num_saves() - 1)
					: 0;

				symphas::io::gp::get_plot_fmt<D, S...>(sets);

				fprintf(gnu, "set output " PHASEFIELD_DATA_NAME POSTFIX_ID_FMT "_%d.png'", id, index);
				fprintf(gnu, sets[id], data_loc, index, names[id]);
				fprintf(gnu, "\nunset output\n");
			}

		}
		fclose(gnu);
	}

}


#undef PIXEL_MULTIPLIER
#undef MAX_RESOLUTION


