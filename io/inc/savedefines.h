
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
 * PURPOSE: Defines basic parts of the saving (writing) functionality for
 * arbitrary data sets which are not phase field data. These are typically
 * used for generic data, such as data that was computed about the phase
 * field. The data is accepted in different formats.
 *
 * ***************************************************************************
 */

#pragma once

#include "write.h"
#include "fieldconv.h"
#include "datalib.h"

//! \cond

#ifdef PROC_EXPORTS
#define DLLPROC DLLEXPORT
#else
#define DLLPROC DLLIMPORT
#endif

#define DEFAULT_PTS_NAME "pts-data-output"
#define DEFAULT_ABS_NAME "abs-data-output"
#define DEFAULT_VEC_NAME "vec-data-output"

//! \endcond

#define PLOT_DISPLAY_SIZE_WIN "500,350"		//!< Typical plot display size for Windows.
#define PLOT_DISPLAY_SQSIZE_WIN "500,500"	//!< Typical square plot display size for Windows.
#define PLOT_DISPLAY_SIZE_LATEX "6,4.5"		//!< Typical plot display size for a Latex document.
#define PLOT_DISPLAY_SQSIZE_LATEX "4.5,4.5"	//!< Typical square plot display size for a Latex document.
#define PLOT_OUTPUT_ACCURACY DATA_OUTPUT_ACCURACY
#define PLOT_OUTPUT_ACCURACY_STR STR(PLOT_OUTPUT_ACCURACY)

namespace symphas::io
{

	void write_data_entry(FILE* out, int data);
	void write_data_entry(FILE* out, scalar_t data);
	void write_data_entry(FILE* out, complex_t data);
	template<typename T, size_t N>
	void write_data_entry(FILE* out, T(&data)[N]);
	void write_data_entry(FILE* out, axis_2d_type data);
	void write_data_entry(FILE* out, axis_3d_type data);

	template<typename... Ts, size_t... Is>
	void write_data_entry(FILE* out, std::tuple<Ts...> const& data, std::index_sequence<Is...>);
	template<typename... Ts>
	void write_data_entry(FILE* out, std::tuple<Ts...> const& data);


	//! Default method for saving a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<typename T>
	void save_pts_default(
		std::vector<std::tuple<iter_type, axis_1d_type, T>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0, 
		const char* name = DEFAULT_PTS_NAME);

	//! Default method for saving a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<size_t N>
	void save_pts_default(
		std::vector<std::tuple<iter_type, axis_1d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_PTS_NAME);



	//! Default method for saving a list of real-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<typename T>
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, T>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);

	//! Default method for saving a list of real-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<size_t N>
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);

	//! Default method for saving a list of real-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<typename... Ts>
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, std::tuple<Ts...>>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);

	//! Default method for saving a list of real-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);

	//! Default method for saving a list of iterated 2D points numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, axis_2d_type>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);

	//! Default method for saving a list of iterated 3D points numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	void save_abs_default(
		std::vector<std::pair<axis_1d_type, axis_3d_type>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_ABS_NAME);


	//! Default method for saving a list of vector-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<size_t D>
	void save_vec_default(
		std::vector<std::pair<typename axis_nd<D>::type, scalar_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_VEC_NAME);

	//! Default method for saving a list of vector-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<size_t D>
	void save_vec_default(
		std::vector<std::pair<typename axis_nd<D>::type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_VEC_NAME);


	//! Default method for saving a list of vector-valued numbers to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 * \param abs_name A string name of the data which is used to make the name of the file.
	 */
	template<size_t D, size_t N>
	void save_vec_default(
		std::vector<std::pair<typename axis_nd<D>::type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0,
		const char* name = DEFAULT_VEC_NAME);


	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_pts(
		std::vector<std::tuple<iter_type, axis_1d_type, scalar_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_pts_default(data, dir, index, id);
	}

	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void, size_t N = 2>
	void save_pts(
		std::vector<std::tuple<iter_type, axis_1d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_pts_default(data, dir, index, id);
	}


	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_pts(
		std::vector<std::tuple<iter_type, axis_1d_type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_pts_default(data, dir, index, id);
	}


	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename T, typename S = void>
	void save_pts(Field<iter_type*, std::pair<iter_type, T>*> field,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		std::vector<std::tuple<iter_type, axis_1d_type, T>> data(field.length());
		for (iter_type i = 0; i < field.length(); ++i)
		{
			data[i] = std::make_tuple(field.data_x()[i], field.data_y()[i].first, field.data_y()[i].second);
		}

		save_pts(data, dir, index, id);
	}


	//! Save a list of real-valued numbers to a file.
	/*!
	 * Save a list of real-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_abs(
		std::vector<std::pair<axis_1d_type, scalar_t>> const& data, 
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	template<typename S = void>
	void save_abs(
		std::vector<std::pair<axis_1d_type, complex_t>> const& data, 
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	template<typename S = void, size_t N = 2>
	void save_abs(
		std::vector<std::pair<axis_1d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	template<typename S = void, typename... Ts>
	void save_abs(
		std::vector<std::pair<axis_1d_type, std::tuple<Ts...>>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	template<typename S = void>
	void save_abs(
		std::vector<std::pair<axis_1d_type, std::vector<std::pair<scalar_t, scalar_t>>>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		throw;
	}

	template<typename S = void>
	void save_abs(
		std::vector<std::pair<axis_1d_type, axis_2d_type>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	template<typename S = void>
	void save_abs(
		std::vector<std::pair<axis_1d_type, axis_3d_type>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_abs_default(data, dir, index, id);
	}

	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename T, typename S = void>
	void save_abs(FieldAxis<1, T*> field,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		std::vector<std::pair<axis_1d_type, T>> data(field.length());
		for (iter_type i = 0; i < field.length(); ++i)
		{
			data[i] = std::make_pair(field.data_x()[i], field.data_y()[i]);
		}

		save_abs(data, dir, index, id);
	}

	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename T, typename S = void>
	void save_abs(Field<scalar_t*, T*> field,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		std::vector<std::pair<axis_1d_type, T>> data(field.length());
		for (iter_type i = 0; i < field.length(); ++i)
		{
			data[i] = std::make_pair(field.data_x()[i], field.data_y()[i]);
		}

		save_abs(data, dir, index, id);
	}

	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename T, typename S = void>
	void save_abs(T* field, len_type len,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		std::vector<std::pair<axis_1d_type, T>> data(len);
		for (iter_type i = 0; i < len; ++i)
		{
			data[i] = std::make_pair(i, field[i]);
		}

		save_abs(data, dir, index, id);
	}



	//! Save a list of 1-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 1-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_1d_type, scalar_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<1>(symphas::lib::to_scalar_field<1>(data), dir, index, id);
	}

	//! Save a list of 1-dimensional vector-valued complex numbers to a file.
	/*!
	 * Save a list of 1-dimensional vector-valued complex numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_1d_type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<1>(data, dir, index, id);
	}

	//! Save a list of 2-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 2-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_2d_type, scalar_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<2>(symphas::lib::to_scalar_field<2>(data), dir, index, id);
	}

	//! Save a list of 2-dimensional vector-valued complex numbers to a file.
	/*!
	 * Save a list of 2-dimensional vector-valued complex numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_2d_type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<2>(data, dir, index, id);
	}

	//! Save a list of 3-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 3-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_3d_type, scalar_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<3>(symphas::lib::to_scalar_field<3>(data), dir, index, id);
	}

	//! Save a list of 3-dimensional vector-valued complex numbers to a file.
	/*!
	 * Save a list of 3-dimensional vector-valued complex numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void>
	void save_vec(
		std::vector<std::pair<axis_3d_type, complex_t>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<3>(data, dir, index, id);
	}




	//! Save a list of 1-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 1-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void, size_t N = 2>
	void save_vec(
		std::vector<std::pair<axis_1d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<1>(symphas::lib::to_scalar_field<1>(data), dir, index, id);
	}

	//! Save a list of 2-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 2-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void, size_t N = 2>
	void save_vec(
		std::vector<std::pair<axis_2d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<2>(symphas::lib::to_scalar_field<2>(data), dir, index, id);
	}

	//! Save a list of 3-dimensional vector-valued numbers to a file.
	/*!
	 * Save a list of 3-dimensional vector-valued numbers to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<typename S = void, size_t N = 2>
	void save_vec(
		std::vector<std::pair<axis_3d_type, scalar_t[N]>> const& data,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		save_vec_default<3>(symphas::lib::to_scalar_field<3>(data), dir, index, id);
	}



	//! Save a list of points to a file.
	/*!
	 * Save a list of points to a file.
	 *
	 * \param data The data that is persisted to the file.
	 * \param dir The directory of the file that is created.
	 * \param index The identifying index of the data that is saved, e.g. the time index.
	 * \param id An identifying number corresponding to the type of data, e.g. the phase-field id
	 * from which the data comes from.
	 */
	template<size_t D, typename T, typename S = void>
	void save_vec(FieldAxis<D, T*> field,
		const char* dir = ".", iter_type index = 0, size_t id = 0)
	{
		std::vector<std::pair<axis_nd_t<D>, T>> data(field.length());
		for (iter_type i = 0; i < field.length(); ++i)
		{
			data[i] = std::make_pair(field.data_x()[i], field.data_y()[i]);
		}

		save_vec(data, dir, index, id);
	}



	DLLIO extern const char* adefset[];			//!< Default gnuplot plotting configuration for single values.
	DLLIO extern const char* defset[];			//!< Default gnuplot plotting configuration for a 2d plot.
	DLLIO extern const char* defset2d[];		//!< Default gnuplot plotting configuration for a 2d point field.
	DLLIO extern const char* defset3d[];		//!< Default gnuplot plotting configuration for a 3d point field.
	DLLIO extern const char* defsettup[];		//!< Default gnuplot plotting configuration for a tuple.


}



/*
 *
 *
 *
 *
 *
 *
 *
 *
 * implementation of save functions
 */

template<typename T>
void symphas::io::save_pts_default(std::vector<std::tuple<iter_type, axis_1d_type, T>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [i, x, y] : data)
	{
		fprintf(def, "%.6d %.6E %.6E\n", i, x, y);
	}
	fclose(def);

	T ranges[4];
	symphas::lib::fill_sorted_ranges(data, ranges);

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	const char* sets[] = { adefset[0], adefset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, rangestr);
}

template<size_t N>
void symphas::io::save_pts_default(std::vector<std::tuple<iter_type, axis_1d_type, scalar_t[N]>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [i, x, y] : data)
	{
		fprintf(def, "%.6d %.6E", i, x);
		for (iter_type n = 0; n < N; ++n)
		{
			fprintf(def, " %.6E", y[n]);
		}
		fprintf(def, "\n");
	}
	fclose(def);

	const char* sets[] = { adefset[0], adefset[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, symphas::io::gp::ranges_auto);
}


template<typename T>
void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, T>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6E %.6E\n", x, y);
	}
	fclose(def);


	axis_coord_t ranges[2];
	T extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);
	T dy = std::max(1e-2, std::max((extrema[1] - extrema[0]) * 0.2, std::abs(extrema[0] * 0.1)));

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, ranges[0], ranges[1], extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { adefset[0], adefset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, rangestr);
}

inline void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, complex_t>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);
	scalar_t *rdata = new scalar_t[data.size()];
	
	scalar_t* it = rdata;
	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6f %.6E %.6E %.6E\n", x, std::real(y), std::imag(y), std::abs(y));
		*it++ = std::abs(y);
	}
	fclose(def);

	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(rdata, static_cast<len_type>(data.size()), extrema);
	scalar_t dy = std::max(1e-2, std::max((extrema[1] - extrema[0]) * 0.2, std::abs(extrema[0] * 0.1)));

	char rangestr[BUFFER_LENGTH];
	snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_xy, data.front().first, data.back().first, extrema[0] - dy, extrema[1] + dy);
	const char* sets[] = { adefset[0], adefset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, rangestr);
}

template<size_t N>
void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, scalar_t[N]>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6E", x);
		for (iter_type n = 0; n < N; ++n)
		{
			fprintf(def, " %.6E", y[n]);
		}
		fprintf(def, "\n");

	}
	fclose(def);

	const char* sets[] = { adefset[0], adefset[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, symphas::io::gp::ranges_auto);
}

inline void symphas::io::write_data_entry(FILE* out, int data)
{
	fprintf(out, "%d", data);
}

inline void symphas::io::write_data_entry(FILE* out, scalar_t data)
{
	fprintf(out, "%." PLOT_OUTPUT_ACCURACY_STR "lf", data);
}

inline void symphas::io::write_data_entry(FILE* out, complex_t data)
{
	write_data_entry(out, data.real());
	fprintf(out, " ");
	write_data_entry(out, data.imag());
}

template<typename T, size_t N>
void symphas::io::write_data_entry(FILE* out, T(&data)[N])
{
	for (iter_type i = 0; i < N; ++i)
	{
		write_data_entry(out, data[i]);
		if (i < N - 1)
		{
			fprintf(out, " ");
		}
	}
}

inline void symphas::io::write_data_entry(FILE* out, axis_2d_type data)
{
	write_data_entry(out, data[0]);
	fprintf(out, " ");
	write_data_entry(out, data[1]);
}

inline void symphas::io::write_data_entry(FILE* out, axis_3d_type data)
{
	write_data_entry(out, data[0]);
	fprintf(out, " ");
	write_data_entry(out, data[1]);
	fprintf(out, " ");
	write_data_entry(out, data[2]);
}

template<typename... Ts, size_t... Is>
void symphas::io::write_data_entry(FILE* out, std::tuple<Ts...> const& data, std::index_sequence<Is...>)
{
	((write_data_entry(out, std::get<Is>(data)), fprintf(out, " ")), ...);
}

template<typename... Ts>
void symphas::io::write_data_entry(FILE* out, std::tuple<Ts...> const& data)
{
	write_data_entry(out, data, std::make_index_sequence<sizeof...(Ts)>{});
}

template<typename... Ts>
void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, std::tuple<Ts...>>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%.6E ", x);
		write_data_entry(def, y);
		fprintf(def, "\n");

	}
	fclose(def);

	const char* sets[] = { defsettup[0], defsettup[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, symphas::io::gp::ranges_auto);
}

inline void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, axis_2d_type>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%." PLOT_OUTPUT_ACCURACY_STR "E", x);
		fprintf(def, " %." PLOT_OUTPUT_ACCURACY_STR "E %." PLOT_OUTPUT_ACCURACY_STR "E\n", y[0], y[1]);

	}
	fclose(def);

	const char* sets[] = { defset2d[0], defset2d[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, symphas::io::gp::ranges_auto);
}


inline void symphas::io::save_abs_default(std::vector<std::pair<axis_1d_type, axis_3d_type>> const& data, const char* dir, int index, size_t id, const char* abs_name)
{
	FILE* def = symphas::io::open_data_file(dir, abs_name, index, id, DataFileType::POSTPROC_DATA);

	for (auto const& [x, y] : data)
	{
		fprintf(def, "%." PLOT_OUTPUT_ACCURACY_STR "E", x);
		fprintf(def, " %." PLOT_OUTPUT_ACCURACY_STR "E %." PLOT_OUTPUT_ACCURACY_STR "E % ." PLOT_OUTPUT_ACCURACY_STR "E\n", y[0], y[1], y[2]);
	}
	fclose(def);

	const char* sets[] = { defset3d[0], defset3d[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, abs_name, index, id, symphas::io::gp::ranges_auto);
}



namespace symphas::internal
{
	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<1>, scalar_t>> const& data, FILE* f)
	{
		auto L = symphas::lib::get_dimensions<1>(data);
		auto it = data.begin();
		fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[0].first);
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "% 10E ", it++->second);
		}
		fprintf(f, "\n");
		fclose(f);
	}

	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<2>, scalar_t>> const& data, FILE* f)
	{
		auto [L, M] = symphas::lib::get_dimensions<2>(data)._2();

		fprintf(f, "% 10d ", 0);
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", data[i].first[0]);
		}
		fprintf(f, "\n");

		for (iter_type j = 0; j < M; j++)
		{
			fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[j * L].first[1]);
			for (iter_type i = 0; i < L; i++)
			{
				auto v = data[j * L + i].second;
				fprintf(f, "% 10E ", v);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
		fclose(f);

	}

	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<3>, scalar_t>> const& data, FILE *f)
	{
		auto [L, M, N] = symphas::lib::get_dimensions<3>(data)._3();
		for (iter_type k = 0; k < N; k++)
		{
			fprintf(f, "% 10d ", k);
			for (iter_type i = 0; i < L; i++)
			{
				fprintf(f, "% 13." AXIS_OUTPUT_ACCURACY_STR "f ", data[i + k * L * M].first[0]);
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < M; j++)
			{
				fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[j * L + k * L * M].first[1]);
				for (iter_type i = 0; i < L; i++)
				{
					auto v = data[k * L * M + j * L + i].second;
					fprintf(f, "% 10E ", v);
				}
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}


	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<1>, complex_t>> const& data, FILE* f)
	{
		auto L = symphas::lib::get_dimensions<1>(data);
		auto it = data.begin();
		fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[0].first);
		for (iter_type i = 0; i < L; i++)
		{
			auto v = it++->second;
			fprintf(f, "% 10E,% 10E ", v.real(), v.imag());
		}
		fprintf(f, "\n");
		fclose(f);
	}

	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<2>, complex_t>> const& data, FILE* f)
	{
		auto [L, M] = symphas::lib::get_dimensions<2>(data)._2();
		
		fprintf(f, "% 10d ", 0);
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "% 27." AXIS_OUTPUT_ACCURACY_STR "f ", data[i].first[0]);
		}
		fprintf(f, "\n");

		for (iter_type j = 0; j < M; j++)
		{
			fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[j * L].first[1]);
			for (iter_type i = 0; i < L; i++)
			{
				auto v = data[j * L + i].second;
				fprintf(f, "% 10E,% 10E ", v.real(), v.imag());
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
		fclose(f);

	}

	inline void save_vec_default_print(std::vector<std::pair<axis_nd_t<3>, complex_t>> const& data, FILE* f)
	{
		auto [L, M, N] = symphas::lib::get_dimensions<3>(data)._3();
		for (iter_type k = 0; k < N; k++)
		{
			fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[k * L * M].first[2]);
			for (iter_type i = 0; i < L; i++)
			{
				fprintf(f, "% 27." AXIS_OUTPUT_ACCURACY_STR "f ", data[i + k * L * M].first[0]);
			}
			fprintf(f, "\n");

			for (iter_type j = 0; j < M; j++)
			{
				fprintf(f, "% 10." AXIS_OUTPUT_ACCURACY_STR "f ", data[j * L + k * L * M].first[1]);
				for (iter_type i = 0; i < L; i++)
				{
					auto v = data[k * L * M + j * L + i].second;
					fprintf(f, "% 10E,% 10E ", v.real(), v.imag());
				}
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
}

template<size_t D>
void symphas::io::save_vec_default(std::vector<std::pair<axis_nd_t<D>, scalar_t>> const& data, const char* dir, int index, size_t id, const char* vec_name)
{
	FILE* def = symphas::io::open_data_file(dir, vec_name, index, id, DataFileType::POSTPROC_DATA);
	symphas::internal::save_vec_default_print(data, def);

	axis_coord_t ranges[D * 2];
	scalar_t extrema[2];
	symphas::lib::fill_sorted_ranges(data, ranges, extrema);

	char rangestr[BUFFER_LENGTH];

	if constexpr (D > 1)
	{
		snprintf(rangestr, BUFFER_LENGTH, symphas::io::gp::ranges_axis_2d, ranges[0], ranges[1], ranges[2], ranges[3]);
	}
	const char* sets[] = { defset[0], defset[1] };

	symphas::io::write_postproc_plot_file(sets, 2, dir, vec_name, index, id, rangestr);
}


template<size_t D>
void symphas::io::save_vec_default(std::vector<std::pair<axis_nd_t<D>, complex_t>> const& data, const char* dir, int index, size_t id, const char* vec_name)
{
	FILE* def = symphas::io::open_data_file(dir, vec_name, index, id, DataFileType::POSTPROC_DATA);
	symphas::internal::save_vec_default_print(data, def);

	const char* sets[] = { defset[0], defset[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, vec_name, index, id, symphas::io::gp::ranges_auto);
}


template<size_t D, size_t N>
void symphas::io::save_vec_default(std::vector<std::pair<typename axis_nd<D>::type, scalar_t[N]>> const& data, const char* dir, int index, size_t id, const char* vec_name)
{
	FILE* def = symphas::io::open_data_file(dir, vec_name, index, id, DataFileType::POSTPROC_DATA);
	symphas::internal::save_vec_default_print(data, def);
	const char* sets[] = { defset[0], defset[1] };
	symphas::io::write_postproc_plot_file(sets, 2, dir, vec_name, index, id, symphas::io::gp::ranges_auto);
}



