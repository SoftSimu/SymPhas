
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
 * MODULE:  sol
 * PURPOSE: Elements related to defining a phase field system, including
 * data persistence and representing information. Also contains methods of
 * getting the correct grid data based on different grid implementations.
 *
 * ***************************************************************************
 */

#pragma once

#include <memory>


#include "gridfunctions.h"
#include "boundary.h"
#include "initialconditions.h"

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif



//**************************************************************************************************

#ifdef USING_IO
#include <thread>

//! Responsible for persisting the grid to file on a parallel thread.
/*!
 * Defines a grid which must be populated with values before its data is saved
 * to a file. The data of the grid is saved to file on a parallel thread, using
 * the save utility selected by the program.
 */
template<typename T>
struct WriteParallel
{
	WriteParallel(len_type len) : snapshot{ len }, done{ true }, thr{ nullptr } {}

	void write(symphas::io::write_info w, symphas::grid_info g) const
	{
		finish_last_write();
		done = false;
		thr = new std::thread{ &WriteParallel<T>::write_delegate, snapshot.values, w, g };
	}

	auto& get_snapshot() const
	{
		finish_last_write();
		return snapshot;
	}

	~WriteParallel()
	{
		finish_last_write();
	}

protected:

	WriteParallel() : WriteParallel(0) {}

	void finish_last_write() const
	{
		if (!done)
		{
			thr->join();
			delete thr;
			done = true;
		}
	}

	static void write_delegate(const T* values, symphas::io::write_info w, symphas::grid_info g)
	{
		symphas::io::save_grid_plotting(values, w, g);
	}


	Block<T> snapshot;

public:

	mutable bool done;			//!< Used to indicate whether the writing job is complete.
	mutable std::thread* thr;	//!< Thread used to in writing to file.
};
#endif

//! Manages information related to representing a phase field system.
/*!
 * Manages information related to representing a phase field system. A
 * phase field system is represented by the grid data and a unique number,
 * which is used to differentiate phase fields in the same model.
 */
struct SystemInfo
{
	//! Create an empty system information instance.
	SystemInfo() : info{ grid::dim_list() }, id{ 0 } {}


	//! Create a new system information instance.
	/*!
	 * Create a new system information instance from the given data.
	 * 
	 * \param info The grid information about the phase field system.
	 * \param id The unique phase field identifier, distinguishing it between
	 * phase fields in the same model.
	 */
	SystemInfo(symphas::grid_info const& info, size_t id) : info{ info }, id{ id } {}

	//! Returns the data representing the system.
	auto const& get_info() const
	{
		return info;
	}

	//! Returns the specific ID value associated with the system.
	auto const get_id() const
	{
		return id;
	}

protected:

	symphas::grid_info info;		//!< Contains the data specifying grid parameters.
	size_t id;						//!< The given number associated with the grid.


};

template<typename G>
struct SystemData;

//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot.
 *
 * The snapshot may be written to disk only if IO is enabled,
 * otherwise no write utility would be available. The snapshot and write utility
 * are available through the object WriteParallel, which is conditionally
 * compiled if IO is enabled.
 */
template<typename T, size_t D>
struct SystemData<Grid<T, D>> : Grid<T, D>, SystemInfo
{
	using SystemInfo::id;
	using SystemInfo::info;

	using Grid<T, D>::dims;
	using Grid<T, D>::len;
	using Grid<T, D>::as_grid;

	//! Create a system data instance.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	SystemData(symphas::interval_data_type const& vdata, size_t id) :
		Grid<T, D>{ grid::construct<::Grid, T, D>(vdata) }, 
		SystemInfo{ { vdata }, id } 
	{
		if (info.dimension() > 0)
		{
			for (iter_type i = info.dimension() - 1; i >= D; --i)
			{
				info.intervals.erase(symphas::index_to_axis(i));
			}
		}
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid.
	 */
	operator const Grid<T, D>&() const
	{
		return as_grid();
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid.
	 */
	operator Grid<T, D>&()
	{
		return as_grid();
	}

	//! Get the number of elements of true data.
	/*!
	 * The elements which are considered true data are the interior domain,
	 * that is, all non-boundary elements. Return the count of these elements.
	 */
	auto length() const
	{
		return len;
	}

	//! Copies the system data into the given array.
	/*!
	 * The values of the system data block are copied into a new one. The copy
	 * is performed point-wise for all data points.
	 */
	void persist(T* out) const
	{
		grid::copy(*this, out);
	}

	//! Copies the input data into the system values.
	/*!
	 * The values of the system data block are initialized from the
	 * given values, correctly transcribing all values.
	 */
	void fill(const T* in) const
	{
		grid::fill(in, *this);
	}

protected:

	SystemData() : Grid<T, D>{}, SystemInfo{ { {} }, 0 } {}
};

//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot. Specialization implementing
 * a system with boundaries.
 * 
 * See SystemData<Grid<T, D>>.
 */
template<typename T, size_t D>
struct SystemData<BoundaryGrid<T, D>> : BoundaryGrid<T, D>, SystemInfo
{
	using SystemInfo::id;
	using SystemInfo::info;

	using Grid<T, D>::dims;
	using Grid<T, D>::len;
	using BoundaryGrid<T, D>::as_grid;


	//! Create a system data instance.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	SystemData(symphas::interval_data_type const& vdata, size_t id) :
		BoundaryGrid<T, D>{ grid::construct<::BoundaryGrid, T, D>(vdata) },
		SystemInfo{ { vdata }, id }
	{
		if (!info.intervals.empty())
		{
			for (iter_type i = 0; i < D; ++i)
			{
				Axis ax = symphas::index_to_axis(i);
				auto& interval = info.intervals.at(ax);

				interval.set_count(dims[i] - 2 * BOUNDARY_DEPTH);
				interval.interval_to_domain();
			}
		}
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid with boundaries.
	 */
	operator const BoundaryGrid<T, D>& () const
	{
		return as_grid();
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid with boundaries.
	 */
	operator BoundaryGrid<T, D>& ()
	{
		return as_grid();
	}

	//! Get the number of elements of true data.
	/*!
	 * The elements which are considered true data are the interior domain,
	 * that is, all non-boundary elements. Return the count of these elements.
	 */
	auto length() const
	{
		return grid::length_interior<D>(dims);
	}

	//! Copies the system data into the given array.
	/*!
	 * The values of the system data block are copied into a new one. The copy
	 * is performed point-wise for all interior points, meaning points not on
	 * the boundary.
	 */
	void persist(T* out) const
	{
		grid::copy_interior(*this, out);
	}

	//! Copies the input data into the system values.
	/*!
	 * The values of the system data block are initialized from the
	 * given values, correctly transcribing all values.
	 */
	void fill(const T* in) const
	{
		grid::fill_interior(in, *this, dims);
	}

protected:

	SystemData() : BoundaryGrid<T, D>{}, SystemInfo{ { {} }, 0 } {}
};



//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot. Specialization implementing
 * a system with boundaries.
 *
 * See SystemData<Grid<T, D>>.
 */
template<typename T, size_t D>
struct SystemData<RegionalGrid<T, D>> : RegionalGrid<T, D>, SystemInfo
{
	using SystemInfo::id;
	using SystemInfo::info;

	using RegionalGrid<T, D>::region;
	using RegionalGrid<T, D>::dims;
	using RegionalGrid<T, D>::len;
	using RegionalGrid<T, D>::as_grid;


	//! Create a system data instance.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	SystemData(symphas::interval_data_type const& vdata, size_t id) :
		RegionalGrid<T, D>{ grid::construct<::RegionalGrid, T, D>(vdata) },
		SystemInfo{ { vdata }, id }
	{
		if (!info.intervals.empty())
		{
			for (iter_type i = 0; i < D; ++i)
			{
				Axis ax = symphas::index_to_axis(i);
				auto& interval = info.intervals.at(ax);
				double boundary = region.boundary_size * interval.width();

				interval.set_count(dims[i] - 2 * BOUNDARY_DEPTH);
				interval.set_interval(interval.left() + boundary, interval.right() - boundary);
			}
		}
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid with boundaries.
	 */
	operator const RegionalGrid<T, D>& () const
	{
		return as_grid();
	}

	//! Get the grid representation of this system.
	/*!
	 * Get the grid representation of this system, which will return the
	 * underlying grid with boundaries.
	 */
	operator RegionalGrid<T, D>& ()
	{
		return as_grid();
	}

	//! Get the number of elements of true data.
	/*!
	 * The elements which are considered true data are the interior domain,
	 * that is, all non-boundary elements. Return the count of these elements.
	 */
	auto length() const
	{
		return grid::length_interior<D>(region.dims);
	}

	//! Copies the system data into the given array.
	/*!
	 * The values of the system data block are copied into a new one. The copy
	 * is performed point-wise for all interior points, meaning points not on
	 * the boundary.
	 */
	void persist(T* out) const
	{
		auto interval = grid::get_iterable_domain(*this);
		auto it = symphas::data_iterator_region(as_grid(), interval);
		auto end = it + grid::length<D>(interval);
		while (it < end)
		{
			*out++ = *it++;
		}
	}

	//! Copies the input data into the system values.
	/*!
	 * The values of the system data block are initialized from the
	 * given values, correctly transcribing all values.
	 */
	void fill(const T* in) const
	{
		grid::fill_interior(in, *this, dims);
	}

protected:

	SystemData() : RegionalGrid<T, D>{}, SystemInfo{ { {} }, 0 } {}
};






template<typename G>
struct PersistentSystemData;

#ifdef USING_IO
#include "datalib.h"

//! Maintains grid data parameters as well as a snapshot.
/*!
 * Maintains grid data parameters, and contains implementations used in data
 * persistence.
 *
 * The snapshot may be written to disk only if IO is enabled,
 * otherwise no write utility would be available. The snapshot and write utility
 * are available through the object WriteParallel, which is conditionally
 * compiled if IO is enabled.
 */
template<template<typename, size_t> typename G, typename T, size_t D>
struct PersistentSystemData<G<T, D>> : SystemData<G<T, D>>
{
	using parent_type = SystemData<G<T, D>>;

	using parent_type::len;
	using parent_type::dims;
	using parent_type::id;
	using parent_type::info;
	using parent_type::persist;
	using parent_type::length;

	//! Create a system that can persist data to disk.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	PersistentSystemData(symphas::interval_data_type const& vdata, size_t id) :
		parent_type(vdata, id), writer{ length() } {}

	void persist() const
	{
		persist(writer.get_snapshot().values);
	}

	//! Writes the current snapshot to the disk if IO is enabled.
	/*!
	 * The current snapshot is written to disk if the IO macro keyword is
	 * defined. The snapshot should be updated using persist before writing.
	 */
	void write(const char* dir, iter_type index) const
	{
		symphas::io::write_info w{ dir, index, id, DataFileType::SOLUTION_DATA };
		symphas::grid_info g{ info.intervals };
		writer.write(w, g);
	}

	//! Writes the current snapshot to the disk if IO is enabled.
	/*!
	 * The current snapshot is written to disk if the IO macro keyword is
	 * defined. The snapshot should be updated using persist before writing.
	 * 
	 * \param dir The directory to which the data file is saved.
	 * \param name The name of the file to save to disk.
	 */
	void write(const char* dir, const char* name, iter_type index) const
	{
		char* write_loc = new char[std::strlen(dir) + std::strlen(name) + 2];
		snprintf(write_loc, BUFFER_LENGTH, "%s/%s", dir, name);

		symphas::io::write_info w{ write_loc, index, id, DataFileType::NAMED_DATA };
		symphas::grid_info g{ info.intervals };
		writer.write(w, g);

		delete[] write_loc;
	}

	//! Save the data information to disk.
	/*!
	 * The system data is saved to a disk.
	 * 
	 * \param dir The directory to which the data file is saved.
	 */
	inline void save(const char* dir, iter_type index) const
	{
		persist();
		write(dir, index);
	}

	//! Save the data information to disk with a name.
	/*!
	 * The system data is saved to a disk.
	 * 
	 * \param dir The directory to which the data file is saved.
	 * \param name The name of the file to save to disk.
	 */
	inline void save(const char* dir, const char* name, iter_type index) const
	{
		persist();
		write(dir, name, index);
	}


	//! Get a copy of the snapshot.
	/*!
	 * This updates the member snapshot data with the current system values,
	 * and copies the snapshot n the snapshot and returns it.
	 */
	Block<T> get_snapshot() const
	{
		persist();
		return writer.get_snapshot();
	}

	using field_type = symphas::FieldAxis<D, T*>;
	operator field_type() const
	{
		auto values = std::shared_ptr<T[]>(new T[length()]);
		persist(values.get());

		return {
			std::shared_ptr<axis_nd_t<D>[]>(std::move(symphas::lib::new_system_axis_list<D>(info.intervals))),
			values,
			length() };
	}

	field_type as_field() const
	{
		return *this;
	}

protected:

	PersistentSystemData() : PersistentSystemData{ {}, 0 } {}

	WriteParallel<T> writer;
};


//! Maintains grid data parameters as well as a snapshot.
/*!
 * Maintains grid data parameters, and contains implementations used in data
 * persistence.
 *
 * The snapshot may be written to disk only if IO is enabled,
 * otherwise no write utility would be available. The snapshot and write utility
 * are available through the object WriteParallel, which is conditionally
 * compiled if IO is enabled.
 */
template<typename T, size_t D>
struct PersistentSystemData<RegionalGrid<T, D>> : SystemData<RegionalGrid<T, D>>
{
	using parent_type = SystemData<RegionalGrid<T, D>>;

	using parent_type::len;
	using parent_type::dims;
	using parent_type::id;
	using parent_type::info;
	using parent_type::persist;
	using parent_type::length;
	using parent_type::values;
	using parent_type::region;

	//! Create a system that can persist data to disk.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	PersistentSystemData(symphas::interval_data_type const& vdata, size_t id) :
		parent_type(vdata, id) {}


	void write(symphas::io::write_info w, symphas::grid_info g) const
	{

		for (auto& [axis, interval] : g)
		{
			w.intervals[axis][0] = interval.domain_left();
			w.intervals[axis][1] = interval.domain_right();
		}

		for (auto& [axis, interval] : g)
		{
			interval.domain_to_interval();

			double offset = region.boundary_size * interval.width();
			interval.set_domain(interval.domain_left() - offset, interval.domain_right() + offset);
		}

		g.update_strides();
		symphas::io::save_grid_plotting(values, w, g);
	}

	//void persist() const
	//{
	//	persist(writer.get_snapshot().values);
	//}

	//! Writes the current snapshot to the disk if IO is enabled.
	/*!
	 * The current snapshot is written to disk if the IO macro keyword is
	 * defined. The snapshot should be updated using persist before writing.
	 */
	void write(const char* dir, iter_type index) const
	{
		symphas::io::write_info w{ dir, index, id, DataFileType::SOLUTION_DATA };
		symphas::grid_info g{ info };
		write(w, g);
	}

	//! Writes the current snapshot to the disk if IO is enabled.
	/*!
	 * The current snapshot is written to disk if the IO macro keyword is
	 * defined. The snapshot should be updated using persist before writing.
	 *
	 * \param dir The directory to which the data file is saved.
	 * \param name The name of the file to save to disk.
	 */
	void write(const char* dir, const char* name, iter_type index) const
	{
		char* write_loc = new char[std::strlen(dir) + std::strlen(name) + 2];
		snprintf(write_loc, BUFFER_LENGTH, "%s/%s", dir, name);

		symphas::io::write_info w{ write_loc, index, id, DataFileType::NAMED_DATA };
		symphas::grid_info g{ info };
		write(w, g);

		delete[] write_loc;
	}

	//! Save the data information to disk.
	/*!
	 * The system data is saved to a disk.
	 *
	 * \param dir The directory to which the data file is saved.
	 */
	inline void save(const char* dir, iter_type index) const
	{
		write(dir, index);
	}

	//! Save the data information to disk with a name.
	/*!
	 * The system data is saved to a disk.
	 *
	 * \param dir The directory to which the data file is saved.
	 * \param name The name of the file to save to disk.
	 */
	inline void save(const char* dir, const char* name, iter_type index) const
	{
		write(dir, name, index);
	}


	//! Get a copy of the snapshot.
	/*!
	 * This updates the member snapshot data with the current system values,
	 * and copies the snapshot n the snapshot and returns it.
	 */
	Block<T> get_snapshot() const
	{
		Block<T> snapshot(length());
		persist(snapshot.values);
		return snapshot;
	}

	using field_type = symphas::FieldAxis<D, T*>;
	operator field_type() const
	{
		auto values = std::shared_ptr<T[]>(new T[length()]);
		persist(values.get());

		return {
			std::shared_ptr<axis_nd_t<D>[]>(std::move(symphas::lib::new_system_axis_list<D>(info.intervals))),
			values,
			length() };
	}

	field_type as_field() const
	{
		return *this;
	}

protected:

	PersistentSystemData() : PersistentSystemData{ {}, 0 } {}
};

#include "writedefines.h"

namespace symphas::io
{
	template<typename G>
	void save_grid(PersistentSystemData<G> const& sys);
	template<typename G>
	void save_grid(PersistentSystemData<G> const& sys, const char* file);
}


//! Default argument overload of symphas::io::save_grid(T*, symphas::io::write_info, symphas::grid_info).
template<typename G>
void symphas::io::save_grid(PersistentSystemData<G> const& sys)
{
	auto data = sys.get_snapshot();
	symphas::io::save_grid(
		data.values, symphas::io::write_info{ 0, sys.get_id() }, 
		sys.dims, sys.info.dimension());
}


//! Default argument overload of symphas::io::save_grid(T*, symphas::io::write_info, symphas::grid_info).
template<typename G>
void symphas::io::save_grid(PersistentSystemData<G> const& sys, const char* file)
{
	auto data = sys.get_snapshot();
	symphas::io::save_grid(
		data.values, symphas::io::write_info{ file, 0, sys.get_id(), DataFileType::NAMED_DATA },
		sys.dims, sys.info.dimension());
}



#else

//! Maintains grid data parameters and allows taking a data snapshot.
/*!
 * Maintains grid data parameters, and allows taking a data snapshot.
 */
template<template<typename, size_t> typename G, typename T, size_t D>
struct PersistentSystemData<G<T, D>> : SystemData<G<T, D>>
{
	using parent_type = SystemData<G<T, D>>;

	using parent_type::SystemData;
	using parent_type::info;

	using parent_type::persist;
	using parent_type::length;


	Block<T> get_snapshot() const
	{
		Block<T> data(length());
		persist(data.values);
		return std::move(data);
	}

};


#endif





// *************************************************************************************************



namespace symphas
{
	
	enum ModelModifiers
	{
		PLOT_DEFAULT,
		PLOT_MAX,
		PLOT_MIN,
		PLOT_CONTOURS
	};

	//! Manages a list of time steps and times they are applied at.
	struct time_step_list
	{


		double* t_dts;			//!< The times at which the time steps apply.
		size_t dts_len;			//!< The number of time steps;

		time_step_list(double dt = 1.0, double time = 0) : t_dts{ new double[2] { time, dt } }, dts_len{ 1 } {}
		time_step_list(double* dts, double* t_dts, size_t dts_len) :
			t_dts{ (dts_len > 0) ? new double[dts_len * 2] {} : nullptr }, dts_len{ dts_len }
		{
			for (iter_type i = 0; i < dts_len; ++i)
			{
				this->t_dts[i * 2] = t_dts[i];
				this->t_dts[i * 2 + 1] = dts[i];
			}
		}

		time_step_list(time_step_list const& other) :
			t_dts{ (other.dts_len > 0) ? new double[other.dts_len * 2] {} : nullptr }, dts_len{ other.dts_len }
		{
			std::copy(other.t_dts, other.t_dts + other.dts_len * 2, t_dts);
		}

		time_step_list(time_step_list&& other) : time_step_list(nullptr, nullptr, 0) 
		{
			swap(*this, other);
		}

		time_step_list& operator=(time_step_list other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(time_step_list& first, time_step_list& second)
		{
			using std::swap;
			swap(first.t_dts, second.t_dts);
			swap(first.dts_len, second.dts_len);
		}

		//! Provide the time step used in the solution of the problem.
		void set_time_step(double dt, double time = 0, bool insert = true);

		//! Clear all the time steps and replace it with one time step set at the default value.
		void clear_time_steps(double default_dt = 1.);

		//! Provide the time step used in the solution of the problem.
		void set_time_steps(const double* dts, const double* t_dts, size_t dts_len);

		inline auto get_time_steps() const
		{
			symphas::lib::array_container<double> steps(dts_len);
			for (iter_type i = 0; i < dts_len; ++i)
			{
				steps[i] = t_dts[i * 2 + 1];
			}
			return steps;
		}

		inline auto get_times_of_steps() const
		{
			symphas::lib::array_container<double> times(dts_len);
			for (iter_type i = 0; i < dts_len; ++i)
			{
				times[i] = t_dts[i * 2];
			}
			return times;
		}

		inline size_t get_num_time_steps() const
		{
			return dts_len;
		}

		double get_time_step(double time = 0) const;

		auto begin() const
		{
			return symphas::lib::zip(
				symphas::stride_type(t_dts, 2, dts_len), 
				symphas::stride_type(t_dts + 1, 2, dts_len)).begin();
		}

		auto end() const
		{
			return symphas::lib::zip(
				symphas::stride_type(t_dts, 2, dts_len), 
				symphas::stride_type(t_dts + 1, 2, dts_len)).end();
		}


		~time_step_list()
		{
			delete[] t_dts;
		}
	};

    //! Representation of the problem parameters for a phase field model.
    /*!
     * Represents all the information about a phase field problem. Contains data
     * about the intervals of the axes, all the boundary conditions, as well as
     * the initial conditions. It also has information about time interval, which
     * depending on the solver implementation, is used to initialize the solver.
     *
     * Each data of the problem parameters is stored as an array, and the order
     * of the elements in the array corresponds to the order of initialization
     * of the phase fields in a defined model. The corresponding order is how the
     * phase fields are initialized when passed to a model.
     *
     * The number of systems that the problem parameters manage is passed to the
     * constructor, and cannot be changed.
     */
    struct problem_parameters_type
    {
    protected:
    
        symphas::init_data_type* tdata;			//!< Data for the initial conditions of each system.
        symphas::interval_data_type* vdata;		//!< Data for the intervals of each system.
        symphas::b_data_type* bdata;			//!< Boundary data for each system.
		len_type* num_fields;
    
        //! Number of elements in the data arrays.
        /*!
         * Number of elements in the data arrays which must match the number of
         * systems in the problem these parameters represent.
         */
		size_t len;
		size_t num_fields_len;
    
		symphas::time_step_list dt_list;		//!< The list of time steps;
		size_t* modifiers;						//!< A bit string of modifiers.
    
        problem_parameters_type() : 
            tdata{ nullptr }, vdata{ nullptr }, bdata{ nullptr },
			num_fields{ new len_type[1]{ len_type(0) } },
			len{ 0 }, num_fields_len{ 0 }, modifiers{ nullptr }, time{ TIME_INIT }, index{ params::start_index } {}
    
    public:
    
        double time;					//!< The current time.
        iter_type index;				//!< The current index;
    
    
        //! Initialize the problem parameters corresponding to \p len systems.
        /*!
            * Initialize the problem parameters corresponding to \p len systems.
            *
            * \param len The number of systems that the problem parameters refer to.
            */
        problem_parameters_type(const size_t len) : 
			tdata{ (len > 0) ? new symphas::init_data_type[len]{} : nullptr },
			vdata{ (len > 0) ? new symphas::interval_data_type[len]{} : nullptr }, 
			bdata{ (len > 0) ? new symphas::b_data_type[len]{} : nullptr },
			num_fields{ new len_type[1]{ len_type(len) } },
			len{ len }, num_fields_len{ 1 }, dt_list{},
			modifiers{ new size_t[1]{} }, time { TIME_INIT }, index{ params::start_index } {}
    
		problem_parameters_type(len_type* num_fields, size_t num_fields_len) : 
			problem_parameters_type((num_fields != nullptr && num_fields_len > 0) ? std::reduce(num_fields, num_fields + num_fields_len) : 0)
		{
			if (num_fields_len > 0)
			{
				delete[] this->num_fields;
				this->num_fields = new len_type[num_fields_len]{};
				this->num_fields_len = num_fields_len;
				std::copy(num_fields, num_fields + num_fields_len, this->num_fields);

				modifiers = new size_t[num_fields_len]{};
			}
		}

        problem_parameters_type(problem_parameters_type const& other);
    
        problem_parameters_type(problem_parameters_type&& other) noexcept : problem_parameters_type()
        {
            swap(*this, other);
        }
        problem_parameters_type& operator=(problem_parameters_type other)
        {
            swap(*this, other);
            return *this;
        }

		void set_modifier(ModelModifiers m, iter_type i = 0)
		{
			modifiers[i] |= (1ull << static_cast<size_t>(m));
		}

		void set_modifier(const char* str, iter_type i = 0);

		void unset_modifier(ModelModifiers m, iter_type i = 0)
		{
			modifiers[i] &= ~(1ull << static_cast<size_t>(m));
		}

		bool check_modifier(ModelModifiers m, iter_type i = 0) const
		{
			return ((1ull << static_cast<size_t>(m)) & modifiers[i]);
		}

		void set_modifiers(const char* const* str, size_t n);

		ModelModifiers get_modifier(iter_type i = 0) const
		{
			if (modifiers[i] == 0)
			{
				return ModelModifiers::PLOT_DEFAULT;
			}
			else
			{
				bool flag = false;
				size_t value = modifiers[i];
				int n = -1;

				while (!flag)
				{
					flag = (value & 1ull);
					value >>= 1ull;
					n += 1;
				}
				return static_cast<ModelModifiers>(n);
			}
		}
    
        //! Set the initial data of the problem parameters.
        /*
            * Set the initial data of the problem parameters, to be applied to one
            * or more systems in a phase field problem initialization. If the
            * length of the given list is smaller than the number of systems
            * represented by this object, the remaining elements will be initialized
            * to the first of this list.
            *
            * \param tdata_set The list of initial data to apply. The i-th element
            * in the given list applies to the i-th phase field to the given phase
            * field problem.
            * \param n The length of the list.
            */
        void set_initial_data(const symphas::init_data_type* tdata_set, size_t n);
        inline void set_initial_data(const symphas::init_data_type* tdata_set)
        {
            set_initial_data(tdata_set, len - 1);
        }
        inline void set_initial_data(symphas::init_data_type const& tdata_set, iter_type i)
        {
            if (i < len)
            {
                tdata[i] = tdata_set;
            }
            else
            {
                throw std::out_of_range("given initial data element larger than list length\n");
            }
        }
		void set_initial_data(symphas::init_entry_type const& tentry_set, iter_type i)
		{
			symphas::init_data_type tdata_set{ { Axis::NONE, tentry_set } };
			set_initial_data(tdata_set, i);
		}

        //! Set the interval data of the problem parameters.
        /*
         * Set the interval data of the problem parameters, to be applied to one
         * or more systems in a phase field problem initialization. The intervals
         * must match the system dimension. If the
         * length of the given list is smaller than the number of systems
         * represented by this object, the remaining elements will be initialized
         * to the first of this list.
         *
         * \param vdata_set The list of interval data to apply. The i-th element
         * in the given list applies to the i-th phase field to the given phase
         * field problem.
         * \param n The length of the list.
         */
        void set_interval_data(const symphas::interval_data_type* vdata_set, size_t n);
        inline void set_interval_data(const symphas::interval_data_type* vdata_set)
        {
            set_interval_data(vdata_set, len - 1);
        }
        inline void set_interval_data(symphas::interval_data_type const& vdata_set, iter_type i)
        {
            if (i < len)
            {
                vdata[i] = vdata_set;
            }
            else
            {
                throw std::out_of_range("given interval data element larger than list length\n");
            }
        }
    
    
        //! Set the boundary data of the problem parameters.
        /*
         * Set the boundary data of the problem parameters, to be applied to one
         * or more systems in a phase field problem initialization. The boundary
         * dimension must match the system dimension. If the
         * length of the given list is smaller than the number of systems
         * represented by this object, the remaining elements will be initialized
         * to the first of this list.
         *
         * \param bdata_set The list of boundary data to apply. The i-th element
         * in the given list applies to the i-th phase field to the given phase
         * field problem.
         * \param n The length of the list.
         */
        void set_boundary_data(const symphas::b_data_type* bdata_set, size_t n);
        inline void set_boundary_data(const symphas::b_data_type* bdata_set)
        {
            set_boundary_data(bdata_set, len - 1);
        }
        inline void set_boundary_data(symphas::b_data_type const& bdata_set, iter_type i)
        {
            if (i < len)
            {
                bdata[i] = bdata_set;
            }
            else
            {
                throw std::out_of_range("given boundary data element larger than list length\n");
            }
        }

		//! Set a list of number of fields, when the number of fields is runtime-configuration.
		/*
		 * Set the number of fields of each field type that will be initialized in a model 
		 * When a model has multiple order parameters types that can have their total count set 
		 * from the configuration, the count will be taken from the parameter set here. The order
		 * of the lengths corresponds to the order of the configurable fields. If the
		 * length of the given list is smaller than the number of systems
		 * represented by this object, the remaining elements will be initialized
		 * to the first of this list.
		 *
		 * \param num_fields_set The list of field lengths to set. The i-th element
		 * in the given list applies to the i-th phase-field type for configurable field counts.
		 * \param n The length of the list.
		 */
		void set_num_fields(len_type* num_fields_set, size_t n);
		inline void set_num_fields(len_type* num_fields_set)
		{
			set_num_fields(num_fields_set, num_fields_len - 1);
		}
		inline void set_num_fields(len_type num_fields_set, iter_type i)
		{
			if (i < num_fields_len)
			{
				len_type diff = num_fields_set - num_fields[i];
				num_fields[i] = num_fields_set;
				len += diff;
			}
			else
			{
				throw std::out_of_range("given field length element larger than list length\n");
			}
		}
		void set_num_fields(len_type num_fields_set)
		{
			set_num_fields(num_fields_set, 0);
		}

        //! Provide the time step used in the solution of the problem.
		void set_time_step(double dt, double time = 0);

		//! Provide the time step used in the solution of the problem.
		void set_time_steps(const double* dts, const double* t_dts, size_t dts_len);

		//! Provide the time step used in the solution of the problem.
		inline void set_time_steps(symphas::time_step_list const& list)
		{
			dt_list = list;
		}

		inline void clear_time_steps(double default_dt = 1.);
    
        //! Get the list of initial data.
        inline const symphas::init_data_type* get_initial_data() const
        {
            return tdata;
        }
    
        //! Get the list of interval data.
        inline const symphas::interval_data_type* get_interval_data() const
        {
            return vdata;
        }
    
        //! Get the list of boundary data.
        inline const symphas::b_data_type* get_boundary_data() const
        {
            return bdata;
        }

		//! Get the list of initial data.
		inline symphas::init_data_type* get_initial_data()
		{
			return tdata;
		}

		//! Get the list of interval data.
		inline symphas::interval_data_type* get_interval_data()
		{
			return vdata;
		}

		//! Get the list of boundary data.
		inline symphas::b_data_type* get_boundary_data()
		{
			return bdata;
		}

		inline const len_type* get_num_fields() const
		{
			return num_fields;
		}

		inline const size_t get_num_fields_len() const
		{
			return num_fields_len;
		}

		inline len_type* get_num_fields()
		{
			return num_fields;
		}

		inline auto get_time_steps() const
		{
			return dt_list.get_time_steps();
		}

		inline auto get_times_of_steps() const
		{
			return dt_list.get_times_of_steps();
		}

		inline size_t get_num_time_steps() const
		{
			return dt_list.get_num_time_steps();
		}

		inline const auto& get_time_step_list() const
		{
			return dt_list;
		}

		inline auto& get_time_step_list()
		{
			return dt_list;
		}

		inline symphas::grid_info get_grid_info_entry(iter_type n)
		{
			if (n < len)
			{
				return { vdata[n] };
			}
			else
			{
				return { nullptr, get_dimension() };
			}
		}

		inline decltype(auto) get_dims_entry(iter_type n)
		{
			if (n < len)
			{
				return get_grid_info_entry(n).get_dims();
			}
			else
			{
				return grid::dim_list(nullptr, get_dimension());
			}
		}

		using modifier_list_type = symphas::lib::array_container<symphas::ModelModifiers>;

		inline auto get_modifiers() const
		{
			modifier_list_type modifier_list(num_fields_len);
			for (iter_type i = 0; i < num_fields_len; ++i)
			{
				modifier_list[i] = get_modifier(i);
			}
			return modifier_list;
		}

    
        //! Get the length of the data elements.
        inline const size_t length() const
        {
            return len;
        }
    
        //! Provide the time step used in the solution of the problem.
        inline double get_time_step(double time = 0) const
        {
            return dt_list.get_time_step(time);
        }
    
        //! Get the dimension represented by the parameters.
        /*!
         * Return the size of the given interval, which by is the first interval
         * by default. The size represents the dimension of the problem.
         */
        inline const size_t get_dimension(iter_type i = 0) const
        {
            return vdata[i].size();
        }
    
    
        //! Sets interval length of all intervals data to the given dimensions. 
        /*!
         * Modifies the length of each interval in the intervals data list so that
         * they correspond to the given dimensions. Each interval will still be
         * in the same spatial range, but the number of discrete points will be
         * made to match what is given. The order of the dimensions given
         * must correspond to \f$x\f$, \f$y\f$ and \f$z\f$.
         *
         * This can only be called once all problem_parameters_type::len systems
         * have been initialized, otherwise it will result in an error.
         *
         * \param dims The dimensions to which all systems are set.
         */
        void equalize_discretization(const len_type* dims);
    
        //! Sets interval length of all intervals data to that at the given index. 
        /*!
         * See problem_parameters_type::equalize_discretization(). Takes the
         * intervals data from the list corresponding to index \p i and equalizes
         * all other intervals data to match its number of discrete points.
         *
         * \param i The index of the intervals data used to set the intervals of
         * all other intervals data elements.
         */
        void equalize_discretization(iter_type i);
    
        friend void swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second);
    
        void extend(size_t new_len)
        {
            if (new_len > len)
            {
                problem_parameters_type pp{ new_len };
                pp.set_boundary_data(get_boundary_data(), len);
                pp.set_initial_data(get_initial_data(), len);
                pp.set_interval_data(get_interval_data(), len);
                pp.set_time_steps(dt_list);
                pp.time = time;
                pp.index = index;
                swap(pp, *this);
            }
        }


		auto begin()
		{
			return symphas::lib::zip(
				symphas::stride_type(tdata, length()),
				symphas::stride_type(vdata, length()),
				symphas::stride_type(bdata, length())).end();
		}

		auto end()
		{
			return symphas::lib::zip(
				symphas::stride_type(tdata, length()), 
				symphas::stride_type(vdata, length()), 
				symphas::stride_type(bdata, length())).end();
		}

		auto begin() const
		{
			return symphas::lib::zip(
				symphas::stride_type(tdata, length()),
				symphas::stride_type(vdata, length()),
				symphas::stride_type(bdata, length())).end();
		}

		auto end() const
		{
			return symphas::lib::zip(
				symphas::stride_type(tdata, length()),
				symphas::stride_type(vdata, length()),
				symphas::stride_type(bdata, length())).end();
		}


		static problem_parameters_type init_default(size_t dimension)
		{
			problem_parameters_type problem_parameters(1);
			symphas::b_data_type bdata(dimension);
			symphas::init_data_type tdata(Inside::NONE);
			symphas::interval_data_type vdata(dimension);

			problem_parameters.set_boundary_data(&bdata, 1);
			problem_parameters.set_initial_data(&tdata, 1);
			problem_parameters.set_interval_data(&vdata, 1);

			return problem_parameters;
		}

        ~problem_parameters_type();
    
    };
    
    void swap(problem_parameters_type& first, problem_parameters_type& second);
}

void swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second);

/*
 * system initialization helper functions
 */

#define INITIAL_CONDITION_INVALID_MSG "the given initial condition algorithm is not valid\n"

namespace symphas::internal
{

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata, 
		Grid<T, D>& data, [[maybe_unused]] symphas::grid_info* info, grid::region_interval<D> const& region, [[maybe_unused]] size_t id)
	{
		auto ginfo = InitialConditions<T, D>{ tdata, info->intervals, data.dims }.initialize(data.values, region, id);
		if (!symphas::is_valid(ginfo))
		{
			fprintf(SYMPHAS_ERR, INITIAL_CONDITION_INVALID_MSG);
		}
	}

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata,
		Grid<any_vector_t<T, D>, D>& data, [[maybe_unused]] symphas::grid_info* info, grid::region_interval<D> const& region, [[maybe_unused]] size_t id)
	{
		any_vector_t<T, D>* values = new any_vector_t<T, D>[data.len];
		auto ginfo = InitialConditions<any_vector_t<T, D>, D>{ tdata, info->intervals, data.dims }.initialize(values, region, id);
		if (!symphas::is_valid(ginfo))
		{
			fprintf(SYMPHAS_ERR, INITIAL_CONDITION_INVALID_MSG);
		}

		for (iter_type n = 0; n < D; ++n)
		{
#			pragma omp parallel for
			for (iter_type i = 0; i < data.len; ++i)
			{
				data.axis(symphas::index_to_axis(n))[i] = values[i][n];
			}
		}
	}


	template<typename T, size_t D>
	void update_info_for_regional(RegionalGrid<T, D> const& data, symphas::grid_info const& info, symphas::grid_info* to_update)
	{
		for (const auto& [axis, interval] : info)
		{
			auto boundary = interval.width() * data.region.boundary_size;
			(*to_update)[axis].set_interval(interval.left() + boundary, interval.right() - boundary);
		}
	}

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata,
		RegionalGrid<T, D>& data, [[maybe_unused]] symphas::grid_info* info, grid::region_interval<D> const& region, [[maybe_unused]] size_t id)
	{
		auto ginfo = InitialConditions<T, D>{ tdata, info->intervals, data.dims }.initialize(data.values, region, id);
		if (!symphas::is_valid(ginfo))
		{
			fprintf(SYMPHAS_ERR, INITIAL_CONDITION_INVALID_MSG);
		}
		else
		{
			bool update_region = false;
			for (auto [ax0, interval0] : ginfo)
			{
				auto interval1 = (info->intervals)[ax0];
				if (interval0.left() != interval1.left() || interval0.right() != interval1.right())
				{
					update_region = true;
				}
			}
			if (update_region)
			{
				grid::resize_adjust_region(data, ginfo);
				update_info_for_regional(data, ginfo, info);
			}
		}
	}

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata,
		RegionalGrid<any_vector_t<T, D>, D>& data, [[maybe_unused]] symphas::grid_info* info, grid::region_interval<D> const& region, [[maybe_unused]] size_t id)
	{
		any_vector_t<T, D>* values = new any_vector_t<T, D>[data.len];
		auto ginfo = InitialConditions<any_vector_t<T, D>, D>{ tdata, info->intervals, data.dims }.initialize(values, region, id);
		if (!symphas::is_valid(ginfo))
		{
			fprintf(SYMPHAS_ERR, INITIAL_CONDITION_INVALID_MSG);
		}
		else
		{
			for (iter_type n = 0; n < D; ++n)
			{
#				pragma omp parallel for
				for (iter_type i = 0; i < data.len; ++i)
				{
					data.axis(symphas::index_to_axis(n))[i] = values[i][n];
				}
			}
			bool update_region = false;
			for (auto [ax0, interval0] : ginfo)
			{
				auto interval1 = (info->intervals)[ax0];
				if (interval0.left() != interval1.left() || interval0.right() != interval1.right())
				{
					update_region = true;
				}
			}
			if (update_region)
			{
				grid::resize_adjust_region(data, ginfo);
				update_info_for_regional(data, ginfo, info);
			}
			grid::resize_adjust_region(data, ginfo);
			update_info_for_regional(data, ginfo, info);
		}
	}

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata,
		Grid<T, D>& data, [[maybe_unused]] symphas::grid_info* info, [[maybe_unused]] size_t id)
	{
		populate_tdata(tdata, data, info, grid::region_interval<D>(info->get_dims()), id);
	}

	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata,
		RegionalGrid<T, D>& data, [[maybe_unused]] symphas::grid_info* info, [[maybe_unused]] size_t id)
	{
		populate_tdata(tdata, data, info, grid::region_interval<D>(info->get_dims()), id);
	}
}




