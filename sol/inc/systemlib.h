
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
#include "datalib.h"

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
	mutable bool done;			//!< Used to indicate whether the writing job is complete.
	mutable std::thread* thr;	//!< Thread used to in writing to file.

	WriteParallel(len_type len) : snapshot{ len }, done{ true }, thr{ nullptr } {}

	void write(symphas::io::write_info w, symphas::grid_info g) const
	{
		finish_last_write();
		done = false;
		thr = new std::thread{ &symphas::io::save_grid_plotting<T>, snapshot.values, w, g };
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

	Block<T> snapshot;
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
	SystemInfo() : info{ {} }, id{ 0 } {}


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

	using Grid<T, D>::values;
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
		for (iter_type i = info.dimension() - 1; i >= D; --i)
		{
			info.intervals.erase(symphas::index_to_axis(i));
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
	auto data_len() const
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
		std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par, 
#endif
			values, values + len, out);
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

	using Grid<T, D>::values;
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
		for (iter_type i = 0; i < D; ++i)
		{
			Axis ax = symphas::index_to_axis(i);
			auto& interval = info.intervals.at(ax);

			interval.set_interval_count(
				interval.left(), 
				interval.right(),
				dims[i] - 2 * THICKNESS);
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
	auto data_len() const
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


protected:

	SystemData() : BoundaryGrid<T, D>{}, SystemInfo{ { {} }, 0 } {}
};







template<typename G>
struct PersistentSystemData;

#ifdef USING_IO


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
	using parent_type::data_len;

	//! Create a system that can persist data to disk.
	/*!
	 * Data for a phase field system is generated, resulting in a new
	 * grid instance and information defining that grid.
	 */
	PersistentSystemData(symphas::interval_data_type const& vdata, size_t id) :
		parent_type(vdata, id), writer{ data_len() } {}

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
		symphas::grid_info g{ info };
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
		symphas::grid_info g{ info };
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
		auto values = std::make_shared<T[]>(data_len());
		persist(values.get());

		return {
			std::shared_ptr<axis_nd_t<D>[]>(std::move(symphas::lib::new_system_axis_list<D>(info.intervals))),
			values,
			data_len() };
	}

	field_type as_field() const
	{
		return *this;
	}

protected:

	PersistentSystemData() : PersistentSystemData{ {}, 0 } {}

	WriteParallel<T> writer;
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
	using parent_type::data_len;


	Block<T> get_snapshot() const
	{
		Block<T> data(data_len());
		persist(data.values);
		return std::move(data);
	}

};


#endif





// *************************************************************************************************


namespace symphas
{

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

		//! Number of elements in the data arrays.
		/*!
		 * Number of elements in the data arrays which must match the number of
		 * systems in the problem these parameters represent.
		 */
		size_t len;

		problem_parameters_type() : len{ 0 }, tdata{ nullptr }, vdata{ nullptr }, bdata{ nullptr }, dt{ 1.0 } {}

	public:

		double dt;						//!< The time discretization.


		//! Initialize the problem parameters corresponding to \p len systems.
		/*!
		 * Initialize the problem parameters corresponding to \p len systems.
		 *
		 * \param len The number of systems that the problem parameters refer to.
		 */
		problem_parameters_type(const size_t len) : len{ len }, dt{ 1.0 },
			tdata{ new symphas::init_data_type[len] }, vdata{ new symphas::interval_data_type[len] }, bdata{ new symphas::b_data_type[len] } {}

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



		//! Set the initial data of the problem parameters.
		/*
		 * Set the initial data of the problem parameters, to be applied to one
		 * or more systems in a phase field problem initialization.
		 *
		 * \param tdata_set The list of initial data to apply. The i-th element
		 * in the given list applies to the i-th phase field to the given phase
		 * field problem.
		 * \param n The length of the list.
		 */
		void set_initial_data(symphas::init_data_type* tdata_set, size_t n);
		inline void set_initial_data(symphas::init_data_type* tdata_set)
		{
			set_initial_data(tdata_set, len - 1);
		}

		//! Set the interval data of the problem parameters.
		/*
		 * Set the interval data of the problem parameters, to be applied to one
		 * or more systems in a phase field problem initialization. The intervals
		 * must match the system dimension.
		 *
		 * \param vdata_set The list of interval data to apply. The i-th element
		 * in the given list applies to the i-th phase field to the given phase
		 * field problem.
		 * \param n The length of the list.
		 */
		void set_interval_data(symphas::interval_data_type* vdata_set, size_t n);
		inline void set_interval_data(symphas::interval_data_type* vdata_set)
		{
			set_interval_data(vdata_set, len - 1);
		}


		//! Set the boundary data of the problem parameters.
		/*
		 * Set the boundary data of the problem parameters, to be applied to one
		 * or more systems in a phase field problem initialization. The boundary
		 * dimension must match the system dimension.
		 *
		 * \param bdata_set The list of boundary data to apply. The i-th element
		 * in the given list applies to the i-th phase field to the given phase
		 * field problem.
		 * \param n The length of the list.
		 */
		void set_boundary_data(symphas::b_data_type* bdata_set, size_t n);
		inline void set_boundary_data(symphas::b_data_type* bdata_set)
		{
			set_boundary_data(bdata_set, len - 1);
		}


		//! Provide the time step used in the solution of the problem.
		inline void set_problem_time_step(double dt_set)
		{
			dt = dt_set;
		}


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

		//! Get the length of the data elements.
		inline const size_t length() const
		{
			return len;
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

		friend void swap(problem_parameters_type& first, problem_parameters_type& second);

		~problem_parameters_type();

	};


}

void swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second);

/*
 * system initialization helper functions
 */

namespace symphas::internal
{


	template<typename T, size_t D>
	void populate_tdata(symphas::init_data_type const& tdata, Grid<T, D>& data,
		[[maybe_unused]] symphas::grid_info* info, [[maybe_unused]] size_t id)
	{

		if (!InitialConditions<T, D>{ tdata, data.dims }.initialize(data.values, id))
		{
			fprintf(SYMPHAS_ERR, "the given initial condition algorithm is not valid\n");
		}
	}

}




