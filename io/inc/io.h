
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
 * PURPOSE: Defines the basic parts of the IO functionality.
 *
 * ***************************************************************************
 */

#pragma once

#include "gridinfo.h"

//! \cond
#ifdef IO_EXPORTS
#define DLLIO DLLEXPORT
#else
#define DLLIO DLLIMPORT
#endif
//! \endcond

/* the save directory names
 */

//! Relative name of the directory storing the data to be plotted.
#define DATA_DIR "data"

//! Relative name of the directory storing the plot configuration files.
#define PLOT_DIR "plot"

//! Relative name of the directory storing all the checkpoint information.
#define CHECKPOINT_DIR "checkpoint"

 /* the name of the backup configuration
  */

//! Key used at the end of the backup configuration to indicate a completed simulation.
#define SIMULATION_DONE_KEY "DONE"

//! Name of the backup configuration file, without extension.
#define BACKUP_CONFIG_NAME "configuration"

//! Relative location of the backup configuration, format string.
#define BACKUP_CONFIG_LOC_FMT "%s/" CHECKPOINT_DIR "/" BACKUP_CONFIG_NAME "." CONFIG_EXTENSION



/*
 * name formatting strings
 */

#define DATA_DIR_RELATIVE_PLOT ".."		//!< The directory of data relative to the plot configuration.
#define PHASEFIELD_DATA_NAME "data"	//!< Name used to identify the file of the phase field grid.
#define POSTFIX_ID_FMT "%zd"			//!< Format for printing the ID value of data.
#define OUTPUT_INDEX_WIDTH_STR "9"		//!< The maximum number of digits in the a solution index value.
#define OUTPUT_DATA_DIR "%s/" DATA_DIR	//!< The output data directory.

#define OUTPUT_PLOT_EXTENSION "gp"		//!< The extension given to plot files which are printed.
#define OUTPUT_DATA_EXTENSION "txt"		//!< The extension given to data output files.

 /* the format of the checkpoint datafile name, which is always the same regardless of the
  * type of printer or compilation options
  * the format is simply the defined data name, without any extension
  * if the data is saved in individual files per index, then the index is also prepended
  */
#define OUTPUT_CHECKPOINT_FMT "%s/" CHECKPOINT_DIR "/" PHASEFIELD_DATA_NAME POSTFIX_ID_FMT
#define OUTPUT_CHECKPOINT_INDEX_FMT "%s/" CHECKPOINT_DIR "/" PHASEFIELD_DATA_NAME POSTFIX_ID_FMT "_%d"



/* in writing the output, the data will be printed according to the specified
 * formatting definitions, related tot he number of decimals in the axis output
 * or data output
 * these definitions apply to the number of digits written by both the plotting
 * data outputs and the full grid outputs
 */

//! Number of decimals specified on the axis of the plot.
#define AXIS_OUTPUT_ACCURACY 2

//! Number of decimals specified in the phase field data output.
#define DATA_OUTPUT_ACCURACY 5

//! Stringified #AXIS_OUTPUT_ACCURACY
#define AXIS_OUTPUT_ACCURACY_STR STR(AXIS_OUTPUT_ACCURACY)

//! Stringified #DATA_OUTPUT_ACCURACY
#define DATA_OUTPUT_ACCURACY_STR STR(DATA_OUTPUT_ACCURACY)

//! The compression ratio for printing data with the XDR utility.
constexpr double XDR_COORD_COMPRESSION = 1E5;




//! Definition for the default stop index.
/*!
 * Definition associated with generalized save information
 */
#define DEFAULT_SAVE_STOP -1
#define DEFAULT_SAVE_BASE 0


//! Values representing the I/O types available.
/*!
 * The enumeration contains the types of I/O methods which may be
 * selected through the parameters. Specifically, see params::writer and
 * params::reader.
 */
enum class WriterType 
{ 
	GNU, 
	XDR,
	COLUMN,
	MOVIE,
	CSV
};


namespace symphas
{
	//! Contains all functionality related to I/O.
	/*!
	 * The functionality implemented by SymPhas to write persistent data and
	 * other information including plot configuration is a member of this namespace.
	 */
	namespace io {}
}

namespace params
{
	//! Only print plot configuration. 
	/*! 
	 * Don't do computation, only provide the plot files when printing output.
	 * 
	 * This is `false` by default.
	 */
	DLLIO extern bool plots_only;

	//! Picks the I/O writer.
	/*! 
	 * Controls which writer is chosen as the output.
	 * The data can be written as data to be read by gnuplot (and the associated
	 * gnuplot config files) or another writer based on the enum.
	 * 
	 * This is equal to WriterType::GNU by default.
	 */
	DLLIO extern WriterType writer;

	//! Picks the I/O reader.
	/*! 
	 * The reader can be chosen when loading data from
	 * another file instead of generating with given initial conditions. See
	 * params::writer.
	 * 
	 * This is equal to WriterType::GNU by default.
	 */
	DLLIO extern WriterType reader;

	//! Use the timestamp as a subdirectory in the output.
	/*!
	 * Allows the user to choose whether the current execution should be put
	 * under a directory with the timestamp name.
	 * 
	 * This is `true` by default.
	 */
	DLLIO extern bool use_timestamp;

	//! Indicates whether the output data should go to the same file.
	/*!
	 * The output is given in either separate files or appended to a single file.
	 * 
	 * This is `true` by default.
	 */
	DLLIO extern bool single_output_file;

	//! Indicates whether the input data comes from the same file.
	/*!
	 * The input is given in either separate files or appended to a single file.
	 *
	 * This is `true` by default.
	 */
	DLLIO extern bool single_input_file;

	//! Indicates whether a checkpoint should be made.
	/*!
	 * The checkpoint is restored from the given directory if this string is
	 * set; the backup configuration will be loaded and the original parameters
	 * of the program restored (this means that any parameters passed to the 
	 * current program will be ignored, and a message is given to indicate it).
	 * When the argument passed to checkpoint is parsed, params::start_index 
	 * will also be initialized.
	 * 
	 * This is expected in the format `checkpointdir`,`index`
	 * If the checkpoint option can't be tokenized again on ',' then no index 
	 * was given in which case the index will be taken from the backup 
	 * configuration file. In other words, if only `checkpointdir` is provided,
	 * then the latest index is used.
	 */
	DLLIO extern char* checkpoint;

	//! Defines the number of checkpoints to be made.
	/*!
	 * The num_checkpoints parameter holds the number of checkpoints that should be
	 * saved throughout the execution of the solution
	 */
	DLLIO extern int checkpoint_count;

}

//! Add the I/O parameters to the given map.
/*
 * Add the key value pairs for parsing and assigning io parameters.
 * 
 * \param param_map The map into which to add key/value pairs corresponding to
 * parameter string names and pointers to parameter variables.
 */
bool add_io_params(param_map_type& param_map);


/*
 *
 *
 *
 *
 *
 * datatypes pertaining to represent information about how data
 * should be saved
 */


//! The type of datafile determines how the file is saved.
/*!
 * The method of persisting the file changes based on what kind of data
 * is being saved.
 */
enum class DataFileType
{
	SOLUTION_DATA,		//!< Value representing phase field data which will be plotted.
	POSTPROC_DATA,		//!< Value representing some additional data not directly related to the phase field.
	CHECKPOINT_DATA,	//!< Value representing the backup of the phase field data.
	GENERIC_DATA,		//!< Some generic data which is saved, for example for testing.
	NAMED_DATA			//!< The written data is given a specific name, and is otherwise same functionally as SOLUTION_DATA.
};


//! Identifies how the spacing between save frames is chosen.
/*!
 * When saving data to file, different selections of when the save frame is
 * taken can be chosen. Typically, the interval between saves is always the
 * same, but for slowly evolving systems for example, it may be desirable
 * to increase the time between incrementally.
 */
enum class SaveType
{
	DEFAULT,			//!< The save intervals are equal.
	MULTIPLICATIVE,		//!< The save interval increases by a percentage.
	EXP,				//!< The save interval grows exponentially.
	LIST,				//!< The save intervals are given in a list.
	SPLIT,				//!< Make a given number of equally spaced intervals.
	NONE				//!< Does not represent any type of saving method.
};

struct SaveParams;

namespace symphas::io
{
	//! Get the save type from the key string.
	/*!
	 * Returns the save type corresponding to the given string key.
	 * If the key is not valid, then the SaveType::NONE value is returned.
	 * 
	 * \param str The key string to interpret.
	 */
	SaveType get_save_type(const char* str);

	//! Get the save key string based on the type. 
	/*!
	 * Returns the key string representing the given type. If the parameter
	 * is invalid, then a `nullptr` is returned.
	 * 
	 * \param type The save type to convert to string.
	 */
	const char* get_save_str(SaveType type);

	//! Parses the specification to initialize save information.
	/*!
	 * The given string is parsed to initialize save information. The parameters
	 * which are not identified in the specification are left unchanged, which is
	 * typically the default values.
	 * 
	 * The format can be given in two different ways:
	 * - `N`: The interval of saving will be every `N` iterations.
	 * - `N` `N0`: The interval of saving will be every `N` iterations starting
	 * with the `N0` frame.
	 * 
	 * Any word specified after the above formats will be interpreted as a
	 * modifier of the save routine. The modifier will change the interval
	 * of saving. There is only one modifier:
	 * - `EXP`: This will increase the spacing between save frame \f$i\f$ and
	 * \f$i+1\f$ according to N\f$^{i}\f$. I.e. if `N` is 10, then the space
	 * between the first save frame and the second save frame is 10, and then
	 * between the second and the third is 100, and so on.
	 * - `MULT`: The given value specifies the number of frames which are saved,
	 * and the intervals are chosen so that subsequent intervals grow as a
	 * percentage of the previous interval.
	 * 
	 * It is also possible to explicitly list the indices which are saved, this
	 * is done by providing the indices in a comma delimited list. This will
	 * automatically make the save type be value SaveType::LIST. _The comma
	 * delimited list must not have spaces_, and is not preceded or followed
	 * by any other tokens.
	 *
	 * \param str The save specification.
	 * \param save The save instance to modify based on the save specification.
	 */
	void parse_save_str(const char* str, SaveParams* save);


	//! Display the save parameters as a string.
	/*!
	 * The save parameters are displayed as a string, in a format which is
	 * compatible with symphas::io::parse_save_str().
	 * 
	 * \param save The save parameters which are written.
	 * \param output The string containing the result, which must be
	 * initialized already.
	 * \param buffer_size The maximum size of output.
	 */
	void save_as_str(SaveParams const* save, char* output, size_t buffer_size);



	//! Obtain the index of the next checkpoint.
	iter_type next_checkpoint(iter_type index, SaveParams const& save);

	//! Returns whether the index corresponds to a checkpoint.
	bool at_checkpoint(iter_type index, SaveParams const& save);

	//! Returns whether there is a checkpoint to load.
	/*!
	 * Indicates whether there is a checkpoint to load using the parameters.
	 */
	bool is_checkpoint_set();


	//! Get a list of save points defined by the given configuration.
	std::vector<iter_type> get_save_indices(SaveParams const& main_save, std::vector<SaveParams> const& other_saves = {});

	//! Get the index of the next global stop point.
	/*
	 * A stop point is anywhere that the solution should be interrupted to do
	 * another kind of activity, typically saving the phase field data.
	 *
	 * \param index The index to begin looking from, when finding the next index
	 * that the program should stop at.
	 */
	iter_type next_save_from_list(iter_type index, SaveParams const& main_save, std::vector<SaveParams> const& other_saves);

	//! Get the index of the next global stop point, given all the points.
	iter_type next_save_from_list(iter_type index, std::vector<iter_type> const& saves);

	//! Returns whether the given index is the last index at which to stop.
	/*!
	 * Checks whether the given index is the last index that the program will
	 * stop at. That is, no more solving will occur after this point according
	 * to the global configuration.
	 *
	 * \param index The index to compare to see if its the last index.
	 * \param saves List of savepoints which are compared.
	 */
	bool is_last_save(iter_type index, std::vector<iter_type> const& saves);

}

//! Information about how data should be saved.
/*!
 * Contains all the parameters related to how data should be saved, most
 * importantly the beginning and end of the iteration range, as well as
 * information about the save interval between frames.
 */
struct SaveParams
{

	SaveParams(SaveType type, double base, iter_type start, iter_type stop, bool init_flag) :
		init_flag{ init_flag }, type{ type }, start{ start }, stop{ stop }, base{ base } {}
	SaveParams(SaveType type, double base, iter_type start, iter_type stop) :
		SaveParams(type, base, start, stop, true) {}
	SaveParams(SaveType type, double base, iter_type stop) :
		SaveParams(type, base, INDEX_INIT, stop) {}
	SaveParams(SaveType type, double base) :
		SaveParams(type, base, DEFAULT_SAVE_STOP) {}

	SaveParams(double base, iter_type start, iter_type stop, bool init_flag) :
		SaveParams(SaveType::DEFAULT, base, start, stop, init_flag) {}
	SaveParams(double base, iter_type start, iter_type stop) :
		SaveParams(base, start, stop, true) {}
	SaveParams(double base, iter_type stop) :
		SaveParams(base, INDEX_INIT, stop) {}
	SaveParams(double base) :
		SaveParams(base, DEFAULT_SAVE_STOP) {}


	SaveParams(iter_type base, iter_type start, iter_type stop) :
		SaveParams(static_cast<double>(base), start, stop) {}
	SaveParams(iter_type base, iter_type stop) :
		SaveParams(static_cast<double>(base), stop) {}
	SaveParams(iter_type base) :
		SaveParams(static_cast<double>(base)) {}

	SaveParams() : SaveParams(DEFAULT_SAVE_BASE) {}

	SaveParams(iter_type* indices, size_t count) :
		init_flag{ true }, type{ SaveType::LIST },
		indices{ (count > 0) ? new iter_type[count] : nullptr }, count{ count }
	{
		std::copy(indices, indices + count, this->indices);
	}

	SaveParams(SaveParams const& other) : SaveParams()
	{
		init_flag = other.init_flag;
		type = other.type;

		if (other.type == SaveType::LIST)
		{
			count = other.count;
			indices = new iter_type[other.count];
			std::copy(other.indices, other.indices + other.count, indices);
		}
		else
		{
			start = other.start;
			stop = other.stop;
			base = other.base;
		}
	}

	SaveParams(SaveParams&& other) noexcept : SaveParams()
	{
		swap(*this, other);
	}

	SaveParams& operator=(SaveParams other)
	{
		swap(*this, other);
		return *this;
	}


	//! Change the size of the list of save indices.
	/*!
	 * Change the size of the list of save indices.
	 *
	 * \param count The new size of the save indices.
	 */
	void set_indices(size_t count)
	{
		if (type == SaveType::LIST)
		{
			delete[] this->indices;
		}
		else
		{
			type = SaveType::LIST;
		}
		this->count = count;
		this->indices = new iter_type[count];
	}

	//! Set a new list of save indices.
	/*!
	 * Set a new list of save indices.
	 *
	 * \param count The new size of the save indices.
	 */
	void set_indices(iter_type* indices, size_t count)
	{
		set_indices(count);
		std::copy(indices, indices + count, this->indices);
	}

	//! Change the parameters, depending on the type.
	/*!
	 * If the type is SaveType::LIST, then \p base will be interpreted as the number of
	 * save points to put between \p start and \p stop, where start will be
	 * omitted if SaveParams::init_flag is false. 
	 * 
	 * If the type is not SaveType::LIST, then the parameters will be 
	 * interpreted as according to the respective save algorithm. 
	 * 
	 * Therefore, set_init_flag()
	 * has to be called with the desired value before the parameters are set.
	 */
	void set_params(SaveType type, double base, iter_type start, iter_type stop);

	void set_params(SaveType type, double base, iter_type stop)
	{
		set_params(type, base, get_start(), stop);
	}

	void set_params(SaveType type, double base)
	{
		set_params(type, base, get_start(), get_stop());
	}

	void set_params(double base, iter_type start, iter_type stop);

	void set_params(double base, iter_type stop)
	{
		set_params(base, get_start(), stop);
	}
	
	void set_params(double base)
	{
		set_params(base, get_start(), get_stop());
	}

	//! Change the flag controlling if the very first index is saved.
	void set_init_flag(bool flag)
	{
		init_flag = flag;
	}

	//! Returns whether the very first index is saved.
	bool get_init_flag() const
	{
		return init_flag;
	}

	//! Return the first index that will be saved.
	iter_type get_start() const
	{
		if (type == SaveType::LIST)
		{
			return indices[0];
		}
		else
		{
			return start;
		}
	}

	//! Return the final index that will be saved.
	iter_type get_stop() const
	{
		if (type == SaveType::LIST)
		{
			return indices[count - 1];
		}
		else
		{
			return stop;
		}
	}

	//! Set the starting index for the save method.
	/*!
	 * Set the starting index for the save method.
	 */
	void set_start(iter_type start);

	//! Set the stop index for the save method.
	/*!
	 * Set the stop index for the save method.
	 */
	void set_stop(iter_type stop);

	//! Returns the number of saves that will be performed.
	iter_type num_saves() const;

	//! Get the index of the next global stop point.
	/*!
	 * Get the index of the next global stop point. The point after the
	 * given index is returned.
	 *
	 * \param current The index before the next save point.
	 */
	iter_type next_save(iter_type current) const;


	// the save state can be the following
	// 1. index 0 with saveinit = true (save the first index)
	// 2. index > 0 with index = save.start
	// 3. we want to save the current index but require index > 0 so no overlapping with 1.

	//! Returns true if the first index needs to be saved.
	/*!
	 * Returns true if the save object needs to be saved on the very first
	 * global save index, which corresponds to the value at params::start_index.
	 *
	 * \param index Return true if this is the first index.
	 */
	bool index_zero_save(iter_type index) const;

	//! Check whether the index is the very first save that is made. 
	/*!
	 * This test checks whether the state of the given save object will correspond
	 * to being saved for the first time.
	 *
	 * \param index Returns true if the given index is the first index to
	 * save.
	 */
	bool first_save(iter_type index) const;

	//! Check whether the current index corresponds to a save index.
	/*!
	 * Returns true if the given index corresponds to a save point for the save
	 * object.
	 *
	 * \param index Returns true if this value is the same as a save point.
	 */
	bool current_save(iter_type index) const;

	//! Returns whether the given index is the last index at which to stop.
	/*!
	 * Checks whether the given index is the last index that the program will
	 * stop at. That is, no more solving will occur after this point according
	 * to the global configuration.
	 *
	 * \param index The index to compare to see if its the last index.
	 * \param save Save parameter object which is compared.
	 */
	bool is_last_save(iter_type index) const;

	~SaveParams()
	{
		if (type == SaveType::LIST)
		{
			delete[] indices;
		}
	}

	friend void swap(SaveParams& first, SaveParams& second)
	{
		using std::swap;
		swap(first.init_flag, second.init_flag);
		swap(first.type, second.type);

		swap(first.start, second.start);
		swap(first.stop, second.stop);
		swap(first.base, second.base);
	}

	friend void symphas::io::parse_save_str(const char* value, SaveParams* save);
	friend void symphas::io::save_as_str(SaveParams const* save, char* output, size_t buffer_size);

protected:

	bool init_flag;				//!< Flag for whether to save the first time index.
	SaveType type;				//!< Changes how the interval is chosen.


	union
	{
		struct
		{
			iter_type
				start,			//!< The first save index.
				stop;			//!< The final index the solver will go to.
			double base;		//!< Used in calculating the each save index.
		};

		struct
		{
			iter_type* indices;	//!< List of indices to be saved at.
			size_t count;		//!< The number of indices in the list.
		};
	};

};

