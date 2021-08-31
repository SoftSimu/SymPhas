
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
 * MODULE:  conf
 * PURPOSE: Defines full configuration functionality.
 *
 * ***************************************************************************
 */

#pragma once


#ifdef USING_PROC
#include "confdata.h"
#else
#include "confsystem.h"
#endif

struct Conf;


namespace symphas::conf
{

	//! Returns a new Conf instance from the given file.
	/*!
	 * The given file is parsed so that each of the entries are put in a
	 * key value pair and appended into a list, which is used to initialize
	 * the configuration. The title is also parsed from the file. Any parameters
	 * that are also specified in the file are then parsed and overwrite the
	 *
	 * \param file The name of the file that configuration properties are read
	 * from.
	 * \param param_map The parameter keys which may be set by this configuration.
	 */
	Conf make_config(const char* file, param_map_type const& param_map);

	//! Returns a new Conf instance from the given file.
	/*!
	 * See symphas::conf::make_config(const char*, param_map_type const&). This
	 * overload does not initialize any parameters that are part of the
	 * configuration file.
	 */
	Conf make_config(const char* file);


	//! Access the global configuration.
	/*!
	 * The global configuration first needs to be initialized using the function
	 * setup_global_config(const char*).
	 */
	const Conf& config();

	//! Initializes the globally accessible configuration.
	/*!
	 * The configuration properties accessible through the whole program is
	 * created, and then accessed using the function config().
	 *
	 * \param file The name of the file that configuration properties are read
	 * from.
	 * \param param_map The key-value list of parameters that can be used to
	 * parse arguments in the configuration.
	 */
	void setup_global_config(const char* file, param_map_type param_map = {});

	//! Initializes the globally accessible configuration.
	/*!
	 * The configuration properties accessible through the whole program is
	 * created, and then accessed using the function config(). The global
	 * configuration is initialized by copy from the given configuration.
	 *
	 * \param configuration The configuration that is used as the global
	 * configuration.
	 */
	void setup_global_config(Conf const& configuration);

	//! Get a list of save points defined by the given configuration.
	/*!
	 * A list of indices is constructed using the given configuration, which
	 * represent the points at which the solution should be paused in order
	 * to persist solution data to disk.
	 *
	 * \param config The configuration which defines the save parameters.
	 */
	std::vector<iter_type> get_save_indices(Conf const& config = symphas::conf::config());

	//! Get the index of the next global stop point.
	/*
	 * A stop point is anywhere that the solution should be interrupted to do
	 * another kind of activity, typically saving the phase field data.
	 *
	 * \param index The index to begin looking from, when finding the next index
	 * that the program should stop at.
	 */
	iter_type next_save(iter_type index, Conf const& config = symphas::conf::config());

	//! Returns whether the given index is the last index at which to stop.
	/*!
	 * Checks whether the given index is the last index that the program will
	 * stop at. That is, no more solving will occur after this point according
	 * to the global configuration.
	 *
	 * \param index The index to compare to see if its the last index.
	 * \param config The configuration that the save parameters are taken from.
	 */
	bool is_last_save(iter_type index, Conf const& config = symphas::conf::config());

	//! Returns the number of saves that will be performed.
	iter_type num_saves(Conf const& config = symphas::conf::config());

	//! Restore the checkpoint from the given file.
	/*!
	 * \param param_map The list of parameter keys to use when reading the backup
	 * checkpoint.
	 * \param dir The directory from which the checkpoint is retrieved. This
	 * directory must have the folder "checkpoint".
	 * \param index The solution index retrieved from the checkpoint.
	 */
	Conf restore_checkpoint(param_map_type param_map, const char* dir, int index = -1);


	//! See restore_checkpoint().
	/*!
	 * Uses information given by the variable params::checkpoint.
	 *
	 * \param param_map The list of parameter keys to use when reading the backup
	 * checkpoint.
	 */
	Conf restore_checkpoint(param_map_type param_map);
}



//! Specifies and manages all configurable options of SymPhas.
/*!
 * Instantiated using a list of string pairs, where each pair is a key-value
 * pair. The keys do not necessarily have to match any configuration property.
 * 
 * The title of the configuration is separately passed, and this typically
 * corresponds to the name of the phase field problem being configured.
 */
struct Conf :
	SystemConf
#ifdef USING_PROC
	, DataConf
#endif
{
	Conf(std::vector<std::pair<std::string, std::string>> options, const char* title, const char* dir = "");
	Conf(const char* file);
	Conf(symphas::problem_parameters_type const& parameters, const char* title, const char* dir = "");

	//! Set the simulation for this configuration as complete.
	void set_done();

	//! Check if the configuration is complete.
	bool is_done() const;

	//! Output the configuration to the provided directory.
	void write(const char* savedir) const;

protected:

	//! Contains the indices for which a solution has already been computed.
	/*!
	 * The solution indices typically record the phase field data for which
	 * there is a backup, by writing the index prepended by 
	 * #CONFIG_INDEX_PREFIX.
	 */
	std::vector<iter_type> computed_indices;

	//! If true, the solution from this configuration has been completed.
	bool sim_done;


	//! Add computed indices to the list.
	/*!
	 * Add computed indices to the list.
	 * The list is always maintained sorted.
	 */
	void append_computed_index(std::vector<iter_type> indices);

	//! Add a computed index to the list.
	/*!
	 * Add a computed index to the list.
	 * The list is always maintained sorted.
	 */
	void append_computed_index(iter_type index);

	//! Return the precomputed indices.
	auto get_computed_indices() const;

	//! Return the last index in the list.
	/*!
	 * Since the list is sorted, this is always the
	 * largest index.
	 */
	iter_type get_last_computed_index() const;

	friend Conf symphas::conf::make_config(const char* file, param_map_type const& param_map);
	friend Conf symphas::conf::restore_checkpoint(param_map_type param_map, const char* dir, int index);
};

