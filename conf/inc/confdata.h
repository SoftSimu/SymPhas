
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
 * PURPOSE: Defines the configuration functionality for including
 * postprocessing jobs.
 *
 * ***************************************************************************
 */

#pragma once

#include "data.h"
#include "confsystem.h"



struct DataConf
{

	DataConf(std::vector<std::pair<std::string, std::string>> = {});

	//! Parses the configuration that specifies a post-processing job.
	/*!
	 * Initializes the post-processing job based on the configuration in the 
	 * form of a string. The string is parsed the appropriate configuration
	 * is selected, which also includes the save interval.
	 * 
	 * \param v The configuration of the post-processing job.
	 * \param d The post-processing job parameters that are initialized.
	 * \param default_proc The default value of the type of the post-processing
	 * job that is performed.
	 */
	void parse_proc_job(const char* v, DataParams* d, ProcessType default_proc);



	void write(const char* savedir, const char* name = BACKUP_CONFIG_NAME) const;
};


