
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
 * MODULE:  lib
 * PURPOSE: A timer class to measure and print elapsed time of code
 * segments. It is clock time and not CPU time.
 *
 * ***************************************************************************
 */


#pragma once

#include <chrono>
#include <cstring>

#include "definitions.h"

namespace symphas
{
	//! Measures the time in a given block of code.
	/*!
	 * A timer will be initialized on object creation and be stopped on object
	 * destruction. Subsequently, the time length will be printed to the desired
	 * stream, or the log (#SYMPHAS_LOG) if no stream is provided.
	 * 
	 * This measures the clock time, not the CPU time.
	 *
	 * Format arguments may be provided in order to customize the timer final
	 * message.
	 * 
	 * \tparam Ts... The types of the format arguments.
	 */
	template<typename... Ts>
	struct Time
	{
		std::chrono::time_point<std::chrono::steady_clock> start_s;
		const char* desc;
		FILE* out;

		Time(FILE* out, const char* str, Ts... args) : start_s{ std::chrono::steady_clock::now() }, desc{ get_desc(str, args...) }, out{ out } {}
		Time(const char* str, Ts... args) : Time(SYMPHAS_LOG, str, args...) {}

		double current_duration()
		{
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_s);
			return duration.count() / 1000.0;
		}

		~Time()
		{
			fprintf(out, "time for %-35s: %lf s\n", desc, current_duration());
			delete[] desc;
		}

	protected:

		static const int DESC_LEN = 1024;

		const char* get_desc(const char* fmt, Ts... args)
		{
			char* desc0 = new char[DESC_LEN] {};
			snprintf(desc0, DESC_LEN, fmt, args...);
			return desc0;
		}
	};

	//! Measures the time in a given block of code.
	/*!
	 * See ::Time. No parameters are provided for formatting, so no output
	 * is generated upon object destruction.
	 */
	template<>
	struct Time<void>
	{
		std::chrono::time_point<std::chrono::steady_clock> start_s;

		Time() : start_s{ std::chrono::steady_clock::now() } {}

		double current_duration()
		{
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_s);
			return duration.count() / 1000.0;
		}
	};

	template<typename... Ts>
	Time(FILE*, const char*, Ts...)->Time<Ts...>;

	template<typename... Ts>
	Time(const char*, Ts...)->Time<Ts...>;

	Time()->Time<void>;
}

