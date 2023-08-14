
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

#ifdef PRINT_TIMINGS_INLINE

#define TIME_THIS_CONTEXT_ONCE(NAME) \
symphas::TimerReport NAME ## _timer(SYMPHAS_LOG, #NAME);
#define TIME_THIS_CONTEXT_LIFETIME(NAME) \
static symphas::TimerReport NAME ## _timer_report(SYMPHAS_LOG, #NAME); symphas::TimerContext NAME ## _timer(&NAME ## _timer_report);

#else

#define TIME_THIS_CONTEXT_ONCE(NAME, ...) 
#define TIME_THIS_CONTEXT_LIFETIME(NAME, ...) 

#endif


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
	struct Time
	{
		const char* desc;
		FILE* out;
		std::chrono::time_point<std::chrono::steady_clock> start_s;

		template<typename... Ts>
		Time(FILE* out, const char* str, Ts&&... args) : desc{ get_desc(str, std::forward<Ts>(args)...) }, out{ out }, start_s{ std::chrono::steady_clock::now() } {}
		template<typename... Ts>
		Time(const char* str, Ts&&... args) : Time(SYMPHAS_LOG, str, std::forward<Ts>(args)...) {}
		Time() : Time(nullptr) {}

		double current_duration()
		{
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_s);
			return duration.count() / 1000.0;
		}

		~Time()
		{
			if (desc != nullptr)
			{
				fprintf(out, "time for %-38s: %12.6lf s\n", desc, current_duration());
				delete[] desc;
			}
		}

	protected:

		static const int DESC_LEN = 1024;

		template<typename... Ts>
		const char* get_desc(const char* fmt, Ts&&... args)
		{
			if (fmt == nullptr)
			{
				return nullptr;
			}
			else
			{
				char* desc0 = new char[DESC_LEN] {};
				snprintf(desc0, DESC_LEN, fmt, std::forward<Ts>(args)...);
				return desc0;
			}
		}
	};


	struct TimerReport
	{
		const char* desc;
		FILE* out;
		std::chrono::nanoseconds count_s;
		int num;

		template<typename... Ts>
		TimerReport(FILE* out, const char* str, Ts&&... args) : desc{ get_desc(str, std::forward<Ts>(args)...) }, out{ out }, count_s{ 0 }, num{ 0 } {}

		template<typename... Ts>
		TimerReport(const char* str, Ts&&... args) : TimerReport(SYMPHAS_LOG, str, std::forward<Ts>(args)...) {}

		TimerReport() : TimerReport("<unspecified>") {}

		void operator+=(std::chrono::nanoseconds const& add_count)
		{
			count_s += add_count;
			++num;
		}

		double current_duration()
		{
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(count_s);
			return duration.count() / 1000.0;
		}

		~TimerReport()
		{
			if (desc != nullptr)
			{
				fprintf(out, "report '%-48s': %12.6lf s; N=%d\n", desc, current_duration(), num);
				delete[] desc;
			}
		}

	protected:

		static const int DESC_LEN = 1024;

		template<typename... Ts>
		const char* get_desc(const char* fmt, Ts&&... args)
		{
			if (fmt == nullptr)
			{
				return nullptr;
			}
			else
			{
				char* desc0 = new char[DESC_LEN] {};
				snprintf(desc0, DESC_LEN, fmt, std::forward<Ts>(args)...);
				return desc0;
			}
		}
	};


	struct TimerContext
	{
		TimerContext(TimerReport* report) : report{ report }, start_s{ std::chrono::steady_clock::now() } {}

		~TimerContext()
		{
			report->operator+=(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start_s));
		}

		TimerReport* report;
		std::chrono::time_point<std::chrono::steady_clock> start_s;
	};
}

