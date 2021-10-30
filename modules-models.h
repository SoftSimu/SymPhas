
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
 * This file is part of the SymPhas API. It is used when model definitions
 * are included from a separate file.
 *
 * ***************************************************************************
 */

// header used when model headers are included

#pragma once

#ifndef MODEL_INCLUDE_HEADER

#include "modelmacros.h"

#else
#include MODEL_INCLUDE_HEADER

#ifdef USING_MODEL_SELECTION

struct model_select
{
	model_select(size_t dimension, StencilParams stp) :
		dimension{ dimension }, stp{ stp } {}

	size_t dimension;
	StencilParams stp;

	/* returns false if there is no model with the given name
	 */
	template<template<typename> typename AppliedSolver, typename... Ts>
	auto call(const char* name, Ts&& ...args)
	{
		constexpr int last_index = decltype(model_counter(model_count_index<255>{}))::value;
		return model_call_wrapper<last_index - 1>::call<AppliedSolver>(dimension, stp, name, std::forward<Ts>(args)...);
	}

};


// undefine all the parameters

#undef c1 
#undef c2 
#undef c3 
#undef c4 
#undef c5 
#undef c6 
#undef c7 
#undef c8 
#undef c9 
#undef c10
#undef c11
#undef c12
#undef c13
#undef c14
#undef c15
#undef c16

// order parameter names

#undef op
#undef dop
#undef var

#undef dx
#undef lap
#undef bilap
#undef gradlap
#undef grad
#undef param

#undef Ii
#undef lit
#undef one


#undef x
#undef y
#undef z
#undef t

// functions

#undef modulus
#undef conj
#undef Im
#undef Re
#undef Gaussian
#undef smoothing
#undef cross


#undef VECTOR
#undef SCALAR
#undef COMPLEX

#undef DYNAMIC
#undef MODE


#endif

#endif


template<size_t D>
struct init_expr_select
{
	/* returns false if there is no model with the given name
	 */
	template<typename T>
	auto call(const char* name, T* values, len_type const* dims, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
	{
        using namespace symphas::internal;
		constexpr int last_index = decltype(init_expr_counter(init_expr_count_index<255>{}))::value;
		return init_expr_call_wrapper<D, last_index - 1>::template call(name, values, dims, vdata, coeff, num_coeff);
	}
};

template<size_t D, typename T>
bool match_init_expr(const char* initname, T* values, len_type const* dims, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
{
	init_expr_select<D> ie;

	if (ie.call(initname, values, dims, vdata, coeff, num_coeff) < 0)
	{
		fprintf(SYMPHAS_ERR, "unknown initialization equation provided, '%s'\n", initname);
		return false;
	}
	return true;
}





