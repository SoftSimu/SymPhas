
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

#undef A
#undef B
#undef C
#undef D
#undef E
#undef F
#undef G
#undef H

#undef AA
#undef BB
#undef CC
#undef DD
#undef EE
#undef FF
#undef GG
#undef HH

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

#undef i
#undef lit
#undef one

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



