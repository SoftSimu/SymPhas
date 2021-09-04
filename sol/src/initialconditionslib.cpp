
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
 */


#include "initialconditionslib.h"


Inside symphas::in_from_str(const char* type)
{
	auto condition = internal::init_key_map.find(type);
	if (condition != internal::init_key_map.end())
	{
		return condition->second;
	}
	return Inside::NONE;
}

const char* symphas::str_from_in(Inside in)
{
	for (auto& [k, v] : internal::init_key_map)
	{
		if (v == in)
		{
			return k;
		}
	}
	return nullptr;
}

InsideTag symphas::in_tag_from_str(const char* type)
{
	auto tag = internal::init_tag_key_map.find(type);
	if (std::strlen(type) == 0)
	{
		return InsideTag::DEFAULT;
	}

	if (tag != internal::init_tag_key_map.end())
	{
		return tag->second;
	}
	return InsideTag::NONE;
}

const char* symphas::str_from_in_tag(InsideTag tag)
{
	for (auto& [k, v] : internal::init_tag_key_map)
	{
		if (v == tag)
		{
			return k;
		}
	}
	return nullptr;
}

std::map<const char*, Inside, symphas::internal::any_case_comparator> symphas::internal::init_key_map = {
	{"GAUSSIAN", Inside::GAUSSIAN},
	{"GA", Inside::GAUSSIAN},
	{"UNIFORM", Inside::UNIFORM},
	{"UN", Inside::UNIFORM},
	{"CAPPED", Inside::CAPPED},
	{"CA", Inside::CAPPED},
	{"CONSTANT", Inside::CONSTANT},
	{"CO", Inside::CONSTANT},
	{"CIRCLE", Inside::CIRCLE},
	{"CI", Inside::CIRCLE},
	{"SQUARE", Inside::SQUARE},
	{"SQ", Inside::SQUARE},
	{"HEXAGONAL", Inside::HEXAGONAL},
	{"HX", Inside::HEXAGONAL},
	{"CUBIC", Inside::CUBIC},
	{"CU", Inside::CUBIC},
	{"SEEDSSQUARE", Inside::SEEDSSQUARE},
	{"SQUARESEEDS", Inside::SEEDSSQUARE},
	{"SS", Inside::SEEDSSQUARE},
	{"SEEDSCIRCLE", Inside::SEEDSCIRCLE},
	{"CIRCLESEEDS", Inside::SEEDSCIRCLE},
	{"CS", Inside::SEEDSCIRCLE},
	{"VORONOI", Inside::VORONOI},
	{"VO", Inside::VORONOI},
	{"EXPRESSION", Inside::EXPRESSION}
};



std::map<const char*, InsideTag, symphas::internal::any_case_comparator> symphas::internal::init_tag_key_map = {
	{"DEFAULT", InsideTag::DEFAULT},
	{"RANDOM", InsideTag::RANDOM},
	{"RND", InsideTag::RANDOM},
	{"INVERT", InsideTag::INVERT},
	{"INV", InsideTag::INVERT},
	{"VARA", InsideTag::VARA},
	{"A", InsideTag::VARA},
	{"VARB", InsideTag::VARB},
	{"B", InsideTag::VARB},
	{"VARC", InsideTag::VARC},
	{"C", InsideTag::VARC}
};

