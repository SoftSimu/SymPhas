
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
 * PURPOSE: Defines macros commonly used in the implementation of SymPhas.
 *
 * ***************************************************************************
 */



#pragma once


#undef COMMA
#undef EMPTY
#undef EXPAND
#undef UNUSED
#undef SINGLE_ARG
#undef STR_ARR_LEN

/*
 * macro expander for dealing with MSVC
 */


#define COMMA ,
#define EMPTY
#define EXPAND(x) x
#define UNUSED(...) do { (void)(__VA_ARGS__); } while (0);
#define STR_ARR_LEN(STR) (sizeof(STR) / sizeof(char))

#define __SYMPHAS_STR(...) #__VA_ARGS__
#define __SYMPHAS_CONCAT(A, B) A ## B


//! Simply outputs arguments as one whole, verbatim.
/*!
 * Definition workaround, typically used to treat an argument with commas as
 * a single argument when used as an input argument to another definition.
 */
#define SINGLE_ARG(...) __VA_ARGS__

//! Stringify the argument.
#define STR(...) __SYMPHAS_STR(__VA_ARGS__)

//! Concatenate the two arguments.
#define CONCAT(A, B) __SYMPHAS_CONCAT(A, B)

