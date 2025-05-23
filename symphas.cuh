
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
 * This file is part of the SymPhas API. This is the main header file that is
 * configured by the CMake build process.
 *
 * ***************************************************************************
 */

#pragma once

#include "definitions.h"

#ifdef USING_CONF
#include "modules-conf.h"
#elif defined(USING_IO)
#include "modules-io.h"
#else
#include "modules-none.h"
#endif

// import sol and sym headers
#include "modeldefs.h"

// cuda implementations
#include "expressioncuda.cuh"
#include "modelarray.cuh"

// import implementation files
#include "modules-models.h"

// model and solver definitions
// #include "prereq-defs.h"

#undef BUFFER_LENGTH_R4
#undef BUFFER_LENGTH_R2
#undef BUFFER_LENGTH_R1
#undef BUFFER_LENGTH
#undef BUFFER_LENGTH_L2
#undef BUFFER_LENGTH_L3
#undef BUFFER_LENGTH_L4
#undef BUFFER_LENGTH_L5
#undef BUFFER_LENGTH_R
#undef BUFFER_LENGTH_L
#undef LINE_READ_BUFFER
#undef FILECHUNK_READ_BUFFER

#undef ERR_CODE_FILE_OPEN
#undef ERR_MSG_FILE_OPEN
