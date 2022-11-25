
/* ***************************************************************************
 * This file is part of the SymPhas package, containing a framework for
 * implementing solvers for phase-field problems with compile-time symbolic
 * algebra.
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
 * This file includes all of the phase-field crystal model examples.
 *
 * ***************************************************************************
 */

#pragma once

#include "modelmacros.h"

// the default parameters (if not specified) for any phase field model are:
// NONCONSERVED
// MODE(1)


PFC_TYPE(PC, DEFAULTS(
	DYNAMICS(CONSERVED)),
	(SCALAR)
)
LINK_PFC_WITH_NAME(PC, PFC_C)

PFC_TYPE(PNC, NO_DEFAULTS, (SCALAR))
LINK_PFC_WITH_NAME(PNC, PFC_NC)

/*
PFC_TYPE(PC2, DEFAULTS(
	DYNAMICS(CONSERVED)),
	(SCALAR, SCALAR)
)
LINK_PFC_WITH_NAME(PC2, PFC_C2)

PFC_TYPE(PNC2, DEFAULTS(
	DYNAMICS(NONCONSERVED)),
	(SCALAR, SCALAR)
)
LINK_PFC_WITH_NAME(PNC2, PFC_NC2)

PFC_TYPE(PCNC, DEFAULTS(
	DYNAMICS(NONCONSERVED)),
	(SCALAR, SCALAR),
	DYNAMICS_OF(1, CONSERVED)
)
LINK_PFC_WITH_NAME(PCNC, PFC_CNC)
*/







