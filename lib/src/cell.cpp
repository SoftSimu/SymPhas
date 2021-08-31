
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

#include "cell.h"


template struct VectorValue<double, 1>;
template struct VectorValue<double, 2>;
template struct VectorValue<double, 3>;
template struct VectorValue<std::complex<double>, 1>;
template struct VectorValue<std::complex<double>, 2>;
template struct VectorValue<std::complex<double>, 3>;
template struct VectorValue<float, 1>;
template struct VectorValue<float, 2>;
template struct VectorValue<float, 3>;
template struct VectorValue<std::complex<float>, 1>;
template struct VectorValue<std::complex<float>, 2>;
template struct VectorValue<std::complex<float>, 3>;


template<>
inline double distance(VectorValue<double, 1> const& c1, VectorValue<double, 1> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0]));
}

template<>
inline double distance(VectorValue<double, 2> const& c1, VectorValue<double, 2> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0])
		+ (c1.v[1] - c2.v[1]) * (c1.v[1] - c2.v[1]));
}

template<>
inline double distance(VectorValue<double, 3> const& c1, VectorValue<double, 3> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0])
		+ (c1.v[1] - c2.v[1]) * (c1.v[1] - c2.v[1])
		+ (c1.v[2] - c2.v[2]) * (c1.v[2] - c2.v[2]));
}

template<>
inline double abs(VectorValue<double, 1> const& c)
{
	return sqrt(c.v[0] * c.v[0]);
}

template<>
inline double abs(VectorValue<double, 2> const& c)
{
	return sqrt(c.v[0] * c.v[0] + c.v[1] * c.v[1]);
}

template<>
inline double abs(VectorValue<double, 3> const& c)
{
	return sqrt(c.v[0] * c.v[0] + c.v[1] * c.v[1] + c.v[2] * c.v[2]);
}




template<>
inline float distance(VectorValue<float, 1> const& c1, VectorValue<float, 1> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0]));
}

template<>
inline float distance(VectorValue<float, 2> const& c1, VectorValue<float, 2> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0])
		+ (c1.v[1] - c2.v[1]) * (c1.v[1] - c2.v[1]));
}

template<>
inline float distance(VectorValue<float, 3> const& c1, VectorValue<float, 3> const& c2)
{
	return sqrt((c1.v[0] - c2.v[0]) * (c1.v[0] - c2.v[0])
		+ (c1.v[1] - c2.v[1]) * (c1.v[1] - c2.v[1])
		+ (c1.v[2] - c2.v[2]) * (c1.v[2] - c2.v[2]));
}

template<>
inline float abs(VectorValue<float, 1> const& c)
{
	return sqrt(c.v[0] * c.v[0]);
}

template<>
inline float abs(VectorValue<float, 2> const& c)
{
	return sqrt(c.v[0] * c.v[0] + c.v[1] * c.v[1]);
}

template<>
inline float abs(VectorValue<float, 3> const& c)
{
	return sqrt(c.v[0] * c.v[0] + c.v[1] * c.v[1] + c.v[2] * c.v[2]);
}




