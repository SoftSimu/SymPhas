
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


#include "initialconditions.h"



inline bool inquad(size_t n, size_t(&dim)[2])
{
	size_t
		& L = dim[0],
		& M = dim[1];

	return
		(n % L < L / 2) &&
		(n / L < M / 2);
}

inline bool inquad(size_t n, size_t(&dim)[3])
{
	size_t
		& L = dim[0],
		& M = dim[1],
		& N = dim[2];

	return
		n % L < L / 2 &&
		(n / L) % M < M / 2 &&
		n / (L * M) < N / 2;
}


// ********************************************************************************


// ********************************************************************************

template<>
scalar_t InitialConditionsAlg<1, Inside::SQUARE>::operator[](iter_type n) const
{
	double ra = dims[0] / 2.0;
	double aa = ra / init.data.gp[0];
	double x = static_cast<double>(n);

	return (x > ra - aa && x < ra + aa) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double shift_circle_x = offsets.get_delta(0);

	double ra = dims[0] / 2.0;
	double aa = ra / init.data.gp[0] - shift_circle_x;
	double x = static_cast<double>(n);

	return (x > ra - aa && x < ra + aa) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<1, Inside::CUBIC>::operator[](iter_type n) const
{
	double const
		dx = dims[0] / init.data.gp[0],
		aa = dx / (4.0 * init.data.gp[1]);

	iter_type const ci = static_cast<iter_type>(n / dx);
	double const mx = n - dx * ci;

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx },
		p{ 0 }
	};

	for (auto const& u : pts)
	{
		if ((mx - u.x) * (mx - u.x) / (aa * aa) < 1.0)
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<1, Inside::CUBIC, InsideTag::VARA>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		aa = dx / (4.0 * init.data.gp[1]);

	iter_type const
		ci = static_cast<iter_type>(n / dx),
		mx = n - static_cast<iter_type>(dx * ci);

	double const rx = dx / 2;

	if (symphas::internal::is_in_circle_1(mx, rx, aa))
	{
		return IC_INNER_VALUE;
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<1, Inside::CUBIC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double const
		dx = dims[0] / init.data.gp[0],
		aa = dx / (4.0 * init.data.gp[1]);

	iter_type const ci = static_cast<iter_type>(n / dx);
	double const mx = n - dx * ci;

	iter_type const
		i_l = static_cast<iter_type>((ci + 1 == static_cast<iter_type>(init.data.gp[0])) ? 0 : ci + 1) + ci,
		i_r = static_cast<iter_type>((ci + 1 == static_cast<iter_type>(init.data.gp[0])) ? 0 : ci + 1);

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx + offsets.get_delta(i_r) },
		p{ 0 + offsets.get_delta(i_l) }
	};

	for (auto const& u : pts)
	{
		if ((mx - u.x) * (mx - u.x) / (aa * aa) < 1.0)
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}



template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets, &values,
		init.data.gp[2], init.data.gp[3], 0);
}




template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}



template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_square_1, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_1(
		n, dims, init, 
		&symphas::internal::is_in_circle_1, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_1(
		n, dims, init,
		&symphas::internal::is_in_circle_1,
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_circle_1,
		&offsets, &values,
		seed_value, field_value, rnd_offset);
}



template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_1(
		n, dims, init,
		&symphas::internal::is_in_circle_1,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_1(
		n, dims, init,
		&symphas::internal::is_in_circle_1,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_1(
		n, dims, init,
		&symphas::internal::is_in_square_1,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}






// ********************************************************************************

template<>
scalar_t InitialConditionsAlg<2, Inside::SQUARE>::operator[](iter_type n) const
{
	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0;

	double const
		aa = (dims[0]) / (2.0 * init.data.gp[0]),
		bb = (dims[1]) / (2.0 * init.data.gp[1]);

	iter_type const
		x = n % dims[0],
		y = n / dims[0];

	return symphas::internal::is_in_square_2(x, ux, aa, y, uy, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	iter_type
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]);

	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0;

	double const
		aa = (dims[0]) / (2.0 * init.data.gp[0]),
		bb = (dims[1]) / (2.0 * init.data.gp[1]);

	iter_type const
		x = n % dims[0] + shift_x,
		y = n / dims[0] + shift_y;

	return symphas::internal::is_in_square_2(x, ux, aa, y, uy, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<2, Inside::SQUARE, InsideTag::VARA>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 2);
	double u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1];

	iter_type const
		x = n % dims[0],
		y = n / dims[0];

	return symphas::internal::is_in_square_2(x, u, aa, y, u, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 2);
	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]);

	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1];

	iter_type const
		x = n % dims[0] + shift_x,
		y = n / dims[0] + shift_y;

	return symphas::internal::is_in_square_2(x, u, aa, y, u, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<2, Inside::CIRCLE>::operator[](iter_type n) const
{
	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0;

	double const
		aa = (dims[0]) / (2.0 * init.data.gp[0]),
		bb = (dims[1]) / (2.0 * init.data.gp[1]);

	iter_type const
		x = n % dims[0],
		y = n / dims[0];

	return symphas::internal::is_in_circle_2(x, ux, aa, y, uy, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CIRCLE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]);

	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0;

	double const
		aa = (dims[0]) / (2.0 * init.data.gp[0]),
		bb = (dims[1]) / (2.0 * init.data.gp[1]);

	iter_type const
		x = n % dims[0] + shift_x,
		y = n / dims[0] + shift_y;

	return symphas::internal::is_in_circle_2(x, ux, aa, y, uy, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CIRCLE, InsideTag::VARA>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 2);
	double u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1];

	iter_type const
		x = n % dims[0],
		y = n / dims[0];

	return symphas::internal::is_in_circle_2(x, u, aa, y, u, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 2);
	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]);

	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1];

	iter_type const
		x = n % dims[0] + shift_x,
		y = n / dims[0] + shift_y;

	return symphas::internal::is_in_circle_2(x, u, aa, y, u, bb) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CUBIC>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = n / dims[0],
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj);

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx, dy },
		p{ dx, 0 },
		p{ 0, 0 },
		p{ 0, dy }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_2(mx, u.x, aa, my, u.y, bb))
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CUBIC, InsideTag::VARA>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = n / dims[0],
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj);

	double const
		rx = dx / 2,
		ry = dy / 2;

	if (symphas::internal::is_in_circle_2(mx, rx, aa, my, ry, bb))
	{
		return IC_INNER_VALUE;
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::CUBIC, InsideTag::RANDOM>::operator[](iter_type n) const
{

	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = n / dims[0],
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj);

	/* these variables are to map each position to the correct random index
	 */

	iter_type const
		ci_len = static_cast<iter_type>(init.data.gp[0]),
		cj_len = static_cast<iter_type>(init.data.gp[1]),
		ci_next = (ci + 1) % ci_len,
		cj_next = (cj + 1) % cj_len;

	iter_type const
		i_top_l = cj_next * ci_len + ci,
		i_top_r = cj_next * ci_len + ci_next,
		i_bot_l = cj * ci_len + ci,
		i_bot_r = cj * ci_len + ci_next;

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx + offsets.get_delta(i_top_r)[0], dy + offsets.get_delta(i_top_r)[1] },
		p{ dx + offsets.get_delta(i_bot_r)[0], 0 + offsets.get_delta(i_bot_r)[1] },
		p{ 0 + offsets.get_delta(i_bot_l)[0], 0 + offsets.get_delta(i_bot_l)[1] },
		p{ 0 + offsets.get_delta(i_top_l)[0], dy + offsets.get_delta(i_top_l)[1] }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_2(mx, u.x, aa, my, u.y, bb))
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::HEXAGONAL>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = n / dims[0],
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj);

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx / 2, 0 },
		p{ dx / 2, dy },
		p{ 0, dy / 2 },
		p{ dx, dy / 2 }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_2(mx, u.x, aa, my, u.y, bb))
		{
			return IC_INNER_VALUE;
		}
	}

	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<2, Inside::HEXAGONAL, InsideTag::RANDOM>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = n / dims[0],
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj);

	/* these variables are to map each position to the correct random index
	 */
	iter_type const
		ci_len = static_cast<iter_type>(init.data.gp[0]),
		cj_len = static_cast<iter_type>(init.data.gp[1]),
		ci_next = (ci + 1) % ci_len,
		cj_next = (cj + 1) % cj_len;

	iter_type const
		i_bot = cj * ci_len * 2 + ci,
		i_top = cj_next * ci_len * 2 + ci,
		i_left = cj * ci_len * 2 + ci_len + ci,
		i_right = cj * ci_len * 2 + ci_len + ci_next;

	// position of center of each ellipse
	using p = struct { double x, y; };
	p pts[] = {
		p{ dx / 2 + offsets.get_delta(i_bot)[0], 0 + offsets.get_delta(i_bot)[1] },
		p{ dx / 2 + offsets.get_delta(i_top)[0], dy + offsets.get_delta(i_top)[1] },
		p{ 0 + offsets.get_delta(i_left)[0], dy / 2 + offsets.get_delta(i_left)[1] },
		p{ dx + offsets.get_delta(i_right)[0], dy / 2 + offsets.get_delta(i_right)[1] }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_2(mx, u.x, aa, my, u.y, bb))
		{
			return IC_INNER_VALUE;
		}
	}

	return IC_OUTER_VALUE;
}




template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets, &values,
		init.data.gp[2], init.data.gp[3], 0);
}




template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}



template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets, &lengths,
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_square_2, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_2(
		n, dims, init, 
		&symphas::internal::is_in_circle_2, 
		&offsets, &lengths, &values,
		init.data.gp[2], init.data.gp[3], 0);
}


template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_2(
		n, dims, init,
		&symphas::internal::is_in_circle_2,
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_circle_2,
		&offsets, &values,
		seed_value, field_value, rnd_offset);
}



template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_2(
		n, dims, init,
		&symphas::internal::is_in_circle_2,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}



template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_2(
		n, dims, init,
		&symphas::internal::is_in_circle_2,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_2(
		n, dims, init,
		&symphas::internal::is_in_square_2,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}







// ********************************************************************************

template<>
scalar_t InitialConditionsAlg<3, Inside::SQUARE>::operator[](iter_type n) const
{
	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0,
		uz = dims[2] / 2.0;

	double const
		aa = ux / init.data.gp[0],
		bb = uy / init.data.gp[1],
		cc = uz / init.data.gp[2];

	iter_type const
		x = n % dims[0],
		y = (n / dims[0]) % dims[1],
		z = n / (dims[0] * dims[1]);

	return symphas::internal::is_in_square_3(x, ux, aa, y, uy, bb, z, uz, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	 iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]),
		shift_z = static_cast<iter_type>(offsets.get_delta(0)[2]);

	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0,
		uz = dims[2] / 2.0;

	double const
		aa = ux / init.data.gp[0],
		bb = uy / init.data.gp[1],
		cc = uz / init.data.gp[2];

	iter_type const
		x = n % dims[0] + shift_x,
		y = (n / dims[0]) % dims[1] + shift_y,
		z = n / (dims[0] * dims[1]) + shift_z;

	return symphas::internal::is_in_square_3(x, ux, aa, y, uy, bb, z, uz, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<3, Inside::SQUARE, InsideTag::VARA>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 3);
	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1],
		cc = u / init.data.gp[2];

	iter_type const
		x = n % dims[0],
		y = (n / dims[0]) % dims[1],
		z = n / (dims[0] * dims[1]);

	return symphas::internal::is_in_square_3(x, u, aa, y, u, bb, z, u, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<3, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 3);

	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]),
		shift_z = static_cast<iter_type>(offsets.get_delta(0)[2]);

	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1],
		cc = u / init.data.gp[2];

	iter_type const
		x = n % dims[0] + shift_x,
		y = (n / dims[0]) % dims[1] + shift_y,
		z = n / (dims[0] * dims[1]) + shift_z;

	return symphas::internal::is_in_square_3(x, u, aa, y, u, bb, z, u, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<3, Inside::CIRCLE>::operator[](iter_type n) const
{
	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0,
		uz = dims[2] / 2.0;

	double const
		aa = ux / init.data.gp[0],
		bb = uy / init.data.gp[1],
		cc = uz / init.data.gp[2];

	iter_type const
		x = n % dims[0],
		y = (n / dims[0]) % dims[1],
		z = n / (dims[0] * dims[1]);

	return symphas::internal::is_in_circle_3(x, ux, aa, y, uy, bb, z, uz, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::CIRCLE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]),
		shift_z = static_cast<iter_type>(offsets.get_delta(0)[2]);

	double const
		ux = dims[0] / 2.0,
		uy = dims[1] / 2.0,
		uz = dims[2] / 2.0;

	double const
		aa = ux / init.data.gp[0],
		bb = uy / init.data.gp[1],
		cc = uz / init.data.gp[2];

	iter_type const
		x = n % dims[0] + shift_x,
		y = (n / dims[0]) % dims[1] + shift_y,
		z = n / (dims[0] * dims[1]) + shift_z;

	return symphas::internal::is_in_circle_3(x, ux, aa, y, uy, bb, z, uz, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::CIRCLE, InsideTag::VARA>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 3);
	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1],
		cc = u / init.data.gp[2];

	iter_type const
		x = n % dims[0],
		y = (n / dims[0]) % dims[1],
		z = n / (dims[0] * dims[1]);

	return symphas::internal::is_in_circle_3(x, u, aa, y, u, bb, z, u, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	auto d = *std::min_element(dims, dims + 3);

	iter_type const
		shift_x = static_cast<iter_type>(offsets.get_delta(0)[0]),
		shift_y = static_cast<iter_type>(offsets.get_delta(0)[1]),
		shift_z = static_cast<iter_type>(offsets.get_delta(0)[2]);

	double const u = d / 2.0;

	double const
		aa = u / init.data.gp[0],
		bb = u / init.data.gp[1],
		cc = u / init.data.gp[2];

	iter_type const
		x = n % dims[0] + shift_x,
		y = (n / dims[0]) % dims[1] + shift_y,
		z = n / (dims[0] * dims[1]) + shift_z;

	return symphas::internal::is_in_circle_3(x, u, aa, y, u, bb, z, u, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}


template<>
scalar_t InitialConditionsAlg<3, Inside::CUBIC>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1],
		dz = dims[2] / init.data.gp[2];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[3]),
		bb = dy / (4.0 * init.data.gp[3]),
		cc = dz / (4.0 * init.data.gp[3]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = (n / dims[0]) % dims[1],
		cz = n / (dims[0] * dims[1]),
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy),
		ck = static_cast<iter_type>(cz / dz);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj),
		mz = cz - static_cast<iter_type>(dz * ck);

	// position of center of each ellipse

	using p = struct { double x, y, z; };
	p pts[] = {
		p{ 0, 0, 0 },
		p{ dx, 0, 0 },
		p{ 0, 0, dz },
		p{ dx, 0, dz },
		p{ 0, dy, 0 },
		p{ dx, dy, 0 },
		p{ 0, dy, dz },
		p{ dx, dy, dz }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_3(mx, u.x, aa, my, u.y, bb, mz, u.z, cc))
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::CUBIC, InsideTag::VARA>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1],
		dz = dims[2] / init.data.gp[2];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[2]),
		bb = dy / (4.0 * init.data.gp[2]),
		cc = dy / (4.0 * init.data.gp[2]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = (n / dims[0]) % dims[1],
		cz = n / (dims[0] * dims[1]),
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy),
		ck = static_cast<iter_type>(cz / dz);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj),
		mz = cz - static_cast<iter_type>(dz * ck);

	double const
		rx = dx / 2,
		ry = dy / 2,
		rz = dz / 2;

	return symphas::internal::is_in_circle_3(mx, rx, aa, my, ry, bb, mz, rz, cc) ? IC_INNER_VALUE : IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::CUBIC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	// size of symmetry square
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1],
		dz = dims[2] / init.data.gp[2];

	// dimensions of ellipse
	double const
		aa = dx / (4.0 * init.data.gp[3]),
		bb = dy / (4.0 * init.data.gp[3]),
		cc = dz / (4.0 * init.data.gp[3]);

	// x y cursor position
	iter_type const
		cx = n % dims[0],
		cy = (n / dims[0]) % dims[1],
		cz = n / (dims[0] * dims[1]),
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy),
		ck = static_cast<iter_type>(cz / dz);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj),
		mz = cz - static_cast<iter_type>(dz * ck);

	iter_type const
		ci_len = static_cast<iter_type>(init.data.gp[0]),
		cj_len = static_cast<iter_type>(init.data.gp[1]),
		ck_len = static_cast<iter_type>(init.data.gp[2]),
		ci_next = (ci + 1) % ci_len,
		cj_next = (cj + 1) % cj_len,
		ck_next = (ck + 1) % ck_len;

	iter_type const
		i_top_l_f = ck * ci_len * cj_len + cj_next * ci_len + ci,
		i_top_r_f = ck * ci_len * cj_len + cj_next * ci_len + ci_next,
		i_bot_l_f = ck * ci_len * cj_len + cj * ci_len + ci,
		i_bot_r_f = ck * ci_len * cj_len + cj * ci_len + ci_next,
		i_top_l_b = ck_next * ci_len * cj_len + cj_next * ci_len + ci,
		i_top_r_b = ck_next * ci_len * cj_len + cj_next * ci_len + ci_next,
		i_bot_l_b = ck_next * ci_len * cj_len + cj * ci_len + ci,
		i_bot_r_b = ck_next * ci_len * cj_len + cj * ci_len + ci_next;

	using p = struct { double x, y, z; };
	p pts[] = {
		p{ offsets.get_delta(i_bot_l_f)[0], offsets.get_delta(i_bot_l_f)[1], offsets.get_delta(i_bot_l_f)[2] },
		p{ dx + offsets.get_delta(i_bot_r_f)[0], offsets.get_delta(i_bot_r_f)[1], offsets.get_delta(i_bot_r_f)[2] },
		p{ offsets.get_delta(i_bot_l_b)[0], offsets.get_delta(i_bot_l_b)[1], dz + offsets.get_delta(i_bot_l_b)[2] },
		p{ dx + offsets.get_delta(i_bot_r_b)[0], offsets.get_delta(i_bot_r_b)[1], dz + offsets.get_delta(i_bot_r_b)[2] },
		p{ offsets.get_delta(i_top_l_f)[0], dy + offsets.get_delta(i_top_l_f)[1], offsets.get_delta(i_top_l_f)[2] },
		p{ dx + offsets.get_delta(i_top_r_f)[0], dy + offsets.get_delta(i_top_r_f)[1], offsets.get_delta(i_top_r_f)[2] },
		p{ offsets.get_delta(i_top_l_b)[0], dy + offsets.get_delta(i_top_l_b)[1], dz + offsets.get_delta(i_top_l_b)[2] },
		p{ dx + offsets.get_delta(i_top_r_b)[0], dy + offsets.get_delta(i_top_r_b)[1], dz + offsets.get_delta(i_top_r_b)[2] }
	};

	for (auto const &u : pts)
	{
		if (symphas::internal::is_in_circle_3(mx, u.x, aa, my, u.y, bb, mz, u.z, cc))
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::HEXAGONAL>::operator[](iter_type n) const
{
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1],
		dz = dims[2] / init.data.gp[2];

	double const
		aa = dx / (4.0 * init.data.gp[3]),
		bb = dy / (4.0 * init.data.gp[3]),
		cc = dz / (4.0 * init.data.gp[3]);

	iter_type const
		cx = n % dims[0],
		cy = (n / dims[0]) % dims[1],
		cz = n / (dims[0] * dims[1]),
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy),
		ck = static_cast<iter_type>(cz / dz);

	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj),
		mz = cz - static_cast<iter_type>(dz * ck);

	using p = struct { double x, y, z; };
	p pts[] = {
		p{ dx / 2, 0, dz / 2 },
		p{ dx / 2, dy, dz / 2 },
		p{ 0, dy / 2, dz / 2 },
		p{ dx, dy / 2, dz / 2 },
		p{ dx / 2, dy / 2, 0 },
		p{ dx / 2, dy / 2, dz }
	};

	for (auto const& u : pts)
	{
		if (symphas::internal::is_in_circle_3(mx, u.x, aa, my, u.y, bb, mz, u.z, cc))
		{
			return IC_INNER_VALUE;
		}
	}

	return IC_OUTER_VALUE;
}

template<>
scalar_t InitialConditionsAlg<3, Inside::HEXAGONAL, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double const
		dx = dims[0] / init.data.gp[0],
		dy = dims[1] / init.data.gp[1],
		dz = dims[2] / init.data.gp[2];

	double const
		aa = dx / (4.0 * init.data.gp[3]),
		bb = dy / (4.0 * init.data.gp[3]),
		cc = dz / (4.0 * init.data.gp[3]);

	iter_type const
		cx = n % dims[0],
		cy = (n / dims[0]) % dims[1],
		cz = n / (dims[0] * dims[1]),
		ci = static_cast<iter_type>(cx / dx),
		cj = static_cast<iter_type>(cy / dy),
		ck = static_cast<iter_type>(cz / dz);

	// coordinates local to the symmetry square
	iter_type const
		mx = cx - static_cast<iter_type>(dx * ci),
		my = cy - static_cast<iter_type>(dy * cj),
		mz = cz - static_cast<iter_type>(dz * ck);

	iter_type const
		ci_len = static_cast<iter_type>(init.data.gp[0]),
		cj_len = static_cast<iter_type>(init.data.gp[1]),
		ck_len = static_cast<iter_type>(init.data.gp[2]),
		ci_next = (ci + 1) % ci_len,
		cj_next = (cj + 1) % cj_len,
		ck_next = (ck + 1) % ck_len;

	iter_type const
		i_front = ck * ci_len * cj_len * 3 + cj * ci_len + ci,
		i_back = ck_next * ci_len * cj_len * 3 + cj * ci_len + ci,
		i_bot = (ck * 3 + 1) * ci_len * cj_len + cj * ci_len * 2 + ci,
		i_top = (ck * 3 + 1) * ci_len * cj_len + cj_next * ci_len * 2 + ci,
		i_left = (ck * 3 + 1) * ci_len * cj_len + cj * ci_len * 2 + ci_len + ci,
		i_right = (ck * 3 + 1) * ci_len * cj_len + cj * ci_len * 2 + ci_len + ci_next;

	using p = struct { double x, y, z; };
	p pts[] = {
		p{ dx / 2 + offsets.get_delta(i_bot)[0], offsets.get_delta(i_bot)[1], dz / 2 + offsets.get_delta(i_bot)[2] },
		p{ dx / 2 + offsets.get_delta(i_top)[0], dy + offsets.get_delta(i_top)[1], dz / 2 + offsets.get_delta(i_top)[2] },
		p{ offsets.get_delta(i_left)[0], dy / 2 + offsets.get_delta(i_left)[1], dz / 2 + offsets.get_delta(i_left)[2] },
		p{ dx + offsets.get_delta(i_right)[0], dy / 2 + offsets.get_delta(i_right)[1], dz / 2 + offsets.get_delta(i_right)[2] },
		p{ dx / 2 + offsets.get_delta(i_front)[0], dy / 2 + offsets.get_delta(i_front)[1], offsets.get_delta(i_front)[2] },
		p{ dx / 2 + offsets.get_delta(i_back)[0], dy / 2 + offsets.get_delta(i_back)[1], dz + offsets.get_delta(i_back)[2] }
	};

	for (auto const &u : pts)
	{
		if (symphas::internal::is_in_circle_3(mx, u.x, aa, my, u.y, bb, mz, u.z, cc))
		{
			return IC_INNER_VALUE;
		}
	}
	return IC_OUTER_VALUE;
}






template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE>::operator[](iter_type n) const
{
	return symphas::internal::seeds_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}




template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &lengths, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, &lengths, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &lengths, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_A_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, &lengths, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}



template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &lengths, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARB>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, &lengths, 
		init.data.gp[2], init.data.gp[3]);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &lengths, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return symphas::internal::seeds_B_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets, &lengths, &values, 
		init.data.gp[2], init.data.gp[3], 0);
}




template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_3(
		n, dims, init, 
		&symphas::internal::is_in_circle_3, 
		&offsets,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_3(
		n, dims, init, 
		&symphas::internal::is_in_square_3, 
		&offsets, &values, 
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_rnd_3(
		n, dims, init,
		&symphas::internal::is_in_circle_3,
		&offsets, &values,
		seed_value, field_value, rnd_offset);
}



template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_A_3(
		n, dims, init,
		&symphas::internal::is_in_circle_3,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_A_rnd_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}


template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	return symphas::internal::seeds_B_3(
		n, dims, init,
		&symphas::internal::is_in_circle_3,
		&offsets, &lengths,
		seed_value, field_value);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>::operator[](iter_type n) const
{
	double seed_value = seed_value_dis(gen);
	double field_value = init.data.gp[4];
	double rnd_offset = seed_value - ((symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)) ? init.data.gp[3] : init.data.gp[2]);
	return symphas::internal::seeds_B_rnd_3(
		n, dims, init,
		&symphas::internal::is_in_square_3,
		&offsets, &lengths, &values,
		seed_value, field_value, rnd_offset);
}




scalar_t voronoi_value_1(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<1> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{
	int x = n;
	axis_1d_type p{ static_cast<axis_coord_t>(x) };

	axis_coord_t d0 = symphas::lib::length(
		axis_1d_type{
			static_cast<axis_coord_t>(dims[0]) });

	iter_type index = 0;
	for (iter_type i = 0; i < N; ++i)
	{
		axis_1d_type v{ offsets.get_delta(i) };
		axis_coord_t d = symphas::lib::distance(p, v);

		if (d < d0)
		{
			index = i;
			d0 = d;
		}
	}

	return values.get_offset(index);
}


scalar_t voronoi_value_2(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<2> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{
	int x = n % dims[0];
	int y = n / dims[0];

	axis_2d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y) };

	axis_coord_t d0 = symphas::lib::length(
		axis_2d_type{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]) });

	iter_type index = 0;
	for (iter_type i = 0; i < N; ++i)
	{
		axis_2d_type v{ offsets.get_delta(i)[0], offsets.get_delta(i)[1] };
		axis_coord_t d = symphas::lib::distance(p, v);

		if (d < d0)
		{
			index = i;
			d0 = d;
		}
	}

	return values.get_offset(index);
}

scalar_t voronoi_value_3(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<3> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{
	int x = n % dims[0];
	int y = (n / dims[0]) % dims[1];
	int z = n / (dims[0] * dims[1]);

	axis_3d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y),
		static_cast<axis_coord_t>(z) };

	axis_3d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]),
			static_cast<axis_coord_t>(dims[2]) };

	axis_coord_t d0 = symphas::lib::length(dm);

	iter_type index = 0;
	for (iter_type i = 0; i < N; ++i)
	{
		axis_3d_type v{ offsets.get_delta(i)[0], offsets.get_delta(i)[1], offsets.get_delta(i)[2] };
		axis_coord_t d = symphas::lib::distance(p, v);
		if (d < d0)
		{
			index = i;
			d0 = d;
		}
	}

	return values.get_offset(index);
}


scalar_t voronoi_value_A_1(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<1> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{
	int x = n;
	axis_1d_type p{ static_cast<axis_coord_t>(x) };

	axis_1d_type dm{
			static_cast<axis_coord_t>(dims[0]) };
	axis_coord_t d0 = symphas::lib::length(dm * 2);

	iter_type index = 0;
	for (iter_type i = 0; i < N; ++i)
	{
		axis_1d_type v0{ offsets.get_delta(i) - p };
		axis_1d_type v{
			v0 - dm * std::round(v0 / dm) };
		//axis_1d_type v_l{ offsets.get_delta(i) - dm };
		//axis_1d_type v_r{ offsets.get_delta(i) + dm };
		//
		//for (auto const& v : { v0, v_l, v_r })
		//{
			//axis_coord_t d = symphas::lib::distance(p, v);
		axis_coord_t d = symphas::lib::length(v);
		if (d < d0)
			{
				index = i;
				d0 = d;
			}
		//}
	}

	return values.get_offset(index);

}


scalar_t voronoi_value_A_2(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<2> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{

	int x = n % dims[0];
	int y = n / dims[0];

	axis_2d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y) };

	axis_2d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]) };

	axis_coord_t d0 = symphas::lib::length(
		axis_2d_type{
			dm[0] * 2,
			dm[1] * 2 });

	iter_type index = 0;
	for (iter_type i = 0; i < N; ++i)
	{
		axis_2d_type v0{ offsets.get_delta(i)[0] - p[0], offsets.get_delta(i)[1] - p[1] };

		axis_2d_type v{
			v0[0] - dm[0] * std::round(v0[0] / dm[0]),
			v0[1] - dm[1] * std::round(v0[1] / dm[1]) };

		//axis_2d_type v_t{ v0[0], v0[1] + dm[1] };
		//axis_2d_type v_b{ v0[0], v0[1] - dm[1] };
		//axis_2d_type v_r{ v0[0] + dm[0], v0[1] };
		//axis_2d_type v_l{ v0[0] - dm[0], v0[1] };
		//
		//axis_2d_type v_tr{ v0[0] + dm[0], v0[1] + dm[1] };
		//axis_2d_type v_tl{ v0[0] - dm[0], v0[1] + dm[1] };
		//axis_2d_type v_br{ v0[0] + dm[0], v0[1] - dm[1] };
		//axis_2d_type v_bl{ v0[0] - dm[0], v0[1] - dm[1] };
		//
		//for (auto const& v : { v0, v_t, v_b, v_r, v_l, v_tr, v_tl, v_br, v_bl })
		//{
		axis_coord_t d = symphas::lib::length(v);
		//axis_coord_t d = symphas::lib::distance(p, v);
			if (d < d0)
			{
				index = i;
				d0 = d;
			}
		//}
	}

	return values.get_offset(index);
}

scalar_t voronoi_value_A_3(iter_type n, const len_type* dims, size_t N,
	symphas::internal::RandomDeltas<3> const& offsets, symphas::internal::RandomOffsets<scalar_t, 1> const& values)
{
	int x = n % dims[0];
	int y = (n / dims[0]) % dims[1];
	int z = n / (dims[0] * dims[1]);

	axis_3d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y),
		static_cast<axis_coord_t>(z) };

	axis_3d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]),
			static_cast<axis_coord_t>(dims[2]) };

	axis_coord_t d0 = symphas::lib::length(
		axis_3d_type{
			dm[0] * 16,
			dm[1] * 16,
			dm[2] * 16 });

	iter_type index = 0;

	for (iter_type i = 0; i < N; ++i)
	{
		axis_3d_type v0{ offsets.get_delta(i)[0] - p[0], offsets.get_delta(i)[1] - p[1], offsets.get_delta(i)[2] - p[2]};

		axis_3d_type v{
			v0[0] - dm[0] * std::round(v0[0] / dm[0]),
			v0[1] - dm[1] * std::round(v0[1] / dm[1]),
			v0[2] - dm[2] * std::round(v0[2] / dm[2]) };

		//axis_3d_type v_t{ v0[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_b{ v0[0], v0[1] - dm[1], v0[2] };
		//axis_3d_type v_r{ v0[0] + dm[0], v0[1], v0[2] };
		//axis_3d_type v_l{ v0[0] - dm[0], v0[1], v0[2] };
		//axis_3d_type v_n{ v0[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_f{ v0[0], v0[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_nr{ v0[0] + dm[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_nl{ v0[0] - dm[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_fr{ v0[0] + dm[0], v0[1], v0[2] - dm[2] };
		//axis_3d_type v_fl{ v0[0] - dm[0], v0[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_tr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_tl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_br{ v0[0] + dm[0], v0[1] - dm[1], v0[2] };
		//axis_3d_type v_bl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] };
		//
		//axis_3d_type v_nt{ v0[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nb{ v0[0], v0[1] - dm[1], v0[2] + dm[2] };
		//axis_3d_type v_ft{ v0[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fb{ v0[0], v0[1] - dm[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_ntr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_ntl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nbr{ v0[0] + dm[0], v0[1] - dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nbl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] + dm[2] };
		//
		//axis_3d_type v_ftr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_ftl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fbr{ v0[0] + dm[0], v0[1] - dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fbl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] - dm[2] };
		//
		//for (auto const& v : {
		//	v0, v_n, v_f, v_t, v_b, v_r, v_l, v_tr, v_tl, v_br, v_bl,
		//	v_fr, v_fl, v_nr, v_nl, v_ft, v_fb, v_nt, v_nb,
		//	v_ntr, v_ntl, v_nbr, v_nbl, v_ftr, v_ftl, v_fbr, v_fbl })
		//{
			//axis_coord_t d = symphas::lib::distance(p, v);
			axis_coord_t d = symphas::lib::length(v);
			if (d < d0)
			{
				index = i;
				d0 = d;
			}
		//}
	}

	return values.get_offset(index);
}



template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI>::operator[](iter_type n) const
{
	return voronoi_value_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI>::operator[](iter_type n) const
{
	return voronoi_value_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI>::operator[](iter_type n) const
{
	return voronoi_value_3(n, dims, N, offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_3(n, dims, N, offsets, values);
}



template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return voronoi_value_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return voronoi_value_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return voronoi_value_3(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return voronoi_value_A_3(n, dims, N, offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_3(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_A_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_A_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return voronoi_value_A_3(n, dims, N, offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_3(n, dims, N, offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_A_1(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_A_2(n, dims, N, offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return voronoi_value_A_3(n, dims, N, offsets, values);
}


scalar_t bubble_value_1(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<1> const& offsets, 
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{
	int x = n;
	axis_1d_type p{ static_cast<axis_coord_t>(x) };

	axis_coord_t d0 = symphas::lib::length(
		axis_1d_type{
			static_cast<axis_coord_t>(dims[0]) });

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_1d_type v{ offsets.get_delta(i) };
		axis_coord_t d = symphas::lib::distance(p, v);

		if (d < R)
		{
			return values.get_offset(i);
		}
	}
	return default_value;
}


scalar_t bubble_value_2(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<2> const& offsets,
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{
	int x = n % dims[0];
	int y = n / dims[0];

	axis_2d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y) };

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_2d_type v{ offsets.get_delta(i)[0], offsets.get_delta(i)[1] };
		axis_coord_t d = symphas::lib::distance(p, v);

		if (d < R)
		{
			return values.get_offset(i);
		}
	}
	return default_value;
}

scalar_t bubble_value_3(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<3> const& offsets,
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{
	int x = n % dims[0];
	int y = (n / dims[0]) % dims[1];
	int z = n / (dims[0] * dims[1]);

	axis_3d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y),
		static_cast<axis_coord_t>(z) };

	axis_3d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]),
			static_cast<axis_coord_t>(dims[2]) };

	axis_coord_t d0 = symphas::lib::length(dm);

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_3d_type v{ offsets.get_delta(i)[0], offsets.get_delta(i)[1], offsets.get_delta(i)[2] };
		axis_coord_t d = symphas::lib::distance(p, v);

		if (d < R)
		{
			return values.get_offset(i);
		}
	}

	return default_value;
}


scalar_t bubble_value_A_1(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<1> const& offsets,
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{
	int x = n;
	axis_1d_type p{ static_cast<axis_coord_t>(x) };

	axis_1d_type dm{
			static_cast<axis_coord_t>(dims[0]) };
	axis_coord_t d0 = symphas::lib::length(dm * 2);

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_1d_type v0{ offsets.get_delta(i) - p };
		axis_1d_type v{ v0 - dm * std::round(v0 / dm) };
		//axis_1d_type v_l{ offsets.get_delta(i) - dm };
		//axis_1d_type v_r{ offsets.get_delta(i) + dm };

		//for (auto const& v : { v0, v_l, v_r })
		//{
			//axis_coord_t d = symphas::lib::distance(p, v);
			axis_coord_t d = symphas::lib::length(v);
			if (d < R)
			{
				return values.get_offset(i);
			}
		//}
	}

	return default_value;
}


scalar_t bubble_value_A_2(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<2> const& offsets,
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{

	int x = n % dims[0];
	int y = n / dims[0];

	axis_2d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y) };

	axis_2d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]) };

	axis_coord_t d0 = symphas::lib::length(
		axis_2d_type{
			dm[0] * 2,
			dm[1] * 2 });

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_2d_type v0{ offsets.get_delta(i)[0] - p[0], offsets.get_delta(i)[1] - p[1]};
		axis_2d_type v{
			v0[0] - dm[0] * std::round(v0[0] / dm[0]),
			v0[1] - dm[1] * std::round(v0[1] / dm[1]) };


		//axis_2d_type v_t{ v0[0], v0[1] + dm[1] };
		//axis_2d_type v_b{ v0[0], v0[1] - dm[1] };
		//axis_2d_type v_r{ v0[0] + dm[0], v0[1] };
		//axis_2d_type v_l{ v0[0] - dm[0], v0[1] };
		//
		//axis_2d_type v_tr{ v0[0] + dm[0], v0[1] + dm[1] };
		//axis_2d_type v_tl{ v0[0] - dm[0], v0[1] + dm[1] };
		//axis_2d_type v_br{ v0[0] + dm[0], v0[1] - dm[1] };
		//axis_2d_type v_bl{ v0[0] - dm[0], v0[1] - dm[1] };

		//for (auto const& v : { v0, v_t, v_b, v_r, v_l, v_tr, v_tl, v_br, v_bl })
		//{
			//axis_coord_t d = symphas::lib::distance(p, v);
			axis_coord_t d = symphas::lib::length(v);
			if (d < R)
			{
				return values.get_offset(i);
			}
		//}
	}

	return default_value;
}

scalar_t bubble_value_A_3(iter_type n, const len_type* dims, size_t N, double R,
	scalar_t default_value,
	symphas::internal::RandomDeltas<3> const& offsets,
	symphas::internal::RandomOffsets<scalar_t, 1> const& values,
	int select = -1)
{
	int x = n % dims[0];
	int y = (n / dims[0]) % dims[1];
	int z = n / (dims[0] * dims[1]);

	axis_3d_type p{
		static_cast<axis_coord_t>(x),
		static_cast<axis_coord_t>(y),
		static_cast<axis_coord_t>(z) };

	axis_3d_type dm{
			static_cast<axis_coord_t>(dims[0]),
			static_cast<axis_coord_t>(dims[1]),
			static_cast<axis_coord_t>(dims[2]) };

	axis_coord_t d0 = symphas::lib::length(
		axis_3d_type{
			dm[0] * 16,
			dm[1] * 16,
			dm[2] * 16 });

	int i_start = (select < 0) ? 0 : select;
	int i_end = (select < 0) ? offsets.size() : select + 1;
	for (iter_type i = i_start; i < i_end; ++i)
	{
		axis_3d_type v0{ offsets.get_delta(i)[0] - p[0], offsets.get_delta(i)[1] - p[1], offsets.get_delta(i)[2] - p[2] };

		axis_3d_type v{
			v0[0] - dm[0] * std::round(v0[0] / dm[0]),
			v0[1] - dm[1] * std::round(v0[1] / dm[1]),
			v0[2] - dm[2] * std::round(v0[2] / dm[2]) };

		//axis_3d_type v_t{ v0[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_b{ v0[0], v0[1] - dm[1], v0[2] };
		//axis_3d_type v_r{ v0[0] + dm[0], v0[1], v0[2] };
		//axis_3d_type v_l{ v0[0] - dm[0], v0[1], v0[2] };
		//axis_3d_type v_n{ v0[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_f{ v0[0], v0[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_nr{ v0[0] + dm[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_nl{ v0[0] - dm[0], v0[1], v0[2] + dm[2] };
		//axis_3d_type v_fr{ v0[0] + dm[0], v0[1], v0[2] - dm[2] };
		//axis_3d_type v_fl{ v0[0] - dm[0], v0[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_tr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_tl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] };
		//axis_3d_type v_br{ v0[0] + dm[0], v0[1] - dm[1], v0[2] };
		//axis_3d_type v_bl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] };
		//
		//axis_3d_type v_nt{ v0[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nb{ v0[0], v0[1] - dm[1], v0[2] + dm[2] };
		//axis_3d_type v_ft{ v0[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fb{ v0[0], v0[1] - dm[1], v0[2] - dm[2] };
		//
		//axis_3d_type v_ntr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_ntl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nbr{ v0[0] + dm[0], v0[1] - dm[1], v0[2] + dm[2] };
		//axis_3d_type v_nbl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] + dm[2] };
		//
		//axis_3d_type v_ftr{ v0[0] + dm[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_ftl{ v0[0] - dm[0], v0[1] + dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fbr{ v0[0] + dm[0], v0[1] - dm[1], v0[2] - dm[2] };
		//axis_3d_type v_fbl{ v0[0] - dm[0], v0[1] - dm[1], v0[2] - dm[2] };

		//for (auto const& v : {
		//	v0, v_n, v_f, v_t, v_b, v_r, v_l, v_tr, v_tl, v_br, v_bl,
		//	v_fr, v_fl, v_nr, v_nl, v_ft, v_fb, v_nt, v_nb,
		//	v_ntr, v_ntl, v_nbr, v_nbl, v_ftr, v_ftl, v_fbr, v_fbl })
		//{
			//axis_coord_t d = symphas::lib::distance(p, v);
			axis_coord_t d = symphas::lib::length(v);
			if (d < R)
			{
				return values.get_offset(i);
			}
		//}
	}

	return default_value;
}




template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE>::operator[](iter_type n) const
{
	return bubble_value_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE>::operator[](iter_type n) const
{
	return bubble_value_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE>::operator[](iter_type n) const
{
	return bubble_value_3(n, dims, N, R, init.data.gp[1], offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_3(n, dims, N, R, init.data.gp[1], offsets, values);
}



template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return bubble_value_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return bubble_value_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>::operator[](iter_type n) const
{
	return bubble_value_3(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>::operator[](iter_type n) const
{
	return bubble_value_A_3(n, dims, N, R, init.data.gp[1], offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_3(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_A_1(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_A_2(n, dims, N, R, init.data.gp[1], offsets, values);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>::operator[](iter_type n) const
{
	return bubble_value_A_3(n, dims, N, R, init.data.gp[1], offsets, values);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_1(n, dims, N, R, init.data.gp[1], offsets, values, select);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_2(n, dims, N, R, init.data.gp[1], offsets, values, select);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_3(n, dims, N, R, init.data.gp[1], offsets, values, select);
}


template<>
scalar_t InitialConditionsAlg<1, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_A_1(n, dims, N, R, init.data.gp[1], offsets, values, select);
}

template<>
scalar_t InitialConditionsAlg<2, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_A_2(n, dims, N, R, init.data.gp[1], offsets, values, select);
}

template<>
scalar_t InitialConditionsAlg<3, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>::operator[](iter_type n) const
{
	return bubble_value_A_3(n, dims, N, R, init.data.gp[1], offsets, values, select);
}






