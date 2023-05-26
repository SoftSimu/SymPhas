
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
		return InsideTag::NONE;
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
	{"BUBBLES", Inside::BUBBLE},
	{"BUBBLE", Inside::BUBBLE},
	{"BB", Inside::BUBBLE},
	{"SIN", Inside::SIN},
	{"COS", Inside::COS},
	{"EXPRESSION", Inside::EXPRESSION}
};



std::map<const char*, InsideTag, symphas::internal::any_case_comparator> symphas::internal::init_tag_key_map = {
	{"DEFAULT", InsideTag::DEFAULT},
	{"RANDOM", InsideTag::RANDOM},
	{"RAND", InsideTag::RANDOM},
	{"R", InsideTag::RANDOM},
	{"INVERT", InsideTag::INVERT},
	{"INV", InsideTag::INVERT},
	{"VARA", InsideTag::VARA},
	{"A", InsideTag::VARA},
	{"VARB", InsideTag::VARB},
	{"B", InsideTag::VARB},
	{"VARC", InsideTag::VARC},
	{"C", InsideTag::VARC},
	{"FIXEDSEED", InsideTag::FIXEDSEED},
	{"FIXED", InsideTag::FIXEDSEED},
	{"FIX", InsideTag::FIXEDSEED}
};

namespace symphas::internal
{


	double compute_bubble_overlap(double coverage, double R, double ratio)
	{
		double Ac = R * R * symphas::PI;
		double Rr = sqrt(Ac) / 2;
		return (R - Rr) * ratio;
	}

	double compute_bubble_overlap_range(double coverage, double R, double ratio)
	{
		return compute_bubble_overlap(coverage, R, ratio) / 2;
		//double overlap = compute_bubble_overlap(coverage);
		//return (overlap > -2)
		//	? 0.168887 * std::pow(overlap - compute_bubble_overlap(1), 2) + compute_bubble_overlap(1)
		//	: 1;
	}

	double compute_bubble_R(symphas::interval_data_type const& vdata, double coverage, size_t N)
	{
		N = (N == 0) ? 1 : N;

		len_type A = 1;
		for (auto const& [key, value] : vdata)
		{
			A *= value.get_count();
		}
		return std::sqrt(A * coverage / (symphas::PI * N));
	}


	auto get_bubbles_in_range_A(
		axis_2d_type pos, double Rc, 
		std::vector<axis_2d_type> const& positions, 
		const len_type* dims)
	{
		std::vector<axis_2d_type> prox;
		prox.reserve(positions.size());

		axis_2d_type dm{
				static_cast<axis_coord_t>(dims[0]),
				static_cast<axis_coord_t>(dims[1]) };

		for (const auto& v0 : positions)
		{
			axis_2d_type v_t{ v0[0], v0[1] + dm[1] };
			axis_2d_type v_b{ v0[0], v0[1] - dm[1] };
			axis_2d_type v_r{ v0[0] + dm[0], v0[1] };
			axis_2d_type v_l{ v0[0] - dm[0], v0[1] };

			axis_2d_type v_tr{ v0[0] + dm[0], v0[1] + dm[1] };
			axis_2d_type v_tl{ v0[0] - dm[0], v0[1] + dm[1] };
			axis_2d_type v_br{ v0[0] + dm[0], v0[1] - dm[1] };
			axis_2d_type v_bl{ v0[0] - dm[0], v0[1] - dm[1] };

			auto check = { v0, v_t, v_b, v_r, v_l, v_tr, v_tl, v_br, v_bl };
			for (auto const& v : check)
			{
				if (symphas::lib::distance(pos, v) < Rc)
				{
					//prox.push_back(v);
					prox.insert(prox.end(), check);
				}
			}
		}
		return prox;
	}


	auto get_bubbles_in_range(
		axis_2d_type pos, double Rc,
		std::vector<axis_2d_type> const& positions)
	{
		std::vector<axis_2d_type> prox;
		prox.reserve(positions.size());

		for (const auto& pos0 : positions)
		{
			if (symphas::lib::distance(pos, pos0) < Rc)
			{
				prox.push_back(pos0);
			}
		}
		return prox;
	}

	//! Generate the position of bubbles in a 2d system.
	/*!
	 * Generate the position of bubbles in a 2d system.
	 *
	 * \param N The number of bubbles to generate.
	 * \param R The radius of the bubble.
	 * \param overlap How much bubbles can be separated.
	 * \param x0 The first x position.
	 * \param y0 The first y position.
	 */
	symphas::internal::RandomDeltas<1> get_bubble_positions_1(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, iter_type x0)
	{
		symphas::internal::RandomDeltas<1> positions(N, dims);
		return positions;
	}

	//! Generate the position of bubbles in a 2d system.
	/*!
	 * Generate the position of bubbles in a 2d system.
	 *
	 * \param N The number of bubbles to generate.
	 * \param R The radius of the bubble.
	 * \param max_overlap How much bubbles can be separated.
	 * \param x0 The first x position.
	 * \param y0 The first y position.
	 */
	symphas::internal::RandomDeltas<2> get_bubble_positions_2(
		size_t N, 
		double R, 
		double max_overlap, 
		symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, 
		const len_type* dims, 
		iter_type x0, iter_type y0)
	{
		std::vector<axis_2d_type> positions;
		positions.reserve(N);
		
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, 1);

		axis_2d_type current_pos = {
			static_cast<axis_coord_t>((x0 < 0) ? x0 + dims[0] : x0), 
			static_cast<axis_coord_t>((y0 < 0) ? y0 + dims[1] : y0)};

		// Generate N bubbles, starting at the first random delta.
		for (iter_type n = 0; n < N - 1; ++n)
		{
			double overlap = overlaps.get_offset(n);
			double Rs = 2 * R - overlap;			// the spawn radius from the current bubble
			double Re = 2 * R - max_overlap;		// the exclusion radius of other bubbles

			double RR = -Re * Re + Rs * Rs;
			double _R = 1 / (2 * Rs);

			auto all_nearby = get_bubbles_in_range_A(current_pos, Re + Rs - max_overlap, positions, dims);
			positions.push_back(current_pos);

			std::vector<std::pair<double, double>> exclusions;
			std::vector<std::pair<double, double>> ranges;
			for (axis_2d_type prox_pos : all_nearby)
			{
				axis_2d_type relative_pos = prox_pos - current_pos;
				double r = symphas::lib::distance(prox_pos, current_pos);

				if (r <= abs(Re - Rs))
				{
					exclusions.push_back({ 0, 2 * symphas::PI });
					break;
				}
				else if (abs(RR / r + r) <= 2 * Rs)
				{
					double alpha = atan2(relative_pos[1], relative_pos[0]);
					double u1 = acos(_R * (RR / r + r));
					double u0 = -u1;
					double alph0 = u0 + alpha;
					double alph1 = u1 + alpha;

					if (alph0 < 0 && alph1 < 0)
					{
						alph0 += 2 * symphas::PI;
						alph1 += 2 * symphas::PI;
					}

					exclusions.push_back({ alph0, alph1 });
				}
			}

			double theta;
			if (exclusions.size() > 0)
			{

				sort(exclusions.begin(), exclusions.end(),
					[] (auto a, auto b) { return a.second > b.second; });
				double alph = exclusions.front().second - 2 * symphas::PI;

				sort(exclusions.begin(), exclusions.end(),
					[] (auto a, auto b) { return a.first < b.first; });
				
				double alph00 = exclusions.front().first + 2 * symphas::PI;
				for (auto [alph0, alph1] : exclusions)
				{
					if (alph < alph0)
					{
						ranges.push_back({ alph, std::min(alph0, alph00) });
					}
					alph = std::max(alph, alph1);
				}

				if (ranges.size() == 1)
				{
					auto [beta0, beta1] = ranges.front();
					if (beta0 > beta1) ranges.clear();
				}

				// range testing to validate ranges
				//for (auto [alph0, alph1] : ranges)
				//{
				//	for (auto alphtest : { alph0 + 1e-4, alph1 - 1e-4 })
				//	{
				//		double Xtest = std::cos(alphtest);
				//		double Ytest = std::sin(alphtest);

				//		double x1test = current_pos[0] + Rs * Xtest;
				//		double y1test = current_pos[1] + Rs * Ytest;
				//		axis_2d_type check_pos0{
				//			(x1test <= 0) ? x1test + dims[0] : (x1test >= dims[0]) ? x1test - dims[0] : x1test,
				//			(y1test <= 0) ? y1test + dims[1] : (y1test >= dims[1]) ? y1test - dims[1] : y1test };


				//		auto nearby0 = get_bubbles_in_range_A(check_pos0, Re, positions, dims);
				//		if (nearby0.size() > 0)
				//		{

				//			double pointertest = 0;
				//			double r = dims[0];
				//			double rr = dims[0];
				//			for (auto bubble : nearby0)
				//			{
				//				double dist = symphas::lib::distance(bubble, current_pos);
				//				double dist2 = symphas::lib::distance(bubble, check_pos0);
				//				rr = std::min(rr, dist2);

				//				if (dist < r)
				//				{
				//					r = dist;
				//					pointertest = atan2(bubble[1] - current_pos[1], bubble[0] - current_pos[0]);
				//				}
				//			}

				//			double u1test = acos(_R * (RR / r + r));
				//			double u0test = -u1test;
				//			double alph0test = u0test + pointertest;
				//			double alph1test = u1test + pointertest;

				//			double Xtest0 = std::cos(pointertest);
				//			double Ytest0 = std::sin(pointertest);
				//			symphas::internal::RandomDeltas<2> positions0(n + 2, dims);
				//		}
				//	}

				//}

				// the reason there's overlap is in the case the exclusion bounds are -ve to -ve.
				if (ranges.size() == 0)
				{
					iter_type x = iter_type(current_pos[0] + R);
					iter_type y = iter_type(current_pos[1] + R);
					iter_type d0 = iter_type(R / 4);

					bool flag = false;
					for (iter_type j = 0; j < dims[0] && !flag; ++j, y += d0)
					{
						y = (y >= dims[1]) ? y - dims[1] : y;
						for (iter_type i = 0; i < dims[0] && !flag; ++i, x += d0)
						{
							x = (x >= dims[0]) ? x - dims[0] : x;
							axis_2d_type pos{ double(x), double(y) };

							if (get_bubbles_in_range_A(pos, Re, positions, dims).size() == 0)
							{
								current_pos = pos;
								flag = true;
							}
						}
					}
					if (!flag)
					{
						return positions;
					}
					else
					{
						continue;
					}
				}
				else if (ranges.size() == 1)
				{
					auto [beta0, beta1] = ranges.front();
					theta = (beta1 - beta0) * dis(gen) + beta0;
				}
				else
				{
					double len = 0;
					for (auto [beta0, beta1] : ranges)
					{
						len += (beta1 - beta0);
					}
					double pick = dis(gen) * len;
					for (auto [beta0, beta1] : ranges)
					{
						pick -= (beta1 - beta0);
						if (pick <= 0)
						{
							theta = beta1 + pick;
							break;
						}
					}
				}
			}
			else
			{
				theta = 2 * symphas::PI * dis(gen);
			}

			double X = std::cos(theta);
			double Y = std::sin(theta);

			double x1 = current_pos[0] + Rs * X;
			double y1 = current_pos[1] + Rs * Y;
			current_pos = {
				(x1 <= 0) ? x1 + dims[0] : (x1 >= dims[0]) ? x1 - dims[0] : x1,
				(y1 <= 0) ? y1 + dims[1] : (y1 >= dims[1]) ? y1 - dims[1] : y1 };
		}
		positions.push_back(current_pos);

		//positions.sort();
		return positions;
	}

	symphas::internal::RandomDeltas<3> get_bubble_positions_3(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, iter_type x0, iter_type y0, iter_type z0)
	{
		symphas::internal::RandomDeltas<3> positions(N, dims);
		return positions;
	}


}