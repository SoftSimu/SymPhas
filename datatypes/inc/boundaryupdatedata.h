
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
 * MODULE:  datatypes
 * PURPOSE: Defines how to update the boundaries, depending on the type
 * and dimension.
 *
 * ***************************************************************************
 */

#pragma once

#include "griditer.h"
#include "grid.h"
#include "boundary.h"


namespace symphas::internal
{
	// Gives the origin and size of the boundary in the global domain. "Boundary" refers to the 
	// ghost cells themselves, which will be updated using real data from the interior.
	// The origin member represents the position in the region that the boundary starts, 
	// and the dimension is the size of the region.
	template<size_t D, Side... s>
	struct boundary_interval;

	// left edge
	template<>
	struct boundary_interval<1, Side::LEFT>
	{
		boundary_interval(const len_type(&dims)[1], len_type boundary_size) :
			origin{ 0 },
			dims{ boundary_size },
			stride{ 1 }
		{}

		iter_type origin[1];
		len_type dims[1];
		len_type stride[1];
	};

	// left edge
	template<>
	struct boundary_interval<1, Side::RIGHT>
	{
		boundary_interval(const len_type(&dims)[1], len_type boundary_size) :
			origin{ dims[0] - 1 },
			dims{ -boundary_size },
			stride{ -1 }
		{}

		iter_type origin[1];
		len_type dims[1];
		len_type stride[1];
	};

	//***********************************************************************
	// 2D
	//***********************************************************************

	//***********************************************************************
	// Fully periodic
	//***********************************************************************

	// left edge
	template<>
	struct boundary_interval<2, Side::LEFT, Side::LEFT>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, boundary_size },
			dims{ boundary_size, dims[1] - boundary_size * 2 },
			stride{ dims[1] - boundary_size * 2, 1 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// top left corner
	template<>
	struct boundary_interval<2, Side::LEFT, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, dims[1] - boundary_size },
			dims{ boundary_size, boundary_size },
			stride{ 1, boundary_size }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// bottom left corner
	template<>
	struct boundary_interval<2, Side::LEFT, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, 0 },
			dims{ boundary_size, boundary_size },
			stride{ 1, boundary_size }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};


	// right edge
	template<>
	struct boundary_interval<2, Side::RIGHT, Side::RIGHT>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ dims[0] - 1, boundary_size },
			dims{ -boundary_size, dims[1] - boundary_size * 2 },
			stride{ -dims[1] + boundary_size * 2, 1 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// top right corner
	template<>
	struct boundary_interval<2, Side::RIGHT, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ dims[0] - 1, dims[1] - boundary_size },
			dims{ -boundary_size, boundary_size },
			stride{ -1, boundary_size }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// bottom right corner
	template<>
	struct boundary_interval<2, Side::RIGHT, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ dims[0] - 1, 0 },
			dims{ -boundary_size, boundary_size },
			stride{ -1, boundary_size }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// top edge
	template<>
	struct boundary_interval<2, Side::TOP, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ boundary_size, dims[1] - 1 },
			dims{ dims[0] - boundary_size * 2, -boundary_size },
			stride{ 1, -dims[0] + boundary_size * 2 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	// bottom edge
	template<>
	struct boundary_interval<2, Side::BOTTOM, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ boundary_size, 0 },
			dims{ dims[0] - boundary_size * 2, boundary_size },
			stride{ 1, dims[0] - boundary_size * 2 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};


	//***********************************************************************
	// Periodic in one direction
	//***********************************************************************

	template<>
	struct boundary_interval<2, Side::LEFT>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, 0 },
			dims{ boundary_size, dims[1] },
			stride{ dims[1], 1 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	template<>
	struct boundary_interval<2, Side::RIGHT>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ dims[0] - 1, 0 },
			dims{ -boundary_size, dims[1] },
			stride{ -dims[1], 1 }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	template<>
	struct boundary_interval<2, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, dims[0] - 1 },
			dims{ dims[0], -boundary_size },
			stride{ 1, -dims[0] }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	template<>
	struct boundary_interval<2, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[2], len_type boundary_size) :
			origin{ 0, 0 },
			dims{ dims[0], boundary_size },
			stride{ 1, dims[0] }
		{}

		iter_type origin[2];
		len_type dims[2];
		len_type stride[2];
	};

	//***********************************************************************
	// 3D
	//***********************************************************************


	//***********************************************************************
	// Fully periodic
	//***********************************************************************


	template<>
	struct boundary_interval<3, Side::LEFT, Side::LEFT, Side::LEFT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, boundary_size, boundary_size },
			dims{ boundary_size, dims[1] - boundary_size * 2, dims[2] - boundary_size * 2 },
			stride{ (dims[1] - boundary_size * 2) * (dims[2] - boundary_size * 2), (dims[1] - boundary_size * 2), 1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::LEFT, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, boundary_size, 0 },
			dims{ boundary_size, dims[1] - boundary_size * 2, boundary_size },
			stride{ (dims[1] - boundary_size * 2) * boundary_size, (dims[1] - boundary_size * 2), 1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::LEFT, Side::BACK>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, boundary_size, dims[2] - 1 },
			dims{ boundary_size, dims[1] - boundary_size * 2, -boundary_size },
			stride{ (dims[1] - boundary_size * 2) * boundary_size, (dims[1] - boundary_size * 2), -1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::LEFT, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, dims[1] - 1, boundary_size },
			dims{ boundary_size, -boundary_size, dims[2] - boundary_size * 2 },
			stride{ (dims[2] - boundary_size * 2) * boundary_size, -(dims[2] - boundary_size * 2), 1}
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::LEFT, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, boundary_size },
			dims{ boundary_size, boundary_size, dims[2] - boundary_size * 2 },
			stride{ (dims[2] - boundary_size * 2) * boundary_size, (dims[2] - boundary_size * 2), 1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::TOP, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, dims[1] - 1, 0 },
			dims{ boundary_size, -boundary_size, boundary_size },
			stride{ 1, -boundary_size, boundary_size * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::TOP, Side::BACK>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, dims[1] - 1, dims[2] - 1 },
			dims{ boundary_size, -boundary_size, -boundary_size },
			stride{ 1, -boundary_size, -boundary_size * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::BOTTOM, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, 0 },
			dims{ boundary_size, boundary_size, boundary_size },
			stride{ 1, boundary_size, -boundary_size * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::LEFT, Side::BOTTOM, Side::BACK>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, dims[2] - 1 },
			dims{ boundary_size, boundary_size, -boundary_size },
			stride{ 1, boundary_size, -boundary_size * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::RIGHT, Side::RIGHT, Side::RIGHT>
		: boundary_interval<3, Side::LEFT, Side::LEFT, Side::LEFT>
	{
		using parent_type = boundary_interval<3, Side::LEFT, Side::LEFT, Side::LEFT>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) 
		{
			parent_type::origin[0] = dims[0] - 1;
			parent_type::dims[0] = -parent_type::dims[0];
			parent_type::stride[0] = -parent_type::stride[0];
		} 
	};

	template<Side side0>
	struct boundary_interval<3, Side::RIGHT, Side::RIGHT, side0>
		: boundary_interval<3, Side::LEFT, Side::LEFT, side0>
	{
		using parent_type = boundary_interval<3, Side::LEFT, Side::LEFT, side0>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			parent_type::origin[0] = dims[0] - 1;
			parent_type::dims[0] = -parent_type::dims[0];
			parent_type::stride[0] = -parent_type::stride[0];
		}
	};

	template<Side side0, Side side1>
	struct boundary_interval<3, Side::RIGHT, side0, side1>
		: boundary_interval<3, Side::LEFT, side0, side1>
	{
		using parent_type = boundary_interval<3, Side::LEFT, side0, side1>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			parent_type::origin[0] = dims[0] - 1;
			parent_type::dims[0] = -parent_type::dims[0];
			parent_type::stride[0] = -parent_type::stride[0];
		}
	};

	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ boundary_size, 0, boundary_size },
			dims{ dims[0] - boundary_size * 2, boundary_size, dims[2] - boundary_size * 2 },
			stride{ 1, (dims[0] - boundary_size * 2) * (dims[2] - boundary_size * 2), (dims[0] - boundary_size * 2) }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};


	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ boundary_size, 0, 0 },
			dims{ (dims[0] - boundary_size * 2), boundary_size, boundary_size },
			stride{ 1, (dims[0] - boundary_size * 2), (dims[0] - boundary_size * 2) * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BACK>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ boundary_size, 0, dims[2] - 1 },
			dims{ (dims[0] - boundary_size * 2), boundary_size, -boundary_size },
			stride{ 1, (dims[0] - boundary_size * 2), -(dims[0] - boundary_size * 2) * boundary_size }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::TOP, Side::TOP, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ boundary_size, dims[1] - 1, boundary_size },
			dims{ dims[0] - boundary_size * 2, -boundary_size, dims[2] - boundary_size * 2 },
			stride{ 1, -(dims[0] - boundary_size * 2) * (dims[2] - boundary_size * 2), (dims[0] - boundary_size * 2) }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<Side side0>
	struct boundary_interval<3, Side::TOP, Side::TOP, side0>
		: boundary_interval<3, Side::BOTTOM, Side::BOTTOM, side0>
	{
		using parent_type = boundary_interval<3, Side::BOTTOM, Side::BOTTOM, side0>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			parent_type::origin[1] = dims[1] - 1;
			parent_type::dims[1] = -parent_type::dims[1];
			parent_type::stride[1] = -parent_type::stride[1];
		}
	};
	template<>
	struct boundary_interval<3, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ boundary_size, boundary_size, 0 },
			dims{ dims[0] - boundary_size * 2, dims[1] - boundary_size * 2, boundary_size },
			stride{ 1, dims[0] - boundary_size * 2, (dims[0] - boundary_size * 2) * (dims[1] - boundary_size * 2) }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::BACK, Side::BACK, Side::BACK>
		: boundary_interval<3, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		using parent_type = boundary_interval<3, Side::FRONT, Side::FRONT, Side::FRONT>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			parent_type::origin[2] = dims[2] - 1;
			parent_type::dims[2] = -parent_type::dims[2];
			parent_type::stride[2] = -parent_type::stride[2];
		}
	};


	//***********************************************************************
	// Periodic in one direction
	//***********************************************************************

	template<>
	struct boundary_interval<3, Side::LEFT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, 0 },
			dims{ boundary_size, dims[1], dims[2] },
			stride{ dims[1] * dims[2], dims[1], 1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::RIGHT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ dims[0] - 1, 0, 0 },
			dims{ -boundary_size, dims[1], dims[2] },
			stride{ -dims[1] * dims[2], dims[1], 1 }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::TOP>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, dims[1] - 1, 0 },
			dims{ dims[0], -boundary_size, dims[2] },
			stride{ 1, -dims[0] * dims[2], dims[0] }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::BOTTOM>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, 0 },
			dims{ dims[0], boundary_size, dims[2] },
			stride{ 1, dims[0] * dims[2], dims[0] }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::FRONT>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, 0 },
			dims{ dims[0], dims[1], boundary_size },
			stride{ 1, dims[0], dims[0] * dims[1] }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};

	template<>
	struct boundary_interval<3, Side::BACK>
	{
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			origin{ 0, 0, dims[2] - 1 },
			dims{ dims[0], dims[1], -boundary_size },
			stride{ 1, dims[0], -dims[0] * dims[1] }
		{}

		iter_type origin[3];
		len_type dims[3];
		len_type stride[3];
	};


	//***********************************************************************
	// Periodic in x and y
	//***********************************************************************

	template<Side side0>
	struct boundary_interval<3, Side::LEFT, side0>
		: boundary_interval<3, Side::LEFT, Side::LEFT, side0>
	{
		using parent_type = boundary_interval<3, Side::LEFT, Side::LEFT, side0>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<Side side0>
	struct boundary_interval<3, Side::RIGHT, side0>
		: boundary_interval<3, Side::RIGHT, Side::RIGHT, side0>
	{
		using parent_type = boundary_interval<3, Side::RIGHT, Side::RIGHT, side0>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::BOTTOM>
		: boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>
	{
		using parent_type = boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::TOP, Side::TOP>
		: boundary_interval<3, Side::TOP, Side::TOP, Side::TOP>
	{
		using parent_type = boundary_interval<3, Side::TOP, Side::TOP, Side::TOP>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	//***********************************************************************
	// Periodic in y and z
	//***********************************************************************


	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::FRONT>
		: boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::FRONT>
	{
		using parent_type = boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::FRONT>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::BOTTOM, Side::BACK>
		: boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BACK>
	{
		using parent_type = boundary_interval<3, Side::BOTTOM, Side::BOTTOM, Side::BACK>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::TOP, Side::FRONT>
		: boundary_interval<3, Side::TOP, Side::TOP, Side::FRONT>
	{
		using parent_type = boundary_interval<3, Side::TOP, Side::TOP, Side::FRONT>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::TOP, Side::BACK>
		: boundary_interval<3, Side::TOP, Side::TOP, Side::BACK>
	{
		using parent_type = boundary_interval<3, Side::TOP, Side::TOP, Side::BACK>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::FRONT, Side::FRONT>
		: boundary_interval<3, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		using parent_type = boundary_interval<3, Side::FRONT, Side::FRONT, Side::FRONT>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct boundary_interval<3, Side::BACK, Side::BACK>
		: boundary_interval<3, Side::BACK, Side::BACK, Side::BACK>
	{
		using parent_type = boundary_interval<3, Side::BACK, Side::BACK, Side::BACK>;
		boundary_interval(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};


	template<size_t D>
	struct offset_data
	{
		template<typename... Ts>
		offset_data(iter_type arg0, Ts... args) : offset{ arg0, args... } {}

		using arr_ref_t = iter_type[D];

		operator arr_ref_t& () { return offset; }
		operator const arr_ref_t& () const { return offset; }

		arr_ref_t& data() { return offset; }

		iter_type operator[](iter_type i) { return offset[i]; }
		iter_type offset[D];

	};

	// and the dimension is the size of the region.
	template<size_t D, BoundaryType type, Side... sides>
	struct periodic_offset;

	//***********************************************************************
	// 1D
	//***********************************************************************

	template<>
	struct periodic_offset<1, BoundaryType::PERIODIC, Side::LEFT> : offset_data<1>
	{
		periodic_offset(const len_type(&dims)[1], len_type boundary_size) :
			offset_data<1>{ dims[0] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<1, BoundaryType::PERIODIC, Side::RIGHT> 
		: periodic_offset<1, BoundaryType::PERIODIC, Side::LEFT>
	{
		using parent_type = periodic_offset<1, BoundaryType::PERIODIC, Side::LEFT>;
		periodic_offset(const len_type(&dims)[1], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};


	//***********************************************************************
	// 2D
	//***********************************************************************


	//***********************************************************************
	// Fully periodic
	//***********************************************************************


	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT> : offset_data<2>
	{
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			offset_data<2>{ dims[0] - boundary_size * 2, 0 } {}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::TOP> : offset_data<2>
	{
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			offset_data<2>{ dims[0] - boundary_size * 2, -dims[0] + boundary_size * 2 } {}

	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::TOP>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::TOP>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[1] = -offset[1];
		}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT>;

		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT, Side::BOTTOM>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT, Side::TOP>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::TOP>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::TOP>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM> : offset_data<2>
	{
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			offset_data<2>{ 0, dims[1] - boundary_size * 2 } {}

	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::TOP, Side::TOP>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[1] = -offset[1];
		}
	};


	//***********************************************************************
	// Periodic along one axis
	//***********************************************************************

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT>;

		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT>;

		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<2, BoundaryType::PERIODIC, Side::TOP>
		: periodic_offset<2, BoundaryType::PERIODIC, Side::TOP, Side::TOP>
	{
		using parent_type = periodic_offset<2, BoundaryType::PERIODIC, Side::TOP, Side::TOP>;
		periodic_offset(const len_type(&dims)[2], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};


	//***********************************************************************
	// 3D
	//***********************************************************************


	//***********************************************************************
	// Fully periodic
	//***********************************************************************


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, 0, 0 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::FRONT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, 0, dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BACK> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, 0, -dims[2] + boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::TOP> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, -dims[1] + boundary_size * 2, 0 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BOTTOM> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, dims[1] - boundary_size * 2, 0 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::TOP, Side::BACK> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, -dims[1] + boundary_size * 2, -dims[2] + boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::TOP, Side::FRONT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, -dims[1] + boundary_size * 2, dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM, Side::BACK> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, dims[1] - boundary_size * 2, -dims[2] + boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM, Side::FRONT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ dims[0] - boundary_size * 2, dims[1] - boundary_size * 2, dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::RIGHT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>;
		using parent_type::offset;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<Side side0>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, side0>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, side0>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, side0>;
		using parent_type::offset;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<Side side0, Side side1>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, side0, side1> 
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, side0, side1>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, side0, side1>;
		using parent_type::offset;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[0] = -offset[0];
		}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::TOP> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ 0, -dims[1] + boundary_size * 2, 0 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::FRONT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ 0, -dims[1] + boundary_size * 2, dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::BACK> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ 0, -dims[1] - boundary_size * 2, -dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ 0, dims[1] - boundary_size * 2, 0 } {}
	};

	template<Side side0>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, side0>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, side0>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, side0>;
		using parent_type::offset;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[2] = -offset[2];
		}
	};


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT> : offset_data<3>
	{
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			offset_data<3>{ 0, 0, dims[2] - boundary_size * 2 } {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>;
		using parent_type::offset;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size)
		{
			offset[2] = -offset[2];
		}
	};

	//***********************************************************************
	// Periodic only along one axis
	//***********************************************************************


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::RIGHT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::RIGHT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::TOP>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::TOP>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	//***********************************************************************
	// Periodic along both x and z
	//***********************************************************************


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::LEFT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::RIGHT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::RIGHT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BACK, Side::BACK, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::FRONT, Side::FRONT, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	//***********************************************************************
	// Periodic along both x and y
	//***********************************************************************


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::TOP>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::TOP>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::TOP>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::BOTTOM>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BOTTOM>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::LEFT, Side::LEFT, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::TOP>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::TOP>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::TOP>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::BOTTOM>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::BOTTOM>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::RIGHT, Side::RIGHT, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::TOP>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::TOP>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	//***********************************************************************
	// Periodic along both y and z
	//***********************************************************************


	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::TOP, Side::TOP, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::FRONT>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::FRONT>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::FRONT>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};

	template<>
	struct periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BACK>
		: periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BACK>
	{
		using parent_type = periodic_offset<3, BoundaryType::PERIODIC, Side::BOTTOM, Side::BOTTOM, Side::BACK>;
		periodic_offset(const len_type(&dims)[3], len_type boundary_size) :
			parent_type(dims, boundary_size) {}
	};


	



	//***********************************************************************
	//
	//***********************************************************************



	template<Side... sides, typename T, size_t D>
	void get_overlap(iter_type(&overlaps)[D][2], RegionalGrid<T, D> const& grid, const iter_type(&origin)[D], const len_type(&dims)[D])
	{
		iter_type adjusted_origin[D]{};
		for (iter_type i = 0; i < D; ++i)
		{
			adjusted_origin[i] = ((symphas::axis_to_index(symphas::axis_of_side(sides)) == i) || ...) ? grid.region.origin[i] : 0;
		}

		for (iter_type i = 0; i < D; ++i)
		{
			overlaps[i][0] = std::max(adjusted_origin[i], origin[i]);
			overlaps[i][1] = std::min(adjusted_origin[i] + grid.region.dims[i], origin[i] + dims[i]);
		}
	}

	template<Side... sides, typename T, size_t D>
	void get_overlap_adjusted(iter_type(&overlaps)[D][2], RegionalGrid<T, D> const& grid, const iter_type(&origin)[D], const len_type(&dims)[D])
	{
		iter_type adjusted_origin[D];
		iter_type adjusted_dims[D];
		for (iter_type i = 0; i < D; ++i)
		{
			if (dims[i] < 0)
			{
				adjusted_origin[i] = origin[i] + dims[i] + 1;
				adjusted_dims[i] = -dims[i];
			}
			else
			{
				adjusted_origin[i] = origin[i];
				adjusted_dims[i] = dims[i];
			}
		}
		return get_overlap<sides...>(overlaps, grid, adjusted_origin, adjusted_dims);
	}

	template<typename T, size_t D, Side... sides>
	void get_overlap(iter_type(&overlaps)[D][2], RegionalGrid<T, D> const& grid, boundary_interval<D, sides...> const& interval)
	{
		get_overlap_adjusted<sides...>(overlaps, grid, interval.origin, interval.dims);
	}


	template<size_t D, size_t... Is>
	bool is_overlapping(const len_type(&overlaps)[D][2], std::index_sequence<Is...>)
	{
		return ((overlaps[Is][1] - overlaps[Is][0] > 0) && ...);
	}

	template<size_t D>
	bool is_overlapping(const len_type(&overlaps)[D][2])
	{
		return is_overlapping(overlaps, std::make_index_sequence<D>{});
	}

	template<size_t D>
	void get_dims_and_origin(iter_type(&origin)[D], len_type(&dims)[D], const len_type(&overlaps)[D][2])
	{
		for (iter_type i = 0; i < D; ++i)
		{
			origin[i] = overlaps[i][0];
			dims[i] = overlaps[i][1] - overlaps[i][0];
		}
	}

	/* 2 dimensional boundaries
	 *
	 */

	template<Side... sides, typename T, size_t D>
	iter_type get_regional_index(iter_type(&pos)[D], RegionalGrid<T, D>& grid)
	{
		for (iter_type i = 0; i < D; ++i)
		{
			if (((symphas::axis_to_index(symphas::axis_of_side(sides)) == i) || ...))
			{
				pos[i] = (pos[i] >= grid.region.origin[i]) ? pos[i] - grid.region.origin[i] : pos[i] - grid.region.origin[i] + grid.dims[i];
			}
		}
		return grid::index_from_position(pos, grid.region.stride);
	}

	template<Side... sides, typename T, size_t D>
	bool side_non_wrapping(RegionalGrid<T, D>& grid)
	{
		bool flag = true;
		for (iter_type i : { symphas::axis_to_index(symphas::axis_of_side(sides))... })
		{
			if (grid.region.origin[i] + grid.region.dims[i] > grid.dims[i])
			{
				flag = false;
			}
		}
		return flag;
	}

	template<typename T, size_t D, Side... sides>
	void regional_update_boundary(symphas::lib::side_list<sides...>, RegionalGrid<T, D>& grid)
	{
		// to update the LEFT face. 
		// first see if the grid region overlaps the left boundary.
		boundary_interval<D, sides...> interval(grid.dims, grid.region.boundary_size);

		iter_type overlaps[D][2];
		get_overlap(overlaps, grid, interval);

		if (is_overlapping(overlaps) && side_non_wrapping<sides...>(grid))
		{
			get_dims_and_origin(interval.origin, interval.dims, overlaps);

			// shift origin to region to copy the boundary from.
			periodic_offset<D, BoundaryType::PERIODIC, sides...> offset(grid.dims, grid.region.boundary_size);

			for (iter_type i = 0; i < D; ++i)
			{
				interval.origin[i] += offset[i];
			}

			get_overlap(overlaps, grid, interval);
			if (is_overlapping(overlaps))
			{
				get_dims_and_origin(interval.origin, interval.dims, overlaps);

				iter_type m = grid::index_from_position(offset.data(), grid.region.stride);
				for (iter_type n = 0; n < grid::length<D>(interval.dims); ++n)
				{
					iter_type pos[D]{};
					grid::get_grid_position_offset(pos, interval.dims, interval.origin, n);
					iter_type index = get_regional_index<sides...>(pos, grid);
					grid.values[index - m] = grid.values[index];
				}
			}
		}
	}

	inline iter_type get_boundary_offset(symphas::lib::side_list<Side::LEFT>, const iter_type(&pos)[2], const iter_type(&origin)[2])
	{
		return std::abs(pos[1] - origin[1]);
	}

	inline iter_type get_boundary_offset(symphas::lib::side_list<Side::RIGHT>, const iter_type(&pos)[2], const iter_type(&origin)[2])
	{
		return std::abs(pos[1] - origin[1]);
	}

	inline iter_type get_boundary_offset(symphas::lib::side_list<Side::TOP>, const iter_type(&pos)[2], const iter_type(&origin)[2])
	{
		return std::abs(pos[0] - origin[0]);
	}

	inline iter_type get_boundary_offset(symphas::lib::side_list<Side::BOTTOM>, const iter_type(&pos)[2], const iter_type(&origin)[2])
	{
		return std::abs(pos[0] - origin[0]);
	}



	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::LEFT>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[2] - origin[2]), std::abs(pos[1] - origin[1]) };
	}

	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::RIGHT>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[2] - origin[2]), std::abs(pos[1] - origin[1]) };
	}

	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::TOP>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[0] - origin[0]), std::abs(pos[2] - origin[2]) };
	}

	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::BOTTOM>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[0] - origin[0]), std::abs(pos[2] - origin[2]) };
	}

	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::FRONT>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[0] - origin[0]), std::abs(pos[1] - origin[1]) };
	}

	inline std::pair<iter_type, iter_type> get_boundary_offset(symphas::lib::side_list<Side::BACK>, const iter_type(&pos)[3], const iter_type(&origin)[3])
	{
		return { std::abs(pos[0] - origin[0]), std::abs(pos[1] - origin[1]) };
	}

	template<Side side, size_t D>
	decltype(auto) get_boundary_offset(const iter_type(&pos)[D], const iter_type(&origin)[D])
	{
		return get_boundary_offset(symphas::lib::side_list<side>{}, pos, origin);
	}

	template<Side side0, Side... sides, typename T>
	void set_boundary_position(double(&pos0)[2], const grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT>* bd, const iter_type(&pos)[1], const iter_type(&origin)[1], double v)
	{
		pos0[0] = v;
		pos0[1] = 0;
	}

	template<Side side0, Side... sides, typename T>
	void set_boundary_position(double(&pos0)[2], const grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT>* bd, const iter_type(&pos)[2], const iter_type(&origin)[2], const double(&v)[2])
	{
		pos0[0] = v[0] + get_boundary_offset<side0>(pos, origin) * bd->h;
		pos0[1] = 0;
	}

	template<Side side0, Side... sides, typename T>
	void set_boundary_position(double(&pos0)[2], const grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT>* bd, const iter_type(&pos)[3], const iter_type(&origin)[3], const double(&v)[4])
	{
		auto [_0, _1] = get_boundary_offset<side0>(pos, origin);
		pos0[0] = v[0] + _0 * bd->h[0];
		pos0[1] = v[2] + _1 * bd->h[1];
	}

	template<typename T, size_t Dm1>
	void apply_boundary_function(const grid::BoundaryApplied<T, Dm1, BoundaryType::DEFAULT>* bd, carry_value<T> const& value, const double(&pos)[2], double time)
	{
		bd->update(*value.value, pos[0], pos[1], time);
	}

	template<typename T, size_t D, Side... sides>
	void regional_update_boundary(symphas::lib::side_list<sides...>, const grid::Boundary<T, D - 1>* b, RegionalGrid<T, D>& grid, double time)
	{
		// to update the LEFT face. 
		// first see if the grid region overlaps the left boundary.
		boundary_interval<D, sides...> interval(grid.dims, grid.region.boundary_size);

		iter_type overlaps[D][2];
		get_overlap(overlaps, grid, interval);

		if (is_overlapping(overlaps))
		{
			iter_type origin[D];
			len_type dims[D];
			get_dims_and_origin(origin, dims, overlaps);

			// shift origin to region to copy the boundary from.
			periodic_offset<D, BoundaryType::PERIODIC, sides...> offset(grid.dims, grid.region.boundary_size);

			for (iter_type i = 0; i < D; ++i)
			{
				origin[i] += offset[i];
			}

			get_overlap(overlaps, grid, origin, dims);
			if (is_overlapping(overlaps))
			{
				get_dims_and_origin(origin, dims, overlaps);

				iter_type stride[D]{};
				grid::get_stride(stride, grid.dims);

				iter_type m = grid::index_from_position(offset.data(), stride);
				
				auto* bd = static_cast<grid::BoundaryApplied<T, D - 1, BoundaryType::DEFAULT> const*>(b);
				
				for (iter_type n = 0; n < grid::length<D>(dims); ++n)
				{
					iter_type pos[D]{};
					grid::get_grid_position_offset(pos, dims, origin, n);
					iter_type index = grid::index_from_position(pos, stride);

					double pos0[2]{};
					set_boundary_position<sides...>(pos0, bd, pos, origin, bd->v);
					apply_boundary_function(bd, grid[index - m], pos0, time);
				}

			}
		}
	}
}

