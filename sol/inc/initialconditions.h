
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
 * MODULE:  sol
 * PURPOSE: Defines the functionality around applying intial conditions
 * to a phase field system.
 *
 * ***************************************************************************
 */

#pragma once

#include <algorithm>
#include <complex>
#include <random>
#include <ctime>
#include <cstring>


#ifdef USING_IO
#include "write.h"
#include "read.h"
#endif

#include "initialconditionslib.h"
#include "initialconditionsexpression.h"


#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif


 /*!
  * \defgroup initial Initial Conditions and File Loading
  * Initial conditions for a system are quite flexible; they can be defined from an algorithm, a file or even a symbolic expression.
  * 
  * \section Overview
  * 
  * Initial conditions can be either taken from the pool of built-in algorithms,
  * be user-generated, or be loaded in from a file.
  * 
  * The built in algorithms are individually listed and described on the
  * main page, and the specifics are given with their respective
  * documentations. The documentation for particular initial conditions can
  * be associated with an initial condition algorithm type and tag
  * by finding the object `InitialConditionsAlg<D, TYPE, TAGS...>`, where
  * `TYPE` is the initial condition algorithm value from enum Inside, and
  * `TAGS` are a list of InsideTag values. The list of initial condition
  * types and their tags are:
  * 
  * | Value | Description | Variations |
  * |------------|--------------------------------------|-------------|
  * |`GAUSSIAN`| 		Uses a Gaussian distribution to randomly assign values.|  *None* |
  * |`UNIFORM`| 			Uses a uniform distribution to randomly assign values.| *None* |
  * |`CAPPED`| 			Values are assigned to be either the minimum or maximum parameter.| *None* |
  * |`CONSTANT`|			All values are assigned to the same parameter.| *None* |
  * |`CIRCLE`| 			Value are assigned to be the shape of a circle.| `RANDOM`, `A`, `A+RANDOM` |
  * |`HEXAGONAL`| 		Values are assigned into circles arranged in a hexagonal pattern.| `RANDOM` |
  * |`CUBIC`| 			Values are assigned into circles arranged in a cubic pattern.| `RANDOM` |
  * |\ref initsquare |			Values are assigned to be the shape of a square.| `RANDOM`, `A`, `A+RANDOM` |
  * |`SQUARESEEDS`| 		Values are put into randomly arranged squares.| `RANDOM`, `A`, `A+RANDOM`, `B`, `B+RANDOM`  |
  * |`CIRCLESEEDS`|		Values are put into randomly arranged circles.| `RANDOM`, `A`, `A+RANDOM`, `B`, `B+RANDOM` |
  * |`FILE`| 			Values are read in from a file.| *None* |
  * |`CHECKPOINT`|		Values are read in from a checkpoint.| *None* |
  * |`NONE`|				Represents no initial condition.| *None* |
  * 
  * |Value | Description |
  * |--------|---------------|
  * | `DEFAULT`|	The default initial generation algorithm is chosen.|
  * | `RANDOM`|		The generation is modified to include some kind of randomness.|
  * | `VARA`|		The A variation is chosen (different from the default variation).|
  * | `VARB`|		The B variation is chosen.|
  * | `INVERT`|		The interior and outer values are switched in the generation. <br> (In some algorithms, the generated value is multiplied by -1.)|
  * | `NONE`|	Represents no tag.|
  *
  * If an initial condition is generated with <span style="color:violet">`Inside`</span>`::NONE` 
  * and <span style="color:violet">`InsideTag`</span>`::NONE`, then no values are populated. If an 
  * initial condition is combined with an invalid modifier, the program reports an error and no 
  * values will be populated. The default constructed  <span style="color:teal">`init_entry_type`</span> 
  * is initialized with <span style="color:violet">`Inside`</span>`::NONE` and 
  * <span style="color:violet">`InsideTag`</span>`::NONE`, so it is valid but will not generate 
  * any initial conditions.
  *
  * 
  * \section Loading From a File
  * 
  * Loading from a file has been as streamlined as possible. There are two
  * possible ways to load from a file: 1) loading from a checkpoint, and 2) 
  * loading from a named file. Typically, we will want to load from a named file.
  * 
  * To load from a named file, the file **must** be in checkpoint format.
  * When we are using the default WriterType::GNU option, we can load from text files.
  * The checkpoint file uses the following format for text files:
  * 
  * ```
  * dim dimension0 dimension1 x0 x1 y0 y1 index
  * ```
  * 
  * For three dimensions, another dimension and interval would be included. This
  * information only appears at the top of a file. The grid data for the given
  * index `index` then follows where each row is on a new line. For every additional 
  * grid block that is also included in the file, only the index is needed. 
  * The index of the next data block will not be separated from the prevoius
  * block by any empty new lines. 
  *  
  * Now, create the initial conditions parameter and the system can be initialized.
  * 
  * ```
  * symphas::init_entry_type tdata{ Inside::FILE, { index, "path/to/file/filename" } };
  * System<scalar_t, 2> sys(tdata, vdata);
  * ```
  * 
  * The index has to correspond to an index available in the file. Also,
  * care has to be paid attention to ensure that the dimension and the given
  * intervals `vdata` corresponds with the data being loaded. 
  * 
  * @{
  */


template<size_t D>
struct InitialConditionsData
{
	//! Create information for initial conditions.
	/*!
	 * Using the given initial condition data \p tdata, the initial
	 * conditions data for a system with dimensions \p dims is generated.
	 *
	 * \param tdata The initial conditions information.
	 * \param dims The dimensions of the system.
	 */
	InitialConditionsData(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		init{ init }, vdata{ vdata }, dims{ 0 }
	{
		std::copy(dims, dims + D, this->dims);
	}

	virtual scalar_t operator[](iter_type) const = 0;

	virtual operator bool() const
	{
		return true;
	}

	virtual ~InitialConditionsData() {};

	symphas::init_entry_type init;		//!< Initial condition parameters.
	symphas::interval_data_type vdata;	//!< Interval parameters.
	len_type dims[D];					//!< Dimensions of the system.

};


//! Used to generate initial conditions data for a system of given dimension.
/*!
 * Given dimensions of a system, a specific routine can be
 * given the index of in (flattened) array of the system, and the initial
 * conditions at that point will be returned. The algorithms
 * are typically designed for Cartesian grids.
 *
 * \tparam D Dimension of the system.
 * \tparam in The initial condition this algorithm generates.
 * \tparam tags The tags which modify the initial condition.
 */
template<size_t D, Inside in, InsideTag... tags>
struct InitialConditionsAlg;


namespace symphas::internal
{
	template<InsideTag... tags>
	struct init_tag_values_list {};
	template<Inside... ins>
	struct init_values_list {};

	template<Inside in>
	struct tags_for_init_value
	{
		using type = symphas::lib::types_list<>;
	};

	template<size_t D, Inside in>
	struct make_new_ic_delegate
	{
		InitialConditionsData<D>* operator()(
			symphas::init_entry_type const& init,
			symphas::interval_data_type const& vdata,
			len_type const* origin) const
		{
			if (init.in == in)
			{
				return check_tags_start(init, vdata, origin, typename tags_for_init_value<in>::type{});
			}
			else
			{
				return nullptr;
			}
		}

	protected:

		InitialConditionsData<D>* check_tags(
			symphas::init_entry_type const& init,
			symphas::interval_data_type const& vdata,
			len_type const* origin) const
		{
			return new InitialConditionsAlg<D, in>(init, vdata, grid::interior_dimensions(vdata));
		}

		template<typename... other_tag_types, InsideTag... tags>
		InitialConditionsData<D>* check_tags(
			symphas::init_entry_type const& init,
			symphas::interval_data_type const& vdata,
			len_type const* origin,
			init_tag_values_list<tags...>, other_tag_types... other_tags) const
		{
			size_t check_tag = build_intag(tags...);
			if (init.intag == check_tag)
			{
				return new InitialConditionsAlg<D, in, tags...>(init, vdata, grid::interior_dimensions(vdata));
			}
			else
			{
				return check_tags(init, vdata, origin, other_tags...);
			}
		}

		template<typename... tag_tuple_types>
		InitialConditionsData<D>* check_tags_start(
			symphas::init_entry_type const& init,
			symphas::interval_data_type const& vdata,
			len_type const* origin,
			symphas::lib::types_list<tag_tuple_types...>) const
		{
			return check_tags(init, vdata, origin, tag_tuple_types{}...);
		}

		
	};
}





/***************************************************************************
**SQUARE INITIAL CONDITIONS************************************************/

/*!
 * \addtogroup initial
 * \subsection initsquare Square
 * Generate a square centered at the halfway point of the dimensions.
 * 
 * The side lengths of the square will be scaled with respect to the
 * dimensions of the system. A rectangular system will contain
 * a rectangle in the middle.
 * 
 * The first `D` parameters are factors that divide the side length
 * with respect to the system dimensions. For instance, if the first
 * two parameters are equal to 2, then the generated square will have its
 * side lengths equal to half the dimension length.
 * 
 * In describing the function of the parameters, let \f$ \bar{x}, \bar{y}\f$
 * and \f$\bar{z}\f$ represent the midpoints of their respective intervals,
 * thus \f$(\bar{x}, \bar{y}, \bar{z})\f$ is the center point of the 
 * system. Moreover, let \f$\Delta x, \Delta y\f$ and \f$\Delta z\f$
 * be the position relative to the center of the system, i.e.
 * \f$\Delta x = x - \bar{x} \f$.
 * 
 * - In 2D, two parameters are given: \f$r_x, r_y\f$.
 * The square will have \ref params::init_inside_val "interior" values assigned for the points:
 * \f{gather*}
 *		\Delta x < A \quad \text{and} \quad  \Delta y < B
 * \f} and \ref params::init_outside_val "outer" values otherwise, where
 * \f{gather*}
 *		2A = (x_1 - x_0) / r_x \quad \text{and} \quad
 *		2B = (y_1 - y_0) / r_y
 * \f} and subscripts 0 and 1 indicate the start and end of the interval.
 * Thus, the width and height are equal to half the system
 * dimension, divided by the given parameters.
 * 
 * - In 3D, three parameters are given: \f$r_x, r_y, r_z\f$. The equation
 * of the sphere follows the same recipe as for 2D.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SQUARE> : 
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generate a square randomly offset from the center.
/*!
 * See InitialConditionsAlg<D, Inside::SQUARE>. The offset is chosen by applying
 * #IC_SQUARE_RND_FACTOR to RandomOffsets.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SQUARE, InsideTag::RANDOM> : 
	InitialConditionsAlg<D, Inside::SQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SQUARE>;
	
	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		offsets = symphas::internal::RandomDeltas<D>(1, dims, IC_SQUARE_RND_FACTOR);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomDeltas<D> offsets;		//!< Manages a list of random offsets.
};

//! Generate a square centered at the halfway point of the dimensions.
/*!
 * See InitialConditionsAlg<D, Inside::SQUARE>.
 *
 * This variation will always generate a square with equal sides. The
 * smallest side is chosen to scale the side of the square.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SQUARE, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SQUARE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generate a square randomly offset from the center.
/*!
 * See InitialConditionsAlg::square_A(). The offset is chosen by applying
 * #IC_SQUARE_RND_FACTOR to RandomOffsets.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generate a square randomly offset from the center.
/*!
 * See InitialConditionsAlg<D, Inside::SQUARE, InsideTag::VARA>. The offset is chosen by applying
 * #IC_SQUARE_RND_FACTOR to RandomOffsets.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SQUARE, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM> 
{
	using parent_type = InitialConditionsAlg<D, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


template<>
struct InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA> :
	InitialConditionsAlg<1, Inside::SQUARE>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE>;
	using parent_type::parent_type;
};

template<>
struct InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<1, Inside::SQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


template<>
struct symphas::internal::tags_for_init_value<Inside::SQUARE>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA, InsideTag::RANDOM>>;
};

/***************************************************************************
**CIRCLE INITIAL CONDITIONS************************************************/

//! Generate a circle centered at the halfway point of the dimensions.
/*!
 * The radii of the circle will be scaled with respect to the
 * dimensions of the system. A rectangular system will contain an
 * ellipse in the middle, thus the major and minor radii are aligned
 * with the axis and are the ones which are scaled.
 * 
 * The first `D` parameters are factors that divide the radii
 * with respect to the system dimensions. For instance, if the first
 * two parameters are equal to 2, then the generated circle will have its
 * radii equal to half the dimension length.
 * 
 * In describing the function of the parameters, let \f$ \bar{x}, \bar{y}\f$
 * and \f$\bar{z}\f$ represent the midpoints of their respective intervals,
 * thus \f$(\bar{x}, \bar{y}, \bar{z})\f$ is the center point of the 
 * system. Moreover, let \f$\Delta x, \Delta y\f$ and \f$\Delta z\f$
 * be the position relative to the center of the system, i.e.
 * \f$\Delta x = x - \bar{x} \f$.
 * 
 * - In 2D, two parameters are given: \f$r_x, r_y\f$.
 * The circle will have \ref params::init_inside_val "interior" values assigned for the points:
 * \f{gather*}
 *		1 > (\Delta x/A)^2 + (\Delta y/B)^2
 * \f} and \ref params::init_outside_val "outer" values otherwise, where
 * \f{gather*}
 *		2A = (x_1 - x_0) / r_x \quad \text{and} \quad
 *		2B = (y_1 - y_0) / r_y
 * \f} and subscripts 0 and 1 indicate the start and end of the interval.
 * Thus, the primary and secondary radii are equal to half the system
 * dimension, divided by the given parameters.
 * 
 * - In 3D, three parameters are given: \f$r_x, r_y, r_z\f$. The equation
 * of the sphere follows the same recipe as for 2D.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CIRCLE> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generate a circle randomly offset from the center.
/*!
 * See InitialConditionsAlg<D, Inside::CIRCLE>. The offset is chosen by applying
 * #IC_CIRCLE_RND_FACTOR to RandomOffsets.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::CIRCLE>
{
	using parent_type = InitialConditionsAlg<D, Inside::CIRCLE>;
	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		offsets = symphas::internal::RandomDeltas<D>(1, dims, IC_CIRCLE_RND_FACTOR);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomDeltas<D> offsets;		//!< Manages a list of random offsets.
};

//! Generate a circle centered at the halfway point of the dimensions.
/*!
 * See InitialConditionsAlg<D, Inside::CIRCLE>.
 *
 * This variation will always generate a circle with uniform radii. The
 * smallest side is chosen to scale the radius.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::CIRCLE>
{
	using parent_type = InitialConditionsAlg<D, Inside::CIRCLE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generate a circle randomly offset from the center.
/*!
 * See InitialConditionsAlg::circle_A(). The offset is chosen by applying
 * #IC_CIRCLE_RND_FACTOR to RandomOffsets.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::RANDOM, InsideTag::VARA>
	: InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM> 
{
	using parent_type = InitialConditionsAlg<D, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};



template<>
struct InitialConditionsAlg<1, Inside::CIRCLE> :
	InitialConditionsAlg<1, Inside::SQUARE>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE>;
	using parent_type::parent_type;
};

template<>
struct InitialConditionsAlg<1, Inside::CIRCLE, InsideTag::RANDOM> :
	InitialConditionsAlg<1, Inside::SQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
};

template<>
struct InitialConditionsAlg<1, Inside::CIRCLE, InsideTag::VARA> :
	InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA>;
	using parent_type::parent_type;
};

template<>
struct InitialConditionsAlg<1, Inside::CIRCLE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<1, Inside::SQUARE, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


template<>
struct symphas::internal::tags_for_init_value<Inside::CIRCLE>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA, InsideTag::RANDOM>>;
};

/***************************************************************************
**CUBIC INITIAL CONDITIONS*************************************************/

//! Generates circles arranged in a cubic lattice.
/*!
 * For 2 dimensions, this is a square lattice arrangement. The first
 * `D` parameters represent the number of tiling blocks for the
 * arrangement. In particular, the blocks are chosen by splitting
 * each dimension in the number of parts indicated by the corresponding
 * parameter. The `D + 1` parameter represents the scaling size of the 
 * spot.
 * 
 * In describing the function of the parameters, let \f$ \hat{x}, \hat{y}\f$
 * and \f$\hat{z}\f$ represent the width of the interval, i.e.
 * \f$\hat{x} = (x_1 - x_0)\f$.
 * 
 * - In 2D, three parameters are given: \f$c_x, c_y, R\f$. The first two
 * are used to partition the system, forming subsystems \f$S_{n, m}\f$,
 * \f$n=1,\ldots, c_x\f$ and \f$m=1,\dots, c_y\f$ defined as the set of 
 * points in the interval 
 * \f[
 *		\left[x_0 + (n-1)\frac{\hat{x}}{c_x}, x_0 + n\frac{\hat{x}}{c_x}\right) \bigcup
 *		\left[y_0 + (m-1)\frac{\hat{y}}{c_y}, y_0 + m\frac{\hat{y}}{c_y}\right).
 * \f]
 * A node is then placed at the intersection of four subsystems, 
 * where subsystems are considered periodic with respect to the
 * whole system (the first and last subsystem horizontally and vertically
 * are also considered adjacent).
 * In other words, for \f$n=1,\ldots, c_x\f$ and \f$m=1,\dots, c_y\f$, 
 * nodes are placed at all the points containing coordinates:
 * \f{eqnarray*}{
 *		& x_0 + (n-1)\frac{\hat{x}}{c_x}, \\ 
 *		& x_0 + n\frac{\hat{x}}{c_x}, \\
 *		& y_0 + (m-1)\frac{\hat{y}}{c_y}, \\ 
 *		& y_0 + m\frac{\hat{y}}{c_y}.
 * \f}
 * At each node, points are initialized with 
 * \ref params::init_inside_val "interior" values in an ellipse 
 * defined by the function:
 * \f{gather*}
 *		1 > (\Delta x/A)^2 + (\Delta y/B)^2
 * \f} where \f$\Delta x\f$ and \f$\Delta y\f$ represent coordinates
 * relative the position of the node, and \f$A\f$ and \f$B\f$ are defined:
 * \f{gather*}
 *		4A = \frac{\hat{x}}{Rc_x} \quad \text{and} \quad
 *		4B = \frac{\hat{x}}{Rc_y}.
 * \f}
 * Everywhere else, values are initialized as
 * \ref params::init_outside_val "outer" values.
 * 
 * Thus, the primary and secondary radii of each ellipse are defined as
 * a quarter of the length of the respective dimension of a subsystem, 
 * divided by the third parameter \f$R\f$. 
 * 
 * - In 3D, four parameters are given: \f$c_x, c_y, c_z, R\f$. THe first 
 * three are used to partition the system using the same recipe as in the
 * 2D case. The third parameter then adjusts the size of each circle,
 * also like the 2D case.
 */

template<size_t D>
struct InitialConditionsAlg<D, Inside::CUBIC> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Generates circles randomly offset arranged in a cubic lattice.
/*!
 * See InitialConditionsAlg<D, Inside::CUBIC>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CUBIC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::CUBIC>
{
	using parent_type = InitialConditionsAlg<D, Inside::CUBIC>;
	using parent_type::parent_type;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		size_t prod = 1;
		len_type local_dims[D];
		for (iter_type i = 0; i < D; ++i)
		{
			prod *= static_cast<size_t>(init.data.gp[i]);
			local_dims[i] = static_cast<len_type>(dims[i] / init.data.gp[i]);
		}
		offsets = symphas::internal::RandomDeltas<D>(prod, local_dims, IC_CUBIC_RND_FACTOR);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomDeltas<D> offsets;		//!< Manages a list of random offsets.
};

//! Generates circles randomly offset arranged in a cubic lattice.
/*!
 * See InitialConditionsAlg<D, Inside::CUBIC>.
 * 
 * The algorithm adjusts the placement of the square to be in the center
 * of the tile instead of in the corners. 
 * In other words, for \f$n=1,\ldots, c_x\f$ and \f$m=1,\dots, c_y\f$, 
 * nodes are placed at all the points containing coordinates:
 * \f{align*}
 *		& x_0 + (n-1)\frac{\hat{x}}{2c_x} \\
 *		& y_0 + (m-1)\frac{\hat{y}}{2c_y}
 * \f}
 *
 * An implication is that this
 * algorithm with the first `D` parameters set to 1 is similar to
 * the InitialConditionsAlg<D, Inside::CIRCLE> function, with a different radius.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::CUBIC, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::CUBIC>
{
	using parent_type = InitialConditionsAlg<D, Inside::CUBIC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


template<>
struct symphas::internal::tags_for_init_value<Inside::CUBIC>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::RANDOM>>;
};

/***************************************************************************
**HEXAGONAL INITIAL CONDITIONS*********************************************/


//! Generates circles arranged in a hexagonal lattice.
/*!
 * For 3 dimensions, this is a tetrahedral lattice arrangement. The first
 * `D` parameters represent the number of tiling blocks for the
 * arrangement. In particular, the blocks are chosen by splitting
 * each dimension in the number of parts indicated by the corresponding
 * parameter. The `D + 1` parameter represents the scaling size of the
 * spot.
 * 
 * In describing the function of the parameters, let \f$ \hat{x}, \hat{y}\f$
 * and \f$\hat{z}\f$ represent the width of the interval, i.e.
 * \f$\hat{x} = (x_1 - x_0)\f$.
 * 
 * - In 2D, three parameters are given: \f$c_x, c_y, R\f$. The first two
 * are used to partition the system, forming subsystems \f$S_{n, m}\f$,
 * \f$n=1,\ldots, c_x\f$ and \f$m=1,\dots, c_y\f$ defined as the set of 
 * points in the interval 
 * \f[
 *		\left[x_0 + (n-1)\frac{\hat{x}}{c_x}, x_0 + n\frac{\hat{x}}{c_x}\right) \bigcup
 *		\left[y_0 + (m-1)\frac{\hat{y}}{c_y}, y_0 + m\frac{\hat{y}}{c_y}\right).
 * \f]
 * A node is then placed at the center of the edge of two adjacent 
 * subsystems, where subsystems are considered periodic with respect to the
 * whole system (the first and last subsystem horizontally and vertically
 * are also considered adjacent).
 * In other words, for \f$n=1,\ldots, c_x\f$ and \f$m=1,\dots, c_y\f$, 
 * nodes are placed at the middle points of all intervals.
 * At each node, points are initialized with 
 * \ref params::init_inside_val "interior" values in an ellipse 
 * defined by the function:
 * \f{gather*}
 *		1 > (\Delta x/A)^2 + (\Delta y/B)^2
 * \f} where \f$\Delta x\f$ and \f$\Delta y\f$ represent coordinates
 * relative the position of the node, and \f$A\f$ and \f$B\f$ are defined:
 * \f{gather*}
 *		4A = \frac{\hat{x}}{Rc_x} \quad \text{and} \quad
 *		4B = \frac{\hat{x}}{Rc_y}.
 * \f}
 * Everywhere else, values are initialized as
 * \ref params::init_outside_val "outer" values.
 * 
 * Thus, the primary and secondary radii of each ellipse are defined as
 * a quarter of the length of the respective dimension of a subsystem, 
 * divided by the third parameter \f$R\f$. 
 * 
 * - In 3D, four parameters are given: \f$c_x, c_y, c_z, R\f$. THe first 
 * three are used to partition the system using the same recipe as in the
 * 2D case. The third parameter then adjusts the size of each circle,
 * also like the 2D case.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::HEXAGONAL> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Generates circles randomly offset arranged in a hexagonal lattice.
/*!
 * See InitialConditionsAlg<D, Inside::HEXAGONAL>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::HEXAGONAL, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::HEXAGONAL>
{
	using parent_type = InitialConditionsAlg<D, Inside::HEXAGONAL>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		size_t prod = D;
		len_type local_dims[D];
		for (iter_type i = 0; i < D; ++i)
		{
			prod *= static_cast<size_t>(init.data.gp[i]);
			local_dims[i] = static_cast<len_type>(dims[i] / init.data.gp[i]);
		}
		offsets = symphas::internal::RandomDeltas<D>(prod, local_dims, IC_HEX_RND_FACTOR);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomDeltas<D> offsets;		//!< Manages a list of random offsets.
};



template<>
struct InitialConditionsAlg<1, Inside::HEXAGONAL> :
	InitialConditionsAlg<1, Inside::CUBIC>
{
	using parent_type = InitialConditionsAlg<1, Inside::CUBIC>;
	using parent_type::parent_type;
};

template<>
struct InitialConditionsAlg<1, Inside::HEXAGONAL, InsideTag::RANDOM> :
	InitialConditionsAlg<1, Inside::CUBIC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<1, Inside::CUBIC, InsideTag::RANDOM>;
	using parent_type::parent_type;
};




template<>
struct symphas::internal::tags_for_init_value<Inside::HEXAGONAL>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::RANDOM>>;
};

/***************************************************************************
**SEED INITIAL CONDITIONS**************************************************/

//! Randomly places squares throughout the system.
/*!
 * The square seeds initial condition takes four parameters. The first
 * parameter determines the number of seeds which are generated.
 * The side lengths of the seed is chosen by multiplying the second  
 * parameter by #IC_SEED_SCALE_FACTOR times the system dimensions. The last
 * two parameters choose the initial values inside the seed and outside,
 * respectively.
 * 
 * - In 2 dimensions, each seed is defined by the points in the set:
 * \f[ 2\Delta x < \mathcal{S} B \hat{x} \cup 
 * 2\Delta y < \mathcal{S} B \hat{y} \f]
 * where \f$(\Delta x, \Delta y)\f$ is the randomly selected position around
 * which the seed is generated, \f$\mathcal{S}\f$ is a scaling factor of the 
 * system axis length and \f$\hat{x} = x_1 - x_0\f$ is the system axis 
 * length for the \f$x\f$ axis (equivalently for \f$\hat{y}\f$). These points
 * are assigned the value \f$C\f$, and all other points are assigned \f$D\f$.
 * 
 * - In 3 dimensions, the same recipe is followed.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		size_t n = static_cast<size_t>(init.data.gp[0]);

		len_type center_dims[D];
		double center_delta[D];
		double scaled_dims[D];
		for (iter_type i = 0; i < D; ++i)
		{
			center_dims[i] = static_cast<len_type>(dims[i] / 2);
			center_delta[i] = dims[i] / 2;
			scaled_dims[i] = dims[i] * IC_SEED_SIZE_FACTOR;
		}

		offsets = symphas::internal::RandomDeltas<D>(n, center_dims);
		lengths = symphas::internal::RandomOffsets<len_type, D>(n, scaled_dims);
		offsets.add_to_all(center_delta);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<len_type, D> lengths;		//!< Manages a list of random sizes.
};


//! Randomly places squares which are randomly sized throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSSQUARE> by adding 
 * randomness to the value generated inside each seed.
 * 
 * The value of each seed is chosen from the uniform distribution:
 * \f[\mathcal{X} \sim \mathcal{D}(C - 2C\delta, C) \f]
 * where \f$\delta\f$ is a scaling factor equal to #IC_SEED_RND_FACTOR. 
 * The mean is then \f$C - C\delta\f$. Thus, in this algorithm, the
 * third parameter determines the maximum value for the uniform
 * distribution. All other points are still assigned \f$D\f$.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE>;
	
	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims)
	{
		size_t n = static_cast<size_t>(init.data.gp[0]);
		scalar_t in_value = (symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT)
			? init.data.gp[3]
			: init.data.gp[2]);
		scalar_t value_rng[] = { in_value * (1 - IC_SEED_RND_FACTOR) };

		values = symphas::internal::RandomOffsets<scalar_t, 1>(n, value_rng, IC_SEED_RND_FACTOR);
	}

	scalar_t operator[](iter_type) const override;

	symphas::internal::RandomOffsets<scalar_t, 1> values;			//!< Manages a list of random values.
};


//! Randomly places squares throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSSQUARE>
 * by adding randomness to the seed sizes.
 * 
 * Instead of using a constant \f$\mathcal{S}\f$ (From InitialConditionsAlg<D, Inside::SEEDSSQUARE>), 
 * the algorithm chooses a random scaling factors for each side
 * derived from the corresponding system axis length. Without loss of 
 * generality, let \f$\hat{\mathcal{S}}_x\f$ be the scaling factor for the 
 * \f$x\f$-axis; it is taken from the uniform distribution:
 * \f[ \mathcal{X} \sim \mathcal{D}(\mathcal{S} - \delta_x, 
 * \mathcal{S} + \delta_x)\f].
 * Here, \f$\delta_x = \eta\hat{x}\f$, with value of \f$\eta\f$ given by
 * #IC_SEED_SIZE_FACTOR. 
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Randomly places squares which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA> and
 * InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


//! Randomly places squares throughout the system.
/*!
 * See InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA>.
 *
 * Each seed is chosen as a perfect square, rather than a rectangle.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Randomly places squares which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB> and
 * InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


//! Randomly places squares throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSSQUARE>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * The circular seeds initial condition takes five parameters. The first
 * parameter determines the number of seeds which are generated.
 * The axial radii of the seed is chosen by multiplying the second
 * parameter by #IC_SEED_SCALE_FACTOR times the system dimensions.
 * The next two parameters choose the initial values inside the seed, such
 * that the overall third parameter is the minimum of the uniform distribution
 * and the overall fourth parameter is the maximum of the uniform distribution.
 * Finally, the fifth parameter chooses the value of the outside.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		gen{ std::random_device{}() },
		seed_value_dis{
			(init.data.gp[2] < init.data.gp[3]) ? init.data.gp[2] : init.data.gp[3],
			(init.data.gp[3] > init.data.gp[2]) ? init.data.gp[3] : init.data.gp[2] } {}


	mutable std::mt19937 gen;
	mutable std::uniform_real_distribution<scalar_t> seed_value_dis;

};


//! Randomly places squares which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>;

	scalar_t operator[](iter_type) const override;

	InitialConditionsAlg(
		symphas::init_entry_type const& tdata,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) : 
		parent_type(tdata, vdata, dims),
		gen{ std::random_device{}() },
		seed_value_dis{
			(tdata.data.gp[2] < tdata.data.gp[3]) ? tdata.data.gp[2] : tdata.data.gp[3],
			(tdata.data.gp[3] > tdata.data.gp[2]) ? tdata.data.gp[3] : tdata.data.gp[2] } {}


	mutable std::mt19937 gen;
	mutable std::uniform_real_distribution<scalar_t> seed_value_dis;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
};



//! Randomly places squares throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * Instead of using a constant \f$\mathcal{S}\f$ (From InitialConditionsAlg<D, Inside::SEEDSSQUARE>),
 * the algorithm chooses a random scaling factors for each radius
 * derived from the corresponding system axis length. Without loss of
 * generality, let \f$\hat{\mathcal{S}}_x\f$ be the scaling factor for the
 * \f$x\f$-axis; it is taken from the uniform distribution:
 * \f[ \mathcal{X} \sim \mathcal{D}(\mathcal{S} - \delta_x,
 * \mathcal{S} + \delta_x)\f].
 * Here, \f$\delta_x = \eta\hat{x}\f$, with value of \f$\eta\f$ given by
 * #IC_SEED_SIZE_FACTOR.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places squares which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARC, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARA, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};




//! Randomly places circles throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * Each seed is chosen as a perfect circle, rather than an ellipse.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARC, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM, InsideTag::VARB, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARB, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};







template<>
struct symphas::internal::tags_for_init_value<Inside::SEEDSSQUARE>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::VARA, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARB>,
		init_tag_values_list<InsideTag::VARB, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARC, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA, InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARA, InsideTag::VARC, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARB, InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARB, InsideTag::VARC, InsideTag::RANDOM>>;
};


//! Randomly places circles throughout the system.
/*!
 * The circular seeds initial condition takes four parameters. The first
 * parameter determines the number of seeds which are generated.
 * The axial radii of the seed is chosen by multiplying the second
 * parameter by #IC_SEED_SCALE_FACTOR times the system dimensions. The last
 * two parameters choose the initial values inside the seed and outside,
 * respectively.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSCIRCLE> by adding 
 * randomness to the value generated inside each seed.
 * 
 * The value of each seed is chosen from the uniform distribution:
 * \f[\mathcal{X} \sim \mathcal{D}(C - 2C\delta, C) \f]
 * where \f$\delta\f$ is a scaling factor equal to #IC_SEED_RND_FACTOR. 
 * The mean is then \f$C - C\delta\f$. Thus, in this algorithm, the
 * third parameter determines the maximum value for the uniform
 * distribution. All other points are still assigned \f$D\f$.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSCIRCLE>
 * by adding randomness to the seed sizes.
 * 
 * Instead of using a constant \f$\mathcal{S}\f$ (From InitialConditionsAlg<D, Inside::SEEDSCIRCLE>), 
 * the algorithm chooses a random scaling factors for each radius
 * derived from the corresponding system axis length. Without loss of 
 * generality, let \f$\hat{\mathcal{S}}_x\f$ be the scaling factor for the 
 * \f$x\f$-axis; it is taken from the uniform distribution:
 * \f[ \mathcal{X} \sim \mathcal{D}(\mathcal{S} - \delta_x, 
 * \mathcal{S} + \delta_x)\f].
 * Here, \f$\delta_x = \eta\hat{x}\f$, with value of \f$\eta\f$ given by
 * #IC_SEED_SIZE_FACTOR.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA> and
 * InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};

//! Randomly places circles throughout the system.
/*!
 * See InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA>.
 *
 * Each seed is chosen as a perfect circle, rather than an ellipse.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB> and
 * InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};



//! Randomly places circles throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSCIRCLE>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * The circular seeds initial condition takes five parameters. The first
 * parameter determines the number of seeds which are generated.
 * The axial radii of the seed is chosen by multiplying the second
 * parameter by #IC_SEED_SCALE_FACTOR times the system dimensions.
 * The next two parameters choose the initial values inside the seed, such
 * that the overall third parameter is the minimum of the uniform distribution
 * and the overall fourth parameter is the maximum of the uniform distribution.
 * Finally, the fifth parameter chooses the value of the outside.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSSQUARE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
};



//! Randomly places circles throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * Instead of using a constant \f$\mathcal{S}\f$ (From InitialConditionsAlg<D, Inside::SEEDSCIRCLE>), 
 * the algorithm chooses a random scaling factors for each radius
 * derived from the corresponding system axis length. Without loss of 
 * generality, let \f$\hat{\mathcal{S}}_x\f$ be the scaling factor for the 
 * \f$x\f$-axis; it is taken from the uniform distribution:
 * \f[ \mathcal{X} \sim \mathcal{D}(\mathcal{S} - \delta_x, 
 * \mathcal{S} + \delta_x)\f].
 * Here, \f$\delta_x = \eta\hat{x}\f$, with value of \f$\eta\f$ given by
 * #IC_SEED_SIZE_FACTOR.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARC, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARA, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARA, InsideTag::RANDOM>;
	using parent_type::parent_type;
};




//! Randomly places circles throughout the system.
/*!
 * This algorithm modifies InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB>
 * by populating the inside of the seed with random values from the
 * uniform distribution.
 *
 * Each seed is chosen as a perfect circle, rather than an ellipse.
 *
 * If the invert tag is also provided, then the outside will have random
 * values from the uniform distribution, and the inside of the seeds will
 * have constant values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


//! Randomly places circles which are randomly sized throughout the system.
/*!
 * This algorithm combines InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC> and
 * InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM>.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM>;
	using parent_type::parent_type;
	scalar_t operator[](iter_type) const override;
};


template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::VARC, InsideTag::RANDOM> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARC, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::RANDOM, InsideTag::VARB, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARB, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};
template<size_t D>
struct InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::SEEDSCIRCLE, InsideTag::VARC, InsideTag::VARB, InsideTag::RANDOM>;
	using parent_type::parent_type;
};




template<>
struct symphas::internal::tags_for_init_value<Inside::SEEDSCIRCLE>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::VARA, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARB>,
		init_tag_values_list<InsideTag::VARB, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARC, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARA, InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARA, InsideTag::VARC, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::VARB, InsideTag::VARC>,
		init_tag_values_list<InsideTag::VARB, InsideTag::VARC, InsideTag::RANDOM>>;
};



//! Generates values for a Voronoi Diagram. 
/*!
 * Given the number of Voronoi points as the first initial condition parameter, 
 * this algorithm will generate a 
 * corresponding number of random values in the uniform range between the
 * second and third initial conditions parameters for each crystal. 
 * That is, given parameters `N`, `A` and `B`, (specified in that order), 
 * `N` is the number of crystals generated, `A` is the lowest value a crystal
 * can be populated with, and `B` is the greatest value a crystal can
 * be populated with.
 * 
 * The routine simply finds the center of the closest Voronoi crystal for each
 * point that must be generated, and fill it with the value of that cystal.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ static_cast<size_t>(init.data.gp[0]) }
	{
		scalar_t center_delta[D];
		for (iter_type i = 0; i < D; ++i)
		{
			center_delta[i] = dims[i] / 2;
		}

		offsets = symphas::internal::RandomDeltas<D>(N, dims, 0.5);
		offsets.add_to_all(center_delta);
		offsets.sort();

		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		values = symphas::internal::RandomOffsets<scalar_t, 1>(N, value_rng, (init.data.gp[2] - init.data.gp[1]) / 2);
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

#define BUBBLE_R_RATIO (1.0 / 16)


//! Fills the system with equally sized semi-overlapping bubbles. 
/*!
 * Fill the bubble with N points, where N is the first parameter of the initial conditions.
 * The second and third parameters are the interior and exterior values of the bubbles.
 * The fourth parameter is the radius of the bubbles, computed based on the area of placing the all
 * the bubbles. Typically, setting this value to 1 will mean that the desired number of bubbles
 * cannot be placed.
 * The fifth parameter is the packing ratio, which represents how tightly together the bubbles
 * should be placed. When this value is zero, there is no overlap between bubbles. 
 * For packing ratios close to zero,
 * the final number of bubbles may not be reached, but the algorithm will attempt to 
 * place as many bubbles as possible. 
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ size_t(init.data.gp[0]) }, R{ symphas::internal::compute_bubble_R(vdata, init.data.gp[3], N) }
	{
		// The amount of overlap controls how close bubbles should be to each other.
		// For high packing ratios, the overlap will be positive because bubble can be generated
		// on top of each other. For small packing ratios, the overlap will be ngative to give
		// each bubble more space.
		scalar_t overlap = symphas::internal::compute_bubble_overlap(init.data.gp[3], R, init.data.gp[4]);

		// The change in overlap varies from .2R to 1R as overlap itself goes from .2R to -2R, after which it is always R
		double overlap_eps = symphas::internal::compute_bubble_overlap_range(init.data.gp[3], R, init.data.gp[4]);

		// when overlap is negative, then the episolon needs to be added in order for the offsets to be computed correctly.
		symphas::internal::RandomOffsets<scalar_t, 1> overlaps(N, overlap, overlap_eps);
		symphas::internal::RandomOffsets<iter_type, D> start(1, dims);

		offsets = symphas::internal::get_bubble_positions<D>(N, R, overlap + overlap_eps, overlaps, dims, start.get_offset(0));
		offsets.sort();

		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		values = symphas::internal::RandomOffsets<scalar_t, 1>(N, value_rng, (init.data.gp[2] - init.data.gp[1]) / 2);
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	double R;
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SIN> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	using parent_type::init;
	using parent_type::dims;
	using parent_type::vdata;

	scalar_t operator[](iter_type n) const override
	{
		return symphas::internal::get_function_value(n, &symphas::math::sin<scalar_t>, dims, vdata, init.data.gp);
	}
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::SIN, InsideTag::VARA> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	using parent_type::init;
	using parent_type::dims;
	using parent_type::vdata;

	scalar_t operator[](iter_type n) const override
	{
		return symphas::internal::get_function_value_A(n, &symphas::math::sin<scalar_t>, dims, vdata, init.data.gp);
	}
};


template<>
struct symphas::internal::tags_for_init_value<Inside::SIN>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>>;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::COS> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	using parent_type::init;
	using parent_type::dims;
	using parent_type::vdata;

	scalar_t operator[](iter_type n) const override
	{
		return symphas::internal::get_function_value(n, &symphas::math::cos<scalar_t>, dims, vdata, init.data.gp);
	}
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::COS, InsideTag::VARA> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	using parent_type::init;
	using parent_type::dims;
	using parent_type::vdata;

	scalar_t operator[](iter_type n) const override
	{
		return symphas::internal::get_function_value_A(n, &symphas::math::cos<scalar_t>, dims, vdata, init.data.gp);
	}
};

template<>
struct symphas::internal::tags_for_init_value<Inside::COS>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>>;
};

//! Generates values for a Voronoi Diagram with periodic boundaries. 
/*!
 * Given the number of Voronoi points as the first initial condition parameter,
 * this algorithm will generate a
 * corresponding number of random values in the uniform range between the
 * second and third initial conditions parameters for each crystal.
 * Each crystal will be periodic around the system.
 * 
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::VORONOI>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI>;
	using parent_type::parent_type;

	scalar_t operator[](iter_type) const override;
};

//! Generates values for a Voronoi Diagram which, once set, are always the same. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED> 
	: InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims), voronoi{ get_voronoi(init, vdata, dims) } {}


	scalar_t operator[](iter_type n) const override
	{
		return voronoi[n];
	}

	auto get_voronoi(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims)
	{
		static InitialConditionsAlg<D, Inside::VORONOI> v(init, vdata, dims);
		return v;
	}

	InitialConditionsAlg<D, Inside::VORONOI> voronoi;
	
};

//! Generates values for a Voronoi Diagram which, once set, are always the same. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same values.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::VARA>
	: InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims), voronoi{ get_voronoi(init, vdata, dims) } {}


	scalar_t operator[](iter_type n) const override
	{
		return voronoi[n];
	}

	auto get_voronoi(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims)
	{
		static InitialConditionsAlg<D, Inside::VORONOI> v(init, vdata, dims);
		return v;
	}

	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::VARA> voronoi;

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::VARA, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::VARA>;
	using parent_type::parent_type;
};

//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ static_cast<size_t>(init.data.gp[0]) }
	{
		static auto offsets0 = get_offsets(dims);
		offsets = offsets0;

		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		values = symphas::internal::RandomOffsets<scalar_t, 1>(N, value_rng, (init.data.gp[1] - init.data.gp[2]) / 2);
	}

	auto get_offsets(len_type const* dims)
	{
		scalar_t center_delta[D];
		for (iter_type i = 0; i < D; ++i)
		{
			center_delta[i] = dims[i] / 2;
		}

		symphas::internal::RandomDeltas<D> offsets0(N, dims, 0.5);
		offsets0.add_to_all(center_delta);
		offsets0.sort();
		return offsets0;
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::RANDOM, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>;
	using parent_type::parent_type;
};

//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ static_cast<size_t>(init.data.gp[0]) }
	{
		static auto offsets0 = get_offsets(dims);
		offsets = offsets0;

		values = get_values(init);
	}

	auto get_offsets(len_type const* dims) const
	{
		scalar_t center_delta[D];
		for (iter_type i = 0; i < D; ++i)
		{
			center_delta[i] = dims[i] / 2;
		}

		symphas::internal::RandomDeltas<D> offsets0(N, dims, 0.5);
		offsets0.add_to_all(center_delta);
		offsets0.sort();
		return offsets0;
	}

	auto get_values(symphas::init_entry_type const& init) const
	{
		static int I = 0;
		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		static symphas::internal::RandomOffsets<scalar_t, 1> values0(N, value_rng, (init.data.gp[1] - init.data.gp[2]) / 2);
		
		double value_rng0[] = { init.data.gp[1] };
		symphas::internal::RandomOffsets<scalar_t, 1> values1(N, value_rng0, 0);
		values1.set_offset(I, values0.get_offset(I));

		I += 1;
		return values1;
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::VARB, InsideTag::RANDOM, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};



template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARA, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARB, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};


//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ static_cast<size_t>(init.data.gp[0]) }
	{
		static auto offsets0 = get_offsets(dims);
		offsets = offsets0;

		values = get_values(init);
	}

	auto get_offsets(len_type const* dims) const
	{
		scalar_t center_delta[D];
		for (iter_type i = 0; i < D; ++i)
		{
			center_delta[i] = dims[i] / 2;
		}

		symphas::internal::RandomDeltas<D> offsets0(N, dims, 0.5);
		offsets0.add_to_all(center_delta);
		offsets0.sort();
		return offsets0;
	}

	auto get_values(symphas::init_entry_type const& init) const
	{
		static int I = 0;

		double value_rng0[] = { init.data.gp[1] };
		symphas::internal::RandomOffsets<scalar_t, 1> values1(N, value_rng0, 0);
		values1.set_offset(I, init.data.gp[2]);

		I += 1;
		return values1;
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>;
	using parent_type::parent_type;
};



//! Generates values for a Voronoi Diagram where the crystal always has the same shape. 
/*!
 * If generating multiple fields, but the voronoi crystal should have the exact same
 * parameters, then this variation will always return the same crystal, but the values
 * of the crystal will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::VORONOI, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};



template<>
struct symphas::internal::tags_for_init_value<Inside::VORONOI>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::FIXEDSEED>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>>;
};





//! Generates values for a bubble arrangement with periodic boundaries. 
/*!
 * Given the number of bubbles as the first initial condition parameter,
 * this algorithm will generate a
 * corresponding number of random values in the uniform range between the
 * second and third initial conditions parameters for each arrangement.
 * Each arrangement will be periodic around the system.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::BUBBLE>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE>;
	using parent_type::parent_type;

	scalar_t operator[](iter_type) const override;
};

//! Fills the system with equally sized semi-overlapping bubbles. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same values.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED>
	: InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims), bubble{ get_bubble(init, vdata, dims) } {}


	scalar_t operator[](iter_type n) const override
	{
		return bubble[n];
	}

	auto get_bubble(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims)
	{
		static InitialConditionsAlg<D, Inside::BUBBLE> v(init, vdata, dims);
		return v;
	}

	InitialConditionsAlg<D, Inside::BUBBLE> bubble;

};

//! Fills the system with equally sized semi-overlapping bubbles. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same values.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::VARA>
	: InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims), bubble{ get_bubble(init, vdata, dims) } {}


	scalar_t operator[](iter_type n) const override
	{
		return bubble[n];
	}

	auto get_bubble(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims)
	{
		static InitialConditionsAlg<D, Inside::BUBBLE> v(init, vdata, dims);
		return v;
	}

	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::VARA> bubble;

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::VARA, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::VARA>;
	using parent_type::parent_type;
};

//! Fills the system with equally sized semi-overlapping bubbles. 
/*!
 * Fill the bubble with N points, where N is the first parameter of the initial conditions.
 * The thid parameter of the initial conditions is the packing ratio. For high packing ratios,
 * the value may not be reached, but the algorithm will attempt to place as many bubbles
 * as possible.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ size_t(init.data.gp[0]) }, R{ symphas::internal::compute_bubble_R(vdata, init.data.gp[3], N) }
	{
		// The amount of overlap controls how close bubbles should be to each other.
		// For high packing ratios, the overlap will be positive because bubble can be generated
		// on top of each other. For small packing ratios, the overlap will be ngative to give
		// each bubble more space.
		scalar_t overlap = symphas::internal::compute_bubble_overlap(init.data.gp[3], R, init.data.gp[4]);

		// The change in overlap varies from .2R to 1R as overlap itself goes from .2R to -2R, after which it is always R
		double overlap_eps = symphas::internal::compute_bubble_overlap_range(init.data.gp[3], R, init.data.gp[4]);

		// when overlap is negative, then the episolon needs to be added in order for the offsets to be computed correctly.
		static symphas::internal::RandomOffsets<scalar_t, 1> overlaps(N, overlap, overlap_eps);
		static symphas::internal::RandomOffsets<iter_type, D> start(1, dims);
		static auto offsets0 = symphas::internal::get_bubble_positions<D>(N, R, overlap + overlap_eps, overlaps, dims, start.get_offset(0));

		offsets = offsets0;

		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		values = symphas::internal::RandomOffsets<scalar_t, 1>(N, value_rng, (init.data.gp[2] - init.data.gp[1]) / 2);;
	}

	scalar_t operator[](iter_type) const override;

	size_t N;													//!< Number of regions.
	double R;
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::RANDOM, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>;
	using parent_type::parent_type;
};


//! Fills the system with equally sized semi-overlapping bubbles. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same arrangement, but the values
 * of the arrangement will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>;
	using parent_type::parent_type;
};

//! Generates values for a bubble arrangement where the arrangement always has the same shape. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same arrangement, but the values
 * of the arrangement will be different.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ size_t(init.data.gp[0]) }, R{ symphas::internal::compute_bubble_R(vdata, init.data.gp[3], N) }
	{
		// The amount of overlap controls how close bubbles should be to each other.
		// For high packing ratios, the overlap will be positive because bubble can be generated
		// on top of each other. For small packing ratios, the overlap will be ngative to give
		// each bubble more space.
		scalar_t overlap = symphas::internal::compute_bubble_overlap(init.data.gp[3], R, init.data.gp[4]);

		// The change in overlap varies from .2R to 1R as overlap itself goes from .2R to -2R, after which it is always R
		double overlap_eps = symphas::internal::compute_bubble_overlap_range(init.data.gp[3], R, init.data.gp[4]);

		// when overlap is negative, then the episolon needs to be added in order for the offsets to be computed correctly.
		static symphas::internal::RandomOffsets<scalar_t, 1> overlaps(N, overlap, overlap_eps);
		static symphas::internal::RandomOffsets<iter_type, D> start(1, dims);
		static auto offsets0 = symphas::internal::get_bubble_positions<D>(N, R, overlap + overlap_eps, overlaps, dims, start.get_offset(0));

		offsets = offsets0;

		double value_rng[] = { (init.data.gp[2] + init.data.gp[3]) / 2 };
		values = get_values(init);
	}

	scalar_t operator[](iter_type) const override;


	auto get_values(symphas::init_entry_type const& init) const
	{
		static int I = 0;
		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		static auto values0 = symphas::internal::RandomOffsets<scalar_t, 1>(N, value_rng, (init.data.gp[2] - init.data.gp[1]) / 2);

		double value_rng0[] = { init.data.gp[1] };
		symphas::internal::RandomOffsets<scalar_t, 1> values1(N, value_rng0, 0);
		values1.set_offset(I, values0.get_offset(I));

		I += 1;
		return values1;
	}

	size_t N;													//!< Number of regions.
	double R;
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::VARB, InsideTag::RANDOM, InsideTag::FIXEDSEED> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

//! Generates values for a bubble arrangement where the arrangement always has the same shape. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same arrangement, but the values
 * of the arrangement will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};



template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARA, InsideTag::VARB> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARB, InsideTag::VARA> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>;
	using parent_type::parent_type;
};


//! Generates values for a bubble arrangement where the arrangement always has the same shape. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same arrangement, but the values
 * of the arrangement will be different.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC> :
	InitialConditionsData<D>
{

	using parent_type = InitialConditionsData<D>;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		N{ size_t(init.data.gp[0]) }, R{ symphas::internal::compute_bubble_R(vdata, init.data.gp[3], N) },
		select{ get_select() }
	{
		// The amount of overlap controls how close bubbles should be to each other.
		// For high packing ratios, the overlap will be positive because bubble can be generated
		// on top of each other. For small packing ratios, the overlap will be ngative to give
		// each bubble more space.
		scalar_t overlap = symphas::internal::compute_bubble_overlap(init.data.gp[3], R, init.data.gp[4]);

		// The change in overlap varies from .2R to 1R as overlap itself goes from .2R to -2R, after which it is always R
		double overlap_eps = symphas::internal::compute_bubble_overlap_range(init.data.gp[3], R, init.data.gp[4]);

		// when overlap is negative, then the episolon needs to be added in order for the offsets to be computed correctly.
		static symphas::internal::RandomOffsets<scalar_t, 1> overlaps(N, overlap, overlap_eps);
		static symphas::internal::RandomOffsets<iter_type, D> start(1, dims);
		static auto offsets0 = symphas::internal::get_bubble_positions<D>(N, R, overlap + overlap_eps, overlaps, dims, start.get_offset(0));

		offsets = offsets0;

		double value_rng[] = { (init.data.gp[1] + init.data.gp[2]) / 2 };
		values = get_values(init);
	}

	scalar_t operator[](iter_type) const override;

	auto get_select()
	{
		static int I = 0;
		I %= N;
		return I++;
	}

	auto get_values(symphas::init_entry_type const& init) const
	{
		double value_rng0[] = { init.data.gp[1] };
		symphas::internal::RandomOffsets<scalar_t, 1> values1(N, value_rng0, 0);
		values1.set_offset(select, init.data.gp[2]);
		return values1;
	}

	size_t N;													//!< Number of regions.
	double R;
	symphas::internal::RandomDeltas<D> offsets;					//!< Manages a list of random offsets.
	symphas::internal::RandomOffsets<scalar_t, 1> values;		//!< Manages a list of random values.
	int select;

};

template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::RANDOM, InsideTag::FIXEDSEED, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>;
	using parent_type::parent_type;
};



//! Generates values for a bubble arrangement where the arrangement always has the same shape. 
/*!
 * If generating multiple fields, but the bubble arrangement should have the exact same
 * parameters, then this variation will always return the same arrangement, but the values
 * of the arrangement will be different.
 *
 * This variation will apply periodic boundaries when determining the closest
 * point.
 */
template<size_t D>
struct InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC> :
	InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>
{
	using parent_type = InitialConditionsAlg<D, Inside::BUBBLE, InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>;
	using parent_type::parent_type;
	using parent_type::N;
	using parent_type::offsets;
	using parent_type::values;

	scalar_t operator[](iter_type) const override;

};



template<>
struct symphas::internal::tags_for_init_value<Inside::BUBBLE>
{
	using type = symphas::lib::types_list<
		init_tag_values_list<InsideTag::VARA>,
		init_tag_values_list<InsideTag::FIXEDSEED>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARB>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARC>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARB>,
		init_tag_values_list<InsideTag::FIXEDSEED, InsideTag::RANDOM, InsideTag::VARA, InsideTag::VARC>>;
};



/***************************************************************************
**RANDOM INITIAL CONDITIONS************************************************/


template<size_t D>
struct InitialConditionsAlg<D, Inside::GAUSSIAN> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::init;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		gen{ std::random_device{}() },
		dis{ init.data.gp[0], init.data.gp[1] }
	{}

	scalar_t operator[](iter_type) const override
	{
		return dis(gen);
	}

	mutable std::mt19937 gen;
	mutable std::normal_distribution<scalar_t> dis;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::UNIFORM> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::init;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		gen{ std::random_device{}() },
		dis{ 0, 1 }, a{ init.data.gp[0] }, b{ init.data.gp[1] }
	{}

	scalar_t operator[](iter_type) const override
	{
		return (a + dis(gen) * (b - a));
	}

	mutable std::mt19937 gen;
	mutable std::uniform_real_distribution<scalar_t> dis;
	double a;
	double b;
};




/***************************************************************************
**OTHER INITIAL CONDITIONS*************************************************/

template<size_t D>
struct InitialConditionsAlg<D, Inside::CONSTANT> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;
	using parent_type::init;

	scalar_t operator[](iter_type) const override
	{
		return init.data.gp[0];
	}

	scalar_t value;
};

template<size_t D>
struct InitialConditionsAlg<D, Inside::CAPPED> :
	InitialConditionsAlg<D, Inside::UNIFORM>
{
	using parent_type = InitialConditionsAlg<D, Inside::UNIFORM>;
	using parent_type::parent_type;
	using parent_type::init;

	InitialConditionsAlg(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* dims
	) :
		parent_type(init, vdata, dims),
		mean{ (init.data.gp[0] + init.data.gp[1]) / 2 }
	{}

	scalar_t operator[](iter_type n) const override
	{
		return (parent_type::operator[](n) < mean)
			? init.data.gp[0]
			: init.data.gp[1];
	}

	scalar_t mean;

};


template<size_t D>
struct InitialConditionsAlg<D, Inside::NONE, InsideTag::NONE> :
	InitialConditionsData<D>
{
	using parent_type = InitialConditionsData<D>;
	using parent_type::parent_type;

	scalar_t operator[](iter_type n) const override
	{
		return 0;
	}

	operator bool() const override
	{
		return false;
	}
};






template<typename T, size_t D>
struct InitialConditions
{
	InitialConditions(
		symphas::init_data_type const& tdata,
		symphas::interval_data_type const& vdata,
		len_type const* origin
	) :
		data{ get_data(tdata, vdata, origin) }
	{}

	~InitialConditions()
	{
		for (auto& [key, entry] : data)
		{
			delete entry;
		}
	}

	auto begin(Axis ax, T* values, grid::region_interval<D> const& interval) const
	{
		return symphas::internal::ic_iterator<T, D>(ax, values, *(data.at(ax)), interval);
	}

	auto end(Axis ax, T* values, grid::region_interval<D> const& interval) const
	{
		return symphas::internal::ic_iterator<T, D>(ax, values, *(data.at(ax)), interval, grid::length<D>(interval));
	}

	bool initialize(T* values, grid::region_interval<D> const& interval, size_t id = 0) const
	{
		bool initialized = true;
		if (data.find(Axis::NONE) != data.end())
		{
			initialized = initialize(Axis::NONE, values, interval, id);
		}

		for (auto const& [key, entry] : data)
		{
			if (key != Axis::NONE)
			{
				initialized = initialized && initialize(key, values, interval, id);
			}
		}
		return initialized;
	}

	bool initialize(Axis ax, T* values, grid::region_interval<D> const& interval, size_t id = 0) const
	{
		// If the initial conditions are not set, then nothing is done.
		if (data.at(ax)->init.in == Inside::NONE)
		{
			if (symphas::internal::tag_bit_compare(data.at(ax)->init.intag, InsideTag::NONE))
			{
				return true;
			}
			else if (symphas::internal::tag_bit_compare(data.at(ax)->init.intag, InsideTag::DEFAULT))
			{
				data.at(ax)->init.f_init->initialize(values, interval);
				return true;
			}
			else
			{
				return false;
			}
		}
#ifdef USING_IO
		else if (data.at(ax)->init.in == Inside::FILE || data.at(ax)->init.in == Inside::CHECKPOINT)
		{
			symphas::io::read_info rinfo{ 
				data.at(ax)->init.file.get_index(), id,
				data.at(ax)->init.file.get_name(), data.at(ax)->init.in == Inside::CHECKPOINT };

			symphas::grid_info ginfo = symphas::io::read_header(rinfo);
			for (iter_type n = 0; n < D; ++n)
			{
				auto file_dims = ginfo.get_dims();
				if (file_dims[n] != interval.dims[n])
				{
					Axis ax0 = symphas::index_to_axis(n);
					fprintf(SYMPHAS_WARN, "%c-dimension read from file header (%d) is inconsistent with "
						"system dimension that is being initialized (%d)\n",
						(ax0 == Axis::X) ? 'x' : (ax0 == Axis::Y) ? 'y' : (ax0 == Axis::Z) ? 'z' : '0',
						file_dims[n], interval.dims[n]);

					return false;
				}
			}

			int read_index = symphas::io::read_grid(values, rinfo);
			if (data.at(ax)->init.file.get_index() != read_index)
			{
				fprintf(SYMPHAS_ERR,
					"system initialization requires the loaded datafile to contain the given index '%d'\n",
					data.at(ax)->init.file.get_index());
				return false;
			}

		}
#endif
		else if (data.at(ax)->init.in == Inside::EXPRESSION)
		{
			return match_init_expr<D>(data.at(ax)->init.expr_data.get_name(), ax, values,
				interval, data.at(ax)->vdata, data.at(ax)->init.expr_data.get_coeff(), data.at(ax)->init.expr_data.get_num_coeff());
		}
		else
		{
			if (*this)
			{
				auto it = symphas::data_iterator_region(values, interval);
				std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
					std::execution::par_unseq,
#endif
					begin(ax, values, interval), end(ax, values, interval), it);

			}
			else
			{
				return false;
			}
		}

		return true;
	}

	operator bool() const
	{
		for (auto const& [key, entry] : data)
		{
			if (!entry) return false;
		}
		return true;
	}


protected:
	
	std::map<Axis, InitialConditionsData<D>*> data;

	auto get_data(
		symphas::init_data_type const& tdata,
		symphas::interval_data_type const& vdata,
		len_type const* origin)
	{
		std::map<Axis, InitialConditionsData<D>*> data0;
		for (auto& [key, entry] : tdata)
		{
			InitialConditionsData<D>* init;
			if (entry.in == Inside::NONE)
			{
				init = next_ic(entry, vdata, origin,
					symphas::internal::init_values_list<>{});
			}
			else
			{
				init = next_ic(entry, vdata, origin,
					symphas::internal::init_values_list<
					ALL_INSIDE_GEN_VALUES
					>{});
			}
			data0[key] = init;
		}
		return data0;
	}

	template<Inside in>
	using icd_type = symphas::internal::make_new_ic_delegate<D, in>;

	InitialConditionsData<D>* next_ic(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* origin,
		symphas::internal::init_values_list<>)
	{
		return new InitialConditionsAlg<D, Inside::NONE, InsideTag::NONE>(init, vdata, grid::interior_dimensions(vdata));
	}

	template<Inside in0, Inside... ins>
	InitialConditionsData<D>* next_ic(
		symphas::init_entry_type const& init,
		symphas::interval_data_type const& vdata,
		len_type const* origin,
		symphas::internal::init_values_list<in0, ins...>)
	{
		InitialConditionsData<D>* ic = icd_type<in0>{}(init, vdata, origin);
		if (!ic)
		{
			return next_ic(init, vdata, origin, symphas::internal::init_values_list<ins...>{});
		}
		else
		{
			return ic;
		}
	}

};

//! @}
