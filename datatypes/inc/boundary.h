
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
 * PURPOSE: Defines the boundary class used for the finite difference
 * implementation of solvers, i.e. when central differencing is used for
 * all derivative approximations.
 *
 * ***************************************************************************
 */


#pragma once


#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4250)
#endif


#include <cassert>
#include <map>
#include <cmath>
#include <random>
#include <cstring>
#include <algorithm>

#include "spslib.h"
#include "grid.h"

//! \cond

#define NUM_BOUNDARY_CONSTANTS 5	//!< Maximum number of parameters for a boundary.

#define GAUSSIAN_COEF 1.0			//!< The default value used in the Gaussian condition.
#define GAUSSIAN_STD 1.0			//!< The default value of the standard deviation used in the Gaussian condition.
#define GAUSSIAN_B 1.0				//!< The default value of the offset of the mean of the Gaussian condition.
#define SIN_COEF_A 1.0				//!< The first default coefficient for the sine condition.
#define SIN_COEF_B 1.0				//!< The second default coefficient for the cosine condition.
#define SIN_K 1.0					//!< The default frequency of the cosine condition.
#define LINEAR_COEF 1.0				//!< The default slope of the linear condition.
#define LINEAR_B 0.0				//!< The default value of the offset of the linear condition.
#define CONSTANT_A 0.0				//!< The default value of the constant value of the constant condition.
#define UNIFORM_A 0					//!< The default starting value of the uniform distribution condition.
#define UNIFORM_B 1					//!< The default end value of the uniform distribution condition.
#define TIME_SCALE 1.0				//!< The default multiplier on the time in time-dependent conditions.
//! \endcond



/*!
 * \defgroup Boundaries Boundary Objects
 * @{
 */


//! Value representing the boundary type.
/*!
 * The boundary type determines the update strategy used to update the
 * boundaries.
 */
enum class BoundaryType
{
	DEFAULT,			//!< Default boundaries, usually Dirichlet and Neumann conditions.
	OPEN,				//!< Boundary representing an open or free space on the other side.
	NONE,				//!< Not a boundary type, for a "no-op" boundary or for error checking.
	PERIODIC,			//!< Periodic boundaries, which tile the system.
	PERIODIC0,			//!< Same as PERIODIC when used as a boundary parameter. In generating boundaries, represents sides that are only opposite.
	PERIODIC3XZ,		//!< Same as PERIODIC when used as a boundary parameter. Periodic in x and z (but not y).
	PERIODIC3XY,		//!< Same as PERIODIC when used as a boundary parameter. Periodic in x and y (but not z).
	PERIODIC3YZ			//!< Same as PERIODIC when used as a boundary parameter. Periodic in y and z (but not x).
};


//! Value representing a modifier of the boundary.
/*!
 * Up to two modifiers may be applied to a single boundary, but only specific 
 * combinations are valid. Typically, any boundary which is generated from
 * a function will also support a time dependence, i.e. it is combinable
 * with the `BoundaryTag::TIME` value.
 * 
 * The modifiers define the algorithm used to update the boundaries.
 * 
 * Refer to grid::BoundaryDefaultTagged and associated specializations
 * for the supported boundary implementations.
 */
enum class BoundaryTag
{
	GAUSSIAN,		//!< Uses a Gaussian distribution to generate random values at each point.
	RANDOM,			//!< Uses a uniform distribution to generate random values at each point.
	TRIG,			//!< Evaluates a trigonometric function of the position on the boundary.
	CONSTANT,		//!< A constant value is set all points on the boundary.
	LINEAR,			//!< Values distributed according to a linear equation.
	TIME,			//!< Modifies an algorithm to include time dependence.
	NONE			//!< No modifier, values are never updated.
};



namespace symphas
{

	//! Manages boundary data along one boundary.
	/*! 
	 * Contain the boundary information for configuring the phase field 
	 * problem.
	 */
	struct b_element_type
	{
		BoundaryType type;				//!< The boundary type.
		BoundaryTag tag[2];				//!< The tag associated with the main type.
		double *params;					//!< Parameters controlling the boundary behaviour.
		int argc;						//!< The number of parameters used.

		//! Construct a new boundary data element.
		/*!
		 * Using the boundary type and a list of the boundary tags, create
		 * a new boundary data element. The given parameters are stored.
		 * 
		 * \param type The type of boundary condition.
		 * \param tags The modifiers of the boundary condition.
		 * \param params A list of parameters used in the boundary algorithm.
		 * \param argc The number of parameters.
		 */
		b_element_type(BoundaryType type, std::initializer_list<BoundaryTag> tags, const double* params, int argc);

		//! Construct a new boundary data element.
		/*!
		 * Using the boundary type and a list of the boundary tags, create
		 * a new boundary data element. The given parameters are stored.
		 *
		 * \param type The type of boundary condition.
		 * \param tags An array of 2 modifiers of the boundary condition.
		 * \param params A list of parameters used in the boundary algorithm.
		 * \param argc The number of parameters.
		 */
		b_element_type(BoundaryType type, const BoundaryTag tag[2], const double* params, int argc) :
			b_element_type(type, { tag[0], tag[1] }, params, argc) {}

		//! Construct a new boundary data element.
		/*!
		 * Using the boundary type and a list of the boundary tags, create
		 * a new boundary data element. No parameters are provided.
		 *
		 * \param type The type of boundary condition.
		 * \param tags The modifiers of the boundary condition.
		 */
		b_element_type(BoundaryType type, std::initializer_list<BoundaryTag> tags) :
			b_element_type(type, tags, nullptr, 0) {}

		//! Construct a new boundary data element.
		/*!
		 * Using the boundary type and a list of the boundary tags, create
		 * a new boundary data element. No parameters are provided.
		 *
		 * \param type The type of boundary condition.
		 * \param tags An array of 2 modifiers of the boundary condition.
		 */
		b_element_type(BoundaryType type, BoundaryTag tag[2]) :
			b_element_type(type, { tag[0], tag[1] }) {}

		//! Construct a new boundary data element.
		/*!
		 * Using the boundary type and the parameters, create
		 * a new boundary data element. No modifiers are provided.
		 *
		 * \param type The type of boundary condition.
		 * \param params A list of parameters used in the boundary algorithm.
		 * \param argc The number of parameters.
		 */
		b_element_type(BoundaryType type, const double* params, int argc) :
			b_element_type(type, { BoundaryTag::NONE, BoundaryTag::NONE }, params, argc) {}

		//! Construct a new boundary data element.
		/*!
		 * Only the boundary type is provided to generate this boundary data.
		 *
		 * \param type The type of boundary condition.
		 */
		b_element_type(BoundaryType type) :
			b_element_type(type, { BoundaryTag::NONE, BoundaryTag::NONE }, nullptr, 0) {}

		//! Construct a new boundary data element.
		/*!
		 * No information about the boundary is provided, so a BoundaryType::DEFAULT
		 * boundary value is assumed.
		 */
		b_element_type() : b_element_type(BoundaryType::DEFAULT, { BoundaryTag::NONE, BoundaryTag::NONE }, nullptr, 0) {}
		
		b_element_type(b_element_type const& other);
		b_element_type(b_element_type&& other) noexcept;

		b_element_type& operator=(b_element_type other);

		//! Set the boundary value parameter at the given index.
		/*!
		 * Set the boundary value parameter at the given index. If the
		 * index is beyond the range of the current list, the list will be
		 * extended, and any unset parameters will be set to 0.
		 */
		void set_parameter(double value, iter_type n);

		friend void swap(b_element_type& first, b_element_type& second)
		{
			using std::swap;
			swap(first.type, second.type);
			swap(first.tag, second.tag);
			swap(first.params, second.params);
			swap(first.argc, first.argc);
		}

		~b_element_type();
	};

	//! The container representing boundary information about a system.
	/*!
	 * The data container which represents information about boundary 
	 * conditions in the context of the problem parameters is represented 
	 * by a map that relates the enum Side to a b_element_type.
	 */
	using b_data_type = std::map<Side, symphas::b_element_type>;



	//! From the given string, get the corresponding initial condition.
	/*!
	 * The initial condition is of type ::Inside. The given string is a key
	 * representing one of the possible initial conditions, which it will
	 * return.
	 */
	BoundaryType boundary_from_str(const char* type);

	//! From the given initial condition, find its corresponding key.
	/*!
	 * The initial condition is of type ::Inside. One of the strings that
	 * represents the initial condition is returned.
	 */
	const char* str_from_boundary(BoundaryType in);


	//! From the given string, get the corresponding initial condition.
	/*!
	 * The initial condition is of type ::Inside. The given string is a key
	 * representing one of the possible initial conditions, which it will
	 * return.
	 */
	BoundaryTag boundary_tag_from_str(const char* type);

	//! From the given initial condition, find its corresponding key.
	/*!
	 * The initial condition is of type ::Inside. One of the strings that
	 * represents the initial condition is returned.
	 */
	const char* str_from_boundary_tag(BoundaryTag tag);

}

namespace grid
{


	//! A boundary object, used as a feature of a Grid.
	/*!
	 * A boundary is a feature of a phase field system. Functionally, the boundary
	 * object supplements the grid by maintaining the interval information for the
	 * grid boundary which it corresponds to, and allowing the boundary points
	 * to be updated.
	 */
	template<typename T, size_t D>
	struct Boundary
	{
		//! Generate a copy of the specialized boundary instance.
		virtual Boundary<T, D>* new_copy() const = 0;

		//! Get the parameters used to construct this boundary.
		/*!
		 * Returns a new parameters object that can be used to construct a new boundary
		 * object with identical properties.
		 */
		virtual symphas::b_element_type get_parameters() const = 0;

		virtual ~Boundary() {}
	};

	//! A boundary object, used as a feature of a Grid.
	/*!
	 * A boundary is a feature of a phase field system. Functionally, the boundary
	 * object supplements the grid by maintaining the interval information for the
	 * grid boundary which it corresponds to, and allowing the boundary points
	 * to be updated.
	 */
	template<typename T, size_t D, BoundaryType type>
	struct BoundaryApplied;

	// **************************************************************************************


	 //! Periodic boundary.
	 /*!
	  * Implements functions that allow the grid boundary to be updated to
	  * reflect periodic boundaries.
	  */
	template<typename T, size_t D>
	struct BoundaryApplied<T, D, BoundaryType::PERIODIC> : Boundary<T, D>
	{
		using Boundary<T, D>::Boundary;
		Boundary<T, D>* new_copy() const
		{
			return new BoundaryApplied<T, D, BoundaryType::PERIODIC>(*this);
		}

		//! Get the parameters used to construct this boundary.
		/*!
		 * Returns a new parameters object that can be used to construct a new boundary
		 * object with identical properties.
		 */
		symphas::b_element_type get_parameters() const
		{
			return symphas::b_element_type(BoundaryType::PERIODIC);
		}

	};



	// **************************************************************************************

	//! A default point boundary for a 1 dimensional system.
	/*!
	 * A default boundary has 5 constants that are used in constructing the
	 * boundary value. An update function is supplied
	 */
	template<typename T>
	struct BoundaryApplied<T, 0, BoundaryType::DEFAULT> : Boundary<T, 0>
	{
		double
			_A, //!< First parameter.
			_B, //!< Second parameter.
			_C, //!< Third parameter.
			_D, //!< Fourth parameter.
			_E; //!< Fifth parameter.

		double v;	//!< The position of the boundary.

		//! Create a new default boundary with the given parameters.
		BoundaryApplied(double* constants, size_t num_constants) : Boundary<T, 0>()
		{
			double* values[] = { &_A, &_B, &_C, &_D, &_E };
			constexpr size_t len_values = sizeof(values) / sizeof(double*);
			for (iter_type i = 0; i < std::min(len_values, num_constants); ++i)
			{
				*values[i] = constants[i];
			}
			for (size_t i = std::min(len_values, num_constants); i < len_values; ++i)
			{
				*values[i] = 0;
			}
		}

		//! Create a new default boundary with no parameters.
		BoundaryApplied() : BoundaryApplied(nullptr, 0) {}

		//! Initialize the interval data of the boundary.
		void init(double);

		void fill_constants(double* constants)
		{
			constants[0] = _A;
			constants[1] = _B;
			constants[2] = _C;
			constants[3] = _D;
			constants[4] = _E;
		}

		//! Update the given value which is at a position on the boundary.
		/*!
		 * Given the position of the point, an update algorithm is chosen based
		 * on the default boundary specialization to update the given boundary
		 * point.
		 *
		 * \param val Reference to a point on the boundary which is updated.
		 * \param time The solution time at which the boundary is updated.
		 */
		virtual void update(T& val, axis_coord_t, axis_coord_t, double time) const = 0;

		template<typename vector_type = T, typename T0 = typename vector_element_type<vector_type>::type>
		void update(multi_value<1, T0> val, axis_coord_t, axis_coord_t, double time) const
		{
			any_vector_t<T0, 1> vector = val;
			this->update(vector, 0, 0, time);
			val = vector;
		}

		void update(carry_value<T> val, axis_coord_t, axis_coord_t, double time) const
		{
			if (!val.clear)
			{
				update(*val.value, 0, 0, time);
			}
		}

		virtual ~BoundaryApplied() {}

	};

	//! A default edge boundary, for a 2 dimensional system.
	template<typename T>
	struct BoundaryApplied<T, 1, BoundaryType::DEFAULT> : Boundary<T, 1>
	{
		double
			_A, //!< First parameter.
			_B, //!< Second parameter.
			_C, //!< Third parameter.
			_D, //!< Fourth parameter.
			_E; //!< Fifth parameter.

		double v[2];	//!< The interval of the boundary.
		double h;		//!< The spatial separation along the boundary.


		//! Create a new default boundary with the given parameters.
		BoundaryApplied(double* constants, size_t num_constants) : Boundary<T, 1>()
		{
			double* values[] = { &_A, &_B, &_C, &_D, &_E };
			constexpr size_t len_values = sizeof(values) / sizeof(double*);
			for (iter_type i = 0; i < std::min(len_values, num_constants); ++i)
			{
				*values[i] = constants[i];
			}
			for (size_t i = std::min(len_values, num_constants); i < len_values; ++i)
			{
				*values[i] = 0;
			}
		}

		//! Create a new default boundary with no parameters.
		BoundaryApplied() : BoundaryApplied(nullptr, 0) {}

		//! Initialize the interval data of the boundary.
		void init(double*, double);

		void fill_constants(double* constants)
		{
			constants[0] = _A;
			constants[1] = _B;
			constants[2] = _C;
			constants[3] = _D;
			constants[4] = _E;
		}

		//! Update the given value which is at a position on the boundary.
		/*!
		 * Given the position of the point, an update algorithm is chosen based
		 * on the default boundary specialization to update the given boundary
		 * point.
		 *
		 * \param val Reference to a point on the boundary which is updated.
		 * \param x The \f$x\f$ position of the point in the system.
		 * \param time The solution time at which the boundary is updated.
		 */
		virtual void update(T& val, axis_coord_t x, axis_coord_t, double time) const = 0;

		template<typename vector_type = T, typename T0 = typename vector_element_type<vector_type>::type>
		void update(multi_value<2, T0> val, axis_coord_t x, axis_coord_t, double time) const
		{
			any_vector_t<T0, 2> vector = val;
			this->update(vector, x, 0, time);
			val = vector;
		}

		void update(carry_value<T> val, axis_coord_t x, axis_coord_t, double time) const
		{
			if (!val.clear)
			{
				update(*val.value, x, 0, time);
			}
		}

		virtual ~BoundaryApplied() {}

	};

	//! A default plane boundary, for a 3 dimensional system.
	template<typename T>
	struct BoundaryApplied<T, 2, BoundaryType::DEFAULT> : Boundary<T, 2>
	{
		double
			_A, //!< First parameter.
			_B, //!< Second parameter.
			_C, //!< Third parameter.
			_D, //!< Fourth parameter.
			_E; //!< Fifth parameter.

		double v[4];	//!< The interval of the boundary.
		double h[2];	//!< The axial spacings of the boundary.


		//! Create a new default boundary with the given parameters.
		BoundaryApplied(double* constants, size_t num_constants) : Boundary<T, 2>()
		{
			double* values[] = { &_A, &_B, &_C, &_D, &_E };
			constexpr size_t len_values = sizeof(values) / sizeof(double*);
			for (iter_type i = 0; i < std::min(len_values, num_constants); ++i)
			{
				*values[i] = constants[i];
			}
			for (size_t i = std::min(len_values, num_constants); i < len_values; ++i)
			{
				*values[i] = 0;
			}
		}

		//! Create a new default boundary with no parameters.
		BoundaryApplied() : BoundaryApplied(nullptr, 0) {}

		//! Initialize the interval data of the boundary.
		void init(double*, double*);

		//! Update the given value which is at a position on the boundary.
		/*!
		 * Given the position of the point, an update algorithm is chosen based
		 * on the default boundary specialization to update the given boundary
		 * point.
		 *
		 * \param val Reference to a point on the boundary which is updated.
		 * \param x The \f$x\f$ position of the point in the system.
		 * \param y The \f$y\f$ position of the point in the system.
		 * \param time The solution time at which the boundary is updated.
		 */
		virtual void update(T& val, axis_coord_t x, axis_coord_t y, double time) const = 0;

		template<typename vector_type = T, typename T0 = typename vector_element_type<vector_type>::type>
		void update(multi_value<3, T0> val, axis_coord_t x, axis_coord_t y, double time) const
		{
			any_vector_t<T0, 3> vector = val;
			this->update(vector, x, y, time);
			val = vector;
		}

		void update(carry_value<T> val, axis_coord_t x, axis_coord_t y, double time) const
		{
			if (!val.clear)
			{
				update(*val.value, x, y, time);
			}
		}

		virtual ~BoundaryApplied() {}

	};


	//! A default boundary which has tags applied to it.
	/*!
	 * A default boundary which has tags applied to it. Specializations
	 * define the algorithm used.
	 */
	template<typename T, size_t D, BoundaryTag...>
	struct BoundaryDefaultTagged;

	namespace internal
	{
		template<typename T, size_t D, BoundaryTag... tags>
		symphas::b_element_type get_parameters(BoundaryDefaultTagged<T, D, tags...> const* boundary)
		{
			double* params = new double[NUM_BOUNDARY_CONSTANTS] {0};
			params[0] = boundary->_A;
			params[1] = boundary->_B;
			params[2] = boundary->_C;
			params[3] = boundary->_D;
			params[4] = boundary->_E;

			symphas::b_element_type bdata(BoundaryType::DEFAULT, { tags... }, params, NUM_BOUNDARY_CONSTANTS);

			delete[] params;

			return bdata;
		}
	}

	//! Functional object aiding in the implementation of specialized boundary.
	template<template<typename, size_t> typename B, typename T, size_t D>
	struct BoundarySupportClass
	{
		//! Generate a copy of the specialized boundary instance.
		/*!
		 * Implemented by each specialization of a default boundary in order
		 * to return a copy with the correct parameters.
		 */
		Boundary<T, D>* new_copy() const
		{
			return new B<T, D>(*cast());
		}

		void set_defaults()
		{
			(*static_cast<B<T, D>*>(this)).set_defaults();
		}

		symphas::b_element_type get_parameters() const
		{
			return grid::internal::get_parameters(cast());
		}

		const B<T, D>* cast() const
		{
			return static_cast<B<T, D> const*>(this);
		}
	};


	
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::RANDOM>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::GAUSSIAN>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::LINEAR>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::TRIG>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::TRIG, BoundaryTag::TIME>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::LINEAR, BoundaryTag::TIME>;
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::CONSTANT>;


	/*
	 * Alias for the boundary algorith implementations.
	 */
	template<typename T, size_t D>
	using BoundaryRandom = BoundaryDefaultTagged<T, D, BoundaryTag::RANDOM>;
	template<typename T, size_t D>
	using BoundaryGaussian = BoundaryDefaultTagged<T, D, BoundaryTag::GAUSSIAN>;
	template<typename T, size_t D>
	using BoundaryLinear = BoundaryDefaultTagged<T, D, BoundaryTag::LINEAR>;
	template<typename T, size_t D>
	using BoundaryTrig = BoundaryDefaultTagged<T, D, BoundaryTag::TRIG>;
	template<typename T, size_t D>
	using BoundaryTrigTime = BoundaryDefaultTagged<T, D, BoundaryTag::TRIG, BoundaryTag::TIME>;
	template<typename T, size_t D>
	using BoundaryLinearTime = BoundaryDefaultTagged<T, D, BoundaryTag::LINEAR, BoundaryTag::TIME>;
	template<typename T, size_t D>
	using BoundaryConstant = BoundaryDefaultTagged<T, D, BoundaryTag::CONSTANT>;



	/* default boundary specializations
	 */

	 //! A specialization of the default boundary applying trigonmetric functions.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::TRIG> : 
		BoundaryApplied<T, D, BoundaryType::DEFAULT>, 
		BoundarySupportClass<BoundaryTrig, T, D>
	{
		using parent_type = BoundaryApplied<T, D, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryTrig, T, D>;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;

		using support_type::new_copy;


		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}

		//! Initialize a boundary value to the result of the Gaussian function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A \sin(C(x + y)) + B\cos(C(x + y))\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t y, double) const
		{
			val = (_A * std::sin(_C * (x + y)) + _B * std::cos(_C * (x + y)))
				* symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = SIN_COEF_A;
			_B = SIN_COEF_B;
			_C = SIN_K;
		}

		Boundary<T, D>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! A specialization of the default boundary applying trigonmetric functions.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::TRIG, BoundaryTag::TIME> : 
		BoundaryApplied<T, D, BoundaryType::DEFAULT>, 
		BoundarySupportClass<BoundaryTrigTime, T, D>
	{
		using parent_type = BoundaryApplied<T, D, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryTrigTime, T, D>;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;
		using parent_type::_D;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the Gaussian function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A \sin(C(x + y) + Dt) + B\cos(C(x + y) + Dt)\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t y, double time) const
		{
			val = (_A * std::sin(_C * (x + y) + _D * time)
				+ _B * std::cos(_C * (x + y) + _D * time))
				* symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = SIN_COEF_A;
			_B = SIN_COEF_B;
			_C = SIN_K;
			_D = TIME_SCALE;
		}

		Boundary<T, D>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};


	//! A specialization of the default boundary that applies the Gaussian function.
	template<typename T>
	struct BoundaryDefaultTagged<T, 0, BoundaryTag::GAUSSIAN> :
		BoundaryApplied<T, 0, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryGaussian, T, 0>
	{
		using parent_type = BoundaryApplied<T, 0, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryGaussian, T, 0>;
		using parent_type::v;
		using parent_type::_A;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the Gaussian function.
		/*!
		 * Uses only parameter \f$A\f$.
		 * For each value, returns \f$A\f$.
		 */
		void update(T& val, axis_coord_t, axis_coord_t, double) const
		{
			val = _A * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = GAUSSIAN_COEF;
		}

		Boundary<T, 0>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! Specialization based on BoundaryGaussian.
	template<typename T>
	struct BoundaryDefaultTagged<T, 1, BoundaryTag::GAUSSIAN> :
		BoundaryApplied<T, 1, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryGaussian, T, 1>
	{
		using parent_type = BoundaryApplied<T, 1, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryGaussian, T, 1>;
		using parent_type::v;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the Gaussian function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A \exp((-(x-B)^2 / 2C^2))\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t, double) const
		{
			double xx = x - (v[1] + v[0]) / 2.0;
			val = _A * std::exp(-(std::pow(xx - _B, 2.0)) / std::pow(_C, 2.0) / 2.0)
				* symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = GAUSSIAN_COEF;
			_B = GAUSSIAN_B;
			_C = GAUSSIAN_STD;
		}

		Boundary<T, 1>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! See BoundaryGaussian.
	template<typename T>
	struct BoundaryDefaultTagged<T, 2, BoundaryTag::GAUSSIAN> :
		BoundaryApplied<T, 2, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryGaussian, T, 2>
	{
		using parent_type = BoundaryApplied<T, 2, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryGaussian, T, 2>;
		using parent_type::v;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;
		using parent_type::_D;

		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the Gaussian function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A \exp((-(x-B)^2 - (y-C)^2) / 2D^2)\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t y, double) const
		{
			double xx = x - (v[1] + v[0]) / 2.0;
			double yy = y - (v[3] + v[3]) / 2.0;
			val = _A * std::exp(
				-(std::pow(x - _B, 2.0) + std::pow(y - _C, 2.0)) / std::pow(_D, 2.0) / 2.0)
				* symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = GAUSSIAN_COEF;
			_B = GAUSSIAN_B;
			_C = GAUSSIAN_B;
			_D = GAUSSIAN_STD;
		}

		Boundary<T, 2>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! A specialization of the default boundary applying a constant.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::CONSTANT> :
		BoundaryApplied<T, D, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryConstant, T, D>
	{
		using parent_type = BoundaryApplied<T, D, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryConstant, T, D>;
		using parent_type::_A;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to a constant value.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A \f].
		 */
		void update(T& val, axis_coord_t, axis_coord_t, double) const
		{
			val = _A * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_B;
		}

		Boundary<T, D>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! Specialization based on BoundaryLinear.
	template<typename T>
	struct BoundaryDefaultTagged<T, 0, BoundaryTag::LINEAR> :
		BoundaryApplied<T, 0, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinear, T, 0>
	{
		using parent_type = BoundaryApplied<T, 0, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinear, T, 0>;
		using parent_type::_A;
		using parent_type::_B;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = A\f].
		 */
		void update(T& val, axis_coord_t, axis_coord_t, double) const
		{
			val = _A * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
		}

		Boundary<T, 0>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! Specialization based on BoundaryLinear.
	template<typename T>
	struct BoundaryDefaultTagged<T, 1, BoundaryTag::LINEAR> :
		BoundaryApplied<T, 1, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinear, T, 1>
	{
		using parent_type = BoundaryApplied<T, 1, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinear, T, 1>;
		using parent_type::_A;
		using parent_type::_B;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = Ax + B\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t, double) const
		{
			val = (_A * x + _B) * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
			_B = LINEAR_B;
		}

		Boundary<T, 1>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! Specialization based on BoundaryLinear.
	template<typename T>
	struct BoundaryDefaultTagged<T, 2, BoundaryTag::LINEAR> :
		BoundaryApplied<T, 2, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinear, T, 2>
	{
		using parent_type = BoundaryApplied<T, 2, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinear, T, 2>;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}

		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = Ax + By + C\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t y, double) const
		{
			val = (_A * x + _B * y + _C) * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
			_B = LINEAR_COEF;
			_C = LINEAR_B;
		}

		Boundary<T, 2>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};



	//! A specialization of the default boundary applying a linear function.
	template<typename T>
	struct BoundaryDefaultTagged<T, 0, BoundaryTag::LINEAR, BoundaryTag::TIME> :
		BoundaryApplied<T, 0, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinearTime, T, 0>
	{
		using parent_type = BoundaryApplied<T, 0, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinearTime, T, 0>;
		using parent_type::_A;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = At\f].
		 */
		void update(T& val, axis_coord_t, axis_coord_t, double time) const
		{
			val = _A * time * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
		}

		Boundary<T, 0>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};


	//! A specialization of the default boundary applying a linear function.
	template<typename T>
	struct BoundaryDefaultTagged<T, 1, BoundaryTag::LINEAR, BoundaryTag::TIME> :
		BoundaryApplied<T, 1, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinearTime, T, 1>
	{
		using parent_type = BoundaryApplied<T, 1, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinearTime, T, 1>;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}

		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = Ct(Ax + B)\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t, double time) const
		{
			val = _C * time * (_A * x + _B) * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
			_B = LINEAR_B;
			_C = TIME_SCALE;
		}

		Boundary<T, 1>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! A specialization of the default boundary applying a linear function.
	template<typename T>
	struct BoundaryDefaultTagged<T, 2, BoundaryTag::LINEAR, BoundaryTag::TIME> :
		BoundaryApplied<T, 2, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryLinearTime, T, 2>
	{
		using parent_type = BoundaryApplied<T, 2, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryLinearTime, T, 2>;
		using parent_type::_A;
		using parent_type::_B;
		using parent_type::_C;
		using parent_type::_D;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type() {}


		//! Initialize a boundary value to the result of the linear function.
		/*!
		 * Updates the given value \f$V\f$ using the equation:
		 * \f[V = Dt(Ax + By + C)\f].
		 */
		void update(T& val, axis_coord_t x, axis_coord_t y, double time) const
		{
			val = _D * time * (_A * x + _B * y + _C) * symphas::lib::get_identity<T>();
		}

		void set_defaults()
		{
			_A = LINEAR_COEF;
			_B = LINEAR_COEF;
			_C = LINEAR_B;
			_D = TIME_SCALE;
		}

		Boundary<T, 2>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}
	};

	//! A specialization of the default boundary applying a random value.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::RANDOM> :
		BoundaryApplied<T, D, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryRandom, T, D>
	{
		using parent_type = BoundaryApplied<T, D, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryRandom, T, D>;
		using parent_type::_A;
		using parent_type::_B;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged(double* constants, size_t num_constants) :
			parent_type(constants, num_constants), support_type(),
			mt_eng{ std::random_device{}() }, prob_dist(_A, _B), th(-symphas::PI, symphas::PI) {}

		//! Initialize a boundary value a random value.
		/*!
		 * Updates the given value \f$V\f$ using the uniform distribution:
		 * \f[\mathcal{X} \sim \mathcal{U}(A, B)\f].
		 */
		void update(T& val, axis_coord_t, axis_coord_t, double) const
		{
			update(val);
		}

		void set_defaults()
		{
			_A = UNIFORM_A;
			_B = UNIFORM_B;
		}

		Boundary<T, D>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}

	protected:
		mutable std::mt19937 mt_eng;
		mutable std::uniform_real_distribution<scalar_t> prob_dist;
		mutable std::uniform_real_distribution<scalar_t> th;


		template<typename S>
		void update(S& val) const
		{
			val = prob_dist(mt_eng);
		}

		void update(vector_t<1>& val) const
		{
			val = { prob_dist(mt_eng) };
		}

		void update(vector_t<2>& val) const
		{

			scalar_t
				a = th(mt_eng),
				r = prob_dist(mt_eng);
			val = { r * cos(a), r * sin(a) };
		}

		void update(vector_t<3>& val) const
		{

			scalar_t
				a = th(mt_eng),
				b = th(mt_eng),
				r = prob_dist(mt_eng);
			val = { r * sin(a) * cos(b), r * sin(a) * sin(b), r * cos(a) };

		}

	};


	//! A specialization of the default boundary applying a random value.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::NONE>;
	template<typename T, size_t D>
	using BoundaryNone = BoundaryDefaultTagged<T, D, BoundaryTag::NONE>;

	//! A specialization of the default boundary applying a random value.
	template<typename T, size_t D>
	struct BoundaryDefaultTagged<T, D, BoundaryTag::NONE> :
		BoundaryApplied<T, D, BoundaryType::DEFAULT>,
		BoundarySupportClass<BoundaryNone, T, D>
	{
		using parent_type = BoundaryApplied<T, D, BoundaryType::DEFAULT>;
		using support_type = BoundarySupportClass<BoundaryNone, T, D>;

		//! Create a new boundary with the given parameters.
		BoundaryDefaultTagged() : parent_type(nullptr, 0), support_type() {}

		void update(T& val, axis_coord_t, axis_coord_t, double) const {}
		void set_defaults() {}

		Boundary<T, D>* new_copy() const
		{
			return support_type::new_copy();
		}

		symphas::b_element_type get_parameters() const
		{
			return support_type::get_parameters();
		}

	};

	/* implementation of boundary iniitialization functions for each specialization
	 */


	template<typename T>
	void BoundaryApplied<T, 0, BoundaryType::DEFAULT>::init(double _v)
	{
		v = _v;
	}

	template<typename T>
	void BoundaryApplied<T, 1, BoundaryType::DEFAULT>::init(double* _v, double _h)
	{
		std::copy(_v, _v + 2, v);
		h = _h;
	}

	template<typename T>
	void BoundaryApplied<T, 2, BoundaryType::DEFAULT>::init(double* _v, double* _h)
	{
		std::copy(_v, _v + 4, v);
		std::copy(_h, _h + 2, h);
	}
}

//! @}

/*
 *
 *
 *
 *
 *
 * implementation
 *
 */

namespace symphas::internal
{

	extern std::map<const char*, BoundaryType, any_case_comparator> boundary_key_map;
	extern std::map<const char*, BoundaryTag, any_case_comparator> boundary_tag_key_map;





	//! Creates a new boundary based on the boundary data.
	template<typename T, size_t D>
	grid::Boundary<T, D>* new_boundary(symphas::b_element_type const& data)
	{

		grid::Boundary<T, D>* b = nullptr;

		
		if (data.tag[1] == BoundaryTag::TIME)
		{
			switch (data.type)
			{
			case BoundaryType::DEFAULT:
			{
				switch (data.tag[0])
				{
				case BoundaryTag::TRIG:
				{
					b = new grid::BoundaryTrigTime<T, D>(data.params, data.argc);
					break;
				}
				case BoundaryTag::LINEAR:
				{
					b = new grid::BoundaryLinearTime<T, D>(data.params, data.argc);
					break;
				}
				default:
					fprintf(SYMPHAS_WARN, "initializing boundary encountered the wrong tag\n");
				}
				break;
			}
			default:
				fprintf(SYMPHAS_WARN, "boundary condition does not support time\n");
			}
		}
		else
		{
			switch (data.type)
			{
			case BoundaryType::PERIODIC:
			{
				b = new grid::BoundaryApplied<T, D, BoundaryType::PERIODIC>{};
				break;
			}
			case BoundaryType::DEFAULT:
			{
				switch (data.tag[0])
				{
				case BoundaryTag::GAUSSIAN:
				{
					b = new grid::BoundaryGaussian<T, D>(data.params, data.argc);
					break;
				}
				case BoundaryTag::TRIG:
				{
					b = new grid::BoundaryTrig<T, D>(data.params, data.argc);
					break;
				}
				case BoundaryTag::CONSTANT:
				{
					b = new grid::BoundaryConstant<T, D>(data.params, data.argc);
					break;
				}
				case BoundaryTag::LINEAR:
				{
					b = new grid::BoundaryLinear<T, D>(data.params, data.argc);
					break;
				}
				case BoundaryTag::RANDOM:
				{
					b = new grid::BoundaryRandom<T, D>(data.params, data.argc);
					break;
				}
				default:
					fprintf(SYMPHAS_WARN, "initializing boundary encountered the wrong tag\n");
				}
				break;
			}
			case BoundaryType::NONE:
				break;
			default:
				fprintf(SYMPHAS_WARN, "boundary condition is not supported\n");
			}
		}

		return b;
	}
}


namespace symphas
{



	//! Helper functions for setting up the boundaries.
	/*!
	 * Sets up a list of boundaries of the prescribed dimension `D+1`.
	 *
	 * A boundary is created and depending on the type, boundary values will
	 * also be initialized. When described here, the boundary values are
	 * denoted as capital letters, \f$A, B, C, D\f$ and \f$E\f$. The order of
	 * the letter corresponds with its position in the array of values in
	 * the boundary data. This is order is the same as parsed from the
	 * configuration.
	 *
	 * What follows is a description of how the boundaries are initialized and
	 * what parameters different boundary types should expect.
	 *
	 * _Periodic boundary types_ (`PERIODIC`):
	 * - Does not take values or modifiers.
	 *
	 * _Default boundary types_ (`DEFAULT`): There are a number of
	 * boundary tags applicable to this type, which simply specializes the type
	 * of default boundary that it is. The list of these and their descriptions:
	 * - `GAUSSIAN` (\ref BoundaryGaussian "spec")
	 *    - Takes 4 values for a 3d system or 3 for a 2d system.
	 *       - In 2D: \f$A \exp((-(x-B)^2 / 2C^2))\f$
	 *       - In 3D: \f$A \exp((-(x-B)^2 - (y-C)^2) / 2D^2)\f$
	 *       - In both cases, the peak is designed to be in the center.
	 *    - Does not support any modifiers.
	 * - `TRIG` (\ref BoundaryTrig "spec")
	 *    - Takes 3 values.
	 *       - In 2D: \f$A \sin(Cx) + B\cos(Cx)\f$
	 *       - In 3D: \f$A \sin(C(x + y)) + B\cos(C(x + y))\f$
	 *    - Supports `TIME` modifier, (\ref BoundaryTrigTime "spec")
	 * adding one constant.
	 *       - In 2D: \f$A \sin(Cx + Dt) + B\cos(Cx + Dt)\f$
	 *       - In 3D: \f$A \sin(C(x + y) + Dt) + B\cos(C(x + y) + Dt)\f$
	 * - `CONSTANT` (\ref BoundaryConstant "spec")
	 *    - Takes 1 value. In any dimension, it returns this value.
	 *    - Does not support any modifiers.
	 * - `LINEAR` (\ref BoundaryLinear<T, 1> "spec")
	 *    - Takes 2 values in 2D and 3 in 3D.
	 *       - In 2D: \f$Ax + B\f$
	 *       - In 3D: \f$Ax + By + C\f$
	 *    - Supports `TIME` modifier (\ref BoundaryLinearTime "spec"),
	 * in which case it always takes 3 values.
	 *       - In 2D: \f$Ax + B + Ct\f$
	 *       - In 3D: \f$A(x + y) + B + Ct\f$
	 * - `RANDOM` (\ref BoundaryRandom "spec")
	 *    - Takes 2 values.
	 *       - In any dimension, the values are generated from the uniform
	 * distribution in the interval \f$[A, B]\f$.
	 *    - Does not support any modifiers.
	 *
	 *
	 *
	 * \param data The boundary data used to initialize the boundary.
	 *
	 * \tparam D The dimension of the boundaries, _not the system_. The dimension
	 * of the boundaries is always one minus the system dimension.
	 */
	template<size_t D>
	struct boundary_setup
	{
	protected:
		//! Typedef for referring to the default boundary of dimension `D`.
		template<typename T>
		using b_default_type = grid::BoundaryApplied<T, D, BoundaryType::DEFAULT>;

	public:
		//! Initialize the list of boundaries from the given boundary data.
		/*!
		 * The array of boundaries is initialized using the given boundary data.
		 * The number of pointers to boundaries in the list must be equal to
		 * twice the system dimension (which is the number of boundaries in a `D+1`-
		 * dimensional system).
		 *
		 * \param boundaries The list of boundaries which are initialized. The
		 * length of this array is twice the dimension `D`.
		 * \param bdata The boundary parameter data.
		 * \param intervals The interval parameter data.
		 *
		 * \tparam T The system value type.
		 */
		template<typename T>
		void operator()(grid::Boundary<T, D>* boundaries[], symphas::b_data_type const& bdata, symphas::interval_data_type const& intervals);
	};



	template<>
	template<typename T>
	void boundary_setup<0>::operator()(grid::Boundary<T, 0>* boundaries[], symphas::b_data_type const& bdata, symphas::interval_data_type const& intervals)
	{
		// there are only 2 boundaries for a 1d system, interval_index = 0 for left and interval_index = 1 for right
		// the boundaries correspond to only one point

		for (auto side : { Side::LEFT, Side::RIGHT })
		{
			if (bdata.find(side) != bdata.end())
			{
				boundaries[symphas::side_to_index(side)] = symphas::internal::new_boundary<T, 0>(bdata.at(side));

				if (bdata.at(side).type == BoundaryType::DEFAULT)
				{
					double interval = 0;
					if (side == Side::LEFT)
					{
						interval = DOMAIN_X0;
					}
					if (side == Side::RIGHT)
					{
						interval = DOMAIN_Xn;
					}

					static_cast<b_default_type<T>*>(boundaries[symphas::side_to_index(side)])->init(interval);
				}
			}
		}
	}

	template<>
	template<typename T>
	void boundary_setup<1>::operator()(grid::Boundary<T, 1>* boundaries[], symphas::b_data_type const& bdata, symphas::interval_data_type const& intervals)
	{
		// each boundary condition is an edge, and is specified by an interval [a, b]
		// the intervals are assigned to each boundary relative to the global coordinate system by
		// rotating the system clockwise so that the boundary for which the interval is constructed is 
		// at the top of the system, then the interval is selected to begin at the left


		for (auto side : { Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM })
		{
			if (bdata.find(side) != bdata.end())
			{
				boundaries[symphas::side_to_index(side)] = symphas::internal::new_boundary<T, 1>(bdata.at(side));

				if (bdata.at(side).type == BoundaryType::DEFAULT)
				{
					double interval[2]{ 0 };
					double h = 0;

					if (side == Side::TOP)
					{
						interval[0] = DOMAIN_X0;
						interval[1] = DOMAIN_Xn;
						h = INTERVAL_Xh;
					}
					else if (side == Side::BOTTOM)
					{
						interval[0] = DOMAIN_X0;
						interval[1] = DOMAIN_Xn;
						h = INTERVAL_Xh;
					}
					else if (side == Side::RIGHT)
					{
						interval[0] = DOMAIN_Y0;
						interval[1] = DOMAIN_Yn;
						h = INTERVAL_Yh;
					}
					else if (side == Side::LEFT)
					{
						interval[0] = DOMAIN_Y0;
						interval[1] = DOMAIN_Yn;
						h = INTERVAL_Yh;
					}

					static_cast<b_default_type<T>*>(boundaries[symphas::side_to_index(side)])->init(interval, h);
				}
			}
		}
	}

	template<>
	template<typename T>
	void boundary_setup<2>::operator()(grid::Boundary<T, 2>* boundaries[], symphas::b_data_type const& bdata, symphas::interval_data_type const& intervals)
	{
		// each boundary condition is a surface, and is specified by an interval [x0, x1], [y0, y1]
		// note that the order the axes are mentioned in this specification also indicates the local x and y
		// the first mentioned axis is the local x axis, the second menteiond axis is the local y axis
		// the interval order is determined by considering the orientation of the face after `rotating' the cube
		// to get to it, for the back face, the rotation is around the y axis, twice

		for (auto side : { Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM, Side::FRONT, Side::BACK })
		{
			if (bdata.find(side) != bdata.end())
			{
				boundaries[symphas::side_to_index(side)] = symphas::internal::new_boundary<T, 2>(bdata.at(side));

				if (bdata.at(side).type == BoundaryType::DEFAULT)
				{
					double interval[4]{ 0 };
					double h[2]{ 0 };

					if (side == Side::RIGHT)
					{
						interval[0] = DOMAIN_Z0;
						interval[1] = DOMAIN_Zn;
						interval[2] = DOMAIN_Y0;
						interval[3] = DOMAIN_Yn;

						h[0] = INTERVAL_Zh;
						h[1] = INTERVAL_Yh;
					}
					else if (side == Side::LEFT)
					{
						interval[0] = DOMAIN_Z0;
						interval[1] = DOMAIN_Zn;
						interval[2] = DOMAIN_Y0;
						interval[3] = DOMAIN_Yn;

						h[0] = INTERVAL_Zh;
						h[1] = INTERVAL_Yh;
					}
					else if (side == Side::TOP)
					{
						interval[0] = DOMAIN_X0;
						interval[1] = DOMAIN_Xn;
						interval[2] = DOMAIN_Z0;
						interval[3] = DOMAIN_Zn;

						h[0] = INTERVAL_Xh;
						h[1] = INTERVAL_Zh;
					}
					else if (side == Side::BOTTOM)
					{
						interval[0] = DOMAIN_X0;
						interval[1] = DOMAIN_Xn;
						interval[2] = DOMAIN_Z0;
						interval[3] = DOMAIN_Zn;

						h[0] = INTERVAL_Xh;
						h[1] = INTERVAL_Zh;
					}
					else if (side == Side::FRONT)
					{
						interval[0] = DOMAIN_X0;
						interval[1] = DOMAIN_Xn;
						interval[2] = DOMAIN_Y0;
						interval[3] = DOMAIN_Yn;

						h[0] = INTERVAL_Xh;
						h[1] = INTERVAL_Yh;
					}
					else if (side == Side::BACK)
					{
						interval[0] = DOMAIN_Xn;
						interval[1] = DOMAIN_X0;
						interval[2] = DOMAIN_Yn;
						interval[3] = DOMAIN_Y0;

						h[0] = INTERVAL_Xh;
						h[1] = INTERVAL_Yh;
					}


					static_cast<b_default_type<T>*>(boundaries[symphas::side_to_index(side)])->init(interval, h);
				}
			}
		}
	}

}





#undef GAUSSIAN_COEF
#undef GAUSSIAN_STD
#undef GAUSSIAN_B
#undef SIN_COEF_A
#undef SIN_COEF_B
#undef SIN_K
#undef LINEAR_COEF
#undef LINEAR_B
#undef CONSTANT_A
#undef UNIFORM_A
#undef UNIFORM_B
#undef TIME_SCALE









#ifdef _MSC_VER
#pragma warning(pop)
#endif



