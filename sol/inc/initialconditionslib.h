
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
 * PURPOSE: Contains functions and elements that are used by the initial
 * conditions in populating phase field values.
 *
 * ***************************************************************************
 */

#pragma once


#include <random>

#include "params.h"
#include "gridfunctions.h"


 //! The width of a distribution when using randomness.
#define IC_RND_STRENGTH params::init_rand_val
//! The percentage by which a square may be adjusted when randomly offset.
#define IC_SQUARE_RND_FACTOR 0.15
//! The percentage by which a circle may be adjusted when randomly offset.
#define IC_CIRCLE_RND_FACTOR IC_SQUARE_RND_FACTOR
//! The percentage by which spots may differ for randomly offset values.
#define IC_HEX_RND_FACTOR 0.07
//! The percentage by which spots may differ for randomly offset values.
#define IC_CUBIC_RND_FACTOR IC_HEX_RND_FACTOR
//! The percentage by which seeds may differ for randomly offset values.
#define IC_SEED_RND_FACTOR 0.20
//! The percentage by which voronoi cells may differ for randomly offset values.
#define IC_VORONOI_RND_FACTOR IC_SEED_RND_FACTOR


//! The default diameter of the seed is 5% of the given system dimension.
#define IC_SEED_SCALE_FACTOR 0.05


//! Returns the outer value to an initial condition.
#define IC_OUTER_VALUE (symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT) \
	? params::init_inside_val \
	: params::init_outside_val)
//! Returns the interior value to an initial condition.
#define IC_INNER_VALUE (symphas::internal::tag_bit_compare(init.intag, InsideTag::INVERT) \
	? params::init_outside_val \
	: params::init_inside_val)


//! The random size of seeds is determined by this factor.
#define IC_SEED_SIZE_FACTOR (IC_SEED_SCALE_FACTOR * 0.1)




//! All possible values for the initial conditions.
/*!
 * The specification of each is defined in InitialConditions.
 */
enum class Inside
{
	GAUSSIAN, 			//!< Uses a Gaussian distribution to randomly assign values.
	UNIFORM, 			//!< Uses a uniform distribution to randomly assign values.
	CAPPED, 			//!< Values are assigned to be either the minimum or maximum parameter.
	CONSTANT,			//!< All values are assigned to the same parameter.
	CIRCLE, 			//!< Value are assigned to be the shape of a circle.
	HEXAGONAL, 			//!< Values are assigned into circles arranged in a hexagonal pattern.
	CUBIC, 				//!< Values are assigned into circles arranged in a cubic pattern.
	SQUARE,				//!< Values are assigned to be the shape of a square.
	SEEDSSQUARE, 		//!< Values are put into randomly arranged squares.
	SEEDSCIRCLE,		//!< Values are put into randomly arranged circles.
	VORONOI,			//!< Generate values in a Voronoi diagram.
	BUBBLE,				//!< Generate random semi-overlapping circles.
	SPIRALHEX,			//!< Generate hexagonal pattern but one circle per field, spirally ordered.
	SIN,				//!< Generate values from the sin function.
	COS,				//!< Generate values from the sin function.
	FILE, 				//!< Values are read in from a file.
	CHECKPOINT,			//!< Values are read in from a checkpoint.
	EXPRESSION,			//!< An equation that is specified similar to the model definitions.
	NONE				//!< Represents no initial condition.
};

#define ALL_INSIDE_GEN_VALUES \
Inside::CAPPED,					\
Inside::CIRCLE,					\
Inside::SQUARE,					\
Inside::CONSTANT,				\
Inside::CUBIC,					\
Inside::HEXAGONAL,				\
Inside::GAUSSIAN,				\
Inside::UNIFORM,				\
Inside::SEEDSCIRCLE,			\
Inside::SEEDSSQUARE,			\
Inside::VORONOI,				\
Inside::BUBBLE,					\
Inside::SPIRALHEX,				\
Inside::SIN,					\
Inside::COS

//! Modifies the initial condition generation algorithm.
/*!
 * The values represented here can each be applied to modify the generation
 * algorithm, their value represents the bit which they assign to 1.
 */
enum class InsideTag
{
	DEFAULT,	//!< The default initial generation algorithm is chosen.
	RANDOM,		//!< The generation is modified to include some kind of randomness.
	FIXEDSEED,	//!< Randomness will not vary throughout the field; not supported by all algorithms.
	VARA,		//!< The A variation is chosen (different from the default variation).
	VARB,		//!< The B variation is chosen.
	VARC,		//!< The C variation is chosen.
	INVERT,		//!< The interior and outer values are switched in the generation.
	NONE		//!< Represents no tag.
};



template<size_t D>
struct InitialConditionsData;


namespace symphas::internal
{

	template<typename T>
	struct value_fill;

	template<>
	struct value_fill<scalar_t>
	{
		scalar_t operator()(Axis ax, scalar_t const& current, scalar_t value) const
		{
			return value;
		}
	};

	template<>
	struct value_fill<complex_t>
	{
		value_fill() :
			gen{ std::random_device{}() },
			th{ -symphas::PI, symphas::PI } {}

		complex_t operator()(Axis ax, complex_t const& current, scalar_t value) const
		{
			using namespace std;
			switch (ax)
			{
			case Axis::NONE:
			{
				scalar_t
					a = th(gen),
					r = std::abs(value);
				return { r * cos(a), r * sin(a) };
			}
			case Axis::X:
			{
				return { value, current.imag() };
			}
			case Axis::Y:
			{
				return { current.real(), value };
			}
			case Axis::T:
			{
				scalar_t m = abs(current);
				return { m * cos(PI * value), m * sin(PI * value) };
			}
			case Axis::R:
			{
				scalar_t r = abs(current) / value;
				return { current.real() / r, current.imag() / r };
			}
			default:
			{
				return current;
			}
			}
		}

		mutable std::mt19937 gen;
		mutable std::uniform_real_distribution<scalar_t> th;
	};


	template<>
	struct value_fill<vector_t<1>>
	{
		vector_t<1> operator()(Axis ax, vector_t<1> const& current, scalar_t value) const
		{
			return { value };
		}
	};

	template<>
	struct value_fill<vector_t<2>> : value_fill<complex_t>
	{
		using value_fill<complex_t>::value_fill;

		vector_t<2> operator()(Axis ax, vector_t<2> const& current, scalar_t value) const
		{
			auto c = value_fill<complex_t>::operator()(ax, complex_t(current[0], current[1]), value);
			return { c.real(), c.imag() };
		}
	};

	template<>
	struct value_fill<vector_t<3>> : value_fill<complex_t>
	{
		using value_fill<complex_t>::value_fill;

		vector_t<3> operator()(Axis ax, vector_t<3> const& current, scalar_t value) const
		{
			using namespace std;
			switch (ax)
			{
			case Axis::NONE:
			{
				scalar_t
					a = th(gen),
					b = th(gen),
					r = std::abs(value);
				return { r * sin(a) * cos(b), r * sin(a) * sin(b), r * cos(a) };
			}
			case Axis::X:
			{
				return { value, current[1], current[2] };
			}
			case Axis::Y:
			{
				return { current[0], value, current[2] };
			}
			case Axis::Z:
			{
				return { current[0], current[1], value };
			}
			case Axis::T:
			{
				scalar_t m = abs(current);
				scalar_t th = PI * value;
				scalar_t phi = acos(current[2] / m);
				return { m * sin(phi) * cos(th), m * sin(phi) * sin(th), m * cos(phi) };
			}
			case Axis::S:
			{
				scalar_t m = abs(current);
				scalar_t th = atan(current[1] / current[0]);
				scalar_t phi = PI * value;
				return { m * sin(phi) * cos(th), m * sin(phi) * sin(th), m * cos(phi) };
			}
			case Axis::R:
			{
				scalar_t r = abs(current) / value;
				return { current[0] / r, current[1] / r, current[2] / r };
			}
			default:
			{
				return current;
			}
			}
		}
	};


	//! Convert an InsideTag value to a size_t value.
	inline size_t tagtoul(InsideTag t)
	{
		return static_cast<size_t>(t);
	}

	extern std::map<const char*, Inside, symphas::internal::any_case_comparator> init_key_map;
	extern std::map<const char*, InsideTag, symphas::internal::any_case_comparator> init_tag_key_map;


	template<typename T, size_t D>
	using ic_iterator_difference_type = symphas::internal::iterator_region_difference_type<T, D>;
	//struct ic_iterator_difference_type : 
	//{
	//	using parent_type = symphas::internal::iterator_region_difference_type<T, D>;
	//	using parent_type::parent_type;
	//
	//
	//};

	template<typename T, size_t D, typename E = InitialConditionsData<D>>
	struct ic_iterator : symphas::iterator_type_impl<ic_iterator<T, D, E>,
		std::random_access_iterator_tag,
		T,
		T&,
		std::ptrdiff_t,
		T*,
		int
	>
	{
		using difference_type = std::ptrdiff_t;

		//! Create an iterator starting at the given position.
		/*!
		 * Create an iterator over an expression starting at the given
		 * position. The expression is explicitly given.
		 *
		 * \param e The expression for this iterator.
		 * \param pos The index of the underlying data in the expression
		 * which is the first index in the iterator.
		 */
		explicit ic_iterator(Axis ax, T* values, E const& e, grid::region_interval<D> const& interval, difference_type pos = 0)
			: ptr{ values, interval, pos }, init{ static_cast<E const*>(&e) }, ax{ ax } {}

		explicit ic_iterator(Axis ax, T* values, E const& e, grid::region_extent<D> const& region, difference_type pos = 0)
			: ptr{ values, region, pos }, init{ static_cast<E const*>(&e) }, ax{ ax } {}

		explicit ic_iterator(Axis ax, T* values, E const& e, symphas::interval_data_type const& interval, difference_type pos = 0)
			: ic_iterator(ax, values, e, grid::region_interval<D>(interval), pos) {}


		//! Dereference the iterator.
		inline T operator*() const
		{
			return value_fill<T>{}(ax, ptr.ptr[*ptr], (*init)[ptr.pos]);
		}

		//! Dereference past the iterator.
		inline T operator[](difference_type given_pos) const
		{
			return value_fill<T>{}(ax, ptr.ptr[ptr[given_pos]], (*init)[ptr.pos + given_pos]);
		}

		//! Member access of the iterated expression.
		inline E* operator->() const
		{
			return init;
		}

		symphas::internal::ic_iterator_difference_type<T, D> ptr;
		E const* init;			//!< Pointer to the initial conditions generator.
		Axis ax;				//!< Axis of the entries to fill. If NONE, all entries are initialized.
	};



}


namespace symphas
{

	//! Creates a tag index from the given ::InsideTag. \relatesalso init_entry_type const& InsideTag Inside
	/*!
	 * The tag is created by bit-shifting the number of bits equal to the index
	 * of the InsideTag. This represents the tag value associated with
	 * determining how the initial condition algorithm functions.
	 *
	 * \param tag The ::InsideTag which generates the index value.
	 */
	inline size_t build_intag(InsideTag tag)
	{
		if (tag == InsideTag::NONE)
		{
			return 0ull;
		}
		else
		{
			return (1ull << symphas::internal::tagtoul(tag));
		}
	}

	template<typename... tag_types>
	inline size_t build_intag(tag_types... tags)
	{
		return (0ull | ... | build_intag(tags));
	}

	//! Creates a tag index from the given ::InsideTag and index.
	/*
	 * The tag is created by bit-shifting the number of bits equal to the index
	 * of the ::InsideTag, and then combining that with the existing index of
	 * the tag. The result represents the tag value associated with
	 * determining how the initial condition algorithm functions.
	 *
	 * \param index The existing index value of the tag.
	 * \param tag The ::InsideTag which generates the index value.
	 */
	inline size_t build_intag(size_t index, InsideTag tag)
	{
		return index | build_intag(tag);
	}


	namespace internal
	{
		inline bool tag_bit_compare(size_t a, InsideTag b)
		{
			return (a & build_intag(b)) > 0;
		}
	}




	//! Contains an array of all possible parameters for initial conditions.
	/*!
	 * Manages the list of parameters which are used by the initial conditions.
	 * Each algorithm has different methods of using the parameters, so the
	 * parameters are set based on the algorithm.
	 */
	struct init_data_parameters
	{
		double* gp;			//!< The generation parameters for the algorithm generating the values.
		size_t N;			//!< The number of parameters.

		template<typename... Ts>
		init_data_parameters(double param0, Ts... params) : gp{ new double[sizeof...(Ts) + 1] {param0, (double)params...} }, N{ sizeof...(Ts) + 1 } {}

		init_data_parameters(size_t N = NUM_INIT_CONSTANTS) : gp{ (N > 0) ? new double[N] {0} : nullptr }, N{ N } {}
		init_data_parameters(init_data_parameters const& other) : init_data_parameters(other.N)
		{
			std::copy(other.gp, other.gp + other.N, gp);
		}
		init_data_parameters(init_data_parameters&& other) noexcept : init_data_parameters(0)
		{
			swap(*this, other);
		}
		init_data_parameters& operator=(init_data_parameters other)
		{
			swap(*this, other);
			return *this;
		}

		~init_data_parameters()
		{
			delete[] gp;
		}

		friend void swap(init_data_parameters& first, init_data_parameters& second)
		{
			using std::swap;
			swap(first.gp, second.gp);
			swap(first.N, second.N);
		}

		//! Creates the list of parameters which are all set to 1.
		/*!
		 * Creates the list of parameters which are all set to 1.
		 */
		static init_data_parameters one(size_t N = NUM_INIT_CONSTANTS)
		{
			init_data_parameters tp(N);
			for (iter_type i = 0; i < N; ++i)
			{
				tp.gp[i] = 1;
			}
			return tp;
		}


	protected:

		friend struct init_entry_type;

		void clear()
		{
			gp = nullptr;
			N = 0;
		}
	};

	//! Contains data representing how to read data as initial conditions.
	struct init_data_read
	{
		init_data_read(const char* name, iter_type index) :
			index{ index }, name{ (name && std::strlen(name) > 0) ? new char[std::strlen(name) + 1] : nullptr }
		{
			if (this->name)
			{
				std::strcpy(this->name, name);
			}
		}


		init_data_read(const char* name) : init_data_read(name, 0) {}
		init_data_read() : init_data_read("", 0) {}

		init_data_read(init_data_read const& other) : init_data_read(other.name, other.index) {}
		init_data_read(init_data_read&& other) : init_data_read() { swap(*this, other); }
		init_data_read& operator=(init_data_read other) { swap(*this, other); return *this; }


		friend void swap(init_data_read& first, init_data_read& second)
		{
			using std::swap;
			swap(first.index, second.index);
			swap(first.name, second.name);
		}

		//! Get the index that the file object refers to.
		/*!
		 * The index is associated with the data that is read from the file.
		 * When this object is used to read data at the given file, the index
		 * might be used to get to the correct index data in the file, or it
		 * might be compared to the index found in the file.
		 */
		iter_type get_index() const
		{
			return index;
		}

		//! Get the name of the file that will be accessed. 
		/*!
		 * Get the name of the file that will be accessed. If this is referring
		 * to checkpoint data, then only the name of the solution directory,
		 * that is the one that the `checkpoint` folder is, will be named.
		 */
		const char* get_name() const
		{
			return name;
		}

		//! Set the index of the read data.
		/*!
		 * Set the index of the read data.
		 */
		void set_index(iter_type index)
		{
			this->index = index;
		}

		~init_data_read()
		{
			delete[] name;
		}

	protected:

		iter_type index;					//!< The index to retrieve from the file.
		char* name;							//!< The name of the file to retrieve.

		friend struct init_entry_type;

		void clear()
		{
			name = nullptr;
			index = 0;
		}
	};


	//! Generates initial conditions using a non-SymPhas algorithm. 
	/*!
	 * Initial conditions can be generated using a non-SymPhas algorithm,
	 * by passing a functor of type `F` which accepts the flattened index
	 * position, the list of dimensions of the system for which
	 * the initial conditions are being generated, and the dimension of the
	 * system (the length of the dimensions array).
	 *
	 * \tparam F The functor type which is used to generate the initial
	 * conditions.
	 */
	template<typename F>
	struct init_data_functor;


	//! Abstract class for generating initial conditions.
	/*!
	 * Initial conditions can be generated using a non-SymPhas algorithm.
	 * In order to type erase the functor, this specialization serves as the
	 * base class, implementing an initialization function for each of the
	 * possible order parameter types.
	 *
	 * \tparam F The functor type which is used to generate the initial
	 * conditions.
	 */
	template<>
	struct init_data_functor<void>
	{
		//! Initializes a list of scalar values.
		virtual void initialize(scalar_t* values, grid::region_interval<1> const& interval) const = 0;

		//! Initializes a list of scalar values.
		virtual void initialize(scalar_t* values, grid::region_interval<2> const& interval) const = 0;

		//! Initializes a list of scalar values.
		virtual void initialize(scalar_t* values, grid::region_interval<3> const& interval) const = 0;

		//! Initializes a list of complex values.
		virtual void initialize(complex_t* values, grid::region_interval<1> const& interval) const = 0;

		//! Initializes a list of complex values.
		virtual void initialize(complex_t* values, grid::region_interval<2> const& interval) const = 0;

		//! Initializes a list of complex values.
		virtual void initialize(complex_t* values, grid::region_interval<3> const& interval) const = 0;

		//! Initializes a list of 1 dimensional vector values.
		virtual void initialize(vector_t<1>* values, grid::region_interval<1> const& interval) const = 0;

		//! Initializes a list of 2 dimensional vector values.
		virtual void initialize(vector_t<2>* values, grid::region_interval<2> const& interval) const = 0;

		//! Initializes a list of 3 dimensional vector values.
		virtual void initialize(vector_t<3>* values, grid::region_interval<3> const& interval) const = 0;


		//! Initializes a list of scalar values.
		void initialize(scalar_t* values, const len_type* dims, size_t dimension)
		{
			if (dimension == 1)
			{
				initialize(values, grid::region_interval<1>((len_type[1]) { dims[0] }));
			}
			else if (dimension == 2)
			{
				initialize(values, grid::region_interval<2>((len_type[2]) { dims[0], dims[1] }));
			}
			else if (dimension == 3)
			{
				initialize(values, grid::region_interval<3>((len_type[3]) { dims[0], dims[1], dims[2] }));
			}
			else
			{
				throw;
			}
		}


		//! Generate a copy of the functor.
		virtual init_data_functor<void>* make_copy() const = 0;

		virtual ~init_data_functor() {}

	};


	template<typename init_data_functor_specialized>
	struct init_data_functor_impl : init_data_functor<void>
	{
		//! Initializes a list of scalar values.
		void initialize(scalar_t* values, grid::region_interval<1> const& interval) const override
		{
			cast().initialize(values, interval, 1);
		}

		//! Initializes a list of scalar values.
		void initialize(scalar_t* values, grid::region_interval<2> const& interval) const override
		{
			cast().initialize(values, interval, 2);
		}

		//! Initializes a list of scalar values.
		void initialize(scalar_t* values, grid::region_interval<3> const& interval) const override
		{
			cast().initialize(values, interval, 3);
		}

		//! Initializes a list of complex values.
		void initialize(complex_t* values, grid::region_interval<1> const& interval) const override
		{
			cast().initialize(values, interval, 1);
		}

		//! Initializes a list of complex values.
		void initialize(complex_t* values, grid::region_interval<2> const& interval) const override
		{
			cast().initialize(values, interval, 2);
		}

		//! Initializes a list of complex values.
		void initialize(complex_t* values, grid::region_interval<3> const& interval) const override
		{
			cast().initialize(values, interval, 3);
		}

		//! Initializes a list of 1 dimensional vector values.
		void initialize(vector_t<1>* values, grid::region_interval<1> const& interval) const override
		{
			cast().initialize(values, interval, 1);
		}

		//! Initializes a list of 2 dimensional vector values.
		void initialize(vector_t<2>* values, grid::region_interval<2> const& interval) const override
		{
			cast().initialize(values, interval, 2);
		}

		//! Initializes a list of 3 dimensional vector values.
		void initialize(vector_t<3>* values, grid::region_interval<3> const& interval) const override
		{
			cast().initialize(values, interval, 3);
		}

		const init_data_functor_specialized& cast() const
		{
			return *static_cast<init_data_functor_specialized const*>(this);
		}

		init_data_functor<void>* make_copy() const override
		{
			return new init_data_functor_specialized(cast());
		}
	};

	template<typename F>
	struct init_data_functor : init_data_functor_impl<init_data_functor<F>>
	{
		using parent_type = init_data_functor<void>;
		using ret_type = std::invoke_result_t<F, iter_type, len_type const*, size_t>;

		//! Create the initial condition algorithm based on a functor.
		/*!
		 * Create the initial condition algorithm based on a functor.
		 */
		init_data_functor(F f) : f{ f } {}

		template<size_t D>
		void initialize(ret_type* values, grid::region_interval<D> const& interval)
		{
			auto it = symphas::data_iterator_region(values, interval);

#			pragma omp parallel for
			for (iter_type n = 0; n < grid::length<D>(interval); ++n)
			{
				*it++ = f(n, interval.dims, D);
			}
		}

		F f;
	};


	template<typename T>
	struct init_data_functor<Block<T>> : init_data_functor_impl<init_data_functor<Block<T>>>
	{
		using parent_type = init_data_functor<void>;

		//! Create the initial condition algorithm based on a functor.
		/*!
		 * Create the initial condition algorithm based on a functor.
		 */
		init_data_functor(Block<T> const& source) : source{ &source } {}

		template<size_t D>
		void initialize(T* values, grid::region_interval<D> const& interval)
		{
			auto it = symphas::data_iterator_region(values, interval);
			std::copy(symphas::data_iterator(source), symphas::data_iterator(source, grid::length<D>(interval)), it);
		}

		Block<T> const* source;
	};


	template<size_t N, typename T>
	struct init_data_functor<MultiBlock<N, T>> : init_data_functor_impl<init_data_functor<MultiBlock<N, T>>>
	{
		using parent_type = init_data_functor<void>;

		//! Create the initial condition algorithm based on a functor.
		/*!
		 * Create the initial condition algorithm based on a functor.
		 */
		init_data_functor(MultiBlock<N, T> const& source) : source{ &source } {}

		template<size_t D>
		void initialize(T* values, grid::region_interval<D> const& interval)
		{
			auto it = symphas::data_iterator_region(values, interval);
			std::copy(symphas::data_iterator(source), symphas::data_iterator(source, grid::length<D>(interval)), it);
		}

		MultiBlock<N, T> const* source;
	};

	template<typename T, size_t D>
	init_data_functor(Grid<T, D>) -> init_data_functor<Block<T>>;
	template<typename T, size_t D>
	init_data_functor(BoundaryGrid<T, D>) -> init_data_functor<Block<T>>;
	template<typename T, size_t D>
	init_data_functor(RegionalGrid<T, D>) -> init_data_functor<Block<T>>;
	template<typename T, size_t D>
	init_data_functor(Grid<any_vector_t<T, D>, D>) -> init_data_functor<MultiBlock<D, T>>;
	template<typename T, size_t D>
	init_data_functor(BoundaryGrid<any_vector_t<T, D>, D>) -> init_data_functor<MultiBlock<D, T>>;
	template<typename T, size_t D>
	init_data_functor(RegionalGrid<any_vector_t<T, D>, D>) -> init_data_functor<MultiBlock<D, T>>;
	template<typename T>
	init_data_functor(Block<T>) -> init_data_functor<Block<T>>;

	//! Stores information about a initial condition expression.
	/*!
	 * Initial conditions can be defined using a symbolic algebra expression,
	 * in terms of the axis variables (\f$x\f$, \f$y\f$ and \f$z\f$), as well as the
	 * coefficients which are passed as part of the configuration. Without
	 * using a configuration, the coefficients can be initialized with a member
	 * function.
	 */
	struct init_data_expr : init_data_read
	{
		init_data_expr(const char* name, const double* coeff, size_t num_coeff) :
			init_data_read(name), coeff{ (num_coeff > 0) ? new double[num_coeff] : nullptr }, num_coeff{ num_coeff }
		{
			std::copy(coeff, coeff + num_coeff, this->coeff);
		}

		init_data_expr(const char* name) : init_data_expr(name, nullptr, 0) {}
		init_data_expr() : init_data_read(), coeff{ nullptr }, num_coeff{ 0 } {}

		init_data_expr(init_data_expr const& other) : init_data_expr(other.name, other.coeff, other.num_coeff) {}
		init_data_expr(init_data_expr&& other) : init_data_expr() { swap(*this, other); }
		init_data_expr& operator=(init_data_expr other) { swap(*this, other); return *this; }

		//! Set the coefficients for the initial expression.
		/*!
		 * Set the coefficients to be used in the initial expression.
		 *
		 * \param new_coeff The list of new coefficients to store in this
		 * object.
		 * \param new_num_coeff The length of the new list of coefficients.
		 */
		void set_coeff(const double* new_coeff, size_t new_num_coeff)
		{
			if (num_coeff != new_num_coeff)
			{
				delete[] coeff;
				coeff = new double[new_num_coeff];
			}

			std::copy(new_coeff, new_coeff + new_num_coeff, coeff);
			num_coeff = new_num_coeff;
		}

		//! Return the list of coefficients in the initial expression.
		const double* get_coeff() const
		{
			return coeff;
		}

		//! Return the number of coefficients in the initial expression.
		size_t get_num_coeff() const
		{
			return num_coeff;
		}

		friend void swap(init_data_expr& first, init_data_expr& second)
		{
			using std::swap;
			swap(static_cast<init_data_read&>(first), static_cast<init_data_read&>(second));
			swap(first.coeff, second.coeff);
			swap(first.num_coeff, second.num_coeff);
		}

	protected:

		double* coeff;
		size_t num_coeff;

		friend struct init_entry_type;

		void clear()
		{
			coeff = nullptr;
			name = nullptr;
			num_coeff = 0;
		}

	};

	struct init_entry_data_type
	{
		init_entry_data_type(init_data_functor<void>* f_init) : f_init{ f_init }, data{  }, file{  }, expr_data{  } {}
		init_entry_data_type(init_data_parameters data) : f_init{ nullptr }, data{ data }, file{  }, expr_data{  } {}
		init_entry_data_type(init_data_read file) : f_init{ nullptr }, data{  }, file{ file }, expr_data{  } {}
		init_entry_data_type(init_data_expr expr_data) : f_init{ nullptr }, data{  }, file{  }, expr_data{ expr_data } {}
		init_entry_data_type() : f_init{ nullptr }, data{  }, file{  }, expr_data{  } {}

		init_data_functor<void>* f_init;	//!< Separate functor to generate the initial conditions.
		init_data_parameters data;			//!< The parameters used by the initial condition algorithm.
		init_data_read file;				//!< The file from which the field is populated.
		init_data_expr expr_data;			//!< The expression that will populate the field.

		~init_entry_data_type()
		{
			delete f_init;
		}

	};

	//template<typename F>
	//init_data_functor(F)->init_data_functor<F>;

	//! Collection of data representing the initial conditions.
	/*!
	 * This type represents the initial conditions of a phase field system,
	 * or can represent the initial conditions of a grid in general.
	 *
	 * The initial conditions can also be sourced from a file, in which case
	 * the file name should conform to *SymPhas* standards, as the file
	 * name can't be set directly. See ::init_data_read.
	 */
	struct init_entry_type : init_entry_data_type
	{

		using init_entry_data_type::f_init;
		using init_entry_data_type::data;
		using init_entry_data_type::file;
		using init_entry_data_type::expr_data;

		//! Create initial conditions data with no parameters.
		/*!
		 * The default initial condition is disabled.
		 */
		init_entry_type() :
			init_entry_data_type(), in{ Inside::NONE }, intag{ symphas::build_intag(InsideTag::NONE) } {}

		//! Create initial conditions data.
		/*!
		 * Initial conditions data is set by the given initial conditions
		 * value that represents the algorithm, and the provided modifier
		 * value. The modifier that is passed must be appropriately constructed
		 * using symphas::build_intag().
		 * The parameters to the initial conditions are also passed.
		 *
		 * \param in A value representing initial condition generation
		 * algorithm.
		 * \param intag The modifiers to the initial conditions algorithm.
		 * \param data Parameters used in generating the initial conditions.
		 */
		init_entry_type(Inside in, size_t intag, init_data_parameters data) :
			init_entry_data_type(data), in{ in }, intag{ intag } {}

		//! Create initial conditions data.
		/*!
		 * Initial conditions data is set by the given initial conditions
		 * value that represents the algorithm, and the provided modifier
		 * value. The parameters to the initial conditions are also passed.
		 *
		 * \param in A value representing initial condition generation
		 * algorithm.
		 * \param intag The modifiers to the initial conditions algorithm.
		 * \param data Parameters used in generating the initial conditions.
		 */
		init_entry_type(Inside in, InsideTag intag, init_data_parameters data) :
			init_entry_type(in, symphas::build_intag(intag), data) {}

		//! Create initial conditions data.
		/*!
		 * Initial conditions data is set by the given initial conditions
		 * value that represents the algorithm. No modifier is provided and
		 * is set to InsideTag::DEFAULT. The parameters to the initial
		 * conditions are also passed.
		 *
		 * \param in A value representing initial condition generation
		 * algorithm.
		 * \param data Parameters used in generating the initial conditions.
		 */
		init_entry_type(Inside in, init_data_parameters data) :
			init_entry_type(in, InsideTag::DEFAULT, data) {}

		//! Create initial conditions data.
		/*!
		 * Initial conditions data is set by the given initial conditions
		 * value. No parameters are passed and the parameters of the initial
		 * condition algorithm are all set to 1.
		 *
		 * \param in A value representing initial condition generation
		 * algorithm.
		 */
		init_entry_type(Inside in) :
			init_entry_type(in, InsideTag::DEFAULT, init_data_parameters::one()) {}

		//! Create initial conditions data from a file.
		/*!
		 * Initial conditions data is sourced from a file. A modifier value
		 * is not provided as it is not applicable.
		 *
		 * \param in A value representing initial condition generation
		 * algorithm.
		 * \param file Information about the file.
		 */
		init_entry_type(Inside in, init_data_read file) :
			init_entry_data_type(file), in{ in }, intag{ symphas::build_intag(InsideTag::DEFAULT) } {}

		//! Create initial conditions data from a file.
		/*!
		 * Initial conditions data is sourced from a file. The file type
		 * is assumed to be Inside::FILE.
		 *
		 * \param file Information about the file.
		 */
		init_entry_type(init_data_read file) :
			init_entry_type(Inside::FILE, file) {}

		//! Create initial conditions data from a file.
		/*!
		 * Initial conditions data is sourced from a file. The file type
		 * is assumed to be Inside::FILE.
		 *
		 * \param file Information about the file.
		 */
		init_entry_type(init_data_expr expr_data) :
			init_entry_data_type(expr_data), in{ Inside::EXPRESSION }, intag{ symphas::build_intag(InsideTag::DEFAULT) } {}

		//! Create initial conditions using a given functor.
		/*!
		 * Initial conditions data is generated using a functor from the
		 * given object. The underlying functor must be a function of only
		 * the flattened index position and the dimensions.
		 *
		 * \param f The functor object which will be used to generate
		 * the initial conditions.
		 */
		template<typename F>
		init_entry_type(init_data_functor<F> const& f) :
			init_entry_data_type(f.make_copy()), in{ Inside::NONE }, intag{ symphas::build_intag(InsideTag::DEFAULT) } {}

		template<typename F, typename = std::invoke_result_t<F, iter_type, len_type const*, size_t>>
		init_entry_type(F&& f) : init_entry_type(init_data_functor{ std::forward<F>(f) }) {}


		init_entry_type(init_entry_type&& other) noexcept : init_entry_type()
		{
			swap(*this, other);
		}

		init_entry_type(init_entry_type const& other) : init_entry_type(other.in, other.intag, *static_cast<init_entry_type const*>(&other)) {}

		init_entry_type operator=(init_entry_type other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(init_entry_type& first, init_entry_type& second)
		{
			using std::swap;
			swap(*static_cast<init_entry_data_type*>(&first), *static_cast<init_entry_data_type*>(&second));
			swap(first.in, second.in);
			swap(first.intag, second.intag);
		}

		Inside in;							//!< Type of interior random generation, min/max values.
		size_t intag;						//!< Modifies the algorithm generating interior values.


	protected:

		init_entry_type(Inside in, size_t intag, init_entry_data_type other) : init_entry_data_type(other), in{ in }, intag{ intag } {}

	};

	void swap(init_entry_type& first, init_entry_type& second);

	struct init_data_type : std::map<Axis, init_entry_type>
	{
		using parent_type = std::map<Axis, init_entry_type>;
		using parent_type::parent_type;

		init_data_type(init_entry_type const& tdata) : parent_type({ { Axis::NONE, tdata } }) {}
		template<typename... Ts>
		init_data_type(Inside in, Ts&&... args) : init_data_type(init_entry_type(in, std::forward<Ts>(args)...)) {}
		init_data_type(Inside in, symphas::init_data_parameters const& args) : init_data_type(init_entry_type(in, args)) {}
		init_data_type() : init_data_type(init_entry_type(Inside::NONE)) {}
	};

	//! From the given string, get the corresponding initial condition.
	/*!
	 * The initial condition is of type ::Inside. The given string is a key
	 * representing one of the possible initial conditions, which it will
	 * return.
	 */
	Inside in_from_str(const char* type);

	//! From the given initial condition, find its corresponding key.
	/*!
	 * The initial condition is of type ::Inside. One of the strings that
	 * represents the initial condition is returned.
	 */
	const char* str_from_in(Inside in);

	//! From the given string, get the corresponding initial tag.
	/*!
	 * The initial tag is of type ::InsideTag. The given string is a key
	 * representing one of the possible initial tags, which it will
	 * return.
	 */
	InsideTag in_tag_from_str(const char* type);

	//! From the given initial tag, find its corresponding key.
	/*!
	 * The initial tag is of type ::InsideTag. One of the strings that
	 * represents the initial tag is returned.
	 */
	const char* str_from_in_tag(InsideTag in);

}

inline void swap(symphas::init_entry_type& first, symphas::init_entry_type& second)
{
	symphas::swap(first, second);
}


// ********************************************************************************

namespace symphas::internal
{

	//! Generates a list of values distributed around a mean.
	/*!
	 * Generates a list of `D`-dimensional values where the values are uniformly
	 * distributed around the given `D`-dimensional mean, and the
	 * width of the distribution is scaled.
	 *
	 * \tparam T The type of the value which is randomly chosen.
	 * \tparam D The dimension of elements in the list.
	 */
	template<typename T, size_t D>
	struct RandomOffsets
	{
		RandomOffsets() {}

		//! Generate a given number of random offsets.
		/*!
		 * Create \p n random offsets, each with the given means.
		 *
		 * \param n The number of offsets to create.
		 * \param means The mean of each offset.
		 * \param s The offset from the means, which is multiplied by
		 * #IC_RND_STRENGTH.
		 */
		template<typename S>
		RandomOffsets(size_t n, S const* means, double s);


		//! Generate a given number of random offsets.
		/*!
		 * Create \p n random offsets. Each offset has a mean around 0, and the
		 * given scale indicates the range around which to place the offset.
		 *
		 * \param n The number of offsets to create.
		 * \param scales The offset from 0, which is multiplied by
		 * #IC_RND_STRENGTH.
		 */
		template<typename S>
		RandomOffsets(size_t n, S const* scales);

		//! Get the element at the chosen index.
		const T* get_offset(size_t i) const
		{
			return offsets[i].data();
		}

		void set_offset(size_t i, const T* values) const
		{
			for (iter_type n = 0; n < D; ++n)
			{
				offsets[i].data()[n] = values[n];
			}
		}


	private:
		std::vector<std::array<T, D>> offsets;
	};

	//! Specialization based on symphas::internal::RandomOffets.
	template<typename T>
	struct RandomOffsets<T, 1>
	{
		RandomOffsets() {}
		template<typename S>
		RandomOffsets(size_t n, S const* means, double s);

		template<typename S>
		RandomOffsets(size_t n, S const* scales);

		template<typename S, typename = std::enable_if_t<(!std::is_pointer_v<S> && !std::is_array_v<S>), int>>
		RandomOffsets(size_t n, S const& mean, double s) : RandomOffsets(n, &mean, s) {}

		auto get_offset(size_t i) const
		{
			return offsets[i];
		}

		void set_offset(size_t i, T const& value)
		{
			offsets[i] = value;
		}

	protected:
		std::vector<T> offsets;
	};


	template<>
	template<typename S>
	inline RandomOffsets<double, 1>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = means[0] + -s * IC_RND_STRENGTH;
		double x_max = means[0] + s * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(x_min, x_max);

		for (auto& e : offsets)
		{
			e = dis(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<double, 1>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = -scales[0] * IC_RND_STRENGTH;
		double x_max = scales[0] * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(x_min, x_max);

		for (auto& e : offsets)
		{
			e = dis(gen);
		}
	}


	template<>
	template<typename S>
	inline RandomOffsets<double, 2>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = means[0] + -s * IC_RND_STRENGTH;
		double x_max = means[0] + s * IC_RND_STRENGTH;
		double y_min = means[1] + -s * IC_RND_STRENGTH;
		double y_max = means[1] + s * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(x_min, x_max);
		std::uniform_real_distribution<> disy(y_min, y_max);

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<double, 2>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = -scales[0] * IC_RND_STRENGTH;
		double x_max = scales[0] * IC_RND_STRENGTH;
		double y_min = -scales[1] * IC_RND_STRENGTH;
		double y_max = scales[1] * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(x_min, x_max);
		std::uniform_real_distribution<> disy(y_min, y_max);

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
		}
	}


	template<>
	template<typename S>
	inline RandomOffsets<double, 3>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = means[0] + -s * IC_RND_STRENGTH;
		double y_min = means[1] + -s * IC_RND_STRENGTH;
		double z_min = means[2] + -s * IC_RND_STRENGTH;
		double x_max = means[0] + s * IC_RND_STRENGTH;
		double y_max = means[1] + s * IC_RND_STRENGTH;
		double z_max = means[2] + s * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		if (z_min > z_max)
		{
			std::swap(z_min, z_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(x_min, x_max);
		std::uniform_real_distribution<> disy(y_min, y_max);
		std::uniform_real_distribution<> disz(z_min, z_max);

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
			e[2] = disz(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<double, 3>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = -scales[0] * IC_RND_STRENGTH;
		double y_min = -scales[1] * IC_RND_STRENGTH;
		double z_min = -scales[2] * IC_RND_STRENGTH;
		double x_max = scales[0] * IC_RND_STRENGTH;
		double y_max = scales[1] * IC_RND_STRENGTH;
		double z_max = scales[2] * IC_RND_STRENGTH;

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		if (z_min > z_max)
		{
			std::swap(z_min, z_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(x_min, x_max);
		std::uniform_real_distribution<> disy(y_min, y_max);
		std::uniform_real_distribution<> disz(z_min, z_max);

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
			e[2] = disz(gen);
		}
	}



	template<>
	template<typename S>
	inline RandomOffsets<len_type, 1>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = std::round(means[0] + -s * IC_RND_STRENGTH);
		double x_max = std::round(means[0] + s * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));

		for (auto& e : offsets)
		{
			e = dis(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<len_type, 1>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = std::round(-scales[0] * IC_RND_STRENGTH);
		double x_max = std::round(scales[0] * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));

		for (auto& e : offsets)
		{
			e = dis(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<len_type, 2>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = std::round(means[0] + -s * IC_RND_STRENGTH);
		double y_min = std::round(means[1] + -s * IC_RND_STRENGTH);
		double x_max = std::round(means[0] + s * IC_RND_STRENGTH);
		double y_max = std::round(means[1] + s * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_int_distribution<> disx(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));
		std::uniform_int_distribution<> disy(static_cast<iter_type>(y_min), static_cast<iter_type>(y_max));

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<len_type, 2>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = std::round(-scales[0] * IC_RND_STRENGTH);
		double y_min = std::round(-scales[1] * IC_RND_STRENGTH);
		double x_max = std::round(scales[0] * IC_RND_STRENGTH);
		double y_max = std::round(scales[1] * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_int_distribution<> disx(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));
		std::uniform_int_distribution<> disy(static_cast<iter_type>(y_min), static_cast<iter_type>(y_max));

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
		}
	}

	template<>
	template<typename S>
	inline RandomOffsets<len_type, 3>::RandomOffsets(size_t n, S const* means, double s)
	{
		offsets.resize(n);

		double x_min = std::round(means[0] + -s * IC_RND_STRENGTH);
		double y_min = std::round(means[1] + -s * IC_RND_STRENGTH);
		double z_min = std::round(means[2] + -s * IC_RND_STRENGTH);
		double x_max = std::round(means[0] + s * IC_RND_STRENGTH);
		double y_max = std::round(means[1] + s * IC_RND_STRENGTH);
		double z_max = std::round(means[2] + s * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		if (z_min > z_max)
		{
			std::swap(z_min, z_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_int_distribution<> disx(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));
		std::uniform_int_distribution<> disy(static_cast<iter_type>(y_min), static_cast<iter_type>(y_max));
		std::uniform_int_distribution<> disz(static_cast<iter_type>(z_min), static_cast<iter_type>(z_max));

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
			e[2] = disz(gen);
		}
	}


	template<>
	template<typename S>
	inline RandomOffsets<len_type, 3>::RandomOffsets(size_t n, S const* scales)
	{
		offsets.resize(n);

		double x_min = std::round(-scales[0] * IC_RND_STRENGTH);
		double y_min = std::round(-scales[1] * IC_RND_STRENGTH);
		double z_min = std::round(-scales[2] * IC_RND_STRENGTH);
		double x_max = std::round(scales[0] * IC_RND_STRENGTH);
		double y_max = std::round(scales[1] * IC_RND_STRENGTH);
		double z_max = std::round(scales[2] * IC_RND_STRENGTH);

		if (x_min > x_max)
		{
			std::swap(x_min, x_max);
		}

		if (y_min > y_max)
		{
			std::swap(y_min, y_max);
		}

		if (z_min > z_max)
		{
			std::swap(z_min, z_max);
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_int_distribution<> disx(static_cast<iter_type>(x_min), static_cast<iter_type>(x_max));
		std::uniform_int_distribution<> disy(static_cast<iter_type>(y_min), static_cast<iter_type>(y_max));
		std::uniform_int_distribution<> disz(static_cast<iter_type>(z_min), static_cast<iter_type>(z_max));

		for (auto& e : offsets)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
			e[2] = disz(gen);
		}
	}








	//! Generates a list of values distributed around 0.
	/*!
	 * Generates a list of `D`-dimensional values where the values are uniformly
	 * distributed around 0, such that the width of the distribution is equal to
	 * `D`-dimensional dims and each one scaled.
	 */
	template<size_t D>
	struct RandomDeltas
	{
		RandomDeltas() {}

		//! Create a given number of point deltas.
		/*!
		 * Create \p n point deltas of the prescribed dimension, where each
		 * point delta is chosen from a range given by the parameter \p dims.
		 *
		 * \param n The number of point deltas to generate.
		 * \param dims The ranges with which each axis of the point is spread.
		 * \param s The strength applied to the spread.
		 */
		RandomDeltas(size_t n, len_type const* dims, double s = 1);

		//! Create the delts from an existing list.
		/*!
		 * \param points The list of points to copy in.
		 */
		RandomDeltas(std::vector<std::array<double, D>> const& points) : points{ points } {}

		//! Return the element at the chosen index..
		std::array<double, D> get_delta(size_t i) const
		{
			return points[i];
		}

		//! Return the list of all elements.
		auto get_deltas() const
		{
			return points;
		}

		auto set_delta(size_t i, std::array<double, D> const& value)
		{
			for (iter_type n = 0; n < D; ++n)
			{
				points[i][n] = value[n];
			}
		}

		//! Add a value to each of the random elements.
		void add_to_all(const double* value)
		{
			for (auto& e : points)
			{
				for (iter_type i = 0; i < D; ++i)
				{
					e[i] += value[i];
				}
			}
		}

		size_t size() const
		{
			return points.size();
		}

		void sort();

	private:
		std::vector<std::array<double, D>> points;
	};


	//! Specialization based on symphas::internal::RandomDeltas.
	template<>
	struct RandomDeltas<1>
	{
		RandomDeltas() {}
		RandomDeltas(size_t n, len_type const* dims, double s = 1);

		auto get_delta(size_t i) const
		{
			return points[i];
		}

		auto get_deltas() const
		{
			return points;
		}

		auto set_delta(size_t i, double const& value)
		{
			points[i] = value;
		}

		void add_to_all(const double* value)
		{
			for (auto& e : points)
			{
				e += *value;
			}
		}

		size_t size() const
		{
			return points.size();
		}

		void sort()
		{
			std::sort(points.begin(), points.end());
		}

	private:
		std::vector<double> points;
	};


	template<>
	inline void RandomDeltas<2>::sort()
	{
		std::sort(points.begin(), points.end(),
			[&] (auto a, auto b) {
				return (a[0] == b[0])
					? (a[1] < b[1])
					: (a[0] < b[0]); });
	}

	template<>
	inline void RandomDeltas<3>::sort()
	{
		std::sort(points.begin(), points.end(),
			[&] (auto a, auto b) {
				return (a[0] == b[0] && a[1] == b[1])
					? (a[2] < b[2])
					: ((a[0] == b[0])
						? (a[1] < b[1])
						: (a[0] < b[0])); });
	}

	inline RandomDeltas<1>::RandomDeltas(size_t n, len_type const* dims, double s)
	{
		points.resize(n);

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> disx(dims[0] * -s, dims[0] * s);

		for (auto& e : points)
		{
			e = disx(gen);
		}
	}

	template<>
	inline RandomDeltas<2>::RandomDeltas(size_t n, len_type const* dims, double s)
	{
		points.resize(n);

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(dims[0] * -s, dims[0] * s);
		std::uniform_real_distribution<> disy(dims[1] * -s, dims[1] * s);

		for (auto& e : points)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
		}
	}

	template<>
	inline RandomDeltas<3>::RandomDeltas(size_t n, len_type const* dims, double s)
	{
		points.resize(n);

		std::random_device rd;
		std::mt19937 gen(rd());

		std::uniform_real_distribution<> disx(dims[0] * -s, dims[0] * s);
		std::uniform_real_distribution<> disy(dims[1] * -s, dims[1] * s);
		std::uniform_real_distribution<> disz(dims[2] * -s, dims[2] * s);

		for (auto& e : points)
		{
			e[0] = disx(gen);
			e[1] = disy(gen);
			e[2] = disz(gen);
		}
	}

	inline bool is_in_square_1(iter_type cx, double ux, double aa)
	{
		return std::abs(cx - ux) < aa;
	}

	inline bool is_in_circle_1(iter_type cx, double ux, double aa)
	{
		return is_in_square_1(cx, ux, aa);
	}

	inline bool is_in_square_2(iter_type cx, double ux, double aa, iter_type cy, double uy, double bb)
	{
		return std::abs(cx - ux) < aa
			&& std::abs(cy - uy) < bb;
	}

	inline bool is_in_circle_2(iter_type cx, double ux, double aa, iter_type cy, double uy, double bb)
	{
		return (cx - ux) * (cx - ux) / (aa * aa)
			+ (cy - uy) * (cy - uy) / (bb * bb) < 1.0;
	}

	inline bool is_in_square_3(iter_type cx, double ux, double aa, iter_type cy, double uy, double bb, iter_type cz, double uz, double cc)
	{
		return std::abs(cx - ux) < aa
			&& std::abs(cy - uy) < bb
			&& std::abs(cz - uz) < cc;
	}

	inline bool is_in_circle_3(iter_type cx, double ux, double aa, iter_type cy, double uy, double bb, iter_type cz, double uz, double cc)
	{
		return (cx - ux) * (cx - ux) / (aa * aa)
			+ (cy - uy) * (cy - uy) / (bb * bb)
			+ (cz - uz) * (cz - uz) / (cc * cc) < 1.0;
	}

	inline bool is_in_square_1_dd(double ux, double aa)
	{
		return is_in_square_1(0, ux, aa);
	}

	inline bool is_in_circle_1_dd(double ux, double aa)
	{
		return is_in_circle_1(0, ux, aa);
	}

	inline bool is_in_square_2_dd(double ux, double aa, double uy, double bb)
	{
		return is_in_square_2(0, ux, 0, aa, uy, bb);
	}

	inline bool is_in_circle_2_dd(double ux, double aa, double uy, double bb)
	{
		return is_in_circle_2(0, ux, 0, aa, uy, bb);
	}

	inline bool is_in_square_3_dd(double ux, double aa, double uy, double bb, double uz, double cc)
	{
		return is_in_square_3(0, ux, aa, 0, uy, bb, 0, uz, cc);
	}

	inline bool is_in_circle_3_dd(double ux, double aa, double uy, double bb, double uz, double cc)
	{
		return is_in_circle_3(0, ux, aa, 0, uy, bb, 0, uz, cc);
	}




	inline bool is_in_square_A_1(double dx, double aa, const len_type(&dims)[1])
	{
		dx = dx - dims[0] * std::round(dx / dims[0]);
		return std::abs(dx) < aa;
	}

	inline bool is_in_circle_A_1(double dx, double aa, const len_type(&dims)[1])
	{
		return is_in_square_A_1(dx, aa, dims);
	}

	inline bool is_in_square_A_2(double dx, double aa, double dy, double bb, const len_type(&dims)[2])
	{
		dx = dx - dims[0] * std::round(dx / dims[0]);
		dy = dy - dims[1] * std::round(dy / dims[1]);
		return std::abs(dx) < aa
			&& std::abs(dy) < bb;
	}

	inline bool is_in_circle_A_2(double dx, double aa, double dy, double bb, const len_type(&dims)[2])
	{
		dx = dx - dims[0] * std::round(dx / dims[0]);
		dy = dy - dims[1] * std::round(dy / dims[1]);
		return (dx) * (dx) / (aa * aa)
			+ (dy) * (dy) / (bb * bb) < 1.0;
	}

	inline bool is_in_square_A_3(double dx, double aa, double dy, double bb, double dz, double cc, const len_type(&dims)[3])
	{
		dx = dx - dims[0] * std::round(dx / dims[0]);
		dy = dy - dims[1] * std::round(dy / dims[1]);
		dz = dz - dims[2] * std::round(dz / dims[2]);
		return std::abs(dx) < aa
			&& std::abs(dy) < bb
			&& std::abs(dz) < cc;
	}

	inline bool is_in_circle_A_3(double dx, double aa, double dy, double bb, double dz, double cc, const len_type(&dims)[3])
	{
		dx = dx - dims[0] * std::round(dx / dims[0]);
		dy = dy - dims[1] * std::round(dy / dims[1]);
		dz = dz - dims[2] * std::round(dz / dims[2]);
		return (dx) * (dx) / (aa * aa)
			+ (dy) * (dy) / (bb * bb)
			+ (dz) * (dz) / (cc * cc) < 1.0;
	}




	template<typename Pred>
	double seeds_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);
		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		iter_type const cx = n;

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i);
			if (p(cx, ux, aa))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_rnd_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);
		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		iter_type const cx = n;

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i);
			if (p(cx, ux, aa))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_A_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, RandomOffsets<len_type, 1> const* scales, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);
		iter_type const cx = n;

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i);
			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i);

			if (p(cx, ux, aa))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}

	template<typename Pred>
	double seeds_A_rnd_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, RandomOffsets<len_type, 1> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);
		iter_type const cx = n;

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i);
			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i);

			if (p(cx, ux, aa))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_B_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, RandomOffsets<len_type, 1> const* scales, double seed_val, double field_val)
	{
		return seeds_A_1(n, dims, init_data, p, offsets, scales, seed_val, field_val);
	}


	template<typename Pred>
	double seeds_B_rnd_1(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<1> const* offsets, RandomOffsets<len_type, 1> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		return seeds_A_rnd_1(n, dims, init_data, p, offsets, scales, values, rnd_offset, seed_val, field_val);
	}




	// ********************************************************************************

	template<typename Pred>
	double seeds_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = n / dims[0];

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];

			if (p(cx, ux, aa, cy, uy, bb))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_rnd_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = n / dims[0];

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];

			if (p(cx, ux, aa, cy, uy, bb))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_A_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, RandomOffsets<len_type, 2> const* scales, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = n / dims[0];

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];

			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[0];
			double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[1];

			if (p(cx, ux, aa, cy, uy, bb))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}

	template<typename Pred>
	double seeds_A_rnd_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, RandomOffsets<len_type, 2> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = n / dims[0];

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];

			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[0];
			double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[1];

			if (p(cx, ux, aa, cy, uy, bb))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_B_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, RandomOffsets<len_type, 2> const* scales, double seed_val, double field_val)
	{
		auto d = *std::min_element(dims, dims + 2);
		len_type square_dims[] = { d, d };
		return seeds_A_2(n, square_dims, init_data, p, offsets, scales, seed_val, field_val);
	}


	template<typename Pred>
	double seeds_B_rnd_2(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<2> const* offsets, RandomOffsets<len_type, 2> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		auto d = *std::min_element(dims, dims + 2);
		len_type square_dims[] = { d, d };
		return seeds_A_rnd_2(n, square_dims, init_data, p, offsets, scales, values, rnd_offset, seed_val, field_val);
	}






	// ********************************************************************************

	template<typename Pred>
	double seeds_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double cc = dims[2] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = (n / dims[0]) % dims[1],
			cz = n / (dims[0] * dims[1]);

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];
			double uz = offsets->get_delta(i)[2];

			if (p(cx, ux, aa, cy, uy, bb, cz, uz, cc))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_rnd_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;
		double cc = dims[2] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR;

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = (n / dims[0]) % dims[1],
			cz = n / (dims[0] * dims[1]);

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];
			double uz = offsets->get_delta(i)[2];

			if (p(cx, ux, aa, cy, uy, bb, cz, uz, cc))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_A_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, RandomOffsets<len_type, 3> const* scales, double seed_val, double field_val)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = (n / dims[0]) % dims[1],
			cz = n / (dims[0] * dims[1]);

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];
			double uz = offsets->get_delta(i)[2];

			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[0];
			double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[1];
			double cc = dims[2] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[2];

			if (p(cx, ux, aa, cy, uy, bb, cz, uz, cc))
			{
				return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? field_val : seed_val;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}

	template<typename Pred>
	double seeds_A_rnd_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, RandomOffsets<len_type, 3> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		size_t num_points = static_cast<size_t>(init_data.data.gp[0]);

		// x y cursor position
		iter_type const
			cx = n % dims[0],
			cy = (n / dims[0]) % dims[1],
			cz = n / (dims[0] * dims[1]);

		for (iter_type i = 0; i < num_points; ++i)
		{
			double ux = offsets->get_delta(i)[0];
			double uy = offsets->get_delta(i)[1];
			double uz = offsets->get_delta(i)[2];

			double aa = dims[0] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[0];
			double bb = dims[1] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[1];
			double cc = dims[2] * init_data.data.gp[1] * IC_SEED_SCALE_FACTOR + scales->get_offset(i)[2];

			if (p(cx, ux, aa, cy, uy, bb, cz, uz, cc))
			{
				return values->get_offset(i) + rnd_offset;
			}
		}
		return (tag_bit_compare(init_data.intag, InsideTag::INVERT)) ? seed_val : field_val;
	}


	template<typename Pred>
	double seeds_B_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, RandomOffsets<len_type, 3> const* scales, double seed_val, double field_val)
	{
		auto d = *std::min_element(dims, dims + 3);
		len_type square_dims[] = { d, d };
		return seeds_A_3(n, square_dims, init_data, p, offsets, scales, seed_val, field_val);
	}


	template<typename Pred>
	double seeds_B_rnd_3(iter_type n, len_type const* dims, symphas::init_entry_type const& init_data, Pred p,
		RandomDeltas<3> const* offsets, RandomOffsets<len_type, 3> const* scales, RandomOffsets<double, 1> const* values,
		double seed_val, double field_val, double rnd_offset)
	{
		auto d = *std::min_element(dims, dims + 3);
		len_type square_dims[] = { d, d };
		return seeds_A_rnd_3(n, square_dims, init_data, p, offsets, scales, values, rnd_offset, seed_val, field_val);
	}



	template<typename F>
	auto get_function_value(iter_type n, F f, const len_type(&dims)[1], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();

		return params[0] * f(x * params[1]);
	}

	template<typename F>
	auto get_function_value(iter_type n, F f, const len_type(&dims)[2], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();
		double y = (n / dims[0]) * vdata[Axis::Y].width() + vdata[Axis::Y].left();

		return params[0] * f(x * params[1] + y * params[2]);
	}

	template<typename F>
	auto get_function_value(iter_type n, F f, const len_type(&dims)[3], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();
		double y = ((n / dims[0]) % dims[1]) * vdata[Axis::Y].width() + vdata[Axis::Y].left();
		double z = (n / (dims[1] * dims[2])) * vdata[Axis::Z].width() + vdata[Axis::Z].left();

		return params[0] * f(x * params[1] + y * params[2] + z * params[3]);
	}

	template<typename F>
	auto get_function_value_A(iter_type n, F f, const len_type(&dims)[1], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();

		return params[0] * f(x * params[1]);
	}

	template<typename F>
	auto get_function_value_A(iter_type n, F f, const len_type(&dims)[2], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();
		double y = (n / dims[0]) * vdata[Axis::Y].width() + vdata[Axis::Y].left();

		return params[0] * (f(x * params[1]) + f(y * params[2]));
	}

	template<typename F>
	auto get_function_value_A(iter_type n, F f, const len_type(&dims)[3], symphas::interval_data_type vdata, const double* params, double zero_at_left = false)
	{
		double x = (n % dims[0]) * vdata[Axis::X].width() + vdata[Axis::X].left();
		double y = ((n / dims[0]) % dims[1]) * vdata[Axis::Y].width() + vdata[Axis::Y].left();
		double z = (n / (dims[1] * dims[2])) * vdata[Axis::Z].width() + vdata[Axis::Z].left();

		return params[0] * (f(x * params[1]) + f(y * params[2]) + f(z * params[3]));
	}


	//! Compute the radius of a bubble.
	double compute_bubble_R(symphas::interval_data_type const& vdata, double coverage, size_t N);

	//! Compute overlap given coverage.
	/*!
	 * When the coverage is 1, the overlap will be about 20% of the radius.
	 * If the coverage is 0.5, the overlap will be about -125% of the radius (meaning
	 * bubbles will be generated so they are away from each other).
	 *
	 * \param coverage How full the system should be, from 0 to 1.
	 */
	double compute_bubble_overlap(double coverage, double R, double ratio);

	//! Compute how much overlap should change by when randomly generated.
	/*!
	 * When the coverage is 1, the overlap will vary equal to its magnitude. When
	 * coverage is small, bubbles are far away from each other and should stay far,
	 * and the overlap range is fixed to the radius of a bubble.
	 *
	 * \param coverage How full the system should be, from 0 to 1.
	 */
	double compute_bubble_overlap_range(double coverage, double R, double ratio);


	//! Generate the position of bubbles in a 1d system.
	/*!
	 * Generate the position of bubbles in a 1d system.
	 *
	 * \param N The number of bubbles to generate.
	 * \param R The radius of the bubble.
	 * \param overlap How much bubbles can overlap.
	 * \param x0 The first x position.
	 */
	symphas::internal::RandomDeltas<1> get_bubble_positions_1(
		size_t N,
		double R,
		double max_overlap,
		symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps,
		const len_type* dims,
		iter_type x0);

	//! Generate the position of bubbles in a 2d system.
	/*!
	 * Generate the position of bubbles in a 2d system.
	 *
	 * \param N The number of bubbles to generate.
	 * \param R The radius of the bubble.
	 * \param overlap How much bubbles can overlap.
	 * \param x0 The first x position.
	 * \param y0 The first y position.
	 */
	symphas::internal::RandomDeltas<2> get_bubble_positions_2(
		size_t N,
		double R,
		double max_overlap,
		symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps,
		const len_type* dims,
		iter_type x0, iter_type y0);

	//! Generate the position of bubbles in a 3d system.
	/*!
	 * Generate the position of bubbles in a 3d system.
	 *
	 * \param N The number of bubbles to generate.
	 * \param R The radius of the bubble.
	 * \param overlap How much bubbles can overlap.
	 * \param x0 The first x position.
	 * \param y0 The first y position.
	 * \param z0 The first z position.
	 */
	symphas::internal::RandomDeltas<3> get_bubble_positions_3(
		size_t N,
		double R,
		double max_overlap,
		symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps,
		const len_type* dims,
		iter_type x0, iter_type y0, iter_type z0);


	template<size_t D>
	auto get_bubble_positions(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, iter_type pos0);

	template<size_t D>
	auto get_bubble_positions(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, const iter_type* pos0);


	template<>
	inline auto get_bubble_positions<1>(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, iter_type pos0)
	{
		return get_bubble_positions_1(N, R, max_overlap, overlaps, dims, pos0);
	}

	template<>
	inline auto get_bubble_positions<1>(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, const iter_type* pos0)
	{
		return get_bubble_positions_1(N, R, max_overlap, overlaps, dims, pos0[0]);
	}

	template<>
	inline auto get_bubble_positions<2>(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, const iter_type* pos0)
	{
		return get_bubble_positions_2(N, R, max_overlap, overlaps, dims, pos0[0], pos0[1]);
	}

	template<>
	inline auto get_bubble_positions<3>(size_t N, double R, double max_overlap, symphas::internal::RandomOffsets<scalar_t, 1> const& overlaps, const len_type* dims, const iter_type* pos0)
	{
		return get_bubble_positions_3(N, R, max_overlap, overlaps, dims, pos0[0], pos0[1], pos0[2]);
	}


	template<size_t D>
	void get_num_hex_tiles(len_type(&num_tiles)[D], size_t N, len_type const* dims);
	template<> void get_num_hex_tiles<1>(len_type(&num_tiles)[1], size_t N, len_type const* dims);
	template<> void get_num_hex_tiles<2>(len_type(&num_tiles)[2], size_t N, len_type const* dims);
	template<> void get_num_hex_tiles<3>(len_type(&num_tiles)[3], size_t N, len_type const* dims);

	template<size_t D>
	symphas::internal::RandomDeltas<D> get_hex_positions(size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<1> get_hex_positions<1>(size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<2> get_hex_positions<2>(size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<3> get_hex_positions<3>(size_t N, len_type const* dims);

	template<size_t D>
	symphas::internal::RandomDeltas<D> to_spiral_order(symphas::internal::RandomDeltas<D> const& positions, size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<1> to_spiral_order<1>(symphas::internal::RandomDeltas<1> const& positions, size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<2> to_spiral_order<2>(symphas::internal::RandomDeltas<2> const& positions, size_t N, len_type const* dims);
	template<> symphas::internal::RandomDeltas<3> to_spiral_order<3>(symphas::internal::RandomDeltas<3> const& positions, size_t N, len_type const* dims);
}

