
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
 * PURPOSE: Defines the model class, which associates phase fields with
 * the numerical solver and acts as the interface between the solver
 * and model.
 *
 * ***************************************************************************
 */

#pragma once


//! \cond
#define SN sizeof...(S)
//! \endcond


#include "solver.h"

//! Used to obtain names of phase fields for a given model.
/*!
 * Used for getting the name of the phase field at the given index, typically
 * for printing to output or for using the name in expressions.
 * 
 * Specialization of this class by model type is used to associate different
 * phase field names than the default ones. All models use default phase field
 * names unless specialized.
 */
template<template<typename> typename M>
struct model_field_name
{
	//! Returns the name of the \f$i\f$-th phase field.
	/*!
	 * Get the name of the phase field at the given index for the prescribed
	 * model.
	 * 
	 * \param i The index of the phase field.
	 */
	const char* operator()(int i)
	{
		constexpr size_t MAX_NAME_COUNT = sizeof(ORDER_PARAMETER_NAMES) / sizeof(*ORDER_PARAMETER_NAMES);

		if (i < MAX_NAME_COUNT)
		{
			return ORDER_PARAMETER_NAMES[i];
		}
		else
		{
			static size_t EXTRA_NAME_COUNT = 0;
			static std::vector<char*> EXTRA_NAMES;

			if (i - MAX_NAME_COUNT < EXTRA_NAME_COUNT)
			{
				return EXTRA_NAMES[i - MAX_NAME_COUNT];
			}
			else
			{
				char* name = new char[BUFFER_LENGTH];
				sprintf(name, ORDER_PARAMETER_NAME_EXTRA_FMT, i);
				EXTRA_NAMES.push_back(name);
				return EXTRA_NAMES.back();
			}
		}
	}
};



// *****************************************************************************************


#ifdef VTK_ON
#include "colormap.h"
#include <thread>
#endif



// *****************************************************************************************


//! Represents a phase field problem.
/*!
 * A base type for the class representing a phase field problem. Contains all 
 * information of the phase field systems. The dynamical equations are not
 * included in this implementation, but are added through specialization of this
 * class.
 * 
 * Manages the solver and acts as the interface to the solver to manage
 * the time evolution of the phase fields. Does not contain the interface to
 * the `equation` step; this is left to the specialization.
 * 
 * \tparam D The dimension size of the phase field problem, hence the dimension
 * of all systems in the model.
 * \tparam Sp The stencil type.
 * \tparam S... The phase field types of the model.
 */
template<size_t D, typename Sp, typename... S>
struct Model
{

	//! The type of the system storing the phase fields, used by the solver.
	template<typename Ty>
	using SolverSystemApplied = typename symphas::solver_system_type<Sp>::template type<Ty, D>;


protected:

	Model()
		: _s{ construct_systems({}, {}, {}, 0) }, solver{ Sp::make_solver() }, coeff{ nullptr }, num_coeff{ 0 }, 
		index{ params::start_index }, time{ 0 }
#ifdef VTK_ON
		, viz_update{ nullptr }
#endif 
	{}

public:

	//! Get the type of the `N`-th phase field.
	/*!
	 * Get the type of the `N` System from the list of systems. The system
	 * represents the phase field.
	 *
	 * \tparam N The phase field index to get the type.
	 */
	template<size_t N>
	using type_of_S = typename std::tuple_element<N, std::tuple<S...>>::type;


	//! Create a new model.
	/*!
	 * A new model is created using the problem parameters. The model
	 * coefficients for the dynamical equations are also stored.
	 *
	 * \param coeff The list of coefficients for this model.
	 * \param num_coeff The number of coefficients provided.
	 * \param parameters The problem parameters for the phase field problem
	 * represented by the model.
	 */
	Model(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		_s{ construct_systems(parameters.get_initial_data(), parameters.get_interval_data(), parameters.get_boundary_data(), parameters.length()) },
		solver{ Sp::make_solver(get_updated_parameters(parameters)) }, coeff{ (num_coeff > 0) ? new double[num_coeff] : nullptr },
		num_coeff{ num_coeff }, index{ parameters.index }, time{ parameters.time }
#ifdef VTK_ON
		, viz_update{ nullptr }
#endif
	{
		std::copy(coeff, coeff + num_coeff, this->coeff);
		visualize();
	}

	//! Create a new model without coefficients.
	/*!
	 * The coefficient list is assumed to be empty, so no coefficients are
	 * saved.
	 * 
	 * \param parameters The problem parameters for the phase field problem
	 * represented by the model.
	 */
	Model(symphas::problem_parameters_type const& parameters) : Model(nullptr, 0, parameters) {}

	Model(Model<D, Sp, S...> const& other) : 
		_s{ other._s }, solver{ other.solver }, coeff{ (other.num_coeff > 0) ? new double[other.num_coeff] : nullptr }, 
		num_coeff{ other.num_coeff }, index{ other.index }, time{ other.time }
#ifdef VTK_ON
		, viz_update{ nullptr }
#endif
	{
		std::copy(other.coeff, other.coeff + other.num_coeff, coeff);
		visualize();
	}

	Model(Model<D, Sp, S...>&& other) noexcept : Model()
	{
		swap(*this, other);
	}

	Model<D, Sp, S...>& operator=(Model<D, Sp, S...> other)
	{
		swap(*this, other);
		return *this;
	}

//	template<typename Sp0, typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
//	Model(Model<D, Sp0, S...> const& other) :
//		_s{ construct_systems(other.systems_tuple()) }, solver{ Sp::make_solver(generate_parameters()) }, coeff{ (other.get_num_coeff() > 0) ? new double[other.get_num_coeff()] : nullptr },
//		num_coeff{ other.get_num_coeff() }, index{ other.index() }, time{ other.get_time()}
//#ifdef VTK_ON
//		, viz_update{ nullptr }
//#endif
//	{
//		std::copy(other.get_coeff(), other.get_coeff() + other.get_num_coeff(), coeff);
//		update_systems(time);
//		visualize();
//	}
//
//	template<typename Sp0, typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
//	Model(Model<D, Sp0, S...>&& other) noexcept : Model()
//	{
//		swap(*this, other);
//	}
//
//	template<typename Sp0, typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
//	Model<D, Sp0, S...>& operator=(Model<D, Sp0, S...> other)
//	{
//		swap(*this, other);
//		return *this;
//	}

	~Model()
	{
		devisualize();
		delete[] coeff;
	}

	friend void swap(Model<D, Sp, S...>& first, Model<D, Sp, S...>& second)
	{
		using std::swap;


		swap(first.solver, second.solver);
		swap(first._s, second._s);
		swap(first.coeff, second.coeff);
		swap(first.num_coeff, second.num_coeff);
		swap(first.index, second.index);
		swap(first.time, second.time);

#ifdef VTK_ON
		swap(first.viz_thread, second.viz_thread);
		swap(first.viz_update, second.viz_update);
#endif
	}

	template<typename Sp0, typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
	friend void swap(Model<D, Sp, S...>& first, Model<D, Sp0, S...>& second)
	{
		using std::swap;

		swap(first._s, second._s);
		swap(first.coeff, second.coeff);
		swap(first.num_coeff, second.num_coeff);
		swap(first.index, second.index);
		swap(first.time, second.time);

#ifdef VTK_ON
		swap(first.viz_thread, second.viz_thread);
		swap(first.viz_update, second.viz_update);
#endif
	}

	//! Returns the list of coefficients.
	/*!
	 * Returns a pointer to the array of coefficients used in the dynamic equations.
	 */
	const double* get_coeff() const
	{
		return coeff;
	}

	//! Get the current solution index.
	/*!
	 * At each complete iteration of the solver, the index is incremented to
	 * track the total number of solver iterations performed. Returns this
	 * value.
	 */
	iter_type get_index() const
	{
		return index;
	}

	//! Return the number of coefficients saved by the model.
	/*!
	 * Returns the number of coefficients that are used by the model in
	 * the dynamic equations.
	 */
	auto get_num_coeff() const
	{
		return num_coeff;
	}

	//! Returns the current simulation time of the model.
	/*!
	 * Returns the current simulation time of the model.
	 */
	auto get_time() const
	{
		return time;
	}

	void set_index(iter_type index)
	{
		this->index = index;
	}

	void set_time(double time)
	{
		this->time = time;
	}

	//! Returns a symbolic variable for the time value of the model.
	/*!
	 * Returns a symbolic variable for the time value of the model. The
	 * variable maintains the time value as a pointer to the time value managed
	 * by the model.
	 */
	auto get_time_var()
	{
		return expr::make_term(TimeValue{ &time });
	}

	//! Get the number of fields of the given types.
	/*!
	 * Return the number of fields which correspond to the given types for this 
	 * model.
	 *
	 * \tparam S0 The type which is counted from the list of all types.
	 */
	template<typename S0>
	static constexpr size_t num_fields()
	{
		return ((static_cast<size_t>(std::is_same<S0, S>::value) + ...));
	}

	//! Get the number of fields of the given types.
	/*!
	 * Return the number of fields which correspond to the given types for this 
	 * model. The number returned is a sum of the number of times each type 
	 * appears in the list.
	 * 
	 * \tparam S0 The type which is counted from the list of all types.
	 * \tparam S1 The second type which is counted.
	 */
	template<typename S0, typename S1, typename... Ss>
	static constexpr size_t num_fields()
	{
		return ((static_cast<size_t>(std::is_same<S0, S>::value) + ...)) + num_fields<S1, Ss...>();
	}


	//! Get the index of the `I`-th given type.
	/*!
	 * Get the index of the selected type from the list of types for this model.
	 * The index is the position in the list, starting from 0. Determines the
	 * index of the type after appearing `I` times.
	 * 
	 * \tparam Type The type to find in the list and return its position index.
	 * \tparam I Chose the `I`-th appearance of the chosen type.
	 */
	template<typename Type, size_t I = 0>
	static constexpr int index_of_type = symphas::lib::nth_index_of_type<Type, I, S...>;

	//! Execute a function for all fields of the given type.
	/*!
	 * Execute a function for all fields of the given type.
	 * The function expects as parameters the axes values, the data values and 
	 * the index of the corresponding field, followed by desired arguments.
	 * 
	 * \param f The function to execute on the field of the given type.
	 * \param args The arguments to pass to the function.
	 * 
	 * \tparam Type The field type that is selected to be passed to the
	 * function.
	 */
	template<typename Type, typename F, typename... Args, typename std::enable_if<(num_fields<Type>() > 0), int>::type = 0>
	void do_for_field_type(F f, Args&& ...args) const
	{
		do_for_field_type<Type, num_fields<Type>() - 1>(f, std::forward<Args>(args)...);
	}

	//! Execute a function for the `N`-th field.
	/*!
	 * Execute a function for all fields of the given type.
	 * The function expects as parameters the axes values, the data values and
	 * the index of the corresponding field, followed by desired arguments.
	 *
	 * \param f The function to execute on the field of the given type.
	 * \param args The arguments to pass to the function.
	 * 
	 * \tparam N the index of the field that is passed to the function.
	 */
	template<size_t N, typename F, typename... Args, typename std::enable_if<(N <= sizeof...(S)), int>::type = 0>
	void do_for_field(F f, Args&& ... args) const
	{
		auto&& data = std::get<N>(_s).get_snapshot();
		f(data.values, data.len, std::forward<Args>(args)...);
	}

	//! Copy the values of the `N`-th field into the given array. 
	/*!
	 * Copies the values from the field into the given array, ensuring that
	 * only the true phase field values are copied. I.e. this means that for
	 * systems which use boundaries, the boundary domain is not included.
	 * 
	 * \param into The array into which to copy the values.
	 * 
	 * \tparam N The index of the phase field which is copied from.
	 */
	template<size_t N>
	void copy_field_values(type_of_S<N>* into) const
	{
		std::get<N>(_s).persist(into);
	}

	//! Copy the values of the `I`-th field of the given type.
	/*!
	 * Copy values from the field corresponding to the given type. The `I`-th
	 * field of that type is chosen from the list. Depending on the system type,
	 * the boundary domain will not be copied.
	 * 
	 * \param[in] into The array into which the field is copied into.
	 * 
	 * \tparam Type is the type to match and copy into.
	 */
	template<typename Type, size_t I, typename std::enable_if_t<(num_fields<Type>() > I), int> = 0>
	void copy_field_type_values([[maybe_unused]] Type* into) const
	{
		copy_field_values<index_of_type<Type, I>>(into);
	}

	//! Copy the values of the `I`-th field of the given type.
	/*!
	 * Copy values from the field corresponding to the given type. The `I`-th
	 * field of that type is chosen from the list. Depending on the system type,
	 * the boundary domain will not be copied.
	 *
	 * \param[in] into The array into which the field is copied into.
	 *
	 * \tparam Type is the type to match and copy into.
	 */
	template<typename Type, size_t I, typename std::enable_if_t<(num_fields<Type>() == 0), int> = 0>
	void copy_field_type_values([[maybe_unused]] Type* into) const {}

	//! Return the values of the `N`-th field as a Grid. 
	/*!
	 * Copies the values from the field and returns a new Grid, ensuring that
	 * only the true phase field values are copied. I.e. this means that for
	 * systems which use boundaries, the boundary domain is not included.
	 *
	 * \tparam N The index of the phase field which is copied from.
	 */
	template<size_t N>
	auto get_field() const
	{
		Grid<type_of_S<N>, D> out(std::get<N>(_s).get_info().get_dims().get());
		std::get<N>(_s).persist(out.values);
		return out;
	}

	//! Persists all the systems managed by the model.
	/*!
	 * Persists the phase field systems to disk, saving the entire array.
	 * This function has no effect without the **io** module.
	 * 
	 * \param dir The directory into which to store the resulting data file.
	 */
	void save_systems(const char* dir) const
	{
		save_systems(dir, std::make_index_sequence<SN>{});
	}

	//! Persists all the systems managed by the model.
	/*!
	 * Persists the phase field systems to disk, saving the entire array.
	 * This function has no effect without the **io** module.
	 *
	 * \param dir The directory into which to store the resulting data file.
	 * \param name The suffix name, to which the system id number is appended.
	 */
	void save_systems(const char* dir, const char* name) const
	{
		save_systems(dir, name, std::make_index_sequence<SN>{});
	}

	//! Updates each of the phase field systems.
	/*!
	 * Updates the systems by directly calling their `update` function.
	 * 
	 * \param time The current time of the solution.
	 */
	void update_systems(double time)
	{
		this->time = time;
		update_systems(time, std::make_index_sequence<SN>{});
	}
	
	//! Advances to the next solution iteration.
	/*!
	 * Calls the `step` function of the solver, stepping forward the solution
	 * iteration. This also increments Moodel::index.
	 * 
	 * \param dt The time increment.
	 */
	void step(double dt)
	{
		++index;
		step(std::make_index_sequence<SN>{}, dt);

#ifdef VTK_ON
		if (viz_update 
			&& params::viz_interval 
			&& index % params::viz_interval == 0)
		{
			viz_update->update();
		}
#endif
	}

	//! Returns the system at the given index.
	/*!
	 * Returns the reference to the system at the given index from the list of
	 * phase field systems of the model.
	 * 
	 * \tparam I The index of the phase field to return.
	 */
	template<size_t I>
	const auto& system() const
	{
		return std::get<I>(_s);
	}

	//! Returns the system at the given index.
	/*!
	 * Returns the reference to the system at the given index from the list of
	 * phase field systems of the model.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	template<size_t I>
	auto& system()
	{
		using chosen_sys_type = typename std::tuple_element<I, decltype(_s)>::type;
		return const_cast<chosen_sys_type&>(static_cast<const Model<D, Sp, S...>&>(*this).system<I>());
	}

	//! Returns the underlying grid of the system at the given index.
	/*!
	 * Returns a reference to the grid of the system at the given index.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	template<size_t I>
	auto& grid()
	{
		return system<I>().as_grid();
	}

	//! Returns the underlying grid of the system at the given index.
	/*!
	 * Returns a reference to the grid of the system at the given index.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	template<size_t I>
	const auto& grid() const
	{
		return system<I>().as_grid();
	}

	//! Returns the tuple of all model systems.
	/*!
	 * Directly returns the list of all systems of this model.
	 */
	const auto& systems_tuple() const
	{
		return _s;
	}

	//! Print information about this model to the given output.
	/*!
	 * The information about this model is written to the given output.
	 * Data about this model includes the number of coefficients and the number
	 * of fields.
	 * 
	 * Since the name is not directly associated with the model, the name is
	 * not provided in this information.
	 * 
	 * \param out The output file or stream to which to print the
	 * information.
	 */
	void print_info(FILE* out) const
	{
		fprintf(out, OUTPUT_BANNER);
		fprintf(out, "Phase field problem of %zd system%s:\n", SN, (SN > 1) ? "s" : "");
		print_all_dimensions(out);
		print_all_coeff(out);
		fprintf(out, OUTPUT_BANNER);
	}

	auto generate_parameters() const
	{
		symphas::problem_parameters_type parameters(SN);
		populate_parameters(parameters, std::make_index_sequence<SN>{});

		parameters.time = get_time();;
		parameters.index = index;

		return parameters;
	}

protected:
	std::tuple<SolverSystemApplied<S>...> _s;	//! Container managing pointers to phase field data.
	Sp solver;									//! Solver for determining the phase field solution.

	//! Coefficients in the equations of motion.
	/*!
	 * List of coefficients used in the equations of motion. If
	 * the number of coefficients is insufficient, then the default coefficient
	 * value #DEFAULT_COEFF_VALUE should be used.
	 */
	double *coeff;

	//! Number of supplied coefficients.
	/*!
	 * The length of the coefficient list, coeff.
	 */
	size_t num_coeff;

	//! The latest index of simulation.
	/*
	 * Internally tracks the current index in the simulation. It is 
	 * incremented for every call to Model::step(), and can be accessed
	 * (but not modified) through the function Model::index().
	 */
	iter_type index;

	//! The current simulation time.
	double time;

#ifdef VTK_ON

	std::thread viz_thread;
	ColourPlotUpdater* viz_update;
	

	void visualize()
	{
		if constexpr (num_fields<scalar_t>() > 0)
		{
			if (params::viz_interval > 0)
			{
				auto& frame = grid<index_of_type<scalar_t>>();
				viz_thread = std::thread([&] ()
				{
					ColourPlot2d viz{};
					viz.init(frame.values, frame.dims, index, viz_update);
				});
			}
		}
	}

	void devisualize() 
	{
		if constexpr (num_fields<scalar_t>() > 0)
		{
			if (params::viz_interval > 0)
			{
				viz_thread.join();
			}
		}
	}


#else
	void visualize() {}
	void devisualize() {}
#endif

	template<size_t... Is>
	void construct_interval_data(symphas::interval_data_type (&intervals)[SN], std::index_sequence<Is...>) const
	{
		(((intervals[Is] = system<Is>().get_info().intervals), ...));
	}

	auto get_updated_parameters(symphas::problem_parameters_type parameters) const
	{
		parameters.extend(SN);
		symphas::interval_data_type intervals[SN];
		construct_interval_data(intervals, std::make_index_sequence<SN>{});

		parameters.set_interval_data(intervals, SN);
		return parameters;
	}

	template<size_t... Is>
	auto populate_parameters(symphas::problem_parameters_type &parameters, std::index_sequence<Is...>) const
	{
		((fill_interval_data<Is>(system<Is>(), parameters),
			fill_boundary_data<Is>(system<Is>(), parameters),
			fill_initial_data<Is>(system<Is>(), parameters)), ...);
	}

	template<size_t I, typename T0>
	void fill_interval_data(PhaseFieldSystem<BoundaryGrid, T0, D> const& system, symphas::problem_parameters_type& parameters) const
	{
		auto intervals = system.get_info().intervals;
		for (iter_type i = 0; i < D; ++i)
		{
			Axis ax = symphas::index_to_axis(i);
			auto& interval = intervals.at(ax);

			interval.set_interval_count(
				interval.left(),
				interval.right(),
				system.dims[i]);
		}

		parameters.set_interval_data(intervals, I);
	}

	template<size_t I, template<typename, size_t> typename G, typename T0>
	void fill_interval_data(PhaseFieldSystem<G, T0, D> const& system, symphas::problem_parameters_type& parameters) const
	{
		parameters.set_interval_data(system.get_info().intervals, I);
	}

	template<size_t I, template<typename, size_t> typename G, typename T0>
	void fill_boundary_data(PhaseFieldSystem<G, T0, D> const& system, symphas::problem_parameters_type& parameters) const
	{
		symphas::b_data_type bdata;

		for (iter_type i = 0; i < D * 2; ++i)
		{
			Side side = symphas::index_to_side(i);
			bdata[side] = symphas::b_element_type(BoundaryType::PERIODIC);
		}

		parameters.set_boundary_data(bdata, I);
	}

	template<size_t I, typename T0>
	void fill_boundary_data(PhaseFieldSystem<BoundaryGrid, T0, D> const& system, symphas::problem_parameters_type& parameters) const
	{
		symphas::b_data_type bdata;

		for (iter_type i = 0; i < D * 2; ++i)
		{
			Side side = symphas::index_to_side(i);
			BoundaryType type = system.types[i];
			bdata[side] = system.boundaries[i]->get_parameters();
		}

		parameters.set_boundary_data(bdata, I);
	}

	template<size_t I, template<typename, size_t> typename G, typename T0>
	void fill_initial_data(PhaseFieldSystem<G, T0, D> const& system, symphas::problem_parameters_type& parameters) const
	{
		parameters.set_initial_data(init_data_functor(system.as_grid()), I);
	}

	template<size_t... Is>
	auto construct_systems(const symphas::init_data_type* tdata, const symphas::interval_data_type *vdata, const symphas::b_data_type *bdata, std::index_sequence<Is...>) const
	{
		return std::make_tuple(SolverSystemApplied<S>(tdata[Is], vdata[Is], bdata[Is], Is)...);
	}

	auto construct_systems(const symphas::init_data_type* tdata, const symphas::interval_data_type* vdata, const symphas::b_data_type* bdata, size_t data_len) const
	{
		if (data_len < SN)
		{
			symphas::init_data_type tdata_extend[SN];
			symphas::interval_data_type vdata_extend[SN];
			symphas::b_data_type bdata_extend[SN];

			for (size_t i = 0; i < data_len; ++i)
			{
				tdata_extend[i] = tdata[i];
				vdata_extend[i] = vdata[i];
				bdata_extend[i] = bdata[i];
			}
			for (size_t i = data_len; i < SN; ++i)
			{ 
				tdata_extend[i] = tdata_extend[0];
				vdata_extend[i] = vdata_extend[0];
				bdata_extend[i] = bdata_extend[0];
			}
			return construct_systems(tdata_extend, vdata_extend, bdata_extend, std::make_index_sequence<SN>{});
		}
		else
		{
			return construct_systems(tdata, vdata, bdata, std::make_index_sequence<SN>{});
		}
	}

	template<size_t I, typename other_sys_type,
		typename std::enable_if_t<!std::is_same<other_sys_type, std::tuple_element_t<I, std::tuple<SolverSystemApplied<S>...>>>::value, int> = 0>
	SolverSystemApplied<std::tuple_element_t<I, std::tuple<S...>>> construct_system(other_sys_type const& system) const
	{
		using value_type = std::tuple_element_t<I, std::tuple<S...>>;

		symphas::b_data_type boundaries{};
		for (iter_type i = 0; i < D * 2; ++i)
		{
			boundaries[symphas::index_to_side(i)] = BoundaryType::PERIODIC;
		}

		bool extend_boundary = params::extend_boundary;
		params::extend_boundary = true;

		auto new_system = SolverSystemApplied<std::tuple_element_t<I, std::tuple<S...>>>(
			symphas::init_data_type{}, system.info.intervals, boundaries, system.id);

		params::extend_boundary = extend_boundary;

		value_type* swap_arr = new value_type[system.length()];
		system.persist(swap_arr);
		new_system.fill(swap_arr);

		delete[] swap_arr;

		return new_system;
	}

	template<size_t I, typename other_sys_type,
		typename std::enable_if_t<std::is_same<other_sys_type, std::tuple_element_t<I, std::tuple<SolverSystemApplied<S>...>>>::value, int> = 0>
	SolverSystemApplied<std::tuple_element_t<I, std::tuple<S...>>> construct_system(other_sys_type const& system) const
	{
		return SolverSystemApplied<std::tuple_element_t<I, std::tuple<S...>>>(system);
	}

	template<typename... Ss, size_t... Is>
	std::tuple<SolverSystemApplied<S>...> construct_systems(std::tuple<Ss...> const& systems, std::index_sequence<Is...>) const
	{
		return std::make_tuple(construct_system<Is>(std::get<Is>(systems))...);
	}

	template<typename... Ss>
	auto construct_systems(std::tuple<Ss...> const& systems) const
	{
		return construct_systems(systems, std::make_index_sequence<sizeof...(Ss)>{});
	}

	template<size_t... Is>
	void update_systems(double time, std::index_sequence<Is...>)
	{
		((..., system<Is>().update(index, time)));
	}

	template<size_t... Is>
	void save_systems(const char* dir, std::index_sequence<Is...>) const
	{
		((..., std::get<Is>(_s).save(dir, index)));
	}

	template<size_t... Is>
	void save_systems(const char* dir, const char* name, std::index_sequence<Is...>) const
	{
		char* names[SN];
		for (iter_type i = 0; i < SN; ++i)
		{
			names[i] = new char[std::strlen(name) + symphas::lib::num_digits(i) + 1];
			snprintf(names[i], BUFFER_LENGTH, "%d%s", i, name);
		}

		((..., std::get<Is>(_s).save(dir, names[Is], index)));

		for (iter_type i = 0; i < SN; ++i)
		{
			delete[] names[i];
		}
	}

	template<size_t... Is>
	void step(std::index_sequence<Is...>, double dt)
	{
		((..., solver.step(system<Is>(), dt)));
	}


	/*
	 * for the I-th type field in the whole list of fields, execute
	 * the given function
	 */

	template<typename Type, size_t I, typename F, typename... Args, typename std::enable_if_t<(I > 0), int> = 0>
	void do_for_field_type(F &&f, Args&& ... args) const
	{
		do_for_field_type<Type, I - 1, Args...>(f, std::forward<Args>(args)...);
	}

	template<typename Type, size_t I, typename F, typename... Args, typename std::enable_if_t<(I == 0), int> = 0>
	void do_for_field_type(F &&f, Args&& ... args) const
	{
		do_for_field<index_of_type<Type, I>>(f, std::forward<Args>(args)...);
	}

	template<size_t I>
	void print_grid_dimensions(FILE* out) const
	{
		fprintf(out, "\tsystem<%zd> ... ", I);
		for (iter_type d = 0; d < D; ++d)
		{
			fprintf(out, "%d", grid<I>().dims[d]);
			if (d < D - 1)
			{
				fprintf(out, " x ");
			}
		}
		fprintf(out, "\n");
	}

	template<size_t... Is>
	void print_all_dimensions(FILE* out, std::index_sequence<Is...>) const
	{
		((..., print_grid_dimensions<Is>(out)));
	}

	void print_all_dimensions(FILE* out) const
	{
		print_all_dimensions(out, std::make_index_sequence<SN>{});
	}

	void print_all_coeff(FILE* out) const
	{
		if (num_coeff > 0)
		{
			fprintf(out, "Number of coefficients provided: %zd.\n", num_coeff);
			fprintf(out, "List of those coefficients:\n");
			for (iter_type i = 0; i < num_coeff; ++i)
			{
				fprintf(out, "\tc%-8d ... % .4E\n", i + 1, coeff[i]);
			}
		}
		else
		{
			fprintf(out, "No coefficients provided.\n");
		}
		fprintf(out, "Default coefficient value: % .2lf.\n", DEFAULT_COEFF_VALUE);
	}
}; 




//! Get the dimension of the phase field problem.
/*!
 * Get the dimension of the phase field problem represented by the given model.
 * The given type is assumed to have parent type model.
 * 
 * \tparam M The model representing the phase field problem
 */
template<typename M>
struct model_dimension
{

protected:

	template<size_t D>
	struct dim_wrap
	{
		static const size_t value = D;
	};

	template<size_t D, typename Sp, typename... S>
	static constexpr auto _pack_dim(Model<D, Sp, S...>*)
	{
		return dim_wrap<D>{};
	}

	static constexpr auto pack_dim(M* m)
	{
		return _pack_dim(m);
	}

	using T = typename std::invoke_result_t<decltype(&model_dimension<M>::pack_dim), M*>;

public:

	static const size_t value = T::value;
};

//! Obtain the number of phase fields of a model.
/*!
 * The number of phase field systems of the given model type is determined.
 * The given type is assumed to have parent class Model and only counts the
 * phase fields which are defined there.
 * 
 * \tparam M The model type that from which is counted the phase fields
 */
template<typename M>
struct model_num_parameters
{

protected:

	template<size_t N>
	struct count_wrap
	{
		static const size_t value = N;
	};

	template<size_t D, typename Sp, typename... S>
	static constexpr auto _pack_count(Model<D, Sp, S...>*)
	{
		return count_wrap<sizeof...(S)>{};
	}

	static constexpr auto pack_count(M* m)
	{
		return _pack_count(m);
	}

	using T = typename std::invoke_result_t<decltype(&model_num_parameters<M>::pack_count), M*>;

public:

	static const size_t value = T::value;
};


//! Obtain the solver type used by the given model.
/*!
 * Uses type traits to obtain the solver type used by the given model. The 
 * type provided is assumed to have parent type Model.
 * 
 * \tparam M The model type.
 */
template<typename M>
struct model_solver
{

protected:

	template<typename Sp>
	struct solver_wrap
	{
		using type = Sp;
	};

	template<size_t D, typename Sp, typename... S>
	static solver_wrap<Sp> _get_type(Model<D, Sp, S...>)
	{
		return {};
	}

	static auto get_type(M m)
	{
		return _get_type(m);
	}

public:

	using type = typename std::invoke_result_t<decltype(&model_solver<M>::get_type), M>::type;
};

template<typename M>
using model_solver_t = typename model_solver<M>::type;

//! Obtain the solver type used by the given model.
/*!
 * Uses type traits to obtain the solver type used by the given model. The
 * type provided is assumed to have parent type Model.
 *
 * \tparam M The model type.
 */
template<typename M>
struct model_base
{

protected:

	template<size_t D, typename Sp, typename... S>
	static Model<D, Sp, S...> _get_type(Model<D, Sp, S...> m)
	{
		return m;
	}

	static auto get_type(M m)
	{
		return _get_type(m);
	}

public:

	using type = typename std::invoke_result_t<decltype(&model_base<M>::get_type), M>;
};

template<typename M>
using model_base_t = typename model_base<M>::type;



//! Obtain the phase-field value type corresponding to the given index. 
/*!
 * Wrapper for a model's member type, Model::type_of_S, which provides
 * the type of the field at the given index.
 * 
 * \tparam M The model type.
 * \tparam N The index of the field to get the type.
 */
template<typename M, size_t N, typename = typename M::template type_of_S<0>>
struct model_field
{
	using type = typename M::template type_of_S<N>;
};

template<typename M, size_t N>
using model_field_t = typename model_field<M, N>::type;




#undef SN

