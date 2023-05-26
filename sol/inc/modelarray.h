
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
 * PURPOSE: Defines the model used for phase field crystal problems.
 *
 * ***************************************************************************
 */

#pragma once

#include <utility>
#include <array>

#include "solver.h"





namespace symphas::internal
{

	template<typename T, typename S = void>
	struct parameterized_type;

	template<typename M, typename T>
	struct model_field_parameters;

	template<>
	struct parameterized_type<void, void>
	{
		parameterized_type() : parameters{ nullptr }, len{ 0 } {}

		template<typename T0, typename... Ts>
		parameterized_type(T0 parameter0, Ts... parameters) :
			parameters{ (sizeof...(Ts) + 1 > 0) ? new double[sizeof...(Ts) + 1] {} : nullptr }, len{ sizeof...(Ts) + 1 }
		{
			populate_parameters(std::make_index_sequence<sizeof...(Ts) + 1>{}, parameter0, parameters...);
		}

		template<typename M, typename T>
		parameterized_type(model_field_parameters<M, T>&& other) : 
			parameterized_type(static_cast<parameterized_type<void, void>&&>(other)) {}

		template<typename M, typename T>
		parameterized_type(model_field_parameters<M, T> const& other) : 
			parameterized_type(static_cast<parameterized_type<void, void> const&>(other)) {}

		parameterized_type(parameterized_type<void, void> const& other) :
			parameters{ (other.len > 0) ? new double[other.len] {} : nullptr }, len{ other.len }
		{
			std::copy(other.parameters, other.parameters + other.len, parameters);
		}

		parameterized_type(parameterized_type<void, void>&& other) : parameterized_type()
		{
			swap(*this, other);
		}

		parameterized_type<void, void>& operator=(parameterized_type<void, void> other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(parameterized_type<void, void>& first, parameterized_type<void, void>& second)
		{
			using std::swap;
			swap(first.parameters, second.parameters);
			swap(first.len, second.len);
		}

		const parameterized_type<void, void>& operator()() const
		{
			return *this;
		}

		parameterized_type<void, void>& operator()()
		{
			return *this;
		}

		double operator[](iter_type i) const
		{
			if (i >= len)
			{
				return DEFAULT_COEFF_VALUE;
			}
			else
			{
				return parameters[i];
			}
		}

		double* parameters;
		len_type len;


		~parameterized_type()
		{
			delete[] parameters;
		}

	protected:

		auto as_parameter(double parameter)
		{
			return parameter;
		}

		template<size_t I0, size_t... Is, typename T0, typename... Ts>
		void populate_parameters(std::index_sequence<I0, Is...>, T0&& parameter0, Ts&&... parameters)
		{
			this->parameters[I0] = as_parameter(std::forward<T0>(parameter0));
			((this->parameters[Is] = as_parameter(std::forward<Ts>(parameters))), ...);
		}

		inline void populate_parameters(std::index_sequence<>) {}
	};

	template<typename T, typename S>
	struct parameterized_type<parameterized_type<T, S>, void> : parameterized_type<T, S>
	{
		using parent_type = parameterized_type<T, S>;
		using parent_type::parent_type;
	};

	template<typename S>
	struct parameterized_type<S, void> : parameterized_type<void, void>
	{
		template<typename M>
		parameterized_type(M const& model) :
			parameterized_type<void>(model_field_parameters<M, S>(model)) {}

		parameterized_type() : parameterized_type<void>() {}

		const parameterized_type<void, void>& operator()() const
		{
			return *this;
		}

		parameterized_type<void, void>& operator()()
		{
			return *this;
		}
	};


	template<typename M, typename T>
	struct model_field_parameters : model_field_parameters<M, void>, parameterized_type<void, void>
	{
		using parent_type = model_field_parameters<M, void>;
		model_field_parameters(M const& m) : parent_type{ m }, parameterized_type<void, void>{ m->param(0) } {}
	};

	template<typename M>
	struct model_field_parameters<M, void> 
	{
		model_field_parameters(M const& m) : model{ &m } {}

		template<size_t I = 0>
		auto op() const
		{
			return model->template op<I>();
		}

		template<size_t I = 0>
		auto dop() const
		{
			return model->template dop<I>();
		}

		template<int N, int P>
		auto op_i(expr::symbols::i_<N, P>) const
		{
			return model->op_i(expr::symbols::i_<N, P>{});
		}

		auto param(size_t I) const
		{
			return model->param(I);
		}

		template<size_t N>
		auto param_matrix() const
		{
			return model->template param_matrix<N>();
		}

		template<expr::NoiseType nt, typename T, typename... T0s>
		auto make_noise(T0s&& ...args) const
		{
			return model->template make_noise<nt, T>(std::forward<T0s>(args)...);
		}

		template<expr::NoiseType nt, typename G, typename... T0s>
		auto make_noise(OpTerm<OpIdentity, G> const& term, T0s&& ...args) const
		{
			return model->template make_noise<nt>(term, std::forward<T0s>(args)...);
		}
		const M* model;
	};

	template<typename T>
	struct non_parameterized_type_impl
	{
		using type = T;
	};

	template<typename T, typename S>
	struct non_parameterized_type_impl<parameterized_type<T, S>>
	{
		using type = S;
	};

	template<typename T>
	using non_parameterized_type = typename non_parameterized_type_impl<T>::type;

}

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
				snprintf(name, BUFFER_LENGTH, ORDER_PARAMETER_NAME_EXTRA_FMT, i);
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



template<size_t D, typename Sp, typename S, typename... Ts>
using ArrayModel = Model<D, Sp, symphas::internal::field_array_t<S>, Ts...>;

template<size_t D, typename Sp, typename... S>
struct Model;


namespace symphas::internal
{
	template<typename T>
	struct field_array_t;
	
	template<typename S>
	constexpr bool not_model_array_type = true;
	template<typename S>
	constexpr bool not_model_array_type<field_array_t<S>> = false;
	template<typename T, typename S>
	constexpr bool not_model_array_type<parameterized_type<T, S>> = not_model_array_type<S>;
	
	template<ModelModifiers modifier = ModelModifiers::PLOT_DEFAULT>
	struct modifier_save
	{
		template<typename S>
		void operator()(S* _s, len_type len, iter_type index, const char* dir) const
		{
			if (dir)
			{
				for (iter_type i = 0; i < len; ++i)
				{
					_s[i].save(dir, index);
				}
			}
		}

		template<typename S>
		auto& get_output(S* _s) const
		{
			return _s[0];
		}
	};

	template<>
	struct modifier_save<ModelModifiers::PLOT_MAX>
	{
		template<typename S>
		void operator()(S* _s, len_type len, iter_type index, const char* dir) const
		{
			S& output = get_output(_s);
			if (dir)
			{
#				pragma omp parallel for
				for (iter_type n = 0; n < output.len; ++n)
				{
					output[n] = _s[0][n];
				}

				for (iter_type i = 1; i < len; ++i)
				{
#					pragma omp parallel for
					for (iter_type n = 0; n < output.len; ++n)
					{
						output[n] = std::max(output[n], _s[i][n]);
					}
				}
				output.save(dir, index);
			}
		}

		template<typename S>
		auto& get_output(S* _s) const
		{
			static S ss(_s[0]);
			return ss;
		}
	};



	template<size_t D, typename Sp, typename S, typename... Ts>
	auto num_fields_as_literal(ArrayModel<D, Sp, S, Ts...> const& model)
	{
		return expr::make_literal(model.len);
	}

	template<size_t D, typename Sp, typename... S>
	auto num_fields_as_literal(Model<D, Sp, S...> const& model)
	{
		return expr::val<sizeof...(S)>;
	}

}

namespace symphas
{

	template<size_t D, typename Sp, typename S, typename... Ts>
	auto model_num_fields(ArrayModel<D, Sp, S, Ts...> const& model)
	{
		return model.len;
	}

	template<size_t D, typename Sp, typename... S>
	auto model_num_fields(Model<D, Sp, S...> const& model)
	{
		return sizeof...(S);
	}
}


//! A representation of the phase field crystal model.
/*!
 * A representation of the phase field crystal model. Implements the phase
 * field crystal model equation, and allows 
 * 
 * \tparam PFC The phase field crystal parameters specialized object.
 * \tparam D The dimension of the phase field crystal problem.
 * \tparam Sp The solver type for numerically solving the phase field crystal
 * problem. 
 * \tparam S... The types of the order parameters.
 */
template<size_t D, typename Sp, typename S, typename... Ts>
struct Model<D, Sp, symphas::internal::field_array_t<S>, Ts...>
{

	//! Get the type of the phase-fields, takes an argument for comptability with other Model.
	/*!
	 * Get the type of the `N`th system, but since they are all the same type, will return
	 * the type directly regardless of the value of `N`.
	 *
	 * \tparam N The phase field index to get the type.
	 */
	template<size_t N = 0>
	using type_of_S = symphas::internal::non_parameterized_type<S>;


	//! The type of the system storing the phase fields, used by the solver.
	using SolverSystemApplied = typename symphas::solver_system_type<Sp>::template type<type_of_S<>, D>;



protected:

	Model()
		: len{ 0 }, _s{ construct_systems({}, {}, {}, 0, 0) }, solver{ Sp::make_solver() }, coeff{ nullptr }, num_coeff{ 0 },
		index{ params::start_index }, time{ 0 }, plot_type{ symphas::ModelModifiers::PLOT_DEFAULT }
#ifdef VTK_ON
		, viz_update{ nullptr }
#endif 
	{}

public:


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
		len{ static_cast<len_type>(parameters.length()) },
		_s{ construct_systems(parameters.get_initial_data(), parameters.get_interval_data(), parameters.get_boundary_data(), parameters.length(), parameters.length()) },
		solver{ Sp::make_solver(get_updated_parameters(parameters)) }, coeff{ (num_coeff > 0) ? new double[num_coeff] : nullptr },
		num_coeff{ num_coeff }, index{ parameters.index }, time{ parameters.time }, plot_type{ parameters.get_modifier() }
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

	Model(ArrayModel<D, Sp, S, Ts...> const& other) :
		len{ other.len },
		_s{ construct_systems(other._s, other.len) }, solver{other.solver}, coeff{(other.num_coeff > 0) ? new double[other.num_coeff] : nullptr},
		num_coeff{ other.num_coeff }, index{ other.index }, time{ other.time }, plot_type{ other.plot_type }
#ifdef VTK_ON
		, viz_update{ nullptr }
#endif
	{
		std::copy(other._s, other._s + other.len, _s);
		std::copy(other.coeff, other.coeff + other.num_coeff, coeff);
		visualize();
	}

	Model(ArrayModel<D, Sp, S, Ts...>&& other) noexcept : Model()
	{
		swap(*this, other);
	}

	ArrayModel<D, Sp, S, Ts...>& operator=(ArrayModel<D, Sp, S, Ts...> other)
	{
		swap(*this, other);
		return *this;
	}


	~Model()
	{
		for (int i = len - 1; i >= 0; --i)
		{
			_s[i].~SolverSystemApplied();
		}
		operator delete[](_s);

		devisualize();
		delete[] coeff;
	}

	friend void swap(ArrayModel<D, Sp, S, Ts...>& first, ArrayModel<D, Sp, S, Ts...>& second)
	{
		using std::swap;


		swap(first.solver, second.solver);
		swap(first._s, second._s);
		swap(first.coeff, second.coeff);
		swap(first.len, second.len);
		swap(first.num_coeff, second.num_coeff);
		swap(first.index, second.index);
		swap(first.time, second.time);
		swap(first.plot_type, second.timeplot_type);

#ifdef VTK_ON
		swap(first.viz_thread, second.viz_thread);
		swap(first.viz_update, second.viz_update);
#endif
	}

	template<typename Sp0, typename = std::enable_if_t<!std::is_same<Sp, Sp0>::value, int>>
	friend void swap(ArrayModel<D, Sp, S, Ts...>& first, ArrayModel<D, Sp0, S, Ts...>& second)
	{
		using std::swap;

		swap(first._s, second._s);
		swap(first.coeff, second.coeff);
		swap(first.len, second.len);
		swap(first.num_coeff, second.num_coeff);
		swap(first.index, second.index);
		swap(first.time, second.time);
		swap(first.plot_type, second.timeplot_type);

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
	auto get_time_var() const
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
	size_t num_fields() const
	{
		if constexpr (std::is_same<S0, type_of_S<>>::value)
		{
			return len;
		}
		else
		{
			return 0;
		}
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
	size_t num_fields() const
	{
		return num_fields<S0>() + num_fields<S1, Ss...>();
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
	//template<typename Type, size_t I = 0>
	//int index_of_type = (std::is_same<Type, S>::value) ? I : -1;

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
	template<typename Type, typename F, typename... Args, typename std::enable_if<std::is_same<Type, type_of_S<>>::value, int>::type = 0>
	void do_for_field_type(F f, Args&& ...args) const
	{
		do_for_field_type<Type, 0>(f, std::forward<Args>(args)...);
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
	template<size_t N, typename F, typename... Args>
	void do_for_field(F f, Args&& ... args) const
	{
		auto&& data = _s[N].get_snapshot();
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
	void copy_field_values(type_of_S<>* into) const
	{
		_s[N].persist(into);
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
	template<typename Type, size_t I, typename std::enable_if_t<std::is_same<Type, type_of_S<>>::value, int> = 0>
	void copy_field_type_values(Type* into) const
	{
		copy_field_values<I>(into);
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
	template<typename Type, size_t I, typename std::enable_if_t<!std::is_same<Type, type_of_S<>>::value, int> = 0>
	void copy_field_type_values(Type*) const {}

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
		Grid<type_of_S<>, D> out(_s[N].get_info().get_dims().get());
		_s[N].persist(out.values);
		return out;
	}

	//! Persists all the systems managed by the model.
	/*!
	 * Persists the phase field systems to disk, saving the entire array.
	 * This function has no effect without the **io** module.
	 *
	 * \param dir The directory into which to store the resulting data file.
	 */
	void save_systems(const char* dir = nullptr) const
	{
		if (plot_type == symphas::ModelModifiers::PLOT_MAX)
		{
			symphas::internal::modifier_save<symphas::ModelModifiers::PLOT_MAX>{}(_s, len, index, dir);
		}
		else
		{
			symphas::internal::modifier_save<>{}(_s, len, index, dir);
		}
		iter_type n = len - sizeof...(Ts);
		for (iter_type i = 0; i < sizeof...(Ts); ++i)
		{
			_s[i + n].save(dir, i + 1);
		}
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
		if (plot_type == ModelModifiers::PLOT_MAX)
		{
			symphas::internal::modifier_save<ModelModifiers::PLOT_MAX>{}(_s, len, index, dir);
		}
		else
		{
			symphas::internal::modifier_save<>{}(_s, len, index, dir);
		}
		iter_type n = len - sizeof...(Ts);
		for (iter_type i = 0; i < sizeof...(Ts); ++i)
		{
			_s[i + n].save(dir, index);
		}
		//char* name = new char[std::strlen(name) + symphas::lib::num_digits(len) + 1];
		//for (iter_type i = 0; i < len; ++i)
		//{
		//	snprintf(name, BUFFER_LENGTH, "%d%s", i, name);
		//	_s[i].save(dir, name, index)));
		//}
		//delete[] name;
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
		for (iter_type i = 0; i < len; ++i)
		{
			_s[i].update(index, time);
		}
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
		solver.dt = dt;
		for (iter_type i = 0; i < len; ++i)
		{
			solver.step(system(i));
		}

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
		return _s[I];
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
		return const_cast<SolverSystemApplied&>(static_cast<const ArrayModel<D, Sp, S, Ts...>&>(*this).system<I>());
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

	//! Returns the system at the given index.
	/*!
	 * Returns the reference to the system at the given index from the list of
	 * phase field systems of the model.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	const auto& system(iter_type i) const
	{
		return _s[i];
	}

	//! Returns the system at the given index.
	/*!
	 * Returns the reference to the system at the given index from the list of
	 * phase field systems of the model.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	auto& system(iter_type i)
	{
		return _s[i];
	}

	//! Returns the underlying grid of the system at the given index.
	/*!
	 * Returns a reference to the grid of the system at the given index.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	auto& grid(iter_type i)
	{
		return system(i).as_grid();
	}

	//! Returns the underlying grid of the system at the given index.
	/*!
	 * Returns a reference to the grid of the system at the given index.
	 *
	 * \tparam I The index of the phase field to return.
	 */
	const auto& grid(iter_type i) const
	{
		return system(i).as_grid();
	}


	//! Returns the tuple of all model systems.
	/*!
	 * Directly returns the list of all systems of this model.
	 */
	decltype(auto) systems_tuple() const
	{
		return std::make_pair(_s, len);
	}

	//! Returns the tuple of all model systems.
	/*!
	 * Directly returns the list of all systems of this model.
	 */
	decltype(auto) systems_tuple()
	{
		return std::make_pair(_s, len);
	}

	//! Returns the tuple of all model systems.
	/*!
	 * Directly returns the list of all systems of this model.
	 */
	const auto& systems() const
	{
		return _s;
	}

	//! Returns the tuple of all model systems.
	/*!
	 * Directly returns the list of all systems of this model.
	 */
	auto& systems()
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
		fprintf(out, "Phase field problem of %d system%s:\n", len, (len > 1) ? "s" : "");
		print_all_dimensions(out);
		print_all_coeff(out);
		fprintf(out, OUTPUT_BANNER);
	}

	auto generate_parameters() const
	{
		symphas::problem_parameters_type parameters(len);
		populate_parameters(parameters);

		parameters.time = get_time();;
		parameters.index = index;

		return parameters;
	}

	len_type len;						//!< Number of fields.

protected:
	SolverSystemApplied* _s;			//!< Container managing pointers to phase field data.

	Sp solver;							//! Solver for determining the phase field solution.

	//! Coefficients in the equations of motion.
	/*!
	 * List of coefficients used in the equations of motion. If
	 * the number of coefficients is insufficient, then the default coefficient
	 * value #DEFAULT_COEFF_VALUE should be used.
	 */
	double* coeff;

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

	symphas::ModelModifiers plot_type;

#ifdef VTK_ON

	std::thread viz_thread;
	ColourPlotUpdater* viz_update;


	void visualize()
	{
		if constexpr (std::is_same<type_of_S<>, scalar_t>::value)
		{
			S* ss;
			if (plot_type == symphas::ModelModifiers::PLOT_MAX)
			{
				ss = &symphas::internal::modifier_save<symphas::ModelModifiers::PLOT_MAX>{}.get_output(_s);
			}
			else
			{
				ss = &symphas::internal::modifier_save<>{}.get_output(_s);
			}

			if (params::viz_interval > 0)
			{
				viz_thread = std::thread([] (auto* viz_grid, int* index, auto** viz_update)
					{
						ColourPlot2d viz{};
						viz.init(viz_grid->values, viz_grid->dims, *index, *viz_update);
					}, &(ss->as_grid()), & index, & viz_update);
			}
		}
	}

	void devisualize()
	{
		if constexpr (std::is_same<type_of_S<>, scalar_t>::value)
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

	auto get_updated_parameters(symphas::problem_parameters_type parameters) const
	{
		parameters.extend(len);
		for (iter_type i = 0; i < len; ++i)
		{
			parameters.set_interval_data(system(i).get_info().intervals, i);
		}
		return parameters;
	}

	template<size_t... Is>
	auto populate_parameters(symphas::problem_parameters_type& parameters) const
	{
		for (iter_type i = 0; i < len; ++i)
		{
			fill_interval_data(system(i), parameters, i);
			fill_boundary_data(system(i), parameters, i);
			fill_initial_data(system(i), parameters, i);
		}
	}

	void fill_interval_data(PhaseFieldSystem<BoundaryGrid, type_of_S<>, D> const& system, symphas::problem_parameters_type& parameters, iter_type n) const
	{
		auto intervals = system.get_info().intervals;
		for (iter_type i = 0; i < D; ++i)
		{
			Axis ax = symphas::index_to_axis(i);
			auto& interval = intervals.at(ax);

			interval.set_count(
				interval.left(),
				interval.right(),
				system.dims[i]);
		}

		parameters.set_interval_data(intervals, n);
	}

	template<template<typename, size_t> typename G>
	void fill_interval_data(PhaseFieldSystem<G, type_of_S<>, D> const& system, symphas::problem_parameters_type& parameters, iter_type n) const
	{
		parameters.set_interval_data(system.get_info().intervals, n);
	}

	template<template<typename, size_t> typename G, typename T0>
	void fill_boundary_data(PhaseFieldSystem<G, T0, D> const& system, symphas::problem_parameters_type& parameters, iter_type n) const
	{
		symphas::b_data_type bdata;

		for (iter_type i = 0; i < D * 2; ++i)
		{
			Side side = symphas::index_to_side(i);
			bdata[side] = symphas::b_element_type(BoundaryType::PERIODIC);
		}

		parameters.set_boundary_data(bdata, n);
	}

	template<typename T0>
	void fill_boundary_data(PhaseFieldSystem<BoundaryGrid, T0, D> const& system, symphas::problem_parameters_type& parameters, iter_type n) const
	{
		symphas::b_data_type bdata;

		for (iter_type i = 0; i < D * 2; ++i)
		{
			Side side = symphas::index_to_side(i);
			BoundaryType type = system.types[i];
			bdata[side] = system.boundaries[i]->get_parameters();
		}

		parameters.set_boundary_data(bdata, n);
	}

	template<template<typename, size_t> typename G, typename T0>
	void fill_initial_data(PhaseFieldSystem<G, T0, D> const& system, symphas::problem_parameters_type& parameters, iter_type n) const
	{
		parameters.set_initial_data(init_data_functor(system.as_grid()), n);
	}

	auto construct_systems(const symphas::init_data_type* tdata, const symphas::interval_data_type* vdata, const symphas::b_data_type* bdata, size_t data_len, size_t len) const
	{
		void* systems = operator new[](len * sizeof(SolverSystemApplied));
		SolverSystemApplied* ptr = static_cast<SolverSystemApplied*>(systems);
		
		for (iter_type i = 0; i < std::min(data_len, len); ++i)
		{
			new(&ptr[i]) SolverSystemApplied(tdata[i], vdata[i], bdata[i], i);
		}
		for (iter_type i = std::min(data_len, len); i < len; ++i)
		{
			new(&ptr[i]) SolverSystemApplied(tdata[0], vdata[0], bdata[0], i);
		}

		return ptr;
	}

	template<typename other_sys_type, typename std::enable_if_t<!std::is_same<other_sys_type, SolverSystemApplied>::value, int> = 0>
	SolverSystemApplied construct_system(other_sys_type const& system, iter_type n) const
	{
		symphas::b_data_type boundaries{};
		for (iter_type i = 0; i < D * 2; ++i)
		{
			boundaries[symphas::index_to_side(i)] = BoundaryType::PERIODIC;
		}

		bool extend_boundary = params::extend_boundary;
		params::extend_boundary = true;

		auto new_system = SolverSystemApplied(
			symphas::init_data_type{}, system.info.intervals, boundaries, system.id);

		params::extend_boundary = extend_boundary;

		type_of_S<>* swap_arr = new type_of_S<>[system.length()];
		system.persist(swap_arr);
		new_system.fill(swap_arr);

		delete[] swap_arr;

		return new_system;
	}


	template<typename S0>
	SolverSystemApplied* construct_systems(S0* const& systems, len_type len) const
	{
		if (len > 0)
		{
			void* new_systems = operator new[](len * sizeof(SolverSystemApplied));
			SolverSystemApplied* ptr = static_cast<SolverSystemApplied*>(new_systems);

			for (iter_type i = 0; i < len; ++i)
			{
				new(&ptr[i]) SolverSystemApplied(systems[i]);
			}
			return ptr;
		}
		else
		{
			return nullptr;
		}
	}

	/*
	 * for the I-th type field in the whole list of fields, execute
	 * the given function
	 */

	template<typename Type, size_t I, typename F, typename... Args, typename std::enable_if_t<(I > 0), int> = 0>
	void do_for_field_type(F&& f, Args&& ... args) const
	{
		do_for_field_type<Type, I - 1, Args...>(f, std::forward<Args>(args)...);
	}

	template<typename Type, size_t I, typename F, typename... Args, typename std::enable_if_t<(I == 0), int> = 0>
	void do_for_field_type(F&& f, Args&& ... args) const
	{
		do_for_field<I>(f, std::forward<Args>(args)...);
	}

	template<size_t... Is>
	void print_all_dimensions(FILE* out, std::index_sequence<Is...>) const
	{
		((..., print_grid_dimensions<Is>(out)));
	}

	void print_all_dimensions(FILE* out) const
	{
		fprintf(out, "\tsystem<*> ... ");
		for (iter_type d = 0; d < D; ++d)
		{
			fprintf(out, "%d", grid(0).dims[d]);
			if (d < D - 1)
			{
				fprintf(out, " x ");
			}
		}
		fprintf(out, "\n");
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










