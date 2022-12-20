
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
 * PURPOSE: Defines the solver object and workflow. The Solver object
 * is used in creating specialized solvers.
 *
 * ***************************************************************************
 */

#pragma once

#include <tuple>

#include "expressiontransforms.h"
#include "systemlib.h"

template<typename Sp>
struct Solver;




namespace symphas::internal
{
	using solver_unsupported_equation = int;

	template<typename S, typename T>
	struct solver_supported_type_match;

	//! Allows a system to be used by the solver.
	/*! See solver_supported_type.
	 */
	template<typename T, typename S>
	struct solver_supported_system;

	//! Allows a system to be used by the solver.
	/*! 
	 * See solver_supported_type.
	 */
	template<typename T, typename... S>
	struct solver_supported_systems;


	template<typename T>
	struct solver_system_type_match;

	template<typename T>
	struct provisional_system_type_match;


}

namespace symphas
{

	//! Allows a type to be used by the solver.
	/*!
	 * A model which implements a type not compatible with a solver may not
	 * correctly compile, so compile time constant information must be
	 * available whether to proceed with building that model. If a type
	 * is not supported, the set of equations will not be processed.
	 */
	template<typename T, typename Ty>
	struct solver_supported_type
	{
		static const bool value = symphas::internal::solver_supported_type_match<typename T::id_type, Ty>::value;
	};

	//! Associates a solver system type to a solver specialization.
	/*!
	 * When a new solver is defined, a different solver system type may be
	 * desired to support its functionality.
	 *
	 * Solver systems are typically used to store temporary information about
	 * the solution in order to time evolve the phase fields on the next
	 * iteration, and we may wish to associate a particular solver system
	 * implementation with the solver we are using.
	 */
	template<typename T>
	struct solver_system_type
	{
		template<typename Ty, size_t D>
		using type = typename symphas::internal::solver_system_type_match<typename T::id_type>::template type<Ty, D>;
	};

	//! Associates a provisional system type to a solver specialization.
	/*!
	 * When a new solver is defined, a different provisional system type may be
	 * desired to support its functionality.
	 *
	 * Provisional systems are typically used to store the output of provisional
	 * variables, and at the same time, we may with to add additional information
	 * about provisional variables in order to support their implementation in
	 * the solver.
	 */
	template<typename T>
	struct provisional_system_type
	{
		template<typename Ty, size_t D>
		using type = typename symphas::internal::provisional_system_type_match<typename T::id_type>::template type<Ty, D>;
	};
}

/*!
 * \defgroup solver Solver Definition
 * Every new solver inherits from Solver through CRTP, and has to implement a 
 * set of functions in order to have a complete work flow.

 * Firstly, there are the following functions which allow the derivatives to be 
 * computed, and these are _optionally_ implemented, as Solver
 * implements these already:
 *
 * ```cpp
 * template<size_t O, typename T>
 * auto applied_generalized_derivative(Block<T> const& e, iter_type n) const
 * 
 * template<typename T>
 * auto applied_laplacian(Block<T> const& e, iter_type n) const
 * 
 * template<typename T>
 * auto applied_bilaplacian(Block<T> const& e, iter_type n) const
 * 
 * template<typename T>
 * auto applied_gradlaplacian(Block<T> const& e, iter_type n) const
 * 
 * template<typename T>
 * auto applied_gradient(Block<T> const& e, iter_type n) const
 * ```
 *
 *
 * For the core solver functionality, three functions must be implemented,
 * `form_expr_one`, `equation` and `step`. These are all `const` functions.
 * 
 * The first function will be given an equation of motion and will be required
 * to create an object from this equation that will be used to update the
 * order parameter in the next two functions. 
 * The signature of this function is:
 * 
 * ```cpp
 * template<size_t En, typename SS, typename S, typename E>
 * auto form_expr_one(SS&&, std::pair<S, E>&& e) const
 * ```
 * 
 * The template parameter `En` is the index of the order parameter for which
 * the equation of motion was provided.
 * 
 * The first argument is the list of all systems, including provisional
 * systems which have been concatenated to the list of phase field systems.
 * The second argument will always been an equation of motion that is
 * provided as an `std::pair` type, where the first entry
 * is the phase field system, and the second entry is the expression
 * representing the RHS of the equation of motion. 
 * 
 * The phase field system type
 * is ::SolverSystem, specialized on the order parameter type and model 
 * dimension. This can be another type, such as the specialized systems
 * SolverSystemFD or SolverSystemSpectral that are used for ::SolverFT and ::SolverSP,
 * respectively, if it is added through #ASSOCIATE_SOLVER_SYSTEM_TYPE. It is
 * also possible to associate a custom system, but this has to be
 * included in the *SymPhas* installation beforehand.
 * 
 * The `form_expr_one` function will be run for each equation of motion that
 * is provided, therefore an output will be produced for each order parameter.
 * 
 * If more information beyond the full system list and equation of motion of the
 * current order parameter is required for the solver
 * implementation, then the function `form_expr` may be overloaded as well,
 * and delegate to `form_expr_one`. It has the signature:
 * 
 * ```cpp
 * template<size_t En, typename SS, typename... Es>
 * decltype(auto) form_expr(SS&& s, Es&& ...es) const
 * ```
 * 
 * where `es...` are all the equations of motion.
 * 
 * The next required function is `equation`, and 
 * will take the result of the `form_expr_one`
 * function and use it to update data, possibly stored in the phase field
 * system, in order to prepare for the next function, `step`. The 
 * implementation in this function is not restricted. 
 * The signature of `equation` is:
 * 
 * ```cpp
 * template<typename S, typename E>
 * void equation(std::pair<S, E>& r) const
 * ```
 * 
 * The last function will step the solution through time, that is, it has
 * the effect of updating the phase field systems with values of the
 * next time index. The function signature is:
 * 
 * ```cpp
 * template<typename S>
 * void step(S& sys) const
 * ```
 * 
 * The type `S` is the phase field system, depending on the type that
 * has been associated with this solver as described above. The other parameter
 * is the time step, which is optionally used. 
 *
 * Lastly, a function is required in order to construct the specialized solver
 * objects using general program parameters, called by the model class in order
 * to create the solver:
 *
 * ```cpp
 * static auto&& make_solver(problem_parameters_type const& parameters)
 * ```
 *
 * @{
 * 
 */

//! Associate the solver system type to a specialized solver.
/*!
 * The solver system will contain the phase field data and be used by
 * the solver in the time evolution of the phase fields.
 */
#define ASSOCIATE_SOLVER_SYSTEM_TYPE(SOLVER_NAME, SYSTEM_TYPE) \
template<> struct symphas::internal::solver_system_type_match<solver_id_type_ ## SOLVER_NAME> \
{ template<typename Ty, size_t D> using type = SYSTEM_TYPE<Ty, D>; };

//! Associate the provisional system type to a specialized solver.
/*!
 * The provisional system is used to store variables, and a particular
 * implementation may be desired which supports the solver in some way.
 */
#define ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SOLVER_NAME, SYSTEM_TYPE) \
template<> struct symphas::internal::provisional_system_type_match<solver_id_type_ ## SOLVER_NAME> \
{ template<typename Ty, size_t D> using type = SYSTEM_TYPE<Ty, D>; };

//! Adds one type to the type of solvable phase field variable.
/*!
 * Some solvers may not be compatible with all types of values of phase fields,
 * such as vector types. This command is used to add a particular type
 * to the allowed phase field value-types which can be solved. Attempting
 * to solve a phase field of value-type that is not allowed will immediately
 * finish the simulation.
 */
#define SYMPHAS_SOLVER_SUPPORTED_TYPE(SOLVER_NAME, TYPE) \
template<> \
struct symphas::internal::solver_supported_type_match<solver_id_type_ ## SOLVER_NAME, TYPE> \
{ \
	static const bool value = true; \
};

//! Add all value-types to the types that can be solved by the solver.
/*!
 * See #SYMPHAS_SOLVER_SUPPORTED_TYPE.
 */
#define SYMPHAS_SOLVER_ALL_SUPPORTED(SOLVER_NAME) \
template<> \
struct symphas::internal::solver_supported_type_match<solver_id_type_ ## SOLVER_NAME, void> \
{ \
	static const bool value = true; \
};


//! Define a new solver of the given name.
/*!
 * Definition of a new solver, which is a template class with a template type
 * that is typically the finite difference stencil. The new solver is given
 * the name provided.
 * 
 * To define a custom solver, see \ref solver.
 * 
 * \param NAME The name of the new solver.
 */
#define NEW_SOLVER_WITH_STENCIL(NAME) \
struct solver_id_type_ ## NAME; \
template<typename stencil_t> \
struct NAME : Solver<NAME<stencil_t>>, stencil_t \
{ \
	using this_type = NAME<stencil_t>; \
	using parent_type = Solver<this_type>; \
	using parent_type::parent_type; \
	using parent_type::generalized_derivative; \
	using parent_type::laplacian; \
	using parent_type::bilaplacian; \
	using parent_type::gradlaplacian; \
	using parent_type::gradient; \
	using stencil_t::applied_generalized_derivative; \
	using stencil_t::applied_laplacian; \
	using stencil_t::applied_bilaplacian; \
	using stencil_t::applied_gradlaplacian; \
	using stencil_t::applied_gradient; \
	using parent_type::make_solver; \
	using id_type = solver_id_type_ ## NAME; \
	using stencil_t::stencil_t; \
	using parent_type::dt; \
	stencil_t const& stencil() const { return (*static_cast<stencil_t const*>(this)); } \
	stencil_t& stencil() { return (*static_cast<stencil_t*>(this)); }

#define NEW_SOLVER(NAME, ...) \
struct solver_id_type_ ## NAME; \
struct NAME : Solver<NAME>, NoStencil, __VA_ARGS__ \
{ \
	using parent_type = Solver<NAME>; \
	using parent_type::generalized_derivative; \
	using parent_type::laplacian; \
	using parent_type::bilaplacian; \
	using parent_type::gradlaplacian; \
	using parent_type::gradient; \
	using NoStencil::applied_generalized_derivative; \
	using NoStencil::applied_laplacian; \
	using NoStencil::applied_bilaplacian; \
	using NoStencil::applied_gradlaplacian; \
	using NoStencil::applied_gradient; \
	using parent_type::make_solver; \
	using id_type = solver_id_type_ ## NAME;



struct NoStencil
{
	template<Axis ax, size_t O, typename G>
	auto applied_generalized_derivative(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<Axis ax, size_t O, typename G>
	auto applied_generalized_directional_derivative(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<size_t... Os, typename G>
	auto applied_generalized_mixed_derivative(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<typename G>
	auto applied_laplacian(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<typename G>
	auto applied_bilaplacian(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<Axis ax, typename G>
	auto applied_gradlaplacian(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

	template<Axis ax, typename G>
	auto applied_gradient(G&& e, iter_type n) const
	{
		return expr::BaseData<G>::get(std::forward<G>(e), n);
	}

};

//! @}



//! Base class for the solver functionality.
/*!
 * Defines the structure and function of a solver. Specialized solvers use
 * this class through CRTP inheritance.  
 * 
 * To define a custom solver, see \ref solver.
 *
 * \param Sp The specialized solver.
 */
template<typename Sp>
struct Solver
{

	double dt;

	Solver(double dt = 0) : dt{ dt } {}

	//! The object which computes the derivative.
	template<Axis ax, size_t O>
	struct derivative;

	//! The object which computes the directional derivative.
	template<Axis ax, size_t O>
	struct directional_derivative;

	//! The object which computes the directional derivative.
	template<size_t... Os>
	struct mixed_derivative;

	template<Axis ax, size_t O, typename G>
	decltype(auto) generalized_derivative(G&& e, iter_type n) const
	{
		return cast().template applied_generalized_derivative<ax, O>(std::forward<G>(e), n);
	}

	template<Axis ax, size_t O, typename G>
	decltype(auto) generalized_directional_derivative(G&& e, iter_type n) const
	{
		return cast().template applied_generalized_directional_derivative<ax, O>(std::forward<G>(e), n);
	}

	template<size_t... Os, typename G>
	decltype(auto) generalized_mixed_derivative(G&& e, iter_type n) const
	{
		return cast().template applied_generalized_mixed_derivative<Os...>(std::forward<G>(e), n);
	}

	template<typename G>
	decltype(auto) laplacian(G&& e, iter_type n) const
	{
		return cast().applied_laplacian(std::forward<G>(e), n);
	}

	template<typename G>
	decltype(auto) bilaplacian(G&& e, iter_type n) const
	{
		return cast().applied_bilaplacian(std::forward<G>(e), n);
	}

	template<Axis ax, typename G>
	decltype(auto) gradlaplacian(G&& e, iter_type n) const
	{
		return cast().template applied_gradlaplacian<ax>(std::forward<G>(e), n);
	}

	template<Axis ax, typename G>
	decltype(auto) gradient(G&& e, iter_type n) const
	{
		return cast().template applied_gradient<ax>(std::forward<G>(e), n);
	}





	//! Move forward one solution iteration.
	/*!
	 * Computes the phase field data for the next time step and updates
	 * the given object, which is typically the phase field system. Depending
	 * on the specific solver implementation, the time step provided here will
	 * not be used, only the time step provided with the solver construction
	 * will be used.
	 * 
	 * \param ss The system which is stepped forward.
	 */
	template<typename S>
	void step(S&& ss) const
	{
		cast().template step(std::forward<S>(ss));
	}

	//! Evaluate the equations.
	/*!
	 * The equations are evaluated, typically updating the data in the 
	 * specialized phase field systems that are defined with the solver. The
	 * data that is provided is not necessarily an expression, and comes
	 * from what was computed in the form_expr() function.
	 * 
	 * \param rs A list of objects that are used in computing part of the
	 * solution at this solution stage.
	 */
	template<typename... Rs>
	void equations(std::tuple<Rs...>& rs) const
	{
		equations_apply(rs, std::make_index_sequence<std::tuple_size<std::tuple<Rs...>>::value>{});
	}

	//! Evaluate the equations.
	/*!
	 * The equations are evaluated, typically updating the data in the
	 * specialized phase field systems that are defined with the solver. The
	 * data that is provided is not necessarily an expression, and comes
	 * from what was computed in the form_expr() function.
	 *
	 * \param rs A list of objects that are used in computing part of the
	 * solution at this solution stage.
	 */
	template<typename... Rs>
	void equations(std::tuple<Rs...>&& rs) const
	{
		equations_apply(rs, std::make_index_sequence<std::tuple_size<std::tuple<Rs...>>::value>{});
	}

	//! Evaluate the equations.
	/*!
	 * The equations are evaluated, typically updating the data in the
	 * specialized phase field systems that are defined with the solver. The
	 * data that is provided is not necessarily an expression, and comes
	 * from what was computed in the form_expr() function.
	 *
	 * \param rs... A list of objects that are used in computing part of the
	 * solution at this solution stage.
	 */
	template<typename... Rs>
	void equations(Rs&& ... rs) const
	{
		((..., cast().equation(std::forward<Rs>(rs))));
	}

	//! Evaluate the equations.
	/*!
	 * For unsupported types, OpVoid will be returned for the equation that
	 * needs to be processed. When OpVoid is given, then no action is taken
	 * (since it is equivalently zero).
	 */
	void equations(symphas::internal::solver_unsupported_equation) {}


	//! Evaluate the given list of elements, typically equation/data pairs.
	/*!
	 * Each equation from the list is evaluated, and the result is
	 * stored by the corresponding left hand side data.
	 */
	template<typename... Rs>
	inline void evaluate(std::tuple<Rs...>& rs) const
	{
		cast().evaluate_apply(rs, std::make_index_sequence<sizeof...(Rs)>{});
	}

	//! Evaluate the given list of elements, typically equation/data pairs.
	/*!
	 * Each equation from the list is evaluated, and the result is
	 * stored by the corresponding left hand side data.
	 */
	template<typename... Rs>
	inline void evaluate(std::tuple<Rs...>&& rs) const
	{
		cast().evaluate_apply(rs, std::make_index_sequence<sizeof...(Rs)>{});
	}

	//! Evaluate the given list of elements, typically equation/data pairs.
	/*!
	 * Each equation from the list is evaluated, and the result is
	 * stored by the corresponding left hand side data.
	 */
	template<typename... Rs>
	inline void evaluate(Rs&&... rs) const
	{
		((..., cast().evaluate_one(std::forward<Rs>(rs))));
	}


	//! Evaluate the given list of elements, typically equation/data pairs.
	/*!
	 * Each equation from the list is evaluated, and the result is
	 * stored by the corresponding left hand side data.
	 */
	template<typename R>
	inline void equation(R&& r) const
	{
		cast().equation(std::forward<R>(r));
	}


	//! Create the equation which is evaluated by the solver iteratively.
	/*!
	 * Solvers need to parse the expression and return another expression
	 * or tuple of expressions that is amenable to their solution method.
	 * The result is then ingested by equations().
	 *
	 * \param systems The list of data fields used in the phase field problem.
	 * \param es The list of equations.
	 *
	 * \tparam En The number of equation systems, the rest are supporting
	 * provisional systems.
	 */
	template<size_t En, typename... S, typename... Es, 
		typename std::enable_if_t<symphas::internal::solver_supported_systems<Sp, S...>::value, int> = 0>
	decltype(auto) form_expr_all(std::tuple<S...> const& systems, Es&& ...es) const
	{
		return cast().template form_expr<En>(systems, std::forward<Es>(es)...);
	}

	//! No equation is created for incompatible types.
	/*!
	 * For order parameter types types which are not compatible with the 
	 * solver, a special type is returned. When given to the equations(), 
	 * there is no effect. This means that models can be defined with any
	 * solver, and there will not be a compiler error if a type is not 
	 * supported.
	 */
	template<size_t En, typename... S, typename... Es, 
		typename std::enable_if_t<!symphas::internal::solver_supported_systems<Sp, S...>::value, int> = 0>
	decltype(auto) form_expr_all(std::tuple<S...> const&, Es&& ...) const
	{
		return symphas::internal::solver_unsupported_equation{};
	}

	//! Create the equation which is evaluated by the solver iteratively.
	/*!
	 * Solvers need to parse the expression and return another expression
	 * or tuple of expressions that is amenable to their solution method.
	 * The result is then ingested by equations().
	 *
	 * \param systems The list of data fields used in the phase field problem.
	 * \param es The list of equations.
	 *
	 * \tparam En The number of equation systems, the rest are supporting
	 * provisional systems.
	 */
	template<size_t En, typename... S, typename... Es, 
		typename std::enable_if_t<symphas::internal::solver_supported_systems<Sp, S...>::value, int> = 0>
	decltype(auto) form_expr_all(std::tuple<S const&...>&& systems, Es&& ...es) const
	{
		return cast().template form_expr<En>(systems, std::forward<Es>(es)...);
	}

	//! Create the equation which is evaluated by the solver iteratively.
	/*!
	 * For order parameter types types which are not compatible with the 
	 * solver, a special type is returned. When given to the equations(), 
	 * there is no effect. This means that models can be defined with any
	 * solver, and there will not be a compiler error if a type is not 
	 * supported.
	 */
	template<size_t En, typename... S, typename... Es, 
		typename std::enable_if_t<!symphas::internal::solver_supported_systems<Sp, S...>::value, int> = 0>
	decltype(auto) form_expr_all(std::tuple<S const&...>&& systems, Es&& ...) const
	{
		return symphas::internal::solver_unsupported_equation{};
	}

	//! Generate a specialized solver instance from problem parameters.
	/*!
	 * Takes the problem parameters and passes
	 * the correct information to the solver constructor in order to create 
	 * the specialized solver instance.
	 * 
	 * \param parameters The problem parameters used to construct the solver.
	 */
	static decltype(auto) make_solver(symphas::problem_parameters_type const& parameters)
	{
		return Sp::make_solver(parameters);
	}

	//! Generate a default specialized solver.
	/*!
	 * In many cases, it is desirable for this function to be implemented 
	 * in conjunction with make_solver(symphas::problem_parameters_type const&),
	 * to avoid any bugs which may occur due to undefined parameters.
	 */
	static decltype(auto) make_solver()
	{
		return Sp::make_solver({0});
	}

	//! Cast the solver to the specialized type.
	/*!
	 * Cast the solver to the specialized type.
	 */
	Sp const& cast() const
	{
		return *static_cast<Sp const*>(this);
	}

	//! Cast the solver to the specialized type.
	/*!
	 * Cast the solver to the specialized type.
	 */
	Sp& cast()
	{
		return *static_cast<Sp*>(this);
	}

	friend void swap(Solver<Sp>& first, Solver<Sp>& second)
	{
		using std::swap;
		swap(*static_cast<Sp*>(&first), *static_cast<Sp*>(&second));
	}

	//! Indirection used when the expression may be invalid.
	template<size_t En, typename SS, typename... Es>
	decltype(auto) form_expr(SS&& s, Es&& ...es) const
	{
		return std::make_tuple(cast().template form_expr_one<En>(std::forward<SS>(s), std::forward<Es>(es))...);
	}

	//! Applied implementation of evaluating a list of elements.
	/*!
	 * forwards to Solver::equation, but needs a different name so that
	 * it avoids disambiguation when using std::apply to forward tuples.
	 */
	template<typename... Rs, size_t... Is>
	inline void equations_apply(std::tuple<Rs...>& rs, std::index_sequence<Is...>) const
	{
		((..., equations(std::get<Is>(rs))));
	}


	//! Applied implementation of evaluating a list of elements.
	/*!
	 * See Solver::_equations.
	 */
	template<typename... Rs, size_t... Is>
	inline void evaluate_apply(std::tuple<Rs...>& rs, std::index_sequence<Is...>) const
	{
		((..., cast().evaluate_one(std::get<Is>(rs))));
	}
	template<typename... Rs, size_t... Is>
	inline void evaluate_apply(std::tuple<Rs...>&& rs, std::index_sequence<Is...>) const
	{
		((..., cast().evaluate_one(std::get<Is>(rs))));
	}


	//! Evaluate the result of the equation into the given data.
	/*!
	 * Evaluate the equation into the given data.
	 *
	 * \param r The pair consisting of the data and equation.
	 */
	template<typename G, typename E, typename std::enable_if_t<expr::has_state<E>::value, int> = 0>
	void evaluate_one(std::pair<G, E>& r) const
	{
		auto&& [grid, equation] = r;
		expr::prune::update(equation);
		expr::result(equation, expr::BaseData<G>::get(grid));
	}

	template<typename G, typename E, typename std::enable_if_t<!expr::has_state<E>::value, int> = 0>
	void evaluate_one(std::pair<G, E>& r) const
	{
		auto&& [grid, equation] = r;
		expr::result(equation, expr::BaseData<G>::get(grid));
	}

	//! Evaluate the result of the equation into the given data.
	/*!
	 * Evaluate the equation into the given data.
	 *
	 * \param r The pair consisting of the data and equation.
	 */
	template<typename G, typename E, typename std::enable_if_t<expr::has_state<E>::value, int> = 0>
	void evaluate_one(std::pair<G, E>&& r) const
	{
		auto& [grid, equation] = r;
		expr::prune::update(equation);
		expr::result(equation, expr::BaseData<G>::get(grid));
	}

	template<typename G, typename E, typename std::enable_if_t<!expr::has_state<E>::value, int> = 0>
	void evaluate_one(std::pair<G, E>&& r) const
	{
		auto& [grid, equation] = r;
		expr::result(equation, expr::BaseData<G>::get(grid));
	}

};



template<Axis ax, size_t O, typename Sp, typename G>
auto apply_directional_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.template generalized_directional_derivative<ax, O>(std::forward<G>(e), n);
}

template<Axis ax, size_t O, typename Sp, typename G, typename std::enable_if_t<(O == 1), int> = 0>
auto apply_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.template gradient<ax>(std::forward<G>(e), n);
}

template<Axis ax, size_t O, typename Sp, typename G, typename std::enable_if_t<(O == 2), int> = 0>
auto apply_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.laplacian(std::forward<G>(e), n);
}

template<Axis ax, size_t O, typename Sp, typename G, typename std::enable_if_t<(O == 3), int> = 0>
auto apply_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.template gradlaplacian<ax>(std::forward<G>(e), n);
}

template<Axis ax, size_t O, typename Sp, typename G, typename std::enable_if_t<(O == 4), int> = 0>
auto apply_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.bilaplacian(std::forward<G>(e), n);
}

template<Axis ax, size_t O, typename Sp, typename G, typename std::enable_if_t<(O > 4), int> = 0>
auto apply_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.template generalized_derivative<ax, O>(std::forward<G>(e), n);
}

template<size_t... Os, typename Sp, typename G>
auto apply_mixed_derivative(Sp const& sp, G&& e, iter_type n)
{
	return sp.template generalized_mixed_derivative<Os...>(std::forward<G>(e), n);
}


template<typename Sp>
template<Axis ax, size_t O>
struct Solver<Sp>::derivative
{
	template<typename G>
	auto operator()(Sp const& sp, G&& e, iter_type n) const
	{
		return apply_derivative<ax, O>(sp, std::forward<G>(e), n);
	}
	static const size_t order = O;
	static const Axis axis = ax;
	static const bool is_directional = (O == 1);
	static const bool is_mixed = false;
};


template<typename Sp>
template<Axis ax, size_t O>
struct Solver<Sp>::directional_derivative
{
	template<typename G>
	auto operator()(Sp const& sp, G&& e, iter_type n) const
	{
		return apply_directional_derivative<ax, O>(sp, std::forward<G>(e), n);
	}
	static const size_t order = O;
	static const Axis axis = ax;
	static const bool is_directional = true;
	static const bool is_mixed = false;
};

namespace symphas::internal
{
	template<typename A, typename L>
	struct max_order_axis;

	template<Axis ax1, size_t O1>
	struct max_order_axis<symphas::lib::axis_list<ax1>, std::index_sequence<O1>>
	{
		static const Axis value = ax1;
	};

	template<Axis ax1, Axis ax2, size_t O1, size_t O2>
	struct max_order_axis<symphas::lib::axis_list<ax1, ax2>, std::index_sequence<O1, O2>>
	{
		static const Axis value = (O1 >= O2) ? ax1 : ax2;
	};

	template<Axis ax1, Axis ax2, Axis ax3, size_t O1, size_t O2, size_t O3>
	struct max_order_axis<symphas::lib::axis_list<ax1, ax2, ax3>, std::index_sequence<O1, O2, O3>>
	{
		static const Axis value = (O1 >= O2 && O1 >= O3) ? ax1 : (O2 >= O1 && O2 >= O3) ? ax2 : ax3;
	};
}

template<typename Sp>
template<size_t... Os>
struct Solver<Sp>::mixed_derivative
{
	template<typename G>
	auto operator()(Sp const& sp, G&& e, iter_type n) const
	{
		return apply_mixed_derivative<Os...>(sp, std::forward<G>(e), n);
	}
	static const size_t order = symphas::math::sum<Os...>();
	static const Axis axis = symphas::internal::max_order_axis<
		decltype(symphas::lib::make_axis_list<sizeof...(Os)>()), std::index_sequence<Os...>>::value;
	static const bool is_directional = true;
	static const bool is_mixed = true;
};

#include "provisionalsystemgroup.h"



template<typename T, typename Ty>
struct symphas::internal::solver_supported_type_match
{
	static const bool value = false;
};

template<typename T, typename S>
struct symphas::internal::solver_supported_system
{

protected:

	template<typename Ty>
	struct supported_wrap
	{
		static const bool value = solver_supported_type<T, void>::value || solver_supported_type<T, Ty>::value;
	};

	template<size_t N, typename Ty>
	static constexpr auto is_supported(MultiBlock<N, Ty>)
	{
		return supported_wrap<any_vector_t<Ty, N>>{};
	}

	template<typename Ty>
	static constexpr auto is_supported(Block<Ty>)
	{
		return supported_wrap<Ty>{};
	}

	static constexpr auto get_value(S s)
	{
		return is_supported(s);
	}

	using wrap_type = typename std::invoke_result_t<decltype(&solver_supported_system<T, S>::get_value), S>;

public:

	static const size_t value = wrap_type::value;
};



template<typename T, typename... S>
struct symphas::internal::solver_supported_systems
{
	static const size_t value = ((solver_supported_system<T, S>::value && ...));
};


template<typename T>
struct symphas::internal::solver_system_type_match
{
	template<typename Ty, size_t D>
	using type = SolverSystem<Ty, D>;
};

template<typename T>
struct symphas::internal::provisional_system_type_match
{
	template<typename Ty, size_t D>
	using type = ProvisionalSystem<Ty, D>;
};

