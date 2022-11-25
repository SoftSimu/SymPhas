
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
 * PURPOSE: Defines the base type and construction of a specialized model.
 *
 * ***************************************************************************
 */

#pragma once

#include "model.h"
#include "expressionaggregates.h"
#include "provisionalsystemgroup.h"
#include "expressionsubstitutables.h"

namespace expr
{
	namespace
	{
		template<Axis ax, typename G, typename Sp>
		using grad_term_t = OpFuncDerivative<typename Solver<Sp>::template derivative<ax, 1>, OpIdentity, OpTerm<OpIdentity, G>, Sp>;

		template<typename E, typename Sp, typename G, Axis... axs, size_t... Is>
		auto euler_lagrange_deriv(OpExpression<E> const& e, Sp const& solver, 
			symphas::lib::types_list<G, symphas::lib::axis_list<axs...>>,
			std::index_sequence<Is...>)
		{
			return expr::apply_operators(
				expr::make_operator_derivative<1>(solver) * (
					expr::apply_operators(
						expr::make_column_vector<Is, sizeof...(Is)>() * expr::make_derivative<1, grad_term_t<axs, G, Sp>>(*static_cast<E const*>(&e))) + ...));
		}

		template<typename G, typename E, typename Sp>
		auto euler_lagrange_deriv(OpExpression<E> const& e, Sp const& solver)
		{
			constexpr size_t dimension = expr::grid_dim<G>::value;

			return euler_lagrange_deriv(*static_cast<E const*>(&e), solver, 
				symphas::lib::make_axis_list<dimension, G>(), std::make_index_sequence<dimension>{});
		}
	}

	template<typename G, typename Sp>
	auto landau_fe(OpTerm<OpIdentity, G> const& term, Sp const& solver)
	{
		return -expr::make_fraction<1, 2>() * expr::pow<2>(term) + expr::make_fraction<1, 4>() * expr::pow<4>(term)
			+ expr::make_fraction<1, 2>() * expr::pow<2>(expr::make_operator_derivative<1>(solver)(term));
	}

	template<typename G, typename Sp>
	auto doublewell_fe(OpTerm<OpIdentity, G> const& term, Sp const& solver)
	{
		return expr::make_fraction<1, 4>() * expr::pow<2>(expr::pow<2>(term) - expr::symbols::one)
			+ expr::make_fraction<1, 2>() * expr::pow<2>(expr::make_operator_derivative<1>(solver)(term));
	}

	template<typename G, typename E, typename Sp>
	auto euler_lagrange_apply(OpTerm<OpIdentity, G> const& term, OpExpression<E> const& e, Sp const& solver)
	{
		auto interface_term = euler_lagrange_deriv<G>(*static_cast<E const*>(&e), solver);
		auto bulk_term = expr::apply_operators(expr::make_derivative<1, G>(*static_cast<E const*>(&e)));
		return bulk_term - interface_term;
	}
}

namespace symphas::internal
{

	/* whether the model exhibits conserved or nonconserved EVOLUTION
	 */
	enum class DynamicType { CONSERVED, NONCONSERVED };

	template<DynamicType... dynamics>
	struct dynamics_list {};

	template<typename parent_model>
	struct MakeEquation : parent_model
	{
		using parent_model::parent_model;
		using parent_model::_s;
		using parent_model::solver;

		template<typename... As>
		auto make_equations(std::tuple<As...> const& eqns) const
		{
			return make_equations(eqns, std::make_index_sequence<sizeof...(As)>{});
		}

		template<typename... As, size_t... Is>
		auto make_equations(std::tuple<As...> const& eqns, std::index_sequence<Is...>) const
		{
			((..., expr::printf(std::get<Is>(eqns).second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_model>::value>(_s, std::get<Is>(eqns)...);
		}

		template<typename... A>
		auto make_equations(A const& ...a) const
		{
			((..., expr::printf(a.second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_model>::value>(_s, a...);
		}

	};

	template<typename parent_trait>
	struct MakeEquationProvisional : parent_trait
	{
		using parent_trait::parent_trait;
		using parent_trait::_s;
		using parent_trait::temp;
		using parent_trait::solver;

		template<typename... As>
		auto make_equations(std::tuple<As...> const& eqns) const
		{
			return make_equations(eqns, std::make_index_sequence<sizeof...(As)>{});
		}

		template<typename... As, size_t... Is>
		auto make_equations(std::tuple<As...> const& eqns, std::index_sequence<Is...>) const
		{
			((..., expr::printf(std::get<Is>(eqns).second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_trait>::value>(forward_systems(), std::get<Is>(eqns)...);
		}

		template<typename... A>
		auto make_equations(A const& ...a) const
		{
			((..., expr::printf(a.second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_trait>::value>(forward_systems(), a...);
		}

	protected:

		template<typename... Ts1, typename... Ts2, size_t... Is1, size_t... Is2>
		decltype(auto) forward_systems(
			std::tuple<Ts1...> const& a, std::tuple<Ts2...> const& b,
			std::index_sequence<Is1...>, std::index_sequence<Is2...>) const
		{
			return std::forward_as_tuple(std::get<Is1>(a)..., std::get<Is2>(b)...);
		}

		template<typename... Ts1, typename... Ts2>
		decltype(auto) forward_systems(std::tuple<Ts1...> const& a, std::tuple<Ts2...> const& b) const
		{
			return forward_systems(a, b, std::make_index_sequence<sizeof...(Ts1)>{}, std::make_index_sequence<sizeof...(Ts2)>{});
		}

		decltype(auto) forward_systems() const
		{
			return forward_systems(_s, temp._s);
		}


	};

}


//! Enclosing class for the generalized phase field problem representation.
/*!
 * An applied model uses the base ::Model class and defines equations of motion
 * in order to implement a phase field problem.
 * 
 * \tparam D The dimension of the phase field problem.
 * \tparam Sp The solver type used to solve the problem.
 */
template<size_t D, typename Sp>
struct ModelApplied
{
	//! Encapsulates order parameter types of the phase field problem.
	/*!
	 * Encapsulates order parameter types of the phase field problem.
	 * 
	 * \tparam S... The order parameter types.
	 */
	template<typename... S>
	struct OpTypes
	{
		//! Specialized model representing the generalized phase field problem.
		/*!
		 * The equation of motion is created on object
		 * initialization, and the model implements update and equation functions.
		 */
		template<template<typename> typename Eq,
			typename eq_type = Eq<symphas::internal::MakeEquation<Model<D, Sp, S...>>>>
		struct Specialized;


		//! Encapsulates specialized model with provisional equations.
		/*!
		 * Encapsulates specialized model with provisional equations.
		 * 
		 * \tparam P... The provisional variable types.
		 */
		template<typename... P>
		struct ProvTypes
		{

			//! Specialized model that uses provisional equations.
			/*!
			 * Both equation and provisional expression objects are created on object
			 * initialization, and the model implements update and equation functions.
			 */
			template<template<typename> typename Eq, template<typename, typename...> typename Pr,
				typename pr_type = Pr<Model<D, Sp, S...>, P...>, 
				typename eq_type = Eq<symphas::internal::MakeEquationProvisional<pr_type>>>
			struct Specialized;
		};
	};
};

// ****************************************************************************************

template<size_t D, typename Sp>
template<typename... S>
template<template<typename> typename Eq, typename eq_type>
struct ModelApplied<D, Sp>::OpTypes<S...>::Specialized : eq_type
{
	using M = Model<D, Sp, S...>;
	using M::solver;
	
	template<typename Sp0>
	using eq_type_solver = Eq<symphas::internal::MakeEquation<Model<D, Sp0, S...>>>;

	template<template<template<typename> typename, typename> typename SpecializedModel, typename Sp0>
	using impl_type = SpecializedModel<Eq, eq_type_solver<Sp0>>;
	
	using parent_type = eq_type;
	using parent_type::parent_type;

	using eqs = typename std::invoke_result_t<decltype(&parent_type::make_equations), parent_type>;
	eqs equations;

	Specialized(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_type(coeff, num_coeff, parameters), equations{ parent_type::make_equations() } {}
	Specialized(symphas::problem_parameters_type const& parameters) : Specialized(nullptr, 0, parameters) {}

	impl_type<ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized, Sp>& operator=(
		impl_type<ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized, Sp> other)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(this), *static_cast<parent_type*>(&other));
		swap(equations, other.equations);
	}



	void update(double time)
	{
		M::update_systems(time);
	}

	void equation()
	{
		M::solver.equations(equations);
	}

protected:

	Specialized() : parent_type(), equations{ parent_type::make_equations() } {}

};

template<size_t D, typename Sp>
template<typename... S>
template<typename... P>
template<template<typename> typename Eq, template<typename, typename...> typename Pr, typename pr_type, typename eq_type>
struct ModelApplied<D, Sp>::OpTypes<S...>::ProvTypes<P...>::Specialized : eq_type
{
	using M = Model<D, Sp, S...>;
	using M::solver;



	template<typename Sp0>
	using pr_type_solver = Pr<Model<D, Sp0, S...>, P...>;
	template<typename Sp0>
	using eq_type_solver = Eq<symphas::internal::MakeEquationProvisional<pr_type_solver<Sp0>>>;

	template<template<template<typename> typename, template<typename, typename...> typename, typename, typename> typename SpecializedModel, typename Sp0>
	using impl_type = SpecializedModel<Eq, Pr, pr_type_solver<Sp0>, eq_type_solver<Sp0>>;

	using parent_type = eq_type;
	using parent_type::parent_type;
	using parent_type::temp;

	using eqs = typename std::invoke_result_t<decltype(&parent_type::make_equations), parent_type>;
	using prs = typename std::invoke_result_t<decltype(&parent_type::make_provisionals), parent_type>;


	prs provisionals;
	eqs equations;

	Specialized(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_type(coeff, num_coeff, parameters), provisionals{ parent_type::make_provisionals() }, equations{ parent_type::make_equations() } {}
	Specialized(symphas::problem_parameters_type const& parameters) : Specialized(nullptr, 0, parameters) {}

	impl_type<ModelApplied<D, Sp>::template OpTypes<S...>::template ProvTypes<P...>::template Specialized, Sp>& operator=(
		impl_type<ModelApplied<D, Sp>::template OpTypes<S...>::template ProvTypes<P...>::template Specialized, Sp> other)
	{
		using std::swap;
		
		swap(*static_cast<parent_type*>(this), *static_cast<parent_type*>(&other));
		swap(equations, other.equations);
		swap(provisionals, other.provisionals);
	}



	void update(double time)
	{
		M::update_systems(time);
		M::solver.evaluate(provisionals);
		temp.update_systems(parent_type::lastindex, time);
	}

	void equation()
	{
		M::solver.equations(equations);
	}


	//! Returns the underlying grid of the provisional variable.
	/*!
	 * Returns a reference to the grid for the provisional variable
	 * at the given index.
	 *
	 * \tparam I The index of the provisional variable to return.
	 */
	template<size_t I>
	auto& provisional()
	{
		return temp.template grid<I>();
	}

	//! Returns the underlying grid of the provisional variable.
	/*!
	 * Returns a reference to the grid for the provisional variable 
	 * at the given index.
	 *
	 * \tparam I The index of the provisional variable to return.
	 */
	template<size_t I>
	const auto& provisional() const
	{
		return temp.template grid<I>();
	}

protected:

	Specialized() : parent_type(), provisionals{ parent_type::make_provisionals() }, equations{ parent_type::make_equations() } {}

};

// ****************************************************************************************


template<symphas::internal::DynamicType dynamics>
struct apply_dynamics;

template<>
struct apply_dynamics<symphas::internal::DynamicType::CONSERVED>
{
	template<typename U_D, typename U, typename F, typename Sp>
	auto operator()(U_D const& dop, U const& op, F const& fe, Sp const& solver)
	{
		return (dop = expr::make_operator_derivative<2>(solver) * expr::euler_lagrange_apply(op, fe, solver));
	}
};

template<>
struct apply_dynamics<symphas::internal::DynamicType::NONCONSERVED>
{
	template<typename U_D, typename U, typename F, typename Sp>
	auto operator()(U_D const& dop, U const& op, F const& fe, Sp const& solver)
	{
		return (dop = -expr::euler_lagrange_apply(op, fe, solver));
	}
};

//! Provides usage for phase field variables in the equations of motion.
/*! 
 * This class uses CRTP to complete the specialized model.
 * 
 * Provides functions for referring to the order parameter.
 * 
 * \tparam enclosing_type The specialized TraitEquation object (hence the CRTP
 * pattern).
 * \tparam parent_trait Depending on whether the specialized model uses
 * provisional equations or not, this is either TraitProvisional or 
 * MakeEquation.
 */
template<template<typename> typename enclosing_type, typename parent_trait>
struct TraitEquation : parent_trait
{
	using parent_trait::parent_trait;
	using parent_trait::solver;
	using seq_t = std::make_index_sequence<model_num_parameters<parent_trait>::value>;

	template<size_t I>
	auto op()
	{
		return expr::make_term<I>(
			NamedData(
				parent_trait::template grid<I>(),
				model_field_name<enclosing_type>{}(I)
			));
	}

	template<size_t I>
	auto dop()
	{
		return OpLHS(expr::as_variable<I>(parent_trait::template system<I>()));
	}

	template<size_t I>
	auto param()
	{
		if (I < parent_trait::num_coeff)
		{
			return expr::make_literal(parent_trait::coeff[I]);
		}
		else
		{
			return expr::make_literal(DEFAULT_COEFF_VALUE);
		}
	}

	template<template<typename> typename other_enclosing_type, typename other_parent_trait>
	explicit operator other_enclosing_type<other_parent_trait>()
	{
		return *static_cast<model_base_t<parent_trait>*>(this);
	}

	template<template<typename> typename other_enclosing_type, typename other_parent_trait>
	explicit operator const other_enclosing_type<other_parent_trait>() const
	{
		return *static_cast<model_base_t<parent_trait> const*>(this);
	}

	template<symphas::internal::DynamicType... dynamics, typename E>
	auto generate_equations(OpExpression<E> const& e)
	{
		using dynamics_list = symphas::lib::types_list<symphas::internal::dynamics_list<dynamics>...>;
		return generate_equations(*static_cast<E const*>(&e), seq_t{}, dynamics_list{});
	}

	template<symphas::internal::DynamicType... dynamics, typename... As>
	auto generate_equations(OpAdd<As...> const& sums)
	{
		using dynamics_list = symphas::lib::types_list<symphas::internal::dynamics_list<dynamics>...>;
		return generate_equations(expand_sum(sums, std::make_index_sequence<sizeof...(As)>{}), seq_t{}, dynamics_list{});
	}

	template<symphas::internal::DynamicType... dynamics, typename I0, typename E, typename... Is>
	auto generate_equations(OpCompound<expr::internal::CompoundOp::ADD, I0, E, Is...> const& sum)
	{
		using dynamics_list = symphas::lib::types_list<symphas::internal::dynamics_list<dynamics>...>;
		return generate_equations(expand_sum(sum, seq_t{}), seq_t{}, dynamics_list{});
	}

protected:

	template<typename E, size_t... Ns, symphas::internal::DynamicType... dynamics>
	auto generate_equations(E const& fe, std::index_sequence<Ns...>, symphas::lib::types_list<symphas::internal::dynamics_list<dynamics>...>)
	{
		return std::make_tuple(apply_dynamics<dynamics>{}(dop<Ns>(), op<Ns>(), fe, solver)...);
	}

	template<typename E>
	auto expand_sum(OpExpression<E> const& e)
	{
		return (*static_cast<E const*>(&e));
	}

	template<typename I0, typename E, typename... Is, size_t... Ns>
	auto expand_sum(OpCompound<expr::internal::CompoundOp::ADD, I0, E, Is...> const& sum, std::index_sequence<Ns...>)
	{
		return sum.template expand<0, Is...>(op<Ns>()...);
	}

	template<typename I0, typename E, typename... Is>
	auto expand_sum(OpCompound<expr::internal::CompoundOp::ADD, I0, E, Is...> const& sum)
	{
		return expand_sum(sum, seq_t{});
	}

	template<typename... As, size_t... Is>
	auto expand_sum(OpAdd<As...> const& sums, std::index_sequence<Is...>)
	{
		return (expand_sum(expr::get<Is>(sums)) + ...);
	}


};


//! Manages the provisional variables.
/*!
 * This class uses CRTP to complete the specialized model in the same way as
 * TraitEquation.
 * 
 * A state variable is added for the provisional equations.
 * 
 * Provisional variables come with some exceptions:
 * 
 * - For solvers that implement boundary conditions, the boundary conditions
 * of provisional variables are not well defined. This is because typically,
 * the boundaries of the provisional variables are never used. If they are
 * used, (like in a nonlocal operation or finite difference stencil) then
 * the result is not well defined.
 * - Taking the derivative of provisional variables is numerically not 
 * well defined, as it would evaluate the evaluated values, rather than the
 * expression that defines the provisional variable.
 * - Expressions of provisional variables are executed sequentially, so
 * while using the 2nd provisional variable in the expression of the first
 * would be undefined, it is appropriate to use the first provisional
 * variable in the expression of the second.
 * 
 * \tparam enclosing_type The specialized TraitEquation object (hence the CRTP
 * pattern).
 * \tparam parent_model This is the base Model.
 * \tparam P... The types of the provisional variables.
 */
template<template<typename> typename enclosing_type, typename parent_model, typename... P>
struct TraitProvisional : TraitEquation<enclosing_type, parent_model>
{
	using parent_trait = TraitEquation<enclosing_type, parent_model>;

	//! Creates the provisional equation container.
	TraitProvisional(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_trait(coeff, num_coeff, parameters), temp{ parameters.get_interval_data()[0], parameters.get_boundary_data()[0] } {}

	using solver_type = typename model_solver<parent_model>::type;

	template<typename Ty, size_t D>
	using ProvisionalSystemApplied = typename symphas::provisional_system_type<solver_type>::template type<Ty, D>;

	ProvisionalSystemGroup<ProvisionalSystemApplied, model_dimension<parent_model>::value, P...> temp;

	//! Method for using the provisional variable. 
	/*!
	 * Provisional variables are used to store the result of expressions, where
	 * they serve as an intermediate result before being used in the equations
	 * of motion. 
	 * 
	 * Provisional variables are given a display name that corresponds to their
	 * equation.
	 */
	template<size_t I>
	auto var()
	{
#ifdef PRINTABLE_EQUATIONS
		std::ostringstream ss;
		ss << "var" << I;
		return expr::make_term<model_num_parameters<parent_model>::value + I>(
			NamedData(
				temp.template grid<I>(),
				ss.str()
			));
#else
		return expr::make_term<model_num_parameters<parent_model>::value + I>(temp.template grid<I>());
#endif
	}

protected:

	template<typename... A>
	auto make_provisionals(A&&... a) const
	{
		return std::make_tuple(std::forward<A>(a)...);
	}
};


// ****************************************************************************************

/*!
 * \addtogroup modelmacros
 * @{
 */


//! \cond
#define PROVISIONAL_TRAIT_FORWARD_DECL template<typename parent_trait, typename... P> struct TraitProvisionalModel;
#define EQUATION_TRAIT_FORWARD_DECL template<typename parent_trait> struct TraitEquationModel;
//! \endcond


//! Defines a TraitProvisional child class used to define provisional variables.
/*!
 * The equations for the provisional variables are provided as arguments, and
 * a new class of type TraitProvisionalModel is generated. Since this name is
 * not unique, this macro should be invoked within a namespace.
 * 
 * \param ... The equations of the provisional variables.
 */
#define PROVISIONAL_TRAIT_DEFINITION(...) \
template<typename parent_model, typename... P> \
struct TraitProvisionalModel : TraitProvisional<TraitEquationModel, parent_model, P...> \
{ \
	using parent_type = TraitProvisional<TraitEquationModel, parent_model, P...>; \
	using parent_type::solver; \
	using parent_type::parent_type; \
	auto make_provisionals() { return parent_type::template make_provisionals(__VA_ARGS__); } \
};

//! Defines a TraitEquation child class used to define dynamical equations.
/*!
 * This macro only defines the beginning of the class, but the definition
 * is incomplete without specifying #EQUATION_TRAIT_DEFINITION afterwards.
 * This is only to provide a preamble section that can be
 * referenced when the equations of motion are provided. Since this name is
 * not unique, this macro should be invoked within a namespace.
 *
 * \param ... The preamble which goes before the equations of motion are
 * specified.
 */
#define EQUATION_TRAIT_PREAMBLE(...) \
using namespace expr; \
template<typename parent_trait> \
struct TraitEquationModel : TraitEquation<TraitEquationModel, parent_trait> \
{ \
	using parent_type = TraitEquation<TraitEquationModel, parent_trait>; \
	using parent_type::solver; \
	using parent_type::parent_type; \
	using parent_type::generate_equations; \
	auto make_equations() { \
		using namespace std; using namespace expr; using namespace expr::symbols; \
		constexpr size_t D = model_dimension<parent_type>::value; \
		UNUSED(D) \
		__VA_ARGS__

//! Defines a TraitEquation child class used to define dynamical equations.
/*!
 * The equations of motion for the phase fields are provided as arguments, and
 * a new class of type TraitEquationModel is generated. Since this name is
 * not unique, this macro should be invoked within a namespace.
 *
 * \param ... The equations of the phase fields.
 */
#define EQUATION_TRAIT_DEFINITION(...) \
return parent_type::template make_equations(__VA_ARGS__); } \
};


 //! Select custom names of phase fields for a given model.
 /*!
  * The given model will have the phase fields have names different from
  * the default ones. Each name in the list has to be given in quotes, as
  * a typical C++ string. The entire list must be surrounded in parentheses.
  *
  * If changing the names of a multi-phase field model, all the names
  * have to be provided. Moreover, the order of the names corresponds to the 
  * order of the phase fields specified to the model.
  *
  * \param MODEL_NAME The name of the model that the names apply to.
  * \param NAMES The parentheses surrounded list of names that the phase
  * fields will be called.
  */
#define DEFINE_MODEL_FIELD_NAMES(MODEL_NAME, NAMES) \
template<> \
struct model_field_name<model_ ## MODEL_NAME::TraitEquationModel> \
{ \
	const char* operator()(int index) { \
		static const char* names[]{SINGLE_ARG NAMES}; \
		return names[index]; \
	} \
};

//! @}



