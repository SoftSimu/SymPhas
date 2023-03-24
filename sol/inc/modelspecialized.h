
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
#include "symboliceval.h"

namespace expr
{
	namespace
	{
		template<Axis ax, typename G, typename Sp>
		using grad_term_t = OpDerivative<typename Solver<Sp>::template derivative<ax, 1>, OpIdentity, OpTerm<OpIdentity, G>, Sp>;

		template<typename E, typename Sp, typename G, Axis... axs, size_t... Is>
		auto euler_lagrange_deriv(E const& e, Sp const& solver, 
			symphas::lib::types_list<G, symphas::lib::axis_list<axs...>>,
			std::index_sequence<Is...>)
		{
			return expr::apply_operators(
				expr::make_operator_derivative<1>(solver) * (
					expr::apply_operators(
						expr::make_column_vector<Is, sizeof...(Is)>() * expr::make_derivative<1, grad_term_t<axs, G, Sp>>(e)) + ...));
		}

		template<typename G, typename E, typename Sp>
		auto euler_lagrange_deriv(E const& e, Sp const& solver)
		{
			constexpr size_t dimension = expr::grid_dim<G>::value;

			return euler_lagrange_deriv(e, solver, 
				symphas::lib::make_axis_list<dimension, G>(), std::make_index_sequence<dimension>{});
		}
	}

	//! -c1/2 * op^2 + c2/4 * op^4 + |grad(op)|^2
	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto landau_fe(OpTerm<OpIdentity, G> const& term, Sp const& solver, coeff_t1 const& c1 = coeff_t1{}, coeff_t2 const& c2 = coeff_t2{})
	{
		return -c1 * expr::make_fraction<1, 2>() * expr::pow<2>(term) + c2 * expr::make_fraction<1, 4>() * expr::pow<4>(term)
			+ expr::make_fraction<1, 2>() * expr::pow<2>(expr::make_operator_derivative<1>(solver)(term));
	}

	//! c2/4 * (op^2 - c1)^2 + |grad(op)|^2
	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto doublewell_fe(OpTerm<OpIdentity, G> const& term, Sp const& solver, coeff_t1 const& c1 = coeff_t1{}, coeff_t2 const& c2 = coeff_t2{})
	{
		return c2 * expr::make_fraction<1, 4>() * expr::pow<2>(expr::pow<2>(term) - c1)
			+ expr::make_fraction<1, 2>() * expr::pow<2>(expr::make_operator_derivative<1>(solver)(term));
	}

	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto cellular_fe(OpTerm<OpIdentity, G> const& term, Sp const& solver, coeff_t1 const& c1 = coeff_t1{})
	{
		return expr::make_fraction<1, 4>() * expr::make_integer<30>() / (c1 * c1) * term * term * expr::pow<2>(expr::symbols::one - term)
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
	enum class DynamicType { NONE, CONSERVED, NONCONSERVED, HYDRODYNAMIC };

	template<DynamicType dynamic0, DynamicType... dynamics>
	struct dynamics_list {};

	template<DynamicType dynamic0>
	struct dynamics_list<dynamic0> 
	{
		static const DynamicType value = dynamic0;
	};

	template<DynamicType dynamic0, DynamicType dynamic1, DynamicType... dynamics>
	struct dynamics_list<dynamic0, dynamic1, dynamics...> {};

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
			((..., expr::printe(std::get<Is>(eqns).second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_model>::value>(_s, std::get<Is>(eqns)...);
		}

		template<typename... A>
		auto make_equations(A const& ...a) const
		{
			((..., expr::printe(a.second, "given equation")));
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
			((..., expr::printe(std::get<Is>(eqns).second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_trait>::value>(forward_systems(), std::get<Is>(eqns)...);
		}

		template<typename... A>
		auto make_equations(A const& ...a) const
		{
			((..., expr::printe(a.second, "given equation")));
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


namespace symphas::internal
{
	template<size_t Dm, typename Sp, typename... Ts>
	struct expand_types_to_model
	{
		using type = typename ModelApplied<Dm, Sp>::template OpTypes<Ts...>;
	};

	template<size_t Dm, typename Sp, typename... Ts>
	struct expand_types_to_model<Dm, Sp, symphas::lib::types_list<Ts...>>
	{
		using type = typename ModelApplied<Dm, Sp>::template OpTypes<Ts...>;
	};
}


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


namespace symphas::internal
{
	using namespace expr::symbols;

	template<DynamicType dynamic>
	using dynamics_key_t = dynamics_list<dynamic>;

	template<DynamicType dynamic>
	using all_dynamics_key_t = symphas::lib::types_list<dynamics_key_t<dynamic>>;


	template<symphas::internal::DynamicType dynamics>
	struct apply_dynamics;

	template<>
	struct apply_dynamics<symphas::internal::DynamicType::CONSERVED>
	{
		template<typename U_D, typename F, typename Sp>
		auto operator()(U_D const& dop, F const& dfe, Sp const& solver)
		{
			return (dop = expr::apply_operators(expr::make_operator_derivative<2>(solver) * dfe));
		}
	};

	template<>
	struct apply_dynamics<symphas::internal::DynamicType::NONCONSERVED>
	{
		template<typename U_D, typename F, typename Sp>
		auto operator()(U_D const& dop, F const& dfe, Sp const& solver)
		{
			return (dop = -expr::apply_operators(dfe));
		}
	};

	template<>
	struct apply_dynamics<symphas::internal::DynamicType::NONE>
	{
		template<typename U_D, typename F, typename Sp>
		auto operator()(U_D const& dop, F const& dfe, Sp const& solver)
		{
			return (dop = dfe);
		}
	};

	// Specifies the dynamics for all fields.
	template<size_t N, typename I, typename E>
	struct special_dynamics
	{
		special_dynamics(OpExpression<E> const& e) : e{ *static_cast<E const*>(&e) } {}
		special_dynamics(I) {}

		E e;
	};

	// Specifies the dynamics for a single field.
	template<size_t N, symphas::internal::DynamicType dynamic>
	struct special_dynamics<N, void, dynamics_key_t<dynamic>> {};

	// Specifies the dynamics for a single field.
	template<size_t N, typename E>
	struct special_dynamics<N, void, E>
	{
		special_dynamics(OpExpression<E> const& e) : e{ *static_cast<E const*>(&e) } {}
		E e;
	};

	// Specifies the dynamics for all fields.
	template<size_t N, typename I>
	struct special_dynamics<N, I, void>
	{
		special_dynamics(I) {}
		special_dynamics(std::index_sequence<N>, I) {}

		template<typename E0>
		auto operator()(OpExpression<E0> const& e)
		{
			return special_dynamics<N, I, E0>(*static_cast<E0 const*>(&e));
		}

		template<symphas::internal::DynamicType dynamic>
		auto operator()(dynamics_key_t<dynamic>)
		{
			return special_dynamics<N, void, dynamics_key_t<dynamic>>();
		}
	};

	// Specifies the dynamics for a single field.
	template<size_t N>
	struct special_dynamics<N, void, void>
	{
		special_dynamics(std::index_sequence<N>) {}

		template<typename E0>
		auto operator()(OpExpression<E0> const& e)
		{
			return special_dynamics<N, void, E0>(*static_cast<E0 const*>(&e));
		}

		template<symphas::internal::DynamicType dynamic>
		auto operator()(dynamics_key_t<dynamic>)
		{
			return special_dynamics<N, void, dynamics_key_t<dynamic>>();
		}
	};


	template<typename I>
	special_dynamics(I) -> special_dynamics<0, I, void>;
	template<size_t N, typename I>
	special_dynamics(std::index_sequence<N>, I) -> special_dynamics<N, I, void>;
	template<size_t N>
	special_dynamics(std::index_sequence<N>) -> special_dynamics<N, void, void>;





	template<int N>
	using func_deriv_i = i_<N, 0>;

	template<typename I>
	constexpr auto func_deriv_sub_var(I)
	{
		return make_sum_variable(I{}, func_deriv_i<-1>{});
	}

	template<int N>
	constexpr auto func_deriv_sub_var()
	{
		return make_sum_variable(func_deriv_i<N>{});
	}


	template<int N, typename E0, typename... Es>
	auto substitute_for_dynamics(
		symphas::lib::types_list<>,
		OpExpression<E0> const& dyn0,
		std::tuple<Es...> const& dfs)
	{
		return *static_cast<E0 const*>(&dyn0);
	}

	template<int N, int N0, int P0, int... Ns, int... Ps, typename... Vs, typename E0, typename... Es>
	auto substitute_for_dynamics(
		symphas::lib::types_list<
			expr::symbols::v_id_type<expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>,
			Vs...>,
		OpExpression<E0> const& dyn0,
		std::tuple<Es...> const& dfs)
	{
		using vvt = expr::symbols::v_id_type<expr::symbols::i_<N0, P0>, expr::symbols::i_<Ns, Ps>...>;
		auto dyn = expr::transform::swap_grid<vvt>(*static_cast<E0 const*>(&dyn0), std::get<size_t(N + N0)>(dfs));
		return substitute_for_dynamics<N>(symphas::lib::types_list<Vs...>{}, dyn, dfs);
	}


	template<typename S>
	struct apply_special_dynamics;

	template<size_t N, DynamicType dynamic>
	struct apply_special_dynamics<special_dynamics<N, void, dynamics_key_t<dynamic>>>
		: apply_dynamics<dynamic>
	{
		apply_special_dynamics(special_dynamics<N, void, dynamics_key_t<dynamic>>)
			: apply_dynamics<dynamic>() {}

		template<typename... U_Ds, typename... Us, typename... Fs, typename Sp>
		auto operator()(std::tuple<U_Ds...> const& dops, std::tuple<Us...> const& ops, std::tuple<Fs...> const& dfes, Sp const& solver)
		{
			return apply_dynamics<dynamic>::operator()(std::get<N>(dops), std::get<N>(dfes), solver);
		}
	};
	
	template<size_t N, typename I0, typename E>
	struct apply_special_dynamics<special_dynamics<N, I0, E>>
	{
		apply_special_dynamics(special_dynamics<N, I0, E> const& e) : e{ e.e } {}

		template<typename... U_Ds, typename... Us, typename... Fs, typename Sp>
		auto operator()(std::tuple<U_Ds...> const& dops, std::tuple<Us...> const& ops, std::tuple<Fs...> const& dfes, Sp const& solver)
		{
			using df_v_types = typename select_v_i_<func_deriv_i<-1>, expr::op_types_t<E>>::type;

			auto with_fes = substitute_for_dynamics<N>(df_v_types{}, e, dfes);
			auto ev = substitute_ops(ops, with_fes);
			return (std::get<N>(dops) = expr::apply_operators(ev));
		}
		
		template<typename... Us, typename E0>
		auto substitute_ops(std::tuple<Us...> const& ops, OpExpression<E0> const& with_fes)
		{
			using v_types = typename select_v_<expr::op_types_t<E0>>::type;
			return substitute_for_dynamics<N>(v_types{}, *static_cast<E0 const*>(&with_fes), ops);
		}

		E e;
	};

	template<size_t N, typename E>
	struct apply_special_dynamics<special_dynamics<N, void, E>>
	{
		apply_special_dynamics(special_dynamics<N, void, E> const& e) : e{ e.e } {}

		template<typename... U_Ds, typename... Us, typename... Fs, typename Sp>
		auto operator()(std::tuple<U_Ds...> const& dops, std::tuple<Us...> const& ops, std::tuple<Fs...> const& dfes, Sp const& solver)
		{
			using v_types = typename select_v_<expr::op_types_t<E>>::type;
			auto ev = substitute_for_dynamics<0>(v_types{}, e, dfes);
			return (std::get<N>(dops) = expr::apply_operators(ev));
		}

		E e;
	};

	template<typename S>
	apply_special_dynamics(S) -> apply_special_dynamics<S>;

	template<DynamicType dynamic, size_t I>
	constexpr DynamicType dynamic_i = dynamic;

	template<typename Dyn0, typename Dyn1, typename... Us>
	auto replace_dynamics_if_none(Dyn0 const& dynamic, Dyn1 const&, std::tuple<Us...> const&)
	{
		return dynamic;
	}

	template<size_t N, size_t N0, typename I, typename E, typename... Us>
	auto replace_dynamics_if_none(
		special_dynamics<N, void, dynamics_key_t<DynamicType::NONE>>,
		special_dynamics<N0, I, E> const& dynamic, 
		std::tuple<Us...> const& ops)
	{
		return special_dynamics(std::index_sequence<N>{}, I{})(dynamic.e);
	}

	template<size_t N, size_t N0, typename I, DynamicType dynamic0, typename... Us>
	auto replace_dynamics_if_none(
		special_dynamics<N, void, dynamics_key_t<DynamicType::NONE>>,
		special_dynamics<N0, I, dynamics_key_t<dynamic0>> const& dynamic, 
		std::tuple<Us...> const& ops)
	{
		return special_dynamics(std::index_sequence<N>{}, I{})(dynamics_key_t<dynamic0>{});
	}

	template<size_t... Is, typename... Dyns, size_t N, typename E, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<Dyns...> const& dynamics,
		special_dynamics<N, void, E> const& dynamic0, std::tuple<Us...> const& ops)
	{
		return std::tuple_cat(
			symphas::lib::get_tuple_lt<N>(dynamics),
			std::make_tuple(special_dynamics(std::index_sequence<N>{})(dynamic0.e)),
			symphas::lib::get_tuple_ge<N + 1>(dynamics));
	}

	template<size_t... Is, typename... Dyns, size_t N, typename I, typename E, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<Dyns...> const& dynamics,
		special_dynamics<N, I, E> const& dynamic0, std::tuple<Us...> const& ops)
	{
		return std::make_tuple(replace_dynamics_if_none(std::get<Is>(dynamics), dynamic0, ops)...);
	}


	template<size_t... Is, typename Dyn0, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<Dyn0> const& dynamic, std::tuple<Us...> const& ops)
	{
		auto dl = std::make_tuple(special_dynamics(std::index_sequence<Is>{})(dynamics_key_t<DynamicType::NONE>{})...);
		return interpret_dynamics(std::index_sequence<Is...>{}, dl, std::get<0>(dynamic), ops);
	}


	template<size_t... Is, typename... DynDs, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<DynDs...> const& all_dynamics,
		std::tuple<> const& dynamics, std::tuple<Us...> const& ops)
	{
		return all_dynamics;
	}

	template<size_t... Is, typename... DynDs, typename Dyn0, typename... Dyns, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<DynDs...> const& all_dynamics,
		std::tuple<Dyn0, Dyns...> const& dynamics, std::tuple<Us...> const& ops)
	{
		auto next = interpret_dynamics(std::index_sequence<Is...>{}, all_dynamics, std::get<0>(dynamics), ops);
		return interpret_dynamics(std::index_sequence<Is...>{}, next, symphas::lib::get_tuple_ge<1>(dynamics), ops);
	}

	template<size_t... Is, typename Dyn0, typename Dyn1, typename... Dyns, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<Dyn0, Dyn1, Dyns...> const& dynamics, std::tuple<Us...> const& ops)
	{
		auto first = interpret_dynamics(std::index_sequence<Is...>{}, std::make_tuple(std::get<0>(dynamics)), ops);
		return interpret_dynamics(std::index_sequence<Is...>{}, first, symphas::lib::get_tuple_ge<1>(dynamics), ops);
	}

	template<size_t... Is, DynamicType dynamic, typename... Us>
	auto interpret_dynamics(std::index_sequence<Is...>, std::tuple<all_dynamics_key_t<dynamic>>, std::tuple<Us...> const& ops)
	{
		return dynamics_list<dynamic_i<dynamic, Is>...>{};
	}

}

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
	static const size_t Dm = model_dimension<parent_trait>::value;

	static const size_t Sn = model_num_parameters<parent_trait>::value;
	using seq_t = std::make_index_sequence<Sn>;

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

	auto param(size_t I)
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

	auto param(size_t N, size_t I)
	{
		return param(N * Sn + I);
	}

    template<NoiseType nt, typename T, typename... Ts>
    auto make_noise(Ts&& ...args) const
    {
        return expr::make_noise<nt, T, Dm>(
            parent_trait::template system<0>().dims,
            parent_trait::template system<0>().get_info().get_widths().get(),
            &solver.dt, std::forward<Ts>(args)...);
    }

    template<NoiseType nt, size_t Z, typename G, typename... Ts>
    auto make_noise(OpTerm<OpIdentity, Variable<Z, G>>, Ts&& ...args) const
    {
        using T = model_field_t<parent_trait, Z>;
        return expr::make_noise<nt, T, Dm>(
            parent_trait::template system<0>().dims,
            parent_trait::template system<0>().get_info().get_widths().get(),
            &solver.dt, std::forward<Ts>(args)...);
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

	template<symphas::internal::DynamicType dynamic0, symphas::internal::DynamicType... dynamics, typename E>
	auto generate_equations(std::tuple<symphas::internal::dynamics_key_t<dynamic0>, symphas::internal::dynamics_key_t<dynamics>...>, OpExpression<E> const& e)
	{
		using dynamics_list = symphas::internal::dynamics_list<dynamic0, dynamics...>;
		return generate_equations_apply(*static_cast<E const*>(&e), seq_t{}, dynamics_list{});
	}

	template<typename... special_dynamics, typename E, size_t... Ns>
	auto generate_equations(std::tuple<special_dynamics...> const& dynamics, OpExpression<E> const& e, std::index_sequence<Ns...>)
	{
		return generate_equations_apply(
			*static_cast<E const*>(&e), 
			symphas::internal::interpret_dynamics(seq_t{}, dynamics, std::make_tuple(op<Ns>()...)));
	}

	template<typename... special_dynamics, typename E>
	auto generate_equations(std::tuple<special_dynamics...> const& dynamics, OpExpression<E> const& e)
	{
		return generate_equations(dynamics, *static_cast<E const*>(&e), seq_t{});
	}

protected:

	template<typename E, symphas::internal::DynamicType... dynamics>
	auto generate_equations_apply(
		OpExpression<E> const& e,
		symphas::internal::dynamics_list<dynamics...> const& ed)
	{
		return generate_equations(std::make_tuple(symphas::internal::dynamics_key_t<dynamics>{}...), *static_cast<E const*>(&e));
	}

	template<typename E, size_t... Ns, typename... Is, typename... Es>
	auto generate_equations_apply(
		OpExpression<E> const& e,
		std::tuple<symphas::internal::special_dynamics<Ns, Is, Es>...> const& ed)
	{
		auto dfes = generate_equations(
			std::make_tuple(symphas::internal::all_dynamics_key_t<symphas::internal::DynamicType::NONE>{}),
			*static_cast<E const*>(&e));

		return std::make_tuple(symphas::internal::apply_special_dynamics(std::get<Ns>(ed))(
			std::make_tuple(dop<Ns>()...), 
			std::make_tuple(op<Ns>()...), 
			std::make_tuple(std::get<1>(std::get<Ns>(dfes))...),
			solver)...);
	}

	template<typename E, size_t... Ns, symphas::internal::DynamicType... dynamics>
	auto generate_equations_apply(
		OpExpression<E> const& fe, 
		std::index_sequence<Ns...>, 
		symphas::internal::dynamics_list<dynamics...>)
	{
		return std::make_tuple(
			symphas::internal::apply_dynamics<dynamics>{}(
				dop<Ns>(),
				expr::euler_lagrange_apply(op<Ns>(), *static_cast<E const*>(&fe), solver),
				solver)...);
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
	static const size_t Dm = model_dimension<parent_type>::value; \
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

#define DEFINE_MODEL_FIELD_NAMES_FORMAT(MODEL_NAME, FORMAT) \
template<> \
struct model_field_name<model_ ## MODEL_NAME::TraitEquationModel> \
{ \
	const char* operator()(int index) { \
		static char** names = build_names(); \
		return names[index]; \
	} \
	char** build_names() \
	{ \
		constexpr size_t len = 100; \
		char** names = new char*[len]; \
		char buffer[100]; \
		for (size_t i = 1; i <= len; ++i) { \
			sprintf(buffer, FORMAT, int(i)); \
			names[i - 1] = new char[std::strlen(buffer) + 1]; \
			std::strcpy(names[i - 1], buffer); \
		} \
		return names; \
	} \
};

//! @}



namespace symphas::internal
{
	template<typename S, typename... Xs>
	struct series_index_selection;

	template<typename, size_t N>
	constexpr size_t ivalue = N;

	template<typename... Ss, typename... Xs>
	struct series_index_selection<std::tuple<Ss...>, Xs...>
	{
		series_index_selection(std::tuple<Ss...> const& systems, Xs... xs) : systems{ &const_cast<std::tuple<Ss...>&>(systems) }, limits{ xs... } {}

		template<typename E, size_t... Ls, size_t... Ns>
		auto substitute_systems(OpExpression<E> const& e, std::index_sequence<Ls...>, std::index_sequence<Ns...>)
		{
			return sum(*static_cast<E const*>(&e))
				.template select<Ls...>(Xs{}...)
				(std::make_tuple(Term(expr::as_variable<Ns>(std::get<Ns>(*systems).as_grid()))...));
		}

		template<typename E>
		auto operator()(OpExpression<E> const& e)
		{
			return substitute_systems(
				*static_cast<E const*>(&e), 
				symphas::lib::seq_repeating_value_t<sizeof...(Xs), size_t, sizeof...(Ss)>{},
				std::make_index_sequence<sizeof...(Ss)>{});
		}


		std::tuple<Ss...> *systems;
		std::tuple<Xs...> limits;
	};

	template<typename... Ss, typename... Xs>
	series_index_selection(std::tuple<Ss...>, Xs...) -> series_index_selection<std::tuple<Ss...>, Xs...>;
}


struct ExpressionStats
{
protected:

	template<typename T, size_t D>
	static auto max(any_vector_t<T, D> const& v1, any_vector_t<T, D> const& v2)
	{
		using std::max;
		any_vector_t<T, D> vm;
		for (iter_type i = 0; i < D; ++i)
		{
			vm[i] = max(v1[i], v2[i]);
		}
		return vm;
	}

	template<typename T, size_t D>
	static auto min(any_vector_t<T, D> const& v1, any_vector_t<T, D> const& v2)
	{
		using std::min;
		any_vector_t<T, D> vm;
		for (iter_type i = 0; i < D; ++i)
		{
			vm[i] = min(v1[i], v2[i]);
		}
		return vm;
	}

public:

	template<typename E>
	auto mean(OpExpression<E> const& e) const
	{
		expr::eval_type_t<E> sum;
		expr::result(OpVoid{}, sum);

		auto sumop = expr::make_term(sum);
		
		len_type len = expr::data_length(*static_cast<E const*>(&e));
		for (iter_type i = 0; i < len; ++i)
		{
			sum += (*static_cast<E const*>(&e)).eval(i);
		}

		sum *= (1. / len);
		return expr::make_literal(sum);
	}

	template<typename E>
	auto max(OpExpression<E> const& e) const
	{
		using std::max;

		expr::eval_type_t<E> result;
		expr::result((*static_cast<E const*>(&e)).eval(0), result);

		len_type len = expr::data_length(*static_cast<E const*>(&e));
		for (iter_type i = 1; i < len; ++i)
		{
			result = max((*static_cast<E const*>(&e)).eval(i), result);
		}

		return expr::make_literal(result);
	}

	template<typename E>
	auto min(OpExpression<E> const& e) const
	{
		using std::max;

		expr::eval_type_t<E> result;
		expr::result((*static_cast<E const*>(&e)).eval(0), result);

		len_type len = expr::data_length(*static_cast<E const*>(&e));
		for (iter_type i = 1; i < len; ++i)
		{
			result = min((*static_cast<E const*>(&e)).eval(i), result);
		}

		return expr::make_literal(result);
	}
};
