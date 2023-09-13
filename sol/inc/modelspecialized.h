
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
#include "provisionalsystemgroup.h"

namespace expr
{
	//! -c1/2 * op^2 + c2/4 * op^4 + |grad(op)|^2
	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto landau_fe(Sp const& solver, OpTerm<OpIdentity, G> const& term, coeff_t1 const& c1 = coeff_t1{}, coeff_t2 const& c2 = coeff_t2{})
	{
		return -c1 * expr::make_fraction<1, 2>() * expr::pow<2>(term) + c2 * expr::make_fraction<1, 4>() * expr::pow<4>(term)
			+ OpChain(OpOperatorChain(expr::make_fraction<1, 2>(), OpIdentity{}), expr::pow<2>(expr::make_operator_derivative<1>(solver)(term)));
	}

	//! c2/4 * (op^2 - c1)^2 + |grad(op)|^2
	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto doublewell_fe(Sp const& solver, OpTerm<OpIdentity, G> const& term, coeff_t1 const& c1 = coeff_t1{}, coeff_t2 const& c2 = coeff_t2{})
	{
		return c2 * expr::make_fraction<1, 4>() * expr::pow<2>(expr::pow<2>(term) - c1)
			+ OpChain(OpOperatorChain(expr::make_fraction<1, 2>(), OpIdentity{}), expr::pow<2>(expr::make_operator_derivative<1>(solver)(term)));
	}

	template<typename G, typename Sp, typename coeff_t1 = OpIdentity, typename coeff_t2 = OpIdentity>
	auto cellular_fe(Sp const& solver, OpTerm<OpIdentity, G> const& term, coeff_t1 const& c1 = coeff_t1{})
	{
		return expr::make_integer<30>() / (c1 * c1) * term * term * expr::pow<2>(expr::symbols::one - term)
			+ OpChain(OpOperatorChain(OpIdentity{}, OpIdentity{}), expr::pow<2>(expr::make_operator_derivative<1>(solver)(term)));
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
		using parent_model::solver;

		template<typename... As>
		auto make_equations(As const& ...as) const
		{
			((..., expr::printe(as.second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_model>::value>(parent_model::systems_tuple(), as...);
		}

	};


	template<typename parent_trait>
	struct MakeEquationProvisional : parent_trait
	{
		using parent_trait::parent_trait;
		using parent_trait::temp;
		using parent_trait::solver;

		template<typename... As>
		auto make_equations(As const& ...as) const
		{
			((..., expr::printe(as.second, "given equation")));
			return solver.template form_expr_all<model_num_parameters<parent_trait>::value>(forward_systems(), as...);
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
			return forward_systems(parent_trait::systems_tuple(), temp._s);
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

	template<typename... Ts>
	struct ArrayType
	{
		//! Specialized model representing the generalized phase field problem.
		/*!
		 * The equation of motion is created on object
		 * initialization, and the model implements update and equation functions.
		 */
		template<template<typename> typename Eq,
			typename eq_type = Eq<symphas::internal::MakeEquation<ArrayModel<D, Sp, Ts...>>>>
		struct Specialized;
	};
};


namespace symphas::internal
{
	template<typename T, size_t Dm, typename Sp, typename... Ts>
	struct expand_types_to_model;

	template<typename... field_types, size_t Dm, typename Sp>
	struct expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp>
	{
		using type = typename ModelApplied<Dm, Sp>::template OpTypes<field_types...>;
	};

	template<typename... field_types, size_t Dm, typename Sp>
	struct expand_types_to_model<symphas::lib::types_list<field_array_t<void>, field_types...>, Dm, Sp>
	{
		using type = typename ModelApplied<Dm, Sp>::template ArrayType<field_types...>;
	};

	template<typename... field_types, size_t Dm, typename Sp, typename... Ts>
	struct expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp, symphas::internal::parameterized::VECTOR, Ts...>
	{
		using type = typename expand_types_to_model<
			symphas::lib::types_list<field_types..., symphas::internal::parameterized::VECTOR_D<Dm>>, 
			Dm, Sp, Ts...>::type;
	};

	template<typename... field_types, size_t Dm, typename Sp, typename T0, typename... Ts>
	struct expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp, T0, Ts...>
	{
		using type = typename expand_types_to_model<symphas::lib::types_list<field_types..., T0>, Dm, Sp, Ts...>::type;
	};

	template<typename... field_types, size_t Dm, typename Sp, typename T0, typename... Ts>
	struct expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp, field_array_t<T0>, Ts...>
	{
		using type = typename expand_types_to_model<symphas::lib::types_list<field_array_t<void>, field_types..., field_array_t<T0>>, Dm, Sp, Ts...>::type;
	};

	template<typename... field_types, size_t Dm, typename Sp, typename T0, typename... Ts>
	struct expand_types_to_model<symphas::lib::types_list<field_array_t<void>, field_types...>, Dm, Sp, field_array_t<T0>, Ts...>
	{
		using type = typename expand_types_to_model<symphas::lib::types_list<field_array_t<void>, field_types..., field_array_t<T0>>, Dm, Sp, Ts...>::type;
	};

	template<typename... field_types, size_t Dm, typename Sp, typename... Ts, typename... Rest>
	struct expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp, symphas::lib::types_list<Ts...>, Rest...>
	{
		using type = typename expand_types_to_model<symphas::lib::types_list<field_types...>, Dm, Sp, Ts..., Rest...>::type;
	};
	
	template<size_t Dm, typename Sp, typename... Ts>
	using expand_types_to_model_t = typename expand_types_to_model<symphas::lib::types_list<>, Dm, Sp, Ts...>::type;
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
	//using parent_type::parent_type;

    using this_type = typename ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized<Eq, eq_type>;
	using eqs = typename std::invoke_result_t<decltype(&parent_type::make_equations), parent_type>;
	eqs equations;

	Specialized(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_type(coeff, num_coeff, parameters), equations{ parent_type::make_equations() } {}
	Specialized(symphas::problem_parameters_type const& parameters) : Specialized(nullptr, 0, parameters) {}

	/*ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized<Eq, eq_type>& operator=(
		ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized<Eq, eq_type> other)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(this), *static_cast<parent_type*>(&other));
		swap(equations, other.equations);
	}*/



	void update(double time)
	{
		M::update_systems(time);
	}

	void equation()
	{
		equation(std::make_index_sequence<sizeof...(S)>{});
	}

protected:

	template<size_t... Is>
	void equation(std::index_sequence<Is...>)
	{
		(..., M::solver.equation(std::get<Is>(equations)));
	}

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
	//using parent_type::parent_type;
	using parent_type::temp;

	using eqs = typename std::invoke_result_t<decltype(&parent_type::make_equations), parent_type>;
	using prs = typename std::invoke_result_t<decltype(&parent_type::make_provisionals), parent_type>;

    using this_type = typename ModelApplied<D, Sp>::template OpTypes<S...>::template ProvTypes<P...>::template Specialized<Eq, Pr, pr_type, eq_type>;

	prs provisionals;
	eqs equations;

	Specialized(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_type(coeff, num_coeff, parameters), provisionals{ parent_type::make_provisionals() }, equations{ parent_type::make_equations() } {}
	Specialized(symphas::problem_parameters_type const& parameters) : Specialized(nullptr, 0, parameters) {}

	/*this_type& operator=(this_type other)
	{
		using std::swap;
		
		swap(*static_cast<parent_type*>(this), *static_cast<parent_type*>(&other));
		swap(equations, other.equations);
		swap(provisionals, other.provisionals);
	}*/



	void update(double time)
	{
		M::update_systems(time);
		M::solver.evaluate(provisionals);
		temp.update_systems(parent_type::lastindex, time);
	}

	void equation()
	{
		equation(std::make_index_sequence<sizeof...(S)>{});
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

	template<size_t... Is>
	void equation(std::index_sequence<Is...>)
	{
		(..., M::solver.equation(std::get<Is>(equations)));
	}

	Specialized() : parent_type(), provisionals{ parent_type::make_provisionals() }, equations{ parent_type::make_equations() } {}

};


template<size_t D, typename Sp>
template<typename... Ts>
template<template<typename> typename Eq, typename eq_type>
struct ModelApplied<D, Sp>::ArrayType<Ts...>::Specialized : eq_type
{
	using M = ArrayModel<D, Sp, Ts...>;
	using M::solver;

	template<typename Sp0>
	using eq_type_solver = Eq<symphas::internal::MakeEquation<ArrayModel<D, Sp0, Ts...>>>;

	template<template<template<typename> typename, typename> typename SpecializedModel, typename Sp0>
	using impl_type = SpecializedModel<Eq, eq_type_solver<Sp0>>;

	using parent_type = eq_type;
	//using parent_type::parent_type;

    using this_type = typename ModelApplied<D, Sp>::template ArrayType<Ts...>::template Specialized<Eq, eq_type>;

	using eqs = std::invoke_result_t<decltype(&parent_type::make_equations), parent_type>;
	eqs equations;

	Specialized(double const* coeff, size_t num_coeff, symphas::problem_parameters_type const& parameters) :
		parent_type(coeff, num_coeff, parameters), equations{ parent_type::make_equations() } {}
	Specialized(symphas::problem_parameters_type const& parameters) : Specialized(nullptr, 0, parameters) {}

    /*
	this_type& operator=(this_type other)
	{
		using std::swap;

		swap(*static_cast<parent_type*>(this), *static_cast<parent_type*>(&other));
		swap(equations, other.equations);
	}*/

	void update(double time)
	{
		M::update_systems(time);
	}

	void equation()
	{
#		ifndef DEBUG
#		pragma omp parallel for
#		endif
		for (iter_type i = 0; i < parent_type::len; ++i)
		{
			M::solver.equation(equations[i]);
		}
	}

	~Specialized() { parent_type::delete_equations(equations); }

protected:

	Specialized() : parent_type(), equations{ parent_type::make_equations() } {}

};


// ****************************************************************************************

ADD_EXPR_TYPE_SYMBOL(diff_F)
ADD_EXPR_TYPE_SYMBOL_TEMPLATE(diff_F_i, (typename I), (I))

namespace symphas::internal
{


	template<int N, int P>
	constexpr auto dFE_var(expr::symbols::i_<N, P>)
	{
		return expr::make_term(expr::symbols::diff_F_i<expr::symbols::i_<N, P>>{});
	}

	template<size_t N>
	constexpr auto dFE_var()
	{
		return expr::make_term<N>(expr::symbols::diff_F_symbol{});
	}

	template<typename F, typename List>
	struct select_dFE_impl;

	template<typename... dFE_ts>
	struct select_dFE_impl<symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<dFE_ts...>;
	};

	template<typename... dFE_ts, size_t N, typename... Rest>
	struct select_dFE_impl<symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<Variable<N, expr::symbols::diff_F_symbol>, Rest...>>
	{
		using type = typename select_dFE_impl<symphas::lib::types_list<dFE_ts..., Variable<N, expr::symbols::diff_F_symbol>>, symphas::lib::types_list<Rest...>>::type;
	};

	template<typename... dFE_ts, typename T, typename... Rest>
	struct select_dFE_impl<symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_dFE_impl<symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<typename List>
	using select_dFE_t = typename select_dFE_impl<symphas::lib::types_list<>, List>::type;


	template<typename I, typename F, typename List>
	struct select_dFE_i_impl;

	template<int N, int P, typename... dFE_ts>
	struct select_dFE_i_impl<expr::symbols::i_<N, P>, symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<>>
	{
		using type = symphas::lib::types_list<dFE_ts...>;
	};

	template<int N, int P, typename... dFE_ts, int P0, typename... Rest>
	struct select_dFE_i_impl<expr::symbols::i_<N, P>, 
		symphas::lib::types_list<dFE_ts...>, 
		symphas::lib::types_list<expr::symbols::diff_F_i_symbol<expr::symbols::i_<N, P0>>, Rest...>>
	{
		using type = typename select_dFE_i_impl<
			expr::symbols::i_<N, P>, 
			symphas::lib::types_list<dFE_ts..., expr::symbols::diff_F_i_symbol<expr::symbols::i_<N, P0>>>, 
			symphas::lib::types_list<Rest...>>::type;
	};

	template<int N, int P, typename... dFE_ts, typename T, typename... Rest>
	struct select_dFE_i_impl<expr::symbols::i_<N, P>, symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<T, Rest...>>
	{
		using type = typename select_dFE_i_impl<expr::symbols::i_<N, P>, symphas::lib::types_list<dFE_ts...>, symphas::lib::types_list<Rest...>>::type;
	};

	template<typename I, typename List>
	using select_dFE_i_t = typename select_dFE_i_impl<I, symphas::lib::types_list<>, List>::type;


	using namespace expr::symbols;

	template<DynamicType dynamic>
	using dynamics_key_t = dynamics_list<dynamic>;

	//template<DynamicType dynamic>
	//using all_dynamics_key_t = symphas::lib::types_list<dynamics_key_t<dynamic>>;


	template<symphas::internal::DynamicType dynamics>
	struct apply_dynamics;

	template<>
	struct apply_dynamics<symphas::internal::DynamicType::CONSERVED>
	{
		template<typename U_D, typename F, typename Sp>
		auto operator()(U_D const& dop, F const& dfe, Sp const& solver)
		{
			return (dop = expr::make_operator_derivative<2>(solver) * dfe);
		}
	};

	template<>
	struct apply_dynamics<symphas::internal::DynamicType::NONCONSERVED>
	{
		template<typename U_D, typename F, typename Sp>
		auto operator()(U_D const& dop, F const& dfe, Sp const& solver)
		{
			return (dop = -dfe);
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
		special_dynamics(OpOperator<E> const& e) : e{ *static_cast<E const*>(&e) } {}
		special_dynamics(I) {}
		special_dynamics() {}

		E e;
	};

	// Specifies the dynamics for a single field.
	template<size_t N, symphas::internal::DynamicType dynamic>
	struct special_dynamics<N, void, dynamics_key_t<dynamic>> 
	{
	};

	// Specifies the dynamics for a single field.
	template<size_t N, typename E>
	struct special_dynamics<N, void, E>
	{
		special_dynamics(OpExpression<E> const& e) : e{ *static_cast<E const*>(&e) } {}
		special_dynamics(OpOperator<E> const& e) : e{ *static_cast<E const*>(&e) } {}
		special_dynamics() : e{} {}

		E e;
	};

	// Specifies the dynamics for all fields.
	template<size_t N, typename I>
	struct special_dynamics<N, I, void>
	{
		special_dynamics(I) {}
		special_dynamics(std::index_sequence<N>, I) {}
		special_dynamics() {}

		template<typename E0>
		auto operator()(OpExpression<E0> const& e)
		{
			return special_dynamics<N, I, E0>(*static_cast<E0 const*>(&e));
		}

		template<typename E0>
		auto operator()(OpOperator<E0> const& e)
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
		special_dynamics() {}

		template<typename E0>
		auto operator()(OpExpression<E0> const& e)
		{
			return special_dynamics<N, void, E0>(*static_cast<E0 const*>(&e));
		}

		template<typename E0>
		auto operator()(OpOperator<E0> const& e)
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





	template<typename T>
	struct special_dynamics_select;

	template<typename T>
	struct special_dynamics_select<T const> : special_dynamics_select<T> {};

	template<int N, int P>
	struct special_dynamics_select<expr::symbols::i_<N, P>>
	{
		template<int, typename I>
		auto select(I)
		{
			return special_dynamics(I{});
		}

	};

	template<>
	struct special_dynamics_select<int>
	{
		template<int N>
		auto select(size_t)
		{
			return special_dynamics(std::index_sequence<N - 1>{});
		}
	};



	template<size_t N, size_t D>
	constexpr bool is_index_divisible(OpFractionLiteral<N, D>)
	{
		return (D == 1);
	}

	template<size_t N, size_t D>
	constexpr bool is_index_divisible(OpNegFractionLiteral<N, D>)
	{
		return (D == 1);
	}

	inline constexpr bool is_index_divisible(OpIdentity)
	{
		return true;
	}

	inline constexpr bool is_index_divisible(OpNegIdentity)
	{
		return true;
	}

	template<size_t N, int N0, int P0>
	constexpr bool is_index_divisible(expr::symbols::i_<N0, P0>)
	{
		return true;
	}

	template<size_t N, typename index_t, typename E>
	constexpr bool is_index_divisible(symphas::lib::types_list<index_t>, OpExpression<E> const& e)
	{
		auto [lhs, rhs] = separate_index(index_t{}, *static_cast<E const*>(&e) - expr::val<N>);
		return is_index_divisible(rhs / lhs);
	}

	template<size_t N, typename E>
	constexpr bool is_index_divisible(OpExpression<E> const& e)
	{
		using index_t = symphas::internal::select_unique_i_<expr::op_types_t<E>>;
		return is_index_divisible<N>(index_t{}, *static_cast<E const*>(&e));
	}

	template<size_t N, typename I>
	struct dynamics_index_compare
	{
		static const bool flag = is_index_divisible<N>(I{});
	};

	template<size_t N, int N0, int P0, int M>
	struct dynamics_index_compare<N, expr::symbols::index_eq_N<expr::symbols::i_<N0, P0>, M>>
	{
		static const bool flag = (N >= M);
	};

	template<size_t N, int N0, int P0, int M0, int M1>
	struct dynamics_index_compare<N, expr::symbols::index_NtN<expr::symbols::i_<N0, P0>, M0, M1>>
	{
		static const bool flag = (N >= M0 && N <= M1);
	};

	template<size_t N, typename L, typename R>
	struct dynamics_rule_compare
	{
		static const bool value = true;		//<! `true` means L rule has greater priority than R.
	};

	template<size_t N, typename R>
	struct dynamics_rule_compare<N, void, R>
	{
		static const bool value = false;
	};

	template<size_t N, size_t N0, typename E, typename R>
	struct dynamics_rule_compare<N, special_dynamics<N0, void, E>, R>
	{
		static const bool value = (N == N0);
	};

	template<size_t N, typename L, typename E>
	struct dynamics_rule_compare<N, L, special_dynamics<N, void, E>>
	{
		static const bool value = false;
	};

	template<size_t N, typename E>
	struct dynamics_rule_compare<N, void, special_dynamics<N, void, E>>
	{
		static const bool value = false;
	};

	template<size_t N, typename L, size_t N0, typename I, typename E>
	struct dynamics_rule_compare<N, L, special_dynamics<N0, I, E>>
	{
		static const bool value = !dynamics_index_compare<N, I>::value;
	};

	template<size_t N, size_t N0, typename I, typename E>
	struct dynamics_rule_compare<N, void, special_dynamics<N0, I, E>>
	{
		static const bool value = false;
	};

	template<size_t N, size_t N0, typename I, typename E, typename R>
	struct dynamics_rule_compare<N, special_dynamics<N0, I, E>, R>
	{
		static const bool value = dynamics_index_compare<N, I>::value;
	};

	template<size_t N, size_t N0, typename I, typename E0, typename E1>
	struct dynamics_rule_compare<N, special_dynamics<N0, I, E0>, special_dynamics<N, void, E1>>
	{
		static const bool value = false;
	};

	template<size_t N, size_t N0, size_t N1, typename I, typename E0, typename E1>
	struct dynamics_rule_compare<N, special_dynamics<N0, void, E0>, special_dynamics<N1, I, E1>>
	{
		static const bool value = (N0 != N);
	};

	template<size_t N, size_t N0, size_t N1, typename E0, typename E1>
	struct dynamics_rule_compare<N, special_dynamics<N0, void, E0>, special_dynamics<N1, void, E1>>
	{
		static const bool value = (N0 == N);
	};
	
	template<size_t N, size_t N0, typename E0, typename E1>
	struct dynamics_rule_compare<N, special_dynamics<N0, void, E0>, special_dynamics<N, void, E1>>
	{
		static const bool value = false;
	};

	template<size_t N, size_t N0, size_t N1, typename I0, typename I1, typename E0, typename E1>
	struct dynamics_rule_compare<N, special_dynamics<N0, I0, E0>, special_dynamics<N1, I1, E1>>
	{
		static const bool value = dynamics_index_compare<N, I0>::value;
	};

	//template<size_t N, size_t N0, typename I0, typename E0, typename E1>
	//struct dynamics_rule_compare<N, special_dynamics<N0, I0, E0>, special_dynamics<N, I0, E1>>
	//{
	//	static const bool value = dynamics_index_compare<N, I0>::value;
	//};

	//template<size_t N, size_t N0, size_t N1, typename E0, typename E1>
	//struct dynamics_rule_compare<N, special_dynamics<N0, void, E0>, special_dynamics<N1, void, E1>>
	//{
	//	static const bool value = (N1 == N);
	//};

	template<size_t N, size_t P, typename... Ts>
	struct dynamics_rule_search_impl;

	template<size_t N, size_t P, typename R>
	struct dynamics_rule_search_impl<N, P, R, std::index_sequence<>, symphas::lib::types_list<>>
	{
		static const size_t value = P;
	};

	template<size_t N, size_t P, typename R, size_t I0, size_t... Is, typename T, typename... Rest>
	struct dynamics_rule_search_impl<N, P, R, std::integer_sequence<bool, false>, std::index_sequence<I0, Is...>, symphas::lib::types_list<T, Rest...>>
	{
		static const size_t value = dynamics_rule_search_impl<N, I0, T, 
			std::index_sequence<Is...>, symphas::lib::types_list<Rest...>>::value;
	};

	template<size_t N, size_t P, typename R, size_t I0, size_t... Is, typename T, typename... Rest>
	struct dynamics_rule_search_impl<N, P, R, std::integer_sequence<bool, true>, std::index_sequence<I0, Is...>, symphas::lib::types_list<T, Rest...>>
	{
		static const size_t value = dynamics_rule_search_impl<N, P, R, 
			std::index_sequence<Is...>, symphas::lib::types_list<Rest...>>::value;
	};

	template<size_t N, size_t P, typename R, size_t I0, size_t... Is, typename T, typename... Rest>
	struct dynamics_rule_search_impl<N, P, R, std::index_sequence<I0, Is...>, symphas::lib::types_list<T, Rest...>>
	{
	protected:
		static const bool flag = dynamics_rule_compare<N, R, T>::value;

	public:
		static const size_t value = dynamics_rule_search_impl<N, P, R, 
			std::integer_sequence<bool, flag>, std::index_sequence<I0, Is...>, symphas::lib::types_list<T, Rest...>>::value;
	};


	template<size_t N, typename T>
	constexpr size_t dynamics_rule_N = dynamics_rule_search_impl<N, 0, void, std::make_index_sequence<symphas::lib::types_list_size<T>::value>, T>::value;


	template<typename E0, typename F, typename... S>
	auto substitute_for_dynamics(
		symphas::lib::types_list<>,
		E0 const& dyn0,
		F const& fe,
		std::tuple<S...> const& ops)
	{
		return dyn0;
	}

	template<size_t N, size_t... Ns, typename E0, typename F, typename... S>
	auto substitute_for_dynamics(
		symphas::lib::types_list<
			Variable<N, expr::symbols::diff_F_symbol>,
			Variable<Ns, expr::symbols::diff_F_symbol>...>,
		E0 const& dyn0,
		F const& fe,
		std::tuple<S...> const& ops)
	{
		return expr::transform::swap_grid<Variable<N, expr::symbols::diff_F_symbol>, Variable<Ns, expr::symbols::diff_F_symbol>...>
			(dyn0,
				expr::make_functional_derivative(fe, std::get<0>(ops)),
				expr::make_functional_derivative(fe, std::get<symphas::lib::index_of_value<size_t, Ns, Ns...> + 1>(ops))...);
	}
	
	template<size_t NN, typename E0>
	auto substitute_for_dynamics(
		symphas::lib::types_list<>,
		E0 const& dyn0)
	{
		return dyn0;
	}


	template<size_t NN, int N, int... Ps, typename E0>
	auto substitute_for_dynamics(
		symphas::lib::types_list<expr::symbols::diff_F_i_symbol<expr::symbols::i_<N, Ps>>...>,
		E0 const& dyn0)
	{
		return expr::transform::swap_grid<expr::symbols::diff_F_i_symbol<expr::symbols::i_<N, Ps>>...>
			(dyn0, expr::make_term<NN + size_t(Ps)>(expr::symbols::diff_F_symbol{})...);
	}

	template<size_t N, int N0, int... Q0s, int... P0s, typename E0, typename model_t>
	auto substitute_ops(
		symphas::lib::types_list<expr::symbols::i_<N0, Q0s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<N0, P0s>>...>,
		E0 const& eq,
		model_t const& model)
	{
		if constexpr (model_num_parameters<model_t>::value > 0)
		{
			auto e0 = expr::transform::swap_grid<OpCoeffSwap<expr::symbols::i_<N0, Q0s>>...>
				(eq, expr::val<size_t(Q0s) + 1>...);
			return expr::transform::swap_grid<expr::symbols::i_<N0, Q0s>..., expr::symbols::v_id_type<expr::symbols::i_<N0, P0s>>...>
				(e0, expr::val<size_t(Q0s) + 1>..., model.template op<N + size_t(P0s)>()...);
		}
		else
		{
			auto e0 = expr::transform::swap_grid<OpCoeffSwap<expr::symbols::i_<N0, Q0s>>...>
				(eq, model.index);
			return expr::transform::swap_grid<expr::symbols::i_<N0, Q0s>..., expr::symbols::v_id_type<expr::symbols::i_<N0, P0s>>...>
				(e0, model.index + expr::val<size_t(Q0s) + 1>..., model.template op<N + size_t(P0s)>()...);
		}
	}

	template<size_t N, typename E0, typename model_t>
	auto substitute_ops(
		symphas::lib::types_list<>,
		symphas::lib::types_list<>,
		E0 const& eq,
		model_t const& model)
	{
		return eq;
	}

	template<size_t N, typename I0, typename model_t, typename E0>
	auto substitute_ops(model_t const& model, E0 const& with_fes)
	{
		using v_types = select_v_i_<I0, expr::op_types_t<E0>>;
		using i_types = select_all_i_<I0, expr::op_types_t<OpAdd<E0, OpTerm<OpIdentity, I0>>>>;
		return substitute_ops<N>(i_types{}, v_types{}, with_fes, model);
	}

	template<typename model_t, size_t N, size_t... Ns>
	auto all_ops(model_t const& model, symphas::lib::types_list<
		Variable<N, expr::symbols::diff_F_symbol>,
		Variable<Ns, expr::symbols::diff_F_symbol>...>)
	{
		return std::make_tuple(model.template op<N>(), model.template op<Ns>()...);
	}

	template<typename model_t>
	auto all_ops(model_t const& model, symphas::lib::types_list<>)
	{
		return std::make_tuple();
	}

	template<typename S>
	struct apply_special_dynamics;

	template<size_t N, DynamicType dynamic>
	struct apply_special_dynamics<special_dynamics<N, void, dynamics_key_t<dynamic>>>
		: apply_dynamics<dynamic>
	{
		apply_special_dynamics(special_dynamics<N, void, dynamics_key_t<dynamic>>)
			: apply_dynamics<dynamic>() {}
		apply_special_dynamics(dynamics_key_t<dynamic>) {}

		template<size_t N0>
		auto select() const
		{
			return apply_special_dynamics<special_dynamics<N0, void, dynamics_key_t<dynamic>>>(dynamics_key_t<dynamic>{});
		}

		template<typename U_D, typename model_t, typename F, typename Sp>
		auto operator()(U_D const& dop, model_t const& model, OpExpression<F> const& fe, Sp const& solver)
		{
			using df_v_type = Variable<N, expr::symbols::diff_F_symbol>;

			auto evolution = substitute_for_dynamics(
				symphas::lib::types_list<df_v_type>{}, OpTerm<OpIdentity, df_v_type>{},
				*static_cast<F const*>(&fe),
				all_ops(model, symphas::lib::types_list<df_v_type>{}));

			//return apply_dynamics<symphas::internal::DynamicType::NONCONSERVED>{}(dop, evolution, solver);
			return apply_dynamics<dynamic>{}(dop, evolution, solver);
		}
	};
	
	template<size_t N, typename I0, typename E>
	struct apply_special_dynamics<special_dynamics<N, I0, E>>
	{
		apply_special_dynamics(E const& e) : e{ e } {}
		apply_special_dynamics(special_dynamics<N, I0, E> const& e) : e{ e.e } {}

		template<size_t N0>
		auto select() const
		{
			return apply_special_dynamics<special_dynamics<N0, I0, E>>(e);
		}

		template<typename U_D, typename model_t, typename F, typename Sp>
		auto operator()(U_D const& dop, model_t const& model, OpExpression<F> const& fe, Sp const& solver)
		{
			using dfi_v_types = select_dFE_i_t<I0, expr::op_types_t<E>>;
			auto subbed_dfs = substitute_for_dynamics<N>(dfi_v_types{}, e);
			auto subbed_ops = substitute_ops<N, I0>(model, subbed_dfs);

			using df_v_types = select_dFE_t<expr::op_types_t<decltype(subbed_dfs)>>;
			auto evolution = substitute_for_dynamics(df_v_types{}, subbed_ops, 
				*static_cast<F const*>(&fe), 
				all_ops(model, df_v_types{}));

			return (dop = evolution);
		}

		

		E e;
	};

	template<size_t N, typename E>
	struct apply_special_dynamics<special_dynamics<N, void, E>>
	{
		apply_special_dynamics(E const& e) : e{ e } {}
		apply_special_dynamics(special_dynamics<N, void, E> const& e) : e{ e.e } {}

		template<size_t N0>
		auto select() const
		{
			return apply_special_dynamics<special_dynamics<N0, void, E>>(e);
		}

		template<typename U_D, typename model_t, typename F, typename Sp>
		auto operator()(U_D const& dop, model_t const& model, OpExpression<F> const& fe, Sp const& solver)
		{
			using df_v_types = select_dFE_t<expr::op_types_t<E>>;

			auto evolution = substitute_for_dynamics(
				df_v_types{}, e,
				*static_cast<F const*>(&fe),
				all_ops(model, df_v_types{}));

			return (dop = evolution);
		}

		E e;
	};

	template<typename S>
	apply_special_dynamics(S) -> apply_special_dynamics<S>;
	template<DynamicType dynamic>
	apply_special_dynamics(dynamics_key_t<dynamic>) -> apply_special_dynamics<special_dynamics<0, void, dynamics_key_t<dynamic>>>;

	template<DynamicType dynamic, size_t I>
	constexpr DynamicType dynamic_i = dynamic;

	
	template<typename M, size_t D, typename Sp, typename... Ts>
	auto get_parameterized_values(M const& model0, Model<D, Sp, Ts...> const& model)
	{
		return new parameterized_type<void, void>[sizeof...(Ts)] { parameterized_type<Ts>(model0)()... };
	}

	template<typename M, size_t D, typename Sp, typename... Ts>
	auto get_parameterized_values(M const& model0, ArrayModel<D, Sp, Ts...> const& model)
	{
		return new parameterized_type<void, void>[sizeof...(Ts)] {
			parameterized_type<Ts>(model0)()... };
	}

	struct param_matrix_factory
	{
		param_matrix_factory() : coeff{ nullptr }, num_fields{ 0 }, num_coeff{ 0 } {}

		template<typename M>
		param_matrix_factory(M const& model, len_type num_fields, len_type num_coeff) :
			coeff{ generate_coeff_matrix(model, model, num_coeff, num_fields) }, 
			num_fields{ num_fields }, num_coeff{ num_coeff } {}

		param_matrix_factory(param_matrix_factory const& other) :
			coeff{ (other.num_fields * other.num_coeff > 0) ? new double[other.num_fields * other.num_coeff] : nullptr},
			num_fields{ other.num_fields }, num_coeff{ other.num_coeff } 
		{
			std::copy(other.coeff, other.coeff + other.num_fields * other.num_coeff, coeff);
		}

		param_matrix_factory(param_matrix_factory&& other) : param_matrix_factory()
		{
			swap(*this, other);
		}

		param_matrix_factory& operator=(param_matrix_factory other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(param_matrix_factory& first, param_matrix_factory& second)
		{
			using std::swap;
			swap(first.coeff, second.coeff);
			swap(first.num_coeff, second.num_coeff);
			swap(first.num_fields, second.num_fields);
		}

		auto operator()(iter_type n) const
		{
			return expr::make_coeff(coeff + (n - 1), num_fields, num_coeff);
		}

	protected:

		template<typename M, size_t D, typename Sp, typename... Ts>
		double* generate_coeff_matrix(
			M const& model0, Model<D, Sp, Ts...> const& model,
			len_type num_coeff, len_type num_fields)
		{
			double* arr = new double[num_fields * num_coeff];
			auto* parameterized_values = get_parameterized_values(model0, model);

			for (iter_type i = 0; i < num_fields; ++i)
			{
				for (iter_type n = 0; n < num_coeff; ++n)
				{
					iter_type index = i * num_coeff + n;
					arr[index] = parameterized_values[i][n];
				}
			}
			delete[] parameterized_values;
			return arr;
		}

		template<typename M, size_t D, typename Sp, typename... Ts>
		double* generate_coeff_matrix(
			M const& model0, ArrayModel<D, Sp, Ts...> const& model,
			len_type num_coeff, len_type num_fields)
		{
			double* arr = new double[num_fields * num_coeff];
			auto* parameterized_values = get_parameterized_values(model0, model);

			//constexpr size_t N = M::num_array_types;
			iter_type offset = 0;
			for (iter_type type_ind = 0; type_ind < sizeof...(Ts); ++type_ind)
			{
				for (iter_type i = 0; i < model.num_fields(type_ind); ++i)
				{
					for (iter_type n = 0; n < num_coeff; ++n)
					{
						iter_type index = i * num_coeff + n;
						arr[index + offset] = parameterized_values[type_ind][n];
					}
				}
				offset += model.num_fields(type_ind) * num_coeff;
			}

			//for (iter_type i = model.len - sizeof...(Ts); i < num_fields; ++i)
			//{
			//	for (iter_type n = 0; n < num_coeff; ++n)
			//	{
			//		iter_type index = i * num_coeff + n;
			//		arr[index] = parameterized_values[i - (model.len - sizeof...(Ts)) + 1][n];
			//	}
			//}

			delete[] parameterized_values;
			return arr;
		}

	public:

		double* coeff;
		len_type num_fields;
		len_type num_coeff;
	};


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

	using names_t = model_field_name<enclosing_type>;

	template<size_t I>
	auto op() const
	{
		return const_cast<TraitEquation<enclosing_type, parent_trait>*>(this)->template _op<I>();
	}

	template<size_t I>
	auto dop() const
	{
		return const_cast<TraitEquation<enclosing_type, parent_trait>*>(this)->template _dop<I>();
	}

	template<int N, int P>
	auto op_n(expr::symbols::i_<N, P>) const
	{
		return expr::symbols::v_<expr::symbols::i_<N, P>>{};
	}

	auto param(size_t I) const
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

	template<size_t N>
	auto param_matrix() const
	{
		static symphas::internal::param_matrix_factory coeff_matrix(*this, Sn, N);
		return coeff_matrix;
	}

	auto c(size_t I) const
	{
		return param(I - 1);
	}

    template<expr::NoiseType nt, typename T, typename... T0s>
    auto make_noise(T0s&& ...args) const
    {
        return symphas::internal::parameterized::NOISE<nt, symphas::internal::parameterized::dimensionalized_t<Dm, T>, Dm>(
			parent_trait::template system<0>().info,
			grid::get_data_domain(parent_trait::template system<0>().as_grid()), 
			&solver.dt)(std::forward<T0s>(args)...);
    }

    template<expr::NoiseType nt, size_t Z, typename G, typename... T0s>
    auto make_noise(OpTerm<OpIdentity, Variable<Z, G>>, T0s&& ...args) const
    {
        using T = model_field_t<parent_trait, Z>;
        return symphas::internal::parameterized::NOISE<nt, symphas::internal::parameterized::dimensionalized_t<Dm, T>, Dm>(
			parent_trait::template system<Z>().info,
			grid::get_data_domain(parent_trait::template system<Z>().as_grid()), 
			&solver.dt)(std::forward<T0s>(args)...);
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

	template<typename... dynamics_ts, typename E>
	auto generate_equations(std::tuple<dynamics_ts...> const& dynamics, OpExpression<E> const& e) const
	{
		expr::printe(*static_cast<E const*>(&e), "free energy");
		return generate_equations_apply(*static_cast<E const*>(&e), dynamics, seq_t{});
	}

	template<typename... As>
	auto make_equations(As&& ...as) const
	{
		((..., expr::printe(as.second, "given equation")));
		return solver.template form_expr_all<model_num_parameters<parent_trait>::value>(parent_trait::systems_tuple(), as...);
	}

protected:


	template<size_t I>
	auto _op()
	{
		return expr::make_term<I>(
			NamedData(
				parent_trait::template grid<I>(),
				names_t{}(I)
			));
	}

	template<size_t I>
	auto _dop()
	{
		return OpLHS(expr::as_variable<I>(parent_trait::template system<I>()));
	}


	template<typename E, typename... dynamics_ts, size_t... Is>
	auto generate_equations_apply(
		OpExpression<E> const& e,
		std::tuple<dynamics_ts...> const& dynamics,
		std::index_sequence<Is...>) const
	{
		using namespace symphas::internal;
		auto dys = std::make_tuple(dynamics_ts{}...);

		return parent_trait::make_equations(
			apply_special_dynamics(
				std::get<dynamics_rule_N<Is, symphas::lib::types_list<dynamics_ts...>>>(dynamics)).template select<Is>()(
					dop<Is>(),
					*this,
					*static_cast<E const*>(&e),
					solver)...);
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
template<template<typename> typename enclosing_type, size_t D, typename Sp, typename... Ts>
struct TraitEquation<enclosing_type, symphas::internal::MakeEquation<ArrayModel<D, Sp, Ts...>>> : 
	symphas::internal::MakeEquation<ArrayModel<D, Sp, Ts...>>
{
	using parent_trait = symphas::internal::MakeEquation<ArrayModel<D, Sp, Ts...>>;
	//using parent_trait::parent_trait;
	using parent_trait::solver;
	using parent_trait::len;

	using names_t = model_field_name<enclosing_type>;
	iter_type index_value;
	DynamicIndex index;

	template<typename... T0s>
	TraitEquation(T0s&&... args) : parent_trait(std::forward<T0s>(args)...), index_value{ 0 }, index{ index_value, 0, len - 1 } {}

	template<size_t I = 0>
	auto op() const
	{
		return const_cast<TraitEquation<enclosing_type, parent_trait>*>(this)->template _op<I>();
	}

	template<size_t I = 0>
	auto dop() const
	{
		return const_cast<TraitEquation<enclosing_type, parent_trait>*>(this)->template _dop<I>();
	}

	template<int N, int P>
	auto op_n(expr::symbols::i_<N, P>) const
	{
		return expr::symbols::v_<expr::symbols::i_<N, P>>{};
	}

	auto param(size_t I) const
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

	template<size_t N>
	auto param_matrix() const
	{
		static symphas::internal::param_matrix_factory coeff_matrix(*this, parent_trait::len, N);
		return coeff_matrix;
	}

	auto c(size_t I) const
	{
		return param(I - 1);
	}

	template<expr::NoiseType nt, typename T, typename... T0s>
	auto make_noise(T0s&& ...args) const
	{
		return symphas::internal::parameterized::NOISE<nt, T, D>(
			parent_trait::template system<0>().info, &solver.dt)
			(std::forward<T0s>(args)...);
	}

	template<expr::NoiseType nt, size_t Z, typename G, typename... T0s>
	auto make_noise(OpTerm<OpIdentity, Variable<Z, G>>, T0s&& ...args) const
	{
		using T = model_field_t<parent_trait, Z>;
		return symphas::internal::parameterized::NOISE<nt, T, D>(
			parent_trait::template system<Z>().info, &solver.dt)
			(std::forward<T0s>(args)...);
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

	template<typename E, typename dynamics_t>
	auto generate_equations(
		std::tuple<dynamics_t> const& dynamics,
		OpExpression<E> const& e) const
	{
		using namespace symphas::internal;

		expr::printe(*static_cast<E const*>(&e), "free energy");

		return make_equations(
			apply_special_dynamics(std::get<0>(dynamics))(
				dop(),
				*this,
				*static_cast<E const*>(&e),
				solver));
	}

	template<typename Dd, typename Ee>
	auto make_equations(std::pair<Dd, Ee> const& a) const
	{
		using this_type = TraitEquation<enclosing_type, symphas::internal::MakeEquation<ArrayModel<D, Sp, Ts...>>>;
		using scheme_type = std::invoke_result_t<decltype(&this_type::_make_equations<Dd, Ee>), this_type, std::pair<Dd, Ee>, iter_type>;

		expr::printe(std::get<1>(a), "given equation");

		void* ptr = operator new[](parent_trait::len * sizeof(scheme_type));
		scheme_type* eqns = static_cast<scheme_type*>(ptr);

		for (iter_type i = 0; i < parent_trait::len; ++i)
		{
			new(&eqns[i]) scheme_type(_make_equations(a, i));
		}
		return eqns;
	}


protected:

	template<typename scheme_type>
	void delete_equations(scheme_type* eqns)
	{
		for (int i = parent_trait::len - 1; i >= 0; --i)
		{
			eqns[i].~scheme_type();
		}
		operator delete[](eqns);
	}

	template<typename Dd, typename Ee>
	auto _make_equations(std::pair<Dd, Ee> a, iter_type i) const
	{
		auto&& [d, e] = a;
		expr::fix_index(e, index == i);
		expr::fix_index(d, index == i);

		return solver.template form_expr_one<model_num_parameters<parent_trait>::value>(parent_trait::systems_tuple(), a);
	}

	template<size_t I = 0>
	auto _op()
	{
		return expr::make_term_dynamic(index,
			NamedData(parent_trait::systems(), names_t{}())
			);
	}

	template<size_t I = 0>
	auto _dop()
	{
		return OpLHS(DynamicVariable(index, parent_trait::systems()));
	}

	//template<typename L, typename E>
	//auto update_placeholders(std::pair<L, E> const& equation) const
	//{
	//	using swap_t = GridSymbol<symphas::internal::non_parameterized_type<S>, D>;
	//	auto&& [l, e] = equation;
	//	return std::make_pair(l, expr::transform::swap_grid<swap_t, expr::symbols::placeholder_N_symbol>(e, op(), index));
	//}

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
	using parent_type::c; \
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
	using parent_type::c; \
	using names_t = typename parent_type::names_t; \
	static const size_t Dm = model_dimension<parent_type>::value; \
	auto make_equations() const { \
		using namespace std; using namespace expr; using namespace expr::symbols; using namespace std::complex_literals; \
		auto [x, y, z] = expr::make_coords<Dm>(DIMENSIONS_OF(0), INTERVALS(0)); \
		auto t = parent_type::get_time_var(); \
		using symphas::internal::parameterized::INT; \
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

 //! Defines a TraitEquation child class used to define dynamical equations.
 /*!
  * The equations of motion for the phase fields are provided as arguments, and
  * a new class of type TraitEquationModel is generated. Since this name is
  * not unique, this macro should be invoked within a namespace.
  *
  * \param ... The equations of the phase fields.
  */
#define EQUATION_FE_TRAIT_DEFINITION(SELECTED_DYNAMICS, ...) \
return parent_type::template generate_equations(std::make_tuple(SINGLE_ARG SELECTED_DYNAMICS), __VA_ARGS__); } \
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
template<> struct model_field_name<model_ ## MODEL_NAME::TraitEquationModel> : model_field_name_builder<model_ ## MODEL_NAME::TraitEquationModel> {}; \
template<> struct model_field_name_format<model_ ## MODEL_NAME::TraitEquationModel> \
{ static const char* value[]; }; \
inline const char* model_field_name_format<model_ ## MODEL_NAME::TraitEquationModel>::value[] = { SINGLE_ARG NAMES };

#define DEFINE_MODEL_FIELD_NAMES_FORMAT(MODEL_NAME, FORMAT) \
template<> struct model_field_name<model_ ## MODEL_NAME::TraitEquationModel> : model_field_name_builder<model_ ## MODEL_NAME::TraitEquationModel> {}; \
template<> struct model_field_name_format<model_ ## MODEL_NAME::TraitEquationModel> \
{ static const char* value; }; \
inline const char* model_field_name_format<model_ ## MODEL_NAME::TraitEquationModel>::value = FORMAT;

//! @}



namespace symphas::internal
{
	template<typename names_t, typename S, typename... Xs>
	struct series_index_selection;

	template<typename, size_t N>
	constexpr size_t ivalue = N;

	template<typename names_t, typename... Ss, typename... Xs>
	struct series_index_selection<names_t, std::tuple<Ss...>, Xs...>
	{
		series_index_selection(names_t, std::tuple<Ss...> const& systems, Xs... xs) : 
			systems{ &const_cast<std::tuple<Ss...>&>(systems) }, limits{ xs... } {}

		template<typename E, size_t... Ls, size_t... Ns>
		auto substitute_systems(OpExpression<E> const& e, std::index_sequence<Ls...>, std::index_sequence<Ns...>)
		{
			auto list = std::make_tuple(
				Term(expr::as_variable<Ns>(
					NamedData(SymbolicData(std::get<Ns>(*systems).as_grid()), names_t{}(Ns))
				))...);
			
			return expr::sum(*static_cast<E const*>(&e))
				.template select<Ls...>(Xs{}...)
				(([&] (auto) { return list; })(Ls)...);
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

	template<typename names_t, typename S, typename... Xs>
	struct series_index_selection<names_t, symphas::internal::field_array_t<S>, Xs...>
	{
		series_index_selection(names_t, std::pair<S*, len_type> const& systems, Xs... xs) :
			systems{ const_cast<S*>(systems.first) }, len{ systems.second }, limits{ xs... } {}

		template<typename E, size_t... Ls>
		auto substitute_systems(OpExpression<E> const& e, std::index_sequence<Ls...>)
		{
			auto list = expr::array_arg(len, NamedData(systems, names_t{}()));

			return expr::sum(*static_cast<E const*>(&e))
				.template select<Ls...>(Xs{}...)
				(([&] (auto) { return list; })(Ls)...);
		}

		template<typename E>
		auto operator()(OpExpression<E> const& e)
		{
			return substitute_systems(
				*static_cast<E const*>(&e),
				symphas::lib::seq_repeating_value_t<sizeof...(Xs), size_t, 0>{});
		}


		S* systems;
		len_type len;
		std::tuple<Xs...> limits;
	};

	template<typename names_t, typename... Ss, typename... Xs>
	series_index_selection(names_t, std::tuple<Ss...>, Xs...) -> series_index_selection<names_t, std::tuple<Ss...>, Xs...>;

	template<typename names_t, typename S, typename... Xs>
	series_index_selection(names_t, std::pair<S*, len_type>, Xs...) -> series_index_selection<names_t, symphas::internal::field_array_t<S>, Xs...>;
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
	auto sum(OpExpression<E> const& e) const
	{
		expr::eval_type_t<E> sum;
		expr::result(OpVoid{}, sum);

		auto sumop = expr::make_term(sum);

		len_type len = expr::data_length(*static_cast<E const*>(&e));
		for (iter_type i = 0; i < len; ++i)
		{
			sum += (*static_cast<E const*>(&e)).eval(i);
		}

		return expr::make_literal(sum);
	}

	template<typename E>
	auto mean(OpExpression<E> const& e) const
	{
		len_type len = expr::data_length(*static_cast<E const*>(&e));
		return sum(*static_cast<E const*>(&e)) * (1. / len);
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
