
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
 * MODULE:  expr
 * PURPOSE: Implements derivatives of the symbolic algebra library.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionaggregates.h"
#include "expressionoperators.h"


template<typename G>
struct SymbolicDerivative
{
	template<Axis ax, size_t O, typename T>
	decltype(auto) generalized_derivative(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}

	template<Axis ax, size_t O, typename T>
	decltype(auto) generalized_directional_derivative(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}

	template<typename T>
	decltype(auto) laplacian(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}

	template<typename T>
	decltype(auto) bilaplacian(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}

	template<Axis ax, typename T>
	decltype(auto) gradlaplacian(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}

	template<Axis ax, typename T>
	decltype(auto) gradient(T&& e, iter_type n) const
	{
		return expr::symbols::Symbol{};
	}
};


template<typename G>
struct SymbolicDerivative<DynamicVariable<G>> : SymbolicDerivative<G>
{
	SymbolicDerivative() : index{} {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index } {}

	DynamicIndex index;
};

template<typename G>
struct SymbolicDerivative<expr::variational_t<DynamicVariable<G>>> : SymbolicDerivative<G>
{
	SymbolicDerivative() : index{} {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index } {}

	DynamicIndex index;
};

template<typename G>
struct SymbolicDerivative<expr::variational_t<DynamicVariable<NamedData<G>>>> : SymbolicDerivative<G>
{

#ifdef PRINTABLE_EQUATIONS

	SymbolicDerivative() : index{}, name{ 0, "" } {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index }, name{ 0, "" } {}
	SymbolicDerivative(DynamicVariable<NamedData<G>> const& var) : index{ var }, name{ 0, var.data.name } {}

	const char* get_name() const
	{
		return name.name;
	}

	char* get_name()
	{
		return name.name;
	}
	
	NamedData<void*> name;

#else

	SymbolicDerivative() : index{} {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index } {}
	SymbolicDerivative(DynamicVariable<NamedData<G>> const& var) : index{ var } {}

#endif

	DynamicIndex index;
};

template<size_t Z, typename G>
struct SymbolicDerivative<expr::variational_t<Variable<Z, NamedData<G>>>> : SymbolicDerivative<G>
{
#ifdef PRINTABLE_EQUATIONS

	SymbolicDerivative() : index{}, name{ 0, "" } {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index }, name{ 0, "" } {}
	SymbolicDerivative(Variable<Z, NamedData<G>> const& var) : index{ var }, name{ 0, var.data.name } {}

	NamedData<void*> name;

	const char* get_name() const
	{
		return name.name;
	}

	char* get_name()
	{
		return name.name;
	}

#else

	SymbolicDerivative() : index{} {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index } {}
	SymbolicDerivative(Variable<Z, NamedData<G>> const& var) : index{ var } {}

#endif

};

template<size_t Z, typename G>
struct SymbolicDerivative<expr::variational_t<Variable<Z, G>>> : SymbolicDerivative<G>
{
	SymbolicDerivative() {}
	SymbolicDerivative(Variable<Z, G> const& var) {}
};

template<typename G>
struct SymbolicDerivative<expr::variational_t<G>> : SymbolicDerivative<G>
{
	SymbolicDerivative() {}
	SymbolicDerivative(G const& var) {}
};

template<typename G, size_t D>
struct SymbolicDerivative<expr::variational_t<GridSymbol<G, D>>> : SymbolicDerivative<G>
{
	SymbolicDerivative() : index{} {}
	SymbolicDerivative(DynamicIndex const& index) : index{ index } {}

	DynamicIndex index;
};

template<typename G>
SymbolicDerivative(DynamicVariable<G>) -> SymbolicDerivative<expr::variational_t<DynamicVariable<G>>>;
template<size_t Z, typename G>
SymbolicDerivative(Variable<Z, G>) -> SymbolicDerivative<expr::variational_t<Variable<Z, G>>>;

template<typename V, typename E, typename G>
using OpFunctionalDerivative = OpDerivative<std::index_sequence<1>, V, E, SymbolicFunctionalDerivative<G>>;

namespace expr
{
    template<typename T>
    constexpr bool is_functional_derivative = false;

    template<typename G>
    constexpr bool is_functional_derivative<
        SymbolicFunctionalDerivative<G>> = true;
}

//! \cond
//! The solver type used by the derivative expression.
template<typename Sp>
using solver_op_type = Sp const&;
//! \endcond

// *************************************************************************************

namespace symphas::internal
{

	//! Implementation of derivative expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct derivative expressions. Wraps the template deduction necessary
	 * to initialize a derivative expression.
	 */
	template<typename Dd>
	struct make_derivative
	{
		//! Constructs the derivative with the identity coefficient.
		template<typename A, typename B>
		static auto get(A&& a, B&& b);

		//! Constructs the derivative applied to an expression.
		template<typename V, typename E, typename Sp>
		static auto get(V const& v, OpExpression<E> const& e, solver_op_type<Sp> solver);

		template<typename V, typename E, typename Sp>
		static auto get(V const& v, OpOperator<E> const& e, solver_op_type<Sp> solver);

		//! Constructs the derivative applied to a variable.
		template<typename V, typename S, typename G, typename Sp>
		static auto get(V const& v, OpTerm<S, G> const& e, solver_op_type<Sp> solver);

		//! Constructs the derivative applied to a constant.
		template<typename V, typename coeff_t, typename Sp, typename std::enable_if_t<(expr::is_coeff<coeff_t> || expr::is_identity<coeff_t>), int> = 0>
		static auto get(V const& v, coeff_t, solver_op_type<Sp> solver);

		template<typename V, typename S, typename G, typename G0>
		static auto get(V const& v, OpTerm<S, G> const& e, SymbolicDerivative<G0> s);

		template<typename V, typename E, typename G0>
		static auto get(V const& v, OpExpression<E> const& e, SymbolicDerivative<G0> s);

		template<typename V, typename E, typename G0>
		static auto get(V const& v, OpOperator<E> const& e, SymbolicDerivative<G0> s);


		//! Constructs the derivative using a grid instead of an expression.
		/*!
		 * Used for the derivative specialization for the OpTerm.
		 */
		template<typename V, typename G, typename Sp>
		static auto get_g(V const& v, G g, solver_op_type<Sp> solver);


		// If passed an OpLiteral, uses its value rather than the object.
		template<typename V, typename E1, typename E2>
		static auto get(OpLiteral<V> const& v, E1&& e1, E2&& e2)
		{
			return get(v.value, std::forward<E1>(e1), std::forward<E2>(e2));
		}

		// If passed an OpLiteral, uses its value rather than the object.
		template<typename V, typename G, typename Sp>
		static auto get_g(OpLiteral<V> const& v, G g, solver_op_type<Sp> solver)
		{
			return get_g(v.value, g, solver);
		}
	};


	//! Implementation of generalized derivative expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct generalized derivative expressions. Wraps the template 
	 * deduction necessary to initialize a derivative expression. As opposed to
	 * the usual derivative, there is no expression associated with generalized
	 * derivatives.
	 * 
	 * \tparam O Order of the generalized derivative.
	 */
	template<size_t O>
	struct make_operator_derivative
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver);

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver);

	};

	template<>
	struct make_operator_derivative<0>
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver)
		{
			return v;
		}

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver)
		{
			return get(OpIdentity{}, solver);
		}
	};

	//! Implementation of generalized derivative expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct generalized derivative expressions. Wraps the template
	 * deduction necessary to initialize a derivative expression. As opposed to
	 * the usual derivative, there is no expression associated with generalized
	 * derivatives.
	 *
	 * \tparam O Order of the generalized derivative.
	 */
	template<Axis ax, size_t O>
	struct make_operator_directional_derivative
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver);

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver);

		// If passed an OpLiteral, uses its value rather than the object.
		template<typename V, typename E1>
		static auto get(OpLiteral<V> const& v, E1&& e1)
		{
			return get(v.value, std::forward<E1>(e1));
		}
	};

	template<Axis ax>
	struct make_operator_directional_derivative<ax, 0>
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver)
		{
			return v;
		}

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver)
		{
			return OpIdentity{};
		}
	};


	//! Implementation of generalized mixed derivative expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct generalized derivative expressions. Wraps the template
	 * deduction necessary to initialize a derivative expression. As opposed to
	 * the usual derivative, there is no expression associated with generalized
	 * derivatives.
	 *
	 * \tparam O Order of the generalized derivative.
	 */
	template<Axis ax, size_t O>
	struct make_operator_mixed_directional_derivative : make_operator_directional_derivative<ax, O> {};


	//! Implementation of generalized mixed derivative expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct generalized derivative expressions. Wraps the template
	 * deduction necessary to initialize a derivative expression. As opposed to
	 * the usual derivative, there is no expression associated with generalized
	 * derivatives.
	 *
	 * \tparam O Order of the generalized derivative.
	 */
	template<size_t... Os>
	struct make_operator_mixed_derivative
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver);

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver);

		// If passed an OpLiteral, uses its value rather than the object.
		template<typename V, typename E1>
		static auto get(OpLiteral<V> const& v, E1&& e1)
		{
			return get(v.value, std::forward<E1>(e1));
		}
	};

	template<>
	struct make_operator_mixed_derivative<0, 0>
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver)
		{
			return v;
		}

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver)
		{
			return OpIdentity{};
		}
	};

	template<>
	struct make_operator_mixed_derivative<0, 0, 0>
	{
		template<typename V, typename Sp>
		static auto get(V const& v, solver_op_type<Sp> solver)
		{
			return v;
		}

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver)
		{
			return OpIdentity{};
		}
	};



	//! Alias for constructing nth derivative. 
	template<Axis ax, size_t O, typename Sp>
	using nth_directional_derivative_apply = make_derivative<typename Solver<Sp>::template directional_derivative<ax, O>>;

	//! Alias for constructing nth derivative. 
	template<typename Sp, size_t... Os>
	using nth_mixed_derivative_apply = make_derivative<typename Solver<Sp>::template mixed_derivative<Os...>>;

	//! Alias for constructing nth derivative. 
	template<Axis ax, size_t O, typename Sp>
	using nth_derivative_apply = make_derivative<typename Solver<Sp>::template derivative<ax, O>>;

	////! Alias for constructing gradient.
	//template<Axis ax, typename Sp>
	//using gradient_apply = nth_directional_derivative_apply<ax, 1, Sp>;

	////! Alias for constructing laplacian.
	//template<Axis ax, typename Sp>
	//using laplacian_apply = nth_derivative_apply<ax, 2, Sp>;

	////! Alias for constructing gradlaplacian.
	//template<Axis ax, typename Sp>
	//using gradlaplacian_apply = nth_derivative_apply<ax, 3, Sp>;

	////! Alias for constructing bilaplacian.
	//template<Axis ax, typename Sp>
	//using bilaplacian_apply = nth_derivative_apply<ax, 4, Sp>;


	template<typename E>
	struct setup_result_data
	{
		E operator()(grid::dim_list const& dims)
		{
			return { dims };
		}
	};

	template<>
	struct setup_result_data<expr::symbols::Symbol>
	{
		expr::symbols::Symbol operator()(grid::dim_list const& dims)
		{
			return {};
		}
	};

}


namespace expr
{


	//! Create a derivative expression with the given term.
	/*!
	 * Create a derivative of the given expression using the given solver
	 * to numerically approximate the derivative.
	 * 
	 * \param a The expression being differentiated.
	 * \param solver The solver which approximates the derivative.
	 */
	template<typename Dd, typename A, typename B>
	auto make_derivative(A&& a, B&& solver)
	{
		return symphas::internal::make_derivative<Dd>::template get(std::forward<A>(a), std::forward<B>(solver));
	}

	//! Create a derivative expression with the given term.
	/*!
	 * Create a derivative of the given expression using the given solver
	 * to numerically approximate the derivative.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The expression being differentiated.
	 * \param solver The solver which approximates the derivative.
	 */
	template<typename Dd, typename V, typename A, typename B>
	auto make_derivative(V&& v, A&& a, B&& solver)
	{
		return symphas::internal::make_derivative<Dd>::template get(std::forward<V>(v), std::forward<A>(a), std::forward<B>(solver));
	}

	//! Create a non-directional derivative expression with the given term.
	/*!
	 * Create a derivative of the given expression using the given solver
	 * to numerically approximate the derivative. The directive will be correctly constructed
	 * for the given axis. If the order is even, then the axis is irrelevant.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The expression being differentiated.
	 * \param solver The solver which approximates the derivative.
	 */
	template<Axis ax, size_t O, typename V, typename A, typename Sp>
	auto make_derivative(V&& v, A&& a, solver_op_type<Sp> solver)
	{
		return symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(std::forward<V>(v), std::forward<A>(a), solver);
	}

	//! Create a non-directional derivative expression with the given term.
	/*!
	 * Create a non-directional derivative of the given expression using the given solver
	 * to numerically approximate the derivative. The directive will be correctly constructed
	 * for the given axis. If the order is even, then the axis is irrelevant.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The expression being differentiated.
	 * \param solver The solver which approximates the derivative.
	 */
	template<Axis ax, size_t O, typename A, typename Sp>
	auto make_derivative(A&& a, solver_op_type<Sp> solver)
	{
		return symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(std::forward<A>(a), solver);
	}

	//! Create a symbolic derivative expression with the given term, of the given order.
	/*!
	 * Create a derivative of the given expression using the given derivative order, and
	 * with respect to the given term.
	 *
	 * \param a The expression being differentiated.
	 * 
	 * \tparam O The order of the derivative.
	 * \tparam G The variable against which this expression is differentiated.
	 */
	template<size_t O, typename G, typename A>
	auto make_derivative(A&& a)
	{
		return make_derivative<std::index_sequence<O>>(std::forward<A>(a), SymbolicDerivative<G>{});
	}

	//! Create a symbolic derivative expression with the given term, of the given order.
	/*!
	 * Create a derivative of the given expression using the given derivative order, and
	 * with respect to the given term.
	 *
	 * \param a The expression being differentiated.
	 *
	 * \tparam O The order of the derivative.
	 * \tparam G The variable against which this expression is differentiated.
	 */
	template<size_t O, typename G, typename A, typename B>
	auto make_derivative(A&& a, B&& b)
	{
		return make_derivative<std::index_sequence<O>>(OpIdentity{}, std::forward<A>(a), SymbolicDerivative<G>(std::forward<B>(b)));
	}

	//! Create a symbolic derivative expression with the given term, of the given order.
	/*!
	 * Create a derivative of the given expression using the given derivative order, and
	 * with respect to the given term.
	 *
	 * \param a The expression being differentiated.
	 *
	 * \tparam O The order of the derivative.
	 * \tparam Z The index of the variable against which to differentiate.
	 */
	template<size_t O, size_t Z, typename A>
	auto make_derivative(A&& a)
	{
		return make_derivative<std::index_sequence<O>>(std::forward<A>(a), SymbolicDerivative<Variable<Z>>{});
	}

	//! Create a symbolic derivative expression with the given term, of the given order.
	/*!
	 * Create a derivative of the given expression using the given derivative order, and
	 * with respect to the given term.
	 *
	 * \param a The expression being differentiated.
	 *
	 * \tparam O The order of the derivative.
	 * \tparam Z The index of the variable against which to differentiate.
	 */
	template<size_t O, size_t Z, typename V, typename A>
	auto make_derivative(V&& v, A&& a)
	{
		return make_derivative<std::index_sequence<O>>(std::forward<V>(v), std::forward<A>(a), SymbolicDerivative<Variable<Z>>{});
	}

	template<size_t O, typename E, typename A>
	auto make_symbolic_derivative(OpExpression<E> const& e, A const& symbol)
	{
		return expr::make_derivative<O, A>(*static_cast<E const*>(&e));
	}

	template<size_t O, typename E, typename A>
	auto make_symbolic_derivative(OpOperator<E> const& e, A const& symbol)
	{
		return expr::make_derivative<O, A>(*static_cast<E const*>(&e));
	}

	template<size_t O, typename E, size_t Z, typename G>
	auto make_symbolic_derivative(OpExpression<E> const& e, OpTerm<OpIdentity, Variable<Z, G>> const& symbol)
	{
		return expr::make_derivative<O, Z>(*static_cast<E const*>(&e));
	}

	template<size_t O, typename E, size_t Z, typename G>
	auto make_symbolic_derivative(OpOperator<E> const& e, OpTerm<OpIdentity, Variable<Z, G>> const& symbol)
	{
		return expr::make_derivative<O, Z>(*static_cast<E const*>(&e));
	}

	template<size_t Z, typename V, typename E>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e)
	{
		return OpFunctionalDerivative<V, E, Variable<Z>>(value, *static_cast<E const*>(&e));
	}

	template<size_t Z, typename V, typename E>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e)
	{
		return OpFunctionalDerivative<V, E, Variable<Z>>(value, *static_cast<E const*>(&e));
	}

	template<size_t Z, typename E>
	auto make_functional_derivative(E&& e)
	{
		return make_functional_derivative<Z>(OpIdentity{}, std::forward<E>(e));
	}

	template<typename V, typename E, size_t Z, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, Variable<Z, G>)
	{
		return OpFunctionalDerivative<V, E, Variable<Z, G>>(value, *static_cast<E const*>(&e));
	}

	template<typename V, typename E, size_t Z, typename G>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e, Variable<Z, G>)
	{
		return OpFunctionalDerivative<V, E, Variable<Z, G>>(value, *static_cast<E const*>(&e));
	}

	template<typename E, size_t Z, typename G>
	auto make_functional_derivative(OpExpression<E> const& e, Variable<Z, G>)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), Variable<Z, G>{});
	}

	template<typename E, size_t Z, typename G>
	auto make_functional_derivative(OpOperator<E> const& e, Variable<Z, G>)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), Variable<Z, G>{});
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, DynamicVariable<G> const& symbol)
	{
		return OpFunctionalDerivative<V, E, DynamicVariable<G>>(value, *static_cast<E const*>(&e), symbol);
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e, DynamicVariable<G> const& symbol)
	{
		return OpFunctionalDerivative<V, E, DynamicVariable<G>>(value, *static_cast<E const*>(&e), symbol);
	}

	template<typename E, typename G>
	auto make_functional_derivative(OpExpression<E> const& e, DynamicVariable<G> const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), symbol);
	}

	template<typename E, typename G>
	auto make_functional_derivative(OpOperator<E> const& e, DynamicVariable<G> const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), symbol);
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, G const& symbol)
	{
		return OpFunctionalDerivative<V, E, G>(value, *static_cast<E const*>(&e));
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e, G const& symbol)
	{
		return OpFunctionalDerivative<V, E, G>(value, *static_cast<E const*>(&e));
	}

	template<typename E, typename G>
	auto make_functional_derivative(OpExpression<E> const& e, G const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), G{});
	}

	template<typename E, typename G>
	auto make_functional_derivative(OpOperator<E> const& e, G const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, *static_cast<E const*>(&e), G{});
	}

	template<typename V, typename E, size_t Z, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, OpTerm<OpIdentity, Variable<Z, G>> const& symbol)
	{
		return OpFunctionalDerivative<V, E, Variable<Z, G>>(value, *static_cast<E const*>(&e), SymbolicDerivative(expr::get<1>(symbol).data()));
	}

	template<typename E, size_t Z, typename G>
	auto make_functional_derivative(E&& e, OpTerm<OpIdentity, Variable<Z, G>> const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, std::forward<E>(e), symbol);
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, OpTerm<OpIdentity, DynamicVariable<G>> const& symbol)
	{
		return OpFunctionalDerivative<V, E, DynamicVariable<G>>(value, *static_cast<E const*>(&e), SymbolicDerivative(expr::get<1>(symbol).data()));
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e, OpTerm<OpIdentity, DynamicVariable<G>> const& symbol)
	{
		return OpFunctionalDerivative<V, E, DynamicVariable<G>>(value, *static_cast<E const*>(&e), SymbolicDerivative(expr::get<1>(symbol).data()));
	}

	template<typename E, typename G>
	auto make_functional_derivative(E&& e, OpTerm<OpIdentity, DynamicVariable<G>> const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, std::forward<E>(e), symbol);
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpExpression<E> const& e, OpTerm<OpIdentity, G> const& symbol)
	{
		return OpFunctionalDerivative<V, E, G>(value, *static_cast<E const*>(&e));
	}

	template<typename V, typename E, typename G>
	auto make_functional_derivative(V const& value, OpOperator<E> const& e, OpTerm<OpIdentity, G> const& symbol)
	{
		return OpFunctionalDerivative<V, E, G>(value, *static_cast<E const*>(&e));
	}

	template<typename E, typename G>
	auto make_functional_derivative(E&& e, OpTerm<OpIdentity, G> const& symbol)
	{
		return make_functional_derivative(OpIdentity{}, std::forward<E>(e), symbol);
	}


	//! Create a derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 * 
	 * \tparam O The order of the derivative to generate.
	 */
	template<size_t O, typename A>
	auto make_operator_derivative(A&& solver)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<A>(solver));
	}

	//! Create a derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The solver which approximates the derivative.
	 *
	 * \tparam O The order of the derivative to generate.
	 */
	template<size_t O, typename V, typename A>
	auto make_operator_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<V>(v), std::forward<A>(solver));
	}

	//! Create a symbolic derivative operator.
	/*!
	 * Create a derivative of the given order for symbolic differentiation.
	 *
	 * \tparam O The order of the derivative to generate.
	 * \tparam G The variable against which this expression is differentiated.
	 */
	template<size_t O, typename G>
	auto make_operator_derivative()
	{
		return symphas::internal::make_operator_derivative<O>::template get(SymbolicDerivative<G>{});
	}

	//! Create a symbolic derivative operator.
	/*!
	 * Create a derivative of the given order for symbolic differentiation.
	 *
	 * \param v The coefficient to the derivative.
	 *
	 * \tparam O The order of the derivative to generate.
	 * \tparam G The variable against which this expression is differentiated.
	 */
	template<size_t O, typename G, typename V>
	auto make_operator_derivative(V&& v)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<V>(v), SymbolicDerivative<G>{});
	}

	//! Create a symbolic derivative operator.
	/*!
	 * Create a derivative of the given order for symbolic differentiation.
	 *
	 * \tparam O The order of the derivative to generate.
	 * \tparam Z The index of the variable against which to differentiate.
	 */
	template<size_t O, size_t Z>
	auto make_operator_derivative()
	{
		return symphas::internal::make_operator_derivative<O>::template get(SymbolicDerivative<Variable<Z>>{});
	}

	//! Create a symbolic derivative operator.
	/*!
	 * Create a derivative of the given order for symbolic differentiation.
	 *
	 * \param v The coefficient to the derivative.
	 *
	 * \tparam O The order of the derivative to generate.
	 * \tparam Z The index of the variable against which to differentiate.
	 */
	template<size_t O, size_t Z, typename V>
	auto make_operator_derivative(V&& v)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<V>(v), SymbolicDerivative<Variable<Z>>{});
	}

	//! Create a directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax, size_t O, typename A>
	auto make_operator_directional_derivative(A&& solver)
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(std::forward<A>(solver));
	}

	//! Create a directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax, size_t O, typename V, typename A>
	auto make_operator_directional_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(std::forward<V>(v), std::forward<A>(solver));
	}


	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, Axis ax2, size_t O2, Axis ax3, size_t O3, typename A,
		typename std::enable_if_t<(ax1 != ax2 && ax2 != ax3 && ax1 != ax3), int> = 0>
	auto make_operator_mixed_derivative(A&& solver)
	{
		return symphas::internal::make_operator_mixed_derivative<
			(ax1 == Axis::X) ? O1 : (ax2 == Axis::X) ? O2 : (ax3 == Axis::X) ? O3 : 0,
			(ax1 == Axis::Y) ? O1 : (ax2 == Axis::Y) ? O2 : (ax3 == Axis::Y) ? O3 : 0,
			(ax1 == Axis::Z) ? O1 : (ax2 == Axis::Z) ? O2 : (ax3 == Axis::Z) ? O3 : 0>
			::template get(std::forward<A>(solver));
	}

	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, Axis ax2, size_t O2, typename A,
		typename std::enable_if_t<(ax1 != ax2), int> = 0>
	auto make_operator_mixed_derivative(A&& solver)
	{
		return symphas::internal::make_operator_mixed_derivative<
			(ax1 == Axis::X) ? O1 : (ax2 == Axis::X) ? O2 : 0,
			(ax1 == Axis::Y) ? O1 : (ax2 == Axis::Y) ? O2 : 0,
			(ax1 == Axis::Z) ? O1 : (ax2 == Axis::Z) ? O2 : 0>
			::template get(std::forward<A>(solver));
	}

	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, typename A>
	auto make_operator_mixed_derivative(A&& solver)
	{
		return make_operator_directional_derivative<ax1, O1>(std::forward<A>(solver));
	}
	
	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, Axis ax2, size_t O2, Axis ax3, size_t O3, typename V, typename A,
		typename std::enable_if_t<(ax1 != ax2 && ax2 != ax3 && ax1 != ax3), int> = 0>
	auto make_operator_mixed_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_mixed_derivative<
			(ax1 == Axis::X) ? O1 : (ax2 == Axis::X) ? O2 : (ax3 == Axis::X) ? O3 : 0,
			(ax1 == Axis::Y) ? O1 : (ax2 == Axis::Y) ? O2 : (ax3 == Axis::Y) ? O3 : 0,
			(ax1 == Axis::Z) ? O1 : (ax2 == Axis::Z) ? O2 : (ax3 == Axis::Z) ? O3 : 0>
			::template get(std::forward<V>(v), std::forward<A>(solver));
	}

	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, Axis ax2, size_t O2, typename V, typename A,
		typename std::enable_if_t<(ax1 != ax2), int> = 0>
	auto make_operator_mixed_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_mixed_derivative<
			(ax1 == Axis::X) ? O1 : (ax2 == Axis::X) ? O2 : 0,
			(ax1 == Axis::Y) ? O1 : (ax2 == Axis::Y) ? O2 : 0,
			(ax1 == Axis::Z) ? O1 : (ax2 == Axis::Z) ? O2 : 0>
			::template get(std::forward<V>(v), std::forward<A>(solver));
	}

	//! Create a mixed directional derivative operator.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate.
	 */
	template<Axis ax1, size_t O1, typename V, typename A>
	auto make_operator_mixed_derivative(V&& v, A&& solver)
	{
		return make_operator_directional_derivative<ax1, O1>(std::forward<V>(v), std::forward<A>(solver));
	}


	template<size_t O, typename Sp, Axis... axs>
	auto break_up_derivative(solver_op_type<Sp> solver, symphas::lib::axis_list<axs...>)
	{
		if constexpr (O == 0)
		{
			return OpIdentity{};
		}
		else if constexpr (O == 1)
		{
			return ((expr::make_column_vector<static_cast<size_t>(axs), sizeof...(axs)>() 
				* make_operator_directional_derivative<axs, 1>(solver)) + ...);
		}
		else if constexpr (O == 2)
		{
			return (make_operator_directional_derivative<axs, 2>(solver) + ...);
		}
		else
		{
			return break_up_derivative<2>(solver, symphas::lib::axis_list<axs...>{}) * break_up_derivative<O - 2>(solver, symphas::lib::axis_list<axs...>{});
		}
	}

	template<size_t O, size_t D, typename Sp>
	auto break_up_derivative(solver_op_type<Sp> solver)
	{
		return break_up_derivative<O>(solver, symphas::lib::make_axis_list<D>());
	}
}

// *************************************************************************************


//! Concrete derivative expression.
/*!
 * If the derivative operator is applied to something that's not a variable or 
 * a linear variable, it needs to evaluate everything inside.
 * This object uses a grid to update the values on a prune call.
 *
 * The implementation of the derivative classes are typically done with respect 
 * to an estimation using finite differences, where the order is chosen 
 * statically and by the given differentiation object type, `Dd`.
 * 
 * Additionally, derivatives require that the call to `eval' does not call 
 * an index which would be considered on the "boundary", since this methodology 
 * in general will apply finite difference approximations directly by the
 * stencil, which visit direct memory neighbours.
 * 
 * \tparam Dd An object with an operator that can apply the derivative.
 * \tparam V Type of the coefficient.
 * \tparam E Expression type that the derivative applies to.
 * \tparam Sp The solver type.
 */
template<typename Dd, typename V, typename E, typename Sp>
struct OpDerivative : OpExpression<OpDerivative<Dd, V, E, Sp>>
{
	using result_grid = expr::storage_t<E>;

	static const size_t order = Dd::order;		//!< The order of this derivative.
	static const Axis axis = Dd::axis;			//!< The axis of this derivative.
	static const bool is_directional = Dd::is_directional;	//!< Whether the derivative is directional.


	OpDerivative() : grid{ 0 }, value{ V{} }, solver{}, e{} {}

	//! Generate a derivative of the given expression.
	/*!
	 * Create a derivative expression representing applying a derivative
	 * to the given expression, where the order depends on the differentiation
	 * object.
	 */
	OpDerivative(V value, E const& e, solver_op_type<Sp> solver) : 
		grid{ symphas::internal::setup_result_data<result_grid>{}(expr::data_dimensions(e)) }, 
		value{value}, solver{solver}, e{e} { /*update();*/ }

	OpDerivative(E const& e, solver_op_type<Sp> solver) : 
		grid{ expr::data_dimensions(e) }, value{ OpIdentity{} }, solver{ solver }, e{ e } { /*update();*/ }

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * Dd{}(solver, grid, n);
	}

	void update()
	{
		expr::prune::update(e);
		expr::result(e, grid);
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<Dd>::get(-value, e, solver);
	}

#ifdef PRINTABLE_EQUATIONS
	using print_type = decltype(symphas::internal::select_print_deriv<Sp>(Dd{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length() + e.print_length();
	}

#endif
	
    template<typename Dd0, typename V0, typename E0, typename Sp0>
	friend auto const& expr::get_enclosed_expression(OpDerivative<Dd0, V0, E0, Sp0> const&);
    template<typename Dd0, typename V0, typename E0, typename Sp0>
	friend auto& expr::get_enclosed_expression(OpDerivative<Dd0, V0, E0, Sp0>&);
    template<typename Dd0, typename V0, typename E0, typename Sp0>
	friend auto const& expr::get_result_data(OpDerivative<Dd0, V0, E0, Sp0> const&);
    template<typename Dd0, typename V0, typename E0, typename Sp0>
	friend auto& expr::get_result_data(OpDerivative<Dd0, V0, E0, Sp0>&);


protected:

	result_grid grid;				//!< Grid storing the intermediate values.

public:

	V value;						//!< Value multiplying the result of this derivative.
	Sp solver;						//!< Solver that applies the derivative function.

protected:

	E e;							//!< Expression object specifying grid values.


};

//! Specialization of an OpTerm of the concrete derivative expression.
/*!
 * \tparam Dd An object with an operator that can apply the derivative.
 * \tparam V Type of the coefficient.
 * \tparam G Grid type of the linear variable.
 * \tparam Sp The solver type.
 */
template<typename Dd, typename V, typename G, typename Sp>
struct OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> : OpExpression<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>>, Dd
{
	using E = OpTerm<OpIdentity, G>;
	using result_grid = expr::storage_t<E>;

	using Dd::order;
	using Dd::axis;
	using Dd::is_directional;

	template<typename V0, typename V1, typename std::enable_if_t<std::is_convertible<mul_result_t<V0, V1>, V>::value, int> = 0>
	OpDerivative(V0 value, OpTerm<V1, G> const& e, solver_op_type<Sp> solver) : value{ expr::coeff(e) * value }, solver{ solver }, data{ expr::data(e) } {}
	OpDerivative(V value, G data, solver_op_type<Sp> solver) : value{ value }, solver{ solver }, data{ data } {}

	OpDerivative() : value{ V{} }, solver{ Sp::make_solver() }, data{} {}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * Dd{}(solver, expr::BaseData<G>::get(data), n);
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<Dd>::get_g(-value, data, solver);
	}

	void update()
	{
		expr::prune::update(data);
	}


#ifdef PRINTABLE_EQUATIONS
	using print_type = decltype(symphas::internal::select_print_deriv<Sp>(Dd{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out, expr::get_op_name(data));
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n, expr::get_op_name(data));
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length() + std::strlen(expr::get_op_name(data));
	}

#endif

    template<typename Dd0, typename V0, typename G0, typename Sp0>
	friend decltype(auto) expr::get_enclosed_expression(OpDerivative<Dd0, V0, OpTerm<OpIdentity, G0>, Sp0> const&);
    template<typename Dd0, typename V0, typename G0, typename Sp0>
	friend decltype(auto) expr::get_enclosed_expression(OpDerivative<Dd0, V0, OpTerm<OpIdentity, G0>, Sp0>&);

public:

	V value;						//!< Value multiplying the result of this derivative.
	Sp solver;						//!< Solver that applies the derivative function.

protected:

	G data;


};



//! Symbolic derivative expression.
/*!
 * Represents the symbolic differentiation with respect to a variable defined by the
 * template parameter `G`. This is a specialization of the base derivative object which applies
 * a concrete solver to an expression. This object cannot evaluate any derivatives, and is only
 * used to represent differentiation in order to perform symbolic differentiation of terms.
 * 
 * \tparam O The order of the derivative.
 * \tparam G Defines the variable against which the expression is differentiated. For example,
 * this can be an std::index_sequence type.
 * \tparam V Type of the coefficient.
 * \tparam E Expression type that the derivative applies to.
 */
template<size_t O, typename G, typename V, typename E>
struct OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>> : 
	OpExpression<OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>>
{
	static const size_t order = O;			//!< The order of this derivative.

	OpDerivative() : value{ V{} }, e{}, solver{} {}

	//! Generate a derivative of the given expression.
	/*!
	 * Create a derivative expression representing applying a derivative
	 * to the given expression, where the order depends on the differentiation
	 * object.
	 */
	OpDerivative(V value, E const& e, SymbolicDerivative<G> const& solver) :
		value{ value }, solver{ solver }, e{ e } {}

	OpDerivative(E const& e, SymbolicDerivative<G> const& solver) :
		value{ OpIdentity{} }, solver{ solver }, e{ e } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<std::index_sequence<O>>::get(-value, e, solver);
	}

#ifdef PRINTABLE_EQUATIONS
	using print_type = decltype(symphas::internal::select_print_deriv(std::index_sequence<O>{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out, G{});
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n, G{});
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length(G{}) + e.print_length();
	}

#endif

    template<size_t O0, typename V0, typename E0, typename G0>
	friend auto const& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>> const&);
    template<size_t O0, typename V0, typename E0, typename G0>
	friend auto& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>>&);

public:

	V value;							//!< Value multiplying the result of this derivative.
	SymbolicDerivative<G> solver;		//!< Ensures this object is valid in the symbolic algebra.

protected:

	E e;								//!< Expression object specifying grid values.

};


//! Specialization of the symbolic derivative.
/*!
 * 
 * \tparam O The order of the derivative.
 * \tparam G Defines the variable against which the expression is differentiated. For example,
 * this can be an std::index_sequence type.
 * \tparam V Type of the coefficient.
 * \tparam E Expression type that the derivative applies to.
 */
template<size_t O, typename G1, typename G2, typename V>
struct OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G1>, SymbolicDerivative<G2>> : 
	OpExpression<OpDerivative<std::index_sequence<O>, V, OpTerm<OpIdentity, G1>, SymbolicDerivative<G2>>>
{
	using E = OpTerm<OpIdentity, G1>;

	static const size_t order = O; //!< The order of this derivative.

	OpDerivative(V value, OpTerm<OpIdentity, G1> e, SymbolicDerivative<G2> const& solver) : value{ value }, solver{ solver }, e{ e } {}
	OpDerivative() : value{ V{} }, e{}, solver{} {}

	template<typename V0, typename V1, typename std::enable_if_t<std::is_convertible<mul_result_t<V0, V1>, V>::value, int> = 0>
	OpDerivative(V0 value, OpTerm<V1, G1> const& e, SymbolicDerivative<G2> const& solver) : value{ e.term * value }, solver{ solver }, e{ e } {}
	
	OpDerivative(V value, G1 data, SymbolicDerivative<G2> const& solver) : value{ value }, solver{ solver }, e{ OpTerm<OpIdentity, G1>(OpIdentity{}, Term(data)) } {}
	
	template<typename V0, typename V1, typename std::enable_if_t<std::is_convertible<mul_result_t<V0, V1>, V>::value, int> = 0>
	OpDerivative(V0 value, OpTerm<V1, G1> const& e) : value{ e.value * value }, solver{}, e{ e } {}
	
	OpDerivative(V value, G1 data) : value{ value }, solver{}, e{ OpTerm<OpIdentity, G1>(OpIdentity{}, Term(data)) } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<std::index_sequence<O>>::get(-value, e, solver);
	}


#ifdef PRINTABLE_EQUATIONS
	using print_type = decltype(symphas::internal::select_print_deriv(std::index_sequence<O>{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out, expr::get_op_name(e.data), G2{});
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n, expr::get_op_name(e.data), G2{});
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length(G2{}) + std::strlen(expr::get_op_name(e.data));
	}

#endif

    template<size_t O0, typename V0, typename G10, typename G20>
	friend auto const& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, OpTerm<OpIdentity, G10>, SymbolicDerivative<G20>> const&);
    template<size_t O0, typename V0, typename G10, typename G20>
	friend auto& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, OpTerm<OpIdentity, G10>, SymbolicDerivative<G20>>&);

public:

	V value;						//!< Value multiplying the result of this derivative.
	SymbolicDerivative<G2> solver;	//!< Solver that applies the derivative function.

protected:

	E e;							//!< Expression object specifying grid values.

};

//! Symbolic derivative expression.
/*!
 * Represents the symbolic differentiation with respect to a variable defined by the
 * template parameter `G`. This is a specialization of the base derivative object which applies
 * a concrete solver to an expression. This object cannot evaluate any derivatives, and is only
 * used to represent differentiation in order to perform symbolic differentiation of terms.
 *
 * \tparam O The order of the derivative.
 * \tparam G Defines the variable against which the expression is differentiated. For example,
 * this can be an std::index_sequence type.
 * \tparam V Type of the coefficient.
 * \tparam E Expression type that the derivative applies to.
 */
template<typename G, typename V, typename E>
struct OpDerivative<std::index_sequence<1>, V, E, SymbolicFunctionalDerivative<G>> :
	OpExpression<OpDerivative<std::index_sequence<1>, V, E, SymbolicFunctionalDerivative<G>>>
{
	static const size_t order = 1;			//!< The order of this derivative.

	using solver_t = std::invoke_result_t<decltype(expr::get_solver<E>), E>;

	OpDerivative() : value{ V{} }, solver{}, implicit_solver{}, e{} {}

	//! Generate a derivative of the given expression.
	/*!
	 * Create a derivative expression representing applying a derivative
	 * to the given expression, where the order depends on the differentiation
	 * object.
	 */
	OpDerivative(V value, E const& e, SymbolicFunctionalDerivative<G> const& solver = SymbolicFunctionalDerivative<G>{}) :
		value{ value }, solver{ solver }, implicit_solver{ expr::get_solver(e) }, e{ e } {}

	OpDerivative(E const& e, SymbolicFunctionalDerivative<G> const& solver = SymbolicFunctionalDerivative<G>{}) :
		value{ OpIdentity{} }, solver{ solver }, implicit_solver{ expr::get_solver(e) }, e{ e } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<std::index_sequence<1>>::get(-value, e, solver);
	}


#ifdef PRINTABLE_EQUATIONS

	using print_type = decltype(symphas::internal::select_print_deriv(SymbolicFunctionalDerivative<G>{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out, solver);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n, solver);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print_length() const
	{

		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length(solver) + e.print_length();
	}

#endif

	template<size_t O0, typename V0, typename E0, typename G0>
	friend auto const& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>> const&);
	template<size_t O0, typename V0, typename E0, typename G0>
	friend auto& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>>&);

public:

	V value;									//!< Value multiplying the result of this derivative.
	SymbolicFunctionalDerivative<G> solver;		//!< Ensures this object is valid in the symbolic algebra.
	solver_t implicit_solver;					//!< The implicit solver, derived from the given expression.

protected:

	E e;										//!< Expression object specifying grid values.

};

//! Symbolic derivative expression.
/*!
 * Represents the symbolic differentiation with respect to a variable defined by the
 * template parameter `G`. This is a specialization of the base derivative object which applies
 * a concrete solver to an expression. This object cannot evaluate any derivatives, and is only
 * used to represent differentiation in order to perform symbolic differentiation of terms.
 *
 * \tparam O The order of the derivative.
 * \tparam G Defines the variable against which the expression is differentiated. For example,
 * this can be an std::index_sequence type.
 * \tparam V Type of the coefficient.
 * \tparam E Expression type that the derivative applies to.
 */
template<typename V, typename G1, typename G2>
struct OpDerivative<std::index_sequence<1>, V, OpTerm<OpIdentity, G1>, SymbolicFunctionalDerivative<G2>> :
	OpExpression<OpDerivative<std::index_sequence<1>, V, OpTerm<OpIdentity, G1>, SymbolicFunctionalDerivative<G2>>>
{
	using E = OpTerm<OpIdentity, G1>;
	static const size_t order = 1;			//!< The order of this derivative.

	using solver_t = std::invoke_result_t<decltype(expr::get_solver<E>), E>;

	OpDerivative() : value{ V{} }, solver{}, implicit_solver{}, e{} {}

	//! Generate a derivative of the given expression.
	/*!
	 * Create a derivative expression representing applying a derivative
	 * to the given expression, where the order depends on the differentiation
	 * object.
	 */
	OpDerivative(V value, E const& e, SymbolicFunctionalDerivative<G2> const& solver = SymbolicFunctionalDerivative<G2>{}) :
		value{ value }, solver{ solver }, implicit_solver{ expr::get_solver(e) }, e{ e } {}

	OpDerivative(E const& e, SymbolicFunctionalDerivative<G2> const& solver = SymbolicFunctionalDerivative<G2>{}) :
		value{ OpIdentity{} }, solver{ solver }, implicit_solver{ expr::get_solver(e) }, e{ e } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<std::index_sequence<1>>::get(-value, e, solver);
	}

#ifdef PRINTABLE_EQUATIONS
	using print_type = decltype(symphas::internal::select_print_deriv(expr::variational_t<G2>{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ print_type::print_length() + e.print_length();
	}

#endif

	template<size_t O0, typename V0, typename E0, typename G0>
	friend auto const& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>> const&);
	template<size_t O0, typename V0, typename E0, typename G0>
	friend auto& expr::get_enclosed_expression(OpDerivative<std::index_sequence<O0>, V0, E0, SymbolicDerivative<G0>>&);

public:

	V value;									//!< Value multiplying the result of this derivative.
	SymbolicFunctionalDerivative<G2> solver;	//!< Ensures this object is valid in the symbolic algebra.
	solver_t implicit_solver;					//!< The implicit solver, derived from the given expression.

protected:

	E e;										//!< Expression object specifying grid values.

};


/* multiplication of a derivative object by a literal
 */
template<typename coeff_t, typename Dd2, typename V2, typename E2, typename Sp,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpDerivative<Dd2, V2, E2, Sp> const& b)
{
	return symphas::internal::make_derivative<Dd2>::template get(value * b.value, expr::get_enclosed_expression(b), b.solver);
}

template<typename coeff_t, typename Dd2, typename tensor_t, typename E2, typename Sp,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpDerivative<Dd2, tensor_t, E2, Sp> const& b)
{
	return (value * b.value) * symphas::internal::make_derivative<Dd2>::template get(expr::get_enclosed_expression(b), b.solver);
}



//
///* multiplication of a derivative object that owns only a grid by a literal
// */
//template<typename S1, typename Dd2, typename V2, typename G2, typename T2>
//auto operator*(OpLiteral<S1> const& a, OpDerivative<Dd2, V2, OpTerm<OpIdentity, G2>, T2> const& b)
//{
//	return symphas::internal::make_derivative<Dd2>::template get_g(a.value * b.value, b.data, b.solver);
//}
//
//template<typename fraction_t, typename coeff_t, typename Dd2, typename G2, typename Sp2, 
//	typename = std::enable_if_t<(expr::is_fraction<fraction_t> && (expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>)), int>>
//auto operator*(fraction_t, OpDerivative<Dd2, coeff_t, OpTerm<OpIdentity, G2>, Sp2> const& b)
//{
//	return symphas::internal::make_derivative<Dd2>::template get_g(fraction_t{} * coeff_t{}, b.data, b.solver);
//}
//
//template<typename fraction_t, typename coeff_t, typename Dd2, typename G2, typename G0,
//	typename = std::enable_if_t<(expr::is_fraction<fraction_t> && (expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>)), int>>
//auto operator*(fraction_t, OpDerivative<Dd2, coeff_t, OpTerm<OpIdentity, G2>, SymbolicDerivative<G0>> const& b)
//{
//	auto&& e = expr::get_enclosed_expression(b);
//	return symphas::internal::make_derivative<Dd2>::template get(fraction_t{} * coeff_t{}, e, b.solver);
//}


// ******************************************************************************************




// ******************************************************************************************



//! Generalized derivative object.
/*! 
 * A specialized operation that will apply a derivative.
 * Acts simply as an operator; it is not applied to any expressions, but can
 * be distributed to expressions with the multiplication operation.
 */
template<size_t O, typename V, typename Sp>
struct OpOperatorDerivative : OpOperator<OpOperatorDerivative<O, V, Sp>>
{
	using parent_type = OpOperator<OpOperatorDerivative<O, V, Sp>>;
	using parent_type::operator*;
	using parent_type::operator-;
	using parent_type::operator+;

	OpOperatorDerivative() : value{ V{} }, solver{} {}

	static_assert(O <= DERIV_MAX_ORDER);
	OpOperatorDerivative(V value, solver_op_type<Sp> solver) : value{ value }, solver{ solver } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_operator_derivative<O>::template get(
			-value, solver);
	}

	//! Apply to an expression.
	/*!
	 * Apply the generalized derivative to an expression to produce a
	 * concrete derivative.
	 */
	template<typename E>
	auto apply_impl(E const& e) const
	{
		return _apply_impl(e);
	}

	template<typename E>
	auto operator*(OpExpression<E> const& e) const;

	//template<size_t O1, typename V1>
	//auto operator*(OpOperatorDerivative<O1, V1, Sp> const& other)
	//{
	//	return expr::coeff(other) * value * expr::make_operator_derivative<O + O1>(solver);
	//}

	//template<typename E>
	//auto operator*(OpOperator<E> const& e)
	//{
	//	return OpOperator<OpOperatorDerivative<O, V, Sp>>::operator*(*static_cast<E const*>(&e));
	//}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out, solver);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out + n, solver);
		return n;
	}


	size_t print_length() const
	{
		return expr::coeff_print_length(value) + symphas::internal::print_deriv<O>::print_length();
	}

#endif

	V value;
	Sp solver;

protected:

	template<typename E>
	auto _apply_impl(OpExpression<E> const& e) const;


	template<typename E>
	auto _apply_impl(OpOperator<E> const& e) const
	{
		return OpOperatorChain(*this, *static_cast<E const*>(&e));
	}
};


template<size_t O, typename V1, typename V2, typename Sp>
auto operator+(OpOperatorDerivative<O, V1, Sp> const& a, OpOperatorDerivative<O, V2, Sp> const& b)
{
	return symphas::internal::make_operator_derivative<O>::template get(a.value + b.value, a.solver);
}

template<size_t O, typename V1, typename V2, typename Sp>
auto operator-(OpOperatorDerivative<O, V1, Sp> const& a, OpOperatorDerivative<O, V2, Sp> const& b)
{
	return a + (-b);
}

template<size_t O1, size_t O2, typename V1, typename V2, typename Sp>
auto operator*(OpOperatorDerivative<O1, V1, Sp> const& a, OpOperatorDerivative<O2, V2, Sp> const& b)
{
	return symphas::internal::make_operator_derivative<O1 + O2>::template get(a.value * b.value, a.solver);
}


template<typename coeff_t, size_t O2, typename V2, typename Sp2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorDerivative<O2, V2, Sp2> const& b)
{
	return symphas::internal::make_operator_derivative<O2>::template get(value * b.value, b.solver);
}

template<typename coeff_t, typename tensor_t, size_t O2, typename Sp2,
	typename std::enable_if_t<(expr::is_tensor<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorDerivative<O2, tensor_t, Sp2> const& b)
{
	return (value * b.value) * symphas::internal::make_operator_derivative<O2>::template get(b.solver);
}




//! Generalized derivative object.
/*!
 * A specialized operation that will apply a directional derivative.
 * Acts simply as an operator; it is not applied to any expressions, but can
 * be distributed to expressions with the multiplication operation.
 */
template<Axis ax, size_t O, typename V, typename Sp>
struct OpOperatorDirectionalDerivative : OpOperator<OpOperatorDirectionalDerivative<ax, O, V, Sp>>
{
	using parent_type = OpOperator<OpOperatorDirectionalDerivative<ax, O, V, Sp>>;
	using parent_type::operator*;
	using parent_type::operator-;
	using parent_type::operator+;

	OpOperatorDirectionalDerivative() : value{ V{} }, solver{} {}

	static_assert(O <= DERIV_MAX_ORDER);
	OpOperatorDirectionalDerivative(V value, solver_op_type<Sp> solver) : value{ value }, solver{ solver } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(
			-value, solver);
	}

	template<typename E>
	auto apply_impl(E const& e) const
	{
		return _apply_impl(e);
	}


	template<typename E>
	auto operator*(OpExpression<E> const& e) const;


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O, ax, true>::print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O, ax, true>::print(out + n);
		return n;
	}


	size_t print_length() const
	{
		return expr::coeff_print_length(value) + symphas::internal::print_deriv<O, ax, true>::print_length();
	}

#endif

	V value;
	Sp solver;


	//! Apply to an expression.
	/*!
	 * Apply the generalized derivative to an expression to produce a
	 * concrete derivative.
	 */
	template<typename E, typename std::enable_if_t<!(expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0), int> = 0>
	auto _apply_impl(OpExpression<E> const& e) const;

	//! Apply to a Symbol.
	/*!
	 * Apply the generalized derivative to an expression to produce a Chain operator.
	 */
	template<typename E, typename std::enable_if_t<(expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0), int> = 0>
	auto _apply_impl(OpExpression<E> const& e) const
	{
		return symphas::internal::nth_directional_derivative_apply<ax, O, Sp>::template get(value, *static_cast<E const*>(&e), solver);
	}

	template<typename E>
	auto _apply_impl(OpOperator<E> const& e) const
	{
		return OpOperatorChain(*this, *static_cast<E const*>(&e));
	}

};


template<Axis ax, size_t O, typename V1, typename V2, typename Sp>
auto operator+(OpOperatorDirectionalDerivative<ax, O, V1, Sp> const& a, 
	OpOperatorDirectionalDerivative<ax, O, V2, Sp> const& b)
{
	return symphas::internal::make_operator_directional_derivative<ax, O>::template get(a.value + b.value, a.solver);
}

template<Axis ax, size_t O, typename V1, typename V2, typename Sp>
auto operator-(OpOperatorDirectionalDerivative<ax, O, V1, Sp> const& a, 
	OpOperatorDirectionalDerivative<ax, O, V2, Sp> const& b)
{
	return a + (-b);
}

template<Axis ax, size_t O1, size_t O2, typename V1, typename V2, typename Sp>
auto operator*(OpOperatorDirectionalDerivative<ax, O1, V1, Sp> const& a,
	OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& b)
{
	return symphas::internal::make_operator_directional_derivative<ax, O1 + O2>::template get(a.value * b.value, a.solver);
}


template<typename coeff_t, typename V, Axis ax2, size_t O2, typename Sp2, 
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorDirectionalDerivative<ax2, O2, V, Sp2> const& b)
{
	return symphas::internal::make_operator_directional_derivative<ax2, O2>::template get(value * expr::coeff(b), b.solver);
}

template<typename coeff_t, typename tensor_t, Axis ax2, size_t O2, typename Sp2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorDirectionalDerivative<ax2, O2, tensor_t, Sp2> const& b)
{
	return (value * expr::coeff(b)) * symphas::internal::make_operator_directional_derivative<ax2, O2>::template get(OpIdentity{}, b.solver);
}





//! Generalized derivative object.
/*!
 * A specialized operation that will apply a directional derivative.
 * Acts simply as an operator; it is not applied to any expressions, but can
 * be distributed to expressions with the multiplication operation.
 */
template<typename V, typename Sp, size_t... Os>
struct OpOperatorMixedDerivative : OpOperator<OpOperatorMixedDerivative<V, Sp, Os...>>
{
	using parent_type = OpOperator<OpOperatorMixedDerivative<V, Sp, Os...>>;
	using parent_type::operator*;
	using parent_type::operator-;
	using parent_type::operator+;

	OpOperatorMixedDerivative() : value{ V{} }, solver{} {}

	OpOperatorMixedDerivative(V value, solver_op_type<Sp> solver) : value{ value }, solver{ solver } {}

	inline auto eval(iter_type n = 0) const
	{
		return value * expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_operator_mixed_derivative<Os...>::template get(
			-value, solver);
	}

	template<typename E>
	auto apply_impl(E const& e) const
	{
		return _apply_impl(e);
	}

	template<typename E>
	auto operator*(OpExpression<E> const& e) const;

#ifdef PRINTABLE_EQUATIONS
    using print_type = decltype(symphas::internal::select_print_deriv<Sp>(typename Sp::template mixed_derivative<Os...>{}));

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += print_type::print(out + n);
		return n;
	}


	size_t print_length() const
	{
		return expr::coeff_print_length(value) + print_type::print_length();
	}

#endif

	V value;
	Sp solver;

protected:

	//! Apply to an expression.
	/*!
	 * Apply the generalized derivative to an expression to produce a
	 * concrete derivative.
	 */
	template<typename E, typename std::enable_if_t<!(expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0), int> = 0>
	auto _apply_impl(OpExpression<E> const& e) const;

	//! Apply to a Symbol.
	/*!
	 * Apply the generalized derivative to an expression to produce a Chain operator.
	 */
	template<typename E, typename std::enable_if_t<(expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0), int> = 0>
	auto _apply_impl(OpExpression<E> const& e) const
	{
		return symphas::internal::nth_mixed_derivative_apply<Sp, Os...>::template get(value, *static_cast<E const*>(&e), solver);
	}


	template<typename E>
	auto _apply_impl(OpOperator<E> const& e) const
	{
		return OpOperatorChain(*this, *static_cast<E const*>(&e));
	}

};




template<typename V1, typename V2, typename Sp, size_t... Os>
auto operator+(OpOperatorMixedDerivative<V1, Sp, Os...> const& a,
	OpOperatorMixedDerivative<V2, Sp, Os...> const& b)
{
	return symphas::internal::make_operator_mixed_derivative<Os...>::template get(a.value + b.value, a.solver);
}

template<typename V1, typename V2, typename Sp, size_t... Os>
auto operator-(OpOperatorMixedDerivative<V1, Sp, Os...> const& a,
	OpOperatorMixedDerivative<V1, Sp, Os...> const& b)
{
	return a + (-b);
}

template<typename V1, typename V2, typename Sp, size_t... O1s, size_t... O2s>
auto operator*(OpOperatorMixedDerivative<V1, Sp, O1s...> const& a,
	OpOperatorMixedDerivative<V1, Sp, O2s...> const& b)
{
	return symphas::internal::make_operator_mixed_derivative<(O1s + O2s)...>::template get(a.value * b.value, a.solver);
}


template<typename coeff_t, typename V, typename Sp2, size_t... Os,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorMixedDerivative<V, Sp2, Os...> const& b)
{
	return symphas::internal::make_operator_mixed_derivative<Os...>::template get(value * expr::coeff(b), b.solver);
}

template<typename coeff_t, typename tensor_t, typename Sp2, size_t... Os,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpOperatorMixedDerivative<tensor_t, Sp2, Os...> const& b)
{
	return (value * expr::coeff(b)) * symphas::internal::make_operator_mixed_derivative<Os...>::template get(OpIdentity{}, b.solver);
}

template<Axis ax, size_t... Os, size_t O2, typename V1, typename V2, typename Sp, Axis... axs>
auto mul_mixed_directional(OpOperatorMixedDerivative<V1, Sp, Os...> const& a,
	OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& b, symphas::internal::axis_list<axs...>)
{
	return symphas::internal::make_operator_mixed_derivative<((ax == axs) ? Os + O2 : Os)...>
		::template get(a.value * b.value, a.solver);
}

template<Axis ax, size_t... Os, size_t O2, typename V1, typename V2, typename Sp, Axis... axs>
auto mul_mixed_directional(OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& a,
	OpOperatorMixedDerivative<V1, Sp, Os...> const& b, symphas::internal::axis_list<axs...>)
{
	return symphas::internal::make_operator_mixed_derivative<((ax == axs) ? Os + O2 : Os)...>
		::template get(a.value * b.value, a.solver);
}

template<Axis ax, size_t... Os, size_t O2, typename V1, typename V2, typename Sp>
auto operator*(OpOperatorMixedDerivative<V1, Sp, Os...> const& a,
	OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& b)
{
	return mul_mixed_directional(a, b, symphas::lib::make_axis_list<sizeof...(Os)>());
}

template<Axis ax, size_t... Os, size_t O2, typename V1, typename V2, typename Sp>
auto operator*(OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& a,
	OpOperatorMixedDerivative<V1, Sp, Os...> const& b)
{
	return mul_mixed_directional(a, b, symphas::lib::make_axis_list<sizeof...(Os)>());
}


// ******************************************************************************************
// Rules for applying directional derivatives.
// ******************************************************************************************


// Applying the sum of directional derivatives of the same (even) order is the same as 
// applying the equivalent derivative. I.e. if `O` is two, then this will be the laplacian.
template<
	Axis ax1, Axis ax2, Axis ax3, typename V, typename Sp, typename E,
	typename std::enable_if_t<
		(expr::grid_dim<E>::value <= 3 && expr::grid_dim<E>::value > 0 
			&& (ax1 != ax2 && ax1 != ax3 && ax2 != ax3) && (expr::is_identity<V> || expr::is_fraction<V>)),
	int> = 0>
auto operator*(
	OpOperatorCombination<
		OpOperatorCombination<
			OpOperatorDirectionalDerivative<ax2, 2, V, Sp>,
			OpOperatorDirectionalDerivative<ax3, 2, V, Sp>>,
		OpOperatorDirectionalDerivative<ax1, 2, V, Sp>
	> const& op,
	OpExpression<E> const& e)
{
	return symphas::internal::make_operator_derivative<2>::template get(V{}, op.g.solver) * (*static_cast<E const*>(&e));
}

template<
	Axis ax1, Axis ax2, typename V, typename Sp, typename E,
	typename std::enable_if_t<
		(expr::grid_dim<E>::value <= 2 && expr::grid_dim<E>::value > 0 
			&& (ax1 != ax2) && (ax1 != Axis::Z && ax2 != Axis::Z) 
			&& (expr::is_identity<V> || expr::is_fraction<V>)),
	int> = 0>
auto operator*(
	OpOperatorCombination<
		OpOperatorDirectionalDerivative<ax1, 2, V, Sp>,
		OpOperatorDirectionalDerivative<ax2, 2, V, Sp>
	> const& op,
	OpExpression<E> const& e)
{
	return symphas::internal::make_operator_derivative<2>::template get(V{}, op.f.solver) * (*static_cast<E const*>(&e));
}

template<
	size_t O, typename V, typename Sp, typename E,
	typename std::enable_if_t<(expr::grid_dim<E>::value == 1), int> = 0>
auto operator*(
	OpOperatorDirectionalDerivative<Axis::X, O, V, Sp> const& op,
	OpExpression<E> const& e)
{
	return symphas::internal::make_operator_derivative<O>::template get(op.value, op.solver) * (*static_cast<E const*>(&e));
}


// Summing together directional derivatives will transform into a combination
template<
	typename V, typename Sp, typename E,
	typename Dd1, typename Dd2, Axis ax1 = Dd1::axis, Axis ax2 = Dd2::axis,
	typename std::enable_if_t<
	(expr::grid_dim<E>::value == 2
		&& (ax1 != ax2) && (expr::is_identity<V> || expr::is_fraction<V>)
		&& (ax1 != Axis::Z && ax2 != Axis::Z)
		&& Dd1::is_directional && Dd1::order == 2 && Dd2::is_directional && Dd2::order == 2),
	int> = 0>
auto operator+(
	OpDerivative<Dd1, V, E, Sp> const& dop0,
	OpDerivative<Dd2, V, E, Sp> const& dop1)
{
	return expr::make_operator_derivative<2>(expr::coeff(dop0), dop0.solver) * expr::get_enclosed_expression(dop0);
}


// Summing together directional derivatives will transform into a combination
template<
	typename V, typename Sp, typename E,
	typename Dd1, typename Dd2, Axis ax1 = Dd1::axis, Axis ax2 = Dd2::axis,
	typename std::enable_if_t<
	(expr::grid_dim<E>::value == 3
		&& (ax1 != ax2) && (expr::is_identity<V> || expr::is_fraction<V>)
		&& Dd1::is_directional && Dd1::order == 2 && Dd2::is_directional && Dd2::order == 2),
	int> = 0>
auto operator+(
	OpDerivative<Dd1, V, E, Sp> const& dop0,
	OpDerivative<Dd2, V, E, Sp> const& dop1)
{
	return (expr::make_operator_directional_derivative<ax1, 2>(dop0.solver)
		+ expr::make_operator_directional_derivative<ax2, 2>(dop0.solver))
		(expr::coeff(dop0) * expr::get_enclosed_expression(dop0));
}

// Summing together directional derivatives will transform into a combination
template<
	typename V, typename Sp, typename E,
	typename Dd1, typename Dd2, Axis ax1 = Dd1::axis, Axis ax2 = Dd2::axis,
	typename std::enable_if_t<
	(expr::grid_dim<E>::value == 3
		&& (ax1 != ax2) && (expr::is_identity<V> || expr::is_fraction<V>)
		&& Dd1::is_directional && Dd1::order == 2 && Dd2::is_directional && Dd2::order == 2),
	int> = 0>
auto operator-(
	OpDerivative<Dd1, V, E, Sp> const& dop0,
	OpDerivative<Dd2, decltype(-std::declval<V>()), E, Sp> const& dop1)
{
	return (expr::make_operator_directional_derivative<ax1, 2>(dop0.solver)
		+ expr::make_operator_directional_derivative<ax2, 2>(dop0.solver))
		(expr::coeff(dop0) * expr::get_enclosed_expression(dop0));
}

// Summing together directional derivatives will transform into a combination
template<
	Axis ax1, Axis ax2, typename V, typename Sp, typename E,
	typename Dd, Axis ax3 = Dd::axis,
	typename std::enable_if_t<
	(expr::grid_dim<E>::value == 3
		&& (ax1 != ax2 && ax1 != ax3 && ax2 != ax3) && (expr::is_identity<V> || expr::is_fraction<V>)
		&& Dd::is_directional && Dd::order == 2),
	int> = 0>
auto operator+(
	OpCombination<
		OpOperatorDirectionalDerivative<ax1, 2, V, Sp>,
		OpOperatorDirectionalDerivative<ax2, 2, V, Sp>,
	E> const& combination,
	OpDerivative<Dd, V, E, Sp> const& dop)
{
	return symphas::internal::make_operator_derivative<2>::template get(V{}, dop.solver) * expr::get_enclosed_expression(dop);
}

// Summing together directional derivatives will transform into a combination
template<
	Axis ax1, Axis ax2, typename V, typename Sp, typename E,
	typename Dd, Axis ax3 = Dd::axis,
	typename std::enable_if_t<
	(expr::grid_dim<E>::value == 3
		&& (ax1 != ax2 && ax1 != ax3 && ax2 != ax3) && (expr::is_identity<V> || expr::is_fraction<V>)
		&& Dd::is_directional && Dd::order == 2),
	int> = 0>
auto operator-(
	OpCombination<
		OpOperatorDirectionalDerivative<ax1, 2, V, Sp>,
		OpOperatorDirectionalDerivative<ax2, 2, V, Sp>,
	E> const& combination,
	OpDerivative<Dd, decltype(-std::declval<V>()), E, Sp> const& dop)
{
	return symphas::internal::make_operator_derivative<2>::template get(V{}, dop.solver) * expr::get_enclosed_expression(dop);
}





namespace symphas::internal
{

	template<size_t R, Axis ax, size_t O, typename Sp>
	auto apply_directional_derivative(OpVoid, solver_op_type<Sp> solver)
	{
		return OpVoid{};
	}

	template<size_t R, Axis ax, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank, size_t N = expr::grid_dim<E>::value - R,
		typename std::enable_if_t<(R0 == 0), int> = 0>
	auto apply_directional_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		constexpr size_t D = expr::grid_dim<E>::value;
		return symphas::internal::nth_directional_derivative_apply<ax, O, Sp>::template get(*static_cast<E const*>(&e), solver);
	}

	template<size_t R, Axis ax, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank, size_t R1 = expr::eval_type<E>::template rank_<1>,
		typename std::enable_if_t<(R0 > 0 && R1 == 1), int> = 0>
	auto apply_directional_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		constexpr size_t D = expr::eval_type<E>::rank;
		constexpr size_t Q = expr::eval_type<E>::template rank_<1>;

		auto rtensor = expr::make_row_vector<R - 1, D>();
		auto ctensor = expr::make_column_vector<R - 1, D>();

		auto d = ctensor * apply_directional_derivative<0, ax, O>(rtensor * (*static_cast<E const*>(&e)), solver);
		if constexpr (R > 1)
		{
			return d + apply_directional_derivative<R - 1, ax, O>((*static_cast<E const*>(&e)), solver);
		}
		else
		{
			return d;
		}
	}

	template<size_t R, Axis ax, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank, size_t R1 = expr::eval_type<E>::template rank_<1>,
		typename std::enable_if_t<(R0 > 0 && R1 > 1), int> = 0>
	auto apply_directional_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		auto rtensor = expr::make_row_vector<R - 1, R0>();
		auto ctensor = expr::make_column_vector<R - 1, R1>();

		auto d = apply_directional_derivative<R0, ax, O>((*static_cast<E const*>(&e) * ctensor), solver) * rtensor;
		if constexpr (R > 1)
		{
			return d + apply_directional_derivative<R - 1, ax, O>((*static_cast<E const*>(&e)), solver);
		}
		else
		{
			return d;
		}
	}

	template<size_t R, size_t... Os, typename Sp>
	auto apply_directional_mixed(OpVoid, solver_op_type<Sp> solver)
	{
		return OpVoid{};
	}

	template<size_t R, size_t... Os, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank, size_t N = expr::grid_dim<E>::value - R,
		typename std::enable_if_t<(R0 == 0), int> = 0>
	auto apply_directional_mixed(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		return symphas::internal::nth_mixed_derivative_apply<Sp, Os...>::template get(*static_cast<E const*>(&e), solver);
	}

	template<size_t R, size_t... Os, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
		typename std::enable_if_t<(R0 > 0), int> = 0>
	auto apply_directional_mixed(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		constexpr size_t D = expr::eval_type<E>::rank;
		auto rtensor = expr::make_row_vector<R - 1, D>();
		auto ctensor = expr::make_column_vector<R - 1, D>();

		auto d = ctensor * apply_directional_mixed<0, Os...>(rtensor * (*static_cast<E const*>(&e)), solver);
		if constexpr (R > 1)
		{
			return d + apply_directional_mixed<R - 1, Os...>((*static_cast<E const*>(&e)), solver);
		}
		else
		{
			return d;
		}
	}


	template<size_t R, size_t O, size_t R0, size_t order_parity = O % 2>
	struct apply_derivative;

	template<size_t R, size_t O>
	struct apply_derivative<R, O, 0, 0>
	{
		template<typename E, typename Sp>
		auto operator()(OpVoid, solver_op_type<Sp> solver)
		{
			return OpVoid{};
		}

		template<typename E, typename Sp>
		auto operator()(OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			using deriv_type = symphas::internal::nth_derivative_apply<Axis::X, O, Sp>;
			return deriv_type::template get(*static_cast<E const*>(&e), solver);
		}
	};

	template<size_t R, size_t O, size_t R0>
	struct apply_derivative<R, O, R0, 0>
	{
		template<typename E, typename Sp>
		auto operator()(OpVoid, solver_op_type<Sp> solver)
		{
			return OpVoid{};
		}

		template<typename E, typename Sp>
		auto operator()(OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			constexpr size_t D = expr::eval_type<E>::rank;
			constexpr Axis ax = (R == 1) ? Axis::X : (R == 2) ? Axis::Y : Axis::Z;
			auto rtensor = expr::make_row_vector<R - 1, D>();
			auto ctensor = expr::make_column_vector<R - 1, D>();

			auto d = ctensor * symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(rtensor * (*static_cast<E const*>(&e)), solver);
			if constexpr (R > 1)
			{
				return d + apply_derivative<R - 1, O, R0, 0>{}((*static_cast<E const*>(&e)), solver);
			}
			else
			{
				return d;
			}
		}
	};

	template<size_t R, size_t O>
	struct apply_derivative<R, O, 0, 1>
	{
		template<typename E, typename Sp>
		auto operator()(OpVoid, solver_op_type<Sp> solver)
		{
			return OpVoid{};
		}

		template<typename E, typename Sp>
		auto operator()(OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			constexpr size_t N = expr::grid_dim<E>::value - R;

			constexpr size_t D = expr::grid_dim<E>::value;
			constexpr Axis ax = (N == 1) ? Axis::X : (N == 2) ? Axis::Y : Axis::Z;
			auto ctensor = expr::make_column_vector<N - 1, D>();

			auto d = ctensor * symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(*static_cast<E const*>(&e), solver);
			if constexpr (N > 1)
			{
				return d + apply_derivative<R + 1, O, 0, 1>{}((*static_cast<E const*>(&e)), solver);
			}
			else
			{
				return d;
			}
		}
	};

	template<size_t R, size_t O, size_t R0>
	struct apply_derivative<R, O, R0, 1>
	{
		template<typename E, typename Sp>
		auto operator()(OpVoid, solver_op_type<Sp> solver)
		{
			return OpVoid{};
		}

		template<typename E, typename Sp>
		auto operator()(OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			constexpr size_t D = expr::eval_type<E>::rank;
			constexpr Axis ax = (R == 1) ? Axis::X : (R == 2) ? Axis::Y : Axis::Z;
			auto rtensor = expr::make_row_vector<R - 1, D>();

			auto d = apply_derivative<0, O, 0, 1>{}(rtensor * (*static_cast<E const*>(&e)), solver) * rtensor;
			if constexpr (R > 1)
			{
				return d + apply_derivative<R - 1, O, R0, 1>{}((*static_cast<E const*>(&e)), solver);
			}
			else
			{
				return d;
			}
		}
	};

	//template<size_t R, size_t O, typename Sp>
	//auto apply_derivative(OpVoid, solver_op_type<Sp> solver)
	//{
	//	return OpVoid{};
	//}

	//template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
	//	typename std::enable_if_t<(R0 == 0 && O % 2 == 0), int> = 0>
	//auto apply_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	//{
	//	return symphas::internal::nth_derivative_apply<Axis::X, O, Sp>::template get(*static_cast<E const*>(&e), solver);
	//}

	//template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
	//	typename std::enable_if_t<(R0 > 0 && O % 2 == 0), int> = 0>
	//auto apply_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	//{
	//	constexpr size_t D = expr::eval_type<E>::rank;
	//	constexpr Axis ax = (R == 1) ? Axis::X : (R == 2) ? Axis::Y : Axis::Z;
	//	auto rtensor = expr::make_row_vector<R - 1, D>();
	//	auto ctensor = expr::make_column_vector<R - 1, D>();

	//	auto d = ctensor * symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(rtensor * (*static_cast<E const*>(&e)), solver);
	//	if constexpr (R > 1)
	//	{
	//		return d + apply_derivative<R - 1, O>((*static_cast<E const*>(&e)), solver);
	//	}
	//	else
	//	{
	//		return d;
	//	}
	//}

	//template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank, size_t N = expr::grid_dim<E>::value - R,
	//	typename std::enable_if_t<(R0 == 0 && O % 2 == 1), int> = 0>
	//auto apply_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	//{
	//	constexpr size_t D = expr::grid_dim<E>::value;
	//	constexpr Axis ax = (N == 1) ? Axis::X : (N == 2) ? Axis::Y : Axis::Z;
	//	auto ctensor = expr::make_column_vector<N - 1, D>();

	//	auto d = ctensor * symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(*static_cast<E const*>(&e), solver);
	//	if constexpr (N > 1)
	//	{
	//		return d + apply_derivative<R + 1, O>((*static_cast<E const*>(&e)), solver);
	//	}
	//	else
	//	{
	//		return d;
	//	}
	//}

	//template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
	//	typename std::enable_if_t<(R0 > 0 && O % 2 == 1), int> = 0>
	//auto apply_derivative(OpExpression<E> const& e, solver_op_type<Sp> solver)
	//{
	//	constexpr size_t D = expr::eval_type<E>::rank;
	//	constexpr Axis ax = (R == 1) ? Axis::X : (R == 2) ? Axis::Y : Axis::Z;
	//	auto rtensor = expr::make_row_vector<R - 1, D>();

	//	auto d = apply_derivative<0, O>(rtensor * (*static_cast<E const*>(&e)), solver) * rtensor;
	//	if constexpr (R > 1)
	//	{
	//		return d + apply_derivative<R - 1, O>((*static_cast<E const*>(&e)), solver);
	//	}
	//	else
	//	{
	//		return d;
	//	}
	//}



	template<size_t R, size_t O, typename Sp>
	auto apply_derivative_dot(OpVoid, solver_op_type<Sp> solver)
	{
		return OpVoid{};
	}

	template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
		typename std::enable_if_t<(R0 > 0 && O % 2 == 1), int> = 0>
	auto apply_derivative_dot(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		constexpr size_t D = expr::eval_type<E>::rank;
		constexpr Axis ax = (R == 1) ? Axis::X : (R == 2) ? Axis::Y : Axis::Z;
		auto rtensor = expr::make_row_vector<R - 1, D>();

		auto d = symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(rtensor * (*static_cast<E const*>(&e)), solver);
		if constexpr (R > 1)
		{
			return d + apply_derivative_dot<R - 1, O>((*static_cast<E const*>(&e)), solver);
		}
		else
		{
			return d;
		}
	}

	template<size_t R, size_t O, typename E, typename Sp, size_t R0 = expr::eval_type<E>::rank,
		typename std::enable_if_t<(R0 == 0 || O % 2 == 0), int> = 0>
	auto apply_derivative_dot(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		return apply_derivative<R, O, R0>{}(*static_cast<E const*>(&e), solver);
	}

	template<size_t O, typename E, typename Sp>
	auto apply_derivative_dot(OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		if constexpr (expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0)
		{
			return expr::make_mul(expr::make_operator_derivative<O>(solver), *static_cast<E const*>(&e));
		}
		else
		{
			constexpr size_t R = expr::eval_type<E>::rank;
			if constexpr (R > 0)
			{
				return apply_derivative_dot<R, O>(*static_cast<E const*>(&e), solver);
			}
			else
			{
				return apply_derivative<0, O, 0>{}(*static_cast<E const*>(&e), solver);
			}
		}
	}

	template<size_t O, typename E, typename G>
	auto apply_derivative_dot(OpExpression<E> const& e, SymbolicDerivative<G> solver)
	{
		return expr::make_derivative<O, G>(*static_cast<E const*>(&e));
	}

	template<size_t O>
	struct initialize_derivative_order
	{
		template<typename V, typename E, typename Sp, size_t R = expr::eval_type<E>::rank,
			typename std::enable_if_t<!(expr::is_coeff<E> || expr::is_identity<E>), int> = 0>
		auto operator()(V const& v, OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			if constexpr (expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0)
			{
				return OpOperatorChain(expr::make_operator_derivative<O>(v, solver),
					*static_cast<E const*>(&e));
			}
			else
			{
				return v * apply_derivative<R, O, R>{}(*static_cast<E const*>(&e), solver);
			}
		}
		
		template<typename V, typename E, typename Sp, size_t R = expr::eval_type<E>::rank,
			typename std::enable_if_t<(expr::is_coeff<E> || expr::is_identity<E>), int> = 0>
		auto operator()(V const& v, OpExpression<E> const& e, solver_op_type<Sp> solver)
		{
			return OpVoid{};
		}


		template<typename V, typename E, typename G>
		auto operator()(V const& v, OpExpression<E> const& e, SymbolicDerivative<G>)
		{
			return v * expr::make_derivative<O, G>(*static_cast<E const*>(&e));
		}

		template<typename V, typename E, typename G>
		auto operator()(V const& v, OpOperator<E> const& e, SymbolicDerivative<G>)
		{
			return v * expr::make_derivative<O, G>(*static_cast<E const*>(&e));
		}
	};

	template<size_t O>
	template<typename V, typename Sp>
	inline auto make_operator_derivative<O>::get(V const& v, solver_op_type<Sp> solver)
	{
		return OpOperatorDerivative<O, V, Sp>(v, solver);
	}

	template<size_t O>
	template<typename Sp>
	inline auto make_operator_derivative<O>::get(solver_op_type<Sp> solver)
	{
		return OpOperatorDerivative<O, OpIdentity, Sp>(OpIdentity{}, solver);
	}

	template<Axis ax, size_t O>
	template<typename V, typename Sp>
	inline auto make_operator_directional_derivative<ax, O>::get(V const& v, solver_op_type<Sp> solver)
	{
		return OpOperatorDirectionalDerivative<ax, O, V, Sp>(v, solver);
	}

	template<Axis ax, size_t O>
	template<typename Sp>
	inline auto make_operator_directional_derivative<ax, O>::get(solver_op_type<Sp> solver)
	{
		return OpOperatorDirectionalDerivative<ax, O, OpIdentity, Sp>(OpIdentity{}, solver);
	}

	template<size_t... Os>
	template<typename V, typename Sp>
	inline auto make_operator_mixed_derivative<Os...>::get(V const& v, solver_op_type<Sp> solver)
	{
		return OpOperatorMixedDerivative<V, Sp, Os...>(v, solver);
	}

	template<size_t... Os>
	template<typename Sp>
	inline auto make_operator_mixed_derivative<Os...>::get(solver_op_type<Sp> solver)
	{
		return OpOperatorMixedDerivative<OpIdentity, Sp, Os...>(OpIdentity{}, solver);
	}

	template<typename Dd>
	template<typename V, typename E, typename Sp>
	inline auto make_derivative<Dd>::get(V const& v, OpExpression<E> const& e, solver_op_type<Sp> solver)
	{
		//if constexpr (expr::eval_type_t<E>::rank > 0)
		//{
		//	auto d = expr::make_derivative<Dd>(expr::symbols::Symbol{}, solver);
		//	auto [op, _] = expr::split::separate_operator(d);
		//	return op(*static_cast<E const*>(&e));
		//}
		//else
		{
			return OpDerivative<Dd, V, E, Sp>(v, *static_cast<const E*>(&e), solver);
		}
	}

	template<typename Dd>
	template<typename V, typename E, typename Sp>
	inline auto make_derivative<Dd>::get(V const& v, OpOperator<E> const& e, solver_op_type<Sp> solver)
	{
		return OpDerivative<Dd, V, E, Sp>(v, *static_cast<const E*>(&e), solver);
	}

	template<typename Dd>
	template<typename V, typename S, typename G, typename G0>
	inline auto make_derivative<Dd>::get(V const& v, OpTerm<S, G> const& e, SymbolicDerivative<G0> s)
	{
		auto coeff = v * expr::coeff(e);
		return OpDerivative<Dd, decltype(coeff), OpTerm<OpIdentity, G>, SymbolicDerivative<G0>>(coeff, OpTerms(OpIdentity{}, expr::terms_after_first(e)), s);
	}

	template<typename Dd>
	template<typename V, typename E, typename G0>
	inline auto make_derivative<Dd>::get(V const& v, OpExpression<E> const& e, SymbolicDerivative<G0> s)
	{
		return OpDerivative<Dd, V, E, SymbolicDerivative<G0>>(v, *static_cast<E const*>(&e), s);
	}

	template<typename Dd>
	template<typename V, typename E, typename G0>
	inline auto make_derivative<Dd>::get(V const& v, OpOperator<E> const& e, SymbolicDerivative<G0> s)
	{
		return OpDerivative<Dd, V, E, SymbolicDerivative<G0>>(v, *static_cast<E const*>(&e), s);
	}

	template<typename Dd>
	template<typename V, typename S, typename G, typename Sp>
	inline auto make_derivative<Dd>::get(V const& v, OpTerm<S, G> const& e, solver_op_type<Sp> solver)
	{
		return make_derivative<Dd>::get_g(v * expr::coeff(e), expr::data(e), solver);
	}

	template<typename Dd>
	template<typename V, typename G, typename Sp>
	inline auto make_derivative<Dd>::get_g(V const& v, G g, solver_op_type<Sp> solver)
	{
		return OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>(v, g, solver);
	}

	template<typename Dd>
	template<typename A, typename B>
	inline auto make_derivative<Dd>::get(A&& a, B&& b)
	{
		return get(OpIdentity{}, std::forward<A>(a), std::forward<B>(b));
	}

	template<typename Dd>
	template<typename V, typename coeff_t, typename Sp, typename std::enable_if_t<(expr::is_coeff<coeff_t> || expr::is_identity<coeff_t>), int>>
	inline auto make_derivative<Dd>::get(V const&, coeff_t, solver_op_type<Sp>)
	{
		return OpVoid{};
	}
}

//
//template<size_t O, typename V, typename Sp, typename E, size_t R = expr::eval_type<E>::rank,
//	typename std::enable_if_t<(R > 0), int> = 0>
//auto operator*(E const& e, OpOperatorDerivative<O, V, Sp> const& d)
//{
//	return expr::dot(e, expr::break_up_derivative<O, R>(d.solver));
//}

//template<size_t O, typename V, typename Sp, typename E>
//auto operator*(OpOperatorDerivative<O, V, Sp> const& d, E const& e)
//{
//	if constexpr (expr::is_coeff<E>)
//	{
//		return e * d;
//	}
//	else if constexpr (std::is_convertible_v<E, OpOperator<E>>)
//	{
//		return d.operator*(e);
//	}
//	else
//	{
//		return expr::coeff(d) * symphas::internal::apply_derivative_dot<O>(e, d.solver);
//	}
//}



template<size_t O, typename V, typename Sp>
template<typename E>
auto OpOperatorDerivative<O, V, Sp>::operator*(OpExpression<E> const& e) const
{
	if constexpr (expr::is_coeff<E>)
	{
		return *static_cast<E const*>(&e) * (*this);
	}
	else
	{
		return value * symphas::internal::apply_derivative_dot<O>(e, solver);
	}
}


template<Axis ax, size_t O, typename V, typename Sp>
template<typename E>
auto OpOperatorDirectionalDerivative<ax, O, V, Sp>::operator*(OpExpression<E> const& e) const
{
	if constexpr (expr::is_coeff<E>)
	{
		return *static_cast<E const*>(&e) * *this;
	}
	//else if constexpr (std::is_convertible_v<E, OpOperator<E>>)
	//{
	//	return d.operator*(e);
	//}
	else
	{
		return apply_impl(e);
		//return d.operator*(e);
	}
}

template<typename V, typename Sp, size_t... Os>
template<typename E>
auto OpOperatorMixedDerivative<V, Sp, Os...>::operator*(OpExpression<E> const& e) const
{
	if constexpr (expr::is_coeff<E>)
	{
		return *static_cast<E const*>(&e) * (*this);
	}
	else
	{
		return apply_impl(e);
		//return value * symphas::internal::apply_mixed_derivative<O>(e, solver);
	}
}
//
//template<size_t O, typename V, typename Sp, typename T>
//auto operator*(OpOperatorDerivative<O, V, Sp> const& d, OpLiteral<T> const& e)
//{
//	return e * d;
//}
//
//template<Axis ax, size_t O, typename V, typename Sp, typename E>
//auto operator*(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& d, E const& e)
//{
//	if constexpr (expr::is_coeff<E>)
//	{
//		return e * d;
//	}
//	else if constexpr (std::is_convertible_v<E, OpOperator<E>>)
//	{
//		return d.operator*(e);
//	}
//	else
//	{
//		return d.operator*(e);
//	}
//}
//
//template<Axis ax, size_t O, typename V, typename Sp, typename T>
//auto operator*(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& d, OpLiteral<T> const& e)
//{
//	return e * d;
//}
//
//template<typename V, typename Sp, size_t... Os, typename E>
//auto operator*(OpOperatorMixedDerivative<V, Sp, Os...> const& d, E const& e)
//{
//	if constexpr (expr::is_coeff<E>)
//	{
//		return e * d;
//	}
//	else if constexpr (std::is_convertible_v<E, OpOperator<E>>)
//	{
//		return d.operator*(e);
//	}
//	else
//	{
//		return d.operator*(e);
//	}
//}
//
//template<typename V, typename Sp, size_t... Os, typename T>
//auto operator*(OpOperatorMixedDerivative<V, Sp, Os...> const& d, OpLiteral<T> const& e)
//{
//	return e * d;
//}




template<size_t O, typename V, typename Sp>
template<typename E>
auto OpOperatorDerivative<O, V, Sp>::_apply_impl(OpExpression<E> const& e) const
{
	return symphas::internal::initialize_derivative_order<O>{}(value, *static_cast<E const*>(&e), solver);
}

template<Axis ax, size_t O, typename V, typename Sp>
template<typename E, typename std::enable_if_t<!(expr::is_symbol<expr::eval_type_t<E>>&& expr::grid_dim<E>::value == 0), int>>
auto OpOperatorDirectionalDerivative<ax, O, V, Sp>::_apply_impl(OpExpression<E> const& e) const
{
	constexpr size_t R = expr::eval_type<E>::rank;
	return value * symphas::internal::apply_directional_derivative<R, ax, O>(*static_cast<E const*>(&e), solver);
}

template<typename V, typename Sp, size_t... Os>
template<typename E, typename std::enable_if_t<!(expr::is_symbol<expr::eval_type_t<E>> && expr::grid_dim<E>::value == 0), int>>
auto OpOperatorMixedDerivative<V, Sp, Os...>::_apply_impl(OpExpression<E> const& e) const
{
	constexpr size_t R = expr::eval_type<E>::rank;
	return value * symphas::internal::apply_directional_mixed<R, Os...>(*static_cast<E const*>(&e), solver);
}



namespace expr
{


	//! Function wrappers to return the laplacian concrete derivative.
	template<typename E, typename Sp>
	auto make_laplacian(E&& e, solver_op_type<Sp> solver)
	{
		return symphas::internal::initialize_derivative_order<2>{}(OpIdentity{}, std::forward<E>(e), solver);
		//	return symphas::internal::laplacian_apply<Axis::X, Sp>::template get();
	}

	//! Function wrappers to return the bilaplacian concrete derivative.
	template<typename E, typename Sp>
	auto make_bilaplacian(E&& e, solver_op_type<Sp> solver)
	{
		return symphas::internal::initialize_derivative_order<4>{}(OpIdentity{}, std::forward<E>(e), solver);
		//return symphas::internal::bilaplacian_apply<Axis::X, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the gradlaplacian concrete derivative.
	template<Axis ax, typename E, typename Sp>
	auto make_gradlaplacian(E&& e, solver_op_type<Sp> solver)
	{
		return symphas::internal::initialize_derivative_order<3>{}(OpIdentity{}, std::forward<E>(e), solver);
		//return symphas::internal::gradlaplacian_apply<ax, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the gradient concrete derivative.
	template<Axis ax, typename E, typename Sp>
	auto make_gradient(E&& e, solver_op_type<Sp> solver)
	{
		return symphas::internal::initialize_derivative_order<1>{}(OpIdentity{}, std::forward<E>(e), solver);
		//return symphas::internal::gradient_apply<ax, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}
}


namespace symphas::math
{
	//! Apply the curl of a vector.
	/*!
	 * Compound assignment operator to calculate the conventional cross product and assign it
	 * to the calling instance.
	 *
	 * Assignment operator to compute the cross product, which is only applicable to the
	 * 3-dimensional VectorValue template specialization. Operator is called by the left hand
	 * instance, taking data from the right hand VectorValue to compute the cross product and
	 * assigning the result to the left hand side.
	 *
	 * \param rhs The VectorValue instance where data is taken from to compute the cross product
	 * the data from the left hand instance.
	 */
	template<size_t O, typename V, typename Sp, typename E, 
		typename = std::enable_if_t<(expr::eval_type<E>::dimension == 3), int>>
	auto cross(OpOperatorDerivative<O, V, Sp> const& d_op, OpExpression<E> const& e)
	{
		any_vector_t<scalar_t, 3> result;
		auto lhs0 = expr::make_operator_directional_derivative<Axis::X, 1>(d_op.solver);
		auto lhs1 = expr::make_operator_directional_derivative<Axis::X, 1>(d_op.solver);
		auto lhs2 = expr::make_operator_directional_derivative<Axis::X, 1>(d_op.solver);

		//result[0] = (lhs[1] * rhs[2] - lhs[2] * rhs[1]);
		//result[1] = (lhs[2] * rhs[0] - lhs[0] * rhs[2]);
		//result[2] = (lhs[0] * rhs[1] - lhs[1] * rhs[0]);

		return result;
	}

}


namespace expr
{

	template<typename E>
	int _get_solver(OpExpression<E> const& e);
	template<typename E>
	int _get_solver(OpOperator<E> const& e);
	int _get_solver(OpAddList<> const& e);
	template<typename E0, typename... Es>
	auto _get_solver(OpAddList<E0, Es...> const& e);
	template<typename... Es>
	auto _get_solver(OpAdd<Es...> const& e);
	template<typename A, typename B>
	auto _get_solver(OpBinaryMul<A, B> const& e);
	template<typename A, typename B>
	auto _get_solver(OpBinaryDiv<A, B> const& e);
	template<typename Dd, typename V, typename E, typename Sp>
	auto _get_solver(OpDerivative<Dd, V, E, Sp> const& e);
	template<size_t O, typename V, typename Sp>
	auto _get_solver(OpOperatorDerivative<O, V, Sp> const& e);
	template<Axis ax, size_t O, typename V, typename Sp>
	auto _get_solver(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e);
	template<typename V, typename Sp, size_t... Os>
	auto _get_solver(OpOperatorMixedDerivative<V, Sp, Os...> const& e);
	template<typename A1, typename A2>
	auto _get_solver(OpOperatorChain<A1, A2> const& e);
	template<typename A1, typename A2, typename E>
	auto _get_solver(OpChain<A1, A2, E> const& e);
	template<typename A1, typename A2>
	auto _get_solver(OpOperatorCombination<A1, A2> const& e);
	template<typename A1, typename A2, typename E>
	auto _get_solver(OpCombination<A1, A2, E> const& e);
	template<typename V, typename sub_t, typename F>
	auto _get_solver(OpSymbolicEval<V, sub_t, F> const& e);
	template<typename V, typename E, typename T>
	auto _get_solver(OpIntegral<V, E, T> const& e);

	template<typename Sp, typename E>
	Sp const& _get_solver(Sp const& solver, OpExpression<E> const& e)
	{
		return solver;
	}

	template<typename E>
	auto _get_solver(int, OpExpression<E> const& e)
	{
		return _get_solver(*static_cast<E const*>(&e));
	}

	template<typename Sp, typename... Es>
	Sp const& _get_solver(Sp const& solver, OpAddList<Es...> const& e)
	{
		return solver;
	}

	template<typename... Es>
	auto _get_solver(int, OpAddList<Es...> const& e)
	{
		return _get_solver(e);
	}


	template<typename E>
	int _get_solver(OpExpression<E> const& e)
	{
		return 0;
	}

	template<typename E>
	int _get_solver(OpOperator<E> const& e)
	{
		return 0;
	}

	inline int _get_solver(OpAddList<> const& e)
	{
		return 0;
	}

	template<typename E0, typename... Es>
	auto _get_solver(OpAddList<E0, Es...> const& e)
	{
		auto e0 = expr::get<0>(e);
		return _get_solver(_get_solver(e0), expr::terms_after_first(e));
	}

	template<typename... Es>
	auto _get_solver(OpAdd<Es...> const& e)
	{
		return _get_solver(*static_cast<OpAddList<Es...> const*>(&e));
	}

	template<typename A, typename B>
	auto _get_solver(OpBinaryMul<A, B> const& e)
	{
		return _get_solver(_get_solver(e.a), e.b);
	}

	template<typename A, typename B>
	auto _get_solver(OpBinaryDiv<A, B> const& e)
	{
		return _get_solver(_get_solver(e.a), e.b);
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto _get_solver(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return _get_solver(e.solver, expr::get_enclosed_expression(e));
	}

	template<size_t O, typename V, typename Sp>
	auto _get_solver(OpOperatorDerivative<O, V, Sp> const& e)
	{
		return e.solver;
	}

	template<Axis ax, size_t O, typename V, typename Sp>
	auto _get_solver(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e)
	{
		return e.solver;
	}

	template<typename V, typename Sp, size_t... Os>
	auto _get_solver(OpOperatorMixedDerivative<V, Sp, Os...> const& e)
	{
		return e.solver;
	}

	template<typename A1, typename A2>
	auto _get_solver(OpOperatorChain<A1, A2> const& e)
	{
		return _get_solver(_get_solver(e.f), e.g);
	}

	template<typename A1, typename A2, typename E>
	auto _get_solver(OpChain<A1, A2, E> const& e)
	{
		return _get_solver(_get_solver(e.combination), e.e);
	}

	template<typename A1, typename A2>
	auto _get_solver(OpOperatorCombination<A1, A2> const& e)
	{
		return _get_solver(_get_solver(e.f), e.g);
	}

	template<typename A1, typename A2, typename E>
	auto _get_solver(OpCombination<A1, A2, E> const& e)
	{
		return _get_solver(_get_solver(e.combination), e.e);
	}

	template<typename V, typename sub_t, typename F>
	auto _get_solver(OpSymbolicEval<V, sub_t, F> const& e)
	{
		return _get_solver(e.f.e);
	}

	template<typename V, typename E, typename T>
	auto _get_solver(OpIntegral<V, E, T> const& e)
	{
 		return _get_solver(expr::get_enclosed_expression(e));
	}

	template<typename E>
	auto get_solver(E const& e)
	{
		return _get_solver(e);
	}
}



