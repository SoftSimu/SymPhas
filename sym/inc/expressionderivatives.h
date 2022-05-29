
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
		static auto get(V v, OpExpression<E> const& e, solver_op_type<Sp> solver);

		//! Constructs the derivative applied to a variable.
		template<typename V, typename S, typename G, typename Sp>
		static auto get(V v, OpLVariable<S, G> const& e, solver_op_type<Sp> solver);

		//! Constructs the derivative applied to a constant.
		template<typename V, typename S, typename Sp>
		static auto get(V v, OpLiteral<S> const& e, solver_op_type<Sp> solver);
		//! Constructs the derivative applied to a constant.
		template<typename V, typename Sp>
		static auto get(V v, OpVoid const e, solver_op_type<Sp> solver);
		//! Constructs the derivative applied to a constant.
		template<typename V, typename Sp>
		static auto get(V v, OpIdentity const e, solver_op_type<Sp> solver);
		//! Constructs the derivative applied to a constant.
		template<typename V, typename Sp>
		static auto get(V v, OpNegIdentity const e, solver_op_type<Sp> solver);


		//! Constructs the derivative using a grid instead of an expression.
		/*!
		 * Used for the derivative specialization for the OpLVariable.
		 */
		template<typename V, typename G, typename Sp>
		static auto get_g(V v, G g, solver_op_type<Sp> solver);
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
		static auto get(V v, solver_op_type<Sp> t);

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver);


		//! Constructs a concrete derivative object.
		/*!
		 * A derivative of the prescribed order is applied to the given expression,
		 * using the accompanying solver type to generate a derivative of the 
		 * expression.
		 * 
		 * \param v The coefficient of the derivative.
		 * \param e The expression the derivative is applied to.
		 * \param solver The solver applying the derivative.
		 */
		template<typename V, typename E, typename Sp>
		static decltype(auto) apply(V v, OpExpression<E> const& e, solver_op_type<Sp> solver);

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
		static auto get(V v, solver_op_type<Sp> t);

		template<typename Sp>
		static auto get(solver_op_type<Sp> solver);


		//! Constructs a concrete derivative object.
		/*!
		 * A derivative of the prescribed order is applied to the given expression,
		 * using the accompanying solver type to generate a derivative of the
		 * expression.
		 *
		 * \param v The coefficient of the derivative.
		 * \param e The expression the derivative is applied to.
		 * \param solver The solver applying the derivative.
		 */
		template<typename V, typename E, typename Sp>
		static decltype(auto) apply(V v, OpExpression<E> const& e, solver_op_type<Sp> solver);

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


	//! Create a derivative operator of the prescribed dimension.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 * 
	 * \tparam O The order of the derivative to generate
	 */
	template<size_t O, typename A>
	auto make_operator_derivative(A&& solver)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<A>(solver));
	}

	//! Create a derivative operator of the prescribed dimension.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The solver which approximates the derivative.
	 *
	 * \tparam O The order of the derivative to generate
	 */
	template<size_t O, typename V, typename A>
	auto make_operator_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_derivative<O>::template get(std::forward<V>(v), std::forward<A>(solver));
	}

	//! Create a directional derivative operator of the prescribed dimension.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param solver The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate
	 */
	template<Axis ax, size_t O, typename A>
	auto make_operator_directional_derivative(A&& solver)
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(std::forward<A>(solver));
	}

	//! Create a directional derivative operator of the prescribed dimension.
	/*!
	 * Create a derivative of the given order, using the given solver
	 * to apply the numerical approximation.
	 *
	 * \param v The coefficient to the derivative.
	 * \param a The solver which approximates the derivative.
	 *
	 * \tparam ax The axis along which the directional derivative is applied.
	 * \tparam O The order of the derivative to generate
	 */
	template<Axis ax, size_t O, typename V, typename A>
	auto make_operator_directional_derivative(V&& v, A&& solver)
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(std::forward<V>(v), std::forward<A>(solver));

	}
}


template<typename E>
using deriv_working_grid = typename expr::grid_type<E>::type;

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
struct OpFuncDerivative : OpExpression<OpFuncDerivative<Dd, V, E, Sp>>
{
	using result_grid = deriv_working_grid<E>;	
	static const size_t order = Dd::order; //!< The order of this derivative.

	//! Generate a derivative of the given expression.
	/*!
	 * Create a derivative expression represesnting applying a derivative
	 * to the given expression, where the order depends on the differentiation
	 * object.
	 */
	OpFuncDerivative(V value, E const& e, solver_op_type<Sp> solver) : 
		grid{ expr::property::data_dimensions(e) }, value{ value }, e{ e }, solver{ solver } { /*update();*/ }

	OpFuncDerivative(E const& e, solver_op_type<Sp> solver) : 
		grid{ expr::property::data_dimensions(e) }, value{ OpIdentity{} }, e{ e }, solver{ solver } { /*update();*/ }

	inline auto eval(iter_type n) const
	{
		return value * Dd{}(solver, grid, n);
	}

	void update()
	{
		expr::result(e, grid);
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<Dd>::get(-value, e, solver);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<order>::print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<order>::print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ symphas::internal::print_deriv<order>::print_length() + e.print_length();
	}

#endif
	
	friend struct expr::compound_get;


protected:

	result_grid grid;				//!< Grid storing the intermediate values.

public:

	V value;						//!< Value multiplying the result of this derivative.

protected:

	E e;							//!< Expression object specifying grid values.

public:

	solver_op_type<Sp> solver;		//!< Solver that applies the derivative function.

};

//! Specialization of an OpLVariable of the concrete derivative expression.
/*!
 * \tparam f The derivative function, typically a member function of the solver.
 * \tparam V Type of the coefficient.
 * \tparam S Coefficient of the .
 * \tparam Sp The solver type.
 */
template<typename Dd, typename V, typename G, typename Sp>
struct OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp> : OpExpression<OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>>
{
	using E = OpLVariable<OpIdentity, G>;
	using result_grid = deriv_working_grid<E>;
	static const size_t order = Dd::order; //!< The order of this derivative.

	template<typename V0, typename V1, typename std::enable_if<std::is_convertible<mul_result_t<V0, V1>, V>::value, int>::type = 0>
	OpFuncDerivative(V0 value, OpLVariable<V1, G> const& e, solver_op_type<Sp> solver) : data{ e.data }, value{ e.value * value }, solver{ solver } {}
	OpFuncDerivative(V value, G data, solver_op_type<Sp> solver) : data{ data }, value{ value }, solver{ solver } {}

	inline auto eval(iter_type n) const
	{
		return value * Dd{}(solver, expr::BaseData<G>::get(data), n);
	}

	auto operator-() const
	{
		return symphas::internal::make_derivative<Dd>::get_g(-value, data, solver);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<order>::print(out);
		n += fprintf(out, SYEX_DERIV_APPLIED_EXPR_FMT_A "%s" SYEX_DERIV_APPLIED_EXPR_FMT_B, expr::get_op_name(data));
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<order>::print(out + n);
		n += sprintf(out + n, SYEX_DERIV_APPLIED_EXPR_FMT_A "%s" SYEX_DERIV_APPLIED_EXPR_FMT_B, expr::get_op_name(data));
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + STR_ARR_LEN(SYEX_DERIV_APPLIED_EXPR_FMT_A SYEX_DERIV_APPLIED_EXPR_FMT_B) - 1
			+ symphas::internal::print_deriv<order>::print_length() + std::strlen(expr::get_op_name(data));
	}

#endif

	friend struct expr::compound_get;

	template<typename S1, typename Dd2, typename V2, typename G2, typename T2>
	friend auto operator*(OpLiteral<S1> const& a, OpFuncDerivative<Dd2, V2, OpLVariable<OpIdentity, G2>, T2> const& b);

protected:

	G data;

public:

	V value;						//!< Value multiplying the result of this derivative.
	solver_op_type<Sp> solver;		//!< Solver that applies the derivative function.

};


/* multiplication of a derivative object by a literal
 */
template<typename S1, typename Dd2, typename V2, typename E2, typename T2>
auto operator*(OpLiteral<S1> const& a, OpFuncDerivative<Dd2, V2, E2, T2> const& b)
{
	return symphas::internal::make_derivative<Dd2>::template get(a.value * b.value, expr::compound_get::template expr(b), b.solver);
}


/* multiplication of a derivative object that owns only a grid by a literal
 */
template<typename S1, typename Dd2, typename V2, typename G2, typename T2>
auto operator*(OpLiteral<S1> const& a, OpFuncDerivative<Dd2, V2, OpLVariable<OpIdentity, G2>, T2> const& b)
{
	return symphas::internal::make_derivative<Dd2>::template get_g(a.value * b.value, b.data, b.solver);
}


// ******************************************************************************************

namespace expr
{

	//! Alias for constructing nth derivative. 
	template<Axis ax, size_t O, typename Sp>
	using nth_directional_derivative_apply = symphas::internal::make_derivative<typename Solver<Sp>::template directional_derivative<ax, O>>;

	//! Alias for constructing nth derivative. 
	template<size_t O, typename Sp>
	using nth_derivative_apply = symphas::internal::make_derivative<typename Solver<Sp>::template derivative<O>>;

	//! Alias for constructing gradient.
	template<typename Sp>
	using gradient_apply = symphas::internal::make_derivative<typename Solver<Sp>::template derivative<1>>;

	//! Alias for constructing laplacian.
	template<typename Sp>
	using laplacian_apply = symphas::internal::make_derivative<typename Solver<Sp>::template derivative<2>>;

	//! Alias for constructing gradlaplacian.
	template<typename Sp>
	using gradlaplacian_apply = symphas::internal::make_derivative<typename Solver<Sp>::template derivative<3>>;

	//! Alias for constructing bilaplacian.
	template<typename Sp>
	using bilaplacian_apply = symphas::internal::make_derivative<typename Solver<Sp>::template derivative<4>>;




	// ******************************************************************************************

	//! Function wrappers to return the laplacian concrete derivative.
	template<typename E, typename Sp>
	auto gradient_x(E&& e, solver_op_type<Sp> solver)
	{
		return expr::nth_directional_derivative_apply<Axis::X, 1, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the laplacian concrete derivative.
	template<typename E, typename Sp>
	auto gradient_y(E&& e, solver_op_type<Sp> solver)
	{
		return expr::nth_directional_derivative_apply<Axis::Y, 1, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the laplacian concrete derivative.
	template<typename E, typename Sp>
	auto gradient_z(E&& e, solver_op_type<Sp> solver)
	{
		return expr::nth_directional_derivative_apply<Axis::Z, 1, Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the laplacian concrete derivative.
	template<typename E, typename Sp>
	auto laplacian(E&& e, solver_op_type<Sp> solver)
	{
		return expr::laplacian_apply<Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the bilaplacian concrete derivative.
	template<typename E, typename Sp>
	auto bilaplacian(E&& e, solver_op_type<Sp> solver)
	{
		return expr::bilaplacian_apply<Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the gradlaplacian concrete derivative.
	template<typename E, typename Sp>
	auto gradlaplacian(E&& e, solver_op_type<Sp> solver)
	{
		return expr::gradlaplacian_apply<Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

	//! Function wrappers to return the gradient concrete derivative.
	template<typename E, typename Sp>
	auto gradient(E&& e, solver_op_type<Sp> solver)
	{
		return expr::gradient_apply<Sp>::template get(OpIdentity{}, std::forward<E>(e), solver);
	}

}



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
	static_assert(O <= DERIV_MAX_ORDER);
	OpOperatorDerivative(V value, solver_op_type<Sp> solver) : value{ value }, solver{ solver } {}

	inline auto eval(iter_type) const {}

	auto operator-() const
	{
		return symphas::internal::make_operator_derivative<O>::template get(
			-value, solver);
	}

	template<typename VV>
	auto operator+(OpOperatorDerivative<O, VV, Sp> const& a) const
	{
		return symphas::internal::make_operator_derivative<O>::template get(
			value + a.value, solver);
	}

	template<typename VV>
	auto operator-(OpOperatorDerivative<O, VV, Sp> const& a) const
	{
		return symphas::internal::make_operator_derivative<O>::template get(
			value - a.value, solver);
	}

	template<size_t O2, typename V2>
	auto operator*(OpOperatorDerivative<O2, V2, Sp> const& a) const
	{
		return symphas::internal::make_operator_derivative<O + O2>::template get(
			value * a.value, solver);
	}

	template<typename Dd, typename V2, typename E, std::enable_if_t<expr::is_directional_derivative<Dd>::value>>
	auto operator*(OpFuncDerivative<Dd, V2, E, Sp> const& d) const
	{
		return symphas::internal::make_operator_derivative<O + Dd::order>::template get(
			value * d.value, solver) * expr::compound_get::template expr(d);
	}

	//! No distribution property, operator is directly applied.
	template<typename E1, typename E2>
	auto operator*(OpBinaryAdd<E1, E2> const& e) const
	{
		return apply(e);
		//return symphas::internal::make_operator_derivative<O>::template get(value, solver) * e.a
		//	+ symphas::internal::make_operator_derivative<O>::template get(value, solver) * e.b;
	}

	//! No distribution property, operator is directly applied.
	template<typename E1, typename E2>
	auto operator*(OpBinarySub<E1, E2> const& e) const
	{
		return apply(e);
		//return symphas::internal::make_operator_derivative<O>::template get(value, solver) * e.a
		//	- symphas::internal::make_operator_derivative<O>::template get(value, solver) * e.b;
	}

	//! Apply to an expression.
	/*!
	 * Apply the generalized derivative to an expression to produce a
	 * concrete derivative.
	 */
	template<typename E>
	auto apply(OpExpression<E> const& e) const
	{
		return symphas::internal::make_operator_derivative<O>::template apply(value, *static_cast<E const*>(&e), solver);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out + n);
		return n;
	}


	size_t print_length() const
	{
		return expr::coeff_print_length(value) + symphas::internal::print_deriv<O>::print_length();
	}

#endif

	V value;
	solver_op_type<Sp> solver;
};


template<typename S1, size_t O2, typename V2, typename T2>
auto operator*(OpLiteral<S1> const& a, OpOperatorDerivative<O2, V2, T2> const& b)
{
	return symphas::internal::make_operator_derivative<O2>::template get(b.value * a.value, b.solver);
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
	static_assert(O <= DERIV_MAX_ORDER);
	OpOperatorDirectionalDerivative(V value, solver_op_type<Sp> solver) : value{ value }, solver{ solver } {}

	inline auto eval(iter_type) const {}

	auto operator-() const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(
			-value, solver);
	}

	template<typename VV>
	auto operator+(OpOperatorDirectionalDerivative<ax, O, VV, Sp> const& a) const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(
			value + a.value, solver);
	}

	template<typename VV>
	auto operator-(OpOperatorDirectionalDerivative<ax, O, VV, Sp> const& a) const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template get(
			value - a.value, solver);
	}

	template<size_t O2, typename V2>
	auto operator*(OpOperatorDirectionalDerivative<ax, O2, V2, Sp> const& a) const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O + O2>::template get(
			value * a.value, solver);
	}

	template<typename Dd, typename V2, typename E, std::enable_if_t<expr::is_directional_derivative<Dd>::value>>
	auto operator*(OpFuncDerivative<Dd, V2, E, Sp> const& d) const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O + Dd::order>::template get(
			value * d.value, solver) * expr::compound_get::template expr(d);
	}

	//! No distribution property, operator is directly applied.
	template<typename E1, typename E2>
	auto operator*(OpBinaryAdd<E1, E2> const& e) const
	{
		return apply(e);
		//return symphas::internal::make_operator_directional_derivative<O>::template get(value, solver) * e.a
		//	+ symphas::internal::make_operator_directional_derivative<O>::template get(value, solver) * e.b;
	}

	//! No distribution property, operator is directly applied.
	template<typename E1, typename E2>
	auto operator*(OpBinarySub<E1, E2> const& e) const
	{
		return apply(e);
		//return symphas::internal::make_operator_directional_derivative<O>::template get(value, solver) * e.a
		//	- symphas::internal::make_operator_directional_derivative<O>::template get(value, solver) * e.b;
	}

	//! Apply to an expression.
	/*!
	 * Apply the generalized derivative to an expression to produce a
	 * concrete derivative.
	 */
	template<typename E>
	auto apply(OpExpression<E> const& e) const
	{
		return symphas::internal::make_operator_directional_derivative<ax, O>::template apply(value, *static_cast<E const*>(&e), solver);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_deriv<O>::print(out + n);
		return n;
	}


	size_t print_length() const
	{
		return expr::coeff_print_length(value) + symphas::internal::print_deriv<O>::print_length();
	}

#endif

	V value;
	solver_op_type<Sp> solver;
};



template<typename S1, Axis ax2, size_t O2, typename V2, typename T2>
auto operator*(OpLiteral<S1> const& a, OpOperatorDirectionalDerivative<ax2, O2, V2, T2> const& b)
{
	return symphas::internal::make_operator_directional_derivative<ax2, O2>::template get(b.value * a.value, b.solver);
}





namespace symphas::internal
{

	template<size_t O>
	struct initialize_derivative_order
	{
		template<typename V, typename E, typename Sp>
		auto operator()(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
		{
			return expr::nth_derivative_apply<O, Sp>::template get(v, *static_cast<E const*>(&e), t);
		}
	};

	template<>
	struct initialize_derivative_order<1>
	{
		template<typename V, typename E, typename Sp>
		auto operator()(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
		{
			return expr::gradient_apply<Sp>::template get(v, *static_cast<E const*>(&e), t);
		}
	};


	template<>
	struct initialize_derivative_order<2>
	{
		template<typename V, typename E, typename Sp>
		auto operator()(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
		{
			return expr::laplacian_apply<Sp>::template get(v, *static_cast<E const*>(&e), t);
		}
	};


	template<>
	struct initialize_derivative_order<3>
	{
		template<typename V, typename E, typename Sp>
		auto operator()(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
		{
			return expr::gradlaplacian_apply<Sp>::template get(v, *static_cast<E const*>(&e), t);
		}
	};

	template<>
	struct initialize_derivative_order<4>
	{
		template<typename V, typename E, typename Sp>
		auto operator()(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
		{
			return expr::bilaplacian_apply<Sp>::template get(v, *static_cast<E const*>(&e), t);
		}
	};



	template<size_t O>
	template<typename V, typename Sp>
	inline auto make_operator_derivative<O>::get(V v, solver_op_type<Sp> t)
	{
		return OpOperatorDerivative<O, V, Sp>(v, t);
	}

	template<size_t O>
	template<typename Sp>
	inline auto make_operator_derivative<O>::get(solver_op_type<Sp> t)
	{
		return OpOperatorDerivative<O, OpIdentity, Sp>(OpIdentity{}, t);
	}


	template<Axis ax, size_t O>
	template<typename V, typename Sp>
	inline auto make_operator_directional_derivative<ax, O>::get(V v, solver_op_type<Sp> t)
	{
		return OpOperatorDirectionalDerivative<ax, O, V, Sp>(v, t);
	}

	template<Axis ax, size_t O>
	template<typename Sp>
	inline auto make_operator_directional_derivative<ax, O>::get(solver_op_type<Sp> t)
	{
		return OpOperatorDirectionalDerivative<ax, O, OpIdentity, Sp>(OpIdentity{}, t);
	}


	template<typename Dd>
	template<typename V, typename E, typename Sp>
	inline auto make_derivative<Dd>::get(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
	{
		return OpFuncDerivative<Dd, V, E, Sp>(v, *static_cast<const E*>(&e), t);
	}

	template<typename Dd>
	template<typename V, typename S, typename G, typename Sp>
	inline auto make_derivative<Dd>::get(V v, OpLVariable<S, G> const& e, solver_op_type<Sp> t)
	{
		return make_derivative<Dd>::get_g(v * e.value, e.data, t);
	}

	template<typename Dd>
	template<typename V, typename G, typename Sp>
	inline auto make_derivative<Dd>::get_g(V v, G g, solver_op_type<Sp> t)
	{
		return OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>(v, g, t);
	}

	template<typename Dd>
	template<typename A, typename B>
	inline auto make_derivative<Dd>::get(A&& a, B&& b)
	{
		return get(OpIdentity{}, std::forward<A>(a), std::forward<B>(b));
	}

	template<typename Dd>
	template<typename V, typename S, typename Sp>
	inline auto make_derivative<Dd>::get(V, OpLiteral<S> const&, solver_op_type<Sp>)
	{
		return OpVoid{};
	}
	template<typename Dd>
	template<typename V, typename Sp>
	inline auto make_derivative<Dd>::get(V, OpVoid const, solver_op_type<Sp>)
	{
		return OpVoid{};
	}
	template<typename Dd>
	template<typename V, typename Sp>
	inline auto make_derivative<Dd>::get(V, OpIdentity const, solver_op_type<Sp>)
	{
		return OpVoid{};
	}
	template<typename Dd>
	template<typename V, typename Sp>
	inline auto make_derivative<Dd>::get(V, OpNegIdentity const, solver_op_type<Sp>)
	{
		return OpVoid{};
	}


	template<size_t O>
	template<typename V, typename E, typename Sp>
	inline decltype(auto) make_operator_derivative<O>::apply(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
	{
		return initialize_derivative_order<O>{}(v, *static_cast<E const*>(&e), t);
	}

	template<Axis ax, size_t O>
	template<typename V, typename E, typename Sp>
	inline decltype(auto) make_operator_directional_derivative<ax, O>::apply(V v, OpExpression<E> const& e, solver_op_type<Sp> t)
	{
		return expr::nth_directional_derivative_apply<ax, O, Sp>::template get(v, *static_cast<E const*>(&e), t);
	}
}






// directional derivative rules -> must be based on the dimension.








