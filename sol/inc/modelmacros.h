	
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
 * PURPOSE: Defines macros used in specifying models
 * with the equations of motion.
 *
 * ***************************************************************************
 */


#pragma once


/*
 *
 *
 *
 */

#include "modelpfc.h"
#include "modelvirtual.h"
#include "stencilincludes.h"
#include "expressiontypeincludes.h"


//! Copy construct the model with a new solver.
/*!
 * A new model is constructed with using instead the given solver. The phase-field
 * data from the given model is copied over as well.
 *
 * \tparam M The model type that from which is counted the phase fields
 */
template<typename Sp>
struct model_swap_solver
{

protected:

	template<template<size_t, typename> typename M, size_t D, typename Sp0>
	static constexpr auto with_new_solver(M<D, Sp0> model)
	{
		auto parameters = model.generate_parameters();
		return M<D, Sp>(model.get_coeff(), model.get_num_coeff(), parameters);
	}

	template<
		template<template<typename> typename, typename> typename SpecializedModel,
		template<typename> typename Eq,
		size_t D, typename Sp0, typename... S>
	static constexpr auto with_new_solver(SpecializedModel<Eq, Eq<symphas::internal::MakeEquation<Model<D, Sp0, S...>>>> const& model)
	{
		using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized<Eq>;
		auto parameters = model.generate_parameters();
		return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
	}

	template<
		template<template<typename> typename, typename> typename SpecializedModel,
		template<size_t, typename> typename PFC,
		size_t D, typename Sp0, typename... S>
	static constexpr auto with_new_solver(SpecializedModel<symphas::internal::MakeEquation, ModelPFCEquation<PFC, D, Sp0, S...>> const& model)
	{
		using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<S...>::template Specialized<symphas::internal::MakeEquation, ModelPFCEquation<PFC, D, Sp, S...>>;
		auto parameters = model.generate_parameters();
		return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
	}

	template<
		template<template<typename> typename, template<typename, typename...> typename, typename, typename> typename Specialized,
		template<typename> typename Eq,
		template<typename, typename...> typename Pr,
		size_t D, typename Sp0, typename... S, typename... P>
	static constexpr auto with_new_solver(
		Specialized<Eq, Pr, Pr<Model<D, Sp0, S...>, P...>, Eq<symphas::internal::MakeEquationProvisional<Pr<Model<D, Sp0, S...>, P...>>>> const& model)
	{
		using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<S...>::template ProvTypes<P...>::template Specialized<Eq, Pr>;
		auto parameters = model.generate_parameters();
		return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
	}


public:

	template<typename M>
	auto operator()(M const& model)
	{
		return with_new_solver(model);
	}
};




 /*!
  * \defgroup modelmacros Macros in Model Definitions
  * @{
  */

#define INVALID_MODEL -1
#define MAX_DEFINED_MODELS 64

 //! \cond

namespace symphas::internal
{

	/*!
	 * Primary object using recursive inheritance in order to
	 * define a member variable and return one of incremented value.
	 *
	 * A higher order index is given on the next usage by overloading
	 * the function for which the terminating point is instantiated. below
	 * it is never defined, because only the return type matters
	 */
	template<int N>
	struct model_count_index : model_count_index<N - 1>
	{
		static const int value = N + 1;
	};

	//! Base specialization to terminate the recursive inheritance.
	/*!
	 * Base specialization which terminates the recursive inheritance for
	 * incrementing model indices.
	 */
	template<>
	struct model_count_index<0>
	{
		static const int value = 1;
	};

	/* overloading this function on a new parameter and updating the return
	 * type will provide a new index on the next usage by decltype
	 */
	constexpr model_count_index<0> model_counter(model_count_index<0>);

}


//! Used in assigning unique names to models for indexing.
/*!
 * The naming format of a model index is defined. Each index name has 
 * to be different and requires a prefix.
 */ 
#define MODEL_INDEX_NAME(PREFIX_NAME) __ ## PREFIX_NAME ## _index

//! Iterates to the next model index for compile-time constant model indexing.
/*!
 * Convenience definition for setting the next index and providing
 * the PREFIX argument which names it.
 * Importantly, there cannot be more than #MAX_DEFINED_MODELS models defined 
 * because after that, the counter will no longer increment.
 */
#define NEXT_MODEL_INDEX(PREFIX_NAME) \
namespace symphas::internal { \
constexpr int MODEL_INDEX_NAME(PREFIX_NAME) = decltype(model_counter(model_count_index<MAX_DEFINED_MODELS>{}))::value; \
constexpr model_count_index<MODEL_INDEX_NAME(PREFIX_NAME)> model_counter(model_count_index<MODEL_INDEX_NAME(PREFIX_NAME)>); }


//! \endcond





#ifdef MODEL_APPLY_CALL

#define USING_MODEL_SELECTION

/* 
 * in order to link a defined model to a string, so that an arbitrary
 * function can be executed, a compile time index counter is maintained
 * 
 * the index counter is used to specialize a template class and define
 * a member function which acts as a selection and wrapper for providing
 * the model type to the user defined function
 */

template<int N>
struct model_call_wrapper
{
	template<template<typename, size_t> typename AppliedSolver, typename... Ts>
	static int call(size_t, StencilParams, const char*, Ts&& ...)
	{
		return INVALID_MODEL;
	}

	template<template<size_t> typename AppliedSolver, typename... Ts>
	static int call(size_t, const char*, Ts&& ...)
	{
		return INVALID_MODEL;
	}
};



// Calls a desired function with the model.
/*!
 * Sets up being able to call models with the ::model_select class. See this
 * class for more details.
 *
 * Select the model type using the solver and the model type, and the type name of the 
 * specialized model along with the string name are passed to a function defined in the
 * build process, so that the model can be run. If the solver uses FD stencils, parameters
 * are also passed so that the stencil coefficients are chosen. In this way
 * however, a model is instantiated for each combination of stencils.
 *
 * This requires that one definition is specified: #MODEL_APPLY_CALL.
 *
 * The definition #MODEL_APPLY_CALL, which is a function that needs to be
 * a template function of at least one type, the model type, which the selection
 * tree will then call if the name corresponds to the given model name by
 * checking all names until the corresponding model name is found, moreover, the given
 * function may only return an integer value which should not equal #INVALID_MODEL.
 *
 * The parameters to the #MODEL_APPLY_CALL method are given after the selection
 * parameters.
 *
 * At the individual function level, the macro definition assembles a
 * ternary tree in order to run the model on the correct solver type based on the
 * stencil and dimension.
 *
 * The given function #MODEL_APPLY_CALL may have multiple template parameters, but
 * the first must always be the model type, meaning that specifying the first parameter
 * may not necessarily fully qualify the method MODEL_APPLY_CALL.
 */
#define MODEL_WRAPPER_FUNC(NAME, GIVEN_NAME, MODEL, SOLVER) \
template<> \
struct model_call_wrapper<symphas::internal::MODEL_INDEX_NAME(NAME)> \
{ \
	template<template<typename, size_t> typename AppliedSolver, typename... Ts> \
	static int call(size_t dimension, StencilParams stp, const char* name, Ts&& ...args) \
	{ \
		if (std::strcmp(name, #GIVEN_NAME) == 0) { return ModelSelectStencil<MODEL, SOLVER>{ dimension, stp }(std::forward<Ts>(args)...); } \
		return model_call_wrapper<symphas::internal::MODEL_INDEX_NAME(NAME) - 1>::call<AppliedSolver>(dimension, stp, name, std::forward<Ts>(args)...); \
	} \
	template<template<size_t> typename AppliedSolver, typename... Ts> \
	static int call(size_t dimension, const char* name, Ts&& ...args) \
	{ \
		if (std::strcmp(name, #GIVEN_NAME) == 0) { return ModelSelect<MODEL, SOLVER>{ dimension }(std::forward<Ts>(args)...); } \
		return model_call_wrapper<symphas::internal::MODEL_INDEX_NAME(NAME) - 1>::call<AppliedSolver>(dimension, name, std::forward<Ts>(args)...); \
	} \
};



// **************************************************************************************

namespace symphas::internal
{

	template<
		template<size_t, typename> typename Model,
		template<typename, size_t> typename Solver,
		size_t N, size_t D, size_t O, size_t... Ps,
		typename... Ts
	>
	auto run_model_call(std::index_sequence<N, D, O, Ps...>, Ts&& ...args);

	template<
		template<size_t, typename> typename Model, template<size_t> typename Solver, 
		size_t N, size_t D, typename... Ts
	>
	auto run_model_call(std::index_sequence<N, D>, Ts&& ...args);

	template<template<size_t, typename> typename M, size_t D>
	struct allowed_model_dimensions { static const bool value = true; };

}

template<
	template<size_t, typename> typename Model, 
	template<typename, size_t = 0> typename Solver>
struct ModelSelectStencil
{

protected:


	template<typename>
	struct StencilFromSeq
	{
#if defined(ALL_STENCILS) && defined(GENERATE_UNDEFINED_STENCILS_ON)
		template<typename T0, typename... Ts>
		auto operator()(T0 const&, Ts&& ...)
		{
			fprintf(SYMPHAS_WARN, "the provided stencil point values are not implemented\n");
			return INVALID_MODEL;
		}

		template<size_t N, size_t D, size_t O, typename... Ts>
		auto operator()(std::index_sequence<N, D, O>, Ts&&... args)
		{
			fprintf(SYMPHAS_WARN, "all stencils will be auto-generated\n");
			return symphas::internal::run_model_call<Model, Solver>(std::index_sequence<N, D, O>{}, std::forward<Ts>(args)...);
		}
#else

		template<typename... Ts>
		auto operator()(Ts&& ...)
		{
			fprintf(SYMPHAS_WARN, "the provided stencil point values are not implemented\n");
			return INVALID_MODEL;
		}
#endif
	};

	template<size_t N, size_t D, size_t O, size_t... Ps>
	struct StencilFromSeq<std::index_sequence<N, D, O, Ps...>>
	{
		template<typename... Ts>
		auto operator()(Ts&& ...args)
		{
			return symphas::internal::run_model_call<Model, Solver>(std::index_sequence<N, D, O, Ps...>{}, std::forward<Ts>(args)...);
		}
	};

	size_t parameters[6];

	// Run model from stencil parameters.

	template<size_t N, size_t D, typename... Ts>
	auto search_ord(std::index_sequence<>, Ts&& ...) const
	{
		fprintf(SYMPHAS_WARN, "the provided order of accuracy '%zd' is invalid for constructing the model\n", parameters[2]);
		return INVALID_MODEL;
	}

	template<size_t N, size_t D, size_t O, size_t... Os, typename... Ts>
	auto search_ord(std::index_sequence<O, Os...>, Ts&& ...args) const
	{
		if (parameters[2] == O)
		{
			return symphas::lib::internal::search<StencilFromSeq>(
				parameters, symphas::internal::cross_list_t<N, D, O>{}, std::forward<Ts>(args)...);
		}
		else
		{
			return search_ord<N, D>(std::index_sequence<Os...>{}, std::forward<Ts>(args)...);
		}
	}

	template<size_t N, typename... Ts>
	auto search_dim(symphas::lib::types_list<>, Ts&& ...) const
	{
		fprintf(SYMPHAS_WARN, "the provided dimension value '%zd' is invalid for constructing the model\n", parameters[1]);
		return INVALID_MODEL;
	}

	template<size_t N, size_t D, size_t... Ns, typename... Seqs, typename... Ts>
	auto search_dim(symphas::lib::types_list<symphas::lib::types_list<std::index_sequence<D>, std::index_sequence<Ns...>>, Seqs...>, Ts&& ...args) const
	{
		if constexpr (symphas::internal::allowed_model_dimensions<Model, D>::value)
		{
			if (parameters[1] == D)
			{
				return search_ord<N, D>(std::index_sequence<Ns...>{}, std::forward<Ts>(args)...);
			}
		}
		return search_dim<N>(symphas::lib::types_list<Seqs...>{}, std::forward<Ts>(args)...);
	}

	template<typename... Ts>
	auto search_type(std::index_sequence<>, Ts&& ...) const
	{
		fprintf(SYMPHAS_WARN, "the provided solver variation '%zd' is invalid\n", parameters[0]);
		return INVALID_MODEL;
	}

	template<size_t N, size_t... Ns, typename... Seqs, typename... Ts>
	auto search_type(std::index_sequence<N, Ns...>, Ts&& ...args) const
	{
		if (parameters[0] == N)
		{
			return search_dim<N>(symphas::internal::dim_ord_list_t<AVAILABLE_DIMENSIONS>{}, std::forward<Ts>(args)...);
		}
		return search_type(std::index_sequence<Ns...>{}, std::forward<Ts>(args)...);
	}

public:

	ModelSelectStencil(size_t dimension, StencilParams stp)
		: parameters{ stp.type, dimension, stp.ord, stp.ptl, stp.ptg, stp.ptb } {}

	template<typename... Ts>
	auto operator()(Ts&& ...args) const
	{
		using namespace symphas::internal;
		return search_type(std::make_index_sequence<decltype(solver_counter(solver_count_index<solver_id_t<Solver<void>>, SOLVER_MAX_VARIATIONS>{}))::value - 1>{}, std::forward<Ts>(args)...);
	}
};

template<template<size_t, typename> typename Model, template<size_t N = 0> typename Solver>
struct ModelSelect
{
protected:

	template<size_t N, typename... Ts>
	auto search_dim(symphas::lib::types_list<>, Ts&& ...) const
	{
		fprintf(SYMPHAS_WARN, "the provided dimension value '%zd' is invalid for constructing the model\n", dimension);
		return INVALID_MODEL;
	}

	template<size_t N, size_t D, size_t... Ns, typename... Seqs, typename... Ts>
	auto search_dim(symphas::lib::types_list<symphas::lib::types_list<std::index_sequence<D>, std::index_sequence<Ns...>>, Seqs...>, Ts&& ...args) const
	{
		if constexpr (symphas::internal::allowed_model_dimensions<Model, D>::value)
		{
			if (dimension == D)
			{
				return symphas::internal::run_model_call<Model, Solver>(std::index_sequence<N, D>{}, std::forward<Ts>(args)...);
			}
		}
		return search_dim<N>(symphas::lib::types_list<Seqs...>{}, std::forward<Ts>(args)...);
	}

	template<typename... Ts>
	auto search_type(std::index_sequence<>, Ts&& ...) const
	{
		fprintf(SYMPHAS_WARN, "the provided solver variation '%zd' is invalid\n", type);
		return INVALID_MODEL;
	}

	template<size_t N, size_t... Ns, typename... Seqs, typename... Ts>
	auto search_type(std::index_sequence<N, Ns...>, Ts&& ...args) const
	{
		if (type == N)
		{
			return search_dim<N>(symphas::internal::dim_ord_list_t<AVAILABLE_DIMENSIONS>{}, std::forward<Ts>(args)...);
		}
		return search_type(std::index_sequence<Ns...>{}, std::forward<Ts>(args)...);
	}

	size_t dimension;
	size_t type;

public:

	ModelSelect(size_t dimension, size_t type) : dimension{ dimension }, type{ type } {}

	template<typename... Ts>
	auto operator()(Ts&& ...args) const
	{
		using namespace symphas::internal;
		return search_type(std::make_index_sequence<decltype(solver_counter(solver_count_index<solver_id_t<Solver<>>, SOLVER_MAX_VARIATIONS>{}))::value - 1>{}, std::forward<Ts>(args)...);
	}
};

#else

#endif


#define PARAMETERS_FROM_MODEL(PARAMETERIZED_TYPE, PARAMETERS) \
template<typename M> struct symphas::internal::model_field_parameters<M, PARAMETERIZED_TYPE> : \
	model_field_parameters<M, void>, parameterized_type<void, void> { \
	using parent_trait = model_field_parameters<M, void>; \
	using parent_type = parent_trait; \
	model_field_parameters(M const& m) : parent_type{ m }, parameterized_type<void, void>{ SINGLE_ARG PARAMETERS } {} \
}; \

#define PARAMETERIZED_TYPE(NAME, TYPE, PARAMETERS) \
namespace symphas::internal::parameterized { struct type_ ## NAME ## _t; } \
PARAMETERS_FROM_MODEL(symphas::internal::parameterized::type_ ## NAME ## _t, (SINGLE_ARG PARAMETERS)) \
template<> struct symphas::internal::parameterized_type<symphas::internal::parameterized::type_ ## NAME ## _t, symphas::internal::parameterized::TYPE> : \
	symphas::internal::parameterized_type<void, void> { \
	parameterized_type() : parameterized_type<void, void>() {} \
	template<typename M> parameterized_type(M const& model) : \
		parameterized_type<void, void>{ model_field_parameters<M, parameterized::type_ ## NAME ## _t>(model) } {} \
}; \
namespace symphas::internal { namespace parameterized { using NAME = symphas::internal::parameterized_type<type_ ## NAME ## _t, symphas::internal::parameterized::TYPE>; } } 





// **************************************************************************************

//! The beginning of a definition for a model of a phase field problem.
/*!
 * Begins the definition for a Model, an encapsulation of a phase field problem. The
 * implementation of a new phase field problem starts with this definition. The name
 * and types of the phase field order parameters are provided. The specification of
 * each new type is a new order parameter. It is important to note that _the order of 
 * the provided types defines the order of the variables specified in the equations of 
 * motion_.
 * 
 * Moreover, the name must not conflict with a previously used model name, it is a unique 
 * identifier. This is because the underlying implementation creates a namespace based on 
 * the identifier. Moreover, the identifier is also used in a new typedef which aliases the 
 * full model type. As a result, the identifier cannot have any spaces or punctuation, it
 * can only use characters as part of a valid C++ name, with the exception it may begin 
 * with a number.
 * 
 * The model type is referred to using the typedef created according to the format `model_<NAME>_t`.
 * For example a model named `1cp` is aliased by `model_1cp_t`. Importantly, this means that
 * the name is case sensitive.
 * 
 * The rest of the model definition is then provided as the third argument. 
 * 
 * The phase field definition supports the specification of the intermediate
 * variables. This is done using the definition:
 *
 * - #PROVISIONAL_DEF(...)
 * 
 * Intermediate variables are throughout SymPhas called provisional variables.
 * The elements associated with the intermediate variables in the solution are named
 * accordingly, including this definition. 
 * 
 * If there are no intermediate variables used in the model definition, then
 * this can be omitted. 
 * 
 * Immediately following these definitions, the equations of motion are provided using
 * the definition:
 * 
 * - #EVOLUTION(...)
 * 
 * \param NAME The name of the model.
 * \param TYPES The order parameter types of the phase field, given in a bracketed list.
 * The types **must** be provided in a list surrounded by brackets, i.e. `(SCALAR, SCALAR)`.
 * \param ... The rest of the phase field definition.
 */
#define MODEL(NAME, TYPES, ...) \
namespace model_ ## NAME { \
template<typename T, size_t D> struct allowed_model_dimensions { static const bool value = true; }; \
EQUATION_TRAIT_FORWARD_DECL \
using namespace symphas::internal; \
using namespace symphas::internal::parameterized; \
template<size_t Dm, typename Sp> \
using OpTypes = symphas::internal::expand_types_to_model_t<Dm, Sp, SINGLE_ARG TYPES>; \
template<typename> struct using_provisional { template<size_t Dm, typename Sp> using type = typename OpTypes<Dm, Sp>::template Specialized<TraitEquationModel>; }; \
	__VA_ARGS__ } \
template<size_t Dm, typename Sp> \
using model_ ## NAME ## _t = model_ ## NAME ::SpecializedModel<Dm, Sp>; \
template<size_t D> \
struct symphas::internal::allowed_model_dimensions<model_ ## NAME ## _t, D> \
{ static const bool value = model_ ## NAME::allowed_model_dimensions<void, D>::value; };


//! Definition for defining the intermediate variables.
/*!
 * Intermediate variables used in a phase field problem are called provisional variables.
 * The elements associated with the intermediate variables in the solution are named
 * accordingly, including this definition. This specifies the equations of the provisional
 * variables used in the phase field problem specification.
 * 
 * Provisional variables are also grids, and have their dimensions inherited 
 * from the phase field systems. Provisional variables are useful to precompute
 * complicated expressions before it is evaluated in an equation of motion.
 * Since they do not directly participate in the EVOLUTION of the
 * system, they can also be used as virtual variables. Provisional
 * variables are always evaluated before the dynamical equations are evaluated.
 * 
 * Some important remarks about using provisional variables in general:
 * 
 * - Provisional variables can be defined in terms of other provisional
 * variables, but they are always *sequentially evaluated*. This means that if
 * the first provisional variable is defined in terms of the second one, 
 * the first provisional variable may be evaluated to undefined values. One
 * particular implication of this is that values of previous iterations may be
 * saved in provisional variables. For instance, if we defined provisional
 * variables `A` and `B` such that `B` copies the phase field, and `A` copies
 * `B`, then `B` will always be equal to the phase field data of the previous 
 * solution index, and `A` will be equal to the phase field data of two
 * solution indices previous, except on the first iteration where it will be
 * zero because provisional variables are initialized to zero.
 * 
 * - The derivative of a provisional variable can be taken, but this should
 * typically not be done. While the specific problems that can arise depend
 * on the solver, take the example of differentiating the provisional variable,
 * and then differentiating it again in the equation of motion. The provisional
 * variable would be completely evaluated, using the solver to approximate
 * the derivative to obtain its values, but differentiating it again would
 * differentiate the approximation, thereby spoiling the numerical approximation
 * of the actual equations of motion. 
 * 
 * 
 * This definition is used to specify the types of the intermediate variables and
 * their definitions. The definition is provided just like a dynamical equation,
 * but here, the form is `var(N) = <equation>` where `N` is the index of the
 * provisional variable from the 1-indexed list and `var` is the macro
 * used to refer to a provisional variable.
 * 
 * > **NOTE**: The provisional variable indexing begins at 1, so that `var(1)`
 * > refers to the first provisional variable.
 * 
 * \param TYPES The order parameter types of the intermediate variables, given in 
 * a bracketed list. As for the primary model definition, the types **must** be 
 * provided in a list surrounded by brackets, i.e. `(SCALAR, SCALAR)`.
 * \param ... The equations defining computation of the intermediate variables.
 * The order of these equations is determined by the order the types are provided.
 */
#define PROVISIONAL_DEF(TYPES, ...) \
PROVISIONAL_TRAIT_FORWARD_DECL \
template<> struct using_provisional<void> { \
template<size_t Dm, typename Sp> using type = \
	typename OpTypes<Dm, Sp>::template ProvTypes<SINGLE_ARG TYPES>::template Specialized<TraitEquationModel, TraitProvisionalModel>; }; \
PROVISIONAL_TRAIT_DEFINITION(__VA_ARGS__)

//! Definition for providing phase field equations of motion.
/*!
 * The phase field equations are specified by providing equations in a comma
 * delimited list. Each equation is specified with the format 
 * `dop(N) = <equation>`, where `N` is the index of the phase field for which
 * the dynamical equation is being defined, `dop` represents the time
 * derivative of the order parameter of this field, and `<equation>` is the
 * equation of the dynamical equation written in any way. The term `dop` can
 * only appear on the left hand side of an equation.
 * 
 * > **NOTE**: The phase fields are indexed beginning with 1, meaning that
 * > `dop(1)` refers to the first phase field.
 *
 * The order parameter of the `N`-th field is referred to as `op(N)` in the
 * equation.
 * 
 * *A short explanation of some symbolic algebra features available in
 * writing equations of motion:*
 * 
 * The following derivatives of expressions are
 * supported:
 * - #grad(E)
 * - #lap(E)
 * - #gradlap(E)
 * - #bilap(E)
 * 
 * where `E` is the expression that is being differentiated.
 * 
 * Derivatives of higher orders can also be specified using `dx(O)`, where `O`
 * is the desired order. The usage of `dx` is different from the derivatives
 * above, as in order to apply it, it must be multiplied by an expression. This
 * applies the derivative to the expression using symbolic algebra rules.
 * 
 * When `O`\f$\le 4\f$, then it is equivalent to the 
 * functions above when multiplied by an expression. 
 * 
 * A unique aspect of using `dx` is that it can be added to another 
 * `dx` of a different order to create a new operator. I.e.:
 * 
 * ```cpp
 * dop(1) = (dx(2) + dx(4)) * op(1)
 * ```
 * results in:
 * 
 * \f{align*}
 * \frac{d\psi}{dt} &= \left(\nabla^2 + \nabla^4\right)\psi \\
 * &= \nabla^2 \psi + \nabla^4 \psi\,.
 * \f}
 * 
 * To use the model coefficients, there is the macro `param(N)`, where `N` 
 * refers to the `N`-th coefficient in the that is given to the model, as a
 * 1-indexed list. That is, the first coefficient is `param(1)`.
 * 
 * Constants are used with `lit(v)`, where `v` is the value of the constant.
 * 
 * Imaginary numbers are also supported, used with `im`, and support:
 * - `modulus(Z)`, the modulus of the number.
 * - `conj(Z)`, the conjugate of the number.
 * - `Re(Z)`, the real part of the complex number.
 * - `Im(Z)`, the imaginary part of the complex number.
 * 
 * where `Z` is the complex number. This usage applies to the order parameter
 * as well, not just coefficients.
 * 
 * It is also possible to use the Gaussian kernel in the equations of motion,
 * which is a grid of the Fourier wavenumbers. See #Gaussian and 
 * OpConvolution.
 * 
 * Also, where the number 1 (i.e. the multiplicative identity) should be used
 * in an expression explicitly, it is recommended to use #one, to maximize
 * the advantage of the compile-time symbolic algebra.
 * 
 * It may be desired to write some preliminary expressions or work that can
 * be used to write the equations of motion. In this case,
 * a preamble section can be specified. The preamble is a code block that 
 * immediately precedes the line that processes the equations of motion.
 * Therefore, everything within the preamble is only available to the equations 
 * of motion.
 * 
 * If no preamble is required and only the equations of motion need to be
 * specified, #EVOLUTION should be used. 
 * 
 * \param PREAMBLE Setup work in order to write the equations of motion. If
 * this part is not required, use the macro #EVOLUTION.
 * \param ... The list of the equations of motion.
 */
#define EVOLUTION_PREAMBLE(PREAMBLE, ...) \
template<size_t Dm, typename Sp> \
using SpecializedModel = typename using_provisional<void>::template type<Dm, Sp>; \
EQUATION_TRAIT_PREAMBLE(SINGLE_ARG PREAMBLE) \
EQUATION_TRAIT_DEFINITION(__VA_ARGS__)

//! See #EVOLUTION_PREAMBLE.
/*!
 * Shortcut macro for when no preamble is needed. Refer to #EVOLUTION_PREAMBLE
 * for full information on defining the equations of motion of a model.
 */
#define EVOLUTION(...) \
EVOLUTION_PREAMBLE((), __VA_ARGS__)

#define FREE_ENERGY_PREAMBLE(PREAMBLE, SELECTED_DYNAMICS, ...) \
template<size_t Dm, typename Sp> \
using SpecializedModel = typename using_provisional<void>::template type<Dm, Sp>; \
EQUATION_TRAIT_PREAMBLE(SINGLE_ARG PREAMBLE) \
EQUATION_FE_TRAIT_DEFINITION((SINGLE_ARG SELECTED_DYNAMICS), __VA_ARGS__)
//EVOLUTION_PREAMBLE(SINGLE_ARG PREAMBLE, parent_type::template generate_equations(std::make_tuple(SINGLE_ARG SELECTED_DYNAMICS), __VA_ARGS__))

#define FREE_ENERGY(SELECTED_DYNAMICS, ...) \
FREE_ENERGY_PREAMBLE((), (SINGLE_ARG SELECTED_DYNAMICS), __VA_ARGS__)


#ifdef USING_MODEL_SELECTION

//! Associate model definition with string.
/*!
 * Used for associating a model with a string when building the model selection tree.
 * The specified model name is linked with the provided string, which will be compared with
 * a selection string when searching for correct model. 
 * 
 * Linking models is entirely optional. For instance, it may be desired to link only some 
 * of the defined models. Any number of models can be linked, and the same model can be 
 * linked multiple times to provide additional names.
 * 
 * \param NAME The name of the model to link.
 * \param GIVEN_NAME The name given to the model when searching by string. This parameter
 * is not passed as a string because it is stringified by the macro.
 */
#define LINK_WITH_NAME(NAME, GIVEN_NAME) \
NEXT_MODEL_INDEX(NAME) \
MODEL_WRAPPER_FUNC(NAME, GIVEN_NAME, model_ ## NAME::SpecializedModel, AppliedSolver)

#else

//! Without model selection, linking is an empty call.
#define LINK_WITH_NAME(NAME, GIVEN_NAME)

#endif

// **************************************************************************************



// **************************************************************************************

//! The beginning of a definition for a Specialized Model.
/*!
 * Begins the definition for a Model of the phase field crystal variation, an encapsulation 
 * of a phase field crystal problem. The implementation of a new phase field problem 
 * starts with this definition. The name 
 * and types of the phase field order parameters are provided. The specification of 
 * each new type is a new order parameter. It is important to note that _the order of 
 * the provided types defines the order of the variables specified in the equations of 
 * motion_.
 *
 * Moreover, the name must not conflict with a previously used model name, it is a unique
 * identifier. The names from the PFC definitions will never conflict with the names
 * of the usual Model definitions. Moreover, the identifier is also used in a new typedef which aliases the
 * full model type. As a result, the identifier cannot have any spaces or punctuation, it
 * can only use characters as part of a valid C++ name, with the exception it may begin
 * with a number. 
 *
 * The model type is referred to using the typedef created according to the 
 * format `modelpfc_NAME_t`. For example a model named `1conserved` is aliased 
 * by `modelpfc_1conserved_t`. Importantly, this means that 
 * the name is case sensitive.
 *
 * The rest of the model definition is then provided as the third argument.
 *
 * The rest of the phase field crystal definition requires specifying first the
 * default behaviour of the model, followed by the order parameter types,
 * followed then by specialized behaviour of particular order parameters. Specifically,
 * this entails by first providing:
 *
 * - #DYNAMICS(...) or #DEFAULT_MODE(...); or
 * - #NO_DEFAULTS (if no default is changed)
 *
 * Immediately following these definitions, the types of the order parameters
 * are provided in a bracketd list. Finally, the behaviour of individual phase
 * fields may be modified using:
 *
 * - #MODE_OF(...); or
 * - #DYNAMICS_OF(...)
 *
 * \param NAME The name of the model.
 * \param DEFAULTS Specifies a different set of default behaviours to modify the phase 
 * field crystal model that is originally defined. See #DEFAULTS(...).
 * \param TYPES The order parameter types of the phase field, given in a bracketed list.
 * The types **must** be provided in a list surrounded by brackets, i.e. `(SCALAR, SCALAR)`.
 * \param ... The rest of the phase field crystal definition.
 */
#define PFC_TYPE(NAME, DEFAULTS, TYPES, ...) \
namespace modelpfc_ ## NAME { \
using namespace symphas::internal::parameterized; \
struct PFCParameters : PFCParametersDefault<modelpfc_ ## NAME::PFCParameters> \
{ \
	using parent_type = PFCParametersDefault<modelpfc_ ## NAME::PFCParameters>; \
	template<size_t N> static constexpr symphas::internal::DynamicType dynamic_val_apply() { return parent_type::template dynamic_val_apply<N>(); } \
	template<size_t N> static constexpr size_t mode_val_apply() { return parent_type::template mode_val_apply<N>(); } \
	DEFAULTS \
}; \
template<size_t Dm, typename Sp> \
struct ModelPFCApplied : ModelPFCEquation<modelpfc_ ## NAME::ModelPFCApplied, Dm, Sp, SINGLE_ARG TYPES>, PFCParameters \
{ \
	using parent_type = ModelPFCEquation<modelpfc_ ## NAME::ModelPFCApplied, Dm, Sp, SINGLE_ARG TYPES>; \
	using parent_type::parent_type; \
}; \
template<size_t Dm, typename Sp> \
using ModelPFCSpecialized = GeneralizedPFCModel<ModelPFCApplied, Dm, Sp, SINGLE_ARG TYPES>; \
__VA_ARGS__ } \
template<size_t Dm, typename Sp> \
using model_ ## NAME ## _t = modelpfc_ ## NAME ::ModelPFCSpecialized<Dm, Sp>;


//! Specify different default behaviours.
/*!
 * Provide a list of the different default behaviours of the phase field crystal
 * model, which will apply to all fields.
 */
#define DEFAULTS(...) __VA_ARGS__
#define NO_DEFAULTS

//! Sets a given dynamic type to be the default.
/*!
 * Changes the EVOLUTION of all the fields, unless the field is
 * explicitly changed in the definition in #DYNAMICS_OF.
 * 
 * Used with the #DEFAULTS definition.
 * 
 * \param TYPE The dynamic type which dictates the default behaviour of
 * the phase field crystal problem.
 */
#define DYNAMICS(TYPE) \
const static symphas::internal::DynamicType DEFAULT_DYNAMIC = decltype(TYPE)::value;

//! Sets the default mode of the phase field crystal problem.
/*!
 * The mode refers roughly to the how the densities are distributed. Typically
 * the mode is only of order 2, resulting in a tetrahedral pattern.
 * 
 * Changes the order of the mode for all fields, unless explicitly changed
 * in #MODE_OF. The mode should always be greater than or equal to 2.
 * 
 * Used with the #DEFAULTS definition.
 * 
 * \param N The numeric value of the phase field crystal mode.
 */
#define MODE(N) \
const static size_t DEFAULT_MODE_N = N;

//! Indicates that a given field follows conserved EVOLUTION.
#define CONSERVED symphas::internal::dynamics_key_t<symphas::internal::DynamicType::CONSERVED>{}

//! Indicates that a given field follows non-conserved EVOLUTION.
#define NONCONSERVED symphas::internal::dynamics_key_t<symphas::internal::DynamicType::NONCONSERVED>{}


//! Modify the EVOLUTION of the given field.
/*!
 * Specializes a particular field with the given non-default behaviour 
 * (although a default characteristic can still be provided with undesired 
 * no side effects).
 * 
 * \param FIELD_N The field index which is specialized.
 * \param TYPE The dynamic type which dictates the new EVOLUTION of the field.
 */
#define DYNAMICS_OF(FIELD_N, TYPE) template<> constexpr symphas::internal::DynamicType PFCParameters::dynamic_val_apply<FIELD_N - 1>() { return decltype(TYPE)::value; }

//! Modify the mode of a field.
/*!
 * The mode refers roughly to the how the densities are distributed. Typically
 * the mode is only of order 2, resulting in a tetrahedral pattern.
 * 
 * Changes the order of the mode, which should always be greater than or equal
 * to 2.
 * 
 * \param FIELD_N The field index which is specialized.
 * \param N The new mode value of the field.
 */
#define MODE_OF(FIELD_N, N) template<> constexpr size_t PFCParameters::mode_val_apply<FIELD_N - 1>() { return N; }

#define EQUATION_OF(FIELD_N) symphas::internal::special_dynamics_select<decltype(FIELD_N)>{}.select<int(FIELD_N)>(FIELD_N)

#define ALL_CONSERVED(I) EQUATION_OF(I)(CONSERVED)
#define ALL_NONCONSERVED(I) EQUATION_OF(I)(NONCONSERVED)


#ifdef USING_MODEL_SELECTION

//! Associate model definition with string.
/*!
 * See #LINK_WITH_NAME.
 */
#define LINK_PFC_WITH_NAME(NAME, GIVEN_NAME) \
NEXT_MODEL_INDEX(NAME) \
MODEL_WRAPPER_FUNC(NAME, GIVEN_NAME, modelpfc_ ## NAME::ModelPFCSpecialized, AppliedSolver)

#else

//! Without model selection, linking is an empty call.
#define LINK_PFC_WITH_NAME(NAME, GIVEN_NAME)

#endif



// **************************************************************************************



#define RESTRICT_DIMENSIONS(...) \
template<size_t D> \
struct allowed_model_dimensions<void, D> \
{ static const bool value = symphas::lib::is_value_in_seq<size_t, D, std::index_sequence<__VA_ARGS__>>::value; };

// **************************************************************************************

namespace symphas::internal
{
	using namespace parameterized;

	template<int N, typename T>
	struct model_repeating_type
	{
		using type = list_repeating_type_t<N, T>;
	};

	template<typename T>
	struct model_repeating_type<CONFIGURATION, T>
	{
		using type = symphas::internal::field_array_t<T>;
	};

	template<int N, typename T>
	using model_repeating_type_t = typename model_repeating_type<N, T>::type;
}

//! Multiple phase fields of the same type.
/*!
 * The given order parameter of `TYPE` will be repeated `N` times.
 */
#define MANY(TYPE, N) symphas::internal::model_repeating_type_t<N, TYPE>

 //! Multiple real valued order parameter.
 /*!
  * `N` order parameters will be of type ::scalar_t. Using this definition
  * will repeat the #SCALAR type `N` times, thereby creating `N` order parameters
  * of this type.
  */
#define SCALARS(N) MANY(SCALAR, N)

 //! Multiple real valued order parameter.
 /*!
  * `N` order parameters will be of type ::complex_t. Using this definition
  * will repeat the #COMPLEX type `N` times, thereby creating `N` order parameters
  * of this type.
  */
#define COMPLEXES(N) MANY(COMPLEX, N)

 //! Multiple real valued order parameter.
 /*!
  * `N` order parameters will be of type ::vector_t. Using this definition
  * will repeat the #VECTOR type `N` times, thereby creating `N` order parameters
  * of this type.
  */
#define VECTORS(N) MANY(VECTOR, N)

#define INT INT

// **************************************************************************************

//! The `i`-th order parameter.
/*!
 * The variable which represents the `i`-th order parameter. This is used
 * in an equation of motion. It is indexed from 1, meaning that `op(1)` is the
 * first order parameter.
 */
#define op(N) parent_type::template op<N - 1>()

//! The time derivative of the `i`-th order parameter.
/*!
 * The time derivative of the `i`-th parameter, which always appears on the
 * left hand side of the equals sign in an equation of motion.
 * It is indexed from 1, meaning that `dop(1)` is the derivative of the
 * first order parameter.
 */
#define dop(N) parent_type::template dop<N - 1>()

//! The `i`-th provisional variable.
/*!
 * The `i`-th provisional variable. It is indexed from 1, 
 * meaning that `var(1)` is the first provisional variable.
 */
#define var(N) parent_type::template var<N - 1>()

//! The differential operator of order `O`.
/*!
 * The differential operator of order `O`, which is applied to expressions
 * through multiplication.
 * 
 * \param ORDER The order of the derivative.
 */
#define Diff(ORDER) expr::make_operator_derivative<ORDER>(solver)

 //! The directional derivative along x of first order.
 /*!
  * The first order directional derivative along x, which is applied to 
  * expressions through multiplication.
  */
#define gradx expr::make_operator_directional_derivative<Axis::X, 1>(solver)

//! The directional derivative along x of first order.
/*!
 * The first order directional derivative along x, which is applied to
 * expressions through multiplication.
 */
#define grady expr::make_operator_directional_derivative<Axis::Y, 1>(solver)

//! The directional derivative along x of first order.
/*!
 * The first order directional derivative along x, which is applied to
 * expressions through multiplication.
 */
#define gradz expr::make_operator_directional_derivative<Axis::Z, 1>(solver)


//! The laplacian of an expression.
/*!
 * The laplacian is the derivative of 2nd order. Applies the laplacian to 
 * the given expression.
 * 
 * \param EXPR The expression which is differentiated.
 */
#define lap Diff(2)

//! The laplacian of an expression.
/*!
 * The bilaplacian is the derivative of 4th order. Applies the bilaplacian to
 * the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define bilap Diff(4)

 //! The gradient of an expression.
 /*!
  * The gradient applies the first derivative to all spatial components of the
  * given expression.
  *
  * \param EXPR The expression which is differentiated.
  */
#define grad Diff(1)

//! The laplacian of an expression.
/*!
 * The gradient of the laplacian. Apply the gradlaplacian to the given 
 * expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define gradlap Diff(3)

//! The `N`-th coefficient of the model.
/*!
 * Refer to the value of the `N`-th parameter which was passed to the
 * phase field model. It is indexed from 1, 
 * meaning that `param(1)` is the first coefficient.
 * 
 * \param N The `N`-th parameter/coefficient.
 */
#define param(N) parent_type::param(N - 1)

//! Imposes a matrix shape with `N` columns to the coefficient array.
/*!
 * Allows accessing coefficients with an index, either with a number
 * or a series index. This is generally for associating each field
 * with a sequence of coefficients, which can change for each field.
 * 
 * \param N The number of columns, i.e. number of unique coefficients, 
 * associated with each field.
 */
#define param_matrix(N) parent_type::template param_matrix<N>()

#define diff_N(E, VAR, N) expr::make_symbolic_derivative<N>(E, VAR)
#define diff(E, VAR) diff_N(E, VAR, 1)

#define NOISE(VARIETY, TYPE, ...) parent_type::template make_noise<expr::NoiseType:: VARIETY, TYPE>(__VA_ARGS__)
#define WHITE_NOISE(TYPE, ...) NOISE(WHITE, TYPE, __VA_ARGS__)
#define POISSON_NOISE(INTENSITY, LAMBDA, ...) NOISE(POISSON, SCALAR, INTENSITY, LAMBDA, __VA_ARGS__)

#define _nW(TYPE, ...) WHITE_NOISE(TYPE, __VA_ARGS__)
#define _nP(INTENSITY, LAMBDA, ...) POISSON_NOISE(INTENSITY, LAMBDA, __VA_ARGS__)

// **************************************************************************************


//! Get the area of the system.
#define AREA expr::make_literal(parent_type::template system<0>().get_info().area())

//! Get the length of the interval along the given axis.
#define LENGTH(AXIS) expr::make_literal(parent_type::template system<0>().get_info()[Axis::AXIS].length())

//! Get the total number of fields of the model.
#define NUM_FIELDS symphas::internal::num_fields_as_literal(*this)


//! Get the dimensions of the system.
#define DIMENSIONS_OF(N) parent_type::template system<N>().dims

//! Get the intervals of the axis in the system.
#define INTERVALS_OF(N, AXIS) parent_type::template system<N>().get_info().at(AXIS)

//! Get the intervals of the system.
#define INTERVALS(N) parent_type::template system<N>().get_info().intervals


// **************************************************************************************


//! Standard Gaussian with spatial distance equal to unity.
#define Gaussian GaussianSmoothing<Dm>{ \
	parent_type::template system<0>().dims, 1.0, 1.0 }

//! Gaussian with spatial width based on the i-th field.
#define Gaussian_n(i) GaussianSmoothing<Dm>{ \
	parent_type::template system<i>().dims, \
	parent_type::template system<i>().get_info().get_widths().get(), \
	1.0, 1.0 }

//! Apply smoothing to the expression `Z`.
#define smoothing(Z) expr::make_convolution(Gaussian, Z)

//! Smoothing with Gaussian based on the spatial width of the i-th field.
#define smoothing_n(i, Z) expr::make_convolution(Gaussian_F(i), Z)

//! Use an integer in an expression as a constant value.
/*!
 * The integer is generated as a compile-time constant value.
 */
#define integer(N) expr::make_integer<N>()


//! The unit vector, which can be defined with one or two (for 3D) angles.
/*!
 * The correct unit vector will be chosen according to the dimension. In two dimensions,
 * only the first direction will be used, and in three, the second direction will be used.
 * If only one direction is provided, then in 3D, the second direction is assumed to be 0.
 */
#define e(...) expr::make_unit_vector<Dm>(__VA_ARGS__)

//! Apply the power N to an expression E.
/*!
 * Construct an expression representing an expression E to the power of the given value, N.
 */
#define power(E, N) expr::pow<N>(E)



#define c1 param(1)		 //!< Coefficient at index 1 in the list of coefficients.
#define c2 param(2)		 //!< Coefficient at index 2 in the list of coefficients.
#define c3 param(3)		 //!< Coefficient at index 3 in the list of coefficients.
#define c4 param(4)		 //!< Coefficient at index 4 in the list of coefficients.
#define c5 param(5)		 //!< Coefficient at index 5 in the list of coefficients.
#define c6 param(6)		 //!< Coefficient at index 6 in the list of coefficients.
#define c7 param(7)		 //!< Coefficient at index 7 in the list of coefficients.
#define c8 param(8)		 //!< Coefficient at index 8 in the list of coefficients.
#define c9  param(9)	 //!< Coefficient at index 9 in the list of coefficients.
#define c10 param(10)	 //!< Coefficient at index 10 in the list of coefficients.
#define c11 param(11)	 //!< Coefficient at index 11 in the list of coefficients.
#define c12 param(12)	 //!< Coefficient at index 12 in the list of coefficients.
#define c13 param(13)	 //!< Coefficient at index 13 in the list of coefficients.
#define c14 param(14)	 //!< Coefficient at index 14 in the list of coefficients.
#define c15 param(15)	 //!< Coefficient at index 15 in the list of coefficients.
#define c16 param(16)	 //!< Coefficient at index 16 in the list of coefficients.

#define C(N) param_matrix(N)	//!< Initialization of a coefficient matrix.

// free energy parameters

#define FE_FUNCTION(NAME, ...) NAME(solver, __VA_ARGS__)
#define LANDAU_FE(...) FE_FUNCTION(landau_fe, __VA_ARGS__)
#define DOUBLE_WELL_FE(...) FE_FUNCTION(doublewell_fe, __VA_ARGS__)
#define CELLULAR_FE(...) FE_FUNCTION(cellular_fe, __VA_ARGS__)


// **************************************************************************************

#define SUM(...) symphas::internal::series_index_selection(names_t{}, parent_type::systems_tuple(), __VA_ARGS__)
#define op_(In, P) parent_type::op_n(In)
#define index_(I) (i_<I, 0>{})

#define op_ii op_(ii, 0)
#define op_jj op_(jj, 0)
#define op_kk op_(kk, 0)

//#define SUM_INDEX(I, P) (i_<I, P>{})

//#define ii_(P) SUM_INDEX(0, P)
//#define ii ii_(0)
//#define jj_(P) SUM_INDEX(1, P)
//#define jj jj_(0)
//#define kk_(P) SUM_INDEX(2, P)
//#define kk kk_(0)

//#define op_(In, P) SUM_VARIABLE(In, P)
//#define op_ii_(P) op_(ii_(P), P)
//#define op_jj_(P) op_(jj_(P), P)
//#define op_jj op_jj_(0)
//#define op_kk_(P) op_(kk_(P), P)
//#define op_kk op_kk_(0)


#define DF(N) symphas::internal::dFE_var<N - 1>()
#define DF_(In) symphas::internal::dFE_var(In)
#define DF_ii DF_(ii_(0))
#define DF_jj DF_(jj_(0))

 // **************************************************************************************

 //! @}



