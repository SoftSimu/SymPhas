
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

#include "modelspecialized.h"
#include "modelpfc.h"
#include "stencilincludes.h"
#include "expressionfunctions.h"



 /*!
  * \defgroup modelmacros Macros in Model Definitions
  * @{
  */

#define INVALID_MODEL -1
#define MAX_DEFINED_MODELS 256

 //! \cond

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
constexpr int MODEL_INDEX_NAME(PREFIX_NAME) = decltype(model_counter(model_count_index<MAX_DEFINED_MODELS>{}))::value; \
constexpr model_count_index<MODEL_INDEX_NAME(PREFIX_NAME)> model_counter(model_count_index<MODEL_INDEX_NAME(PREFIX_NAME)>);


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

//! \cond
template<int N>
struct model_call_wrapper
{
	template<template<typename> typename AppliedSolver, typename... Ts>
	static int call(size_t, StencilParams, const char*, Ts&& ...)
	{
		return INVALID_MODEL;
	}
};


//! Implements the association between a model and a name.
#define MODEL_WRAPPER_FUNC(NAME, IMPL) \
template<> \
struct model_call_wrapper<MODEL_INDEX_NAME(NAME)> \
{ \
	template<template<typename> typename AppliedSolver, typename... Ts> \
	static int call(size_t dimension, StencilParams stp, const char* name, Ts&& ...args) \
	{ \
		IMPL \
		return model_call_wrapper<MODEL_INDEX_NAME(NAME) - 1>::call<AppliedSolver>(dimension, stp, name, std::forward<Ts>(args)...); \
	} \
};
//! \endcond


// Calls a desired function with the model.
/*!
 * Sets up being able to call models with the ::model_select class. See this
 * class for more details.
 * 
 * Select the model type using the first three provided parameters,
 * and the remaining generic parameters are passed to a function defined in the
 * build process, so that executing the model. The parameters apply to the 
 * finite different stencils; the stencil sizes are chosen. In this way
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
#define RUNMODELNNN(SOLVER, MODEL, Ll, Bb, Gg, DIM, ORD) \
MODEL_APPLY_CALL <MODEL<DIM, SOLVER<Stencil ## DIM ## d ## ORD ## h<Ll, Bb, Gg>>>> (std::forward<Ts>(args)...)

#if !defined(DEBUG) && defined(ALL_STENCILS)

  // 2d 2h
#define RUNMODEL2D2HNN(SOLVER, MODEL, L, B) ( \
(stp.ptg == 6) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 6, 2, 2)) : \
(stp.ptg == 8) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 8, 2, 2)) : \
(stp.ptg == 12) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 12, 2, 2)) : \
(stp.ptg == 16) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 16, 2, 2)) : INVALID_MODEL ) 

#define RUNMODEL2D2HN(SOLVER, MODEL, L) ( \
(stp.ptb == 13) ? (RUNMODEL2D2HNN(SOLVER, MODEL, L, 13)) : \
(stp.ptb == 17) ? (RUNMODEL2D2HNN(SOLVER, MODEL, L, 17)) : \
(stp.ptb == 21) ? (RUNMODEL2D2HNN(SOLVER, MODEL, L, 21)) : INVALID_MODEL )

#define RUNMODEL2D2H(SOLVER, MODEL) ( \
(stp.ptl == 5) ? (RUNMODEL2D2HN(SOLVER, MODEL, 5)) : \
(stp.ptl == 9) ? (RUNMODEL2D2HN(SOLVER, MODEL, 9)) : 0 )

// 3d 2h
#define RUNMODEL3D2HNN(SOLVER, MODEL, L, B) ( \
(stp.ptg == 10) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 10, 3, 2)) : \
(stp.ptg == 12) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 12, 3, 2)) : \
(stp.ptg == 28) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 28, 3, 2)) : \
(stp.ptg == 36) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 36, 3, 2)) : \
(stp.ptg == 40) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 40, 3, 2)) : INVALID_MODEL )

#define RUNMODEL3D2HN(SOLVER, MODEL, L) ( \
(stp.ptb == 21) ? (RUNMODEL3D2HNN(SOLVER, MODEL, L, 21)) : \
(stp.ptb == 25) ? (RUNMODEL3D2HNN(SOLVER, MODEL, L, 25)) : \
(stp.ptb == 41) ? (RUNMODEL3D2HNN(SOLVER, MODEL, L, 41)) : \
(stp.ptb == 52) ? (RUNMODEL3D2HNN(SOLVER, MODEL, L, 52)) : \
(stp.ptb == 57) ? (RUNMODEL3D2HNN(SOLVER, MODEL, L, 57)) : INVALID_MODEL )

#define RUNMODEL3D2H(SOLVER, MODEL) ( \
(stp.ptl == 7) ? (RUNMODEL3D2HN(SOLVER, MODEL, 7)) : \
(stp.ptl == 15) ? (RUNMODEL3D2HN(SOLVER, MODEL, 15)) : \
(stp.ptl == 19) ? (RUNMODEL3D2HN(SOLVER, MODEL, 19)) : \
(stp.ptl == 21) ? (RUNMODEL3D2HN(SOLVER, MODEL, 21)) : \
(stp.ptl == 27) ? (RUNMODEL3D2HN(SOLVER, MODEL, 27)) : INVALID_MODEL )

// 2d 4h
#define RUNMODEL2D4HNN(SOLVER, MODEL, L, B) ( \
(stp.ptg == 14) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 14, 2, 4)) : \
(stp.ptg == 18) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 18, 2, 4)) : \
(stp.ptg == 26) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 26, 2, 4)) : \
(stp.ptg == 30) ? (RUNMODELNNN(SOLVER, MODEL, L, B, 30, 2, 4)) : INVALID_MODEL ) 

#define RUNMODEL2D4HN(SOLVER, MODEL, L) ( \
(stp.ptb == 21) ? (RUNMODEL2D4HNN(SOLVER, MODEL, L, 21)) : \
(stp.ptb == 25) ? (RUNMODEL2D4HNN(SOLVER, MODEL, L, 25)) : \
(stp.ptb == 33) ? (RUNMODEL2D4HNN(SOLVER, MODEL, L, 33)) : \
(stp.ptb == 37) ? (RUNMODEL2D4HNN(SOLVER, MODEL, L, 37)) : INVALID_MODEL )

#define RUNMODEL2D4H(SOLVER, MODEL) ( \
(stp.ptl == 9) ? (RUNMODEL2D4HN(SOLVER, MODEL, 9)) : \
(stp.ptl == 17) ? (RUNMODEL2D4HN(SOLVER, MODEL, 17)) : \
(stp.ptl == 21) ? (RUNMODEL2D4HN(SOLVER, MODEL, 21)) : INVALID_MODEL )

// 4h
#define RUNMODEL4H(SOLVER, MODEL_NAMESPACE, MODEL) ( \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 2>::value) \
	if (dimension == 2) return RUNMODEL2D4H(SOLVER, MODEL); \
return INVALID_MODEL;

// 2h
#define RUNMODEL2H(SOLVER, MODEL_NAMESPACE, MODEL) ( \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 2>::value) \
	if (dimension == 2) return RUNMODEL2D2H(SOLVER, MODEL); \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 3>::value) \
	if (dimension == 3) return RUNMODEL3D2H(SOLVER, MODEL); \
return INVALID_MODEL;

// full
#define RUNMODEL(SOLVER, MODEL_NAMESPACE, MODEL) ( \
if (stp.ord == 2) (RUNMODEL2H(SOLVER, MODEL_NAMESPACE, MODEL)); \
else if (stp.ord == 4) (RUNMODEL4H(SOLVER, MODEL_NAMESPACE, MODEL)) \
else return INVALID_MODEL;

#else


#define RUNMODEL(SOLVER, MODEL_NAMESPACE, MODEL) \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 1>::value) \
	if (dimension == 1) return RUNMODELNNN(SOLVER, MODEL, 3, 5, 4, 1, 2); \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 2>::value) \
	if (dimension == 2) return RUNMODELNNN(SOLVER, MODEL, 9, 13, 6, 2, 2); \
if constexpr (MODEL_NAMESPACE::allowed_model_dimensions<void, 3>::value) \
	if (dimension == 3) return RUNMODELNNN(SOLVER, MODEL, 15, 41, 10, 3, 2); \
return INVALID_MODEL;

#endif

#else

#endif


// **************************************************************************************



#define RESTRICT_DIMENSIONS(...) \
template<size_t D> \
struct allowed_model_dimensions<void, D> \
{ static const bool value = symphas::lib::value_in_seq<D, std::index_sequence<__VA_ARGS__>>::value; };

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
 * - #MODEL_DEF(...)
 * 
 * \param NAME The name of the model.
 * \param TYPES The order parameter types of the phase field, given in a bracketed list.
 * The types **must** be provided in a list surrounded by brackets, i.e. `(SCALAR, SCALAR)`.
 * \param ... The rest of the phase field definition.
 */
#define MODEL(NAME, TYPES, ...) \
namespace model_ ## NAME { \
template<typename, size_t D> \
struct allowed_model_dimensions { static const bool value = true; }; \
EQUATION_TRAIT_FORWARD_DECL \
template<size_t Dm, typename Sp> \
using OpTypes = typename ModelApplied<Dm, Sp>::template OpTypes<SINGLE_ARG TYPES>; \
template<typename> struct using_provisional { template<size_t Dm, typename Sp> using type = typename OpTypes<Dm, Sp>::template Specialized<TraitEquationModel>; }; \
__VA_ARGS__ } \
template<size_t Dm, typename Sp> \
using model_ ## NAME ## _t = model_ ## NAME ::SpecializedModel<Dm, Sp>;



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
 * Since they do not directly participate in the dynamics of the
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
 * OpFuncConvolution.
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
 * specified, #MODEL_DEF should be used. 
 * 
 * \param PREAMBLE Setup work in order to write the equations of motion. If
 * this part is not required, use the macro #MODEL_DEF.
 * \param ... The list of the equations of motion.
 */
#define MODEL_PREAMBLE_DEF(PREAMBLE, ...) \
template<size_t Dm, typename Sp> \
using SpecializedModel = typename using_provisional<void>::template type<Dm, Sp>; \
EQUATION_TRAIT_PREAMBLE(SINGLE_ARG PREAMBLE) \
EQUATION_TRAIT_DEFINITION(__VA_ARGS__)

//! Same as #MODEL_PREAMBLE_DEF.
#define MODEL_DEF_PREAMBLE(PREAMBLE, ...) MODEL_PREAMBLE_DEF((SINGLE_ARG PREAMBLE), __VA_ARGS__)

//! See #MODEL_PREAMBLE_DEF.
/*!
 * Shortcut macro for when no preamble is needed. Refer to #MODEL_PREAMBLE_DEF
 * for full information on defining the equations of motion of a model.
 */
#define MODEL_DEF(...) \
MODEL_PREAMBLE_DEF((), __VA_ARGS__)



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
MODEL_WRAPPER_FUNC(NAME, \
if (std::strcmp(name, #GIVEN_NAME) == 0) { RUNMODEL(AppliedSolver, model_ ## NAME, model_ ## NAME::SpecializedModel); })

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
 * - #DEFAULT_DYNAMIC(...) or #DEFAULT_MODE(...); or
 * - #NO_DEFAULTS (if no default is changed)
 *
 * Immediately following these definitions, the types of the order parameters
 * are provided in a bracketd list. Finally, the behaviour of individual phase
 * fields may be modified using:
 *
 * - #SPECIALIZED_DEF(...)
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
template<typename, size_t D> \
struct allowed_model_dimensions { static const bool value = true; }; \
struct PFCParameters : PFCParametersDefault<modelpfc_ ## NAME::PFCParameters> \
{ \
	using parent_type = PFCParametersDefault<modelpfc_ ## NAME::PFCParameters>; \
	template<size_t N> static constexpr DynamicType dynamic_val_apply() { return parent_type::template dynamic_val_apply<N>(); } \
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
using modelpfc_ ## NAME ## _t = modelpfc_ ## NAME ::ModelPFCSpecialized<Dm, Sp>;


//! Specify different default behaviours.
/*!
 * Provide a list of the different default behaviours of the phase field crystal
 * model, which will apply to all fields.
 */
#define DEFAULTS(...) __VA_ARGS__
#define NO_DEFAULTS

//! Sets a given dynamic type to be the default.
/*!
 * Changes the dynamics of all the fields, unless the field is
 * explicitly changed in the definition in #SPECIALIZED_DEF.
 * 
 * Used with the #DEFAULTS definition.
 * 
 * \param TYPE The dynamic type which dictates the default behaviour of
 * the phase field crystal problem.
 */
#define DEFAULT_DYNAMIC(TYPE) \
const static DynamicType DEFAULT_DYNAMIC = TYPE;

//! Sets the default mode of the phase field crystal problem.
/*!
 * The mode refers roughly to the how the densities are distributed. Typically
 * the mode is only of order 2, resulting in a tetrahedral pattern.
 * 
 * Changes the order of the mode for all fields, unless explicitly changed
 * in #SPECIALIZED_DEF. The mode should always be greater than or equal to 2.
 * 
 * Used with the #DEFAULTS definition.
 * 
 * \param N The numeric value of the phase field crystal mode.
 */
#define DEFAULT_MODE(N) \
const static size_t DEFAULT_MODE_N = N;

//! Indicates that a given field follows conserved dynamics.
#define PFC_CONSERVED DynamicType::CONSERVED

//! Indicates that a given field follows non-conserved dynamics.
#define PFC_NONCONSERVED DynamicType::NONCONSERVED



//! Modifies individual phase fields in the phase field crystal problem.
/*! 
 * The commands #DYNAMIC and #MODE are provided as a list of arguments
 * to this definition, one item for each of the phase field and characteristic
 * combinations which should be modified.
 * 
 * \param ... The list of modifications to the default behaviour of specific
 * fields.
 */
#define SPECIALIZED_DEF(...) __VA_ARGS__

//! Modify the dynamics of the given field.
/*!
 * Specializes a particular field with the given non-default behaviour 
 * (although a default characteristic can still be provided with undesired 
 * no side effects).
 * 
 * \param FIELD_N The field index which is specialized.
 * \param TYPE The dynamic type which dictates the new dynamics of the field.
 */
#define DYNAMIC(FIELD_N, TYPE) template<> constexpr DynamicType PFCParameters::dynamic_val_apply<FIELD_N - 1>() { return TYPE; }

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
#define MODE(FIELD_N, N) template<> constexpr size_t PFCParameters::mode_val_apply<FIELD_N - 1>() { return N; }


#ifdef USING_MODEL_SELECTION

//! Associate model definition with string.
/*!
 * See #LINK_WITH_NAME.
 */
#define LINK_PFC_WITH_NAME(NAME, GIVEN_NAME) \
NEXT_MODEL_INDEX(NAME) \
MODEL_WRAPPER_FUNC(NAME, \
if (std::strcmp(name, #GIVEN_NAME) == 0) { RUNMODEL(AppliedSolver, modelpfc_ ## NAME, modelpfc_ ## NAME::ModelPFCSpecialized); } )

#else

//! Without model selection, linking is an empty call.
#define LINK_PFC_WITH_NAME(NAME, GIVEN_NAME)

#endif



#define SYS_DIMS(N) parent_type::template system<N>().dims
#define SYS_INTERVAL(N, AXIS) parent_type::template system<N>().get_info().intervals


// **************************************************************************************


//! Real valued order parameter.
/*!
 * The order parameter will be of type ::scalar_t.
 */
#define SCALAR scalar_t
//! Complex value order parameter.
/*!
 * The order parameter will be of type ::complex_t.
 */
#define COMPLEX complex_t
//! Vector valued order parameter of real elements.
/*!
 * The order parameter will a vector type, specified by ::vector_t.
 */
#define VECTOR vector_t<Dm>


 // **************************************************************************************

//! The `i`-th order parameter.
/*!
 * The variable which represents the `i`-th order parameter. This is used
 * in an equation of motion. It is indexed from 1, meaning that `op(1)` is the
 * first order parameter.
 */
#define op(i) parent_type::template op<i - 1>()

//! The time derivative of the `i`-th order parameter.
/*!
 * The time derivative of the `i`-th parameter, which always appears on the
 * left hand side of the equals sign in an equation of motion.
 * It is indexed from 1, meaning that `dop(1)` is the derivative of the
 * first order parameter.
 */
#define dop(i) parent_type::template dop<i - 1>()

//! The `i`-th provisional variable.
/*!
 * The `i`-th provisional variable. It is indexed from 1, 
 * meaning that `var(1)` is the first provisional variable.
 */
#define var(i) parent_type::template var<i - 1>()

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
#define lap(EXPR) expr::laplacian(EXPR, solver)

//! The laplacian of an expression.
/*!
 * The bilaplacian is the derivative of 4th order. Applies the bilaplacian to
 * the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define bilap(EXPR) expr::bilaplacian(EXPR, solver)

//! The laplacian of an expression.
/*!
 * The gradient of the laplacian. Apply the gradlaplacian to the given 
 * expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define gradlap(EXPR) expr::gradlaplacian(EXPR, solver)

//! The gradient of an expression.
/*!
 * The gradient applies the first derivative to all spatial components of the
 * given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define grad(EXPR) expr::gradient(EXPR, solver)

//! The `N`-th coefficient of the model.
/*!
 * Refer to the value of the `N`-th parameter which was passed to the
 * phase field model. It is indexed from 1, 
 * meaning that `param(1)` is the first coefficient.
 * 
 * \param N The `N`-th parameter/coefficient.
 */
#define param(N) parent_type::template param<N - 1>()


// **************************************************************************************



//! The modulus of the given expression.
/*!
 * The modulus of the expression is computed, interpreting it as a complex
 * number. The expression can also be real.
 */
#define modulus(EXPR) expr::modulus(EXPR)

//! Take the conjugate of a complex number resulting from an expression.
/*!
 * Take the conjugate of the complex number that is computed by
 * evaluating the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define conj(EXPR) expr::conj(EXPR)

//! Take the imaginary part of a complex number resulting from an expression.
/*!
 * Take the imaginary part of the complex number that is computed by
 * evaluating the given expression.
 * 
 * \param EXPR The expression which is differentiated.
 */
#define Im(EXPR) expr::imag(EXPR)

//! Take the imaginary part of a complex number resulting from an expression.
/*!
 * Take the real part of the complex number that is computed by
 * evaluating the given expression.
 *
 * \param EXPR The expression which is differentiated.
 */
#define Re(EXPR) expr::real(EXPR)

//! Standard Gaussian with spatial distance equal to unity.
#define Gaussian GaussianSmoothing<model_dimension<parent_type>::value>{ \
	parent_type::template system<0>().dims, 1.0, 1.0 }

//! Gaussian with spatial width based on the i-th field.
#define Gaussian_F(i) GaussianSmoothing<model_dimension<parent_type>::value>{ \
	parent_type::template system<i>().dims, parent_type::template system<i>().get_info().get_widths().get(), 1.0, 1.0 }

//! Apply smoothing to the expression `Z`.
#define smoothing(Z) expr::make_convolution(Gaussian, Z)

//! Smoothing with Gaussian based on the spatial width of the i-th field.
#define smoothing_F(i, Z) expr::make_convolution(Gaussian_F(i), Z)

//! Use a number in an expression.
/*!
 * Create a number with the given value to use in an expression. 
 * This is a constant.
 * 
 * \param V The value of the number.
 */
#define lit(V) expr::make_literal(V)

//! The number 1 for an expression.
/*!
 * Uses the multiplicative identity in an expression, useful in particular
 * for ensuring that symbolic algebra rules are followed.
 */
#define one OpIdentity{}


//! The imaginary number for an expression.
/*!
 * An imaginary number literal that can be used in an expression.
 */
#define Ii lit(complex_t(0, 1))


#define x expr::make_var<Axis::X, D>(SYS_DIMS(0), SYS_INTERVAL(0, X))
#define y expr::make_var<Axis::Y, D>(SYS_DIMS(0), SYS_INTERVAL(0, Y))
#define z expr::make_var<Axis::Z, D>(SYS_DIMS(0), SYS_INTERVAL(0, Z))
#define t parent_type::get_time_var()


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

 //! @}



