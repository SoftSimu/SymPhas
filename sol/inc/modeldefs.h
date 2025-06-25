
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
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

#include "expressiontypeincludes.h"
#include "modelpfc.h"
#include "modelvirtual.h"
#include "stencilincludes.h"

//! Copy construct the model with a new solver.
/*!
 * A new model is constructed with using instead the given solver. The
 * phase-field data from the given model is copied over as well.
 *
 * \tparam M The model type that from which is counted the phase fields
 */
template <typename Sp>
struct model_swap_solver {
 protected:
  template <template <size_t, typename> typename M, size_t D, typename Sp0>
  static constexpr auto with_new_solver(M<D, Sp0> model) {
    auto parameters = model.generate_parameters();
    return M<D, Sp>(model.get_coeff(), model.get_num_coeff(), parameters);
  }

  template <template <template <typename> typename,
                      typename> typename SpecializedModel,
            template <typename> typename Eq, size_t D, typename Sp0,
            typename... S>
  static constexpr auto with_new_solver(
      SpecializedModel<
          Eq, Eq<symphas::internal::MakeEquation<Model<D, Sp0, S...>>>> const&
          model) {
    using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<
        S...>::template Specialized<Eq>;
    auto parameters = model.generate_parameters();
    return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
  }

  template <template <template <typename> typename,
                      typename> typename SpecializedModel,
            template <size_t, typename> typename PFC, size_t D, typename Sp0,
            typename... S>
  static constexpr auto with_new_solver(
      SpecializedModel<symphas::internal::MakeEquation,
                       ModelPFCEquation<PFC, D, Sp0, S...>> const& model) {
    using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<
        S...>::template Specialized<symphas::internal::MakeEquation,
                                    ModelPFCEquation<PFC, D, Sp, S...>>;
    auto parameters = model.generate_parameters();
    return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
  }

  template <
      template <template <typename> typename, typename> typename Specialized,
      template <typename> typename Pr, size_t D, typename Sp0, typename... S,
      typename... P>
  static constexpr auto with_new_solver(
      Specialized<Pr, Pr<symphas::internal::MakeEquationProvisional<
                          Model<D, Sp0, S...>, P...>>> const& model) {
    using model_Sp = typename ModelApplied<D, Sp>::template OpTypes<
        S...>::template ProvTypes<P...>::template Specialized<Pr>;
    auto parameters = model.generate_parameters();
    return model_Sp(model.get_coeff(), model.get_num_coeff(), parameters);
  }

 public:
  template <typename M>
  auto operator()(M const& model) {
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

namespace symphas::internal {

/*!
 * Primary object using recursive inheritance in order to
 * define a member variable and return one of incremented value.
 *
 * A higher order index is given on the next usage by overloading
 * the function for which the terminating point is instantiated. below
 * it is never defined, because only the return type matters
 */
template <int N>
struct model_count_index : model_count_index<N - 1> {
  static const int value = N + 1;
};

//! Base specialization to terminate the recursive inheritance.
/*!
 * Base specialization which terminates the recursive inheritance for
 * incrementing model indices.
 */
template <>
struct model_count_index<0> {
  static const int value = 1;
};

/* overloading this function on a new parameter and updating the return
 * type will provide a new index on the next usage by decltype
 */
constexpr model_count_index<0> model_counter(model_count_index<0>);

}  // namespace symphas::internal

//! Used in assigning unique names to models for indexing.
/*!
 * The naming format of a model index is defined. Each index name has
 * to be different and requires a prefix.
 */
#define MODEL_INDEX_NAME(PREFIX_NAME, GIVEN_NAME) \
  __##PREFIX_NAME##_##GIVEN_NAME##_index

//! Iterates to the next model index for compile-time constant model indexing.
/*!
 * Convenience definition for setting the next index and providing
 * the PREFIX argument which names it.
 * Importantly, there cannot be more than #MAX_DEFINED_MODELS models defined
 * because after that, the counter will no longer increment.
 */
#define NEXT_MODEL_INDEX(PREFIX_NAME, GIVEN_NAME)                              \
  namespace symphas::internal {                                                \
  constexpr int MODEL_INDEX_NAME(PREFIX_NAME, GIVEN_NAME) =                    \
      decltype(model_counter(model_count_index<MAX_DEFINED_MODELS>{}))::value; \
  constexpr model_count_index<MODEL_INDEX_NAME(PREFIX_NAME, GIVEN_NAME)>       \
      model_counter(                                                           \
          model_count_index<MODEL_INDEX_NAME(PREFIX_NAME, GIVEN_NAME)>);       \
  }

//! \endcond

namespace symphas::internal {

template <template <size_t, typename> typename M, size_t D>
struct allowed_model_dimensions {
  static const bool value = true;
};

}  // namespace symphas::internal

/*
 * in order to link a defined model to a string, so that an arbitrary
 * function can be executed, a compile time index counter is maintained
 *
 * the index counter is used to specialize a template class and define
 * a member function which acts as a selection and wrapper for providing
 * the model type to the user defined function
 */

template <template <typename> typename model_apply_type, int N>
struct model_call_wrapper {
  template <template <typename, size_t> typename AppliedSolver, typename... Ts>
  static int call(size_t, StencilParams, const char*, Ts&&...) {
    return INVALID_MODEL;
  }

  template <template <size_t> typename AppliedSolver, typename... Ts>
  static int call(size_t, const char*, Ts&&...) {
    return INVALID_MODEL;
  }
};
