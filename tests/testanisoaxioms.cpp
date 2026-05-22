// Compile-time axiom checks.  See testanisoaxioms.h for the rationale.
//
// Each `static_assert` below encodes one shape axiom from the algebra:
//
//   S1.  rank(scalar Term)             == 0
//   S2.  rank(vector Term in D dims)   == D
//   S3.  rank(a + b)                   == rank(a)              [iff equal]
//   S4.  rank(scalar * vector)         == D                    [scalar promotes]
//   S5.  rank(dot(v, w))               == 0
//   S6.  rank(grad(scalar))            == D
//   S7.  rank(div(vector))             == 0
//   S8.  rank(lap(e))                  == rank(e)
//   S9.  rank(pow<N>(scalar))          == 0
//
// C1 (rank query terminates without a hard-fail) is implicit: if the
// static_assert itself doesn't compile, C1 is violated for that shape.
//
// Several blocks are gated behind SYMPHAS_TEST_AXIOMS_STRICT.  Those
// are the axioms known-broken at HEAD (May 2026).  Default-off so a
// fresh checkout still builds; flip on to surface failures.

#include "testanisoaxioms.h"

#include <cstdio>
#include <type_traits>

namespace {

using Sp = SolverFT<Stencil2d2h<>>;

// ----- Build representative expression types ----------------------------
// We don't need values here -- only types -- but we still need the
// solver type wired up so derivative nodes are well-formed.

template <size_t D>
using scalar_grid_t = Grid<double, D>;
template <size_t D>
using vector_grid_t = Grid<VectorValue<double, D>, D>;

template <size_t D>
using scalar_term_t = decltype(expr::make_term(std::declval<scalar_grid_t<D>&>()));
template <size_t D>
using vector_term_t = decltype(expr::make_term(std::declval<vector_grid_t<D>&>()));

// Two-D versions for the bulk of the suite.
using S2 = scalar_term_t<2>;
using V2 = vector_term_t<2>;

constexpr size_t S2_rank = expr::eval_type<S2>::rank;
constexpr size_t V2_rank = expr::eval_type<V2>::rank;

// ----- S1: scalar leaf is rank 0 ---------------------------------------
static_assert(S2_rank == 0,
              "S1 violated: scalar Grid term must report rank 0");

// ----- S2: vector leaf is rank D ---------------------------------------
static_assert(V2_rank == 2,
              "S2 violated: vector Grid<VectorValue<T,D>,D> term must "
              "report rank D");

// ----- S3: rank-matched addition preserves rank -----------------------
using S2_plus_S2 = decltype(std::declval<S2>() + std::declval<S2>());
static_assert(expr::eval_type<S2_plus_S2>::rank == 0,
              "S3 violated: scalar + scalar must report rank 0");

using V2_plus_V2 = decltype(std::declval<V2>() + std::declval<V2>());
static_assert(expr::eval_type<V2_plus_V2>::rank == 2,
              "S3 violated: vector + vector must report rank D");

// ----- S4: scalar-times-anything keeps the other side's rank ----------
using lit_t = decltype(expr::make_literal(1.0));
using lit_times_S2 = decltype(std::declval<lit_t>() * std::declval<S2>());
static_assert(expr::eval_type<lit_times_S2>::rank == 0,
              "S4 violated: literal * scalar must be rank 0");

using lit_times_V2 = decltype(std::declval<lit_t>() * std::declval<V2>());
static_assert(expr::eval_type<lit_times_V2>::rank == 2,
              "S4 violated: literal * vector must be rank D (scalar "
              "should not change the vector's rank)");

// ----- S5: dot of two vectors is a scalar -----------------------------
using dot_VV = decltype(expr::dot(std::declval<V2>(), std::declval<V2>()));
static_assert(expr::eval_type<dot_VV>::rank == 0,
              "S5 violated: dot(vector, vector) must report rank 0");

// ----- S8: laplacian preserves rank (scalar branch) -------------------
// We need a solver instance to build a laplacian; use a default-constructed
// one purely for type deduction.
using lap_S2 = decltype(expr::make_laplacian(
    std::declval<S2>(), std::declval<solver_op_type<Sp>>()));
static_assert(expr::eval_type<lap_S2>::rank == 0,
              "S8 violated: lap(scalar) must report rank 0");

// ----- S9: pow<N>(scalar) is a scalar ---------------------------------
using pow2_S2 = decltype(expr::pow<2>(std::declval<S2>()));
static_assert(expr::eval_type<pow2_S2>::rank == 0,
              "S9 violated: pow<2>(scalar) must report rank 0");

using pow3_S2 = decltype(expr::pow<3>(std::declval<S2>()));
static_assert(expr::eval_type<pow3_S2>::rank == 0,
              "S9 violated: pow<3>(scalar) must report rank 0");

// ----- KNOWN-FAILING axioms -------------------------------------------
// The following are believed to be broken at HEAD (May 2026).  They are
// the targets of the Bug 4 / eval_type refactor.  Building with
// SYMPHAS_TEST_AXIOMS_STRICT surfaces them as compile failures.

#ifdef SYMPHAS_TEST_AXIOMS_STRICT

// ----- Closure of compositions used in Liam's model -------------------
// dot(M, grad psi) -- the exact shape AnisoT3c builds.  Building this
// type itself currently FAILS in non-STRICT compilation because the
// SFINAE-gated operator* between an OpOperatorDerivative and an
// OpTerms triggers expr::eval_type<OpOperatorDerivative<...>>::rank_<1>,
// which hits __is_complete_or_unbounded.  That is exactly Evidence-1
// of the eval_type architectural defect (operators.h:2967-style enable_if
// asking eval-as-shape).
//
// Until the declared-shape refactor lands these compositions can only
// be exercised under STRICT.  They will start passing one node family
// at a time as shape rules are declared.

using grad_op_t = decltype(expr::make_operator_derivative<1>(
    std::declval<solver_op_type<Sp>>()));
using grad_psi_t = decltype(std::declval<grad_op_t>() * std::declval<S2>());
// S6: grad(scalar) is rank D
static_assert(expr::eval_type<grad_psi_t>::rank == 2,
              "S6 violated: grad(scalar) must report rank D");

using dot_M_grad_psi_t =
    decltype(expr::dot(std::declval<V2>(), std::declval<grad_psi_t>()));
// S5 again, on a composite
static_assert(expr::eval_type<dot_M_grad_psi_t>::rank == 0,
              "S5 violated (composite): dot(M, grad psi) must report rank 0");

// S7: div(vector) must be a scalar.  Today divergence_of returns an
// expression whose rank is mis-reported under deep contexts because
// the un-applied operator form's eval is shape-deduced rather than
// shape-declared.  See review/_anisofmpfc_plan.md "I2".
using div_M_t = decltype(expr::divergence_of(
    std::declval<V2>(), std::declval<solver_op_type<Sp>>()));
static_assert(expr::eval_type<div_M_t>::rank == 0,
              "S7 violated: div(vector) must report rank 0");

// Bug 4: OpTensor<T,0,0,1,1> is algebraically a scalar but evaluates to
// a 1x1 matrix.  Lock the *expected* (math) answer here so when the
// refactor lands, the test starts passing.
using tensor_1x1_t = decltype(expr::make_tensor<0, 0, 1, 1>(
    expr::make_literal(1.0)));
static_assert(expr::eval_type<tensor_1x1_t>::rank == 0,
              "C2/C3 violated: OpTensor<T,0,0,1,1> denotes a scalar; "
              "its declared rank should be 0 (currently reports rank>0 "
              "because eval returns a 1x1 matrix).");

// Bug 4 composite: lap(power(dot(M, grad psi), 3) * div(M)) is rank 0.
// This is exactly M1 from the axiom list.  Building this static_assert
// is itself the C1 test for the model-target shape: if the rank query
// even completes, C1 holds for this expression.
using inner_pow_t = decltype(expr::pow<3>(std::declval<dot_M_grad_psi_t>()));
using inner_mul_t = decltype(std::declval<inner_pow_t>() * std::declval<div_M_t>());
using model_M1_t = decltype(expr::make_laplacian(
    std::declval<inner_mul_t>(), std::declval<solver_op_type<Sp>>()));
static_assert(expr::eval_type<model_M1_t>::rank == 0,
              "M1 violated: lap(pow<3>(dot(M,grad psi)) * div(M)) must "
              "report rank 0 (the user's model target shape).");

#endif  // SYMPHAS_TEST_AXIOMS_STRICT

}  // namespace

void testanisoaxioms() {
    // Nothing to run at runtime -- the file's value is its static_asserts.
    // Print a tiny banner so the suite output shows it was compiled.
    std::printf("--- testanisoaxioms: compile-time axiom checks OK ---\n");
}
