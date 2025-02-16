
#include "testexpression.h"

#include "timer.h"

using namespace symphas;

// **********************************************************************************

void testexpressioneval() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op0 = expr::make_term(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  printf(
      "simple test of checking evaluation of expression at a single point\n");

  grid0[0] = 2.;
  auto e2 = op0 * op0;                                  // 4
  auto e3 = expr::make_literal(2) * (op0 + op0 * op0);  // 12
  printf("%lf %lf\n", e2.eval(0), e3.eval(0));

  iter_type offset = solver.dims[0];
  grid0[offset + 0] = -1;
  grid0[offset + 1] = 2;
  grid0[offset + 2] = 3;

  printf(
      "more indepth test using the laplacian in 1d at the first (inner) grid "
      "point\n");

  // a test of whether the correct value is determined for some expressions
  auto l0 = expr::make_laplacian(op0, solver);

  expr::prune::update(l0);
  auto test1 = l0.eval(offset + 1);
  auto test2 = (-l0).eval(offset + 1);
  auto test3 = -l0.eval(offset + 1);
  printf("%lf %lf %lf\n", test1, test2, test3);

  auto test4 = (op0 * op0 * op0).eval(offset + 1);
  auto test5 = (op0 + op0 * op0).eval(offset + 1);
  auto test6 =
      ((op0 + expr::make_literal(2)) * (op0 + op0)).eval(grid0.dims[0] + 2);
  printf("%lf %lf %lf\n", test4, test5, test6);

  auto nl1 = op0 * op0 + expr::make_literal(2) * op0;            // = 8
  auto nl2 = (expr::make_operator_derivative<2>(solver) * op0);  // = -6
  auto nl3 = nl2 * op0;                                          // = -12
  printf("%lf %lf %lf\n", nl1.eval(offset + 1), nl2.eval(offset + 1),
         nl3.eval(offset + 1));
}

// **********************************************************************************

void testexpressionspeed() {
  len_type len = 1'000;

  srand((unsigned int)time(NULL));
  double constant = rand() / RAND_MAX;
  int n = 2;
  int test_len = static_cast<int>(len), test_run = 100'000'000;

  Grid<double, 1> arr({test_len, test_len});
  MyEquation<double> eq{constant, n};

  printf(
      "testing the difference in evaluation between different methods; test "
      "case "
      "involves a random number (simulating a long expression that takes "
      "longer to "
      "evaluate along with a derivative. Computed %d times.\n",
      test_run);

  MyWrapper weq(eq, &MyEquation<double>::eval, arr);
  {
    Time t0("direct evaluation");
    for (int i = 0; i < test_run; ++i) {
      arr.values[n + i % (test_len - n)] +=
          rand() / RAND_MAX +
          neighbour_average(arr.values + n + i % (test_len - n), n) + constant;
    }
  }
  {
    Time t0("evaluation using wrapping object");
    for (int i = 0; i < test_run; ++i) {
      arr.values[n + i % (test_len - n)] += weq.eval(n + i % (test_len - n));
    }
  }
  {
    Time t0("evaluation using function");
    for (int i = 0; i < test_run; ++i) {
      arr.values[n + i % (test_len - n)] +=
          eq.eval(arr, n + i % (test_len - n));
    }
  }
}

// **********************************************************************************

void testexpressionmodel() {
  Grid<double, 2> grid0({100, 100});
  auto psi = expr::make_term(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  // use the coefficients from the configuration
  double a = symphas::conf::config().model_settings.get_coeff(0);
  double b = symphas::conf::config().model_settings.get_coeff(1);

  // auto vop = Variable<2, ref<Grid<double, 2>>>(op);
  auto vop = expr::as_variable<2>(grid0);
  auto var_op = expr::make_term<1>(OpIdentity{}, grid0);
  auto mvar = var_op * var_op * var_op;

  // check the generated expression type from applying derivatives together with
  // other terms typical of an expression
  auto pfc0 = -(expr::make_literal(a) + expr::make_literal(a) * psi +
                expr::make_literal(b) * psi * psi) *
              psi;
  auto pfc1 = (expr::make_literal(b) * psi +
               expr::make_literal(2.) * expr::make_literal(b) *
                   expr::make_laplacian(psi, solver) +
               expr::make_bilaplacian(psi, solver));
}

void testexpressionmodelspeed() {
  Grid<double, 2> grid0({100, 100});
  Grid<double, 2> grid1({100, 100});
  grid::fill_random(grid0);
  grid::fill_random(grid1);
  auto op = expr::make_term(grid0);
  auto dop = OpLHS(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  // use the coefficients from the configuration
  double a = symphas::conf::config().model_settings.get_coeff(0);
  double b = symphas::conf::config().model_settings.get_coeff(1);

  int iter_cap = 100'000'000;
  double* arr = new double[iter_cap];
  double* arr2 = new double[iter_cap];
  double* arr3 = new double[iter_cap];

  Grid<double, 1> g(
      {iter_cap - BOUNDARY_DEPTH * 2, iter_cap - BOUNDARY_DEPTH * 2});
  Grid<double, 1> r(
      {iter_cap - BOUNDARY_DEPTH * 2, iter_cap - BOUNDARY_DEPTH * 2});
  auto arr_t = expr::make_term(1.0, g);

  srand((unsigned int)time(NULL));
  for (int i = 0; i < iter_cap; ++i) {
    arr[i] = rand() / 1000.;
    g.values[i] = arr[i];
  }

  printf(
      "compares how long it takes to directly compute a simple expression on "
      "an array %d "
      "cells long, with how long it takes to compute the same thing inside an "
      "expression object.\n",
      iter_cap);
  {
    Time t0("regular double array eval");
    for (int i = 0; i < iter_cap; ++i) {
      // modela
      // dpsi = (A - lit(4.) * B * psi * psi) * psi
      arr2[i] = (arr[i] + arr[i] * arr[i]) * arr[i];
    }
  }

  auto e = (OpLHS(r) = (arr_t + arr_t * arr_t) * arr_t);
  auto e2 = (arr_t + arr_t * arr_t) * arr_t;

  {
    Time t0("expression template eval");
    for (int i = 0; i < iter_cap; ++i) {
      arr3[i] = std::get<1>(e).eval(i);
    }
  }

  delete[] arr;
  delete[] arr2;
  delete[] arr3;

  printf(
      "tests how long it takes to evaluate model A using different methods.\n");

  iter_type steps = symphas::conf::config().simulation_settings.save.get_stop();
  Stencil2d2h<5, 13, 6> stencil(grid0.dims, 1.);

  {
    /* the speed of this 'raw form' is the same as the speed if e.eval is called
     */
    Time t0("raw form");
    for (int step = 0; step < steps; ++step) {
      iter_type i =
          BOUNDARY_DEPTH + (grid0.dims[0] + BOUNDARY_DEPTH) * BOUNDARY_DEPTH;
      for (iter_type m = 0, mm = grid0.dims[1]; m < mm; ++m) {
        for (iter_type l = 0, ll = grid0.dims[0]; l < ll; ++l) {
          grid1.values[i] = stencil.laplacian(grid0.values + i) +
                            (a - 4. * b * grid0.values[i] * grid0.values[i]) *
                                grid0.values[i];
          ++i;
        }
        i += BOUNDARY_DEPTH + BOUNDARY_DEPTH;
      }
    }
  }

  {
    /* the speed of this 'raw form' is the same as the speed if e.eval is called
     */

    auto l0 = expr::make_laplacian(op, solver);
    auto r0 =
        (OpLiteral<double>(a) - OpLiteral(4.) * OpLiteral(b) * op * op) * op;
    auto e0 = l0 + r0;

    Time t0("direct expression evaluation");
    for (int step = 0; step < steps; ++step) {
      iter_type i =
          BOUNDARY_DEPTH + (grid0.dims[0] + BOUNDARY_DEPTH) * BOUNDARY_DEPTH;
      for (iter_type m = 0, mm = grid0.dims[1]; m < mm; ++m) {
        for (iter_type l = 0, ll = grid0.dims[0]; l < ll; ++l) {
          grid1.values[i] = e0.eval(i);
          ++i;
        }
        i += BOUNDARY_DEPTH + BOUNDARY_DEPTH;
      }
    }
  }
  {
    /* this is extremely (EXTREMELY) slow when using the following in place of
     * the laplacian: OpFunction laplacian2(op, [&](auto a, auto n) { return
     * solver.laplacian(a, n); }); note the lack of reference on the 'a'
     * parameter
     */
    OpFunction laplacian2(
        op, [&](auto& a, auto n) { return solver.laplacian(a, n); });
    auto r0 = (dop = laplacian2 + (expr::make_literal(a) -
                                   expr::make_literal(4.) *
                                       expr::make_literal(b) * op * op) *
                                      op);
    Time t0("grid iteration evaluation");
    for (int step = 0; step < steps; ++step) {
      expr::result(r0);
    }
  }

#ifdef SOLVERFT
  SolverFT<Stencil2d2h<5, 13, 6>> solver2(op.dims);
  {
    Time t0("solver equation evaluation");
    for (int step = 0; step < 1000; ++step) {
      solver2.evaluate(expr::make_term(dop) =
                           expr::make_laplacian(expr::make_term(op), solver2) +
                           (expr::make_literal(a) -
                            expr::make_literal(4.) * expr::make_literal(b) *
                                expr::make_term(op) * expr::make_term(op)) *
                               expr::make_term(op));
    }
  }
#endif
}

// **********************************************************************************

void testexpressionoperator() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  // test the combination of different operator expressions with variables and
  // if they correctly apply the generalized operator objects
  auto d = expr::make_laplacian(op, solver) * expr::make_literal(3.);
  auto y = expr::make_literal(1) + expr::make_operator_derivative<2>(solver);
  auto y2 = y * y;
  auto y22 = expr::make_operator_derivative<2>(solver) * y;
  auto applyme = (expr::make_literal(1) + y);
  auto applymeop = (applyme * applyme) * op;
  auto applymeop2 = applyme * (applyme * op);

  // check if the dimension is correctly returned for this type of expression
  constexpr int dimension = expr::grid_dim<decltype(applymeop)>::dimension;

  // test correct printing of a complicated operator expression

  auto testdexp = op * ((y * y) * op);
  expr::printe(applymeop);
  expr::printe(testdexp);

  /* test if the distribution of the derivative operator is correctly applied
   */

  // distribution of derivative onto another derivative
  auto deriv_exp1 =
      expr::make_operator_derivative<2>(solver) *
      (expr::make_literal(2) + expr::make_operator_derivative<2>(solver));
  auto deriv_exp2 =
      (expr::make_literal(2) + expr::make_operator_derivative<2>(solver)) *
      expr::make_operator_derivative<2>(solver);

  // distribution of derivative operator onto a variable
  auto dv1 =
      expr::make_operator_derivative<2>(solver) * (expr::make_literal(2) * op);
  auto dv2 = expr::make_laplacian(op, solver);

  // a variable times a derivative
  auto m0 = op * expr::make_laplacian(op, solver);
  auto m1 = expr::make_literal(2) * m0;

  iter_type offset = grid0.dims[0] * BOUNDARY_DEPTH + BOUNDARY_DEPTH;

  // check the operator object
  auto opcc1 = OpOperatorCombination(op, expr::make_literal(2));
  auto opcc2 = opcc1 * op;
  printf("%lf\n", opcc2.eval(offset + 1));

  auto ovar_2 = expr::make_term<1>(grid0);
  auto ovar_4 = expr::make_term<2>(grid0);

  // test the various combinations of derivatives into operators
  auto opcc3 = opcc1 * opcc2;
  auto opcc4 = expr::make_operator_derivative<2>(solver) +
               expr::make_operator_derivative<4>(solver);
  auto opcc5 = expr::make_operator_derivative<2>(solver) + opcc4;
  auto opcc6 =
      expr::make_operator_derivative<4>(solver) + expr::make_literal(4.0);
  auto opcc7 = expr::make_operator_derivative<2>(solver) * opcc2;
  auto opcc8 = expr::apply_operators(opcc1 * (ovar_2 + ovar_4));

  // check whether the grid can be interpreted from the expression
  typename expr::storage_type<decltype(op)>::type gg(grid0.dims);
}

void testexpressionoperatorspeed() {
  Grid<double, 2> grid0({100, 100});
  Grid<double, 2> grid1({100, 100});
  grid::fill_random(grid0);
  grid::fill_random(grid1);
  auto op = expr::make_term(grid0);
  auto dop = expr::make_term(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  // use the coefficients from the configuration
  double a = symphas::conf::config().model_settings.get_coeff(0);
  double b = symphas::conf::config().model_settings.get_coeff(1);

  printf(
      "test the evaluation speed of operators evaluated on expressions, "
      "typical of "
      "the PFC case with operator combinations except simplified");

  auto opce1 = expr::make_operator_derivative<2>(solver) * op +
               (expr::make_literal(a) -
                expr::make_literal(4.) * expr::make_literal(b) * op * op) *
                   op;
  {
    Time t0("direct evaluation, no combination");
    expr::result(opce1, grid1);
  }

  auto opco1 =
      (expr::make_operator_derivative<2>(solver) + expr::make_literal(a)) -
      OpLiteral(4.) * OpLiteral(b) * op * op;
  auto opco2 = opco1 * op;
  {
    Time t0("operator combination");
    expr::result(expr::apply_operators(opco2), grid1);
  }
}

void testexpressionstate() {
  Grid<double, 2> grid0({100, 100});
  Grid<double, 2> grid1({100, 100});
  grid::fill_random(grid0);
  grid::fill_random(grid1);
  auto op = expr::make_term(grid0);
  auto dop = OpLHS(expr::make_term(grid1));
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  auto applyme =
      (expr::make_literal(1) + expr::make_operator_derivative<2>(solver));
  auto applymeop = (applyme * applyme) * op;

  // verifies that expressions which have a state are identified as such
  auto state_pair = (dop = expr::make_laplacian(op + op, solver));
  constexpr bool state1 =
      expr::has_state<decltype(expr::make_laplacian(op + op, solver) * op *
                               op)>::value;
  constexpr bool state2 = expr::has_state<decltype(state_pair)>::value;
  constexpr bool state3 = expr::has_state<decltype(std::make_tuple(
      dop = applymeop, dop = op + op, state_pair))>::value;
}

void testexpressioncollect() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto ovar_2 = expr::make_term<2>(grid0);
  auto ovar_4 = expr::make_term<4>(grid0);
  auto ovar_5 = expr::make_term<5>(grid0);

  // test of collecting like terms
  auto lt1 = op + op * op;
  auto lt2 = lt1 + lt1;
  auto lt3 = op + op;
  auto lt4 = expr::make_literal(3.4) * ovar_2 + ovar_2;
  auto adr = ovar_4 + ovar_2 + op + ovar_2 + ovar_4 + ovar_5;
}

void testexpressionconvolution() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);
  double h[2]{1, 1};

  // testing the correct behaviour of convolution
  auto cc1 =
      expr::make_convolution::get(GaussianSmoothing<2>{grid0.dims, h, 1}, op);
  auto cc2 =
      expr::make_literal(1.0) *
      expr::make_convolution::get(GaussianSmoothing<2>{grid0.dims, h, 1}, cc1);
  auto cc3 =
      expr::make_convolution::get(GaussianSmoothing<2>{grid0.dims, h, 1}, op) -
      op;
}

void testexpressionvariable() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);

  auto ovar_2 = expr::make_term<2>(grid0);

  // tests whether the nonlinear variables are correctly combined as like terms

  auto ovar_test_add2 =
      expr::make_literal(2.) * ovar_2 * ovar_2 + ovar_2 * ovar_2;
  using GG = Grid<double, 2>;  // decltype(expr::get<1>(ovar_2).data().get());
  constexpr auto is_nlc1 = expr::satisfies_combination<
      Variable<2, GG>, symphas::lib::types_list<Variable<2, GG>>,
      symphas::lib::types_list<Variable<2, GG>, Variable<2, GG>>>;
  constexpr auto is_nlc2 = expr::satisfies_combination<
      K<2>, symphas::lib::types_list<Variable<2, GG>>,
      symphas::lib::types_list<K<2>, Variable<2, GG>>>;
  constexpr auto is_nl1 = expr::is_combinable<WaveVectorData<2, 2>>;
}

void testexpressionswap() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto ovar_2 = expr::make_term<2>(grid0);
  auto ovar_4 = expr::make_term<4>(grid0);

  // the variable to swap in is a new one of variable id 3
  auto to_swp = expr::as_variable<3>(grid0);

  // check that the swap works correctly
  auto org_1 = ovar_4 * ovar_2 * ovar_2;
  auto swp_1 = expr::transform::swap_grid<2>(org_1, to_swp);

  expr::printe(org_1);
  expr::printe(swp_1);
}

void testexpressionsort() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  auto ovar_2 = expr::make_term<2>(grid0);

  // a selection of the greatest hits
  auto l0 = expr::make_laplacian(op, solver);
  auto y = expr::make_literal(1) + expr::make_operator_derivative<2>(solver);
  auto y2 = y * y;
  auto y22 = expr::make_operator_derivative<2>(solver) * y;
  auto applyme = (expr::make_literal(1) + y);

  // the result of determining the sort index for derivatives
  constexpr auto sortderiv0 = expr::derivative_index<0, OpIdentity>::value;
  constexpr auto sortderiv1 =
      expr::derivative_index<0, decltype(ovar_2)>::value;
  constexpr auto sortderiv2 = expr::derivative_index<0, decltype(l0)>::value;
  constexpr auto sortderiv3 =
      expr::derivative_index<0, decltype(applyme * op)>::value;
  constexpr auto sortderiv4 =
      expr::derivative_index<0, decltype(y2 * op)>::value;
  constexpr auto sortderiv5 =
      expr::derivative_index<0, decltype(y22 * op)>::value;
  constexpr auto sortderiv6 =
      expr::derivative_index<0, decltype(y * op)>::value;

  // this test establishes that the smallest over the given bound is returned
  // when the given bound is over 5 in this case, then the largest dervative
  // will be returned
  constexpr auto sortderiv7 =
      expr::derivative_index<0, decltype(y22 * op + y * op)>::value;
  constexpr auto sortderiv8 =
      expr::derivative_index<200, decltype(y22 * op + y * op)>::value;
  constexpr auto sortderiv9 =
      expr::derivative_index<2, decltype(y22 * op + op)>::value;
  constexpr auto sortderivA =
      expr::derivative_index<2, decltype(y22 * op + y * op + op)>::value;
  constexpr auto aaa = fixed_min<fixed_max<2, 1>, fixed_max<2, 4>>;
}

void testexpressiondivision() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto ovar_2 = expr::make_term<2>(grid0);
  auto ovar_4 = expr::make_term<4>(grid0);
  auto ovar_5 = expr::make_term<5>(grid0);

  // testing division
  auto dive1 = ovar_2 / ovar_4;
  auto dive2 = op / op;
  auto dive3 = op / (expr::make_literal(2) * op);
  auto divinv = expr::make_literal(2) / dive1;  // invert division
  auto divr1 = (ovar_2 * ovar_4) / ovar_2;      // simple cancellation
  auto divr2 =
      (ovar_4 / ovar_2) * ovar_2;  // simple cancellation by multiplication
  auto divr3 = ovar_2 * ovar_4 / ovar_5 / ovar_4;
  auto divr4 = ovar_4 * ovar_2 / ovar_5 / ovar_4;
  auto divr5 = (ovar_2 * ovar_4 / ovar_2 / ovar_4) * (ovar_2 * ovar_4);
  auto divr6 =
      (ovar_2 * ovar_4 * ovar_5 * ovar_2) / (ovar_5 * ovar_4 * ovar_2 * ovar_2);

  auto ee1 =
      (ovar_2 * ovar_2 * ovar_2) / (ovar_2 + ovar_2);  // this one should divide
  auto ee2 = (op * op * op) / (op + op);               // this should not divide
}

void testexpressionexponential() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto exp1 = expr::exp(op);
  auto exp2 = expr::exp(op * op);
}

// distinct from the evaluate function for expressions; tests the expression
// object which acts as its own simple evaluate block
void testexpressionevaluation() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto evo = expr::make_fourier_map(op);
  constexpr bool evob = expr::has_state<decltype(evo)>::value;
}

void testexpressionkgrid() {
  len_type dims[]{100, 100};
  double h[]{1, 1};
  K<2> kgrid1(dims, h);
  auto constexpr order1 = order_K_type<K<2>>::value;
}
#include "solversp.h"
#include "spectrallib.h"
void testderivativefactor() {
  Grid<double, 2> grid0({100, 100});
  Grid<VectorValue<double, 2>, 2> grid1({100, 100});

  len_type dims[2]{100, 100};
  double h[2]{1., 1.};

  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);
  auto opv = expr::make_term(grid1);
  SolverFT<Stencil2d2h<5, 6, 13>> solver(grid0.dims, 1.);

  auto [d_l, e_l] =
      expr::split::separate_operator(expr::make_laplacian(op, solver));

  auto e_simple1 = expr::make_laplacian(op, solver);
  auto e_simple2 = expr::make_laplacian(op * op, solver);
  auto e_simple3 = expr::make_laplacian(opv, solver);
  auto [d_simple1, ee_simple1] = expr::split::separate_operator(e_simple1);
  auto [d_simple2, ee_simple2] = expr::split::separate_operator(e_simple2);
  auto [d_simple3, ee_simple3] =
      expr::split::separate_operator(expr::get<0>(e_simple3));
  auto [d_simple4, ee_simple4] =
      expr::split::separate_operator(expr::get<1>(e_simple3));
  auto cee_simple3 =
      expr::coeff(ee_simple4) * expr::make_term(k_grid_type<2, 2>(dims, h)) +
      expr::coeff(ee_simple3) * expr::make_term(k_grid_type<2, 2>(dims, h));
  constexpr bool flagaa = std::is_same<
      decltype(cee_simple3),
      OpTerms<OpAdd<OpTensor<OpIdentity, 0, 2>, OpTensor<OpIdentity, 1, 2>>,
              Term<WaveVectorData<2, 2>>>>::value;
  constexpr size_t R = expr::eval_type<
      OpTerms<OpAdd<OpTensor<OpIdentity, 0, 2>, OpTensor<OpIdentity, 1, 2>>,
              Term<WaveVectorData<2, 2>>>>::rank;

  using namespace expr::symbols;

  constexpr auto sss = decltype(e_simple1)::axis;
  auto cee3 = expr::coeff(cee_simple3.term);
  constexpr bool flagaab = expr::is_tensor<decltype(cee3)>;
  constexpr bool flagaac =
      expr::is_coeff<decltype(expr::make_row_vector<0, R>())>;
  auto r0 = expr::exp(
      (expr::make_row_vector<1, R>() * expr::coeff(cee_simple3)) *
      OpTerms(OpIdentity{}, expr::terms_after_first(cee_simple3))) /*.eval(0)*/;
  /*auto r1 = expr::make_column_vector<1, R>() *
  expr::exp(expr::make_literal(2.) * ((expr::make_row_vector<1, R>() *
  expr::coeff(cee_simple3)) * OpTerms(OpIdentity{},
  expr::terms_after_first(cee_simple3)))); auto r2 = expr::make_column_vector<0,
  R>() * expr::exp(expr::make_literal(2.) * (expr::make_row_vector<0, R>() *
  cee_simple3));*/

  auto eq = (3_c + 4_c) * 5_c;

  /*auto te_impl = TE{};
  auto te_applied = expr::apply_operators(te_impl);*/

  // auto r = solver_sp::form_A_op<2>(cee_simple3, 1., dims,
  // std::index_sequence<0, 1>{});

  auto ed1 = expr::make_bilaplacian(op, solver);
  // auto ed2 = expr::make_operator_derivative<2>(solver) * ed1;

  // automatically factor by the greatest derivative order (the smallest common
  // derivative)
  auto [dd1, ee1] = expr::split::separate_operator(ed1);

  // try to factor by a derivative larger than the one available
  auto [dd1_e, ee1_e] = expr::split::separate_operator(ed1 + e_simple2);
  constexpr auto ee1i = expr::derivative_index<2, decltype(ee1)>::value;

  // checking the order
  constexpr auto sortderiv1 =
      expr::derivative_index<0, decltype(e_simple1)>::value;
  // constexpr auto sortderiv2 = expr::derivative_index<0,
  // decltype(ed2)>::value;

  // a far more complicated operator factoring
  auto comb1 = expr::make_operator_derivative<2>(solver) + OpIdentity{};
  auto chain1 = expr::make_operator_derivative<2>(solver) * comb1;
  auto chain2 = expr::make_operator_derivative<2>(solver)(comb1);
  auto comb2 = expr::make_operator_derivative<2>(solver) + chain1;
  auto [dd2, ee2] = expr::split::separate_operator(chain1 * (op * op));
}

void testexpressiondistribution() {
  Grid<double, 2> grid0({100, 100});
  grid::fill_random(grid0);
  auto op = expr::make_term(grid0);

  auto one = OpIdentity{};
  auto c = OpLiteral{5.1};
  auto e1 = -op * op + c * (one - op);
  auto e2 = op - op * op + c * (one - op);
}

void testexpressionnoisetype() {
  // using TE = OpAdd<
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 0, 2>, Term<NamedData<Grid<double,
  //         2>>>>, OpBinaryMul<OpTerms<OpIdentity, Term<NamedData<Grid<double,
  //         2>>>>,
  //                     OpSymbolicEval<OpTensor<OpIdentity, 0, 0, 1, 2>,
  //                                    NoiseData<expr::NoiseType::WHITE,
  //                                              VectorValue<double, 2>, 2>,
  //                                    SymbolicFunction<OpLiteral<int>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 1, 2>, Term<NamedData<Grid<double,
  //         2>>>>, OpBinaryMul<OpTerms<OpIdentity, Term<NamedData<Grid<double,
  //         2>>>>,
  //                     OpSymbolicEval<OpTensor<OpIdentity, 0, 1, 1, 2>,
  //                                    NoiseData<expr::NoiseType::WHITE,
  //                                              VectorValue<double, 2>, 2>,
  //                                    SymbolicFunction<OpLiteral<int>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 0, 0, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 0, 0, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 0, 0, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 0, 0, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 1, 1, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 1, 1, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 1, 1, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>,
  //     OpBinaryMul<
  //         OpTerms<OpTensor<OpIdentity, 1, 1, 2, 2>,
  //                 Term<NamedData<Grid<double, 2>>>>,
  //         OpMap<
  //             MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //             OpBinaryMul<
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<
  //                           OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //                 OpMap<MapGridInverseFourier<std::complex<double>, double,
  //                 2>,
  //                       OpIdentity,
  //                       OpTerms<OpLiteral<double>,
  //                               Term<NamedData<Grid<std::complex<double>,
  //                               2>>>, Term<Variable<
  //                                   2,
  //                                   NamedData<GridData<std::complex<double>,
  //                                                         2>>>>>>>>>>;

  // auto noise = OpBinaryMul<
  //     OpTerms<OpTensor<OpIdentity, 1, 2>, Term<NamedData<Grid<double, 2>>>>,
  //     OpBinaryMul<OpTerms<OpIdentity, Term<NamedData<Grid<double, 2>>>>,
  //                 OpSymbolicEval<OpTensor<OpIdentity, 0, 1, 1, 2>,
  //                                NoiseData<expr::NoiseType::WHITE,
  //                                          VectorValue<double, 2>, 2>,
  //                                SymbolicFunction<OpLiteral<int>>>>>{};

  // auto mul1 = OpBinaryMul<
  //     OpTerms<OpTensor<OpIdentity, 0, 0, 2, 2>,
  //             Term<NamedData<Grid<double, 2>>>>,
  //     OpMap<
  //         MapGridFourier<double, std::complex<double>, 2>, OpIdentity,
  //         OpBinaryMul<
  //             OpMap<MapGridInverseFourier<std::complex<double>, double, 2>,
  //                   OpIdentity,
  //                   OpTerms<OpLiteral<double>,
  //                           Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                           Term<Variable<2, NamedData<GridData<
  //                                                std::complex<double>,
  //                                                2>>>>>>,
  //             OpMap<MapGridInverseFourier<std::complex<double>, double, 2>,
  //                   OpIdentity,
  //                   OpTerms<
  //                       OpLiteral<double>,
  //                       Term<NamedData<Grid<std::complex<double>, 2>>>,
  //                       Term<Variable<2, NamedData<GridData<
  //                                            std::complex<double>,
  //                                            2>>>>>>>>>{};

  // auto eval_noise = noise.eval(0);
  // auto eval_mul = mul1.eval(0);
}

// ************************************************************************************************

size_t k_len = 20;
size_t t_len = 1000;

complex_t*** read_data() {
  const char dir_name[] = "C:/Users/Zirconix/Dropbox/FidelityCalc/";
  const char f_name_u[] = "udataL20vp0001.csv";
  const char f_name_v[] = "vdataL20vp0001.csv";
  const char* f_name[] = {f_name_u, f_name_v};

  double a, b;
  complex_t*** all_data = new complex_t**[2];

  for (int n = 0; n < 2; ++n) {
    char* open_name =
        new char[std::strlen(dir_name) + std::strlen(f_name[n]) + 1];
    sprintf(open_name, "%s%s", dir_name, f_name[n]);
    FILE* f = fopen(open_name, "r");
    delete[] open_name;

    complex_t** data = new complex_t*[t_len];
    all_data[n] = data;

    for (int t = 0; t < t_len; ++t)  // for each row
    {
      complex_t* data_row = new complex_t[k_len];
      data[t] = data_row;
    }

    char line[11'000'000];
    for (int k = 0; k < k_len; ++k)  // for each row
    {
      fgets(line, 11'000'000, f);
      char* tok = std::strtok(line, ",");

      // read the data at each column
      for (int t = 0; t < t_len; ++t) {
        sscanf(tok, "%lf%lf", &a, &b);
        data[t][k] = complex_t{a, b};
        tok = std::strtok(NULL, ",");
      }
    }

    fclose(f);
  }

  return all_data;
}

auto new_data() {
  int L = static_cast<int>(k_len);
  double hi = 0.5;

  double* thetak_put = new double[L];
  complex_t* u_data_put = new complex_t[L];
  complex_t* v_data_put = new complex_t[L];

  double ks, ak, bk;

  for (iter_type i = 0, v = -L + 1; v < L; ++i, v += 2) {
    ks = v * symphas::PI / L;
    ak = hi - cos(ks);
    bk = sin(ks);
    thetak_put[i] = 0.5 * atan(bk / ak);

    u_data_put[i] = cos(thetak_put[i]);
    v_data_put[i] = sin(thetak_put[i]);
  }

  return std::make_tuple(u_data_put, v_data_put, thetak_put);
}

template <size_t N>
auto get_data() {
  auto data = new_data();
  return std::get<N>(data);
}

complex_t* get_u_data() { return get_data<0>(); }

complex_t* get_v_data() { return get_data<1>(); }

double* get_theta_data() { return get_data<2>(); }

complex_t* u_data = get_u_data();
complex_t* v_data = get_v_data();
double* thetak = get_theta_data();

Number_uk uk;
Number_vmk vmk;
Number_umk umk;
Number_vk vk;
Number_udk udk;
Number_vdmk vdmk;
Number_udmk udmk;
Number_vdk vdk;
Eta_dk etadk;
Eta_mk etamk;
Eta_dmk etadmk;
Eta_k etak;

Number_vk_0 vk0;
Number_umk_0 umk0;
Number_udmk_0 udmk0;
Number_vdk_0 vdk0;

// constexpr char* ORDER_PARAMETER_NAMES[] = {
//	"u_{k}(t')", "v_{-k}(t)", "u_{-k}(t')", "v_{k}(t)",
//	"\\bar u_{k}(t)", "\\bar v_{-k}(t')", "\\bar u_{-k}(t)", "\\bar
// v_{k}(t')", 	"v_k", "u_{-k}", "\\bar u_{-k}", "\\bar v_{k}",
//	"\\eta^\\dagger_{k}", "\\eta_{-k}", "\\eta^\\dagger_{-k}", "\\eta_{k}"
//};

// constexpr char* ORDER_PARAMETER_NAMES[] = {
//	"utable[[p, tp + dt]]", "vtable[[mp, tp]]", "utable[[mp, tp + dt]]",
//"vtable[[p, tp]]", 	"Conjugate[utable[[p, tp]]]", "Conjugate[vtable[[mp, tp
//+ dt]]]", "Conjugate[utable[[mp, tp]]]", "Conjugate[vtable[[p, tp + dt]]]",
//	"vki[[p]]", "uki[[mp]]", "Conjugate[uki[[mp]]]", "Conjugate[vki[[p]]]"
//};

void test12varexpression() {
  complex_t*** data = read_data();

  auto op_uk = expr::make_term(uk);
  auto op_vmk = expr::make_term(vmk);
  auto op_umk = expr::make_term(umk);
  auto op_vk = expr::make_term(vk);
  auto op_udk = expr::make_term(udk);
  auto op_vdmk = expr::make_term(vdmk);
  auto op_udmk = expr::make_term(udmk);
  auto op_vdk = expr::make_term(vdk);
  auto op_etadk = expr::make_term(etadk);
  auto op_etamk = expr::make_term(etamk);
  auto op_etadmk = expr::make_term(etadmk);
  auto op_etak = expr::make_term(etak);

  auto op_thetak = expr::make_term(thetak);

  // adds the variables to the lists in the correct order

  expr::get_op_name(uk);
  expr::get_op_name(vmk);
  expr::get_op_name(umk);
  expr::get_op_name(vk);
  expr::get_op_name(udk);
  expr::get_op_name(vdmk);
  expr::get_op_name(udmk);
  expr::get_op_name(vdk);

  expr::get_op_name(vk0);
  expr::get_op_name(umk0);
  expr::get_op_name(udmk0);
  expr::get_op_name(vdk0);

  expr::get_op_name(etadk);
  expr::get_op_name(etamk);
  expr::get_op_name(etadmk);
  expr::get_op_name(etak);

  auto ddk = op_udk * op_etadk + op_vmk * op_etamk;
  auto ddmk = op_udmk * op_etadmk + op_vk * op_etak;
  auto dk = op_uk * op_etak + op_vdmk * op_etadmk;
  auto dmk = op_umk * op_etamk + op_vdk * op_etadk;

  auto second_moment_0 = eta_identity_prune(ddk * ddmk);
  auto second_moment_1 = eta_identity_prune(dmk * dk);
  auto fourth_moment = dmk * dk * ddk * ddmk;

  expr::printe(fourth_moment);

  using expr::cos;
  using expr::sin;

  // auto fourth_moment = (sin(op_thetak) * dmk * dk) * (sin(op_thetak) * ddk *
  // ddmk); auto e_right = cos(op_thetak) + sin(op_thetak) * ddk * ddmk; auto
  // e_left = cos(op_thetak) + sin(op_thetak) * dmk * dk;
  //  auto expectation = e_left * e_right;

  auto expectation =
      cos(op_thetak) * cos(op_thetak) +
      cos(op_thetak) * sin(op_thetak) * (second_moment_0 + second_moment_1) +
      sin(op_thetak) * sin(op_thetak) * fourth_moment;

  for (iter_type t = 0; t < t_len; ++t) {
    u_data = data[0][t];
    v_data = data[1][t];

    uk.update_data();
    vmk.update_data();
    umk.update_data();
    vk.update_data();
    udk.update_data();
    vdmk.update_data();
    udmk.update_data();
    vdk.update_data();

    complex_t product = 1;
    for (iter_type i = 0; i < k_len; ++i) {
      complex_t value = expectation.eval(i);
      product *= value;
    }
    printf("<E|E> = %.4lf + %.4lfi ~ %.24lf\n", real(product), imag(product),
           std::abs(product) * std::abs(product));
  }
}

void testccliketypes() {
  using vv = Variable<0, NamedData<std::complex<double>*>>;
  using v1 = Variable<1, NamedData<std::complex<double>*>>;
  using G = OpTerm<OpLiteral<double>, vv>;
  using Gg = OpTerm<OpLiteral<double>, v1>;
  using G2 = mul_result_t<G, G>;
  using G21 = mul_result_t<G2, Gg>;

  using cct = symphas::lib::combine_collect_like_types<G, G>;
  using cc_type = typename cct::type;
  using cc_count = typename cct::count_seq;
  using opgt = typename expr::op_types<G21>::type;
}

void testcoefftype() {
  using vv = Variable<0, NamedData<std::complex<double>*>>;
  using v1 = Variable<1, NamedData<std::complex<double>*>>;
  using Gneg = OpTerm<OpNegIdentity, vv>;
  using Gpos = OpTerm<OpIdentity, vv>;

  constexpr bool is_neg2 = expr::has_nmi_coeff<Gpos>;
  constexpr bool is_neg1 = expr::has_nmi_coeff<Gneg>;
  constexpr bool is_pos1 = expr::has_pmi_coeff<Gpos>;
  constexpr bool is_pos2 = expr::has_val_coeff<Gneg>;
}

void testorganizederivative() {
  len_type dims[]{100, 100, 100};
  SolverFT<Stencil3d2h<7, 10, 21>> solver(dims, 1.);
  {
    Grid<double, 3> grid0(dims);
    grid::fill_random(grid0);
    auto op0 = expr::make_term(grid0);
    SolverFT<Stencil3d2h<7, 10, 21>> solver(grid0.dims, 1.);

    auto opx2 = expr::make_operator_directional_derivative<Axis::X, 2>(solver);
    auto opy2 = expr::make_operator_directional_derivative<Axis::Y, 2>(solver);
    auto opz2 = expr::make_operator_directional_derivative<Axis::Z, 2>(solver);

    auto eq = opx2 * op0 + opy2 * op0 + opz2 * op0;
  }

  {
    Grid<any_vector_t<double, 3>, 3> grid0(dims);
    grid::fill_random(grid0);
    auto op0 = expr::make_term(grid0);

    auto opx2 = expr::make_operator_directional_derivative<Axis::X, 2>(solver);
    auto opy2 = expr::make_operator_directional_derivative<Axis::Y, 2>(solver);
    auto opz2 = expr::make_operator_directional_derivative<Axis::Z, 2>(solver);

    auto xx = expr::get<0>(opx2 * op0);
    auto yx = expr::get<0>(opy2 * op0);
    auto zx = expr::get<0>(opz2 * op0);

    auto xxx = xx + yx + zx;

    auto eq = opx2 * op0 + opy2 * op0 + opz2 * op0;
  }
  {
    auto lap = expr::make_operator_derivative<2>(solver);

    Grid<double, 3> grid0(dims);
    Grid<any_vector_t<double, 3>, 3> grid1(dims);
    grid::fill_random(grid0);
    auto op0 = expr::make_term<1>(grid0);
    auto op1 = expr::make_term<1>(grid1);
    auto opd = expr::make_operator_derivative<1>(solver)(op0);

    // auto list = std::make_tuple(
    //     Term(expr::as_variable<1>(NamedData(SymbolicData(grid0)))));
    // using namespace expr::symbols;
    // using ii = i_<0>;
    // using v = v_<ii>;

    // auto lfe =
    //     expr::sum(expr::landau_fe(solver, v{})).template
    //     select<0>(ii{})(list);

    // expr::sum_over(lfe) / ii{};?
    auto lfe = expr::landau_fe(solver, op0);

    auto intg0 = expr::make_integral(lfe, op0) + op1 * op1;
    auto fe_derv0 = expr::make_functional_derivative(intg0, op0);
    auto intg1 = expr::make_integral(lfe, op0);
    auto fe_derv1 = expr::make_functional_derivative(intg1, op0);
    auto intg2 = expr::make_integral(lfe + op1 * op1, op0);
    auto fe_derv2 = expr::make_functional_derivative(intg2, op0);

    auto ee0 = expr::apply_operators(fe_derv0);
    auto ee1 = expr::apply_operators(fe_derv1);
    auto ee2 = expr::apply_operators(fe_derv2);
  }
}
