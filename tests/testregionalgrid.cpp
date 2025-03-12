
#include "testregionalgrid.h"

#include "symphas.h"

void testregionalgrid() {
  len_type side_length_small = 16;
  auto test_dimensions = grid::dim_list{side_length_small, side_length_small};
  RegionalGrid<double, 2> grid(test_dimensions);
  for (iter_type x = 0; x < 10; ++x) {
    for (iter_type y = 0; y < 10; ++y) {
      if (x > 1 && x < 5 && y > 3 && y < 6) {
        grid(x, y) = 1.;
      }
    }
  }

  grid(1, 2) = 1.;
  grid(1, 3) = 1.;
  grid(2, 3) = 1.;
  symphas::io::print_grid(grid);

  iter_type new_origin[2]{4, 4};
  grid.adjust(new_origin);
  symphas::io::print_grid(grid);

  len_type intervals[2][2]{{8, 8 + side_length_small},
                           {8, 8 + side_length_small}};

  grid::resize_adjust_region(grid, intervals);
  symphas::io::print_grid(grid);

  new_origin[0] = 10;
  new_origin[1] = 10;
  grid.adjust(new_origin);
  symphas::io::print_grid(grid);

  grid::resize_adjust_region(grid, 1.);
  symphas::io::print_grid(grid);
  grid::resize_adjust_region(grid, 1.);
  symphas::io::print_grid(grid);
  grid::adjust_region(grid, 1.);
  symphas::io::print_grid(grid);

  // grid.adjust((iter_type[2]) { 0, 0 });
  // symphas::io::print_grid(grid);

  grid.adjust((iter_type[2]){2, 2});
  symphas::io::print_grid(grid);
  grid.adjust((iter_type[2]){0, 0}, grid.dims);
  symphas::io::print_grid(grid);
  symphas::io::print_region(grid);

  auto interval0 = grid::get_iterable_domain(grid);
  auto interval1 = grid::get_iterable_domain(grid);

  grid::region_interval<2> remove_region0(grid.dims);
  for (iter_type i = 0; i < 2; ++i) {
    interval0[i][1] -= 1;
    remove_region0[i][0] = grid.region.boundary_size + 2;
    remove_region0[i][1] = grid.dims[i] - grid.region.boundary_size - 2;
  }

  grid::region_interval_multiple<2> regions0(grid.dims,
                                             grid.region.boundary_size);

  regions0 += interval0;
  symphas::io::print_region(regions0);
  printf("\n");
  symphas::io::print_region(remove_region0);
  printf("\n");

  regions0 /= remove_region0;
  symphas::io::print_region(regions0);
  printf("\n");

  symphas::io::print_region(remove_region0 - regions0);
  printf("\n");

  grid::region_interval<2> interval_diff(grid.dims);
  interval_diff[1][0] = 4;
  interval_diff[1][1] = 7;
  symphas::io::print_region(remove_region0 - regions0 - interval_diff);
  printf("\n");

  auto test_diff_region =
      remove_region0 - (remove_region0 - regions0 - interval_diff);
  symphas::io::print_region(test_diff_region);
  printf("\n");
  grid::region_interval_multiple<2> regions2(grid.dims,
                                             grid.region.boundary_size);
  regions2 += interval0;
  symphas::io::print_region(regions2 / test_diff_region);
  printf("\n");
  symphas::io::print_region(regions2 + test_diff_region);
  printf("\n");

  grid::region_interval<2> middle_interval(grid.dims);
  for (iter_type i = 0; i < 2; ++i) {
    middle_interval[i][0] = 5;
    middle_interval[i][1] = 12;
  }
  symphas::io::print_region(middle_interval);
  printf("\n");
  symphas::io::print_region(middle_interval / (regions2 + test_diff_region));
  printf("\n");

  grid::region_interval<2> left_interval(grid.dims);
  left_interval[0][0] = 0;
  left_interval[0][1] = 7;
  symphas::io::print_region(left_interval);
  printf("\n");
  symphas::io::print_region(left_interval / (regions2 + test_diff_region));
  printf("\n");

  for (auto region : regions0) {
    printf("%d ", region[0][0]);
    region[0][0] += 1;
  }
  printf("\n");

  grid.adjust((iter_type[2]){8, 8}, (iter_type[2]){0, 0});
  symphas::io::print_grid(grid);
  symphas::io::print_region(grid);

  grid::region_interval<2> remove_region1(grid.dims);
  for (iter_type i = 0; i < 2; ++i) {
    remove_region1[i][0] = grid.region.boundary_size + 2;
    remove_region1[i][1] = grid.dims[i] - grid.region.boundary_size - 2;
  }

  grid::region_interval_multiple<2> regions1(grid.dims,
                                             grid.region.boundary_size);
  regions1 += interval1;
  regions1 /= remove_region1;

  for (grid::region_interval<2> region : regions1) {
    printf("%d ", region[0][0]);
  }
  printf("\n");

  iter_type region_dims[2];
  for (iter_type i = 0; i < 2; ++i) {
    region_dims[i] = grid.dims[i] - grid.region.boundary_size * 2;
  }
  Grid<double, 2> region(region_dims);
  grid::region_index_list<2> list0(grid.dims, grid.region.boundary_size);
  for (iter_type n = 0; n < region.len; ++n) {
    region[n] = grid[list0[n]];
  }
  symphas::io::print_grid(region);
  symphas::io::print_region(grid);

  grid(0, 0) = 0;
  grid(1, 2) = 0;
  symphas::io::print_grid(grid);
  grid::resize_adjust_region(grid, 1.);
  symphas::io::print_region(grid);
  symphas::io::print_grid(grid);

  grid::region_interval<2> test_region_combine;
  for (iter_type i = 0; i < 2; ++i) {
    test_region_combine.dims[i] = 100;
    test_region_combine[i][0] = 80;
    test_region_combine[i][1] = 115;
  }

  grid::region_interval_multiple<2> combined_multiple(test_region_combine.dims,
                                                      3);
  combined_multiple += test_region_combine;
  auto test_region_combine2 = +combined_multiple;
  grid::region_interval_multiple<2> combined_multiple2(test_region_combine.dims,
                                                       3);
  combined_multiple += test_region_combine2;

  auto term = expr::make_term<0>(grid);
  auto s = &expr::get<1>(term).data();

  len_type repeat = 40;
  len_type side_length = 2000;
  auto grid_dimensions = grid::dim_list{side_length, side_length};

  // using E = OpTerms<OpLiteral<double>,
  // Term<expr::symbols::v_<expr::symbols::i_<0, 0>>, 1>>; using G = typename
  // expr::base_data_type<E>::type; constexpr bool flag_combinable =
  // expr::is_combinable<E>;

  // symphas::io::print_grid(grid);
  // symphas::io::print_region(grid);

  grid = RegionalGrid<double, 2>((iter_type[2]){256, 256});
  grid.adjust((iter_type[2]){250, 239}, (iter_type[2]){136, 142});
  auto term01 = expr::make_term(grid);
  auto interval01 = expr::iterable_domain(2. * term01 * term01);

  grid = RegionalGrid<double, 2>(grid_dimensions);
  RegionalGrid<double, 2> frid(grid_dimensions);
  RegionalGrid<double, 2> result_grid(grid_dimensions);

  auto serm = expr::make_term<1>(frid);
  auto expr = term * term * term * serm + serm * serm * term;

  double *values0 = new double[side_length * side_length]{};
  double *values1 = new double[side_length * side_length]{};
  double *values2 = new double[side_length * side_length]{};
  double *values3 = new double[side_length * side_length]{};
  auto vterm0 = expr::make_term<0>(values0);
  auto vterm1 = expr::make_term<1>(values1);
  auto exprv = vterm0 * vterm0 * vterm0 * vterm1 + vterm1 * vterm1 * vterm0;

  auto interval = expr::iterable_domain(expr);
  symphas::data_iterator_region ittest(grid, interval);
  auto last = ittest + 1900;
  grid::fill_random(values0, side_length * side_length);
  grid::fill_random(values1, side_length * side_length);

  SolverFT<Stencil2d2h<5, 6, 13>> solver{grid.dims, 1.};
  auto dx2 = expr::make_operator_directional_derivative<Axis::X, 2>(solver);
  auto dy2 = expr::make_operator_directional_derivative<Axis::Y, 2>(solver);
  auto term00 = 0.35 * term;
  auto lapp = (dx2 * term) + (dy2 * term);

  using ii = expr::symbols::i_<0, 0>;
  using v_ii = expr::symbols::v_<ii>;
  using namespace expr::symbols;

  double *arr = new double[10]{};
  arr[1] = 2;
  arr[0] = 3;

  auto cond = ii{} != 2_n;
  auto sum = expr::sum(v_ii{}).select(ii{} != 2_n)(expr::array_arg(10, arr));

  // test reduce iterator:
  iter_type ind = 0;
  DynamicIndex dyn(ind);
  auto op = expr::make_term_dynamic(dyn, arr);
  auto expr_reduce_it = sum + op + op * op;
  double res;

  expr::prune::update(sum);
  expr::result_by_term<expr::matches_series>(sum, res);
  expr::prune::update(expr_reduce_it);
  expr::result_by_term<expr::matches_series>(expr_reduce_it, res);

  using namespace symphas::internal;
  auto test_reduce_it =
      *(symphas::reduce_iterator<
          expr::expression_iterator_group, 2,
          OpTerms<
              OpCoeff<OpLiteral<double>, DynamicIndex>,
              Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 2>>,
          OpTerms<
              OpCoeff<double, DynamicIndex>,
              Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>,
          OpBinaryMul<OpIntegral<OpIdentity,
                                 OpTerms<OpLiteral<double>,
                                         Term<DynamicVariable<NamedData<
                                                  SolverSystemFD<double, 2> *>>,
                                              2>>,
                                 expr::variational_t<symphas::grid_info>>,
                      OpTerms<OpNegFractionLiteral<2, 1>,
                              Term<DynamicVariable<
                                       NamedData<SolverSystemFD<double, 2> *>>,
                                   1>>>,
          OpTerms<
              OpCoeff<double, DynamicIndex>,
              Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 3>>,
          OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                              0>::derivative<Axis::X, 2>,
                       OpCoeff<OpLiteral<double>, DynamicIndex>,
                       OpTerms<OpIdentity,
                               Term<DynamicVariable<
                                        NamedData<SolverSystemFD<double, 2> *>>,
                                    1>>,
                       SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
          OpBinaryMul<
              OpIntegral<
                  OpIdentity,
                  OpAdd<
                      OpBinaryMul<
                          OpTerms<OpCoeff<OpLiteral<double>, DynamicIndex>,
                                  Term<DynamicVariable<NamedData<
                                           SolverSystemFD<double, 2> *>>,
                                       1>>,
                          OpBinaryMul<
                              OpDerivative<
                                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                                         0>::derivative<Axis::X, 1>,
                                  OpIdentity,
                                  OpTerms<
                                      OpIdentity,
                                      Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
                              OpSymbolicEval<
                                  OpIdentity,
                                  SymbolicSeries<
                                      ReduceOp<
                                          symphas::internal::SeriesOp::ADD>,
                                      Substitution<SymbolicDataArray<NamedData<
                                          SolverSystemFD<double, 2> *>>>,
                                      symphas::lib::types_list<
                                          OpTerms<
                                              OpIdentity,
                                              Term<expr::symbols::v_id_type<
                                                       expr::symbols::i_<1, 0>>,
                                                   2>>,
                                          symphas::lib::types_list<
                                              expr::symbols::i_<1, 0>>,
                                          symphas::lib::types_list<
                                              expr::series_limits<
                                                  OpAdd<OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpLiteral<double>>,
                                                  int>>,
                                          symphas::lib::types_list<
                                              expr::symbols::v_id_type<
                                                  expr::symbols::i_<1, 0>>>>>,
                                  SymbolicFunction<
                                      OpTerms<
                                          OpIdentity,
                                          Term<
                                              Variable<
                                                  0, std::reference_wrapper<
                                                         NamedData<SymbolicData<
                                                             SolverSystemFD<
                                                                 double, 2>>>>>,
                                              2>>,
                                      Variable<0,
                                               NamedData<SymbolicData<
                                                   SolverSystemFD<double, 2>>>>,
                                      Variable<1, OpLiteral<int>>>>>>,
                      OpBinaryMul<
                          OpTerms<OpCoeff<OpLiteral<double>, DynamicIndex>,
                                  Term<DynamicVariable<NamedData<
                                           SolverSystemFD<double, 2> *>>,
                                       1>>,
                          OpBinaryMul<
                              OpDerivative<
                                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                                         0>::derivative<Axis::X, 1>,
                                  OpIdentity,
                                  OpTerms<
                                      OpIdentity,
                                      Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
                              OpSymbolicEval<
                                  OpIdentity,
                                  SymbolicSeries<
                                      symphas::internal::ReduceOp<
                                          SeriesOp::ADD>,
                                      Substitution<SymbolicDataArray<NamedData<
                                          SolverSystemFD<double, 2> *>>>,
                                      symphas::lib::types_list<
                                          OpTerms<
                                              OpIdentity,
                                              Term<expr::symbols::v_id_type<
                                                       expr::symbols::i_<1, 0>>,
                                                   2>>,
                                          symphas::lib::types_list<
                                              expr::symbols::i_<1, 0>>,
                                          symphas::lib::types_list<
                                              expr::series_limits<
                                                  int,
                                                  OpAdd<OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpLiteral<double>>>>,
                                          symphas::lib::types_list<
                                              expr::symbols::v_id_type<
                                                  expr::symbols::i_<1, 0>>>>>,
                                  SymbolicFunction<
                                      OpTerms<
                                          OpIdentity,
                                          Term<
                                              Variable<
                                                  0, std::reference_wrapper<
                                                         NamedData<SymbolicData<
                                                             SolverSystemFD<
                                                                 double, 2>>>>>,
                                              2>>,
                                      Variable<0,
                                               NamedData<SymbolicData<
                                                   SolverSystemFD<double, 2>>>>,
                                      Variable<1, OpLiteral<int>>>>>>>,
                  expr::variational_t<symphas::grid_info>>,
              OpDerivative<
                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                         0>::derivative<Axis::X, 1>,
                  OpIdentity,
                  OpTerms<OpIdentity, Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>>,
          OpBinaryMul<
              OpIntegral<
                  OpIdentity,
                  OpAdd<
                      OpBinaryMul<
                          OpTerms<OpCoeff<OpLiteral<double>, DynamicIndex>,
                                  Term<DynamicVariable<NamedData<
                                           SolverSystemFD<double, 2> *>>,
                                       1>>,
                          OpBinaryMul<
                              OpDerivative<
                                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                                         0>::derivative<Axis::Y, 1>,
                                  OpIdentity,
                                  OpTerms<
                                      OpIdentity,
                                      Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
                              OpSymbolicEval<
                                  OpIdentity,
                                  SymbolicSeries<
                                      symphas::internal::ReduceOp<
                                          SeriesOp::ADD>,
                                      Substitution<SymbolicDataArray<NamedData<
                                          SolverSystemFD<double, 2> *>>>,
                                      symphas::lib::types_list<
                                          OpTerms<
                                              OpIdentity,
                                              Term<expr::symbols::v_id_type<
                                                       expr::symbols::i_<1, 0>>,
                                                   2>>,
                                          symphas::lib::types_list<
                                              expr::symbols::i_<1, 0>>,
                                          symphas::lib::types_list<
                                              expr::series_limits<
                                                  OpAdd<OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpLiteral<double>>,
                                                  int>>,
                                          symphas::lib::types_list<
                                              expr::symbols::v_id_type<
                                                  expr::symbols::i_<1, 0>>>>>,
                                  SymbolicFunction<
                                      OpTerms<
                                          OpIdentity,
                                          Term<
                                              Variable<
                                                  0, std::reference_wrapper<
                                                         NamedData<SymbolicData<
                                                             SolverSystemFD<
                                                                 double, 2>>>>>,
                                              2>>,
                                      Variable<0,
                                               NamedData<SymbolicData<
                                                   SolverSystemFD<double, 2>>>>,
                                      Variable<1, OpLiteral<int>>>>>>,
                      OpBinaryMul<
                          OpTerms<OpCoeff<OpLiteral<double>, DynamicIndex>,
                                  Term<DynamicVariable<NamedData<
                                           SolverSystemFD<double, 2> *>>,
                                       1>>,
                          OpBinaryMul<
                              OpDerivative<
                                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                                         0>::derivative<Axis::Y, 1>,
                                  OpIdentity,
                                  OpTerms<
                                      OpIdentity,
                                      Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
                              OpSymbolicEval<
                                  OpIdentity,
                                  SymbolicSeries<
                                      symphas::internal::ReduceOp<
                                          SeriesOp::ADD>,
                                      Substitution<SymbolicDataArray<NamedData<
                                          SolverSystemFD<double, 2> *>>>,
                                      symphas::lib::types_list<
                                          OpTerms<
                                              OpIdentity,
                                              Term<expr::symbols::v_id_type<
                                                       expr::symbols::i_<1, 0>>,
                                                   2>>,
                                          symphas::lib::types_list<
                                              expr::symbols::i_<1, 0>>,
                                          symphas::lib::types_list<
                                              expr::series_limits<
                                                  int,
                                                  OpAdd<OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpBinaryMul<
                                                            OpLiteral<double>,
                                                            DynamicIndex>,
                                                        OpLiteral<double>>>>,
                                          symphas::lib::types_list<
                                              expr::symbols::v_id_type<
                                                  expr::symbols::i_<1, 0>>>>>,
                                  SymbolicFunction<
                                      OpTerms<
                                          OpIdentity,
                                          Term<
                                              Variable<
                                                  0, std::reference_wrapper<
                                                         NamedData<SymbolicData<
                                                             SolverSystemFD<
                                                                 double, 2>>>>>,
                                              2>>,
                                      Variable<0,
                                               NamedData<SymbolicData<
                                                   SolverSystemFD<double, 2>>>>,
                                      Variable<1, OpLiteral<int>>>>>>>,
                  expr::variational_t<symphas::grid_info>>,
              OpDerivative<
                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                         0>::derivative<Axis::Y, 1>,
                  OpIdentity,
                  OpTerms<OpIdentity, Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>>,
          OpBinaryMul<
              OpFunctionApply<
                  &symphas::math::cos<double>, double,
                  OpSymbolicEval<
                      OpLiteral<double>,
                      SymbolicListIndex<
                          OpAdd<DynamicIndex, DynamicIndex, OpIdentity>,
                          expr::symbols::i_<0, 0>>,
                      SymbolicFunction<OpSymbolicEval<
                          OpIdentity,
                          NoiseData<expr::NoiseType::POISSON, double, 2, Grid>,
                          SymbolicFunction<
                              OpLiteral<double>,
                              Variable<2,
                                       OpTerms<OpIdentity, Term<TimeValue, 1>>>,
                              Variable<3, expr::symbols::i_<0, 0>>>>>>>,
              OpDerivative<
                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                         0>::derivative<Axis::X, 1>,
                  OpIdentity,
                  OpTerms<OpIdentity, Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>>,
          OpBinaryMul<
              OpFunctionApply<
                  &symphas::math::sin<double>, double,
                  OpSymbolicEval<
                      OpLiteral<double>,
                      SymbolicListIndex<
                          OpAdd<DynamicIndex, DynamicIndex, OpIdentity>,
                          expr::symbols::i_<0, 0>>,
                      SymbolicFunction<OpSymbolicEval<
                          OpIdentity,
                          NoiseData<expr::NoiseType::POISSON, double, 2, Grid>,
                          SymbolicFunction<
                              OpLiteral<double>,
                              Variable<2,
                                       OpTerms<OpIdentity, Term<TimeValue, 1>>>,
                              Variable<3, expr::symbols::i_<0, 0>>>>>>>,
              OpDerivative<
                  Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>,
                         0>::derivative<Axis::Y, 1>,
                  OpIdentity,
                  OpTerms<OpIdentity, Term<DynamicVariable<NamedData<
                                               SolverSystemFD<double, 2> *>>,
                                           1>>,
                  SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>{});

  auto start =
      symphas::reduce_iterator(op.begin(symphas::it_grp, interval),
                               (op * op).begin(symphas::it_grp, interval));
  auto it = expr::forward_value{}(*(++std::_Get_unwrapped(start)));
  auto end = symphas::reduce_iterator(op.end(symphas::it_grp, interval));

  // grid::copy(values0, grid);
  // grid::copy(values1, frid);
  grid::fill_random(grid);
  grid::fill_random(frid);
  // expr::result(expr, grid);
  printf("starting\n");

  iter_type index = (frid.region.boundary_size + side_length / 2) +
                    side_length * (frid.region.boundary_size + side_length / 2);

  auto start3 = std::chrono::high_resolution_clock::now();

  iter_type n;
  ;
  for (iter_type r = 0; r < repeat; ++r) {
    n = frid.region.boundary_size + side_length * frid.region.boundary_size;
    for (iter_type y = frid.region.boundary_size;
         y < side_length - frid.region.boundary_size; ++y) {
      for (iter_type x = frid.region.boundary_size;
           x < side_length - frid.region.boundary_size; ++x) {
        values2[n] =
            exprv[n];  // values0[n] * values0[n] * values0[n] * values1[n] +
                       // values1[n] * values1[n] * values0[n];
        ++n;
      }
      n += frid.region.boundary_size * 2;
    }
  }

  auto stop3 = std::chrono::high_resolution_clock::now();
  auto duration3 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3)
          .count();
  printf("%.4E seconds, %lf\n", duration3 / 1e6, values2[index]);

  auto start0 = std::chrono::high_resolution_clock::now();

  // std::copy(
  //	exprv.begin(interval),
  //	exprv.end(interval), values3_it);

  for (iter_type r = 0; r < repeat; ++r) {
    symphas::data_iterator_region values3_it(values3, interval);
    for (auto expr_it = exprv.begin(interval); expr_it < exprv.end(interval);
         ++expr_it, ++values3_it) {
      *values3_it = *expr_it;
    }
  }

  auto stop0 = std::chrono::high_resolution_clock::now();
  auto duration0 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0)
          .count();
  printf("%.4E seconds, %lf\n", duration0 / 1e6, values3[index]);
  printf("increase: %.2lf times\n", (double)duration0 / duration3);

  bool flag = true;
  n = frid.region.boundary_size + side_length * frid.region.boundary_size;
  for (iter_type y = frid.region.boundary_size;
       y < side_length - frid.region.boundary_size; ++y) {
    for (iter_type x = frid.region.boundary_size;
         x < side_length - frid.region.boundary_size; ++x) {
      flag = flag && (values3[n] == values2[n]);
      ++n;
    }
    n += frid.region.boundary_size * 2;
  }
  printf("all same? %s\n", (flag) ? "True" : "False");

  auto start7 = std::chrono::high_resolution_clock::now();

  grid::region_index_list<2> list(interval);
  for (iter_type r = 0; r < repeat; ++r) {
    symphas::data_iterator_selection gvalues3_it(values3, list);
    for (auto expr_it = expr::expression_iterator_selection(exprv, list);
         expr_it < expr::expression_iterator_selection(
                       exprv, list,
                       (side_length - grid.region.boundary_size * 2) *
                           (side_length - grid.region.boundary_size * 2));
         ++gvalues3_it, ++expr_it) {
      *gvalues3_it = *expr_it;
    }
  }

  auto stop7 = std::chrono::high_resolution_clock::now();
  auto duration7 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop7 - start7)
          .count();
  printf("%.4E seconds, %lf\n", duration7 / 1e6, values3[index]);
  printf("increase selection: %.2lf times\n", (double)duration7 / duration3);

  flag = true;
  n = frid.region.boundary_size + side_length * frid.region.boundary_size;
  for (iter_type y = frid.region.boundary_size;
       y < side_length - frid.region.boundary_size; ++y) {
    for (iter_type x = frid.region.boundary_size;
         x < side_length - frid.region.boundary_size; ++x) {
      flag = flag && (values3[n] == values2[n]);
      ++n;
    }
    n += frid.region.boundary_size * 2;
  }
  printf("all same? %s\n", (flag) ? "True" : "False");

  auto start6 = std::chrono::high_resolution_clock::now();

  // std::copy(
  //	exprv.begin(interval),
  //	exprv.end(interval), values3_it);
  for (iter_type r = 0; r < repeat; ++r) {
    symphas::data_iterator_group gvalues3_it(result_grid, interval);
    for (auto expr_it = expr::expression_iterator_group(expr, interval);
         expr_it <
         expr::expression_iterator_group(
             expr, interval, side_length - grid.region.boundary_size * 2);
         ++gvalues3_it, ++expr_it) {
      *gvalues3_it = *expr_it;
    }
  }

  auto stop6 = std::chrono::high_resolution_clock::now();
  auto duration6 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop6 - start6)
          .count();
  printf("%.4E seconds, %lf\n", duration6 / 1e6, values3[index]);
  printf("increase group: %.2lf times\n", (double)duration6 / duration3);

  // flag = true;
  // n = frid.region.boundary_size + side_length * frid.region.boundary_size;
  // for (iter_type y = frid.region.boundary_size; y < side_length -
  // frid.region.boundary_size; ++y)
  //{
  //	for (iter_type x = frid.region.boundary_size; x < side_length -
  // frid.region.boundary_size; ++x)
  //	{
  //		flag = flag && (values3[n] == values2[n]);
  //		++n;
  //	}
  //	n += frid.region.boundary_size * 2;
  // }
  // printf("all same? %s\n", (flag) ? "True" : "False");

  auto start1 = std::chrono::high_resolution_clock::now();
  for (iter_type r = 0; r < repeat; ++r) {
    expr::result(expr, result_grid);
  }
  auto stop1 = std::chrono::high_resolution_clock::now();
  auto duration1 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1)
          .count();

  printf("%.4E seconds, %lf\n", duration1 / 1e6, (double)result_grid[index]);
  printf("increase expr: %.2lf times\n", (double)duration1 / duration3);

  flag = true;
  n = frid.region.boundary_size + side_length * frid.region.boundary_size;
  for (iter_type y = frid.region.boundary_size;
       y < side_length - frid.region.boundary_size; ++y) {
    for (iter_type x = frid.region.boundary_size;
         x < side_length - frid.region.boundary_size; ++x) {
      if (result_grid(x - frid.region.boundary_size,
                      y - frid.region.boundary_size) != result_grid[n]) {
        flag = false;
        break;
      }
      ++n;
    }
    n += frid.region.boundary_size * 2;
  }
  printf("all same? %s\n", (flag) ? "True" : "False");

  auto start2 = std::chrono::high_resolution_clock::now();
  for (iter_type r = 0; r < repeat; ++r) {
    n = frid.region.boundary_size + side_length * frid.region.boundary_size;
    for (iter_type y = frid.region.boundary_size;
         y < side_length - frid.region.boundary_size; ++y) {
      for (iter_type x = frid.region.boundary_size;
           x < side_length - frid.region.boundary_size; ++x) {
        result_grid[n] = expr[n];
        ++n;
      }
      n += frid.region.boundary_size * 2;
    }
  }
  auto stop2 = std::chrono::high_resolution_clock::now();
  auto duration2 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2)
          .count();
  printf("%.4E seconds\n", duration2 / 1e6);
  printf("increase direct: %.2lf times\n", (double)duration2 / duration3);

  auto start4 = std::chrono::high_resolution_clock::now();
  for (iter_type r = 0; r < repeat; ++r) {
    symphas::data_iterator_region data_it(result_grid, interval);
    for (auto expr_it = expr.begin(interval); expr_it < expr.end(interval);) {
      *data_it++ = *expr_it++;
    }
  }
  auto stop4 = std::chrono::high_resolution_clock::now();
  auto duration4 =
      std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4)
          .count();
  printf("%.4E seconds\n", duration4 / 1e6);
  printf("increase region: %.2lf times\n", (double)duration4 / duration3);
}
