
#include "testoperator.h"


void testoperator()
{
	using add_t = OpAdd<OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::directional_derivative<Axis::Y, 2>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>>, OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>>, OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::directional_derivative<Axis::X, 2>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>>, OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>>, OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>>, OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>>, OpTerms<OpLiteral<double>, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 2>>, OpTerms<OpLiteral<double>, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>, 4>>>;
	auto addit = add_t{};
	auto term = Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>{};
	auto sym = substitute_arg<0>(NamedData<SymbolicData<BoundaryGrid<double, 2>>>{});
	auto swapped = expr::transform::swap_grid<Variable<0>>(add_t{}, sym);
	constexpr bool f = std::is_same<add_t, decltype(swapped)>::value;

	static symphas::internal::GeneratedStencilApply stencil{ expr::get_central_space_stencil<6, 2, 2>() };

	auto solver = SolverFT<Stencil2d2h< 5, 6, 13>>::make_solver(symphas::problem_parameters_type(0));
	auto op = expr::make_operator_derivative<2>(solver);
	auto opg = expr::make_operator_derivative<1>(solver);
	Grid<vector_t<2>, 2> grid({ 10, 10 });
	Grid<double, 2> grid2({ 10, 10 });
	auto term = expr::make_term<1>(grid);
	auto serm = expr::make_term<2>(grid2);
	auto term2 = term * term;
	auto termp = expr::pow<2>(term);
	using ev2 = expr::eval_type_t<decltype(term2)>;
	using ev2p = expr::eval_type_t<decltype(termp)>;
	auto test = (term * op) * term;
	auto app = expr::apply_operators(test);
	
	auto flo = op * (expr::pow<2>(1. + op) * serm);
	flo.print(stdout);
	printf("\n");
	auto floa = expr::apply_operators(flo);
	floa.print(stdout);
	printf("\n");

	auto flo1 = expr::pow<2>(term * opg) * term;
	flo1.print(stdout);
	printf("\n");
	auto flo1a = expr::apply_operators(flo1);
	flo1a.print(stdout);
	printf("\n");

	auto flo2 = expr::transpose(term) * opg(term);
	flo2.print(stdout);
	printf("\n");
	auto flo2a = expr::apply_operators(flo2);
	flo2a.print(stdout);
	printf("\n");


	auto div_flo3 = 1. / (1. + serm);
	auto chain_flo3 = term * opg;
	auto flo3 = (div_flo3 * chain_flo3) * term;
	flo3.print(stdout);
	printf("\n");
	auto flo3a = expr::apply_operators(flo3);
	flo3a.print(stdout);
	printf("\n");

	auto term_tensor = OpTensor<OpLiteral<double>, 0, 2>{} *(term * term);
	auto deriv_tensor = expr::make_derivative< SolverFT<Stencil2d2h< 5, 6, 13>>::derivative<Axis::X, 1>>(term_tensor, solver);

	auto e = expr::make_literal(2.) * expr::make_unit_vector<2>(2 * expr::symbols::Pi * (NoiseData<expr::NoiseType::POISSON, scalar_t, 2>(grid::dim_list(), 1.)(OpIdentity{})));

	auto e0 = expr::get<0>(e);

	using types = OpAdd<
		OpIntegral<OpIdentity, 
			OpAdd<
				OpBinaryMul<
					OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, 
					OpBinaryMul<
						OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpTensor<OpIdentity, 1, 2>, OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, 
						OpSymbolicEval<OpIdentity, SymbolicSeries<symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>, Substitution<SymbolicDataArray<SolverSystemFD<double, 2>>>, symphas::lib::types_list<OpTerms<OpIdentity, Term<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>, 2>>, symphas::lib::types_list<expr::symbols::i_<1, 0>>, symphas::lib::types_list<expr::series_limits<OpAdd<OpBinaryMul<OpLiteral<double>, DynamicIndex>, OpLiteral<double>>, int>>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>>>>, SymbolicFunction<OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 2>>, Variable<0, NamedData<SymbolicData<SolverSystemFD<double, 2>>>>, Variable<1, OpLiteral<int>>>>>>, 
				OpBinaryMul<
					OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, 
					OpBinaryMul<
						OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::Y, 1>, OpTensor<OpIdentity, 1, 2>, OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, 
						OpSymbolicEval<OpIdentity, SymbolicSeries<symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>, Substitution<SymbolicDataArray<SolverSystemFD<double, 2>>>, symphas::lib::types_list<OpTerms<OpIdentity, Term<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>, 2>>, symphas::lib::types_list<expr::symbols::i_<1, 0>>, symphas::lib::types_list<expr::series_limits<int, OpBinaryMul<OpLiteral<double>, DynamicIndex>>>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>>>>, SymbolicFunction<OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 2>>, Variable<0, NamedData<SymbolicData<SolverSystemFD<double, 2>>>>, Variable<1, OpLiteral<int>>>>>>, 
				OpBinaryMul<
					OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, 
					OpBinaryMul<
						OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpTensor<OpIdentity, 0, 2>, OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, 
						SolverFT<Stencil2d2h<9, 6, 13>>>, OpSymbolicEval<OpIdentity, SymbolicSeries<symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>, Substitution<SymbolicDataArray<SolverSystemFD<double, 2>>>, symphas::lib::types_list<OpTerms<OpIdentity, Term<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>, 2>>, symphas::lib::types_list<expr::symbols::i_<1, 0>>, symphas::lib::types_list<expr::series_limits<int, OpBinaryMul<OpLiteral<double>, DynamicIndex>>>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>>>>, SymbolicFunction<OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 2>>, Variable<0, NamedData<SymbolicData<SolverSystemFD<double, 2>>>>, Variable<1, OpLiteral<int>>>>>>, 
				OpBinaryMul<
					OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>,	
					OpBinaryMul<OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>>>::derivative<Axis::X, 1>, OpTensor<OpIdentity, 0, 2>, OpTerms<OpIdentity, Term<DynamicVariable<NamedData<SolverSystemFD<double, 2> *>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>>>, 
				OpSymbolicEval<OpIdentity, SymbolicSeries<symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>, Substitution<SymbolicDataArray<SolverSystemFD<double, 2>>>, symphas::lib::types_list<OpTerms<OpIdentity, Term<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>, 2>>, symphas::lib::types_list<expr::symbols::i_<1, 0>>, symphas::lib::types_list<expr::series_limits<OpAdd<OpBinaryMul<OpLiteral<double>, DynamicIndex>, OpLiteral<double>>, int>>, symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<1, 0>>>>>, SymbolicFunction<OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 2>>, Variable<0, NamedData<SymbolicData<SolverSystemFD<double, 2>>>>, Variable<1, OpLiteral<int>>>>>>>, expr::variational_t<void>>, 
		OpFunctionApply<&symphas::math::cos, OpIdentity, OpSymbolicEval<OpLiteral<double>, SymbolicListIndex<OpAdd<DynamicIndex, OpIdentity>, expr::symbols::i_<0, 0>>, SymbolicFunction<OpSymbolicEval<OpIdentity, NoiseData<expr::NoiseType::POISSON, double, 2>, SymbolicFunction<OpTerms<OpIdentity, Term<TimeValue, 1>>, Variable<2, expr::symbols::i_<0, 0>>>>>>>, 
		OpFunctionApply<&symphas::math::sin, OpIdentity, OpSymbolicEval<OpLiteral<double>, SymbolicListIndex<OpAdd<DynamicIndex, OpIdentity>, expr::symbols::i_<0, 0>>, SymbolicFunction<OpSymbolicEval<OpIdentity, NoiseData<expr::NoiseType::POISSON, double, 2>, SymbolicFunction<OpTerms<OpIdentity, Term<TimeValue, 1>>, Variable<2, expr::symbols::i_<0, 0>>>>>>>>;


	using types1 = OpAdd<
		OpDerivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<1, 1>, OpTensor<OpFractionLiteral<2, 1>, 0, 1, 2, 2>, OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::Y, std::reference_wrapper<Grid<VectorValue<double, 2>, 2>>>>, 1>>, SolverFT<Stencil2d2h<5, 6, 13>>>, 
		OpDerivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<0, 2>, OpTensor<OpFractionLiteral<2, 1>, 1, 1, 2, 2>, OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::Y, std::reference_wrapper<Grid<VectorValue<double, 2>, 2>>>>, 1>>, SolverFT<Stencil2d2h<5, 6, 13>>>, 
		OpDerivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<0, 2>, OpTensor<OpFractionLiteral<2, 1>, 1, 0, 2, 2>, OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::X, std::reference_wrapper<Grid<VectorValue<double, 2>, 2>>>>, 1>>, SolverFT<Stencil2d2h<5, 6, 13>>>, 
		OpDerivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<1, 1>, OpTensor<OpFractionLiteral<2, 1>, 0, 0, 2, 2>, OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::X, std::reference_wrapper<Grid<VectorValue<double, 2>, 2>>>>, 1>>, SolverFT<Stencil2d2h<5, 6, 13>>>>;

	using evt0 = expr::eval_type_t< types1>;

	symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_mixed_stencil<2>(std::index_sequence<0, 2>{}));


	auto l1 = evt1{};
	auto grid = Grid < VectorValue<double, 2>, 2>({ 10, 10 });
	auto term = expr::make_term<0>(grid);
	auto termx = expr::make_row_vector<1, 2>() * term;
	auto l2 = OpTensor<OpFractionLiteral<2, 1>, 1, 1, 2, 2>{} *expr::make_derivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<0, 2>>(termx, SolverFT<Stencil2d2h<5, 6, 13>>::make_solver(symphas::problem_parameters_type{ 0 }));
	auto l2_ = l2.eval(0);

	(types1{} + types2{}).print(stdout);
	printf("\n");
	types3{}.print(stdout);
	printf("\n");



	constexpr size_t N = symphas::lib::index_of_value<size_t, 2, 0, 1, 2, 3>;

}

