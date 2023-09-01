
#include "testoperator.h"
#include "solver.h"
#include "stencilincludes.h"
#include "conf.h"
#include "modelspecialized.h"

NEW_SOLVER_WITH_STENCIL(SolverFT)

/*
 * the parameter is a tuple with the first parameter being a Variable object with ref
 * base type, and the ref is a reference to the system; like this, the variable index is
 * packaged and the equality operator is deferred to the oplvariable
 */
	template<typename S, typename E>
	inline void equation(std::pair<S, E>& r) const
	{
		auto& [sys, equation] = r;
		expr::prune::update(equation);
		expr::result(equation, sys.get().dframe);
	}

	static auto make_solver()
	{
		len_type dims[2]{ 20,20 };
		double h = 1.0;
		return this_type{ &dims[0], h};
	}
};

void testoperator()
{

	//auto term = Variable<0, std::reference_wrapper<NamedData<SymbolicData<BoundaryGrid<double, 2>>>>>{};
	//auto sym = expr::substitute_arg<0>(NamedData<SymbolicData<BoundaryGrid<double, 2>>>{});

	//static symphas::internal::GeneratedStencilApply stencil{ expr::get_central_space_stencil<6, 2, 2>() };

	len_type dims[2]{ 20,20 };
	double h = 1.0;
	SolverFT<Stencil2d2h< 5, 6, 13>> solver(&dims[0], h);

	Grid<vector_t<2>, 2> grid({ 10, 10 });
	Grid<double, 2> grid2({ 10, 10 });
	 
	auto op = expr::make_operator_derivative<2>(solver);
	auto opg = expr::make_operator_derivative<1>(solver);
	auto term = expr::make_term<1>(grid);
	auto serm = expr::make_term<2>(grid2);
	auto term2 = term * term;
	auto termp = expr::pow<2>(term);
	using ev2 = expr::eval_type_t<decltype(term2)>;
	using ev2p = expr::eval_type_t<decltype(termp)>;
	
	/*

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
	*/

	using namespace expr::symbols;

	// an interesting functionality to generate an interval.
	auto constant_range1 = 2_n --> 8.4_n | 1.6_n;
	auto constant_range2 = 2_n --> 8.4_n & 4_n;
	auto nn = 2_c;

	// "shortcut" symphas way to make derivatives
	auto d = d_(solver) / dx ^ 2_n;
	auto functional = dF_(x{});

	// getting the non-directional component
	auto op2_2 = d_(solver) ^ 2_n;

	// testing automatic symbolic derivative of functions
	auto f = (fn(x{}, y{}) = (x{} ^ 2_n) + y{});
	auto f_prime = D(f, 2_n);
	auto f_prime_x = f_prime(x{});

	auto dxx = d_(solver) / dx;
	auto dyy = d_(solver) / dy;
	auto basic_deriv0 = d_(solver)(serm);
	auto basic_deriv0_2 = basic_deriv0 * basic_deriv0;
	
	auto term_tensor = OpTensor<OpLiteral<double>, 0, 2>{} *(term * term);
	auto deriv_tensor = expr::make_derivative<SolverFT<Stencil2d2h< 5, 6, 13>>::derivative<Axis::X, 1>>(term_tensor, solver);
	auto cell_fe_grad = expr::apply_operators(expr::apply_operators(expr::make_operator_derivative<1>(solver)(serm)) * serm);
	auto cell_fe_grad_a = cell_fe_grad.a;
	auto cell_fe_grad_b = cell_fe_grad.b;

	double* dt = new double(0.05);
	double* h_ptr = &h;

	auto termin = expr::make_term(dt);
	auto termin2 = termin * termin;
	auto add = termin2 + termin;
	auto testswapexpr0 = expr::transform::swap_expression<decltype(termin2)>(add, OpIdentity{});
	constexpr bool ff = expr::satisfies<decltype(termin), expr::matches_with<expr::matches_with<decltype(termin2)>>>;
	
	
	auto noisetest = expr::make_literal(2.) * expr::make_unit_vector<2>(2 * expr::symbols::Pi * (NoiseData<expr::NoiseType::POISSON, scalar_t, 2>(grid::dim_list(), h_ptr, dt)(OpIdentity{})));
	auto noisetest0 = expr::get<0>(noisetest);

	//symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_mixed_stencil<2>(std::index_sequence<0, 2>{}));

	auto termx = expr::make_row_vector<1, 2>() * term;
	auto l2 = OpTensor<OpFractionLiteral<2, 1>, 1, 1, 2, 2>{} *expr::make_derivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<0, 2>>(termx, solver);
	auto l2_ = l2.eval(0);

	//using EE = OpBinaryMul<
	//		OpOperatorChain<
	//			OpIdentity, 
	//			OpAdd<
	//				OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::Y, 1>, OpTensor<OpIdentity, 1, 2>, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//				OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 1>, OpTensor<OpIdentity, 0, 2>, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>, 
	//		OpAdd<
	//			OpOperatorCombination<
	//				OpOperatorChain<
	//					OpOperatorDirectionalDerivative<Axis::Y, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//					OpOperatorChain<
	//						OpIdentity, 
	//						OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::Y, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>, 
	//				OpOperatorChain<
	//					OpOperatorDirectionalDerivative<Axis::X, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//					OpOperatorChain<
	//						OpIdentity, 
	//						OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>>>>;


	//using EE0 = OpAdd<
	//	OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 4>, OpNegFractionLiteral<1, 2>, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//	OpOperatorChain<
	//		OpOperatorDerivative<2, OpNegIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//		OpOperatorCombination<
	//			OpOperatorChain<
	//				OpOperatorDirectionalDerivative<Axis::Y, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//				OpOperatorChain<
	//					OpIdentity, 
	//					OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::Y, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>, 
	//			OpOperatorChain<
	//				OpOperatorDirectionalDerivative<Axis::X, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//				OpOperatorChain<OpIdentity, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>>>, 
	//	OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 2>, OpNegIdentity, OpTerms<OpNegIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 3>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//	OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 2>, OpNegIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//	OpBinaryMul<
	//		OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::Y, 1>, OpNegIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//		OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::Y, NamedData<std::reference_wrapper<BoundaryGrid<VectorValue<double, 2>, 2>>>>>, 1>>>, 
	//	OpBinaryMul<
	//		OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 1>, OpNegIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, 
	//		OpTerms<OpIdentity, Term<Variable<1, VectorComponent<Axis::X, NamedData<std::reference_wrapper<BoundaryGrid<VectorValue<double, 2>, 2>>>>>, 1>>>>;

	using EE = OpOperatorChain<
		OpOperatorDerivative<2, OpNegIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
		OpOperatorCombination<
		OpOperatorChain<
		OpOperatorDirectionalDerivative<Axis::Y, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
		OpOperatorChain<
		OpIdentity,
		OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::Y, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>,
		OpOperatorChain<
		OpOperatorDirectionalDerivative<Axis::X, 1, OpFractionLiteral<1, 2>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>,
		OpOperatorChain<OpIdentity, OpDerivative<Solver<SolverFT<Stencil2d2h<9, 6, 13>, 0>, 0>::derivative<Axis::X, 1>, OpIdentity, OpTerms<OpIdentity, Term<Variable<0, NamedData<std::reference_wrapper<BoundaryGrid<double, 2>>>>, 1>>, SolverFT<Stencil2d2h<9, 6, 13>, 0>>>>>>;

	auto test_expand_chain = EE{};
	auto test_expand_chain_applied = expr::apply_operators(test_expand_chain).g;

	constexpr bool flagflag = std::is_same_v<decltype(test_expand_chain), decltype(test_expand_chain_applied)>;

	constexpr size_t N = symphas::lib::index_of_value<size_t, 2, 0, 1, 2, 3>;


	Grid<double, 2>* data = new Grid<double, 2>[5] { grid2, grid2, grid2, grid2, grid2 };
	iter_type ind = 0;
	DynamicIndex index(ind);

	auto dynamic_term = expr::make_term_dynamic(index, Variable<0, Grid<double, 2>*>(data));
	auto df = expr::make_operator_functional_derivative(dynamic_term);
	auto free_energy = expr::cellular_fe(solver, dynamic_term);

	auto ss = dynamic_term.eval(0);

	//auto functional = expr::make_domain_integral(2. * free_energy);

	//auto motion = expr::apply_operators(df(functional));


	/*
	auto test_eval_chain = OpChain<OpCoeff<double, expr::symbols::i_<0, 0>>, OpIdentity, 
		OpBinaryMul<
			OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>, 
			OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>
		>>{};
	auto test_eval_chain_eval = expr::eval_type_t<decltype(test_eval_chain)>{};

	auto test_eval_mul = OpBinaryMul<
		OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>,
		OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>>{};
	auto test_eval_mul_eval = expr::eval_type_t<decltype(test_eval_mul)>{};

	auto test_eval_iden_mul = OpIdentity{}(OpBinaryMul<
		OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>,
		OpOperatorChain<OpOperatorDerivative<1, OpIdentity, SolverFT<Stencil2d2h<9, 6, 13>, 0>>, OpTerms<OpIdentity, Term<Variable<0, std::reference_wrapper<NamedData<SymbolicData<SolverSystemFD<double, 2>>>>>, 1>>>>{});
	auto test_eval_iden_mul_eval = expr::eval_type_t<decltype(test_eval_iden_mul)>{};

	auto a = test_eval_iden_mul.a;

	auto test_coeff_eval = OpCoeff<double, expr::symbols::i_<0, 0>>{}.eval(0);
	*/
}

