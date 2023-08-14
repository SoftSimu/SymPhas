
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
	 
	//auto op = expr::make_operator_derivative<2>(solver);
	//auto opg = expr::make_operator_derivative<1>(solver);
	auto term = expr::make_term<1>(grid);
	//auto serm = expr::make_term<2>(grid2);
	//auto term2 = term * term;
	//auto termp = expr::pow<2>(term);
	//using ev2 = expr::eval_type_t<decltype(term2)>;
	//using ev2p = expr::eval_type_t<decltype(termp)>;
	//auto test = (term * op) * term;
	//auto app = expr::apply_operators(test);
	//
	//auto flo = op * (expr::pow<2>(1. + op) * serm);
	//flo.print(stdout);
	//printf("\n");
	//auto floa = expr::apply_operators(flo);
	//floa.print(stdout);
	//printf("\n");
	//
	//auto flo1 = expr::pow<2>(term * opg) * term;
	//flo1.print(stdout);
	//printf("\n");
	//auto flo1a = expr::apply_operators(flo1);
	//flo1a.print(stdout);
	//printf("\n");
	//
	//auto flo2 = expr::transpose(term) * opg(term);
	//flo2.print(stdout);
	//printf("\n");
	//auto flo2a = expr::apply_operators(flo2);
	//flo2a.print(stdout);
	//printf("\n");
	//
	//
	//auto div_flo3 = 1. / (1. + serm);
	//auto chain_flo3 = term * opg;
	//auto flo3 = (div_flo3 * chain_flo3) * term;
	//flo3.print(stdout);
	//printf("\n");
	//auto flo3a = expr::apply_operators(flo3);
	//flo3a.print(stdout);
	//printf("\n");

	auto term_tensor = OpTensor<OpLiteral<double>, 0, 2>{} *(term * term);
	auto deriv_tensor = expr::make_derivative< SolverFT<Stencil2d2h< 5, 6, 13>>::derivative<Axis::X, 1>>(term_tensor, solver);

	double* dt = new double(0.05);
	double* h_ptr = &h;

	auto e = expr::make_literal(2.) * expr::make_unit_vector<2>(2 * expr::symbols::Pi * (NoiseData<expr::NoiseType::POISSON, scalar_t, 2>(grid::dim_list(), h_ptr, dt)(OpIdentity{})));

	auto e0 = expr::get<0>(e);


	//symphas::internal::GeneratedStencilApply stencil(expr::get_central_space_mixed_stencil<2>(std::index_sequence<0, 2>{}));


	auto termx = expr::make_row_vector<1, 2>() * term;
	auto l2 = OpTensor<OpFractionLiteral<2, 1>, 1, 1, 2, 2>{} *expr::make_derivative<Solver<SolverFT<Stencil2d2h<5, 6, 13>>>::mixed_derivative<0, 2>>(termx, solver);
	auto l2_ = l2.eval(0);


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
}

