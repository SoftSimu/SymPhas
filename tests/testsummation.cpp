
#include "testsummation.h"
#include "symphas.h"

void testsummation()
{
	using namespace expr;
	using namespace expr::symbols;
	using ii = i_<0, 0>;
	using jj = i_<1, 0>;
	using v_ii = v_<ii>;
	using v_jj = v_<jj>;


	len_type N = 100;
	len_type sum_n = 6;
	iter_type dims[]{ N, N };


	Grid<double, 2> grid(dims);
	grid::fill_random(grid);

	printf("testing simple sum\n");


	double** values_list = new double* [sum_n];
	for (iter_type i = 0; i < sum_n; ++i)
	{
		values_list[i] = new double[N] {};
		grid::fill_random(values_list[i], N);
	}

	double* sum_result_1 = new double[N] {};
	for (iter_type i = 0; i < sum_n; ++i)
	{
		for (iter_type n = 0; n < N; ++n) 
		{
			sum_result_1[n] += values_list[i][n];
		}
	}

	int ind0;
	DynamicIndex index(&ind0, 0, sum_n - 1);
	auto term = expr::make_term_dynamic(index, values_list);

	double* sum_result_2 = new double[N] {};
	for (iter_type i = 0; i < sum_n; ++i)
	{
		ind0 = i;
		for (iter_type n = 0; n < N; ++n)
		{
			sum_result_2[n] += term[n];
		}
	}

	auto sum_expr_0 = expr::sum<ii>(v_ii{})(expr::series_limits(1, sum_n))(values_list);
	expr::prune::update(sum_expr_0);
	double* sum_result_3 = new double[N] {};
	for (iter_type n = 0; n < N; ++n)
	{
		sum_result_3[n] = sum_expr_0[n];
	}


	bool result_12 = true;
	bool result_23 = true;
	bool result_13 = true;
	for (iter_type n = 0; n < N; ++n)
	{
		result_12 = result_12 && (sum_result_1[n] == sum_result_2[n]);
		result_23 = result_23 && (std::abs(sum_result_2[n] - sum_result_3[n]) < 1e-8);
		result_13 = result_13 && (std::abs(sum_result_1[n] - sum_result_3[n]) < 1e-8);
	}

	printf("12: %s, 23: %s, 13: %s\n",
		(result_12) ? "true" : "false",
		(result_23) ? "true" : "false",
		(result_13) ? "true" : "false");

	printf("testing iterated sums\n");

	for (iter_type n = 0; n < N; ++n)
	{
		sum_result_1[n] = 0;
		sum_result_2[n] = 0;
		sum_result_3[n] = 0;
	}

	for (iter_type i = 0; i < sum_n; ++i)
	{
		for (iter_type j = 0; j < sum_n; ++j)
		{
			if (j != i)
			{
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_1[n] += values_list[i][n] * values_list[j][n];
				}
			}
		}
	}

	int ind1;
	DynamicIndex index2(&ind1, 0, sum_n - 1);
	auto termi = expr::make_term_dynamic(index, values_list);
	auto termj = expr::make_term_dynamic(index2, values_list);
	for (iter_type i = 0; i < sum_n; ++i)
	{
		for (iter_type j = 0; j < sum_n; ++j)
		{
			if (j != i)
			{
				ind0 = i;
				ind1 = j;
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_2[n] += termi[n] * termj[n];
				}
			}
		}
	}

	auto sum_expr_1 = expr::sum(v_ii{} *v_jj{}).select(ii{}, jj{} != ii{})(expr::as_array_arg(values_list, sum_n), expr::as_array_arg(values_list, sum_n));
	expr::prune::update(sum_expr_1);
	for (iter_type n = 0; n < N; ++n)
	{
		sum_result_3[n] = sum_expr_1[n];
	}

	result_12 = true;
	result_23 = true;
	result_13 = true;
	for (iter_type n = 0; n < N; ++n)
	{
		result_12 = result_12 && (sum_result_1[n] == sum_result_2[n]);
		result_23 = result_23 && (std::abs(sum_result_2[n] - sum_result_3[n]) < 1e-8);
		result_13 = result_13 && (std::abs(sum_result_1[n] - sum_result_3[n]) < 1e-8);
	}

	printf("12: %s, 23: %s, 13: %s\n",
		(result_12) ? "true" : "false",
		(result_23) ? "true" : "false",
		(result_13) ? "true" : "false");

	printf("testing iterated sums with dynamic limits\n");

	double** sum_result_11 = new double*[sum_n];
	double** sum_result_22 = new double*[sum_n];
	double** sum_result_33 = new double*[sum_n];
	for (iter_type nn = 0; nn < sum_n; ++nn)
	{
		sum_result_11[nn] = new double[N] {};
		sum_result_22[nn] = new double[N] {};
		sum_result_33[nn] = new double[N] {};
	}

	for (iter_type nn = 0; nn < sum_n; ++nn)
	{
		for (iter_type j = 0; j < sum_n; ++j)
		{
			if (j != nn)
			{
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_11[nn][n] += 2 * values_list[nn][n] * values_list[j][n] * values_list[j][n];
				}
			}
		}

		for (iter_type i = 0; i < sum_n; ++i)
		{
			if (i != nn)
			{
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_11[nn][n] += 2 * values_list[nn][n] * values_list[i][n] * values_list[i][n];
				}
			}
		}
	}


	for (iter_type nn = 0; nn < sum_n; ++nn)
	{
		ind0 = nn;
		for (iter_type j = 0; j < sum_n; ++j)
		{
			if (j != nn)
			{
				ind1 = j;
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_22[nn][n] += 2 * termi[n] * termj[n] * termj[n];
				}
			}
		}

		ind1 = nn;
		for (iter_type i = 0; i < sum_n; ++i)
		{
			if (i != nn)
			{
				ind0 = i;
				for (iter_type n = 0; n < N; ++n)
				{
					sum_result_22[nn][n] += 2 * termj[n] * termi[n] * termi[n];
				}
			}
		}
	}

	using G = DynamicVariable<Grid<double, 2>>;

	Grid<double, 2>* grid_list = new Grid<double, 2>[6]{
		Grid<double, 2>{dims},
		Grid<double, 2>{dims},
		Grid<double, 2>{dims},
		Grid<double, 2>{dims},
		Grid<double, 2>{dims},
		Grid<double, 2>{dims}};

	for (iter_type i = 0; i < sum_n; ++i)
	{
		grid::copy(values_list[i], grid_list[i]);
	}

	auto sum_expr_20 = expr::sum(v_ii{} ).select(ii{})(expr::as_array_arg(grid_list, sum_n)).data.e;
	auto sum_expr_2 = expr::sum(v_ii{} *v_jj{} *v_ii{} *v_jj{}).select(ii{}, jj{} != ii{})(expr::as_array_arg(grid_list, sum_n), expr::as_array_arg(grid_list, sum_n));
	auto df = expr::make_operator_derivative<1>(SymbolicFunctionalDerivative<G>{index});
	auto sum_evaluate_2 = expr::apply_operators(df(expr::make_domain_integral(sum_expr_2)));
	for (iter_type nn = 0; nn < sum_n; ++nn)
	{
		ind0 = nn;
		expr::prune::update(sum_evaluate_2);
		for (iter_type n = 0; n < N; ++n)
		{
			sum_result_33[nn][n] = sum_evaluate_2[n];
		}
	}

	for (iter_type i = 0; i < sum_n; ++i)
	{
		result_12 = true;
		result_23 = true;
		result_13 = true;
		for (iter_type n = 0; n < N; ++n)
		{
			result_12 = result_12 && (sum_result_11[i][n] == sum_result_22[i][n]);
			result_23 = result_23 && (std::abs(sum_result_22[i][n] - sum_result_33[i][n]) < 1e-8);
			result_13 = result_13 && (std::abs(sum_result_11[i][n] - sum_result_33[i][n]) < 1e-8);
		}

		printf("%d -> 12: %s, 23: %s, 13: %s\n",
			i,
			(result_12) ? "true" : "false",
			(result_23) ? "true" : "false",
			(result_13) ? "true" : "false");
	}
}








