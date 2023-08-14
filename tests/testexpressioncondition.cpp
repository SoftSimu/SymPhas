
#include "testexpressioncondition.h"
#include "symphas.h"

void testexpressioncondition()
{
	using namespace expr;
	using namespace expr::symbols;
	using namespace expr::split;
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


	int ind0;
	DynamicIndex index(&ind0, 0, sum_n - 1);
	auto term = expr::make_term_dynamic(index, values_list);

	auto sum_expr_0 = expr::sum<ii>(v_ii{})(expr::series_limits(1, sum_n))(values_list);


	auto split0 = filter<matches_series>(sum_expr_0);
	auto split1 = filter<matches_series>(one);
	auto split2 = filter<matches_series, matches_mul>(sum_expr_0 + one + expr::make_mul(one, one));
}








