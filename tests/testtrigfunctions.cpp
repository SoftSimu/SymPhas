
#include "testtrigfunctions.h"
#include "expressiontypeincludes.h"

#define TEST_FUNCTION(FUNC, LEFT, RIGHT) { \
		constexpr double WIDTH = (RIGHT - LEFT) / TEST_POINTS; \
		auto xx = expr::make_varx(TEST_POINTS, LEFT, RIGHT); \
		(xx * xx * expr::FUNC (xx)).print(stdout); printf("\n"); \
		for (iter_type i = 0; i < TEST_POINTS; ++i) \
		{ \
			if (expr::FUNC(xx).eval(i) != std::FUNC(xx.eval(i))) \
				printf(#FUNC "(%lf) = %lf vs " #FUNC "(%lf) = %lf\n", xx.eval(i), expr::FUNC(xx).eval(i), xx.eval(i), std::FUNC(LEFT + i * WIDTH));\
		} }

bool testtrigfunctions()
{

	// test functions

	
	constexpr size_t TEST_POINTS = 100;


	using namespace symphas;
		TEST_FUNCTION(cos, 0, 2 * PI);
		TEST_FUNCTION(sin, 0, 2 * PI);
		TEST_FUNCTION(tan, 0, 2 * PI);
		TEST_FUNCTION(acos, -1, 1);
		TEST_FUNCTION(asin, -1, 1);
		TEST_FUNCTION(atan, -1, 1);
		TEST_FUNCTION(cosh, 0, 2 * PI);
		TEST_FUNCTION(sinh, 0, 2 * PI);
		TEST_FUNCTION(tanh, 0, 2 * PI);
		TEST_FUNCTION(acosh, -1, 1);
		TEST_FUNCTION(asinh, -1, 1);
		TEST_FUNCTION(atanh, -1, 1);
		TEST_FUNCTION(sqrt, 0, 100);

		
}
