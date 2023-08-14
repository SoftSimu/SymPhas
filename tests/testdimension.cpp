
#include "testdimension.h"


void testdimension()
{
	srand((unsigned int)time(NULL));

	iter_type
		L = 10,
		M = 15,
		N = 8;
	
	std::vector<std::pair<axis_2d_type, double>> test_2d;

	// i'm filling it in a weird way so that i can be sure it gives me the right dimensions
	for (iter_type i = 0; i < L; ++i)
	{
		for (iter_type j = 0; j < M; ++j)
		{
			test_2d.emplace_back(axis_2d_type{ i / 2.0, j / 2.0 }, rand() / RAND_MAX);
		}
	}


	std::vector<std::pair<axis_3d_type, double>> test_3d;

	// i'm filling it in a weird way so that i can be sure it gives me the right dimensions
	for (iter_type i = 0; i < L; ++i)
	{
		for (iter_type k = 0; k < N; ++k)
		{
			for (iter_type j = 0; j < M; ++j)
			{
				test_3d.emplace_back(axis_3d_type{ i / 2.0, j / 2.0, k / 2.0 }, rand() / RAND_MAX);
			}
		}
	}

	auto[L2, M2] = symphas::lib::get_dimensions<2>(test_2d)._2();
	auto[L3, M3, N3] = symphas::lib::get_dimensions<3>(test_3d)._3();

	printf("Original values: %d %d %d\n", L, M, N);
	printf("2d test values: %d %d\n", L2, M2);
	printf("3d test values: %d %d %d\n", L3, M3, N3);
}

