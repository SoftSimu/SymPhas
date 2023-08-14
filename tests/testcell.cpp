
#include "testcell.h"


void testcell()
{
	any_vector_t<double, 2> z{ 1, 1 };
	any_vector_t<double, 2> u{ 5, 2 };
	any_vector_t<double, 2> v{ 3, -1.5 };


	double a = 3;

	/*
	 * dot product, interchanging operands
	 */
	auto uv = u * v;
	auto vu = v * u;
	printf("dot product...\n");
	printf("(%f, %f) * (%f, %f) = %f\n", u.v[0], u.v[1], v.v[0], v.v[1], uv);
	printf("(%f, %f) * (%f, %f) = %f\n", v.v[0], v.v[1], u.v[0], u.v[1], vu);

	/*
	 * scalar multiplication with scalar on left and right sides
	 */
	auto au = a * u;
	auto ua = u * a;
	printf("scalar multiplication...\n");
	printf("(%f, %f) * %f = (%f, %f)\n", u.v[0], u.v[1], a, ua.v[0], ua.v[1]);
	printf("%f * (%f, %f) = (%f, %f)\n", a, u.v[0], u.v[1], au.v[0], au.v[1]);

	// identity
	auto vz = u * z;
	//auto uu = u * (-u);
	printf("identity...\n");
	printf("(%f, %f) * (%f, %f) = %f\n", u.v[0], u.v[1], z.v[0], z.v[1], vz);
	//printf("(%f, %f) * (%f, %f) = (%f, %f)\n", u.v[0], u.v[1], (-u).v[0], (-u).v[1], uu[0], uu[1]);


	// zero scalar
	auto zero = u * 0.;
	printf("zero scalar...\n");
	printf("(%f, %f) * 0 = (%f, %f)\n", u.v[0], u.v[1], zero.v[0], zero.v[1]);
}