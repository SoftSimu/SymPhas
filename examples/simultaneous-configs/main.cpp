#include <iostream>

#ifdef MODEL_INCLUDE_HEADER
#include "simulation.h"

#else
#include "symphas.h"

#define dpsi dop(1)
#define psi op(1)
MODEL(EX, (SCALAR),
	MODEL_DEF(
		dpsi = lap(psi) + (c1 - lit(4.) * c2 * psi * psi) * psi)
)
//LINK_WITH_NAME(MA, MODELA)
#endif

#ifdef UNIT_TEST

#include "testnextsavedefines.h"
#include "testcell.h"
#include "testdimension.h"
#include "testexpression.h"

#endif

struct full_fill
{
	scalar_t constant;
	full_fill(scalar_t constant) : constant{ constant } {}
	scalar_t operator()(iter_type, len_type const*, size_t)
	{
		return constant;
	}
};


struct half_fill
{
	scalar_t constant;
	half_fill(scalar_t constant) : constant{ constant } {}
	scalar_t operator()(iter_type n, len_type const* dims, size_t dimension)
	{
		if (dimension == 2)
		{
			axis_coord_t x = n % dims[0];
			axis_coord_t y = n / dims[0];
			return (x < dims[0] / 2) ? constant : 0;
		}
		else if (dimension == 3)
		{
			axis_coord_t x = n % dims[0];
			axis_coord_t y = (n / dims[0]) % dims[1];
			axis_coord_t z = n / (dims[0] * dims[1]);
			return (x < dims[0] / 2) ? constant : 0;
		}
		return 0;
	}
};


struct quarter_fill
{
	scalar_t constant;
	quarter_fill(scalar_t constant) : constant{ constant } {}
	scalar_t operator()(iter_type n, len_type const* dims, size_t dimension)
	{
		if (dimension == 2)
		{
			axis_coord_t x = n % dims[0];
			axis_coord_t y = n / dims[0];
			return (
				x < dims[0] / 2 
				&& y < dims[1] / 2) 
				? constant : 0;
		}
		else if (dimension == 3)
		{
			axis_coord_t x = n % dims[0];
			axis_coord_t y = (n / dims[0]) % dims[1];
			axis_coord_t z = n / (dims[0] * dims[1]);
			return (
				x < dims[0] / 2
				&& y < dims[1] / 2
				&& z < dims[2] / 2)
				? constant : 0;
		}
		return 0;
	}
};

int main(int argc, char* argv[])
{
	//MPI_Init(&argc, &argv);

#ifdef UNIT_TEST

	testcell();
	testdimension();
	testnextsave();
	testexpression();

#else

#	if !defined(MODEL_INCLUDE_HEADER) || !defined(USING_MODEL_SELECTION)

	UNUSED(argc);
	UNUSED(argv);

	double dt = 0.1;
	symphas::problem_parameters_type pp{ 1 };

	symphas::b_data_type bdata;
	symphas::init_data_type tdata{ Inside::UNIFORM, { -1, 1 } };
	symphas::interval_data_type vdata;
	symphas::interval_element_type interval;
	interval.set_interval(1, 64, 0.5);

	bdata[Side::LEFT] = BoundaryType::PERIODIC;
	bdata[Side::RIGHT] = BoundaryType::PERIODIC;
	bdata[Side::TOP] = BoundaryType::PERIODIC;
	bdata[Side::BOTTOM] = BoundaryType::PERIODIC;

	vdata[Axis::X] = interval;
	vdata[Axis::Y] = interval;

	pp.set_boundary_data(&bdata);
	pp.set_initial_data(&tdata);
	pp.set_interval_data(&vdata);
	pp.set_problem_time_step(dt);

	model_EX_t<2, SOLVER_TYPE<Stencil2d2h<5, 9, 6>>> model{ pp };
	symphas::find_solution(model, dt, 50);
	auto pfdata = model.grid<0>();

#		ifdef USING_IO
	symphas::io::save_grid(pfdata.values, pfdata.dims, 2);
#		endif

#	else


	Time t("entire simulation");

	if (argc == 1)
	{
		fprintf(SYMPHAS_ERR, "the first argument specifying the configuration must be provided\n");
		exit(102);
	}

	

	symphas::init(argv[1], argv + 2, argc - 2);


#	ifdef USING_CONF
	simulate::initiate(symphas::conf::config().get_model_name(), symphas::conf::config().get_coeff_list(), symphas::conf::config().coeff_count());
#	else
	//simulate::initiate("MODELA", nullptr, 0);
#	endif
#	endif

#endif

}


