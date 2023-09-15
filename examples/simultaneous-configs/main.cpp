#include <iostream>

#include "symphas.h"

#if (defined(MODEL_INCLUDE_HEADER) && defined(USING_MODEL_SELECTION))
#include "simulation.h"

#else

#define dpsi dop(1)
#define psi op(1)
MODEL(EX, (SCALAR),
	EVOLUTION(
		dpsi = lap(psi) + (c1 - lit(4.) * c2 * psi * psi) * psi)
)
//LINK_WITH_NAME(MA, MODELA)
#endif


int main(int argc, char* argv[])
{

#ifdef UNIT_TEST

	testcell();
	testdimension();
	testnextsave();
	testexpression();

#else

#	if (!defined(MODEL_INCLUDE_HEADER) || !defined(USING_MODEL_SELECTION))

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
	pp.set_time_step(dt);

#		ifdef SOLVER_INCLUDE_HEADER
	model_EX_t<2, SolverFT<Stencil2d2h<5, 9, 6>>> model{ pp };
	symphas::find_solution(model, dt, 50);
#		else
	model_EX_t<2, GenericSolver<>> model{ pp };
#		endif

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
	simulate::initiate(symphas::conf::config().get_model_name(), symphas::conf::config().get_coeff_list(), symphas::conf::config().get_coeff_len());
#	else
	simulate::initiate("MODELA", nullptr, 0);
#	endif
#	endif

#endif

}


