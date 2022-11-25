#include "symphas.h"

#define psi op(1)
#define dpsi dop(1)
#define A c1
#define B c2

MODEL(EX, (SCALAR),
	EVOLUTION(dpsi = lap(psi) + (A - lit (4.) * B * psi * psi) * psi)
)


int  main(int argc , char* argv [])
{
	UNUSED(argc);
	UNUSED(argv);

	double dt = 0.1;
	symphas::problem_parameters_type pp{ 1 };

	symphas::b_data_type bdata;
	symphas::init_data_type tdata{ Inside::UNIFORM, { -1, 1 } };
	symphas::interval_data_type vdata;
	symphas::interval_element_type interval;
	interval.set_interval_count(0, 128, 128);

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


	model_EX_t<2, SolverSP> model{ pp };
	symphas::find_solution(model, dt, 50);
	auto pfdata = model.grid<0>();
}
