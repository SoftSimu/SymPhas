#include "symphas.h"

#define psi op(1)
#define dpsi dop(1)

MODEL(EX, (SCALAR),
	EVOLUTION(dpsi = lap(psi) + (c(1) - c(2) * psi * psi) * psi)
)


int  main(int argc , char* argv [])
{
	UNUSED(argc);
	UNUSED(argv);

	double dt = 0.1;
	symphas::problem_parameters_type pp{ 1 };
	symphas::interval_element_type interval;
	interval.set_count(0, 80, 128);

	symphas::b_data_type bdata(2, BoundaryType::PERIODIC);
	symphas::interval_data_type vdata(2, interval);
	symphas::init_data_type tdata(Inside::UNIFORM, { -1, 1 });

	bdata[Side::TOP] = BoundaryTag::CONSTANT << 1;
	bdata[Side::BOTTOM] = BoundaryTag::CONSTANT << 1;

	pp.set_boundary_data(&bdata);
	pp.set_initial_data(&tdata);
	pp.set_interval_data(&vdata);
	pp.set_time_step(dt);


	model_EX_t<2, SolverFT<Stencil2d2h<>>> model{ pp };
	symphas::find_solution(model, dt, 50);
	auto pfdata = model.grid<0>();
}
