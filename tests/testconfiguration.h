#pragma once

#include "symphas.h"

void testconfiguration()
{

	double dt = 0.1;
	symphas::problem_parameters_type pp{ 1 };
	auto f = [](iter_type, len_type const*, size_t) -> scalar_t { return 1; };
	symphas::init_data_functor init_functor{ f };
	
	
	
	symphas::b_data_type bdata;
	symphas::init_data_type tdata1{ { Axis::NONE, f } };
	symphas::init_data_type tdata{ { Axis::NONE, Inside::UNIFORM, { -1, 1 } } };
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
	model_MA_t<2, SOLVER_TYPE<Stencil2d2h<5, 9, 6>>> model{ pp };
	auto data = model.grid<0>();

	auto solver = SOLVER_TYPE<Stencil2d2h<5, 9, 6>>::make_solver( pp);
	auto op = expr::make_term(Grid<double, 2>({ 10, 10 }));
	auto d = expr::make_operator_derivative<2>(solver);
	auto combination = (expr::make_literal(3.)^2) + d;
	auto ec = combination * (op * op);
	
	len_type dims[] = { 10, 10 };
	double h[] = { 1, 1 };

}


