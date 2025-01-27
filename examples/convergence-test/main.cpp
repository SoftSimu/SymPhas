#include <iostream>

#include "symphas.h"



int main(int argc, char* argv[])
{
	symphas::init(argv[1], argv + 2, argc - 2);


	using solver = SolverFT<Stencil2d2h<9, 13, 6>>;
	using model = model_AC_MMS_t<2, solver>;

	char convergence_f_name[] = "convergence_results_dt.csv";
	FILE* f = fopen(convergence_f_name, "w");
	fprintf(f, "dt, N, L2\n");
	fclose(f);


	// for running changing time

	double dt = symphas::conf::config().dt;
	double mdt = 0.6;


	Conf c(symphas::conf::config());
	auto pp = c.get_problem_parameters();
	auto init_data = *pp.get_initial_data();

	double h = pp.get_interval(Axis::X, 0).width();

	symphas::interval_data_type intervals;
	intervals[Axis::X].set_interval(0, 1 - h, h);
	intervals[Axis::Y].set_interval(0, 0.5, h);
	c.set_intervals(intervals, 0);

	symphas::init_data_type compare_init(symphas::init_data_expr("MMS_AC_FINAL", init_data.expr_data.get_coeff(), init_data.expr_data.get_num_coeff()));
	pp.set_initial_data(compare_init, 0);

	char buffer[100];
	sprintf(buffer, "acmms/h=%1.4lf_analytic", h);
	c.set_directory(buffer);

	c.save = SaveParams(0, 0);
	model model_final(c.get_coeff_list(), c.coeff_count(), pp);
	symphas::find_solution(model_final, c);
	symphas::io::write_plot_config(model_final, c.get_result_dir(), c.save);

	auto data_final = model_final.get_field<0>();

	for (iter_type i = 0; i < 10; ++i)
	{
		Conf c(symphas::conf::config());
		c.dt = dt;

		sprintf(buffer, "acmms/h=%1.4lf_dt=%1.4lf_numeric", h, dt);
		c.set_directory(buffer);

		auto init_data = *pp.get_initial_data();
		symphas::init_data_type first_init(symphas::init_data_expr("MMS_AC_INIT", init_data.expr_data.get_coeff(), init_data.expr_data.get_num_coeff()));
		pp.set_initial_data(first_init, 0);

		iter_type N = 8. / dt;
		c.save = SaveParams(N, N);


		model model_simulate(c.get_coeff_list(), c.coeff_count(), pp);
		symphas::find_solution(model_simulate, c);
		symphas::io::write_plot_config(model_simulate, c.get_result_dir());


		auto data_simulate = model_simulate.get_field<0>();

		double norm = 0;
		auto dx = model_simulate.system<0>().get_info().at(Axis::X).width();
		for (iter_type n = 0; n < data_simulate.len; ++n)
		{
			auto diff = data_simulate[n] - data_final[n];
			norm += diff * diff * dx * dx;
		}


		norm = std::sqrt(norm);

		f = fopen(convergence_f_name, "a");
		fprintf(f, "%E, %E, %d, %E\n", h, dt, N, norm);
		fclose(f);

		dt *= mdt;
	}


	double mh = 0.8;
	h = 0.01;

	char convergenceh_f_name[] = "convergence_results_h.csv";
	f = fopen(convergenceh_f_name, "w");
	fprintf(f, "h, sysh, x, y, L2\n");
	fclose(f);

	for (iter_type i = 0; i < 10; ++i)
	{
		Conf c(symphas::conf::config());
		c.set_width(h);

		intervals[Axis::X].set_interval(0, 1 - h, h);
		intervals[Axis::Y].set_interval(0, 0.5, h);
		c.set_intervals(intervals, 0);

		auto pp = c.get_problem_parameters();

		char buffer[100];
		sprintf(buffer, "acmms/h=%1.4lf_numeric", h);
		c.set_directory(buffer);

		auto init_data = *pp.get_initial_data();
		symphas::init_data_type first_init(symphas::init_data_expr("MMS_AC_INIT", init_data.expr_data.get_coeff(), init_data.expr_data.get_num_coeff()));
		pp.set_initial_data(first_init, 0);

		model model_simulate(c.get_coeff_list(), c.coeff_count(), pp);
		symphas::find_solution(model_simulate, c);
		symphas::io::write_plot_config(model_simulate, c.get_result_dir());

		symphas::init_data_type compare_init(symphas::init_data_expr("MMS_AC_FINAL", init_data.expr_data.get_coeff(), init_data.expr_data.get_num_coeff()));
		pp.set_initial_data(compare_init, 0);

		sprintf(buffer, "acmms/h=%1.4lf_analytic", h);
		c.set_directory(buffer);

		c.save = SaveParams(0, 0);
		model model_final(c.get_coeff_list(), c.coeff_count(), pp);
		symphas::find_solution(model_final, c);
		symphas::io::write_plot_config(model_final, c.get_result_dir(), c.save);


		auto data_simulate = model_simulate.get_field<0>();
		auto data_final = model_final.get_field<0>();

		double norm = 0;
		auto dx = model_simulate.system<0>().get_info().at(Axis::X).width();
		for (iter_type n = 0; n < data_simulate.len; ++n)
		{
			auto diff = data_simulate[n] - data_final[n];
			norm += diff * diff * dx * dx;
		}

		norm = std::sqrt(norm);

		f = fopen(convergenceh_f_name, "a");
		fprintf(f, "%E, %E, %d, %d, %E\n", h, model_simulate.system<0>().get_info().at(Axis::X).width(), data_simulate.dims[0], data_simulate.dims[1], norm);
		fclose(f);

		h *= mh;
	}


}


