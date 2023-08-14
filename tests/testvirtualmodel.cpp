
#include "testvirtualmodel.h"

#include "symphas.h"

void testvirtualmodel()
{
	// These two files must exist and be in checkpoint format.
	const char* first_file = "tests\\data0";
	const char* second_file = "tests\\data02";

	auto pp = symphas::conf::config().get_problem_parameters();

	symphas::init_data_type tdataf{ { Axis::NONE, { Inside::FILE, symphas::init_data_read{ first_file, 0 } } } };
	pp.set_initial_data(&tdataf);

	SaveParams pyma_save{ 100, 1000 };
	ModelVirtual<2, scalar_t> pyma("load_dir", pyma_save);
	Conf c(symphas::conf::config());
	c.set_directory("pyma");
	c.save = pyma_save;
	symphas::find_solution(pyma, c);


	symphas::init_data_type tdata{ { Axis::NONE, { Inside::FILE, { second_file, 1000 } } } };
	System<scalar_t, 2> sys(tdata, pp.get_interval_data()[0]);

	auto mypy_end = sys.as_field();
	auto sim_end = pyma.template system<0>().as_field();
	symphas::io::save_vec(symphas::lib::abs(mypy_end - sim_end));

}

