
#include "testdft.h"

#include "symphas.h"

void testdft()
{

	symphas::init_data_type tdata{ { Axis::NONE, { Inside::UNIFORM, symphas::init_data_parameters{ -1, 1 } } }};
	symphas::interval_data_type vdata;
	symphas::interval_element_type interval;
	interval.set_count(1, 16);
	vdata[Axis::X] = interval;
	interval.set_count(1, 22);
	vdata[Axis::Y] = interval;
	interval.set_count(1, 13);
	vdata[Axis::Z] = interval;

	System<scalar_t, 3> sys(tdata, vdata);
	System<complex_t, 3> sysc(tdata, vdata);
	System<complex_t, 3> sysc1(tdata, vdata);

	auto f = sysc.as_field();
	auto dft0 = symphas::lib::fourier_transform(f);

	symphas::dft::long_dft(sysc.values, sysc1.values, sysc.dims[0], sysc.dims[1]);
	symphas::dft::long_dft(sysc1.values, sysc.values, sysc.dims[0], sysc.dims[1], true);

	auto dft1 = sysc.as_field();

	symphas::io::save_vec(f, "original");
	symphas::io::save_vec((1.0 / sysc.len) * dft1, "long");
	//symphas::io::save_vec((1.0 / sysc.len) * dft0);
	symphas::io::save_vec((1.0 / sysc.len) * symphas::lib::fourier_transform(dft0, true));



}

