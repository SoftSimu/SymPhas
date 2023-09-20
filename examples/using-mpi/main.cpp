#include <iostream>

#include "symphas.h"
#include "solverft.h"
#include "modelcellnomot.h"


int main(int argc, char* argv[])
{
	symphas::init(argv[1], argv + 2, argc - 2);
	MPI_Init(&argc, &argv);
	
	model_CELL_MIGRATION_NO_MOTILITY_t<2, SolverFT<Stencil2d2h<>>> model(
		symphas::conf::config().get_coeff_list(), 
		symphas::conf::config().get_coeff_len(), 
		symphas::conf::config().get_problem_parameters());

	symphas::io::write_plot_config(model);
	symphas::find_solution(model);
}


