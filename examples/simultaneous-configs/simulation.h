#pragma once

#include "timer.h"
#include "symphas.h"


 // **************************************************************************************

using namespace symphas;

namespace simulate
{

#ifdef USING_CONF
	template<typename M>
	int simulate(double const *coeff, size_t num_coeff)
	{
		auto pp = symphas::conf::config().get_problem_parameters();

		if (params::plots_only)
		{
			M model(coeff, num_coeff, pp);
			symphas::io::write_plot_config(model);
		}
		else
		{
			/*
			 * execute a single model and collect data if there is only a single run
			 */
			if (symphas::conf::config().runs == 1)
			{
				M model(coeff, num_coeff, pp);
				symphas::io::write_plot_config(model);
				find_solution(model);
			}

			/*
			 * assemble a vector of duplicated models if there are
			 * multiple runs
			 */
			else
			{
				len_type runs = static_cast<len_type>(symphas::conf::config().runs);
				std::vector<M> models(runs, { coeff, num_coeff, pp });
				symphas::io::write_plot_config(models.front());
				find_solution(models.data(), runs);
			}

		}
		return 1;
	}

#else

	template<typename M>
	int simulate(double const* coeff, size_t num_coeff)
	{
		problem_parameters_type pp{ 1 };
		M model(coeff, num_coeff, pp);
		find_solution(model, 0.05, 100);
		return 1;
	}



#endif


	inline void initiate(const char* modelname, double const *coeff, size_t num_coeff)
	{
#ifdef USING_MODEL_SELECTION
#ifdef USING_CONF
		model_select m{ symphas::conf::config().dimension, symphas::conf::config().stp };
#else
		model_select m{ 2, StencilParams{} };
#endif

		if (m.call<SOLVER_TYPE>(modelname, coeff, num_coeff) == INVALID_MODEL)
		{
			fprintf(SYMPHAS_ERR, "Unknown model provided, '%s'\n", modelname);
			exit(101);
		}

#ifdef PRINT_TIMINGS
		print_timings(SYMPHAS_LOG);
#endif
#endif

	}

}

