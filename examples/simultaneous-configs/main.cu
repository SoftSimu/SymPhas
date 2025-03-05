
#include "symphas.cuh"

#ifndef _MSC_VER
#include <unistd.h>
#endif

using namespace symphas;

template <typename M>
struct Simulation {
#ifdef USING_CONF
  int simulate(double const *coeff, size_t num_coeff) {
    auto pp = symphas::conf::config().get_problem_parameters();

    if (params::plots_only) {
      M model(coeff, num_coeff, pp);
      symphas::io::write_plot_config(model);
    } else {
      /*
       * execute a single model and collect data if there is only a single run
       */
      if (symphas::conf::config().runs == 1) {
        M model(coeff, num_coeff, pp);
        symphas::io::write_plot_config(model);
        find_solution(model);
      }

      /*
       * assemble a vector of duplicated models if there are
       * multiple runs
       */
      else {
        len_type runs = static_cast<len_type>(symphas::conf::config().runs);
        std::vector<M> models(runs, {coeff, num_coeff, pp});
        symphas::io::write_plot_config(models.front());
        find_solution(models.data(), runs);
      }
    }
    return 1;
  }

#else

  int simulate(double const *coeff, size_t num_coeff) {
    problem_parameters_type pp{1};
    M model(coeff, num_coeff, pp);
    find_solution(model, 0.05, 100);
    return 1;
  }

#endif

  template <typename... Ts>
  auto operator()(Ts &&...ts) {
    return simulate(std::forward<Ts>(ts)...);
  }
};

inline void initiate(const char *modelname, double const *coeff,
                     size_t num_coeff) {
#ifdef USING_CONF
  model_select<Simulation> m{
      symphas::conf::config().simulation_settings.dimension,
      symphas::conf::config().simulation_settings.stp};
#else
  model_select<Simulation> m{2, StencilParams{2, 9, 6, 13}};
#endif

  if (m.call<SolverFT>(modelname, coeff, num_coeff) == INVALID_MODEL) {
    fprintf(SYMPHAS_ERR, "Unknown model provided, '%s'\n", modelname);
    exit(101);
  }

#ifdef PRINT_TIMINGS
  print_timings(SYMPHAS_LOG);
#endif
}

int main(int argc, char *argv[]) {
#ifdef MODEL_INCLUDE_HEADER

  Time t("entire simulation");

  symphas::init(argv[1], argv + 2, argc - 2);

#ifdef USING_CONF
  initiate(symphas::conf::config().model_settings.model,
           symphas::conf::config().model_settings.coeff,
           symphas::conf::config().model_settings.coeff_len);
#else
  initiate("MODELA", nullptr, 0);
#endif

  symphas::finalize();
#else
  printf("Nothing to do.\n");
#endif
}
