
#include "symphas.h"
#include <chrono>
#include <cstdio>

#ifdef MODEL_INCLUDE_HEADER
#include "simulation.h"
#endif

#ifndef _MSC_VER
#include <unistd.h>
#endif

namespace {
struct PhaseTimer {
  using clk = std::chrono::high_resolution_clock;
  clk::time_point t0;
  PhaseTimer() : t0(clk::now()) {}
  double elapsed_ms() const {
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
  }
};
}

int main(int argc, char* argv[]) {
  PhaseTimer t_main;
#ifdef MODEL_INCLUDE_HEADER

  Time t("entire simulation");

  PhaseTimer t_init;
  symphas::init(argv[1], argv + 2, argc - 2);
  fprintf(stderr, "[PHASE rank=? symphas_init_ms=%.1f]\n", t_init.elapsed_ms());

#ifdef USING_CONF
  PhaseTimer t_initiate;
  initiate(symphas::conf::config().model_settings.model,
           symphas::conf::config().model_settings.coeff,
           symphas::conf::config().model_settings.coeff_len);
  fprintf(stderr, "[PHASE rank=? initiate_ms=%.1f]\n",
          t_initiate.elapsed_ms());
#else
  initiate("MODELA", nullptr, 0);
#endif

  PhaseTimer t_finalize;
  symphas::finalize();
  fprintf(stderr, "[PHASE rank=? finalize_ms=%.1f]\n",
          t_finalize.elapsed_ms());
  fprintf(stderr, "[PHASE rank=? total_main_ms=%.1f]\n",
          t_main.elapsed_ms());

#else
  printf("Nothing to do.\n");
#endif
}

