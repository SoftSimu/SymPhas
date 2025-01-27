
#include "expressioncuda.cuh"
#include "simulation.h"
#include "symphas.h"

#ifndef _MSC_VER
#include <unistd.h>
#endif

int main(int argc, char* argv[]) {
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
