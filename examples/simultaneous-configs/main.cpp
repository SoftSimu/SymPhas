#include <iostream>

#include "symphas.h"

#if (defined(MODEL_INCLUDE_HEADER) && defined(USING_MODEL_SELECTION))
#include "simulation.h"
#endif


int main(int argc, char* argv[])
{
#if (defined(MODEL_INCLUDE_HEADER) || !defined(USING_MODEL_SELECTION))

    symphas::Time t("entire simulation");
	symphas::init(argv[1], argv + 2, argc - 2);

#	ifdef USING_CONF
	simulate::initiate(symphas::conf::config().get_model_name(), symphas::conf::config().get_coeff_list(), symphas::conf::config().get_coeff_len());
#	else
	simulate::initiate("MODELA", nullptr, 0);
#	endif

#endif


}


