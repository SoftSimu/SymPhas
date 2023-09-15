#include "symphas.h"


int  main(int argc , char* argv [])
{
	UNUSED(argc);
	UNUSED(argv);

	// configure program-wide parameters

	using namespace params;

	// Configure single input and single output.
	PARAMS += SINGLE_INPUT << true, SINGLE_OUTPUT << true;

	// Configure that the writer will output in CSV format.
	PARAMS += WRITER << symphas::IOType::CSV;

	// Configure to turn on VTK and set the interval of visualization one update every 11 frames.
	PARAMS += VIZ_INTERVAL << 11;

	symphas::init();
}

