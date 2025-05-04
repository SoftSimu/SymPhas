#pragma once

#include "procinclude.h"
#include "collectgroup.h"



// **************************************************************************************************

namespace symphas
{
	//! Functions that execute the data collection.
	namespace proc {}
}

/*
 * in the template aliases below, each collector that should be recorded
 * for the particular field type is added
 */

using CollectorGroupScalar = CollectorGroup<scalar_t, COLLECT_SCALAR>;
using CollectorGroupComplex = CollectorGroup<complex_t, COLLECT_COMPLEX>;




namespace symphas::proc
{
	template<typename M>
	bool collect_data_avg(const M* models, const char* dir, size_t runs)
	{
		bool bc = CollectorGroupComplex::aggregate_start(models, dir, runs);
		bool bs = CollectorGroupScalar::aggregate_start(models, dir, runs);
		return bc || bs;
	}
	
	template<typename M>
	bool collect_data(M const& model, const char* dir)
	{
		bool bc = CollectorGroupComplex::collect_start(model, dir);
		bool bs = CollectorGroupScalar::collect_start(model, dir);
		return bc || bs;
	
	}
}


#undef SCALAR
#undef COMPLEX
#undef VECTOR




	