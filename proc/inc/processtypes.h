#pragma once

#include "io.h"
#include "datalib.h"

/* enums for the data processing utility
 */

enum class ProcessType
{
	NO, SCALAR, VECTOR, POINT, DYNAMIC
};

enum class ProcessTag
{
	NONE, DFT
};

struct DataParams
{
	DataParams() : type{ ProcessType::NO }, tag{ ProcessTag::NONE }, save{ SaveParams{} } {}

	ProcessType type;			// Flag for collecting structure factor.
	ProcessTag tag;				// A tag on how to collect the data.
	SaveParams save;			// Save type for this data.
};

namespace symphas::internal
{
	void parse_proc_job(const char* v, DataParams* d, ProcessType default_proc, char split_char);
	void parse_proc_spec(const char* str, DataParams* d, ProcessType default_proc);
	void parse_proc(const char* proc_name, char split_char);
	void update_data_stop(iter_type stop);

	inline std::map<std::string, std::vector<DataParams>*, symphas::internal::any_case_comparator> data_map;
	inline std::map<std::string, ProcessType, symphas::internal::any_case_comparator> default_procs_map;

	std::map<std::string, std::vector<DataParams>*, symphas::internal::any_case_comparator> get_data_map();
	std::map<std::string, ProcessType, symphas::internal::any_case_comparator> get_default_procs_map();

}

/*
 * alias templates for the vector and scalar types emitted by the processing functions
 */

template<typename Y, size_t D>
using vector_data = symphas::FieldAxis<D, Y*>;

template<typename Y>
using scalar_data = symphas::Field<scalar_t*, Y*>;

template<typename Y>
using point_data = symphas::Field<scalar_t, Y>;


template<typename F>
using dynamic_data = symphas::Field<iter_type*, F**>;

