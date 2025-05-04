#pragma once

#include "data.h"
#include "savedefines.h"


/*
 * the collection process is used for computing static information about the system and
 * immediately saving them to a file
 */
template<typename X, typename Y, typename S, typename Dt>
struct ProcessStatic : Process<ProcessStatic, X, Y, S, Dt>
{
	using parent_type = Process<::ProcessStatic, X, Y, S, Dt>;

	// inherit the constructor
	using parent_type::parent_type;



	template<size_t I, typename R, typename M = int>
	void collect_one(point_data<R> const& data, M&&, const char* dir, iter_type index, iter_type)
	{
		symphas::io::save_pts<S>(symphas::lib::combine_data(data.data_x(), data.data_y(), data.length()), dir, index, I);
	}

	template<size_t I, typename R, typename M = int>
	void collect_one(point_data<R*> const& data, M&&, const char* dir, iter_type index, iter_type)
	{
		symphas::io::save_pts<S>(symphas::lib::combine_data(data.data_x(), data.y, data.length()), dir, index, I);
	}

	template<size_t I, typename R, typename M = int>
	void collect_one(scalar_data<R> const& data, M&&, const char* dir, iter_type index, iter_type)
	{
		symphas::io::save_abs<S>(symphas::lib::combine_data(data.x, data.y, data.length()), dir, index, I);
	}

	template<size_t I, typename R, size_t D, typename M = int>
	void collect_one(vector_data<R, D> const& data, M&&, const char* dir, iter_type index, iter_type)
	{
		symphas::io::save_vec<S>(symphas::lib::combine_data(data.x, data.y, data.length()), dir, index, I);
	}


};



