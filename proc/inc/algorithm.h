
#pragma once

#include <functional>
#include <numeric>

#include "processtypes.h"



namespace symphas::internal
{
	inline axis_2d_type& operator+=(axis_2d_type& a, axis_2d_type const& b)
	{
		a[0] += b[0];
		a[1] += b[1];
		return a;
	}

	inline axis_3d_type& operator+=(axis_3d_type& a, axis_3d_type const& b)
	{
		a[0] += b[0];
		a[1] += b[1];
		a[2] += b[2];
		return a;
	}

	inline axis_2d_type& operator*=(axis_2d_type& a, scalar_t const& b)
	{
		a[0] *= b;
		a[1] *= b;
		return a;
	}

	inline axis_3d_type& operator*=(axis_3d_type& a, scalar_t const& b)
	{
		a[0] *= b;
		a[1] *= b;
		a[2] *= b;
		return a;
	}

	template<typename T>
	static void aggregate_add_to(T& data0, T const& data)
	{
		data0 += data;
	}


	template<typename T>
	static void aggregate_add_to(std::vector<T>& data0, std::vector<T> const& data);
	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X*, Y*>& data0, symphas::Field<X*, Y*> const& data);
	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X, Y*>& data0, symphas::Field<X, Y*> const& data);
	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X, Y>& data0, symphas::Field<X, Y> const& data);
	template<size_t D, typename Y>
	static void aggregate_add_to(symphas::FieldAxis<D, Y>& data0, symphas::FieldAxis<D, Y> const& data);
	template<typename T1, typename T2>
	static void aggregate_add_to(std::pair<T1, T2>& data0, std::pair<T1, T2> const& data);
	template<typename... Ts, size_t... Ns>
	static void aggregate_add_to(std::tuple<Ts...>& data0, std::tuple<Ts...> const& data, std::index_sequence<Ns...>);
	template<typename... Ts>
	static void aggregate_add_to(std::tuple<Ts...>& data0, std::tuple<Ts...> const& data);

	template<typename T>
	static void aggregate_add_to(std::vector<T>& data0, std::vector<T> const& data)
	{
		for (iter_type i = 0; i < data0.size(); ++i)
		{
			aggregate_add_to(data0[i], data[i]);
		}
	}

	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X*, Y*>& data0, symphas::Field<X*, Y*> const& data)
	{
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_add_to(data0.data_x()[i], data.data_x()[i]);
			aggregate_add_to(data0.data_y()[i], data.data_y()[i]);
		}
	}

	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X, Y*>& data0, symphas::Field<X, Y*> const& data)
	{
		aggregate_add_to(data0.data_x(), data.data_x());
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_add_to(data0.data_y()[i], data.data_y()[i]);
		}
	}

	template<typename X, typename Y>
	static void aggregate_add_to(symphas::Field<X, Y>& data0, symphas::Field<X, Y> const& data)
	{
		aggregate_add_to(data0.data_x(), data.data_x());
		aggregate_add_to(data0.data_y(), data.data_y());
	}

	template<size_t D, typename Y>
	static void aggregate_add_to(symphas::FieldAxis<D, Y>& data0, symphas::FieldAxis<D, Y> const& data)
	{
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_add_to(data0.data_y()[i], data.data_y()[i]);
		}
	}

	template<typename T1, typename T2>
	static void aggregate_add_to(std::pair<T1, T2>& data0, std::pair<T1, T2> const& data)
	{
		aggregate_add_to(data0.first, data.first);
		aggregate_add_to(data0.second, data.second);
	}

	template<typename... Ts, size_t... Ns>
	static void aggregate_add_to(std::tuple<Ts...>& data0, std::tuple<Ts...> const& data, std::index_sequence<Ns...>)
	{
		((aggregate_add_to(std::get<Ns>(data0), std::get<Ns>(data)), ...));
	}

	template<typename... Ts>
	static void aggregate_add_to(std::tuple<Ts...>& data0, std::tuple<Ts...> const& data)
	{
		aggregate_add_to(data0, data, std::make_index_sequence<sizeof...(Ts)>{});
	}

	

	template<typename T>
	static void aggregate_avg(std::vector<T>& data0, size_t count);
	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X*, Y*>& data0, size_t count);
	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X, Y*>& data0, size_t count);
	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X, Y>& data0, size_t count);
	template<size_t D, typename Y>
	static void aggregate_avg(symphas::FieldAxis<D, Y>& data0, size_t count);
	template<typename T1, typename T2>
	static void aggregate_avg(std::pair<T1, T2>& data0, size_t count);
	template<typename... Ts, size_t... Ns>
	static void aggregate_avg(std::tuple<Ts...>& data0, size_t count, std::index_sequence<Ns...>);
	template<typename... Ts>
	static void aggregate_avg(std::tuple<Ts...>& data0, size_t count);


	template<typename T>
	static void aggregate_avg(T& data0, size_t count)
	{
		data0 *= (1.0 / count);
	}

	template<typename T>
	static void aggregate_avg(Block<T>& data0, size_t count)
	{
		std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			data0.values, data0.values + data0.len, 
			[&] (auto& value) { aggregate_avg(value, count); });
	}

	template<typename T>
	static void aggregate_avg(std::vector<T>& data0, size_t count)
	{
		std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			data0.begin(), data0.end(), 
			[&] (auto& value) { aggregate_avg(value, count); });
	}

	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X*, Y*>& data0, size_t count)
	{
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_avg(data0.x[i], count);
			aggregate_avg(data0.y[i], count);
		}
	}

	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X, Y*>& data0, size_t count)
	{
		aggregate_avg(data0.data_x(), count);
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_avg(data0.y[i], count);
		}
	}

	template<typename X, typename Y>
	static void aggregate_avg(symphas::Field<X, Y>& data0, size_t count)
	{
		aggregate_avg(data0.data_x(), count);
		aggregate_avg(data0.data_y(), count);
	}

	template<size_t D, typename Y>
	static void aggregate_avg(symphas::FieldAxis<D, Y>& data0, size_t count)
	{
		for (iter_type i = 0; i < data0.length(); ++i)
		{
			aggregate_avg(data0.y[i], count);
		}
	}

	template<typename T1, typename T2>
	static void aggregate_avg(std::pair<T1, T2>& data0, size_t count)
	{
		aggregate_avg(data0.first, count);
		aggregate_avg(data0.second, count);
	}

	template<typename... Ts, size_t... Ns>
	static void aggregate_avg(std::tuple<Ts...>& data0, size_t count, std::index_sequence<Ns...>)
	{
		((aggregate_avg(std::get<Ns>(data0), count), ...));
	}

	template<typename... Ts>
	static void aggregate_avg(std::tuple<Ts...>& data0, size_t count)
	{
		aggregate_avg(data0, count, std::make_index_sequence<sizeof...(Ts)>{});
	}

}

/*
 * this class is used to conveniently maintain aggregate data
 * it does not own the pointer to the data, it relies on it to exist prior.
 * the raw_y values are then aggregated to the given pointer and
 * averaged...
 *
 * in other words: this is a REDUCE operation on the first input!
 */

 //template<typename R>
struct Aggr
{
	/* the given data is summed with the existing pointer based on
	 * the output of the given function, then the arithmetic average computed
	 * the arithmetic average is wrt the given number of (runs + 1), as
	 * this accounts for the provided value
	 */
	template<typename F, typename Y>
	static auto aggregate(F f, Y** raw_y, size_t runs)
	{
		auto data0 = f(raw_y[0]);
		for (iter_type i = 1; i < runs; ++i)
		{
			symphas::internal::aggregate_add_to(data0, f(raw_y[i]));
		}

		symphas::internal::aggregate_avg(data0, runs);
		return data0;
	}
};


// *****************************************************************************

//! Provides compile-time flag for using fourier_transform() function.
/*!
 *
 */
template<typename F>
struct ft_valid
{

	template<typename> struct bool_type { using type = std::true_type; };

	template<typename F0 = F, typename = decltype(symphas::lib::fourier_transform(std::declval<F0>()))>
	static constexpr auto test(int)
	{
		return true;
	}

	static constexpr bool test(...)
	{
		return false;
	}

	static const bool value = test(0);
};

template<typename X, typename Y, typename std::enable_if_t<ft_valid<symphas::Field<X, Y>>::value, int> = 0>
auto fourier_transform(symphas::Field<X, Y> const& field)
{
	return symphas::lib::fourier_transform(field);
}

template<size_t D, typename Y, typename std::enable_if_t<ft_valid<symphas::FieldAxis<D, Y*>>::value, int> = 0>
auto fourier_transform(symphas::FieldAxis<D, Y*> const& field)
{
	return symphas::lib::fourier_transform(field);
}

template<typename X, typename Y, typename std::enable_if_t<!ft_valid<symphas::Field<X, Y>>::value, int> = 0>
auto fourier_transform(symphas::Field<X, Y> const& field)
{
	return field;
}

template<size_t D, typename Y, typename std::enable_if_t<!ft_valid<symphas::FieldAxis<D, Y*>>::value, int> = 0>
auto fourier_transform(symphas::FieldAxis<D, Y*> const& field)
{
	return field;
}

// *****************************************************************************

/* type traits for analyzing the algorithm given to the process function
 * the goal is to determine the Y type from the given algorithm which generates a field
 * in particular, a problem was encountered when an algorithm is checked by types which
 * has a void return type (isn't one which is defined), and a compilation error would
 * be triggered because the other type traits aren't able to infer the Y axis type
 */

template<typename F>
struct alg_Y
{
public:
	using type = symphas::field_y_t<F>;
};


template<>
struct alg_Y<void>
{
	using type = void;
};


template<typename F>
struct alg_FT
{
protected:

	static auto apply(F f)
	{
		return fourier_transform(f);
	}

public:
	
	using type = typename alg_Y<std::invoke_result_t<decltype(&alg_FT<F>::apply), F>>::type;
};



// *****************************************************************************


/*
 * the process class generalizes the ability of the program to record data
 * it has the ability to either collect immediately (on every time step) or
 * collect immediately for many systems (thermal average)
 */
template<template<typename, typename, typename, typename> typename P,
	typename X, typename Y, typename S, typename Dt>
	struct Process
{

private:

	template<ProcessType type, typename M>
	using alg_return_type = typename Dt::template alg_return_type<type, X, Y, M>;
	template<ProcessType type, typename M>
	using alg_Y_type = typename alg_Y<alg_return_type<type, M>>::type;
	template<ProcessType type, typename M>
	using alg_FT_type = typename alg_FT<alg_return_type<type, M>>::type;
	

public:

	std::vector<DataParams> datas;
	Process(std::vector<DataParams> datas) : datas{ datas } 
	{}

	inline bool is_scheduled(iter_type index)
	{
		auto pred = [&](DataParams data) {
			return
				(data.type != ProcessType::NO)
				&& (data.save.current_save(index)); };
		
		return std::transform_reduce(
			datas.begin(), datas.end(), false, std::logical_or<bool>{}, pred);
	}

	/* specific collect methods for each data type that can be emitted
	 * a single point (x, y)
	 * regular line plot (linear x-axis, any y-values)
	 * a field (x-vector axis, any y-values)
	 */

	template<size_t I, typename R, size_t D, typename M>
	void collect_one(vector_data<R, D> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		auto sorted = data;
		static_cast<P<X, Y, S, Dt>&>(*this).template collect_one<I, R, D>(sorted.sort(), model, dir, index, data_index);
	}

	template<size_t I, typename R, typename M>
	void collect_one(point_data<R> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		static_cast<P<X, Y, S, Dt>&>(*this).template collect_one<I, R>(data, model, dir, index, data_index);
	}

	template<size_t I, typename R, typename M>
	void collect_one(scalar_data<R> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		static_cast<P<X, Y, S, Dt>&>(*this).template collect_one<I, R>(data, model, dir, index, data_index);
	}


	/*
	 * collection function executes the particular collect functions based on
	 * the process type and process tag
	 */
	template<size_t I, typename M>
	void collect(X data_x, Y data_y, len_type data_len, M const& model, const char* dir, iter_type index, iter_type data_index)
	{

		auto& data = datas[data_index];
		if constexpr (Dt::template has_alg<ProcessType::POINT, X, Y>())
		{
			if (data.type == ProcessType::POINT)
			{
				if (data.tag == ProcessTag::DFT)
				{
					using R = alg_Y_type<ProcessType::POINT, M>;
					if constexpr (std::is_pointer<R>::value)
					{
						collect_one<I>(
							fourier_transform(
								Dt::template get_data<ProcessType::POINT, X, Y>(data_x, data_y, data_len, model)), 
							model, dir, index, data_index);
					}
					else
					{
						fprintf(SYMPHAS_ERR, "attempting to perform an invalid transformation by taking the Fourier transform of"
							"a single point; this operation cannot be performed on the result of this algorithm\n");
						exit(1);
					}
				}
				else
				{
					collect_one<I>(
						Dt::template get_data<ProcessType::POINT, X, Y>(data_x, data_y, data_len, model), 
						model, dir, index, data_index);
				}
			}
		}
		if constexpr (Dt::template has_alg<ProcessType::SCALAR, X, Y>())
		{
			if (data.type == ProcessType::SCALAR)
			{
				if (data.tag == ProcessTag::DFT)
				{
					collect_one<I>(
						fourier_transform(
							Dt::template get_data<ProcessType::SCALAR, X, Y>(data_x, data_y, data_len, model)), 
						model, dir, index, data_index);
				}
				else
				{
					collect_one<I>(
						Dt::template get_data<ProcessType::SCALAR, X, Y>(data_x, data_y, data_len, model), 
						model, dir, index, data_index);
				}
			}
		}
		if constexpr (Dt::template has_alg<ProcessType::VECTOR, X, Y>())
		{
			if (data.type == ProcessType::VECTOR)
			{
				if (data.tag == ProcessTag::DFT)
				{
					collect_one<I>(
						fourier_transform(
							Dt::template get_data<ProcessType::VECTOR, X, Y>(data_x, data_y, data_len, model)),
						model, dir, index, data_index);
				}
				else
				{
					collect_one<I>(
						Dt::template get_data<ProcessType::VECTOR, X, Y>(data_x, data_y, data_len, model),
						model, dir, index, data_index);
				}
			}
		}
	}

	/*
	 * the aggregation is given an array of y values and the number of systems to aggregate over
	 * and then executes the given specialized aggregation function based on the parameters
	 */
	template<size_t I, typename M>
	void aggregate(X data_x, Y* data_ys, len_type data_len, M const& model, const char* dir, iter_type index, size_t runs, iter_type data_index)
	{
		auto& data = datas[data_index];
		if constexpr (Dt::template has_alg<ProcessType::POINT, X, Y>())
		{
			using R = alg_Y_type<ProcessType::POINT, M>;
			if (data.type == ProcessType::POINT)
			{
				if (data.tag == ProcessTag::DFT)
				{
					if constexpr (std::is_pointer<R>::value)
					{
						auto f = [&](auto data_y) {
							return fourier_transform(
								Dt::template get_data<ProcessType::POINT, X, Y>(data_x, data_y, data_len, model)); };
						collect_one<I>(Aggr/*<Rft>*/::aggregate(f, data_ys, runs), model, dir, index, data_index);
					}
					else
					{
						fprintf(SYMPHAS_ERR, "attempting to perform an invalid transformation by taking the Fourier transform of"
							"a single point; this operation cannot be performed on the result of this algorithm\n");
						exit(1);
					}
				}
				else
				{
					auto f = [&](auto data_y) {
						return Dt::template get_data<ProcessType::POINT, X, Y>(data_x, data_y, data_len, model); };
					collect_one<I>(Aggr::aggregate(f, data_ys, runs), model, dir, index, data_index);
				}
			}
		}
		if constexpr (Dt::template has_alg<ProcessType::SCALAR, X, Y>())
		{
			if (data.type == ProcessType::SCALAR)
			{
				if (data.tag == ProcessTag::DFT)
				{
					auto f = [&](auto data_y) {
						return fourier_transform(
							Dt::template get_data<ProcessType::SCALAR, X, Y>(data_x, data_y, data_len, model)); };
					collect_one<I>(Aggr::aggregate(f, data_ys, runs), model, dir, index, data_index);
				}
				else
				{
					auto f = [&] (auto data_y) {
						return Dt::template get_data<ProcessType::SCALAR, X, Y>(data_x, data_y, data_len, model); };
					collect_one<I>(Aggr::aggregate(f, data_ys, runs), model, dir, index, data_index);
				}
			}
		}
		if constexpr (Dt::template has_alg<ProcessType::VECTOR, X, Y>())
		{
			if (data.type == ProcessType::VECTOR)
			{
				if (data.tag == ProcessTag::DFT)
				{
					auto f = [&](auto data_y) {
						return fourier_transform(
							Dt::template get_data<ProcessType::VECTOR, X, Y>(data_x, data_y, data_len, model)); };
					collect_one<I>(Aggr::aggregate(f, data_ys, runs), model, dir, index, data_index);
				}
				else
				{
					auto f = [&](auto data_y) {
						return Dt::template get_data<ProcessType::VECTOR, X, Y>(data_x, data_y, data_len, model); };
					collect_one<I>(Aggr::aggregate(f, data_ys, runs), model, dir, index, data_index);
				}
			}
		}
	}


	template<size_t I, typename M>
	void collect(X data_x, Y data_y, len_type data_len, M const& model, const char* dir, iter_type index)
	{
		for (iter_type i = 0; i < datas.size(); ++i)
		{
			collect<I>(data_x, data_y, data_len, model, dir, index, i);
		}
	}
	template<size_t I, typename M>
	void aggregate(X data_x, Y* data_ys, len_type data_len, M const& model, const char* dir, iter_type index, size_t runs)
	{
		for (iter_type i = 0; i < datas.size(); ++i)
		{
			aggregate<I>(data_x, data_ys, data_len, model, dir, index, runs, i);
		}
	}
};



