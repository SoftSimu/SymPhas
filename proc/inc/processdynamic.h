#pragma once

#include "data.h"
#include "savedefines.h"
#include "../modules-models.h"


template<size_t I, typename F, typename P, typename M>
struct PersistentData;

template<typename X, typename Y, typename S, typename Dt>
struct ProcessDynamic;

/*
 * PersistentData class maintains data, it is meant to be used in
 * a static context for data collection
 */
template<size_t I, typename F, typename X, typename Y, typename S, typename Dt, typename M>
struct PersistentData<I, F, ProcessDynamic<X, Y, S, Dt>, M>
{
	using P = ProcessDynamic<X, Y, S, Dt>;

	PersistentData() : keep{}, process{ *((P*)0) }, model{ *((M*)0) }, dir{ nullptr }, stop{ 0 } {}

	PersistentData(P const& process, M const& model, const char* dir, iter_type stop) :
		keep{}, process{ process }, model{ model }, 
		dir{ (std::strlen(dir) > 0) ? new char[std::strlen(dir + 1)] : nullptr }, stop { stop } 
	{
		std::strcpy(this->dir, dir);
	}

	PersistentData(PersistentData<I, F, P, M> const& other) : PersistentData(other.process, other.model, other.dir, other.stop) {}
	PersistentData(PersistentData<I, F, P, M>&& other) : PersistentData() 
	{
		swap(*this, other);
	}

	PersistentData<I, F, P, M>& operator=(PersistentData<I, F, P, M> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(PersistentData<I, F, P, M>& first, PersistentData<I, F, P, M>& second)
	{
		using std::swap;
		swap(first.keep, second.keep);
		swap(first.process, second.process);
		swap(first.model, second.model);
		swap(first.dir, second.dir);
		swap(first.stop, second.stop);
	}

	void append(iter_type index, F data)
	{
		if (index <= stop)
		{
			keep.emplace_back(index, data);
		}
	}

	dynamic_data<F> get_data()
	{
		std::shared_ptr<iter_type[]> x{ new iter_type[keep.size()] };
		std::shared_ptr<F* []> y{ new F * [keep.size()] };

		for (iter_type i = 0; i < keep.size(); ++i)
		{
			x[i] = keep[i].first;
			y[i] = &(keep[i].second);
		}
		return { x, y, static_cast<len_type>(keep.size()) };
	}


	~PersistentData()
	{
		process.template output<I>(Dt::template get_data<ProcessType::DYNAMIC>(get_data(), model), dir, stop);
	}

	std::vector<std::pair<iter_type, F>> keep;

	P const& process;
	M const& model;
	char* dir;
	iter_type stop;

};


template<size_t I, typename F, typename P, typename M>
auto get_persistent_data(P const& process, M const& model, const char* dir, iter_type stop)
{
	return PersistentData<I, F, P, M>(process, model, dir, stop);
}

namespace symphas::internal
{
	template<typename T1, typename T2>
	void assign_data(T1& to, T2 const& from)
	{
		to = from;
	}

	template<typename T1, typename T2, size_t N>
	void assign_data(T1(&to)[N], T2 const (&from)[N])
	{
		std::copy(from, from + N, to);
	}

	template<typename... Ts, size_t... Is>
	void assign_data(std::tuple<Ts*...>(&to), std::tuple<Ts...> const* from, len_type length, std::index_sequence<Is...>)
	{
		for (iter_type i = 0; i < length; ++i)
		{
			(assign_data(std::get<Is>(to)[i], std::get<Is>(from[i])), ...);
		}
	}

	template<typename... Ts>
	void assign_data(std::tuple<Ts*...> (&to), std::tuple<Ts...> const* from, len_type length)
	{
		assign_data(to, from, length, std::make_index_sequence<sizeof...(Ts)>{});
	}
}

/*
 * the dynamic process is used for computing dynamic information about the system.
 * this is accomplished by storing a static variable on each call to the collection
 * functions, and then building it until it needs to be emitted at the end
 */
template<typename X, typename Y, typename S, typename Dt>
struct ProcessDynamic : Process<ProcessDynamic, X, Y, S, Dt>
{
	using parent_type = Process<::ProcessDynamic, X, Y, S, Dt>;

	// inherit the constructor
	using parent_type::parent_type;




	template<size_t I, typename R>
	static void output(point_data<R> const& data, const char* dir, iter_type index)
	{
		symphas::io::save_pts<S>(symphas::lib::combine_data(data.data_x(), data.data_y(), data.length()), dir, index, I);
	}

	template<size_t I, typename R>
	static void output(point_data<R*> const& data, const char* dir, iter_type index)
	{
		symphas::io::save_pts<S>(symphas::lib::combine_data(data.data_x(), data.y, data.length()), dir, index, I);
	}

	template<size_t I, typename R>
	static void output(scalar_data<R> const& data, const char* dir, iter_type index)
	{
		symphas::io::save_abs<S>(symphas::lib::combine_data(data.x, data.y, data.length()), dir, index, I);
	}

	template<size_t I, typename T1, typename T2>
	static void output(scalar_data<std::pair<T1, T2>> const& data, const char* dir, iter_type index)
	{
		auto* data_first = new T1[data.length()];
		auto* data_second = new T2[data.length()];
		for (iter_type i = 0; i < data.length(); ++i)
		{
			symphas::internal::assign_data(data_first[i], data.data_y()[i].first);
			symphas::internal::assign_data(data_second[i], data.data_y()[i].second);
		}
		output<I>(scalar_data<T1>{ data.data_x(), data_first, data.length() }, dir, index);
		output<I>(scalar_data<T2>{ data.data_x(), data_second, data.length() }, dir, index);
	}
	/*template<size_t I, typename... Ts, size_t... Ns>
	static void fill_and_save(scalar_data<std::tuple<Ts...>> const& data, const char* dir, iter_type index, std::tuple<Ts...>& tuple, std::index_sequence<Ns...>)
	{
		for (iter_type i = 0; i < data.length(); ++i)
		{
			((symphas::internal::assign_data(std::get<Ns>(tuple)[i], std::get<Ns>(data.data_y()[i])), ...));
		}
		((output<I>(scalar_data<Ts>{ data.data_x(), std::get<Ns>(tuple), data.length() }, dir, index), ...));
	}

	template<size_t I, typename... Ts>
	static void output(scalar_data<std::tuple<Ts...>> const& data, const char* dir, iter_type index)
	{
		auto tuple = std::make_tuple(new Ts[data.length()]...);
		fill_and_save<I>(data, tuple, dir, index, std::make_index_sequence<sizeof...(Ts)>{});
	}*/

	template<size_t I, typename R, size_t D>
	static void output(vector_data<R, D> const& data, const char* dir, iter_type index)
	{
		symphas::io::save_vec<S>(symphas::lib::combine_data(data.x, data.y, data.length()), dir, index, I);
	}

	template<size_t I, typename T1, typename T2, size_t D>
	static void output(vector_data<std::pair<T1, T2>, D> const& data, const char* dir, iter_type index)
	{
		auto* data_first = new T1[data.length()];
		auto* data_second = new T2[data.length()];
		for (iter_type i = 0; i < data.length(); ++i)
		{
			symphas::internal::assign_data(data_first[i], data.data_y()[i].first);
			symphas::internal::assign_data(data_second[i], data.data_y()[i].second);
		}
		output<I>(vector_data<T1, D>{ data.data_x(), data_first, data.length() }, dir, index);
		output<I>(vector_data<T2, D>{ data.data_x(), data_second, data.length() }, dir, index);
	}

	template<size_t I, typename... Ts, size_t D, size_t... Is>
	static void output(vector_data<std::tuple<Ts...>, D> const& data, const char* dir, iter_type index, std::index_sequence<Is...>)
	{
		std::tuple<Ts*...> data_expand;
		((std::get<Is>(data_expand) = new Ts[data.length()]), ...);
		symphas::internal::assign_data(data_expand, data.data_y().get(), data.length());
		
		(output<I>(vector_data<Ts, D>{ data.data_x(), std::get<Is>(data_expand), data.length() }, dir, index), ...);

	}

	template<size_t I, typename T1, typename... Ts, size_t D>
	static void output(vector_data<std::tuple<T1, Ts...>, D> const& data, const char* dir, iter_type index)
	{
		output<I>(data, dir, index, std::make_index_sequence<sizeof...(Ts) + 1>{});
	}

	template<size_t I, typename R>
	static void output(std::vector<R> const& datas, const char* dir, iter_type index)
	{
		iter_type i = 0;
		for (const auto& data : datas)
		{
			output<I>(data, dir, i++);
		}
	}


	// ************************************************************************************


	template<size_t I, typename R, typename M>
	void collect_one(point_data<R> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		static auto keep = get_persistent_data<I, point_data<R>>(*this, model, dir, parent_type::datas[data_index].save.get_stop());
		keep.append(index, data);

		//if (index == parent_type::datas[data_index].save.get_stop())
		//{
		//	output<I>(Dt::template get_data<ProcessType::DYNAMIC>(keep.get_data(), model), dir, index);
		//}

	}

	template<size_t I, typename R, typename M>
	void collect_one(scalar_data<R> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		static auto keep = get_persistent_data<scalar_data<R>>(*this, model, dir, parent_type::datas[data_index].save.get_stop());
		keep.append(index, data);
		//static PersistentData<scalar_data<R>> keep;
		//keep.append(index, data);
		//
		//if (index == parent_type::datas[data_index].save.get_stop())
		//{
		//	output<I>(Dt::template get_data<ProcessType::DYNAMIC>(keep.get_data(), model), dir, index);
		//}
	}

	template<size_t I, typename R, size_t D, typename M>
	void collect_one(vector_data<R, D> const& data, M const& model, const char* dir, iter_type index, iter_type data_index)
	{
		static auto keep = get_persistent_data<vector_data<R, D>>(*this, model, dir, parent_type::datas[data_index].save.get_stop());
		keep.append(index, data);
		//static PersistentData<vector_data<R, D>> keep;
		//keep.append(index, data);
		//
		//if (index == parent_type::datas[data_index].save.get_stop())
		//{
		//	output<I>(Dt::template get_data<ProcessType::DYNAMIC>(keep.get_data(), model), dir, index);
		//}
	}





};



