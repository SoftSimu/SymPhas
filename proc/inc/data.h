#pragma once

#include <type_traits>

#include "datamacros.h"
#include "algorithm.h"


// ****************************************************************************************

// inheriting from these structs activates the corresponding operations
// takes template parameters for the ProcessType, the specialized algorithm which is based on the
// CRTP, and then the X and Y axis types


template<ProcessType type, template<typename, typename> typename SA, typename X, typename Y>
struct Algorithm
{
	static inline const ProcessType value = type;
	template<typename M>
	static auto get_data(X data_x, Y data_y, len_type len, M const& model)
	{
		return SA<X, Y>::get_data(data_x, data_y, len, model);
	}
};

template<ProcessType type, template <typename, typename> typename SA, typename Y>
struct Algorithm<type, SA, axis_1d_type*, Y*>
{
	static inline const ProcessType value = type;
	template<typename M>
	static auto get_data(const axis_1d_type* data_x, Y* data_y, len_type len, len_type l, M const& model)
	{
		return SA<axis_1d_type*, Y*>::get_data(data_x, data_y, len, l, model);
	}
};

template<ProcessType type, template <typename, typename> typename SA, typename Y>
struct Algorithm<type, SA, axis_2d_type*, Y*>
{
	static inline const ProcessType value = type;
	template<typename M>
	static auto get_data(const axis_2d_type* data_x, Y* data_y, len_type len, len_type l, len_type m, M const& model)
	{
		return SA<axis_2d_type*, Y*>::get_data(data_x, data_y, len, l, m, model);
	}
};

template<ProcessType type, template <typename, typename> typename SA, typename Y>
struct Algorithm<type, SA, axis_3d_type*, Y*>
{
	static inline const ProcessType value = type;
	template<typename M>
	static auto get_data(const axis_3d_type* data_x, Y* data_y, len_type len, len_type l, len_type m, iter_type n, M const& model)
	{
		return SA<axis_3d_type*, Y*>::get_data(data_x, data_y, len, l, m, n, model);
	}
};

template<template <typename, typename> typename SA, typename X, typename Y>
using AlgorithmVectorType = Algorithm<ProcessType::VECTOR, SA, X, Y>;
template<template <typename, typename> typename SA, typename X, typename Y>
using AlgorithmScalarType = Algorithm<ProcessType::SCALAR, SA, X, Y>;
template<template <typename, typename> typename SA, typename X, typename Y>
using AlgorithmPointType = Algorithm<ProcessType::POINT, SA, X, Y>;
template<template <typename, typename> typename SA, typename X, typename Y>
using AlgorithmDynamicType = Algorithm<ProcessType::DYNAMIC, SA, X, Y>;




// ****************************************************************************************

/*!
 * Data assumes that the x and y values given to it correspond by index,
 * and that the data is sorted. The object does NOT own the x and y arrays!
 * If given two arrays, it does not copy the data and also does not change it.
 */
template<template<typename, typename> typename... Algorithms>
struct Data
{
	template<ProcessType type, typename X, typename Y>
	static constexpr bool has_alg()
	{
		return ((... || (Algorithms<X, Y>::value == type)));
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(X data_x, Y data_y, len_type len, M const& model)
	{
		return get_data_apply<type, X, Y, M, Algorithms...>(data_x, data_y, len, model);
	}


	/* overloads based on being given a field as the argument to process data
	 */

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X, Y> &field, M const& model)
	{
		return get_data<type, X, Y>(field.data_x(), field.data_y(), field.length(), model);
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X, Y*>& field, M const& model)
	{
		return get_data<type, X, Y*>(field.data_x(), field.y, field.length(), model);
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X*, Y*>& field, M const& model)
	{
		return get_data<type, X*, Y*>(field.x, field.y, field.length(), model);
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X, Y>&& field, M const& model)
	{
		return get_data<type, X, Y>(field.data_x(), field.data_y(), field.length(), model);
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X, Y*>&& field, M const& model)
	{
		return get_data<type, X, Y*>(field.data_x(), field.y, field.length(), model);
	}

	template<ProcessType type, typename X, typename Y, typename M>
	static auto get_data(symphas::Field<X*, Y*>&& field, M const& model)
	{
		return get_data<type, X*, Y*>(field.x, field.y, field.length(), model);
	}

protected:

	template<ProcessType type, typename X, typename Y, typename M, 
		template<typename, typename> typename A, template<typename, typename> typename... As>
	static auto get_data_apply(X data_x, Y data_y, len_type len, M const& model)
	{
		if constexpr (A<X, Y>::value == type)
		{
			return call_algorithm<A>(data_x, data_y, len, model);
		}
		else if constexpr (sizeof...(As) > 0)
		{
			return get_data_apply<type, X, Y, M, As...>(data_x, data_y, len, model);
		}
		else
		{
			UNUSED(data_x, data_y, len);
		}
	}

	template<template<typename, typename> typename A, typename X, typename Y, typename M>
	static auto call_algorithm(X data_x, Y data_y, len_type len, M const& model)
	{
		return A<X, Y>::get_data(data_x, data_y, len, model);
	}

	template<template<typename, typename> typename A, typename Y, typename M>
	static auto call_algorithm(axis_1d_type* data_x, Y* data_y, len_type len, M const& model)
	{
		auto l = symphas::lib::get_dimensions<1>(data_x, len);
		return A<axis_1d_type*, Y*>::get_data(data_x, data_y, len, l, model);
	}

	template<template<typename, typename> typename A, typename Y, typename M>
	static auto call_algorithm(axis_2d_type* data_x, Y* data_y, len_type len, M const& model)
	{
		auto [l, m] = symphas::lib::get_dimensions<2>(data_x, len)._2();
		return A<axis_2d_type*, Y*>::get_data(data_x, data_y, len, l, m, model);
	}

	template<template<typename, typename> typename A, typename Y, typename M>
	static auto call_algorithm(axis_3d_type* data_x, Y* data_y, len_type len, M const& model)
	{
		auto [l, m, n] = symphas::lib::get_dimensions<3>(data_x, len)._3();
		return A<axis_3d_type*, Y*>::get_data(data_x, data_y, len, l, m, n, model);
	}


public:

	template<ProcessType type, typename X, typename Y, typename M>
	using alg_return_type = typename std::invoke_result_t<
		decltype(&Data<Algorithms...>::template get_data_apply<type, X, Y, M, Algorithms...>), 
		X, Y, len_type, M>;


};

// ****************************************************************************************


namespace symphas::lib
{

	template<typename... Ts>
	auto to_field_data(scalar_t* xs, std::tuple<Ts...> *ys, len_type len)
	{
		return scalar_data<std::tuple<Ts...>>{ xs, ys, len };
	}

	template<size_t N>
	auto to_field_data(scalar_t* xs, scalar_t (*&ys)[N], len_type len)
	{
		return scalar_data<scalar_t[N]>{ xs, ys, len };
	}

	inline auto to_field_data(scalar_t* xs, axis_nd_t<1>* ys, len_type len)
	{
		return scalar_data<axis_nd_t<1>>{ xs, ys, len };
	}

	inline auto to_field_data(scalar_t* xs, axis_nd_t<2>* ys, len_type len)
	{
		return scalar_data<axis_nd_t<2>>{ xs, ys, len };
	}

	inline auto to_field_data(scalar_t* xs, axis_nd_t<3>* ys, len_type len)
	{
		return scalar_data<axis_nd_t<3>>{ xs, ys, len };
	}

	inline auto to_field_data(scalar_t x, axis_nd_t<1> y)
	{
		return point_data<axis_nd_t<1>>{ x, y };
	}

	inline auto to_field_data(scalar_t x, axis_nd_t<2> y)
	{
		return point_data<axis_nd_t<2>>{ x, y };
	}

	inline auto to_field_data(scalar_t x, axis_nd_t<3> y)
	{
		return point_data<axis_nd_t<3>>{ x, y };
	}
}




