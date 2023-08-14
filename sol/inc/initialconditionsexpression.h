
#pragma once

#include "initialconditionslib.h"
#include "expressionlib.h"
#include "expressions.h"


#define INITIAL_CONDITION_EQUATION(NAME, DIMENSIONS, EQUATION) \
DEFINE_INIT_EXPR_EQUATION(NAME, DIMENSIONS, EQUATION) \
NEXT_INIT_EXPR_INDEX(NAME) \
INIT_EXPR_WRAPPER_FUNC(INIT_EXPR_INDEX_NAME(NAME), DIMENSIONS, \
if (std::strcmp(name, #NAME) == 0) { return symphas::internal::run_init_expr<D, InitExpression_ ## NAME<D>>(ax, values, interval, vdata, coeff, num_coeff); })



//! \cond


#define DEFINE_INIT_EXPR_EQUATION(NAME, DIMENSIONS, EQUATION)				\
template<size_t D>															\
struct InitExpression_ ## NAME : InitExpression<D>							\
{																			\
	using parent_type = InitExpression<D>;									\
	using parent_type::parent_type;											\
	auto get_equation() { using namespace expr; using namespace std; return EQUATION; }	\
};


//! Used in assigning unique names to models for indexing.
/*!
 * The naming format of a model index is defined. Each index name has
 * to be different and requires a prefix.
 */
#define INIT_EXPR_INDEX_NAME(PREFIX_NAME) __INIT_EXPR_ ## PREFIX_NAME ## _index

 //! Iterates to the next model index for compile-time constant model indexing.
 /*!
  * Convenience definition for setting the next index and providing
  * the PREFIX argument which names it.
  * Importantly, there cannot be more than #MAX_DEFINED_MODELS models defined
  * because after that, the counter will no longer increment.
  */
#define NEXT_INIT_EXPR_INDEX(PREFIX_NAME) \
constexpr int INIT_EXPR_INDEX_NAME(PREFIX_NAME) = decltype(symphas::internal::init_expr_counter(symphas::internal::init_expr_count_index<255>{}))::value; \
namespace symphas::internal { \
constexpr init_expr_count_index<INIT_EXPR_INDEX_NAME(PREFIX_NAME)> init_expr_counter(init_expr_count_index<INIT_EXPR_INDEX_NAME(PREFIX_NAME)>); }


//! Implements the association between a model and a name.
#define INIT_EXPR_WRAPPER_FUNC(N, DIMENSIONS, ...) \
template<size_t D> \
struct symphas::internal::init_expr_call_wrapper<D, N> \
{ \
	template<typename T> \
	static int call(const char* name, Axis ax, T* values, grid::region_interval<D> const& interval, \
		symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff) \
	{ \
		if constexpr (symphas::lib::is_value_in_seq<size_t, D, std::index_sequence<SINGLE_ARG DIMENSIONS>>::value) \
		{ __VA_ARGS__ } \
		return init_expr_call_wrapper<D, N - 1>::call(name, ax, values, interval, vdata, coeff, num_coeff); \
	} \
};


//! Generates initial conditions using a defined equation. 
/*!
 * Initial conditions can be generated based on a predefined symbolic
 * algebra expression. This expression is
 * defined in a way similar to the model definitions. This class
 * defines that expression, and then this will be created similarly
 * to the model selection routine. The model will select
 * the expression-based initial condition based on the configuration.
 *
 * The expression can refer to the parameters of the model, or
 * parameters of the initial condition routine itself, which would
 * be passed by the configuration.
 *
 * \tparam F The functor type which is used to generate the initial
 * conditions.
 */

template<size_t D>
struct InitExpression
{
	struct AxisData
	{
		AxisData() : dims{}, info{ nullptr, 0 } {}
		template<size_t... Is>
		AxisData(const len_type(&dims)[D], symphas::grid_info const& info, std::index_sequence<Is...>) : dims{ dims[Is]... }, info{info} {}
		AxisData(const len_type(&dims)[D], symphas::grid_info const& info) : AxisData(dims, info, std::make_index_sequence<D>{}) {}

		auto get_info() const
		{
			return info;
		}

		len_type dims[D];
		symphas::grid_info info;
	};

	InitExpression(grid::region_interval<D> const& interval, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
		: coeff{ (num_coeff > 0) ? new double[num_coeff] : nullptr }, num_coeff{ num_coeff }, data{ interval.dims, vdata }
	{
		std::copy(coeff, coeff + num_coeff, this->coeff);
	}

	InitExpression(InitExpression const& other) : InitExpression(other.coeff, other.num_coeff) {}
	InitExpression(InitExpression&& other) : InitExpression() { swap(*this, other); }
	InitExpression& operator=(InitExpression other) { swap(*this, other); return *this; }

	friend void swap(InitExpression& first, InitExpression& second)
	{
		using std::swap;
		swap(first.coeff, second.coeff);
		swap(first.num_coeff, second.num_coeff);
	}

	template<size_t>
	auto system()
	{
		return data;
	}

	auto param(size_t I)
	{
		if (I < num_coeff)
		{
			return expr::make_literal(coeff[I]);
		}
		else
		{
			return expr::make_literal(DEFAULT_COEFF_VALUE);
		}
	}

	double* coeff;
	size_t num_coeff;
	AxisData data;

protected:

	InitExpression() : coeff{ nullptr }, num_coeff{ 0 }, data{} {}

};

template<size_t... Ds>
struct InitExpressionBuild;

template<size_t D0>
struct InitExpressionBuild<D0> : InitExpression<D0>
{
	using parent_type = InitExpression<D0>;
	using parent_type::parent_type;
};

template<size_t D0, size_t D1, size_t... Ds>
struct InitExpressionBuild<D0, D1, Ds...> : InitExpression<D0>, InitExpressionBuild<D1, Ds...>
{
	using parent_type = InitExpression<D0>;
	using parent_type::parent_type;
};


namespace symphas::internal
{

	//! Recursive inheritance defining an incrementing value.
	/*!
	 * Primary object using recursive inheritance in order to
	 * define a member variable and return one of incremented value.
	 *
	 * A higher order index is given on the next usage by overloading
	 * the function for which the terminating point is instantiated. below
	 * it is never defined, because only the return type matters.
	 */
	template<int N>
	struct init_expr_count_index : init_expr_count_index<N - 1>
	{
		static const int value = N + 1;
	};

	//! Base specialization to terminate the recursive inheritance.
	/*!
	 * Base specialization which terminates the recursive inheritance for
	 * incrementing model indices.
	 */
	template<>
	struct init_expr_count_index<0>
	{
		static const int value = 1;
	};

	/* overloading this function on a new parameter and updating the return
	 * type will provide a new index on the next usage by decltype
	 */
	constexpr init_expr_count_index<0> init_expr_counter(init_expr_count_index<0>);


	template<size_t D, int N>
	struct init_expr_call_wrapper
	{
		template<typename T>
		static int call(
			const char* name, Axis ax, T* values, grid::region_interval<D> const& interval,
			symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
		{
			return -1;
		}
	};

	template<size_t D, typename EI, typename T>
	bool run_init_expr(
		Axis ax, T* values, grid::region_interval<D> const& interval,
		symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
	{
		EI init_expr{ interval, vdata, coeff, num_coeff };
		auto expr = init_expr.get_equation();

		auto begin = symphas::internal::ic_iterator<T, D, decltype(expr)>(ax, values, expr, interval);
		auto end = symphas::internal::ic_iterator<T, D, decltype(expr)>(ax, values, expr, interval, grid::length<D>(interval));

		auto it = symphas::data_iterator_region(values, interval);

		std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			begin, end, it);

		return true;
	}


}

template<size_t D, typename T>
bool match_init_expr(const char* initname, Axis ax, T* values,
	grid::region_interval<D> const& interval, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff);


//! \endcond


