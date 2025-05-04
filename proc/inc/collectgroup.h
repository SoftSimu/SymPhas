#pragma once

#include "datalib.h"
#include "model.h"



// **************************************************************************************************

namespace symphas::internal
{

	template<size_t N, typename Y, typename M>
	Y** make_data_ys(M const* models, size_t runs)
	{
		Y** data_ys = new Y * [runs];
		for (iter_type i = 0; i < runs; ++i)
		{
			auto& b = models[i].template system<N>();
			data_ys[i] = new Y[b.length()];
		}
		return data_ys;

	}
	template<size_t N, typename Y, typename M>
	void populate_data_ys(Y** data_ys, M const* models, size_t runs)
	{
		for (iter_type i = 0; i < runs; ++i)
		{
			auto& b = models[i].template system<N>();
			b.persist(data_ys[i]);
		}

	}

}



/*
 * facilitates the actual process step by running the specialized
 * processing objects (Collector) over each field corresponding to type T
 */
template<typename Y, template<typename, typename> typename... Collectors>
struct CollectorGroup
{

protected:

	template<size_t I, template<typename, typename> typename Collector, size_t D, typename Sp, typename... S>
	static bool collect_combine(Model<D, Sp, S...> const& model, const char* dir)
	{
		using M = Model<D, Sp, S...>;
		constexpr size_t N = M::template index_of_type<Y, I>;

		static Collector<axis_nd_t<D>*, Y*> c;

		if (c.is_scheduled(model.get_index()))
		{
			static auto* data_x = symphas::lib::new_system_axis_list<D>(model.template system<N>().get_info().intervals);
			model.template do_for_field<N>(&collect_one<N, D, Collector, M>, model, dir, model.get_index(), c, data_x);
		}

		if constexpr (I + 1 < M::template num_fields<Y>())
		{
			return collect_combine<I + 1, Collector>(model, dir) || c.is_scheduled(model.get_index());
		}
		else
		{
			return c.is_scheduled(model.get_index());
		}
	}

	template<size_t I, template<typename, typename> typename Collector, size_t D, typename Sp, typename S, typename... Ts>
	static bool collect_combine(ArrayModel<D, Sp, S, Ts...> const& model, const char* dir)
	{
		using M = ArrayModel<D, Sp, S, Ts...>;
		constexpr size_t N = M::template index_of_type<Y, I>;

		static Collector<axis_nd_t<D>*, Y*> c;

		if (c.is_scheduled(model.get_index()))
		{
			static auto* data_x = symphas::lib::new_system_axis_list<D>(model.template system<N>().get_info().intervals);
			model.template do_for_field<N>(&collect_one<N, D, Collector, M>, model, dir, model.get_index(), c, data_x);
		}

		if constexpr (I + 1 < sizeof...(Ts) + 1)
		{
			return collect_combine<I + 1, Collector>(model, dir) || c.is_scheduled(model.get_index());
		}
		else
		{
			return c.is_scheduled(model.get_index());
		}
	}

	/* performs a collection procedure for the given set of values
	 */
	template<size_t N, size_t D, template<typename, typename> typename Collector, typename M, typename X = axis_nd_t<D>>
	static void collect_one(Y* data_y, len_type data_len, M const& model, const char* dir, 
		iter_type index, Collector<X*, Y*> c, X* data_x)
	{
		c.template collect<N>(data_x, data_y, data_len, model, dir, index);
	}


	/* in order to perform the aggregation between models, construct an array of field
	 * values to pass to the collector
	 */
	template<size_t I, template<typename, typename> typename Collector, size_t D, typename Sp, typename... S>
	static bool aggregate_combine(Model<D, Sp, S...> const* models, const char* dir, size_t runs)
	{
		using M = Model<D, Sp, S...>;
		constexpr size_t N = M::template index_of_type<Y, I>;

		M const& model0 = *models;
		static Collector<axis_nd_t<D>*, Y*> c;

		if (c.is_scheduled(model0.get_index()))
		{
			static Y** data_ys = symphas::internal::make_data_ys<N, Y>(models, runs);
			symphas::internal::populate_data_ys<N, Y>(data_ys, models, runs);

			static auto* data_x = symphas::lib::new_system_axis_list<D>(model0.template system<N>().get_info().intervals);
			aggregate_one<N, D, Collector, M>(data_ys, model0.template system<N>().length(), model0, dir, model0.get_index(), runs, c, data_x);
		}

		if constexpr (I + 1 < M::template num_fields<Y>())
		{
			return aggregate_combine<I + 1, Collector>(models, dir, runs) || c.is_scheduled(model0.get_index());
		}
		else
		{
			return c.is_scheduled(model0.get_index());
		}
	}

	/* in order to perform the aggregation between models, construct an array of field
	 * values to pass to the collector
	 */
	template<size_t I, template<typename, typename> typename Collector, size_t D, typename Sp, typename S, typename... Ts>
	static bool aggregate_combine(ArrayModel<D, Sp, S, Ts...> const* models, const char* dir, size_t runs)
	{
		using M = ArrayModel<D, Sp, S, Ts...>;
		constexpr size_t N = M::template index_of_type<Y, I>;

		M const& model0 = *models;
		static Collector<axis_nd_t<D>*, Y*> c;

		if (c.is_scheduled(model0.get_index()))
		{
			static Y** data_ys = symphas::internal::make_data_ys<N, Y>(models, runs);
			symphas::internal::populate_data_ys<N, Y>(data_ys, models, runs);

			static auto* data_x = symphas::lib::new_system_axis_list<D>(model0.template system<N>().get_info().intervals);
			aggregate_one<N, D, Collector, M>(data_ys, model0.template system<N>().length(), model0, dir, model0.get_index(), runs, c, data_x);
		}

		if constexpr (I + 1 < sizeof...(Ts) + 1)
		{
			return aggregate_combine<I + 1, Collector>(models, dir, runs) || c.is_scheduled(model0.get_index());
		}
		else
		{
			return c.is_scheduled(model0.get_index());
		}
	}

	/* apply the aggregation with the given axes and array of yvalues
	 */
	template<size_t N, size_t D, template<typename, typename> typename Collector, typename M, typename X = axis_nd_t<D>>
	static void aggregate_one(Y** data_ys, len_type data_len, M const& model, const char* dir, 
		iter_type index, size_t runs, Collector<X*, Y*> c, X* data_x)
	{
		c.template aggregate<N>(data_x, data_ys, data_len, model, dir, index, runs);
	}

public:


	template<typename M>
	static bool collect_start([[maybe_unused]] M const& model, [[maybe_unused]] const char* dir)
	{
		if constexpr (M::template index_of_type<Y> >= 0)
		{
			int value = ((... + static_cast<iter_type>(collect_combine<0, Collectors>(model, dir))));
			return value > 0;
		}
		else
		{
			return false;
		}
	}

	template<typename M>
	static bool aggregate_start([[maybe_unused]] M const* models, [[maybe_unused]] const char* dir, [[maybe_unused]] size_t runs)
	{
		if constexpr (M::template index_of_type<Y> >= 0)
		{
			int value = ((... + static_cast<iter_type>(aggregate_combine<0, Collectors>(models, dir, runs))));
			return value > 0;
		}
		else
		{
			return false;
		}
	}
};


