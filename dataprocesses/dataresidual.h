#pragma once

#include <execution>

#include "processdynamic.h"
#include "modelspecialized.h"

//! Finite difference solver.
/*!
 * The boundaries are updated after the provisional variables are computed
 * from the corresponding equations. The boundary data is typically the same
 * as the first field that is given, but in practise it does not matter
 * for the numerical results unless finite difference approximations are
 * applied extensively to the provisional variables. This is not recommended.
 */
NEW_SOLVER_WITH_STENCIL(SolverRS)


public:

	SolverRS(len_type* dims, double h, double dt) : stencil_t(dims, h), parent_type(dt) {}


	template<typename S>
	void step(S& sys, double dt) const {}

	template<typename S, typename E>
	inline void equation(std::pair<S, E>& r) const
	{
		auto& [sys, equation] = r;
		expr::prune::update(equation);
		expr::result_interior(equation, sys.get().dframe);
	}

	template<size_t En, typename SS, typename S, typename E>
	auto form_expr_one(SS&&, std::pair<S, E> const& e) const
	{
		auto [sys, equation] = e;
		auto eq_ft = expr::apply_operators(equation);

		return std::make_pair(sys, eq_ft);
	}

	static auto make_solver(symphas::problem_parameters_type const& parameters)
	{
		double h = parameters.get_interval_data()[0].at(Axis::X).width();
		size_t dim = parameters.get_dimension();
		double dt = parameters.get_time_step();

		/* the dimensions of the problem are taken from the first system
		 * since the solver assumes they are homogeneous
		 */
		len_type* dims = new len_type[dim];
		for (iter_type i = 0; i < dim; ++i)
		{
			Axis side = symphas::index_to_axis(i);
			dims[i] = parameters.get_interval_data()[0].at(side).get_count() + 2 * THICKNESS;
		}

		auto s = this_type{ dims, h, dt };
		delete[] dims;
		return s;
	}


	// given the grid/equation pair, evaluate the equation into the grid
	// element by element
	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>& r) const
	{
		auto& [grid, equation] = r;
		if constexpr (expr::has_state<E>::value)
		{
			expr::prune::update(equation);
		}

		expr::result_interior(equation, expr::BaseData<G>::get(grid));
	}

	template<typename G, typename E>
	void evaluate_one(std::pair<G, E>&& r) const
	{
		auto& [grid, equation] = r;
		if constexpr (expr::has_state<E>::value)
		{
			expr::prune::update(equation);
		}

		expr::result_interior(equation, expr::BaseData<G>::get(grid));
	}
};


ASSOCIATE_SOLVER_SYSTEM_TYPE(SolverRS, SolverSystemFD)
ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverRS, ProvisionalSystemFD)
SYMPHAS_SOLVER_ALL_SUPPORTED(SolverRS)




DEFINE_ALGORITHM_POINT(Residual)
DEFINE_ALGORITHM_DYNAMIC(Residual)
DEFINE_DATA(Residual, POINT, ALG_POINT, ALG_DYNAMIC)

//! Diffusion only.
MODEL(RES_LAP, (SCALAR),
	EVOLUTION(
		dop(1) = lap(op(1)))
)

//! lap of cubic term.
MODEL(RES_LAP3, (SCALAR),
	EVOLUTION(
		dop(1) = lap(op(1) * op(1) * op(1)))
)

//! Diffusion only.
MODEL(RES_BILAP, (SCALAR),
	EVOLUTION(
		dop(1) = bilap(op(1)))
)

//! Diffusion only.
MODEL(RES_QULAP, (SCALAR),
	EVOLUTION(
		dop(1) = Diff(6) * op(1))
)

ALGORITHM_POINT_1D_DECLARATION(Residual)
{
	len_type length = model.template system<0>().length();
	std::tuple<std::vector<Y>, std::vector<Y>> ys{ length, length };
	if constexpr (model_num_parameters<model_type>::value == 1)
	{
		auto m = model_swap_solver<SolverRS<typename SelfSelectingStencil<1, 2>::template Points<3, 4, 5>>>{}(model);
		m.equation();

		m.template system<0>().persist(std::get<0>(ys).data());
		grid::copy_interior(m.template system<0>().dframe, std::get<1>(ys).data());

		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}
	else
	{
		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}
}

ALGORITHM_POINT_2D_DECLARATION(Residual)
{
	len_type length = model.template system<0>().length();
	std::tuple<std::vector<Y>, std::vector<Y>> ys{ length, length };
	if constexpr (model_num_parameters<model_type>::value == 1)
	{
		//auto m = model_swap_solver<SolverRS<typename SelfSelectingStencil<D, 2>::template Points<9, 6, 13>>>{}(model);

		auto parameters = model.generate_parameters();
		auto* coeffs = model.get_coeff();
		size_t num_coeffs = model.get_num_coeff();

		//auto m = model_swap_solver<SolverRS<typename SelfSelectingStencil<2, 2>::template Points<9, 6, 13>>>{}(model);

		model_RES_LAP_t<2, SolverRS<typename SelfSelectingStencil<2, 2>::template Points<9, 6, 13>>> m(coeffs, num_coeffs, parameters);
		//model_RES_LAP3_t<2, SolverRS<typename SelfSelectingStencil<2, 2>::template Points<9, 6, 13>>> m1(coeffs, num_coeffs, parameters);
		//model_RES_BILAP_t<2, SolverRS<typename SelfSelectingStencil<2, 2>::template Points<9, 6, 13>>> m2(coeffs, num_coeffs, parameters);
		m.equation();
		//m1.equation();
		//m2.equation();

		len_type length = m.template system<0>().length();

		model.template system<0>().persist(std::get<0>(ys).data());
		grid::copy_interior(m.template system<0>().dframe, std::get<1>(ys).data());
		//grid::copy_interior(m1.template system<0>().dframe, std::get<2>(ys).data());
		//grid::copy_interior(m2.template system<0>().dframe, std::get<3>(ys).data());

		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}
	else
	{
		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}
}

ALGORITHM_POINT_3D_DECLARATION(Residual)
{
	len_type length = model.template system<0>().length();
	std::tuple<std::vector<Y>, std::vector<Y>> ys{ length, length };
	if constexpr (model_num_parameters<model_type>::value == 1)
	{
		constexpr size_t D = model_dimension<model_type>::value;
		auto m = model_swap_solver<SolverRS<typename SelfSelectingStencil<D, 2>::template Points<7, 10, 21>>>{}(model);
		m.equation();

		len_type length = m.template system<0>().length();

		m.template system<0>().persist(std::get<0>(ys).data());
		grid::copy_interior(m.template system<0>().dframe, std::get<1>(ys).data());

		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}
	else
	{
		return point_data<std::tuple<std::vector<Y>, std::vector<Y>>>{ model.get_time(), ys, length };
	}

}


ALGORITHM_DYNAMIC_DECLARATION(Residual)
{
	using Y0 = model_field_t<model_type, 0>;

	size_t num_residuals = len - 1;
	std::vector<vector_data<std::tuple<Y0, Y0>, D>> result;
	
	for (iter_type i = 0; i < num_residuals; ++i)
	{

		iter_type pos = static_cast<iter_type>(std::ceil((static_cast<double>(len) / num_residuals) * (i + 1)) - 1);
		len_type len_data = data_y[pos]->length();

		std::tuple<Y0, Y0>* ys = new std::tuple<Y0, Y0>[len_data];
		auto* xs = symphas::lib::new_system_axis_list<D>(model.template system<0>().get_info().intervals);

		double dt = data_y[pos]->data_x() - data_y[pos - 1]->data_x();

		for (iter_type i = 0; i < len_data; ++i)
		{
			auto psi_t = (std::get<0>(data_y[pos]->data_y())[i] - std::get<0>(data_y[pos - 1]->data_y())[i]) / dt;
			auto psi_x = std::get<1>(data_y[pos - 1]->data_y())[i];
			//auto psi_x3 = std::get<2>(data_y[pos - 1]->data_y())[i];
			//auto psi_2x = std::get<3>(data_y[pos - 1]->data_y())[i];
			std::get<0>(ys[i]) = psi_t;
			std::get<1>(ys[i]) = psi_x;
			//std::get<2>(ys[i]) = psi_x3;
			//std::get<3>(ys[i]) = psi_2x;
		}
		result.emplace_back(xs, ys, len_data);
	}

	return result;
}


