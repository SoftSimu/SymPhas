#include "symphas.h"
#include "solverinclude.h"

// model definition: note that power is NOT denoted using ^
MODEL(A, (SCALAR), EVOLUTION(dop(1) = lap(op(1)) - (POWER(op(1), 2) - 1_n) * op(1)))
MODEL(C, (SCALARS(2)), EVOLUTION(dop(1) = lap(op(1)) - (POWER(op(1), 2) - 1_n) * op(1), dop(2) = lap(op(2)) - (POWER(op(2), 2) - 1_n) * op(2)))
//MODEL(AF, (SCALARS(1)), FREE_ENERGY(dop(1) = lap(op(1)) - (POWER(op(1), 2) - 1_n) * op(1), dop(2) = lap(op(2)) - (POWER(op(2), 2) - 1_n) * op(2)))

	    MODEL(AFE, (SCALAR), 
			FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1)))))
INITIAL_CONDITION_EQUATION(INIT, (2, 3), sin(x) + cos(y))
MODEL(MH_FE, (SCALAR, VECTOR),
	FREE_ENERGY((
			EQUATION_OF(1)(-lap(-DF(1)) - grad(op(1)) * DF(2)), 
			EQUATION_OF(2)(lap(DF(2)) + grad(op(1)) * -DF(1))
			),
		INT(LANDAU_FE(op(1)) + _2 * POWER(op(2), 2)))
)
int main()
{
        
    auto init = [](auto index, const auto* dims, auto dimension) 
    {
        int pos[2]{};
        grid::get_grid_position(pos, dims, index);
        return std::sin(pos[0]) * std::cos(pos[1]);    
    };

    symphas::interval_element_type interval(1000);
    symphas::b_data_type bdata(2, BoundaryType::PERIODIC);
	symphas::init_data_type tdata(Inside::UNIFORM, { -1, 1 });

    symphas::problem_parameters_type parameters(1);

    using namespace symphas;

    auto tdata1 = Inside::UNIFORM <<= { -1, 1 };
    auto tdat02 = Inside::CIRCLE / InsideTag::RANDOM <<= { -1, 1 };
    auto tdata2 = Inside::EXPRESSION << "INIT";
    auto tdat3a = Inside::CIRCLE / InsideTag::VARA <<= {1, 1};
    auto tdat4a = Inside::LAMBDA << init;
    symphas::init_data_functor lambd(init);
    symphas::init_data_type tdata00(lambd);

    //another way to make an interval:
    auto interval0 = BoundaryType::PERIODIC || 100_h / 0.1_dh || BoundaryType::PERIODIC;
    auto grid_intervals = interval0 * interval0;

    // set the parameters using 'problem spec notation'
    parameters[0] << interval*interval << bdata << tdata;

    

    // the model typename is derived from the model definition, the given name is inserted into 'model_*_t'. 
    model_C_t<2, SolverFT<Stencil2d2h<>>> model{ parameters };
    model_AFE_t<2, SolverSP<>> modela{ parameters };
    symphas::find_solution(modela, 0.1, 50);
	auto field = model.get_field<0>();
	symphas::io::save_grid(field);
}