#include "symphas.h"

// Various model definitions demonstrating different approaches

// Basic evolution equation approach
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

// Multi-field model
MODEL(TwoField, (SCALARS(2)),
    EVOLUTION(
        dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1),
        dop(2) = lap(op(2)) - (power(op(2), 2) - 1_n) * op(2)
    )
)

// Free energy approach
MODEL(AllenCahn_FE, (SCALAR), 
    FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1))))
)

// Mixed field types with custom dynamics
MODEL(MixedFields, (SCALAR, VECTOR),
    FREE_ENERGY(
        (EQUATION_OF(1)(-lap(-DF(1)) - grad(op(1)) * DF(2)),
         EQUATION_OF(2)(lap(DF(2)) + grad(op(1)) * -DF(1))),
        INT(LANDAU_FE(op(1)) + 2_n * power(op(2), 2))
    )
)

// Custom initial condition function 
INITIAL_CONDITION_EQUATION(SINCOS_INIT, (2, 3), sin(x) + cos(y))

int main() {
    using namespace symphas;
    
    std::cout << "=== SymPhas Comprehensive Demo ===" << std::endl;
    std::cout << "This example demonstrates various SymPhas features:" << std::endl;
    std::cout << "- Multiple model definition approaches" << std::endl;
    std::cout << "- Different initial condition methods" << std::endl;
    std::cout << "- Problem Spec Notation variations" << std::endl;
    std::cout << "- Field analysis and I/O" << std::endl << std::endl;
    
    // Lambda-based initial condition
    auto lambda_init = [](auto index, const auto* dims, auto dimension) {
        int pos[2]{};
        grid::get_grid_position(pos, dims, index);
        return std::sin(pos[0] * 0.1) * std::cos(pos[1] * 0.1);
    };
    
    // Different ways to set up initial conditions
    auto uniform_ic = Inside::UNIFORM <<= {-1, 1};
    auto circle_ic = Inside::CIRCLE / InsideTag::RANDOM <<= {-1, 1};
    auto expression_ic = Inside::EXPRESSION << "SINCOS_INIT";
    auto lambda_ic = Inside::LAMBDA << lambda_init;
    
    // Different ways to set up grids using Problem Spec Notation
    auto interval1 = BoundaryType::PERIODIC || 50_h / 0.5_dh;
    auto grid1 = interval1 * interval1;
    
    auto interval2 = 50_h / 0.5_dh * 50_h / 0.5_dh;
    auto grid2 = interval2 << (Side::LEFT / BoundaryType::PERIODIC)
                          << (Side::RIGHT / BoundaryType::PERIODIC);
    
    // Set up parameters for different models
    auto params1 = grid1 << uniform_ic;
    auto params2 = grid1 << expression_ic;
    
    // Example 1: Basic Allen-Cahn with uniform initial conditions
    std::cout << "Running Allen-Cahn model with uniform IC..." << std::endl;
    model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model1{params1};
    find_solution(model1, 0.01, 1000);
    
    auto field1 = model1.get_field<0>();
    auto min1 = *std::min_element(field1.values, field1.values + field1.len);
    auto max1 = *std::max_element(field1.values, field1.values + field1.len);
    std::cout << "  Result: Field range [" << min1 << ", " << max1 << "]" << std::endl;
    io::save_grid(field1, "allencahn");
    
    // Example 2: Free energy formulation with expression-based IC
    std::cout << "Running free energy model with expression IC..." << std::endl;
    model_AllenCahn_FE_t<2, SolverFT<Stencil2d2h<>>> model2{params2};
    find_solution(model2, 0.01, 1000);
    
    auto field2 = model2.get_field<0>();
    auto min2 = *std::min_element(field2.values, field2.values + field2.len);
    auto max2 = *std::max_element(field2.values, field2.values + field2.len);
    std::cout << "  Result: Field range [" << min2 << ", " << max2 << "]" << std::endl;
    io::save_grid(field2, "freeenergy");
    
    // Example 3: Multi-field system
    std::cout << "Running two-field system..." << std::endl;
    auto multifield_params = grid1 << uniform_ic;
    model_TwoField_t<2, SolverFT<Stencil2d2h<>>> model3{multifield_params};
    find_solution(model3, 0.01, 1000);
    
    auto field3a = model3.get_field<0>();
    auto field3b = model3.get_field<1>();
    std::cout << "  Field 1 range: [" 
              << *std::min_element(field3a.values, field3a.values + field3a.len)
              << ", " 
              << *std::max_element(field3a.values, field3a.values + field3a.len)
              << "]" << std::endl;
    std::cout << "  Field 2 range: [" 
              << *std::min_element(field3b.values, field3b.values + field3b.len)
              << ", " 
              << *std::max_element(field3b.values, field3b.values + field3b.len)
              << "]" << std::endl;
    io::save_grid(field3a, "twofield_1");
    io::save_grid(field3b, "twofield_2");
    
    std::cout << std::endl << "=== Demo Complete ===" << std::endl;
    std::cout << "Check the 'data/' directory for output files." << std::endl;
    
    return 0;
}
