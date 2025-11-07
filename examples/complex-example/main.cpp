/*!
 * \file main.cpp
 * \brief Complex SymPhas example demonstrating coupled phase fields.
 * 
 * This example shows:
 * - Multiple coupled phase fields with cross-coupling
 * - Custom initial conditions
 * - Periodic boundary conditions
 * - Multi-field dynamics
 */

#include "solverft.h"
#include "symphas.h"

// Define a coupled two-field model
// Field 1: Cahn-Hilliard type dynamics with coupling
// Field 2: Allen-Cahn type dynamics with coupling
MODEL(CoupledFields, (SCALARS(2)),
    EVOLUTION(
        // Field 1: Conservative dynamics with coupling to field 2
        dop(1) = lap(lap(op(1)) - power(op(1), 3) + op(1) + c(1) * op(2)),
        // Field 2: Non-conservative dynamics with coupling to field 1  
        dop(2) = lap(op(2)) - power(op(2), 3) + op(2) + c(2) * op(1)
    )
)

// Custom initial condition functions
INITIAL_CONDITION_EQUATION(SPIRAL_PATTERN, (1, 1), 
    0.1 * tanh(5 * (0.3 - sqrt(x*x + y*y))) * cos(4 * atan2(y, x))
)

INITIAL_CONDITION_EQUATION(WAVE_PATTERN, (1, 1),
    0.05 * exp(-(x*x + y*y) / 0.1) * sin(2 * PI * x) * cos(2 * PI * y) 
)

int main() {
    using namespace symphas;
    
    std::cout << "=== Complex Coupled Fields Example ===" << std::endl;
    
    // Set up 2D grid: 32x32 points, domain size 1x1, periodic boundaries
    auto interval = BoundaryType::PERIODIC || 32_h / 1.0_dh;
    auto grid = interval * interval;
    
    // Set up initial conditions for both fields
    auto params = grid << (
        Inside::VORONOI <<= {6},
        Inside::CIRCLE <<= {1}
    );
    
    // Set coupling parameters
    // c(1): coupling strength from field 2 to field 1
    // c(2): coupling strength from field 1 to field 2
    double param_values[]{1.0, 0.5};  // Asymmetric coupling
    
    // Create the coupled model with basic solver (no FFT required)
    model_CoupledFields_t<2, SolverFT<Stencil2d2h<>>> model{param_values, 2,
                                                            params};
    
    std::cout << "Running coupled field simulation..." << std::endl;
    std::cout << "Grid: 32x32, Domain: 1x1, Time step: 0.001" << std::endl;
    std::cout << "Coupling parameters: c1=" << param_values[0]
              << ", c2=" << param_values[1] << std::endl;
    
    // Run simulation: 100 steps with dt=0.001
    find_solution(model, 0.001, 100);
    
    // Print final statistics for both fields
    auto field1 = model.get_field<0>();
    auto field2 = model.get_field<1>();
    
    auto [min1, max1] = std::minmax_element(field1.values, field1.values + field1.len);
    auto [min2, max2] = std::minmax_element(field2.values, field2.values + field2.len);
    
    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "Field 1 range: [" << *min1 << ", " << *max1 << "]" << std::endl;
    std::cout << "Field 2 range: [" << *min2 << ", " << *max2 << "]" << std::endl;
    
    // Calculate and display coupling effects
    double field1_mean = 0.0, field2_mean = 0.0;
    for (iter_type i = 0; i < field1.len; ++i) {
        field1_mean += field1.values[i];
        field2_mean += field2.values[i];
    }
    field1_mean /= field1.len;
    field2_mean /= field2.len;
    
    std::cout << "Field 1 mean: " << field1_mean << std::endl;
    std::cout << "Field 2 mean: " << field2_mean << std::endl;
    std::cout << "\n✓ Coupled field dynamics successfully demonstrated!" << std::endl;
    
    return 0;
}
