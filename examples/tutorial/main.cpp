#include "symphas.h"

// Define the Allen-Cahn model
MODEL(AllenCahn, (SCALAR), 
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

int main() {
    using namespace symphas;
    
    // Set up a 2D grid: 64x64 with spacing 0.5 and periodic boundaries
    auto interval = BoundaryType::PERIODIC || 32_h / 0.5_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::UNIFORM <<= {-0.1, 0.1});
    
    // Create model with Forward Euler solver
    model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Run simulation for 1000 steps with dt=0.01
    find_solution(model, 0.01, 1000);
    
    // Print final field statistics
    auto data = model.get_field<0>();
    std::cout << "Simulation completed. Field range: " 
              << *std::min_element(data.values, data.values + data.len) 
              << " to " 
              << *std::max_element(data.values, data.values + data.len) 
              << std::endl;
    
    return 0;
}
