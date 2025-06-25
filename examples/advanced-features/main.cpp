#include "symphas.h"
#include "solverft.h"

// Model A using free energy formulation
MODEL(ModelA_FreeEnergy, (SCALAR), 
    FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1))))
)

int main() {
    using namespace symphas;
    
    // Set up a larger 2D grid: 200x200 with spacing 0.5
    auto interval = BoundaryType::PERIODIC || 100_h / 0.5_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::UNIFORM <<= {-1, 1});
    
    // Create model using free energy formulation
    model_ModelA_FreeEnergy_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Run simulation for 5000 steps with dt=0.01
    std::cout << "Starting Allen-Cahn simulation using free energy formulation..." << std::endl;
    find_solution(model, 0.01, 5000);
    
    // Get and analyze the final field
    auto field = model.get_field<0>();
    
    // Calculate some statistics
    auto min_val = *std::min_element(field.values, field.values + field.len);
    auto max_val = *std::max_element(field.values, field.values + field.len);
    
    // Calculate average
    double sum = 0.0;
    for (len_type i = 0; i < field.len; ++i) {
        sum += field.values[i];
    }
    double average = sum / field.len;
    
    std::cout << "\n=== Simulation Results ===" << std::endl;
    std::cout << "Grid size: " << field.dims[0] << "x" << field.dims[1] << std::endl;
    std::cout << "Field range: [" << min_val << ", " << max_val << "]" << std::endl;
    std::cout << "Average value: " << average << std::endl;
    
    // Save the field data to file
    io::save_grid(field);
    std::cout << "Field data saved to 'data/data_0.txt'" << std::endl;
    
    return 0;
}