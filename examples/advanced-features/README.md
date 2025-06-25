# Advanced Features Example

This directory demonstrates advanced SymPhas features with focus on **free energy formulation** and **I/O capabilities**.

## What this example demonstrates:

- **Free Energy Models**: Defining models using free energy functionals instead of direct evolution equations
- **I/O Integration**: Saving simulation data to files using the `io` module
- **Field Analysis**: Computing and displaying field statistics
- **Larger Scale**: Simulation on a 200×200 grid
- **Clean Output**: Professional result presentation

## The Model

This example uses the Allen-Cahn model defined through its free energy:

```cpp
MODEL(ModelA_FreeEnergy, (SCALAR), 
    FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1))))
)
```

This is equivalent to the evolution equation:
```
∂φ/∂t = ∇²φ - (φ² - 1)φ
```

But demonstrates how SymPhas can automatically derive evolution equations from free energy functionals.

## Building and Running

From the SymPhas root directory:

```bash
# Create build directory
mkdir build && cd build

# Configure (include solver headers for built-in solvers)
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DMAIN_FILE="examples/advanced-features/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run
./symphas_impl
```

**Note**: The `SOLVER_INCLUDE_HEADER_DIR` and `SOLVER_INCLUDE_HEADER_NAME` parameters are required because this example uses built-in solver objects. The `SOLVER_INCLUDE_HEADER_NAME` automatically includes the necessary solver definitions (e.g., for `SolverFT`, `SolverSP`, etc.).

## Expected Output

```
Starting Allen-Cahn simulation using free energy formulation...
[simulation progress...]

=== Simulation Results ===
Grid size: 200x200
Field range: [-0.XXXX, 0.XXXX]
Average value: X.XXXX
Field data saved to 'data/data_0.txt'
```

## Key Features Demonstrated

1. **Free Energy Approach**: Using `FREE_ENERGY` macro with `LANDAU_FE` predefined functional
2. **Larger Grid**: 200×200 grid showing scalability
3. **Field Analysis**: Computing min, max, and average values
4. **Professional I/O**: Using `io::save_grid()` for data output
5. **User Feedback**: Clear progress and result reporting

## Next Steps

- **Comprehensive Demo**: Check out `examples/comprehensive-demo/` for multiple models and advanced features
- **Custom Models**: Experiment with your own free energy functionals
- **Visualization**: Use the saved data with plotting tools like gnuplot or Python
- **Parameter Studies**: Modify grid size, time step, or simulation length
