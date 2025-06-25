# Simple SymPhas Tutorial

This # Configure (include solver headers for built-in solvers)
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DMAIN_FILE="examples/tutorial/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run
./symphas_impl
```

**Note**: The `SOLVER_INCLUDE_HEADER_DIR` and `SOLVER_INCLUDE_HEADER_NAME` parameters are required because this example uses built-in solver objects (`SolverFT<Stencil2d2h<>>`). The `SOLVER_INCLUDE_HEADER_NAME` automatically includes the necessary solver definitions.ains a minimal example of using SymPhas to simulate the Allen-Cahn equation.

## What this example demonstrates:

- **Minimal Setup**: Everything is defined in a single `main.cpp` file
- **Model Definition**: Allen-Cahn equation using the `MODEL` macro
- **Grid Setup**: 64×64 2D grid with periodic boundaries
- **Problem Spec Notation**: Setting up grid dimensions, spacing, and initial conditions
- **Solver Usage**: Using the explicit finite difference solver (`SolverFT`)
- **Basic Output**: Printing field statistics after simulation

## The Allen-Cahn Model

The Allen-Cahn equation is one of the simplest phase-field models:

```
∂φ/∂t = ∇²φ - (φ² - 1)φ
```

Where:
- `φ` is the order parameter (phase field)
- The `∇²φ` term represents diffusion
- The `-(φ² - 1)φ` term represents a double-well potential that drives phase separation

## Building and Running

From the SymPhas root directory:

```bash
# Create build directory
mkdir build && cd build

# Configure (include solver headers for SolverFT)
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DMAIN_FILE="examples/tutorial/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run
./symphas_impl
```

**Note**: The `SOLVER_INCLUDE_HEADER_DIR` parameters are required because this example uses `SolverFT` (Forward Time integrator).

## Expected Output

You should see output similar to:
```
-------------------------------------------------------------
Phase field problem of 1 system:
        system<0> ... 64 x 64
No coefficients provided.
Default coefficient value:  1.00.
-------------------------------------------------------------
[  0.0%]          0        +     0.000
Progress      Index   Runtime (seconds)
[100.0%]       1000        +     X.XXX
Simulation completed. Field range: -X.XXX to X.XXX
```

## Next Steps

After successfully running this tutorial:

1. **Try the complex example**: Check out `examples/complex-example/` for a more advanced tutorial with I/O, configuration files, and more features
2. **Experiment with parameters**: Modify the grid size, time step, or number of iterations
3. **Add output**: Include `io::save_grid(data)` to save field data to files
4. **Try different models**: Check the "Creating Your Own Models" section in the documentation
