# SymPhas MPI Example

This example demonstrates how to use SymPhas with MPI (Message Passing Interface) for parallel computing across multiple processes.

## Features Demonstrated

- **MPI Integration**: Running SymPhas simulations across multiple processes
- **Cell Migration Model**: Demonstrates the built-in cell migration model without motility
- **Parallel Execution**: Distributing computation across multiple CPU cores/nodes
- **SolverFT with MPI**: Using the Forward Time integrator in a distributed environment

## Requirements

- **MPI Library**: OpenMPI, MPICH, or Intel MPI
- **SymPhas with MPI Support**: Built with `-DUSE_MPI=ON`

## Building and Running

From the SymPhas root directory:

```bash
# Create build directory
mkdir build && cd build

# Configure with MPI support
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DUSE_MPI=ON \
    -DMAIN_FILE="examples/using-mpi/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run with MPI (example: 4 processes)
mpirun -np 4 ./symphas_impl
```

**Note**: The `SOLVER_INCLUDE_HEADER_DIR` and `SOLVER_INCLUDE_HEADER_NAME` parameters are required because this example uses built-in solver objects. The `SOLVER_INCLUDE_HEADER_NAME` automatically includes the necessary solver definitions when using solvers like `SolverFT`.

## Model Description

This example uses the predefined `CELL_MIGRATION_NO_MOTILITY` model, which simulates:
- **Cell density fields** representing cellular populations
- **Diffusion and chemotaxis** without directed motility
- **Interaction dynamics** between different cell types

## MPI Parallelization

SymPhas automatically handles:
- **Domain decomposition** across MPI processes
- **Boundary communication** between neighboring processes
- **Load balancing** for optimal performance
- **Collective operations** for global reductions

## Expected Output

With MPI, you'll see output from multiple processes, showing:
- Process rank identification
- Domain partitioning information
- Simulation progress from each process
- Synchronized timing information

## Performance Tips

- Use powers of 2 for the number of processes for optimal domain decomposition
- Ensure the grid size is divisible by the number of processes
- Monitor communication overhead vs computation time
- Use appropriate MPI binding strategies for your hardware

## Troubleshooting

- **MPI not found**: Install MPI development libraries and rebuild SymPhas
- **Communication errors**: Check firewall settings and MPI configuration
- **Load imbalance**: Ensure grid dimensions are compatible with process count
- **Memory issues**: Reduce grid size or increase available memory per process
