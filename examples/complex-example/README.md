# Complex SymPhas Example

This example demonstrates coupled phase field dynamics using SymPhas, showing how two fields can interact through coupling terms in their evolution equations.

## Features Demonstrated

- **Coupled Multi-Field System**: Two interacting phase fields with asymmetric coupling
- **Custom Initial Conditions**: Mathematical patterns defined using INITIAL_CONDITION_EQUATION
- **Conservative vs Non-Conservative Dynamics**: Field 1 uses Cahn-Hilliard (4th-order), Field 2 uses Allen-Cahn (2nd-order)
- **Periodic Boundary Conditions**: Seamless domain wrapping
- **Parameter Control**: Adjustable coupling strengths

## Model Description

The example implements two coupled phase fields:

**Field 1 (Conservative)**: Cahn-Hilliard dynamics with coupling
```
∂φ₁/∂t = ∇²(∇²φ₁ - φ₁³ + φ₁ + c₁φ₂)
```

**Field 2 (Non-Conservative)**: Allen-Cahn dynamics with coupling  
```
∂φ₂/∂t = ∇²φ₂ - φ₂³ + φ₂ + c₂φ₁
```

### Initial Conditions

- **Field 1**: Spiral pattern - tanh profile with angular modulation
- **Field 2**: Wave pattern - Gaussian envelope with sinusoidal modulation

### Parameters

- `c(1) = 1.0`: Coupling strength from field 2 to field 1
- `c(2) = 0.5`: Coupling strength from field 1 to field 2 (asymmetric coupling)

## Building and Running

Since CMakeLists.txt files are not needed, build from the SymPhas root directory:

1. **Build SymPhas with this example:**
   ```bash
   cd /path/to/symphas
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release \
        -DMAIN_FILE=../examples/complex-example/main.cpp \
        -DSOLVER_INCLUDE_HEADER_DIR=../examples/solvers \
        -DSOLVER_INCLUDE_HEADER_NAME=solverinclude.h
   make
   ```

   **Note**: The `SOLVER_INCLUDE_HEADER_DIR` and `SOLVER_INCLUDE_HEADER_NAME` parameters are required because this example uses built-in solver objects. The `SOLVER_INCLUDE_HEADER_NAME` automatically includes the necessary solver definitions when using solvers like `SolverFT`.

2. **Run the simulation:**
   ```bash
   ./symphas
   ```

## Understanding the Output

The simulation runs for 2000 time steps with dt=0.01, printing:

- Grid and domain information
- Coupling parameter values  
- Final field ranges and mean values
- Evolution statistics

## Key Learning Points

1. **Model Definition**: Use `MODEL(name, (field_types), EVOLUTION(...))` syntax
2. **Multi-Field Coupling**: Fields can reference each other using `op(1)`, `op(2)`, etc.
3. **Initial Conditions**: Define patterns using `INITIAL_CONDITION_EQUATION`
4. **Parameter Usage**: Access model parameters with `c(1)`, `c(2)`, etc.
5. **Conservative vs Non-Conservative**: Different equation orders for different physics

## Experimentation Ideas

Try modifying:
- **Coupling strengths**: Change the `param_values` to see different interactions
- **Grid resolution**: Adjust `128_h` to change spatial resolution
- **Domain size**: Modify `4.0_dh` to change physical domain size
- **Time stepping**: Change dt or number of steps
- **Initial conditions**: Create new patterns in the INITIAL_CONDITION_EQUATION definitions

## Physical Interpretation

This model could represent:
- **Phase separation with defects**: Field 1 as composition, Field 2 as structural defects
- **Reaction-diffusion systems**: Chemical species with different diffusivities
- **Multi-physics coupling**: Different physical phenomena influencing each other

The asymmetric coupling (c₁ ≠ c₂) creates interesting non-equilibrium dynamics where the fields affect each other differently.
