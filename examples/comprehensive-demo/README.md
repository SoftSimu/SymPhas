# Comprehensive SymPhas Demo

This directory contains a comprehensive demonstration of SymPhas capabilities, showcasing multiple models, initial condition methods, and advanced features in a single example.

## What this example demonstrates:

- **Multiple Model Approaches**: Evolution equations, free energy functionals, and custom dynamics
- **Various Initial Conditions**: Uniform, lambda functions, expressions, and geometric shapes
- **Problem Spec Notation**: Different ways to set up grids and boundary conditions
- **Multi-Field Systems**: Both single and multi-field models
- **Mixed Field Types**: Scalar and vector fields together
- **Comprehensive I/O**: Multiple output files with descriptive names
- **Modular Headers**: Using separate solver include files

## Models Included

### 1. Allen-Cahn (Evolution Approach)
```cpp
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)
```

### 2. Two-Field System
```cpp
MODEL(TwoField, (SCALARS(2)),
    EVOLUTION(
        dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1),
        dop(2) = lap(op(2)) - (power(op(2), 2) - 1_n) * op(2)
    )
)
```

### 3. Free Energy Formulation
```cpp
MODEL(AllenCahn_FE, (SCALAR), 
    FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1))))
)
```

### 4. Mixed Field Types with Custom Dynamics
```cpp
MODEL(MixedFields, (SCALAR, VECTOR),
    FREE_ENERGY(
        (EQUATION_OF(1)(-lap(-DF(1)) - grad(op(1)) * DF(2)),
         EQUATION_OF(2)(lap(DF(2)) + grad(op(1)) * -DF(1))),
        INT(LANDAU_FE(op(1)) + 2_n * power(op(2), 2))
    )
)
```

## Initial Condition Methods Demonstrated

1. **Uniform Random**: `Inside::UNIFORM <<= {-1, 1}`
2. **Lambda Functions**: Custom initialization functions
3. **Expression-Based**: Using `INITIAL_CONDITION_EQUATION` macro
4. **Geometric**: Circle-based initialization

## Building and Running

From the SymPhas root directory:

```bash
# Create build directory
mkdir build && cd build

# Configure with separate solver header
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DMAIN_FILE="examples/comprehensive-demo/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="examples/comprehensive-demo" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run
./symphas_impl
```

## Expected Output

```
=== SymPhas Comprehensive Demo ===
This example demonstrates various SymPhas features:
- Multiple model definition approaches
- Different initial condition methods
- Problem Spec Notation variations
- Field analysis and I/O

Running Allen-Cahn model with uniform IC...
  Result: Field range [-X.XXX, X.XXX]
Running free energy model with expression IC...
  Result: Field range [-X.XXX, X.XXX]
Running two-field system...
  Field 1 range: [-X.XXX, X.XXX]
  Field 2 range: [-X.XXX, X.XXX]

=== Demo Complete ===
Check the 'data/' directory for output files.
```

## Output Files Generated

- `data/allencahn_0.txt`: Basic Allen-Cahn simulation
- `data/freeenergy_0.txt`: Free energy formulation result
- `data/twofield_1_0.txt`: First field of two-field system
- `data/twofield_2_1.txt`: Second field of two-field system

## Learning Progression

1. **Simple Tutorial** (`examples/tutorial/`): Basic concepts
2. **Advanced Features** (`examples/advanced-features/`): Free energy and I/O
3. **This Example**: Comprehensive feature showcase
4. **Custom Development**: Build your own models and applications

## Key Concepts Illustrated

- **Model Definition Flexibility**: Multiple approaches for the same physics
- **Problem Spec Notation**: Intuitive domain setup language
- **Extensibility**: Easy to add new models and features  
- **Professional Output**: Organized data management and user feedback
- **Modular Design**: Separation of concerns for maintainability
