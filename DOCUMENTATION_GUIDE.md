# SymPhas Complete Guide & User Manual

## Table of Contents

1. [Overview and Introduction](#overview-and-introduction)
   - [What is SymPhas?](#what-is-symphas)
   - [Key Features](#key-features)
   - [Applications](#applications)

2. [Installation Guide](#installation-guide)
   - [Overview](#overview)
   - [Prerequisites and Dependencies](#prerequisites-and-dependencies)
   - [Platform-Specific Prerequisites](#platform-specific-prerequisites)
   - [Installation Process](#installation-process)
   - [Configuration Options](#configuration-options)
   - [Troubleshooting](#troubleshooting)
   - [Installation Verification](#installation-verification)

3. [Getting Started Tutorial](#getting-started-tutorial)
   - [Prerequisites](#prerequisites-1)
   - [CMake Configuration for Projects](#cmake-configuration-for-projects)
   - [Tutorial: Your First SymPhas Simulation](#tutorial-your-first-symphas-simulation)
   - [More Advanced Examples](#more-advanced-examples)
   - [Building Your Own Project](#building-your-own-project)
   - [Next Steps](#next-steps-1)

4. [Build Configuration Options](#build-configuration-options)
   - [Required Configuration Parameters](#required-configuration-parameters)
   - [Core Module Configuration](#core-module-configuration)
   - [Advanced Solver and Dimension Configuration](#advanced-solver-and-dimension-configuration)
   - [Stencil Override System](#stencil-override-system)
   - [Parallel Computing Configuration](#parallel-computing-configuration)
   - [Output and Debugging Configuration](#output-and-debugging-configuration)
   - [Complete Configuration Examples](#complete-configuration-examples)
   - [Troubleshooting Configuration Issues](#troubleshooting-configuration-issues)

4. [Architecture and Design](#architecture-and-design)
   - [System Architecture](#system-architecture)
   - [Design Principles](#design-principles)

5. [Module Documentation](#module-documentation)
   - [Core Modules (Required)](#core-modules-required)
     - [lib Module](#1-lib-module)
     - [datatypes Module](#2-datatypes-module)
     - [sym Module](#3-sym-module)
     - [sol Module](#4-sol-module)
   - [Optional Modules](#optional-modules)
     - [io Module (Optional)](#5-io-module-optional)
     - [conf Module (Optional)](#6-conf-module-optional)

6. [Input/Output (IO) System](#inputoutput-io-system)
   - [Overview](#overview-2)
   - [Output Format Selection](#output-format-selection)
   - [Available Output Formats](#available-output-formats)
   - [Setting Output Format](#setting-output-format)
   - [SaveParams Configuration](#saveparams-configuration)
   - [Checkpoint System](#checkpoint-system)
   - [Program Spec Notation for IO](#program-spec-notation-for-io)
   - [File Organization](#file-organization)
   - [Practical Examples](#practical-examples-1)

7. [Core Components](#core-components)
   - [Grid System](#grid-system)
   - [Symbolic Algebra System](#symbolic-algebra-system)
   - [Model Definition System](#model-definition-system)
   - [Solver System](#solver-system)

8. [Symbols and Functions Reference](#symbols-and-functions-reference)
   - [Coefficients and Parameters](#coefficients-and-parameters)
   - [Order Parameters and DerivaWtives](#order-parameters-and-derivatives)
   - [Supported Field Types](#supported-field-types)
   - [Constants and Literals](#constants-and-literals)
   - [Time Symbol](#time-symbol)
   - [Mathematical Operators and Functions](#mathematical-operators-and-functions)
   - [Differential Operators](#differential-operators)
   - [Vector and Tensor Operations](#vector-and-tensor-operations)
   - [Summation and Indexing](#summation-and-indexing)
   - [Noise and Randomness](#noise-and-randomness)
   - [Special Macros and Utilities](#special-macros-and-utilities)

9. [Creating Your Own Models](#creating-your-own-models)
   - [Model Definition Approaches](#model-definition-approaches)
   - [1. Direct Evolution Equations](#1-direct-evolution-equations)
   - [2. Free Energy Formulation](#2-free-energy-formulation)
   - [3. Evolution with Preamble](#3-evolution-with-preamble)
   - [4. Provisional Variables](#4-provisional-variables)
     - [What are Provisional Variables?](#what-are-provisional-variables)
     - [Defining Provisional Variables](#defining-provisional-variables)
     - [Complete Working Example](#complete-working-example)
     - [Multi-Variable Provisional Examples](#multi-variable-provisional-examples)
     - [Provisional Variable Persistence and Analysis](#provisional-variable-persistence-and-analysis)
     - [Saving Provisional Variables](#saving-provisional-variables)
     - [Practical Applications](#practical-applications)
     - [Provisional vs. Preamble Variables](#provisional-vs-preamble-variables)
   - [5. Model Linking and Naming](#5-model-linking-and-naming)
     - [How LINK_WITH_NAME Works](#how-link_with_name-works)
     - [Basic Model Linking](#basic-model-linking)
     - [Multiple Aliases for the Same Model](#multiple-aliases-for-the-same-model)
     - [Using Linked Names in Configuration Files](#using-linked-names-in-configuration-files)
     - [Runtime Model Selection](#runtime-model-selection)
     - [Benefits of Model Linking](#benefits-of-model-linking)
     - [Best Practices for Model Linking](#best-practices-for-model-linking)

10. [Boundary Conditions](#boundary-conditions)
   - [Overview of Boundary Types](#overview-of-boundary-types)
   - [DEFAULT Boundary Type with Tags](#default-boundary-type-with-tags)
   - [Common Use Cases and Examples](#common-use-cases-and-examples)
   - [Advanced Features](#advanced-features)

11. [JSON Configuration System](#json-configuration-system)
   - [JSON Configuration Overview](#json-configuration-overview)
   - [Complete JSON Structure](#complete-json-structure)
   - [Detailed Section Documentation](#detailed-section-documentation)
   - [Enhanced Save Parameter Types](#enhanced-save-parameter-types)
   - [Multi-Phase Time Stepping](#multi-phase-time-stepping)
   - [Working with Configuration Files](#working-with-configuration-files)

12. [Problem Spec Notation](#problem-spec-notation)
   - [Overview](#overview-1)
   - [Grid Dimensions and Discretization](#grid-dimensions-and-discretization)
   - [Boundary Conditions](#boundary-conditions-1)
   - [Initial Conditions](#initial-conditions)
   - [Complete Parameter Assembly](#complete-parameter-assembly)
   - [Multi-Field Systems](#multi-field-systems)
   - [Practical Examples](#practical-examples)
   - [Best Practices](#best-practices)
   - [Common Pitfalls](#common-pitfalls)

13. [Defining Your Own Solver](#defining-your-own-solver)
   - [Solver Architecture Overview](#solver-architecture-overview)
   - [Basic Solver Structure](#basic-solver-structure)
   - [Required Functions](#required-functions)
   - [Solver Implementation Examples](#solver-implementation-examples)
   - [Advanced Features](#advanced-features-1)
   - [Integration with SymPhas](#integration-with-symphas)
   - [Best Practices and Guidelines](#best-practices-and-guidelines)
   - [Common Pitfalls](#common-pitfalls-1)

14. [Additional Resources](#additional-resources)
   - [Provisional Variable Saving](#provisional-variable-saving)
   - [API Reference](#api-reference)
   - [Developer's Guide](#developers-guide)
   - [Example Projects](#example-projects)

---

## Overview and Introduction

### What is SymPhas?

**SymPhas** (Symbolic Phase Field Solver) is a comprehensive C++ framework for implementing solvers for phase-field problems with compile-time symbolic algebra. It provides a **Domain Specific Language (DSL)** implemented through C++ macros that allows you to write mathematical equations in natural syntax, which is then compiled into high-performance C++ code.

The framework supports multiple field types (**scalar**, **complex**, and **vector** fields) and provides a high-performance, modular approach to solving multi-phase-field differential equations across various scientific and engineering domains.

### Key Features

- **Domain Specific Language (DSL)**: Write mathematical equations in natural syntax using C++ macros
- **Multiple Field Types**: Support for scalar, complex, and vector fields
- **Compile-time Symbolic Algebra**: Expression trees are fully formulated at compile time, eliminating runtime branching
- **High Performance**: Optimized runtimes with linear scaling for any model type
- **Parallelization Support**: Built-in MPI and CUDA support for distributed and GPU computing
- **Modular Architecture**: Separation of concerns with well-defined module interfaces
- **Extensive Documentation**: Comprehensive API documentation and examples
- **Third-party Development**: Framework designed for easy extension and customization

### Applications

SymPhas supports any multi-phase-field problem that can be formulated field-theoretically:
- Dendrite growth simulations
- Reaction-diffusion systems
- Biological systems modeling
- Phase-field crystal problems
- Turing pattern formation
- Gray-Scott reaction systems

---

## Installation Guide

### Quick Start

For users who want to get SymPhas running immediately:

**1. Install Dependencies:**
```bash
# Ubuntu/Debian
sudo apt install build-essential cmake libfftw3-dev

# macOS (with Homebrew)
brew install cmake fftw

# Windows: Use vcpkg or Visual Studio package manager
```

**2. Clone and Build:**
```bash
git clone https://github.com/SoftSimu/SymPhas.git
cd SymPhas
mkdir build && cd build

# Basic build with examples
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMAIN_FILE="../examples/tutorial/main.cpp"

make
./symphas_impl
```


> **Need more details?** Continue reading for comprehensive installation instructions, troubleshooting, and advanced configuration options.

This quick start guide continues in [Getting Started Tutorial](#getting-started-tutorial).

---

### Overview

SymPhas is a cross-platform C++ library that uses CMake for configuration and building. This comprehensive guide covers all aspects of installation, from dependencies to platform-specific considerations and troubleshooting.

### Prerequisites and Dependencies

#### Core Requirements

- **C++17 compatible compiler** (see supported versions below)
- **CMake 3.22 or higher** (3.14 minimum, but 3.22+ recommended)

#### Supported Compilers

| Compiler | Minimum Version | Status | Notes |
|----------|----------------|---------|-------|
| GCC | 7.5 | ✅ Fully Supported | Recommended: GCC 10+ |
| Clang | 11.0.1 | ✅ Fully Supported | Good performance |
| MSVC | Any | ❌ Not Typically Compatible | Template compilation issues |

> **Important**: SymPhas uses C++17 features extensively. Ensure your compiler fully supports the C++17 standard.

#### Required Dependencies

- **FFTW3** (version 3.3.7+)
  - **Purpose**: Fourier Transform solvers (`SolverFT`, spectral methods)
  - **Note**: Can be disabled, but eliminates Fourier transform capabilities

#### Optional Dependencies

| Package | Version | Purpose | Installation Priority |
|---------|---------|---------|---------------------|
| **MPI** | Any recent | Parallel computing | High (for large simulations) |
| **CUDA Toolkit** | 10.0+ | GPU acceleration | Medium (performance boost) |
| **VTK** | 9.0+ | Real-time visualization | Medium (convenience) |
| **xdrfile** | 2.1.2+ | Binary XDR output (Gromacs format) | Low (specialized output) |

### Platform-Specific Prerequisites

#### Linux Prerequisites

**Ubuntu/Debian:**
```bash
# Core build tools
sudo apt update
sudo apt install build-essential cmake git

# GCC (if not using system default)
sudo apt install gcc-10 g++-10

# FFTW3 (required)
sudo apt install libfftw3-dev

# Optional dependencies
sudo apt install libopenmpi-dev    # MPI support
sudo apt install libvtk9-dev       # VTK for visualization
```

**RHEL/CentOS/Fedora:**
```bash
# Core build tools
sudo dnf groupinstall "Development Tools"
sudo dnf install cmake git

# FFTW3 (required)
sudo dnf install fftw3-devel

# Optional dependencies
sudo dnf install openmpi-devel vtk-devel
```

**From Source (Advanced):**
If system packages are outdated or unavailable:
```bash
# FFTW3 from source (recommended method)
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j$(nproc)
sudo make install
```

#### Windows Prerequisites

**Required Software:**
1. **Visual Studio 2019+** with "Desktop Development with C++" workload
2. **CMake for Windows** (3.22+)
3. **Git for Windows**

**Dependency Installation:**
- Use **vcpkg** (recommended) or build dependencies from source
- **vcpkg installation:**
  ```cmd
  git clone https://github.com/Microsoft/vcpkg.git
  cd vcpkg
  .\bootstrap-vcpkg.bat
  .\vcpkg install fftw3:x64-windows
  .\vcpkg install openmpi:x64-windows  # Optional
  .\vcpkg install vtk:x64-windows      # Optional
  ```
**How to make CMake and your executable find all vcpkg-installed libraries**

1. *Tell CMake where to find all vcpkg packages*:
When running CMake, always add the following argument (replace the path with your vcpkg location):
    ```
    -DCMAKE_TOOLCHAIN_FILE=C:/Users/<youruser>/vcpkg/scripts/buildsystems/vcpkg.cmake
    ```

    This allows CMake to automatically find all libraries installed by vcpkg (FFTW, OpenMPI, VTK, etc.) without setting each `*_DIR` variable manually.

2. *Ensure Windows can find DLLs at runtime*:
Many libraries installed by vcpkg (like FFTW) are provided as DLLs. To avoid “DLL not found” errors:

    - Add the vcpkg `bin` directory to your system or user `PATH` environment variable:

      ```
      C:\Users\<youruser>\vcpkg\installed\x64-windows\bin
      ```
    - This allows Windows to find `fftw3.dll` and other DLLs automatically, no matter where your executable is located.
    - (Alternative: You can copy DLLs to your executable directory, but this is not recommended for portability.)

3. *Summary*: 

    Use the vcpkg toolchain file for CMake to find all dependencies.
    
    Add the vcpkg `bin` directory to your `PATH` so DLLs are found at runtime.

#### macOS Prerequisites

**Using Homebrew:**
```bash
# Install Xcode command line tools
xcode-select --install

# Install dependencies
brew install cmake git
brew install fftw
brew install open-mpi    # Optional
brew install vtk         # Optional
```

### Installation Process

#### 1. Download SymPhas

```bash
git clone https://github.com/SoftSimu/SymPhas.git
cd SymPhas
git checkout main  # or dev for latest features
```

#### 2. Prepare Build Directory

```bash
mkdir build
cd build
```

#### 3. CMake Configuration

**SymPhas relies on a driver file to execute specific models with chosen solvers**. The driver file is specified with the cmake configuration option `MAIN_FILE`. It is also possible to install SymPhas with only the core header and library files (in order to be linked to from a separate directory).

The installation process varies depending on your intended use case:

**Simulation Installation (With Examples):**
```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install \
      -DCMAKE_BUILD_TYPE=Release \
      -DSOLVER_INCLUDE_HEADER_DIR=../examples/solvers \
      -DMODEL_INCLUDE_HEADER_DIR=../examples/models \
      -DSOLVER_INCLUDE_HEADER_NAME=solverinclude.h \
      -DMODEL_INCLUDE_HEADER_NAME=modelinclude.h \
      -DMAIN_FILE=../examples/simultaneous-configs/main.cpp \
      ..
```
> ℹ️ The driver file is `examples/simultaneous-configs/main.cpp`, which only needs to include the `symphas.h` header as the models and solvers are automatically included from there.

**Basic Installation (Headers and Core Library):**
```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install \
      -DCMAKE_BUILD_TYPE=Release \
      ..
```
> ⚠️ This does not compile with any driver file or models that can be run!

**Windows with vcpkg:**
```cmd
cmake -DCMAKE_INSTALL_PREFIX=C:\path\to\install ^
      -DCMAKE_TOOLCHAIN_FILE=C:\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake ^
      -DCMAKE_BUILD_TYPE=Release ^
      ..
```

#### 4. Build and Install

**Linux/macOS:**
```bash
make -j$(nproc)
make install
```

**Windows:**
```cmd
cmake --build . --config Release --parallel
cmake --install . --config Release
```

### Configuration Options

#### Essential CMake Variables

- **`MAIN_FILE`**: Path to your main driver file (required for executable builds)
- **`SOLVER_INCLUDE_HEADER_DIR`**: Directory with solver definitions (required when using built-in solvers)
- **`SOLVER_INCLUDE_HEADER_NAME`**: Solver include file name (e.g., "solverinclude.h")
- **`MODEL_INCLUDE_HEADER_DIR`**: Directory with model definitions
- **`MODEL_INCLUDE_HEADER_NAME`**: Model include file name

#### Performance Optimization

- **`AVAILABLE_DIMENSIONS`**: Specify dimensions (1, 2, 3) to reduce compile time
  ```bash
  -DAVAILABLE_DIMENSIONS=2  # Only 2D simulations
  ```

- **`AVAILABLE_STENCILS_BASIC_ONLY=ON`**: Use subset of stencils for faster compilation

#### Stencil Configuration

Control which numerical stencils are compiled:

```bash
# Basic stencil control
-DAVAILABLE_STENCIL_POINTS_LAPLACIAN_2D="5,9" \
-DAVAILABLE_STENCIL_POINTS_BILAPLACIAN_2D="13,17"

# Advanced: Override system
-DENABLE_STENCIL_OVERRIDE=ON
```

### Troubleshooting

#### Common Issues and Solutions

**1. FFTW Not Found**
```
Error: Could NOT find FFTW3
```
**Solution:**
```bash
# Specify FFTW location manually
cmake -DFFTW3_DIR=/path/to/fftw3/lib/cmake/fftw3/ ..
```

**2. Compiler Version Issues**
```
Error: C++17 features not supported
```
**Solutions:**
- Update compiler to supported version
- Explicitly set C++17 standard:
  ```bash
  cmake -DCMAKE_CXX_STANDARD=17 ..
  ```

**3. Template Compilation Errors**
```
Error: Template instantiation depth exceeded
```
**Solutions:**
- Use `AVAILABLE_STENCILS_BASIC_ONLY=ON`
- Limit dimensions with `AVAILABLE_DIMENSIONS`
- Increase template depth: `-DCMAKE_CXX_FLAGS="-ftemplate-depth=1000"`

**4. Windows vcpkg Integration**
```
Error: Package not found
```
**Solution:**
```cmd
cmake -DCMAKE_TOOLCHAIN_FILE=C:\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake ..
```

**5. MPI Compilation Issues**
**Solution:**
```bash
# Use MPI compiler wrappers
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
```

#### Performance Considerations

**Compile Time Optimization:**
- Use `AVAILABLE_STENCILS_BASIC_ONLY=ON` for faster builds
- Limit `AVAILABLE_DIMENSIONS` to needed dimensions
- Consider `COMBINE_SHARED=ON` for single shared library

**Runtime Performance:**
- Use `CMAKE_BUILD_TYPE=Release` for production
- Enable compiler optimizations: `-DCMAKE_CXX_FLAGS="-O3 -march=native"`
- Link with optimized FFTW and MPI libraries

#### Debug Builds

For development and debugging:
```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="-g -O0" \
      ..
```

### Installation Verification

After installation, verify SymPhas works correctly:

**1. Check Headers:**
```bash
find /path/to/install -name "symphas.h"
```

**2. Test Compilation:**
Create a minimal test file and compile it against your installation.

**3. Run Examples:**
```bash
cd examples/tutorial
mkdir build && cd build
cmake -DSymPhas_DIR=/path/to/install/lib/cmake/SymPhas ..
make
./tutorial_example
```

---

## Getting Started Tutorial

### Prerequisites

Before starting this tutorial, ensure you have completed the [Installation Guide](#installation-guide) above.

### CMake Configuration for Projects

When building SymPhas projects, you need to configure specific CMake variables:

**Essential Variables:**
- **`MAIN_FILE`**: Path to your main driver file containing `main()` function

**Solver Configuration (required when using built-in solvers):**
- **`SOLVER_INCLUDE_HEADER_DIR`**: Directory containing solver definitions
- **`SOLVER_INCLUDE_HEADER_NAME`**: Solver include file name (e.g., "solverinclude.h")

**Optional Variables:**
- **`MODEL_INCLUDE_HEADER_DIR`**: Directory containing model definitions  
- **`MODEL_INCLUDE_HEADER_NAME`**: Model include file name

> **Note**: When using built-in solvers like `SolverFT` or `SolverSP`, you must specify `SOLVER_INCLUDE_HEADER_DIR` and typically `SOLVER_INCLUDE_HEADER_NAME`, or the build will fail.

### Tutorial: Your First SymPhas Simulation

This tutorial uses the minimal Allen-Cahn example from `examples/tutorial/` to demonstrate basic SymPhas usage.

#### 1. Set Up the Tutorial Example

The tutorial example contains a complete minimal driver file (`main.cpp`):
```cpp
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
```

3. **Build and run the existing tutorial**:
```bash
# Create build directory (if not already exists)
mkdir build && cd build

# Configure with CMake (specify path to your main file)
# Note: SOLVER_INCLUDE_HEADER_DIR is required when using built-in solver objects like SolverFT
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DMAIN_FILE="../examples/tutorial/main.cpp" \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"

# Build
make symphas_impl

# Run the simulation
./symphas_impl
```
>⚠️ The path to the main file and solver directory is configured relative to the path you are configuring from, not the root cmake file (in the root directory of SymPhas).

**Expected Output**:
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
Simulation completed. Field range: -0.XXX to 0.XXX
```

**Understanding this Simple Example**:

- **Model Definition**: Defines Allen-Cahn model directly in main file using `MODEL` macro
- **Grid Setup**: Creates a 64×64 grid with periodic boundaries using [Problem Spec Notation](#problem-spec-notation)
- **Initial Conditions**: Sets random uniform initial values between -0.1 and 0.1
- **Solver**: Uses explicit finite difference solver (`SolverFT`) with 2nd order stencils
- **Output**: Prints basic field statistics after simulation




### More Advanced Examples

For more comprehensive examples that demonstrate different aspects of SymPhas:

**Advanced Features** (`examples/advanced-features/`): Demonstrates:
- **Free Energy Formulation**: Defining models using free energy functionals
- **I/O Integration**: Saving field data to files using the `io` module  
- **Field Analysis**: Computing and displaying simulation statistics
- **Larger Scale**: 200×200 grid simulation

**Comprehensive Demo** (`examples/comprehensive-demo/`): Showcases:
- **Multiple Model Approaches**: Evolution equations, free energy, custom dynamics
- **Various Initial Conditions**: Uniform, lambda functions, expressions, geometric shapes
- **Multi-Field Systems**: Both single and multi-field models with mixed types
- **Modular Design**: Separate header files and organized code structure

### Building Your Own Project

<!-- To create your own SymPhas project, you have two main approaches: -->

<!-- #### Approach 1: Using the SymPhas Build System -->

1. **Create project directory inside SymPhas**:
```bash
mkdir myproject
cd myproject
```

2. **Create your main file** (`main.cpp`):
```cpp
#include "solverinclude.h"
#include "symphas.h"

// Define your model
MODEL(MyModel, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - power(op(1), 3))
)

int main() {
    using namespace symphas;
    
    // Setup parameters
    auto interval = BoundaryType::PERIODIC || 64_h / 0.1_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::UNIFORM <<= {-0.1, 0.1});
    
    // Create and run model
    model_MyModel_t<2, SolverFT<Stencil2d2h<>>> model{params};
    symphas::find_solution(model, 0.01, 1000);
    
    return 0;
}
```

3. **Return to the project root**:
```bash
cd ../
```

4. **Build and run**:
```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMAIN_FILE="../myproject/main.cpp"
make
./myproject
```
> The built-in finite difference solver is linked to the build.

<!-- 
#### Approach 2: Standalone Installation

1. **Install SymPhas system-wide**:
```bash
cd SymPhas
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DUSE_IO=ON \
    -DUSE_CONF=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"
make install
```

2. **Create standalone project**:
```bash
mkdir myproject && cd myproject
# Create main.cpp and CMakeLists.txt as above
```

3. **CMakeLists.txt for standalone**:
```cmake
cmake_minimum_required(VERSION 3.22)
project(MyProject)

find_package(SymPhas REQUIRED)
add_executable(myproject main.cpp)
target_link_libraries(myproject SymPhas::SymPhas)
```
 -->

### Next Steps

After completing this tutorial, you're ready to explore SymPhas further:

1. **Learn More About Models**: Read [Creating Your Own Models](#creating-your-own-models) to understand:
   - Evolution equations vs free energy formulation
   - Multi-field systems
   - Complex field types
   - Provisional variables

2. **Explore Examples**: Check out the comprehensive examples:
   - `examples/advanced-features/` — Free energy models and I/O
   - `examples/comprehensive-demo/` — Multi-field systems  
   - `examples/models/` — Various phase-field models
   - `examples/solvers/` — Different numerical solvers

3. **Configuration Options**: Learn about [Build Configuration Options](#build-configuration-options) for:
   - Performance optimization
   - Different stencil configurations
   - MPI parallelization
   - GPU acceleration

4. **API Reference**: Consult the [Symbols and Functions Reference](#symbols-and-functions-reference) for detailed documentation of all SymPhas functions and operators.

---

- **Explore Examples**: Check `examples/` directory for more complex models
- **Read Model Grammar**: Learn advanced model definition syntax
- **Configuration Files**: Use JSON configuration for complex setups
- **Parallel Computing**: Enable MPI for large-scale simulations
- **Visualization**: Set up VTK for real-time visualization


---


## Build Configuration Options

SymPhas provides extensive compile-time configuration through CMake variables. Understanding these options is crucial for setting up functional builds.

#### Required Configuration Parameters

**Essential for any executable build:**
```bash
# Solver configuration - REQUIRED
-DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers"  # Path to solver directory
-DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h"     # Solver include file name

# Main file - REQUIRED for standalone executables
-DMAIN_FILE="../myproject/main.cpp"                          # Path to main driver file
```

**Optional but recommended:**
```bash
# Model configuration - optional if models defined in main file
-DMODEL_INCLUDE_HEADER_DIR="../examples/models"    # Path to model directory  
-DMODEL_INCLUDE_HEADER_NAME="modelinclude.h"       # Model include file name

# I/O support - highly recommended
-DUSE_IO=ON                                         # Enable file I/O capabilities
```

#### Core Module Configuration

```bash
# Major modules
-DUSE_IO=ON                    # Enable I/O module (reading/writing data files)
-DUSE_CONF=ON                  # Enable configuration module (JSON configs)
-DUSE_VTK=ON                   # Enable VTK visualization support
-DUSE_FFTW3=ON                 # Enable FFTW3 library integration

# Build type and optimization
-DCMAKE_BUILD_TYPE=Release     # Release, Debug, RelWithDebInfo
-DCOMBINE_SHARED=ON            # Create single combined library
```

#### Advanced Solver and Dimension Configuration

```bash
# Available dimensions (affects compile time significantly)
-DAVAILABLE_DIMENSIONS="1,2,3"                     # Compile for 1D, 2D, and 3D
-DAVAILABLE_DIMENSIONS="2"                       # Only 2D (faster compile)

# Stencil accuracy per dimension
-DAVAILABLE_STENCIL_ACCURACY_1D="2,4,6"           # 2nd, 4th, 6th order accuracy in 1D
-DAVAILABLE_STENCIL_ACCURACY_2D="2,4"             # 2nd and 4th order in 2D  
-DAVAILABLE_STENCIL_ACCURACY_3D="2"               # Only 2nd order in 3D (recommended)

# Stencil optimization options
-DAVAILABLE_STENCILS_AUTOGENERATION=ON            # Auto-generate needed stencils (slow compile)
-DAVAILABLE_STENCILS_BASIC_ONLY=ON                # Use only basic stencil set (fast compile)
```

#### Stencil Override System

The stencil override system provides fine-grained control over which finite difference stencils are available at build time through CMake cache variables. This system automatically discovers available stencils and provides individual boolean controls for maximum flexibility.

**Enable the Override System:**
```bash
-DENABLE_STENCIL_OVERRIDE=ON                       # Enable individual stencil selection
```

When enabled, this creates boolean cache variables for each stencil implementation from:
- `stencilh2.h` - Second-order accurate stencils
- `stencilh4.h` - Fourth-order accurate stencils

**Available Stencils by Dimension and Order:**

*1D Second-order (2H):*
- Laplacian: 3 points
- Bilaplacian: 5 points  
- Gradlaplacian: 4 points

*2D Second-order (2H):*
- Laplacian: 5, 9 points
- Bilaplacian: 13, 17, 21 points
- Gradlaplacian: 6, 8, 12, 16 points

*2D Fourth-order (4H):*
- Laplacian: 9, 17, 21 points
- Bilaplacian: 21, 25, 33, 37 points
- Gradlaplacian: 14, 18, 26, 30 points

*3D Second-order (2H):*
- Laplacian: 7, 15, 19, 21, 27 points
- Bilaplacian: 21, 25, 41, 52, 57 points
- Gradlaplacian: 10, 12, 28, 36, 40 points

**Individual Stencil Selection:**
```bash
# Example: Enable only specific 2D stencils
-DENABLE_STENCIL_OVERRIDE=ON \
-DSTENCIL_LAPLACIAN_2D_2H_5PTS=ON \
-DSTENCIL_LAPLACIAN_2D_2H_9PTS=OFF \
-DSTENCIL_BILAPLACIAN_2D_2H_17PTS=ON \
-DSTENCIL_GRADLAPLACIAN_2D_2H_6PTS=ON
```

**CMake GUI Workflow:**
1. Configure with `-DENABLE_STENCIL_OVERRIDE=ON`
2. Open CMake GUI
3. Check/uncheck desired stencil checkboxes
4. Configure and generate

**Benefits:**
- **Fine-grained Control**: Enable only the stencils you need
- **Reduced Binary Size**: Fewer template instantiations
- **Faster Compilation**: Less code to compile

**Integration with JSON Configuration:**
The system works seamlessly with JSON stencil validation:

```json
{
  "simulation": {
    "stencil": {
      "order": 2,
      "ptl": 5,    // Validated against selected Laplacian stencils
      "ptg": 12,   // Validated against selected Gradlaplacian stencils  
      "ptb": 17    // Validated against selected Bilaplacian stencils
    }
  }
}
```

#### Parallel Computing Configuration

```bash
# MPI support
-DUSING_MPI=ON                                     # Enable MPI parallelization
-DMPI_CXX_COMPILER=mpicxx                         # Specify MPI compiler wrapper

# CUDA support (experimental)
-DUSING_CUDA=ON                                    # Enable CUDA GPU support
-DSYMPHAS_CUDA_ARCHITECTURES="70;75;80"          # Target CUDA architectures
```

#### Output and Debugging Configuration

```bash
# Equation printing and formatting
-DPRINTABLE_EQUATIONS=ON                          # Enable equation text output
-DUSE_LATEX_FORMAT=ON                             # Use LaTeX formatting for equations
```

#### Complete Configuration Examples
> The following example configurations assume you wrote a driver file `main.cpp` in the directory `myproject` relative to the SymPhas root path.

**Minimal working configuration:**
```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMAIN_FILE="../myproject/main.cpp"
```

**Full-featured development configuration:**
```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_IO=ON \
    -DUSE_CONF=ON \
    -DUSE_VTK=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMODEL_INCLUDE_HEADER_DIR="../examples/models" \
    -DMODEL_INCLUDE_HEADER_NAME="modelinclude.h" \
    -DMAIN_FILE="../myproject/main.cpp" \
    -DAVAILABLE_DIMENSIONS="2,3" \
    -DAVAILABLE_STENCIL_ACCURACY_2D="2,4" \
    -DAVAILABLE_STENCIL_ACCURACY_3D="2" \
    -DAVAILABLE_STENCILS_BASIC_ONLY=ON \
    -DPRINTABLE_EQUATIONS=ON
```

**High-performance computing configuration:**
```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMAIN_FILE="../myproject/main.cpp" \
    -DAVAILABLE_DIMENSIONS="3" \
    -DAVAILABLE_STENCIL_ACCURACY_3D="2;4" \
```

**Custom stencil selection configuration:**
```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_IO=ON \
    -DUSE_CONF=ON \
    -DSOLVER_INCLUDE_HEADER_DIR="../examples/solvers" \
    -DSOLVER_INCLUDE_HEADER_NAME="solverinclude.h" \
    -DMAIN_FILE="../main.cpp" \
    -DAVAILABLE_DIMENSIONS="2" \
    -DENABLE_STENCIL_OVERRIDE=ON \
    -DSTENCIL_LAPLACIAN_2D_2H_5PTS=ON \
    -DSTENCIL_LAPLACIAN_2D_2H_9PTS=OFF \
    -DSTENCIL_BILAPLACIAN_2D_2H_13PTS=OFF \
    -DSTENCIL_BILAPLACIAN_2D_2H_17PTS=ON \
    -DSTENCIL_GRADLAPLACIAN_2D_2H_6PTS=ON \
    -DSTENCIL_GRADLAPLACIAN_2D_2H_12PTS=ON
```

#### Troubleshooting Configuration Issues

**Common problems and solutions:**
1. **"Could not find solver"**:
   - Verify `SOLVER_INCLUDE_HEADER_DIR` path is correct
   - Ensure `SOLVER_INCLUDE_HEADER_NAME` file exists in specified directory

2. **"Unknown model provided"**:
   - Either define models in main file, or set `MODEL_INCLUDE_HEADER_*` parameters
   - Check that model files contain properly formatted `MODEL()` definitions
   - Check that the model has actually been linked with `LINK_WITH_NAME`

3. **"FFTW3 not found"**:
   - Install FFTW3 library or set `-DUSE_FFTW3=OFF` if not needed
   - Point CMake to FFTW3 installation: `-DFFTW3_DIR=/path/to/fftw3`

4. **Very long compile times**:
   - Use `-DAVAILABLE_STENCILS_BASIC_ONLY=ON`
   - Reduce dimensions: `-DAVAILABLE_DIMENSIONS="2"`
   - Avoid `-DAVAILABLE_STENCILS_AUTOGENERATION=ON`
   - Use stencil override system with specific selections: `-DENABLE_STENCIL_OVERRIDE=ON`

5. **"Invalid stencil point count"**:
   - When using stencil override system, ensure JSON `ptl`/`ptg`/`ptb` values match enabled stencils
   - Check that `ENABLE_STENCIL_OVERRIDE=ON` and required stencil options are enabled
   - Use `cmake-gui` to verify which stencils are currently enabled

6. **"Stencil not found" compilation errors**:
   - Verify that required stencil implementations exist in `stencilh2.h` or `stencilh4.h`
   - Ensure stencil override selections match available implementations
   - Check CMake output for discovered stencil configurations

---

## Architecture and Design

### System Architecture

SymPhas follows a modular, layered architecture with clear separation of concerns:

```
┌─────────────────────────────────────────────────────────┐
│                    Driver Layer                         │
│              (User Applications)                        │
├─────────────────────────────────────────────────────────┤
│                   Model Layer                           │
│        (Phase Field Problem Definitions)                │
├─────────────────────────────────────────────────────────┤
│                  Solver Layer                           │
│           (Numerical Solution Methods)                  │
├─────────────────────────────────────────────────────────┤
│               Symbolic Algebra                          │
│        (Expression Trees & Operations)                  │
├─────────────────────────────────────────────────────────┤
│                 Data Types                              │
│         (Grids, Boundaries, Arrays)                     │
├─────────────────────────────────────────────────────────┤
│                Core Library                             │
│     (Basic Functions & Utilities)                       │
└─────────────────────────────────────────────────────────┘
```

### Design Principles

1. **Template Meta-Programming**: Extensive use of C++ templates for compile-time optimization and typing system prevents errors at runtime
2. **Zero-Cost Abstractions**: High-level abstractions with no runtime performance overhead
3. **Modularity**: Each module has well-defined responsibilities and interfaces
4. **Extensibility**: Framework designed for easy addition of new solvers and models

---

## Module Documentation

### Core Modules (Required)

#### 1. **lib** Module
- **Purpose**: Basic elements and functions used by other modules
- **Namespace**: `symphas`
- **Key Components**:
  - Program parameters management
  - Mathematical constants and utilities

#### 2. **datatypes** Module
- **Purpose**: Basic objects for finite difference numerical solvers
- **Namespace**: `grid`
- **Key Components**:
  - `Grid<T, D>` N-dimensional data grids
  - Boundary condition management
  - Grid iterators and accessors
  - CUDA grid implementations

#### 3. **sym** Module
- **Purpose**: Symbolic algebra elements and functionality
- **Namespace**: `expr`
- **Key Components**:
  - Expression tree construction
  - Symbolic differentiation
  - Expression simplification rules
  - Compile-time expression evaluation
  - Custom operator definitions

#### 4. **sol** Module
- **Purpose**: Consolidates sym and datatypes for solution framework
- **Key Components**:
  - `Solver<Sp, N>` class hierarchy
  - `Model<D, Sp, S...>` class templates
  - Finite difference stencils (`Stencil<O, D>`)
  - Time integration schemes
  - Phase field system management

### Optional Modules

#### 5. **io** Module (Optional)
- **Purpose**: Input and output functionality
- **Namespace**: `symphas::io`
- **Key Components**:
  - File I/O operations
  - Data serialization
  - Checkpoint management
  - Visualization output formats
  - Result export utilities

#### 6. **conf** Module (Optional)
- **Purpose**: Configuration file management
- **Dependencies**: Requires `io` module
- **Key Components**:
  - Configuration file parsing
  - Parameter management
  - Simulation setup from config files
  - Runtime configuration updates

---

## Input/Output (IO) System

### Overview

SymPhas provides a comprehensive input/output (IO) system that handles data persistence, visualization, and checkpoint management for phase-field simulations. The IO system supports multiple output formats, configurable save parameters, and automatic checkpoint creation for long-running simulations.

**Key Features:**
- **Multiple Output Formats**: GNU/gnuplot structured matrix, CSV, and XDR binary formats
- **Flexible Save Parameters**: Control when and how data is saved during simulations
- **Checkpoint System**: Automatic checkpointing with restart capabilities
- **Data Organization**: Automatic directory structure and file naming
- **Field Type Support**: All field types (scalar, complex, vector) with appropriate formatting

### Output Format Selection

SymPhas provides three main output formats optimized for different use cases:

| Format | Description | Use Cases | Requirements |
|--------|-------------|-----------|--------------|
| **GNU** | Structured matrix format | Data visualization with gnuplot, mathematical analysis | [gnuplot](http://www.gnuplot.info/) software (optional) |
| **CSV** | Comma-separated values | Spreadsheet analysis, general data processing | None |
| **XDR** | Binary format (Gromacs-compatible) | Large datasets, space-efficient storage | [libxdrfile](https://github.com/wesbarnett/libxdrfile) library |

**Default Behavior**: SymPhas uses the GNU structured matrix format by default.

### Available Output Formats

#### GNU: Structured Matrix Format

The GNU format outputs data as a structured grid where:
- **First row**: y-axis coordinate values  
- **First column**: x-axis coordinate values
- **Data grid**: Field values at corresponding (x,y) coordinates
- **Frame separation**: Two newlines separate time frames
- **Compatibility**: Can be read directly by gnuplot for visualization
- **Vector fields**: Decomposed into normalized direction components plus magnitude

**Field Type Handling:**
- **Scalar fields**: Direct numerical values in scientific notation
- **Complex fields**: Norm of the complex number in scientific notation
- **Vector fields**: Decomposed into coordinates + normalized direction components + magnitude
  - 1D vectors: `x value` (coordinate + component value)
  - 2D vectors: `x y dx dy magnitude` (coordinates + unit vector components + magnitude)
  - 3D vectors: `x y z dx dy dz magnitude` (coordinates + unit vector components + magnitude)

**Format Rules:**
- First entry in grid is always 0
- Ghost cells are excluded
- Multi-dimensional data is flattened appropriately for the format

**Scalar Field Example:**
```
            0          3.00          4.00          5.00          6.00
         3.00  9.533376E-01  7.855586E-01  9.237062E-01 -2.333816E-02
         4.00 -9.023678E-01  3.934696E-01 -5.117627E-01  2.793115E-01
         5.00 -1.757295E-01 -8.584688E-02 -1.780544E-01  6.065830E-01
         6.00 -8.881749E-01 -4.034086E-01  1.099667E-01 -1.756146E-02

            0          3.00          4.00          5.00          6.00
         3.00  9.534012E-01  7.856234E-01  9.237891E-01 -2.332145E-02
         4.00 -9.022345E-01  3.935123E-01 -5.116789E-01  2.794567E-01
         5.00 -1.756123E-01 -8.583456E-02 -1.779345E-01  6.067123E-01
         6.00 -8.880567E-01 -4.033234E-01  1.100456E-01 -1.754789E-02
```

**2D Vector Field Example:**
```
3.00 3.00 0.70711 0.70711 1.41421
4.00 3.00 0.89443 0.44721 2.23607
5.00 3.00 0.60000 0.80000 2.50000
3.00 4.00 0.70711 0.70711 1.41521
4.00 4.00 0.89443 0.44721 2.23707
5.00 4.00 0.60000 0.80000 2.50100

3.00 3.00 0.70712 0.70712 1.41422
4.00 3.00 0.89444 0.44722 2.23608
5.00 3.00 0.60001 0.80001 2.50001
3.00 4.00 0.70712 0.70712 1.41522
4.00 4.00 0.89444 0.44722 2.23708
5.00 4.00 0.60001 0.80001 2.50101
```
*Note: Format is `x y dx dy magnitude` where (dx, dy) are normalized direction components. One line per grid point, separated by blank lines between time frames.*

#### CSV: Comma-Separated Data

The CSV format provides maximum compatibility with data analysis tools:
- **Standard delimiters**: Comma-separated values
- **Field handling**: Complex numbers as `"real,imag"`, vectors as quoted strings with space-separated components
- **Clean formatting**: Proper quoting and escaping for all data types

**Complex Field Example:**
```csv
"1.23,0.45","2.34,1.56","3.45,2.67"
"4.56,3.78","5.67,4.89","6.78,5.91"
```

**Vector Field Examples:**
```csv
# 2D vectors: "dx dy magnitude" format (normalized direction + magnitude)
"0.7071 0.7071 1.414","0.8944 0.4472 2.236","0.6000 0.8000 2.500"
"0.9487 0.3162 3.162","0.5547 0.8321 3.606","0.7071 0.7071 4.243"

# 3D vectors: "dx dy dz magnitude" format
"0.577 0.577 0.577 1.732","0.267 0.535 0.802 3.742","0.378 0.756 0.535 4.583"
"0.447 0.894 0.000 4.472","0.516 0.688 0.516 5.196","0.408 0.816 0.408 5.745"
```

#### XDR: Binary Format

The XDR format provides space-efficient binary storage:
- **Binary encoding**: Efficient storage for large datasets
- **Gromacs compatibility**: Uses standard XDR format from computational chemistry
- **Header information**: Grid dimensions, coordinates, and metadata
- **Cross-platform**: Consistent byte ordering across systems

### Setting Output Format

#### Using Program Spec Notation

The most direct way to set the output format is using Program Spec Notation in your main driver file:

```cpp
#include "symphas.h"

using namespace symphas;

int main() {
    // Set output format to CSV
    PARAMS += WRITER << IOType::CSV;
    
    // Alternative: Set both reader and writer
    PARAMS += READER_AND_WRITER << IOType::GNU;
    
    // Your simulation code here...
    return 0;
}
```

**Available IO Types:**
- `IOType::GNU` - Structured matrix format (default)
- `IOType::CSV` - Comma-separated values  
- `IOType::XDR` - Binary XDR format
- `IOType::COLUMN` - Column-style plain text
- `IOType::MOVIE` - Movie/animation output format

#### Using Command Line Arguments

When parsing command line parameters with `symphas::init()`, use the `--writer` flag:

```bash
# Set output format to CSV
./my_simulation --writer=csv

# Set output format to XDR binary
./my_simulation --writer=xdr

# Set both reader and writer
./my_simulation --reader-and-writer=gnu
```

**Command Line Options:**
- `--writer=TYPE` or `-w TYPE` - Set output format
- `--reader=TYPE` or `-r TYPE` - Set input format  
- `--reader-and-writer=TYPE` or `-x TYPE` - Set both formats

where `TYPE` is one of: `gnu`, `csv`, `xdr`, `column`, `movie`

### SaveParams Configuration

The `SaveParams` class controls when and how often simulation data is saved. It provides a standardized way to define saving intervals and pass them to solution functions.

#### Basic SaveParams Usage

```cpp
// Constructor: SaveParams(save_interval, stop_iteration)
SaveParams save_params(1000, 100000);

// This means:
// - Save data every 1000 simulation steps
// - Stop simulation after 100,000 iterations
```

#### Advanced SaveParams Configuration

```cpp
// More detailed constructor
SaveParams save_params(
    SaveType::DEFAULT,    // Save type
    1000.0,               // Save interval (as double)
    0,                    // Start index
    100000                // Stop index
);

// Using with different save types
SaveParams list_save(SaveType::LIST, "1000,5000,10000,25000,50000");
SaveParams count_save(SaveType::COUNT, 20, 0, 100000);  // Save 20 times total
```

#### SaveParams Methods

```cpp
SaveParams save(1000, 50000);

// Check if current iteration should save
bool should_save = save.is_save_index(current_iteration);

// Get next save point
iter_type next_save = save.next_save(current_iteration);

// Check if simulation should stop
bool is_done = save.is_last_save(current_iteration);

// Get stop iteration
iter_type stop_iter = save.get_stop();
```

#### Integration with find_solution()

```cpp
// Use SaveParams with simulation functions
SaveParams save(500, 10000);  // Save every 500 steps, stop at 10,000

// Single model simulation
find_solution(model, dt, save);

// Multi-step simulation with specific save directory
find_solution(model, dt_list, save, "output_directory");
```

### Checkpoint System

SymPhas provides an automatic checkpoint system for simulation recovery and restart capabilities.

#### Checkpoint Overview

**Purpose**: Checkpoints save the complete simulation state at specified intervals, allowing:
- **Recovery**: Resume simulations after interruption
- **Analysis**: Access intermediate states for detailed study
- **Workflow**: Break long simulations into manageable segments

**Default Behavior**: 
- Checkpoints are saved 10 times during simulation
- For N total steps, checkpoints occur every N/10 steps
- Stored in plain text format for easy inspection

#### Checkpoint Functions

**Automatic Checkpointing:**
```cpp
// Check and save checkpoint if needed
symphas::save_if_checkpoint(model, save_params, "checkpoint_directory");

// This function:
// 1. Checks if current iteration matches checkpoint schedule
// 2. Saves checkpoint if needed
// 3. Records index for restart information
```

**Manual Checkpointing:**
```cpp
// Save checkpoint for all systems in a model
symphas::checkpoint_systems(model.systems_tuple(), "checkpoint_dir", iteration);

// Save checkpoint for individual system
symphas::checkpoint_system(system, "checkpoint_dir", iteration);
```

#### Checkpoint Configuration

**Command Line Control:**
```bash
# Set number of checkpoints (default: 10)
./simulation --checkpoint-count=20

# Disable checkpoints
./simulation --checkpoint-count=0

# Load from existing checkpoint
./simulation --checkpoint=path/to/checkpoint,index
```

**Programmatic Control:**
```cpp
// Set checkpoint frequency
params::checkpoint_count = 15;  // Save 15 checkpoints total

// Disable checkpointing
params::checkpoint_count = 0;
```

#### Checkpoint File Format

Checkpoint files are stored in the `checkpoint/` directory with the naming format `data<N>` where N is the field index.

**File Structure:**
```
# Line 1: Info line
iteration_index dimension_x dimension_y dimension_z dx dy dz

# Lines 2+: Field values (one per line, organized by x-slice)
field_value_1
field_value_2
...
```

**Loading Checkpoints:**
```bash
# Load most recent checkpoint from directory
./simulation --checkpoint=path/to/checkpoint

# Load specific iteration
./simulation --checkpoint=path/to/checkpoint,5000
```

### Program Spec Notation for IO

SymPhas uses Program Spec Notation for setting various IO-related parameters:

#### Writer/Reader Configuration

```cpp
using namespace symphas;

// Basic format selection
PARAMS += WRITER << IOType::CSV;
PARAMS += READER << IOType::GNU;
PARAMS += READER_AND_WRITER << IOType::XDR;
```

#### Advanced IO Parameters

```cpp
// File handling options
PARAMS += SINGLE_OUTPUT << true;     // All data in single file
PARAMS += SINGLE_INPUT << true;      // Read from single file
PARAMS += USE_TIMESTAMP << true;     // Use timestamp in directory names

// Boundary handling
PARAMS += EXTEND_BOUNDARY << true;   // Include ghost cells in output dimensions
```

#### Available IO Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `WRITER` | IOType | Output format selection | `IOType::GNU` |
| `READER` | IOType | Input format selection | `IOType::GNU` |
| `READER_AND_WRITER` | IOType | Set both reader and writer | N/A |
| `SINGLE_OUTPUT` | bool | Append all data to single file | `true` |
| `SINGLE_INPUT` | bool | Read all data from single file | `true` |
| `USE_TIMESTAMP` | bool | Use timestamp in output directories | `true` |
| `EXTEND_BOUNDARY` | bool | Include ghost cells in grid dimensions | `false` |
| `PLOTS_ONLY` | bool | Generate only plot configs, no computation | `false` |

### File Organization

SymPhas automatically organizes output files in a structured directory hierarchy:

#### Default Directory Structure

```
output/
├── data/                  # Main simulation data
│   ├── data0.txt   # Field data files
│   ├── data0.txt
│   └── ...
├── checkpoint/            # Checkpoint files
│   ├── data0              # Field 0 checkpoint
│   ├── data1              # Field 1 checkpoint
│   └── ...
└── plot/                  # Plot configuration files
    └── data.gp
```

#### File Naming Conventions

**Data Files:**
- Format: `data_XXXXXX.txt` where XXXXXX is the iteration number
- Extension varies by format: `.txt` (GNU), `.csv` (CSV), `.xdr` (XDR)

**Checkpoint Files:**
- Format: `data<N>` where N is the field index (0-based)
- Plain text format regardless of main output format

**Timestamp Directories:**
When `USE_TIMESTAMP` is enabled, output is organized by execution time:
```
output_2024-01-15_14-30-45/
├── data/
├── checkpoint/
└── plots/
```

#### Custom Directory Control

```cpp
// Specify custom output directory
find_solution(model, dt_list, save_params, "custom_output_dir");

// Directory will be created if it doesn't exist
// Full path: custom_output_dir/data/, custom_output_dir/checkpoint/, etc.
```

### Practical Examples

#### Example 1: Basic CSV Output

```cpp
#include "symphas.h"

using namespace symphas;

MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

int main() {
    // Set output format to CSV
    PARAMS += WRITER << IOType::CSV;
    
    // Configure save parameters: save every 100 steps, stop at 5000
    SaveParams save(100, 5000);
    
    // Set up simulation
    auto interval = BoundaryType::PERIODIC || 64_h / 0.1_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::UNIFORM <<= {-0.1, 0.1});
    
    model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Run with automatic checkpointing
    find_solution(model, 0.01, save, "csv_output");
    
    return 0;
}
```

#### Example 2: High-Performance Binary Output

```cpp
#include "symphas.h"

using namespace symphas;

int main() {
    // Configure for large simulation with binary output
    PARAMS += WRITER << IOType::XDR;           // Binary format for efficiency
    PARAMS += SINGLE_OUTPUT << true;           // Single output file
    PARAMS += USE_TIMESTAMP << true;           // Timestamped directories
    PARAMS += EXTEND_BOUNDARY << false;        // Exclude ghost cells
    
    // Save frequently for analysis
    SaveParams save(50, 100000);  // Save every 50 steps
    
    // Large grid simulation setup
    auto interval = BoundaryType::PERIODIC || 256_h / 0.05_dh;
    auto grid = interval * interval * interval;  // 3D simulation
    auto params = grid << (Inside::UNIFORM <<= {-0.2, 0.2});
    
    // High-performance model
    model_AllenCahn_t<3, SolverFT<Stencil3d2h<>>> model{params};
    
    // Run with custom checkpoint count
    params::checkpoint_count = 20;  // 20 checkpoints
    find_solution(model, 0.001, save, "large_simulation");
    
    return 0;
}
```

#### Example 3: Analysis-Focused Setup

```cpp
#include "symphas.h"

using namespace symphas;

int main() {
    // Set up for detailed analysis
    PARAMS += WRITER << IOType::GNU;           // Gnuplot-compatible
    PARAMS += SINGLE_OUTPUT << false;          // Separate files per timestep
    PARAMS += PLOTS_ONLY << false;             // Include computation
    
    // Custom save schedule for detailed analysis
    SaveParams analysis_save(SaveType::LIST, "100,250,500,1000,2500,5000,10000");
    
    // Medium resolution for detailed tracking
    auto interval = BoundaryType::PERIODIC || 128_h / 0.1_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::CIRCLE <<= 20_c * 0.5_c, Inside::UNIFORM <<= {-1.0, 1.0});
    
    model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Run with analysis configuration
    find_solution(model, 0.01, analysis_save, "analysis_output");
    
    return 0;
}
```

#### Example 4: Checkpoint Recovery

```cpp
#include "symphas.h"

using namespace symphas;

int main(int argc, char* argv[]) {
    PARAMS += WRITER << IOType::CSV;
    
    // Set up for continuation run
    SaveParams save(200, 20000);  // Continue for more steps
    
    auto interval = BoundaryType::PERIODIC || 64_h / 0.1_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::CHECKPOINT);  // Load from checkpoint
    
    model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Continue simulation from checkpoint
    find_solution(model, 0.01, save, "continued_run");
    
    return 0;
}
```

These examples demonstrate the flexibility and power of SymPhas's IO system, from basic data output to sophisticated checkpoint management for long-running simulations.

---

## Core Components

### Grid System

The `Grid<T, D>` class is the fundamental data structure for storing phase field data. It inherits from `Block<T>` and provides dimensional information and multi-dimensional access patterns:

```cpp
template<typename T, size_t D>
struct Grid : Block<T> {
    len_type dims[D];       // Dimensions in each direction
    
    // Constructors
    Grid();
    Grid(const len_type* dimensions);
    Grid(grid::dim_list dimensions);
    
    // Multi-dimensional element access via operator()
    template<typename... Ts>
    decltype(auto) operator()(Ts&&... indices) const;
};
```

**Key Features:**
- **Inherits data storage**: Gets `T* values` and `len_type len` from `Block<T>`
- **Dimensional structure**: Maintains shape information in `dims[D]` array
- **Flattened storage**: Data stored as 1D array internally for memory efficiency
- **Template-based**: Generic over data type `T` and dimension `D`
- **Multi-dimensional access**: Uses `operator()` for coordinate-based indexing

**Specialized Grid Types:**
- `BoundaryGrid<T, D>`: Includes ghost cells for boundary conditions
- `RegionalGrid<T, D>`: For sub-domain grids in multi-field simulations

### Symbolic Algebra System

The symbolic algebra system builds expression trees at compile time:

```cpp
// Basic expressions
auto expr1 = op(1);                   // Phase field 1
auto expr2 = lap(op(1));              // Laplacian of phase field 1
auto expr3 = grad(op(1));             // Gradient of phase field 1
auto expr4 = op(1) * op(1) - 1_n;     // Nonlinear term

// Complex expressions
auto evolution = lap(op(1)) - (op(1) * op(1) - 1_n) * op(1);
```
> The above examples can only be used within a model definition.

### Model Definition System

SymPhas provides a **Domain Specific Language (DSL)** implemented through C++ macros that allows you to define phase field models using natural mathematical syntax. This DSL translates mathematical equations into efficient C++ code at compile time.

**Supported Field Types:**
- **SCALAR**: Real-valued scalar fields (e.g., order parameter φ)
- **COMPLEX**: Complex-valued fields (e.g., for pattern formation, Ginzburg-Landau)
- **VECTOR**: Vector-valued fields (e.g., velocity, displacement, multi-component systems)

Models are defined using the macro-based DSL:

```cpp
// Model A: Allen-Cahn equation (scalar field)
MODEL(ModelA, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

// Complex Ginzburg-Landau equation (complex field)
MODEL(CGL, (COMPLEX),
    EVOLUTION(dop(1) = (1_n + c(1) * Ii) * op(1) + 
                      (1_n + c(2) * Ii) * lap(op(1)) - 
                      abs(op(1), 2) * op(1))
)

// Model C: Two-field system (multiple scalar fields)
MODEL(ModelC, (SCALARS(2)),
    EVOLUTION(
        dop(1) = -bilap(op(1)) - lap((c(1) - c(2)*op(1)*op(1))*op(1)),
        dop(2) = lap(op(2)) + (c(3) - c(4)*op(2)*op(2))*op(2)
    )
)

// Navier-Stokes-like system (scalar + vector)
MODEL(FluidFlow, (SCALAR, VECTOR),
    EVOLUTION(
        dop(1) = lap(op(1)) - div(op(1) * op(2)),  // Scalar transport
        dop(2) = c(1) * lap(op(2)) - grad(op(1))   // Vector momentum
    )
)
```

The DSL enables:
- **Type safety**: Field types are checked at compile time
- **Expression optimization**: Mathematical expressions are optimized during compilation
- **Memory management**: Grid operations are automatically vectorized and optimized

### Solver System


Solvers implement numerical time-stepping schemes. The main solvers provided by SymPhas are:

**Available Built-in Solvers:**
- `SolverFT`: Forward-Time Central-Space (explicit Euler) finite difference solver. This is the standard explicit time-stepping method for phase-field models See `solverft.h`. 
- `SolverSP`: Semi-implicit Fourier Spectral solver (requires FFTW). This solver uses spectral methods for high accuracy and efficiency on periodic domains. See `solversp.h` and `spectrallib.h`.

New solvers can also be defined, which is explained in [Defining Your Own Solver](#defining-your-own-solver).


## Symbols and Functions Reference

The SymPhas DSL for defining models includes grammar for using symbolic algebra. This section documents all the key symbols, constants, operators, and functions you can use in model definitions, especially for constructing evolution equations and free energy expressions.

### Coefficients and Parameters

- `c(N)` — The N-th model coefficient (1-indexed). Use for parameterizing models. Example: `c(1) * lap(op(1))`.
- `param(N)` — Alias for `c(N)`; both retrieve the N-th coefficient.
- `param_matrix(N)` — Access a matrix of coefficients for multi-field models.

### Order Parameters and Derivatives

- `op(N)` — The N-th order parameter (field variable). Example: `op(1)`.
- `dop(N)` — Time derivative of the N-th order parameter. Appears on the left-hand side of evolution equations.
- `var(N)` — The N-th provisional variable (for intermediate expressions). Must be defined using `PROVISIONAL_DEF` with assignment operator: `var(1) <= expression`.

**Note on Provisional Variable Persistence**: Provisional variables can be saved to disk alongside phase field data, making them available for post-processing and analysis. See the [Provisional Variable Saving](#provisional-variable-saving) section for details.

### Supported Field Types

- `SCALAR`, `COMPLEX`, `VECTOR` — Macros for real, complex, and vector fields.
- `SCALARS(N)`, `COMPLEXES(N)`, `VECTORS(N)` — Define N fields of the given type.

### Constants and Literals

- `1_n`, `2_n`, `0.5_n`, etc. — Compile-time constant integer and real literals (user-defined literal syntax).
- `pi_n` — The mathematical constant π.
- `e_n` — The mathematical constant e.
- `Ii` — The imaginary unit (√-1) for complex-valued models.
- `lit(x)` — Create a literal constant with value `x`.

### Time Symbol

- `TIME` — The current simulation time (symbolic variable for time-dependent models).

### Mathematical Operators and Functions

- `+`, `-`, `*`, `/` — Standard arithmetic operators.
- `power(E, N)` — Raise expression `E` to the power `N`.
- `abs(x)` — Absolute value.
- `conj(x)` — Complex conjugate.
- `modulus(x)` — Modulus (magnitude) of a complex number.
- `Re(x)`, `Im(x)` — Real and imaginary parts.
- `sin(x)`, `cos(x)`, `tan(x)`, `exp(x)`, `sqrt(x)`, `log(x)`, etc. — Standard math functions (work for real and complex arguments).

### Differential Operators

- `grad(x)` — Gradient of `x`.
- `div(x)` — Divergence of `x`.
- `curl(x)` — Curl of `x`.
- `lap(x)` — Laplacian of `x`.
- `bilap(x)` — Biharmonic (bi-Laplacian) of `x`.
- `diff(N)` — N-th order derivative representing the dot product between N gradient operators (e.g., `diff(3)` for gradient of a laplacian). This term is multiplied with an expression to take its derivative.
- `symDiff(E, VAR, N)` — Symbolic `N`-th derivative of `E` with respect to variable `VAR`.

### Vector and Tensor Operations

- `dot(a, b)` — Dot product.
- `cross(a, b)` — Cross product.

### Summation and Indexing

- `SUM(ii)(...)` — Sum over index `ii`. Most commonly used for specifying phase-field models in terms of their free energy.
- `ALL_NONCONSERVED(ii)` — Specify that all fields indexed by `ii` (in conjunction with `SUM(ii)`) use nonconserved dynamics. Can also be used with literal field indices like `ALL_NONCONSERVED(1)`.
- `ALL_CONSERVED(ii)` — Specify that all fields indexed by `ii` (in conjunction with `SUM(ii)`) use conserved dynamics. Can also be used with literal field indices like `ALL_CONSERVED(1)`.
- `EQUATION_OF(ii)(...)` — Apply dynamics specification to fields indexed by `ii`. Can be used with both symbolic indices and specific field numbers. When both are present, integer indices take precedence.

#### Using ALL_CONSERVED and ALL_NONCONSERVED with SUM

The `ALL_CONSERVED(ii)` and `ALL_NONCONSERVED(ii)` macros are designed to work with `SUM(ii)` to specify dynamics for all fields in a summation. **`SUM(ii)` is most commonly used when specifying phase-field models in terms of their free energy**, allowing you to define models with many fields using compact notation.

**Example Usage:**
```cpp
// Model with 4 scalar fields, all with conserved dynamics
MODEL(MBM, (SCALARS(4)),
    FREE_ENERGY((ALL_CONSERVED(ii)),
        INT(SUM(ii)(LANDAU_FE(op_ii))))
)

// This is equivalent to manually specifying:
// FREE_ENERGY((CONSERVED, CONSERVED, CONSERVED, CONSERVED), ...)
```

**How it works:**
- `ALL_CONSERVED(ii)` uses the symbolic index `ii` to specify that every field indexed by `ii` should use conserved dynamics
- When used with `SUM(ii)`, this applies to all fields that participate in the summation
- `op_ii` is shorthand for `op_(ii, 0)` and represents the field at index `ii` within the sum

**Using `EQUATION_OF` with indexing:**
`EQUATION_OF` can also be used with both symbolic indices and specific field numbers:

```cpp
// Apply custom dynamics to all fields indexed by ii
FREE_ENERGY((EQUATION_OF(ii)(-lap(DF(ii)))), ...)

// Apply custom dynamics to field 2 specifically
FREE_ENERGY((EQUATION_OF(2)(-bilap(DF(2)))), ...)

// Mixed usage - integer index takes precedence over symbolic index
FREE_ENERGY((EQUATION_OF(ii)(-lap(DF(ii))), EQUATION_OF(2)(-bilap(DF(2)))), ...)
// Field 2 will use the bilap equation, all other fields use the lap equation
```

**Alternative with specific field indices:**
```cpp
// Apply conserved dynamics to field 1 specifically
FREE_ENERGY((ALL_CONSERVED(1)), ...)

// Apply nonconserved dynamics to field 2 specifically  
FREE_ENERGY((ALL_NONCONSERVED(2)), ...)
```

### Noise and Randomness

- `NOISE(VARIETY, TYPE, ...)` — Add noise of a given type (e.g., `WHITE_NOISE(SCALAR)`).

### Special Macros and Utilities

- `LANDAU_FE(...)`, `DOUBLE_WELL_FE(...)` — Predefined free energy functionals.
- `e(...)` — Unit vector in a given direction, which takes the rotation as an argument. Providing `e(x)` (or equivalently using `y` or `z`) is a shortcut for rotation in the corresponding direction.
- `e_(...)` — Unit row vector in a given direction.
- `ARRAY(ii)(...)` — Indexed array access.
- `INT(E)` — Domain integral of expression `E`.
- `T(E)` — Transpose of `E`.

---

## Creating Your Own Models

SymPhas provides multiple ways to define phase-field models using the DSL This section covers all the different approaches with detailed examples and explanations.

### Model Definition Approaches

There are several ways to define models in SymPhas, each with different advantages:

1. **Direct Evolution Equations** (`EVOLUTION`) - Define time derivatives directly
2. **Free Energy Formulation** (`FREE_ENERGY`) - Define models in terms of free energy functionals
3. **Preamble Evolution** (`EVOLUTION_PREAMBLE`) - Evolution with variable definitions
4. **Specialized Models** - Advanced models with custom behavior

### 1. Direct Evolution Equations

The most straightforward approach is to define the time evolution equations directly using the `EVOLUTION` macro.

#### Single Field Models

```cpp
// Allen-Cahn equation (Model A)
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

// Complex Ginzburg-Landau equation
MODEL(ComplexGL, (COMPLEX),
    EVOLUTION(dop(1) = (1_n + c(1) * Ii) * op(1) + 
                      (1_n + c(2) * Ii) * lap(op(1)) - 
                      power(modulus(op(1)), 2) * op(1))
)

// Vector Allen-Cahn (for multi-component order parameters)
MODEL(VectorAC, (VECTOR),
    EVOLUTION(dop(1) = lap(op(1)) - c(1) * (dot(op(1), op(1)) - 1_n) * op(1))
)
```

#### Multi-Field Models

```cpp
// Model C: Cahn-Hilliard + Allen-Cahn coupling
MODEL(ModelC, (SCALARS(2)),
    EVOLUTION(
        dop(1) = -bilap(op(1)) - lap((c(1) - c(2)*power(op(1), 2))*op(1)) 
                 - c(3) * lap(op(2)),
        dop(2) = lap(op(2)) + (c(4) - c(5)*power(op(2), 2))*op(2) 
                 + c(6) * op(1)
    )
)

// Mixed field types: scalar + vector + complex
MODEL(MultiPhysics, (SCALAR, VECTOR, COMPLEX),
    EVOLUTION(
        dop(1) = lap(op(1)) - div(op(1) * op(2)),     // Scalar transport
        dop(2) = c(1) * lap(op(2)) - grad(op(1))      // Vector momentum
                 + c(2) * Re(op(3)) * op(2),
        dop(3) = c(3) * lap(op(3)) + c(4) * op(3)     // Complex field
                 - power(modulus(op(3)), 2) * op(3)
    )
)
```

### 2. Free Energy Formulation

Define models by specifying their free energy functional. SymPhas automatically derives the evolution equations using variational principles.

#### Basic Free Energy Models

```cpp
// Allen-Cahn from free energy (Model A)
MODEL(ModelA_FE, (SCALAR),
    FREE_ENERGY((NONCONSERVED), 
        INT(c(1)/2.0_n * dot(grad(op(1)), grad(op(1))) + 
            c(2)/4.0_n * power(op(1), 4) - c(3)/2.0_n * power(op(1), 2)))
)

// Cahn-Hilliard from free energy (Model B)
MODEL(ModelB_FE, (SCALAR),
    FREE_ENERGY((CONSERVED),
        INT(c(1)/2.0_n * dot(grad(op(1)), grad(op(1))) + 
            c(2)/4.0_n * power(op(1), 4) - c(3)/2.0_n * power(op(1), 2)))
)

// Using predefined functionals
MODEL(LandauModel, (SCALAR),
    FREE_ENERGY((NONCONSERVED), INT(LANDAU_FE(op(1))))
)

MODEL(DoubleWellModel, (SCALAR),
    FREE_ENERGY((CONSERVED), INT(DOUBLE_WELL_FE(op(1))))
)
```

#### Multi-Field Free Energy with Summation

```cpp
// Many-field model using summation notation
MODEL(ManyFieldModel, (SCALARS(5)),
    FREE_ENERGY((ALL_CONSERVED(ii)),
        INT(SUM(ii)(c(1)/2.0_n * dot(grad(op_ii), grad(op_ii)) + 
                    c(2)/4.0_n * power(op_ii, 4) - c(3)/2.0_n * power(op_ii, 2))))
)

// Mixed dynamics: some conserved, some non-conserved
MODEL(MixedDynamics, (SCALARS(3)),
    FREE_ENERGY((CONSERVED, NONCONSERVED, CONSERVED),
        INT(c(1)/2.0_n * dot(grad(op(1)), grad(op(1))) +
            c(2)/2.0_n * dot(grad(op(2)), grad(op(2))) +
            c(3)/2.0_n * dot(grad(op(3)), grad(op(3))) +
            c(4) * op(1) * op(2) * op(3)))
)
```

#### Custom Dynamics with `EQUATION_OF`

```cpp
// Custom evolution equations for specific fields
MODEL(CustomDynamics, (SCALARS(3)),
    FREE_ENERGY((EQUATION_OF(1)(-c(1) * bilap(DF(1))),  // Custom for field 1
                 NONCONSERVED,                            // Standard for field 2
                 EQUATION_OF(3)(-c(2) * lap(DF(3)) + c(3) * DF(3))), // Custom for field 3
        INT(SUM(ii)(LANDAU_FE(op_ii))))
)

// Using symbolic indices with custom dynamics
MODEL(SymbolicCustom, (SCALARS(4)),
    FREE_ENERGY((EQUATION_OF(ii)(-c(1) * lap(DF(ii)) + c(2) * bilap(DF(ii))),
                 EQUATION_OF(2)(-c(3) * bilap(DF(2)))),  // Override for field 2
        INT(SUM(ii)(LANDAU_FE(op_ii))))
)
```

### 3. Evolution with Preamble

Use `EVOLUTION_PREAMBLE` to define local variables and more complex expressions before the evolution equations.

```cpp
// Complex Ginzburg-Landau with preamble
MODEL(CGL_Preamble, (COMPLEX),
    EVOLUTION_PREAMBLE(
        (auto alpha = c(1);
         auto beta = c(2);
         auto gamma = c(3);
         auto nonlinear_term = power(modulus(op(1)), 2) * op(1);),
        dop(1) = alpha * lap(op(1)) + beta * op(1) - gamma * nonlinear_term
    )
)

// Multi-field with coupling calculations
MODEL(CoupledFields, (SCALARS(2)),
    EVOLUTION_PREAMBLE(
        (auto D1 = c(1);
         auto D2 = c(2);
         auto coupling_strength = c(3);
         auto field_coupling = coupling_strength * op(1) * op(2);
         auto nonlinear1 = power(op(1), 3) - op(1);
         auto nonlinear2 = power(op(2), 3) - op(2);),
        dop(1) = D1 * lap(op(1)) - nonlinear1 + field_coupling,
        dop(2) = D2 * lap(op(2)) - nonlinear2 - field_coupling
    )
)
```

### 4. Provisional Variables

Provisional variables are intermediate calculations that are computed at each time step and can be used in evolution equations (when they appear in phase-field equations, they are not substituted as expressions). Unlike preamble variables which are local C++ variables, provisional variables are full grid-based fields that persist throughout the simulation and can be saved to disk for analysis.

#### What are Provisional Variables?

Provisional variables serve as "virtual fields" - they store the results of intermediate expressions (like gradients, laplacians, or complex combinations) that are:
- **Computed every time step** before the evolution equations
- **Available throughout the model** for use in multiple equations  
- **Grid-based fields** with the same dimensions as phase fields
- **Persisted to disk** for post-processing and analysis

#### Defining Provisional Variables

```cpp
MODEL(WithProvisional, (SCALAR),
    PROVISIONAL_DEF((VECTOR),
        var(1) <= grad(op(1))
    );
    EVOLUTION(dop(1) = -lap(div(grad(op(1)))) + c(1) * op(1))
)
```

**Key Points:**
- Use `PROVISIONAL_DEF((TYPE1, TYPE2, ...), ...)` to specify provisional variable types
- Define each provisional variable with `var(N) <= expression` (note the `<=` operator which is used to do the assignment)
- Provisional variables are 1-indexed: `var(1)`, `var(2)`, etc.
- Types can be `SCALAR`, `VECTOR`, or `COMPLEX`

> **⚠️ Note: Provisional Field Evaluation**
Provisional fields are **evaluated before the solver step**, not substituted as expressions into phase-field equations. This means:
> - If you compute `var(1) <= grad(op(1))`, then `var(1)` contains the **numerical values** of the gradient
> - Using `var(1)` in equations uses these **pre-computed values**, not the symbolic expression
> - **Never apply the same operation twice**: If `var(1) = grad(op(1))`, don't use `grad(var(1))` in equations - this applies the gradient stencil twice, causing massive numerical instability

#### Complete Working Example

```cpp
// Advanced provisional model with energy density calculation
MODEL(WithProvisional, (SCALARS(2)),
    PROVISIONAL_DEF((VECTOR, SCALAR), 
        var(1) <= grad(op(1)),                    // Pre-compute gradient
        var(2) <= dot(grad(op(1)), grad(op(1)))             // Pre-compute gradient magnitude squared
    );
    EVOLUTION_PREAMBLE(
        (auto energy_density = 
            -c(1) * dot(grad(op(1)), grad(op(1))) 
            + c(2) * power(op(1), 2) 
            + c(2) * power(op(1), 4);),
        dop(1) = -bilap(op(1)) + lap(symDiff(energy_density, op(1), 1)),
        dop(2) = lap(op(2)) + c(3) * op(1)
    )
)
LINK_WITH_NAME(WithProvisional, WITH_PROVISIONAL)
```

**Key Points:**
- Use `PROVISIONAL_DEF((TYPE1, TYPE2, ...), ...)` to specify provisional variable types
- Define each provisional variable with `var(N) <= expression` (note the `<=` operator)
- Provisional variables are 1-indexed: `var(1)`, `var(2)`, etc.
- Types can be `SCALAR`, `VECTOR`, or `COMPLEX`

#### Multi-Variable Provisional Examples

```cpp
// Multiple provisional variables with different types
MODEL(MultiProvisional, (SCALARS(2)),
    PROVISIONAL_DEF((VECTOR, SCALAR, SCALAR),
        var(1) <= grad(op(1)),                    // Gradient (SCALAR → VECTOR)
        var(2) <= div(grad(op(1))),                    // Divergence of gradient (Laplacian)
        var(3) <= power(op(2), 2) - c(1)          // Nonlinear expression
    );
    EVOLUTION(
        dop(1) = -var(2) + c(2) * var(3),        // Use pre-computed Laplacian
        dop(2) = lap(op(2)) + c(3) * var(3)
    )
)

// Provisional variables can depend on each other (sequential evaluation)
MODEL(ChainedProvisional, (SCALAR),
    PROVISIONAL_DEF((VECTOR, SCALAR),
        var(1) <= grad(op(1)),                    // Computed first
        var(2) <= dot(var(1), var(1))             // Uses result of var(1)
    );
    EVOLUTION(dop(1) = lap(op(1)) - var(2) * op(1))
)
```

**⚠️ Numerical Stability Warning:**
```cpp
// ❌ WRONG: This applies gradient twice!
MODEL(BadExample, (SCALAR),
    PROVISIONAL_DEF((VECTOR),
        var(1) <= grad(op(1))     // First gradient application
    );
    EVOLUTION(dop(1) = div(var(1)))  // Second differentiation - UNSTABLE!
)

// ✅ CORRECT: Use the pre-computed values directly
MODEL(GoodExample, (SCALAR),
    PROVISIONAL_DEF((SCALAR),
        var(1) <= div(grad(op(1)))     // Pre-compute div of gradient
    );
    EVOLUTION(dop(1) = var(1))   // Use values, apply different operation
)
```

#### Provisional Variable Persistence and Analysis

**New Feature**: Provisional variables can now be saved to disk alongside phase field data, making them invaluable for:

- **Post-processing**: Analyze intermediate quantities like gradients and energy densities
- **Debugging**: Verify that complex expressions are computed correctly
- **Visualization**: Create animations of how derived quantities evolve
- **Research**: Study the behavior of intermediate expressions over time

#### Saving Provisional Variables


Provisional variables are automatically saved when the system is saved with `model.save_systems`.

```cpp
// Save provisional variables only
model.save_systems("./output");

// Automatic saving during simulation
symphas::find_solution(model, dt_list, save_params, save_points, "./output");
// ↑ Automatically saves provisional variables at each save interval
```

**File Format and Naming Convention:**

Provisional fields are saved to the `provisional` directory at the root output path as `provisional/data<I>_<N>.txt`, where `I` is the provisional field index, and `N` is the solution iteration. The file format is the same as the phase-field output format and is split into individual files per saved iteration. 

#### Practical Applications

```cpp
// Example 1: Phase field crystal model with energy density tracking
MODEL(PFCWithEnergy, (SCALAR),
    PROVISIONAL_DEF((SCALAR),
        var(1) <= power(op(1), 2) + c(1) * power(lap(op(1)), 2)  // Local energy density
    )
    EVOLUTION(dop(1) = -lap(lap(op(1)) + power(op(1), 3) - c(2) * op(1)))
)

// Example 2: Optimized gradient computation for large systems
MODEL(OptimizedGradientModel, (SCALARS(2)),
    PROVISIONAL_DEF((VECTOR, SCALAR),
        var(1) <= grad(op(1)),                    // Compute gradient once
        var(2) <= dot(var(1), var(1))             // Gradient magnitude squared
    )
    EVOLUTION(
        dop(1) = c(1) * div(var(1)) - c(2) * var(2) * op(1),  // Use pre-computed values
        dop(2) = lap(op(2)) + c(3) * var(2)                   // Reuse gradient magnitude
    )
)

// Example 3: Complex field decomposition
MODEL(ComplexFieldAnalysis, (COMPLEX),
    PROVISIONAL_DEF((SCALAR, SCALAR, VECTOR),
        var(1) <= modulus(op(1)),                 // Amplitude
        var(2) <= Re(op(1) * conj(grad(op(1)))), // Phase gradient
        var(3) <= grad(var(1))                   // Amplitude gradient
    )
    EVOLUTION(dop(1) = c(1) * lap(op(1)) - var(1) * (var(1) - c(2)) * op(1))
)

// Example 4: Multi-scale coupling (avoid repeated expensive computations)
MODEL(MultiScaleCoupling, (SCALARS(3)),
    PROVISIONAL_DEF((SCALAR, VECTOR, SCALAR),
        var(1) <= c(1) * op(1) + c(2) * op(2) + c(3) * op(3), // Coupling term
        var(2) <= grad(var(1)),                                // Coupling gradient
        var(3) <= div(var(2))                                  // Coupling divergence
    )
    EVOLUTION(
        dop(1) = lap(op(1)) - var(3) + c(4) * var(1),
        dop(2) = lap(op(2)) + c(5) * var(1),
        dop(3) = lap(op(3)) - c(6) * var(1)
    )
)
```

**Best Practices:**
- **Pre-compute expensive operations**: Gradients, Laplacians, nonlinear combinations that are used in multiple phase-field equations.
- **Reuse across multiple equations**: Compute once, use in multiple evolution equations
- **Avoid double derivatives**: Never apply the same differential operator to provisional variables
- **Sequential dependencies**: Provisional variables can depend on previously defined ones

#### Provisional vs. Preamble Variables

| Feature | Provisional Variables | Preamble Variables |
|---------|----------------------|-------------------|
| **Storage** | Full grid fields | Local C++ variables |
| **Evaluation** | Pre-computed values before solver | Computed as expressions during equation building |
| **Persistence** | Saved to disk | Temporary only |
| **Scope** | Available throughout model | Only in equations block |
| **Performance** | Computed once per time step | Computed per equation evaluation |
| **Memory** | Uses grid memory | Integrated into memory used by equation |
| **Derivatives** | ⚠️ Never apply same operation twice | Safe to use in any expression |
| **Use Case** | Complex field operations, expensive computations | Making your equations more human-readable |
| **Numerical Stability** | Requires careful handling | No special considerations |

**When to Use Provisional Variables:**
- Computing gradients, Laplacians, or other differential operators
- Expensive nonlinear combinations used in multiple equations
- Field quantities you want to analyze/visualize
- Avoiding repeated computations in multi-field systems

**When to Use Preamble Variables:**
- Simple coefficient combinations
- Local mathematical expressions
- Temporary values that don't need persistence
- Cases where you need full symbolic manipulation

### 5. Model Linking and Naming

The `LINK_WITH_NAME` construct is a powerful feature that enables **configuration-driven model selection**. It allows you to create alternative names for models that can be automatically selected at runtime based on configuration files, making your simulations more flexible and user-friendly.

#### How LINK_WITH_NAME Works

The `LINK_WITH_NAME` mechanism enables **configuration-driven model selection** through a compile-time dispatch system:

1. **Model Registration**: When you use `LINK_WITH_NAME(OriginalModel, LINKED_NAME)`, SymPhas registers the string name with the model type
2. **Compile-time Dispatch Table**: Creates template specializations that map string names to model types
3. **Runtime String Matching**: Configuration systems can reference models by their linked names
4. **Automatic Model Selection**: The framework automatically selects and instantiates the correct model type based on the string name

This enables you to change models by editing configuration files without recompiling your code.

#### Basic Model Linking

```cpp
// Define a model
MODEL(MyComplexModel, (SCALARS(2)),
    EVOLUTION(
        dop(1) = lap(op(1)) - power(op(1), 3) + c(1) * op(2),
        dop(2) = c(2) * lap(op(2)) + c(3) * op(1) * op(2)
    )
)

// Link with a simpler name for configuration files
LINK_WITH_NAME(MyComplexModel, SIMPLE_MODEL)

// Now you can use either name in C++ code:
// model_MyComplexModel_t<2, SolverFT<Stencil2d2h<>>> model1{params};
// model_SIMPLE_MODEL_t<2, SolverFT<Stencil2d2h<>>> model2{params};
```

#### Multiple Aliases for the Same Model

```cpp
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

LINK_WITH_NAME(AllenCahn, AC)
LINK_WITH_NAME(AllenCahn, MODEL_A)
LINK_WITH_NAME(AllenCahn, PHASE_SEPARATION)

// All of these are equivalent:
// model_AllenCahn_t<2, SolverFT<Stencil2d2h<>>> model1{params};
// model_AC_t<2, SolverFT<Stencil2d2h<>>> model2{params};
// model_MODEL_A_t<2, SolverFT<Stencil2d2h<>>> model3{params};
// model_PHASE_SEPARATION_t<2, SolverFT<Stencil2d2h<>>> model4{params};
```

#### Using Linked Names in Configuration Files

**JSON Configuration:**
The linked model names can be used directly in JSON configuration files:

```json
{
  "model": {
    "name": "SIMPLE_MODEL",     // Uses the linked name
    "domain": "main_domain",
    "coefficients": {
      "c1": 0.5,
      "c2": 1.0,
      "c3": -0.3
    }
  }
}
```

#### Runtime Model Selection

**Method 1: Using model_select (Recommended)**
```cpp
#include "symphas.h"

// Define a simulation application type
template<typename ModelType>
struct Simulation {
    int simulate(const double* coeff, size_t num_coeff) {
        // Create model instance with coefficients
        symphas::problem_parameters_type params = setup_parameters();
        ModelType model{coeff, num_coeff, params};
        
        // Run simulation
        symphas::find_solution(model, 0.01, 1000);
        return 0;
    }
    
    template<typename... Ts>
    auto operator()(Ts&&... args) {
        return simulate(std::forward<Ts>(args)...);
    }
};

int main() {
    // Create model selector with dimension and stencil parameters
    model_select<Simulation> selector{2, StencilParams{2, 9, 6, 13}};
    
    // Select and run model by string name
    double coeffs[] = {1.0, -1.0, 1.0};
    if (selector.call<SolverFT>("SIMPLE_MODEL", coeffs, 3) == INVALID_MODEL) {
        fprintf(stderr, "Unknown model: SIMPLE_MODEL\n");
        return 1;
    }
    
    return 0;
}
```

**Method 2: Direct Model Instantiation**
```cpp
#include "symphas.h"

int main() {
    // Create problem parameters
    symphas::problem_parameters_type params = setup_parameters();
    
    // Manually instantiate using linked name (both names work identically)
    model_SIMPLE_MODEL_t<2, SolverFT<Stencil2d2h<>>> model{params};
    // This is equivalent to:
    // model_MyComplexModel_t<2, SolverFT<Stencil2d2h<>>> model{params};
    
    // Run simulation
    symphas::find_solution(model, 0.01, 1000);
    return 0;
}
```

**Method 3: Configuration-Driven with model_select**
```cpp
#include "symphas.h"

template<typename ModelType>
struct ConfigDrivenSimulation {
    int simulate() {
        // Access the global configuration (already loaded by symphas::init)
        auto& config = symphas::conf::config();
        
        // Create and run model with config parameters
        ModelType model{config.model_settings.coeff, 
                       config.model_settings.coeff_len, 
                       config.get_problem_parameters()};
        
        // Run simulation with parameters from config
        symphas::find_solution(model, config.simulation_settings.get_time_step_list());
        return 0;
    }
    
    template<typename... Ts>
    auto operator()(Ts&&... args) { return simulate(); }
};

int main() {
    // Initialize SymPhas with configuration file
    symphas::init("simulation.json", nullptr, 0);
    
    // Extract model name from the loaded configuration
    const char* model_name = symphas::conf::config().model_settings.model;
    
    // Create model selector with parameters from config
    auto& config = symphas::conf::config();
    model_select<ConfigDrivenSimulation> selector{
        config.simulation_settings.dimension,
        config.simulation_settings.stp
    };
    
    // Select and run the model by name from configuration
    if (selector.call<SolverFT>(model_name) == INVALID_MODEL) {
        fprintf(stderr, "Error: Unknown model '%s'\n", model_name);
        return 1;
    }
    
    return 0;
}
```

#### Benefits of Model Linking

1. **Configuration-Driven Simulation**: Change models without recompiling by editing JSON/config files
2. **User-Friendly Names**: Use intuitive names instead of complex C++ identifiers
3. **Backward Compatibility**: Support both technical and simplified naming schemes
4. **Runtime Flexibility**: Same executable can run different models based on configuration
5. **Batch Processing**: Easy to create parameter studies with different models

#### Best Practices for Model Linking

**1. Naming Conventions:**
```cpp
// Use descriptive, configuration-friendly names
LINK_WITH_NAME(AllenCahnConserved, ALLEN_CAHN)          // Good: Clear physics meaning
LINK_WITH_NAME(AllenCahnConserved, AC_CONSERVED)       // Good: Standard abbreviation
LINK_WITH_NAME(AllenCahnConserved, model_type_47)      // Poor: No semantic meaning
LINK_WITH_NAME(AllenCahnConserved, MyCrazyModel)       // Poor: Unprofessional
```

**2. Organization Strategies:**
```cpp
// Group related models with consistent prefixes
LINK_WITH_NAME(PhaseField2D, PF_2D)
LINK_WITH_NAME(PhaseField3D, PF_3D)
LINK_WITH_NAME(PhaseFieldNoise, PF_NOISE)

// Alternative physics-based grouping
LINK_WITH_NAME(AllenCahn, SPINODAL_DECOMPOSITION)
LINK_WITH_NAME(CahnHilliard, PHASE_SEPARATION)
LINK_WITH_NAME(SwiftHohenberg, PATTERN_FORMATION)
```

**3. Avoiding Common Pitfalls:**
```cpp
LINK_WITH_NAME(ModelA, SHARED_NAME)
LINK_WITH_NAME(ModelB, SHARED_NAME)  // Compile error!

// ❌ DON'T: Use names that conflict with C++ keywords or SymPhas internals
LINK_WITH_NAME(MyModel, CLASS)      // C++ keyword
LINK_WITH_NAME(MyModel, MODEL)      // SymPhas reserved
LINK_WITH_NAME(MyModel, STENCIL)    // SymPhas reserved

// ✅ DO: Use unique, descriptive names
LINK_WITH_NAME(ModelA, MODEL_A_VARIANT)
LINK_WITH_NAME(ModelB, MODEL_B_VARIANT)
```


**4. Testing Model Links:**
```cpp
// Always test that your linked names work in configuration files
{
  "model": {
    "name": "BINARY_MIXTURE",  // Test each linked name
    "coefficients": {"c1": -1.0, "c2": 1.0}
  }
}
```

#### Technical Summary: LINK_WITH_NAME Implementation

The `LINK_WITH_NAME` system is a compile-time metaprogramming framework that enables string-based model selection while maintaining full type safety and performance.

**Key Benefits:**
- **Configuration-Driven Simulation**: Change models without recompiling by editing JSON/config files
- **User-Friendly Names**: Use intuitive names instead of complex C++ identifiers  
- **Performance**: Zero runtime overhead after model selection
- **Flexibility**: Same executable can run different models based on configuration


### 6. Advanced Model Features

#### Time-Dependent Parameters

```cpp
MODEL(TimeDependent, (SCALAR),
    EVOLUTION_PREAMBLE(
        (auto time_factor = sin(c(1) * t);
         auto mobility = c(2) * (1_n + c(3) * time_factor);),
        dop(1) = mobility * lap(op(1)) - (power(op(1), 2) - 1_n) * op(1)
    )
)
```

#### Noise Integration

```cpp
MODEL(StochasticAC, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1) + 
                      c(1) * WHITE_NOISE(SCALAR))
)

// Multi-field with white noise
MODEL(StochasticMultiField, (SCALARS(2)),
    EVOLUTION(
        dop(1) = lap(op(1)) - power(op(1), 3) + c(1) * WHITE_NOISE(SCALAR),
        dop(2) = c(2) * lap(op(2)) + c(3) * op(1) + c(4) * WHITE_NOISE(SCALAR)
    )
)
```

### 7. Model Organization and Best Practices

#### File Organization

```cpp
// In your model header file (e.g., mymodels.h)
#pragma once
#include "modelmacros.h"

// Group related models together
namespace MyModels {
    
    // Basic phase separation
    MODEL(PhaseSeparation, (SCALAR),
        FREE_ENERGY((CONSERVED), INT(DOUBLE_WELL_FE(op(1))))
    )
    LINK_WITH_NAME(PhaseSeparation, PS)
    
    // With noise
    MODEL(PhaseSeparationNoise, (SCALAR),
        EVOLUTION(dop(1) = -lap(symDiff(DOUBLE_WELL_FE(op(1)), op(1))) + 
                          c(1) * WHITE_NOISE(SCALAR))
    )
    LINK_WITH_NAME(PhaseSeparationNoise, PSN)
    
    // Multi-component
    MODEL(MultiComponentPS, (SCALARS(3)),
        FREE_ENERGY((ALL_CONSERVED(ii)),
            INT(SUM(ii)(DOUBLE_WELL_FE(op_ii)) + 
                c(1) * op(1) * op(2) + c(2) * op(2) * op(3)))
    )
    LINK_WITH_NAME(MultiComponentPS, MCPS)
}
```

#### Parameter Conventions

```cpp
// Document your parameter meanings
MODEL(DocumentedModel, (SCALARS(2)),
    // Parameters:
    // c(1): Diffusion coefficient for field 1
    // c(2): Diffusion coefficient for field 2  
    // c(3): Coupling strength
    // c(4): Nonlinearity strength field 1
    // c(5): Nonlinearity strength field 2
    EVOLUTION(
        dop(1) = c(1) * lap(op(1)) - c(4) * power(op(1), 3) + c(3) * op(2),
        dop(2) = c(2) * lap(op(2)) - c(5) * power(op(2), 3) - c(3) * op(1)
    )
)
```


#### Complete Example: Model Library with Configuration

**Model definitions (models.h):**
```cpp
#pragma once
#include "symphas.h"

// Allen-Cahn model
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)
LINK_WITH_NAME(AllenCahn, ALLEN_CAHN)
LINK_WITH_NAME(AllenCahn, MODEL_A)

// Cahn-Hilliard model  
MODEL(CahnHilliard, (SCALAR),
    EVOLUTION(dop(1) = -lap(lap(op(1)) - power(op(1), 3) + c(1) * op(1)))
)
LINK_WITH_NAME(CahnHilliard, CAHN_HILLIARD)
LINK_WITH_NAME(CahnHilliard, MODEL_B)

// Swift-Hohenberg model
MODEL(SwiftHohenberg, (SCALAR),
    EVOLUTION(dop(1) = c(1) * op(1) - power(1_n + lap(op(1)), 2) * op(1) - power(op(1), 3))
)
LINK_WITH_NAME(SwiftHohenberg, SWIFT_HOHENBERG)
LINK_WITH_NAME(SwiftHohenberg, SH_MODEL)

// Advanced model with provisional variables
MODEL(WithProvisional, (SCALARS(2)),
    PROVISIONAL_DEF((VECTOR, SCALAR), 
        var(1) <= grad(op(1)),                    // Pre-compute gradient
        var(2) <= dot(var(1), var(1))             // Pre-compute gradient magnitude squared
    )
    EVOLUTION_PREAMBLE(
        (auto energy_density = -c(1) * dot(grad(op(1)), grad(op(1))) +
                               c(2) * power(op(1), 2) +
                               c(2) * power(op(1), 4);),
        dop(1) = -bilap(op(1)) + lap(symDiff(energy_density, op(1), 1)),
        dop(2) = lap(op(2)) + c(3) * op(1)
    )
)
LINK_WITH_NAME(WithProvisional, WITH_PROVISIONAL)
```

**Configuration file (simulation.json):**
```json
{
  "definitions": {
    "GRID_SIZE": 128,
    "SAVE_INTERVAL": 100
  },
  "simulation": {
    "dimension": 2,
    "time_steps": [{"dt": 0.01, "iterations": 10000}],
    "save": {"type": "DEFAULT", "interval": "${SAVE_INTERVAL}"}
  },
  "model": {
    "name": "WITH_PROVISIONAL",    // Uses model with provisional variables
    "domain": "main_domain",
    "coefficients": {
      "c1": 0.5,
      "c2": 1.0,
      "c3": -0.3
    }
  },
  "domains": {
    "main_domain": {
      "dimension": 2,
      "intervals": {
        "x": {"start": 0.0, "end": 1.0, "points": "${GRID_SIZE}"},
        "y": {"start": 0.0, "end": 1.0, "points": "${GRID_SIZE}"}
      },
      "boundaries": {
        "x": {"left": {"type": "PERIODIC"}, "right": {"type": "PERIODIC"}},
        "y": {"left": {"type": "PERIODIC"}, "right": {"type": "PERIODIC"}}
      },
      "initial_conditions": [
        {
          "type": "UNIFORM",
          "parameters": [-1, 1]
        },
        {
          "type": "UNIFORM", 
          "parameters": [-0.5, 0.5]
        }
      ]
    }
  }
}
```

**Alternative Configuration for Swift-Hohenberg:**
```json
{
  "model": {
    "name": "SWIFT_HOHENBERG",    // Switch to different model
    "domain": "main_domain",
    "coefficients": {
      "c1": 0.5,
      "c2": 1.0,
      "c3": -0.3
    }
  }
}
```

**Driver program (main.cpp):**
```cpp
#include "symphas.h"
#include "models.h"  // Include your model definitions

// Define application type for running simulations
template<typename ModelType>
struct PhaseFieldSimulation {
    int simulate() {
        // In a real application, you would:
        // 1. Access the loaded configuration via symphas::conf::config()
        // 2. Create problem parameters from config data
        // 3. Run the simulation
        
        printf("Running simulation with model: %s\n", typeid(ModelType).name());
        return 0;
    }
    
    template<typename... Ts>
    auto operator()(Ts&&... args) { 
        return simulate(); 
    }
};

int main() {
    // Model name extracted from configuration file
    const char* model_name = "SWIFT_HOHENBERG";  // From simulation.json
    
    // Create model selector with 2D simulation parameters
    model_select<PhaseFieldSimulation> selector{2, StencilParams{2, 9, 6, 13}};
    
    // Select and run the model by name
    if (selector.call<SolverFT>(model_name) == INVALID_MODEL) {
        fprintf(stderr, "Error: Unknown model '%s'\n", model_name);
        return 1;
    }
    
    printf("Simulation completed successfully\n");
    return 0;
}
```

---


## Boundary Conditions

SymPhas provides a comprehensive boundary condition system that supports various types of boundary behaviors for phase-field simulations. The system is designed to be flexible, allowing you to specify different boundary conditions for each side of your simulation domain.

### Overview of Boundary Types

The boundary condition system supports two main categories:
- **Boundary Types**: Define the fundamental boundary behavior
- **Boundary Tags**: Modify DEFAULT boundaries with specific algorithms

SymPhas supports the following boundary types:

#### 1. `DEFAULT`
- The primary boundary type for non-periodic conditions
- Supports various tags to define specific behavior (Dirichlet-like, Neumann-like, etc.)
- Most flexible boundary type with algorithmic modifications

#### 2. `PERIODIC`
- Implements true periodic (wrap-around) boundaries
- Values at opposite boundaries are enforced to be identical
- No additional parameters or tags supported

#### 3. `OPEN`
- **Currently not implemented** - throws an exception if used
- Reserved for future implementation of absorbing/outflow boundaries
- Defaults to DEFAULT behavior if encountered

### DEFAULT Boundary Type with Tags

The DEFAULT boundary type is the most versatile, supporting various algorithmic modifications through tags:

#### Available Tags

**`CONSTANT`**: Sets constant values at boundary points
- **Parameters**: 1 value (the constant)
- **Behavior**: All boundary points get the same fixed value
- **Example**: Dirichlet-like boundary with fixed temperature

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [0.0]},
      "right": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [1.0]}
    }
  }
}
```

**`LINEAR`**: Applies linear functions across boundary points
- **Parameters**: 
  - 2D: 2 values (Ax + B)
  - 3D: 3 values (Ax + By + C)
- **Behavior**: Values vary linearly along the boundary
- **Supports**: TIME modifier for time-dependent linear changes

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["LINEAR"], "parameters": [1.0, 0.0]}
    }
  }
}
```

**`GAUSSIAN`**: Creates Gaussian distributions on boundaries
- **Parameters**:
  - 2D: 3 values (A * exp(-(x-B)²/(2C²)))
  - 3D: 4 values (A * exp(-(x-B)² - (y-C)²)/(2D²)))
- **Behavior**: Gaussian peak centered on the boundary
- **No modifiers supported**

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["GAUSSIAN"], "parameters": [1.0, 0.5, 0.1]}
    }
  }
}
```

**`TRIG`**: Applies trigonometric functions to boundaries
- **Parameters**: 3 values
  - 2D: A * sin(Cx) + B * cos(Cx)
  - 3D: A * sin(C(x+y)) + B * cos(C(x+y))
- **Supports**: TIME modifier for oscillating boundaries
- **With TIME**: Adds phase term (+ Dt)

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["TRIG"], "parameters": [1.0, 0.0, 6.28]}
    }
  }
}
```

**`RANDOM`**: Generates random values on boundaries
- **Parameters**: 2 values (uniform distribution range [A, B])
- **Behavior**: Each boundary point gets a random value in the specified range
- **No modifiers supported**

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["RANDOM"], "parameters": [0.0, 1.0]}
    }
  }
}
```

**`TIME` (Modifier)**: Can be combined with LINEAR and TRIG for time-dependent boundaries
- **With LINEAR**: Adds time-dependent term
- **With TRIG**: Creates oscillating phase
- **Parameters**: Adds one additional time scaling parameter

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["LINEAR", "TIME"], "parameters": [1.0, 0.0, 0.1]}
    }
  }
}
```

### Common Use Cases and Examples

#### 1. Periodic Boundaries (Bulk Simulations)
Use when simulating infinite or bulk systems:
```json
{
  "boundaries": {
    "x": {
      "left": {"type": "PERIODIC"},
      "right": {"type": "PERIODIC"}
    },
    "y": {
      "left": {"type": "PERIODIC"},
      "right": {"type": "PERIODIC"}
    }
  }
}
```

#### 2. Fixed Value Boundaries (Dirichlet-like)
Use for fixed temperature, concentration, or field values:
```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [0.0]},
      "right": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [1.0]}
    }
  }
}
```

#### 3. Linear Gradient Boundaries
Create linear variations across boundaries:
```json
{
  "boundaries": {
    "y": {
      "left": {"type": "DEFAULT", "tags": ["LINEAR"], "parameters": [1.0, 0.0]},
      "right": {"type": "DEFAULT", "tags": ["LINEAR"], "parameters": [-1.0, 2.0]}
    }
  }
}
```

#### 4. Time-Dependent Boundaries
For oscillating or time-varying boundary conditions:
```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["TRIG", "TIME"], "parameters": [1.0, 0.0, 6.28, 0.1]}
    }
  }
}
```

#### 5. Mixed Boundary Types
Combine different boundary types for complex setups:
```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [0.0]},
      "right": {"type": "DEFAULT", "tags": ["GAUSSIAN"], "parameters": [1.0, 0.5, 0.1]}
    },
    "y": {
      "left": {"type": "PERIODIC"},
      "right": {"type": "PERIODIC"}
    }
  }
}
```

### Advanced Features

#### Multi-Dimensional Linear Functions
For 3D simulations, LINEAR boundaries can be vary across the entire plane:

```json
{
  "boundaries": {
    "z": {
      "left": {"type": "DEFAULT", "tags": ["LINEAR"], "parameters": [1.0, 0.5, 2.0]}
    }
  }
}
```
This creates: f(x,y) = 1.0*x + 0.5*y + 2.0

#### Oscillating Boundaries
Combine TRIG and TIME for dynamic boundary conditions:

```json
{
  "boundaries": {
    "x": {
      "left": {"type": "DEFAULT", "tags": ["TRIG", "TIME"], "parameters": [1.0, 0.0, 6.28, 0.1]}
    }
  }
}
```
This creates: f(x,t) = 1.0 * sin(6.28*x + 0.1*t)

#### Error Handling and Validation

SymPhas provides comprehensive boundary condition validation:

- **Type Validation**: Ensures supported boundary types
- **Parameter Count**: Validates correct number of parameters for each tag
- **Tag Compatibility**: Checks valid tag combinations
- **Dimension Consistency**: Ensures parameter count matches simulation dimension

**Common Error Messages:**
```
Error: Invalid boundary type "INVALID_TYPE"
Supported types: DEFAULT, PERIODIC, OPEN

Error: GAUSSIAN tag requires 3 parameters for 2D, got 2
Parameters: [amplitude, center, width]

Error: TIME modifier only supported with LINEAR and TRIG tags
```

---

## JSON Configuration System

SymPhas provides a modern JSON-based configuration system as an alternative to the legacy configuration format. The JSON system offers structured, readable configuration files with variable substitution and type-safe parameter parsing.

### JSON Configuration Overview

The JSON configuration system is implemented through the `JsonConfManager` class in `conf/inc/confjson.h` and `conf/src/confjson.cpp`. It provides:

- **Structured Configuration**: Organized sections for simulation, model, domain, and naming parameters
- **Variable Substitution**: Define reusable variables with `${VARIABLE_NAME}` syntax
- **Type Safety**: Automatic type checking and conversion
- **Legacy Compatibility**: Support for existing configuration formats
- **Enhanced Save Parameters**: Fine-grained control over data output
- **Multi-Phase Time Stepping**: Support for different time steps in simulation phases
- **Comprehensive Stencil Control**: Integration with the stencil override system

- 
### Complete JSON Structure

A comprehensive JSON configuration consists of these main sections:

```json
{
  "definitions": {
    "DT": 0.01,
    "GRID_SIZE": 128,
    "SIMULATION_TIME": 100.0,
    "SAVE_INTERVAL": 10.0
  },
  "simulation": {
    "dimension": 2,
    "time_steps": [
      {
        "dt": "${DT}",
        "time": "${SIMULATION_TIME}"
      }
    ],
    "save": {
      "type": "DEFAULT",
      "interval": "${SAVE_INTERVAL}",
      "save_initial": true
    },
    "stencil": {
      "order": 2,
      "ptl": 5,
      "ptg": 6,
      "ptb": 13
    },
    "checkpoint": {
      "count": 5
    },
    "solver_variation": 0
  },
  "model": {
    "name": "PHASE_FIELD_CRYSTAL_MODEL",
    "domain": "main_domain",
    "coefficients": {
      "alpha": -0.5,
      "beta": 1.0,
      "gamma": 1.0
    },
    "fields": [
      {
        "count": 1,
        "modifiers": "NONE"
      }
    ]
  },
  "intervals": {
    "unit_interval": {
      "start": 0.0,
      "end": 1.0,
      "points": "${GRID_SIZE}"
    },
    "large_domain": {
      "start": -5.0,
      "end": 5.0,
      "points": 256
    }
  },
  "domains": {
    "main_domain": {
      "dimension": 2,
      "intervals": {
        "x": "unit_interval",
        "y": "unit_interval"
      },
      "boundaries": {
        "x": {
          "left": {"type": "PERIODIC"},
          "right": {"type": "PERIODIC"}
        },
        "y": {
          "left": {"type": "PERIODIC"},
          "right": {"type": "PERIODIC"}
        }
      }
    }
  },
  "naming": {
    "dir": "./results",
    "data": "phase_field",
    "solution": "solution",
    "checkpoint": "checkpoint"
  }
}
```

### Detailed Section Documentation

#### 1. Definitions (Optional)

Define reusable variables that can be referenced throughout the configuration:

```json
{
  "definitions": {
    "DT": 0.01,
    "DIMS": [64, 64, 32],
    "SYSTEM_INTERVAL": [0.0, 1.0],
    "MATERIAL_PARAMS": {
      "diffusion": 0.5,
      "mobility": 1.0
    }
  }
}
```

Variables can be referenced using `${VARIABLE_NAME}` syntax throughout the configuration.

#### 2. Simulation Parameters

Core simulation settings with enhanced save parameter support:

```json
{
  "simulation": {
    "dimension": 2,
    "time_steps": [
      {
        "dt": 0.01,
        "time": 100.0
      }
    ],
    "save": {
      "type": "DEFAULT",
      "interval": 10.0,
      "save_initial": true
    },
    "stencil": {
      "order": 2,
      "ptl": 5,
      "ptg": 5,
      "ptb": 13
    },
    "checkpoint": {
      "count": 5
    },
    "solver_variation": 0
  }
}
```

#### Enhanced Save Parameter Types

The `save` section supports multiple save types that map directly to the SymPhas `SaveParams` structure:

**DEFAULT**: Regular interval-based saving
```json
"save": {
  "type": "DEFAULT",
  "interval": 10.0,
  "save_initial": true
}
```

**MUL**: Multiplicative/adaptive save intervals
```json
"save": {
  "type": "MUL",
  "interval": 6,
  "save_initial": false
}
```

**How MUL works**: You specify the total number of saves you want (`interval`), and SymPhas automatically calculates optimal exponential scaling. For example, with `"interval": 6` over 10,000 time steps, saves occur at: 1, 6, 39, 251, 1584, 10000.

**EXP**: Exponential save intervals based on index
```json
"save": {
  "type": "EXP", 
  "interval": 2.0,
  "save_initial": true
}
```

**How EXP works**: Save spacing = `interval`^`index`. With `"interval": 2.0`, saves occur at: step 1, 2, 4, 8, 16, 32, etc.

**LIST**: Specific list of save indices
```json
"save": {
  "type": "LIST",
  "indices": [10, 50, 100, 200, 500, 1000],
  "save_initial": false
}
```

**SPLIT**: Split save functionality
```json
"save": {
  "type": "SPLIT",
  "interval": 10,
  "save_initial": true
}
```

**How SPLIT works**: You specify the total number of saves (`interval`), and SymPhas splits them evenly across the simulation.

#### Multi-Phase Time Stepping

Support for multiple time stepping approaches with different phases:

**Single time step:**
```json
{
  "simulation": {
    "time_step": 0.01,
    "time": 100.0
  }
}
```

**Multiple time steps with different dt values:**
```json
{
  "simulation": {
    "time_steps": [
      {
        "dt": 0.001,
        "time": 10.0        // Run for 10.0 time units
      },
      {
        "dt": 0.01,
        "iterations": 5000   // Run for exactly 5000 iterations
      }
    ]
  }
}
```

**Mixed time/iterations specification:**
```json
{
  "definitions": {
    "DT_SMALL": 0.005,
    "DT_LARGE": 0.02,
    "BURN_IN_STEPS": 1000,
    "MAIN_SIM_TIME": 50.0
  },
  "simulation": {
    "time_steps": [
      {
        "dt": "${DT_SMALL}",
        "iterations": "${BURN_IN_STEPS}",
        "comment": "Fine time step for burn-in"
      },
      {
        "dt": "${DT_LARGE}", 
        "time": "${MAIN_SIM_TIME}",
        "comment": "Larger time step for main simulation"
      }
    ]
  }
}
```

#### 3. Model Configuration

Model definition that references a named domain:

```json
{
  "model": {
    "name": "PHASE_FIELD_CRYSTAL_MODEL",
    "domain": "main_domain",
    "coefficients": {
      "alpha": -0.5,
      "beta": 1.0,
      "gamma": 1.0
    },
    "fields": [
      {
        "count": 1,
        "modifiers": "NONE"
      }
    ]
  }
}
```

#### 4. Intervals

Named interval definitions that can be reused across domains:

```json
{
  "intervals": {
    "unit_interval": {
      "start": 0.0,
      "end": 1.0,
      "points": 64
    },
    "large_domain": {
      "start": -5.0,
      "end": 5.0,
      "points": 256
    },
    "fine_grid": {
      "start": 0.0,
      "end": 10.0,
      "width": 0.02    // Alternative to points
    }
  }
}
```

#### 5. Domains

Named domain definitions that reference interval objects:

```json
{
  "domains": {
    "main_domain": {
      "dimension": 2,
      "intervals": {
        "x": "unit_interval",
        "y": "unit_interval"
      },
      "boundaries": {
        "x": {
          "left": {"type": "PERIODIC"},
          "right": {"type": "PERIODIC"}
        },
        "y": {
          "left": {"type": "PERIODIC"},
          "right": {"type": "PERIODIC"}
        }
      }
    },
    "custom_domain": {
      "dimension": 2,
      "intervals": {
        "x": "large_domain",
        "y": "unit_interval"
      },
      "boundaries": {
        "x": {
          "left": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [0.0]},
          "right": {"type": "DEFAULT", "tags": ["CONSTANT"], "parameters": [0.0]}
        },
        "y": {
          "left": {"type": "PERIODIC"},
          "right": {"type": "PERIODIC"}
        }
      }
    }
  }
}
```

#### 6. Naming

Output file naming configuration:

```json
{
  "naming": {
    "dir": "./results",
    "data": "phase_field",
    "solution": "solution",
    "checkpoint": "checkpoint"
  }
}
```

### Working with Configuration Files

SymPhas provides powerful configuration capabilities through JSON files, allowing you to specify simulation parameters, initial conditions, boundary conditions, and even model selection without recompiling your code.

#### Loading Configuration from Code

To load and use JSON configuration files directly from your C++ code:

```cpp
#include "symphas.h"
#include "solverinclude.h"

// Define a simulation application that can work with any model type
template<typename ModelType>
struct ConfigDrivenSimulation {
    int simulate() {
        // Access the loaded configuration
        auto& config = symphas::conf::config();
        
        // Get model parameters from configuration
        auto problem_params = config.get_problem_parameters();
        ModelType model{problem_params};
        
        // Run simulation with parameters from config
        symphas::find_solution(model, config.simulation_settings.get_time_step_list());
        
        return 0;
    }
    
    template<typename... Ts>
    auto operator()(Ts&&... args) { return simulate(); }
};

int main() {
    // Load JSON configuration file
    symphas::init("simulation_config.json", nullptr, 0);
    
    // Get configuration settings
    auto& config = symphas::conf::config();
    const char* model_name = config.model_settings.model;
    
    // Create model selector with parameters from configuration
    model_select<ConfigDrivenSimulation> selector{
        config.simulation_settings.dimension,
        config.simulation_settings.stp
    };
    
    // Select and run the model specified in configuration
    if (selector.call<SolverFT>(model_name) == INVALID_MODEL) {
        fprintf(stderr, "Error: Unknown model '%s'\n", model_name);
        return 1;
    }
    return 0;
}
```

#### Sample Configuration File

Create a file named `simulation_config.json`:

```json
{
    "definitions": {
        "GRID_SIZE": 64,
        "SAVE_INTERVAL": 100
    },
    "simulation": {
        "dimension": 2,
        "time_steps": [{"dt": 0.01, "iterations": 10000}],
        "save": {"type": "DEFAULT", "interval": "${SAVE_INTERVAL}"}
    },
    "model": {
        "name": "SWIFT_HOHENBERG",    // Automatically selects SwiftHohenberg model
        "domain": "main_domain",
        "coefficients": {
            "c1": 0.5,
            "c2": 1.0,
            "c3": -0.3
        }
    },
    "domains": {
        "main_domain": {
            "dimension": 2,
            "intervals": {
                "x": {"start": 0.0, "end": 1.0, "points": "${GRID_SIZE}"},
                "y": {"start": 0.0, "end": 1.0, "points": "${GRID_SIZE}"}
            },
            "boundaries": {
                "x": {"left": {"type": "PERIODIC"}, "right": {"type": "PERIODIC"}},
                "y": {"left": {"type": "PERIODIC"}, "right": {"type": "PERIODIC"}}
            }
        }
    }
}
```

#### Configuration-Based Model Selection

When using configuration files, you can select different models without recompiling:

**1. Define multiple models in your code:**
```cpp
// Define multiple models that can be selected via configuration
MODEL(AllenCahn, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - (power(op(1), 2) - 1_n) * op(1))
)

MODEL(CahnHilliard, (SCALAR),
    EVOLUTION(dop(1) = -lap(lap(op(1)) - (power(op(1), 3) - op(1))))
)

MODEL(GrayScott, (SCALAR, SCALAR),
    EVOLUTION(
        dop(1) = D_u * lap(op(1)) - op(1) * power(op(2), 2) + F * (1_n - op(1)),
        dop(2) = D_v * lap(op(2)) + op(1) * power(op(2), 2) - (F + k) * op(2)
    )
)
```

**2. Use configuration-driven model selection:**
```cpp
// Define a simulation runner that works with any model type
template<typename ModelType>
struct ConfigSimulation {
    int run() {
        auto& config = symphas::conf::config();
        
        // Create model from configuration
        auto problem_params = config.get_problem_parameters();
        ModelType model{problem_params};
        
        // Run simulation with parameters from config
        symphas::find_solution(model, config.simulation_settings.get_time_step_list());
        
        return 0;
    }
    
    template<typename... Ts>
    auto operator()(Ts&&... args) { return run(); }
};

int main() {
    using namespace symphas;
    
    // Load configuration
    symphas::init("config.json", nullptr, 0);
    auto& config = symphas::conf::config();
    
    // Get model name from config
    const char* model_name = config.model_settings.model;
    
    // Create model selector with config parameters
    model_select<ConfigSimulation> selector{
        config.simulation_settings.dimension,
        config.simulation_settings.stp
    };
    
    // Run simulation - model is selected based on "model" field in JSON
    if (selector.call<SolverFT>(model_name) == INVALID_MODEL) {
        fprintf(stderr, "Unknown model: %s\n", model_name);
        return 1;
    }
    
    return 0;
}
```

**3. Switch models by changing the JSON file:**

For Allen-Cahn:
```json
{"model": "AllenCahn", ...}
```

For Cahn-Hilliard:
```json
{"model": "CahnHilliard", ...}
```

For Gray-Scott:
```json
{
    "model": "GrayScott",
    "system": [
        {"field_type": "SCALAR", ...},  // u field
        {"field_type": "SCALAR", ...}   // v field
    ],
    ...
}
```

#### Advanced Configuration Features

#### Variable Substitution

Variables defined in the `definitions` section can be used throughout the configuration:

```json
{
  "definitions": {
    "DOMAIN_SIZE": 64,
    "SIMULATION_TIME": 1000.0,
    "DIFFUSION_COEFF": 0.1,
    "REACTION_RATE": 1.0
  },
  "simulation": {
    "time": "${SIMULATION_TIME}",
    "dt": 0.001
  },
  "model": {
    "coefficients": {
      "c1": "${DIFFUSION}",
      "c3": "${MOBILITY}",
      "c7": 2.5,
      "alpha": -0.8,
      "beta": 1.0
    }
  },
  "intervals": {
    "main_interval": {
      "start": 0.0,
      "end": 1.0,
      "points": "${DOMAIN_SIZE}"
    }
  }
}
```

**Usage with C++:**
```cpp
#include "confjson.h"

// Create and parse JSON configuration
JsonConfManager config("config.json");

// Access parsed settings
auto* sim_settings = config.get_simulation_settings();
auto domain_settings = config.get_domain_settings();
auto model_settings = config.get_model_settings();

// Use in simulation setup
// ... simulation initialization using config
```

**Migration from Legacy Format:**

Legacy format:
```
DT 0.01
SAVE_RATE 10.0
DIMENSIONS 64 64
```

Equivalent JSON:
```json
{
  "simulation": {
    "time_step": 0.01,
    "save": {
      "type": "DEFAULT", 
      "interval": 10.0
    }
  },
  "intervals": {
    "main": {
      "start": 0.0,
      "end": 1.0,
      "points": 64
    }
  }
}
```

---






### Using External Coefficient Files

SymPhas supports specifying model parameters and coupling coefficients in a separate file, referenced from your JSON configuration using the `"coefficients_file"` key. This is especially useful for multi-field models or when managing large sets of parameters.


#### How to Use

1. **Create a coefficient file** in the supported format (see below).
2. **Add the `coefficients_file` key** to your JSON config, pointing to the file (relative to the JSON file's location):

```json
{
  "simulation": {
    "model": "PFC2",
    "dimension": 2,
    ...,
    "coefficients_file": "pfc2_coeffs.txt"
  }
}
```

#### File Format
The coefficient file is a plain text file with toml-style section headers and CSV-style rows.

- **CSV-like** (comma or whitespace separated)
- **Section headers**: `[default]`, `[non_coupling]`, `[coupling]`
- **Field count**: `fields: N` at the top
- **Comments**: lines starting with `#`
- **Order**: non-coupling first, then coupling, flattened for the model


For a two-field PFC model, for example:
```
fields: 2

[default]
1.0

[non_coupling]
1.0, 2.0

[coupling]
0.0, 0.5
0.5, 0.0
```

- `[default]` section: Optional section to define the default value of coefficients that are not specified by the file.
- `[non_coupling]` section: List the self-coupling coefficients for each field, separated by commas.
- `[coupling]` section: Each row gives the coupling coefficient from field i to field j.

#### Notes and Best Practices
- The path to the coefficient file is resolved relative to the JSON config file's directory.
- If coefficients_file is specified, it overrides any inline coefficient arrays in the JSON.
- This system is designed for flexibility and easy editing, especially for large or multi-field systems.


## Problem Spec Notation

**Problem Spec Notation** is SymPhas's intuitive domain-specific language for setting up simulation parameters through operator chaining. It provides a natural, mathematical syntax for defining grid dimensions, boundary conditions, and initial conditions without requiring complex constructor calls.

### Overview

Problem Spec Notation enables you to configure three essential components of any phase-field simulation:
- **Grid dimensions and discretization**: Size and spacing of the computational grid
- **Boundary conditions**: How field values behave at domain boundaries  
- **Initial conditions**: Starting field configurations

The notation uses three primary operators:
- **`<<`** and **`<<=`**: Pass data to objects (single value vs. list of values)
- **`/`**: Select specific behaviors or tags 
- **`||`**: Combine boundaries with intervals
- **`*`**: Combine intervals to create multi-dimensional grids

> **⚠️ Requirement**: You must include `using namespace symphas;` before using Problem Spec Notation.

### Grid Dimensions and Discretization

#### Basic Grid Setup

```cpp
using namespace symphas;

// Grid size with literal postfixes
auto grid1d = 100_h;                           // 100 grid points in 1D
auto grid2d = 64_h * 64_h;                     // 64×64 grid
auto grid3d = 32_h * 32_h * 32_h;              // 32×32×32 grid

// With custom discretization (spacing)
auto spaced_interval = 100_h / 0.1_dh;         // 100 points with spacing 0.1
auto spaced_grid = (50_h / 0.5_dh) * (50_h / 0.5_dh);  // 50×50 with spacing 0.5×0.5
```

**Key Points:**
- `_h` suffix specifies grid points (h for "h-spacing")
- `_dh` suffix specifies discretization/spacing 
- Use parentheses when combining spacing to enforce correct order of operations
- Default discretization is 1.0 if not specified

#### Advanced Grid Configuration

```cpp
using namespace symphas;

// Different spacing in each dimension
auto anisotropic = (100_h / 0.1_dh) * (200_h / 0.05_dh);

// Large 3D grid with fine resolution
auto fine_3d = (128_h / 0.01_dh) * (128_h / 0.01_dh) * (64_h / 0.02_dh);
```

### Boundary Conditions

#### Basic Boundary Types

```cpp
using namespace symphas;

// Periodic boundaries (most common)
auto periodic = BoundaryType::PERIODIC;

// Default boundaries (Dirichlet/Neumann conditions)
auto default_bc = BoundaryType::DEFAULT;
```

#### Boundary Tags and Customization

Default boundaries support various tags for specific behaviors:

```cpp
using namespace symphas;

// Constant value boundary
auto constant_bc = BoundaryType::DEFAULT / BoundaryTag::CONSTANT << 1.0;

// Linear boundary condition
auto linear_bc = BoundaryType::DEFAULT / BoundaryTag::LINEAR;

// Time-dependent constant boundary
auto time_bc = BoundaryType::DEFAULT / BoundaryTag::CONSTANT / BoundaryTag::TIME << 2.0;

// Available boundary tags:
// - BoundaryTag::CONSTANT: Fixed value
// - BoundaryTag::LINEAR: Linear variation
// - BoundaryTag::GAUSSIAN: Gaussian profile
// - BoundaryTag::TRIG: Trigonometric (sine/cosine)
// - BoundaryTag::TIME: Time-dependent
// - BoundaryTag::RANDOM: Random values
```

#### Combined Grid and Boundary Setup

```cpp
using namespace symphas;

// Method 1: Using || operator (recommended for uniform boundaries)
auto interval = BoundaryType::PERIODIC || 64_h / 0.5_dh;
auto grid = interval * interval;  // Creates 64×64 grid with periodic boundaries

// Method 2: Interval first, then boundaries
auto base_interval = 50_h / 0.1_dh;
auto grid_with_bc = base_interval * base_interval 
                    << (Side::LEFT / BoundaryType::PERIODIC)
                    << (Side::RIGHT / BoundaryType::PERIODIC);

// Method 3: Mixed boundary conditions
auto mixed_interval = (BoundaryTag::CONSTANT << 0.0) || 100_h / 0.2_dh || (BoundaryTag::CONSTANT << 1.0);
auto mixed_grid = mixed_interval * mixed_interval;
```

**Boundary Side Order:**
When multiplying intervals, boundaries are assigned in this order:
- **First interval**: LEFT and RIGHT sides
- **Second interval**: TOP and BOTTOM sides  
- **Third interval**: FRONT and BACK sides

#### Modifying Specific Boundaries

```cpp
using namespace symphas;

// Start with basic grid
auto base_grid = (BoundaryType::PERIODIC || 64_h / 0.1_dh) * (BoundaryType::PERIODIC || 64_h / 0.1_dh);

// Override specific sides
auto modified_grid = base_grid 
    << (Side::LEFT / (BoundaryTag::CONSTANT << 2.0))
    << (Side::TOP / (BoundaryTag::LINEAR));
```

### Initial Conditions

#### Key-Based Initial Conditions

```cpp
using namespace symphas;

// Uniform random distribution
auto uniform_ic = Inside::UNIFORM <<= {-1.0, 1.0};

// Constant value
auto constant_ic = Inside::CONSTANT << 0.5;

// Geometric shapes
auto circle_ic = Inside::CIRCLE <<= {0.0, 1.0};
auto square_ic = Inside::SQUARE <<= {-0.5, 0.5};

// With variation tags
auto random_circle = Inside::CIRCLE / InsideTag::RANDOM <<= {-1.0, 1.0};
auto varied_square = Inside::SQUARE / InsideTag::VARA <<= {0.0, 1.0};
```

**Available Initial Condition Keys:**
- `Inside::UNIFORM`: Uniform random distribution
- `Inside::CONSTANT`: Constant value throughout domain
- `Inside::CIRCLE`: Circular region
- `Inside::SQUARE`: Square/rectangular region
- `Inside::HEXAGONAL`: Hexagonal pattern
- `Inside::LAMBDA`: Custom lambda function
- `Inside::EXPRESSION`: Expression-based (defined with macros)

**Available Tags:**
- `InsideTag::RANDOM`: Add randomness
- `InsideTag::VARA`, `InsideTag::VARB`: Geometric variations
- `InsideTag::INVERT`: Invert the pattern

#### Expression-Based Initial Conditions

```cpp
// Define the initial condition expression
INITIAL_CONDITION_EQUATION(SINCOS_INIT, (2, 3), sin(x) + cos(y))

// Use in Problem Spec Notation
using namespace symphas;
auto expr_ic = Inside::EXPRESSION << "SINCOS_INIT";
```

#### Lambda-Based Initial Conditions

```cpp
using namespace symphas;

// Define lambda function
auto init_func = [](auto index, const auto* dims, auto dimension) -> double {
    int pos[2]{};
    symphas::lib::position_from_index(pos, index, dims, dimension);
    return std::sin(pos[0] * 0.1) * std::cos(pos[1] * 0.1);
};

// Use in Problem Spec Notation
auto lambda_ic = Inside::LAMBDA << init_func;
```

### Complete Parameter Assembly

```cpp
using namespace symphas;

// Set up grid with boundaries
auto interval = BoundaryType::PERIODIC || 64_h / 0.5_dh;
auto grid = interval * interval;

// Set up initial conditions
auto initial = Inside::UNIFORM <<= {-0.1, 0.1};

// Combine into problem parameters
auto params = grid << initial;

// Use with model
model_MyModel_t<2, SolverFT<Stencil2d2h<>>> model{params};
```

### Multi-Field Systems

For models with multiple fields, provide initial conditions for each:

```cpp
using namespace symphas;

// For SCALARS(2) model
auto ic1 = Inside::UNIFORM <<= {-1.0, 1.0};
auto ic2 = Inside::CIRCLE <<= {0.0, 0.5};

auto params = grid << ic1 << ic2;
```

### Practical Examples

#### Example 1: Basic Periodic Domain
```cpp
using namespace symphas;

// Simple 2D periodic domain
auto interval = BoundaryType::PERIODIC || 64_h / 0.1_dh;
auto grid = interval * interval;
auto params = grid << (Inside::UNIFORM <<= {-0.1, 0.1});
```

#### Example 2: Mixed Boundary Conditions
```cpp
using namespace symphas;

// Fixed boundaries with different values on each side
auto base = 100_h / 0.05_dh;
auto grid = base * base
    << (Side::LEFT / (BoundaryTag::CONSTANT << 0.0))
    << (Side::RIGHT / (BoundaryTag::CONSTANT << 1.0))
    << (Side::TOP / (BoundaryTag::LINEAR))
    << (Side::BOTTOM / (BoundaryTag::CONSTANT << 0.5));

auto params = grid << (Inside::CONSTANT << 0.25);
```

#### Example 3: Complex Initial Conditions
```cpp
using namespace symphas;

INITIAL_CONDITION_EQUATION(WAVE_INIT, (2), sin(2*M_PI*x/L_x) * cos(2*M_PI*y/L_y))

auto grid = (BoundaryType::PERIODIC || 128_h / 0.1_dh) * (BoundaryType::PERIODIC || 128_h / 0.1_dh);
auto params = grid << (Inside::EXPRESSION << "WAVE_INIT");
```

### Best Practices

1. **Always use `using namespace symphas;`** before Problem Spec Notation
2. **Use parentheses** when combining discretization: `(100_h / 0.1_dh)`
3. **Start simple** with periodic boundaries, then add complexity
4. **Test boundary conditions** with simple constant values first
5. **Match initial condition count** to number of fields in your model
6. **Use descriptive names** for expression-based initial conditions

### Common Pitfalls

```cpp
// ❌ Wrong: Missing parentheses
auto bad_grid = 100_h / 0.1_dh * 100_h / 0.1_dh;  // Order of operations issue

// ✅ Correct: Use parentheses  
auto good_grid = (100_h / 0.1_dh) * (100_h / 0.1_dh);

// ❌ Wrong: Forgot namespace
auto interval = BoundaryType::PERIODIC || 64_h / 0.5_dh;  // Won't compile

// ✅ Correct: Include namespace
using namespace symphas;
auto interval = BoundaryType::PERIODIC || 64_h / 0.5_dh;

// ❌ Wrong: Using << instead of <<=
auto ic = Inside::UNIFORM << {-1.0, 1.0};  // Wrong operator

// ✅ Correct: Use <<= for lists
auto ic = Inside::UNIFORM <<= {-1.0, 1.0};
```

## Defining Your Own Solver

SymPhas provides a flexible framework for implementing custom numerical solvers through its Curiously Recurring Template Pattern (CRTP) architecture. This section provides a comprehensive guide to creating your own solver from basic finite difference methods to advanced spectral solvers.

### Solver Architecture Overview

Every solver in SymPhas inherits from the `Solver<Sp, N>` base class using CRTP, where:
- `Sp` is the specialized solver type (your custom solver)
- `N` is the solver system type index (usually 0 unless you need multiple variants)

Solvers are responsible for:
1. **Parsing equations** into solver-specific forms (`form_expr_one`)
2. **Evaluating equations** to compute derivatives (`equation`)
3. **Time stepping** the solution forward (`step`)
4. **Computing spatial derivatives** (optional, can use stencils)

### Basic Solver Structure

There are two main ways to define a solver:

#### Option 1: Solver with Stencil Support
```cpp
START_NEW_SOLVER_WITH_STENCIL(MySolver)

// Implementation functions go here...

END_SOLVER // Don't forget the closing keyword
```

#### Option 2: Solver without Stencils
```cpp
START_NEW_SOLVER(MySolver)

// Implementation functions go here...

END_SOLVER // Don't forget the closing keyword
```

### Required Functions

Every solver must implement these three core functions:

#### 1. `form_expr_one` - Parse and Transform Equations

This function receives an equation of motion and transforms it into a form suitable for your solver:

```cpp
template <size_t En, typename SS, typename S, typename E>
auto form_expr_one(SS&&, std::pair<S, E> const& e) const {
    auto [sys, equation] = e;
    
    // Transform the equation for your solver
    // Examples:
    // - Apply operator optimizations
    // - Split linear/nonlinear terms
    // - Preprocess for time stepping scheme
    
    auto transformed_eq = expr::transform::optimize(
        expr::apply_operators(equation)
    );
    
    return std::make_pair(sys, transformed_eq);
}
```

**Parameters:**
- `En`: Index of the current order parameter (0, 1, 2, ...)
- `systems`: Tuple containing all phase field systems and provisional systems
- `e`: Pair containing the system (`S`) and equation (`E`)

#### 2. `equation` - Evaluate the Transformed Equation

This function evaluates the transformed equation and stores results:

```cpp
// Single equation version
template <typename S, typename E>
void equation(std::pair<S, E>& r) const {
    auto& [sys, eq] = r;
    
    // Evaluate the equation and store in sys.dframe
    expr::result(eq, sys.get().dframe);
}

// Multiple equations version (for parallel processing)
template <typename S, typename E>
void equation(std::pair<S, E>* r, len_type len) const {
    for (iter_type i = 0; i < len; ++i) {
        equation(r[i]);
    }
}
```

#### 3. `step` - Time Integration

This function advances the solution by one time step:

```cpp
template <typename S>
void step(S& sys) const {
    // Update sys.as_grid() using sys.dframe and dt
    // Example: Forward Euler
    expr::result(
        expr::make_term(sys.as_grid()) + expr::make_term(dt, sys.dframe),
        sys.as_grid(),
        expr::iterable_domain(sys.as_grid())
    );
}
```

#### 4. `make_solver` - Static Constructor

This static function creates solver instances from problem parameters:

```cpp
static auto make_solver(symphas::problem_parameters_type const& parameters) {
    if (parameters.length()) {
        double h = parameters.get_interval_data()[0].at(Axis::X).width();
        size_t dim = parameters.get_dimension();
        
        // Get dimensions for stencil solvers
        len_type* dims = new len_type[dim];
        for (iter_type i = 0; i < dim; ++i) {
            Axis side = symphas::index_to_axis(i);
            dims[i] = parameters.get_interval_data()[0].at(side).get_count() 
                     + 2 * BOUNDARY_DEPTH;
        }
        
        auto solver = this_type{dims, h};
        delete[] dims;
        return solver;
    } else {
        return this_type{grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3)};
    }
}
```

### Solver Implementation Examples

#### Example 1: Simple Forward Euler Solver

```cpp
START_NEW_SOLVER_WITH_STENCIL(MyForwardEuler)

// Transform equations (no special processing needed for Forward Euler)
template <size_t En, typename SS, typename S, typename E>
auto form_expr_one(SS&&, std::pair<S, E> const& e) const {
    auto [sys, equation] = e;
    auto optimized_eq = expr::transform::optimize(expr::apply_operators(equation));
    return std::make_pair(sys, optimized_eq);
}

// Evaluate equation into derivative frame
template <typename S, typename E>
void equation(std::pair<S, E>& r) const {
    auto& [sys, eq] = r;
    expr::result(eq, sys.get().dframe);
}

// Forward Euler time step: u^(n+1) = u^n + dt * f(u^n)
template <typename S>
void step(S& sys) const {
    expr::result(
        expr::make_term(sys.as_grid()) + expr::make_term(dt, sys.dframe),
        sys.as_grid(),
        expr::iterable_domain(sys.as_grid())
    );
}

static auto make_solver(symphas::problem_parameters_type const& parameters) {
    if (parameters.length()) {
        double h = parameters.get_interval_data()[0].at(Axis::X).width();
        size_t dim = parameters.get_dimension();
        
        len_type* dims = new len_type[dim];
        for (iter_type i = 0; i < dim; ++i) {
            Axis side = symphas::index_to_axis(i);
            dims[i] = parameters.get_interval_data()[0].at(side).get_count() 
                     + 2 * BOUNDARY_DEPTH;
        }
        
        auto solver = this_type{dims, h};
        delete[] dims;
        return solver;
    } else {
        return this_type{grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3)};
    }
}

END_SOLVER
```

#### Example 2: Runge-Kutta 4th Order Solver

```cpp
START_NEW_SOLVER_WITH_STENCIL(MyRK4)

// Store intermediate stages
struct RK4Data {
    Grid<double, 2> k1, k2, k3, k4, temp;
    
    RK4Data(len_type const* dims) 
        : k1(dims), k2(dims), k3(dims), k4(dims), temp(dims) {}
};

std::unique_ptr<RK4Data> rk_data;

MyRK4(const len_type* dims, double h, double dt = 1.0) 
    : parent_type(dt), rk_data(std::make_unique<RK4Data>(dims)) {}

template <size_t En, typename SS, typename S, typename E>
auto form_expr_one(SS&&, std::pair<S, E> const& e) const {
    auto [sys, equation] = e;
    auto optimized_eq = expr::transform::optimize(expr::apply_operators(equation));
    return std::make_pair(sys, optimized_eq);
}

template <typename S, typename E>
void equation(std::pair<S, E>& r) const {
    auto& [sys, eq] = r;
    // The RK4 method will call this multiple times per step
    expr::result(eq, sys.get().dframe);
}

template <typename S>
void step(S& sys) const {
    auto& u = sys.as_grid();
    auto& k1 = rk_data->k1;
    auto& k2 = rk_data->k2;
    auto& k3 = rk_data->k3;
    auto& k4 = rk_data->k4;
    auto& temp = rk_data->temp;
    
    // k1 = f(u_n) - already computed in equation()
    std::copy(sys.dframe.values, sys.dframe.values + sys.dframe.len, k1.values);
    
    // k2 = f(u_n + dt/2 * k1)
    expr::result(expr::make_term(u) + expr::make_term(dt/2.0, k1), temp);
    // Would need to re-evaluate equation with temp... (simplified here)
    
    // k3 = f(u_n + dt/2 * k2) 
    // k4 = f(u_n + dt * k3)
    
    // Final update: u_{n+1} = u_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    expr::result(
        expr::make_term(u) + expr::make_term(dt/6.0, 
            expr::make_term(k1) + expr::make_term(2.0, k2) + 
            expr::make_term(2.0, k3) + expr::make_term(k4)
        ),
        u
    );
}

static auto make_solver(symphas::problem_parameters_type const& parameters) {
    // Similar implementation to previous examples
    if (parameters.length()) {
        double h = parameters.get_interval_data()[0].at(Axis::X).width();
        size_t dim = parameters.get_dimension();
        
        len_type* dims = new len_type[dim];
        for (iter_type i = 0; i < dim; ++i) {
            Axis side = symphas::index_to_axis(i);
            dims[i] = parameters.get_interval_data()[0].at(side).get_count() 
                     + 2 * BOUNDARY_DEPTH;
        }
        
        auto solver = this_type{dims, h};
        delete[] dims;
        return solver;
    } else {
        return this_type{grid::dim_list(nullptr, 3), grid::h_list(nullptr, 3)};
    }
}

END_SOLVER
```

### Advanced Features

#### Custom Solver Systems

You can associate custom solver systems with your solver for storing additional data:

```cpp
// Define a custom system type
template <typename T, size_t D>
struct MyCustomSolverSystem : SolverSystemFD<T, D> {
    Grid<T, D> auxiliary_data;  // Extra storage
    
    MyCustomSolverSystem(len_type const* dims, double h) 
        : SolverSystemFD<T, D>(dims, h), auxiliary_data(dims) {}
};

// Associate it with your solver
ASSOCIATE_SOLVER_SYSTEM_TYPE(MySolver, MyCustomSolverSystem)
```

#### Supported Types

Specify which field types your solver supports:

```cpp
// Support all types
SYMPHAS_SOLVER_ALL_SUPPORTED(MySolver)

// Or support specific types
SYMPHAS_SOLVER_SUPPORTED_TYPE(MySolver, double)
SYMPHAS_SOLVER_SUPPORTED_TYPE(MySolver, complex_t)
```

#### Custom Derivative Implementations

Override default derivative computations:

```cpp
// In your solver class
template <typename T>
auto applied_laplacian(Block<T> const& e, iter_type n) const {
    // Custom Laplacian implementation
    // Use stencil or custom finite difference formula
    return stencil().applied_laplacian(e, n);
}

template <Axis ax, typename T>
auto applied_gradient(Block<T> const& e, iter_type n) const {
    // Custom gradient implementation
    return stencil().template applied_gradient<ax>(e, n);
}
```

### Integration with SymPhas

#### Using Your Custom Solver

```cpp
// In your main.cpp or model file
#include "my_custom_solver.h"

MODEL(MyModel, (SCALAR),
    EVOLUTION(dop(1) = lap(op(1)) - power(op(1), 3))
)

int main() {
    using namespace symphas;
    
    // Setup parameters
    auto interval = BoundaryType::PERIODIC || 64_h / 0.1_dh;
    auto grid = interval * interval;
    auto params = grid << (Inside::UNIFORM <<= {-0.1, 0.1});
    
    // Use your custom solver
    model_MyModel_t<2, MyForwardEuler<Stencil2d2h<>>> model{params};
    
    // Run simulation
    find_solution(model, 0.01, 1000);
    
    return 0;
}
```

#### CMake Configuration

Include your solver in the solver include file:

```cpp
// In examples/solvers/solverinclude.h (or your custom include)
#pragma once

#include "solverft.h"
#include "solversp.h"
#include "my_custom_solver.h"  // Add your solver
```

### Best Practices and Guidelines

#### 1. **Solver Stability**
- Always consider stability constraints (CFL condition for explicit methods)
- Validate time step sizes in `make_solver`
- Handle boundary conditions appropriately

#### 2. **Performance Optimization**
- Use `expr::transform::optimize()` to optimize expressions
- Leverage parallel processing with OpenMP
- Consider memory access patterns for cache efficiency

#### 3. **Error Handling**
- Validate input parameters in constructors
- Check dimensional consistency
- Provide meaningful error messages

#### 4. **Documentation**
```cpp
//! My custom solver for phase field problems.
/*!
 * This solver implements [describe method] for solving
 * phase field equations. It is suitable for [describe use cases].
 * 
 * Key features:
 * - [Feature 1]
 * - [Feature 2] 
 * 
 * Stability: CFL condition requires dt < h²/(2*D) where D is diffusion coefficient
 */
START_NEW_SOLVER_WITH_STENCIL(MySolver)
// Implementation...
END_SOLVER
```

#### 5. **Testing**
- Test with known analytical solutions
- Verify convergence rates
- Compare with existing solvers on benchmark problems
- Test with different boundary conditions and grid sizes

### Common Pitfalls

1. **Memory Management**: Always clean up dynamically allocated memory
2. **Time Step Handling**: Store `dt` properly and use consistently
3. **Boundary Conditions**: Ensure proper handling of ghost cells
4. **Expression Evaluation**: Use `expr::result()` correctly for storing results
5. **Template Specialization**: Ensure all required template parameters are handled
This comprehensive guide provides the foundation for implementing custom solvers in SymPhas. Start with simple explicit methods and gradually add complexity as needed for your specific phase field problems.

---

## Additional Resources

### Provisional Variable Saving

**New Feature**: Provisional variables can now be saved to disk for post-processing and analysis.

#### Overview
Provisional variables are intermediate calculations (like gradients, laplacians) that are computed during model evolution. Previously, these were only used internally and discarded. Now they can be persisted to disk alongside regular phase field data.

#### Usage

**Save provisional variables only:**
```cpp
model.save_provisional_systems("./output_directory");
model.save_provisional_systems("./output_directory", "custom_prefix_");
```

**Save phase fields and provisional variables together:**
```cpp
model.save_all_systems("./output_directory", true);  // true = include provisional variables
```

**Automatic saving during simulation:**
When using `symphas::find_solution()`, provisional variables are automatically saved with phase field data at each save interval.

#### File Naming
- Default: `provisional_var_1.dat`, `provisional_var_2.dat`, etc.
- Custom: `{prefix}var_1.dat`, `{prefix}var_2.dat`, etc.

#### Benefits
- **Post-processing**: Analyze intermediate quantities like gradients
- **Debugging**: Verify provisional variable calculations
- **Visualization**: Create plots of derived quantities
- **Research**: Study evolution of intermediate expressions

For detailed examples, see `PROVISIONAL_SAVING_USAGE.md`.

### API Reference

For detailed technical documentation of all SymPhas classes, functions, and interfaces, see the comprehensive [API Reference](API_REFERENCE.md). This includes:

- **Core Classes**: Detailed documentation of `Grid`, `Model`, `Solver`, and other fundamental classes
- **Expression System**: Complete reference for symbolic algebra operators and functions
- **Namespace Documentation**: All `symphas`, `expr`, `grid`, and `io` namespace contents
- **Template Parameters**: Detailed explanations of template specializations
- **Function Signatures**: Complete method and function documentation

### Developer's Guide

For advanced customization, extending SymPhas, and contributing to the framework, see the [Developer's Guide](DEVELOPERS_GUIDE.md). This covers:

- **Architecture Details**: Internal framework design and patterns
- **Extension Points**: How to add new solvers, models, and modules
- **Build System**: Advanced CMake configuration and module development
- **Testing Framework**: How to write and run tests
- **Contribution Guidelines**: Code style, review process, and best practices

### Example Projects

Explore the `examples/` directory for comprehensive demonstrations:

- **`examples/tutorial/`**: Basic tutorial with step-by-step introduction
- **`examples/advanced-features/`**: Advanced SymPhas capabilities and techniques
- **`examples/comprehensive-demo/`**: Large-scale simulation examples
- **`examples/using-mpi/`**: Parallel computing with MPI
- **`examples/complex-example/`**: Complex multi-field systems

Each example includes detailed README files with build instructions and explanations.

### Getting Help

- **GitHub Issues**: Report bugs and request features at the [SymPhas repository](https://github.com/SoftSimu/SymPhas)
- **Documentation**: This guide covers most common use cases and configurations
- **Examples**: Check the examples directory for similar use cases
- **API Reference**: Detailed technical documentation for all classes and functions

---

**© SymPhas Documentation Team** | This documentation is part of the SymPhas framework for phase-field simulations.