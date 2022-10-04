# Installation Instructions

## Prerequisites
- [CMake](https://cmake.org) v3.21+
- A C++14 (17 for SYCL builds) compilant C++ compiler

## Optional Dependencies
- A [CUDA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html) compiler, (e.g. `nvcc`, required only if CUDA is enabled)
- A [HIP](https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP-GUIDE.html) compiler, (e.g. `hipcc` or modern `clang++`, required only if HIP is enabled)
- A LLVM SYCL Compiler (e.g. [DPC++](https://github.com/intel/llvm), required only if SYCL enabled)
- Libxc (required only for host API)

## Configuration, Compilation and Installation

ExchCXX provides a CMake build system with automatic dependency management 
(through 
[FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html)).
As such, a simple CMake invocation will often suffice for most purposes
```
cmake -S /path/to/ExchCXX -B /path/to/binary [ExchCXX configure options]
cmake --build /path/to/binary # to build without installation
cmake --build /path/to/binary --target install # to install
``` 


ExchCXX is linkable both as an installed library as well as a CMake subproject via `FetchContent`
```
# ExchCXX Discovery
find_package( exchcxx REQUIRED )
target_link_libraries( my_target PUBLIC exchcxx::exchcxx )
```

```
# ExchCXX as CMake Subproject
include(FetchContent)

# Set ExchCXX CMake options (see below)

# Pull master branch of ExchCXX
FetchContent_Declare( exchcxx 
  GIT_REPOSITORY https://github/com/wavefunction91/ExchCXX.git 
  GIT_TAG master 
)
FetchContent_MakeAvailable( exchcxx )

# Link to target
target_link_libraries( my_target PUBLIC exchcxx::exchcxx )
```

## CMake Variables

| Variable Name                | Description                                                    | Default  |
|-------------------------     |----------------------------------------------------------------|----------|
| `EXCHCXX_ENABLE_TESTS`       | Enable Testing Framework (Catch2)                              | `ON`     |
| `EXCHCXX_ENABLE_BENCHMARK`   | Enable Performance Benchmark                                   | `OFF`    |
| `EXCHCXX_ENABLE_CUDA`        | Enable CUDA XC evaluator                                       | `OFF`    |
| `EXCHCXX_ENABLE_HIP`         | Enable HIP XC evaluator                                        | `OFF`    |
| `EXCHCXX_ENABLE_SYCL`        | Enable SYCL XC evaluator                                       | `OFF`    |

N.B. ExchCXX accelerator bindings are mutally exclusive - this is intentially and will be enforced on configure.
If it would be desirable to support coexistance of these bindings for your application, please open
an [issue](https://github.com/wavefunction91/ExchCXX/issues)

## Build Notes

### Building with CUDA

For CUDA builds, [`CMAKE_CUDA_ARCHITECTURES`](https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html) 
must be set. For example, to build a fat binary for NVIDIA V100 and A100,
```
cmake [source locations] -DCMAKE_CUDA_ARCHITECTURES="70;80" [additional options]
```
ExchCXX supports NVIDIA GPUs with CC >= 60.

### Building with SYCL

SYCL builds require C++17 per the SYCL-2020 standard.
Due to the volitility in SYCL compiler development, ExchCXX does not provide automatic means to determine
the SYCL targets (e.g. NVIDIA-PTX, OpenMP, etc) for the platform of interest. As such, these must be
specified manually in `CMAKE_CXX_FLAGS`. For example, building for a NVIDIA-PTX target with the DPC++ 
compiler may be achieved via
```
cmake [source locations] -DCMAKE_CXX_FLAGS="-fsycl-targets=nvptx64-nvidia-cuda -fsycl" [additional options]
```
