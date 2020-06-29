# Synopsis

ExchCXX is a modern C++ library for the evaluation of exchange-correlation (XC)
functionals required for density functional theory (DFT) calculations. In
addition to providing high-level wrappers for the ubiquitous
[Libxc](https://www.tddft.org/programs/libxc) library for the CPU evaluation of
XC functionals, ExchCXX also provides locally developed implementations of a
small subset of XC functionals which may be evaluated either on the host (CPU)
or device (GPU, FPGA, etc). Currently GPU support is provided through the
[CUDA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html)
framework and is only accessible to NVIDIA GPUs. However, work is underway to
provide portable interfaces to these implementations as to allow for use with
emerging exascale and post-exascale compute resources. 


ExchCXX is a work in progress. Its development has been funded by the U.S.
Department of Energy Exascale Computing Project 
([NWChemEx](https://github.com/NWChemEx-Project)).

# Design Goals
* Provide a modern C++ wrapper around various XC functional libraries (Libxc, XCFun, etc)
* Provide stable, portable and high-performance implementations for the evaluation of XC functionals on various accelerator achitectures (GPUs, FPGAs, etc)

# Feature Requests / Contributing

Although ExchCXX may in principle wrap any XC functional implemented by third-party software,
generating these interfaces must currently be done by hand. As such, if there is a functional
that you would like to see wrapped in either host or device code, please create a [GitHub issue](https://github.com/wavefunction91/ExchCXX/issues)
which details your request. Alternatively, I am more than happy to accept [pull-requests](https://github.com/wavefunction91/ExchCXX/pulls) if
you have made the changes locally and would like to see them integrated in upstream.


# License

ExchCXX is made freely available under the terms of the 3-Clause BSD license. See
LICENSE.txt for details.

# Acknowledgements

The development of ExchCXX is supported by the Exascale Computing Project
(17-SC-20-SC), a collaborative effort of the U.S. Department of Energy Office
of Science and the National Nuclear Security Administration.
