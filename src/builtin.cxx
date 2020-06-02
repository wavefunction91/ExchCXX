#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/builtin.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

#include <iostream>

namespace ExchCXX {

XCKernel builtin_kernel_factory( XCKernel::Kernel kernel, XCKernel::Spin polar) {

    return XCKernel(
      std::make_unique< detail::BuiltinKernelInterface >( kernel, polar )
    );

}



};
