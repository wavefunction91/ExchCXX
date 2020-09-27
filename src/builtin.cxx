#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/builtin/interface.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX {

XCKernel builtin_kernel_factory( Kernel kernel, Spin polar) {

    return XCKernel(
      std::make_unique< detail::BuiltinKernelInterface >( kernel, polar )
    );

}



}
