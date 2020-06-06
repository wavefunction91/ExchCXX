#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif


#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

XCKernel libxc_kernel_factory(const XCKernel::Kernel, const XCKernel::Spin );
XCKernel builtin_kernel_factory( XCKernel::Kernel, XCKernel::Spin );

static inline XCKernel kernel_factory( 
  XCKernel::Backend backend, XCKernel::Kernel kern, XCKernel::Spin polar
) {

  if( backend == XCKernel::Backend::libxc )
    return libxc_kernel_factory( kern, polar );
  else
    return builtin_kernel_factory( kern, polar );

}

} // namespace ExchCXX

