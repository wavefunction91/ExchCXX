#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif


#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

XCKernel libxc_kernel_factory(const Kernel, const Spin );
XCKernel builtin_kernel_factory( Kernel, Spin );

static inline XCKernel kernel_factory( 
  Backend backend, Kernel kern, Spin polar
) {

  if( backend == Backend::libxc )
    return libxc_kernel_factory( kern, polar );
  else
    return builtin_kernel_factory( kern, polar );

}

} // namespace ExchCXX

