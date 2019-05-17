#ifndef __INCLUDED_FACTORY_XC_KERNEL_HPP__
#define __INCLUDED_FACTORY_XC_KERNEL_HPP__

#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

XCKernel libxc_kernel_factory(const std::string&, const XCKernel::Spin );

}; // namespace ExchCXX

#endif
