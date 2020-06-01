#pragma once
#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/libxc.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

#include <unordered_map>
namespace ExchCXX {

static std::unordered_map< XCKernel::Spin, int > libxc_polar_map {
  { XCKernel::Spin::Polarized,   XC_POLARIZED   },
  { XCKernel::Spin::Unpolarized, XC_UNPOLARIZED }
};

}
