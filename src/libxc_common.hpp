#pragma once
#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/libxc.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

#include <unordered_map>
namespace ExchCXX {

static std::unordered_map< Spin, int > libxc_polar_map {
  { Spin::Polarized,   XC_POLARIZED   },
  { Spin::Unpolarized, XC_UNPOLARIZED }
};

}
