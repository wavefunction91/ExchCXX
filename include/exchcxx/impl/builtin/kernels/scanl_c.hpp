#pragma once

#include <exchcxx/impl/builtin/kernels/scan_c.hpp>
#include <exchcxx/impl/builtin/kernels/pc07_k.hpp>
#include <exchcxx/impl/builtin/kernels/pc07opt_k.hpp>

#include <exchcxx/impl/builtin/kernels/deorbitalized.hpp>

namespace ExchCXX {


template <>
struct kernel_traits<BuiltinSCANL_C> : 
  public kernel_traits<Deorbitalized<BuiltinSCAN_C, BuiltinPC07OPT_K>> {

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = 1e-20;

};

struct BuiltinSCANL_C : detail::BuiltinKernelImpl< BuiltinSCANL_C > {

  BuiltinSCANL_C( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSCANL_C >(p) { }
  
  virtual ~BuiltinSCANL_C() = default;

};

}
