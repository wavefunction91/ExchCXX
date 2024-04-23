#pragma once

#include <exchcxx/impl/builtin/kernels/scan_x.hpp>
#include <exchcxx/impl/builtin/kernels/pc07_k.hpp>
#include <exchcxx/impl/builtin/kernels/pc07opt_k.hpp>

#include <exchcxx/impl/builtin/kernels/deorbitalized.hpp>

namespace ExchCXX {


template <>
struct kernel_traits<BuiltinSCANL_X> : 
  public kernel_traits<Deorbitalized<BuiltinSCAN_X, BuiltinPC07OPT_K>> {

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = 1e-20;

};

struct BuiltinSCANL_X : detail::BuiltinKernelImpl< BuiltinSCANL_X > {

  BuiltinSCANL_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSCANL_X >(p) { }
  
  virtual ~BuiltinSCANL_X() = default;

};

}
