#pragma once

#include <exchcxx/impl/builtin/kernels/r2scan_c.hpp>
#include <exchcxx/impl/builtin/kernels/pc07_k.hpp>
#include <exchcxx/impl/builtin/kernels/pc07opt_k.hpp>

#include <exchcxx/impl/builtin/kernels/deorbitalized.hpp>

namespace ExchCXX {


template <>
struct kernel_traits<BuiltinR2SCANL_C> : 
  public kernel_traits<Deorbitalized<BuiltinR2SCAN_C, BuiltinPC07OPT_K>> {

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = 1e-20;

};

struct BuiltinR2SCANL_C : detail::BuiltinKernelImpl< BuiltinR2SCANL_C > {

  BuiltinR2SCANL_C( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinR2SCANL_C >(p) { }
  
  virtual ~BuiltinR2SCANL_C() = default;

};

}
