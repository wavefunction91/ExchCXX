#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/pbe_x.hpp>
#include <exchcxx/impl/builtin/kernels/pbe_c.hpp>

namespace ExchCXX {

template <>
struct kernel_traits<BuiltinPBE0> {

  static constexpr bool is_hyb  = true;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.25;
  static constexpr double dens_tol  = 1e-32; // not used, delegates to PBE_C/X

  using pbe_x_traits = kernel_traits<BuiltinPBE_X>;
  using pbe_c_traits = kernel_traits<BuiltinPBE_C>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    pbe_x_traits::eval_exc_unpolar( rho, sigma, eps );
    double eps_x = eps;

    pbe_c_traits::eval_exc_unpolar( rho, sigma, eps );

    eps = (1. - exx_coeff) * eps_x + eps;
  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

    pbe_x_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );
    double eps_x    = eps;
    double vrho_x   = vrho;
    double vsigma_x = vsigma;


    pbe_c_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );

    eps    = (1. - exx_coeff) * eps_x    + eps;
    vrho   = (1. - exx_coeff) * vrho_x   + vrho;
    vsigma = (1. - exx_coeff) * vsigma_x + vsigma;

  }

};






struct BuiltinPBE0 : detail::BuiltinKernelImpl< BuiltinPBE0 > {

  BuiltinPBE0( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE0 >(p) { }
  
  virtual ~BuiltinPBE0() = default;

};

}
