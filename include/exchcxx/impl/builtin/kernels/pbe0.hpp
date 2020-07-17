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

  using pbe_x_traits = kernel_traits<BuiltinPBE_X>;
  using pbe_c_traits = kernel_traits<BuiltinPBE_C>;

  static constexpr bool is_hyb  = true;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.25;


  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    pbe_x_traits::eval_exc_unpolar( rho, sigma, eps );
    double eps_x = eps;

    pbe_c_traits::eval_exc_unpolar( rho, sigma, eps );

    eps = (1. - exx_coeff) * eps_x + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps ) {

    pbe_x_traits::eval_exc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, eps );
    double eps_x = eps;

    pbe_c_traits::eval_exc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, eps );

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

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps, double& vrho_a,
      double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    pbe_x_traits::eval_exc_vxc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb,
      eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );
    double eps_x    = eps;
    double vrho_x_a   = vrho_a;
    double vrho_x_b   = vrho_b;
    double vsigma_x_aa = vsigma_aa;
    double vsigma_x_ab = vsigma_ab;
    double vsigma_x_bb = vsigma_bb;


    pbe_c_traits::eval_exc_vxc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb,
      eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );

    eps       = (1. - exx_coeff) * eps_x       + eps;
    vrho_a    = (1. - exx_coeff) * vrho_x_a    + vrho_a;
    vrho_b    = (1. - exx_coeff) * vrho_x_b    + vrho_b;
    vsigma_aa = (1. - exx_coeff) * vsigma_x_aa + vsigma_aa;
    vsigma_ab = (1. - exx_coeff) * vsigma_x_ab + vsigma_ab;
    vsigma_bb = (1. - exx_coeff) * vsigma_x_bb + vsigma_bb;

  }

};






struct BuiltinPBE0 : detail::BuiltinKernelImpl< BuiltinPBE0 > {

  BuiltinPBE0( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE0 >(p) { }
  
  virtual ~BuiltinPBE0() = default;

};

}
