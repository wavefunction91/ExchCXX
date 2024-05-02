#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/slater_exchange.hpp>
#include <exchcxx/impl/builtin/kernels/b88.hpp>
#include <exchcxx/impl/builtin/kernels/vwn_rpa.hpp>
#include <exchcxx/impl/builtin/kernels/lyp.hpp>

namespace ExchCXX {

template <>
struct kernel_traits<BuiltinB3LYP> {

  static constexpr bool is_hyb  = true;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool is_epc  = false;
  static constexpr double exx_coeff = 0.20;

  static constexpr double b3lyp_ax = 0.72; // % GGA X
  static constexpr double b3lyp_ac = 0.81; // % GGA C
  
  static constexpr double b3lyp_slater_coeff = 1. - exx_coeff - b3lyp_ax;
  static constexpr double b3lyp_vwn_coeff    = 1. - b3lyp_ac;

  using slater_traits = kernel_traits<BuiltinSlaterExchange>;
  using b88_traits    = kernel_traits<BuiltinB88>;
  using vwn_traits    = kernel_traits<BuiltinVWN_RPA>;
  using lyp_traits    = kernel_traits<BuiltinLYP>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    slater_traits::eval_exc_unpolar( rho, eps );
    double slater_eps = eps;
  
    b88_traits::eval_exc_unpolar( rho, sigma, eps );
    double b88_eps = eps;

    vwn_traits::eval_exc_unpolar( rho, eps );
    double vwn_eps = eps;

    lyp_traits::eval_exc_unpolar( rho, sigma, eps );
    double lyp_eps = eps;

    eps = b3lyp_slater_coeff * slater_eps + b3lyp_ax * b88_eps + 
          b3lyp_vwn_coeff    * vwn_eps    + b3lyp_ac * lyp_eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps ) {

    slater_traits::eval_exc_polar( rho_a, rho_b, eps );
    double slater_eps = eps;
  
    b88_traits::eval_exc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, eps );
    double b88_eps = eps;

    vwn_traits::eval_exc_polar( rho_a, rho_b, eps );
    double vwn_eps = eps;

    lyp_traits::eval_exc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, eps );
    double lyp_eps = eps;

    eps = b3lyp_slater_coeff * slater_eps + b3lyp_ax * b88_eps + 
          b3lyp_vwn_coeff    * vwn_eps    + b3lyp_ac * lyp_eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

    slater_traits::eval_exc_vxc_unpolar( rho, eps, vrho );
    double slater_eps  = eps;
    double slater_vrho = vrho;
  
    b88_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );
    double b88_eps    = eps;
    double b88_vrho   = vrho;
    double b88_vsigma = vsigma;

    vwn_traits::eval_exc_vxc_unpolar( rho, eps, vrho );
    double vwn_eps  = eps;
    double vwn_vrho = vrho;

    lyp_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );
    double lyp_eps    = eps;
    double lyp_vrho   = vrho;
    double lyp_vsigma = vsigma;

    eps = b3lyp_slater_coeff * slater_eps + b3lyp_ax * b88_eps + 
          b3lyp_vwn_coeff    * vwn_eps    + b3lyp_ac * lyp_eps;

    vrho = b3lyp_slater_coeff * slater_vrho + b3lyp_ax * b88_vrho + 
           b3lyp_vwn_coeff    * vwn_vrho    + b3lyp_ac * lyp_vrho;

    vsigma = b3lyp_ax * b88_vsigma +  b3lyp_ac * lyp_vsigma;

  }



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, double sigma_ab, 
      double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa,
      double& vsigma_ab, double& vsigma_bb ) {

    slater_traits::eval_exc_vxc_polar( rho_a, rho_b, eps, vrho_a, vrho_b );
    double slater_eps  = eps;
    double slater_vrho_a = vrho_a;
    double slater_vrho_b = vrho_b;
  
    b88_traits::eval_exc_vxc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, 
      eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );
    double b88_eps    = eps;
    double b88_vrho_a   = vrho_a;
    double b88_vrho_b   = vrho_b;
    double b88_vsigma_aa = vsigma_aa;
    double b88_vsigma_ab = vsigma_ab;
    double b88_vsigma_bb = vsigma_bb;

    vwn_traits::eval_exc_vxc_polar( rho_a, rho_b, eps, vrho_a, vrho_b );
    double vwn_eps  = eps;
    double vwn_vrho_a = vrho_a;
    double vwn_vrho_b = vrho_b;

    lyp_traits::eval_exc_vxc_polar( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, 
      eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );
    double lyp_eps    = eps;
    double lyp_vrho_a   = vrho_a;
    double lyp_vrho_b   = vrho_b;
    double lyp_vsigma_aa = vsigma_aa;
    double lyp_vsigma_ab = vsigma_ab;
    double lyp_vsigma_bb = vsigma_bb;

    eps = b3lyp_slater_coeff * slater_eps + b3lyp_ax * b88_eps + 
          b3lyp_vwn_coeff    * vwn_eps    + b3lyp_ac * lyp_eps;

    vrho_a = b3lyp_slater_coeff * slater_vrho_a + b3lyp_ax * b88_vrho_a + 
             b3lyp_vwn_coeff    * vwn_vrho_a    + b3lyp_ac * lyp_vrho_a;
    vrho_b = b3lyp_slater_coeff * slater_vrho_b + b3lyp_ax * b88_vrho_b + 
             b3lyp_vwn_coeff    * vwn_vrho_b    + b3lyp_ac * lyp_vrho_b;

    vsigma_aa = b3lyp_ax * b88_vsigma_aa +  b3lyp_ac * lyp_vsigma_aa;
    vsigma_ab = b3lyp_ax * b88_vsigma_ab +  b3lyp_ac * lyp_vsigma_ab;
    vsigma_bb = b3lyp_ax * b88_vsigma_bb +  b3lyp_ac * lyp_vsigma_bb;

  }

};






struct BuiltinB3LYP : detail::BuiltinKernelImpl< BuiltinB3LYP > {

  BuiltinB3LYP( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinB3LYP >(p) { }
  
  virtual ~BuiltinB3LYP() = default;

};

}
