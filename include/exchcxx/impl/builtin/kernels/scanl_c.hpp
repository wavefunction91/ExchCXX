#pragma once

#include <exchcxx/impl/builtin/kernels/scan_c.hpp>
#include <exchcxx/impl/builtin/kernels/pc07_k.hpp>
#include <exchcxx/impl/builtin/kernels/pc07opt_k.hpp>

#include <exchcxx/impl/builtin/kernels/deorbitalized.hpp>

namespace ExchCXX {

#if 0
template <>
struct kernel_traits<BuiltinSCANL_C> { 

  using base_traits = kernel_traits<Deorbitalized<BuiltinSCAN_C, BuiltinPC07OPT_K>>;

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = true;
  static constexpr bool needs_laplacian = true;
  static constexpr double exx_coeff = 0.0;
  
  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps ) {
    base_traits::eval_exc_unpolar_impl(rho, sigma, lapl, tau, eps);
  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau ) {

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps ) {

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b ) {


  }
};
#else

template <>
struct kernel_traits<BuiltinSCANL_C> : 
  public kernel_traits<Deorbitalized<BuiltinSCAN_C, BuiltinPC07OPT_K>>,
  public mgga_screening_interface<BuiltinSCANL_C> {

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = 1e-20;

};
#endif

struct BuiltinSCANL_C : detail::BuiltinKernelImpl< BuiltinSCANL_C > {

  BuiltinSCANL_C( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSCANL_C >(p) { }
  
  virtual ~BuiltinSCANL_C() = default;

};

}
