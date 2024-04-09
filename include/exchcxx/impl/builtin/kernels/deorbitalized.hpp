#pragma once 
#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>

namespace ExchCXX {


template <typename XCEF, typename KEDF>
struct kernel_traits<Deorbitalized<XCEF,KEDF>> {

  using xc_traits = kernel_traits<XCEF>;
  using ke_traits = kernel_traits<KEDF>;

  static constexpr bool is_hyb  = xc_traits::is_hyb  or ke_traits::is_hyb;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = true;
  static constexpr bool needs_laplacian = true;
  static constexpr bool is_kedf = false;
  static constexpr double exx_coeff = xc_traits::exx_coeff + ke_traits::exx_coeff;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps ) {

    double TAU;
    ke_traits::eval_exc_unpolar_impl(rho, sigma, lapl, tau, TAU);
    xc_traits::eval_exc_unpolar_impl(rho, sigma, lapl, TAU, eps);

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau ) {

    double TAU, vrho_k, vsigma_k, vlapl_k, vtau_k;
    ke_traits::eval_exc_vxc_unpolar_impl(rho, sigma, lapl, tau, TAU, vrho_k, vsigma_k, vlapl_k, vtau_k);
    xc_traits::eval_exc_vxc_unpolar_impl(rho, sigma, lapl, TAU, eps, vrho, vsigma, vlapl, vtau);

    vrho   += vtau * vrho_k;
    vsigma += vtau * vsigma_k;
    vlapl   = vtau * vlapl_k;
    vtau    = 0.0;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps ) {

    double TAU_A, TAU_B;
    ke_traits::eval_exc_polar_impl(rho_a, 0.0, sigma_aa, 0.0, 0.0, lapl_a, 0.0, 0.0, 0.0, TAU_A);
    ke_traits::eval_exc_polar_impl(rho_b, 0.0, sigma_bb, 0.0, 0.0, lapl_b, 0.0, 0.0, 0.0, TAU_B);

    xc_traits::eval_exc_polar_impl(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, TAU_A, TAU_B, eps);

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b ) {

    double TAU_A, TAU_B, vrho_a_k, vrho_b_k, vsigma_aa_k, vsigma_bb_k, vlapl_a_k, vlapl_b_k, vtau_k, dummy;
    ke_traits::eval_exc_vxc_polar_impl(rho_a, 0.0, sigma_aa, 0.0, 0.0, lapl_a, 0.0, 0.0, 0.0, TAU_A, vrho_a_k, dummy, vsigma_aa_k, dummy, dummy, vlapl_a_k, dummy, dummy, dummy);
    ke_traits::eval_exc_vxc_polar_impl(rho_b, 0.0, sigma_bb, 0.0, 0.0, lapl_b, 0.0, 0.0, 0.0, TAU_B, vrho_b_k, dummy, vsigma_bb_k, dummy, dummy, vlapl_b_k, dummy, dummy, dummy);

    xc_traits::eval_exc_vxc_polar_impl(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, TAU_A, TAU_B, eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vlapl_a, vlapl_b, vtau_a, vtau_b);

    vrho_a    += vtau_a * vrho_a_k;
    vrho_b    += vtau_b * vrho_b_k;
    vsigma_aa += vtau_a * vsigma_aa_k;
    //vsigma_ab += .....;
    vsigma_bb += vtau_b * vsigma_bb_k;
    vlapl_a    = vtau_a * vlapl_a_k;
    vlapl_b    = vtau_b * vlapl_b_k;
    vtau_a     = 0.0;
    vtau_b     = 0.0;
  }


};

template <typename XCEF, typename KEDF>
struct Deorbitalized : detail::BuiltinKernelImpl< Deorbitalized<XCEF, KEDF> > {

  Deorbitalized( Spin p ) :
    detail::BuiltinKernelImpl< Deorbitalized<XCEF, KEDF> >(p) { }

  virtual ~Deorbitalized() = default;

};

}
