/**
 * ExchCXX 
 *
 * Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). 
 *
 * Portions Copyright (c) Microsoft Corporation.
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#pragma once 
#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>
#include <exchcxx/util/unused.hpp>

namespace ExchCXX {

template <typename XCEF, typename KEDF>
struct kernel_traits<Deorbitalized<XCEF,KEDF>> {

  using xc_traits = kernel_traits<XCEF>;
  using ke_traits = kernel_traits<KEDF>;

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = true;
  static constexpr bool needs_laplacian = true;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;
  static constexpr double exx_coeff = xc_traits::exx_coeff + ke_traits::exx_coeff;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double lapl, double tau, double& eps ) {

    double TAU;
    ke_traits::eval_exc_unpolar(rho, sigma, lapl, 0.0, TAU);

    TAU = TAU * rho;
    xc_traits::eval_exc_unpolar(rho, sigma, lapl, TAU, eps);

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau ) {

    double TAU, vrho_k, vsigma_k, vlapl_k, dummy;
    ke_traits::eval_exc_vxc_unpolar(rho, sigma, lapl, 0.0, TAU, vrho_k, vsigma_k, vlapl_k, dummy);


    TAU = TAU * rho;
    xc_traits::eval_exc_vxc_unpolar(rho, sigma, lapl, TAU, eps, vrho, vsigma, vlapl, vtau);

    vrho   += vtau * vrho_k;
    vsigma += vtau * vsigma_k;
    vlapl   = vtau * vlapl_k;
    vtau    = 0.0;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps ) {

    double TAU_A, TAU_B;
    ke_traits::eval_exc_polar(rho_a, 0.0, sigma_aa, 0.0, 0.0, lapl_a, 0.0, 0.0, 0.0, TAU_A);
    ke_traits::eval_exc_polar(rho_b, 0.0, sigma_bb, 0.0, 0.0, lapl_b, 0.0, 0.0, 0.0, TAU_B);

    TAU_A *= rho_a;
    TAU_B *= rho_b;
    xc_traits::eval_exc_polar(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, TAU_A, TAU_B, eps);

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b ) {

    double TAU_A, TAU_B, vrho_a_k, vrho_b_k, vsigma_aa_k, vsigma_bb_k, vlapl_a_k, vlapl_b_k, vtau_k, dummy;
    ke_traits::eval_exc_vxc_polar(rho_a, 0.0, sigma_aa, 0.0, 0.0, lapl_a, 0.0, 0.0, 0.0, TAU_A, vrho_a_k, dummy, vsigma_aa_k, dummy, dummy, vlapl_a_k, dummy, dummy, dummy);
    ke_traits::eval_exc_vxc_polar(rho_b, 0.0, sigma_bb, 0.0, 0.0, lapl_b, 0.0, 0.0, 0.0, TAU_B, vrho_b_k, dummy, vsigma_bb_k, dummy, dummy, vlapl_b_k, dummy, dummy, dummy);

    TAU_A *= rho_a;
    TAU_B *= rho_b;

    xc_traits::eval_exc_vxc_polar(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, TAU_A, TAU_B, eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vlapl_a, vlapl_b, vtau_a, vtau_b);

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

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar( double rho, double sigma, double lapl, double tau, double& vrho, double& vsigma, double& vlapl, double& vtau,
                         double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau,
                         double& v2sigma2, double& v2sigmalapl, double& v2sigmatau,
                         double& v2lapl2, double& v2lapltau, double& v2tau2 ) {
    #if defined(__CUDACC__) || defined(__HIPCC__)
    printf("eval_vxc_fxc_unpolar not implemented for deorbitalized kernels\n");
    #elif defined(__SYCL_DEVICE_ONLY__) || defined(EXCHCXX_ENABLE_SYCL)
    sycl::ext::oneapi::experimental::printf("eval_vxc_fxc_unpolar not implemented for deorbitalized kernels\n");
    #else
    unused(rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);
    throw std::runtime_error("eval_vxc_fxc_unpolar not implemented for deorbitalized kernels");
    #endif
  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar( double rho_a, double rho_b,
                       double sigma_aa, double sigma_ab, double sigma_bb,
                       double lapl_a, double lapl_b, double tau_a, double tau_b,
                       double& vrho_a, double& vrho_b,
                       double& vsigma_aa, double& vsigma_ab, double& vsigma_bb,
                       double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b,
                       double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb,
                       double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb,
                       double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb,
                       double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b,
                       double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b,
                       double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb,
                       double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb,
                       double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, 
                       double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b,
                       double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, 
                       double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b,
                       double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb,
                       double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b,
                       double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb ) {
    #if defined(__CUDACC__) || defined(__HIPCC__)
    printf("eval_vxc_fxc_polar not implemented for deorbitalized kernels\n");
    #elif defined(__SYCL_DEVICE_ONLY__) || defined(EXCHCXX_ENABLE_SYCL)
    sycl::ext::oneapi::experimental::printf("eval_vxc_fxc_polar not implemented for deorbitalized kernels\n");
    #else
      unused(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, tau_a, tau_b, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vlapl_a, vlapl_b, vtau_a, vtau_b, v2rho2_aa, v2rho2_ab, v2rho2_bb, v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb, v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb, v2rholapl_a_a, v2rholapl_a_b, v2rholapl_b_a, v2rholapl_b_b, v2rhotau_a_a, v2rhotau_a_b, v2rhotau_b_a, v2rhotau_b_b, v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb, v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb, v2sigmalapl_aa_a, v2sigmalapl_aa_b, v2sigmalapl_ab_a, v2sigmalapl_ab_b, v2sigmalapl_bb_a, v2sigmalapl_bb_b, v2sigmatau_aa_a, v2sigmatau_aa_b, v2sigmatau_ab_a, v2sigmatau_ab_b, v2sigmatau_bb_a, v2sigmatau_bb_b, v2lapl2_aa, v2lapl2_ab, v2lapl2_bb, v2lapltau_a_a, v2lapltau_a_b, v2lapltau_b_a, v2lapltau_b_b, v2tau2_aa, v2tau2_ab, v2tau2_bb);
      throw std::runtime_error("eval_vxc_fxc_polar not implemented for deorbitalized kernels");
    #endif

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar( double rho, double sigma, double lapl, double tau,
                     double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, 
                     double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, 
                     double& v2lapl2, double& v2lapltau, double& v2tau2 ) {
    #if defined(__CUDACC__) || defined(__HIPCC__)
    printf("eval_fxc_unpolar not implemented for deorbitalized kernels\n");
    #elif defined(__SYCL_DEVICE_ONLY__) || defined(EXCHCXX_ENABLE_SYCL)
    sycl::ext::oneapi::experimental::printf("eval_fxc_unpolar not implemented for deorbitalized kernels\n");
    #else
    unused(rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);
    throw std::runtime_error("eval_fxc_unpolar not implemented for deorbitalized kernels");
    #endif
  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar( double rho_a, double rho_b, 
                   double sigma_aa, double sigma_ab, double sigma_bb, 
                   double lapl_a, double lapl_b, double tau_a, double tau_b,
                   double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb,
                   double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb,
                   double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb,
                   double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b,
                   double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b,
                   double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb,
                   double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb,
                   double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, 
                   double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b,
                   double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, 
                   double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b,
                   double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb,
                   double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b,
                   double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb ) {
    #if defined(__CUDACC__) || defined(__HIPCC__)
    printf("eval_fxc_polar not implemented for deorbitalized kernels\n");
    #elif defined(__SYCL_DEVICE_ONLY__) || defined(EXCHCXX_ENABLE_SYCL)
    sycl::ext::oneapi::experimental::printf("eval_fxc_polar not implemented for deorbitalized kernels\n");
    #else
    unused(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, tau_a, tau_b, v2rho2_aa, v2rho2_ab, v2rho2_bb, v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb, v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb, v2rholapl_a_a, v2rholapl_a_b, v2rholapl_b_a, v2rholapl_b_b, v2rhotau_a_a, v2rhotau_a_b, v2rhotau_b_a, v2rhotau_b_b, v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb, v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb, v2sigmalapl_aa_a, v2sigmalapl_aa_b, v2sigmalapl_ab_a, v2sigmalapl_ab_b, v2sigmalapl_bb_a, v2sigmalapl_bb_b, v2sigmatau_aa_a, v2sigmatau_aa_b, v2sigmatau_ab_a, v2sigmatau_ab_b, v2sigmatau_bb_a, v2sigmatau_bb_b, v2lapl2_aa, v2lapl2_ab, v2lapl2_bb, v2lapltau_a_a, v2lapltau_a_b, v2lapltau_b_a, v2lapltau_b_b, v2tau2_aa, v2tau2_ab, v2tau2_bb);
    throw std::runtime_error("eval_fxc_polar not implemented for deorbitalized kernels");
    #endif
  }

};

template <typename XCEF, typename KEDF>
struct Deorbitalized : detail::BuiltinKernelImpl< Deorbitalized<XCEF, KEDF> > {

  Deorbitalized( Spin p ) :
    detail::BuiltinKernelImpl< Deorbitalized<XCEF, KEDF> >(p) { }

  virtual ~Deorbitalized() = default;

};

}
