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

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {

template <typename KernelType>
struct lda_screening_interface {

  using traits = kernel_traits<KernelType>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double& eps ) {

    if( rho <= traits::dens_tol ) {
      eps = 0.;
    } else {
      traits::eval_exc_unpolar_impl( rho, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double& eps ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    if( rho_s <= traits::dens_tol ) {
      eps = 0.;
    } else {
      rho_a = safe_max(rho_a, traits::dens_tol);
      rho_b = safe_max(rho_b, traits::dens_tol);
      traits::eval_exc_polar_impl( rho_a, rho_b, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double& eps, double& vrho ) {

    if( rho <= traits::dens_tol ) {
      eps  = 0.;
      vrho = 0.;
    } else {
      traits::eval_exc_vxc_unpolar_impl( rho, eps, vrho );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double& eps, double& vrho_a,
      double& vrho_b ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    if( rho_s <= traits::dens_tol ) {
      eps    = 0.;
      vrho_a = 0.;
      vrho_b = 0.;
    } else {
      rho_a = safe_max(rho_a, traits::dens_tol);
      rho_b = safe_max(rho_b, traits::dens_tol);
      traits::eval_exc_vxc_polar_impl( rho_a, rho_b, eps, vrho_a, vrho_b );
    }

  }

BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar( double rho, double& v2rho2 ) {

    if( rho <= traits::dens_tol ) {
      v2rho2 = 0.;
    } else {
      traits::eval_fxc_unpolar_impl( rho, v2rho2 );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar( double rho_a, double rho_b, double& v2rho2_aa,
                    double& v2rho2_ab, double& v2rho2_bb ) {

    const double rho_s = rho_a + rho_b;

    if( rho_s <= traits::dens_tol ) {
      v2rho2_aa = 0.;
      v2rho2_ab = 0.;
      v2rho2_bb = 0.;
    } else {
      rho_a = safe_max(rho_a, traits::dens_tol);
      rho_b = safe_max(rho_b, traits::dens_tol);
      traits::eval_fxc_polar_impl( rho_a, rho_b, v2rho2_aa, v2rho2_ab, v2rho2_bb );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar( double rho, double& vrho, double& v2rho2 ) {

    if( rho <= traits::dens_tol ) {
      vrho = 0.;
      v2rho2 = 0.;
    } else {
      traits::eval_vxc_fxc_unpolar_impl( rho, vrho, v2rho2 );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar( double rho_a, double rho_b, double& vrho_a,
                        double& vrho_b, double& v2rho2_aa,
                        double& v2rho2_ab, double& v2rho2_bb ) {

    const double rho_s = rho_a + rho_b;

    if( rho_s <= traits::dens_tol ) {
      vrho_a = 0.;
      vrho_b = 0.;
      v2rho2_aa = 0.;
      v2rho2_ab = 0.;
      v2rho2_bb = 0.;
    } else {
      rho_a = safe_max(rho_a, traits::dens_tol);
      rho_b = safe_max(rho_b, traits::dens_tol);
      traits::eval_vxc_fxc_polar_impl( rho_a, rho_b, vrho_a, vrho_b, 
                                       v2rho2_aa, v2rho2_ab, v2rho2_bb );
    }

  }

};


template <typename KernelType>
struct gga_screening_interface {

  using traits = kernel_traits<KernelType>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    if( rho <= traits::dens_tol ) {
      eps = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      traits::eval_exc_unpolar_impl( rho, sigma, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    if( rho_s <= traits::dens_tol ) {
      eps = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);
      traits::eval_exc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab,
        sigma_bb, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

    if( rho <= traits::dens_tol ) {
      eps    = 0.;
      vrho   = 0.;
      vsigma = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      traits::eval_exc_vxc_unpolar_impl( rho, sigma, eps, vrho, vsigma );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps, double& vrho_a,
      double& vrho_b, double& vsigma_aa, double& vsigma_ab, 
      double& vsigma_bb ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    eps       = 0.;
    vrho_a    = 0.;
    vrho_b    = 0.;
    vsigma_aa = 0.;
    vsigma_ab = 0.;
    vsigma_bb = 0.;

    if( rho_s > traits::dens_tol ) {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);
      traits::eval_exc_vxc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab, 
        sigma_bb, eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );
    }

  }

BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar( double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    if( rho <= traits::dens_tol ) {
      v2rho2 = 0.;
      v2rhosigma = 0.;
      v2sigma2 = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      traits::eval_fxc_unpolar_impl( rho, sigma, v2rho2, v2rhosigma, v2sigma2 );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb,
                   double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb,
                   double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb,
                   double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb,
                   double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb,
                   double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    const double rho_s = rho_a + rho_b;

    if( rho_s <= traits::dens_tol ) {
      v2rho2_aa = 0.; v2rho2_ab = 0.; v2rho2_bb = 0.;
      v2rhosigma_a_aa = 0.; v2rhosigma_a_ab = 0.; v2rhosigma_a_bb = 0.;
      v2rhosigma_b_aa = 0.; v2rhosigma_b_ab = 0.; v2rhosigma_b_bb = 0.;
      v2sigma2_aa_aa = 0.; v2sigma2_aa_ab = 0.; v2sigma2_aa_bb = 0.;
      v2sigma2_ab_ab = 0.; v2sigma2_ab_bb = 0.; v2sigma2_bb_bb = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);
      traits::eval_fxc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb,
                                  v2rho2_aa, v2rho2_ab, v2rho2_bb,
                                  v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb,
                                  v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb,
                                  v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb,
                                  v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar( double rho, double sigma, double& vrho, double& vsigma,
                         double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    if( rho <= traits::dens_tol ) {
      vrho = 0.; vsigma = 0.;
      v2rho2 = 0.; v2rhosigma = 0.; v2sigma2 = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      traits::eval_vxc_fxc_unpolar_impl( rho, sigma, vrho, vsigma,
                                        v2rho2, v2rhosigma, v2sigma2 );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb,
                       double& vrho_a, double& vrho_b, 
                       double& vsigma_aa, double& vsigma_ab, double& vsigma_bb,
                       double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb,
                       double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb,
                       double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb,
                       double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb,
                       double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    const double rho_s = rho_a + rho_b;

    vrho_a = 0.; vrho_b = 0.;
    vsigma_aa = 0.; vsigma_ab = 0.; vsigma_bb = 0.;
    v2rho2_aa = 0.; v2rho2_ab = 0.; v2rho2_bb = 0.;
    v2rhosigma_a_aa = 0.; v2rhosigma_a_ab = 0.; v2rhosigma_a_bb = 0.;
    v2rhosigma_b_aa = 0.; v2rhosigma_b_ab = 0.; v2rhosigma_b_bb = 0.;
    v2sigma2_aa_aa = 0.; v2sigma2_aa_ab = 0.; v2sigma2_aa_bb = 0.;
    v2sigma2_ab_ab = 0.; v2sigma2_ab_bb = 0.; v2sigma2_bb_bb = 0.;

    if( rho_s > traits::dens_tol ) {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);
      traits::eval_vxc_fxc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb,
                                      vrho_a, vrho_b, 
                                      vsigma_aa, vsigma_ab, vsigma_bb,
                                      v2rho2_aa, v2rho2_ab, v2rho2_bb,
                                      v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb,
                                      v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb,
                                      v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb,
                                      v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb );
    }

  }

};

template <typename KernelType>
struct mgga_screening_interface {

  using traits = kernel_traits<KernelType>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double lapl, double tau, double& eps ) {

    
    if( rho <= traits::dens_tol ) {
      eps = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      tau   = safe_max(tau, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma = enforce_fermi_hole_curvature(sigma, rho, tau);
      }

      traits::eval_exc_unpolar_impl(rho, sigma, lapl, tau, eps) ;
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double lapl_a, double lapl_b,
      double tau_a, double tau_b, double& eps ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    if( rho_s <= traits::dens_tol ) {
      eps = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      tau_a = safe_max(tau_a, traits::tau_tol);
      tau_b = safe_max(tau_b, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma_aa = enforce_fermi_hole_curvature(sigma_aa, rho_a, tau_a);
        sigma_bb = enforce_fermi_hole_curvature(sigma_bb, rho_b, tau_b);
      }
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);
    
      traits::eval_exc_polar_impl( rho_a, rho_b, 
        sigma_aa, sigma_ab, sigma_bb, lapl_a, lapl_b, 
        safe_max(traits::tau_tol, tau_a), 
        safe_max(traits::tau_tol, tau_b), 
        eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double lapl, double tau, 
    double& eps, double& vrho, double& vsigma, double& vlapl, double &vtau ) {

    if( rho <= traits::dens_tol ) {
      eps    = 0.;
      vrho   = 0.;
      vsigma = 0.;
      vlapl  = 0.;
      vtau   = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      tau   = safe_max(tau, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma = enforce_fermi_hole_curvature(sigma, rho, tau);
      }
      traits::eval_exc_vxc_unpolar_impl(rho, sigma, lapl, tau, 
        eps, vrho, vsigma, vlapl, vtau );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double lapl_a, double lapl_b,
      double tau_a, double tau_b, double& eps, double& vrho_a,
      double& vrho_b, double& vsigma_aa, double& vsigma_ab, 
      double& vsigma_bb, double& vlapl_a, double& vlapl_b,
      double& vtau_a, double& vtau_b ) {

    const double rho_s = rho_a + rho_b;
    //const double rho_z = rho_a - rho_b;

    eps       = 0.;
    vrho_a    = 0.;
    vrho_b    = 0.;
    vsigma_aa = 0.;
    vsigma_ab = 0.;
    vsigma_bb = 0.;
    vlapl_a   = 0.;
    vlapl_b   = 0.;
    vtau_a    = 0.;
    vtau_b    = 0.;

    if( rho_s > traits::dens_tol ) {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      tau_a = safe_max(tau_a, traits::tau_tol);
      tau_b = safe_max(tau_b, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma_aa = enforce_fermi_hole_curvature(sigma_aa, rho_a, tau_a);
        sigma_bb = enforce_fermi_hole_curvature(sigma_bb, rho_b, tau_b);
      }
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);

      traits::eval_exc_vxc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab, 
        sigma_bb, lapl_a, lapl_b, 
        safe_max(traits::tau_tol, tau_a), 
        safe_max(traits::tau_tol, tau_b), 
        eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vlapl_a, vlapl_b, 
        vtau_a, vtau_b );
    }

  }

BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar( double rho, double sigma, double lapl, double tau, 
                     double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, 
                     double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, 
                     double& v2lapl2, double& v2lapltau, double& v2tau2 ) {

    if( rho <= traits::dens_tol ) {
      v2rho2 = 0.; v2rhosigma = 0.; v2rholapl = 0.; v2rhotau = 0.;
      v2sigma2 = 0.; v2sigmalapl = 0.; v2sigmatau = 0.;
      v2lapl2 = 0.; v2lapltau = 0.; v2tau2 = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      tau   = safe_max(tau, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma = enforce_fermi_hole_curvature(sigma, rho, tau);
      }
      traits::eval_fxc_unpolar_impl(rho, sigma, lapl, tau, 
        v2rho2, v2rhosigma, v2rholapl, v2rhotau, 
        v2sigma2, v2sigmalapl, v2sigmatau, 
        v2lapl2, v2lapltau, v2tau2 );
    }

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

    const double rho_s = rho_a + rho_b;

    // Initialize all values to zero
    v2rho2_aa = 0.; v2rho2_ab = 0.; v2rho2_bb = 0.;
    v2rhosigma_a_aa = 0.; v2rhosigma_a_ab = 0.; v2rhosigma_a_bb = 0.;
    v2rhosigma_b_aa = 0.; v2rhosigma_b_ab = 0.; v2rhosigma_b_bb = 0.;
    v2rholapl_a_a = 0.; v2rholapl_a_b = 0.; v2rholapl_b_a = 0.; v2rholapl_b_b = 0.;
    v2rhotau_a_a = 0.; v2rhotau_a_b = 0.; v2rhotau_b_a = 0.; v2rhotau_b_b = 0.;
    v2sigma2_aa_aa = 0.; v2sigma2_aa_ab = 0.; v2sigma2_aa_bb = 0.;
    v2sigma2_ab_ab = 0.; v2sigma2_ab_bb = 0.; v2sigma2_bb_bb = 0.;
    v2sigmalapl_aa_a = 0.; v2sigmalapl_aa_b = 0.; v2sigmalapl_ab_a = 0.; 
    v2sigmalapl_ab_b = 0.; v2sigmalapl_bb_a = 0.; v2sigmalapl_bb_b = 0.;
    v2sigmatau_aa_a = 0.; v2sigmatau_aa_b = 0.; v2sigmatau_ab_a = 0.; 
    v2sigmatau_ab_b = 0.; v2sigmatau_bb_a = 0.; v2sigmatau_bb_b = 0.;
    v2lapl2_aa = 0.; v2lapl2_ab = 0.; v2lapl2_bb = 0.;
    v2lapltau_a_a = 0.; v2lapltau_a_b = 0.; v2lapltau_b_a = 0.; v2lapltau_b_b = 0.;
    v2tau2_aa = 0.; v2tau2_ab = 0.; v2tau2_bb = 0.;

    if( rho_s > traits::dens_tol ) {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      tau_a = safe_max(tau_a, traits::tau_tol);
      tau_b = safe_max(tau_b, traits::tau_tol);
      if constexpr (not traits::is_kedf) {
        sigma_aa = enforce_fermi_hole_curvature(sigma_aa, rho_a, tau_a);
        sigma_bb = enforce_fermi_hole_curvature(sigma_bb, rho_b, tau_b);
      }
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);

      traits::eval_fxc_polar_impl( rho_a, rho_b, 
        sigma_aa, sigma_ab, sigma_bb, 
        lapl_a, lapl_b, 
        safe_max(traits::tau_tol, tau_a), 
        safe_max(traits::tau_tol, tau_b), 
        v2rho2_aa, v2rho2_ab, v2rho2_bb,
        v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb,
        v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb,
        v2rholapl_a_a, v2rholapl_a_b, v2rholapl_b_a, v2rholapl_b_b,
        v2rhotau_a_a, v2rhotau_a_b, v2rhotau_b_a, v2rhotau_b_b,
        v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb,
        v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb,
        v2sigmalapl_aa_a, v2sigmalapl_aa_b, v2sigmalapl_ab_a, 
        v2sigmalapl_ab_b, v2sigmalapl_bb_a, v2sigmalapl_bb_b,
        v2sigmatau_aa_a, v2sigmatau_aa_b, v2sigmatau_ab_a, 
        v2sigmatau_ab_b, v2sigmatau_bb_a, v2sigmatau_bb_b,
        v2lapl2_aa, v2lapl2_ab, v2lapl2_bb,
        v2lapltau_a_a, v2lapltau_a_b, v2lapltau_b_a, v2lapltau_b_b,
        v2tau2_aa, v2tau2_ab, v2tau2_bb );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar( double rho, double sigma, double lapl, double tau,
                         double& vrho, double& vsigma, double& vlapl, double& vtau,
                         double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau,
                         double& v2sigma2, double& v2sigmalapl, double& v2sigmatau,
                         double& v2lapl2, double& v2lapltau, double& v2tau2 ) {

    if( rho <= traits::dens_tol ) {
      vrho = 0.; vsigma = 0.; vlapl = 0.; vtau = 0.;
      v2rho2 = 0.; v2rhosigma = 0.; v2rholapl = 0.; v2rhotau = 0.;
      v2sigma2 = 0.; v2sigmalapl = 0.; v2sigmatau = 0.;
      v2lapl2 = 0.; v2lapltau = 0.; v2tau2 = 0.;
    } else {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      sigma = safe_max(sigma, sigma_tol_sq); 
      tau = safe_max(tau, traits::tau_tol);
      
      if constexpr (not traits::is_kedf) {
        sigma = enforce_fermi_hole_curvature(sigma, rho, tau);
      }
      
      traits::eval_vxc_fxc_unpolar_impl(rho, sigma, lapl, tau,
                                       vrho, vsigma, vlapl, vtau,
                                       v2rho2, v2rhosigma, v2rholapl, v2rhotau,
                                       v2sigma2, v2sigmalapl, v2sigmatau,
                                       v2lapl2, v2lapltau, v2tau2);
    }

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

    const double rho_s = rho_a + rho_b;

    // Initialize all values to zero
    vrho_a = 0.; vrho_b = 0.;
    vsigma_aa = 0.; vsigma_ab = 0.; vsigma_bb = 0.;
    vlapl_a = 0.; vlapl_b = 0.; vtau_a = 0.; vtau_b = 0.;
    v2rho2_aa = 0.; v2rho2_ab = 0.; v2rho2_bb = 0.;
    v2rhosigma_a_aa = 0.; v2rhosigma_a_ab = 0.; v2rhosigma_a_bb = 0.;
    v2rhosigma_b_aa = 0.; v2rhosigma_b_ab = 0.; v2rhosigma_b_bb = 0.;
    v2rholapl_a_a = 0.; v2rholapl_a_b = 0.; v2rholapl_b_a = 0.; v2rholapl_b_b = 0.;
    v2rhotau_a_a = 0.; v2rhotau_a_b = 0.; v2rhotau_b_a = 0.; v2rhotau_b_b = 0.;
    v2sigma2_aa_aa = 0.; v2sigma2_aa_ab = 0.; v2sigma2_aa_bb = 0.;
    v2sigma2_ab_ab = 0.; v2sigma2_ab_bb = 0.; v2sigma2_bb_bb = 0.;
    v2sigmalapl_aa_a = 0.; v2sigmalapl_aa_b = 0.; v2sigmalapl_ab_a = 0.; 
    v2sigmalapl_ab_b = 0.; v2sigmalapl_bb_a = 0.; v2sigmalapl_bb_b = 0.;
    v2sigmatau_aa_a = 0.; v2sigmatau_aa_b = 0.; v2sigmatau_ab_a = 0.; 
    v2sigmatau_ab_b = 0.; v2sigmatau_bb_a = 0.; v2sigmatau_bb_b = 0.;
    v2lapl2_aa = 0.; v2lapl2_ab = 0.; v2lapl2_bb = 0.;
    v2lapltau_a_a = 0.; v2lapltau_a_b = 0.; v2lapltau_b_a = 0.; v2lapltau_b_b = 0.;
    v2tau2_aa = 0.; v2tau2_ab = 0.; v2tau2_bb = 0.;

    if( rho_s > traits::dens_tol ) {
      constexpr auto sigma_tol_sq = traits::sigma_tol * traits::sigma_tol;
      rho_a    = safe_max(rho_a, traits::dens_tol);
      rho_b    = safe_max(rho_b, traits::dens_tol);
      sigma_aa = safe_max(sigma_aa, sigma_tol_sq); 
      sigma_bb = safe_max(sigma_bb, sigma_tol_sq); 
      tau_a = safe_max(tau_a, traits::tau_tol);
      tau_b = safe_max(tau_b, traits::tau_tol);
      
      if constexpr (not traits::is_kedf) {
        sigma_aa = enforce_fermi_hole_curvature(sigma_aa, rho_a, tau_a);
        sigma_bb = enforce_fermi_hole_curvature(sigma_bb, rho_b, tau_b);
      }
      
      sigma_ab = enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb);

      traits::eval_vxc_fxc_polar_impl( rho_a, rho_b,
        sigma_aa, sigma_ab, sigma_bb,
        lapl_a, lapl_b,
        safe_max(traits::tau_tol, tau_a),
        safe_max(traits::tau_tol, tau_b),
        vrho_a, vrho_b,
        vsigma_aa, vsigma_ab, vsigma_bb,
        vlapl_a, vlapl_b, vtau_a, vtau_b,
        v2rho2_aa, v2rho2_ab, v2rho2_bb,
        v2rhosigma_a_aa, v2rhosigma_a_ab, v2rhosigma_a_bb,
        v2rhosigma_b_aa, v2rhosigma_b_ab, v2rhosigma_b_bb,
        v2rholapl_a_a, v2rholapl_a_b, v2rholapl_b_a, v2rholapl_b_b,
        v2rhotau_a_a, v2rhotau_a_b, v2rhotau_b_a, v2rhotau_b_b,
        v2sigma2_aa_aa, v2sigma2_aa_ab, v2sigma2_aa_bb,
        v2sigma2_ab_ab, v2sigma2_ab_bb, v2sigma2_bb_bb,
        v2sigmalapl_aa_a, v2sigmalapl_aa_b, v2sigmalapl_ab_a,
        v2sigmalapl_ab_b, v2sigmalapl_bb_a, v2sigmalapl_bb_b,
        v2sigmatau_aa_a, v2sigmatau_aa_b, v2sigmatau_ab_a,
        v2sigmatau_ab_b, v2sigmatau_bb_a, v2sigmatau_bb_b,
        v2lapl2_aa, v2lapl2_ab, v2lapl2_bb,
        v2lapltau_a_a, v2lapltau_a_b, v2lapltau_b_a, v2lapltau_b_b,
        v2tau2_aa, v2tau2_ab, v2tau2_bb );
    }

  }

};

}
