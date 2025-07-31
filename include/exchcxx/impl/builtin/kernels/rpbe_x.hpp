/**
 * ExchCXX 
 *
 * Copyright (c) Microsoft Corporation.
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

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinRPBE_X > :
  public gga_screening_interface< BuiltinRPBE_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;


  static constexpr double rpbe_kappa = 0.8040;
  static constexpr double rpbe_mu = 0.2195149727645171;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = rpbe_mu * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t34 = 0.1e1 / rpbe_kappa;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t39 = safe_math::exp( -t21 * t25 * t29 * t33 * t34 / 0.24e2 );
    const double t42 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t39 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );


    eps = 0.2e1 * t46;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = rpbe_mu * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t34 = 0.1e1 / rpbe_kappa;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t39 = safe_math::exp( -t21 * t25 * t29 * t33 * t34 / 0.24e2 );
    const double t42 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t39 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );
    const double t52 = t30 * rho;
    const double t55 = t17 / t18 / t52;
    const double t59 = t29 * t39;
    const double t60 = t20 * t25 * t59;
    const double t64 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t6 * t55 * rpbe_mu * t60 / 0.24e2 );
    const double t72 = t25 * t28 * t39;
    const double t73 = t21 * t72;
    const double t76 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t73 / 0.64e2 );


    eps = 0.2e1 * t46;
    vrho = 0.2e1 * rho * t64 + 0.2e1 * t46;
    vsigma = 0.2e1 * rho * t76;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = rpbe_mu * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t34 = 0.1e1 / rpbe_kappa;
    constexpr double t96 = rpbe_mu * rpbe_mu;
    constexpr double t99 = t20 * t20;
    constexpr double t102 = t99 / t23 / t22;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t39 = safe_math::exp( -t21 * t25 * t29 * t33 * t34 / 0.24e2 );
    const double t42 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t39 );
    const double t52 = t30 * rho;
    const double t55 = t17 / t18 / t52;
    const double t59 = t29 * t39;
    const double t60 = t20 * t25 * t59;
    const double t64 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t6 * t55 * rpbe_mu * t60 / 0.24e2 );
    const double t72 = t25 * t28 * t39;
    const double t73 = t21 * t72;
    const double t76 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t73 / 0.64e2 );
    const double t85 = t30 * t30;
    const double t88 = t17 / t18 / t85;
    const double t93 = t85 * t52;
    const double t98 = t6 * t17 / t93 * t96;
    const double t103 = sigma * sigma;
    const double t106 = t27 * t34 * t39;
    const double t107 = t102 * t103 * t106;
    const double t111 = piecewise_functor_3( t2, 0.0, t6 * t17 / t31 / rho * t42 / 0.12e2 - t6 * t88 * rpbe_mu * t60 / 0.8e1 + t98 * t107 / 0.108e3 );
    const double t117 = t85 * t30;
    const double t121 = t6 * t17 / t117 * t96;
    const double t125 = t102 * t27 * sigma * t34 * t39;
    const double t129 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.192e3 * t6 * t55 * t73 - t121 * t125 / 0.288e3 );
    const double t132 = t85 * rho;
    const double t137 = t102 * t106;
    const double t140 = piecewise_functor_3( t2, 0.0, t6 * t17 / t132 * t96 * t137 / 0.768e3 );


    v2rho2 = 0.2e1 * rho * t111 + 0.4e1 * t64;
    v2rhosigma = 0.2e1 * rho * t129 + 0.2e1 * t76;
    v2sigma2 = 0.2e1 * rho * t140;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double sigma, double& vrho, double& vsigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = rpbe_mu * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t34 = 0.1e1 / rpbe_kappa;
    constexpr double t96 = rpbe_mu * rpbe_mu;
    constexpr double t99 = t20 * t20;
    constexpr double t102 = t99 / t23 / t22;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t39 = safe_math::exp( -t21 * t25 * t29 * t33 * t34 / 0.24e2 );
    const double t42 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t39 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );
    const double t52 = t30 * rho;
    const double t55 = t17 / t18 / t52;
    const double t59 = t29 * t39;
    const double t60 = t20 * t25 * t59;
    const double t64 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t6 * t55 * rpbe_mu * t60 / 0.24e2 );
    const double t72 = t25 * t28 * t39;
    const double t73 = t21 * t72;
    const double t76 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t73 / 0.64e2 );
    const double t85 = t30 * t30;
    const double t88 = t17 / t18 / t85;
    const double t93 = t85 * t52;
    const double t98 = t6 * t17 / t93 * t96;
    const double t103 = sigma * sigma;
    const double t106 = t27 * t34 * t39;
    const double t107 = t102 * t103 * t106;
    const double t111 = piecewise_functor_3( t2, 0.0, t6 * t17 / t31 / rho * t42 / 0.12e2 - t6 * t88 * rpbe_mu * t60 / 0.8e1 + t98 * t107 / 0.108e3 );
    const double t117 = t85 * t30;
    const double t121 = t6 * t17 / t117 * t96;
    const double t125 = t102 * t27 * sigma * t34 * t39;
    const double t129 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.192e3 * t6 * t55 * t73 - t121 * t125 / 0.288e3 );
    const double t132 = t85 * rho;
    const double t137 = t102 * t106;
    const double t140 = piecewise_functor_3( t2, 0.0, t6 * t17 / t132 * t96 * t137 / 0.768e3 );


    vrho = 0.2e1 * rho * t64 + 0.2e1 * t46;
    vsigma = 0.2e1 * rho * t76;
    v2rho2 = 0.2e1 * rho * t111 + 0.4e1 * t64;
    v2rhosigma = 0.2e1 * rho * t129 + 0.2e1 * t76;
    v2sigma2 = 0.2e1 * rho * t140;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = rpbe_mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t41 = 0.1e1 / rpbe_kappa;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t45 = safe_math::exp( -t34 * sigma_aa * t39 * t41 / 0.24e2 );
    const double t48 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t45 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = rho_b * rho_b;
    const double t64 = safe_math::cbrt( rho_b );
    const double t65 = t64 * t64;
    const double t67 = 0.1e1 / t65 / t63;
    const double t72 = safe_math::exp( -t34 * sigma_bb * t67 * t41 / 0.24e2 );
    const double t75 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t72 );
    const double t79 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t75 );


    eps = t52 + t79;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = rpbe_mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t41 = 0.1e1 / rpbe_kappa;
    constexpr double t100 = t28 * t33;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t45 = safe_math::exp( -t34 * sigma_aa * t39 * t41 / 0.24e2 );
    const double t48 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t45 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = rho_b * rho_b;
    const double t64 = safe_math::cbrt( rho_b );
    const double t65 = t64 * t64;
    const double t67 = 0.1e1 / t65 / t63;
    const double t72 = safe_math::exp( -t34 * sigma_bb * t67 * t41 / 0.24e2 );
    const double t75 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t72 );
    const double t79 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t75 );
    const double t80 = t6 * t6;
    const double t81 = 0.1e1 / t80;
    const double t82 = t16 * t81;
    const double t84 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t82 );
    const double t87 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t84 );
    const double t88 = t87 * t26;
    const double t92 = t26 * t26;
    const double t93 = 0.1e1 / t92;
    const double t94 = t25 * t93;
    const double t97 = t5 * t94 * t48 / 0.8e1;
    const double t99 = t5 * t27 * rpbe_mu;
    const double t101 = t35 * rho_a;
    const double t103 = 0.1e1 / t37 / t101;
    const double t106 = t100 * sigma_aa * t103 * t45;
    const double t110 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t88 * t48 - t97 + t99 * t106 / 0.24e2 );
    const double t111 = t54 * t81;
    const double t113 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t111 );
    const double t116 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t113 );
    const double t117 = t116 * t26;
    const double t121 = t61 * t93;
    const double t124 = t5 * t121 * t75 / 0.8e1;
    const double t126 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t75 - t124 );
    const double t130 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t82 );
    const double t133 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t130 );
    const double t134 = t133 * t26;
    const double t139 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t48 - t97 );
    const double t141 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t111 );
    const double t144 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t141 );
    const double t145 = t144 * t26;
    const double t150 = t5 * t62 * rpbe_mu;
    const double t151 = t63 * rho_b;
    const double t153 = 0.1e1 / t65 / t151;
    const double t156 = t100 * sigma_bb * t153 * t72;
    const double t160 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t145 * t75 - t124 + t150 * t156 / 0.24e2 );
    const double t163 = t5 * t27;
    const double t166 = t29 * t33 * t39 * t45;
    const double t169 = piecewise_functor_3( t1, 0.0, -t163 * t166 / 0.64e2 );
    const double t170 = t5 * t62;
    const double t173 = t29 * t33 * t67 * t72;
    const double t176 = piecewise_functor_3( t53, 0.0, -t170 * t173 / 0.64e2 );


    eps = t52 + t79;
    vrho_a = t52 + t79 + t6 * ( t110 + t126 );
    vrho_b = t52 + t79 + t6 * ( t139 + t160 );
    vsigma_aa = t6 * t169;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t176;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = rpbe_mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t41 = 0.1e1 / rpbe_kappa;
    constexpr double t100 = t28 * t33;
    constexpr double t224 = rpbe_mu * rpbe_mu;
    constexpr double t227 = t28 * t28;
    constexpr double t230 = t227 / t31 / t30;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t45 = safe_math::exp( -t34 * sigma_aa * t39 * t41 / 0.24e2 );
    const double t48 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t45 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = rho_b * rho_b;
    const double t64 = safe_math::cbrt( rho_b );
    const double t65 = t64 * t64;
    const double t67 = 0.1e1 / t65 / t63;
    const double t72 = safe_math::exp( -t34 * sigma_bb * t67 * t41 / 0.24e2 );
    const double t75 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t72 );
    const double t80 = t6 * t6;
    const double t81 = 0.1e1 / t80;
    const double t82 = t16 * t81;
    const double t84 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t82 );
    const double t87 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t84 );
    const double t88 = t87 * t26;
    const double t92 = t26 * t26;
    const double t93 = 0.1e1 / t92;
    const double t94 = t25 * t93;
    const double t97 = t5 * t94 * t48 / 0.8e1;
    const double t99 = t5 * t27 * rpbe_mu;
    const double t101 = t35 * rho_a;
    const double t103 = 0.1e1 / t37 / t101;
    const double t106 = t100 * sigma_aa * t103 * t45;
    const double t110 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t88 * t48 - t97 + t99 * t106 / 0.24e2 );
    const double t111 = t54 * t81;
    const double t113 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t111 );
    const double t116 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t113 );
    const double t117 = t116 * t26;
    const double t121 = t61 * t93;
    const double t124 = t5 * t121 * t75 / 0.8e1;
    const double t126 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t75 - t124 );
    const double t130 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t82 );
    const double t133 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t130 );
    const double t134 = t133 * t26;
    const double t139 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t48 - t97 );
    const double t141 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t111 );
    const double t144 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t141 );
    const double t145 = t144 * t26;
    const double t150 = t5 * t62 * rpbe_mu;
    const double t151 = t63 * rho_b;
    const double t153 = 0.1e1 / t65 / t151;
    const double t156 = t100 * sigma_bb * t153 * t72;
    const double t160 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t145 * t75 - t124 + t150 * t156 / 0.24e2 );
    const double t163 = t5 * t27;
    const double t166 = t29 * t33 * t39 * t45;
    const double t169 = piecewise_functor_3( t1, 0.0, -t163 * t166 / 0.64e2 );
    const double t170 = t5 * t62;
    const double t173 = t29 * t33 * t67 * t72;
    const double t176 = piecewise_functor_3( t53, 0.0, -t170 * t173 / 0.64e2 );
    const double t179 = t23 * t23;
    const double t180 = 0.1e1 / t179;
    const double t181 = t84 * t84;
    const double t184 = t80 * t6;
    const double t185 = 0.1e1 / t184;
    const double t186 = t16 * t185;
    const double t189 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t81 + 0.2e1 * t186 );
    const double t193 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t180 * t181 + 0.4e1 / 0.3e1 * t23 * t189 );
    const double t194 = t193 * t26;
    const double t198 = t87 * t93;
    const double t200 = t5 * t198 * t48;
    const double t203 = t5 * t88 * rpbe_mu;
    const double t207 = 0.1e1 / t92 / t6;
    const double t208 = t25 * t207;
    const double t211 = t5 * t208 * t48 / 0.12e2;
    const double t213 = t5 * t94 * rpbe_mu;
    const double t214 = t213 * t106;
    const double t216 = t35 * t35;
    const double t218 = 0.1e1 / t37 / t216;
    const double t221 = t100 * sigma_aa * t218 * t45;
    const double t226 = t5 * t27 * t224;
    const double t231 = sigma_aa * sigma_aa;
    const double t232 = t230 * t231;
    const double t235 = 0.1e1 / t36 / t216 / t101;
    const double t237 = t235 * t41 * t45;
    const double t238 = t232 * t237;
    const double t242 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t194 * t48 - t200 / 0.4e1 + t203 * t106 / 0.12e2 + t211 + t214 / 0.36e2 - 0.11e2 / 0.72e2 * t99 * t221 + t226 * t238 / 0.216e3 );
    const double t243 = t59 * t59;
    const double t244 = 0.1e1 / t243;
    const double t245 = t113 * t113;
    const double t248 = t54 * t185;
    const double t251 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t81 + 0.2e1 * t248 );
    const double t255 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t244 * t245 + 0.4e1 / 0.3e1 * t59 * t251 );
    const double t256 = t255 * t26;
    const double t260 = t116 * t93;
    const double t262 = t5 * t260 * t75;
    const double t264 = t61 * t207;
    const double t267 = t5 * t264 * t75 / 0.12e2;
    const double t269 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t256 * t75 - t262 / 0.4e1 + t267 );
    const double t285 = t133 * t93;
    const double t287 = t5 * t285 * t48;
    const double t310 = t144 * t93;
    const double t312 = t5 * t310 * t75;
    const double t320 = t5 * t121 * rpbe_mu;
    const double t321 = t320 * t156;
    const double t329 = t130 * t130;
    const double t334 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t81 + 0.2e1 * t186 );
    const double t338 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t180 * t329 + 0.4e1 / 0.3e1 * t23 * t334 );
    const double t339 = t338 * t26;
    const double t345 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t339 * t48 - t287 / 0.4e1 + t211 );
    const double t346 = t141 * t141;
    const double t351 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t81 + 0.2e1 * t248 );
    const double t355 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t244 * t346 + 0.4e1 / 0.3e1 * t59 * t351 );
    const double t356 = t355 * t26;
    const double t362 = t5 * t145 * rpbe_mu;
    const double t366 = t63 * t63;
    const double t368 = 0.1e1 / t65 / t366;
    const double t371 = t100 * sigma_bb * t368 * t72;
    const double t375 = t5 * t62 * t224;
    const double t376 = sigma_bb * sigma_bb;
    const double t377 = t230 * t376;
    const double t380 = 0.1e1 / t64 / t366 / t151;
    const double t382 = t380 * t41 * t72;
    const double t383 = t377 * t382;
    const double t387 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t356 * t75 - t312 / 0.4e1 + t362 * t156 / 0.12e2 + t267 + t321 / 0.36e2 - 0.11e2 / 0.72e2 * t150 * t371 + t375 * t383 / 0.216e3 );
    const double t390 = t5 * t88;
    const double t393 = t5 * t94;
    const double t395 = t393 * t166 / 0.192e3;
    const double t398 = t29 * t33 * t103 * t45;
    const double t401 = t216 * t35;
    const double t403 = 0.1e1 / t36 / t401;
    const double t406 = sigma_aa * t41 * t45;
    const double t407 = t230 * t403 * t406;
    const double t411 = piecewise_functor_3( t1, 0.0, -t390 * t166 / 0.64e2 - t395 + t163 * t398 / 0.24e2 - t226 * t407 / 0.576e3 );
    const double t416 = t5 * t121;
    const double t418 = t416 * t173 / 0.192e3;
    const double t428 = t5 * t145;
    const double t433 = t29 * t33 * t153 * t72;
    const double t436 = t366 * t63;
    const double t438 = 0.1e1 / t64 / t436;
    const double t441 = sigma_bb * t41 * t72;
    const double t442 = t230 * t438 * t441;
    const double t446 = piecewise_functor_3( t53, 0.0, -t428 * t173 / 0.64e2 - t418 + t170 * t433 / 0.24e2 - t375 * t442 / 0.576e3 );
    const double t448 = t216 * rho_a;
    const double t453 = t230 / t36 / t448 * t41 * t45;
    const double t456 = piecewise_functor_3( t1, 0.0, t226 * t453 / 0.1536e4 );
    const double t457 = t366 * rho_b;
    const double t462 = t230 / t64 / t457 * t41 * t72;
    const double t465 = piecewise_functor_3( t53, 0.0, t375 * t462 / 0.1536e4 );


    v2rho2_aa = 0.2e1 * t110 + 0.2e1 * t126 + t6 * ( t242 + t269 );
    v2rho2_bb = 0.2e1 * t139 + 0.2e1 * t160 + t6 * ( t345 + t387 );
    v2rhosigma_a_aa = t6 * t411 + t169;
    v2rhosigma_b_bb = t6 * t446 + t176;
    v2sigma2_aa_aa = t6 * t456;
    v2sigma2_bb_bb = t6 * t465;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;



  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = rpbe_mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t41 = 0.1e1 / rpbe_kappa;
    constexpr double t100 = t28 * t33;
    constexpr double t224 = rpbe_mu * rpbe_mu;
    constexpr double t227 = t28 * t28;
    constexpr double t230 = t227 / t31 / t30;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t45 = safe_math::exp( -t34 * sigma_aa * t39 * t41 / 0.24e2 );
    const double t48 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t45 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = rho_b * rho_b;
    const double t64 = safe_math::cbrt( rho_b );
    const double t65 = t64 * t64;
    const double t67 = 0.1e1 / t65 / t63;
    const double t72 = safe_math::exp( -t34 * sigma_bb * t67 * t41 / 0.24e2 );
    const double t75 = 0.1e1 + rpbe_kappa * ( 0.1e1 - t72 );
    const double t79 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t75 );
    const double t80 = t6 * t6;
    const double t81 = 0.1e1 / t80;
    const double t82 = t16 * t81;
    const double t84 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t82 );
    const double t87 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t84 );
    const double t88 = t87 * t26;
    const double t92 = t26 * t26;
    const double t93 = 0.1e1 / t92;
    const double t94 = t25 * t93;
    const double t97 = t5 * t94 * t48 / 0.8e1;
    const double t99 = t5 * t27 * rpbe_mu;
    const double t101 = t35 * rho_a;
    const double t103 = 0.1e1 / t37 / t101;
    const double t106 = t100 * sigma_aa * t103 * t45;
    const double t110 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t88 * t48 - t97 + t99 * t106 / 0.24e2 );
    const double t111 = t54 * t81;
    const double t113 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t111 );
    const double t116 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t113 );
    const double t117 = t116 * t26;
    const double t121 = t61 * t93;
    const double t124 = t5 * t121 * t75 / 0.8e1;
    const double t126 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t75 - t124 );
    const double t130 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t82 );
    const double t133 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t130 );
    const double t134 = t133 * t26;
    const double t139 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t48 - t97 );
    const double t141 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t111 );
    const double t144 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t141 );
    const double t145 = t144 * t26;
    const double t150 = t5 * t62 * rpbe_mu;
    const double t151 = t63 * rho_b;
    const double t153 = 0.1e1 / t65 / t151;
    const double t156 = t100 * sigma_bb * t153 * t72;
    const double t160 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t145 * t75 - t124 + t150 * t156 / 0.24e2 );
    const double t163 = t5 * t27;
    const double t166 = t29 * t33 * t39 * t45;
    const double t169 = piecewise_functor_3( t1, 0.0, -t163 * t166 / 0.64e2 );
    const double t170 = t5 * t62;
    const double t173 = t29 * t33 * t67 * t72;
    const double t176 = piecewise_functor_3( t53, 0.0, -t170 * t173 / 0.64e2 );
    const double t179 = t23 * t23;
    const double t180 = 0.1e1 / t179;
    const double t181 = t84 * t84;
    const double t184 = t80 * t6;
    const double t185 = 0.1e1 / t184;
    const double t186 = t16 * t185;
    const double t189 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t81 + 0.2e1 * t186 );
    const double t193 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t180 * t181 + 0.4e1 / 0.3e1 * t23 * t189 );
    const double t194 = t193 * t26;
    const double t198 = t87 * t93;
    const double t200 = t5 * t198 * t48;
    const double t203 = t5 * t88 * rpbe_mu;
    const double t207 = 0.1e1 / t92 / t6;
    const double t208 = t25 * t207;
    const double t211 = t5 * t208 * t48 / 0.12e2;
    const double t213 = t5 * t94 * rpbe_mu;
    const double t214 = t213 * t106;
    const double t216 = t35 * t35;
    const double t218 = 0.1e1 / t37 / t216;
    const double t221 = t100 * sigma_aa * t218 * t45;
    const double t226 = t5 * t27 * t224;
    const double t231 = sigma_aa * sigma_aa;
    const double t232 = t230 * t231;
    const double t235 = 0.1e1 / t36 / t216 / t101;
    const double t237 = t235 * t41 * t45;
    const double t238 = t232 * t237;
    const double t242 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t194 * t48 - t200 / 0.4e1 + t203 * t106 / 0.12e2 + t211 + t214 / 0.36e2 - 0.11e2 / 0.72e2 * t99 * t221 + t226 * t238 / 0.216e3 );
    const double t243 = t59 * t59;
    const double t244 = 0.1e1 / t243;
    const double t245 = t113 * t113;
    const double t248 = t54 * t185;
    const double t251 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t81 + 0.2e1 * t248 );
    const double t255 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t244 * t245 + 0.4e1 / 0.3e1 * t59 * t251 );
    const double t256 = t255 * t26;
    const double t260 = t116 * t93;
    const double t262 = t5 * t260 * t75;
    const double t264 = t61 * t207;
    const double t267 = t5 * t264 * t75 / 0.12e2;
    const double t269 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t256 * t75 - t262 / 0.4e1 + t267 );
    const double t285 = t133 * t93;
    const double t287 = t5 * t285 * t48;
    const double t310 = t144 * t93;
    const double t312 = t5 * t310 * t75;
    const double t320 = t5 * t121 * rpbe_mu;
    const double t321 = t320 * t156;
    const double t329 = t130 * t130;
    const double t334 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t81 + 0.2e1 * t186 );
    const double t338 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t180 * t329 + 0.4e1 / 0.3e1 * t23 * t334 );
    const double t339 = t338 * t26;
    const double t345 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t339 * t48 - t287 / 0.4e1 + t211 );
    const double t346 = t141 * t141;
    const double t351 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t81 + 0.2e1 * t248 );
    const double t355 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t244 * t346 + 0.4e1 / 0.3e1 * t59 * t351 );
    const double t356 = t355 * t26;
    const double t362 = t5 * t145 * rpbe_mu;
    const double t366 = t63 * t63;
    const double t368 = 0.1e1 / t65 / t366;
    const double t371 = t100 * sigma_bb * t368 * t72;
    const double t375 = t5 * t62 * t224;
    const double t376 = sigma_bb * sigma_bb;
    const double t377 = t230 * t376;
    const double t380 = 0.1e1 / t64 / t366 / t151;
    const double t382 = t380 * t41 * t72;
    const double t383 = t377 * t382;
    const double t387 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t356 * t75 - t312 / 0.4e1 + t362 * t156 / 0.12e2 + t267 + t321 / 0.36e2 - 0.11e2 / 0.72e2 * t150 * t371 + t375 * t383 / 0.216e3 );
    const double t390 = t5 * t88;
    const double t393 = t5 * t94;
    const double t395 = t393 * t166 / 0.192e3;
    const double t398 = t29 * t33 * t103 * t45;
    const double t401 = t216 * t35;
    const double t403 = 0.1e1 / t36 / t401;
    const double t406 = sigma_aa * t41 * t45;
    const double t407 = t230 * t403 * t406;
    const double t411 = piecewise_functor_3( t1, 0.0, -t390 * t166 / 0.64e2 - t395 + t163 * t398 / 0.24e2 - t226 * t407 / 0.576e3 );
    const double t416 = t5 * t121;
    const double t418 = t416 * t173 / 0.192e3;
    const double t428 = t5 * t145;
    const double t433 = t29 * t33 * t153 * t72;
    const double t436 = t366 * t63;
    const double t438 = 0.1e1 / t64 / t436;
    const double t441 = sigma_bb * t41 * t72;
    const double t442 = t230 * t438 * t441;
    const double t446 = piecewise_functor_3( t53, 0.0, -t428 * t173 / 0.64e2 - t418 + t170 * t433 / 0.24e2 - t375 * t442 / 0.576e3 );
    const double t448 = t216 * rho_a;
    const double t453 = t230 / t36 / t448 * t41 * t45;
    const double t456 = piecewise_functor_3( t1, 0.0, t226 * t453 / 0.1536e4 );
    const double t457 = t366 * rho_b;
    const double t462 = t230 / t64 / t457 * t41 * t72;
    const double t465 = piecewise_functor_3( t53, 0.0, t375 * t462 / 0.1536e4 );


    vrho_a = t52 + t79 + t6 * ( t110 + t126 );
    vrho_b = t52 + t79 + t6 * ( t139 + t160 );
    vsigma_aa = t6 * t169;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t176;
    v2rho2_aa = 0.2e1 * t110 + 0.2e1 * t126 + t6 * ( t242 + t269 );
    v2rho2_bb = 0.2e1 * t139 + 0.2e1 * t160 + t6 * ( t345 + t387 );
    v2rhosigma_a_aa = t6 * t411 + t169;
    v2rhosigma_b_bb = t6 * t446 + t176;
    v2sigma2_aa_aa = t6 * t456;
    v2sigma2_bb_bb = t6 * t465;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;



  }


};

struct BuiltinRPBE_X : detail::BuiltinKernelImpl< BuiltinRPBE_X > {

  BuiltinRPBE_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinRPBE_X >(p) { }
  
  virtual ~BuiltinRPBE_X() = default;

};



} // namespace ExchCXX
