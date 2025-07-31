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

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinRevPBE_X > :
  public gga_screening_interface< BuiltinRevPBE_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-32;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 2.1544346900318956e-43;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;


  static constexpr double kappa = 1.245;
  static constexpr double mu =  0.2195149727645171;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
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
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t56 = kappa * kappa;
    constexpr double t78 = t20 * t25 * t28;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );
    const double t52 = t30 * rho;
    const double t58 = t6 * t17 / t18 / t52 * t56;
    const double t59 = t37 * t37;
    const double t61 = 0.1e1 / t59 * mu;
    const double t64 = t25 * sigma * t28;
    const double t65 = t61 * t20 * t64;
    const double t69 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t58 * t65 / 0.24e2 );
    const double t79 = t61 * t78;
    const double t82 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t56 * t79 / 0.64e2 );


    eps = 0.2e1 * t46;
    vrho = 0.2e1 * rho * t69 + 0.2e1 * t46;
    vsigma = 0.2e1 * rho * t82;

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
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t56 = kappa * kappa;
    constexpr double t78 = t20 * t25 * t28;
    constexpr double t106 = mu * mu;
    constexpr double t108 = t20 * t20;
    constexpr double t111 = 0.1e1 / t23 / t22;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
    const double t52 = t30 * rho;
    const double t58 = t6 * t17 / t18 / t52 * t56;
    const double t59 = t37 * t37;
    const double t61 = 0.1e1 / t59 * mu;
    const double t64 = t25 * sigma * t28;
    const double t65 = t61 * t20 * t64;
    const double t69 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t58 * t65 / 0.24e2 );
    const double t79 = t61 * t78;
    const double t82 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t56 * t79 / 0.64e2 );
    const double t91 = t30 * t30;
    const double t96 = t6 * t17 / t18 / t91 * t56;
    const double t99 = t91 * t52;
    const double t103 = t6 * t17 / t99 * t56;
    const double t107 = 0.1e1 / t59 / t37 * t106;
    const double t109 = t107 * t108;
    const double t112 = sigma * sigma;
    const double t115 = t109 * t111 * t112 * t27;
    const double t119 = piecewise_functor_3( t2, 0.0, t6 * t17 / t31 / rho * t42 / 0.12e2 - t96 * t65 / 0.8e1 + t103 * t115 / 0.54e2 );
    const double t124 = t91 * t30;
    const double t128 = t6 * t17 / t124 * t56;
    const double t131 = t109 * t111 * t27 * sigma;
    const double t135 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.192e3 * t58 * t79 - t128 * t131 / 0.144e3 );
    const double t138 = t91 * rho;
    const double t145 = t107 * t108 * t111 * t27;
    const double t148 = piecewise_functor_3( t2, 0.0, t6 * t17 / t138 * t56 * t145 / 0.384e3 );


    v2rho2 = 0.2e1 * rho * t119 + 0.4e1 * t69;
    v2rhosigma = 0.2e1 * rho * t135 + 0.2e1 * t82;
    v2sigma2 = 0.2e1 * rho * t148;

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
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t56 = kappa * kappa;
    constexpr double t78 = t20 * t25 * t28;
    constexpr double t106 = mu * mu;
    constexpr double t108 = t20 * t20;
    constexpr double t111 = 0.1e1 / t23 / t22;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );
    const double t52 = t30 * rho;
    const double t58 = t6 * t17 / t18 / t52 * t56;
    const double t59 = t37 * t37;
    const double t61 = 0.1e1 / t59 * mu;
    const double t64 = t25 * sigma * t28;
    const double t65 = t61 * t20 * t64;
    const double t69 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t58 * t65 / 0.24e2 );
    const double t79 = t61 * t78;
    const double t82 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t56 * t79 / 0.64e2 );
    const double t91 = t30 * t30;
    const double t96 = t6 * t17 / t18 / t91 * t56;
    const double t99 = t91 * t52;
    const double t103 = t6 * t17 / t99 * t56;
    const double t107 = 0.1e1 / t59 / t37 * t106;
    const double t109 = t107 * t108;
    const double t112 = sigma * sigma;
    const double t115 = t109 * t111 * t112 * t27;
    const double t119 = piecewise_functor_3( t2, 0.0, t6 * t17 / t31 / rho * t42 / 0.12e2 - t96 * t65 / 0.8e1 + t103 * t115 / 0.54e2 );
    const double t124 = t91 * t30;
    const double t128 = t6 * t17 / t124 * t56;
    const double t131 = t109 * t111 * t27 * sigma;
    const double t135 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.192e3 * t58 * t79 - t128 * t131 / 0.144e3 );
    const double t138 = t91 * rho;
    const double t145 = t107 * t108 * t111 * t27;
    const double t148 = piecewise_functor_3( t2, 0.0, t6 * t17 / t138 * t56 * t145 / 0.384e3 );


    vrho = 0.2e1 * rho * t69 + 0.2e1 * t46;
    vsigma = 0.2e1 * rho * t82;
    v2rho2 = 0.2e1 * rho * t119 + 0.4e1 * t69;
    v2rhosigma = 0.2e1 * rho * t135 + 0.2e1 * t82;
    v2sigma2 = 0.2e1 * rho * t148;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t81 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t77 );


    eps = t52 + t81;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t100 = kappa * kappa;
    constexpr double t171 = t28 * t33;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t81 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t77 );
    const double t82 = t6 * t6;
    const double t83 = 0.1e1 / t82;
    const double t84 = t16 * t83;
    const double t86 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t84 );
    const double t89 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t86 );
    const double t90 = t89 * t26;
    const double t94 = t26 * t26;
    const double t95 = 0.1e1 / t94;
    const double t96 = t25 * t95;
    const double t99 = t5 * t96 * t48 / 0.8e1;
    const double t101 = t27 * t100;
    const double t102 = t5 * t101;
    const double t103 = t43 * t43;
    const double t105 = 0.1e1 / t103 * mu;
    const double t106 = t105 * t28;
    const double t107 = t35 * rho_a;
    const double t109 = 0.1e1 / t37 / t107;
    const double t111 = t106 * t34 * t109;
    const double t115 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t90 * t48 - t99 + t102 * t111 / 0.24e2 );
    const double t116 = t54 * t83;
    const double t118 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t116 );
    const double t121 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t118 );
    const double t122 = t121 * t26;
    const double t126 = t61 * t95;
    const double t129 = t5 * t126 * t77 / 0.8e1;
    const double t131 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t122 * t77 - t129 );
    const double t135 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t84 );
    const double t138 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t135 );
    const double t139 = t138 * t26;
    const double t144 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t139 * t48 - t99 );
    const double t146 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t116 );
    const double t149 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t146 );
    const double t150 = t149 * t26;
    const double t154 = t62 * t100;
    const double t155 = t5 * t154;
    const double t156 = t72 * t72;
    const double t158 = 0.1e1 / t156 * mu;
    const double t159 = t158 * t28;
    const double t160 = t64 * rho_b;
    const double t162 = 0.1e1 / t66 / t160;
    const double t164 = t159 * t63 * t162;
    const double t168 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t150 * t77 - t129 + t155 * t164 / 0.24e2 );
    const double t173 = t105 * t171 * t39;
    const double t176 = piecewise_functor_3( t1, 0.0, -t102 * t173 / 0.64e2 );
    const double t178 = t158 * t171 * t68;
    const double t181 = piecewise_functor_3( t53, 0.0, -t155 * t178 / 0.64e2 );


    eps = t52 + t81;
    vrho_a = t52 + t81 + t6 * ( t115 + t131 );
    vrho_b = t52 + t81 + t6 * ( t144 + t168 );
    vsigma_aa = t6 * t176;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t181;

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
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t100 = kappa * kappa;
    constexpr double t171 = t28 * t33;
    constexpr double t223 = mu * mu;
    constexpr double t225 = t28 * t28;
    constexpr double t228 = 0.1e1 / t31 / t30;
    constexpr double t442 = t225 * t228;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t82 = t6 * t6;
    const double t83 = 0.1e1 / t82;
    const double t84 = t16 * t83;
    const double t86 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t84 );
    const double t89 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t86 );
    const double t90 = t89 * t26;
    const double t94 = t26 * t26;
    const double t95 = 0.1e1 / t94;
    const double t96 = t25 * t95;
    const double t99 = t5 * t96 * t48 / 0.8e1;
    const double t101 = t27 * t100;
    const double t102 = t5 * t101;
    const double t103 = t43 * t43;
    const double t105 = 0.1e1 / t103 * mu;
    const double t106 = t105 * t28;
    const double t107 = t35 * rho_a;
    const double t109 = 0.1e1 / t37 / t107;
    const double t111 = t106 * t34 * t109;
    const double t115 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t90 * t48 - t99 + t102 * t111 / 0.24e2 );
    const double t116 = t54 * t83;
    const double t118 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t116 );
    const double t121 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t118 );
    const double t122 = t121 * t26;
    const double t126 = t61 * t95;
    const double t129 = t5 * t126 * t77 / 0.8e1;
    const double t131 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t122 * t77 - t129 );
    const double t135 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t84 );
    const double t138 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t135 );
    const double t139 = t138 * t26;
    const double t144 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t139 * t48 - t99 );
    const double t146 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t116 );
    const double t149 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t146 );
    const double t150 = t149 * t26;
    const double t154 = t62 * t100;
    const double t155 = t5 * t154;
    const double t156 = t72 * t72;
    const double t158 = 0.1e1 / t156 * mu;
    const double t159 = t158 * t28;
    const double t160 = t64 * rho_b;
    const double t162 = 0.1e1 / t66 / t160;
    const double t164 = t159 * t63 * t162;
    const double t168 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t150 * t77 - t129 + t155 * t164 / 0.24e2 );
    const double t173 = t105 * t171 * t39;
    const double t176 = piecewise_functor_3( t1, 0.0, -t102 * t173 / 0.64e2 );
    const double t178 = t158 * t171 * t68;
    const double t181 = piecewise_functor_3( t53, 0.0, -t155 * t178 / 0.64e2 );
    const double t184 = t23 * t23;
    const double t185 = 0.1e1 / t184;
    const double t186 = t86 * t86;
    const double t189 = t82 * t6;
    const double t190 = 0.1e1 / t189;
    const double t191 = t16 * t190;
    const double t194 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t83 + 0.2e1 * t191 );
    const double t198 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t185 * t186 + 0.4e1 / 0.3e1 * t23 * t194 );
    const double t199 = t198 * t26;
    const double t203 = t89 * t95;
    const double t205 = t5 * t203 * t48;
    const double t208 = t5 * t90 * t100;
    const double t212 = 0.1e1 / t94 / t6;
    const double t213 = t25 * t212;
    const double t216 = t5 * t213 * t48 / 0.12e2;
    const double t218 = t5 * t96 * t100;
    const double t219 = t218 * t111;
    const double t224 = 0.1e1 / t103 / t43 * t223;
    const double t226 = t224 * t225;
    const double t229 = sigma_aa * sigma_aa;
    const double t230 = t228 * t229;
    const double t231 = t35 * t35;
    const double t234 = 0.1e1 / t36 / t231 / t107;
    const double t236 = t226 * t230 * t234;
    const double t240 = 0.1e1 / t37 / t231;
    const double t242 = t106 * t34 * t240;
    const double t246 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t199 * t48 - t205 / 0.4e1 + t208 * t111 / 0.12e2 + t216 + t219 / 0.36e2 + t102 * t236 / 0.108e3 - 0.11e2 / 0.72e2 * t102 * t242 );
    const double t247 = t59 * t59;
    const double t248 = 0.1e1 / t247;
    const double t249 = t118 * t118;
    const double t252 = t54 * t190;
    const double t255 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t83 + 0.2e1 * t252 );
    const double t259 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t248 * t249 + 0.4e1 / 0.3e1 * t59 * t255 );
    const double t260 = t259 * t26;
    const double t264 = t121 * t95;
    const double t266 = t5 * t264 * t77;
    const double t268 = t61 * t212;
    const double t271 = t5 * t268 * t77 / 0.12e2;
    const double t273 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t260 * t77 - t266 / 0.4e1 + t271 );
    const double t289 = t138 * t95;
    const double t291 = t5 * t289 * t48;
    const double t314 = t149 * t95;
    const double t316 = t5 * t314 * t77;
    const double t324 = t5 * t126 * t100;
    const double t325 = t324 * t164;
    const double t333 = t135 * t135;
    const double t338 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t83 + 0.2e1 * t191 );
    const double t342 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t185 * t333 + 0.4e1 / 0.3e1 * t23 * t338 );
    const double t343 = t342 * t26;
    const double t349 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t343 * t48 - t291 / 0.4e1 + t216 );
    const double t350 = t146 * t146;
    const double t355 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t83 + 0.2e1 * t252 );
    const double t359 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t248 * t350 + 0.4e1 / 0.3e1 * t59 * t355 );
    const double t360 = t359 * t26;
    const double t366 = t5 * t150 * t100;
    const double t372 = 0.1e1 / t156 / t72 * t223;
    const double t373 = t372 * t225;
    const double t374 = sigma_bb * sigma_bb;
    const double t375 = t228 * t374;
    const double t376 = t64 * t64;
    const double t379 = 0.1e1 / t65 / t376 / t160;
    const double t381 = t373 * t375 * t379;
    const double t385 = 0.1e1 / t66 / t376;
    const double t387 = t159 * t63 * t385;
    const double t391 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t360 * t77 - t316 / 0.4e1 + t366 * t164 / 0.12e2 + t271 + t325 / 0.36e2 + t155 * t381 / 0.108e3 - 0.11e2 / 0.72e2 * t155 * t387 );
    const double t397 = t218 * t173 / 0.192e3;
    const double t398 = t231 * t35;
    const double t400 = 0.1e1 / t36 / t398;
    const double t403 = t226 * t228 * t400 * sigma_aa;
    const double t407 = t105 * t171 * t109;
    const double t411 = piecewise_functor_3( t1, 0.0, -t208 * t173 / 0.64e2 - t397 - t102 * t403 / 0.288e3 + t102 * t407 / 0.24e2 );
    const double t416 = t324 * t178 / 0.192e3;
    const double t427 = t376 * t64;
    const double t429 = 0.1e1 / t65 / t427;
    const double t432 = t373 * t228 * t429 * sigma_bb;
    const double t436 = t158 * t171 * t162;
    const double t440 = piecewise_functor_3( t53, 0.0, -t366 * t178 / 0.64e2 - t416 - t155 * t432 / 0.288e3 + t155 * t436 / 0.24e2 );
    const double t443 = t231 * rho_a;
    const double t447 = t224 * t442 / t36 / t443;
    const double t450 = piecewise_functor_3( t1, 0.0, t102 * t447 / 0.768e3 );
    const double t451 = t376 * rho_b;
    const double t455 = t372 * t442 / t65 / t451;
    const double t458 = piecewise_functor_3( t53, 0.0, t155 * t455 / 0.768e3 );


    v2rho2_aa = 0.2e1 * t115 + 0.2e1 * t131 + t6 * ( t246 + t273 );
    v2rho2_bb = 0.2e1 * t144 + 0.2e1 * t168 + t6 * ( t349 + t391 );
    v2rhosigma_a_aa = t6 * t411 + t176;
    v2rhosigma_b_bb = t6 * t440 + t181;
    v2sigma2_aa_aa = t6 * t450;
    v2sigma2_bb_bb = t6 * t458;
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
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t100 = kappa * kappa;
    constexpr double t171 = t28 * t33;
    constexpr double t223 = mu * mu;
    constexpr double t225 = t28 * t28;
    constexpr double t228 = 0.1e1 / t31 / t30;
    constexpr double t442 = t225 * t228;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t81 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t77 );
    const double t82 = t6 * t6;
    const double t83 = 0.1e1 / t82;
    const double t84 = t16 * t83;
    const double t86 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t84 );
    const double t89 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t86 );
    const double t90 = t89 * t26;
    const double t94 = t26 * t26;
    const double t95 = 0.1e1 / t94;
    const double t96 = t25 * t95;
    const double t99 = t5 * t96 * t48 / 0.8e1;
    const double t101 = t27 * t100;
    const double t102 = t5 * t101;
    const double t103 = t43 * t43;
    const double t105 = 0.1e1 / t103 * mu;
    const double t106 = t105 * t28;
    const double t107 = t35 * rho_a;
    const double t109 = 0.1e1 / t37 / t107;
    const double t111 = t106 * t34 * t109;
    const double t115 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t90 * t48 - t99 + t102 * t111 / 0.24e2 );
    const double t116 = t54 * t83;
    const double t118 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t116 );
    const double t121 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t118 );
    const double t122 = t121 * t26;
    const double t126 = t61 * t95;
    const double t129 = t5 * t126 * t77 / 0.8e1;
    const double t131 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t122 * t77 - t129 );
    const double t135 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t84 );
    const double t138 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t135 );
    const double t139 = t138 * t26;
    const double t144 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t139 * t48 - t99 );
    const double t146 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t116 );
    const double t149 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t146 );
    const double t150 = t149 * t26;
    const double t154 = t62 * t100;
    const double t155 = t5 * t154;
    const double t156 = t72 * t72;
    const double t158 = 0.1e1 / t156 * mu;
    const double t159 = t158 * t28;
    const double t160 = t64 * rho_b;
    const double t162 = 0.1e1 / t66 / t160;
    const double t164 = t159 * t63 * t162;
    const double t168 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t150 * t77 - t129 + t155 * t164 / 0.24e2 );
    const double t173 = t105 * t171 * t39;
    const double t176 = piecewise_functor_3( t1, 0.0, -t102 * t173 / 0.64e2 );
    const double t178 = t158 * t171 * t68;
    const double t181 = piecewise_functor_3( t53, 0.0, -t155 * t178 / 0.64e2 );
    const double t184 = t23 * t23;
    const double t185 = 0.1e1 / t184;
    const double t186 = t86 * t86;
    const double t189 = t82 * t6;
    const double t190 = 0.1e1 / t189;
    const double t191 = t16 * t190;
    const double t194 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t83 + 0.2e1 * t191 );
    const double t198 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t185 * t186 + 0.4e1 / 0.3e1 * t23 * t194 );
    const double t199 = t198 * t26;
    const double t203 = t89 * t95;
    const double t205 = t5 * t203 * t48;
    const double t208 = t5 * t90 * t100;
    const double t212 = 0.1e1 / t94 / t6;
    const double t213 = t25 * t212;
    const double t216 = t5 * t213 * t48 / 0.12e2;
    const double t218 = t5 * t96 * t100;
    const double t219 = t218 * t111;
    const double t224 = 0.1e1 / t103 / t43 * t223;
    const double t226 = t224 * t225;
    const double t229 = sigma_aa * sigma_aa;
    const double t230 = t228 * t229;
    const double t231 = t35 * t35;
    const double t234 = 0.1e1 / t36 / t231 / t107;
    const double t236 = t226 * t230 * t234;
    const double t240 = 0.1e1 / t37 / t231;
    const double t242 = t106 * t34 * t240;
    const double t246 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t199 * t48 - t205 / 0.4e1 + t208 * t111 / 0.12e2 + t216 + t219 / 0.36e2 + t102 * t236 / 0.108e3 - 0.11e2 / 0.72e2 * t102 * t242 );
    const double t247 = t59 * t59;
    const double t248 = 0.1e1 / t247;
    const double t249 = t118 * t118;
    const double t252 = t54 * t190;
    const double t255 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t83 + 0.2e1 * t252 );
    const double t259 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t248 * t249 + 0.4e1 / 0.3e1 * t59 * t255 );
    const double t260 = t259 * t26;
    const double t264 = t121 * t95;
    const double t266 = t5 * t264 * t77;
    const double t268 = t61 * t212;
    const double t271 = t5 * t268 * t77 / 0.12e2;
    const double t273 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t260 * t77 - t266 / 0.4e1 + t271 );
    const double t289 = t138 * t95;
    const double t291 = t5 * t289 * t48;
    const double t314 = t149 * t95;
    const double t316 = t5 * t314 * t77;
    const double t324 = t5 * t126 * t100;
    const double t325 = t324 * t164;
    const double t333 = t135 * t135;
    const double t338 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t83 + 0.2e1 * t191 );
    const double t342 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t185 * t333 + 0.4e1 / 0.3e1 * t23 * t338 );
    const double t343 = t342 * t26;
    const double t349 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t343 * t48 - t291 / 0.4e1 + t216 );
    const double t350 = t146 * t146;
    const double t355 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t83 + 0.2e1 * t252 );
    const double t359 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t248 * t350 + 0.4e1 / 0.3e1 * t59 * t355 );
    const double t360 = t359 * t26;
    const double t366 = t5 * t150 * t100;
    const double t372 = 0.1e1 / t156 / t72 * t223;
    const double t373 = t372 * t225;
    const double t374 = sigma_bb * sigma_bb;
    const double t375 = t228 * t374;
    const double t376 = t64 * t64;
    const double t379 = 0.1e1 / t65 / t376 / t160;
    const double t381 = t373 * t375 * t379;
    const double t385 = 0.1e1 / t66 / t376;
    const double t387 = t159 * t63 * t385;
    const double t391 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t360 * t77 - t316 / 0.4e1 + t366 * t164 / 0.12e2 + t271 + t325 / 0.36e2 + t155 * t381 / 0.108e3 - 0.11e2 / 0.72e2 * t155 * t387 );
    const double t397 = t218 * t173 / 0.192e3;
    const double t398 = t231 * t35;
    const double t400 = 0.1e1 / t36 / t398;
    const double t403 = t226 * t228 * t400 * sigma_aa;
    const double t407 = t105 * t171 * t109;
    const double t411 = piecewise_functor_3( t1, 0.0, -t208 * t173 / 0.64e2 - t397 - t102 * t403 / 0.288e3 + t102 * t407 / 0.24e2 );
    const double t416 = t324 * t178 / 0.192e3;
    const double t427 = t376 * t64;
    const double t429 = 0.1e1 / t65 / t427;
    const double t432 = t373 * t228 * t429 * sigma_bb;
    const double t436 = t158 * t171 * t162;
    const double t440 = piecewise_functor_3( t53, 0.0, -t366 * t178 / 0.64e2 - t416 - t155 * t432 / 0.288e3 + t155 * t436 / 0.24e2 );
    const double t443 = t231 * rho_a;
    const double t447 = t224 * t442 / t36 / t443;
    const double t450 = piecewise_functor_3( t1, 0.0, t102 * t447 / 0.768e3 );
    const double t451 = t376 * rho_b;
    const double t455 = t372 * t442 / t65 / t451;
    const double t458 = piecewise_functor_3( t53, 0.0, t155 * t455 / 0.768e3 );


    vrho_a = t52 + t81 + t6 * ( t115 + t131 );
    vrho_b = t52 + t81 + t6 * ( t144 + t168 );
    vsigma_aa = t6 * t176;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t181;
    v2rho2_aa = 0.2e1 * t115 + 0.2e1 * t131 + t6 * ( t246 + t273 );
    v2rho2_bb = 0.2e1 * t144 + 0.2e1 * t168 + t6 * ( t349 + t391 );
    v2rhosigma_a_aa = t6 * t411 + t176;
    v2rhosigma_b_bb = t6 * t440 + t181;
    v2sigma2_aa_aa = t6 * t450;
    v2sigma2_bb_bb = t6 * t458;
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

struct BuiltinRevPBE_X : detail::BuiltinKernelImpl< BuiltinRevPBE_X > {

  BuiltinRevPBE_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinRevPBE_X >(p) { }
  
  virtual ~BuiltinRevPBE_X() = default;

};



} // namespace ExchCXX
