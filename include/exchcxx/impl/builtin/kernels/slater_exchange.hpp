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
struct kernel_traits< BuiltinSlaterExchange > :
  public lda_screening_interface< BuiltinSlaterExchange > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;


  static constexpr double alpha = 1.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t15 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t10 * t11 );
    const double t16 = alpha * t15;


    eps = 0.2e1 * t16;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t15 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t10 * t11 );
    const double t16 = alpha * t15;
    const double t17 = rho * alpha;
    const double t18 = t11 * t11;
    const double t23 = piecewise_functor_3( t2, 0.0, -t6 * t10 / t18 / 0.8e1 );


    eps = 0.2e1 * t16;
    vrho = 0.2e1 * t17 * t23 + 0.2e1 * t16;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double& v2rho2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t17 = rho * alpha;
    const double t18 = t11 * t11;
    const double t23 = piecewise_functor_3( t2, 0.0, -t6 * t10 / t18 / 0.8e1 );
    const double t33 = piecewise_functor_3( t2, 0.0, t6 * t10 / t18 / rho / 0.12e2 );


    v2rho2 = 0.2e1 * t17 * t33 + 0.4e1 * alpha * t23;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double& vrho, double& v2rho2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t15 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t10 * t11 );
    const double t16 = alpha * t15;
    const double t17 = rho * alpha;
    const double t18 = t11 * t11;
    const double t23 = piecewise_functor_3( t2, 0.0, -t6 * t10 / t18 / 0.8e1 );
    const double t33 = piecewise_functor_3( t2, 0.0, t6 * t10 / t18 / rho / 0.12e2 );


    vrho = 0.2e1 * t17 * t23 + 0.2e1 * t16;
    v2rho2 = 0.2e1 * t17 * t33 + 0.4e1 * alpha * t23;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t24 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t19 * t20 );
    const double t25 = alpha * t24;
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t39 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t35 * t20 );
    const double t40 = alpha * t39;


    eps = t25 + t40;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t24 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t19 * t20 );
    const double t25 = alpha * t24;
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t39 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t35 * t20 );
    const double t40 = alpha * t39;
    const double t41 = t13 * t7;
    const double t44 = t6 * t6;
    const double t45 = 0.1e1 / t44;
    const double t48 = 0.2e1 * t14 * t45 * t15;
    const double t49 = t15 * t15;
    const double t50 = 0.1e1 / t49;
    const double t51 = t7 * t50;
    const double t53 = -rho_a * t45 + t7;
    const double t58 = piecewise_functor_3( t10, 0.0, 0.2e1 * t41 * t15 - t48 + 0.2e1 / 0.3e1 * t14 * t51 * t53 );
    const double t62 = t20 * t20;
    const double t63 = 0.1e1 / t62;
    const double t66 = t5 * t19 * t63 / 0.8e1;
    const double t68 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t58 * t20 - t66 );
    const double t69 = alpha * t68;
    const double t72 = 0.2e1 * t30 * t45 * t31;
    const double t73 = rho_b * rho_b;
    const double t74 = t13 * t73;
    const double t75 = t44 * t6;
    const double t76 = 0.1e1 / t75;
    const double t77 = t31 * t31;
    const double t78 = 0.1e1 / t77;
    const double t79 = t76 * t78;
    const double t83 = piecewise_functor_3( t29, 0.0, -t72 - 0.2e1 / 0.3e1 * t74 * t79 );
    const double t89 = t5 * t35 * t63 / 0.8e1;
    const double t91 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t83 * t20 - t89 );
    const double t92 = alpha * t91;
    const double t95 = rho_a * rho_a;
    const double t96 = t13 * t95;
    const double t97 = t76 * t50;
    const double t101 = piecewise_functor_3( t10, 0.0, -t48 - 0.2e1 / 0.3e1 * t96 * t97 );
    const double t106 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t101 * t20 - t66 );
    const double t107 = alpha * t106;
    const double t110 = t7 * t78;
    const double t112 = -rho_b * t45 + t7;
    const double t117 = piecewise_functor_3( t29, 0.0, 0.2e1 * t41 * t31 - t72 + 0.2e1 / 0.3e1 * t30 * t110 * t112 );
    const double t122 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t20 - t89 );
    const double t123 = alpha * t122;


    eps = t25 + t40;
    vrho_a = t25 + t40 + t6 * ( t69 + t92 );
    vrho_b = t25 + t40 + t6 * ( t107 + t123 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t41 = t13 * t7;
    const double t44 = t6 * t6;
    const double t45 = 0.1e1 / t44;
    const double t48 = 0.2e1 * t14 * t45 * t15;
    const double t49 = t15 * t15;
    const double t50 = 0.1e1 / t49;
    const double t51 = t7 * t50;
    const double t53 = -rho_a * t45 + t7;
    const double t58 = piecewise_functor_3( t10, 0.0, 0.2e1 * t41 * t15 - t48 + 0.2e1 / 0.3e1 * t14 * t51 * t53 );
    const double t62 = t20 * t20;
    const double t63 = 0.1e1 / t62;
    const double t66 = t5 * t19 * t63 / 0.8e1;
    const double t68 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t58 * t20 - t66 );
    const double t69 = alpha * t68;
    const double t72 = 0.2e1 * t30 * t45 * t31;
    const double t73 = rho_b * rho_b;
    const double t74 = t13 * t73;
    const double t75 = t44 * t6;
    const double t76 = 0.1e1 / t75;
    const double t77 = t31 * t31;
    const double t78 = 0.1e1 / t77;
    const double t79 = t76 * t78;
    const double t83 = piecewise_functor_3( t29, 0.0, -t72 - 0.2e1 / 0.3e1 * t74 * t79 );
    const double t89 = t5 * t35 * t63 / 0.8e1;
    const double t91 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t83 * t20 - t89 );
    const double t92 = alpha * t91;
    const double t95 = rho_a * rho_a;
    const double t96 = t13 * t95;
    const double t97 = t76 * t50;
    const double t101 = piecewise_functor_3( t10, 0.0, -t48 - 0.2e1 / 0.3e1 * t96 * t97 );
    const double t106 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t101 * t20 - t66 );
    const double t107 = alpha * t106;
    const double t110 = t7 * t78;
    const double t112 = -rho_b * t45 + t7;
    const double t117 = piecewise_functor_3( t29, 0.0, 0.2e1 * t41 * t31 - t72 + 0.2e1 / 0.3e1 * t30 * t110 * t112 );
    const double t122 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t20 - t89 );
    const double t123 = alpha * t122;
    const double t128 = t13 * t45;
    const double t129 = t128 * t15;
    const double t131 = t50 * t53;
    const double t136 = 0.4e1 * t14 * t76 * t15;
    const double t137 = t45 * t50;
    const double t139 = t14 * t137 * t53;
    const double t142 = 0.1e1 / t49 / t8;
    const double t143 = t7 * t142;
    const double t144 = t53 * t53;
    const double t150 = 0.2e1 * rho_a * t76 - 0.2e1 * t45;
    const double t155 = piecewise_functor_3( t10, 0.0, -0.4e1 * t129 + 0.4e1 / 0.3e1 * t41 * t131 + t136 - 0.4e1 / 0.3e1 * t139 - 0.4e1 / 0.9e1 * t14 * t143 * t144 + 0.2e1 / 0.3e1 * t14 * t51 * t150 );
    const double t160 = t5 * t58 * t63;
    const double t163 = 0.1e1 / t62 / t6;
    const double t166 = t5 * t19 * t163 / 0.12e2;
    const double t168 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t155 * t20 - t160 / 0.4e1 + t166 );
    const double t169 = alpha * t168;
    const double t172 = 0.4e1 * t30 * t76 * t31;
    const double t173 = t44 * t44;
    const double t174 = 0.1e1 / t173;
    const double t175 = t174 * t78;
    const double t176 = t74 * t175;
    const double t178 = t73 * rho_b;
    const double t179 = t13 * t178;
    const double t181 = 0.1e1 / t173 / t6;
    const double t183 = 0.1e1 / t77 / t27;
    const double t184 = t181 * t183;
    const double t188 = piecewise_functor_3( t29, 0.0, t172 + 0.8e1 / 0.3e1 * t176 - 0.4e1 / 0.9e1 * t179 * t184 );
    const double t193 = t5 * t83 * t63;
    const double t197 = t5 * t35 * t163 / 0.12e2;
    const double t199 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t188 * t20 - t193 / 0.4e1 + t197 );
    const double t200 = alpha * t199;
    const double t207 = t174 * t50;
    const double t208 = t96 * t207;
    const double t220 = t5 * t101 * t63;
    const double t226 = t128 * t31;
    const double t233 = t45 * t78;
    const double t235 = t30 * t233 * t112;
    const double t241 = rho_b * t76;
    const double t253 = t5 * t117 * t63;
    const double t264 = t95 * rho_a;
    const double t265 = t13 * t264;
    const double t266 = t181 * t142;
    const double t270 = piecewise_functor_3( t10, 0.0, t136 + 0.8e1 / 0.3e1 * t208 - 0.4e1 / 0.9e1 * t265 * t266 );
    const double t276 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t270 * t20 - t220 / 0.4e1 + t166 );
    const double t277 = alpha * t276;
    const double t279 = t78 * t112;
    const double t283 = t7 * t183;
    const double t284 = t112 * t112;
    const double t289 = -0.2e1 * t45 + 0.2e1 * t241;
    const double t294 = piecewise_functor_3( t29, 0.0, -0.4e1 * t226 + 0.4e1 / 0.3e1 * t41 * t279 + t172 - 0.4e1 / 0.3e1 * t235 - 0.4e1 / 0.9e1 * t30 * t283 * t284 + 0.2e1 / 0.3e1 * t30 * t110 * t289 );
    const double t300 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t294 * t20 - t253 / 0.4e1 + t197 );
    const double t301 = alpha * t300;


    v2rho2_aa = 0.2e1 * t69 + 0.2e1 * t92 + t6 * ( t169 + t200 );
    v2rho2_bb = 0.2e1 * t107 + 0.2e1 * t123 + t6 * ( t277 + t301 );
    v2rho2_ab = 0.0;



  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t24 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t19 * t20 );
    const double t25 = alpha * t24;
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t39 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t35 * t20 );
    const double t40 = alpha * t39;
    const double t41 = t13 * t7;
    const double t44 = t6 * t6;
    const double t45 = 0.1e1 / t44;
    const double t48 = 0.2e1 * t14 * t45 * t15;
    const double t49 = t15 * t15;
    const double t50 = 0.1e1 / t49;
    const double t51 = t7 * t50;
    const double t53 = -rho_a * t45 + t7;
    const double t58 = piecewise_functor_3( t10, 0.0, 0.2e1 * t41 * t15 - t48 + 0.2e1 / 0.3e1 * t14 * t51 * t53 );
    const double t62 = t20 * t20;
    const double t63 = 0.1e1 / t62;
    const double t66 = t5 * t19 * t63 / 0.8e1;
    const double t68 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t58 * t20 - t66 );
    const double t69 = alpha * t68;
    const double t72 = 0.2e1 * t30 * t45 * t31;
    const double t73 = rho_b * rho_b;
    const double t74 = t13 * t73;
    const double t75 = t44 * t6;
    const double t76 = 0.1e1 / t75;
    const double t77 = t31 * t31;
    const double t78 = 0.1e1 / t77;
    const double t79 = t76 * t78;
    const double t83 = piecewise_functor_3( t29, 0.0, -t72 - 0.2e1 / 0.3e1 * t74 * t79 );
    const double t89 = t5 * t35 * t63 / 0.8e1;
    const double t91 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t83 * t20 - t89 );
    const double t92 = alpha * t91;
    const double t95 = rho_a * rho_a;
    const double t96 = t13 * t95;
    const double t97 = t76 * t50;
    const double t101 = piecewise_functor_3( t10, 0.0, -t48 - 0.2e1 / 0.3e1 * t96 * t97 );
    const double t106 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t101 * t20 - t66 );
    const double t107 = alpha * t106;
    const double t110 = t7 * t78;
    const double t112 = -rho_b * t45 + t7;
    const double t117 = piecewise_functor_3( t29, 0.0, 0.2e1 * t41 * t31 - t72 + 0.2e1 / 0.3e1 * t30 * t110 * t112 );
    const double t122 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t20 - t89 );
    const double t123 = alpha * t122;
    const double t128 = t13 * t45;
    const double t129 = t128 * t15;
    const double t131 = t50 * t53;
    const double t136 = 0.4e1 * t14 * t76 * t15;
    const double t137 = t45 * t50;
    const double t139 = t14 * t137 * t53;
    const double t142 = 0.1e1 / t49 / t8;
    const double t143 = t7 * t142;
    const double t144 = t53 * t53;
    const double t150 = 0.2e1 * rho_a * t76 - 0.2e1 * t45;
    const double t155 = piecewise_functor_3( t10, 0.0, -0.4e1 * t129 + 0.4e1 / 0.3e1 * t41 * t131 + t136 - 0.4e1 / 0.3e1 * t139 - 0.4e1 / 0.9e1 * t14 * t143 * t144 + 0.2e1 / 0.3e1 * t14 * t51 * t150 );
    const double t160 = t5 * t58 * t63;
    const double t163 = 0.1e1 / t62 / t6;
    const double t166 = t5 * t19 * t163 / 0.12e2;
    const double t168 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t155 * t20 - t160 / 0.4e1 + t166 );
    const double t169 = alpha * t168;
    const double t172 = 0.4e1 * t30 * t76 * t31;
    const double t173 = t44 * t44;
    const double t174 = 0.1e1 / t173;
    const double t175 = t174 * t78;
    const double t176 = t74 * t175;
    const double t178 = t73 * rho_b;
    const double t179 = t13 * t178;
    const double t181 = 0.1e1 / t173 / t6;
    const double t183 = 0.1e1 / t77 / t27;
    const double t184 = t181 * t183;
    const double t188 = piecewise_functor_3( t29, 0.0, t172 + 0.8e1 / 0.3e1 * t176 - 0.4e1 / 0.9e1 * t179 * t184 );
    const double t193 = t5 * t83 * t63;
    const double t197 = t5 * t35 * t163 / 0.12e2;
    const double t199 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t188 * t20 - t193 / 0.4e1 + t197 );
    const double t200 = alpha * t199;
    const double t207 = t174 * t50;
    const double t208 = t96 * t207;
    const double t220 = t5 * t101 * t63;
    const double t226 = t128 * t31;
    const double t233 = t45 * t78;
    const double t235 = t30 * t233 * t112;
    const double t241 = rho_b * t76;
    const double t253 = t5 * t117 * t63;
    const double t264 = t95 * rho_a;
    const double t265 = t13 * t264;
    const double t266 = t181 * t142;
    const double t270 = piecewise_functor_3( t10, 0.0, t136 + 0.8e1 / 0.3e1 * t208 - 0.4e1 / 0.9e1 * t265 * t266 );
    const double t276 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t270 * t20 - t220 / 0.4e1 + t166 );
    const double t277 = alpha * t276;
    const double t279 = t78 * t112;
    const double t283 = t7 * t183;
    const double t284 = t112 * t112;
    const double t289 = -0.2e1 * t45 + 0.2e1 * t241;
    const double t294 = piecewise_functor_3( t29, 0.0, -0.4e1 * t226 + 0.4e1 / 0.3e1 * t41 * t279 + t172 - 0.4e1 / 0.3e1 * t235 - 0.4e1 / 0.9e1 * t30 * t283 * t284 + 0.2e1 / 0.3e1 * t30 * t110 * t289 );
    const double t300 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t294 * t20 - t253 / 0.4e1 + t197 );
    const double t301 = alpha * t300;


    vrho_a = t25 + t40 + t6 * ( t69 + t92 );
    vrho_b = t25 + t40 + t6 * ( t107 + t123 );
    v2rho2_aa = 0.2e1 * t69 + 0.2e1 * t92 + t6 * ( t169 + t200 );
    v2rho2_bb = 0.2e1 * t107 + 0.2e1 * t123 + t6 * ( t277 + t301 );
    v2rho2_ab = 0.0;



  }


};

struct BuiltinSlaterExchange : detail::BuiltinKernelImpl< BuiltinSlaterExchange > {

  BuiltinSlaterExchange( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSlaterExchange >(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

};



} // namespace ExchCXX
