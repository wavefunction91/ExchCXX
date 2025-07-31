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
struct kernel_traits< BuiltinOPTX_X > :
  public gga_screening_interface< BuiltinOPTX_X > {

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


  static constexpr double a = 1.05151;
  static constexpr double b = 1.43169/constants::X_FACTOR_C;
  static constexpr double gamma = 0.006;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t24 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = gamma * gamma;
    constexpr double t21 = b * t20;
    constexpr double t32 = t24 * t24;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t22 = sigma * sigma;
    const double t23 = t21 * t22;
    const double t25 = rho * rho;
    const double t26 = t25 * t25;
    const double t27 = t26 * rho;
    const double t33 = t18 * t18;
    const double t35 = 0.1e1 / t33 / t25;
    const double t38 = gamma * sigma * t32 * t35 + 0.1e1;
    const double t39 = t38 * t38;
    const double t40 = 0.1e1 / t39;
    const double t41 = t24 / t18 / t27 * t40;
    const double t44 = 0.2e1 * t23 * t41 + a;
    const double t48 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t44 );


    eps = 0.2e1 * t48;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t24 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = gamma * gamma;
    constexpr double t21 = b * t20;
    constexpr double t32 = t24 * t24;
    constexpr double t62 = b * t20 * gamma;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t22 = sigma * sigma;
    const double t23 = t21 * t22;
    const double t25 = rho * rho;
    const double t26 = t25 * t25;
    const double t27 = t26 * rho;
    const double t33 = t18 * t18;
    const double t35 = 0.1e1 / t33 / t25;
    const double t38 = gamma * sigma * t32 * t35 + 0.1e1;
    const double t39 = t38 * t38;
    const double t40 = 0.1e1 / t39;
    const double t41 = t24 / t18 / t27 * t40;
    const double t44 = 0.2e1 * t23 * t41 + a;
    const double t48 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t44 );
    const double t50 = t17 / t33;
    const double t54 = t26 * t25;
    const double t58 = t24 / t18 / t54 * t40;
    const double t63 = t22 * sigma;
    const double t64 = t26 * t26;
    const double t65 = t64 * rho;
    const double t66 = 0.1e1 / t65;
    const double t69 = 0.1e1 / t39 / t38;
    const double t73 = -0.32e2 / 0.3e1 * t23 * t58 + 0.64e2 / 0.3e1 * t62 * t63 * t66 * t69;
    const double t78 = piecewise_functor_3( t2, 0.0, -t6 * t50 * t44 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t73 );
    const double t81 = t21 * sigma;
    const double t84 = 0.1e1 / t64;
    const double t89 = -0.8e1 * t62 * t22 * t84 * t69 + 0.4e1 * t81 * t41;
    const double t93 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t89 );


    eps = 0.2e1 * t48;
    vrho = 0.2e1 * rho * t78 + 0.2e1 * t48;
    vsigma = 0.2e1 * rho * t93;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t24 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = gamma * gamma;
    constexpr double t21 = b * t20;
    constexpr double t32 = t24 * t24;
    constexpr double t62 = b * t20 * gamma;
    constexpr double t119 = t20 * t20;
    constexpr double t120 = b * t119;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t22 = sigma * sigma;
    const double t23 = t21 * t22;
    const double t25 = rho * rho;
    const double t26 = t25 * t25;
    const double t27 = t26 * rho;
    const double t33 = t18 * t18;
    const double t35 = 0.1e1 / t33 / t25;
    const double t38 = gamma * sigma * t32 * t35 + 0.1e1;
    const double t39 = t38 * t38;
    const double t40 = 0.1e1 / t39;
    const double t41 = t24 / t18 / t27 * t40;
    const double t44 = 0.2e1 * t23 * t41 + a;
    const double t50 = t17 / t33;
    const double t54 = t26 * t25;
    const double t58 = t24 / t18 / t54 * t40;
    const double t63 = t22 * sigma;
    const double t64 = t26 * t26;
    const double t65 = t64 * rho;
    const double t66 = 0.1e1 / t65;
    const double t69 = 0.1e1 / t39 / t38;
    const double t73 = -0.32e2 / 0.3e1 * t23 * t58 + 0.64e2 / 0.3e1 * t62 * t63 * t66 * t69;
    const double t78 = piecewise_functor_3( t2, 0.0, -t6 * t50 * t44 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t73 );
    const double t81 = t21 * sigma;
    const double t84 = 0.1e1 / t64;
    const double t89 = -0.8e1 * t62 * t22 * t84 * t69 + 0.4e1 * t81 * t41;
    const double t93 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t89 );
    const double t98 = t17 / t33 / rho;
    const double t105 = t25 * rho;
    const double t106 = t26 * t105;
    const double t110 = t24 / t18 / t106 * t40;
    const double t113 = t64 * t25;
    const double t114 = 0.1e1 / t113;
    const double t121 = t22 * t22;
    const double t122 = t120 * t121;
    const double t123 = t64 * t26;
    const double t125 = 0.1e1 / t33 / t123;
    const double t126 = t39 * t39;
    const double t127 = 0.1e1 / t126;
    const double t129 = t125 * t127 * t32;
    const double t132 = 0.608e3 / 0.9e1 * t23 * t110 - 0.2752e4 / 0.9e1 * t62 * t63 * t114 * t69 + 0.512e3 / 0.3e1 * t122 * t129;
    const double t137 = piecewise_functor_3( t2, 0.0, t6 * t98 * t44 / 0.12e2 - t6 * t50 * t73 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t132 );
    const double t149 = t120 * t63;
    const double t150 = t64 * t105;
    const double t152 = 0.1e1 / t33 / t150;
    const double t154 = t152 * t127 * t32;
    const double t157 = -0.64e2 / 0.3e1 * t81 * t58 + 0.32e3 / 0.3e1 * t62 * t22 * t66 * t69 - 0.64e2 * t149 * t154;
    const double t162 = piecewise_functor_3( t2, 0.0, -t6 * t50 * t89 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t157 );
    const double t171 = t120 * t22;
    const double t175 = 0.1e1 / t33 / t113 * t127 * t32;
    const double t178 = -0.32e2 * t62 * sigma * t84 * t69 + 0.24e2 * t171 * t175 + 0.4e1 * t21 * t41;
    const double t182 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t178 );


    v2rho2 = 0.2e1 * rho * t137 + 0.4e1 * t78;
    v2rhosigma = 0.2e1 * rho * t162 + 0.2e1 * t93;
    v2sigma2 = 0.2e1 * rho * t182;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double sigma, double& vrho, double& vsigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t24 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = gamma * gamma;
    constexpr double t21 = b * t20;
    constexpr double t32 = t24 * t24;
    constexpr double t62 = b * t20 * gamma;
    constexpr double t119 = t20 * t20;
    constexpr double t120 = b * t119;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t22 = sigma * sigma;
    const double t23 = t21 * t22;
    const double t25 = rho * rho;
    const double t26 = t25 * t25;
    const double t27 = t26 * rho;
    const double t33 = t18 * t18;
    const double t35 = 0.1e1 / t33 / t25;
    const double t38 = gamma * sigma * t32 * t35 + 0.1e1;
    const double t39 = t38 * t38;
    const double t40 = 0.1e1 / t39;
    const double t41 = t24 / t18 / t27 * t40;
    const double t44 = 0.2e1 * t23 * t41 + a;
    const double t48 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t44 );
    const double t50 = t17 / t33;
    const double t54 = t26 * t25;
    const double t58 = t24 / t18 / t54 * t40;
    const double t63 = t22 * sigma;
    const double t64 = t26 * t26;
    const double t65 = t64 * rho;
    const double t66 = 0.1e1 / t65;
    const double t69 = 0.1e1 / t39 / t38;
    const double t73 = -0.32e2 / 0.3e1 * t23 * t58 + 0.64e2 / 0.3e1 * t62 * t63 * t66 * t69;
    const double t78 = piecewise_functor_3( t2, 0.0, -t6 * t50 * t44 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t73 );
    const double t81 = t21 * sigma;
    const double t84 = 0.1e1 / t64;
    const double t89 = -0.8e1 * t62 * t22 * t84 * t69 + 0.4e1 * t81 * t41;
    const double t93 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t89 );
    const double t98 = t17 / t33 / rho;
    const double t105 = t25 * rho;
    const double t106 = t26 * t105;
    const double t110 = t24 / t18 / t106 * t40;
    const double t113 = t64 * t25;
    const double t114 = 0.1e1 / t113;
    const double t121 = t22 * t22;
    const double t122 = t120 * t121;
    const double t123 = t64 * t26;
    const double t125 = 0.1e1 / t33 / t123;
    const double t126 = t39 * t39;
    const double t127 = 0.1e1 / t126;
    const double t129 = t125 * t127 * t32;
    const double t132 = 0.608e3 / 0.9e1 * t23 * t110 - 0.2752e4 / 0.9e1 * t62 * t63 * t114 * t69 + 0.512e3 / 0.3e1 * t122 * t129;
    const double t137 = piecewise_functor_3( t2, 0.0, t6 * t98 * t44 / 0.12e2 - t6 * t50 * t73 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t132 );
    const double t149 = t120 * t63;
    const double t150 = t64 * t105;
    const double t152 = 0.1e1 / t33 / t150;
    const double t154 = t152 * t127 * t32;
    const double t157 = -0.64e2 / 0.3e1 * t81 * t58 + 0.32e3 / 0.3e1 * t62 * t22 * t66 * t69 - 0.64e2 * t149 * t154;
    const double t162 = piecewise_functor_3( t2, 0.0, -t6 * t50 * t89 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t157 );
    const double t171 = t120 * t22;
    const double t175 = 0.1e1 / t33 / t113 * t127 * t32;
    const double t178 = -0.32e2 * t62 * sigma * t84 * t69 + 0.24e2 * t171 * t175 + 0.4e1 * t21 * t41;
    const double t182 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t178 );


    vrho = 0.2e1 * rho * t78 + 0.2e1 * t48;
    vsigma = 0.2e1 * rho * t93;
    v2rho2 = 0.2e1 * rho * t137 + 0.4e1 * t78;
    v2rhosigma = 0.2e1 * rho * t162 + 0.2e1 * t93;
    v2sigma2 = 0.2e1 * rho * t182;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = gamma * gamma;
    constexpr double t29 = b * t28;


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
    const double t30 = sigma_aa * sigma_aa;
    const double t31 = rho_a * rho_a;
    const double t32 = t31 * t31;
    const double t33 = t32 * rho_a;
    const double t34 = safe_math::cbrt( rho_a );
    const double t36 = 0.1e1 / t34 / t33;
    const double t39 = t34 * t34;
    const double t43 = 0.1e1 + gamma * sigma_aa / t39 / t31;
    const double t44 = t43 * t43;
    const double t45 = 0.1e1 / t44;
    const double t48 = t29 * t30 * t36 * t45 + a;
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = sigma_bb * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = t64 * t64;
    const double t66 = t65 * rho_b;
    const double t67 = safe_math::cbrt( rho_b );
    const double t69 = 0.1e1 / t67 / t66;
    const double t72 = t67 * t67;
    const double t76 = 0.1e1 + gamma * sigma_bb / t72 / t64;
    const double t77 = t76 * t76;
    const double t78 = 0.1e1 / t77;
    const double t81 = t29 * t63 * t69 * t78 + a;
    const double t85 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t81 );


    eps = t52 + t85;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = gamma * gamma;
    constexpr double t29 = b * t28;
    constexpr double t111 = b * t28 * gamma;


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
    const double t30 = sigma_aa * sigma_aa;
    const double t31 = rho_a * rho_a;
    const double t32 = t31 * t31;
    const double t33 = t32 * rho_a;
    const double t34 = safe_math::cbrt( rho_a );
    const double t36 = 0.1e1 / t34 / t33;
    const double t39 = t34 * t34;
    const double t43 = 0.1e1 + gamma * sigma_aa / t39 / t31;
    const double t44 = t43 * t43;
    const double t45 = 0.1e1 / t44;
    const double t48 = t29 * t30 * t36 * t45 + a;
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = sigma_bb * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = t64 * t64;
    const double t66 = t65 * rho_b;
    const double t67 = safe_math::cbrt( rho_b );
    const double t69 = 0.1e1 / t67 / t66;
    const double t72 = t67 * t67;
    const double t76 = 0.1e1 + gamma * sigma_bb / t72 / t64;
    const double t77 = t76 * t76;
    const double t78 = 0.1e1 / t77;
    const double t81 = t29 * t63 * t69 * t78 + a;
    const double t85 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t81 );
    const double t86 = t6 * t6;
    const double t87 = 0.1e1 / t86;
    const double t88 = t16 * t87;
    const double t90 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t88 );
    const double t93 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t90 );
    const double t94 = t93 * t26;
    const double t98 = t26 * t26;
    const double t99 = 0.1e1 / t98;
    const double t100 = t25 * t99;
    const double t103 = t5 * t100 * t48 / 0.8e1;
    const double t104 = t32 * t31;
    const double t106 = 0.1e1 / t34 / t104;
    const double t112 = t30 * sigma_aa;
    const double t113 = t32 * t32;
    const double t114 = t113 * rho_a;
    const double t115 = 0.1e1 / t114;
    const double t118 = 0.1e1 / t44 / t43;
    const double t122 = -0.16e2 / 0.3e1 * t29 * t30 * t106 * t45 + 0.16e2 / 0.3e1 * t111 * t112 * t115 * t118;
    const double t127 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t94 * t48 - t103 - 0.3e1 / 0.8e1 * t5 * t27 * t122 );
    const double t128 = t54 * t87;
    const double t130 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t128 );
    const double t133 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t130 );
    const double t134 = t133 * t26;
    const double t138 = t61 * t99;
    const double t141 = t5 * t138 * t81 / 0.8e1;
    const double t143 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t81 - t141 );
    const double t147 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t88 );
    const double t150 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t147 );
    const double t151 = t150 * t26;
    const double t156 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t151 * t48 - t103 );
    const double t158 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t128 );
    const double t161 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t158 );
    const double t162 = t161 * t26;
    const double t166 = t65 * t64;
    const double t168 = 0.1e1 / t67 / t166;
    const double t172 = t63 * sigma_bb;
    const double t173 = t65 * t65;
    const double t174 = t173 * rho_b;
    const double t175 = 0.1e1 / t174;
    const double t178 = 0.1e1 / t77 / t76;
    const double t182 = 0.16e2 / 0.3e1 * t111 * t172 * t175 * t178 - 0.16e2 / 0.3e1 * t29 * t63 * t168 * t78;
    const double t187 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t162 * t81 - t141 - 0.3e1 / 0.8e1 * t5 * t62 * t182 );
    const double t193 = 0.1e1 / t113;
    const double t198 = -0.2e1 * t111 * t30 * t193 * t118 + 0.2e1 * t29 * sigma_aa * t36 * t45;
    const double t202 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t198 );
    const double t206 = 0.1e1 / t173;
    const double t211 = -0.2e1 * t111 * t63 * t206 * t178 + 0.2e1 * t29 * sigma_bb * t69 * t78;
    const double t215 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t211 );


    eps = t52 + t85;
    vrho_a = t52 + t85 + t6 * ( t127 + t143 );
    vrho_b = t52 + t85 + t6 * ( t156 + t187 );
    vsigma_aa = t6 * t202;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t215;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = gamma * gamma;
    constexpr double t29 = b * t28;
    constexpr double t111 = b * t28 * gamma;
    constexpr double t267 = t28 * t28;
    constexpr double t268 = b * t267;


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
    const double t30 = sigma_aa * sigma_aa;
    const double t31 = rho_a * rho_a;
    const double t32 = t31 * t31;
    const double t33 = t32 * rho_a;
    const double t34 = safe_math::cbrt( rho_a );
    const double t36 = 0.1e1 / t34 / t33;
    const double t39 = t34 * t34;
    const double t43 = 0.1e1 + gamma * sigma_aa / t39 / t31;
    const double t44 = t43 * t43;
    const double t45 = 0.1e1 / t44;
    const double t48 = t29 * t30 * t36 * t45 + a;
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = sigma_bb * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = t64 * t64;
    const double t66 = t65 * rho_b;
    const double t67 = safe_math::cbrt( rho_b );
    const double t69 = 0.1e1 / t67 / t66;
    const double t72 = t67 * t67;
    const double t76 = 0.1e1 + gamma * sigma_bb / t72 / t64;
    const double t77 = t76 * t76;
    const double t78 = 0.1e1 / t77;
    const double t81 = t29 * t63 * t69 * t78 + a;
    const double t86 = t6 * t6;
    const double t87 = 0.1e1 / t86;
    const double t88 = t16 * t87;
    const double t90 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t88 );
    const double t93 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t90 );
    const double t94 = t93 * t26;
    const double t98 = t26 * t26;
    const double t99 = 0.1e1 / t98;
    const double t100 = t25 * t99;
    const double t103 = t5 * t100 * t48 / 0.8e1;
    const double t104 = t32 * t31;
    const double t106 = 0.1e1 / t34 / t104;
    const double t112 = t30 * sigma_aa;
    const double t113 = t32 * t32;
    const double t114 = t113 * rho_a;
    const double t115 = 0.1e1 / t114;
    const double t118 = 0.1e1 / t44 / t43;
    const double t122 = -0.16e2 / 0.3e1 * t29 * t30 * t106 * t45 + 0.16e2 / 0.3e1 * t111 * t112 * t115 * t118;
    const double t127 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t94 * t48 - t103 - 0.3e1 / 0.8e1 * t5 * t27 * t122 );
    const double t128 = t54 * t87;
    const double t130 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t128 );
    const double t133 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t130 );
    const double t134 = t133 * t26;
    const double t138 = t61 * t99;
    const double t141 = t5 * t138 * t81 / 0.8e1;
    const double t143 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t81 - t141 );
    const double t147 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t88 );
    const double t150 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t147 );
    const double t151 = t150 * t26;
    const double t156 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t151 * t48 - t103 );
    const double t158 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t128 );
    const double t161 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t158 );
    const double t162 = t161 * t26;
    const double t166 = t65 * t64;
    const double t168 = 0.1e1 / t67 / t166;
    const double t172 = t63 * sigma_bb;
    const double t173 = t65 * t65;
    const double t174 = t173 * rho_b;
    const double t175 = 0.1e1 / t174;
    const double t178 = 0.1e1 / t77 / t76;
    const double t182 = 0.16e2 / 0.3e1 * t111 * t172 * t175 * t178 - 0.16e2 / 0.3e1 * t29 * t63 * t168 * t78;
    const double t187 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t162 * t81 - t141 - 0.3e1 / 0.8e1 * t5 * t62 * t182 );
    const double t193 = 0.1e1 / t113;
    const double t198 = -0.2e1 * t111 * t30 * t193 * t118 + 0.2e1 * t29 * sigma_aa * t36 * t45;
    const double t202 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t198 );
    const double t206 = 0.1e1 / t173;
    const double t211 = -0.2e1 * t111 * t63 * t206 * t178 + 0.2e1 * t29 * sigma_bb * t69 * t78;
    const double t215 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t211 );
    const double t218 = t23 * t23;
    const double t219 = 0.1e1 / t218;
    const double t220 = t90 * t90;
    const double t223 = t86 * t6;
    const double t224 = 0.1e1 / t223;
    const double t225 = t16 * t224;
    const double t228 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t87 + 0.2e1 * t225 );
    const double t232 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t219 * t220 + 0.4e1 / 0.3e1 * t23 * t228 );
    const double t233 = t232 * t26;
    const double t237 = t93 * t99;
    const double t239 = t5 * t237 * t48;
    const double t245 = 0.1e1 / t98 / t6;
    const double t246 = t25 * t245;
    const double t249 = t5 * t246 * t48 / 0.12e2;
    const double t251 = t5 * t100 * t122;
    const double t253 = t31 * rho_a;
    const double t254 = t32 * t253;
    const double t256 = 0.1e1 / t34 / t254;
    const double t261 = t113 * t31;
    const double t262 = 0.1e1 / t261;
    const double t269 = t30 * t30;
    const double t270 = t113 * t32;
    const double t272 = 0.1e1 / t39 / t270;
    const double t274 = t44 * t44;
    const double t275 = 0.1e1 / t274;
    const double t279 = 0.304e3 / 0.9e1 * t29 * t30 * t256 * t45 - 0.688e3 / 0.9e1 * t111 * t112 * t262 * t118 + 0.128e3 / 0.3e1 * t268 * t269 * t272 * t275;
    const double t284 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t233 * t48 - t239 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t94 * t122 + t249 - t251 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t279 );
    const double t285 = t59 * t59;
    const double t286 = 0.1e1 / t285;
    const double t287 = t130 * t130;
    const double t290 = t54 * t224;
    const double t293 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t87 + 0.2e1 * t290 );
    const double t297 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t286 * t287 + 0.4e1 / 0.3e1 * t59 * t293 );
    const double t298 = t297 * t26;
    const double t302 = t133 * t99;
    const double t304 = t5 * t302 * t81;
    const double t306 = t61 * t245;
    const double t309 = t5 * t306 * t81 / 0.12e2;
    const double t311 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t298 * t81 - t304 / 0.4e1 + t309 );
    const double t327 = t150 * t99;
    const double t329 = t5 * t327 * t48;
    const double t351 = t161 * t99;
    const double t353 = t5 * t351 * t81;
    const double t360 = t5 * t138 * t182;
    const double t368 = t147 * t147;
    const double t373 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t87 + 0.2e1 * t225 );
    const double t377 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t219 * t368 + 0.4e1 / 0.3e1 * t23 * t373 );
    const double t378 = t377 * t26;
    const double t384 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t378 * t48 - t329 / 0.4e1 + t249 );
    const double t385 = t158 * t158;
    const double t390 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t87 + 0.2e1 * t290 );
    const double t394 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t286 * t385 + 0.4e1 / 0.3e1 * t59 * t390 );
    const double t395 = t394 * t26;
    const double t404 = t64 * rho_b;
    const double t405 = t65 * t404;
    const double t407 = 0.1e1 / t67 / t405;
    const double t412 = t173 * t64;
    const double t413 = 0.1e1 / t412;
    const double t418 = t63 * t63;
    const double t419 = t173 * t65;
    const double t421 = 0.1e1 / t72 / t419;
    const double t423 = t77 * t77;
    const double t424 = 0.1e1 / t423;
    const double t428 = 0.304e3 / 0.9e1 * t29 * t63 * t407 * t78 - 0.688e3 / 0.9e1 * t111 * t172 * t413 * t178 + 0.128e3 / 0.3e1 * t268 * t418 * t421 * t424;
    const double t433 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t395 * t81 - t353 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t162 * t182 + t309 - t360 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t62 * t428 );
    const double t441 = t5 * t100 * t198 / 0.8e1;
    const double t450 = t113 * t253;
    const double t452 = 0.1e1 / t39 / t450;
    const double t457 = -0.32e2 / 0.3e1 * t29 * sigma_aa * t106 * t45 + 0.8e2 / 0.3e1 * t111 * t30 * t115 * t118 - 0.16e2 * t268 * t112 * t452 * t275;
    const double t462 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t94 * t198 - t441 - 0.3e1 / 0.8e1 * t5 * t27 * t457 );
    const double t469 = t5 * t138 * t211 / 0.8e1;
    const double t490 = t173 * t404;
    const double t492 = 0.1e1 / t72 / t490;
    const double t497 = -0.32e2 / 0.3e1 * t29 * sigma_bb * t168 * t78 + 0.8e2 / 0.3e1 * t111 * t63 * t175 * t178 - 0.16e2 * t268 * t172 * t492 * t424;
    const double t502 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t162 * t211 - t469 - 0.3e1 / 0.8e1 * t5 * t62 * t497 );
    const double t512 = 0.1e1 / t39 / t261;
    const double t517 = -0.8e1 * t111 * sigma_aa * t193 * t118 + 0.6e1 * t268 * t30 * t512 * t275 + 0.2e1 * t29 * t36 * t45;
    const double t521 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t517 );
    const double t530 = 0.1e1 / t72 / t412;
    const double t535 = -0.8e1 * t111 * sigma_bb * t206 * t178 + 0.6e1 * t268 * t63 * t530 * t424 + 0.2e1 * t29 * t69 * t78;
    const double t539 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t535 );


    v2rho2_aa = 0.2e1 * t127 + 0.2e1 * t143 + t6 * ( t284 + t311 );
    v2rho2_bb = 0.2e1 * t156 + 0.2e1 * t187 + t6 * ( t384 + t433 );
    v2rhosigma_a_aa = t6 * t462 + t202;
    v2rhosigma_b_bb = t6 * t502 + t215;
    v2sigma2_aa_aa = t6 * t521;
    v2sigma2_bb_bb = t6 * t539;
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
    constexpr double t5 = t2 / t3;
    constexpr double t28 = gamma * gamma;
    constexpr double t29 = b * t28;
    constexpr double t111 = b * t28 * gamma;
    constexpr double t267 = t28 * t28;
    constexpr double t268 = b * t267;


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
    const double t30 = sigma_aa * sigma_aa;
    const double t31 = rho_a * rho_a;
    const double t32 = t31 * t31;
    const double t33 = t32 * rho_a;
    const double t34 = safe_math::cbrt( rho_a );
    const double t36 = 0.1e1 / t34 / t33;
    const double t39 = t34 * t34;
    const double t43 = 0.1e1 + gamma * sigma_aa / t39 / t31;
    const double t44 = t43 * t43;
    const double t45 = 0.1e1 / t44;
    const double t48 = t29 * t30 * t36 * t45 + a;
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = sigma_bb * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = t64 * t64;
    const double t66 = t65 * rho_b;
    const double t67 = safe_math::cbrt( rho_b );
    const double t69 = 0.1e1 / t67 / t66;
    const double t72 = t67 * t67;
    const double t76 = 0.1e1 + gamma * sigma_bb / t72 / t64;
    const double t77 = t76 * t76;
    const double t78 = 0.1e1 / t77;
    const double t81 = t29 * t63 * t69 * t78 + a;
    const double t85 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t81 );
    const double t86 = t6 * t6;
    const double t87 = 0.1e1 / t86;
    const double t88 = t16 * t87;
    const double t90 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t88 );
    const double t93 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t90 );
    const double t94 = t93 * t26;
    const double t98 = t26 * t26;
    const double t99 = 0.1e1 / t98;
    const double t100 = t25 * t99;
    const double t103 = t5 * t100 * t48 / 0.8e1;
    const double t104 = t32 * t31;
    const double t106 = 0.1e1 / t34 / t104;
    const double t112 = t30 * sigma_aa;
    const double t113 = t32 * t32;
    const double t114 = t113 * rho_a;
    const double t115 = 0.1e1 / t114;
    const double t118 = 0.1e1 / t44 / t43;
    const double t122 = -0.16e2 / 0.3e1 * t29 * t30 * t106 * t45 + 0.16e2 / 0.3e1 * t111 * t112 * t115 * t118;
    const double t127 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t94 * t48 - t103 - 0.3e1 / 0.8e1 * t5 * t27 * t122 );
    const double t128 = t54 * t87;
    const double t130 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t128 );
    const double t133 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t130 );
    const double t134 = t133 * t26;
    const double t138 = t61 * t99;
    const double t141 = t5 * t138 * t81 / 0.8e1;
    const double t143 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t134 * t81 - t141 );
    const double t147 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t88 );
    const double t150 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t147 );
    const double t151 = t150 * t26;
    const double t156 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t151 * t48 - t103 );
    const double t158 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t128 );
    const double t161 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t158 );
    const double t162 = t161 * t26;
    const double t166 = t65 * t64;
    const double t168 = 0.1e1 / t67 / t166;
    const double t172 = t63 * sigma_bb;
    const double t173 = t65 * t65;
    const double t174 = t173 * rho_b;
    const double t175 = 0.1e1 / t174;
    const double t178 = 0.1e1 / t77 / t76;
    const double t182 = 0.16e2 / 0.3e1 * t111 * t172 * t175 * t178 - 0.16e2 / 0.3e1 * t29 * t63 * t168 * t78;
    const double t187 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t162 * t81 - t141 - 0.3e1 / 0.8e1 * t5 * t62 * t182 );
    const double t193 = 0.1e1 / t113;
    const double t198 = -0.2e1 * t111 * t30 * t193 * t118 + 0.2e1 * t29 * sigma_aa * t36 * t45;
    const double t202 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t198 );
    const double t206 = 0.1e1 / t173;
    const double t211 = -0.2e1 * t111 * t63 * t206 * t178 + 0.2e1 * t29 * sigma_bb * t69 * t78;
    const double t215 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t211 );
    const double t218 = t23 * t23;
    const double t219 = 0.1e1 / t218;
    const double t220 = t90 * t90;
    const double t223 = t86 * t6;
    const double t224 = 0.1e1 / t223;
    const double t225 = t16 * t224;
    const double t228 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t87 + 0.2e1 * t225 );
    const double t232 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t219 * t220 + 0.4e1 / 0.3e1 * t23 * t228 );
    const double t233 = t232 * t26;
    const double t237 = t93 * t99;
    const double t239 = t5 * t237 * t48;
    const double t245 = 0.1e1 / t98 / t6;
    const double t246 = t25 * t245;
    const double t249 = t5 * t246 * t48 / 0.12e2;
    const double t251 = t5 * t100 * t122;
    const double t253 = t31 * rho_a;
    const double t254 = t32 * t253;
    const double t256 = 0.1e1 / t34 / t254;
    const double t261 = t113 * t31;
    const double t262 = 0.1e1 / t261;
    const double t269 = t30 * t30;
    const double t270 = t113 * t32;
    const double t272 = 0.1e1 / t39 / t270;
    const double t274 = t44 * t44;
    const double t275 = 0.1e1 / t274;
    const double t279 = 0.304e3 / 0.9e1 * t29 * t30 * t256 * t45 - 0.688e3 / 0.9e1 * t111 * t112 * t262 * t118 + 0.128e3 / 0.3e1 * t268 * t269 * t272 * t275;
    const double t284 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t233 * t48 - t239 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t94 * t122 + t249 - t251 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t279 );
    const double t285 = t59 * t59;
    const double t286 = 0.1e1 / t285;
    const double t287 = t130 * t130;
    const double t290 = t54 * t224;
    const double t293 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t87 + 0.2e1 * t290 );
    const double t297 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t286 * t287 + 0.4e1 / 0.3e1 * t59 * t293 );
    const double t298 = t297 * t26;
    const double t302 = t133 * t99;
    const double t304 = t5 * t302 * t81;
    const double t306 = t61 * t245;
    const double t309 = t5 * t306 * t81 / 0.12e2;
    const double t311 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t298 * t81 - t304 / 0.4e1 + t309 );
    const double t327 = t150 * t99;
    const double t329 = t5 * t327 * t48;
    const double t351 = t161 * t99;
    const double t353 = t5 * t351 * t81;
    const double t360 = t5 * t138 * t182;
    const double t368 = t147 * t147;
    const double t373 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t87 + 0.2e1 * t225 );
    const double t377 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t219 * t368 + 0.4e1 / 0.3e1 * t23 * t373 );
    const double t378 = t377 * t26;
    const double t384 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t378 * t48 - t329 / 0.4e1 + t249 );
    const double t385 = t158 * t158;
    const double t390 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t87 + 0.2e1 * t290 );
    const double t394 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t286 * t385 + 0.4e1 / 0.3e1 * t59 * t390 );
    const double t395 = t394 * t26;
    const double t404 = t64 * rho_b;
    const double t405 = t65 * t404;
    const double t407 = 0.1e1 / t67 / t405;
    const double t412 = t173 * t64;
    const double t413 = 0.1e1 / t412;
    const double t418 = t63 * t63;
    const double t419 = t173 * t65;
    const double t421 = 0.1e1 / t72 / t419;
    const double t423 = t77 * t77;
    const double t424 = 0.1e1 / t423;
    const double t428 = 0.304e3 / 0.9e1 * t29 * t63 * t407 * t78 - 0.688e3 / 0.9e1 * t111 * t172 * t413 * t178 + 0.128e3 / 0.3e1 * t268 * t418 * t421 * t424;
    const double t433 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t395 * t81 - t353 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t162 * t182 + t309 - t360 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t62 * t428 );
    const double t441 = t5 * t100 * t198 / 0.8e1;
    const double t450 = t113 * t253;
    const double t452 = 0.1e1 / t39 / t450;
    const double t457 = -0.32e2 / 0.3e1 * t29 * sigma_aa * t106 * t45 + 0.8e2 / 0.3e1 * t111 * t30 * t115 * t118 - 0.16e2 * t268 * t112 * t452 * t275;
    const double t462 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t94 * t198 - t441 - 0.3e1 / 0.8e1 * t5 * t27 * t457 );
    const double t469 = t5 * t138 * t211 / 0.8e1;
    const double t490 = t173 * t404;
    const double t492 = 0.1e1 / t72 / t490;
    const double t497 = -0.32e2 / 0.3e1 * t29 * sigma_bb * t168 * t78 + 0.8e2 / 0.3e1 * t111 * t63 * t175 * t178 - 0.16e2 * t268 * t172 * t492 * t424;
    const double t502 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t162 * t211 - t469 - 0.3e1 / 0.8e1 * t5 * t62 * t497 );
    const double t512 = 0.1e1 / t39 / t261;
    const double t517 = -0.8e1 * t111 * sigma_aa * t193 * t118 + 0.6e1 * t268 * t30 * t512 * t275 + 0.2e1 * t29 * t36 * t45;
    const double t521 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t517 );
    const double t530 = 0.1e1 / t72 / t412;
    const double t535 = -0.8e1 * t111 * sigma_bb * t206 * t178 + 0.6e1 * t268 * t63 * t530 * t424 + 0.2e1 * t29 * t69 * t78;
    const double t539 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t535 );


    vrho_a = t52 + t85 + t6 * ( t127 + t143 );
    vrho_b = t52 + t85 + t6 * ( t156 + t187 );
    vsigma_aa = t6 * t202;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t215;
    v2rho2_aa = 0.2e1 * t127 + 0.2e1 * t143 + t6 * ( t284 + t311 );
    v2rho2_bb = 0.2e1 * t156 + 0.2e1 * t187 + t6 * ( t384 + t433 );
    v2rhosigma_a_aa = t6 * t462 + t202;
    v2rhosigma_b_bb = t6 * t502 + t215;
    v2sigma2_aa_aa = t6 * t521;
    v2sigma2_bb_bb = t6 * t539;
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

struct BuiltinOPTX_X : detail::BuiltinKernelImpl< BuiltinOPTX_X > {

  BuiltinOPTX_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinOPTX_X >(p) { }
  
  virtual ~BuiltinOPTX_X() = default;

};



} // namespace ExchCXX
