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
struct kernel_traits< BuiltinEPC18_2 > :
  public lda_screening_interface< BuiltinEPC18_2 > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = true;

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;


  static constexpr double a = 3.90;
  static constexpr double b = 0.50;
  static constexpr double c = 0.06;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;


    eps = piecewise_functor_3( t2, 0.0, -t9 * t18 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;
    const double t22 = t17 * t17;
    const double t23 = 0.1e1 / t22;
    const double t27 = 0.32e2 * t13 * rho - 0.4e1 * t10;
    const double t32 = piecewise_functor_3( t2, 0.0, t9 * t23 * t27 / 0.4e1 - t8 * t18 / 0.4e1 );


    eps = piecewise_functor_3( t2, 0.0, -t9 * t18 / 0.4e1 );
    vrho = rho * t32 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double& v2rho2 ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;
    const double t22 = t17 * t17;
    const double t23 = 0.1e1 / t22;
    const double t27 = 0.32e2 * t13 * rho - 0.4e1 * t10;
    const double t32 = piecewise_functor_3( t2, 0.0, t9 * t23 * t27 / 0.4e1 - t8 * t18 / 0.4e1 );
    const double t39 = 0.1e1 / t22 / t17;
    const double t40 = t27 * t27;
    const double t44 = t8 * t8;
    const double t45 = t44 * rho;
    const double t50 = piecewise_functor_3( t2, 0.0, t8 * t23 * t27 / 0.2e1 - t9 * t39 * t40 / 0.2e1 + 0.8e1 * t45 * t23 * c );


    v2rho2 = rho * t50 + 0.2e1 * t32;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double& vrho, double& v2rho2 ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;
    const double tzk0 = piecewise_functor_3( t2, 0.0, -t9 * t18 / 0.4e1 );
    const double t22 = t17 * t17;
    const double t23 = 0.1e1 / t22;
    const double t27 = 0.32e2 * t13 * rho - 0.4e1 * t10;
    const double t32 = piecewise_functor_3( t2, 0.0, t9 * t23 * t27 / 0.4e1 - t8 * t18 / 0.4e1 );
    const double t39 = 0.1e1 / t22 / t17;
    const double t40 = t27 * t27;
    const double t44 = t8 * t8;
    const double t45 = t44 * rho;
    const double t50 = piecewise_functor_3( t2, 0.0, t8 * t23 * t27 / 0.2e1 - t9 * t39 * t40 / 0.2e1 + 0.8e1 * t45 * t23 * c );


    vrho = rho * t32 + tzk0;
    v2rho2 = rho * t50 + 0.2e1 * t32;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;


    eps = piecewise_functor_3( t3, 0.0, -t18 * t40 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;
    const double t43 = t4 * t4;
    const double t44 = 0.1e1 / t43;
    const double t45 = t14 * t44;
    const double t47 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t45 );
    const double t48 = t47 * t4;
    const double t50 = t17 * t22;
    const double t51 = t50 * t39;
    const double t52 = t19 * t44;
    const double t54 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t52 );
    const double t55 = t54 * t39;
    const double t57 = t38 * t38;
    const double t58 = 0.1e1 / t57;
    const double t59 = t22 * t58;
    const double t60 = b * t32;
    const double t61 = t25 * t25;
    const double t63 = t24 / t61;
    const double t64 = t48 + 0.1e1 + t16;
    const double t66 = t28 * t28;
    const double t68 = t24 / t66;
    const double t70 = t54 * t4 + t21 + 0.1e1;
    const double t73 = t63 * t64 / 0.6e1 + t68 * t70 / 0.6e1;
    const double t77 = c * t35 * t31;
    const double t80 = -0.3e1 * t60 * t73 + 0.6e1 * t77 * t73;
    const double t81 = t59 * t80;
    const double t85 = piecewise_functor_3( t3, 0.0, -t18 * t55 / 0.4e1 + t18 * t81 / 0.4e1 - t48 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t88 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t45 );
    const double t89 = t88 * t4;
    const double t92 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t52 );
    const double t93 = t92 * t39;
    const double t95 = t89 + 0.1e1 + t16;
    const double t98 = t92 * t4 + t21 + 0.1e1;
    const double t101 = t63 * t95 / 0.6e1 + t68 * t98 / 0.6e1;
    const double t106 = -0.3e1 * t60 * t101 + 0.6e1 * t77 * t101;
    const double t107 = t59 * t106;
    const double t111 = piecewise_functor_3( t3, 0.0, t18 * t107 / 0.4e1 - t18 * t93 / 0.4e1 - t89 * t40 / 0.4e1 - t51 / 0.4e1 );


    eps = piecewise_functor_3( t3, 0.0, -t18 * t40 / 0.4e1 );
    vrho_a = t4 * t85 + eps;
    vrho_b = t4 * t111 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;
    const double t43 = t4 * t4;
    const double t44 = 0.1e1 / t43;
    const double t45 = t14 * t44;
    const double t47 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t45 );
    const double t48 = t47 * t4;
    const double t50 = t17 * t22;
    const double t51 = t50 * t39;
    const double t52 = t19 * t44;
    const double t54 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t52 );
    const double t55 = t54 * t39;
    const double t57 = t38 * t38;
    const double t58 = 0.1e1 / t57;
    const double t59 = t22 * t58;
    const double t60 = b * t32;
    const double t61 = t25 * t25;
    const double t63 = t24 / t61;
    const double t64 = t48 + 0.1e1 + t16;
    const double t66 = t28 * t28;
    const double t68 = t24 / t66;
    const double t70 = t54 * t4 + t21 + 0.1e1;
    const double t73 = t63 * t64 / 0.6e1 + t68 * t70 / 0.6e1;
    const double t77 = c * t35 * t31;
    const double t80 = -0.3e1 * t60 * t73 + 0.6e1 * t77 * t73;
    const double t81 = t59 * t80;
    const double t85 = piecewise_functor_3( t3, 0.0, -t18 * t55 / 0.4e1 + t18 * t81 / 0.4e1 - t48 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t88 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t45 );
    const double t89 = t88 * t4;
    const double t92 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t52 );
    const double t93 = t92 * t39;
    const double t95 = t89 + 0.1e1 + t16;
    const double t98 = t92 * t4 + t21 + 0.1e1;
    const double t101 = t63 * t95 / 0.6e1 + t68 * t98 / 0.6e1;
    const double t106 = -0.3e1 * t60 * t101 + 0.6e1 * t77 * t101;
    const double t107 = t59 * t106;
    const double t111 = piecewise_functor_3( t3, 0.0, t18 * t107 / 0.4e1 - t18 * t93 / 0.4e1 - t89 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t114 = t43 * t4;
    const double t115 = 0.1e1 / t114;
    const double t116 = t14 * t115;
    const double t119 = piecewise_functor_5( t8, 0.0, t12, 0.0, -0.2e1 * t44 + 0.2e1 * t116 );
    const double t120 = t119 * t4;
    const double t123 = t47 * t22;
    const double t124 = t123 * t39;
    const double t130 = t17 * t54;
    const double t131 = t130 * t39;
    const double t133 = t58 * t80;
    const double t134 = t50 * t133;
    const double t136 = t19 * t115;
    const double t139 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t44 + 0.2e1 * t136 );
    const double t140 = t139 * t39;
    const double t143 = t54 * t58;
    const double t144 = t143 * t80;
    const double t148 = 0.1e1 / t57 / t38;
    const double t149 = t22 * t148;
    const double t150 = t80 * t80;
    const double t151 = t149 * t150;
    const double t154 = b * t31;
    const double t155 = t73 * t73;
    const double t160 = t24 / t61 / t18;
    const double t161 = t64 * t64;
    const double t165 = t120 + 0.2e1 * t47;
    const double t170 = t24 / t66 / t27;
    const double t171 = t70 * t70;
    const double t176 = t139 * t4 + 0.2e1 * t54;
    const double t179 = -t160 * t161 / 0.9e1 + t63 * t165 / 0.6e1 - t170 * t171 / 0.9e1 + t68 * t176 / 0.6e1;
    const double t182 = c * t35;
    const double t187 = -0.6e1 * t154 * t155 + 0.3e2 * t182 * t155 - 0.3e1 * t60 * t179 + 0.6e1 * t77 * t179;
    const double t188 = t59 * t187;
    const double t192 = piecewise_functor_3( t3, 0.0, -t120 * t40 / 0.4e1 - t124 / 0.2e1 - t48 * t55 / 0.2e1 + t48 * t81 / 0.2e1 - t131 / 0.2e1 + t134 / 0.2e1 - t18 * t140 / 0.4e1 + t18 * t144 / 0.2e1 - t18 * t151 / 0.2e1 + t18 * t188 / 0.4e1 );
    const double t195 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t116 );
    const double t196 = t195 * t4;
    const double t199 = t88 * t22;
    const double t200 = t199 * t39;
    const double t211 = t17 * t92;
    const double t212 = t211 * t39;
    const double t215 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t136 );
    const double t216 = t215 * t39;
    const double t219 = t92 * t58;
    const double t220 = t219 * t80;
    const double t225 = t58 * t106;
    const double t226 = t50 * t225;
    const double t228 = t143 * t106;
    const double t231 = t18 * t22;
    const double t232 = t148 * t106;
    const double t233 = t232 * t80;
    const double t236 = t101 * t73;
    const double t239 = t95 * t64;
    const double t242 = t196 + t88 + t47;
    const double t245 = t98 * t70;
    const double t249 = t215 * t4 + t54 + t92;
    const double t252 = -t160 * t239 / 0.9e1 + t63 * t242 / 0.6e1 - t170 * t245 / 0.9e1 + t68 * t249 / 0.6e1;
    const double t259 = -0.6e1 * t154 * t236 + 0.3e2 * t182 * t236 - 0.3e1 * t60 * t252 + 0.6e1 * t77 * t252;
    const double t260 = t59 * t259;
    const double t263 = -t196 * t40 / 0.4e1 - t200 / 0.4e1 - t89 * t55 / 0.4e1 + t89 * t81 / 0.4e1 - t124 / 0.4e1 - t131 / 0.4e1 + t134 / 0.4e1 - t48 * t93 / 0.4e1 - t212 / 0.4e1 - t18 * t216 / 0.4e1 + t18 * t220 / 0.4e1 + t48 * t107 / 0.4e1 + t226 / 0.4e1 + t18 * t228 / 0.4e1 - t231 * t233 / 0.2e1 + t18 * t260 / 0.4e1;
    const double t264 = piecewise_functor_3( t3, 0.0, t263 );
    const double t269 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t44 + 0.2e1 * t116 );
    const double t270 = t269 * t4;
    const double t282 = piecewise_functor_5( t12, 0.0, t8, 0.0, -0.2e1 * t44 + 0.2e1 * t136 );
    const double t283 = t282 * t39;
    const double t286 = t219 * t106;
    const double t289 = t106 * t106;
    const double t290 = t149 * t289;
    const double t293 = t101 * t101;
    const double t296 = t95 * t95;
    const double t300 = t270 + 0.2e1 * t88;
    const double t303 = t98 * t98;
    const double t308 = t282 * t4 + 0.2e1 * t92;
    const double t311 = -t160 * t296 / 0.9e1 + t63 * t300 / 0.6e1 - t170 * t303 / 0.9e1 + t68 * t308 / 0.6e1;
    const double t318 = -0.6e1 * t154 * t293 + 0.3e2 * t182 * t293 - 0.3e1 * t60 * t311 + 0.6e1 * t77 * t311;
    const double t319 = t59 * t318;
    const double t323 = piecewise_functor_3( t3, 0.0, -t270 * t40 / 0.4e1 - t200 / 0.2e1 - t89 * t93 / 0.2e1 + t89 * t107 / 0.2e1 - t212 / 0.2e1 + t226 / 0.2e1 - t18 * t283 / 0.4e1 + t18 * t286 / 0.2e1 - t18 * t290 / 0.2e1 + t18 * t319 / 0.4e1 );


    v2rho2_aa = t4 * t192 + 0.2e1 * t85;
    v2rho2_ab = t4 * t264 + t111 + t85;
    v2rho2_bb = t4 * t323 + 0.2e1 * t111;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;
    const double tzk0 = piecewise_functor_3( t3, 0.0, -t18 * t40 / 0.4e1 );
    const double t43 = t4 * t4;
    const double t44 = 0.1e1 / t43;
    const double t45 = t14 * t44;
    const double t47 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t45 );
    const double t48 = t47 * t4;
    const double t50 = t17 * t22;
    const double t51 = t50 * t39;
    const double t52 = t19 * t44;
    const double t54 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t52 );
    const double t55 = t54 * t39;
    const double t57 = t38 * t38;
    const double t58 = 0.1e1 / t57;
    const double t59 = t22 * t58;
    const double t60 = b * t32;
    const double t61 = t25 * t25;
    const double t63 = t24 / t61;
    const double t64 = t48 + 0.1e1 + t16;
    const double t66 = t28 * t28;
    const double t68 = t24 / t66;
    const double t70 = t54 * t4 + t21 + 0.1e1;
    const double t73 = t63 * t64 / 0.6e1 + t68 * t70 / 0.6e1;
    const double t77 = c * t35 * t31;
    const double t80 = -0.3e1 * t60 * t73 + 0.6e1 * t77 * t73;
    const double t81 = t59 * t80;
    const double t85 = piecewise_functor_3( t3, 0.0, -t18 * t55 / 0.4e1 + t18 * t81 / 0.4e1 - t48 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t88 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t45 );
    const double t89 = t88 * t4;
    const double t92 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t52 );
    const double t93 = t92 * t39;
    const double t95 = t89 + 0.1e1 + t16;
    const double t98 = t92 * t4 + t21 + 0.1e1;
    const double t101 = t63 * t95 / 0.6e1 + t68 * t98 / 0.6e1;
    const double t106 = -0.3e1 * t60 * t101 + 0.6e1 * t77 * t101;
    const double t107 = t59 * t106;
    const double t111 = piecewise_functor_3( t3, 0.0, t18 * t107 / 0.4e1 - t18 * t93 / 0.4e1 - t89 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t114 = t43 * t4;
    const double t115 = 0.1e1 / t114;
    const double t116 = t14 * t115;
    const double t119 = piecewise_functor_5( t8, 0.0, t12, 0.0, -0.2e1 * t44 + 0.2e1 * t116 );
    const double t120 = t119 * t4;
    const double t123 = t47 * t22;
    const double t124 = t123 * t39;
    const double t130 = t17 * t54;
    const double t131 = t130 * t39;
    const double t133 = t58 * t80;
    const double t134 = t50 * t133;
    const double t136 = t19 * t115;
    const double t139 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t44 + 0.2e1 * t136 );
    const double t140 = t139 * t39;
    const double t143 = t54 * t58;
    const double t144 = t143 * t80;
    const double t148 = 0.1e1 / t57 / t38;
    const double t149 = t22 * t148;
    const double t150 = t80 * t80;
    const double t151 = t149 * t150;
    const double t154 = b * t31;
    const double t155 = t73 * t73;
    const double t160 = t24 / t61 / t18;
    const double t161 = t64 * t64;
    const double t165 = t120 + 0.2e1 * t47;
    const double t170 = t24 / t66 / t27;
    const double t171 = t70 * t70;
    const double t176 = t139 * t4 + 0.2e1 * t54;
    const double t179 = -t160 * t161 / 0.9e1 + t63 * t165 / 0.6e1 - t170 * t171 / 0.9e1 + t68 * t176 / 0.6e1;
    const double t182 = c * t35;
    const double t187 = -0.6e1 * t154 * t155 + 0.3e2 * t182 * t155 - 0.3e1 * t60 * t179 + 0.6e1 * t77 * t179;
    const double t188 = t59 * t187;
    const double t192 = piecewise_functor_3( t3, 0.0, -t120 * t40 / 0.4e1 - t124 / 0.2e1 - t48 * t55 / 0.2e1 + t48 * t81 / 0.2e1 - t131 / 0.2e1 + t134 / 0.2e1 - t18 * t140 / 0.4e1 + t18 * t144 / 0.2e1 - t18 * t151 / 0.2e1 + t18 * t188 / 0.4e1 );
    const double t195 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t116 );
    const double t196 = t195 * t4;
    const double t199 = t88 * t22;
    const double t200 = t199 * t39;
    const double t211 = t17 * t92;
    const double t212 = t211 * t39;
    const double t215 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t136 );
    const double t216 = t215 * t39;
    const double t219 = t92 * t58;
    const double t220 = t219 * t80;
    const double t225 = t58 * t106;
    const double t226 = t50 * t225;
    const double t228 = t143 * t106;
    const double t231 = t18 * t22;
    const double t232 = t148 * t106;
    const double t233 = t232 * t80;
    const double t236 = t101 * t73;
    const double t239 = t95 * t64;
    const double t242 = t196 + t88 + t47;
    const double t245 = t98 * t70;
    const double t249 = t215 * t4 + t54 + t92;
    const double t252 = -t160 * t239 / 0.9e1 + t63 * t242 / 0.6e1 - t170 * t245 / 0.9e1 + t68 * t249 / 0.6e1;
    const double t259 = -0.6e1 * t154 * t236 + 0.3e2 * t182 * t236 - 0.3e1 * t60 * t252 + 0.6e1 * t77 * t252;
    const double t260 = t59 * t259;
    const double t263 = -t196 * t40 / 0.4e1 - t200 / 0.4e1 - t89 * t55 / 0.4e1 + t89 * t81 / 0.4e1 - t124 / 0.4e1 - t131 / 0.4e1 + t134 / 0.4e1 - t48 * t93 / 0.4e1 - t212 / 0.4e1 - t18 * t216 / 0.4e1 + t18 * t220 / 0.4e1 + t48 * t107 / 0.4e1 + t226 / 0.4e1 + t18 * t228 / 0.4e1 - t231 * t233 / 0.2e1 + t18 * t260 / 0.4e1;
    const double t264 = piecewise_functor_3( t3, 0.0, t263 );
    const double t269 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t44 + 0.2e1 * t116 );
    const double t270 = t269 * t4;
    const double t282 = piecewise_functor_5( t12, 0.0, t8, 0.0, -0.2e1 * t44 + 0.2e1 * t136 );
    const double t283 = t282 * t39;
    const double t286 = t219 * t106;
    const double t289 = t106 * t106;
    const double t290 = t149 * t289;
    const double t293 = t101 * t101;
    const double t296 = t95 * t95;
    const double t300 = t270 + 0.2e1 * t88;
    const double t303 = t98 * t98;
    const double t308 = t282 * t4 + 0.2e1 * t92;
    const double t311 = -t160 * t296 / 0.9e1 + t63 * t300 / 0.6e1 - t170 * t303 / 0.9e1 + t68 * t308 / 0.6e1;
    const double t318 = -0.6e1 * t154 * t293 + 0.3e2 * t182 * t293 - 0.3e1 * t60 * t311 + 0.6e1 * t77 * t311;
    const double t319 = t59 * t318;
    const double t323 = piecewise_functor_3( t3, 0.0, -t270 * t40 / 0.4e1 - t200 / 0.2e1 - t89 * t93 / 0.2e1 + t89 * t107 / 0.2e1 - t212 / 0.2e1 + t226 / 0.2e1 - t18 * t283 / 0.4e1 + t18 * t286 / 0.2e1 - t18 * t290 / 0.2e1 + t18 * t319 / 0.4e1 );


    vrho_a = t4 * t85 + tzk0;
    vrho_b = t4 * t111 + tzk0;
    v2rho2_aa = t4 * t192 + 0.2e1 * t85;
    v2rho2_ab = t4 * t264 + t111 + t85;
    v2rho2_bb = t4 * t323 + 0.2e1 * t111;

  }


};

struct BuiltinEPC18_2 : detail::BuiltinKernelImpl< BuiltinEPC18_2 > {

  BuiltinEPC18_2( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinEPC18_2 >(p) { }
  
  virtual ~BuiltinEPC18_2() = default;

};



} // namespace ExchCXX
