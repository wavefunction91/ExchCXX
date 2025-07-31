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
struct kernel_traits< BuiltinEPC17_2 > :
  public lda_screening_interface< BuiltinEPC17_2 > {

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


  static constexpr double a = 2.35;
  static constexpr double b = 2.40;
  static constexpr double c = 6.60;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;


    eps = piecewise_functor_3( t2, 0.0, -t9 * t19 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;
    const double t23 = t18 * t18;
    const double t24 = 0.1e1 / t23;
    const double t26 = b / t12;
    const double t30 = t15 * rho / 0.2e1 - t26 * t9 / 0.2e1;
    const double t35 = piecewise_functor_3( t2, 0.0, t9 * t24 * t30 / 0.4e1 - t8 * t19 / 0.4e1 );


    eps = piecewise_functor_3( t2, 0.0, -t9 * t19 / 0.4e1 );
    vrho = rho * t35 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double& v2rho2 ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;
    const double t23 = t18 * t18;
    const double t24 = 0.1e1 / t23;
    const double t26 = b / t12;
    const double t30 = t15 * rho / 0.2e1 - t26 * t9 / 0.2e1;
    const double t35 = piecewise_functor_3( t2, 0.0, t9 * t24 * t30 / 0.4e1 - t8 * t19 / 0.4e1 );
    const double t38 = t8 * t24;
    const double t42 = 0.1e1 / t23 / t18;
    const double t43 = t30 * t30;
    const double t49 = b / t12 / t11;
    const double t50 = t8 * t8;
    const double t55 = t49 * t50 * t10 / 0.2e1 - t26 * t8 / 0.2e1 + t15 / 0.2e1;
    const double t60 = piecewise_functor_3( t2, 0.0, t38 * t30 / 0.2e1 - t9 * t42 * t43 / 0.2e1 + t9 * t24 * t55 / 0.4e1 );


    v2rho2 = rho * t60 + 0.2e1 * t35;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double& vrho, double& v2rho2 ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;
    const double tzk0 = piecewise_functor_3( t2, 0.0, -t9 * t19 / 0.4e1 );
    const double t23 = t18 * t18;
    const double t24 = 0.1e1 / t23;
    const double t26 = b / t12;
    const double t30 = t15 * rho / 0.2e1 - t26 * t9 / 0.2e1;
    const double t35 = piecewise_functor_3( t2, 0.0, t9 * t24 * t30 / 0.4e1 - t8 * t19 / 0.4e1 );
    const double t38 = t8 * t24;
    const double t42 = 0.1e1 / t23 / t18;
    const double t43 = t30 * t30;
    const double t49 = b / t12 / t11;
    const double t50 = t8 * t8;
    const double t55 = t49 * t50 * t10 / 0.2e1 - t26 * t8 / 0.2e1 + t15 / 0.2e1;
    const double t60 = piecewise_functor_3( t2, 0.0, t38 * t30 / 0.2e1 - t9 * t42 * t43 / 0.2e1 + t9 * t24 * t55 / 0.4e1 );


    vrho = rho * t35 + tzk0;
    v2rho2 = rho * t60 + 0.2e1 * t35;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {



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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;


    eps = piecewise_functor_3( t3, 0.0, -t18 * t35 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {



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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;
    const double t38 = 0.1e1 / t23;
    const double t39 = t14 * t38;
    const double t41 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t39 );
    const double t42 = t41 * t4;
    const double t44 = t17 * t22;
    const double t45 = t44 * t34;
    const double t46 = t19 * t38;
    const double t48 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t46 );
    const double t49 = t48 * t34;
    const double t51 = t33 * t33;
    const double t52 = 0.1e1 / t51;
    const double t53 = t22 * t52;
    const double t55 = b / t26;
    const double t56 = t41 * t23;
    const double t58 = t18 * t22;
    const double t59 = 0.2e1 * t58;
    const double t61 = t56 * t22 + t24 * t48 + t59;
    const double t64 = c * t41;
    const double t67 = t4 * t22;
    const double t69 = t29 * t67 / 0.2e1;
    const double t70 = t23 * t48;
    const double t73 = -t55 * t61 / 0.4e1 + t64 * t30 / 0.4e1 + t69 + t29 * t70 / 0.4e1;
    const double t74 = t53 * t73;
    const double t78 = piecewise_functor_3( t3, 0.0, -t18 * t49 / 0.4e1 + t18 * t74 / 0.4e1 - t42 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t81 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t39 );
    const double t82 = t81 * t4;
    const double t85 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t46 );
    const double t86 = t85 * t34;
    const double t88 = t81 * t23;
    const double t91 = t88 * t22 + t24 * t85 + t59;
    const double t94 = c * t81;
    const double t97 = t23 * t85;
    const double t100 = -t55 * t91 / 0.4e1 + t94 * t30 / 0.4e1 + t69 + t29 * t97 / 0.4e1;
    const double t101 = t53 * t100;
    const double t105 = piecewise_functor_3( t3, 0.0, t18 * t101 / 0.4e1 - t18 * t86 / 0.4e1 - t82 * t35 / 0.4e1 - t45 / 0.4e1 );


    eps = piecewise_functor_3( t3, 0.0, -t18 * t35 / 0.4e1 );
    vrho_a = t4 * t78 + eps;
    vrho_b = t4 * t105 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {



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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;
    const double t38 = 0.1e1 / t23;
    const double t39 = t14 * t38;
    const double t41 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t39 );
    const double t42 = t41 * t4;
    const double t44 = t17 * t22;
    const double t45 = t44 * t34;
    const double t46 = t19 * t38;
    const double t48 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t46 );
    const double t49 = t48 * t34;
    const double t51 = t33 * t33;
    const double t52 = 0.1e1 / t51;
    const double t53 = t22 * t52;
    const double t55 = b / t26;
    const double t56 = t41 * t23;
    const double t58 = t18 * t22;
    const double t59 = 0.2e1 * t58;
    const double t61 = t56 * t22 + t24 * t48 + t59;
    const double t64 = c * t41;
    const double t67 = t4 * t22;
    const double t69 = t29 * t67 / 0.2e1;
    const double t70 = t23 * t48;
    const double t73 = -t55 * t61 / 0.4e1 + t64 * t30 / 0.4e1 + t69 + t29 * t70 / 0.4e1;
    const double t74 = t53 * t73;
    const double t78 = piecewise_functor_3( t3, 0.0, -t18 * t49 / 0.4e1 + t18 * t74 / 0.4e1 - t42 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t81 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t39 );
    const double t82 = t81 * t4;
    const double t85 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t46 );
    const double t86 = t85 * t34;
    const double t88 = t81 * t23;
    const double t91 = t88 * t22 + t24 * t85 + t59;
    const double t94 = c * t81;
    const double t97 = t23 * t85;
    const double t100 = -t55 * t91 / 0.4e1 + t94 * t30 / 0.4e1 + t69 + t29 * t97 / 0.4e1;
    const double t101 = t53 * t100;
    const double t105 = piecewise_functor_3( t3, 0.0, t18 * t101 / 0.4e1 - t18 * t86 / 0.4e1 - t82 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t109 = 0.1e1 / t23 / t4;
    const double t110 = t14 * t109;
    const double t113 = piecewise_functor_5( t8, 0.0, t12, 0.0, -0.2e1 * t38 + 0.2e1 * t110 );
    const double t114 = t113 * t4;
    const double t117 = t41 * t22;
    const double t118 = t117 * t34;
    const double t124 = t17 * t48;
    const double t125 = t124 * t34;
    const double t127 = t52 * t73;
    const double t128 = t44 * t127;
    const double t130 = t19 * t109;
    const double t133 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t38 + 0.2e1 * t130 );
    const double t134 = t133 * t34;
    const double t137 = t48 * t52;
    const double t138 = t137 * t73;
    const double t142 = 0.1e1 / t51 / t33;
    const double t143 = t22 * t142;
    const double t144 = t73 * t73;
    const double t145 = t143 * t144;
    const double t150 = b / t26 / t25;
    const double t151 = t61 * t61;
    const double t154 = t113 * t23;
    const double t156 = t42 * t22;
    const double t160 = 0.2e1 * t44;
    const double t161 = t18 * t48;
    const double t164 = t24 * t133 + t154 * t22 + 0.2e1 * t56 * t48 + 0.4e1 * t156 + t160 + 0.4e1 * t161;
    const double t167 = c * t113;
    const double t170 = t64 * t67;
    const double t174 = t29 * t22 / 0.2e1;
    const double t175 = t4 * t48;
    const double t176 = t29 * t175;
    const double t177 = t23 * t133;
    const double t180 = t150 * t151 / 0.8e1 - t55 * t164 / 0.4e1 + t167 * t30 / 0.4e1 + t170 + t64 * t70 / 0.2e1 + t174 + t176 + t29 * t177 / 0.4e1;
    const double t181 = t53 * t180;
    const double t185 = piecewise_functor_3( t3, 0.0, -t114 * t35 / 0.4e1 - t118 / 0.2e1 - t42 * t49 / 0.2e1 + t42 * t74 / 0.2e1 - t125 / 0.2e1 + t128 / 0.2e1 - t18 * t134 / 0.4e1 + t18 * t138 / 0.2e1 - t18 * t145 / 0.2e1 + t18 * t181 / 0.4e1 );
    const double t188 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t110 );
    const double t189 = t188 * t4;
    const double t192 = t81 * t22;
    const double t193 = t192 * t34;
    const double t204 = t17 * t85;
    const double t205 = t204 * t34;
    const double t208 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t130 );
    const double t209 = t208 * t34;
    const double t212 = t85 * t52;
    const double t213 = t212 * t73;
    const double t218 = t52 * t100;
    const double t219 = t44 * t218;
    const double t221 = t137 * t100;
    const double t224 = t142 * t100;
    const double t225 = t224 * t73;
    const double t228 = t91 * t61;
    const double t231 = t188 * t23;
    const double t233 = t82 * t22;
    const double t239 = t18 * t85;
    const double t242 = t24 * t208 + t231 * t22 + t88 * t48 + t56 * t85 + 0.2e1 * t156 + t160 + 0.2e1 * t161 + 0.2e1 * t233 + 0.2e1 * t239;
    const double t245 = c * t188;
    const double t248 = t94 * t67;
    const double t256 = t4 * t85;
    const double t257 = t29 * t256;
    const double t259 = t23 * t208;
    const double t262 = t150 * t228 / 0.8e1 - t55 * t242 / 0.4e1 + t245 * t30 / 0.4e1 + t248 / 0.2e1 + t94 * t70 / 0.4e1 + t170 / 0.2e1 + t174 + t176 / 0.2e1 + t64 * t97 / 0.4e1 + t257 / 0.2e1 + t29 * t259 / 0.4e1;
    const double t263 = t53 * t262;
    const double t266 = -t189 * t35 / 0.4e1 - t193 / 0.4e1 - t82 * t49 / 0.4e1 + t82 * t74 / 0.4e1 - t118 / 0.4e1 - t125 / 0.4e1 + t128 / 0.4e1 - t42 * t86 / 0.4e1 - t205 / 0.4e1 - t18 * t209 / 0.4e1 + t18 * t213 / 0.4e1 + t42 * t101 / 0.4e1 + t219 / 0.4e1 + t18 * t221 / 0.4e1 - t58 * t225 / 0.2e1 + t18 * t263 / 0.4e1;
    const double t267 = piecewise_functor_3( t3, 0.0, t266 );
    const double t272 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t38 + 0.2e1 * t110 );
    const double t273 = t272 * t4;
    const double t285 = piecewise_functor_5( t12, 0.0, t8, 0.0, -0.2e1 * t38 + 0.2e1 * t130 );
    const double t286 = t285 * t34;
    const double t289 = t212 * t100;
    const double t292 = t100 * t100;
    const double t293 = t143 * t292;
    const double t296 = t91 * t91;
    const double t299 = t272 * t23;
    const double t306 = t299 * t22 + t24 * t285 + 0.2e1 * t88 * t85 + t160 + 0.4e1 * t233 + 0.4e1 * t239;
    const double t309 = c * t272;
    const double t314 = t23 * t285;
    const double t317 = t150 * t296 / 0.8e1 - t55 * t306 / 0.4e1 + t309 * t30 / 0.4e1 + t248 + t94 * t97 / 0.2e1 + t174 + t257 + t29 * t314 / 0.4e1;
    const double t318 = t53 * t317;
    const double t322 = piecewise_functor_3( t3, 0.0, -t273 * t35 / 0.4e1 - t193 / 0.2e1 - t82 * t86 / 0.2e1 + t82 * t101 / 0.2e1 - t205 / 0.2e1 + t219 / 0.2e1 - t18 * t286 / 0.4e1 + t18 * t289 / 0.2e1 - t18 * t293 / 0.2e1 + t18 * t318 / 0.4e1 );


    v2rho2_aa = t4 * t185 + 0.2e1 * t78;
    v2rho2_ab = t4 * t267 + t105 + t78;
    v2rho2_bb = t4 * t322 + 0.2e1 * t105;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {



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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;
    const double tzk0 = piecewise_functor_3( t3, 0.0, -t18 * t35 / 0.4e1 );
    const double t38 = 0.1e1 / t23;
    const double t39 = t14 * t38;
    const double t41 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t39 );
    const double t42 = t41 * t4;
    const double t44 = t17 * t22;
    const double t45 = t44 * t34;
    const double t46 = t19 * t38;
    const double t48 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t46 );
    const double t49 = t48 * t34;
    const double t51 = t33 * t33;
    const double t52 = 0.1e1 / t51;
    const double t53 = t22 * t52;
    const double t55 = b / t26;
    const double t56 = t41 * t23;
    const double t58 = t18 * t22;
    const double t59 = 0.2e1 * t58;
    const double t61 = t56 * t22 + t24 * t48 + t59;
    const double t64 = c * t41;
    const double t67 = t4 * t22;
    const double t69 = t29 * t67 / 0.2e1;
    const double t70 = t23 * t48;
    const double t73 = -t55 * t61 / 0.4e1 + t64 * t30 / 0.4e1 + t69 + t29 * t70 / 0.4e1;
    const double t74 = t53 * t73;
    const double t78 = piecewise_functor_3( t3, 0.0, -t18 * t49 / 0.4e1 + t18 * t74 / 0.4e1 - t42 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t81 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t39 );
    const double t82 = t81 * t4;
    const double t85 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t46 );
    const double t86 = t85 * t34;
    const double t88 = t81 * t23;
    const double t91 = t88 * t22 + t24 * t85 + t59;
    const double t94 = c * t81;
    const double t97 = t23 * t85;
    const double t100 = -t55 * t91 / 0.4e1 + t94 * t30 / 0.4e1 + t69 + t29 * t97 / 0.4e1;
    const double t101 = t53 * t100;
    const double t105 = piecewise_functor_3( t3, 0.0, t18 * t101 / 0.4e1 - t18 * t86 / 0.4e1 - t82 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t109 = 0.1e1 / t23 / t4;
    const double t110 = t14 * t109;
    const double t113 = piecewise_functor_5( t8, 0.0, t12, 0.0, -0.2e1 * t38 + 0.2e1 * t110 );
    const double t114 = t113 * t4;
    const double t117 = t41 * t22;
    const double t118 = t117 * t34;
    const double t124 = t17 * t48;
    const double t125 = t124 * t34;
    const double t127 = t52 * t73;
    const double t128 = t44 * t127;
    const double t130 = t19 * t109;
    const double t133 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t38 + 0.2e1 * t130 );
    const double t134 = t133 * t34;
    const double t137 = t48 * t52;
    const double t138 = t137 * t73;
    const double t142 = 0.1e1 / t51 / t33;
    const double t143 = t22 * t142;
    const double t144 = t73 * t73;
    const double t145 = t143 * t144;
    const double t150 = b / t26 / t25;
    const double t151 = t61 * t61;
    const double t154 = t113 * t23;
    const double t156 = t42 * t22;
    const double t160 = 0.2e1 * t44;
    const double t161 = t18 * t48;
    const double t164 = t24 * t133 + t154 * t22 + 0.2e1 * t56 * t48 + 0.4e1 * t156 + t160 + 0.4e1 * t161;
    const double t167 = c * t113;
    const double t170 = t64 * t67;
    const double t174 = t29 * t22 / 0.2e1;
    const double t175 = t4 * t48;
    const double t176 = t29 * t175;
    const double t177 = t23 * t133;
    const double t180 = t150 * t151 / 0.8e1 - t55 * t164 / 0.4e1 + t167 * t30 / 0.4e1 + t170 + t64 * t70 / 0.2e1 + t174 + t176 + t29 * t177 / 0.4e1;
    const double t181 = t53 * t180;
    const double t185 = piecewise_functor_3( t3, 0.0, -t114 * t35 / 0.4e1 - t118 / 0.2e1 - t42 * t49 / 0.2e1 + t42 * t74 / 0.2e1 - t125 / 0.2e1 + t128 / 0.2e1 - t18 * t134 / 0.4e1 + t18 * t138 / 0.2e1 - t18 * t145 / 0.2e1 + t18 * t181 / 0.4e1 );
    const double t188 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t110 );
    const double t189 = t188 * t4;
    const double t192 = t81 * t22;
    const double t193 = t192 * t34;
    const double t204 = t17 * t85;
    const double t205 = t204 * t34;
    const double t208 = piecewise_functor_5( t12, 0.0, t8, 0.0, 0.2e1 * t130 );
    const double t209 = t208 * t34;
    const double t212 = t85 * t52;
    const double t213 = t212 * t73;
    const double t218 = t52 * t100;
    const double t219 = t44 * t218;
    const double t221 = t137 * t100;
    const double t224 = t142 * t100;
    const double t225 = t224 * t73;
    const double t228 = t91 * t61;
    const double t231 = t188 * t23;
    const double t233 = t82 * t22;
    const double t239 = t18 * t85;
    const double t242 = t24 * t208 + t231 * t22 + t88 * t48 + t56 * t85 + 0.2e1 * t156 + t160 + 0.2e1 * t161 + 0.2e1 * t233 + 0.2e1 * t239;
    const double t245 = c * t188;
    const double t248 = t94 * t67;
    const double t256 = t4 * t85;
    const double t257 = t29 * t256;
    const double t259 = t23 * t208;
    const double t262 = t150 * t228 / 0.8e1 - t55 * t242 / 0.4e1 + t245 * t30 / 0.4e1 + t248 / 0.2e1 + t94 * t70 / 0.4e1 + t170 / 0.2e1 + t174 + t176 / 0.2e1 + t64 * t97 / 0.4e1 + t257 / 0.2e1 + t29 * t259 / 0.4e1;
    const double t263 = t53 * t262;
    const double t266 = -t189 * t35 / 0.4e1 - t193 / 0.4e1 - t82 * t49 / 0.4e1 + t82 * t74 / 0.4e1 - t118 / 0.4e1 - t125 / 0.4e1 + t128 / 0.4e1 - t42 * t86 / 0.4e1 - t205 / 0.4e1 - t18 * t209 / 0.4e1 + t18 * t213 / 0.4e1 + t42 * t101 / 0.4e1 + t219 / 0.4e1 + t18 * t221 / 0.4e1 - t58 * t225 / 0.2e1 + t18 * t263 / 0.4e1;
    const double t267 = piecewise_functor_3( t3, 0.0, t266 );
    const double t272 = piecewise_functor_5( t8, 0.0, t12, 0.0, 0.2e1 * t38 + 0.2e1 * t110 );
    const double t273 = t272 * t4;
    const double t285 = piecewise_functor_5( t12, 0.0, t8, 0.0, -0.2e1 * t38 + 0.2e1 * t130 );
    const double t286 = t285 * t34;
    const double t289 = t212 * t100;
    const double t292 = t100 * t100;
    const double t293 = t143 * t292;
    const double t296 = t91 * t91;
    const double t299 = t272 * t23;
    const double t306 = t299 * t22 + t24 * t285 + 0.2e1 * t88 * t85 + t160 + 0.4e1 * t233 + 0.4e1 * t239;
    const double t309 = c * t272;
    const double t314 = t23 * t285;
    const double t317 = t150 * t296 / 0.8e1 - t55 * t306 / 0.4e1 + t309 * t30 / 0.4e1 + t248 + t94 * t97 / 0.2e1 + t174 + t257 + t29 * t314 / 0.4e1;
    const double t318 = t53 * t317;
    const double t322 = piecewise_functor_3( t3, 0.0, -t273 * t35 / 0.4e1 - t193 / 0.2e1 - t82 * t86 / 0.2e1 + t82 * t101 / 0.2e1 - t205 / 0.2e1 + t219 / 0.2e1 - t18 * t286 / 0.4e1 + t18 * t289 / 0.2e1 - t18 * t293 / 0.2e1 + t18 * t318 / 0.4e1 );


    vrho_a = t4 * t78 + tzk0;
    vrho_b = t4 * t105 + tzk0;
    v2rho2_aa = t4 * t185 + 0.2e1 * t78;
    v2rho2_ab = t4 * t267 + t105 + t78;
    v2rho2_bb = t4 * t322 + 0.2e1 * t105;

  }


};

struct BuiltinEPC17_2 : detail::BuiltinKernelImpl< BuiltinEPC17_2 > {

  BuiltinEPC17_2( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinEPC17_2 >(p) { }
  
  virtual ~BuiltinEPC17_2() = default;

};



} // namespace ExchCXX
