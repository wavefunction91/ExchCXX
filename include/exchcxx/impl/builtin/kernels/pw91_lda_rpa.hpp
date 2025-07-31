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
struct kernel_traits< BuiltinPW91_LDA_RPA > :
  public lda_screening_interface< BuiltinPW91_LDA_RPA > {

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


  static constexpr double pp_0 = 0.75;
  static constexpr double pp_1 = 0.75;
  static constexpr double pp_2 = 1.0;
  static constexpr double a_0 = 0.031091;
  static constexpr double a_1 = 0.015545;
  static constexpr double a_2 = 0.016887;
  static constexpr double alpha1_0 = 0.082477;
  static constexpr double alpha1_1 = 0.035374;
  static constexpr double alpha1_2 = 0.028829;
  static constexpr double beta1_0 = 5.1486;
  static constexpr double beta1_1 = 6.4869;
  static constexpr double beta1_2 = 10.357;
  static constexpr double beta2_0 = 1.6483;
  static constexpr double beta2_1 = 1.3083;
  static constexpr double beta2_2 = 3.6231;
  static constexpr double beta3_0 = 0.23647;
  static constexpr double beta3_1 = 0.15180;
  static constexpr double beta3_2 = 0.47990;
  static constexpr double beta4_0 = 0.20614;
  static constexpr double beta4_1 = 0.082349;
  static constexpr double beta4_2 = 0.12279;
  static constexpr double fz20 = 1.709921;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t52 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_0;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_0 * t3;
    constexpr double t29 = beta3_0;
    constexpr double t36 = pp_0 + 0.1e1;
    constexpr double t57 = a_2;
    constexpr double t59 = alpha1_2;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = 0.1e1 / t57;
    constexpr double t65 = beta1_2;
    constexpr double t69 = beta2_2 * t3;
    constexpr double t72 = beta3_2;
    constexpr double t77 = pp_2 + 0.1e1;
    constexpr double t87 = 0.1e1 / fz20;


    const double t10 = safe_math::cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t15 = 0.1e1 + t4 * t12 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = safe_math::sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t34 = t21 / 0.4e1;
    const double t37 = safe_math::pow( t34, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.125e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = safe_math::log( t43 );
    const double t45 = t1 * t15 * t44;
    const double t47 = safe_math::cbrt( zeta_tol );
    const double t49 = piecewise_functor_3( 0.1e1 <= zeta_tol, t47 * zeta_tol, 1.0 );
    const double t56 = ( 0.2e1 * t49 - 0.2e1 ) / ( 0.2e1 * t52 - 0.2e1 );
    const double t63 = 0.1e1 + t60 * t12 / 0.4e1;
    const double t78 = safe_math::pow( t34, t77 );
    const double t79 = beta4_2 * t78;
    const double t80 = t65 * t22 / 0.2e1 + t69 * t12 / 0.4e1 + 0.125e0 * t72 * t30 + t79;
    const double t84 = 0.1e1 + t64 / t80 / 0.2e1;
    const double t85 = safe_math::log( t84 );
    const double t89 = t56 * t57 * t63 * t85 * t87;


    eps = -0.2e1 * t45 + 0.2e1 * t89;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t52 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_0;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_0 * t3;
    constexpr double t29 = beta3_0;
    constexpr double t36 = pp_0 + 0.1e1;
    constexpr double t57 = a_2;
    constexpr double t59 = alpha1_2;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = 0.1e1 / t57;
    constexpr double t65 = beta1_2;
    constexpr double t69 = beta2_2 * t3;
    constexpr double t72 = beta3_2;
    constexpr double t77 = pp_2 + 0.1e1;
    constexpr double t87 = 0.1e1 / fz20;
    constexpr double t94 = t1 * t2 * t3;


    const double t10 = safe_math::cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t15 = 0.1e1 + t4 * t12 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = safe_math::sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t34 = t21 / 0.4e1;
    const double t37 = safe_math::pow( t34, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.125e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = safe_math::log( t43 );
    const double t45 = t1 * t15 * t44;
    const double t47 = safe_math::cbrt( zeta_tol );
    const double t49 = piecewise_functor_3( 0.1e1 <= zeta_tol, t47 * zeta_tol, 1.0 );
    const double t56 = ( 0.2e1 * t49 - 0.2e1 ) / ( 0.2e1 * t52 - 0.2e1 );
    const double t63 = 0.1e1 + t60 * t12 / 0.4e1;
    const double t78 = safe_math::pow( t34, t77 );
    const double t79 = beta4_2 * t78;
    const double t80 = t65 * t22 / 0.2e1 + t69 * t12 / 0.4e1 + 0.125e0 * t72 * t30 + t79;
    const double t84 = 0.1e1 + t64 / t80 / 0.2e1;
    const double t85 = safe_math::log( t84 );
    const double t89 = t56 * t57 * t63 * t85 * t87;
    const double t96 = 0.1e1 / t10 / rho;
    const double t99 = t94 * t9 * t96 * t44;
    const double t101 = t39 * t39;
    const double t102 = 0.1e1 / t101;
    const double t103 = t15 * t102;
    const double t104 = 0.1e1 / t22;
    const double t106 = t18 * t104 * t3;
    const double t107 = t9 * t96;
    const double t112 = safe_math::sqrt( t21 );
    const double t114 = t29 * t112 * t3;
    const double t117 = 0.1e1 / rho;
    const double t121 = -t106 * t107 / 0.12e2 - t26 * t107 / 0.12e2 - 0.625e-1 * t114 * t107 - t38 * t36 * t117 / 0.3e1;
    const double t122 = 0.1e1 / t43;
    const double t123 = t121 * t122;
    const double t124 = t103 * t123;
    const double t127 = t56 * t57 * t59 * t3;
    const double t131 = t127 * t9 * t96 * t85 * t87;
    const double t133 = t56 * t63;
    const double t134 = t80 * t80;
    const double t135 = 0.1e1 / t134;
    const double t137 = t65 * t104 * t3;
    const double t143 = t72 * t112 * t3;
    const double t149 = -t137 * t107 / 0.12e2 - t69 * t107 / 0.12e2 - 0.625e-1 * t143 * t107 - t79 * t77 * t117 / 0.3e1;
    const double t151 = 0.1e1 / t84;
    const double t152 = t151 * t87;
    const double t154 = t133 * t135 * t149 * t152;


    eps = -0.2e1 * t45 + 0.2e1 * t89;
    vrho = -0.2e1 * t45 + 0.2e1 * t89 + rho * ( t99 / 0.6e1 + t124 - t131 / 0.6e1 - t154 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double& v2rho2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t52 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_0;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_0 * t3;
    constexpr double t29 = beta3_0;
    constexpr double t36 = pp_0 + 0.1e1;
    constexpr double t57 = a_2;
    constexpr double t59 = alpha1_2;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = 0.1e1 / t57;
    constexpr double t65 = beta1_2;
    constexpr double t69 = beta2_2 * t3;
    constexpr double t72 = beta3_2;
    constexpr double t77 = pp_2 + 0.1e1;
    constexpr double t87 = 0.1e1 / fz20;
    constexpr double t94 = t1 * t2 * t3;
    constexpr double t168 = t4 * t9;
    constexpr double t183 = t3 * t3;
    constexpr double t185 = t6 * t6;
    constexpr double t186 = t185 * t7;
    constexpr double t205 = t36 * t36;
    constexpr double t260 = t77 * t77;
    constexpr double t278 = t87 * t64;


    const double t10 = safe_math::cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t15 = 0.1e1 + t4 * t12 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = safe_math::sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t34 = t21 / 0.4e1;
    const double t37 = safe_math::pow( t34, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.125e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = safe_math::log( t43 );
    const double t47 = safe_math::cbrt( zeta_tol );
    const double t49 = piecewise_functor_3( 0.1e1 <= zeta_tol, t47 * zeta_tol, 1.0 );
    const double t56 = ( 0.2e1 * t49 - 0.2e1 ) / ( 0.2e1 * t52 - 0.2e1 );
    const double t63 = 0.1e1 + t60 * t12 / 0.4e1;
    const double t78 = safe_math::pow( t34, t77 );
    const double t79 = beta4_2 * t78;
    const double t80 = t65 * t22 / 0.2e1 + t69 * t12 / 0.4e1 + 0.125e0 * t72 * t30 + t79;
    const double t84 = 0.1e1 + t64 / t80 / 0.2e1;
    const double t85 = safe_math::log( t84 );
    const double t96 = 0.1e1 / t10 / rho;
    const double t99 = t94 * t9 * t96 * t44;
    const double t101 = t39 * t39;
    const double t102 = 0.1e1 / t101;
    const double t103 = t15 * t102;
    const double t104 = 0.1e1 / t22;
    const double t106 = t18 * t104 * t3;
    const double t107 = t9 * t96;
    const double t112 = safe_math::sqrt( t21 );
    const double t114 = t29 * t112 * t3;
    const double t117 = 0.1e1 / rho;
    const double t121 = -t106 * t107 / 0.12e2 - t26 * t107 / 0.12e2 - 0.625e-1 * t114 * t107 - t38 * t36 * t117 / 0.3e1;
    const double t122 = 0.1e1 / t43;
    const double t123 = t121 * t122;
    const double t124 = t103 * t123;
    const double t127 = t56 * t57 * t59 * t3;
    const double t131 = t127 * t9 * t96 * t85 * t87;
    const double t133 = t56 * t63;
    const double t134 = t80 * t80;
    const double t135 = 0.1e1 / t134;
    const double t137 = t65 * t104 * t3;
    const double t143 = t72 * t112 * t3;
    const double t149 = -t137 * t107 / 0.12e2 - t69 * t107 / 0.12e2 - 0.625e-1 * t143 * t107 - t79 * t77 * t117 / 0.3e1;
    const double t151 = 0.1e1 / t84;
    const double t152 = t151 * t87;
    const double t154 = t133 * t135 * t149 * t152;
    const double t161 = rho * rho;
    const double t163 = 0.1e1 / t10 / t161;
    const double t166 = t94 * t9 * t163 * t44;
    const double t169 = t96 * t102;
    const double t171 = t168 * t169 * t123;
    const double t173 = t101 * t39;
    const double t174 = 0.1e1 / t173;
    const double t175 = t15 * t174;
    const double t176 = t121 * t121;
    const double t177 = t176 * t122;
    const double t178 = t175 * t177;
    const double t181 = 0.1e1 / t22 / t21;
    const double t184 = t18 * t181 * t183;
    const double t187 = t10 * t10;
    const double t190 = t186 / t187 / t161;
    const double t193 = t9 * t163;
    const double t198 = 0.1e1/safe_math::sqrt( t21 );
    const double t200 = t29 * t198 * t183;
    const double t206 = 0.1e1 / t161;
    const double t213 = -t184 * t190 / 0.18e2 + t106 * t193 / 0.9e1 + t26 * t193 / 0.9e1 + 0.41666666666666666666e-1 * t200 * t190 + 0.83333333333333333333e-1 * t114 * t193 + t38 * t205 * t206 / 0.9e1 + t38 * t36 * t206 / 0.3e1;
    const double t214 = t213 * t122;
    const double t215 = t103 * t214;
    const double t216 = t101 * t101;
    const double t217 = 0.1e1 / t216;
    const double t218 = t15 * t217;
    const double t219 = t43 * t43;
    const double t220 = 0.1e1 / t219;
    const double t222 = t176 * t220 * t17;
    const double t223 = t218 * t222;
    const double t228 = t127 * t9 * t163 * t85 * t87;
    const double t231 = t56 * t60 * t6;
    const double t232 = t8 * t96;
    const double t233 = t232 * t135;
    const double t234 = t149 * t151;
    const double t235 = t234 * t87;
    const double t237 = t231 * t233 * t235;
    const double t239 = t134 * t80;
    const double t240 = 0.1e1 / t239;
    const double t241 = t149 * t149;
    const double t244 = t133 * t240 * t241 * t152;
    const double t247 = t65 * t181 * t183;
    const double t255 = t72 * t198 * t183;
    const double t267 = -t247 * t190 / 0.18e2 + t137 * t193 / 0.9e1 + t69 * t193 / 0.9e1 + 0.41666666666666666666e-1 * t255 * t190 + 0.83333333333333333333e-1 * t143 * t193 + t79 * t260 * t206 / 0.9e1 + t79 * t77 * t206 / 0.3e1;
    const double t270 = t133 * t135 * t267 * t152;
    const double t271 = t134 * t134;
    const double t272 = 0.1e1 / t271;
    const double t274 = t56 * t63 * t272;
    const double t275 = t84 * t84;
    const double t276 = 0.1e1 / t275;
    const double t277 = t241 * t276;
    const double t280 = t274 * t277 * t278;


    v2rho2 = t99 / 0.3e1 + 0.2e1 * t124 - t131 / 0.3e1 - 0.2e1 * t154 + rho * ( -0.2e1 / 0.9e1 * t166 - t171 / 0.6e1 - 0.2e1 * t178 + t215 + t223 / 0.2e1 + 0.2e1 / 0.9e1 * t228 + t237 / 0.6e1 + 0.2e1 * t244 - t270 - t280 / 0.2e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double& vrho, double& v2rho2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t52 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_0;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_0 * t3;
    constexpr double t29 = beta3_0;
    constexpr double t36 = pp_0 + 0.1e1;
    constexpr double t57 = a_2;
    constexpr double t59 = alpha1_2;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = 0.1e1 / t57;
    constexpr double t65 = beta1_2;
    constexpr double t69 = beta2_2 * t3;
    constexpr double t72 = beta3_2;
    constexpr double t77 = pp_2 + 0.1e1;
    constexpr double t87 = 0.1e1 / fz20;
    constexpr double t94 = t1 * t2 * t3;
    constexpr double t168 = t4 * t9;
    constexpr double t183 = t3 * t3;
    constexpr double t185 = t6 * t6;
    constexpr double t186 = t185 * t7;
    constexpr double t205 = t36 * t36;
    constexpr double t260 = t77 * t77;
    constexpr double t278 = t87 * t64;


    const double t10 = safe_math::cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t15 = 0.1e1 + t4 * t12 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = safe_math::sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t34 = t21 / 0.4e1;
    const double t37 = safe_math::pow( t34, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.125e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = safe_math::log( t43 );
    const double t45 = t1 * t15 * t44;
    const double t47 = safe_math::cbrt( zeta_tol );
    const double t49 = piecewise_functor_3( 0.1e1 <= zeta_tol, t47 * zeta_tol, 1.0 );
    const double t56 = ( 0.2e1 * t49 - 0.2e1 ) / ( 0.2e1 * t52 - 0.2e1 );
    const double t63 = 0.1e1 + t60 * t12 / 0.4e1;
    const double t78 = safe_math::pow( t34, t77 );
    const double t79 = beta4_2 * t78;
    const double t80 = t65 * t22 / 0.2e1 + t69 * t12 / 0.4e1 + 0.125e0 * t72 * t30 + t79;
    const double t84 = 0.1e1 + t64 / t80 / 0.2e1;
    const double t85 = safe_math::log( t84 );
    const double t89 = t56 * t57 * t63 * t85 * t87;
    const double t96 = 0.1e1 / t10 / rho;
    const double t99 = t94 * t9 * t96 * t44;
    const double t101 = t39 * t39;
    const double t102 = 0.1e1 / t101;
    const double t103 = t15 * t102;
    const double t104 = 0.1e1 / t22;
    const double t106 = t18 * t104 * t3;
    const double t107 = t9 * t96;
    const double t112 = safe_math::sqrt( t21 );
    const double t114 = t29 * t112 * t3;
    const double t117 = 0.1e1 / rho;
    const double t121 = -t106 * t107 / 0.12e2 - t26 * t107 / 0.12e2 - 0.625e-1 * t114 * t107 - t38 * t36 * t117 / 0.3e1;
    const double t122 = 0.1e1 / t43;
    const double t123 = t121 * t122;
    const double t124 = t103 * t123;
    const double t127 = t56 * t57 * t59 * t3;
    const double t131 = t127 * t9 * t96 * t85 * t87;
    const double t133 = t56 * t63;
    const double t134 = t80 * t80;
    const double t135 = 0.1e1 / t134;
    const double t137 = t65 * t104 * t3;
    const double t143 = t72 * t112 * t3;
    const double t149 = -t137 * t107 / 0.12e2 - t69 * t107 / 0.12e2 - 0.625e-1 * t143 * t107 - t79 * t77 * t117 / 0.3e1;
    const double t151 = 0.1e1 / t84;
    const double t152 = t151 * t87;
    const double t154 = t133 * t135 * t149 * t152;
    const double t161 = rho * rho;
    const double t163 = 0.1e1 / t10 / t161;
    const double t166 = t94 * t9 * t163 * t44;
    const double t169 = t96 * t102;
    const double t171 = t168 * t169 * t123;
    const double t173 = t101 * t39;
    const double t174 = 0.1e1 / t173;
    const double t175 = t15 * t174;
    const double t176 = t121 * t121;
    const double t177 = t176 * t122;
    const double t178 = t175 * t177;
    const double t181 = 0.1e1 / t22 / t21;
    const double t184 = t18 * t181 * t183;
    const double t187 = t10 * t10;
    const double t190 = t186 / t187 / t161;
    const double t193 = t9 * t163;
    const double t198 = 0.1e1/safe_math::sqrt( t21 );
    const double t200 = t29 * t198 * t183;
    const double t206 = 0.1e1 / t161;
    const double t213 = -t184 * t190 / 0.18e2 + t106 * t193 / 0.9e1 + t26 * t193 / 0.9e1 + 0.41666666666666666666e-1 * t200 * t190 + 0.83333333333333333333e-1 * t114 * t193 + t38 * t205 * t206 / 0.9e1 + t38 * t36 * t206 / 0.3e1;
    const double t214 = t213 * t122;
    const double t215 = t103 * t214;
    const double t216 = t101 * t101;
    const double t217 = 0.1e1 / t216;
    const double t218 = t15 * t217;
    const double t219 = t43 * t43;
    const double t220 = 0.1e1 / t219;
    const double t222 = t176 * t220 * t17;
    const double t223 = t218 * t222;
    const double t228 = t127 * t9 * t163 * t85 * t87;
    const double t231 = t56 * t60 * t6;
    const double t232 = t8 * t96;
    const double t233 = t232 * t135;
    const double t234 = t149 * t151;
    const double t235 = t234 * t87;
    const double t237 = t231 * t233 * t235;
    const double t239 = t134 * t80;
    const double t240 = 0.1e1 / t239;
    const double t241 = t149 * t149;
    const double t244 = t133 * t240 * t241 * t152;
    const double t247 = t65 * t181 * t183;
    const double t255 = t72 * t198 * t183;
    const double t267 = -t247 * t190 / 0.18e2 + t137 * t193 / 0.9e1 + t69 * t193 / 0.9e1 + 0.41666666666666666666e-1 * t255 * t190 + 0.83333333333333333333e-1 * t143 * t193 + t79 * t260 * t206 / 0.9e1 + t79 * t77 * t206 / 0.3e1;
    const double t270 = t133 * t135 * t267 * t152;
    const double t271 = t134 * t134;
    const double t272 = 0.1e1 / t271;
    const double t274 = t56 * t63 * t272;
    const double t275 = t84 * t84;
    const double t276 = 0.1e1 / t275;
    const double t277 = t241 * t276;
    const double t280 = t274 * t277 * t278;


    vrho = -0.2e1 * t45 + 0.2e1 * t89 + rho * ( t99 / 0.6e1 + t124 - t131 / 0.6e1 - t154 );
    v2rho2 = t99 / 0.3e1 + 0.2e1 * t124 - t131 / 0.3e1 - 0.2e1 * t154 + rho * ( -0.2e1 / 0.9e1 * t166 - t171 / 0.6e1 - 0.2e1 * t178 + t215 + t223 / 0.2e1 + 0.2e1 / 0.9e1 * t228 + t237 / 0.6e1 + 0.2e1 * t244 - t270 - t280 / 0.2e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t70 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t18 = 0.1e1 / t1;
    constexpr double t19 = beta1_0;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = beta2_0 * t3;
    constexpr double t30 = beta3_0;
    constexpr double t37 = pp_0 + 0.1e1;
    constexpr double t75 = a_1;
    constexpr double t76 = alpha1_1;
    constexpr double t77 = t76 * t3;
    constexpr double t82 = 0.1e1 / t75;
    constexpr double t83 = beta1_1;
    constexpr double t87 = beta2_1 * t3;
    constexpr double t90 = beta3_1;
    constexpr double t95 = pp_1 + 0.1e1;
    constexpr double t105 = a_2;
    constexpr double t106 = alpha1_2;
    constexpr double t107 = t106 * t3;
    constexpr double t112 = 0.1e1 / t105;
    constexpr double t113 = beta1_2;
    constexpr double t117 = beta2_2 * t3;
    constexpr double t120 = beta3_2;
    constexpr double t125 = pp_2 + 0.1e1;
    constexpr double t134 = 0.1e1 / fz20;


    const double t10 = rho_a + rho_b;
    const double t11 = safe_math::cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = safe_math::sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = safe_math::pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.125e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = safe_math::log( t44 );
    const double t46 = t1 * t16 * t45;
    const double t47 = 0.2e1 * t46;
    const double t48 = rho_a - rho_b;
    const double t49 = t48 * t48;
    const double t50 = t49 * t49;
    const double t51 = t10 * t10;
    const double t52 = t51 * t51;
    const double t53 = 0.1e1 / t52;
    const double t54 = t50 * t53;
    const double t55 = 0.1e1 / t10;
    const double t56 = t48 * t55;
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( zeta_tol );
    const double t60 = t59 * zeta_tol;
    const double t61 = safe_math::cbrt( t57 );
    const double t63 = piecewise_functor_3( t58, t60, t61 * t57 );
    const double t64 = 0.1e1 - t56;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t60, t66 * t64 );
    const double t69 = t63 + t68 - 0.2e1;
    const double t73 = 0.1e1 / ( 0.2e1 * t70 - 0.2e1 );
    const double t74 = t69 * t73;
    const double t80 = 0.1e1 + t77 * t13 / 0.4e1;
    const double t96 = safe_math::pow( t35, t95 );
    const double t97 = beta4_1 * t96;
    const double t98 = t83 * t23 / 0.2e1 + t87 * t13 / 0.4e1 + 0.125e0 * t90 * t31 + t97;
    const double t102 = 0.1e1 + t82 / t98 / 0.2e1;
    const double t103 = safe_math::log( t102 );
    const double t110 = 0.1e1 + t107 * t13 / 0.4e1;
    const double t126 = safe_math::pow( t35, t125 );
    const double t127 = beta4_2 * t126;
    const double t128 = t113 * t23 / 0.2e1 + t117 * t13 / 0.4e1 + 0.125e0 * t120 * t31 + t127;
    const double t132 = 0.1e1 + t112 / t128 / 0.2e1;
    const double t133 = safe_math::log( t132 );
    const double t135 = t133 * t134;
    const double t138 = -0.2e1 * t75 * t80 * t103 - 0.2e1 * t105 * t110 * t135 + 0.2e1 * t46;
    const double t139 = t74 * t138;
    const double t140 = t54 * t139;
    const double t143 = t110 * t133 * t134;
    const double t145 = 0.2e1 * t74 * t105 * t143;


    eps = -t47 + t140 + t145;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t70 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t18 = 0.1e1 / t1;
    constexpr double t19 = beta1_0;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = beta2_0 * t3;
    constexpr double t30 = beta3_0;
    constexpr double t37 = pp_0 + 0.1e1;
    constexpr double t75 = a_1;
    constexpr double t76 = alpha1_1;
    constexpr double t77 = t76 * t3;
    constexpr double t82 = 0.1e1 / t75;
    constexpr double t83 = beta1_1;
    constexpr double t87 = beta2_1 * t3;
    constexpr double t90 = beta3_1;
    constexpr double t95 = pp_1 + 0.1e1;
    constexpr double t105 = a_2;
    constexpr double t106 = alpha1_2;
    constexpr double t107 = t106 * t3;
    constexpr double t112 = 0.1e1 / t105;
    constexpr double t113 = beta1_2;
    constexpr double t117 = beta2_2 * t3;
    constexpr double t120 = beta3_2;
    constexpr double t125 = pp_2 + 0.1e1;
    constexpr double t134 = 0.1e1 / fz20;
    constexpr double t147 = t1 * t2 * t3;
    constexpr double t201 = t75 * t76 * t3;
    constexpr double t226 = t105 * t106;
    constexpr double t227 = t226 * t20;
    constexpr double t259 = t226 * t3;


    const double t10 = rho_a + rho_b;
    const double t11 = safe_math::cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = safe_math::sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = safe_math::pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.125e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = safe_math::log( t44 );
    const double t46 = t1 * t16 * t45;
    const double t47 = 0.2e1 * t46;
    const double t48 = rho_a - rho_b;
    const double t49 = t48 * t48;
    const double t50 = t49 * t49;
    const double t51 = t10 * t10;
    const double t52 = t51 * t51;
    const double t53 = 0.1e1 / t52;
    const double t54 = t50 * t53;
    const double t55 = 0.1e1 / t10;
    const double t56 = t48 * t55;
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( zeta_tol );
    const double t60 = t59 * zeta_tol;
    const double t61 = safe_math::cbrt( t57 );
    const double t63 = piecewise_functor_3( t58, t60, t61 * t57 );
    const double t64 = 0.1e1 - t56;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t60, t66 * t64 );
    const double t69 = t63 + t68 - 0.2e1;
    const double t73 = 0.1e1 / ( 0.2e1 * t70 - 0.2e1 );
    const double t74 = t69 * t73;
    const double t80 = 0.1e1 + t77 * t13 / 0.4e1;
    const double t96 = safe_math::pow( t35, t95 );
    const double t97 = beta4_1 * t96;
    const double t98 = t83 * t23 / 0.2e1 + t87 * t13 / 0.4e1 + 0.125e0 * t90 * t31 + t97;
    const double t102 = 0.1e1 + t82 / t98 / 0.2e1;
    const double t103 = safe_math::log( t102 );
    const double t110 = 0.1e1 + t107 * t13 / 0.4e1;
    const double t126 = safe_math::pow( t35, t125 );
    const double t127 = beta4_2 * t126;
    const double t128 = t113 * t23 / 0.2e1 + t117 * t13 / 0.4e1 + 0.125e0 * t120 * t31 + t127;
    const double t132 = 0.1e1 + t112 / t128 / 0.2e1;
    const double t133 = safe_math::log( t132 );
    const double t135 = t133 * t134;
    const double t138 = -0.2e1 * t75 * t80 * t103 - 0.2e1 * t105 * t110 * t135 + 0.2e1 * t46;
    const double t139 = t74 * t138;
    const double t140 = t54 * t139;
    const double t143 = t110 * t133 * t134;
    const double t145 = 0.2e1 * t74 * t105 * t143;
    const double t149 = 0.1e1 / t11 / t10;
    const double t152 = t147 * t9 * t149 * t45;
    const double t153 = t152 / 0.6e1;
    const double t154 = t40 * t40;
    const double t155 = 0.1e1 / t154;
    const double t156 = t16 * t155;
    const double t157 = 0.1e1 / t23;
    const double t159 = t19 * t157 * t3;
    const double t160 = t9 * t149;
    const double t165 = safe_math::sqrt( t22 );
    const double t167 = t30 * t165 * t3;
    const double t173 = -t159 * t160 / 0.12e2 - t27 * t160 / 0.12e2 - 0.625e-1 * t167 * t160 - t39 * t37 * t55 / 0.3e1;
    const double t174 = 0.1e1 / t44;
    const double t175 = t173 * t174;
    const double t176 = t156 * t175;
    const double t177 = t49 * t48;
    const double t178 = t177 * t53;
    const double t179 = t178 * t139;
    const double t180 = 0.4e1 * t179;
    const double t181 = t52 * t10;
    const double t182 = 0.1e1 / t181;
    const double t183 = t50 * t182;
    const double t184 = t183 * t139;
    const double t185 = 0.4e1 * t184;
    const double t186 = 0.1e1 / t51;
    const double t187 = t48 * t186;
    const double t188 = t55 - t187;
    const double t191 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t188 );
    const double t192 = -t188;
    const double t195 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t192 );
    const double t197 = ( t191 + t195 ) * t73;
    const double t198 = t197 * t138;
    const double t199 = t54 * t198;
    const double t206 = t98 * t98;
    const double t207 = 0.1e1 / t206;
    const double t208 = t80 * t207;
    const double t210 = t83 * t157 * t3;
    const double t216 = t90 * t165 * t3;
    const double t222 = -t210 * t160 / 0.12e2 - t87 * t160 / 0.12e2 - 0.625e-1 * t216 * t160 - t97 * t95 * t55 / 0.3e1;
    const double t223 = 0.1e1 / t102;
    const double t224 = t222 * t223;
    const double t228 = t8 * t149;
    const double t232 = t128 * t128;
    const double t233 = 0.1e1 / t232;
    const double t234 = t110 * t233;
    const double t236 = t113 * t157 * t3;
    const double t242 = t120 * t165 * t3;
    const double t248 = -t236 * t160 / 0.12e2 - t117 * t160 / 0.12e2 - 0.625e-1 * t242 * t160 - t127 * t125 * t55 / 0.3e1;
    const double t249 = 0.1e1 / t132;
    const double t251 = t248 * t249 * t134;
    const double t253 = t201 * t9 * t149 * t103 / 0.6e1 + t208 * t224 - t153 - t176 + t227 * t228 * t135 / 0.6e1 + t234 * t251;
    const double t254 = t74 * t253;
    const double t255 = t54 * t254;
    const double t257 = t197 * t105 * t143;
    const double t258 = 0.2e1 * t257;
    const double t260 = t74 * t259;
    const double t263 = t9 * t149 * t133 * t134;
    const double t264 = t260 * t263;
    const double t265 = t264 / 0.6e1;
    const double t266 = t74 * t110;
    const double t268 = t249 * t134;
    const double t269 = t233 * t248 * t268;
    const double t270 = t266 * t269;
    const double t273 = -t55 - t187;
    const double t276 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t273 );
    const double t277 = -t273;
    const double t280 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t277 );
    const double t282 = ( t276 + t280 ) * t73;
    const double t283 = t282 * t138;
    const double t284 = t54 * t283;
    const double t286 = t282 * t105 * t143;
    const double t287 = 0.2e1 * t286;


    eps = -t47 + t140 + t145;
    vrho_a = -t47 + t140 + t145 + t10 * ( t153 + t176 + t180 - t185 + t199 + t255 + t258 - t265 - t270 );
    vrho_b = -t47 + t140 + t145 + t10 * ( t153 + t176 - t180 - t185 + t284 + t255 + t287 - t265 - t270 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t70 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t18 = 0.1e1 / t1;
    constexpr double t19 = beta1_0;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = beta2_0 * t3;
    constexpr double t30 = beta3_0;
    constexpr double t37 = pp_0 + 0.1e1;
    constexpr double t75 = a_1;
    constexpr double t76 = alpha1_1;
    constexpr double t77 = t76 * t3;
    constexpr double t82 = 0.1e1 / t75;
    constexpr double t83 = beta1_1;
    constexpr double t87 = beta2_1 * t3;
    constexpr double t90 = beta3_1;
    constexpr double t95 = pp_1 + 0.1e1;
    constexpr double t105 = a_2;
    constexpr double t106 = alpha1_2;
    constexpr double t107 = t106 * t3;
    constexpr double t112 = 0.1e1 / t105;
    constexpr double t113 = beta1_2;
    constexpr double t117 = beta2_2 * t3;
    constexpr double t120 = beta3_2;
    constexpr double t125 = pp_2 + 0.1e1;
    constexpr double t134 = 0.1e1 / fz20;
    constexpr double t147 = t1 * t2 * t3;
    constexpr double t201 = t75 * t76 * t3;
    constexpr double t226 = t105 * t106;
    constexpr double t227 = t226 * t20;
    constexpr double t259 = t226 * t3;
    constexpr double t302 = t3 * t3;
    constexpr double t304 = t6 * t6;
    constexpr double t305 = t304 * t7;
    constexpr double t326 = t37 * t37;
    constexpr double t364 = t125 * t125;
    constexpr double t435 = t77 * t9;
    constexpr double t461 = t95 * t95;
    constexpr double t484 = t4 * t9;
    constexpr double t493 = t107 * t9;
    constexpr double t528 = t134 * t112;
    constexpr double t543 = t107 * t6;


    const double t10 = rho_a + rho_b;
    const double t11 = safe_math::cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = safe_math::sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = safe_math::pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.125e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = safe_math::log( t44 );
    const double t46 = t1 * t16 * t45;
    const double t48 = rho_a - rho_b;
    const double t49 = t48 * t48;
    const double t50 = t49 * t49;
    const double t51 = t10 * t10;
    const double t52 = t51 * t51;
    const double t53 = 0.1e1 / t52;
    const double t54 = t50 * t53;
    const double t55 = 0.1e1 / t10;
    const double t56 = t48 * t55;
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( zeta_tol );
    const double t60 = t59 * zeta_tol;
    const double t61 = safe_math::cbrt( t57 );
    const double t63 = piecewise_functor_3( t58, t60, t61 * t57 );
    const double t64 = 0.1e1 - t56;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t60, t66 * t64 );
    const double t69 = t63 + t68 - 0.2e1;
    const double t73 = 0.1e1 / ( 0.2e1 * t70 - 0.2e1 );
    const double t74 = t69 * t73;
    const double t80 = 0.1e1 + t77 * t13 / 0.4e1;
    const double t96 = safe_math::pow( t35, t95 );
    const double t97 = beta4_1 * t96;
    const double t98 = t83 * t23 / 0.2e1 + t87 * t13 / 0.4e1 + 0.125e0 * t90 * t31 + t97;
    const double t102 = 0.1e1 + t82 / t98 / 0.2e1;
    const double t103 = safe_math::log( t102 );
    const double t110 = 0.1e1 + t107 * t13 / 0.4e1;
    const double t126 = safe_math::pow( t35, t125 );
    const double t127 = beta4_2 * t126;
    const double t128 = t113 * t23 / 0.2e1 + t117 * t13 / 0.4e1 + 0.125e0 * t120 * t31 + t127;
    const double t132 = 0.1e1 + t112 / t128 / 0.2e1;
    const double t133 = safe_math::log( t132 );
    const double t135 = t133 * t134;
    const double t138 = -0.2e1 * t75 * t80 * t103 - 0.2e1 * t105 * t110 * t135 + 0.2e1 * t46;
    const double t139 = t74 * t138;
    const double t143 = t110 * t133 * t134;
    const double t149 = 0.1e1 / t11 / t10;
    const double t152 = t147 * t9 * t149 * t45;
    const double t153 = t152 / 0.6e1;
    const double t154 = t40 * t40;
    const double t155 = 0.1e1 / t154;
    const double t156 = t16 * t155;
    const double t157 = 0.1e1 / t23;
    const double t159 = t19 * t157 * t3;
    const double t160 = t9 * t149;
    const double t165 = safe_math::sqrt( t22 );
    const double t167 = t30 * t165 * t3;
    const double t173 = -t159 * t160 / 0.12e2 - t27 * t160 / 0.12e2 - 0.625e-1 * t167 * t160 - t39 * t37 * t55 / 0.3e1;
    const double t174 = 0.1e1 / t44;
    const double t175 = t173 * t174;
    const double t176 = t156 * t175;
    const double t177 = t49 * t48;
    const double t178 = t177 * t53;
    const double t179 = t178 * t139;
    const double t181 = t52 * t10;
    const double t182 = 0.1e1 / t181;
    const double t183 = t50 * t182;
    const double t184 = t183 * t139;
    const double t186 = 0.1e1 / t51;
    const double t187 = t48 * t186;
    const double t188 = t55 - t187;
    const double t191 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t188 );
    const double t192 = -t188;
    const double t195 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t192 );
    const double t197 = ( t191 + t195 ) * t73;
    const double t198 = t197 * t138;
    const double t199 = t54 * t198;
    const double t206 = t98 * t98;
    const double t207 = 0.1e1 / t206;
    const double t208 = t80 * t207;
    const double t210 = t83 * t157 * t3;
    const double t216 = t90 * t165 * t3;
    const double t222 = -t210 * t160 / 0.12e2 - t87 * t160 / 0.12e2 - 0.625e-1 * t216 * t160 - t97 * t95 * t55 / 0.3e1;
    const double t223 = 0.1e1 / t102;
    const double t224 = t222 * t223;
    const double t228 = t8 * t149;
    const double t232 = t128 * t128;
    const double t233 = 0.1e1 / t232;
    const double t234 = t110 * t233;
    const double t236 = t113 * t157 * t3;
    const double t242 = t120 * t165 * t3;
    const double t248 = -t236 * t160 / 0.12e2 - t117 * t160 / 0.12e2 - 0.625e-1 * t242 * t160 - t127 * t125 * t55 / 0.3e1;
    const double t249 = 0.1e1 / t132;
    const double t251 = t248 * t249 * t134;
    const double t253 = t201 * t9 * t149 * t103 / 0.6e1 + t208 * t224 - t153 - t176 + t227 * t228 * t135 / 0.6e1 + t234 * t251;
    const double t254 = t74 * t253;
    const double t255 = t54 * t254;
    const double t257 = t197 * t105 * t143;
    const double t258 = 0.2e1 * t257;
    const double t260 = t74 * t259;
    const double t263 = t9 * t149 * t133 * t134;
    const double t264 = t260 * t263;
    const double t266 = t74 * t110;
    const double t268 = t249 * t134;
    const double t269 = t233 * t248 * t268;
    const double t270 = t266 * t269;
    const double t273 = -t55 - t187;
    const double t276 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t273 );
    const double t277 = -t273;
    const double t280 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t277 );
    const double t282 = ( t276 + t280 ) * t73;
    const double t283 = t282 * t138;
    const double t284 = t54 * t283;
    const double t286 = t282 * t105 * t143;
    const double t287 = 0.2e1 * t286;
    const double t290 = t152 / 0.3e1;
    const double t291 = 0.2e1 * t176;
    const double t292 = 0.8e1 * t179;
    const double t293 = 0.8e1 * t184;
    const double t295 = 0.2e1 * t255;
    const double t297 = t264 / 0.3e1;
    const double t298 = 0.2e1 * t270;
    const double t300 = 0.1e1 / t23 / t22;
    const double t303 = t19 * t300 * t302;
    const double t306 = t11 * t11;
    const double t309 = t305 / t306 / t51;
    const double t313 = 0.1e1 / t11 / t51;
    const double t314 = t9 * t313;
    const double t319 = 0.1e1/safe_math::sqrt( t22 );
    const double t321 = t30 * t319 * t302;
    const double t333 = -t303 * t309 / 0.18e2 + t159 * t314 / 0.9e1 + t27 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t321 * t309 + 0.83333333333333333333e-1 * t167 * t314 + t39 * t326 * t186 / 0.9e1 + t39 * t37 * t186 / 0.3e1;
    const double t334 = t333 * t174;
    const double t335 = t156 * t334;
    const double t336 = t49 * t53;
    const double t337 = t336 * t139;
    const double t338 = 0.12e2 * t337;
    const double t339 = t177 * t182;
    const double t340 = t339 * t139;
    const double t341 = 0.32e2 * t340;
    const double t343 = 0.1e1 / t52 / t51;
    const double t344 = t50 * t343;
    const double t345 = t344 * t139;
    const double t346 = 0.2e2 * t345;
    const double t347 = t197 * t110;
    const double t348 = t347 * t269;
    const double t349 = 0.2e1 * t348;
    const double t351 = t113 * t300 * t302;
    const double t359 = t120 * t319 * t302;
    const double t371 = -t351 * t309 / 0.18e2 + t236 * t314 / 0.9e1 + t117 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t359 * t309 + 0.83333333333333333333e-1 * t242 * t314 + t127 * t364 * t186 / 0.9e1 + t127 * t125 * t186 / 0.3e1;
    const double t373 = t233 * t371 * t268;
    const double t374 = t266 * t373;
    const double t375 = t154 * t40;
    const double t376 = 0.1e1 / t375;
    const double t377 = t16 * t376;
    const double t378 = t173 * t173;
    const double t379 = t378 * t174;
    const double t380 = t377 * t379;
    const double t381 = 0.2e1 * t380;
    const double t382 = t154 * t154;
    const double t383 = 0.1e1 / t382;
    const double t384 = t16 * t383;
    const double t385 = t44 * t44;
    const double t386 = 0.1e1 / t385;
    const double t388 = t378 * t386 * t18;
    const double t389 = t384 * t388;
    const double t390 = t389 / 0.2e1;
    const double t391 = t178 * t198;
    const double t392 = 0.8e1 * t391;
    const double t393 = t178 * t254;
    const double t394 = 0.8e1 * t393;
    const double t395 = t183 * t198;
    const double t396 = 0.8e1 * t395;
    const double t397 = t335 + t338 - t341 + t346 - t349 - t374 - t381 + t390 + t392 + t394 - t396;
    const double t398 = t183 * t254;
    const double t399 = 0.8e1 * t398;
    const double t400 = t61 * t61;
    const double t401 = 0.1e1 / t400;
    const double t402 = t188 * t188;
    const double t405 = t51 * t10;
    const double t406 = 0.1e1 / t405;
    const double t407 = t48 * t406;
    const double t409 = -0.2e1 * t186 + 0.2e1 * t407;
    const double t413 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t401 * t402 + 0.4e1 / 0.3e1 * t61 * t409 );
    const double t414 = t66 * t66;
    const double t415 = 0.1e1 / t414;
    const double t416 = t192 * t192;
    const double t419 = -t409;
    const double t423 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t415 * t416 + 0.4e1 / 0.3e1 * t66 * t419 );
    const double t425 = ( t413 + t423 ) * t73;
    const double t426 = t425 * t138;
    const double t427 = t54 * t426;
    const double t428 = t197 * t253;
    const double t429 = t54 * t428;
    const double t430 = 0.2e1 * t429;
    const double t436 = t149 * t207;
    const double t440 = t206 * t98;
    const double t441 = 0.1e1 / t440;
    const double t442 = t80 * t441;
    const double t443 = t222 * t222;
    const double t444 = t443 * t223;
    const double t448 = t83 * t300 * t302;
    const double t456 = t90 * t319 * t302;
    const double t468 = -t448 * t309 / 0.18e2 + t210 * t314 / 0.9e1 + t87 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t456 * t309 + 0.83333333333333333333e-1 * t216 * t314 + t97 * t461 * t186 / 0.9e1 + t97 * t95 * t186 / 0.3e1;
    const double t469 = t468 * t223;
    const double t471 = t206 * t206;
    const double t472 = 0.1e1 / t471;
    const double t473 = t80 * t472;
    const double t474 = t102 * t102;
    const double t475 = 0.1e1 / t474;
    const double t477 = t443 * t475 * t82;
    const double t482 = t147 * t9 * t313 * t45;
    const double t483 = 0.2e1 / 0.9e1 * t482;
    const double t485 = t149 * t155;
    const double t487 = t484 * t485 * t175;
    const double t488 = t487 / 0.6e1;
    const double t489 = t8 * t313;
    const double t494 = t149 * t233;
    const double t498 = t232 * t128;
    const double t499 = 0.1e1 / t498;
    const double t500 = t110 * t499;
    const double t501 = t248 * t248;
    const double t502 = t501 * t249;
    const double t503 = t502 * t134;
    const double t506 = t371 * t249;
    const double t507 = t506 * t134;
    const double t509 = t232 * t232;
    const double t510 = 0.1e1 / t509;
    const double t511 = t110 * t510;
    const double t512 = t511 * t501;
    const double t513 = t132 * t132;
    const double t514 = 0.1e1 / t513;
    const double t515 = t514 * t134;
    const double t516 = t515 * t112;
    const double t519 = -0.2e1 / 0.9e1 * t201 * t9 * t313 * t103 - t435 * t436 * t224 / 0.6e1 - 0.2e1 * t442 * t444 + t208 * t469 + t473 * t477 / 0.2e1 + t483 + t488 + t381 - t335 - t390 - 0.2e1 / 0.9e1 * t227 * t489 * t135 - t493 * t494 * t251 / 0.6e1 - 0.2e1 * t500 * t503 + t234 * t507 + t512 * t516 / 0.2e1;
    const double t520 = t74 * t519;
    const double t521 = t54 * t520;
    const double t523 = t499 * t501 * t268;
    const double t524 = t266 * t523;
    const double t525 = 0.2e1 * t524;
    const double t526 = t74 * t511;
    const double t527 = t501 * t514;
    const double t529 = t527 * t528;
    const double t530 = t526 * t529;
    const double t531 = t530 / 0.2e1;
    const double t532 = t197 * t259;
    const double t533 = t532 * t263;
    const double t534 = t533 / 0.3e1;
    const double t536 = t425 * t105 * t143;
    const double t537 = 0.2e1 * t536;
    const double t540 = t9 * t313 * t133 * t134;
    const double t541 = t260 * t540;
    const double t542 = 0.2e1 / 0.9e1 * t541;
    const double t544 = t74 * t543;
    const double t545 = t228 * t233;
    const double t546 = t545 * t251;
    const double t547 = t544 * t546;
    const double t548 = t547 / 0.6e1;
    const double t549 = -t399 + t427 + t430 + t521 - t483 - t488 + t525 - t531 - t534 + t537 + t542 + t548;
    const double t552 = t282 * t110;
    const double t553 = t552 * t269;
    const double t555 = t282 * t259;
    const double t556 = t555 * t263;
    const double t558 = t548 - t553 - t483 + t525 - t348 - t374 - t488 - t531 + t335 - t533 / 0.6e1 + t542 - t556 / 0.6e1 - t338;
    const double t559 = t178 * t283;
    const double t561 = t183 * t283;
    const double t563 = t401 * t273;
    const double t566 = t61 * t48;
    const double t570 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t563 * t188 + 0.8e1 / 0.3e1 * t566 * t406 );
    const double t571 = t415 * t277;
    const double t574 = t66 * t48;
    const double t578 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t571 * t192 - 0.8e1 / 0.3e1 * t574 * t406 );
    const double t580 = ( t570 + t578 ) * t73;
    const double t581 = t580 * t138;
    const double t582 = t54 * t581;
    const double t586 = t580 * t105 * t143;
    const double t588 = t282 * t253;
    const double t589 = t54 * t588;
    const double t590 = t346 + 0.4e1 * t559 - 0.4e1 * t561 + t582 - t381 + t390 - 0.4e1 * t391 - 0.4e1 * t395 - t399 + t429 + t521 + 0.2e1 * t586 + t589;
    const double t595 = t273 * t273;
    const double t599 = 0.2e1 * t186 + 0.2e1 * t407;
    const double t603 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t401 * t595 + 0.4e1 / 0.3e1 * t61 * t599 );
    const double t604 = t277 * t277;
    const double t607 = -t599;
    const double t611 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t415 * t604 + 0.4e1 / 0.3e1 * t66 * t607 );
    const double t613 = ( t603 + t611 ) * t73;
    const double t614 = t613 * t138;
    const double t615 = t54 * t614;
    const double t617 = t613 * t105 * t143;
    const double t618 = 0.2e1 * t617;
    const double t619 = t615 + t618 + t335 + t338 + t341 + t346 - t374 - t381 + t390 - t394 - t399;
    const double t620 = t556 / 0.3e1;
    const double t621 = 0.8e1 * t559;
    const double t622 = 0.8e1 * t561;
    const double t623 = 0.2e1 * t589;
    const double t624 = 0.2e1 * t553;
    const double t625 = t521 - t483 - t488 + t525 - t531 - t620 - t621 - t622 + t623 + t542 + t548 - t624;


    v2rho2_aa = t290 + t291 + t292 - t293 + 0.2e1 * t199 + t295 + 0.4e1 * t257 - t297 - t298 + t10 * ( t397 + t549 );
    v2rho2_ab = t290 + t291 - t293 + t199 + t295 + t258 - t297 - t298 + t284 + t287 + t10 * ( t558 + t590 );
    v2rho2_bb = t290 + t291 - t292 - t293 + 0.2e1 * t284 + t295 + 0.4e1 * t286 - t297 - t298 + t10 * ( t619 + t625 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t70 = constants::m_cbrt_2;
    constexpr double t1 = a_0;
    constexpr double t2 = alpha1_0;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t18 = 0.1e1 / t1;
    constexpr double t19 = beta1_0;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = beta2_0 * t3;
    constexpr double t30 = beta3_0;
    constexpr double t37 = pp_0 + 0.1e1;
    constexpr double t75 = a_1;
    constexpr double t76 = alpha1_1;
    constexpr double t77 = t76 * t3;
    constexpr double t82 = 0.1e1 / t75;
    constexpr double t83 = beta1_1;
    constexpr double t87 = beta2_1 * t3;
    constexpr double t90 = beta3_1;
    constexpr double t95 = pp_1 + 0.1e1;
    constexpr double t105 = a_2;
    constexpr double t106 = alpha1_2;
    constexpr double t107 = t106 * t3;
    constexpr double t112 = 0.1e1 / t105;
    constexpr double t113 = beta1_2;
    constexpr double t117 = beta2_2 * t3;
    constexpr double t120 = beta3_2;
    constexpr double t125 = pp_2 + 0.1e1;
    constexpr double t134 = 0.1e1 / fz20;
    constexpr double t147 = t1 * t2 * t3;
    constexpr double t201 = t75 * t76 * t3;
    constexpr double t226 = t105 * t106;
    constexpr double t227 = t226 * t20;
    constexpr double t259 = t226 * t3;
    constexpr double t302 = t3 * t3;
    constexpr double t304 = t6 * t6;
    constexpr double t305 = t304 * t7;
    constexpr double t326 = t37 * t37;
    constexpr double t364 = t125 * t125;
    constexpr double t435 = t77 * t9;
    constexpr double t461 = t95 * t95;
    constexpr double t484 = t4 * t9;
    constexpr double t493 = t107 * t9;
    constexpr double t528 = t134 * t112;
    constexpr double t543 = t107 * t6;


    const double t10 = rho_a + rho_b;
    const double t11 = safe_math::cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = safe_math::sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = safe_math::pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.125e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = safe_math::log( t44 );
    const double t46 = t1 * t16 * t45;
    const double t47 = 0.2e1 * t46;
    const double t48 = rho_a - rho_b;
    const double t49 = t48 * t48;
    const double t50 = t49 * t49;
    const double t51 = t10 * t10;
    const double t52 = t51 * t51;
    const double t53 = 0.1e1 / t52;
    const double t54 = t50 * t53;
    const double t55 = 0.1e1 / t10;
    const double t56 = t48 * t55;
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( zeta_tol );
    const double t60 = t59 * zeta_tol;
    const double t61 = safe_math::cbrt( t57 );
    const double t63 = piecewise_functor_3( t58, t60, t61 * t57 );
    const double t64 = 0.1e1 - t56;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t60, t66 * t64 );
    const double t69 = t63 + t68 - 0.2e1;
    const double t73 = 0.1e1 / ( 0.2e1 * t70 - 0.2e1 );
    const double t74 = t69 * t73;
    const double t80 = 0.1e1 + t77 * t13 / 0.4e1;
    const double t96 = safe_math::pow( t35, t95 );
    const double t97 = beta4_1 * t96;
    const double t98 = t83 * t23 / 0.2e1 + t87 * t13 / 0.4e1 + 0.125e0 * t90 * t31 + t97;
    const double t102 = 0.1e1 + t82 / t98 / 0.2e1;
    const double t103 = safe_math::log( t102 );
    const double t110 = 0.1e1 + t107 * t13 / 0.4e1;
    const double t126 = safe_math::pow( t35, t125 );
    const double t127 = beta4_2 * t126;
    const double t128 = t113 * t23 / 0.2e1 + t117 * t13 / 0.4e1 + 0.125e0 * t120 * t31 + t127;
    const double t132 = 0.1e1 + t112 / t128 / 0.2e1;
    const double t133 = safe_math::log( t132 );
    const double t135 = t133 * t134;
    const double t138 = -0.2e1 * t75 * t80 * t103 - 0.2e1 * t105 * t110 * t135 + 0.2e1 * t46;
    const double t139 = t74 * t138;
    const double t140 = t54 * t139;
    const double t143 = t110 * t133 * t134;
    const double t145 = 0.2e1 * t74 * t105 * t143;
    const double t149 = 0.1e1 / t11 / t10;
    const double t152 = t147 * t9 * t149 * t45;
    const double t153 = t152 / 0.6e1;
    const double t154 = t40 * t40;
    const double t155 = 0.1e1 / t154;
    const double t156 = t16 * t155;
    const double t157 = 0.1e1 / t23;
    const double t159 = t19 * t157 * t3;
    const double t160 = t9 * t149;
    const double t165 = safe_math::sqrt( t22 );
    const double t167 = t30 * t165 * t3;
    const double t173 = -t159 * t160 / 0.12e2 - t27 * t160 / 0.12e2 - 0.625e-1 * t167 * t160 - t39 * t37 * t55 / 0.3e1;
    const double t174 = 0.1e1 / t44;
    const double t175 = t173 * t174;
    const double t176 = t156 * t175;
    const double t177 = t49 * t48;
    const double t178 = t177 * t53;
    const double t179 = t178 * t139;
    const double t180 = 0.4e1 * t179;
    const double t181 = t52 * t10;
    const double t182 = 0.1e1 / t181;
    const double t183 = t50 * t182;
    const double t184 = t183 * t139;
    const double t185 = 0.4e1 * t184;
    const double t186 = 0.1e1 / t51;
    const double t187 = t48 * t186;
    const double t188 = t55 - t187;
    const double t191 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t188 );
    const double t192 = -t188;
    const double t195 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t192 );
    const double t197 = ( t191 + t195 ) * t73;
    const double t198 = t197 * t138;
    const double t199 = t54 * t198;
    const double t206 = t98 * t98;
    const double t207 = 0.1e1 / t206;
    const double t208 = t80 * t207;
    const double t210 = t83 * t157 * t3;
    const double t216 = t90 * t165 * t3;
    const double t222 = -t210 * t160 / 0.12e2 - t87 * t160 / 0.12e2 - 0.625e-1 * t216 * t160 - t97 * t95 * t55 / 0.3e1;
    const double t223 = 0.1e1 / t102;
    const double t224 = t222 * t223;
    const double t228 = t8 * t149;
    const double t232 = t128 * t128;
    const double t233 = 0.1e1 / t232;
    const double t234 = t110 * t233;
    const double t236 = t113 * t157 * t3;
    const double t242 = t120 * t165 * t3;
    const double t248 = -t236 * t160 / 0.12e2 - t117 * t160 / 0.12e2 - 0.625e-1 * t242 * t160 - t127 * t125 * t55 / 0.3e1;
    const double t249 = 0.1e1 / t132;
    const double t251 = t248 * t249 * t134;
    const double t253 = t201 * t9 * t149 * t103 / 0.6e1 + t208 * t224 - t153 - t176 + t227 * t228 * t135 / 0.6e1 + t234 * t251;
    const double t254 = t74 * t253;
    const double t255 = t54 * t254;
    const double t257 = t197 * t105 * t143;
    const double t258 = 0.2e1 * t257;
    const double t260 = t74 * t259;
    const double t263 = t9 * t149 * t133 * t134;
    const double t264 = t260 * t263;
    const double t265 = t264 / 0.6e1;
    const double t266 = t74 * t110;
    const double t268 = t249 * t134;
    const double t269 = t233 * t248 * t268;
    const double t270 = t266 * t269;
    const double t273 = -t55 - t187;
    const double t276 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t61 * t273 );
    const double t277 = -t273;
    const double t280 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t277 );
    const double t282 = ( t276 + t280 ) * t73;
    const double t283 = t282 * t138;
    const double t284 = t54 * t283;
    const double t286 = t282 * t105 * t143;
    const double t287 = 0.2e1 * t286;
    const double t290 = t152 / 0.3e1;
    const double t291 = 0.2e1 * t176;
    const double t292 = 0.8e1 * t179;
    const double t293 = 0.8e1 * t184;
    const double t295 = 0.2e1 * t255;
    const double t297 = t264 / 0.3e1;
    const double t298 = 0.2e1 * t270;
    const double t300 = 0.1e1 / t23 / t22;
    const double t303 = t19 * t300 * t302;
    const double t306 = t11 * t11;
    const double t309 = t305 / t306 / t51;
    const double t313 = 0.1e1 / t11 / t51;
    const double t314 = t9 * t313;
    const double t319 = 0.1e1/safe_math::sqrt( t22 );
    const double t321 = t30 * t319 * t302;
    const double t333 = -t303 * t309 / 0.18e2 + t159 * t314 / 0.9e1 + t27 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t321 * t309 + 0.83333333333333333333e-1 * t167 * t314 + t39 * t326 * t186 / 0.9e1 + t39 * t37 * t186 / 0.3e1;
    const double t334 = t333 * t174;
    const double t335 = t156 * t334;
    const double t336 = t49 * t53;
    const double t337 = t336 * t139;
    const double t338 = 0.12e2 * t337;
    const double t339 = t177 * t182;
    const double t340 = t339 * t139;
    const double t341 = 0.32e2 * t340;
    const double t343 = 0.1e1 / t52 / t51;
    const double t344 = t50 * t343;
    const double t345 = t344 * t139;
    const double t346 = 0.2e2 * t345;
    const double t347 = t197 * t110;
    const double t348 = t347 * t269;
    const double t349 = 0.2e1 * t348;
    const double t351 = t113 * t300 * t302;
    const double t359 = t120 * t319 * t302;
    const double t371 = -t351 * t309 / 0.18e2 + t236 * t314 / 0.9e1 + t117 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t359 * t309 + 0.83333333333333333333e-1 * t242 * t314 + t127 * t364 * t186 / 0.9e1 + t127 * t125 * t186 / 0.3e1;
    const double t373 = t233 * t371 * t268;
    const double t374 = t266 * t373;
    const double t375 = t154 * t40;
    const double t376 = 0.1e1 / t375;
    const double t377 = t16 * t376;
    const double t378 = t173 * t173;
    const double t379 = t378 * t174;
    const double t380 = t377 * t379;
    const double t381 = 0.2e1 * t380;
    const double t382 = t154 * t154;
    const double t383 = 0.1e1 / t382;
    const double t384 = t16 * t383;
    const double t385 = t44 * t44;
    const double t386 = 0.1e1 / t385;
    const double t388 = t378 * t386 * t18;
    const double t389 = t384 * t388;
    const double t390 = t389 / 0.2e1;
    const double t391 = t178 * t198;
    const double t392 = 0.8e1 * t391;
    const double t393 = t178 * t254;
    const double t394 = 0.8e1 * t393;
    const double t395 = t183 * t198;
    const double t396 = 0.8e1 * t395;
    const double t397 = t335 + t338 - t341 + t346 - t349 - t374 - t381 + t390 + t392 + t394 - t396;
    const double t398 = t183 * t254;
    const double t399 = 0.8e1 * t398;
    const double t400 = t61 * t61;
    const double t401 = 0.1e1 / t400;
    const double t402 = t188 * t188;
    const double t405 = t51 * t10;
    const double t406 = 0.1e1 / t405;
    const double t407 = t48 * t406;
    const double t409 = -0.2e1 * t186 + 0.2e1 * t407;
    const double t413 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t401 * t402 + 0.4e1 / 0.3e1 * t61 * t409 );
    const double t414 = t66 * t66;
    const double t415 = 0.1e1 / t414;
    const double t416 = t192 * t192;
    const double t419 = -t409;
    const double t423 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t415 * t416 + 0.4e1 / 0.3e1 * t66 * t419 );
    const double t425 = ( t413 + t423 ) * t73;
    const double t426 = t425 * t138;
    const double t427 = t54 * t426;
    const double t428 = t197 * t253;
    const double t429 = t54 * t428;
    const double t430 = 0.2e1 * t429;
    const double t436 = t149 * t207;
    const double t440 = t206 * t98;
    const double t441 = 0.1e1 / t440;
    const double t442 = t80 * t441;
    const double t443 = t222 * t222;
    const double t444 = t443 * t223;
    const double t448 = t83 * t300 * t302;
    const double t456 = t90 * t319 * t302;
    const double t468 = -t448 * t309 / 0.18e2 + t210 * t314 / 0.9e1 + t87 * t314 / 0.9e1 + 0.41666666666666666666e-1 * t456 * t309 + 0.83333333333333333333e-1 * t216 * t314 + t97 * t461 * t186 / 0.9e1 + t97 * t95 * t186 / 0.3e1;
    const double t469 = t468 * t223;
    const double t471 = t206 * t206;
    const double t472 = 0.1e1 / t471;
    const double t473 = t80 * t472;
    const double t474 = t102 * t102;
    const double t475 = 0.1e1 / t474;
    const double t477 = t443 * t475 * t82;
    const double t482 = t147 * t9 * t313 * t45;
    const double t483 = 0.2e1 / 0.9e1 * t482;
    const double t485 = t149 * t155;
    const double t487 = t484 * t485 * t175;
    const double t488 = t487 / 0.6e1;
    const double t489 = t8 * t313;
    const double t494 = t149 * t233;
    const double t498 = t232 * t128;
    const double t499 = 0.1e1 / t498;
    const double t500 = t110 * t499;
    const double t501 = t248 * t248;
    const double t502 = t501 * t249;
    const double t503 = t502 * t134;
    const double t506 = t371 * t249;
    const double t507 = t506 * t134;
    const double t509 = t232 * t232;
    const double t510 = 0.1e1 / t509;
    const double t511 = t110 * t510;
    const double t512 = t511 * t501;
    const double t513 = t132 * t132;
    const double t514 = 0.1e1 / t513;
    const double t515 = t514 * t134;
    const double t516 = t515 * t112;
    const double t519 = -0.2e1 / 0.9e1 * t201 * t9 * t313 * t103 - t435 * t436 * t224 / 0.6e1 - 0.2e1 * t442 * t444 + t208 * t469 + t473 * t477 / 0.2e1 + t483 + t488 + t381 - t335 - t390 - 0.2e1 / 0.9e1 * t227 * t489 * t135 - t493 * t494 * t251 / 0.6e1 - 0.2e1 * t500 * t503 + t234 * t507 + t512 * t516 / 0.2e1;
    const double t520 = t74 * t519;
    const double t521 = t54 * t520;
    const double t523 = t499 * t501 * t268;
    const double t524 = t266 * t523;
    const double t525 = 0.2e1 * t524;
    const double t526 = t74 * t511;
    const double t527 = t501 * t514;
    const double t529 = t527 * t528;
    const double t530 = t526 * t529;
    const double t531 = t530 / 0.2e1;
    const double t532 = t197 * t259;
    const double t533 = t532 * t263;
    const double t534 = t533 / 0.3e1;
    const double t536 = t425 * t105 * t143;
    const double t537 = 0.2e1 * t536;
    const double t540 = t9 * t313 * t133 * t134;
    const double t541 = t260 * t540;
    const double t542 = 0.2e1 / 0.9e1 * t541;
    const double t544 = t74 * t543;
    const double t545 = t228 * t233;
    const double t546 = t545 * t251;
    const double t547 = t544 * t546;
    const double t548 = t547 / 0.6e1;
    const double t549 = -t399 + t427 + t430 + t521 - t483 - t488 + t525 - t531 - t534 + t537 + t542 + t548;
    const double t552 = t282 * t110;
    const double t553 = t552 * t269;
    const double t555 = t282 * t259;
    const double t556 = t555 * t263;
    const double t558 = t548 - t553 - t483 + t525 - t348 - t374 - t488 - t531 + t335 - t533 / 0.6e1 + t542 - t556 / 0.6e1 - t338;
    const double t559 = t178 * t283;
    const double t561 = t183 * t283;
    const double t563 = t401 * t273;
    const double t566 = t61 * t48;
    const double t570 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t563 * t188 + 0.8e1 / 0.3e1 * t566 * t406 );
    const double t571 = t415 * t277;
    const double t574 = t66 * t48;
    const double t578 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t571 * t192 - 0.8e1 / 0.3e1 * t574 * t406 );
    const double t580 = ( t570 + t578 ) * t73;
    const double t581 = t580 * t138;
    const double t582 = t54 * t581;
    const double t586 = t580 * t105 * t143;
    const double t588 = t282 * t253;
    const double t589 = t54 * t588;
    const double t590 = t346 + 0.4e1 * t559 - 0.4e1 * t561 + t582 - t381 + t390 - 0.4e1 * t391 - 0.4e1 * t395 - t399 + t429 + t521 + 0.2e1 * t586 + t589;
    const double t595 = t273 * t273;
    const double t599 = 0.2e1 * t186 + 0.2e1 * t407;
    const double t603 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.9e1 * t401 * t595 + 0.4e1 / 0.3e1 * t61 * t599 );
    const double t604 = t277 * t277;
    const double t607 = -t599;
    const double t611 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t415 * t604 + 0.4e1 / 0.3e1 * t66 * t607 );
    const double t613 = ( t603 + t611 ) * t73;
    const double t614 = t613 * t138;
    const double t615 = t54 * t614;
    const double t617 = t613 * t105 * t143;
    const double t618 = 0.2e1 * t617;
    const double t619 = t615 + t618 + t335 + t338 + t341 + t346 - t374 - t381 + t390 - t394 - t399;
    const double t620 = t556 / 0.3e1;
    const double t621 = 0.8e1 * t559;
    const double t622 = 0.8e1 * t561;
    const double t623 = 0.2e1 * t589;
    const double t624 = 0.2e1 * t553;
    const double t625 = t521 - t483 - t488 + t525 - t531 - t620 - t621 - t622 + t623 + t542 + t548 - t624;


    vrho_a = -t47 + t140 + t145 + t10 * ( t153 + t176 + t180 - t185 + t199 + t255 + t258 - t265 - t270 );
    vrho_b = -t47 + t140 + t145 + t10 * ( t153 + t176 - t180 - t185 + t284 + t255 + t287 - t265 - t270 );
    v2rho2_aa = t290 + t291 + t292 - t293 + 0.2e1 * t199 + t295 + 0.4e1 * t257 - t297 - t298 + t10 * ( t397 + t549 );
    v2rho2_ab = t290 + t291 - t293 + t199 + t295 + t258 - t297 - t298 + t284 + t287 + t10 * ( t558 + t590 );
    v2rho2_bb = t290 + t291 - t292 - t293 + 0.2e1 * t284 + t295 + 0.4e1 * t286 - t297 - t298 + t10 * ( t619 + t625 );

  }


};

struct BuiltinPW91_LDA_RPA : detail::BuiltinKernelImpl< BuiltinPW91_LDA_RPA > {

  BuiltinPW91_LDA_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPW91_LDA_RPA >(p) { }
  
  virtual ~BuiltinPW91_LDA_RPA() = default;

};



} // namespace ExchCXX
