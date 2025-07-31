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
struct kernel_traits< BuiltinMPW91_X > :
  public gga_screening_interface< BuiltinMPW91_X > {

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


  static constexpr double a = 6.0*0.00426/constants::X2S;
  static constexpr double b = 1.0/constants::X2S;
  static constexpr double c = 0.00426/(constants::X_FACTOR_C*constants::X2S*constants::X2S);
  static constexpr double d = -(0.00426 - 5.0*0.0003780762333399851)/(constants::X_FACTOR_C*constants::X2S*constants::X2S);
  static constexpr double f = 1.0e-6/(constants::X_FACTOR_C*0.00048120394750740677);
  static constexpr double alpha = 100.0;
  static constexpr double expo = 3.72;

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
    constexpr double t44 = t20 * t20;
    constexpr double t45 = 0.1e1 / t23;
    constexpr double t46 = t44 * t45;


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
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t37 = safe_math::exp( -alpha * t20 * t25 * t34 / 0.24e2 );
    const double t40 = ( d * t37 + c ) * t20;
    const double t41 = t40 * t25;
    const double t47 = safe_math::sqrt( sigma );
    const double t50 = 0.1e1 / t18 / rho;
    const double t51 = t47 * t27 * t50;
    const double t54 = safe_math::pow( t46 * t51 / 0.12e2, expo );
    const double t55 = f * t54;
    const double t56 = t41 * t34 / 0.24e2 - t55;
    const double t57 = t46 * t47;
    const double t63 = safe_math::log( b * t44 * t45 * t51 / 0.12e2 + safe_math::sqrt( square( b * t44 * t45 * t51 / 0.12e2 ) + 0.1e1 ) );
    const double t64 = a * t63;
    const double t65 = t27 * t50 * t64;
    const double t68 = 0.1e1 + t57 * t65 / 0.12e2 + t55;
    const double t69 = 0.1e1 / t68;
    const double t71 = t56 * t69 + 0.1e1;
    const double t75 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t71 );


    eps = 0.2e1 * t75;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

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
    constexpr double t44 = t20 * t20;
    constexpr double t45 = 0.1e1 / t23;
    constexpr double t46 = t44 * t45;
    constexpr double t81 = d * alpha;
    constexpr double t83 = 0.1e1 / t23 / t22;
    constexpr double t84 = t44 * t83;
    constexpr double t85 = t81 * t84;
    constexpr double t117 = t20 * t25;
    constexpr double t120 = b * b;
    constexpr double t150 = t25 * t28;
    constexpr double t164 = t117 * t28;


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
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t37 = safe_math::exp( -alpha * t20 * t25 * t34 / 0.24e2 );
    const double t40 = ( d * t37 + c ) * t20;
    const double t41 = t40 * t25;
    const double t47 = safe_math::sqrt( sigma );
    const double t50 = 0.1e1 / t18 / rho;
    const double t51 = t47 * t27 * t50;
    const double t54 = safe_math::pow( t46 * t51 / 0.12e2, expo );
    const double t55 = f * t54;
    const double t56 = t41 * t34 / 0.24e2 - t55;
    const double t57 = t46 * t47;
    const double t63 = safe_math::log( b * t44 * t45 * t51 / 0.12e2 + safe_math::sqrt( square( b * t44 * t45 * t51 / 0.12e2 ) + 0.1e1 ) );
    const double t64 = a * t63;
    const double t65 = t27 * t50 * t64;
    const double t68 = 0.1e1 + t57 * t65 / 0.12e2 + t55;
    const double t69 = 0.1e1 / t68;
    const double t71 = t56 * t69 + 0.1e1;
    const double t75 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t71 );
    const double t77 = t17 / t31;
    const double t86 = sigma * sigma;
    const double t87 = t86 * t27;
    const double t88 = t30 * t30;
    const double t89 = t88 * t30;
    const double t91 = 0.1e1 / t18 / t89;
    const double t92 = t91 * t37;
    const double t96 = t30 * rho;
    const double t98 = 0.1e1 / t31 / t96;
    const double t102 = 0.1e1 / rho;
    const double t105 = 0.4e1 / 0.3e1 * t55 * expo * t102;
    const double t106 = t85 * t87 * t92 / 0.108e3 - t41 * t29 * t98 / 0.9e1 + t105;
    const double t108 = t68 * t68;
    const double t109 = 0.1e1 / t108;
    const double t110 = t56 * t109;
    const double t114 = t27 / t18 / t30 * t64;
    const double t118 = t117 * t29;
    const double t125 = 0.6e1 * t120 * t20 * t25 * t34 + 0.144e3;
    const double t126 = safe_math::sqrt( t125 );
    const double t128 = b / t126;
    const double t129 = t98 * a * t128;
    const double t132 = -t57 * t114 / 0.9e1 - 0.2e1 / 0.3e1 * t118 * t129 - t105;
    const double t134 = t106 * t69 - t110 * t132;
    const double t139 = piecewise_functor_3( t2, 0.0, -t6 * t77 * t71 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t134 );
    const double t142 = t88 * rho;
    const double t144 = 0.1e1 / t18 / t142;
    const double t145 = t27 * t144;
    const double t146 = t37 * sigma;
    const double t154 = 0.1e1 / sigma;
    const double t157 = t55 * expo * t154 / 0.2e1;
    const double t158 = -t85 * t145 * t146 / 0.288e3 + t40 * t150 * t33 / 0.24e2 - t157;
    const double t161 = t46 / t47;
    const double t166 = t33 * a * t128;
    const double t169 = t161 * t65 / 0.24e2 + t164 * t166 / 0.4e1 + t157;
    const double t171 = -t110 * t169 + t158 * t69;
    const double t175 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t171 );


    eps = 0.2e1 * t75;
    vrho = 0.2e1 * rho * t139 + 0.2e1 * t75;
    vsigma = 0.2e1 * rho * t175;

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
    constexpr double t44 = t20 * t20;
    constexpr double t45 = 0.1e1 / t23;
    constexpr double t46 = t44 * t45;
    constexpr double t81 = d * alpha;
    constexpr double t83 = 0.1e1 / t23 / t22;
    constexpr double t84 = t44 * t83;
    constexpr double t85 = t81 * t84;
    constexpr double t117 = t20 * t25;
    constexpr double t120 = b * b;
    constexpr double t150 = t25 * t28;
    constexpr double t164 = t117 * t28;
    constexpr double t194 = alpha * alpha;
    constexpr double t195 = d * t194;
    constexpr double t196 = t22 * t22;
    constexpr double t197 = 0.1e1 / t196;
    constexpr double t198 = t195 * t197;
    constexpr double t212 = expo * expo;
    constexpr double t243 = t120 * b;
    constexpr double t312 = t81 * t44;
    constexpr double t313 = t83 * t27;
    constexpr double t341 = t84 * t27;


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
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t37 = safe_math::exp( -alpha * t20 * t25 * t34 / 0.24e2 );
    const double t40 = ( d * t37 + c ) * t20;
    const double t41 = t40 * t25;
    const double t47 = safe_math::sqrt( sigma );
    const double t50 = 0.1e1 / t18 / rho;
    const double t51 = t47 * t27 * t50;
    const double t54 = safe_math::pow( t46 * t51 / 0.12e2, expo );
    const double t55 = f * t54;
    const double t56 = t41 * t34 / 0.24e2 - t55;
    const double t57 = t46 * t47;
    const double t63 = safe_math::log( b * t44 * t45 * t51 / 0.12e2 + safe_math::sqrt( square( b * t44 * t45 * t51 / 0.12e2 ) + 0.1e1 ) );
    const double t64 = a * t63;
    const double t65 = t27 * t50 * t64;
    const double t68 = 0.1e1 + t57 * t65 / 0.12e2 + t55;
    const double t69 = 0.1e1 / t68;
    const double t71 = t56 * t69 + 0.1e1;
    const double t77 = t17 / t31;
    const double t86 = sigma * sigma;
    const double t87 = t86 * t27;
    const double t88 = t30 * t30;
    const double t89 = t88 * t30;
    const double t91 = 0.1e1 / t18 / t89;
    const double t92 = t91 * t37;
    const double t96 = t30 * rho;
    const double t98 = 0.1e1 / t31 / t96;
    const double t102 = 0.1e1 / rho;
    const double t105 = 0.4e1 / 0.3e1 * t55 * expo * t102;
    const double t106 = t85 * t87 * t92 / 0.108e3 - t41 * t29 * t98 / 0.9e1 + t105;
    const double t108 = t68 * t68;
    const double t109 = 0.1e1 / t108;
    const double t110 = t56 * t109;
    const double t114 = t27 / t18 / t30 * t64;
    const double t118 = t117 * t29;
    const double t125 = 0.6e1 * t120 * t20 * t25 * t34 + 0.144e3;
    const double t126 = safe_math::sqrt( t125 );
    const double t128 = b / t126;
    const double t129 = t98 * a * t128;
    const double t132 = -t57 * t114 / 0.9e1 - 0.2e1 / 0.3e1 * t118 * t129 - t105;
    const double t134 = t106 * t69 - t110 * t132;
    const double t139 = piecewise_functor_3( t2, 0.0, -t6 * t77 * t71 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t134 );
    const double t142 = t88 * rho;
    const double t144 = 0.1e1 / t18 / t142;
    const double t145 = t27 * t144;
    const double t146 = t37 * sigma;
    const double t154 = 0.1e1 / sigma;
    const double t157 = t55 * expo * t154 / 0.2e1;
    const double t158 = -t85 * t145 * t146 / 0.288e3 + t40 * t150 * t33 / 0.24e2 - t157;
    const double t161 = t46 / t47;
    const double t166 = t33 * a * t128;
    const double t169 = t161 * t65 / 0.24e2 + t164 * t166 / 0.4e1 + t157;
    const double t171 = -t110 * t169 + t158 * t69;
    const double t175 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t171 );
    const double t180 = t17 / t31 / rho;
    const double t187 = t88 * t96;
    const double t189 = 0.1e1 / t18 / t187;
    const double t190 = t189 * t37;
    const double t199 = t86 * sigma;
    const double t200 = t88 * t88;
    const double t201 = t200 * t30;
    const double t202 = 0.1e1 / t201;
    const double t208 = 0.1e1 / t31 / t88;
    const double t213 = 0.1e1 / t30;
    const double t214 = t212 * t213;
    const double t216 = 0.16e2 / 0.9e1 * t55 * t214;
    const double t219 = 0.4e1 / 0.3e1 * t55 * expo * t213;
    const double t220 = -t85 * t87 * t190 / 0.12e2 + t198 * t199 * t202 * t37 / 0.81e2 + 0.11e2 / 0.27e2 * t41 * t29 * t208 - t216 - t219;
    const double t222 = t106 * t109;
    const double t226 = 0.1e1 / t108 / t68;
    const double t227 = t56 * t226;
    const double t228 = t132 * t132;
    const double t234 = t27 / t18 / t96 * t64;
    const double t238 = t208 * a * t128;
    const double t241 = t84 * t87;
    const double t245 = 0.1e1 / t126 / t125;
    const double t246 = t243 * t245;
    const double t247 = t189 * a * t246;
    const double t250 = 0.7e1 / 0.27e2 * t57 * t234 + 0.1e2 / 0.3e1 * t118 * t238 - 0.32e2 / 0.3e1 * t241 * t247 + t216 + t219;
    const double t252 = -t110 * t250 - 0.2e1 * t222 * t132 + t220 * t69 + 0.2e1 * t227 * t228;
    const double t257 = piecewise_functor_3( t2, 0.0, t6 * t180 * t71 / 0.12e2 - t6 * t77 * t134 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t252 );
    const double t263 = t27 * t91;
    const double t267 = t200 * rho;
    const double t268 = 0.1e1 / t267;
    const double t276 = t212 * t102;
    const double t279 = 0.2e1 / 0.3e1 * t55 * t276 * t154;
    const double t280 = t85 * t263 * t146 / 0.36e2 - t198 * t268 * t86 * t37 / 0.216e3 - t40 * t150 * t98 / 0.9e1 + t279;
    const double t282 = t158 * t109;
    const double t285 = t169 * t132;
    const double t294 = a * t243 * t245 * sigma;
    const double t297 = -t161 * t114 / 0.18e2 - t164 * t129 + 0.4e1 * t84 * t263 * t294 - t279;
    const double t299 = -t110 * t297 - t282 * t132 - t222 * t169 + 0.2e1 * t227 * t285 + t280 * t69;
    const double t304 = piecewise_functor_3( t2, 0.0, -t6 * t77 * t171 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t299 );
    const double t307 = 0.1e1 / t200;
    const double t318 = 0.1e1 / t86;
    const double t321 = t55 * t212 * t318 / 0.4e1;
    const double t324 = t55 * expo * t318 / 0.2e1;
    const double t325 = t198 * t307 * t37 * sigma / 0.576e3 - t312 * t313 * t144 * t37 / 0.144e3 - t321 + t324;
    const double t329 = t169 * t169;
    const double t334 = t46 / t47 / sigma;
    const double t338 = t117 * t154 * t28;
    const double t343 = t144 * a * t246;
    const double t346 = -t334 * t65 / 0.48e2 + t338 * t166 / 0.8e1 - 0.3e1 / 0.2e1 * t341 * t343 + t321 - t324;
    const double t348 = -t110 * t346 - 0.2e1 * t282 * t169 + 0.2e1 * t227 * t329 + t325 * t69;
    const double t352 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t348 );


    v2rho2 = 0.2e1 * rho * t257 + 0.4e1 * t139;
    v2rhosigma = 0.2e1 * rho * t304 + 0.2e1 * t175;
    v2sigma2 = 0.2e1 * rho * t352;

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
    constexpr double t44 = t20 * t20;
    constexpr double t45 = 0.1e1 / t23;
    constexpr double t46 = t44 * t45;
    constexpr double t81 = d * alpha;
    constexpr double t83 = 0.1e1 / t23 / t22;
    constexpr double t84 = t44 * t83;
    constexpr double t85 = t81 * t84;
    constexpr double t117 = t20 * t25;
    constexpr double t120 = b * b;
    constexpr double t150 = t25 * t28;
    constexpr double t164 = t117 * t28;
    constexpr double t194 = alpha * alpha;
    constexpr double t195 = d * t194;
    constexpr double t196 = t22 * t22;
    constexpr double t197 = 0.1e1 / t196;
    constexpr double t198 = t195 * t197;
    constexpr double t212 = expo * expo;
    constexpr double t243 = t120 * b;
    constexpr double t312 = t81 * t44;
    constexpr double t313 = t83 * t27;
    constexpr double t341 = t84 * t27;


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
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t37 = safe_math::exp( -alpha * t20 * t25 * t34 / 0.24e2 );
    const double t40 = ( d * t37 + c ) * t20;
    const double t41 = t40 * t25;
    const double t47 = safe_math::sqrt( sigma );
    const double t50 = 0.1e1 / t18 / rho;
    const double t51 = t47 * t27 * t50;
    const double t54 = safe_math::pow( t46 * t51 / 0.12e2, expo );
    const double t55 = f * t54;
    const double t56 = t41 * t34 / 0.24e2 - t55;
    const double t57 = t46 * t47;
    const double t63 = safe_math::log( b * t44 * t45 * t51 / 0.12e2 + safe_math::sqrt( square( b * t44 * t45 * t51 / 0.12e2 ) + 0.1e1 ) );
    const double t64 = a * t63;
    const double t65 = t27 * t50 * t64;
    const double t68 = 0.1e1 + t57 * t65 / 0.12e2 + t55;
    const double t69 = 0.1e1 / t68;
    const double t71 = t56 * t69 + 0.1e1;
    const double t75 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t71 );
    const double t77 = t17 / t31;
    const double t86 = sigma * sigma;
    const double t87 = t86 * t27;
    const double t88 = t30 * t30;
    const double t89 = t88 * t30;
    const double t91 = 0.1e1 / t18 / t89;
    const double t92 = t91 * t37;
    const double t96 = t30 * rho;
    const double t98 = 0.1e1 / t31 / t96;
    const double t102 = 0.1e1 / rho;
    const double t105 = 0.4e1 / 0.3e1 * t55 * expo * t102;
    const double t106 = t85 * t87 * t92 / 0.108e3 - t41 * t29 * t98 / 0.9e1 + t105;
    const double t108 = t68 * t68;
    const double t109 = 0.1e1 / t108;
    const double t110 = t56 * t109;
    const double t114 = t27 / t18 / t30 * t64;
    const double t118 = t117 * t29;
    const double t125 = 0.6e1 * t120 * t20 * t25 * t34 + 0.144e3;
    const double t126 = safe_math::sqrt( t125 );
    const double t128 = b / t126;
    const double t129 = t98 * a * t128;
    const double t132 = -t57 * t114 / 0.9e1 - 0.2e1 / 0.3e1 * t118 * t129 - t105;
    const double t134 = t106 * t69 - t110 * t132;
    const double t139 = piecewise_functor_3( t2, 0.0, -t6 * t77 * t71 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t134 );
    const double t142 = t88 * rho;
    const double t144 = 0.1e1 / t18 / t142;
    const double t145 = t27 * t144;
    const double t146 = t37 * sigma;
    const double t154 = 0.1e1 / sigma;
    const double t157 = t55 * expo * t154 / 0.2e1;
    const double t158 = -t85 * t145 * t146 / 0.288e3 + t40 * t150 * t33 / 0.24e2 - t157;
    const double t161 = t46 / t47;
    const double t166 = t33 * a * t128;
    const double t169 = t161 * t65 / 0.24e2 + t164 * t166 / 0.4e1 + t157;
    const double t171 = -t110 * t169 + t158 * t69;
    const double t175 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t171 );
    const double t180 = t17 / t31 / rho;
    const double t187 = t88 * t96;
    const double t189 = 0.1e1 / t18 / t187;
    const double t190 = t189 * t37;
    const double t199 = t86 * sigma;
    const double t200 = t88 * t88;
    const double t201 = t200 * t30;
    const double t202 = 0.1e1 / t201;
    const double t208 = 0.1e1 / t31 / t88;
    const double t213 = 0.1e1 / t30;
    const double t214 = t212 * t213;
    const double t216 = 0.16e2 / 0.9e1 * t55 * t214;
    const double t219 = 0.4e1 / 0.3e1 * t55 * expo * t213;
    const double t220 = -t85 * t87 * t190 / 0.12e2 + t198 * t199 * t202 * t37 / 0.81e2 + 0.11e2 / 0.27e2 * t41 * t29 * t208 - t216 - t219;
    const double t222 = t106 * t109;
    const double t226 = 0.1e1 / t108 / t68;
    const double t227 = t56 * t226;
    const double t228 = t132 * t132;
    const double t234 = t27 / t18 / t96 * t64;
    const double t238 = t208 * a * t128;
    const double t241 = t84 * t87;
    const double t245 = 0.1e1 / t126 / t125;
    const double t246 = t243 * t245;
    const double t247 = t189 * a * t246;
    const double t250 = 0.7e1 / 0.27e2 * t57 * t234 + 0.1e2 / 0.3e1 * t118 * t238 - 0.32e2 / 0.3e1 * t241 * t247 + t216 + t219;
    const double t252 = -t110 * t250 - 0.2e1 * t222 * t132 + t220 * t69 + 0.2e1 * t227 * t228;
    const double t257 = piecewise_functor_3( t2, 0.0, t6 * t180 * t71 / 0.12e2 - t6 * t77 * t134 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t252 );
    const double t263 = t27 * t91;
    const double t267 = t200 * rho;
    const double t268 = 0.1e1 / t267;
    const double t276 = t212 * t102;
    const double t279 = 0.2e1 / 0.3e1 * t55 * t276 * t154;
    const double t280 = t85 * t263 * t146 / 0.36e2 - t198 * t268 * t86 * t37 / 0.216e3 - t40 * t150 * t98 / 0.9e1 + t279;
    const double t282 = t158 * t109;
    const double t285 = t169 * t132;
    const double t294 = a * t243 * t245 * sigma;
    const double t297 = -t161 * t114 / 0.18e2 - t164 * t129 + 0.4e1 * t84 * t263 * t294 - t279;
    const double t299 = -t110 * t297 - t282 * t132 - t222 * t169 + 0.2e1 * t227 * t285 + t280 * t69;
    const double t304 = piecewise_functor_3( t2, 0.0, -t6 * t77 * t171 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t299 );
    const double t307 = 0.1e1 / t200;
    const double t318 = 0.1e1 / t86;
    const double t321 = t55 * t212 * t318 / 0.4e1;
    const double t324 = t55 * expo * t318 / 0.2e1;
    const double t325 = t198 * t307 * t37 * sigma / 0.576e3 - t312 * t313 * t144 * t37 / 0.144e3 - t321 + t324;
    const double t329 = t169 * t169;
    const double t334 = t46 / t47 / sigma;
    const double t338 = t117 * t154 * t28;
    const double t343 = t144 * a * t246;
    const double t346 = -t334 * t65 / 0.48e2 + t338 * t166 / 0.8e1 - 0.3e1 / 0.2e1 * t341 * t343 + t321 - t324;
    const double t348 = -t110 * t346 - 0.2e1 * t282 * t169 + 0.2e1 * t227 * t329 + t325 * t69;
    const double t352 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t348 );


    vrho = 0.2e1 * rho * t139 + 0.2e1 * t75;
    vsigma = 0.2e1 * rho * t175;
    v2rho2 = 0.2e1 * rho * t257 + 0.4e1 * t139;
    v2rhosigma = 0.2e1 * rho * t304 + 0.2e1 * t175;
    v2sigma2 = 0.2e1 * rho * t352;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = alpha * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t49 = t28 * t28;
    constexpr double t50 = 0.1e1 / t31;
    constexpr double t51 = t49 * t50;
    constexpr double t63 = b * t49;


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
    const double t40 = t34 * t39;
    const double t43 = safe_math::exp( -t29 * t40 / 0.24e2 );
    const double t46 = ( d * t43 + c ) * t28;
    const double t52 = safe_math::sqrt( sigma_aa );
    const double t54 = 0.1e1 / t36 / rho_a;
    const double t58 = safe_math::pow( t51 * t52 * t54 / 0.12e2, expo );
    const double t59 = f * t58;
    const double t60 = t46 * t40 / 0.24e2 - t59;
    const double t61 = t51 * t52;
    const double t68 = safe_math::log( t63 * t50 * t52 * t54 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t52 * t54 / 0.12e2 ) + 0.1e1 ) );
    const double t69 = t54 * a * t68;
    const double t72 = 0.1e1 + t61 * t69 / 0.12e2 + t59;
    const double t73 = 0.1e1 / t72;
    const double t75 = t60 * t73 + 0.1e1;
    const double t79 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t75 );
    const double t80 = rho_b <= dens_tol;
    const double t81 = -t16;
    const double t83 = piecewise_functor_5( t14, t11, t10, t15, t81 * t7 );
    const double t84 = 0.1e1 + t83;
    const double t85 = t84 <= zeta_tol;
    const double t86 = safe_math::cbrt( t84 );
    const double t88 = piecewise_functor_3( t85, t22, t86 * t84 );
    const double t89 = t88 * t26;
    const double t90 = t33 * sigma_bb;
    const double t91 = rho_b * rho_b;
    const double t92 = safe_math::cbrt( rho_b );
    const double t93 = t92 * t92;
    const double t95 = 0.1e1 / t93 / t91;
    const double t96 = t90 * t95;
    const double t99 = safe_math::exp( -t29 * t96 / 0.24e2 );
    const double t102 = ( d * t99 + c ) * t28;
    const double t105 = safe_math::sqrt( sigma_bb );
    const double t107 = 0.1e1 / t92 / rho_b;
    const double t111 = safe_math::pow( t51 * t105 * t107 / 0.12e2, expo );
    const double t112 = f * t111;
    const double t113 = t102 * t96 / 0.24e2 - t112;
    const double t114 = t51 * t105;
    const double t120 = safe_math::log( t63 * t50 * t105 * t107 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t105 * t107 / 0.12e2 ) + 0.1e1 ) );
    const double t121 = t107 * a * t120;
    const double t124 = 0.1e1 + t114 * t121 / 0.12e2 + t112;
    const double t125 = 0.1e1 / t124;
    const double t127 = t113 * t125 + 0.1e1;
    const double t131 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t127 );


    eps = t79 + t131;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = alpha * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t49 = t28 * t28;
    constexpr double t50 = 0.1e1 / t31;
    constexpr double t51 = t49 * t50;
    constexpr double t63 = b * t49;
    constexpr double t151 = d * alpha * t49;
    constexpr double t153 = 0.1e1 / t31 / t30;
    constexpr double t185 = t28 * t33;
    constexpr double t188 = b * b;
    constexpr double t189 = t188 * t28;
    constexpr double t319 = a * b;


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
    const double t40 = t34 * t39;
    const double t43 = safe_math::exp( -t29 * t40 / 0.24e2 );
    const double t46 = ( d * t43 + c ) * t28;
    const double t52 = safe_math::sqrt( sigma_aa );
    const double t54 = 0.1e1 / t36 / rho_a;
    const double t58 = safe_math::pow( t51 * t52 * t54 / 0.12e2, expo );
    const double t59 = f * t58;
    const double t60 = t46 * t40 / 0.24e2 - t59;
    const double t61 = t51 * t52;
    const double t68 = safe_math::log( t63 * t50 * t52 * t54 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t52 * t54 / 0.12e2 ) + 0.1e1 ) );
    const double t69 = t54 * a * t68;
    const double t72 = 0.1e1 + t61 * t69 / 0.12e2 + t59;
    const double t73 = 0.1e1 / t72;
    const double t75 = t60 * t73 + 0.1e1;
    const double t79 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t75 );
    const double t80 = rho_b <= dens_tol;
    const double t81 = -t16;
    const double t83 = piecewise_functor_5( t14, t11, t10, t15, t81 * t7 );
    const double t84 = 0.1e1 + t83;
    const double t85 = t84 <= zeta_tol;
    const double t86 = safe_math::cbrt( t84 );
    const double t88 = piecewise_functor_3( t85, t22, t86 * t84 );
    const double t89 = t88 * t26;
    const double t90 = t33 * sigma_bb;
    const double t91 = rho_b * rho_b;
    const double t92 = safe_math::cbrt( rho_b );
    const double t93 = t92 * t92;
    const double t95 = 0.1e1 / t93 / t91;
    const double t96 = t90 * t95;
    const double t99 = safe_math::exp( -t29 * t96 / 0.24e2 );
    const double t102 = ( d * t99 + c ) * t28;
    const double t105 = safe_math::sqrt( sigma_bb );
    const double t107 = 0.1e1 / t92 / rho_b;
    const double t111 = safe_math::pow( t51 * t105 * t107 / 0.12e2, expo );
    const double t112 = f * t111;
    const double t113 = t102 * t96 / 0.24e2 - t112;
    const double t114 = t51 * t105;
    const double t120 = safe_math::log( t63 * t50 * t105 * t107 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t105 * t107 / 0.12e2 ) + 0.1e1 ) );
    const double t121 = t107 * a * t120;
    const double t124 = 0.1e1 + t114 * t121 / 0.12e2 + t112;
    const double t125 = 0.1e1 / t124;
    const double t127 = t113 * t125 + 0.1e1;
    const double t131 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t127 );
    const double t132 = t6 * t6;
    const double t133 = 0.1e1 / t132;
    const double t134 = t16 * t133;
    const double t136 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t134 );
    const double t139 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t136 );
    const double t140 = t139 * t26;
    const double t144 = t26 * t26;
    const double t145 = 0.1e1 / t144;
    const double t146 = t25 * t145;
    const double t149 = t5 * t146 * t75 / 0.8e1;
    const double t154 = sigma_aa * sigma_aa;
    const double t155 = t153 * t154;
    const double t156 = t35 * t35;
    const double t157 = t156 * t35;
    const double t159 = 0.1e1 / t36 / t157;
    const double t164 = t35 * rho_a;
    const double t166 = 0.1e1 / t37 / t164;
    const double t170 = 0.1e1 / rho_a;
    const double t173 = 0.4e1 / 0.3e1 * t59 * expo * t170;
    const double t174 = t151 * t155 * t159 * t43 / 0.216e3 - t46 * t34 * t166 / 0.9e1 + t173;
    const double t176 = t72 * t72;
    const double t177 = 0.1e1 / t176;
    const double t178 = t60 * t177;
    const double t182 = 0.1e1 / t36 / t35 * a * t68;
    const double t186 = t185 * sigma_aa;
    const double t192 = 0.6e1 * t189 * t40 + 0.144e3;
    const double t193 = safe_math::sqrt( t192 );
    const double t194 = 0.1e1 / t193;
    const double t195 = b * t194;
    const double t196 = t166 * a * t195;
    const double t199 = -t61 * t182 / 0.9e1 - 0.2e1 / 0.3e1 * t186 * t196 - t173;
    const double t201 = t174 * t73 - t178 * t199;
    const double t206 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t140 * t75 - t149 - 0.3e1 / 0.8e1 * t5 * t27 * t201 );
    const double t207 = t81 * t133;
    const double t209 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t207 );
    const double t212 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t209 );
    const double t213 = t212 * t26;
    const double t217 = t88 * t145;
    const double t220 = t5 * t217 * t127 / 0.8e1;
    const double t222 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t213 * t127 - t220 );
    const double t226 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t134 );
    const double t229 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t226 );
    const double t230 = t229 * t26;
    const double t235 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t230 * t75 - t149 );
    const double t237 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t207 );
    const double t240 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t237 );
    const double t241 = t240 * t26;
    const double t245 = sigma_bb * sigma_bb;
    const double t246 = t153 * t245;
    const double t247 = t91 * t91;
    const double t248 = t247 * t91;
    const double t250 = 0.1e1 / t92 / t248;
    const double t255 = t91 * rho_b;
    const double t257 = 0.1e1 / t93 / t255;
    const double t261 = 0.1e1 / rho_b;
    const double t264 = 0.4e1 / 0.3e1 * t112 * expo * t261;
    const double t265 = t151 * t246 * t250 * t99 / 0.216e3 - t102 * t90 * t257 / 0.9e1 + t264;
    const double t267 = t124 * t124;
    const double t268 = 0.1e1 / t267;
    const double t269 = t113 * t268;
    const double t273 = 0.1e1 / t92 / t91 * a * t120;
    const double t276 = t185 * sigma_bb;
    const double t280 = 0.6e1 * t189 * t96 + 0.144e3;
    const double t281 = safe_math::sqrt( t280 );
    const double t282 = 0.1e1 / t281;
    const double t283 = b * t282;
    const double t284 = t257 * a * t283;
    const double t287 = -t114 * t273 / 0.9e1 - 0.2e1 / 0.3e1 * t276 * t284 - t264;
    const double t289 = t265 * t125 - t269 * t287;
    const double t294 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t241 * t127 - t220 - 0.3e1 / 0.8e1 * t5 * t89 * t289 );
    const double t297 = t156 * rho_a;
    const double t299 = 0.1e1 / t36 / t297;
    const double t300 = t153 * t299;
    const double t301 = t43 * sigma_aa;
    const double t308 = 0.1e1 / sigma_aa;
    const double t311 = t59 * expo * t308 / 0.2e1;
    const double t312 = -t151 * t300 * t301 / 0.576e3 + t46 * t33 * t39 / 0.24e2 - t311;
    const double t315 = t51 / t52;
    const double t320 = t319 * t194;
    const double t323 = t315 * t69 / 0.24e2 + t185 * t39 * t320 / 0.4e1 + t311;
    const double t325 = -t178 * t323 + t312 * t73;
    const double t329 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t325 );
    const double t330 = t247 * rho_b;
    const double t332 = 0.1e1 / t92 / t330;
    const double t333 = t153 * t332;
    const double t334 = t99 * sigma_bb;
    const double t341 = 0.1e1 / sigma_bb;
    const double t344 = t112 * expo * t341 / 0.2e1;
    const double t345 = -t151 * t333 * t334 / 0.576e3 + t102 * t33 * t95 / 0.24e2 - t344;
    const double t348 = t51 / t105;
    const double t352 = t319 * t282;
    const double t355 = t348 * t121 / 0.24e2 + t185 * t95 * t352 / 0.4e1 + t344;
    const double t357 = t345 * t125 - t269 * t355;
    const double t361 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t357 );


    eps = t79 + t131;
    vrho_a = t79 + t131 + t6 * ( t206 + t222 );
    vrho_b = t79 + t131 + t6 * ( t235 + t294 );
    vsigma_aa = t6 * t329;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t361;

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
    constexpr double t29 = alpha * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t49 = t28 * t28;
    constexpr double t50 = 0.1e1 / t31;
    constexpr double t51 = t49 * t50;
    constexpr double t63 = b * t49;
    constexpr double t151 = d * alpha * t49;
    constexpr double t153 = 0.1e1 / t31 / t30;
    constexpr double t185 = t28 * t33;
    constexpr double t188 = b * b;
    constexpr double t189 = t188 * t28;
    constexpr double t319 = a * b;
    constexpr double t406 = alpha * alpha;
    constexpr double t407 = d * t406;
    constexpr double t408 = t30 * t30;
    constexpr double t409 = 0.1e1 / t408;
    constexpr double t410 = t407 * t409;
    constexpr double t424 = expo * expo;
    constexpr double t453 = t49 * t153;
    constexpr double t456 = t188 * b;
    constexpr double t692 = a * t456;


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
    const double t40 = t34 * t39;
    const double t43 = safe_math::exp( -t29 * t40 / 0.24e2 );
    const double t46 = ( d * t43 + c ) * t28;
    const double t52 = safe_math::sqrt( sigma_aa );
    const double t54 = 0.1e1 / t36 / rho_a;
    const double t58 = safe_math::pow( t51 * t52 * t54 / 0.12e2, expo );
    const double t59 = f * t58;
    const double t60 = t46 * t40 / 0.24e2 - t59;
    const double t61 = t51 * t52;
    const double t68 = safe_math::log( t63 * t50 * t52 * t54 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t52 * t54 / 0.12e2 ) + 0.1e1 ) );
    const double t69 = t54 * a * t68;
    const double t72 = 0.1e1 + t61 * t69 / 0.12e2 + t59;
    const double t73 = 0.1e1 / t72;
    const double t75 = t60 * t73 + 0.1e1;
    const double t80 = rho_b <= dens_tol;
    const double t81 = -t16;
    const double t83 = piecewise_functor_5( t14, t11, t10, t15, t81 * t7 );
    const double t84 = 0.1e1 + t83;
    const double t85 = t84 <= zeta_tol;
    const double t86 = safe_math::cbrt( t84 );
    const double t88 = piecewise_functor_3( t85, t22, t86 * t84 );
    const double t89 = t88 * t26;
    const double t90 = t33 * sigma_bb;
    const double t91 = rho_b * rho_b;
    const double t92 = safe_math::cbrt( rho_b );
    const double t93 = t92 * t92;
    const double t95 = 0.1e1 / t93 / t91;
    const double t96 = t90 * t95;
    const double t99 = safe_math::exp( -t29 * t96 / 0.24e2 );
    const double t102 = ( d * t99 + c ) * t28;
    const double t105 = safe_math::sqrt( sigma_bb );
    const double t107 = 0.1e1 / t92 / rho_b;
    const double t111 = safe_math::pow( t51 * t105 * t107 / 0.12e2, expo );
    const double t112 = f * t111;
    const double t113 = t102 * t96 / 0.24e2 - t112;
    const double t114 = t51 * t105;
    const double t120 = safe_math::log( t63 * t50 * t105 * t107 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t105 * t107 / 0.12e2 ) + 0.1e1 ) );
    const double t121 = t107 * a * t120;
    const double t124 = 0.1e1 + t114 * t121 / 0.12e2 + t112;
    const double t125 = 0.1e1 / t124;
    const double t127 = t113 * t125 + 0.1e1;
    const double t132 = t6 * t6;
    const double t133 = 0.1e1 / t132;
    const double t134 = t16 * t133;
    const double t136 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t134 );
    const double t139 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t136 );
    const double t140 = t139 * t26;
    const double t144 = t26 * t26;
    const double t145 = 0.1e1 / t144;
    const double t146 = t25 * t145;
    const double t149 = t5 * t146 * t75 / 0.8e1;
    const double t154 = sigma_aa * sigma_aa;
    const double t155 = t153 * t154;
    const double t156 = t35 * t35;
    const double t157 = t156 * t35;
    const double t159 = 0.1e1 / t36 / t157;
    const double t164 = t35 * rho_a;
    const double t166 = 0.1e1 / t37 / t164;
    const double t170 = 0.1e1 / rho_a;
    const double t173 = 0.4e1 / 0.3e1 * t59 * expo * t170;
    const double t174 = t151 * t155 * t159 * t43 / 0.216e3 - t46 * t34 * t166 / 0.9e1 + t173;
    const double t176 = t72 * t72;
    const double t177 = 0.1e1 / t176;
    const double t178 = t60 * t177;
    const double t182 = 0.1e1 / t36 / t35 * a * t68;
    const double t186 = t185 * sigma_aa;
    const double t192 = 0.6e1 * t189 * t40 + 0.144e3;
    const double t193 = safe_math::sqrt( t192 );
    const double t194 = 0.1e1 / t193;
    const double t195 = b * t194;
    const double t196 = t166 * a * t195;
    const double t199 = -t61 * t182 / 0.9e1 - 0.2e1 / 0.3e1 * t186 * t196 - t173;
    const double t201 = t174 * t73 - t178 * t199;
    const double t206 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t140 * t75 - t149 - 0.3e1 / 0.8e1 * t5 * t27 * t201 );
    const double t207 = t81 * t133;
    const double t209 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t207 );
    const double t212 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t209 );
    const double t213 = t212 * t26;
    const double t217 = t88 * t145;
    const double t220 = t5 * t217 * t127 / 0.8e1;
    const double t222 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t213 * t127 - t220 );
    const double t226 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t134 );
    const double t229 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t226 );
    const double t230 = t229 * t26;
    const double t235 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t230 * t75 - t149 );
    const double t237 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t207 );
    const double t240 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t237 );
    const double t241 = t240 * t26;
    const double t245 = sigma_bb * sigma_bb;
    const double t246 = t153 * t245;
    const double t247 = t91 * t91;
    const double t248 = t247 * t91;
    const double t250 = 0.1e1 / t92 / t248;
    const double t255 = t91 * rho_b;
    const double t257 = 0.1e1 / t93 / t255;
    const double t261 = 0.1e1 / rho_b;
    const double t264 = 0.4e1 / 0.3e1 * t112 * expo * t261;
    const double t265 = t151 * t246 * t250 * t99 / 0.216e3 - t102 * t90 * t257 / 0.9e1 + t264;
    const double t267 = t124 * t124;
    const double t268 = 0.1e1 / t267;
    const double t269 = t113 * t268;
    const double t273 = 0.1e1 / t92 / t91 * a * t120;
    const double t276 = t185 * sigma_bb;
    const double t280 = 0.6e1 * t189 * t96 + 0.144e3;
    const double t281 = safe_math::sqrt( t280 );
    const double t282 = 0.1e1 / t281;
    const double t283 = b * t282;
    const double t284 = t257 * a * t283;
    const double t287 = -t114 * t273 / 0.9e1 - 0.2e1 / 0.3e1 * t276 * t284 - t264;
    const double t289 = t265 * t125 - t269 * t287;
    const double t294 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t241 * t127 - t220 - 0.3e1 / 0.8e1 * t5 * t89 * t289 );
    const double t297 = t156 * rho_a;
    const double t299 = 0.1e1 / t36 / t297;
    const double t300 = t153 * t299;
    const double t301 = t43 * sigma_aa;
    const double t308 = 0.1e1 / sigma_aa;
    const double t311 = t59 * expo * t308 / 0.2e1;
    const double t312 = -t151 * t300 * t301 / 0.576e3 + t46 * t33 * t39 / 0.24e2 - t311;
    const double t315 = t51 / t52;
    const double t320 = t319 * t194;
    const double t323 = t315 * t69 / 0.24e2 + t185 * t39 * t320 / 0.4e1 + t311;
    const double t325 = -t178 * t323 + t312 * t73;
    const double t329 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t325 );
    const double t330 = t247 * rho_b;
    const double t332 = 0.1e1 / t92 / t330;
    const double t333 = t153 * t332;
    const double t334 = t99 * sigma_bb;
    const double t341 = 0.1e1 / sigma_bb;
    const double t344 = t112 * expo * t341 / 0.2e1;
    const double t345 = -t151 * t333 * t334 / 0.576e3 + t102 * t33 * t95 / 0.24e2 - t344;
    const double t348 = t51 / t105;
    const double t352 = t319 * t282;
    const double t355 = t348 * t121 / 0.24e2 + t185 * t95 * t352 / 0.4e1 + t344;
    const double t357 = t345 * t125 - t269 * t355;
    const double t361 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t357 );
    const double t364 = t23 * t23;
    const double t365 = 0.1e1 / t364;
    const double t366 = t136 * t136;
    const double t369 = t132 * t6;
    const double t370 = 0.1e1 / t369;
    const double t371 = t16 * t370;
    const double t374 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t133 + 0.2e1 * t371 );
    const double t378 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t365 * t366 + 0.4e1 / 0.3e1 * t23 * t374 );
    const double t379 = t378 * t26;
    const double t383 = t139 * t145;
    const double t385 = t5 * t383 * t75;
    const double t391 = 0.1e1 / t144 / t6;
    const double t392 = t25 * t391;
    const double t395 = t5 * t392 * t75 / 0.12e2;
    const double t397 = t5 * t146 * t201;
    const double t399 = t156 * t164;
    const double t401 = 0.1e1 / t36 / t399;
    const double t411 = t154 * sigma_aa;
    const double t412 = t156 * t156;
    const double t413 = t412 * t35;
    const double t414 = 0.1e1 / t413;
    const double t420 = 0.1e1 / t37 / t156;
    const double t425 = 0.1e1 / t35;
    const double t426 = t424 * t425;
    const double t428 = 0.16e2 / 0.9e1 * t59 * t426;
    const double t431 = 0.4e1 / 0.3e1 * t59 * expo * t425;
    const double t432 = -t151 * t155 * t401 * t43 / 0.24e2 + t410 * t411 * t414 * t43 / 0.324e3 + 0.11e2 / 0.27e2 * t46 * t34 * t420 - t428 - t431;
    const double t434 = t174 * t177;
    const double t438 = 0.1e1 / t176 / t72;
    const double t439 = t60 * t438;
    const double t440 = t199 * t199;
    const double t446 = 0.1e1 / t36 / t164 * a * t68;
    const double t450 = t420 * a * t195;
    const double t454 = t453 * t154;
    const double t458 = 0.1e1 / t193 / t192;
    const double t459 = t456 * t458;
    const double t463 = 0.7e1 / 0.27e2 * t61 * t446 + 0.1e2 / 0.3e1 * t186 * t450 - 0.16e2 / 0.3e1 * t454 * t401 * a * t459 + t428 + t431;
    const double t465 = -t178 * t463 - 0.2e1 * t434 * t199 + t432 * t73 + 0.2e1 * t439 * t440;
    const double t470 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t379 * t75 - t385 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t140 * t201 + t395 - t397 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t465 );
    const double t471 = t86 * t86;
    const double t472 = 0.1e1 / t471;
    const double t473 = t209 * t209;
    const double t476 = t81 * t370;
    const double t479 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t133 + 0.2e1 * t476 );
    const double t483 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.9e1 * t472 * t473 + 0.4e1 / 0.3e1 * t86 * t479 );
    const double t484 = t483 * t26;
    const double t488 = t212 * t145;
    const double t490 = t5 * t488 * t127;
    const double t492 = t88 * t391;
    const double t495 = t5 * t492 * t127 / 0.12e2;
    const double t497 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t484 * t127 - t490 / 0.4e1 + t495 );
    const double t513 = t229 * t145;
    const double t515 = t5 * t513 * t75;
    const double t537 = t240 * t145;
    const double t539 = t5 * t537 * t127;
    const double t546 = t5 * t217 * t289;
    const double t554 = t226 * t226;
    const double t559 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t133 + 0.2e1 * t371 );
    const double t563 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t365 * t554 + 0.4e1 / 0.3e1 * t23 * t559 );
    const double t564 = t563 * t26;
    const double t570 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t564 * t75 - t515 / 0.4e1 + t395 );
    const double t571 = t237 * t237;
    const double t576 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t133 + 0.2e1 * t476 );
    const double t580 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.9e1 * t472 * t571 + 0.4e1 / 0.3e1 * t86 * t576 );
    const double t581 = t580 * t26;
    const double t590 = t247 * t255;
    const double t592 = 0.1e1 / t92 / t590;
    const double t597 = t245 * sigma_bb;
    const double t598 = t247 * t247;
    const double t599 = t598 * t91;
    const double t600 = 0.1e1 / t599;
    const double t606 = 0.1e1 / t93 / t247;
    const double t610 = 0.1e1 / t91;
    const double t611 = t424 * t610;
    const double t613 = 0.16e2 / 0.9e1 * t112 * t611;
    const double t616 = 0.4e1 / 0.3e1 * t112 * expo * t610;
    const double t617 = -t151 * t246 * t592 * t99 / 0.24e2 + t410 * t597 * t600 * t99 / 0.324e3 + 0.11e2 / 0.27e2 * t102 * t90 * t606 - t613 - t616;
    const double t619 = t265 * t268;
    const double t623 = 0.1e1 / t267 / t124;
    const double t624 = t113 * t623;
    const double t625 = t287 * t287;
    const double t631 = 0.1e1 / t92 / t255 * a * t120;
    const double t635 = t606 * a * t283;
    const double t638 = t453 * t245;
    const double t641 = 0.1e1 / t281 / t280;
    const double t642 = t456 * t641;
    const double t646 = 0.7e1 / 0.27e2 * t114 * t631 + 0.1e2 / 0.3e1 * t276 * t635 - 0.16e2 / 0.3e1 * t638 * t592 * a * t642 + t613 + t616;
    const double t648 = t617 * t125 - t269 * t646 - 0.2e1 * t619 * t287 + 0.2e1 * t624 * t625;
    const double t653 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t581 * t127 - t539 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t241 * t289 + t495 - t546 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t89 * t648 );
    const double t661 = t5 * t146 * t325 / 0.8e1;
    const double t662 = t153 * t159;
    const double t666 = t412 * rho_a;
    const double t667 = 0.1e1 / t666;
    const double t675 = t424 * t170;
    const double t678 = 0.2e1 / 0.3e1 * t59 * t675 * t308;
    const double t679 = t151 * t662 * t301 / 0.72e2 - t410 * t667 * t154 * t43 / 0.864e3 - t46 * t33 * t166 / 0.9e1 + t678;
    const double t681 = t312 * t177;
    const double t684 = t323 * t199;
    const double t691 = t453 * t159;
    const double t694 = t692 * t458 * sigma_aa;
    const double t697 = -t315 * t182 / 0.18e2 - t185 * t166 * t320 + 0.2e1 * t691 * t694 - t678;
    const double t699 = -t178 * t697 - t681 * t199 - t434 * t323 + 0.2e1 * t439 * t684 + t679 * t73;
    const double t704 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t140 * t325 - t661 - 0.3e1 / 0.8e1 * t5 * t27 * t699 );
    const double t711 = t5 * t217 * t357 / 0.8e1;
    const double t724 = t153 * t250;
    const double t728 = t598 * rho_b;
    const double t729 = 0.1e1 / t728;
    const double t737 = t424 * t261;
    const double t740 = 0.2e1 / 0.3e1 * t112 * t737 * t341;
    const double t741 = t151 * t724 * t334 / 0.72e2 - t410 * t729 * t245 * t99 / 0.864e3 - t102 * t33 * t257 / 0.9e1 + t740;
    const double t743 = t345 * t268;
    const double t746 = t355 * t287;
    const double t753 = t453 * t250;
    const double t755 = t692 * t641 * sigma_bb;
    const double t758 = -t348 * t273 / 0.18e2 - t185 * t257 * t352 + 0.2e1 * t753 * t755 - t740;
    const double t760 = t741 * t125 - t269 * t758 - t743 * t287 - t619 * t355 + 0.2e1 * t624 * t746;
    const double t765 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t241 * t357 - t711 - 0.3e1 / 0.8e1 * t5 * t89 * t760 );
    const double t767 = 0.1e1 / t412;
    const double t775 = 0.1e1 / t154;
    const double t778 = t59 * t424 * t775 / 0.4e1;
    const double t781 = t59 * expo * t775 / 0.2e1;
    const double t782 = t410 * t767 * t43 * sigma_aa / 0.2304e4 - t151 * t300 * t43 / 0.288e3 - t778 + t781;
    const double t786 = t323 * t323;
    const double t791 = t51 / t52 / sigma_aa;
    const double t794 = t185 * t308;
    const double t796 = t39 * a * t195;
    const double t800 = t692 * t458;
    const double t803 = -t791 * t69 / 0.48e2 + t794 * t796 / 0.8e1 - 0.3e1 / 0.4e1 * t453 * t299 * t800 + t778 - t781;
    const double t805 = -t178 * t803 - 0.2e1 * t681 * t323 + 0.2e1 * t439 * t786 + t782 * t73;
    const double t809 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t805 );
    const double t810 = 0.1e1 / t598;
    const double t818 = 0.1e1 / t245;
    const double t821 = t112 * t424 * t818 / 0.4e1;
    const double t824 = t112 * expo * t818 / 0.2e1;
    const double t825 = t410 * t810 * t99 * sigma_bb / 0.2304e4 - t151 * t333 * t99 / 0.288e3 - t821 + t824;
    const double t829 = t355 * t355;
    const double t834 = t51 / t105 / sigma_bb;
    const double t837 = t185 * t341;
    const double t839 = t95 * a * t283;
    const double t843 = t692 * t641;
    const double t846 = -t834 * t121 / 0.48e2 + t837 * t839 / 0.8e1 - 0.3e1 / 0.4e1 * t453 * t332 * t843 + t821 - t824;
    const double t848 = t825 * t125 - t269 * t846 - 0.2e1 * t743 * t355 + 0.2e1 * t624 * t829;
    const double t852 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t848 );


    v2rho2_aa = 0.2e1 * t206 + 0.2e1 * t222 + t6 * ( t470 + t497 );
    v2rho2_bb = 0.2e1 * t235 + 0.2e1 * t294 + t6 * ( t570 + t653 );
    v2rhosigma_a_aa = t6 * t704 + t329;
    v2rhosigma_b_bb = t6 * t765 + t361;
    v2sigma2_aa_aa = t6 * t809;
    v2sigma2_bb_bb = t6 * t852;
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
    constexpr double t29 = alpha * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t49 = t28 * t28;
    constexpr double t50 = 0.1e1 / t31;
    constexpr double t51 = t49 * t50;
    constexpr double t63 = b * t49;
    constexpr double t151 = d * alpha * t49;
    constexpr double t153 = 0.1e1 / t31 / t30;
    constexpr double t185 = t28 * t33;
    constexpr double t188 = b * b;
    constexpr double t189 = t188 * t28;
    constexpr double t319 = a * b;
    constexpr double t406 = alpha * alpha;
    constexpr double t407 = d * t406;
    constexpr double t408 = t30 * t30;
    constexpr double t409 = 0.1e1 / t408;
    constexpr double t410 = t407 * t409;
    constexpr double t424 = expo * expo;
    constexpr double t453 = t49 * t153;
    constexpr double t456 = t188 * b;
    constexpr double t692 = a * t456;


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
    const double t40 = t34 * t39;
    const double t43 = safe_math::exp( -t29 * t40 / 0.24e2 );
    const double t46 = ( d * t43 + c ) * t28;
    const double t52 = safe_math::sqrt( sigma_aa );
    const double t54 = 0.1e1 / t36 / rho_a;
    const double t58 = safe_math::pow( t51 * t52 * t54 / 0.12e2, expo );
    const double t59 = f * t58;
    const double t60 = t46 * t40 / 0.24e2 - t59;
    const double t61 = t51 * t52;
    const double t68 = safe_math::log( t63 * t50 * t52 * t54 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t52 * t54 / 0.12e2 ) + 0.1e1 ) );
    const double t69 = t54 * a * t68;
    const double t72 = 0.1e1 + t61 * t69 / 0.12e2 + t59;
    const double t73 = 0.1e1 / t72;
    const double t75 = t60 * t73 + 0.1e1;
    const double t79 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t75 );
    const double t80 = rho_b <= dens_tol;
    const double t81 = -t16;
    const double t83 = piecewise_functor_5( t14, t11, t10, t15, t81 * t7 );
    const double t84 = 0.1e1 + t83;
    const double t85 = t84 <= zeta_tol;
    const double t86 = safe_math::cbrt( t84 );
    const double t88 = piecewise_functor_3( t85, t22, t86 * t84 );
    const double t89 = t88 * t26;
    const double t90 = t33 * sigma_bb;
    const double t91 = rho_b * rho_b;
    const double t92 = safe_math::cbrt( rho_b );
    const double t93 = t92 * t92;
    const double t95 = 0.1e1 / t93 / t91;
    const double t96 = t90 * t95;
    const double t99 = safe_math::exp( -t29 * t96 / 0.24e2 );
    const double t102 = ( d * t99 + c ) * t28;
    const double t105 = safe_math::sqrt( sigma_bb );
    const double t107 = 0.1e1 / t92 / rho_b;
    const double t111 = safe_math::pow( t51 * t105 * t107 / 0.12e2, expo );
    const double t112 = f * t111;
    const double t113 = t102 * t96 / 0.24e2 - t112;
    const double t114 = t51 * t105;
    const double t120 = safe_math::log( t63 * t50 * t105 * t107 / 0.12e2 + safe_math::sqrt( square( t63 * t50 * t105 * t107 / 0.12e2 ) + 0.1e1 ) );
    const double t121 = t107 * a * t120;
    const double t124 = 0.1e1 + t114 * t121 / 0.12e2 + t112;
    const double t125 = 0.1e1 / t124;
    const double t127 = t113 * t125 + 0.1e1;
    const double t131 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t127 );
    const double t132 = t6 * t6;
    const double t133 = 0.1e1 / t132;
    const double t134 = t16 * t133;
    const double t136 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t134 );
    const double t139 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t136 );
    const double t140 = t139 * t26;
    const double t144 = t26 * t26;
    const double t145 = 0.1e1 / t144;
    const double t146 = t25 * t145;
    const double t149 = t5 * t146 * t75 / 0.8e1;
    const double t154 = sigma_aa * sigma_aa;
    const double t155 = t153 * t154;
    const double t156 = t35 * t35;
    const double t157 = t156 * t35;
    const double t159 = 0.1e1 / t36 / t157;
    const double t164 = t35 * rho_a;
    const double t166 = 0.1e1 / t37 / t164;
    const double t170 = 0.1e1 / rho_a;
    const double t173 = 0.4e1 / 0.3e1 * t59 * expo * t170;
    const double t174 = t151 * t155 * t159 * t43 / 0.216e3 - t46 * t34 * t166 / 0.9e1 + t173;
    const double t176 = t72 * t72;
    const double t177 = 0.1e1 / t176;
    const double t178 = t60 * t177;
    const double t182 = 0.1e1 / t36 / t35 * a * t68;
    const double t186 = t185 * sigma_aa;
    const double t192 = 0.6e1 * t189 * t40 + 0.144e3;
    const double t193 = safe_math::sqrt( t192 );
    const double t194 = 0.1e1 / t193;
    const double t195 = b * t194;
    const double t196 = t166 * a * t195;
    const double t199 = -t61 * t182 / 0.9e1 - 0.2e1 / 0.3e1 * t186 * t196 - t173;
    const double t201 = t174 * t73 - t178 * t199;
    const double t206 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t140 * t75 - t149 - 0.3e1 / 0.8e1 * t5 * t27 * t201 );
    const double t207 = t81 * t133;
    const double t209 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t207 );
    const double t212 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t209 );
    const double t213 = t212 * t26;
    const double t217 = t88 * t145;
    const double t220 = t5 * t217 * t127 / 0.8e1;
    const double t222 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t213 * t127 - t220 );
    const double t226 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t134 );
    const double t229 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t226 );
    const double t230 = t229 * t26;
    const double t235 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t230 * t75 - t149 );
    const double t237 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t207 );
    const double t240 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.3e1 * t86 * t237 );
    const double t241 = t240 * t26;
    const double t245 = sigma_bb * sigma_bb;
    const double t246 = t153 * t245;
    const double t247 = t91 * t91;
    const double t248 = t247 * t91;
    const double t250 = 0.1e1 / t92 / t248;
    const double t255 = t91 * rho_b;
    const double t257 = 0.1e1 / t93 / t255;
    const double t261 = 0.1e1 / rho_b;
    const double t264 = 0.4e1 / 0.3e1 * t112 * expo * t261;
    const double t265 = t151 * t246 * t250 * t99 / 0.216e3 - t102 * t90 * t257 / 0.9e1 + t264;
    const double t267 = t124 * t124;
    const double t268 = 0.1e1 / t267;
    const double t269 = t113 * t268;
    const double t273 = 0.1e1 / t92 / t91 * a * t120;
    const double t276 = t185 * sigma_bb;
    const double t280 = 0.6e1 * t189 * t96 + 0.144e3;
    const double t281 = safe_math::sqrt( t280 );
    const double t282 = 0.1e1 / t281;
    const double t283 = b * t282;
    const double t284 = t257 * a * t283;
    const double t287 = -t114 * t273 / 0.9e1 - 0.2e1 / 0.3e1 * t276 * t284 - t264;
    const double t289 = t265 * t125 - t269 * t287;
    const double t294 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t241 * t127 - t220 - 0.3e1 / 0.8e1 * t5 * t89 * t289 );
    const double t297 = t156 * rho_a;
    const double t299 = 0.1e1 / t36 / t297;
    const double t300 = t153 * t299;
    const double t301 = t43 * sigma_aa;
    const double t308 = 0.1e1 / sigma_aa;
    const double t311 = t59 * expo * t308 / 0.2e1;
    const double t312 = -t151 * t300 * t301 / 0.576e3 + t46 * t33 * t39 / 0.24e2 - t311;
    const double t315 = t51 / t52;
    const double t320 = t319 * t194;
    const double t323 = t315 * t69 / 0.24e2 + t185 * t39 * t320 / 0.4e1 + t311;
    const double t325 = -t178 * t323 + t312 * t73;
    const double t329 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t325 );
    const double t330 = t247 * rho_b;
    const double t332 = 0.1e1 / t92 / t330;
    const double t333 = t153 * t332;
    const double t334 = t99 * sigma_bb;
    const double t341 = 0.1e1 / sigma_bb;
    const double t344 = t112 * expo * t341 / 0.2e1;
    const double t345 = -t151 * t333 * t334 / 0.576e3 + t102 * t33 * t95 / 0.24e2 - t344;
    const double t348 = t51 / t105;
    const double t352 = t319 * t282;
    const double t355 = t348 * t121 / 0.24e2 + t185 * t95 * t352 / 0.4e1 + t344;
    const double t357 = t345 * t125 - t269 * t355;
    const double t361 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t357 );
    const double t364 = t23 * t23;
    const double t365 = 0.1e1 / t364;
    const double t366 = t136 * t136;
    const double t369 = t132 * t6;
    const double t370 = 0.1e1 / t369;
    const double t371 = t16 * t370;
    const double t374 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t133 + 0.2e1 * t371 );
    const double t378 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t365 * t366 + 0.4e1 / 0.3e1 * t23 * t374 );
    const double t379 = t378 * t26;
    const double t383 = t139 * t145;
    const double t385 = t5 * t383 * t75;
    const double t391 = 0.1e1 / t144 / t6;
    const double t392 = t25 * t391;
    const double t395 = t5 * t392 * t75 / 0.12e2;
    const double t397 = t5 * t146 * t201;
    const double t399 = t156 * t164;
    const double t401 = 0.1e1 / t36 / t399;
    const double t411 = t154 * sigma_aa;
    const double t412 = t156 * t156;
    const double t413 = t412 * t35;
    const double t414 = 0.1e1 / t413;
    const double t420 = 0.1e1 / t37 / t156;
    const double t425 = 0.1e1 / t35;
    const double t426 = t424 * t425;
    const double t428 = 0.16e2 / 0.9e1 * t59 * t426;
    const double t431 = 0.4e1 / 0.3e1 * t59 * expo * t425;
    const double t432 = -t151 * t155 * t401 * t43 / 0.24e2 + t410 * t411 * t414 * t43 / 0.324e3 + 0.11e2 / 0.27e2 * t46 * t34 * t420 - t428 - t431;
    const double t434 = t174 * t177;
    const double t438 = 0.1e1 / t176 / t72;
    const double t439 = t60 * t438;
    const double t440 = t199 * t199;
    const double t446 = 0.1e1 / t36 / t164 * a * t68;
    const double t450 = t420 * a * t195;
    const double t454 = t453 * t154;
    const double t458 = 0.1e1 / t193 / t192;
    const double t459 = t456 * t458;
    const double t463 = 0.7e1 / 0.27e2 * t61 * t446 + 0.1e2 / 0.3e1 * t186 * t450 - 0.16e2 / 0.3e1 * t454 * t401 * a * t459 + t428 + t431;
    const double t465 = -t178 * t463 - 0.2e1 * t434 * t199 + t432 * t73 + 0.2e1 * t439 * t440;
    const double t470 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t379 * t75 - t385 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t140 * t201 + t395 - t397 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t465 );
    const double t471 = t86 * t86;
    const double t472 = 0.1e1 / t471;
    const double t473 = t209 * t209;
    const double t476 = t81 * t370;
    const double t479 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t133 + 0.2e1 * t476 );
    const double t483 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.9e1 * t472 * t473 + 0.4e1 / 0.3e1 * t86 * t479 );
    const double t484 = t483 * t26;
    const double t488 = t212 * t145;
    const double t490 = t5 * t488 * t127;
    const double t492 = t88 * t391;
    const double t495 = t5 * t492 * t127 / 0.12e2;
    const double t497 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t484 * t127 - t490 / 0.4e1 + t495 );
    const double t513 = t229 * t145;
    const double t515 = t5 * t513 * t75;
    const double t537 = t240 * t145;
    const double t539 = t5 * t537 * t127;
    const double t546 = t5 * t217 * t289;
    const double t554 = t226 * t226;
    const double t559 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t133 + 0.2e1 * t371 );
    const double t563 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t365 * t554 + 0.4e1 / 0.3e1 * t23 * t559 );
    const double t564 = t563 * t26;
    const double t570 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t564 * t75 - t515 / 0.4e1 + t395 );
    const double t571 = t237 * t237;
    const double t576 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t133 + 0.2e1 * t476 );
    const double t580 = piecewise_functor_3( t85, 0.0, 0.4e1 / 0.9e1 * t472 * t571 + 0.4e1 / 0.3e1 * t86 * t576 );
    const double t581 = t580 * t26;
    const double t590 = t247 * t255;
    const double t592 = 0.1e1 / t92 / t590;
    const double t597 = t245 * sigma_bb;
    const double t598 = t247 * t247;
    const double t599 = t598 * t91;
    const double t600 = 0.1e1 / t599;
    const double t606 = 0.1e1 / t93 / t247;
    const double t610 = 0.1e1 / t91;
    const double t611 = t424 * t610;
    const double t613 = 0.16e2 / 0.9e1 * t112 * t611;
    const double t616 = 0.4e1 / 0.3e1 * t112 * expo * t610;
    const double t617 = -t151 * t246 * t592 * t99 / 0.24e2 + t410 * t597 * t600 * t99 / 0.324e3 + 0.11e2 / 0.27e2 * t102 * t90 * t606 - t613 - t616;
    const double t619 = t265 * t268;
    const double t623 = 0.1e1 / t267 / t124;
    const double t624 = t113 * t623;
    const double t625 = t287 * t287;
    const double t631 = 0.1e1 / t92 / t255 * a * t120;
    const double t635 = t606 * a * t283;
    const double t638 = t453 * t245;
    const double t641 = 0.1e1 / t281 / t280;
    const double t642 = t456 * t641;
    const double t646 = 0.7e1 / 0.27e2 * t114 * t631 + 0.1e2 / 0.3e1 * t276 * t635 - 0.16e2 / 0.3e1 * t638 * t592 * a * t642 + t613 + t616;
    const double t648 = t617 * t125 - t269 * t646 - 0.2e1 * t619 * t287 + 0.2e1 * t624 * t625;
    const double t653 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t581 * t127 - t539 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t241 * t289 + t495 - t546 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t89 * t648 );
    const double t661 = t5 * t146 * t325 / 0.8e1;
    const double t662 = t153 * t159;
    const double t666 = t412 * rho_a;
    const double t667 = 0.1e1 / t666;
    const double t675 = t424 * t170;
    const double t678 = 0.2e1 / 0.3e1 * t59 * t675 * t308;
    const double t679 = t151 * t662 * t301 / 0.72e2 - t410 * t667 * t154 * t43 / 0.864e3 - t46 * t33 * t166 / 0.9e1 + t678;
    const double t681 = t312 * t177;
    const double t684 = t323 * t199;
    const double t691 = t453 * t159;
    const double t694 = t692 * t458 * sigma_aa;
    const double t697 = -t315 * t182 / 0.18e2 - t185 * t166 * t320 + 0.2e1 * t691 * t694 - t678;
    const double t699 = -t178 * t697 - t681 * t199 - t434 * t323 + 0.2e1 * t439 * t684 + t679 * t73;
    const double t704 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t140 * t325 - t661 - 0.3e1 / 0.8e1 * t5 * t27 * t699 );
    const double t711 = t5 * t217 * t357 / 0.8e1;
    const double t724 = t153 * t250;
    const double t728 = t598 * rho_b;
    const double t729 = 0.1e1 / t728;
    const double t737 = t424 * t261;
    const double t740 = 0.2e1 / 0.3e1 * t112 * t737 * t341;
    const double t741 = t151 * t724 * t334 / 0.72e2 - t410 * t729 * t245 * t99 / 0.864e3 - t102 * t33 * t257 / 0.9e1 + t740;
    const double t743 = t345 * t268;
    const double t746 = t355 * t287;
    const double t753 = t453 * t250;
    const double t755 = t692 * t641 * sigma_bb;
    const double t758 = -t348 * t273 / 0.18e2 - t185 * t257 * t352 + 0.2e1 * t753 * t755 - t740;
    const double t760 = t741 * t125 - t269 * t758 - t743 * t287 - t619 * t355 + 0.2e1 * t624 * t746;
    const double t765 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t241 * t357 - t711 - 0.3e1 / 0.8e1 * t5 * t89 * t760 );
    const double t767 = 0.1e1 / t412;
    const double t775 = 0.1e1 / t154;
    const double t778 = t59 * t424 * t775 / 0.4e1;
    const double t781 = t59 * expo * t775 / 0.2e1;
    const double t782 = t410 * t767 * t43 * sigma_aa / 0.2304e4 - t151 * t300 * t43 / 0.288e3 - t778 + t781;
    const double t786 = t323 * t323;
    const double t791 = t51 / t52 / sigma_aa;
    const double t794 = t185 * t308;
    const double t796 = t39 * a * t195;
    const double t800 = t692 * t458;
    const double t803 = -t791 * t69 / 0.48e2 + t794 * t796 / 0.8e1 - 0.3e1 / 0.4e1 * t453 * t299 * t800 + t778 - t781;
    const double t805 = -t178 * t803 - 0.2e1 * t681 * t323 + 0.2e1 * t439 * t786 + t782 * t73;
    const double t809 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t805 );
    const double t810 = 0.1e1 / t598;
    const double t818 = 0.1e1 / t245;
    const double t821 = t112 * t424 * t818 / 0.4e1;
    const double t824 = t112 * expo * t818 / 0.2e1;
    const double t825 = t410 * t810 * t99 * sigma_bb / 0.2304e4 - t151 * t333 * t99 / 0.288e3 - t821 + t824;
    const double t829 = t355 * t355;
    const double t834 = t51 / t105 / sigma_bb;
    const double t837 = t185 * t341;
    const double t839 = t95 * a * t283;
    const double t843 = t692 * t641;
    const double t846 = -t834 * t121 / 0.48e2 + t837 * t839 / 0.8e1 - 0.3e1 / 0.4e1 * t453 * t332 * t843 + t821 - t824;
    const double t848 = t825 * t125 - t269 * t846 - 0.2e1 * t743 * t355 + 0.2e1 * t624 * t829;
    const double t852 = piecewise_functor_3( t80, 0.0, -0.3e1 / 0.8e1 * t5 * t89 * t848 );


    vrho_a = t79 + t131 + t6 * ( t206 + t222 );
    vrho_b = t79 + t131 + t6 * ( t235 + t294 );
    vsigma_aa = t6 * t329;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t361;
    v2rho2_aa = 0.2e1 * t206 + 0.2e1 * t222 + t6 * ( t470 + t497 );
    v2rho2_bb = 0.2e1 * t235 + 0.2e1 * t294 + t6 * ( t570 + t653 );
    v2rhosigma_a_aa = t6 * t704 + t329;
    v2rhosigma_b_bb = t6 * t765 + t361;
    v2sigma2_aa_aa = t6 * t809;
    v2sigma2_bb_bb = t6 * t852;
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

struct BuiltinMPW91_X : detail::BuiltinKernelImpl< BuiltinMPW91_X > {

  BuiltinMPW91_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinMPW91_X >(p) { }
  
  virtual ~BuiltinMPW91_X() = default;

};



} // namespace ExchCXX
