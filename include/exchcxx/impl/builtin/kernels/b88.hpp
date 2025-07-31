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
struct kernel_traits< BuiltinB88 > :
  public gga_screening_interface< BuiltinB88 > {

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


  static constexpr double beta = 0.0042;
  static constexpr double gamma = 6.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;


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
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t56 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t52 );


    eps = 0.2e1 * t56;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;
    constexpr double t99 = t21 * t24;
    constexpr double t100 = t25 * t29;


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
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t56 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t52 );
    const double t58 = t17 / t32;
    const double t62 = t31 * rho;
    const double t64 = 0.1e1 / t32 / t62;
    const double t65 = t64 * t47;
    const double t69 = t46 * t46;
    const double t70 = 0.1e1 / t69;
    const double t71 = t34 * t70;
    const double t75 = t28 / t18 / t31 * t43;
    const double t77 = t35 * sigma;
    const double t78 = t29 * t64;
    const double t80 = t30 * t34 + 0.1e1;
    const double t81 = safe_math::sqrt( t80 );
    const double t82 = 0.1e1 / t81;
    const double t83 = t78 * t82;
    const double t86 = -0.4e1 / 0.3e1 * t37 * t75 - 0.4e1 / 0.3e1 * t77 * t83;
    const double t91 = -0.16e2 / 0.27e2 * t27 * t30 * t65 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t86;
    const double t96 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t52 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t91 );
    const double t104 = t35 / t36;
    const double t106 = t29 * t34;
    const double t107 = t106 * t82;
    const double t110 = t104 * t44 / 0.2e1 + t35 * t107 / 0.2e1;
    const double t115 = -0.2e1 / 0.9e1 * t27 * t30 * t71 * t110 + 0.2e1 / 0.9e1 * t99 * t100 * t48;
    const double t119 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t115 );


    eps = 0.2e1 * t56;
    vrho = 0.2e1 * rho * t96 + 0.2e1 * t56;
    vsigma = 0.2e1 * rho * t119;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;
    constexpr double t99 = t21 * t24;
    constexpr double t100 = t25 * t29;
    constexpr double t210 = t35 * t28;


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
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t58 = t17 / t32;
    const double t62 = t31 * rho;
    const double t64 = 0.1e1 / t32 / t62;
    const double t65 = t64 * t47;
    const double t69 = t46 * t46;
    const double t70 = 0.1e1 / t69;
    const double t71 = t34 * t70;
    const double t75 = t28 / t18 / t31 * t43;
    const double t77 = t35 * sigma;
    const double t78 = t29 * t64;
    const double t80 = t30 * t34 + 0.1e1;
    const double t81 = safe_math::sqrt( t80 );
    const double t82 = 0.1e1 / t81;
    const double t83 = t78 * t82;
    const double t86 = -0.4e1 / 0.3e1 * t37 * t75 - 0.4e1 / 0.3e1 * t77 * t83;
    const double t91 = -0.16e2 / 0.27e2 * t27 * t30 * t65 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t86;
    const double t96 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t52 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t91 );
    const double t104 = t35 / t36;
    const double t106 = t29 * t34;
    const double t107 = t106 * t82;
    const double t110 = t104 * t44 / 0.2e1 + t35 * t107 / 0.2e1;
    const double t115 = -0.2e1 / 0.9e1 * t27 * t30 * t71 * t110 + 0.2e1 / 0.9e1 * t99 * t100 * t48;
    const double t119 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t115 );
    const double t124 = t17 / t32 / rho;
    const double t131 = t31 * t31;
    const double t133 = 0.1e1 / t32 / t131;
    const double t134 = t133 * t47;
    const double t138 = t64 * t70;
    const double t144 = 0.1e1 / t69 / t46;
    const double t145 = t34 * t144;
    const double t146 = t86 * t86;
    const double t154 = t28 / t18 / t62 * t43;
    const double t157 = t29 * t133;
    const double t158 = t157 * t82;
    const double t161 = sigma * sigma;
    const double t162 = t35 * t161;
    const double t165 = 0.1e1 / t18 / t131 / t62;
    const double t168 = 0.1e1 / t81 / t80;
    const double t169 = t28 * t165 * t168;
    const double t172 = 0.28e2 / 0.9e1 * t37 * t154 + 0.2e2 / 0.3e1 * t77 * t158 - 0.32e2 / 0.9e1 * t162 * t169;
    const double t177 = 0.176e3 / 0.81e2 * t27 * t30 * t134 + 0.32e2 / 0.27e2 * t27 * t30 * t138 * t86 + 0.4e1 / 0.9e1 * t27 * t30 * t145 * t146 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t172;
    const double t182 = piecewise_functor_3( t2, 0.0, t6 * t124 * t52 / 0.12e2 - t6 * t58 * t91 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t177 );
    const double t191 = t70 * t86;
    const double t200 = t21 * t26 * sigma;
    const double t201 = t144 * t110;
    const double t202 = t201 * t86;
    const double t203 = t106 * t202;
    const double t211 = t131 * t31;
    const double t213 = 0.1e1 / t18 / t211;
    const double t218 = -0.2e1 / 0.3e1 * t104 * t75 - 0.2e1 * t35 * t83 + 0.4e1 / 0.3e1 * t210 * t213 * t168 * sigma;
    const double t223 = -0.16e2 / 0.27e2 * t99 * t100 * t65 - 0.2e1 / 0.9e1 * t27 * t106 * t191 + 0.16e2 / 0.27e2 * t27 * t30 * t138 * t110 + 0.4e1 / 0.9e1 * t200 * t203 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t218;
    const double t228 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t115 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t223 );
    const double t231 = t70 * t110;
    const double t235 = t110 * t110;
    const double t242 = t35 / t36 / sigma;
    const double t245 = 0.1e1 / sigma;
    const double t246 = t35 * t245;
    const double t249 = t131 * rho;
    const double t252 = t28 / t18 / t249;
    const double t253 = t252 * t168;
    const double t256 = -t242 * t44 / 0.4e1 + t246 * t107 / 0.4e1 - t35 * t253 / 0.2e1;
    const double t261 = -0.4e1 / 0.9e1 * t27 * t106 * t231 + 0.4e1 / 0.9e1 * t27 * t30 * t145 * t235 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t256;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t261 );


    v2rho2 = 0.2e1 * rho * t182 + 0.4e1 * t96;
    v2rhosigma = 0.2e1 * rho * t228 + 0.2e1 * t119;
    v2sigma2 = 0.2e1 * rho * t265;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double sigma, double& vrho, double& vsigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;
    constexpr double t99 = t21 * t24;
    constexpr double t100 = t25 * t29;
    constexpr double t210 = t35 * t28;


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
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t56 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t52 );
    const double t58 = t17 / t32;
    const double t62 = t31 * rho;
    const double t64 = 0.1e1 / t32 / t62;
    const double t65 = t64 * t47;
    const double t69 = t46 * t46;
    const double t70 = 0.1e1 / t69;
    const double t71 = t34 * t70;
    const double t75 = t28 / t18 / t31 * t43;
    const double t77 = t35 * sigma;
    const double t78 = t29 * t64;
    const double t80 = t30 * t34 + 0.1e1;
    const double t81 = safe_math::sqrt( t80 );
    const double t82 = 0.1e1 / t81;
    const double t83 = t78 * t82;
    const double t86 = -0.4e1 / 0.3e1 * t37 * t75 - 0.4e1 / 0.3e1 * t77 * t83;
    const double t91 = -0.16e2 / 0.27e2 * t27 * t30 * t65 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t86;
    const double t96 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t52 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t91 );
    const double t104 = t35 / t36;
    const double t106 = t29 * t34;
    const double t107 = t106 * t82;
    const double t110 = t104 * t44 / 0.2e1 + t35 * t107 / 0.2e1;
    const double t115 = -0.2e1 / 0.9e1 * t27 * t30 * t71 * t110 + 0.2e1 / 0.9e1 * t99 * t100 * t48;
    const double t119 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t115 );
    const double t124 = t17 / t32 / rho;
    const double t131 = t31 * t31;
    const double t133 = 0.1e1 / t32 / t131;
    const double t134 = t133 * t47;
    const double t138 = t64 * t70;
    const double t144 = 0.1e1 / t69 / t46;
    const double t145 = t34 * t144;
    const double t146 = t86 * t86;
    const double t154 = t28 / t18 / t62 * t43;
    const double t157 = t29 * t133;
    const double t158 = t157 * t82;
    const double t161 = sigma * sigma;
    const double t162 = t35 * t161;
    const double t165 = 0.1e1 / t18 / t131 / t62;
    const double t168 = 0.1e1 / t81 / t80;
    const double t169 = t28 * t165 * t168;
    const double t172 = 0.28e2 / 0.9e1 * t37 * t154 + 0.2e2 / 0.3e1 * t77 * t158 - 0.32e2 / 0.9e1 * t162 * t169;
    const double t177 = 0.176e3 / 0.81e2 * t27 * t30 * t134 + 0.32e2 / 0.27e2 * t27 * t30 * t138 * t86 + 0.4e1 / 0.9e1 * t27 * t30 * t145 * t146 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t172;
    const double t182 = piecewise_functor_3( t2, 0.0, t6 * t124 * t52 / 0.12e2 - t6 * t58 * t91 / 0.4e1 - 0.3e1 / 0.8e1 * t6 * t19 * t177 );
    const double t191 = t70 * t86;
    const double t200 = t21 * t26 * sigma;
    const double t201 = t144 * t110;
    const double t202 = t201 * t86;
    const double t203 = t106 * t202;
    const double t211 = t131 * t31;
    const double t213 = 0.1e1 / t18 / t211;
    const double t218 = -0.2e1 / 0.3e1 * t104 * t75 - 0.2e1 * t35 * t83 + 0.4e1 / 0.3e1 * t210 * t213 * t168 * sigma;
    const double t223 = -0.16e2 / 0.27e2 * t99 * t100 * t65 - 0.2e1 / 0.9e1 * t27 * t106 * t191 + 0.16e2 / 0.27e2 * t27 * t30 * t138 * t110 + 0.4e1 / 0.9e1 * t200 * t203 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t218;
    const double t228 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t115 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t223 );
    const double t231 = t70 * t110;
    const double t235 = t110 * t110;
    const double t242 = t35 / t36 / sigma;
    const double t245 = 0.1e1 / sigma;
    const double t246 = t35 * t245;
    const double t249 = t131 * rho;
    const double t252 = t28 / t18 / t249;
    const double t253 = t252 * t168;
    const double t256 = -t242 * t44 / 0.4e1 + t246 * t107 / 0.4e1 - t35 * t253 / 0.2e1;
    const double t261 = -0.4e1 / 0.9e1 * t27 * t106 * t231 + 0.4e1 / 0.9e1 * t27 * t30 * t145 * t235 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t256;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t261 );


    vrho = 0.2e1 * rho * t96 + 0.2e1 * t56;
    vsigma = 0.2e1 * rho * t119;
    v2rho2 = 0.2e1 * rho * t182 + 0.4e1 * t96;
    v2rhosigma = 0.2e1 * rho * t228 + 0.2e1 * t119;
    v2sigma2 = 0.2e1 * rho * t265;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;


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
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t59 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t55 );
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t93 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t89 );


    eps = t59 + t93;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;
    constexpr double t119 = t32 * t34;
    constexpr double t120 = t29 * t119;


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
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t59 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t55 );
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t93 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t89 );
    const double t94 = t6 * t6;
    const double t95 = 0.1e1 / t94;
    const double t96 = t16 * t95;
    const double t98 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t96 );
    const double t101 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t98 );
    const double t102 = t101 * t26;
    const double t106 = t26 * t26;
    const double t107 = 0.1e1 / t106;
    const double t108 = t25 * t107;
    const double t111 = t5 * t108 * t55 / 0.8e1;
    const double t112 = t36 * rho_a;
    const double t114 = 0.1e1 / t38 / t112;
    const double t121 = sigma_aa * t40;
    const double t122 = t49 * t49;
    const double t123 = 0.1e1 / t122;
    const double t125 = 0.1e1 / t37 / t36;
    const double t129 = sigma_aa * t114;
    const double t130 = t121 + 0.1e1;
    const double t131 = safe_math::sqrt( t130 );
    const double t132 = 0.1e1 / t131;
    const double t136 = -0.4e1 / 0.3e1 * t41 * t42 * t125 * t46 - 0.4e1 / 0.3e1 * t41 * t129 * t132;
    const double t137 = t123 * t136;
    const double t141 = -0.16e2 / 0.27e2 * t33 * t35 * t114 * t50 - 0.2e1 / 0.9e1 * t120 * t121 * t137;
    const double t146 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t55 - t111 - 0.3e1 / 0.8e1 * t5 * t27 * t141 );
    const double t147 = t61 * t95;
    const double t149 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t147 );
    const double t152 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t149 );
    const double t153 = t152 * t26;
    const double t157 = t68 * t107;
    const double t160 = t5 * t157 * t89 / 0.8e1;
    const double t162 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t153 * t89 - t160 );
    const double t166 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t96 );
    const double t169 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t166 );
    const double t170 = t169 * t26;
    const double t175 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t170 * t55 - t111 );
    const double t177 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t147 );
    const double t180 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t177 );
    const double t181 = t180 * t26;
    const double t185 = t71 * rho_b;
    const double t187 = 0.1e1 / t73 / t185;
    const double t192 = sigma_bb * t75;
    const double t193 = t83 * t83;
    const double t194 = 0.1e1 / t193;
    const double t196 = 0.1e1 / t72 / t71;
    const double t200 = sigma_bb * t187;
    const double t201 = t192 + 0.1e1;
    const double t202 = safe_math::sqrt( t201 );
    const double t203 = 0.1e1 / t202;
    const double t207 = -0.4e1 / 0.3e1 * t41 * t76 * t196 * t80 - 0.4e1 / 0.3e1 * t41 * t200 * t203;
    const double t208 = t194 * t207;
    const double t212 = -0.16e2 / 0.27e2 * t33 * t70 * t187 * t84 - 0.2e1 / 0.9e1 * t120 * t192 * t208;
    const double t217 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t89 - t160 - 0.3e1 / 0.8e1 * t5 * t69 * t212 );
    const double t220 = t34 * t40;
    const double t223 = 0.1e1 / t42;
    const double t230 = t41 * t223 * t44 * t46 / 0.2e1 + t41 * t40 * t132 / 0.2e1;
    const double t231 = t123 * t230;
    const double t235 = -0.2e1 / 0.9e1 * t120 * t121 * t231 + 0.2e1 / 0.9e1 * t33 * t220 * t50;
    const double t239 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t235 );
    const double t240 = t34 * t75;
    const double t243 = 0.1e1 / t76;
    const double t250 = t41 * t243 * t78 * t80 / 0.2e1 + t41 * t75 * t203 / 0.2e1;
    const double t251 = t194 * t250;
    const double t255 = -0.2e1 / 0.9e1 * t120 * t192 * t251 + 0.2e1 / 0.9e1 * t33 * t240 * t84;
    const double t259 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t255 );


    eps = t59 + t93;
    vrho_a = t59 + t93 + t6 * ( t146 + t162 );
    vrho_b = t59 + t93 + t6 * ( t175 + t217 );
    vsigma_aa = t6 * t239;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t259;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;
    constexpr double t119 = t32 * t34;
    constexpr double t120 = t29 * t119;


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
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t94 = t6 * t6;
    const double t95 = 0.1e1 / t94;
    const double t96 = t16 * t95;
    const double t98 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t96 );
    const double t101 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t98 );
    const double t102 = t101 * t26;
    const double t106 = t26 * t26;
    const double t107 = 0.1e1 / t106;
    const double t108 = t25 * t107;
    const double t111 = t5 * t108 * t55 / 0.8e1;
    const double t112 = t36 * rho_a;
    const double t114 = 0.1e1 / t38 / t112;
    const double t121 = sigma_aa * t40;
    const double t122 = t49 * t49;
    const double t123 = 0.1e1 / t122;
    const double t125 = 0.1e1 / t37 / t36;
    const double t129 = sigma_aa * t114;
    const double t130 = t121 + 0.1e1;
    const double t131 = safe_math::sqrt( t130 );
    const double t132 = 0.1e1 / t131;
    const double t136 = -0.4e1 / 0.3e1 * t41 * t42 * t125 * t46 - 0.4e1 / 0.3e1 * t41 * t129 * t132;
    const double t137 = t123 * t136;
    const double t141 = -0.16e2 / 0.27e2 * t33 * t35 * t114 * t50 - 0.2e1 / 0.9e1 * t120 * t121 * t137;
    const double t146 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t55 - t111 - 0.3e1 / 0.8e1 * t5 * t27 * t141 );
    const double t147 = t61 * t95;
    const double t149 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t147 );
    const double t152 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t149 );
    const double t153 = t152 * t26;
    const double t157 = t68 * t107;
    const double t160 = t5 * t157 * t89 / 0.8e1;
    const double t162 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t153 * t89 - t160 );
    const double t166 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t96 );
    const double t169 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t166 );
    const double t170 = t169 * t26;
    const double t175 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t170 * t55 - t111 );
    const double t177 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t147 );
    const double t180 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t177 );
    const double t181 = t180 * t26;
    const double t185 = t71 * rho_b;
    const double t187 = 0.1e1 / t73 / t185;
    const double t192 = sigma_bb * t75;
    const double t193 = t83 * t83;
    const double t194 = 0.1e1 / t193;
    const double t196 = 0.1e1 / t72 / t71;
    const double t200 = sigma_bb * t187;
    const double t201 = t192 + 0.1e1;
    const double t202 = safe_math::sqrt( t201 );
    const double t203 = 0.1e1 / t202;
    const double t207 = -0.4e1 / 0.3e1 * t41 * t76 * t196 * t80 - 0.4e1 / 0.3e1 * t41 * t200 * t203;
    const double t208 = t194 * t207;
    const double t212 = -0.16e2 / 0.27e2 * t33 * t70 * t187 * t84 - 0.2e1 / 0.9e1 * t120 * t192 * t208;
    const double t217 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t89 - t160 - 0.3e1 / 0.8e1 * t5 * t69 * t212 );
    const double t220 = t34 * t40;
    const double t223 = 0.1e1 / t42;
    const double t230 = t41 * t223 * t44 * t46 / 0.2e1 + t41 * t40 * t132 / 0.2e1;
    const double t231 = t123 * t230;
    const double t235 = -0.2e1 / 0.9e1 * t120 * t121 * t231 + 0.2e1 / 0.9e1 * t33 * t220 * t50;
    const double t239 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t235 );
    const double t240 = t34 * t75;
    const double t243 = 0.1e1 / t76;
    const double t250 = t41 * t243 * t78 * t80 / 0.2e1 + t41 * t75 * t203 / 0.2e1;
    const double t251 = t194 * t250;
    const double t255 = -0.2e1 / 0.9e1 * t120 * t192 * t251 + 0.2e1 / 0.9e1 * t33 * t240 * t84;
    const double t259 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t255 );
    const double t262 = t23 * t23;
    const double t263 = 0.1e1 / t262;
    const double t264 = t98 * t98;
    const double t267 = t94 * t6;
    const double t268 = 0.1e1 / t267;
    const double t269 = t16 * t268;
    const double t272 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t95 + 0.2e1 * t269 );
    const double t276 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t263 * t264 + 0.4e1 / 0.3e1 * t23 * t272 );
    const double t277 = t276 * t26;
    const double t281 = t101 * t107;
    const double t283 = t5 * t281 * t55;
    const double t289 = 0.1e1 / t106 / t6;
    const double t290 = t25 * t289;
    const double t293 = t5 * t290 * t55 / 0.12e2;
    const double t295 = t5 * t108 * t141;
    const double t297 = t36 * t36;
    const double t299 = 0.1e1 / t38 / t297;
    const double t308 = 0.1e1 / t122 / t49;
    const double t309 = t136 * t136;
    const double t310 = t308 * t309;
    const double t315 = 0.1e1 / t37 / t112;
    const double t320 = sigma_aa * t299;
    const double t324 = sigma_aa * sigma_aa;
    const double t327 = 0.1e1 / t37 / t297 / t112;
    const double t330 = 0.1e1 / t131 / t130;
    const double t334 = 0.28e2 / 0.9e1 * t41 * t42 * t315 * t46 + 0.2e2 / 0.3e1 * t41 * t320 * t132 - 0.16e2 / 0.9e1 * t41 * t324 * t327 * t330;
    const double t335 = t123 * t334;
    const double t339 = 0.176e3 / 0.81e2 * t33 * t35 * t299 * t50 + 0.32e2 / 0.27e2 * t120 * t129 * t137 + 0.4e1 / 0.9e1 * t120 * t121 * t310 - 0.2e1 / 0.9e1 * t120 * t121 * t335;
    const double t344 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t277 * t55 - t283 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t102 * t141 + t293 - t295 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t339 );
    const double t345 = t66 * t66;
    const double t346 = 0.1e1 / t345;
    const double t347 = t149 * t149;
    const double t350 = t61 * t268;
    const double t353 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t95 + 0.2e1 * t350 );
    const double t357 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t346 * t347 + 0.4e1 / 0.3e1 * t66 * t353 );
    const double t358 = t357 * t26;
    const double t362 = t152 * t107;
    const double t364 = t5 * t362 * t89;
    const double t366 = t68 * t289;
    const double t369 = t5 * t366 * t89 / 0.12e2;
    const double t371 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t358 * t89 - t364 / 0.4e1 + t369 );
    const double t387 = t169 * t107;
    const double t389 = t5 * t387 * t55;
    const double t411 = t180 * t107;
    const double t413 = t5 * t411 * t89;
    const double t420 = t5 * t157 * t212;
    const double t428 = t166 * t166;
    const double t433 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t95 + 0.2e1 * t269 );
    const double t437 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t263 * t428 + 0.4e1 / 0.3e1 * t23 * t433 );
    const double t438 = t437 * t26;
    const double t444 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t438 * t55 - t389 / 0.4e1 + t293 );
    const double t445 = t177 * t177;
    const double t450 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t95 + 0.2e1 * t350 );
    const double t454 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t346 * t445 + 0.4e1 / 0.3e1 * t66 * t450 );
    const double t455 = t454 * t26;
    const double t464 = t71 * t71;
    const double t466 = 0.1e1 / t73 / t464;
    const double t475 = 0.1e1 / t193 / t83;
    const double t476 = t207 * t207;
    const double t477 = t475 * t476;
    const double t482 = 0.1e1 / t72 / t185;
    const double t487 = sigma_bb * t466;
    const double t491 = sigma_bb * sigma_bb;
    const double t494 = 0.1e1 / t72 / t464 / t185;
    const double t497 = 0.1e1 / t202 / t201;
    const double t501 = 0.28e2 / 0.9e1 * t41 * t76 * t482 * t80 + 0.2e2 / 0.3e1 * t41 * t487 * t203 - 0.16e2 / 0.9e1 * t41 * t491 * t494 * t497;
    const double t502 = t194 * t501;
    const double t506 = 0.176e3 / 0.81e2 * t33 * t70 * t466 * t84 + 0.32e2 / 0.27e2 * t120 * t200 * t208 + 0.4e1 / 0.9e1 * t120 * t192 * t477 - 0.2e1 / 0.9e1 * t120 * t192 * t502;
    const double t511 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t455 * t89 - t413 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t181 * t212 + t369 - t420 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t69 * t506 );
    const double t519 = t5 * t108 * t235 / 0.8e1;
    const double t520 = t34 * t114;
    const double t530 = t308 * t230;
    const double t531 = t530 * t136;
    const double t542 = t297 * t36;
    const double t544 = 0.1e1 / t37 / t542;
    const double t545 = t544 * t330;
    const double t549 = -0.2e1 / 0.3e1 * t41 * t223 * t125 * t46 - 0.2e1 * t41 * t114 * t132 + 0.2e1 / 0.3e1 * t41 * t545 * sigma_aa;
    const double t550 = t123 * t549;
    const double t554 = -0.16e2 / 0.27e2 * t33 * t520 * t50 - 0.2e1 / 0.9e1 * t33 * t220 * t137 + 0.16e2 / 0.27e2 * t120 * t129 * t231 + 0.4e1 / 0.9e1 * t120 * t121 * t531 - 0.2e1 / 0.9e1 * t120 * t121 * t550;
    const double t559 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t235 - t519 - 0.3e1 / 0.8e1 * t5 * t27 * t554 );
    const double t566 = t5 * t157 * t255 / 0.8e1;
    const double t579 = t34 * t187;
    const double t589 = t475 * t250;
    const double t590 = t589 * t207;
    const double t601 = t464 * t71;
    const double t603 = 0.1e1 / t72 / t601;
    const double t604 = t603 * t497;
    const double t608 = -0.2e1 / 0.3e1 * t41 * t243 * t196 * t80 - 0.2e1 * t41 * t187 * t203 + 0.2e1 / 0.3e1 * t41 * t604 * sigma_bb;
    const double t609 = t194 * t608;
    const double t613 = -0.16e2 / 0.27e2 * t33 * t579 * t84 - 0.2e1 / 0.9e1 * t33 * t240 * t208 + 0.16e2 / 0.27e2 * t120 * t200 * t251 + 0.4e1 / 0.9e1 * t120 * t192 * t590 - 0.2e1 / 0.9e1 * t120 * t192 * t609;
    const double t618 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t255 - t566 - 0.3e1 / 0.8e1 * t5 * t69 * t613 );
    const double t623 = t230 * t230;
    const double t624 = t308 * t623;
    const double t629 = 0.1e1 / t42 / sigma_aa;
    const double t633 = 0.1e1 / sigma_aa;
    const double t637 = t297 * rho_a;
    const double t639 = 0.1e1 / t37 / t637;
    const double t643 = t41 * t633 * t40 * t132 / 0.4e1 - t41 * t629 * t44 * t46 / 0.4e1 - t41 * t639 * t330 / 0.4e1;
    const double t644 = t123 * t643;
    const double t648 = -0.4e1 / 0.9e1 * t33 * t220 * t231 + 0.4e1 / 0.9e1 * t120 * t121 * t624 - 0.2e1 / 0.9e1 * t120 * t121 * t644;
    const double t652 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t648 );
    const double t656 = t250 * t250;
    const double t657 = t475 * t656;
    const double t662 = 0.1e1 / t76 / sigma_bb;
    const double t666 = 0.1e1 / sigma_bb;
    const double t670 = t464 * rho_b;
    const double t672 = 0.1e1 / t72 / t670;
    const double t676 = t41 * t666 * t75 * t203 / 0.4e1 - t41 * t662 * t78 * t80 / 0.4e1 - t41 * t672 * t497 / 0.4e1;
    const double t677 = t194 * t676;
    const double t681 = -0.4e1 / 0.9e1 * t33 * t240 * t251 + 0.4e1 / 0.9e1 * t120 * t192 * t657 - 0.2e1 / 0.9e1 * t120 * t192 * t677;
    const double t685 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t681 );


    v2rho2_aa = 0.2e1 * t146 + 0.2e1 * t162 + t6 * ( t344 + t371 );
    v2rho2_bb = 0.2e1 * t175 + 0.2e1 * t217 + t6 * ( t444 + t511 );
    v2rhosigma_a_aa = t6 * t559 + t239;
    v2rhosigma_b_bb = t6 * t618 + t259;
    v2sigma2_aa_aa = t6 * t652;
    v2sigma2_bb_bb = t6 * t685;
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
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;
    constexpr double t119 = t32 * t34;
    constexpr double t120 = t29 * t119;


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
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t59 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t55 );
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t93 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t89 );
    const double t94 = t6 * t6;
    const double t95 = 0.1e1 / t94;
    const double t96 = t16 * t95;
    const double t98 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t96 );
    const double t101 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t98 );
    const double t102 = t101 * t26;
    const double t106 = t26 * t26;
    const double t107 = 0.1e1 / t106;
    const double t108 = t25 * t107;
    const double t111 = t5 * t108 * t55 / 0.8e1;
    const double t112 = t36 * rho_a;
    const double t114 = 0.1e1 / t38 / t112;
    const double t121 = sigma_aa * t40;
    const double t122 = t49 * t49;
    const double t123 = 0.1e1 / t122;
    const double t125 = 0.1e1 / t37 / t36;
    const double t129 = sigma_aa * t114;
    const double t130 = t121 + 0.1e1;
    const double t131 = safe_math::sqrt( t130 );
    const double t132 = 0.1e1 / t131;
    const double t136 = -0.4e1 / 0.3e1 * t41 * t42 * t125 * t46 - 0.4e1 / 0.3e1 * t41 * t129 * t132;
    const double t137 = t123 * t136;
    const double t141 = -0.16e2 / 0.27e2 * t33 * t35 * t114 * t50 - 0.2e1 / 0.9e1 * t120 * t121 * t137;
    const double t146 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t55 - t111 - 0.3e1 / 0.8e1 * t5 * t27 * t141 );
    const double t147 = t61 * t95;
    const double t149 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t147 );
    const double t152 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t149 );
    const double t153 = t152 * t26;
    const double t157 = t68 * t107;
    const double t160 = t5 * t157 * t89 / 0.8e1;
    const double t162 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t153 * t89 - t160 );
    const double t166 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t96 );
    const double t169 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t166 );
    const double t170 = t169 * t26;
    const double t175 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t170 * t55 - t111 );
    const double t177 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t147 );
    const double t180 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t177 );
    const double t181 = t180 * t26;
    const double t185 = t71 * rho_b;
    const double t187 = 0.1e1 / t73 / t185;
    const double t192 = sigma_bb * t75;
    const double t193 = t83 * t83;
    const double t194 = 0.1e1 / t193;
    const double t196 = 0.1e1 / t72 / t71;
    const double t200 = sigma_bb * t187;
    const double t201 = t192 + 0.1e1;
    const double t202 = safe_math::sqrt( t201 );
    const double t203 = 0.1e1 / t202;
    const double t207 = -0.4e1 / 0.3e1 * t41 * t76 * t196 * t80 - 0.4e1 / 0.3e1 * t41 * t200 * t203;
    const double t208 = t194 * t207;
    const double t212 = -0.16e2 / 0.27e2 * t33 * t70 * t187 * t84 - 0.2e1 / 0.9e1 * t120 * t192 * t208;
    const double t217 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t89 - t160 - 0.3e1 / 0.8e1 * t5 * t69 * t212 );
    const double t220 = t34 * t40;
    const double t223 = 0.1e1 / t42;
    const double t230 = t41 * t223 * t44 * t46 / 0.2e1 + t41 * t40 * t132 / 0.2e1;
    const double t231 = t123 * t230;
    const double t235 = -0.2e1 / 0.9e1 * t120 * t121 * t231 + 0.2e1 / 0.9e1 * t33 * t220 * t50;
    const double t239 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t235 );
    const double t240 = t34 * t75;
    const double t243 = 0.1e1 / t76;
    const double t250 = t41 * t243 * t78 * t80 / 0.2e1 + t41 * t75 * t203 / 0.2e1;
    const double t251 = t194 * t250;
    const double t255 = -0.2e1 / 0.9e1 * t120 * t192 * t251 + 0.2e1 / 0.9e1 * t33 * t240 * t84;
    const double t259 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t255 );
    const double t262 = t23 * t23;
    const double t263 = 0.1e1 / t262;
    const double t264 = t98 * t98;
    const double t267 = t94 * t6;
    const double t268 = 0.1e1 / t267;
    const double t269 = t16 * t268;
    const double t272 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t95 + 0.2e1 * t269 );
    const double t276 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t263 * t264 + 0.4e1 / 0.3e1 * t23 * t272 );
    const double t277 = t276 * t26;
    const double t281 = t101 * t107;
    const double t283 = t5 * t281 * t55;
    const double t289 = 0.1e1 / t106 / t6;
    const double t290 = t25 * t289;
    const double t293 = t5 * t290 * t55 / 0.12e2;
    const double t295 = t5 * t108 * t141;
    const double t297 = t36 * t36;
    const double t299 = 0.1e1 / t38 / t297;
    const double t308 = 0.1e1 / t122 / t49;
    const double t309 = t136 * t136;
    const double t310 = t308 * t309;
    const double t315 = 0.1e1 / t37 / t112;
    const double t320 = sigma_aa * t299;
    const double t324 = sigma_aa * sigma_aa;
    const double t327 = 0.1e1 / t37 / t297 / t112;
    const double t330 = 0.1e1 / t131 / t130;
    const double t334 = 0.28e2 / 0.9e1 * t41 * t42 * t315 * t46 + 0.2e2 / 0.3e1 * t41 * t320 * t132 - 0.16e2 / 0.9e1 * t41 * t324 * t327 * t330;
    const double t335 = t123 * t334;
    const double t339 = 0.176e3 / 0.81e2 * t33 * t35 * t299 * t50 + 0.32e2 / 0.27e2 * t120 * t129 * t137 + 0.4e1 / 0.9e1 * t120 * t121 * t310 - 0.2e1 / 0.9e1 * t120 * t121 * t335;
    const double t344 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t277 * t55 - t283 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t102 * t141 + t293 - t295 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t27 * t339 );
    const double t345 = t66 * t66;
    const double t346 = 0.1e1 / t345;
    const double t347 = t149 * t149;
    const double t350 = t61 * t268;
    const double t353 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t95 + 0.2e1 * t350 );
    const double t357 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t346 * t347 + 0.4e1 / 0.3e1 * t66 * t353 );
    const double t358 = t357 * t26;
    const double t362 = t152 * t107;
    const double t364 = t5 * t362 * t89;
    const double t366 = t68 * t289;
    const double t369 = t5 * t366 * t89 / 0.12e2;
    const double t371 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t358 * t89 - t364 / 0.4e1 + t369 );
    const double t387 = t169 * t107;
    const double t389 = t5 * t387 * t55;
    const double t411 = t180 * t107;
    const double t413 = t5 * t411 * t89;
    const double t420 = t5 * t157 * t212;
    const double t428 = t166 * t166;
    const double t433 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t95 + 0.2e1 * t269 );
    const double t437 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t263 * t428 + 0.4e1 / 0.3e1 * t23 * t433 );
    const double t438 = t437 * t26;
    const double t444 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t438 * t55 - t389 / 0.4e1 + t293 );
    const double t445 = t177 * t177;
    const double t450 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t95 + 0.2e1 * t350 );
    const double t454 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.9e1 * t346 * t445 + 0.4e1 / 0.3e1 * t66 * t450 );
    const double t455 = t454 * t26;
    const double t464 = t71 * t71;
    const double t466 = 0.1e1 / t73 / t464;
    const double t475 = 0.1e1 / t193 / t83;
    const double t476 = t207 * t207;
    const double t477 = t475 * t476;
    const double t482 = 0.1e1 / t72 / t185;
    const double t487 = sigma_bb * t466;
    const double t491 = sigma_bb * sigma_bb;
    const double t494 = 0.1e1 / t72 / t464 / t185;
    const double t497 = 0.1e1 / t202 / t201;
    const double t501 = 0.28e2 / 0.9e1 * t41 * t76 * t482 * t80 + 0.2e2 / 0.3e1 * t41 * t487 * t203 - 0.16e2 / 0.9e1 * t41 * t491 * t494 * t497;
    const double t502 = t194 * t501;
    const double t506 = 0.176e3 / 0.81e2 * t33 * t70 * t466 * t84 + 0.32e2 / 0.27e2 * t120 * t200 * t208 + 0.4e1 / 0.9e1 * t120 * t192 * t477 - 0.2e1 / 0.9e1 * t120 * t192 * t502;
    const double t511 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t455 * t89 - t413 / 0.4e1 - 0.3e1 / 0.4e1 * t5 * t181 * t212 + t369 - t420 / 0.4e1 - 0.3e1 / 0.8e1 * t5 * t69 * t506 );
    const double t519 = t5 * t108 * t235 / 0.8e1;
    const double t520 = t34 * t114;
    const double t530 = t308 * t230;
    const double t531 = t530 * t136;
    const double t542 = t297 * t36;
    const double t544 = 0.1e1 / t37 / t542;
    const double t545 = t544 * t330;
    const double t549 = -0.2e1 / 0.3e1 * t41 * t223 * t125 * t46 - 0.2e1 * t41 * t114 * t132 + 0.2e1 / 0.3e1 * t41 * t545 * sigma_aa;
    const double t550 = t123 * t549;
    const double t554 = -0.16e2 / 0.27e2 * t33 * t520 * t50 - 0.2e1 / 0.9e1 * t33 * t220 * t137 + 0.16e2 / 0.27e2 * t120 * t129 * t231 + 0.4e1 / 0.9e1 * t120 * t121 * t531 - 0.2e1 / 0.9e1 * t120 * t121 * t550;
    const double t559 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t235 - t519 - 0.3e1 / 0.8e1 * t5 * t27 * t554 );
    const double t566 = t5 * t157 * t255 / 0.8e1;
    const double t579 = t34 * t187;
    const double t589 = t475 * t250;
    const double t590 = t589 * t207;
    const double t601 = t464 * t71;
    const double t603 = 0.1e1 / t72 / t601;
    const double t604 = t603 * t497;
    const double t608 = -0.2e1 / 0.3e1 * t41 * t243 * t196 * t80 - 0.2e1 * t41 * t187 * t203 + 0.2e1 / 0.3e1 * t41 * t604 * sigma_bb;
    const double t609 = t194 * t608;
    const double t613 = -0.16e2 / 0.27e2 * t33 * t579 * t84 - 0.2e1 / 0.9e1 * t33 * t240 * t208 + 0.16e2 / 0.27e2 * t120 * t200 * t251 + 0.4e1 / 0.9e1 * t120 * t192 * t590 - 0.2e1 / 0.9e1 * t120 * t192 * t609;
    const double t618 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t255 - t566 - 0.3e1 / 0.8e1 * t5 * t69 * t613 );
    const double t623 = t230 * t230;
    const double t624 = t308 * t623;
    const double t629 = 0.1e1 / t42 / sigma_aa;
    const double t633 = 0.1e1 / sigma_aa;
    const double t637 = t297 * rho_a;
    const double t639 = 0.1e1 / t37 / t637;
    const double t643 = t41 * t633 * t40 * t132 / 0.4e1 - t41 * t629 * t44 * t46 / 0.4e1 - t41 * t639 * t330 / 0.4e1;
    const double t644 = t123 * t643;
    const double t648 = -0.4e1 / 0.9e1 * t33 * t220 * t231 + 0.4e1 / 0.9e1 * t120 * t121 * t624 - 0.2e1 / 0.9e1 * t120 * t121 * t644;
    const double t652 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t648 );
    const double t656 = t250 * t250;
    const double t657 = t475 * t656;
    const double t662 = 0.1e1 / t76 / sigma_bb;
    const double t666 = 0.1e1 / sigma_bb;
    const double t670 = t464 * rho_b;
    const double t672 = 0.1e1 / t72 / t670;
    const double t676 = t41 * t666 * t75 * t203 / 0.4e1 - t41 * t662 * t78 * t80 / 0.4e1 - t41 * t672 * t497 / 0.4e1;
    const double t677 = t194 * t676;
    const double t681 = -0.4e1 / 0.9e1 * t33 * t240 * t251 + 0.4e1 / 0.9e1 * t120 * t192 * t657 - 0.2e1 / 0.9e1 * t120 * t192 * t677;
    const double t685 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t681 );


    vrho_a = t59 + t93 + t6 * ( t146 + t162 );
    vrho_b = t59 + t93 + t6 * ( t175 + t217 );
    vsigma_aa = t6 * t239;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t259;
    v2rho2_aa = 0.2e1 * t146 + 0.2e1 * t162 + t6 * ( t344 + t371 );
    v2rho2_bb = 0.2e1 * t175 + 0.2e1 * t217 + t6 * ( t444 + t511 );
    v2rhosigma_a_aa = t6 * t559 + t239;
    v2rhosigma_b_bb = t6 * t618 + t259;
    v2sigma2_aa_aa = t6 * t652;
    v2sigma2_bb_bb = t6 * t685;
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

struct BuiltinB88 : detail::BuiltinKernelImpl< BuiltinB88 > {

  BuiltinB88( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinB88 >(p) { }
  
  virtual ~BuiltinB88() = default;

};



} // namespace ExchCXX
