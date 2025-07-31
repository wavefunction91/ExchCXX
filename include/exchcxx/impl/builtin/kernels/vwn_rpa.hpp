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
struct kernel_traits< BuiltinVWN_RPA > :
  public lda_screening_interface< BuiltinVWN_RPA > {

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




  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t21 = t12 + 0.13072e2;
    const double t24 = safe_math::atan( 0.44899888641287296627e-1 / t21 );
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t45 = ( 0.310907e-1 * t19 + 0.20521972937837502661e2 * t24 + 0.44313737677495382697e-2 * t30 ) * t44;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t52 = safe_math::log( t4 * t9 * t48 / 0.4e1 );
    const double t54 = t12 + 0.201231e2;
    const double t57 = safe_math::atan( 0.11716852777089929792e1 / t54 );
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t62 = safe_math::log( t60 * t48 );
    const double t66 = ( 0.1554535e-1 * t52 + 0.61881802979060631482e0 * t57 + 0.26673100072733151594e-2 * t62 ) * t38 * t42;


    eps = t45 + t66;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t73 = t4 * t6;
    constexpr double t81 = t3 * t6;
    constexpr double t90 = t1 * t1;
    constexpr double t92 = 0.1e1 / t3;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t21 = t12 + 0.13072e2;
    const double t24 = safe_math::atan( 0.44899888641287296627e-1 / t21 );
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t45 = ( 0.310907e-1 * t19 + 0.20521972937837502661e2 * t24 + 0.44313737677495382697e-2 * t30 ) * t44;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t52 = safe_math::log( t4 * t9 * t48 / 0.4e1 );
    const double t54 = t12 + 0.201231e2;
    const double t57 = safe_math::atan( 0.11716852777089929792e1 / t54 );
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t62 = safe_math::log( t60 * t48 );
    const double t66 = ( 0.1554535e-1 * t52 + 0.61881802979060631482e0 * t57 + 0.26673100072733151594e-2 * t62 ) * t38 * t42;
    const double t68 = 0.1e1 / t7 / rho;
    const double t69 = t6 * t68;
    const double t74 = t14 * t14;
    const double t75 = 0.1e1 / t74;
    const double t76 = t8 * t75;
    const double t77 = t4 * t69;
    const double t78 = t77 / 0.12e2;
    const double t79 = 0.1e1 / t12;
    const double t80 = t79 * t1;
    const double t83 = t80 * t81 * t68;
    const double t85 = -t78 - 0.10893333333333333333e1 * t83;
    const double t93 = ( -t4 * t69 * t15 / 0.12e2 - t73 * t76 * t85 / 0.4e1 ) * t90 * t92;
    const double t94 = t5 * t7;
    const double t95 = t94 * t14;
    const double t98 = t21 * t21;
    const double t99 = 0.1e1 / t98;
    const double t101 = t99 * t79 * t1;
    const double t103 = 0.2016e-2 * t99 + 0.1e1;
    const double t104 = 0.1e1 / t103;
    const double t109 = t27 * t15;
    const double t110 = t109 * t79;
    const double t113 = t28 * t75;
    const double t115 = -t110 * t77 / 0.6e1 - t113 * t85;
    const double t116 = 0.1e1 / t28;
    const double t117 = t115 * t116;
    const double t121 = ( 0.10363566666666666667e-1 * t93 * t95 + 0.15357238326806922974e0 * t101 * t81 * t68 * t104 + 0.44313737677495382697e-2 * t117 * t14 ) * t44;
    const double t125 = t47 * t47;
    const double t126 = 0.1e1 / t125;
    const double t127 = t8 * t126;
    const double t129 = -t78 - 0.1676925e1 * t83;
    const double t135 = ( -t4 * t69 * t48 / 0.12e2 - t73 * t127 * t129 / 0.4e1 ) * t90 * t92;
    const double t136 = t94 * t47;
    const double t139 = t54 * t54;
    const double t140 = 0.1e1 / t139;
    const double t142 = t140 * t79 * t1;
    const double t144 = 0.137284639e1 * t140 + 0.1e1;
    const double t145 = 0.1e1 / t144;
    const double t150 = t59 * t48;
    const double t151 = t150 * t79;
    const double t154 = t60 * t126;
    const double t156 = -t151 * t77 / 0.6e1 - t154 * t129;
    const double t157 = 0.1e1 / t60;
    const double t158 = t156 * t157;
    const double t163 = ( 0.51817833333333333333e-2 * t135 * t136 + 0.12084332918108974175e0 * t142 * t81 * t68 * t145 + 0.26673100072733151594e-2 * t158 * t47 ) * t38 * t42;


    eps = t45 + t66;
    vrho = t45 + t66 + rho * ( t121 + t163 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double& v2rho2 ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t73 = t4 * t6;
    constexpr double t81 = t3 * t6;
    constexpr double t90 = t1 * t1;
    constexpr double t92 = 0.1e1 / t3;
    constexpr double t191 = t3 * t3;
    constexpr double t192 = t191 * t5;
    constexpr double t254 = t90 * t191;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t21 = t12 + 0.13072e2;
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t54 = t12 + 0.201231e2;
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t68 = 0.1e1 / t7 / rho;
    const double t69 = t6 * t68;
    const double t74 = t14 * t14;
    const double t75 = 0.1e1 / t74;
    const double t76 = t8 * t75;
    const double t77 = t4 * t69;
    const double t78 = t77 / 0.12e2;
    const double t79 = 0.1e1 / t12;
    const double t80 = t79 * t1;
    const double t83 = t80 * t81 * t68;
    const double t85 = -t78 - 0.10893333333333333333e1 * t83;
    const double t93 = ( -t4 * t69 * t15 / 0.12e2 - t73 * t76 * t85 / 0.4e1 ) * t90 * t92;
    const double t94 = t5 * t7;
    const double t95 = t94 * t14;
    const double t98 = t21 * t21;
    const double t99 = 0.1e1 / t98;
    const double t101 = t99 * t79 * t1;
    const double t103 = 0.2016e-2 * t99 + 0.1e1;
    const double t104 = 0.1e1 / t103;
    const double t109 = t27 * t15;
    const double t110 = t109 * t79;
    const double t113 = t28 * t75;
    const double t115 = -t110 * t77 / 0.6e1 - t113 * t85;
    const double t116 = 0.1e1 / t28;
    const double t117 = t115 * t116;
    const double t121 = ( 0.10363566666666666667e-1 * t93 * t95 + 0.15357238326806922974e0 * t101 * t81 * t68 * t104 + 0.44313737677495382697e-2 * t117 * t14 ) * t44;
    const double t125 = t47 * t47;
    const double t126 = 0.1e1 / t125;
    const double t127 = t8 * t126;
    const double t129 = -t78 - 0.1676925e1 * t83;
    const double t135 = ( -t4 * t69 * t48 / 0.12e2 - t73 * t127 * t129 / 0.4e1 ) * t90 * t92;
    const double t136 = t94 * t47;
    const double t139 = t54 * t54;
    const double t140 = 0.1e1 / t139;
    const double t142 = t140 * t79 * t1;
    const double t144 = 0.137284639e1 * t140 + 0.1e1;
    const double t145 = 0.1e1 / t144;
    const double t150 = t59 * t48;
    const double t151 = t150 * t79;
    const double t154 = t60 * t126;
    const double t156 = -t151 * t77 / 0.6e1 - t154 * t129;
    const double t157 = 0.1e1 / t60;
    const double t158 = t156 * t157;
    const double t163 = ( 0.51817833333333333333e-2 * t135 * t136 + 0.12084332918108974175e0 * t142 * t81 * t68 * t145 + 0.26673100072733151594e-2 * t158 * t47 ) * t38 * t42;
    const double t168 = rho * rho;
    const double t170 = 0.1e1 / t7 / t168;
    const double t171 = t6 * t170;
    const double t173 = t4 * t171 * t15;
    const double t175 = t68 * t75;
    const double t180 = 0.1e1 / t74 / t14;
    const double t181 = t8 * t180;
    const double t182 = t85 * t85;
    const double t186 = t4 * t171;
    const double t187 = t186 / 0.9e1;
    const double t189 = 0.1e1 / t12 / t10;
    const double t190 = t189 * t90;
    const double t193 = t7 * t7;
    const double t195 = 0.1e1 / t193 / t168;
    const double t197 = t190 * t192 * t195;
    const double t200 = t80 * t81 * t170;
    const double t202 = t187 - 0.7262222222222222222e0 * t197 + 0.14524444444444444444e1 * t200;
    const double t208 = ( t173 / 0.9e1 + t73 * t175 * t85 / 0.6e1 + t73 * t181 * t182 / 0.2e1 - t73 * t76 * t202 / 0.4e1 ) * t90 * t92;
    const double t212 = t5 / t193;
    const double t213 = t212 * t14;
    const double t216 = t94 * t85;
    const double t219 = t98 * t21;
    const double t221 = 0.1e1 / t219 * t1;
    const double t222 = t221 * t3;
    const double t227 = t99 * t189 * t90;
    const double t236 = t98 * t98;
    const double t238 = 0.1e1 / t236 / t21;
    const double t239 = t238 * t1;
    const double t240 = t239 * t3;
    const double t241 = t103 * t103;
    const double t242 = 0.1e1 / t241;
    const double t247 = t27 * t75;
    const double t248 = t247 * t80;
    const double t249 = t68 * t85;
    const double t253 = t109 * t189;
    const double t255 = t5 * t195;
    const double t256 = t254 * t255;
    const double t261 = t28 * t180;
    const double t265 = t173 / 0.72e2 + t248 * t81 * t249 / 0.3e1 - t253 * t256 / 0.9e1 + 0.2e1 / 0.9e1 * t110 * t186 + 0.2e1 * t261 * t182 - t113 * t202;
    const double t266 = t265 * t116;
    const double t270 = 0.1e1 / t28 / t27;
    const double t271 = t115 * t270;
    const double t272 = t14 * t79;
    const double t273 = t271 * t272;
    const double t279 = ( 0.10363566666666666667e-1 * t208 * t95 + 0.34545222222222222223e-2 * t93 * t213 + 0.10363566666666666667e-1 * t93 * t216 + 0.51190794422689743247e-1 * t222 * t171 * t104 + 0.10238158884537948649e0 * t227 * t192 * t195 * t104 - 0.20476317769075897299e0 * t101 * t81 * t170 * t104 - 0.10320064155614252239e-3 * t240 * t171 * t242 + 0.44313737677495382697e-2 * t266 * t14 + 0.73856229462492304495e-3 * t273 * t77 + 0.44313737677495382697e-2 * t117 * t85 ) * t44;
    const double t281 = t4 * t171 * t48;
    const double t283 = t68 * t126;
    const double t288 = 0.1e1 / t125 / t47;
    const double t289 = t8 * t288;
    const double t290 = t129 * t129;
    const double t296 = t187 - 0.111795e1 * t197 + 0.22359e1 * t200;
    const double t302 = ( t281 / 0.9e1 + t73 * t283 * t129 / 0.6e1 + t73 * t289 * t290 / 0.2e1 - t73 * t127 * t296 / 0.4e1 ) * t90 * t92;
    const double t305 = t212 * t47;
    const double t308 = t94 * t129;
    const double t311 = t139 * t54;
    const double t313 = 0.1e1 / t311 * t1;
    const double t314 = t313 * t3;
    const double t319 = t140 * t189 * t90;
    const double t328 = t139 * t139;
    const double t330 = 0.1e1 / t328 / t54;
    const double t331 = t330 * t1;
    const double t332 = t331 * t3;
    const double t333 = t144 * t144;
    const double t334 = 0.1e1 / t333;
    const double t339 = t59 * t126;
    const double t340 = t339 * t80;
    const double t341 = t68 * t129;
    const double t345 = t150 * t189;
    const double t350 = t60 * t288;
    const double t354 = t281 / 0.72e2 + t340 * t81 * t341 / 0.3e1 - t345 * t256 / 0.9e1 + 0.2e1 / 0.9e1 * t151 * t186 + 0.2e1 * t350 * t290 - t154 * t296;
    const double t355 = t354 * t157;
    const double t359 = 0.1e1 / t60 / t59;
    const double t360 = t156 * t359;
    const double t361 = t47 * t79;
    const double t362 = t360 * t361;
    const double t369 = ( 0.51817833333333333333e-2 * t302 * t136 + 0.17272611111111111111e-2 * t135 * t305 + 0.51817833333333333333e-2 * t135 * t308 + 0.40281109727029913917e-1 * t314 * t171 * t145 + 0.80562219454059827833e-1 * t319 * t192 * t195 * t145 - 0.16112443890811965567e0 * t142 * t81 * t170 * t145 - 0.55299776073946902743e-1 * t332 * t171 * t334 + 0.26673100072733151594e-2 * t355 * t47 + 0.4445516678788858599e-3 * t362 * t77 + 0.26673100072733151594e-2 * t158 * t129 ) * t38 * t42;


    v2rho2 = 0.2e1 * t121 + 0.2e1 * t163 + rho * ( t279 + t369 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double& vrho, double& v2rho2 ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t73 = t4 * t6;
    constexpr double t81 = t3 * t6;
    constexpr double t90 = t1 * t1;
    constexpr double t92 = 0.1e1 / t3;
    constexpr double t191 = t3 * t3;
    constexpr double t192 = t191 * t5;
    constexpr double t254 = t90 * t191;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t21 = t12 + 0.13072e2;
    const double t24 = safe_math::atan( 0.44899888641287296627e-1 / t21 );
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t45 = ( 0.310907e-1 * t19 + 0.20521972937837502661e2 * t24 + 0.44313737677495382697e-2 * t30 ) * t44;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t52 = safe_math::log( t4 * t9 * t48 / 0.4e1 );
    const double t54 = t12 + 0.201231e2;
    const double t57 = safe_math::atan( 0.11716852777089929792e1 / t54 );
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t62 = safe_math::log( t60 * t48 );
    const double t66 = ( 0.1554535e-1 * t52 + 0.61881802979060631482e0 * t57 + 0.26673100072733151594e-2 * t62 ) * t38 * t42;
    const double t68 = 0.1e1 / t7 / rho;
    const double t69 = t6 * t68;
    const double t74 = t14 * t14;
    const double t75 = 0.1e1 / t74;
    const double t76 = t8 * t75;
    const double t77 = t4 * t69;
    const double t78 = t77 / 0.12e2;
    const double t79 = 0.1e1 / t12;
    const double t80 = t79 * t1;
    const double t83 = t80 * t81 * t68;
    const double t85 = -t78 - 0.10893333333333333333e1 * t83;
    const double t93 = ( -t4 * t69 * t15 / 0.12e2 - t73 * t76 * t85 / 0.4e1 ) * t90 * t92;
    const double t94 = t5 * t7;
    const double t95 = t94 * t14;
    const double t98 = t21 * t21;
    const double t99 = 0.1e1 / t98;
    const double t101 = t99 * t79 * t1;
    const double t103 = 0.2016e-2 * t99 + 0.1e1;
    const double t104 = 0.1e1 / t103;
    const double t109 = t27 * t15;
    const double t110 = t109 * t79;
    const double t113 = t28 * t75;
    const double t115 = -t110 * t77 / 0.6e1 - t113 * t85;
    const double t116 = 0.1e1 / t28;
    const double t117 = t115 * t116;
    const double t121 = ( 0.10363566666666666667e-1 * t93 * t95 + 0.15357238326806922974e0 * t101 * t81 * t68 * t104 + 0.44313737677495382697e-2 * t117 * t14 ) * t44;
    const double t125 = t47 * t47;
    const double t126 = 0.1e1 / t125;
    const double t127 = t8 * t126;
    const double t129 = -t78 - 0.1676925e1 * t83;
    const double t135 = ( -t4 * t69 * t48 / 0.12e2 - t73 * t127 * t129 / 0.4e1 ) * t90 * t92;
    const double t136 = t94 * t47;
    const double t139 = t54 * t54;
    const double t140 = 0.1e1 / t139;
    const double t142 = t140 * t79 * t1;
    const double t144 = 0.137284639e1 * t140 + 0.1e1;
    const double t145 = 0.1e1 / t144;
    const double t150 = t59 * t48;
    const double t151 = t150 * t79;
    const double t154 = t60 * t126;
    const double t156 = -t151 * t77 / 0.6e1 - t154 * t129;
    const double t157 = 0.1e1 / t60;
    const double t158 = t156 * t157;
    const double t163 = ( 0.51817833333333333333e-2 * t135 * t136 + 0.12084332918108974175e0 * t142 * t81 * t68 * t145 + 0.26673100072733151594e-2 * t158 * t47 ) * t38 * t42;
    const double t168 = rho * rho;
    const double t170 = 0.1e1 / t7 / t168;
    const double t171 = t6 * t170;
    const double t173 = t4 * t171 * t15;
    const double t175 = t68 * t75;
    const double t180 = 0.1e1 / t74 / t14;
    const double t181 = t8 * t180;
    const double t182 = t85 * t85;
    const double t186 = t4 * t171;
    const double t187 = t186 / 0.9e1;
    const double t189 = 0.1e1 / t12 / t10;
    const double t190 = t189 * t90;
    const double t193 = t7 * t7;
    const double t195 = 0.1e1 / t193 / t168;
    const double t197 = t190 * t192 * t195;
    const double t200 = t80 * t81 * t170;
    const double t202 = t187 - 0.7262222222222222222e0 * t197 + 0.14524444444444444444e1 * t200;
    const double t208 = ( t173 / 0.9e1 + t73 * t175 * t85 / 0.6e1 + t73 * t181 * t182 / 0.2e1 - t73 * t76 * t202 / 0.4e1 ) * t90 * t92;
    const double t212 = t5 / t193;
    const double t213 = t212 * t14;
    const double t216 = t94 * t85;
    const double t219 = t98 * t21;
    const double t221 = 0.1e1 / t219 * t1;
    const double t222 = t221 * t3;
    const double t227 = t99 * t189 * t90;
    const double t236 = t98 * t98;
    const double t238 = 0.1e1 / t236 / t21;
    const double t239 = t238 * t1;
    const double t240 = t239 * t3;
    const double t241 = t103 * t103;
    const double t242 = 0.1e1 / t241;
    const double t247 = t27 * t75;
    const double t248 = t247 * t80;
    const double t249 = t68 * t85;
    const double t253 = t109 * t189;
    const double t255 = t5 * t195;
    const double t256 = t254 * t255;
    const double t261 = t28 * t180;
    const double t265 = t173 / 0.72e2 + t248 * t81 * t249 / 0.3e1 - t253 * t256 / 0.9e1 + 0.2e1 / 0.9e1 * t110 * t186 + 0.2e1 * t261 * t182 - t113 * t202;
    const double t266 = t265 * t116;
    const double t270 = 0.1e1 / t28 / t27;
    const double t271 = t115 * t270;
    const double t272 = t14 * t79;
    const double t273 = t271 * t272;
    const double t279 = ( 0.10363566666666666667e-1 * t208 * t95 + 0.34545222222222222223e-2 * t93 * t213 + 0.10363566666666666667e-1 * t93 * t216 + 0.51190794422689743247e-1 * t222 * t171 * t104 + 0.10238158884537948649e0 * t227 * t192 * t195 * t104 - 0.20476317769075897299e0 * t101 * t81 * t170 * t104 - 0.10320064155614252239e-3 * t240 * t171 * t242 + 0.44313737677495382697e-2 * t266 * t14 + 0.73856229462492304495e-3 * t273 * t77 + 0.44313737677495382697e-2 * t117 * t85 ) * t44;
    const double t281 = t4 * t171 * t48;
    const double t283 = t68 * t126;
    const double t288 = 0.1e1 / t125 / t47;
    const double t289 = t8 * t288;
    const double t290 = t129 * t129;
    const double t296 = t187 - 0.111795e1 * t197 + 0.22359e1 * t200;
    const double t302 = ( t281 / 0.9e1 + t73 * t283 * t129 / 0.6e1 + t73 * t289 * t290 / 0.2e1 - t73 * t127 * t296 / 0.4e1 ) * t90 * t92;
    const double t305 = t212 * t47;
    const double t308 = t94 * t129;
    const double t311 = t139 * t54;
    const double t313 = 0.1e1 / t311 * t1;
    const double t314 = t313 * t3;
    const double t319 = t140 * t189 * t90;
    const double t328 = t139 * t139;
    const double t330 = 0.1e1 / t328 / t54;
    const double t331 = t330 * t1;
    const double t332 = t331 * t3;
    const double t333 = t144 * t144;
    const double t334 = 0.1e1 / t333;
    const double t339 = t59 * t126;
    const double t340 = t339 * t80;
    const double t341 = t68 * t129;
    const double t345 = t150 * t189;
    const double t350 = t60 * t288;
    const double t354 = t281 / 0.72e2 + t340 * t81 * t341 / 0.3e1 - t345 * t256 / 0.9e1 + 0.2e1 / 0.9e1 * t151 * t186 + 0.2e1 * t350 * t290 - t154 * t296;
    const double t355 = t354 * t157;
    const double t359 = 0.1e1 / t60 / t59;
    const double t360 = t156 * t359;
    const double t361 = t47 * t79;
    const double t362 = t360 * t361;
    const double t369 = ( 0.51817833333333333333e-2 * t302 * t136 + 0.17272611111111111111e-2 * t135 * t305 + 0.51817833333333333333e-2 * t135 * t308 + 0.40281109727029913917e-1 * t314 * t171 * t145 + 0.80562219454059827833e-1 * t319 * t192 * t195 * t145 - 0.16112443890811965567e0 * t142 * t81 * t170 * t145 - 0.55299776073946902743e-1 * t332 * t171 * t334 + 0.26673100072733151594e-2 * t355 * t47 + 0.4445516678788858599e-3 * t362 * t77 + 0.26673100072733151594e-2 * t158 * t129 ) * t38 * t42;


    vrho = t45 + t66 + rho * ( t121 + t163 );
    v2rho2 = 0.2e1 * t121 + 0.2e1 * t163 + rho * ( t279 + t369 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t56 = t33 * t55;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t77 = t75 * t49 * t53;


    eps = t56 + t77;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t84 = t4 * t6;
    constexpr double t92 = t3 * t6;
    constexpr double t101 = t1 * t1;
    constexpr double t103 = 0.1e1 / t3;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t56 = t33 * t55;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t77 = t75 * t49 * t53;
    const double t79 = 0.1e1 / t8 / t7;
    const double t80 = t6 * t79;
    const double t85 = t15 * t15;
    const double t86 = 0.1e1 / t85;
    const double t87 = t9 * t86;
    const double t88 = t4 * t80;
    const double t89 = t88 / 0.12e2;
    const double t90 = 0.1e1 / t13;
    const double t91 = t90 * t1;
    const double t94 = t91 * t92 * t79;
    const double t96 = -t89 - 0.10893333333333333333e1 * t94;
    const double t104 = ( -t4 * t80 * t16 / 0.12e2 - t84 * t87 * t96 / 0.4e1 ) * t101 * t103;
    const double t105 = t5 * t8;
    const double t106 = t105 * t15;
    const double t109 = t22 * t22;
    const double t110 = 0.1e1 / t109;
    const double t112 = t110 * t90 * t1;
    const double t114 = 0.2016e-2 * t110 + 0.1e1;
    const double t115 = 0.1e1 / t114;
    const double t120 = t28 * t16;
    const double t121 = t120 * t90;
    const double t124 = t29 * t86;
    const double t126 = -t121 * t88 / 0.6e1 - t124 * t96;
    const double t127 = 0.1e1 / t29;
    const double t128 = t126 * t127;
    const double t131 = 0.10363566666666666667e-1 * t104 * t106 + 0.15357238326806922974e0 * t112 * t92 * t79 * t115 + 0.44313737677495382697e-2 * t128 * t15;
    const double t132 = t131 * t55;
    const double t133 = t7 * t7;
    const double t134 = 0.1e1 / t133;
    const double t135 = t34 * t134;
    const double t136 = t35 - t135;
    const double t139 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t136 );
    const double t140 = -t136;
    const double t143 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t140 );
    const double t144 = t139 + t143;
    const double t146 = t33 * t144 * t53;
    const double t150 = t58 * t58;
    const double t151 = 0.1e1 / t150;
    const double t152 = t9 * t151;
    const double t154 = -t89 - 0.1676925e1 * t94;
    const double t160 = ( -t4 * t80 * t59 / 0.12e2 - t84 * t152 * t154 / 0.4e1 ) * t101 * t103;
    const double t161 = t105 * t58;
    const double t164 = t65 * t65;
    const double t165 = 0.1e1 / t164;
    const double t167 = t165 * t90 * t1;
    const double t169 = 0.137284639e1 * t165 + 0.1e1;
    const double t170 = 0.1e1 / t169;
    const double t175 = t70 * t59;
    const double t176 = t175 * t90;
    const double t179 = t71 * t151;
    const double t181 = -t176 * t88 / 0.6e1 - t179 * t154;
    const double t182 = 0.1e1 / t71;
    const double t183 = t181 * t182;
    const double t186 = 0.51817833333333333333e-2 * t160 * t161 + 0.12084332918108974175e0 * t167 * t92 * t79 * t170 + 0.26673100072733151594e-2 * t183 * t58;
    const double t188 = t186 * t49 * t53;
    const double t190 = t75 * t144 * t53;
    const double t193 = -t35 - t135;
    const double t196 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t193 );
    const double t197 = -t193;
    const double t200 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t197 );
    const double t201 = t196 + t200;
    const double t203 = t33 * t201 * t53;
    const double t205 = t75 * t201 * t53;


    eps = t56 + t77;
    vrho_a = t56 + t77 + t7 * ( t132 - t146 + t188 + t190 );
    vrho_b = t56 + t77 + t7 * ( t132 - t203 + t188 + t205 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t84 = t4 * t6;
    constexpr double t92 = t3 * t6;
    constexpr double t101 = t1 * t1;
    constexpr double t103 = 0.1e1 / t3;
    constexpr double t234 = t3 * t3;
    constexpr double t235 = t234 * t5;
    constexpr double t297 = t101 * t234;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t79 = 0.1e1 / t8 / t7;
    const double t80 = t6 * t79;
    const double t85 = t15 * t15;
    const double t86 = 0.1e1 / t85;
    const double t87 = t9 * t86;
    const double t88 = t4 * t80;
    const double t89 = t88 / 0.12e2;
    const double t90 = 0.1e1 / t13;
    const double t91 = t90 * t1;
    const double t94 = t91 * t92 * t79;
    const double t96 = -t89 - 0.10893333333333333333e1 * t94;
    const double t104 = ( -t4 * t80 * t16 / 0.12e2 - t84 * t87 * t96 / 0.4e1 ) * t101 * t103;
    const double t105 = t5 * t8;
    const double t106 = t105 * t15;
    const double t109 = t22 * t22;
    const double t110 = 0.1e1 / t109;
    const double t112 = t110 * t90 * t1;
    const double t114 = 0.2016e-2 * t110 + 0.1e1;
    const double t115 = 0.1e1 / t114;
    const double t120 = t28 * t16;
    const double t121 = t120 * t90;
    const double t124 = t29 * t86;
    const double t126 = -t121 * t88 / 0.6e1 - t124 * t96;
    const double t127 = 0.1e1 / t29;
    const double t128 = t126 * t127;
    const double t131 = 0.10363566666666666667e-1 * t104 * t106 + 0.15357238326806922974e0 * t112 * t92 * t79 * t115 + 0.44313737677495382697e-2 * t128 * t15;
    const double t132 = t131 * t55;
    const double t133 = t7 * t7;
    const double t134 = 0.1e1 / t133;
    const double t135 = t34 * t134;
    const double t136 = t35 - t135;
    const double t139 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t136 );
    const double t140 = -t136;
    const double t143 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t140 );
    const double t144 = t139 + t143;
    const double t146 = t33 * t144 * t53;
    const double t150 = t58 * t58;
    const double t151 = 0.1e1 / t150;
    const double t152 = t9 * t151;
    const double t154 = -t89 - 0.1676925e1 * t94;
    const double t160 = ( -t4 * t80 * t59 / 0.12e2 - t84 * t152 * t154 / 0.4e1 ) * t101 * t103;
    const double t161 = t105 * t58;
    const double t164 = t65 * t65;
    const double t165 = 0.1e1 / t164;
    const double t167 = t165 * t90 * t1;
    const double t169 = 0.137284639e1 * t165 + 0.1e1;
    const double t170 = 0.1e1 / t169;
    const double t175 = t70 * t59;
    const double t176 = t175 * t90;
    const double t179 = t71 * t151;
    const double t181 = -t176 * t88 / 0.6e1 - t179 * t154;
    const double t182 = 0.1e1 / t71;
    const double t183 = t181 * t182;
    const double t186 = 0.51817833333333333333e-2 * t160 * t161 + 0.12084332918108974175e0 * t167 * t92 * t79 * t170 + 0.26673100072733151594e-2 * t183 * t58;
    const double t188 = t186 * t49 * t53;
    const double t190 = t75 * t144 * t53;
    const double t193 = -t35 - t135;
    const double t196 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t193 );
    const double t197 = -t193;
    const double t200 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t197 );
    const double t201 = t196 + t200;
    const double t203 = t33 * t201 * t53;
    const double t205 = t75 * t201 * t53;
    const double t208 = 0.2e1 * t132;
    const double t210 = 0.2e1 * t188;
    const double t213 = 0.1e1 / t8 / t133;
    const double t214 = t6 * t213;
    const double t216 = t4 * t214 * t16;
    const double t218 = t79 * t86;
    const double t223 = 0.1e1 / t85 / t15;
    const double t224 = t9 * t223;
    const double t225 = t96 * t96;
    const double t229 = t4 * t214;
    const double t230 = t229 / 0.9e1;
    const double t232 = 0.1e1 / t13 / t11;
    const double t233 = t232 * t101;
    const double t236 = t8 * t8;
    const double t238 = 0.1e1 / t236 / t133;
    const double t240 = t233 * t235 * t238;
    const double t243 = t91 * t92 * t213;
    const double t245 = t230 - 0.7262222222222222222e0 * t240 + 0.14524444444444444444e1 * t243;
    const double t251 = ( t216 / 0.9e1 + t84 * t218 * t96 / 0.6e1 + t84 * t224 * t225 / 0.2e1 - t84 * t87 * t245 / 0.4e1 ) * t101 * t103;
    const double t255 = t5 / t236;
    const double t256 = t255 * t15;
    const double t259 = t105 * t96;
    const double t262 = t109 * t22;
    const double t263 = 0.1e1 / t262;
    const double t264 = t263 * t1;
    const double t265 = t264 * t3;
    const double t270 = t110 * t232 * t101;
    const double t279 = t109 * t109;
    const double t281 = 0.1e1 / t279 / t22;
    const double t282 = t281 * t1;
    const double t283 = t282 * t3;
    const double t284 = t114 * t114;
    const double t285 = 0.1e1 / t284;
    const double t290 = t28 * t86;
    const double t291 = t290 * t91;
    const double t292 = t79 * t96;
    const double t296 = t120 * t232;
    const double t298 = t5 * t238;
    const double t299 = t297 * t298;
    const double t304 = t29 * t223;
    const double t308 = t216 / 0.72e2 + t291 * t92 * t292 / 0.3e1 - t296 * t299 / 0.9e1 + 0.2e1 / 0.9e1 * t121 * t229 + 0.2e1 * t304 * t225 - t124 * t245;
    const double t309 = t308 * t127;
    const double t313 = 0.1e1 / t29 / t28;
    const double t314 = t126 * t313;
    const double t315 = t15 * t90;
    const double t316 = t314 * t315;
    const double t321 = 0.10363566666666666667e-1 * t251 * t106 + 0.34545222222222222223e-2 * t104 * t256 + 0.10363566666666666667e-1 * t104 * t259 + 0.51190794422689743247e-1 * t265 * t214 * t115 + 0.10238158884537948649e0 * t270 * t235 * t238 * t115 - 0.20476317769075897299e0 * t112 * t92 * t213 * t115 - 0.10320064155614252239e-3 * t283 * t214 * t285 + 0.44313737677495382697e-2 * t309 * t15 + 0.73856229462492304495e-3 * t316 * t88 + 0.44313737677495382697e-2 * t128 * t96;
    const double t322 = t321 * t55;
    const double t324 = t131 * t144 * t53;
    const double t325 = 0.2e1 * t324;
    const double t326 = t41 * t41;
    const double t327 = 0.1e1 / t326;
    const double t328 = t136 * t136;
    const double t331 = t133 * t7;
    const double t332 = 0.1e1 / t331;
    const double t333 = t34 * t332;
    const double t335 = -0.2e1 * t134 + 0.2e1 * t333;
    const double t339 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t327 * t328 + 0.4e1 / 0.3e1 * t41 * t335 );
    const double t340 = t46 * t46;
    const double t341 = 0.1e1 / t340;
    const double t342 = t140 * t140;
    const double t345 = -t335;
    const double t349 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t341 * t342 + 0.4e1 / 0.3e1 * t46 * t345 );
    const double t350 = t339 + t349;
    const double t352 = t33 * t350 * t53;
    const double t354 = t4 * t214 * t59;
    const double t356 = t79 * t151;
    const double t361 = 0.1e1 / t150 / t58;
    const double t362 = t9 * t361;
    const double t363 = t154 * t154;
    const double t369 = t230 - 0.111795e1 * t240 + 0.22359e1 * t243;
    const double t375 = ( t354 / 0.9e1 + t84 * t356 * t154 / 0.6e1 + t84 * t362 * t363 / 0.2e1 - t84 * t152 * t369 / 0.4e1 ) * t101 * t103;
    const double t378 = t255 * t58;
    const double t381 = t105 * t154;
    const double t384 = t164 * t65;
    const double t385 = 0.1e1 / t384;
    const double t386 = t385 * t1;
    const double t387 = t386 * t3;
    const double t392 = t165 * t232 * t101;
    const double t401 = t164 * t164;
    const double t403 = 0.1e1 / t401 / t65;
    const double t404 = t403 * t1;
    const double t405 = t404 * t3;
    const double t406 = t169 * t169;
    const double t407 = 0.1e1 / t406;
    const double t412 = t70 * t151;
    const double t413 = t412 * t91;
    const double t414 = t79 * t154;
    const double t418 = t175 * t232;
    const double t423 = t71 * t361;
    const double t427 = t354 / 0.72e2 + t413 * t92 * t414 / 0.3e1 - t418 * t299 / 0.9e1 + 0.2e1 / 0.9e1 * t176 * t229 + 0.2e1 * t423 * t363 - t179 * t369;
    const double t428 = t427 * t182;
    const double t432 = 0.1e1 / t71 / t70;
    const double t433 = t181 * t432;
    const double t434 = t58 * t90;
    const double t435 = t433 * t434;
    const double t440 = 0.51817833333333333333e-2 * t375 * t161 + 0.17272611111111111111e-2 * t160 * t378 + 0.51817833333333333333e-2 * t160 * t381 + 0.40281109727029913917e-1 * t387 * t214 * t170 + 0.80562219454059827833e-1 * t392 * t235 * t238 * t170 - 0.16112443890811965567e0 * t167 * t92 * t213 * t170 - 0.55299776073946902743e-1 * t405 * t214 * t407 + 0.26673100072733151594e-2 * t428 * t58 + 0.4445516678788858599e-3 * t435 * t88 + 0.26673100072733151594e-2 * t183 * t154;
    const double t442 = t440 * t49 * t53;
    const double t444 = t186 * t144 * t53;
    const double t445 = 0.2e1 * t444;
    const double t447 = t75 * t350 * t53;
    const double t451 = t131 * t201 * t53;
    const double t452 = t327 * t193;
    const double t455 = t41 * t34;
    const double t459 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t452 * t136 + 0.8e1 / 0.3e1 * t455 * t332 );
    const double t460 = t341 * t197;
    const double t463 = t46 * t34;
    const double t467 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t460 * t140 - 0.8e1 / 0.3e1 * t463 * t332 );
    const double t468 = t459 + t467;
    const double t470 = t33 * t468 * t53;
    const double t472 = t186 * t201 * t53;
    const double t474 = t75 * t468 * t53;
    const double t479 = 0.2e1 * t451;
    const double t480 = t193 * t193;
    const double t484 = 0.2e1 * t134 + 0.2e1 * t333;
    const double t488 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t327 * t480 + 0.4e1 / 0.3e1 * t41 * t484 );
    const double t489 = t197 * t197;
    const double t492 = -t484;
    const double t496 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t341 * t489 + 0.4e1 / 0.3e1 * t46 * t492 );
    const double t497 = t488 + t496;
    const double t499 = t33 * t497 * t53;
    const double t500 = 0.2e1 * t472;
    const double t502 = t75 * t497 * t53;


    v2rho2_aa = t208 - 0.2e1 * t146 + t210 + 0.2e1 * t190 + t7 * ( t322 - t325 - t352 + t442 + t445 + t447 );
    v2rho2_ab = t208 - t146 + t210 + t190 - t203 + t205 + t7 * ( t322 - t324 - t451 - t470 + t442 + t444 + t472 + t474 );
    v2rho2_bb = t208 - 0.2e1 * t203 + t210 + 0.2e1 * t205 + t7 * ( t322 - t479 - t499 + t442 + t500 + t502 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double& vrho_a, double& vrho_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t84 = t4 * t6;
    constexpr double t92 = t3 * t6;
    constexpr double t101 = t1 * t1;
    constexpr double t103 = 0.1e1 / t3;
    constexpr double t234 = t3 * t3;
    constexpr double t235 = t234 * t5;
    constexpr double t297 = t101 * t234;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t56 = t33 * t55;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t77 = t75 * t49 * t53;
    const double t79 = 0.1e1 / t8 / t7;
    const double t80 = t6 * t79;
    const double t85 = t15 * t15;
    const double t86 = 0.1e1 / t85;
    const double t87 = t9 * t86;
    const double t88 = t4 * t80;
    const double t89 = t88 / 0.12e2;
    const double t90 = 0.1e1 / t13;
    const double t91 = t90 * t1;
    const double t94 = t91 * t92 * t79;
    const double t96 = -t89 - 0.10893333333333333333e1 * t94;
    const double t104 = ( -t4 * t80 * t16 / 0.12e2 - t84 * t87 * t96 / 0.4e1 ) * t101 * t103;
    const double t105 = t5 * t8;
    const double t106 = t105 * t15;
    const double t109 = t22 * t22;
    const double t110 = 0.1e1 / t109;
    const double t112 = t110 * t90 * t1;
    const double t114 = 0.2016e-2 * t110 + 0.1e1;
    const double t115 = 0.1e1 / t114;
    const double t120 = t28 * t16;
    const double t121 = t120 * t90;
    const double t124 = t29 * t86;
    const double t126 = -t121 * t88 / 0.6e1 - t124 * t96;
    const double t127 = 0.1e1 / t29;
    const double t128 = t126 * t127;
    const double t131 = 0.10363566666666666667e-1 * t104 * t106 + 0.15357238326806922974e0 * t112 * t92 * t79 * t115 + 0.44313737677495382697e-2 * t128 * t15;
    const double t132 = t131 * t55;
    const double t133 = t7 * t7;
    const double t134 = 0.1e1 / t133;
    const double t135 = t34 * t134;
    const double t136 = t35 - t135;
    const double t139 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t136 );
    const double t140 = -t136;
    const double t143 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t140 );
    const double t144 = t139 + t143;
    const double t146 = t33 * t144 * t53;
    const double t150 = t58 * t58;
    const double t151 = 0.1e1 / t150;
    const double t152 = t9 * t151;
    const double t154 = -t89 - 0.1676925e1 * t94;
    const double t160 = ( -t4 * t80 * t59 / 0.12e2 - t84 * t152 * t154 / 0.4e1 ) * t101 * t103;
    const double t161 = t105 * t58;
    const double t164 = t65 * t65;
    const double t165 = 0.1e1 / t164;
    const double t167 = t165 * t90 * t1;
    const double t169 = 0.137284639e1 * t165 + 0.1e1;
    const double t170 = 0.1e1 / t169;
    const double t175 = t70 * t59;
    const double t176 = t175 * t90;
    const double t179 = t71 * t151;
    const double t181 = -t176 * t88 / 0.6e1 - t179 * t154;
    const double t182 = 0.1e1 / t71;
    const double t183 = t181 * t182;
    const double t186 = 0.51817833333333333333e-2 * t160 * t161 + 0.12084332918108974175e0 * t167 * t92 * t79 * t170 + 0.26673100072733151594e-2 * t183 * t58;
    const double t188 = t186 * t49 * t53;
    const double t190 = t75 * t144 * t53;
    const double t193 = -t35 - t135;
    const double t196 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t193 );
    const double t197 = -t193;
    const double t200 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t197 );
    const double t201 = t196 + t200;
    const double t203 = t33 * t201 * t53;
    const double t205 = t75 * t201 * t53;
    const double t208 = 0.2e1 * t132;
    const double t210 = 0.2e1 * t188;
    const double t213 = 0.1e1 / t8 / t133;
    const double t214 = t6 * t213;
    const double t216 = t4 * t214 * t16;
    const double t218 = t79 * t86;
    const double t223 = 0.1e1 / t85 / t15;
    const double t224 = t9 * t223;
    const double t225 = t96 * t96;
    const double t229 = t4 * t214;
    const double t230 = t229 / 0.9e1;
    const double t232 = 0.1e1 / t13 / t11;
    const double t233 = t232 * t101;
    const double t236 = t8 * t8;
    const double t238 = 0.1e1 / t236 / t133;
    const double t240 = t233 * t235 * t238;
    const double t243 = t91 * t92 * t213;
    const double t245 = t230 - 0.7262222222222222222e0 * t240 + 0.14524444444444444444e1 * t243;
    const double t251 = ( t216 / 0.9e1 + t84 * t218 * t96 / 0.6e1 + t84 * t224 * t225 / 0.2e1 - t84 * t87 * t245 / 0.4e1 ) * t101 * t103;
    const double t255 = t5 / t236;
    const double t256 = t255 * t15;
    const double t259 = t105 * t96;
    const double t262 = t109 * t22;
    const double t263 = 0.1e1 / t262;
    const double t264 = t263 * t1;
    const double t265 = t264 * t3;
    const double t270 = t110 * t232 * t101;
    const double t279 = t109 * t109;
    const double t281 = 0.1e1 / t279 / t22;
    const double t282 = t281 * t1;
    const double t283 = t282 * t3;
    const double t284 = t114 * t114;
    const double t285 = 0.1e1 / t284;
    const double t290 = t28 * t86;
    const double t291 = t290 * t91;
    const double t292 = t79 * t96;
    const double t296 = t120 * t232;
    const double t298 = t5 * t238;
    const double t299 = t297 * t298;
    const double t304 = t29 * t223;
    const double t308 = t216 / 0.72e2 + t291 * t92 * t292 / 0.3e1 - t296 * t299 / 0.9e1 + 0.2e1 / 0.9e1 * t121 * t229 + 0.2e1 * t304 * t225 - t124 * t245;
    const double t309 = t308 * t127;
    const double t313 = 0.1e1 / t29 / t28;
    const double t314 = t126 * t313;
    const double t315 = t15 * t90;
    const double t316 = t314 * t315;
    const double t321 = 0.10363566666666666667e-1 * t251 * t106 + 0.34545222222222222223e-2 * t104 * t256 + 0.10363566666666666667e-1 * t104 * t259 + 0.51190794422689743247e-1 * t265 * t214 * t115 + 0.10238158884537948649e0 * t270 * t235 * t238 * t115 - 0.20476317769075897299e0 * t112 * t92 * t213 * t115 - 0.10320064155614252239e-3 * t283 * t214 * t285 + 0.44313737677495382697e-2 * t309 * t15 + 0.73856229462492304495e-3 * t316 * t88 + 0.44313737677495382697e-2 * t128 * t96;
    const double t322 = t321 * t55;
    const double t324 = t131 * t144 * t53;
    const double t325 = 0.2e1 * t324;
    const double t326 = t41 * t41;
    const double t327 = 0.1e1 / t326;
    const double t328 = t136 * t136;
    const double t331 = t133 * t7;
    const double t332 = 0.1e1 / t331;
    const double t333 = t34 * t332;
    const double t335 = -0.2e1 * t134 + 0.2e1 * t333;
    const double t339 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t327 * t328 + 0.4e1 / 0.3e1 * t41 * t335 );
    const double t340 = t46 * t46;
    const double t341 = 0.1e1 / t340;
    const double t342 = t140 * t140;
    const double t345 = -t335;
    const double t349 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t341 * t342 + 0.4e1 / 0.3e1 * t46 * t345 );
    const double t350 = t339 + t349;
    const double t352 = t33 * t350 * t53;
    const double t354 = t4 * t214 * t59;
    const double t356 = t79 * t151;
    const double t361 = 0.1e1 / t150 / t58;
    const double t362 = t9 * t361;
    const double t363 = t154 * t154;
    const double t369 = t230 - 0.111795e1 * t240 + 0.22359e1 * t243;
    const double t375 = ( t354 / 0.9e1 + t84 * t356 * t154 / 0.6e1 + t84 * t362 * t363 / 0.2e1 - t84 * t152 * t369 / 0.4e1 ) * t101 * t103;
    const double t378 = t255 * t58;
    const double t381 = t105 * t154;
    const double t384 = t164 * t65;
    const double t385 = 0.1e1 / t384;
    const double t386 = t385 * t1;
    const double t387 = t386 * t3;
    const double t392 = t165 * t232 * t101;
    const double t401 = t164 * t164;
    const double t403 = 0.1e1 / t401 / t65;
    const double t404 = t403 * t1;
    const double t405 = t404 * t3;
    const double t406 = t169 * t169;
    const double t407 = 0.1e1 / t406;
    const double t412 = t70 * t151;
    const double t413 = t412 * t91;
    const double t414 = t79 * t154;
    const double t418 = t175 * t232;
    const double t423 = t71 * t361;
    const double t427 = t354 / 0.72e2 + t413 * t92 * t414 / 0.3e1 - t418 * t299 / 0.9e1 + 0.2e1 / 0.9e1 * t176 * t229 + 0.2e1 * t423 * t363 - t179 * t369;
    const double t428 = t427 * t182;
    const double t432 = 0.1e1 / t71 / t70;
    const double t433 = t181 * t432;
    const double t434 = t58 * t90;
    const double t435 = t433 * t434;
    const double t440 = 0.51817833333333333333e-2 * t375 * t161 + 0.17272611111111111111e-2 * t160 * t378 + 0.51817833333333333333e-2 * t160 * t381 + 0.40281109727029913917e-1 * t387 * t214 * t170 + 0.80562219454059827833e-1 * t392 * t235 * t238 * t170 - 0.16112443890811965567e0 * t167 * t92 * t213 * t170 - 0.55299776073946902743e-1 * t405 * t214 * t407 + 0.26673100072733151594e-2 * t428 * t58 + 0.4445516678788858599e-3 * t435 * t88 + 0.26673100072733151594e-2 * t183 * t154;
    const double t442 = t440 * t49 * t53;
    const double t444 = t186 * t144 * t53;
    const double t445 = 0.2e1 * t444;
    const double t447 = t75 * t350 * t53;
    const double t451 = t131 * t201 * t53;
    const double t452 = t327 * t193;
    const double t455 = t41 * t34;
    const double t459 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t452 * t136 + 0.8e1 / 0.3e1 * t455 * t332 );
    const double t460 = t341 * t197;
    const double t463 = t46 * t34;
    const double t467 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t460 * t140 - 0.8e1 / 0.3e1 * t463 * t332 );
    const double t468 = t459 + t467;
    const double t470 = t33 * t468 * t53;
    const double t472 = t186 * t201 * t53;
    const double t474 = t75 * t468 * t53;
    const double t479 = 0.2e1 * t451;
    const double t480 = t193 * t193;
    const double t484 = 0.2e1 * t134 + 0.2e1 * t333;
    const double t488 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.9e1 * t327 * t480 + 0.4e1 / 0.3e1 * t41 * t484 );
    const double t489 = t197 * t197;
    const double t492 = -t484;
    const double t496 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.9e1 * t341 * t489 + 0.4e1 / 0.3e1 * t46 * t492 );
    const double t497 = t488 + t496;
    const double t499 = t33 * t497 * t53;
    const double t500 = 0.2e1 * t472;
    const double t502 = t75 * t497 * t53;


    vrho_a = t56 + t77 + t7 * ( t132 - t146 + t188 + t190 );
    vrho_b = t56 + t77 + t7 * ( t132 - t203 + t188 + t205 );
    v2rho2_aa = t208 - 0.2e1 * t146 + t210 + 0.2e1 * t190 + t7 * ( t322 - t325 - t352 + t442 + t445 + t447 );
    v2rho2_ab = t208 - t146 + t210 + t190 - t203 + t205 + t7 * ( t322 - t324 - t451 - t470 + t442 + t444 + t472 + t474 );
    v2rho2_bb = t208 - 0.2e1 * t203 + t210 + 0.2e1 * t205 + t7 * ( t322 - t479 - t499 + t442 + t500 + t502 );

  }


};

struct BuiltinVWN_RPA : detail::BuiltinKernelImpl< BuiltinVWN_RPA > {

  BuiltinVWN_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN_RPA >(p) { }
  
  virtual ~BuiltinVWN_RPA() = default;

};



} // namespace ExchCXX
