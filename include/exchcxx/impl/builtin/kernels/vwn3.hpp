#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinVWN3 > :
  public lda_screening_interface< BuiltinVWN3 > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t90 = constants::m_pi_sq;
    constexpr double t118 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t91 = 0.1e1 / t90;
    constexpr double t119 = t118 - 0.1e1;
    constexpr double t121 = 0.1e1 / t119 / 0.2e1;
    constexpr double t122 = 0.9e1 * t119;
    constexpr double t123 = t121 * t122;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.186372e1 * t12 + 0.129352e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.372744e1;
    const double t24 = safe_math::atan( 0.61519908197590802322e1 / t21 );
    const double t25 = 0.38783294878113014393e-1 * t24;
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.10498e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t31 = 0.96902277115443742139e-3 * t30;
    const double t33 = t11 + 0.353021e1 * t12 + 0.180578e2;
    const double t34 = 0.1e1 / t33;
    const double t38 = safe_math::log( t4 * t9 * t34 / 0.4e1 );
    const double t40 = t12 + 0.706042e1;
    const double t43 = safe_math::atan( 0.473092690956011283e1 / t40 );
    const double t45 = t26 + 0.325e0;
    const double t46 = t45 * t45;
    const double t48 = safe_math::log( t46 * t34 );
    const double t50 = 0.1554535e-1 * t38 + 0.52491393169780936218e-1 * t43 + 0.22478670955426118383e-2 * t48 - t20 - t25 - t31;
    const double t52 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t53 = 0.1e1 / t52;
    const double t57 = safe_math::log( t4 * t9 * t53 / 0.4e1 );
    const double t59 = t12 + 0.201231e2;
    const double t62 = safe_math::atan( 0.11716852777089929792e1 / t59 );
    const double t64 = t26 + 0.743294e0;
    const double t65 = t64 * t64;
    const double t67 = safe_math::log( t65 * t53 );
    const double t70 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t71 = 0.1e1 / t70;
    const double t75 = safe_math::log( t4 * t9 * t71 / 0.4e1 );
    const double t77 = t12 + 0.13072e2;
    const double t80 = safe_math::atan( 0.44899888641287296627e-1 / t77 );
    const double t82 = t26 + 0.409286e0;
    const double t83 = t82 * t82;
    const double t85 = safe_math::log( t83 * t71 );
    const double t87 = 0.1554535e-1 * t57 + 0.61881802979060631482e0 * t62 + 0.26673100072733151594e-2 * t67 - 0.310907e-1 * t75 - 0.20521972937837502661e2 * t80 - 0.44313737677495382697e-2 * t85;
    const double t88 = 0.1e1 / t87;
    const double t92 = t50 * t88 * t91;
    const double t94 = t11 + 0.534175e0 * t12 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = safe_math::log( t4 * t9 * t95 / 0.4e1 );
    const double t100 = t12 + 0.106835e1;
    const double t103 = safe_math::atan( 0.6692072046645941483e1 / t100 );
    const double t105 = t26 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = safe_math::log( t106 * t95 );
    const double t110 = t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108;
    const double t112 = safe_math::cbrt( zeta_tol );
    const double t114 = piecewise_functor_3( 0.1e1 <= zeta_tol, t112 * zeta_tol, 1.0 );
    const double t116 = 0.2e1 * t114 - 0.2e1;
    const double t124 = t110 * t116 * t123;
    const double t126 = t92 * t124 / 0.24e2;


    eps = t20 + t25 + t31 - t126;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t90 = constants::m_pi_sq;
    constexpr double t118 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t91 = 0.1e1 / t90;
    constexpr double t119 = t118 - 0.1e1;
    constexpr double t121 = 0.1e1 / t119 / 0.2e1;
    constexpr double t122 = 0.9e1 * t119;
    constexpr double t123 = t121 * t122;
    constexpr double t133 = t4 * t6;
    constexpr double t141 = t3 * t6;
    constexpr double t150 = t1 * t1;
    constexpr double t152 = 0.1e1 / t3;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.186372e1 * t12 + 0.129352e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.372744e1;
    const double t24 = safe_math::atan( 0.61519908197590802322e1 / t21 );
    const double t25 = 0.38783294878113014393e-1 * t24;
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.10498e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t31 = 0.96902277115443742139e-3 * t30;
    const double t33 = t11 + 0.353021e1 * t12 + 0.180578e2;
    const double t34 = 0.1e1 / t33;
    const double t38 = safe_math::log( t4 * t9 * t34 / 0.4e1 );
    const double t40 = t12 + 0.706042e1;
    const double t43 = safe_math::atan( 0.473092690956011283e1 / t40 );
    const double t45 = t26 + 0.325e0;
    const double t46 = t45 * t45;
    const double t48 = safe_math::log( t46 * t34 );
    const double t50 = 0.1554535e-1 * t38 + 0.52491393169780936218e-1 * t43 + 0.22478670955426118383e-2 * t48 - t20 - t25 - t31;
    const double t52 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t53 = 0.1e1 / t52;
    const double t57 = safe_math::log( t4 * t9 * t53 / 0.4e1 );
    const double t59 = t12 + 0.201231e2;
    const double t62 = safe_math::atan( 0.11716852777089929792e1 / t59 );
    const double t64 = t26 + 0.743294e0;
    const double t65 = t64 * t64;
    const double t67 = safe_math::log( t65 * t53 );
    const double t70 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t71 = 0.1e1 / t70;
    const double t75 = safe_math::log( t4 * t9 * t71 / 0.4e1 );
    const double t77 = t12 + 0.13072e2;
    const double t80 = safe_math::atan( 0.44899888641287296627e-1 / t77 );
    const double t82 = t26 + 0.409286e0;
    const double t83 = t82 * t82;
    const double t85 = safe_math::log( t83 * t71 );
    const double t87 = 0.1554535e-1 * t57 + 0.61881802979060631482e0 * t62 + 0.26673100072733151594e-2 * t67 - 0.310907e-1 * t75 - 0.20521972937837502661e2 * t80 - 0.44313737677495382697e-2 * t85;
    const double t88 = 0.1e1 / t87;
    const double t92 = t50 * t88 * t91;
    const double t94 = t11 + 0.534175e0 * t12 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = safe_math::log( t4 * t9 * t95 / 0.4e1 );
    const double t100 = t12 + 0.106835e1;
    const double t103 = safe_math::atan( 0.6692072046645941483e1 / t100 );
    const double t105 = t26 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = safe_math::log( t106 * t95 );
    const double t110 = t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108;
    const double t112 = safe_math::cbrt( zeta_tol );
    const double t114 = piecewise_functor_3( 0.1e1 <= zeta_tol, t112 * zeta_tol, 1.0 );
    const double t116 = 0.2e1 * t114 - 0.2e1;
    const double t124 = t110 * t116 * t123;
    const double t126 = t92 * t124 / 0.24e2;
    const double t128 = 0.1e1 / t7 / rho;
    const double t129 = t6 * t128;
    const double t134 = t14 * t14;
    const double t135 = 0.1e1 / t134;
    const double t136 = t8 * t135;
    const double t137 = t4 * t129;
    const double t138 = t137 / 0.12e2;
    const double t139 = 0.1e1 / t12;
    const double t140 = t139 * t1;
    const double t143 = t140 * t141 * t128;
    const double t145 = -t138 - 0.31062e0 * t143;
    const double t153 = ( -t4 * t129 * t15 / 0.12e2 - t133 * t136 * t145 / 0.4e1 ) * t150 * t152;
    const double t154 = t5 * t7;
    const double t155 = t154 * t14;
    const double t156 = t153 * t155;
    const double t157 = 0.10363566666666666667e-1 * t156;
    const double t158 = t21 * t21;
    const double t159 = 0.1e1 / t158;
    const double t161 = t159 * t139 * t1;
    const double t163 = 0.378469910464e2 * t159 + 0.1e1;
    const double t164 = 0.1e1 / t163;
    const double t167 = t161 * t141 * t128 * t164;
    const double t168 = 0.39765745675026770179e-1 * t167;
    const double t169 = t27 * t15;
    const double t170 = t169 * t139;
    const double t173 = t28 * t135;
    const double t175 = -t170 * t137 / 0.6e1 - t173 * t145;
    const double t176 = 0.1e1 / t28;
    const double t177 = t175 * t176;
    const double t178 = t177 * t14;
    const double t179 = 0.96902277115443742139e-3 * t178;
    const double t183 = t33 * t33;
    const double t184 = 0.1e1 / t183;
    const double t185 = t8 * t184;
    const double t187 = -t138 - 0.58836833333333333333e0 * t143;
    const double t193 = ( -t4 * t129 * t34 / 0.12e2 - t133 * t185 * t187 / 0.4e1 ) * t150 * t152;
    const double t194 = t154 * t33;
    const double t197 = t40 * t40;
    const double t198 = 0.1e1 / t197;
    const double t200 = t198 * t139 * t1;
    const double t202 = 0.223816694236e2 * t198 + 0.1e1;
    const double t203 = 0.1e1 / t202;
    const double t208 = t45 * t34;
    const double t209 = t208 * t139;
    const double t212 = t46 * t184;
    const double t214 = -t209 * t137 / 0.6e1 - t212 * t187;
    const double t215 = 0.1e1 / t46;
    const double t216 = t214 * t215;
    const double t219 = 0.51817833333333333333e-2 * t193 * t194 + 0.41388824077869423261e-1 * t200 * t141 * t128 * t203 + 0.22478670955426118383e-2 * t216 * t33 - t157 - t168 - t179;
    const double t221 = t219 * t88 * t91;
    const double t222 = t221 * t124;
    const double t224 = t87 * t87;
    const double t225 = 0.1e1 / t224;
    const double t226 = t50 * t225;
    const double t227 = t91 * t110;
    const double t228 = t226 * t227;
    const double t229 = t116 * t121;
    const double t233 = t52 * t52;
    const double t234 = 0.1e1 / t233;
    const double t235 = t8 * t234;
    const double t237 = -t138 - 0.1676925e1 * t143;
    const double t243 = ( -t4 * t129 * t53 / 0.12e2 - t133 * t235 * t237 / 0.4e1 ) * t150 * t152;
    const double t244 = t154 * t52;
    const double t247 = t59 * t59;
    const double t248 = 0.1e1 / t247;
    const double t250 = t248 * t139 * t1;
    const double t252 = 0.137284639e1 * t248 + 0.1e1;
    const double t253 = 0.1e1 / t252;
    const double t258 = t64 * t53;
    const double t259 = t258 * t139;
    const double t262 = t65 * t234;
    const double t264 = -t259 * t137 / 0.6e1 - t262 * t237;
    const double t265 = 0.1e1 / t65;
    const double t266 = t264 * t265;
    const double t272 = t70 * t70;
    const double t273 = 0.1e1 / t272;
    const double t274 = t8 * t273;
    const double t276 = -t138 - 0.10893333333333333333e1 * t143;
    const double t282 = ( -t4 * t129 * t71 / 0.12e2 - t133 * t274 * t276 / 0.4e1 ) * t150 * t152;
    const double t283 = t154 * t70;
    const double t286 = t77 * t77;
    const double t287 = 0.1e1 / t286;
    const double t289 = t287 * t139 * t1;
    const double t291 = 0.2016e-2 * t287 + 0.1e1;
    const double t292 = 0.1e1 / t291;
    const double t297 = t82 * t71;
    const double t298 = t297 * t139;
    const double t301 = t83 * t273;
    const double t303 = -t298 * t137 / 0.6e1 - t301 * t276;
    const double t304 = 0.1e1 / t83;
    const double t305 = t303 * t304;
    const double t308 = 0.51817833333333333333e-2 * t243 * t244 + 0.12084332918108974175e0 * t250 * t141 * t128 * t253 + 0.26673100072733151594e-2 * t266 * t52 - 0.10363566666666666667e-1 * t282 * t283 - 0.15357238326806922974e0 * t289 * t141 * t128 * t292 - 0.44313737677495382697e-2 * t305 * t70;
    const double t309 = t122 * t308;
    const double t310 = t229 * t309;
    const double t311 = t228 * t310;
    const double t316 = t94 * t94;
    const double t317 = 0.1e1 / t316;
    const double t318 = t8 * t317;
    const double t320 = -t138 - 0.89029166666666666667e-1 * t143;
    const double t326 = ( -t4 * t129 * t95 / 0.12e2 - t133 * t318 * t320 / 0.4e1 ) * t150 * t152;
    const double t327 = t154 * t94;
    const double t330 = t100 * t100;
    const double t331 = 0.1e1 / t330;
    const double t333 = t331 * t139 * t1;
    const double t335 = 0.447838282775e2 * t331 + 0.1e1;
    const double t336 = 0.1e1 / t335;
    const double t341 = t105 * t95;
    const double t342 = t341 * t139;
    const double t345 = t106 * t317;
    const double t347 = -t342 * t137 / 0.6e1 - t345 * t320;
    const double t348 = 0.1e1 / t106;
    const double t349 = t347 * t348;
    const double t352 = t326 * t327 / 0.3e1 + 0.36052240899892258525e0 * t333 * t141 * t128 * t336 + 0.21608710360898267022e-1 * t349 * t94;
    const double t354 = t352 * t116 * t123;
    const double t355 = t92 * t354;


    eps = t20 + t25 + t31 - t126;
    vrho = t20 + t25 + t31 - t126 + rho * ( t157 + t168 + t179 - t222 / 0.24e2 + t311 / 0.24e2 - t355 / 0.24e2 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t91 = constants::m_pi_sq;
    constexpr double t129 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t92 = 0.1e1 / t91;
    constexpr double t130 = t129 - 0.1e1;
    constexpr double t132 = 0.1e1 / t130 / 0.2e1;
    constexpr double t141 = 0.9e1 * t130;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.186372e1 * t13 + 0.129352e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t21 = 0.310907e-1 * t20;
    const double t22 = t13 + 0.372744e1;
    const double t25 = safe_math::atan( 0.61519908197590802322e1 / t22 );
    const double t26 = 0.38783294878113014393e-1 * t25;
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.10498e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t32 = 0.96902277115443742139e-3 * t31;
    const double t34 = t12 + 0.353021e1 * t13 + 0.180578e2;
    const double t35 = 0.1e1 / t34;
    const double t39 = safe_math::log( t4 * t10 * t35 / 0.4e1 );
    const double t41 = t13 + 0.706042e1;
    const double t44 = safe_math::atan( 0.473092690956011283e1 / t41 );
    const double t46 = t27 + 0.325e0;
    const double t47 = t46 * t46;
    const double t49 = safe_math::log( t47 * t35 );
    const double t51 = 0.1554535e-1 * t39 + 0.52491393169780936218e-1 * t44 + 0.22478670955426118383e-2 * t49 - t21 - t26 - t32;
    const double t53 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t54 = 0.1e1 / t53;
    const double t58 = safe_math::log( t4 * t10 * t54 / 0.4e1 );
    const double t60 = t13 + 0.201231e2;
    const double t63 = safe_math::atan( 0.11716852777089929792e1 / t60 );
    const double t65 = t27 + 0.743294e0;
    const double t66 = t65 * t65;
    const double t68 = safe_math::log( t66 * t54 );
    const double t71 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t72 = 0.1e1 / t71;
    const double t76 = safe_math::log( t4 * t10 * t72 / 0.4e1 );
    const double t78 = t13 + 0.13072e2;
    const double t81 = safe_math::atan( 0.44899888641287296627e-1 / t78 );
    const double t83 = t27 + 0.409286e0;
    const double t84 = t83 * t83;
    const double t86 = safe_math::log( t84 * t72 );
    const double t88 = 0.1554535e-1 * t58 + 0.61881802979060631482e0 * t63 + 0.26673100072733151594e-2 * t68 - 0.310907e-1 * t76 - 0.20521972937837502661e2 * t81 - 0.44313737677495382697e-2 * t86;
    const double t89 = 0.1e1 / t88;
    const double t90 = t51 * t89;
    const double t94 = t12 + 0.534175e0 * t13 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = safe_math::log( t4 * t10 * t95 / 0.4e1 );
    const double t100 = t13 + 0.106835e1;
    const double t103 = safe_math::atan( 0.6692072046645941483e1 / t100 );
    const double t105 = t27 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = safe_math::log( t106 * t95 );
    const double t111 = t92 * ( t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108 );
    const double t112 = t90 * t111;
    const double t113 = rho_a - rho_b;
    const double t114 = 0.1e1 / t7;
    const double t115 = t113 * t114;
    const double t116 = 0.1e1 + t115;
    const double t117 = t116 <= zeta_tol;
    const double t118 = safe_math::cbrt( zeta_tol );
    const double t119 = t118 * zeta_tol;
    const double t120 = safe_math::cbrt( t116 );
    const double t122 = piecewise_functor_3( t117, t119, t120 * t116 );
    const double t123 = 0.1e1 - t115;
    const double t124 = t123 <= zeta_tol;
    const double t125 = safe_math::cbrt( t123 );
    const double t127 = piecewise_functor_3( t124, t119, t125 * t123 );
    const double t128 = t122 + t127 - 0.2e1;
    const double t133 = t128 * t132;
    const double t134 = t113 * t113;
    const double t135 = t134 * t134;
    const double t136 = t7 * t7;
    const double t137 = t136 * t136;
    const double t138 = 0.1e1 / t137;
    const double t140 = -t135 * t138 + 0.1e1;
    const double t142 = t140 * t141;
    const double t143 = t133 * t142;
    const double t145 = t112 * t143 / 0.24e2;
    const double t146 = t51 * t128;
    const double t147 = t132 * t135;
    const double t148 = t147 * t138;
    const double t149 = t146 * t148;


    eps = t21 + t26 + t32 - t145 + t149;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t91 = constants::m_pi_sq;
    constexpr double t129 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t92 = 0.1e1 / t91;
    constexpr double t130 = t129 - 0.1e1;
    constexpr double t132 = 0.1e1 / t130 / 0.2e1;
    constexpr double t141 = 0.9e1 * t130;
    constexpr double t156 = t4 * t6;
    constexpr double t164 = t3 * t6;
    constexpr double t173 = t1 * t1;
    constexpr double t175 = 0.1e1 / t3;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.186372e1 * t13 + 0.129352e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t21 = 0.310907e-1 * t20;
    const double t22 = t13 + 0.372744e1;
    const double t25 = safe_math::atan( 0.61519908197590802322e1 / t22 );
    const double t26 = 0.38783294878113014393e-1 * t25;
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.10498e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t32 = 0.96902277115443742139e-3 * t31;
    const double t34 = t12 + 0.353021e1 * t13 + 0.180578e2;
    const double t35 = 0.1e1 / t34;
    const double t39 = safe_math::log( t4 * t10 * t35 / 0.4e1 );
    const double t41 = t13 + 0.706042e1;
    const double t44 = safe_math::atan( 0.473092690956011283e1 / t41 );
    const double t46 = t27 + 0.325e0;
    const double t47 = t46 * t46;
    const double t49 = safe_math::log( t47 * t35 );
    const double t51 = 0.1554535e-1 * t39 + 0.52491393169780936218e-1 * t44 + 0.22478670955426118383e-2 * t49 - t21 - t26 - t32;
    const double t53 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t54 = 0.1e1 / t53;
    const double t58 = safe_math::log( t4 * t10 * t54 / 0.4e1 );
    const double t60 = t13 + 0.201231e2;
    const double t63 = safe_math::atan( 0.11716852777089929792e1 / t60 );
    const double t65 = t27 + 0.743294e0;
    const double t66 = t65 * t65;
    const double t68 = safe_math::log( t66 * t54 );
    const double t71 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t72 = 0.1e1 / t71;
    const double t76 = safe_math::log( t4 * t10 * t72 / 0.4e1 );
    const double t78 = t13 + 0.13072e2;
    const double t81 = safe_math::atan( 0.44899888641287296627e-1 / t78 );
    const double t83 = t27 + 0.409286e0;
    const double t84 = t83 * t83;
    const double t86 = safe_math::log( t84 * t72 );
    const double t88 = 0.1554535e-1 * t58 + 0.61881802979060631482e0 * t63 + 0.26673100072733151594e-2 * t68 - 0.310907e-1 * t76 - 0.20521972937837502661e2 * t81 - 0.44313737677495382697e-2 * t86;
    const double t89 = 0.1e1 / t88;
    const double t90 = t51 * t89;
    const double t94 = t12 + 0.534175e0 * t13 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = safe_math::log( t4 * t10 * t95 / 0.4e1 );
    const double t100 = t13 + 0.106835e1;
    const double t103 = safe_math::atan( 0.6692072046645941483e1 / t100 );
    const double t105 = t27 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = safe_math::log( t106 * t95 );
    const double t111 = t92 * ( t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108 );
    const double t112 = t90 * t111;
    const double t113 = rho_a - rho_b;
    const double t114 = 0.1e1 / t7;
    const double t115 = t113 * t114;
    const double t116 = 0.1e1 + t115;
    const double t117 = t116 <= zeta_tol;
    const double t118 = safe_math::cbrt( zeta_tol );
    const double t119 = t118 * zeta_tol;
    const double t120 = safe_math::cbrt( t116 );
    const double t122 = piecewise_functor_3( t117, t119, t120 * t116 );
    const double t123 = 0.1e1 - t115;
    const double t124 = t123 <= zeta_tol;
    const double t125 = safe_math::cbrt( t123 );
    const double t127 = piecewise_functor_3( t124, t119, t125 * t123 );
    const double t128 = t122 + t127 - 0.2e1;
    const double t133 = t128 * t132;
    const double t134 = t113 * t113;
    const double t135 = t134 * t134;
    const double t136 = t7 * t7;
    const double t137 = t136 * t136;
    const double t138 = 0.1e1 / t137;
    const double t140 = -t135 * t138 + 0.1e1;
    const double t142 = t140 * t141;
    const double t143 = t133 * t142;
    const double t145 = t112 * t143 / 0.24e2;
    const double t146 = t51 * t128;
    const double t147 = t132 * t135;
    const double t148 = t147 * t138;
    const double t149 = t146 * t148;
    const double t151 = 0.1e1 / t8 / t7;
    const double t152 = t6 * t151;
    const double t157 = t15 * t15;
    const double t158 = 0.1e1 / t157;
    const double t159 = t9 * t158;
    const double t160 = t4 * t152;
    const double t161 = t160 / 0.12e2;
    const double t162 = 0.1e1 / t13;
    const double t163 = t162 * t1;
    const double t166 = t163 * t164 * t151;
    const double t168 = -t161 - 0.31062e0 * t166;
    const double t176 = ( -t4 * t152 * t16 / 0.12e2 - t156 * t159 * t168 / 0.4e1 ) * t173 * t175;
    const double t177 = t5 * t8;
    const double t178 = t177 * t15;
    const double t179 = t176 * t178;
    const double t180 = 0.10363566666666666667e-1 * t179;
    const double t181 = t22 * t22;
    const double t182 = 0.1e1 / t181;
    const double t184 = t182 * t162 * t1;
    const double t186 = 0.378469910464e2 * t182 + 0.1e1;
    const double t187 = 0.1e1 / t186;
    const double t190 = t184 * t164 * t151 * t187;
    const double t191 = 0.39765745675026770179e-1 * t190;
    const double t192 = t28 * t16;
    const double t193 = t192 * t162;
    const double t196 = t29 * t158;
    const double t198 = -t193 * t160 / 0.6e1 - t196 * t168;
    const double t199 = 0.1e1 / t29;
    const double t200 = t198 * t199;
    const double t201 = t200 * t15;
    const double t202 = 0.96902277115443742139e-3 * t201;
    const double t206 = t34 * t34;
    const double t207 = 0.1e1 / t206;
    const double t208 = t9 * t207;
    const double t210 = -t161 - 0.58836833333333333333e0 * t166;
    const double t216 = ( -t4 * t152 * t35 / 0.12e2 - t156 * t208 * t210 / 0.4e1 ) * t173 * t175;
    const double t217 = t177 * t34;
    const double t220 = t41 * t41;
    const double t221 = 0.1e1 / t220;
    const double t223 = t221 * t162 * t1;
    const double t225 = 0.223816694236e2 * t221 + 0.1e1;
    const double t226 = 0.1e1 / t225;
    const double t231 = t46 * t35;
    const double t232 = t231 * t162;
    const double t235 = t47 * t207;
    const double t237 = -t232 * t160 / 0.6e1 - t235 * t210;
    const double t238 = 0.1e1 / t47;
    const double t239 = t237 * t238;
    const double t242 = 0.51817833333333333333e-2 * t216 * t217 + 0.41388824077869423261e-1 * t223 * t164 * t151 * t226 + 0.22478670955426118383e-2 * t239 * t34 - t180 - t191 - t202;
    const double t243 = t242 * t89;
    const double t244 = t243 * t111;
    const double t245 = t244 * t143;
    const double t246 = t245 / 0.24e2;
    const double t247 = t88 * t88;
    const double t248 = 0.1e1 / t247;
    const double t249 = t51 * t248;
    const double t250 = t249 * t111;
    const double t254 = t53 * t53;
    const double t255 = 0.1e1 / t254;
    const double t256 = t9 * t255;
    const double t258 = -t161 - 0.1676925e1 * t166;
    const double t264 = ( -t4 * t152 * t54 / 0.12e2 - t156 * t256 * t258 / 0.4e1 ) * t173 * t175;
    const double t265 = t177 * t53;
    const double t268 = t60 * t60;
    const double t269 = 0.1e1 / t268;
    const double t271 = t269 * t162 * t1;
    const double t273 = 0.137284639e1 * t269 + 0.1e1;
    const double t274 = 0.1e1 / t273;
    const double t279 = t65 * t54;
    const double t280 = t279 * t162;
    const double t283 = t66 * t255;
    const double t285 = -t280 * t160 / 0.6e1 - t283 * t258;
    const double t286 = 0.1e1 / t66;
    const double t287 = t285 * t286;
    const double t293 = t71 * t71;
    const double t294 = 0.1e1 / t293;
    const double t295 = t9 * t294;
    const double t297 = -t161 - 0.10893333333333333333e1 * t166;
    const double t303 = ( -t4 * t152 * t72 / 0.12e2 - t156 * t295 * t297 / 0.4e1 ) * t173 * t175;
    const double t304 = t177 * t71;
    const double t307 = t78 * t78;
    const double t308 = 0.1e1 / t307;
    const double t310 = t308 * t162 * t1;
    const double t312 = 0.2016e-2 * t308 + 0.1e1;
    const double t313 = 0.1e1 / t312;
    const double t318 = t83 * t72;
    const double t319 = t318 * t162;
    const double t322 = t84 * t294;
    const double t324 = -t319 * t160 / 0.6e1 - t322 * t297;
    const double t325 = 0.1e1 / t84;
    const double t326 = t324 * t325;
    const double t329 = 0.51817833333333333333e-2 * t264 * t265 + 0.12084332918108974175e0 * t271 * t164 * t151 * t274 + 0.26673100072733151594e-2 * t287 * t53 - 0.10363566666666666667e-1 * t303 * t304 - 0.15357238326806922974e0 * t310 * t164 * t151 * t313 - 0.44313737677495382697e-2 * t326 * t71;
    const double t330 = t142 * t329;
    const double t331 = t133 * t330;
    const double t332 = t250 * t331;
    const double t333 = t332 / 0.24e2;
    const double t337 = t94 * t94;
    const double t338 = 0.1e1 / t337;
    const double t339 = t9 * t338;
    const double t341 = -t161 - 0.89029166666666666667e-1 * t166;
    const double t347 = ( -t4 * t152 * t95 / 0.12e2 - t156 * t339 * t341 / 0.4e1 ) * t173 * t175;
    const double t348 = t177 * t94;
    const double t351 = t100 * t100;
    const double t352 = 0.1e1 / t351;
    const double t354 = t352 * t162 * t1;
    const double t356 = 0.447838282775e2 * t352 + 0.1e1;
    const double t357 = 0.1e1 / t356;
    const double t362 = t105 * t95;
    const double t363 = t362 * t162;
    const double t366 = t106 * t338;
    const double t368 = -t363 * t160 / 0.6e1 - t366 * t341;
    const double t369 = 0.1e1 / t106;
    const double t370 = t368 * t369;
    const double t374 = t92 * ( t347 * t348 / 0.3e1 + 0.36052240899892258525e0 * t354 * t164 * t151 * t357 + 0.21608710360898267022e-1 * t370 * t94 );
    const double t375 = t90 * t374;
    const double t376 = t375 * t143;
    const double t377 = t376 / 0.24e2;
    const double t378 = 0.1e1 / t136;
    const double t379 = t113 * t378;
    const double t380 = t114 - t379;
    const double t383 = piecewise_functor_3( t117, 0.0, 0.4e1 / 0.3e1 * t120 * t380 );
    const double t384 = -t380;
    const double t387 = piecewise_functor_3( t124, 0.0, 0.4e1 / 0.3e1 * t125 * t384 );
    const double t388 = t383 + t387;
    const double t389 = t388 * t132;
    const double t390 = t389 * t142;
    const double t391 = t112 * t390;
    const double t392 = t391 / 0.24e2;
    const double t393 = t134 * t113;
    const double t394 = t393 * t138;
    const double t395 = t137 * t7;
    const double t396 = 0.1e1 / t395;
    const double t397 = t135 * t396;
    const double t399 = -0.4e1 * t394 + 0.4e1 * t397;
    const double t400 = t399 * t141;
    const double t401 = t133 * t400;
    const double t402 = t112 * t401;
    const double t403 = t402 / 0.24e2;
    const double t404 = t242 * t128;
    const double t405 = t404 * t148;
    const double t406 = t51 * t388;
    const double t407 = t406 * t148;
    const double t408 = t132 * t393;
    const double t409 = t408 * t138;
    const double t410 = t146 * t409;
    const double t411 = 0.4e1 * t410;
    const double t412 = t147 * t396;
    const double t413 = t146 * t412;
    const double t414 = 0.4e1 * t413;
    const double t415 = t180 + t191 + t202 - t246 + t333 - t377 - t392 - t403 + t405 + t407 + t411 - t414;
    const double t417 = -t114 - t379;
    const double t420 = piecewise_functor_3( t117, 0.0, 0.4e1 / 0.3e1 * t120 * t417 );
    const double t421 = -t417;
    const double t424 = piecewise_functor_3( t124, 0.0, 0.4e1 / 0.3e1 * t125 * t421 );
    const double t425 = t420 + t424;
    const double t426 = t425 * t132;
    const double t427 = t426 * t142;
    const double t428 = t112 * t427;
    const double t429 = t428 / 0.24e2;
    const double t431 = 0.4e1 * t394 + 0.4e1 * t397;
    const double t432 = t431 * t141;
    const double t433 = t133 * t432;
    const double t434 = t112 * t433;
    const double t435 = t434 / 0.24e2;
    const double t436 = t51 * t425;
    const double t437 = t436 * t148;
    const double t438 = t180 + t191 + t202 - t246 + t333 - t377 - t429 - t435 + t405 + t437 - t411 - t414;


    eps = t21 + t26 + t32 - t145 + t149;
    vrho_a = t7 * t415 - t145 + t149 + t21 + t26 + t32;
    vrho_b = t7 * t438 - t145 + t149 + t21 + t26 + t32;

  }


};

struct BuiltinVWN3 : detail::BuiltinKernelImpl< BuiltinVWN3 > {

  BuiltinVWN3( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN3 >(p) { }
  
  virtual ~BuiltinVWN3() = default;

};



} // namespace ExchCXX