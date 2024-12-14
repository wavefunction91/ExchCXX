#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPBE_C > :
  public gga_screening_interface< BuiltinPBE_C > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-12;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000021e-16;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double beta = 0.06672455060314922;
  static constexpr double gamma = 0.031090690869654895034;
  static constexpr double BB = 1.;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t70 = 0.1e1 / t3;
    constexpr double t75 = BB * beta;
    constexpr double t76 = 0.1e1 / gamma;
    constexpr double t91 = t39 * t39;
    constexpr double t96 = 0.1e1 / t19;
    constexpr double t97 = t1 * t96;
    constexpr double t98 = t97 * t6;
    constexpr double t104 = beta * t76;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.53425e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.379785e1 * t13 + 0.8969e0 * t10 + 0.204775e0 * t16 + 0.123235e0 * t24;
    const double t29 = 0.1e1 + 0.16081979498692535067e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t32 = 0.621814e-1 * t12 * t30;
    const double t33 = 0.1e1 <= zeta_tol;
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( t33, t34 * zeta_tol, 1.0 );
    const double t43 = ( 0.2e1 * t36 - 0.2e1 ) / ( 0.2e1 * t39 - 0.2e1 );
    const double t45 = 0.1e1 + 0.278125e-1 * t10;
    const double t50 = 0.51785e1 * t13 + 0.905775e0 * t10 + 0.1100325e0 * t16 + 0.1241775e0 * t24;
    const double t53 = 0.1e1 + 0.29608749977793437516e2 / t50;
    const double t54 = safe_math::log( t53 );
    const double t57 = 0.19751673498613801407e-1 * t43 * t45 * t54;
    const double t58 = t34 * t34;
    const double t59 = piecewise_functor_3( t33, t58, 1.0 );
    const double t60 = t59 * t59;
    const double t61 = t60 * t59;
    const double t62 = gamma * t61;
    const double t63 = rho * rho;
    const double t65 = 0.1e1 / t7 / t63;
    const double t68 = 0.1e1 / t60;
    const double t72 = t68 * t18 * t70 * t5;
    const double t79 = 0.1e1 / t61;
    const double t81 = safe_math::exp( -( -t32 + t57 ) * t76 * t79 );
    const double t82 = t81 - 0.1e1;
    const double t83 = 0.1e1 / t82;
    const double t84 = t76 * t83;
    const double t85 = sigma * sigma;
    const double t87 = t75 * t84 * t85;
    const double t88 = t63 * t63;
    const double t90 = 0.1e1 / t21 / t88;
    const double t92 = t90 * t91;
    const double t93 = t60 * t60;
    const double t94 = 0.1e1 / t93;
    const double t95 = t92 * t94;
    const double t99 = t95 * t98;
    const double t102 = sigma * t65 * t39 * t72 / 0.96e2 + t87 * t99 / 0.3072e4;
    const double t103 = beta * t102;
    const double t107 = t104 * t83 * t102 + 0.1e1;
    const double t108 = 0.1e1 / t107;
    const double t109 = t76 * t108;
    const double t111 = t103 * t109 + 0.1e1;
    const double t112 = safe_math::log( t111 );
    const double t113 = t62 * t112;


    eps = -t32 + t57 + t113;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t70 = 0.1e1 / t3;
    constexpr double t75 = BB * beta;
    constexpr double t76 = 0.1e1 / gamma;
    constexpr double t91 = t39 * t39;
    constexpr double t96 = 0.1e1 / t19;
    constexpr double t97 = t1 * t96;
    constexpr double t98 = t97 * t6;
    constexpr double t104 = beta * t76;
    constexpr double t125 = t3 * t6;
    constexpr double t170 = gamma * gamma;
    constexpr double t171 = 0.1e1 / t170;
    constexpr double t172 = t75 * t171;
    constexpr double t182 = t96 * t6;
    constexpr double t224 = t18 * t70 * t5;
    constexpr double t234 = beta * beta;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.53425e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.379785e1 * t13 + 0.8969e0 * t10 + 0.204775e0 * t16 + 0.123235e0 * t24;
    const double t29 = 0.1e1 + 0.16081979498692535067e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t32 = 0.621814e-1 * t12 * t30;
    const double t33 = 0.1e1 <= zeta_tol;
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( t33, t34 * zeta_tol, 1.0 );
    const double t43 = ( 0.2e1 * t36 - 0.2e1 ) / ( 0.2e1 * t39 - 0.2e1 );
    const double t45 = 0.1e1 + 0.278125e-1 * t10;
    const double t50 = 0.51785e1 * t13 + 0.905775e0 * t10 + 0.1100325e0 * t16 + 0.1241775e0 * t24;
    const double t53 = 0.1e1 + 0.29608749977793437516e2 / t50;
    const double t54 = safe_math::log( t53 );
    const double t57 = 0.19751673498613801407e-1 * t43 * t45 * t54;
    const double t58 = t34 * t34;
    const double t59 = piecewise_functor_3( t33, t58, 1.0 );
    const double t60 = t59 * t59;
    const double t61 = t60 * t59;
    const double t62 = gamma * t61;
    const double t63 = rho * rho;
    const double t65 = 0.1e1 / t7 / t63;
    const double t68 = 0.1e1 / t60;
    const double t72 = t68 * t18 * t70 * t5;
    const double t79 = 0.1e1 / t61;
    const double t81 = safe_math::exp( -( -t32 + t57 ) * t76 * t79 );
    const double t82 = t81 - 0.1e1;
    const double t83 = 0.1e1 / t82;
    const double t84 = t76 * t83;
    const double t85 = sigma * sigma;
    const double t87 = t75 * t84 * t85;
    const double t88 = t63 * t63;
    const double t90 = 0.1e1 / t21 / t88;
    const double t92 = t90 * t91;
    const double t93 = t60 * t60;
    const double t94 = 0.1e1 / t93;
    const double t95 = t92 * t94;
    const double t99 = t95 * t98;
    const double t102 = sigma * t65 * t39 * t72 / 0.96e2 + t87 * t99 / 0.3072e4;
    const double t103 = beta * t102;
    const double t107 = t104 * t83 * t102 + 0.1e1;
    const double t108 = 0.1e1 / t107;
    const double t109 = t76 * t108;
    const double t111 = t103 * t109 + 0.1e1;
    const double t112 = safe_math::log( t111 );
    const double t113 = t62 * t112;
    const double t115 = 0.1e1 / t7 / rho;
    const double t116 = t6 * t115;
    const double t118 = t4 * t116 * t30;
    const double t119 = 0.11073470983333333333e-2 * t118;
    const double t120 = t26 * t26;
    const double t121 = 0.1e1 / t120;
    const double t122 = t12 * t121;
    const double t124 = 0.1e1 / t13 * t1;
    const double t126 = t125 * t115;
    const double t127 = t124 * t126;
    const double t129 = t4 * t116;
    const double t131 = safe_math::sqrt( t10 );
    const double t132 = t131 * t1;
    const double t133 = t132 * t126;
    const double t138 = t20 * t5 / t21 / rho;
    const double t140 = -0.632975e0 * t127 - 0.29896666666666666667e0 * t129 - 0.1023875e0 * t133 - 0.82156666666666666667e-1 * t138;
    const double t141 = 0.1e1 / t29;
    const double t142 = t140 * t141;
    const double t143 = t122 * t142;
    const double t144 = 0.1e1 * t143;
    const double t145 = t43 * t1;
    const double t148 = t145 * t125 * t115 * t54;
    const double t149 = 0.18311447306006545054e-3 * t148;
    const double t150 = t43 * t45;
    const double t151 = t50 * t50;
    const double t152 = 0.1e1 / t151;
    const double t157 = -0.86308333333333333334e0 * t127 - 0.301925e0 * t129 - 0.5501625e-1 * t133 - 0.82785e-1 * t138;
    const double t159 = 0.1e1 / t53;
    const double t160 = t152 * t157 * t159;
    const double t161 = t150 * t160;
    const double t162 = 0.5848223622634646207e0 * t161;
    const double t163 = t63 * rho;
    const double t165 = 0.1e1 / t7 / t163;
    const double t173 = t82 * t82;
    const double t174 = 0.1e1 / t173;
    const double t175 = t174 * t85;
    const double t176 = t175 * t90;
    const double t177 = t172 * t176;
    const double t179 = 0.1e1 / t93 / t61;
    const double t180 = t91 * t179;
    const double t181 = t180 * t1;
    const double t183 = t119 + t144 - t149 - t162;
    const double t184 = t183 * t81;
    const double t185 = t182 * t184;
    const double t186 = t181 * t185;
    const double t189 = t88 * rho;
    const double t191 = 0.1e1 / t21 / t189;
    const double t192 = t191 * t91;
    const double t193 = t192 * t94;
    const double t194 = t193 * t98;
    const double t197 = -0.7e1 / 0.288e3 * sigma * t165 * t39 * t72 + t177 * t186 / 0.3072e4 - 0.7e1 / 0.4608e4 * t87 * t194;
    const double t198 = beta * t197;
    const double t200 = t107 * t107;
    const double t201 = 0.1e1 / t200;
    const double t202 = t76 * t201;
    const double t204 = beta * t171 * t174;
    const double t206 = t79 * t81;
    const double t211 = t204 * t102 * t183 * t206 + t104 * t83 * t197;
    const double t212 = t202 * t211;
    const double t214 = -t103 * t212 + t198 * t109;
    const double t215 = 0.1e1 / t111;
    const double t217 = t62 * t214 * t215;
    const double t220 = rho * gamma;
    const double t228 = t75 * t84 * sigma;
    const double t231 = t65 * t39 * t68 * t224 / 0.96e2 + t228 * t99 / 0.1536e4;
    const double t232 = beta * t231;
    const double t235 = t234 * t102;
    const double t236 = t235 * t171;
    const double t237 = t201 * t83;
    const double t238 = t237 * t231;
    const double t240 = t232 * t109 - t236 * t238;


    eps = -t32 + t57 + t113;
    vrho = -t32 + t57 + t113 + rho * ( t119 + t144 - t149 - t162 + t217 );
    vsigma = t220 * t61 * t240 * t215;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t56 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t19 = t1 * t1;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = t19 * t20;
    constexpr double t110 = 0.1e1 / t3;
    constexpr double t111 = t110 * t5;
    constexpr double t115 = BB * beta;
    constexpr double t116 = 0.1e1 / gamma;
    constexpr double t130 = t56 * t56;
    constexpr double t135 = 0.1e1 / t20;
    constexpr double t136 = t1 * t135;
    constexpr double t137 = t136 * t6;
    constexpr double t143 = beta * t116;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t11 = t4 * t6 / t8;
    const double t13 = 0.1e1 + 0.53425e-1 * t11;
    const double t14 = safe_math::sqrt( t11 );
    const double t17 = pow_3_2( t11 );
    const double t22 = t8 * t8;
    const double t25 = t21 * t5 / t22;
    const double t27 = 0.379785e1 * t14 + 0.8969e0 * t11 + 0.204775e0 * t17 + 0.123235e0 * t25;
    const double t30 = 0.1e1 + 0.16081979498692535067e2 / t27;
    const double t31 = safe_math::log( t30 );
    const double t33 = 0.621814e-1 * t13 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = t34 * t34;
    const double t36 = t35 * t35;
    const double t37 = t7 * t7;
    const double t38 = t37 * t37;
    const double t39 = 0.1e1 / t38;
    const double t40 = t36 * t39;
    const double t41 = 0.1e1 / t7;
    const double t42 = t34 * t41;
    const double t43 = 0.1e1 + t42;
    const double t44 = t43 <= zeta_tol;
    const double t45 = safe_math::cbrt( zeta_tol );
    const double t46 = t45 * zeta_tol;
    const double t47 = safe_math::cbrt( t43 );
    const double t48 = t47 * t43;
    const double t49 = piecewise_functor_3( t44, t46, t48 );
    const double t50 = 0.1e1 - t42;
    const double t51 = t50 <= zeta_tol;
    const double t52 = safe_math::cbrt( t50 );
    const double t53 = t52 * t50;
    const double t54 = piecewise_functor_3( t51, t46, t53 );
    const double t55 = t49 + t54 - 0.2e1;
    const double t59 = 0.1e1 / ( 0.2e1 * t56 - 0.2e1 );
    const double t60 = t55 * t59;
    const double t62 = 0.1e1 + 0.5137e-1 * t11;
    const double t67 = 0.705945e1 * t14 + 0.1549425e1 * t11 + 0.420775e0 * t17 + 0.1562925e0 * t25;
    const double t70 = 0.1e1 + 0.32163958997385070134e2 / t67;
    const double t71 = safe_math::log( t70 );
    const double t75 = 0.1e1 + 0.278125e-1 * t11;
    const double t80 = 0.51785e1 * t14 + 0.905775e0 * t11 + 0.1100325e0 * t17 + 0.1241775e0 * t25;
    const double t83 = 0.1e1 + 0.29608749977793437516e2 / t80;
    const double t84 = safe_math::log( t83 );
    const double t85 = t75 * t84;
    const double t87 = -0.310907e-1 * t62 * t71 + t33 - 0.19751673498613801407e-1 * t85;
    const double t88 = t60 * t87;
    const double t89 = t40 * t88;
    const double t91 = 0.19751673498613801407e-1 * t60 * t85;
    const double t92 = t45 * t45;
    const double t93 = t47 * t47;
    const double t94 = piecewise_functor_3( t44, t92, t93 );
    const double t95 = t52 * t52;
    const double t96 = piecewise_functor_3( t51, t92, t95 );
    const double t98 = t94 / 0.2e1 + t96 / 0.2e1;
    const double t99 = t98 * t98;
    const double t100 = t99 * t98;
    const double t101 = gamma * t100;
    const double t103 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t105 = 0.1e1 / t8 / t37;
    const double t106 = t103 * t105;
    const double t108 = 0.1e1 / t99;
    const double t112 = t108 * t19 * t111;
    const double t118 = ( -t33 + t89 + t91 ) * t116;
    const double t119 = 0.1e1 / t100;
    const double t121 = safe_math::exp( -t118 * t119 );
    const double t122 = t121 - 0.1e1;
    const double t123 = 0.1e1 / t122;
    const double t124 = t116 * t123;
    const double t125 = t103 * t103;
    const double t127 = t115 * t124 * t125;
    const double t129 = 0.1e1 / t22 / t38;
    const double t131 = t129 * t130;
    const double t132 = t99 * t99;
    const double t133 = 0.1e1 / t132;
    const double t134 = t131 * t133;
    const double t138 = t134 * t137;
    const double t141 = t106 * t56 * t112 / 0.96e2 + t127 * t138 / 0.3072e4;
    const double t142 = beta * t141;
    const double t146 = t143 * t123 * t141 + 0.1e1;
    const double t147 = 0.1e1 / t146;
    const double t148 = t116 * t147;
    const double t150 = t142 * t148 + 0.1e1;
    const double t151 = safe_math::log( t150 );
    const double t152 = t101 * t151;


    eps = -t33 + t89 + t91 + t152;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t56 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t19 = t1 * t1;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = t19 * t20;
    constexpr double t110 = 0.1e1 / t3;
    constexpr double t111 = t110 * t5;
    constexpr double t115 = BB * beta;
    constexpr double t116 = 0.1e1 / gamma;
    constexpr double t130 = t56 * t56;
    constexpr double t135 = 0.1e1 / t20;
    constexpr double t136 = t1 * t135;
    constexpr double t137 = t136 * t6;
    constexpr double t143 = beta * t116;
    constexpr double t164 = t3 * t6;
    constexpr double t275 = t19 * t110;
    constexpr double t280 = t115 * t116;
    constexpr double t288 = t135 * t6;
    constexpr double t404 = t275 * t5;
    constexpr double t414 = beta * beta;
    constexpr double t416 = gamma * gamma;
    constexpr double t417 = 0.1e1 / t416;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t11 = t4 * t6 / t8;
    const double t13 = 0.1e1 + 0.53425e-1 * t11;
    const double t14 = safe_math::sqrt( t11 );
    const double t17 = pow_3_2( t11 );
    const double t22 = t8 * t8;
    const double t25 = t21 * t5 / t22;
    const double t27 = 0.379785e1 * t14 + 0.8969e0 * t11 + 0.204775e0 * t17 + 0.123235e0 * t25;
    const double t30 = 0.1e1 + 0.16081979498692535067e2 / t27;
    const double t31 = safe_math::log( t30 );
    const double t33 = 0.621814e-1 * t13 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = t34 * t34;
    const double t36 = t35 * t35;
    const double t37 = t7 * t7;
    const double t38 = t37 * t37;
    const double t39 = 0.1e1 / t38;
    const double t40 = t36 * t39;
    const double t41 = 0.1e1 / t7;
    const double t42 = t34 * t41;
    const double t43 = 0.1e1 + t42;
    const double t44 = t43 <= zeta_tol;
    const double t45 = safe_math::cbrt( zeta_tol );
    const double t46 = t45 * zeta_tol;
    const double t47 = safe_math::cbrt( t43 );
    const double t48 = t47 * t43;
    const double t49 = piecewise_functor_3( t44, t46, t48 );
    const double t50 = 0.1e1 - t42;
    const double t51 = t50 <= zeta_tol;
    const double t52 = safe_math::cbrt( t50 );
    const double t53 = t52 * t50;
    const double t54 = piecewise_functor_3( t51, t46, t53 );
    const double t55 = t49 + t54 - 0.2e1;
    const double t59 = 0.1e1 / ( 0.2e1 * t56 - 0.2e1 );
    const double t60 = t55 * t59;
    const double t62 = 0.1e1 + 0.5137e-1 * t11;
    const double t67 = 0.705945e1 * t14 + 0.1549425e1 * t11 + 0.420775e0 * t17 + 0.1562925e0 * t25;
    const double t70 = 0.1e1 + 0.32163958997385070134e2 / t67;
    const double t71 = safe_math::log( t70 );
    const double t75 = 0.1e1 + 0.278125e-1 * t11;
    const double t80 = 0.51785e1 * t14 + 0.905775e0 * t11 + 0.1100325e0 * t17 + 0.1241775e0 * t25;
    const double t83 = 0.1e1 + 0.29608749977793437516e2 / t80;
    const double t84 = safe_math::log( t83 );
    const double t85 = t75 * t84;
    const double t87 = -0.310907e-1 * t62 * t71 + t33 - 0.19751673498613801407e-1 * t85;
    const double t88 = t60 * t87;
    const double t89 = t40 * t88;
    const double t91 = 0.19751673498613801407e-1 * t60 * t85;
    const double t92 = t45 * t45;
    const double t93 = t47 * t47;
    const double t94 = piecewise_functor_3( t44, t92, t93 );
    const double t95 = t52 * t52;
    const double t96 = piecewise_functor_3( t51, t92, t95 );
    const double t98 = t94 / 0.2e1 + t96 / 0.2e1;
    const double t99 = t98 * t98;
    const double t100 = t99 * t98;
    const double t101 = gamma * t100;
    const double t103 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t105 = 0.1e1 / t8 / t37;
    const double t106 = t103 * t105;
    const double t108 = 0.1e1 / t99;
    const double t112 = t108 * t19 * t111;
    const double t118 = ( -t33 + t89 + t91 ) * t116;
    const double t119 = 0.1e1 / t100;
    const double t121 = safe_math::exp( -t118 * t119 );
    const double t122 = t121 - 0.1e1;
    const double t123 = 0.1e1 / t122;
    const double t124 = t116 * t123;
    const double t125 = t103 * t103;
    const double t127 = t115 * t124 * t125;
    const double t129 = 0.1e1 / t22 / t38;
    const double t131 = t129 * t130;
    const double t132 = t99 * t99;
    const double t133 = 0.1e1 / t132;
    const double t134 = t131 * t133;
    const double t138 = t134 * t137;
    const double t141 = t106 * t56 * t112 / 0.96e2 + t127 * t138 / 0.3072e4;
    const double t142 = beta * t141;
    const double t146 = t143 * t123 * t141 + 0.1e1;
    const double t147 = 0.1e1 / t146;
    const double t148 = t116 * t147;
    const double t150 = t142 * t148 + 0.1e1;
    const double t151 = safe_math::log( t150 );
    const double t152 = t101 * t151;
    const double t154 = 0.1e1 / t8 / t7;
    const double t155 = t6 * t154;
    const double t157 = t4 * t155 * t31;
    const double t158 = 0.11073470983333333333e-2 * t157;
    const double t159 = t27 * t27;
    const double t160 = 0.1e1 / t159;
    const double t161 = t13 * t160;
    const double t163 = 0.1e1 / t14 * t1;
    const double t165 = t164 * t154;
    const double t166 = t163 * t165;
    const double t168 = t4 * t155;
    const double t170 = safe_math::sqrt( t11 );
    const double t171 = t170 * t1;
    const double t172 = t171 * t165;
    const double t177 = t21 * t5 / t22 / t7;
    const double t179 = -0.632975e0 * t166 - 0.29896666666666666667e0 * t168 - 0.1023875e0 * t172 - 0.82156666666666666667e-1 * t177;
    const double t180 = 0.1e1 / t30;
    const double t181 = t179 * t180;
    const double t182 = t161 * t181;
    const double t183 = 0.1e1 * t182;
    const double t184 = t35 * t34;
    const double t185 = t184 * t39;
    const double t186 = t185 * t88;
    const double t187 = 0.4e1 * t186;
    const double t188 = t38 * t7;
    const double t189 = 0.1e1 / t188;
    const double t190 = t36 * t189;
    const double t191 = t190 * t88;
    const double t192 = 0.4e1 * t191;
    const double t193 = 0.1e1 / t37;
    const double t194 = t34 * t193;
    const double t195 = t41 - t194;
    const double t198 = piecewise_functor_3( t44, 0.0, 0.4e1 / 0.3e1 * t47 * t195 );
    const double t199 = -t195;
    const double t202 = piecewise_functor_3( t51, 0.0, 0.4e1 / 0.3e1 * t52 * t199 );
    const double t204 = ( t198 + t202 ) * t59;
    const double t205 = t204 * t87;
    const double t206 = t40 * t205;
    const double t210 = t67 * t67;
    const double t211 = 0.1e1 / t210;
    const double t212 = t62 * t211;
    const double t217 = -0.1176575e1 * t166 - 0.516475e0 * t168 - 0.2103875e0 * t172 - 0.104195e0 * t177;
    const double t218 = 0.1e1 / t70;
    const double t219 = t217 * t218;
    const double t225 = t80 * t80;
    const double t226 = 0.1e1 / t225;
    const double t227 = t75 * t226;
    const double t232 = -0.86308333333333333334e0 * t166 - 0.301925e0 * t168 - 0.5501625e-1 * t172 - 0.82785e-1 * t177;
    const double t233 = 0.1e1 / t83;
    const double t234 = t232 * t233;
    const double t237 = 0.53237641966666666666e-3 * t4 * t155 * t71 + 0.1e1 * t212 * t219 - t158 - t183 + 0.18311447306006545054e-3 * t4 * t155 * t84 + 0.5848223622634646207e0 * t227 * t234;
    const double t238 = t60 * t237;
    const double t239 = t40 * t238;
    const double t240 = t204 * t85;
    const double t241 = 0.19751673498613801407e-1 * t240;
    const double t242 = t60 * t1;
    const double t244 = t164 * t154 * t84;
    const double t245 = t242 * t244;
    const double t246 = 0.18311447306006545054e-3 * t245;
    const double t247 = t60 * t75;
    const double t249 = t226 * t232 * t233;
    const double t250 = t247 * t249;
    const double t251 = 0.5848223622634646207e0 * t250;
    const double t252 = gamma * t99;
    const double t253 = 0.1e1 / t47;
    const double t256 = piecewise_functor_3( t44, 0.0, 0.2e1 / 0.3e1 * t253 * t195 );
    const double t257 = 0.1e1 / t52;
    const double t260 = piecewise_functor_3( t51, 0.0, 0.2e1 / 0.3e1 * t257 * t199 );
    const double t262 = t256 / 0.2e1 + t260 / 0.2e1;
    const double t263 = t151 * t262;
    const double t264 = t252 * t263;
    const double t265 = 0.3e1 * t264;
    const double t266 = t37 * t7;
    const double t268 = 0.1e1 / t8 / t266;
    const double t269 = t103 * t268;
    const double t272 = 0.7e1 / 0.288e3 * t269 * t56 * t112;
    const double t273 = t56 * t119;
    const double t274 = t106 * t273;
    const double t276 = t5 * t262;
    const double t277 = t275 * t276;
    const double t281 = t122 * t122;
    const double t282 = 0.1e1 / t281;
    const double t283 = t282 * t125;
    const double t285 = t280 * t283 * t129;
    const double t286 = t130 * t133;
    const double t287 = t286 * t1;
    const double t290 = ( t158 + t183 + t187 - t192 + t206 + t239 + t241 - t246 - t251 ) * t116;
    const double t292 = t133 * t262;
    const double t295 = 0.3e1 * t118 * t292 - t290 * t119;
    const double t296 = t295 * t121;
    const double t297 = t288 * t296;
    const double t298 = t287 * t297;
    const double t302 = 0.1e1 / t22 / t188;
    const double t303 = t302 * t130;
    const double t304 = t303 * t133;
    const double t305 = t304 * t137;
    const double t307 = 0.7e1 / 0.4608e4 * t127 * t305;
    const double t308 = t123 * t125;
    const double t310 = t280 * t308 * t129;
    const double t312 = 0.1e1 / t132 / t98;
    const double t313 = t130 * t312;
    const double t314 = t313 * t1;
    const double t316 = t314 * t288 * t262;
    const double t319 = -t272 - t274 * t277 / 0.48e2 - t285 * t298 / 0.3072e4 - t307 - t310 * t316 / 0.768e3;
    const double t320 = beta * t319;
    const double t322 = t146 * t146;
    const double t323 = 0.1e1 / t322;
    const double t324 = t116 * t323;
    const double t325 = t143 * t282;
    const double t326 = t141 * t295;
    const double t331 = -t325 * t326 * t121 + t143 * t123 * t319;
    const double t332 = t324 * t331;
    const double t334 = -t142 * t332 + t320 * t148;
    const double t335 = 0.1e1 / t150;
    const double t336 = t334 * t335;
    const double t337 = t101 * t336;
    const double t338 = t158 + t183 + t187 - t192 + t206 + t239 + t241 - t246 - t251 + t265 + t337;
    const double t340 = -t41 - t194;
    const double t343 = piecewise_functor_3( t44, 0.0, 0.4e1 / 0.3e1 * t47 * t340 );
    const double t344 = -t340;
    const double t347 = piecewise_functor_3( t51, 0.0, 0.4e1 / 0.3e1 * t52 * t344 );
    const double t349 = ( t343 + t347 ) * t59;
    const double t350 = t349 * t87;
    const double t351 = t40 * t350;
    const double t352 = t349 * t85;
    const double t353 = 0.19751673498613801407e-1 * t352;
    const double t356 = piecewise_functor_3( t44, 0.0, 0.2e1 / 0.3e1 * t253 * t340 );
    const double t359 = piecewise_functor_3( t51, 0.0, 0.2e1 / 0.3e1 * t257 * t344 );
    const double t361 = t356 / 0.2e1 + t359 / 0.2e1;
    const double t362 = t151 * t361;
    const double t363 = t252 * t362;
    const double t364 = 0.3e1 * t363;
    const double t365 = t5 * t361;
    const double t366 = t275 * t365;
    const double t370 = ( t158 + t183 - t187 - t192 + t351 + t239 + t353 - t246 - t251 ) * t116;
    const double t372 = t133 * t361;
    const double t375 = 0.3e1 * t118 * t372 - t370 * t119;
    const double t376 = t375 * t121;
    const double t377 = t288 * t376;
    const double t378 = t287 * t377;
    const double t382 = t314 * t288 * t361;
    const double t385 = -t272 - t274 * t366 / 0.48e2 - t285 * t378 / 0.3072e4 - t307 - t310 * t382 / 0.768e3;
    const double t386 = beta * t385;
    const double t388 = t141 * t375;
    const double t393 = -t325 * t388 * t121 + t143 * t123 * t385;
    const double t394 = t324 * t393;
    const double t396 = -t142 * t394 + t386 * t148;
    const double t397 = t396 * t335;
    const double t398 = t101 * t397;
    const double t399 = t158 + t183 - t187 - t192 + t351 + t239 + t353 - t246 - t251 + t364 + t398;
    const double t401 = t7 * gamma;
    const double t402 = t105 * t56;
    const double t405 = t402 * t108 * t404;
    const double t408 = t115 * t124 * t103;
    const double t409 = t408 * t138;
    const double t411 = t405 / 0.96e2 + t409 / 0.1536e4;
    const double t412 = beta * t411;
    const double t415 = t414 * t141;
    const double t418 = t415 * t417;
    const double t419 = t323 * t123;
    const double t420 = t419 * t411;
    const double t422 = t412 * t148 - t418 * t420;
    const double t427 = t405 / 0.48e2 + t409 / 0.768e3;
    const double t428 = beta * t427;
    const double t430 = t419 * t427;
    const double t432 = t428 * t148 - t418 * t430;
    const double t433 = t100 * t432;


    eps = -t33 + t89 + t91 + t152;
    vrho_a = t7 * t338 + t152 - t33 + t89 + t91;
    vrho_b = t7 * t399 + t152 - t33 + t89 + t91;
    vsigma_aa = t401 * t100 * t422 * t335;
    vsigma_ab = t401 * t433 * t335;
    vsigma_bb = vsigma_aa;

  }


};

struct BuiltinPBE_C : detail::BuiltinKernelImpl< BuiltinPBE_C > {

  BuiltinPBE_C( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_C >(p) { }
  
  virtual ~BuiltinPBE_C() = default;

};



} // namespace ExchCXX