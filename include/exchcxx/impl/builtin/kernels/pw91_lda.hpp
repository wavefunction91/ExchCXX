#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPW91_LDA > :
  public lda_screening_interface< BuiltinPW91_LDA > {

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

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double pp_0 = 1.;
  static constexpr double pp_1 = 1.;
  static constexpr double pp_2 = 1.;
  static constexpr double a_0 = 0.031091;
  static constexpr double a_1 = 0.015545;
  static constexpr double a_2 = 0.016887;
  static constexpr double alpha1_0 = 0.21370;
  static constexpr double alpha1_1 = 0.20548;
  static constexpr double alpha1_2 = 0.11125;
  static constexpr double beta1_0 = 7.5957;
  static constexpr double beta1_1 = 14.1189;
  static constexpr double beta1_2 = 10.357;
  static constexpr double beta2_0 = 3.5876;
  static constexpr double beta2_1 = 6.1977;
  static constexpr double beta2_2 = 3.6231;
  static constexpr double beta3_0 = 1.6382;
  static constexpr double beta3_1 = 3.3662;
  static constexpr double beta3_2 = 0.88026;
  static constexpr double beta4_0 = 0.49294;
  static constexpr double beta4_1 = 0.62517;
  static constexpr double beta4_2 = 0.49671;
  static constexpr double fz20 = 1.709921;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
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

    (void)(eps);
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
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
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

    (void)(eps);
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


};

struct BuiltinPW91_LDA : detail::BuiltinKernelImpl< BuiltinPW91_LDA > {

  BuiltinPW91_LDA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPW91_LDA >(p) { }
  
  virtual ~BuiltinPW91_LDA() = default;

};



} // namespace ExchCXX