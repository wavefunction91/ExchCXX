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

  static constexpr double dens_tol  = 1e-24;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

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

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
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


    const double t10 = cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t13 = t4 * t12;
    const double t15 = 0.1e1 + t13 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t37 = pow( t21 / 0.4e1, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.12500000000000000000e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = log( t43 );


    eps = -0.2e1 * t1 * t15 * t44;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
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


    const double t10 = cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t13 = t4 * t12;
    const double t15 = 0.1e1 + t13 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t37 = pow( t21 / 0.4e1, t36 );
    const double t38 = beta4_0 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.12500000000000000000e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = log( t43 );
    const double t53 = rho * t15;
    const double t54 = t39 * t39;
    const double t55 = 0.1e1 / t54;
    const double t58 = t18 / t22 * t3;
    const double t60 = 0.1e1 / t10 / rho;
    const double t61 = t9 * t60;
    const double t66 = sqrt( t21 );
    const double t68 = t29 * t66 * t3;
    const double t71 = 0.1e1 / rho;
    const double t75 = -t58 * t61 / 0.12e2 - t26 * t61 / 0.12e2 - 0.62500000000000000000e-1 * t68 * t61 - t38 * t36 * t71 / 0.3e1;
    const double t77 = 0.1e1 / t43;


    eps = -0.2e1 * t1 * t15 * t44;
    vrho = ( -0.2e1 * t1 * t15 * t44 ) + t11 * t1 * t2 * t19 * t8 * t44 / 0.6e1 + t53 * t55 * t75 * t77;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t1 = a_1;
    constexpr double t2 = alpha1_1;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_1;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_1 * t3;
    constexpr double t29 = beta3_1;
    constexpr double t36 = pp_1 + 0.1e1;


    const double t10 = cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t13 = t4 * t12;
    const double t15 = 0.1e1 + t13 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t37 = pow( t21 / 0.4e1, t36 );
    const double t38 = beta4_1 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.12500000000000000000e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = log( t43 );


    eps = -0.2e1 * t1 * t15 * t44;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t1 = a_1;
    constexpr double t2 = alpha1_1;
    constexpr double t4 = t2 * t3;
    constexpr double t8 = t7 * t7;
    constexpr double t9 = t6 * t8;
    constexpr double t17 = 0.1e1 / t1;
    constexpr double t18 = beta1_1;
    constexpr double t19 = t3 * t6;
    constexpr double t26 = beta2_1 * t3;
    constexpr double t29 = beta3_1;
    constexpr double t36 = pp_1 + 0.1e1;


    const double t10 = cbrt( rho );
    const double t11 = 0.1e1 / t10;
    const double t12 = t9 * t11;
    const double t13 = t4 * t12;
    const double t15 = 0.1e1 + t13 / 0.4e1;
    const double t21 = t19 * t8 * t11;
    const double t22 = sqrt( t21 );
    const double t30 = pow_3_2( t21 );
    const double t37 = pow( t21 / 0.4e1, t36 );
    const double t38 = beta4_1 * t37;
    const double t39 = t18 * t22 / 0.2e1 + t26 * t12 / 0.4e1 + 0.12500000000000000000e0 * t29 * t30 + t38;
    const double t43 = 0.1e1 + t17 / t39 / 0.2e1;
    const double t44 = log( t43 );
    const double t53 = rho * t15;
    const double t54 = t39 * t39;
    const double t55 = 0.1e1 / t54;
    const double t58 = t18 / t22 * t3;
    const double t60 = 0.1e1 / t10 / rho;
    const double t61 = t9 * t60;
    const double t66 = sqrt( t21 );
    const double t68 = t29 * t66 * t3;
    const double t71 = 0.1e1 / rho;
    const double t75 = -t58 * t61 / 0.12e2 - t26 * t61 / 0.12e2 - 0.62500000000000000000e-1 * t68 * t61 - t38 * t36 * t71 / 0.3e1;
    const double t77 = 0.1e1 / t43;


    eps = -0.2e1 * t1 * t15 * t44;
    vrho = ( -0.2e1 * t1 * t15 * t44 ) + t11 * t1 * t2 * t19 * t8 * t44 / 0.6e1 + t53 * t55 * t75 * t77;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t64 = constants::m_cbrt_2;
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
    constexpr double t69 = a_1;
    constexpr double t70 = alpha1_1;
    constexpr double t71 = t70 * t3;
    constexpr double t76 = 0.1e1 / t69;
    constexpr double t77 = beta1_1;
    constexpr double t81 = beta2_1 * t3;
    constexpr double t84 = beta3_1;
    constexpr double t89 = pp_1 + 0.1e1;
    constexpr double t99 = a_2;
    constexpr double t100 = alpha1_2;
    constexpr double t101 = t100 * t3;
    constexpr double t106 = 0.1e1 / t99;
    constexpr double t107 = beta1_2;
    constexpr double t111 = beta2_2 * t3;
    constexpr double t114 = beta3_2;
    constexpr double t119 = pp_2 + 0.1e1;
    constexpr double t128 = 0.1e1 / fz20;


    const double t10 = rho_a + rho_b;
    const double t11 = cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.12500000000000000000e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = log( t44 );
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
    const double t58 = cbrt( t57 );
    const double t60 = 0.1e1 - t56;
    const double t61 = cbrt( t60 );
    const double t63 = t58 * t57 + t61 * t60 - 0.2e1;
    const double t67 = 0.1e1 / ( 0.2e1 * t64 - 0.2e1 );
    const double t68 = t63 * t67;
    const double t74 = 0.1e1 + t71 * t13 / 0.4e1;
    const double t90 = pow( t35, t89 );
    const double t91 = beta4_1 * t90;
    const double t92 = t77 * t23 / 0.2e1 + t81 * t13 / 0.4e1 + 0.12500000000000000000e0 * t84 * t31 + t91;
    const double t96 = 0.1e1 + t76 / t92 / 0.2e1;
    const double t97 = log( t96 );
    const double t104 = 0.1e1 + t101 * t13 / 0.4e1;
    const double t120 = pow( t35, t119 );
    const double t121 = beta4_2 * t120;
    const double t122 = t107 * t23 / 0.2e1 + t111 * t13 / 0.4e1 + 0.12500000000000000000e0 * t114 * t31 + t121;
    const double t126 = 0.1e1 + t106 / t122 / 0.2e1;
    const double t127 = log( t126 );
    const double t129 = t127 * t128;
    const double t132 = -0.2e1 * t99 * t104 * t129 - 0.2e1 * t69 * t74 * t97 + 0.2e1 * t46;
    const double t133 = t68 * t132;
    const double t134 = t54 * t133;
    const double t137 = t104 * t127 * t128;
    const double t139 = 0.2e1 * t68 * t99 * t137;


    eps = -t47 + t134 + t139;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_one_ov_pi;
    constexpr double t7 = constants::m_cbrt_4;
    constexpr double t64 = constants::m_cbrt_2;
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
    constexpr double t69 = a_1;
    constexpr double t70 = alpha1_1;
    constexpr double t71 = t70 * t3;
    constexpr double t76 = 0.1e1 / t69;
    constexpr double t77 = beta1_1;
    constexpr double t81 = beta2_1 * t3;
    constexpr double t84 = beta3_1;
    constexpr double t89 = pp_1 + 0.1e1;
    constexpr double t99 = a_2;
    constexpr double t100 = alpha1_2;
    constexpr double t101 = t100 * t3;
    constexpr double t106 = 0.1e1 / t99;
    constexpr double t107 = beta1_2;
    constexpr double t111 = beta2_2 * t3;
    constexpr double t114 = beta3_2;
    constexpr double t119 = pp_2 + 0.1e1;
    constexpr double t128 = 0.1e1 / fz20;
    constexpr double t141 = t1 * t2 * t3;
    constexpr double t192 = t69 * t70 * t3;
    constexpr double t217 = t99 * t100;
    constexpr double t218 = t217 * t20;
    constexpr double t250 = t217 * t3;


    const double t10 = rho_a + rho_b;
    const double t11 = cbrt( t10 );
    const double t12 = 0.1e1 / t11;
    const double t13 = t9 * t12;
    const double t16 = 0.1e1 + t4 * t13 / 0.4e1;
    const double t22 = t20 * t8 * t12;
    const double t23 = sqrt( t22 );
    const double t31 = pow_3_2( t22 );
    const double t35 = t22 / 0.4e1;
    const double t38 = pow( t35, t37 );
    const double t39 = beta4_0 * t38;
    const double t40 = t19 * t23 / 0.2e1 + t27 * t13 / 0.4e1 + 0.12500000000000000000e0 * t30 * t31 + t39;
    const double t44 = 0.1e1 + t18 / t40 / 0.2e1;
    const double t45 = log( t44 );
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
    const double t58 = cbrt( t57 );
    const double t60 = 0.1e1 - t56;
    const double t61 = cbrt( t60 );
    const double t63 = t58 * t57 + t61 * t60 - 0.2e1;
    const double t67 = 0.1e1 / ( 0.2e1 * t64 - 0.2e1 );
    const double t68 = t63 * t67;
    const double t74 = 0.1e1 + t71 * t13 / 0.4e1;
    const double t90 = pow( t35, t89 );
    const double t91 = beta4_1 * t90;
    const double t92 = t77 * t23 / 0.2e1 + t81 * t13 / 0.4e1 + 0.12500000000000000000e0 * t84 * t31 + t91;
    const double t96 = 0.1e1 + t76 / t92 / 0.2e1;
    const double t97 = log( t96 );
    const double t104 = 0.1e1 + t101 * t13 / 0.4e1;
    const double t120 = pow( t35, t119 );
    const double t121 = beta4_2 * t120;
    const double t122 = t107 * t23 / 0.2e1 + t111 * t13 / 0.4e1 + 0.12500000000000000000e0 * t114 * t31 + t121;
    const double t126 = 0.1e1 + t106 / t122 / 0.2e1;
    const double t127 = log( t126 );
    const double t129 = t127 * t128;
    const double t132 = -0.2e1 * t99 * t104 * t129 - 0.2e1 * t69 * t74 * t97 + 0.2e1 * t46;
    const double t133 = t68 * t132;
    const double t134 = t54 * t133;
    const double t137 = t104 * t127 * t128;
    const double t139 = 0.2e1 * t68 * t99 * t137;
    const double t143 = 0.1e1 / t11 / t10;
    const double t146 = t141 * t9 * t143 * t45;
    const double t147 = t146 / 0.6e1;
    const double t148 = t40 * t40;
    const double t149 = 0.1e1 / t148;
    const double t150 = t16 * t149;
    const double t151 = 0.1e1 / t23;
    const double t153 = t19 * t151 * t3;
    const double t154 = t9 * t143;
    const double t159 = sqrt( t22 );
    const double t161 = t30 * t159 * t3;
    const double t167 = -t153 * t154 / 0.12e2 - t27 * t154 / 0.12e2 - 0.62500000000000000000e-1 * t161 * t154 - t39 * t37 * t55 / 0.3e1;
    const double t168 = 0.1e1 / t44;
    const double t169 = t167 * t168;
    const double t170 = t150 * t169;
    const double t171 = t49 * t48;
    const double t172 = t171 * t53;
    const double t173 = t172 * t133;
    const double t174 = 0.4e1 * t173;
    const double t175 = t52 * t10;
    const double t176 = 0.1e1 / t175;
    const double t177 = t50 * t176;
    const double t178 = t177 * t133;
    const double t179 = 0.4e1 * t178;
    const double t180 = 0.1e1 / t51;
    const double t181 = t48 * t180;
    const double t182 = t55 - t181;
    const double t184 = -t182;
    const double t188 = ( 0.4e1 / 0.3e1 * t58 * t182 + 0.4e1 / 0.3e1 * t61 * t184 ) * t67;
    const double t189 = t188 * t132;
    const double t190 = t54 * t189;
    const double t197 = t92 * t92;
    const double t198 = 0.1e1 / t197;
    const double t199 = t74 * t198;
    const double t201 = t77 * t151 * t3;
    const double t207 = t84 * t159 * t3;
    const double t213 = -t201 * t154 / 0.12e2 - t81 * t154 / 0.12e2 - 0.62500000000000000000e-1 * t207 * t154 - t91 * t89 * t55 / 0.3e1;
    const double t214 = 0.1e1 / t96;
    const double t215 = t213 * t214;
    const double t219 = t8 * t143;
    const double t223 = t122 * t122;
    const double t224 = 0.1e1 / t223;
    const double t225 = t104 * t224;
    const double t227 = t107 * t151 * t3;
    const double t233 = t114 * t159 * t3;
    const double t239 = -t227 * t154 / 0.12e2 - t111 * t154 / 0.12e2 - 0.62500000000000000000e-1 * t233 * t154 - t121 * t119 * t55 / 0.3e1;
    const double t240 = 0.1e1 / t126;
    const double t242 = t239 * t240 * t128;
    const double t244 = t192 * t9 * t143 * t97 / 0.6e1 + t199 * t215 - t147 - t170 + t218 * t219 * t129 / 0.6e1 + t225 * t242;
    const double t245 = t68 * t244;
    const double t246 = t54 * t245;
    const double t248 = t188 * t99 * t137;
    const double t249 = 0.2e1 * t248;
    const double t251 = t68 * t250;
    const double t254 = t9 * t143 * t127 * t128;
    const double t255 = t251 * t254;
    const double t256 = t255 / 0.6e1;
    const double t257 = t68 * t104;
    const double t259 = t240 * t128;
    const double t260 = t224 * t239 * t259;
    const double t261 = t257 * t260;
    const double t264 = -t55 - t181;
    const double t266 = -t264;
    const double t270 = ( 0.4e1 / 0.3e1 * t58 * t264 + 0.4e1 / 0.3e1 * t61 * t266 ) * t67;
    const double t271 = t270 * t132;
    const double t272 = t54 * t271;
    const double t274 = t270 * t99 * t137;
    const double t275 = 0.2e1 * t274;


    eps = -t47 + t134 + t139;
    vrho_a = -t47 + t134 + t139 + t10 * ( t147 + t170 + t174 - t179 + t190 + t246 + t249 - t256 - t261 );
    vrho_b = -t47 + t134 + t139 + t10 * ( t147 + t170 - t174 - t179 + t272 + t246 + t275 - t256 - t261 );

  }


};

struct BuiltinPW91_LDA_RPA : detail::BuiltinKernelImpl< BuiltinPW91_LDA_RPA > {

  BuiltinPW91_LDA_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPW91_LDA_RPA >(p) { }
  
  virtual ~BuiltinPW91_LDA_RPA() = default;

};



} // namespace ExchCXX