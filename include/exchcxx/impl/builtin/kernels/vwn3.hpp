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

  static constexpr double dens_tol  = 1e-24;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = sqrt( t10 );
    const double t14 = t10 / 0.4e1 + 0.18637200000000000000e1 * t12 + 0.129352e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.372744e1;
    const double t24 = atan( 0.61519908197590802322e1 / t21 );
    const double t25 = 0.38783294878113014393e-1 * t24;
    const double t27 = t12 / 0.2e1 + 0.10498e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.96902277115443742139e-3 * t30;


    eps = t20 + t25 + t31;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t38 = t4 * t6;
    constexpr double t46 = t3 * t6;
    constexpr double t55 = t1 * t1;
    constexpr double t57 = 0.1e1 / t3;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = sqrt( t10 );
    const double t14 = t10 / 0.4e1 + 0.18637200000000000000e1 * t12 + 0.129352e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.372744e1;
    const double t24 = atan( 0.61519908197590802322e1 / t21 );
    const double t25 = 0.38783294878113014393e-1 * t24;
    const double t27 = t12 / 0.2e1 + 0.10498e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.96902277115443742139e-3 * t30;
    const double t33 = 0.1e1 / t7 / rho;
    const double t34 = t6 * t33;
    const double t39 = t14 * t14;
    const double t40 = 0.1e1 / t39;
    const double t41 = t8 * t40;
    const double t42 = t4 * t34;
    const double t44 = 0.1e1 / t12;
    const double t45 = t44 * t1;
    const double t50 = -t42 / 0.12e2 - 0.31062000000000000000e0 * t45 * t46 * t33;
    const double t58 = ( -t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1 ) * t55 * t57;
    const double t59 = t5 * t7;
    const double t60 = t59 * t14;
    const double t61 = t58 * t60;
    const double t63 = t21 * t21;
    const double t64 = 0.1e1 / t63;
    const double t66 = t64 * t44 * t1;
    const double t68 = 0.37846991046400000000e2 * t64 + 0.1e1;
    const double t69 = 0.1e1 / t68;
    const double t72 = t66 * t46 * t33 * t69;
    const double t74 = t27 * t15;
    const double t75 = t74 * t44;
    const double t78 = t28 * t40;
    const double t80 = -t75 * t42 / 0.6e1 - t78 * t50;
    const double t81 = 0.1e1 / t28;
    const double t82 = t80 * t81;
    const double t83 = t82 * t14;


    eps = t20 + t25 + t31;
    vrho = t20 + t25 + t31 + rho * ( 0.10363566666666666667e-1 * t61 + 0.39765745675026770179e-1 * t72 + 0.96902277115443742139e-3 * t83 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = sqrt( t10 );
    const double t14 = t10 / 0.4e1 + 0.35302100000000000000e1 * t12 + 0.180578e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.1554535e-1 * t19;
    const double t21 = t12 + 0.706042e1;
    const double t24 = atan( 0.47309269095601128300e1 / t21 );
    const double t25 = 0.52491393169780936218e-1 * t24;
    const double t27 = t12 / 0.2e1 + 0.32500e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.22478670955426118383e-2 * t30;


    eps = t20 + t25 + t31;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double& eps, double& vrho ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t38 = t4 * t6;
    constexpr double t46 = t3 * t6;
    constexpr double t55 = t1 * t1;
    constexpr double t57 = 0.1e1 / t3;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = sqrt( t10 );
    const double t14 = t10 / 0.4e1 + 0.35302100000000000000e1 * t12 + 0.180578e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.1554535e-1 * t19;
    const double t21 = t12 + 0.706042e1;
    const double t24 = atan( 0.47309269095601128300e1 / t21 );
    const double t25 = 0.52491393169780936218e-1 * t24;
    const double t27 = t12 / 0.2e1 + 0.32500e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.22478670955426118383e-2 * t30;
    const double t33 = 0.1e1 / t7 / rho;
    const double t34 = t6 * t33;
    const double t39 = t14 * t14;
    const double t40 = 0.1e1 / t39;
    const double t41 = t8 * t40;
    const double t42 = t4 * t34;
    const double t44 = 0.1e1 / t12;
    const double t45 = t44 * t1;
    const double t50 = -t42 / 0.12e2 - 0.58836833333333333333e0 * t45 * t46 * t33;
    const double t58 = ( -t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1 ) * t55 * t57;
    const double t59 = t5 * t7;
    const double t60 = t59 * t14;
    const double t61 = t58 * t60;
    const double t63 = t21 * t21;
    const double t64 = 0.1e1 / t63;
    const double t66 = t64 * t44 * t1;
    const double t68 = 0.22381669423600000000e2 * t64 + 0.1e1;
    const double t69 = 0.1e1 / t68;
    const double t72 = t66 * t46 * t33 * t69;
    const double t74 = t27 * t15;
    const double t75 = t74 * t44;
    const double t78 = t28 * t40;
    const double t80 = -t75 * t42 / 0.6e1 - t78 * t50;
    const double t81 = 0.1e1 / t28;
    const double t82 = t80 * t81;
    const double t83 = t82 * t14;


    eps = t20 + t25 + t31;
    vrho = t20 + t25 + t31 + rho * ( 0.51817833333333333333e-2 * t61 + 0.41388824077869423261e-1 * t72 + 0.22478670955426118383e-2 * t83 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t91 = constants::m_pi_sq;
    constexpr double t123 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t92 = 0.1e1 / t91;
    constexpr double t124 = t123 - 0.1e1;
    constexpr double t126 = 0.1e1 / t124 / 0.2e1;
    constexpr double t135 = 0.9e1 * t124;


    const double t7 = rho_a + rho_b;
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = sqrt( t11 );
    const double t15 = t12 + 0.18637200000000000000e1 * t13 + 0.129352e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = log( t4 * t10 * t16 / 0.4e1 );
    const double t21 = 0.310907e-1 * t20;
    const double t22 = t13 + 0.372744e1;
    const double t25 = atan( 0.61519908197590802322e1 / t22 );
    const double t26 = 0.38783294878113014393e-1 * t25;
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.10498e0;
    const double t29 = t28 * t28;
    const double t31 = log( t29 * t16 );
    const double t32 = 0.96902277115443742139e-3 * t31;
    const double t34 = t12 + 0.35302100000000000000e1 * t13 + 0.180578e2;
    const double t35 = 0.1e1 / t34;
    const double t39 = log( t4 * t10 * t35 / 0.4e1 );
    const double t41 = t13 + 0.706042e1;
    const double t44 = atan( 0.47309269095601128300e1 / t41 );
    const double t46 = t27 + 0.32500e0;
    const double t47 = t46 * t46;
    const double t49 = log( t47 * t35 );
    const double t51 = 0.1554535e-1 * t39 + 0.52491393169780936218e-1 * t44 + 0.22478670955426118383e-2 * t49 - t21 - t26 - t32;
    const double t53 = t12 + 0.10061550000000000000e2 * t13 + 0.101578e3;
    const double t54 = 0.1e1 / t53;
    const double t58 = log( t4 * t10 * t54 / 0.4e1 );
    const double t60 = t13 + 0.201231e2;
    const double t63 = atan( 0.11716852777089929792e1 / t60 );
    const double t65 = t27 + 0.743294e0;
    const double t66 = t65 * t65;
    const double t68 = log( t66 * t54 );
    const double t71 = t12 + 0.65360000000000000000e1 * t13 + 0.427198e2;
    const double t72 = 0.1e1 / t71;
    const double t76 = log( t4 * t10 * t72 / 0.4e1 );
    const double t78 = t13 + 0.130720e2;
    const double t81 = atan( 0.44899888641287296627e-1 / t78 );
    const double t83 = t27 + 0.409286e0;
    const double t84 = t83 * t83;
    const double t86 = log( t84 * t72 );
    const double t88 = 0.1554535e-1 * t58 + 0.61881802979060631482e0 * t63 + 0.26673100072733151594e-2 * t68 - 0.310907e-1 * t76 - 0.20521972937837502661e2 * t81 - 0.44313737677495382697e-2 * t86;
    const double t89 = 0.1e1 / t88;
    const double t90 = t51 * t89;
    const double t94 = t12 + 0.53417500000000000000e0 * t13 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = log( t4 * t10 * t95 / 0.4e1 );
    const double t100 = t13 + 0.106835e1;
    const double t103 = atan( 0.66920720466459414830e1 / t100 );
    const double t105 = t27 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = log( t106 * t95 );
    const double t111 = t92 * ( t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108 );
    const double t112 = t90 * t111;
    const double t113 = rho_a - rho_b;
    const double t114 = 0.1e1 / t7;
    const double t115 = t113 * t114;
    const double t116 = 0.1e1 + t115;
    const double t117 = cbrt( t116 );
    const double t119 = 0.1e1 - t115;
    const double t120 = cbrt( t119 );
    const double t122 = t117 * t116 + t120 * t119 - 0.2e1;
    const double t127 = t122 * t126;
    const double t128 = t113 * t113;
    const double t129 = t128 * t128;
    const double t130 = t7 * t7;
    const double t131 = t130 * t130;
    const double t132 = 0.1e1 / t131;
    const double t134 = -t129 * t132 + 0.1e1;
    const double t136 = t134 * t135;
    const double t137 = t127 * t136;
    const double t139 = t112 * t137 / 0.24e2;
    const double t140 = t51 * t122;
    const double t141 = t126 * t129;
    const double t142 = t141 * t132;
    const double t143 = t140 * t142;


    eps = t21 + t26 + t32 - t139 + t143;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t91 = constants::m_pi_sq;
    constexpr double t123 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t92 = 0.1e1 / t91;
    constexpr double t124 = t123 - 0.1e1;
    constexpr double t126 = 0.1e1 / t124 / 0.2e1;
    constexpr double t135 = 0.9e1 * t124;
    constexpr double t150 = t4 * t6;
    constexpr double t158 = t3 * t6;
    constexpr double t167 = t1 * t1;
    constexpr double t169 = 0.1e1 / t3;


    const double t7 = rho_a + rho_b;
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = sqrt( t11 );
    const double t15 = t12 + 0.18637200000000000000e1 * t13 + 0.129352e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = log( t4 * t10 * t16 / 0.4e1 );
    const double t21 = 0.310907e-1 * t20;
    const double t22 = t13 + 0.372744e1;
    const double t25 = atan( 0.61519908197590802322e1 / t22 );
    const double t26 = 0.38783294878113014393e-1 * t25;
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.10498e0;
    const double t29 = t28 * t28;
    const double t31 = log( t29 * t16 );
    const double t32 = 0.96902277115443742139e-3 * t31;
    const double t34 = t12 + 0.35302100000000000000e1 * t13 + 0.180578e2;
    const double t35 = 0.1e1 / t34;
    const double t39 = log( t4 * t10 * t35 / 0.4e1 );
    const double t41 = t13 + 0.706042e1;
    const double t44 = atan( 0.47309269095601128300e1 / t41 );
    const double t46 = t27 + 0.32500e0;
    const double t47 = t46 * t46;
    const double t49 = log( t47 * t35 );
    const double t51 = 0.1554535e-1 * t39 + 0.52491393169780936218e-1 * t44 + 0.22478670955426118383e-2 * t49 - t21 - t26 - t32;
    const double t53 = t12 + 0.10061550000000000000e2 * t13 + 0.101578e3;
    const double t54 = 0.1e1 / t53;
    const double t58 = log( t4 * t10 * t54 / 0.4e1 );
    const double t60 = t13 + 0.201231e2;
    const double t63 = atan( 0.11716852777089929792e1 / t60 );
    const double t65 = t27 + 0.743294e0;
    const double t66 = t65 * t65;
    const double t68 = log( t66 * t54 );
    const double t71 = t12 + 0.65360000000000000000e1 * t13 + 0.427198e2;
    const double t72 = 0.1e1 / t71;
    const double t76 = log( t4 * t10 * t72 / 0.4e1 );
    const double t78 = t13 + 0.130720e2;
    const double t81 = atan( 0.44899888641287296627e-1 / t78 );
    const double t83 = t27 + 0.409286e0;
    const double t84 = t83 * t83;
    const double t86 = log( t84 * t72 );
    const double t88 = 0.1554535e-1 * t58 + 0.61881802979060631482e0 * t63 + 0.26673100072733151594e-2 * t68 - 0.310907e-1 * t76 - 0.20521972937837502661e2 * t81 - 0.44313737677495382697e-2 * t86;
    const double t89 = 0.1e1 / t88;
    const double t90 = t51 * t89;
    const double t94 = t12 + 0.53417500000000000000e0 * t13 + 0.114813e2;
    const double t95 = 0.1e1 / t94;
    const double t99 = log( t4 * t10 * t95 / 0.4e1 );
    const double t100 = t13 + 0.106835e1;
    const double t103 = atan( 0.66920720466459414830e1 / t100 );
    const double t105 = t27 + 0.228344e0;
    const double t106 = t105 * t105;
    const double t108 = log( t106 * t95 );
    const double t111 = t92 * ( t99 + 0.32323836906055067299e0 * t103 + 0.21608710360898267022e-1 * t108 );
    const double t112 = t90 * t111;
    const double t113 = rho_a - rho_b;
    const double t114 = 0.1e1 / t7;
    const double t115 = t113 * t114;
    const double t116 = 0.1e1 + t115;
    const double t117 = cbrt( t116 );
    const double t119 = 0.1e1 - t115;
    const double t120 = cbrt( t119 );
    const double t122 = t117 * t116 + t120 * t119 - 0.2e1;
    const double t127 = t122 * t126;
    const double t128 = t113 * t113;
    const double t129 = t128 * t128;
    const double t130 = t7 * t7;
    const double t131 = t130 * t130;
    const double t132 = 0.1e1 / t131;
    const double t134 = -t129 * t132 + 0.1e1;
    const double t136 = t134 * t135;
    const double t137 = t127 * t136;
    const double t139 = t112 * t137 / 0.24e2;
    const double t140 = t51 * t122;
    const double t141 = t126 * t129;
    const double t142 = t141 * t132;
    const double t143 = t140 * t142;
    const double t145 = 0.1e1 / t8 / t7;
    const double t146 = t6 * t145;
    const double t151 = t15 * t15;
    const double t152 = 0.1e1 / t151;
    const double t153 = t9 * t152;
    const double t154 = t4 * t146;
    const double t155 = t154 / 0.12e2;
    const double t156 = 0.1e1 / t13;
    const double t157 = t156 * t1;
    const double t160 = t157 * t158 * t145;
    const double t162 = -t155 - 0.31062000000000000000e0 * t160;
    const double t170 = ( -t4 * t146 * t16 / 0.12e2 - t150 * t153 * t162 / 0.4e1 ) * t167 * t169;
    const double t171 = t5 * t8;
    const double t172 = t171 * t15;
    const double t173 = t170 * t172;
    const double t174 = 0.10363566666666666667e-1 * t173;
    const double t175 = t22 * t22;
    const double t176 = 0.1e1 / t175;
    const double t178 = t176 * t156 * t1;
    const double t180 = 0.37846991046400000000e2 * t176 + 0.1e1;
    const double t181 = 0.1e1 / t180;
    const double t184 = t178 * t158 * t145 * t181;
    const double t185 = 0.39765745675026770179e-1 * t184;
    const double t186 = t28 * t16;
    const double t187 = t186 * t156;
    const double t190 = t29 * t152;
    const double t192 = -t187 * t154 / 0.6e1 - t190 * t162;
    const double t193 = 0.1e1 / t29;
    const double t194 = t192 * t193;
    const double t195 = t194 * t15;
    const double t196 = 0.96902277115443742139e-3 * t195;
    const double t200 = t34 * t34;
    const double t201 = 0.1e1 / t200;
    const double t202 = t9 * t201;
    const double t204 = -t155 - 0.58836833333333333333e0 * t160;
    const double t210 = ( -t4 * t146 * t35 / 0.12e2 - t150 * t202 * t204 / 0.4e1 ) * t167 * t169;
    const double t211 = t171 * t34;
    const double t214 = t41 * t41;
    const double t215 = 0.1e1 / t214;
    const double t217 = t215 * t156 * t1;
    const double t219 = 0.22381669423600000000e2 * t215 + 0.1e1;
    const double t220 = 0.1e1 / t219;
    const double t225 = t46 * t35;
    const double t226 = t225 * t156;
    const double t229 = t47 * t201;
    const double t231 = -t226 * t154 / 0.6e1 - t229 * t204;
    const double t232 = 0.1e1 / t47;
    const double t233 = t231 * t232;
    const double t236 = 0.51817833333333333333e-2 * t210 * t211 + 0.41388824077869423261e-1 * t217 * t158 * t145 * t220 + 0.22478670955426118383e-2 * t233 * t34 - t174 - t185 - t196;
    const double t237 = t236 * t89;
    const double t238 = t237 * t111;
    const double t239 = t238 * t137;
    const double t240 = t239 / 0.24e2;
    const double t241 = t88 * t88;
    const double t242 = 0.1e1 / t241;
    const double t243 = t51 * t242;
    const double t244 = t243 * t111;
    const double t248 = t53 * t53;
    const double t249 = 0.1e1 / t248;
    const double t250 = t9 * t249;
    const double t252 = -t155 - 0.16769250000000000000e1 * t160;
    const double t258 = ( -t4 * t146 * t54 / 0.12e2 - t150 * t250 * t252 / 0.4e1 ) * t167 * t169;
    const double t259 = t171 * t53;
    const double t262 = t60 * t60;
    const double t263 = 0.1e1 / t262;
    const double t265 = t263 * t156 * t1;
    const double t267 = 0.13728463900000000000e1 * t263 + 0.1e1;
    const double t268 = 0.1e1 / t267;
    const double t273 = t65 * t54;
    const double t274 = t273 * t156;
    const double t277 = t66 * t249;
    const double t279 = -t274 * t154 / 0.6e1 - t277 * t252;
    const double t280 = 0.1e1 / t66;
    const double t281 = t279 * t280;
    const double t287 = t71 * t71;
    const double t288 = 0.1e1 / t287;
    const double t289 = t9 * t288;
    const double t291 = -t155 - 0.10893333333333333333e1 * t160;
    const double t297 = ( -t4 * t146 * t72 / 0.12e2 - t150 * t289 * t291 / 0.4e1 ) * t167 * t169;
    const double t298 = t171 * t71;
    const double t301 = t78 * t78;
    const double t302 = 0.1e1 / t301;
    const double t304 = t302 * t156 * t1;
    const double t306 = 0.20160000000000000000e-2 * t302 + 0.1e1;
    const double t307 = 0.1e1 / t306;
    const double t312 = t83 * t72;
    const double t313 = t312 * t156;
    const double t316 = t84 * t288;
    const double t318 = -t313 * t154 / 0.6e1 - t316 * t291;
    const double t319 = 0.1e1 / t84;
    const double t320 = t318 * t319;
    const double t323 = 0.51817833333333333333e-2 * t258 * t259 + 0.12084332918108974175e0 * t265 * t158 * t145 * t268 + 0.26673100072733151594e-2 * t281 * t53 - 0.10363566666666666667e-1 * t297 * t298 - 0.15357238326806922974e0 * t304 * t158 * t145 * t307 - 0.44313737677495382697e-2 * t320 * t71;
    const double t324 = t136 * t323;
    const double t325 = t127 * t324;
    const double t326 = t244 * t325;
    const double t327 = t326 / 0.24e2;
    const double t331 = t94 * t94;
    const double t332 = 0.1e1 / t331;
    const double t333 = t9 * t332;
    const double t335 = -t155 - 0.89029166666666666667e-1 * t160;
    const double t341 = ( -t4 * t146 * t95 / 0.12e2 - t150 * t333 * t335 / 0.4e1 ) * t167 * t169;
    const double t342 = t171 * t94;
    const double t345 = t100 * t100;
    const double t346 = 0.1e1 / t345;
    const double t348 = t346 * t156 * t1;
    const double t350 = 0.44783828277500000000e2 * t346 + 0.1e1;
    const double t351 = 0.1e1 / t350;
    const double t356 = t105 * t95;
    const double t357 = t356 * t156;
    const double t360 = t106 * t332;
    const double t362 = -t357 * t154 / 0.6e1 - t360 * t335;
    const double t363 = 0.1e1 / t106;
    const double t364 = t362 * t363;
    const double t368 = t92 * ( t341 * t342 / 0.3e1 + 0.36052240899892258525e0 * t348 * t158 * t145 * t351 + 0.21608710360898267022e-1 * t364 * t94 );
    const double t369 = t90 * t368;
    const double t370 = t369 * t137;
    const double t371 = t370 / 0.24e2;
    const double t372 = 0.1e1 / t130;
    const double t373 = t113 * t372;
    const double t374 = t114 - t373;
    const double t376 = -t374;
    const double t379 = 0.4e1 / 0.3e1 * t117 * t374 + 0.4e1 / 0.3e1 * t120 * t376;
    const double t380 = t379 * t126;
    const double t381 = t380 * t136;
    const double t382 = t112 * t381;
    const double t383 = t382 / 0.24e2;
    const double t384 = t128 * t113;
    const double t385 = t384 * t132;
    const double t386 = t131 * t7;
    const double t387 = 0.1e1 / t386;
    const double t388 = t129 * t387;
    const double t390 = -0.4e1 * t385 + 0.4e1 * t388;
    const double t391 = t390 * t135;
    const double t392 = t127 * t391;
    const double t393 = t112 * t392;
    const double t394 = t393 / 0.24e2;
    const double t395 = t236 * t122;
    const double t396 = t395 * t142;
    const double t397 = t51 * t379;
    const double t398 = t397 * t142;
    const double t399 = t126 * t384;
    const double t400 = t399 * t132;
    const double t401 = t140 * t400;
    const double t402 = 0.4e1 * t401;
    const double t403 = t141 * t387;
    const double t404 = t140 * t403;
    const double t405 = 0.4e1 * t404;
    const double t406 = t174 + t185 + t196 - t240 + t327 - t371 - t383 - t394 + t396 + t398 + t402 - t405;
    const double t408 = -t114 - t373;
    const double t410 = -t408;
    const double t413 = 0.4e1 / 0.3e1 * t117 * t408 + 0.4e1 / 0.3e1 * t120 * t410;
    const double t414 = t413 * t126;
    const double t415 = t414 * t136;
    const double t416 = t112 * t415;
    const double t417 = t416 / 0.24e2;
    const double t419 = 0.4e1 * t385 + 0.4e1 * t388;
    const double t420 = t419 * t135;
    const double t421 = t127 * t420;
    const double t422 = t112 * t421;
    const double t423 = t422 / 0.24e2;
    const double t424 = t51 * t413;
    const double t425 = t424 * t142;
    const double t426 = t174 + t185 + t196 - t240 + t327 - t371 - t417 - t423 + t396 + t425 - t402 - t405;


    eps = t21 + t26 + t32 - t139 + t143;
    vrho_a = t7 * t406 - t139 + t143 + t21 + t26 + t32;
    vrho_b = t7 * t426 - t139 + t143 + t21 + t26 + t32;

  }


};

struct BuiltinVWN3 : detail::BuiltinKernelImpl< BuiltinVWN3 > {

  BuiltinVWN3( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN3 >(p) { }
  
  virtual ~BuiltinVWN3() = default;

};



} // namespace ExchCXX