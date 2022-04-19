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

  static constexpr double dens_tol  = 1e-12;

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
    constexpr double t37 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t41 = t18 / t3 * t5;
    constexpr double t44 = BB * beta;
    constexpr double t45 = 0.1e1 / gamma;
    constexpr double t58 = t37 * t37;
    constexpr double t60 = 0.1e1 / t19;
    constexpr double t61 = t1 * t60;
    constexpr double t62 = t61 * t6;
    constexpr double t68 = beta * t45;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.53425000000000000000e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.37978500000000000000e1 * t13 + 0.89690000000000000000e0 * t10 + 0.20477500000000000000e0 * t16 + 0.12323500000000000000e0 * t24;
    const double t29 = 0.1e1 + 0.16081979498692535067e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t31 = t12 * t30;
    const double t32 = 0.621814e-1 * t31;
    const double t33 = rho * rho;
    const double t35 = 0.1e1 / t7 / t33;
    const double t48 = safe_math::exp( 0.621814e-1 * t31 * t45 );
    const double t49 = t48 - 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t51 = t45 * t50;
    const double t52 = sigma * sigma;
    const double t54 = t44 * t51 * t52;
    const double t55 = t33 * t33;
    const double t57 = 0.1e1 / t21 / t55;
    const double t59 = t57 * t58;
    const double t63 = t59 * t62;
    const double t66 = sigma * t35 * t37 * t41 / 0.96e2 + t54 * t63 / 0.3072e4;
    const double t67 = beta * t66;
    const double t71 = t68 * t50 * t66 + 0.1e1;
    const double t72 = 0.1e1 / t71;
    const double t73 = t45 * t72;
    const double t75 = t67 * t73 + 0.1e1;
    const double t76 = safe_math::log( t75 );
    const double t77 = gamma * t76;


    eps = -t32 + t77;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t37 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t41 = t18 / t3 * t5;
    constexpr double t44 = BB * beta;
    constexpr double t45 = 0.1e1 / gamma;
    constexpr double t58 = t37 * t37;
    constexpr double t60 = 0.1e1 / t19;
    constexpr double t61 = t1 * t60;
    constexpr double t62 = t61 * t6;
    constexpr double t68 = beta * t45;
    constexpr double t89 = t3 * t6;
    constexpr double t116 = t44 * t45;
    constexpr double t122 = t58 * t1;
    constexpr double t123 = t122 * t60;
    constexpr double t124 = t4 * t6;
    constexpr double t177 = beta * beta;
    constexpr double t179 = gamma * gamma;
    constexpr double t180 = 0.1e1 / t179;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.53425000000000000000e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.37978500000000000000e1 * t13 + 0.89690000000000000000e0 * t10 + 0.20477500000000000000e0 * t16 + 0.12323500000000000000e0 * t24;
    const double t29 = 0.1e1 + 0.16081979498692535067e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t31 = t12 * t30;
    const double t32 = 0.621814e-1 * t31;
    const double t33 = rho * rho;
    const double t35 = 0.1e1 / t7 / t33;
    const double t48 = safe_math::exp( 0.621814e-1 * t31 * t45 );
    const double t49 = t48 - 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t51 = t45 * t50;
    const double t52 = sigma * sigma;
    const double t54 = t44 * t51 * t52;
    const double t55 = t33 * t33;
    const double t57 = 0.1e1 / t21 / t55;
    const double t59 = t57 * t58;
    const double t63 = t59 * t62;
    const double t66 = sigma * t35 * t37 * t41 / 0.96e2 + t54 * t63 / 0.3072e4;
    const double t67 = beta * t66;
    const double t71 = t68 * t50 * t66 + 0.1e1;
    const double t72 = 0.1e1 / t71;
    const double t73 = t45 * t72;
    const double t75 = t67 * t73 + 0.1e1;
    const double t76 = safe_math::log( t75 );
    const double t77 = gamma * t76;
    const double t79 = 0.1e1 / t7 / rho;
    const double t80 = t6 * t79;
    const double t82 = t4 * t80 * t30;
    const double t84 = t26 * t26;
    const double t85 = 0.1e1 / t84;
    const double t86 = t12 * t85;
    const double t88 = 0.1e1 / t13 * t1;
    const double t90 = t89 * t79;
    const double t93 = t4 * t80;
    const double t95 = safe_math::sqrt( t10 );
    const double t96 = t95 * t1;
    const double t104 = -0.63297500000000000000e0 * t88 * t90 - 0.29896666666666666667e0 * t93 - 0.10238750000000000000e0 * t96 * t90 - 0.82156666666666666667e-1 * t20 * t5 / t21 / rho;
    const double t105 = 0.1e1 / t29;
    const double t106 = t104 * t105;
    const double t107 = t86 * t106;
    const double t109 = t33 * rho;
    const double t111 = 0.1e1 / t7 / t109;
    const double t117 = t49 * t49;
    const double t118 = 0.1e1 / t117;
    const double t119 = t118 * t52;
    const double t121 = t116 * t119 * t57;
    const double t132 = -0.11073470983333333333e-2 * t124 * t79 * t30 * t45 - 0.10000000000000000000e1 * t86 * t106 * t45;
    const double t133 = t6 * t132;
    const double t135 = t123 * t133 * t48;
    const double t138 = t55 * rho;
    const double t140 = 0.1e1 / t21 / t138;
    const double t141 = t140 * t58;
    const double t142 = t141 * t62;
    const double t145 = -0.7e1 / 0.288e3 * sigma * t111 * t37 * t41 - t121 * t135 / 0.3072e4 - 0.7e1 / 0.4608e4 * t54 * t142;
    const double t146 = beta * t145;
    const double t148 = t71 * t71;
    const double t149 = 0.1e1 / t148;
    const double t150 = t45 * t149;
    const double t151 = t68 * t118;
    const double t152 = t66 * t132;
    const double t157 = t68 * t50 * t145 - t151 * t152 * t48;
    const double t158 = t150 * t157;
    const double t160 = t146 * t73 - t67 * t158;
    const double t162 = 0.1e1 / t75;
    const double t163 = gamma * t160 * t162;
    const double t166 = rho * gamma;
    const double t171 = t44 * t51 * sigma;
    const double t174 = t35 * t37 * t41 / 0.96e2 + t171 * t63 / 0.1536e4;
    const double t175 = beta * t174;
    const double t178 = t177 * t66;
    const double t181 = t178 * t180;
    const double t182 = t149 * t50;
    const double t183 = t182 * t174;
    const double t185 = t175 * t73 - t181 * t183;


    eps = -t32 + t77;
    vrho = -t32 + t77 + rho * ( 0.11073470983333333333e-2 * t82 + 0.10000000000000000000e1 * t107 + t163 );
    vsigma = t166 * t185 * t162;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t37 = 0.1e1 / t3;
    constexpr double t39 = t18 * t37 * t5;
    constexpr double t42 = BB * beta;
    constexpr double t43 = 0.1e1 / gamma;
    constexpr double t56 = 0.1e1 / t19;
    constexpr double t57 = t1 * t56;
    constexpr double t58 = t57 * t6;
    constexpr double t64 = beta * t43;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.51370000000000000000e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.70594500000000000000e1 * t13 + 0.15494250000000000000e1 * t10 + 0.42077500000000000000e0 * t16 + 0.15629250000000000000e0 * t24;
    const double t29 = 0.1e1 + 0.32163958997385070134e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t31 = t12 * t30;
    const double t32 = 0.3109070e-1 * t31;
    const double t33 = rho * rho;
    const double t35 = 0.1e1 / t7 / t33;
    const double t46 = safe_math::exp( 0.62181400000000000000e-1 * t31 * t43 );
    const double t47 = t46 - 0.1e1;
    const double t48 = 0.1e1 / t47;
    const double t50 = t42 * t43 * t48;
    const double t51 = sigma * sigma;
    const double t52 = t33 * t33;
    const double t54 = 0.1e1 / t21 / t52;
    const double t62 = sigma * t35 * t39 / 0.48e2 + t50 * t51 * t54 * t58 / 0.768e3;
    const double t63 = beta * t62;
    const double t67 = t64 * t48 * t62 + 0.1e1;
    const double t68 = 0.1e1 / t67;
    const double t69 = t43 * t68;
    const double t71 = t63 * t69 + 0.1e1;
    const double t72 = safe_math::log( t71 );
    const double t74 = gamma * t72 / 0.2e1;


    eps = -t32 + t74;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = t3 * t3;
    constexpr double t20 = t18 * t19;
    constexpr double t37 = 0.1e1 / t3;
    constexpr double t39 = t18 * t37 * t5;
    constexpr double t42 = BB * beta;
    constexpr double t43 = 0.1e1 / gamma;
    constexpr double t56 = 0.1e1 / t19;
    constexpr double t57 = t1 * t56;
    constexpr double t58 = t57 * t6;
    constexpr double t64 = beta * t43;
    constexpr double t86 = t3 * t6;
    constexpr double t119 = t4 * t6;
    constexpr double t164 = t37 * t5;
    constexpr double t174 = beta * beta;
    constexpr double t176 = gamma * gamma;
    constexpr double t177 = 0.1e1 / t176;


    const double t7 = safe_math::cbrt( rho );
    const double t10 = t4 * t6 / t7;
    const double t12 = 0.1e1 + 0.51370000000000000000e-1 * t10;
    const double t13 = safe_math::sqrt( t10 );
    const double t16 = pow_3_2( t10 );
    const double t21 = t7 * t7;
    const double t24 = t20 * t5 / t21;
    const double t26 = 0.70594500000000000000e1 * t13 + 0.15494250000000000000e1 * t10 + 0.42077500000000000000e0 * t16 + 0.15629250000000000000e0 * t24;
    const double t29 = 0.1e1 + 0.32163958997385070134e2 / t26;
    const double t30 = safe_math::log( t29 );
    const double t31 = t12 * t30;
    const double t32 = 0.3109070e-1 * t31;
    const double t33 = rho * rho;
    const double t35 = 0.1e1 / t7 / t33;
    const double t46 = safe_math::exp( 0.62181400000000000000e-1 * t31 * t43 );
    const double t47 = t46 - 0.1e1;
    const double t48 = 0.1e1 / t47;
    const double t50 = t42 * t43 * t48;
    const double t51 = sigma * sigma;
    const double t52 = t33 * t33;
    const double t54 = 0.1e1 / t21 / t52;
    const double t62 = sigma * t35 * t39 / 0.48e2 + t50 * t51 * t54 * t58 / 0.768e3;
    const double t63 = beta * t62;
    const double t67 = t64 * t48 * t62 + 0.1e1;
    const double t68 = 0.1e1 / t67;
    const double t69 = t43 * t68;
    const double t71 = t63 * t69 + 0.1e1;
    const double t72 = safe_math::log( t71 );
    const double t74 = gamma * t72 / 0.2e1;
    const double t76 = 0.1e1 / t7 / rho;
    const double t77 = t6 * t76;
    const double t79 = t4 * t77 * t30;
    const double t81 = t26 * t26;
    const double t82 = 0.1e1 / t81;
    const double t83 = t12 * t82;
    const double t85 = 0.1e1 / t13 * t1;
    const double t87 = t86 * t76;
    const double t90 = t4 * t77;
    const double t92 = safe_math::sqrt( t10 );
    const double t93 = t92 * t1;
    const double t101 = -0.11765750000000000000e1 * t85 * t87 - 0.51647500000000000000e0 * t90 - 0.21038750000000000000e0 * t93 * t87 - 0.10419500000000000000e0 * t20 * t5 / t21 / rho;
    const double t102 = 0.1e1 / t29;
    const double t103 = t101 * t102;
    const double t104 = t83 * t103;
    const double t106 = t33 * rho;
    const double t108 = 0.1e1 / t7 / t106;
    const double t112 = t47 * t47;
    const double t113 = 0.1e1 / t112;
    const double t114 = t43 * t113;
    const double t116 = t42 * t114 * t51;
    const double t117 = t54 * t1;
    const double t118 = t117 * t56;
    const double t127 = -0.10647528393333333333e-2 * t119 * t76 * t30 * t43 - 0.20000000000000000000e1 * t83 * t103 * t43;
    const double t129 = t6 * t127 * t46;
    const double t130 = t118 * t129;
    const double t133 = t52 * rho;
    const double t135 = 0.1e1 / t21 / t133;
    const double t140 = -0.7e1 / 0.144e3 * sigma * t108 * t39 - t116 * t130 / 0.768e3 - 0.7e1 / 0.1152e4 * t50 * t51 * t135 * t58;
    const double t141 = beta * t140;
    const double t143 = t67 * t67;
    const double t144 = 0.1e1 / t143;
    const double t145 = t43 * t144;
    const double t146 = t64 * t113;
    const double t147 = t62 * t127;
    const double t152 = t64 * t48 * t140 - t146 * t147 * t46;
    const double t153 = t145 * t152;
    const double t155 = t141 * t69 - t63 * t153;
    const double t157 = 0.1e1 / t71;
    const double t158 = gamma * t155 * t157;
    const double t162 = rho * gamma;
    const double t171 = t35 * t18 * t164 / 0.48e2 + t50 * sigma * t54 * t58 / 0.384e3;
    const double t172 = beta * t171;
    const double t175 = t174 * t62;
    const double t178 = t175 * t177;
    const double t179 = t144 * t48;
    const double t180 = t179 * t171;
    const double t182 = t172 * t69 - t178 * t180;


    eps = -t32 + t74;
    vrho = -t32 + t74 + rho * ( 0.53237641966666666666e-3 * t79 + 0.10000000000000000000e1 * t104 + t158 / 0.2e1 );
    vsigma = t162 * t182 * t157 / 0.2e1;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t19 = t1 * t1;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = t19 * t20;
    constexpr double t101 = 0.1e1 / t3;
    constexpr double t102 = t101 * t5;
    constexpr double t106 = BB * beta;
    constexpr double t107 = 0.1e1 / gamma;
    constexpr double t121 = t50 * t50;
    constexpr double t126 = 0.1e1 / t20;
    constexpr double t127 = t1 * t126;
    constexpr double t128 = t127 * t6;
    constexpr double t134 = beta * t107;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t11 = t4 * t6 / t8;
    const double t13 = 0.1e1 + 0.53425000000000000000e-1 * t11;
    const double t14 = safe_math::sqrt( t11 );
    const double t17 = pow_3_2( t11 );
    const double t22 = t8 * t8;
    const double t25 = t21 * t5 / t22;
    const double t27 = 0.37978500000000000000e1 * t14 + 0.89690000000000000000e0 * t11 + 0.20477500000000000000e0 * t17 + 0.12323500000000000000e0 * t25;
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
    const double t44 = safe_math::cbrt( t43 );
    const double t45 = t44 * t43;
    const double t46 = 0.1e1 - t42;
    const double t47 = safe_math::cbrt( t46 );
    const double t48 = t47 * t46;
    const double t49 = t45 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t54 = t49 * t53;
    const double t56 = 0.1e1 + 0.51370000000000000000e-1 * t11;
    const double t61 = 0.70594500000000000000e1 * t14 + 0.15494250000000000000e1 * t11 + 0.42077500000000000000e0 * t17 + 0.15629250000000000000e0 * t25;
    const double t64 = 0.1e1 + 0.32163958997385070134e2 / t61;
    const double t65 = safe_math::log( t64 );
    const double t69 = 0.1e1 + 0.27812500000000000000e-1 * t11;
    const double t74 = 0.51785000000000000000e1 * t14 + 0.90577500000000000000e0 * t11 + 0.11003250000000000000e0 * t17 + 0.12417750000000000000e0 * t25;
    const double t77 = 0.1e1 + 0.29608749977793437516e2 / t74;
    const double t78 = safe_math::log( t77 );
    const double t79 = t69 * t78;
    const double t81 = -0.3109070e-1 * t56 * t65 + t33 - 0.19751673498613801407e-1 * t79;
    const double t82 = t54 * t81;
    const double t83 = t40 * t82;
    const double t85 = 0.19751673498613801407e-1 * t54 * t79;
    const double t86 = t44 * t44;
    const double t87 = t47 * t47;
    const double t89 = t86 / 0.2e1 + t87 / 0.2e1;
    const double t90 = t89 * t89;
    const double t91 = t90 * t89;
    const double t92 = gamma * t91;
    const double t94 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t96 = 0.1e1 / t8 / t37;
    const double t97 = t94 * t96;
    const double t99 = 0.1e1 / t90;
    const double t103 = t99 * t19 * t102;
    const double t109 = ( -t33 + t83 + t85 ) * t107;
    const double t110 = 0.1e1 / t91;
    const double t112 = safe_math::exp( -t109 * t110 );
    const double t113 = t112 - 0.1e1;
    const double t114 = 0.1e1 / t113;
    const double t115 = t107 * t114;
    const double t116 = t94 * t94;
    const double t118 = t106 * t115 * t116;
    const double t120 = 0.1e1 / t22 / t38;
    const double t122 = t120 * t121;
    const double t123 = t90 * t90;
    const double t124 = 0.1e1 / t123;
    const double t125 = t122 * t124;
    const double t129 = t125 * t128;
    const double t132 = t97 * t50 * t103 / 0.96e2 + t118 * t129 / 0.3072e4;
    const double t133 = beta * t132;
    const double t137 = t114 * t132 * t134 + 0.1e1;
    const double t138 = 0.1e1 / t137;
    const double t139 = t107 * t138;
    const double t141 = t133 * t139 + 0.1e1;
    const double t142 = safe_math::log( t141 );
    const double t143 = t92 * t142;


    eps = -t33 + t83 + t85 + t143;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t19 = t1 * t1;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = t19 * t20;
    constexpr double t101 = 0.1e1 / t3;
    constexpr double t102 = t101 * t5;
    constexpr double t106 = BB * beta;
    constexpr double t107 = 0.1e1 / gamma;
    constexpr double t121 = t50 * t50;
    constexpr double t126 = 0.1e1 / t20;
    constexpr double t127 = t1 * t126;
    constexpr double t128 = t127 * t6;
    constexpr double t134 = beta * t107;
    constexpr double t155 = t3 * t6;
    constexpr double t259 = t19 * t101;
    constexpr double t264 = t106 * t107;
    constexpr double t272 = t126 * t6;
    constexpr double t381 = t259 * t5;
    constexpr double t391 = beta * beta;
    constexpr double t393 = gamma * gamma;
    constexpr double t394 = 0.1e1 / t393;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t11 = t4 * t6 / t8;
    const double t13 = 0.1e1 + 0.53425000000000000000e-1 * t11;
    const double t14 = safe_math::sqrt( t11 );
    const double t17 = pow_3_2( t11 );
    const double t22 = t8 * t8;
    const double t25 = t21 * t5 / t22;
    const double t27 = 0.37978500000000000000e1 * t14 + 0.89690000000000000000e0 * t11 + 0.20477500000000000000e0 * t17 + 0.12323500000000000000e0 * t25;
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
    const double t44 = safe_math::cbrt( t43 );
    const double t45 = t44 * t43;
    const double t46 = 0.1e1 - t42;
    const double t47 = safe_math::cbrt( t46 );
    const double t48 = t47 * t46;
    const double t49 = t45 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t54 = t49 * t53;
    const double t56 = 0.1e1 + 0.51370000000000000000e-1 * t11;
    const double t61 = 0.70594500000000000000e1 * t14 + 0.15494250000000000000e1 * t11 + 0.42077500000000000000e0 * t17 + 0.15629250000000000000e0 * t25;
    const double t64 = 0.1e1 + 0.32163958997385070134e2 / t61;
    const double t65 = safe_math::log( t64 );
    const double t69 = 0.1e1 + 0.27812500000000000000e-1 * t11;
    const double t74 = 0.51785000000000000000e1 * t14 + 0.90577500000000000000e0 * t11 + 0.11003250000000000000e0 * t17 + 0.12417750000000000000e0 * t25;
    const double t77 = 0.1e1 + 0.29608749977793437516e2 / t74;
    const double t78 = safe_math::log( t77 );
    const double t79 = t69 * t78;
    const double t81 = -0.3109070e-1 * t56 * t65 + t33 - 0.19751673498613801407e-1 * t79;
    const double t82 = t54 * t81;
    const double t83 = t40 * t82;
    const double t85 = 0.19751673498613801407e-1 * t54 * t79;
    const double t86 = t44 * t44;
    const double t87 = t47 * t47;
    const double t89 = t86 / 0.2e1 + t87 / 0.2e1;
    const double t90 = t89 * t89;
    const double t91 = t90 * t89;
    const double t92 = gamma * t91;
    const double t94 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t96 = 0.1e1 / t8 / t37;
    const double t97 = t94 * t96;
    const double t99 = 0.1e1 / t90;
    const double t103 = t99 * t19 * t102;
    const double t109 = ( -t33 + t83 + t85 ) * t107;
    const double t110 = 0.1e1 / t91;
    const double t112 = safe_math::exp( -t109 * t110 );
    const double t113 = t112 - 0.1e1;
    const double t114 = 0.1e1 / t113;
    const double t115 = t107 * t114;
    const double t116 = t94 * t94;
    const double t118 = t106 * t115 * t116;
    const double t120 = 0.1e1 / t22 / t38;
    const double t122 = t120 * t121;
    const double t123 = t90 * t90;
    const double t124 = 0.1e1 / t123;
    const double t125 = t122 * t124;
    const double t129 = t125 * t128;
    const double t132 = t97 * t50 * t103 / 0.96e2 + t118 * t129 / 0.3072e4;
    const double t133 = beta * t132;
    const double t137 = t114 * t132 * t134 + 0.1e1;
    const double t138 = 0.1e1 / t137;
    const double t139 = t107 * t138;
    const double t141 = t133 * t139 + 0.1e1;
    const double t142 = safe_math::log( t141 );
    const double t143 = t92 * t142;
    const double t145 = 0.1e1 / t8 / t7;
    const double t146 = t6 * t145;
    const double t148 = t4 * t146 * t31;
    const double t149 = 0.11073470983333333333e-2 * t148;
    const double t150 = t27 * t27;
    const double t151 = 0.1e1 / t150;
    const double t152 = t13 * t151;
    const double t154 = 0.1e1 / t14 * t1;
    const double t156 = t155 * t145;
    const double t157 = t154 * t156;
    const double t159 = t4 * t146;
    const double t161 = safe_math::sqrt( t11 );
    const double t162 = t161 * t1;
    const double t163 = t162 * t156;
    const double t168 = t21 * t5 / t22 / t7;
    const double t170 = -0.63297500000000000000e0 * t157 - 0.29896666666666666667e0 * t159 - 0.10238750000000000000e0 * t163 - 0.82156666666666666667e-1 * t168;
    const double t171 = 0.1e1 / t30;
    const double t172 = t170 * t171;
    const double t173 = t152 * t172;
    const double t174 = 0.10000000000000000000e1 * t173;
    const double t175 = t35 * t34;
    const double t176 = t175 * t39;
    const double t177 = t176 * t82;
    const double t178 = 0.4e1 * t177;
    const double t179 = t38 * t7;
    const double t180 = 0.1e1 / t179;
    const double t181 = t36 * t180;
    const double t182 = t181 * t82;
    const double t183 = 0.4e1 * t182;
    const double t184 = 0.1e1 / t37;
    const double t185 = t34 * t184;
    const double t186 = t41 - t185;
    const double t188 = -t186;
    const double t192 = ( 0.4e1 / 0.3e1 * t44 * t186 + 0.4e1 / 0.3e1 * t47 * t188 ) * t53;
    const double t193 = t192 * t81;
    const double t194 = t40 * t193;
    const double t198 = t61 * t61;
    const double t199 = 0.1e1 / t198;
    const double t200 = t56 * t199;
    const double t205 = -0.11765750000000000000e1 * t157 - 0.51647500000000000000e0 * t159 - 0.21038750000000000000e0 * t163 - 0.10419500000000000000e0 * t168;
    const double t206 = 0.1e1 / t64;
    const double t207 = t205 * t206;
    const double t213 = t74 * t74;
    const double t214 = 0.1e1 / t213;
    const double t215 = t69 * t214;
    const double t220 = -0.86308333333333333334e0 * t157 - 0.30192500000000000000e0 * t159 - 0.55016250000000000000e-1 * t163 - 0.82785000000000000000e-1 * t168;
    const double t221 = 0.1e1 / t77;
    const double t222 = t220 * t221;
    const double t225 = 0.53237641966666666666e-3 * t4 * t146 * t65 + 0.10000000000000000000e1 * t200 * t207 - t149 - t174 + 0.18311447306006545054e-3 * t4 * t146 * t78 + 0.58482236226346462070e0 * t215 * t222;
    const double t226 = t54 * t225;
    const double t227 = t40 * t226;
    const double t228 = t192 * t79;
    const double t229 = 0.19751673498613801407e-1 * t228;
    const double t230 = t54 * t1;
    const double t232 = t155 * t145 * t78;
    const double t233 = t230 * t232;
    const double t234 = 0.18311447306006545054e-3 * t233;
    const double t235 = t54 * t69;
    const double t237 = t214 * t220 * t221;
    const double t238 = t235 * t237;
    const double t239 = 0.58482236226346462070e0 * t238;
    const double t240 = gamma * t90;
    const double t241 = 0.1e1 / t44;
    const double t243 = 0.1e1 / t47;
    const double t246 = t241 * t186 / 0.3e1 + t243 * t188 / 0.3e1;
    const double t248 = t240 * t142 * t246;
    const double t249 = 0.3e1 * t248;
    const double t250 = t37 * t7;
    const double t252 = 0.1e1 / t8 / t250;
    const double t253 = t94 * t252;
    const double t256 = 0.7e1 / 0.288e3 * t253 * t50 * t103;
    const double t257 = t50 * t110;
    const double t258 = t97 * t257;
    const double t260 = t5 * t246;
    const double t261 = t259 * t260;
    const double t265 = t113 * t113;
    const double t266 = 0.1e1 / t265;
    const double t267 = t266 * t116;
    const double t269 = t264 * t267 * t120;
    const double t270 = t121 * t124;
    const double t271 = t270 * t1;
    const double t274 = ( t149 + t174 + t178 - t183 + t194 + t227 + t229 - t234 - t239 ) * t107;
    const double t276 = t124 * t246;
    const double t279 = 0.3e1 * t109 * t276 - t110 * t274;
    const double t280 = t279 * t112;
    const double t281 = t272 * t280;
    const double t282 = t271 * t281;
    const double t286 = 0.1e1 / t22 / t179;
    const double t287 = t286 * t121;
    const double t288 = t287 * t124;
    const double t289 = t288 * t128;
    const double t291 = 0.7e1 / 0.4608e4 * t118 * t289;
    const double t292 = t114 * t116;
    const double t294 = t264 * t292 * t120;
    const double t296 = 0.1e1 / t123 / t89;
    const double t297 = t121 * t296;
    const double t298 = t297 * t1;
    const double t300 = t298 * t272 * t246;
    const double t303 = -t256 - t258 * t261 / 0.48e2 - t269 * t282 / 0.3072e4 - t291 - t294 * t300 / 0.768e3;
    const double t304 = beta * t303;
    const double t306 = t137 * t137;
    const double t307 = 0.1e1 / t306;
    const double t308 = t107 * t307;
    const double t309 = t134 * t266;
    const double t310 = t132 * t279;
    const double t315 = -t112 * t309 * t310 + t114 * t134 * t303;
    const double t316 = t308 * t315;
    const double t318 = -t133 * t316 + t139 * t304;
    const double t319 = 0.1e1 / t141;
    const double t320 = t318 * t319;
    const double t321 = t92 * t320;
    const double t322 = t149 + t174 + t178 - t183 + t194 + t227 + t229 - t234 - t239 + t249 + t321;
    const double t324 = -t41 - t185;
    const double t326 = -t324;
    const double t330 = ( 0.4e1 / 0.3e1 * t44 * t324 + 0.4e1 / 0.3e1 * t47 * t326 ) * t53;
    const double t331 = t330 * t81;
    const double t332 = t40 * t331;
    const double t333 = t330 * t79;
    const double t334 = 0.19751673498613801407e-1 * t333;
    const double t338 = t241 * t324 / 0.3e1 + t243 * t326 / 0.3e1;
    const double t339 = t142 * t338;
    const double t340 = t240 * t339;
    const double t341 = 0.3e1 * t340;
    const double t342 = t5 * t338;
    const double t343 = t259 * t342;
    const double t347 = ( t149 + t174 - t178 - t183 + t332 + t227 + t334 - t234 - t239 ) * t107;
    const double t349 = t124 * t338;
    const double t352 = 0.3e1 * t109 * t349 - t110 * t347;
    const double t353 = t352 * t112;
    const double t354 = t272 * t353;
    const double t355 = t271 * t354;
    const double t359 = t298 * t272 * t338;
    const double t362 = -t256 - t258 * t343 / 0.48e2 - t269 * t355 / 0.3072e4 - t291 - t294 * t359 / 0.768e3;
    const double t363 = beta * t362;
    const double t365 = t132 * t352;
    const double t370 = -t112 * t309 * t365 + t114 * t134 * t362;
    const double t371 = t308 * t370;
    const double t373 = -t133 * t371 + t139 * t363;
    const double t374 = t373 * t319;
    const double t375 = t92 * t374;
    const double t376 = t149 + t174 - t178 - t183 + t332 + t227 + t334 - t234 - t239 + t341 + t375;
    const double t378 = t7 * gamma;
    const double t379 = t96 * t50;
    const double t382 = t379 * t99 * t381;
    const double t385 = t106 * t115 * t94;
    const double t386 = t385 * t129;
    const double t388 = t382 / 0.96e2 + t386 / 0.1536e4;
    const double t389 = beta * t388;
    const double t392 = t391 * t132;
    const double t395 = t392 * t394;
    const double t396 = t307 * t114;
    const double t397 = t396 * t388;
    const double t399 = t139 * t389 - t395 * t397;
    const double t404 = t382 / 0.48e2 + t386 / 0.768e3;
    const double t405 = beta * t404;
    const double t407 = t396 * t404;
    const double t409 = t139 * t405 - t395 * t407;
    const double t410 = t91 * t409;


    eps = -t33 + t83 + t85 + t143;
    vrho_a = t322 * t7 + t143 - t33 + t83 + t85;
    vrho_b = t376 * t7 + t143 - t33 + t83 + t85;
    vsigma_aa = t378 * t91 * t399 * t319;
    vsigma_ab = t378 * t410 * t319;
    vsigma_bb = vsigma_aa;

  }


};

struct BuiltinPBE_C : detail::BuiltinKernelImpl< BuiltinPBE_C > {

  BuiltinPBE_C( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_C >(p) { }
  
  virtual ~BuiltinPBE_C() = default;

};



} // namespace ExchCXX