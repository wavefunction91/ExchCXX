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

  static constexpr double dens_tol  = 1e-25;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double beta = 0.0042;
  static constexpr double gamma = 6.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t1 * t3 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t12 = t1 * t1;
    constexpr double t13 = beta * t12;
    constexpr double t14 = 0.1e1 / t3;
    constexpr double t15 = t14 * t5;
    constexpr double t16 = t13 * t15;
    constexpr double t22 = gamma * beta;


    const double t10 = cbrt( rho );
    const double t11 = t9 * t10;
    const double t17 = sigma * t9;
    const double t18 = rho * rho;
    const double t19 = t10 * t10;
    const double t21 = 0.1e1 / t19 / t18;
    const double t23 = sqrt( sigma );
    const double t24 = t22 * t23;
    const double t25 = t10 * rho;
    const double t26 = 0.1e1 / t25;
    const double t30 = log( t23 * t8 * t26 + sqrt( square( t23 * t8 * t26 ) + 0.1e1 ) );
    const double t31 = t8 * t26 * t30;
    const double t33 = 0.10e1 + t24 * t31;
    const double t34 = 0.1e1 / t33;
    const double t35 = t21 * t34;
    const double t39 = 0.10e1 + 0.2e1 / 0.9e1 * t16 * t17 * t35;
    const double t41 = t7 * t11 * t39;


    eps = -0.3e1 / 0.16e2 * t41;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t1 * t3 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t12 = t1 * t1;
    constexpr double t13 = beta * t12;
    constexpr double t14 = 0.1e1 / t3;
    constexpr double t15 = t14 * t5;
    constexpr double t16 = t13 * t15;
    constexpr double t22 = gamma * beta;
    constexpr double t46 = t6 * t9;
    constexpr double t80 = t13 * t14;
    constexpr double t81 = t5 * t9;


    const double t10 = cbrt( rho );
    const double t11 = t9 * t10;
    const double t17 = sigma * t9;
    const double t18 = rho * rho;
    const double t19 = t10 * t10;
    const double t21 = 0.1e1 / t19 / t18;
    const double t23 = sqrt( sigma );
    const double t24 = t22 * t23;
    const double t25 = t10 * rho;
    const double t26 = 0.1e1 / t25;
    const double t30 = log( t23 * t8 * t26 + sqrt( square( t23 * t8 * t26 ) + 0.1e1 ) );
    const double t31 = t8 * t26 * t30;
    const double t33 = 0.10e1 + t24 * t31;
    const double t34 = 0.1e1 / t33;
    const double t35 = t21 * t34;
    const double t39 = 0.10e1 + 0.2e1 / 0.9e1 * t16 * t17 * t35;
    const double t41 = t7 * t11 * t39;
    const double t45 = t25 * t1 * t3;
    const double t47 = t18 * rho;
    const double t49 = 0.1e1 / t19 / t47;
    const double t50 = t49 * t34;
    const double t54 = t33 * t33;
    const double t55 = 0.1e1 / t54;
    const double t56 = t21 * t55;
    const double t60 = t8 / t10 / t18 * t30;
    const double t62 = t22 * sigma;
    const double t63 = t9 * t49;
    const double t65 = t17 * t21 + 0.1e1;
    const double t66 = sqrt( t65 );
    const double t67 = 0.1e1 / t66;
    const double t68 = t63 * t67;
    const double t71 = -0.4e1 / 0.3e1 * t24 * t60 - 0.4e1 / 0.3e1 * t62 * t68;
    const double t76 = -0.16e2 / 0.27e2 * t16 * t17 * t50 - 0.2e1 / 0.9e1 * t16 * t17 * t56 * t71;
    const double t85 = t22 / t23;
    const double t87 = t9 * t21;
    const double t88 = t87 * t67;
    const double t91 = t22 * t88 / 0.2e1 + t85 * t31 / 0.2e1;
    const double t97 = t46 * ( -0.2e1 / 0.9e1 * t16 * t17 * t56 * t91 + 0.2e1 / 0.9e1 * t80 * t81 * t35 );


    eps = -0.3e1 / 0.16e2 * t41;
    vrho = -t41 / 0.4e1 - 0.3e1 / 0.16e2 * t45 * t46 * t76;
    vsigma = -0.3e1 / 0.16e2 * t45 * t97;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double sigma, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t9 = t1 * t1;
    constexpr double t10 = beta * t9;
    constexpr double t11 = 0.1e1 / t3;
    constexpr double t12 = t10 * t11;
    constexpr double t18 = gamma * beta;


    const double t7 = cbrt( rho );
    const double t8 = t6 * t7;
    const double t13 = t5 * sigma;
    const double t14 = rho * rho;
    const double t15 = t7 * t7;
    const double t17 = 0.1e1 / t15 / t14;
    const double t19 = sqrt( sigma );
    const double t20 = t7 * rho;
    const double t21 = 0.1e1 / t20;
    const double t22 = t19 * t21;
    const double t23 = log( t22 + sqrt( t22 * t22 + 0.1e1 ) );
    const double t26 = 0.10e1 + t18 * t22 * t23;
    const double t27 = 0.1e1 / t26;
    const double t32 = 0.10e1 + 0.2e1 / 0.9e1 * t12 * t13 * t17 * t27;
    const double t34 = t4 * t8 * t32;


    eps = -0.3e1 / 0.8e1 * t34;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t9 = t1 * t1;
    constexpr double t10 = beta * t9;
    constexpr double t11 = 0.1e1 / t3;
    constexpr double t12 = t10 * t11;
    constexpr double t18 = gamma * beta;
    constexpr double t38 = t3 * t6;
    constexpr double t46 = t11 * t5;
    constexpr double t47 = t10 * t46;


    const double t7 = cbrt( rho );
    const double t8 = t6 * t7;
    const double t13 = t5 * sigma;
    const double t14 = rho * rho;
    const double t15 = t7 * t7;
    const double t17 = 0.1e1 / t15 / t14;
    const double t19 = sqrt( sigma );
    const double t20 = t7 * rho;
    const double t21 = 0.1e1 / t20;
    const double t22 = t19 * t21;
    const double t23 = log( t22 + sqrt( t22 * t22 + 0.1e1 ) );
    const double t26 = 0.10e1 + t18 * t22 * t23;
    const double t27 = 0.1e1 / t26;
    const double t32 = 0.10e1 + 0.2e1 / 0.9e1 * t12 * t13 * t17 * t27;
    const double t34 = t4 * t8 * t32;
    const double t37 = t20 * t1;
    const double t39 = t14 * rho;
    const double t41 = 0.1e1 / t15 / t39;
    const double t48 = sigma * t17;
    const double t49 = t26 * t26;
    const double t50 = 0.1e1 / t49;
    const double t52 = 0.1e1 / t7 / t14;
    const double t56 = sigma * t41;
    const double t57 = t48 + 0.1e1;
    const double t58 = sqrt( t57 );
    const double t59 = 0.1e1 / t58;
    const double t63 = -0.4e1 / 0.3e1 * t18 * t19 * t52 * t23 - 0.4e1 / 0.3e1 * t18 * t56 * t59;
    const double t64 = t50 * t63;
    const double t68 = -0.16e2 / 0.27e2 * t12 * t13 * t41 * t27 - 0.2e1 / 0.9e1 * t47 * t48 * t64;
    const double t72 = t5 * t17;
    const double t75 = 0.1e1 / t19;
    const double t82 = t18 * t75 * t21 * t23 / 0.2e1 + t18 * t17 * t59 / 0.2e1;
    const double t83 = t50 * t82;
    const double t88 = t38 * ( 0.2e1 / 0.9e1 * t12 * t72 * t27 - 0.2e1 / 0.9e1 * t47 * t48 * t83 );


    eps = -0.3e1 / 0.8e1 * t34;
    vrho = -t34 / 0.2e1 - 0.3e1 / 0.8e1 * t37 * t38 * t68;
    vsigma = -0.3e1 / 0.8e1 * t37 * t88;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = beta * t18;
    constexpr double t20 = 0.1e1 / t3;
    constexpr double t21 = t19 * t20;
    constexpr double t28 = gamma * beta;


    const double t8 = rho_a - rho_b;
    const double t9 = rho_a + rho_b;
    const double t10 = 0.1e1 / t9;
    const double t11 = t8 * t10;
    const double t13 = 0.1e1 / 0.2e1 + t11 / 0.2e1;
    const double t14 = cbrt( t13 );
    const double t15 = t14 * t13;
    const double t16 = cbrt( t9 );
    const double t17 = t15 * t16;
    const double t22 = t5 * sigma_aa;
    const double t23 = rho_a * rho_a;
    const double t24 = cbrt( rho_a );
    const double t25 = t24 * t24;
    const double t27 = 0.1e1 / t25 / t23;
    const double t29 = sqrt( sigma_aa );
    const double t31 = 0.1e1 / t24 / rho_a;
    const double t32 = t29 * t31;
    const double t33 = log( t32 + sqrt( t32 * t32 + 0.1e1 ) );
    const double t36 = 0.10e1 + t28 * t32 * t33;
    const double t37 = 0.1e1 / t36;
    const double t42 = 0.10e1 + 0.2e1 / 0.9e1 * t21 * t22 * t27 * t37;
    const double t44 = t7 * t17 * t42;
    const double t46 = 0.1e1 / 0.2e1 - t11 / 0.2e1;
    const double t47 = cbrt( t46 );
    const double t48 = t47 * t46;
    const double t49 = t48 * t16;
    const double t50 = t5 * sigma_bb;
    const double t51 = rho_b * rho_b;
    const double t52 = cbrt( rho_b );
    const double t53 = t52 * t52;
    const double t55 = 0.1e1 / t53 / t51;
    const double t56 = sqrt( sigma_bb );
    const double t58 = 0.1e1 / t52 / rho_b;
    const double t59 = t56 * t58;
    const double t60 = log( t59 + sqrt( t59 * t59 + 0.1e1 ) );
    const double t63 = 0.10e1 + t28 * t59 * t60;
    const double t64 = 0.1e1 / t63;
    const double t69 = 0.10e1 + 0.2e1 / 0.9e1 * t21 * t50 * t55 * t64;
    const double t71 = t7 * t49 * t69;


    eps = -0.3e1 / 0.8e1 * t44 - 0.3e1 / 0.8e1 * t71;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t18 = t1 * t1;
    constexpr double t19 = beta * t18;
    constexpr double t20 = 0.1e1 / t3;
    constexpr double t21 = t19 * t20;
    constexpr double t28 = gamma * beta;
    constexpr double t98 = t20 * t5;
    constexpr double t99 = t19 * t98;


    const double t8 = rho_a - rho_b;
    const double t9 = rho_a + rho_b;
    const double t10 = 0.1e1 / t9;
    const double t11 = t8 * t10;
    const double t13 = 0.1e1 / 0.2e1 + t11 / 0.2e1;
    const double t14 = cbrt( t13 );
    const double t15 = t14 * t13;
    const double t16 = cbrt( t9 );
    const double t17 = t15 * t16;
    const double t22 = t5 * sigma_aa;
    const double t23 = rho_a * rho_a;
    const double t24 = cbrt( rho_a );
    const double t25 = t24 * t24;
    const double t27 = 0.1e1 / t25 / t23;
    const double t29 = sqrt( sigma_aa );
    const double t31 = 0.1e1 / t24 / rho_a;
    const double t32 = t29 * t31;
    const double t33 = log( t32 + sqrt( t32 * t32 + 0.1e1 ) );
    const double t36 = 0.10e1 + t28 * t32 * t33;
    const double t37 = 0.1e1 / t36;
    const double t42 = 0.10e1 + 0.2e1 / 0.9e1 * t21 * t22 * t27 * t37;
    const double t44 = t7 * t17 * t42;
    const double t46 = 0.1e1 / 0.2e1 - t11 / 0.2e1;
    const double t47 = cbrt( t46 );
    const double t48 = t47 * t46;
    const double t49 = t48 * t16;
    const double t50 = t5 * sigma_bb;
    const double t51 = rho_b * rho_b;
    const double t52 = cbrt( rho_b );
    const double t53 = t52 * t52;
    const double t55 = 0.1e1 / t53 / t51;
    const double t56 = sqrt( sigma_bb );
    const double t58 = 0.1e1 / t52 / rho_b;
    const double t59 = t56 * t58;
    const double t60 = log( t59 + sqrt( t59 * t59 + 0.1e1 ) );
    const double t63 = 0.10e1 + t28 * t59 * t60;
    const double t64 = 0.1e1 / t63;
    const double t69 = 0.10e1 + 0.2e1 / 0.9e1 * t21 * t50 * t55 * t64;
    const double t71 = t7 * t49 * t69;
    const double t73 = 0.3e1 / 0.8e1 * t44;
    const double t74 = 0.3e1 / 0.8e1 * t71;
    const double t75 = t14 * t16;
    const double t76 = t9 * t9;
    const double t77 = 0.1e1 / t76;
    const double t78 = t8 * t77;
    const double t80 = t10 / 0.2e1 - t78 / 0.2e1;
    const double t81 = t42 * t80;
    const double t83 = t7 * t75 * t81;
    const double t84 = t83 / 0.2e1;
    const double t85 = t16 * t16;
    const double t86 = 0.1e1 / t85;
    const double t87 = t15 * t86;
    const double t89 = t7 * t87 * t42;
    const double t90 = t89 / 0.8e1;
    const double t91 = t23 * rho_a;
    const double t93 = 0.1e1 / t25 / t91;
    const double t100 = sigma_aa * t27;
    const double t101 = t36 * t36;
    const double t102 = 0.1e1 / t101;
    const double t104 = 0.1e1 / t24 / t23;
    const double t108 = sigma_aa * t93;
    const double t109 = t100 + 0.1e1;
    const double t110 = sqrt( t109 );
    const double t111 = 0.1e1 / t110;
    const double t115 = -0.4e1 / 0.3e1 * t28 * t29 * t104 * t33 - 0.4e1 / 0.3e1 * t28 * t108 * t111;
    const double t116 = t102 * t115;
    const double t120 = -0.16e2 / 0.27e2 * t21 * t22 * t93 * t37 - 0.2e1 / 0.9e1 * t99 * t100 * t116;
    const double t122 = t7 * t17 * t120;
    const double t123 = 0.3e1 / 0.8e1 * t122;
    const double t124 = t47 * t16;
    const double t125 = -t80;
    const double t126 = t69 * t125;
    const double t128 = t7 * t124 * t126;
    const double t129 = t128 / 0.2e1;
    const double t130 = t48 * t86;
    const double t132 = t7 * t130 * t69;
    const double t133 = t132 / 0.8e1;
    const double t137 = -t10 / 0.2e1 - t78 / 0.2e1;
    const double t138 = t42 * t137;
    const double t140 = t7 * t75 * t138;
    const double t141 = t140 / 0.2e1;
    const double t142 = -t137;
    const double t143 = t69 * t142;
    const double t145 = t7 * t124 * t143;
    const double t146 = t145 / 0.2e1;
    const double t147 = t51 * rho_b;
    const double t149 = 0.1e1 / t53 / t147;
    const double t154 = sigma_bb * t55;
    const double t155 = t63 * t63;
    const double t156 = 0.1e1 / t155;
    const double t158 = 0.1e1 / t52 / t51;
    const double t162 = sigma_bb * t149;
    const double t163 = t154 + 0.1e1;
    const double t164 = sqrt( t163 );
    const double t165 = 0.1e1 / t164;
    const double t169 = -0.4e1 / 0.3e1 * t28 * t56 * t158 * t60 - 0.4e1 / 0.3e1 * t28 * t162 * t165;
    const double t170 = t156 * t169;
    const double t174 = -0.16e2 / 0.27e2 * t21 * t50 * t149 * t64 - 0.2e1 / 0.9e1 * t99 * t154 * t170;
    const double t176 = t7 * t49 * t174;
    const double t177 = 0.3e1 / 0.8e1 * t176;
    const double t181 = t16 * t9 * t1;
    const double t182 = t181 * t3;
    const double t183 = t6 * t15;
    const double t184 = t5 * t27;
    const double t187 = 0.1e1 / t29;
    const double t194 = t28 * t187 * t31 * t33 / 0.2e1 + t28 * t27 * t111 / 0.2e1;
    const double t195 = t102 * t194;
    const double t199 = -0.2e1 / 0.9e1 * t99 * t100 * t195 + 0.2e1 / 0.9e1 * t21 * t184 * t37;
    const double t203 = t6 * t48;
    const double t204 = t5 * t55;
    const double t207 = 0.1e1 / t56;
    const double t214 = t28 * t207 * t58 * t60 / 0.2e1 + t28 * t55 * t165 / 0.2e1;
    const double t215 = t156 * t214;
    const double t219 = -0.2e1 / 0.9e1 * t99 * t154 * t215 + 0.2e1 / 0.9e1 * t21 * t204 * t64;


    eps = -0.3e1 / 0.8e1 * t44 - 0.3e1 / 0.8e1 * t71;
    vrho_a = -t73 - t74 + t9 * ( -t84 - t90 - t123 - t129 - t133 );
    vrho_b = -t73 - t74 + t9 * ( -t141 - t90 - t146 - t133 - t177 );
    vsigma_aa = -0.3e1 / 0.8e1 * t182 * t183 * t199;
    vsigma_ab = 0.0e0;
    vsigma_bb = -0.3e1 / 0.8e1 * t182 * t203 * t219;

  }


};

struct BuiltinB88 : detail::BuiltinKernelImpl< BuiltinB88 > {

  BuiltinB88( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinB88 >(p) { }
  
  virtual ~BuiltinB88() = default;

};



} // namespace ExchCXX