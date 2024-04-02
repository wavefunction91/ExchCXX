#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinLYP > :
  public gga_screening_interface< BuiltinLYP > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;

  static constexpr double dens_tol  = 1e-32;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 2.1544346900318956e-43;
  static constexpr double tau_tol = 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double a = 0.04918;
  static constexpr double b = 0.132;
  static constexpr double c = 0.2533;
  static constexpr double d = 0.349;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t20 = constants::m_cbrt_3;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t48 = constants::m_cbrt_2;
    constexpr double t21 = t20 * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t49 = t48 * t48;


    const double t1 = safe_math::cbrt( rho );
    const double t2 = 0.1e1 / t1;
    const double t4 = d * t2 + 0.1e1;
    const double t5 = 0.1e1 / t4;
    const double t7 = safe_math::exp( -c * t2 );
    const double t8 = b * t7;
    const double t9 = rho * rho;
    const double t10 = t1 * t1;
    const double t12 = 0.1e1 / t10 / t9;
    const double t13 = sigma * t12;
    const double t15 = d * t5 + c;
    const double t16 = t15 * t2;
    const double t18 = -0.1e1 / 0.72e2 - 0.7e1 / 0.72e2 * t16;
    const double t26 = 0.1e1 <= zeta_tol;
    const double t27 = zeta_tol * zeta_tol;
    const double t28 = safe_math::cbrt( zeta_tol );
    const double t29 = t28 * t28;
    const double t31 = piecewise_functor_3( t26, t29 * t27, 1.0 );
    const double t35 = 0.5e1 / 0.2e1 - t16 / 0.18e2;
    const double t36 = t35 * sigma;
    const double t37 = t12 * t31;
    const double t40 = t16 - 0.11e2;
    const double t41 = t40 * sigma;
    const double t44 = piecewise_functor_3( t26, t29 * t27 * zeta_tol, 1.0 );
    const double t45 = t12 * t44;
    const double t50 = sigma * t49;
    const double t53 = piecewise_functor_3( t26, t27, 1.0 );
    const double t54 = t53 * sigma;
    const double t56 = t49 * t12 * t31;
    const double t62 = -t13 * t18 - 0.3e1 / 0.1e2 * t21 * t24 * t31 + t36 * t37 / 0.8e1 + t41 * t45 / 0.144e3 - t48 * ( 0.4e1 / 0.3e1 * t50 * t37 - t54 * t56 / 0.2e1 ) / 0.8e1;


    eps = a * ( t8 * t5 * t62 - t5 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t20 = constants::m_cbrt_3;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t48 = constants::m_cbrt_2;
    constexpr double t21 = t20 * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t49 = t48 * t48;
    constexpr double t74 = b * c;
    constexpr double t91 = d * d;


    const double t1 = safe_math::cbrt( rho );
    const double t2 = 0.1e1 / t1;
    const double t4 = d * t2 + 0.1e1;
    const double t5 = 0.1e1 / t4;
    const double t7 = safe_math::exp( -c * t2 );
    const double t8 = b * t7;
    const double t9 = rho * rho;
    const double t10 = t1 * t1;
    const double t12 = 0.1e1 / t10 / t9;
    const double t13 = sigma * t12;
    const double t15 = d * t5 + c;
    const double t16 = t15 * t2;
    const double t18 = -0.1e1 / 0.72e2 - 0.7e1 / 0.72e2 * t16;
    const double t26 = 0.1e1 <= zeta_tol;
    const double t27 = zeta_tol * zeta_tol;
    const double t28 = safe_math::cbrt( zeta_tol );
    const double t29 = t28 * t28;
    const double t31 = piecewise_functor_3( t26, t29 * t27, 1.0 );
    const double t35 = 0.5e1 / 0.2e1 - t16 / 0.18e2;
    const double t36 = t35 * sigma;
    const double t37 = t12 * t31;
    const double t40 = t16 - 0.11e2;
    const double t41 = t40 * sigma;
    const double t44 = piecewise_functor_3( t26, t29 * t27 * zeta_tol, 1.0 );
    const double t45 = t12 * t44;
    const double t50 = sigma * t49;
    const double t53 = piecewise_functor_3( t26, t27, 1.0 );
    const double t54 = t53 * sigma;
    const double t56 = t49 * t12 * t31;
    const double t62 = -t13 * t18 - 0.3e1 / 0.1e2 * t21 * t24 * t31 + t36 * t37 / 0.8e1 + t41 * t45 / 0.144e3 - t48 * ( 0.4e1 / 0.3e1 * t50 * t37 - t54 * t56 / 0.2e1 ) / 0.8e1;
    const double t66 = rho * a;
    const double t67 = t4 * t4;
    const double t68 = 0.1e1 / t67;
    const double t69 = t68 * d;
    const double t71 = 0.1e1 / t1 / rho;
    const double t75 = t74 * t71;
    const double t76 = t7 * t5;
    const double t77 = t76 * t62;
    const double t80 = t8 * t68;
    const double t81 = t62 * d;
    const double t85 = t9 * rho;
    const double t87 = 0.1e1 / t10 / t85;
    const double t88 = sigma * t87;
    const double t92 = t91 * t68;
    const double t94 = 0.1e1 / t10 / rho;
    const double t97 = t15 * t71 - t92 * t94;
    const double t98 = 0.7e1 / 0.216e3 * t97;
    const double t100 = t97 / 0.54e2;
    const double t101 = t100 * sigma;
    const double t104 = t87 * t31;
    const double t108 = -t97 / 0.3e1;
    const double t109 = t108 * sigma;
    const double t112 = t87 * t44;
    const double t118 = t49 * t87 * t31;
    const double t124 = 0.8e1 / 0.3e1 * t88 * t18 - t13 * t98 + t101 * t37 / 0.8e1 - t36 * t104 / 0.3e1 + t109 * t45 / 0.144e3 - t41 * t112 / 0.54e2 - t48 * ( -0.32e2 / 0.9e1 * t50 * t104 + 0.4e1 / 0.3e1 * t54 * t118 ) / 0.8e1;
    const double t127 = -t69 * t71 / 0.3e1 + t75 * t77 / 0.3e1 + t80 * t81 * t71 / 0.3e1 + t8 * t5 * t124;
    const double t129 = t66 * b;
    const double t138 = t53 * t49;
    const double t144 = -t12 * t18 + t35 * t12 * t31 / 0.8e1 + t40 * t12 * t44 / 0.144e3 - t48 * ( 0.4e1 / 0.3e1 * t56 - t138 * t37 / 0.2e1 ) / 0.8e1;
    const double t145 = t76 * t144;


    eps = a * ( t8 * t5 * t62 - t5 );
    vrho = t66 * t127 + eps;
    vsigma = t129 * t145;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(eps);
    constexpr double t32 = constants::m_cbrt_3;
    constexpr double t35 = constants::m_cbrt_pi_sq;
    constexpr double t62 = constants::m_cbrt_2;
    constexpr double t33 = t32 * t32;
    constexpr double t36 = t35 * t35;
    constexpr double t37 = t33 * t36;


    const double t1 = rho_a - rho_b;
    const double t2 = t1 * t1;
    const double t3 = rho_a + rho_b;
    const double t4 = t3 * t3;
    const double t5 = 0.1e1 / t4;
    const double t7 = -t2 * t5 + 0.1e1;
    const double t8 = safe_math::cbrt( t3 );
    const double t9 = 0.1e1 / t8;
    const double t11 = d * t9 + 0.1e1;
    const double t12 = 0.1e1 / t11;
    const double t15 = safe_math::exp( -c * t9 );
    const double t16 = b * t15;
    const double t18 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t19 = t8 * t8;
    const double t21 = 0.1e1 / t19 / t4;
    const double t22 = t18 * t21;
    const double t24 = d * t12 + c;
    const double t25 = t24 * t9;
    const double t27 = 0.47e2 - 0.7e1 * t25;
    const double t30 = t7 * t27 / 0.72e2 - 0.2e1 / 0.3e1;
    const double t38 = 0.1e1 / t3;
    const double t39 = t1 * t38;
    const double t40 = 0.1e1 + t39;
    const double t41 = t40 <= zeta_tol;
    const double t42 = zeta_tol * zeta_tol;
    const double t43 = safe_math::cbrt( zeta_tol );
    const double t44 = t43 * t43;
    const double t45 = t44 * t42;
    const double t46 = t40 * t40;
    const double t47 = safe_math::cbrt( t40 );
    const double t48 = t47 * t47;
    const double t49 = t48 * t46;
    const double t50 = piecewise_functor_3( t41, t45, t49 );
    const double t51 = 0.1e1 - t39;
    const double t52 = t51 <= zeta_tol;
    const double t53 = t51 * t51;
    const double t54 = safe_math::cbrt( t51 );
    const double t55 = t54 * t54;
    const double t56 = t55 * t53;
    const double t57 = piecewise_functor_3( t52, t45, t56 );
    const double t58 = t50 + t57;
    const double t63 = t62 * t7;
    const double t65 = 0.5e1 / 0.2e1 - t25 / 0.18e2;
    const double t66 = rho_a * rho_a;
    const double t67 = safe_math::cbrt( rho_a );
    const double t68 = t67 * t67;
    const double t70 = 0.1e1 / t68 / t66;
    const double t71 = sigma_aa * t70;
    const double t72 = t71 * t50;
    const double t73 = rho_b * rho_b;
    const double t74 = safe_math::cbrt( rho_b );
    const double t75 = t74 * t74;
    const double t77 = 0.1e1 / t75 / t73;
    const double t78 = sigma_bb * t77;
    const double t79 = t78 * t57;
    const double t80 = t72 + t79;
    const double t81 = t65 * t80;
    const double t84 = t25 - 0.11e2;
    const double t86 = t44 * t42 * zeta_tol;
    const double t89 = piecewise_functor_3( t41, t86, t48 * t46 * t40 );
    const double t93 = piecewise_functor_3( t52, t86, t55 * t53 * t51 );
    const double t95 = t71 * t89 + t78 * t93;
    const double t96 = t84 * t95;
    const double t101 = piecewise_functor_3( t41, t42, t46 );
    const double t102 = t101 * sigma_bb;
    const double t103 = t77 * t57;
    const double t106 = piecewise_functor_3( t52, t42, t53 );
    const double t107 = t106 * sigma_aa;
    const double t108 = t70 * t50;
    const double t114 = -t22 * t30 - 0.3e1 / 0.2e2 * t37 * t7 * t58 + t63 * t81 / 0.32e2 + t63 * t96 / 0.576e3 - t62 * ( 0.2e1 / 0.3e1 * t72 + 0.2e1 / 0.3e1 * t79 - t102 * t103 / 0.4e1 - t107 * t108 / 0.4e1 ) / 0.8e1;


    eps = a * ( t16 * t12 * t114 - t7 * t12 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    constexpr double t32 = constants::m_cbrt_3;
    constexpr double t35 = constants::m_cbrt_pi_sq;
    constexpr double t62 = constants::m_cbrt_2;
    constexpr double t33 = t32 * t32;
    constexpr double t36 = t35 * t35;
    constexpr double t37 = t33 * t36;
    constexpr double t134 = b * c;
    constexpr double t151 = d * d;


    const double t1 = rho_a - rho_b;
    const double t2 = t1 * t1;
    const double t3 = rho_a + rho_b;
    const double t4 = t3 * t3;
    const double t5 = 0.1e1 / t4;
    const double t7 = -t2 * t5 + 0.1e1;
    const double t8 = safe_math::cbrt( t3 );
    const double t9 = 0.1e1 / t8;
    const double t11 = d * t9 + 0.1e1;
    const double t12 = 0.1e1 / t11;
    const double t15 = safe_math::exp( -c * t9 );
    const double t16 = b * t15;
    const double t18 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t19 = t8 * t8;
    const double t21 = 0.1e1 / t19 / t4;
    const double t22 = t18 * t21;
    const double t24 = d * t12 + c;
    const double t25 = t24 * t9;
    const double t27 = 0.47e2 - 0.7e1 * t25;
    const double t30 = t7 * t27 / 0.72e2 - 0.2e1 / 0.3e1;
    const double t38 = 0.1e1 / t3;
    const double t39 = t1 * t38;
    const double t40 = 0.1e1 + t39;
    const double t41 = t40 <= zeta_tol;
    const double t42 = zeta_tol * zeta_tol;
    const double t43 = safe_math::cbrt( zeta_tol );
    const double t44 = t43 * t43;
    const double t45 = t44 * t42;
    const double t46 = t40 * t40;
    const double t47 = safe_math::cbrt( t40 );
    const double t48 = t47 * t47;
    const double t49 = t48 * t46;
    const double t50 = piecewise_functor_3( t41, t45, t49 );
    const double t51 = 0.1e1 - t39;
    const double t52 = t51 <= zeta_tol;
    const double t53 = t51 * t51;
    const double t54 = safe_math::cbrt( t51 );
    const double t55 = t54 * t54;
    const double t56 = t55 * t53;
    const double t57 = piecewise_functor_3( t52, t45, t56 );
    const double t58 = t50 + t57;
    const double t63 = t62 * t7;
    const double t65 = 0.5e1 / 0.2e1 - t25 / 0.18e2;
    const double t66 = rho_a * rho_a;
    const double t67 = safe_math::cbrt( rho_a );
    const double t68 = t67 * t67;
    const double t70 = 0.1e1 / t68 / t66;
    const double t71 = sigma_aa * t70;
    const double t72 = t71 * t50;
    const double t73 = rho_b * rho_b;
    const double t74 = safe_math::cbrt( rho_b );
    const double t75 = t74 * t74;
    const double t77 = 0.1e1 / t75 / t73;
    const double t78 = sigma_bb * t77;
    const double t79 = t78 * t57;
    const double t80 = t72 + t79;
    const double t81 = t65 * t80;
    const double t84 = t25 - 0.11e2;
    const double t86 = t44 * t42 * zeta_tol;
    const double t89 = piecewise_functor_3( t41, t86, t48 * t46 * t40 );
    const double t93 = piecewise_functor_3( t52, t86, t55 * t53 * t51 );
    const double t95 = t71 * t89 + t78 * t93;
    const double t96 = t84 * t95;
    const double t101 = piecewise_functor_3( t41, t42, t46 );
    const double t102 = t101 * sigma_bb;
    const double t103 = t77 * t57;
    const double t106 = piecewise_functor_3( t52, t42, t53 );
    const double t107 = t106 * sigma_aa;
    const double t108 = t70 * t50;
    const double t114 = -t22 * t30 - 0.3e1 / 0.2e2 * t37 * t7 * t58 + t63 * t81 / 0.32e2 + t63 * t96 / 0.576e3 - t62 * ( 0.2e1 / 0.3e1 * t72 + 0.2e1 / 0.3e1 * t79 - t102 * t103 / 0.4e1 - t107 * t108 / 0.4e1 ) / 0.8e1;
    const double t118 = t3 * a;
    const double t119 = t1 * t5;
    const double t120 = t4 * t3;
    const double t121 = 0.1e1 / t120;
    const double t122 = t2 * t121;
    const double t124 = -0.2e1 * t119 + 0.2e1 * t122;
    const double t126 = t11 * t11;
    const double t127 = 0.1e1 / t126;
    const double t128 = t7 * t127;
    const double t130 = 0.1e1 / t8 / t3;
    const double t131 = d * t130;
    const double t133 = t128 * t131 / 0.3e1;
    const double t135 = t134 * t130;
    const double t136 = t15 * t12;
    const double t137 = t136 * t114;
    const double t139 = t135 * t137 / 0.3e1;
    const double t140 = t16 * t127;
    const double t141 = t114 * d;
    const double t144 = t140 * t141 * t130 / 0.3e1;
    const double t146 = 0.1e1 / t19 / t120;
    const double t147 = t18 * t146;
    const double t149 = 0.8e1 / 0.3e1 * t147 * t30;
    const double t152 = t151 * t127;
    const double t154 = 0.1e1 / t19 / t3;
    const double t157 = t24 * t130 - t152 * t154;
    const double t158 = 0.7e1 / 0.3e1 * t157;
    const double t159 = t7 * t158;
    const double t161 = t124 * t27 / 0.72e2 + t159 / 0.72e2;
    const double t166 = t48 * t40;
    const double t167 = t38 - t119;
    const double t168 = t166 * t167;
    const double t170 = piecewise_functor_3( t41, 0.0, 0.8e1 / 0.3e1 * t168 );
    const double t171 = t55 * t51;
    const double t172 = -t167;
    const double t173 = t171 * t172;
    const double t175 = piecewise_functor_3( t52, 0.0, 0.8e1 / 0.3e1 * t173 );
    const double t176 = t170 + t175;
    const double t180 = t62 * t124;
    const double t183 = t157 / 0.54e2;
    const double t184 = t183 * t80;
    const double t186 = t63 * t184 / 0.32e2;
    const double t189 = 0.1e1 / t68 / t66 / rho_a;
    const double t190 = sigma_aa * t189;
    const double t191 = t190 * t50;
    const double t193 = t71 * t170;
    const double t194 = t78 * t175;
    const double t195 = -0.8e1 / 0.3e1 * t191 + t193 + t194;
    const double t196 = t65 * t195;
    const double t202 = -t157 / 0.3e1;
    const double t203 = t202 * t95;
    const double t205 = t63 * t203 / 0.576e3;
    const double t210 = piecewise_functor_3( t41, 0.0, 0.11e2 / 0.3e1 * t49 * t167 );
    const double t214 = piecewise_functor_3( t52, 0.0, 0.11e2 / 0.3e1 * t56 * t172 );
    const double t216 = -0.8e1 / 0.3e1 * t190 * t89 + t71 * t210 + t78 * t214;
    const double t217 = t84 * t216;
    const double t225 = piecewise_functor_3( t41, 0.0, 0.2e1 * t40 * t167 );
    const double t226 = t225 * sigma_bb;
    const double t229 = t77 * t175;
    const double t234 = piecewise_functor_3( t52, 0.0, 0.2e1 * t51 * t172 );
    const double t235 = t234 * sigma_aa;
    const double t238 = t189 * t50;
    const double t241 = t70 * t170;
    const double t247 = t149 - t22 * t161 - 0.3e1 / 0.2e2 * t37 * t124 * t58 - 0.3e1 / 0.2e2 * t37 * t7 * t176 + t180 * t81 / 0.32e2 + t186 + t63 * t196 / 0.32e2 + t180 * t96 / 0.576e3 + t205 + t63 * t217 / 0.576e3 - t62 * ( -0.16e2 / 0.9e1 * t191 + 0.2e1 / 0.3e1 * t193 + 0.2e1 / 0.3e1 * t194 - t226 * t103 / 0.4e1 - t102 * t229 / 0.4e1 - t235 * t108 / 0.4e1 + 0.2e1 / 0.3e1 * t107 * t238 - t107 * t241 / 0.4e1 ) / 0.8e1;
    const double t250 = t16 * t12 * t247 - t124 * t12 - t133 + t139 + t144;
    const double t253 = 0.2e1 * t119 + 0.2e1 * t122;
    const double t257 = t253 * t27 / 0.72e2 + t159 / 0.72e2;
    const double t262 = -t38 - t119;
    const double t263 = t166 * t262;
    const double t265 = piecewise_functor_3( t41, 0.0, 0.8e1 / 0.3e1 * t263 );
    const double t266 = -t262;
    const double t267 = t171 * t266;
    const double t269 = piecewise_functor_3( t52, 0.0, 0.8e1 / 0.3e1 * t267 );
    const double t270 = t265 + t269;
    const double t274 = t62 * t253;
    const double t277 = t71 * t265;
    const double t280 = 0.1e1 / t75 / t73 / rho_b;
    const double t281 = sigma_bb * t280;
    const double t282 = t281 * t57;
    const double t284 = t78 * t269;
    const double t285 = t277 - 0.8e1 / 0.3e1 * t282 + t284;
    const double t286 = t65 * t285;
    const double t293 = piecewise_functor_3( t41, 0.0, 0.11e2 / 0.3e1 * t49 * t262 );
    const double t299 = piecewise_functor_3( t52, 0.0, 0.11e2 / 0.3e1 * t56 * t266 );
    const double t301 = t71 * t293 - 0.8e1 / 0.3e1 * t281 * t93 + t78 * t299;
    const double t302 = t84 * t301;
    const double t310 = piecewise_functor_3( t41, 0.0, 0.2e1 * t40 * t262 );
    const double t311 = t310 * sigma_bb;
    const double t314 = t280 * t57;
    const double t317 = t77 * t269;
    const double t322 = piecewise_functor_3( t52, 0.0, 0.2e1 * t51 * t266 );
    const double t323 = t322 * sigma_aa;
    const double t326 = t70 * t265;
    const double t332 = t149 - t22 * t257 - 0.3e1 / 0.2e2 * t37 * t253 * t58 - 0.3e1 / 0.2e2 * t37 * t7 * t270 + t274 * t81 / 0.32e2 + t186 + t63 * t286 / 0.32e2 + t274 * t96 / 0.576e3 + t205 + t63 * t302 / 0.576e3 - t62 * ( 0.2e1 / 0.3e1 * t277 - 0.16e2 / 0.9e1 * t282 + 0.2e1 / 0.3e1 * t284 - t311 * t103 / 0.4e1 + 0.2e1 / 0.3e1 * t102 * t314 - t102 * t317 / 0.4e1 - t323 * t108 / 0.4e1 - t107 * t326 / 0.4e1 ) / 0.8e1;
    const double t335 = t16 * t12 * t332 - t253 * t12 - t133 + t139 + t144;
    const double t337 = t118 * b;
    const double t338 = t21 * t30;
    const double t339 = t65 * t70;
    const double t340 = t339 * t50;
    const double t343 = t84 * t70;
    const double t344 = t343 * t89;
    const double t348 = t106 * t70;
    const double t354 = -t338 + t63 * t340 / 0.32e2 + t63 * t344 / 0.576e3 - t62 * ( 0.2e1 / 0.3e1 * t108 - t348 * t50 / 0.4e1 ) / 0.8e1;
    const double t355 = t136 * t354;
    const double t356 = t154 * a;
    const double t357 = t356 * b;
    const double t358 = t136 * t30;
    const double t361 = t65 * t77;
    const double t362 = t361 * t57;
    const double t365 = t84 * t77;
    const double t366 = t365 * t93;
    const double t370 = t101 * t77;
    const double t376 = -t338 + t63 * t362 / 0.32e2 + t63 * t366 / 0.576e3 - t62 * ( 0.2e1 / 0.3e1 * t103 - t370 * t57 / 0.4e1 ) / 0.8e1;
    const double t377 = t136 * t376;


    eps = a * ( t16 * t12 * t114 - t7 * t12 );
    vrho_a = t118 * t250 + eps;
    vrho_b = t118 * t335 + eps;
    vsigma_aa = t337 * t355;
    vsigma_ab = -0.2e1 * t357 * t358;
    vsigma_bb = t337 * t377;

  }


};

struct BuiltinLYP : detail::BuiltinKernelImpl< BuiltinLYP > {

  BuiltinLYP( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinLYP >(p) { }
  
  virtual ~BuiltinLYP() = default;

};



} // namespace ExchCXX