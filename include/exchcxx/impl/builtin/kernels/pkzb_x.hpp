#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPKZB_X > :
  public mgga_screening_interface< BuiltinPKZB_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = true;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps ) {

    (void)(lapl);
    (void)(eps);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t71 = piecewise_functor_3( t3, 0.0, -0.3e1 / 0.8e1 * t7 * t18 * t19 * t67 );


    eps = 0.2e1 * t71;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau ) {

    (void)(lapl);
    (void)(eps);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;
    constexpr double t117 = t25 * t28;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t71 = piecewise_functor_3( t3, 0.0, -0.3e1 / 0.8e1 * t7 * t18 * t19 * t67 );
    const double t72 = 0.1e1 / t31;
    const double t77 = t4 * t18;
    const double t78 = t64 * t64;
    const double t79 = 0.1e1 / t78;
    const double t80 = t19 * t79;
    const double t81 = t30 * rho;
    const double t83 = 0.1e1 / t31 / t81;
    const double t84 = t29 * t83;
    const double t85 = t26 * t84;
    const double t91 = -0.5e1 / 0.12e2 * t26 * t37 * t33 + t85 / 0.108e3;
    const double t94 = t91 * t21;
    const double t95 = t94 * t25;
    const double t100 = t57 * t30;
    const double t102 = 0.1e1 / t19 / t100;
    const double t106 = -0.1e2 / 0.729e3 * t85 + 0.292e3 / 0.2025e4 * t44 * t91 - 0.73e2 / 0.972e4 * t95 * t34 + 0.73e2 / 0.3645e4 * t48 * t84 - 0.24436516267640330302e-2 * t54 * t56 * t102;
    const double t111 = piecewise_functor_3( t3, 0.0, -t7 * t18 * t72 * t67 / 0.8e1 - 0.16551095363746320496e0 * t77 * t80 * t106 );
    const double t118 = t117 * t33;
    const double t119 = t47 * t118;
    const double t123 = t54 * t27 * t60 * sigma;
    const double t125 = 0.5e1 / 0.972e3 * t26 * t28 * t33 - 0.146e3 / 0.18225e5 * t119 + 0.96852413827153753492e-3 * t123;
    const double t129 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t125 );
    const double t131 = t117 * t39;
    const double t140 = 0.73e2 / 0.2025e4 * t47 * t131 - 0.73e2 / 0.1944e5 * t54 * t27 / t19 / t57 * sigma;
    const double t144 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t140 );


    eps = 0.2e1 * t71;
    vrho = 0.2e1 * rho * t111 + 0.2e1 * t71;
    vsigma = 0.2e1 * rho * t129;
    vlapl = 0.e0;
    vtau = 0.2e1 * rho * t144;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t76 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t26 * t27 * t72 );
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t124 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t85 * t27 * t120 );


    eps = t76 + t124;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t76 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t26 * t27 * t72 );
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t124 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t85 * t27 * t120 );
    const double t125 = t7 * t7;
    const double t126 = 0.1e1 / t125;
    const double t127 = t17 * t126;
    const double t129 = piecewise_functor_5( t11, 0.0, t15, 0.0, t8 - t127 );
    const double t132 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t129 );
    const double t137 = t27 * t27;
    const double t138 = 0.1e1 / t137;
    const double t142 = t6 * t26 * t138 * t72 / 0.8e1;
    const double t143 = t3 * t26;
    const double t144 = t69 * t69;
    const double t145 = 0.1e1 / t144;
    const double t146 = t27 * t145;
    const double t147 = t35 * rho_a;
    const double t149 = 0.1e1 / t37 / t147;
    const double t151 = t34 * sigma_aa * t149;
    const double t157 = -0.5e1 / 0.12e2 * t34 * tau_a * t39 + t151 / 0.108e3;
    const double t160 = t157 * t29;
    const double t163 = t53 * t149;
    const double t166 = t62 * t35;
    const double t168 = 0.1e1 / t36 / t166;
    const double t172 = -0.1e2 / 0.729e3 * t151 + 0.292e3 / 0.2025e4 * t49 * t157 - 0.73e2 / 0.972e4 * t160 * t54 + 0.73e2 / 0.3645e4 * t52 * t163 - 0.12218258133820165151e-2 * t60 * t61 * t168;
    const double t173 = t146 * t172;
    const double t177 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t132 * t27 * t72 - t142 - 0.16551095363746320496e0 * t143 * t173 );
    const double t178 = t78 * t126;
    const double t180 = piecewise_functor_5( t15, 0.0, t11, 0.0, -t8 - t178 );
    const double t183 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t180 );
    const double t191 = t6 * t85 * t138 * t120 / 0.8e1;
    const double t193 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t183 * t27 * t120 - t191 );
    const double t197 = piecewise_functor_5( t11, 0.0, t15, 0.0, -t8 - t127 );
    const double t200 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t197 );
    const double t206 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t200 * t27 * t72 - t142 );
    const double t208 = piecewise_functor_5( t15, 0.0, t11, 0.0, t8 - t178 );
    const double t211 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t208 );
    const double t216 = t3 * t85;
    const double t217 = t117 * t117;
    const double t218 = 0.1e1 / t217;
    const double t219 = t27 * t218;
    const double t220 = t87 * rho_b;
    const double t222 = 0.1e1 / t89 / t220;
    const double t224 = t34 * sigma_bb * t222;
    const double t230 = -0.5e1 / 0.12e2 * t34 * tau_b * t91 + t224 / 0.108e3;
    const double t233 = t230 * t29;
    const double t236 = t105 * t222;
    const double t239 = t110 * t87;
    const double t241 = 0.1e1 / t88 / t239;
    const double t245 = -0.1e2 / 0.729e3 * t224 + 0.292e3 / 0.2025e4 * t101 * t230 - 0.73e2 / 0.972e4 * t233 * t106 + 0.73e2 / 0.3645e4 * t104 * t236 - 0.12218258133820165151e-2 * t60 * t109 * t241;
    const double t246 = t219 * t245;
    const double t250 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t211 * t27 * t120 - t191 - 0.16551095363746320496e0 * t216 * t246 );
    const double t255 = t33 * t39;
    const double t256 = t52 * t255;
    const double t259 = t60 * t65 * sigma_aa;
    const double t261 = 0.5e1 / 0.972e3 * t34 * t39 - 0.146e3 / 0.18225e5 * t256 + 0.48426206913576876746e-3 * t259;
    const double t262 = t146 * t261;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t262 );
    const double t268 = t33 * t91;
    const double t269 = t104 * t268;
    const double t272 = t60 * t113 * sigma_bb;
    const double t274 = 0.5e1 / 0.972e3 * t34 * t91 - 0.146e3 / 0.18225e5 * t269 + 0.48426206913576876746e-3 * t272;
    const double t275 = t219 * t274;
    const double t278 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t275 );
    const double t279 = t33 * t44;
    const double t283 = 0.1e1 / t36 / t62;
    const double t287 = 0.73e2 / 0.2025e4 * t52 * t279 - 0.73e2 / 0.3888e5 * t60 * t283 * sigma_aa;
    const double t288 = t146 * t287;
    const double t291 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t288 );
    const double t292 = t33 * t96;
    const double t296 = 0.1e1 / t88 / t110;
    const double t300 = 0.73e2 / 0.2025e4 * t104 * t292 - 0.73e2 / 0.3888e5 * t60 * t296 * sigma_bb;
    const double t301 = t219 * t300;
    const double t304 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t301 );


    eps = t76 + t124;
    vrho_a = t76 + t124 + t7 * ( t177 + t193 );
    vrho_b = t76 + t124 + t7 * ( t206 + t250 );
    vsigma_aa = t7 * t265;
    vsigma_ab = 0.e0;
    vsigma_bb = t7 * t278;
    vlapl_a = 0.e0;
    vlapl_b = 0.e0;
    vtau_a = t7 * t291;
    vtau_b = t7 * t304;

  }


};

struct BuiltinPKZB_X : detail::BuiltinKernelImpl< BuiltinPKZB_X > {

  BuiltinPKZB_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPKZB_X >(p) { }
  
  virtual ~BuiltinPKZB_X() = default;

};



} // namespace ExchCXX