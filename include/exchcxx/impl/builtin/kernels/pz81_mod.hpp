#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPZ81_MOD > :
  public lda_screening_interface< BuiltinPZ81_MOD > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double gamma_0 = -0.1423;
  static constexpr double gamma_1 = -0.0843;
  static constexpr double beta1_0 = 1.0529;
  static constexpr double beta1_1 = 1.3981;
  static constexpr double beta2_0 = 0.3334;
  static constexpr double beta2_1 = 0.2611;
  static constexpr double a_0 = 0.0311;
  static constexpr double a_1 = 0.01555;
  static constexpr double b_0 = -0.048;
  static constexpr double b_1 = -0.0269;
  static constexpr double c_0 = 0.0020191519406228;
  static constexpr double c_1 = 0.00069255121311694;
  static constexpr double d_0 = -0.0116320663789130;
  static constexpr double d_1 = -0.00480126353790614;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t76 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_0;
    constexpr double t14 = beta1_0;
    constexpr double t19 = beta2_0 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_0;
    constexpr double t32 = c_0 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_0 * t1;
    constexpr double t43 = gamma_1;
    constexpr double t44 = beta1_1;
    constexpr double t48 = beta2_1 * t1;
    constexpr double t54 = a_1;
    constexpr double t58 = c_1 * t1;
    constexpr double t59 = t58 * t3;
    constexpr double t63 = d_1 * t1;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = safe_math::sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = safe_math::log( t11 );
    const double t34 = t9 * t28;
    const double t42 = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_0 + t33 * t34 / 0.4e1 + t38 * t21 / 0.4e1 );
    const double t51 = 0.1e1 + t44 * t15 / 0.2e1 + t48 * t21 / 0.4e1;
    const double t67 = piecewise_functor_3( t12, t43 / t51, t54 * t28 + b_1 + t59 * t34 / 0.4e1 + t63 * t21 / 0.4e1 );
    const double t70 = safe_math::cbrt( zeta_tol );
    const double t72 = piecewise_functor_3( 0.1e1 <= zeta_tol, t70 * zeta_tol, 1.0 );
    const double t74 = 0.2e1 * t72 - 0.2e1;
    const double t79 = 0.1e1 / ( 0.2e1 * t76 - 0.2e1 );
    const double t80 = ( t67 - t42 ) * t74 * t79;


    eps = t42 + t80;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t76 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_0;
    constexpr double t14 = beta1_0;
    constexpr double t19 = beta2_0 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_0;
    constexpr double t32 = c_0 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_0 * t1;
    constexpr double t43 = gamma_1;
    constexpr double t44 = beta1_1;
    constexpr double t48 = beta2_1 * t1;
    constexpr double t54 = a_1;
    constexpr double t58 = c_1 * t1;
    constexpr double t59 = t58 * t3;
    constexpr double t63 = d_1 * t1;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = safe_math::sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = safe_math::log( t11 );
    const double t34 = t9 * t28;
    const double t42 = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_0 + t33 * t34 / 0.4e1 + t38 * t21 / 0.4e1 );
    const double t51 = 0.1e1 + t44 * t15 / 0.2e1 + t48 * t21 / 0.4e1;
    const double t67 = piecewise_functor_3( t12, t43 / t51, t54 * t28 + b_1 + t59 * t34 / 0.4e1 + t63 * t21 / 0.4e1 );
    const double t70 = safe_math::cbrt( zeta_tol );
    const double t72 = piecewise_functor_3( 0.1e1 <= zeta_tol, t70 * zeta_tol, 1.0 );
    const double t74 = 0.2e1 * t72 - 0.2e1;
    const double t79 = 0.1e1 / ( 0.2e1 * t76 - 0.2e1 );
    const double t80 = ( t67 - t42 ) * t74 * t79;
    const double t81 = t24 * t24;
    const double t83 = t13 / t81;
    const double t84 = 0.1e1 / t15;
    const double t86 = t14 * t84 * t1;
    const double t88 = 0.1e1 / t7 / rho;
    const double t89 = t20 * t88;
    const double t93 = -t19 * t89 / 0.12e2 - t86 * t89 / 0.12e2;
    const double t95 = 0.1e1 / rho;
    const double t99 = t6 * t88 * t28;
    const double t107 = piecewise_functor_3( t12, -t83 * t93, -t27 * t95 / 0.3e1 - t33 * t99 / 0.12e2 - t32 * t89 / 0.12e2 - t38 * t89 / 0.12e2 );
    const double t108 = t51 * t51;
    const double t110 = t43 / t108;
    const double t112 = t44 * t84 * t1;
    const double t116 = -t112 * t89 / 0.12e2 - t48 * t89 / 0.12e2;
    const double t127 = piecewise_functor_3( t12, -t110 * t116, -t54 * t95 / 0.3e1 - t59 * t99 / 0.12e2 - t58 * t89 / 0.12e2 - t63 * t89 / 0.12e2 );
    const double t130 = ( t127 - t107 ) * t74 * t79;


    eps = t42 + t80;
    vrho = t42 + t80 + rho * ( t107 + t130 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t87 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t14 = gamma_0;
    constexpr double t15 = beta1_0;
    constexpr double t20 = beta2_0 * t1;
    constexpr double t21 = t3 * t6;
    constexpr double t28 = a_0;
    constexpr double t33 = c_0 * t1;
    constexpr double t34 = t33 * t3;
    constexpr double t39 = d_0 * t1;
    constexpr double t44 = gamma_1;
    constexpr double t45 = beta1_1;
    constexpr double t49 = beta2_1 * t1;
    constexpr double t55 = a_1;
    constexpr double t59 = c_1 * t1;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = d_1 * t1;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t1 * t3 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = 0.1e1 <= t12;
    const double t16 = safe_math::sqrt( t11 );
    const double t22 = t21 * t9;
    const double t25 = 0.1e1 + t15 * t16 / 0.2e1 + t20 * t22 / 0.4e1;
    const double t29 = safe_math::log( t12 );
    const double t35 = t10 * t29;
    const double t43 = piecewise_functor_3( t13, t14 / t25, t28 * t29 + b_0 + t34 * t35 / 0.4e1 + t39 * t22 / 0.4e1 );
    const double t52 = 0.1e1 + t45 * t16 / 0.2e1 + t49 * t22 / 0.4e1;
    const double t68 = piecewise_functor_3( t13, t44 / t52, t55 * t29 + b_1 + t60 * t35 / 0.4e1 + t64 * t22 / 0.4e1 );
    const double t69 = t68 - t43;
    const double t70 = rho_a - rho_b;
    const double t71 = 0.1e1 / t7;
    const double t72 = t70 * t71;
    const double t73 = 0.1e1 + t72;
    const double t74 = t73 <= zeta_tol;
    const double t75 = safe_math::cbrt( zeta_tol );
    const double t76 = t75 * zeta_tol;
    const double t77 = safe_math::cbrt( t73 );
    const double t79 = piecewise_functor_3( t74, t76, t77 * t73 );
    const double t80 = 0.1e1 - t72;
    const double t81 = t80 <= zeta_tol;
    const double t82 = safe_math::cbrt( t80 );
    const double t84 = piecewise_functor_3( t81, t76, t82 * t80 );
    const double t85 = t79 + t84 - 0.2e1;
    const double t90 = 0.1e1 / ( 0.2e1 * t87 - 0.2e1 );
    const double t91 = t69 * t85 * t90;


    eps = t43 + t91;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t87 = constants::m_cbrt_2;
    constexpr double t6 = t5 * t5;
    constexpr double t14 = gamma_0;
    constexpr double t15 = beta1_0;
    constexpr double t20 = beta2_0 * t1;
    constexpr double t21 = t3 * t6;
    constexpr double t28 = a_0;
    constexpr double t33 = c_0 * t1;
    constexpr double t34 = t33 * t3;
    constexpr double t39 = d_0 * t1;
    constexpr double t44 = gamma_1;
    constexpr double t45 = beta1_1;
    constexpr double t49 = beta2_1 * t1;
    constexpr double t55 = a_1;
    constexpr double t59 = c_1 * t1;
    constexpr double t60 = t59 * t3;
    constexpr double t64 = d_1 * t1;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t1 * t3 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = 0.1e1 <= t12;
    const double t16 = safe_math::sqrt( t11 );
    const double t22 = t21 * t9;
    const double t25 = 0.1e1 + t15 * t16 / 0.2e1 + t20 * t22 / 0.4e1;
    const double t29 = safe_math::log( t12 );
    const double t35 = t10 * t29;
    const double t43 = piecewise_functor_3( t13, t14 / t25, t28 * t29 + b_0 + t34 * t35 / 0.4e1 + t39 * t22 / 0.4e1 );
    const double t52 = 0.1e1 + t45 * t16 / 0.2e1 + t49 * t22 / 0.4e1;
    const double t68 = piecewise_functor_3( t13, t44 / t52, t55 * t29 + b_1 + t60 * t35 / 0.4e1 + t64 * t22 / 0.4e1 );
    const double t69 = t68 - t43;
    const double t70 = rho_a - rho_b;
    const double t71 = 0.1e1 / t7;
    const double t72 = t70 * t71;
    const double t73 = 0.1e1 + t72;
    const double t74 = t73 <= zeta_tol;
    const double t75 = safe_math::cbrt( zeta_tol );
    const double t76 = t75 * zeta_tol;
    const double t77 = safe_math::cbrt( t73 );
    const double t79 = piecewise_functor_3( t74, t76, t77 * t73 );
    const double t80 = 0.1e1 - t72;
    const double t81 = t80 <= zeta_tol;
    const double t82 = safe_math::cbrt( t80 );
    const double t84 = piecewise_functor_3( t81, t76, t82 * t80 );
    const double t85 = t79 + t84 - 0.2e1;
    const double t90 = 0.1e1 / ( 0.2e1 * t87 - 0.2e1 );
    const double t91 = t69 * t85 * t90;
    const double t92 = t25 * t25;
    const double t94 = t14 / t92;
    const double t95 = 0.1e1 / t16;
    const double t97 = t15 * t95 * t1;
    const double t99 = 0.1e1 / t8 / t7;
    const double t100 = t21 * t99;
    const double t104 = -t20 * t100 / 0.12e2 - t97 * t100 / 0.12e2;
    const double t109 = t6 * t99 * t29;
    const double t117 = piecewise_functor_3( t13, -t94 * t104, -t28 * t71 / 0.3e1 - t34 * t109 / 0.12e2 - t33 * t100 / 0.12e2 - t39 * t100 / 0.12e2 );
    const double t118 = t52 * t52;
    const double t120 = t44 / t118;
    const double t122 = t45 * t95 * t1;
    const double t126 = -t122 * t100 / 0.12e2 - t49 * t100 / 0.12e2;
    const double t137 = piecewise_functor_3( t13, -t120 * t126, -t55 * t71 / 0.3e1 - t60 * t109 / 0.12e2 - t59 * t100 / 0.12e2 - t64 * t100 / 0.12e2 );
    const double t138 = t137 - t117;
    const double t140 = t138 * t85 * t90;
    const double t141 = t7 * t7;
    const double t142 = 0.1e1 / t141;
    const double t143 = t70 * t142;
    const double t144 = t71 - t143;
    const double t147 = piecewise_functor_3( t74, 0.0, 0.4e1 / 0.3e1 * t77 * t144 );
    const double t148 = -t144;
    const double t151 = piecewise_functor_3( t81, 0.0, 0.4e1 / 0.3e1 * t82 * t148 );
    const double t152 = t147 + t151;
    const double t154 = t69 * t152 * t90;
    const double t157 = -t71 - t143;
    const double t160 = piecewise_functor_3( t74, 0.0, 0.4e1 / 0.3e1 * t77 * t157 );
    const double t161 = -t157;
    const double t164 = piecewise_functor_3( t81, 0.0, 0.4e1 / 0.3e1 * t82 * t161 );
    const double t165 = t160 + t164;
    const double t167 = t69 * t165 * t90;


    eps = t43 + t91;
    vrho_a = t43 + t91 + t7 * ( t117 + t140 + t154 );
    vrho_b = t43 + t91 + t7 * ( t117 + t140 + t167 );

  }


};

struct BuiltinPZ81_MOD : detail::BuiltinKernelImpl< BuiltinPZ81_MOD > {

  BuiltinPZ81_MOD( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPZ81_MOD >(p) { }
  
  virtual ~BuiltinPZ81_MOD() = default;

};



} // namespace ExchCXX