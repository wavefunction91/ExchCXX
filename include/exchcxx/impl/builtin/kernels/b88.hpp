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
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1e-25;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double beta = 0.0042;
  static constexpr double gamma = 6.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t56 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t52 );


    eps = 0.2e1 * t56;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t23 = constants::m_cbrt_one_ov_pi;
    constexpr double t25 = constants::m_cbrt_4;
    constexpr double t28 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t20 = t3 * t3;
    constexpr double t21 = beta * t20;
    constexpr double t24 = 0.1e1 / t23;
    constexpr double t26 = t24 * t25;
    constexpr double t27 = t21 * t26;
    constexpr double t29 = t28 * t28;
    constexpr double t35 = gamma * beta;
    constexpr double t99 = t21 * t24;
    constexpr double t100 = t25 * t29;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t30 = sigma * t29;
    const double t31 = rho * rho;
    const double t32 = t18 * t18;
    const double t34 = 0.1e1 / t32 / t31;
    const double t36 = safe_math::sqrt( sigma );
    const double t37 = t35 * t36;
    const double t39 = 0.1e1 / t18 / rho;
    const double t43 = safe_math::log( t36 * t28 * t39 + safe_math::sqrt( square( t36 * t28 * t39 ) + 0.1e1 ) );
    const double t44 = t28 * t39 * t43;
    const double t46 = t37 * t44 + 0.1e1;
    const double t47 = 0.1e1 / t46;
    const double t48 = t34 * t47;
    const double t52 = 0.1e1 + 0.2e1 / 0.9e1 * t27 * t30 * t48;
    const double t56 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t52 );
    const double t58 = t17 / t32;
    const double t62 = t31 * rho;
    const double t64 = 0.1e1 / t32 / t62;
    const double t65 = t64 * t47;
    const double t69 = t46 * t46;
    const double t70 = 0.1e1 / t69;
    const double t71 = t34 * t70;
    const double t75 = t28 / t18 / t31 * t43;
    const double t77 = t35 * sigma;
    const double t78 = t29 * t64;
    const double t80 = t30 * t34 + 0.1e1;
    const double t81 = safe_math::sqrt( t80 );
    const double t82 = 0.1e1 / t81;
    const double t83 = t78 * t82;
    const double t86 = -0.4e1 / 0.3e1 * t37 * t75 - 0.4e1 / 0.3e1 * t77 * t83;
    const double t91 = -0.16e2 / 0.27e2 * t27 * t30 * t65 - 0.2e1 / 0.9e1 * t27 * t30 * t71 * t86;
    const double t96 = piecewise_functor_3( t2, 0.0, -t6 * t58 * t52 / 0.8e1 - 0.3e1 / 0.8e1 * t6 * t19 * t91 );
    const double t104 = t35 / t36;
    const double t106 = t29 * t34;
    const double t107 = t106 * t82;
    const double t110 = t104 * t44 / 0.2e1 + t35 * t107 / 0.2e1;
    const double t115 = -0.2e1 / 0.9e1 * t27 * t30 * t71 * t110 + 0.2e1 / 0.9e1 * t99 * t100 * t48;
    const double t119 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t115 );


    eps = 0.2e1 * t56;
    vrho = 0.2e1 * rho * t96 + 0.2e1 * t56;
    vsigma = 0.2e1 * rho * t119;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t59 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t55 );
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t93 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t89 );


    eps = t59 + t93;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t31 = constants::m_cbrt_one_ov_pi;
    constexpr double t34 = constants::m_cbrt_4;
    constexpr double t5 = t2 / t3;
    constexpr double t28 = t2 * t2;
    constexpr double t29 = beta * t28;
    constexpr double t32 = 0.1e1 / t31;
    constexpr double t33 = t29 * t32;
    constexpr double t41 = gamma * beta;
    constexpr double t119 = t32 * t34;
    constexpr double t120 = t29 * t119;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t35 = t34 * sigma_aa;
    const double t36 = rho_a * rho_a;
    const double t37 = safe_math::cbrt( rho_a );
    const double t38 = t37 * t37;
    const double t40 = 0.1e1 / t38 / t36;
    const double t42 = safe_math::sqrt( sigma_aa );
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t45 = t42 * t44;
    const double t46 = safe_math::log( t45 + safe_math::sqrt( t45 * t45 + 0.1e1 ) );
    const double t49 = t41 * t45 * t46 + 0.1e1;
    const double t50 = 0.1e1 / t49;
    const double t55 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t35 * t40 * t50;
    const double t59 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t55 );
    const double t60 = rho_b <= dens_tol;
    const double t61 = -t16;
    const double t63 = piecewise_functor_5( t14, t11, t10, t15, t61 * t7 );
    const double t64 = 0.1e1 + t63;
    const double t65 = t64 <= zeta_tol;
    const double t66 = safe_math::cbrt( t64 );
    const double t68 = piecewise_functor_3( t65, t22, t66 * t64 );
    const double t69 = t68 * t26;
    const double t70 = t34 * sigma_bb;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = safe_math::sqrt( sigma_bb );
    const double t78 = 0.1e1 / t72 / rho_b;
    const double t79 = t76 * t78;
    const double t80 = safe_math::log( t79 + safe_math::sqrt( t79 * t79 + 0.1e1 ) );
    const double t83 = t41 * t79 * t80 + 0.1e1;
    const double t84 = 0.1e1 / t83;
    const double t89 = 0.1e1 + 0.2e1 / 0.9e1 * t33 * t70 * t75 * t84;
    const double t93 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t89 );
    const double t94 = t6 * t6;
    const double t95 = 0.1e1 / t94;
    const double t96 = t16 * t95;
    const double t98 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t96 );
    const double t101 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t98 );
    const double t102 = t101 * t26;
    const double t106 = t26 * t26;
    const double t107 = 0.1e1 / t106;
    const double t108 = t25 * t107;
    const double t111 = t5 * t108 * t55 / 0.8e1;
    const double t112 = t36 * rho_a;
    const double t114 = 0.1e1 / t38 / t112;
    const double t121 = sigma_aa * t40;
    const double t122 = t49 * t49;
    const double t123 = 0.1e1 / t122;
    const double t125 = 0.1e1 / t37 / t36;
    const double t129 = sigma_aa * t114;
    const double t130 = t121 + 0.1e1;
    const double t131 = safe_math::sqrt( t130 );
    const double t132 = 0.1e1 / t131;
    const double t136 = -0.4e1 / 0.3e1 * t41 * t42 * t125 * t46 - 0.4e1 / 0.3e1 * t41 * t129 * t132;
    const double t137 = t123 * t136;
    const double t141 = -0.16e2 / 0.27e2 * t33 * t35 * t114 * t50 - 0.2e1 / 0.9e1 * t120 * t121 * t137;
    const double t146 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t102 * t55 - t111 - 0.3e1 / 0.8e1 * t5 * t27 * t141 );
    const double t147 = t61 * t95;
    const double t149 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t147 );
    const double t152 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t149 );
    const double t153 = t152 * t26;
    const double t157 = t68 * t107;
    const double t160 = t5 * t157 * t89 / 0.8e1;
    const double t162 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t153 * t89 - t160 );
    const double t166 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t96 );
    const double t169 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t166 );
    const double t170 = t169 * t26;
    const double t175 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t170 * t55 - t111 );
    const double t177 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t147 );
    const double t180 = piecewise_functor_3( t65, 0.0, 0.4e1 / 0.3e1 * t66 * t177 );
    const double t181 = t180 * t26;
    const double t185 = t71 * rho_b;
    const double t187 = 0.1e1 / t73 / t185;
    const double t192 = sigma_bb * t75;
    const double t193 = t83 * t83;
    const double t194 = 0.1e1 / t193;
    const double t196 = 0.1e1 / t72 / t71;
    const double t200 = sigma_bb * t187;
    const double t201 = t192 + 0.1e1;
    const double t202 = safe_math::sqrt( t201 );
    const double t203 = 0.1e1 / t202;
    const double t207 = -0.4e1 / 0.3e1 * t41 * t76 * t196 * t80 - 0.4e1 / 0.3e1 * t41 * t200 * t203;
    const double t208 = t194 * t207;
    const double t212 = -0.16e2 / 0.27e2 * t33 * t70 * t187 * t84 - 0.2e1 / 0.9e1 * t120 * t192 * t208;
    const double t217 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t181 * t89 - t160 - 0.3e1 / 0.8e1 * t5 * t69 * t212 );
    const double t220 = t34 * t40;
    const double t223 = 0.1e1 / t42;
    const double t230 = t41 * t223 * t44 * t46 / 0.2e1 + t41 * t40 * t132 / 0.2e1;
    const double t231 = t123 * t230;
    const double t235 = -0.2e1 / 0.9e1 * t120 * t121 * t231 + 0.2e1 / 0.9e1 * t33 * t220 * t50;
    const double t239 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t235 );
    const double t240 = t34 * t75;
    const double t243 = 0.1e1 / t76;
    const double t250 = t41 * t243 * t78 * t80 / 0.2e1 + t41 * t75 * t203 / 0.2e1;
    const double t251 = t194 * t250;
    const double t255 = -0.2e1 / 0.9e1 * t120 * t192 * t251 + 0.2e1 / 0.9e1 * t33 * t240 * t84;
    const double t259 = piecewise_functor_3( t60, 0.0, -0.3e1 / 0.8e1 * t5 * t69 * t255 );


    eps = t59 + t93;
    vrho_a = t59 + t93 + t6 * ( t146 + t162 );
    vrho_b = t59 + t93 + t6 * ( t175 + t217 );
    vsigma_aa = t6 * t239;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t259;

  }


};

struct BuiltinB88 : detail::BuiltinKernelImpl< BuiltinB88 > {

  BuiltinB88( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinB88 >(p) { }
  
  virtual ~BuiltinB88() = default;

};



} // namespace ExchCXX