#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinVWN_RPA > :
  public lda_screening_interface< BuiltinVWN_RPA > {

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



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t21 = t12 + 0.13072e2;
    const double t24 = safe_math::atan( 0.44899888641287296627e-1 / t21 );
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t45 = ( 0.310907e-1 * t19 + 0.20521972937837502661e2 * t24 + 0.44313737677495382697e-2 * t30 ) * t44;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t52 = safe_math::log( t4 * t9 * t48 / 0.4e1 );
    const double t54 = t12 + 0.201231e2;
    const double t57 = safe_math::atan( 0.11716852777089929792e1 / t54 );
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t62 = safe_math::log( t60 * t48 );
    const double t66 = ( 0.1554535e-1 * t52 + 0.61881802979060631482e0 * t57 + 0.26673100072733151594e-2 * t62 ) * t38 * t42;


    eps = t45 + t66;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t39 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t73 = t4 * t6;
    constexpr double t81 = t3 * t6;
    constexpr double t90 = t1 * t1;
    constexpr double t92 = 0.1e1 / t3;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = safe_math::sqrt( t10 );
    const double t14 = t11 + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = safe_math::log( t4 * t9 * t15 / 0.4e1 );
    const double t21 = t12 + 0.13072e2;
    const double t24 = safe_math::atan( 0.44899888641287296627e-1 / t21 );
    const double t26 = t12 / 0.2e1;
    const double t27 = t26 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = safe_math::log( t28 * t15 );
    const double t34 = safe_math::cbrt( zeta_tol );
    const double t36 = piecewise_functor_3( 0.1e1 <= zeta_tol, t34 * zeta_tol, 1.0 );
    const double t38 = 0.2e1 * t36 - 0.2e1;
    const double t42 = 0.1e1 / ( 0.2e1 * t39 - 0.2e1 );
    const double t44 = -t38 * t42 + 0.1e1;
    const double t45 = ( 0.310907e-1 * t19 + 0.20521972937837502661e2 * t24 + 0.44313737677495382697e-2 * t30 ) * t44;
    const double t47 = t11 + 0.1006155e2 * t12 + 0.101578e3;
    const double t48 = 0.1e1 / t47;
    const double t52 = safe_math::log( t4 * t9 * t48 / 0.4e1 );
    const double t54 = t12 + 0.201231e2;
    const double t57 = safe_math::atan( 0.11716852777089929792e1 / t54 );
    const double t59 = t26 + 0.743294e0;
    const double t60 = t59 * t59;
    const double t62 = safe_math::log( t60 * t48 );
    const double t66 = ( 0.1554535e-1 * t52 + 0.61881802979060631482e0 * t57 + 0.26673100072733151594e-2 * t62 ) * t38 * t42;
    const double t68 = 0.1e1 / t7 / rho;
    const double t69 = t6 * t68;
    const double t74 = t14 * t14;
    const double t75 = 0.1e1 / t74;
    const double t76 = t8 * t75;
    const double t77 = t4 * t69;
    const double t78 = t77 / 0.12e2;
    const double t79 = 0.1e1 / t12;
    const double t80 = t79 * t1;
    const double t83 = t80 * t81 * t68;
    const double t85 = -t78 - 0.10893333333333333333e1 * t83;
    const double t93 = ( -t4 * t69 * t15 / 0.12e2 - t73 * t76 * t85 / 0.4e1 ) * t90 * t92;
    const double t94 = t5 * t7;
    const double t95 = t94 * t14;
    const double t98 = t21 * t21;
    const double t99 = 0.1e1 / t98;
    const double t101 = t99 * t79 * t1;
    const double t103 = 0.2016e-2 * t99 + 0.1e1;
    const double t104 = 0.1e1 / t103;
    const double t109 = t27 * t15;
    const double t110 = t109 * t79;
    const double t113 = t28 * t75;
    const double t115 = -t110 * t77 / 0.6e1 - t113 * t85;
    const double t116 = 0.1e1 / t28;
    const double t117 = t115 * t116;
    const double t121 = ( 0.10363566666666666667e-1 * t93 * t95 + 0.15357238326806922974e0 * t101 * t81 * t68 * t104 + 0.44313737677495382697e-2 * t117 * t14 ) * t44;
    const double t125 = t47 * t47;
    const double t126 = 0.1e1 / t125;
    const double t127 = t8 * t126;
    const double t129 = -t78 - 0.1676925e1 * t83;
    const double t135 = ( -t4 * t69 * t48 / 0.12e2 - t73 * t127 * t129 / 0.4e1 ) * t90 * t92;
    const double t136 = t94 * t47;
    const double t139 = t54 * t54;
    const double t140 = 0.1e1 / t139;
    const double t142 = t140 * t79 * t1;
    const double t144 = 0.137284639e1 * t140 + 0.1e1;
    const double t145 = 0.1e1 / t144;
    const double t150 = t59 * t48;
    const double t151 = t150 * t79;
    const double t154 = t60 * t126;
    const double t156 = -t151 * t77 / 0.6e1 - t154 * t129;
    const double t157 = 0.1e1 / t60;
    const double t158 = t156 * t157;
    const double t163 = ( 0.51817833333333333333e-2 * t135 * t136 + 0.12084332918108974175e0 * t142 * t81 * t68 * t145 + 0.26673100072733151594e-2 * t158 * t47 ) * t38 * t42;


    eps = t45 + t66;
    vrho = t45 + t66 + rho * ( t121 + t163 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t56 = t33 * t55;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t77 = t75 * t49 * t53;


    eps = t56 + t77;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t50 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t84 = t4 * t6;
    constexpr double t92 = t3 * t6;
    constexpr double t101 = t1 * t1;
    constexpr double t103 = 0.1e1 / t3;


    const double t7 = rho_a + rho_b;
    const double t8 = safe_math::cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = safe_math::sqrt( t11 );
    const double t15 = t12 + 0.6536e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = safe_math::log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.13072e2;
    const double t25 = safe_math::atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = safe_math::log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = t37 <= zeta_tol;
    const double t39 = safe_math::cbrt( zeta_tol );
    const double t40 = t39 * zeta_tol;
    const double t41 = safe_math::cbrt( t37 );
    const double t43 = piecewise_functor_3( t38, t40, t41 * t37 );
    const double t44 = 0.1e1 - t36;
    const double t45 = t44 <= zeta_tol;
    const double t46 = safe_math::cbrt( t44 );
    const double t48 = piecewise_functor_3( t45, t40, t46 * t44 );
    const double t49 = t43 + t48 - 0.2e1;
    const double t53 = 0.1e1 / ( 0.2e1 * t50 - 0.2e1 );
    const double t55 = -t49 * t53 + 0.1e1;
    const double t56 = t33 * t55;
    const double t58 = t12 + 0.1006155e2 * t13 + 0.101578e3;
    const double t59 = 0.1e1 / t58;
    const double t63 = safe_math::log( t4 * t10 * t59 / 0.4e1 );
    const double t65 = t13 + 0.201231e2;
    const double t68 = safe_math::atan( 0.11716852777089929792e1 / t65 );
    const double t70 = t27 + 0.743294e0;
    const double t71 = t70 * t70;
    const double t73 = safe_math::log( t71 * t59 );
    const double t75 = 0.1554535e-1 * t63 + 0.61881802979060631482e0 * t68 + 0.26673100072733151594e-2 * t73;
    const double t77 = t75 * t49 * t53;
    const double t79 = 0.1e1 / t8 / t7;
    const double t80 = t6 * t79;
    const double t85 = t15 * t15;
    const double t86 = 0.1e1 / t85;
    const double t87 = t9 * t86;
    const double t88 = t4 * t80;
    const double t89 = t88 / 0.12e2;
    const double t90 = 0.1e1 / t13;
    const double t91 = t90 * t1;
    const double t94 = t91 * t92 * t79;
    const double t96 = -t89 - 0.10893333333333333333e1 * t94;
    const double t104 = ( -t4 * t80 * t16 / 0.12e2 - t84 * t87 * t96 / 0.4e1 ) * t101 * t103;
    const double t105 = t5 * t8;
    const double t106 = t105 * t15;
    const double t109 = t22 * t22;
    const double t110 = 0.1e1 / t109;
    const double t112 = t110 * t90 * t1;
    const double t114 = 0.2016e-2 * t110 + 0.1e1;
    const double t115 = 0.1e1 / t114;
    const double t120 = t28 * t16;
    const double t121 = t120 * t90;
    const double t124 = t29 * t86;
    const double t126 = -t121 * t88 / 0.6e1 - t124 * t96;
    const double t127 = 0.1e1 / t29;
    const double t128 = t126 * t127;
    const double t131 = 0.10363566666666666667e-1 * t104 * t106 + 0.15357238326806922974e0 * t112 * t92 * t79 * t115 + 0.44313737677495382697e-2 * t128 * t15;
    const double t132 = t131 * t55;
    const double t133 = t7 * t7;
    const double t134 = 0.1e1 / t133;
    const double t135 = t34 * t134;
    const double t136 = t35 - t135;
    const double t139 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t136 );
    const double t140 = -t136;
    const double t143 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t140 );
    const double t144 = t139 + t143;
    const double t146 = t33 * t144 * t53;
    const double t150 = t58 * t58;
    const double t151 = 0.1e1 / t150;
    const double t152 = t9 * t151;
    const double t154 = -t89 - 0.1676925e1 * t94;
    const double t160 = ( -t4 * t80 * t59 / 0.12e2 - t84 * t152 * t154 / 0.4e1 ) * t101 * t103;
    const double t161 = t105 * t58;
    const double t164 = t65 * t65;
    const double t165 = 0.1e1 / t164;
    const double t167 = t165 * t90 * t1;
    const double t169 = 0.137284639e1 * t165 + 0.1e1;
    const double t170 = 0.1e1 / t169;
    const double t175 = t70 * t59;
    const double t176 = t175 * t90;
    const double t179 = t71 * t151;
    const double t181 = -t176 * t88 / 0.6e1 - t179 * t154;
    const double t182 = 0.1e1 / t71;
    const double t183 = t181 * t182;
    const double t186 = 0.51817833333333333333e-2 * t160 * t161 + 0.12084332918108974175e0 * t167 * t92 * t79 * t170 + 0.26673100072733151594e-2 * t183 * t58;
    const double t188 = t186 * t49 * t53;
    const double t190 = t75 * t144 * t53;
    const double t193 = -t35 - t135;
    const double t196 = piecewise_functor_3( t38, 0.0, 0.4e1 / 0.3e1 * t41 * t193 );
    const double t197 = -t193;
    const double t200 = piecewise_functor_3( t45, 0.0, 0.4e1 / 0.3e1 * t46 * t197 );
    const double t201 = t196 + t200;
    const double t203 = t33 * t201 * t53;
    const double t205 = t75 * t201 * t53;


    eps = t56 + t77;
    vrho_a = t56 + t77 + t7 * ( t132 - t146 + t188 + t190 );
    vrho_b = t56 + t77 + t7 * ( t132 - t203 + t188 + t205 );

  }


};

struct BuiltinVWN_RPA : detail::BuiltinKernelImpl< BuiltinVWN_RPA > {

  BuiltinVWN_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN_RPA >(p) { }
  
  virtual ~BuiltinVWN_RPA() = default;

};



} // namespace ExchCXX