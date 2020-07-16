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
    const double t14 = t10 / 0.4e1 + 0.65360000000000000000e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.130720e2;
    const double t24 = atan( 0.44899888641287296627e-1 / t21 );
    const double t25 = 0.20521972937837502661e2 * t24;
    const double t27 = t12 / 0.2e1 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.44313737677495382697e-2 * t30;


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
    const double t14 = t10 / 0.4e1 + 0.65360000000000000000e1 * t12 + 0.427198e2;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.130720e2;
    const double t24 = atan( 0.44899888641287296627e-1 / t21 );
    const double t25 = 0.20521972937837502661e2 * t24;
    const double t27 = t12 / 0.2e1 + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.44313737677495382697e-2 * t30;
    const double t33 = 0.1e1 / t7 / rho;
    const double t34 = t6 * t33;
    const double t39 = t14 * t14;
    const double t40 = 0.1e1 / t39;
    const double t41 = t8 * t40;
    const double t42 = t4 * t34;
    const double t44 = 0.1e1 / t12;
    const double t45 = t44 * t1;
    const double t50 = -t42 / 0.12e2 - 0.10893333333333333333e1 * t45 * t46 * t33;
    const double t58 = ( -t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1 ) * t55 * t57;
    const double t59 = t5 * t7;
    const double t60 = t59 * t14;
    const double t61 = t58 * t60;
    const double t63 = t21 * t21;
    const double t64 = 0.1e1 / t63;
    const double t66 = t64 * t44 * t1;
    const double t68 = 0.20160000000000000000e-2 * t64 + 0.1e1;
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
    vrho = t20 + t25 + t31 + rho * ( 0.10363566666666666667e-1 * t61 + 0.15357238326806922974e0 * t72 + 0.44313737677495382697e-2 * t83 );

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
    const double t14 = t10 / 0.4e1 + 0.10061550000000000000e2 * t12 + 0.101578e3;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.1554535e-1 * t19;
    const double t21 = t12 + 0.201231e2;
    const double t24 = atan( 0.11716852777089929792e1 / t21 );
    const double t25 = 0.61881802979060631482e0 * t24;
    const double t27 = t12 / 0.2e1 + 0.743294e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.26673100072733151594e-2 * t30;


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
    const double t14 = t10 / 0.4e1 + 0.10061550000000000000e2 * t12 + 0.101578e3;
    const double t15 = 0.1e1 / t14;
    const double t19 = log( t4 * t9 * t15 / 0.4e1 );
    const double t20 = 0.1554535e-1 * t19;
    const double t21 = t12 + 0.201231e2;
    const double t24 = atan( 0.11716852777089929792e1 / t21 );
    const double t25 = 0.61881802979060631482e0 * t24;
    const double t27 = t12 / 0.2e1 + 0.743294e0;
    const double t28 = t27 * t27;
    const double t30 = log( t28 * t15 );
    const double t31 = 0.26673100072733151594e-2 * t30;
    const double t33 = 0.1e1 / t7 / rho;
    const double t34 = t6 * t33;
    const double t39 = t14 * t14;
    const double t40 = 0.1e1 / t39;
    const double t41 = t8 * t40;
    const double t42 = t4 * t34;
    const double t44 = 0.1e1 / t12;
    const double t45 = t44 * t1;
    const double t50 = -t42 / 0.12e2 - 0.16769250000000000000e1 * t45 * t46 * t33;
    const double t58 = ( -t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1 ) * t55 * t57;
    const double t59 = t5 * t7;
    const double t60 = t59 * t14;
    const double t61 = t58 * t60;
    const double t63 = t21 * t21;
    const double t64 = 0.1e1 / t63;
    const double t66 = t64 * t44 * t1;
    const double t68 = 0.13728463900000000000e1 * t64 + 0.1e1;
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
    vrho = t20 + t25 + t31 + rho * ( 0.51817833333333333333e-2 * t61 + 0.12084332918108974175e0 * t72 + 0.26673100072733151594e-2 * t83 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t44 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;


    const double t7 = rho_a + rho_b;
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = sqrt( t11 );
    const double t15 = t12 + 0.65360000000000000000e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.130720e2;
    const double t25 = atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = cbrt( t37 );
    const double t40 = 0.1e1 - t36;
    const double t41 = cbrt( t40 );
    const double t43 = t38 * t37 + t41 * t40 - 0.2e1;
    const double t47 = 0.1e1 / ( 0.2e1 * t44 - 0.2e1 );
    const double t49 = -t43 * t47 + 0.1e1;
    const double t50 = t33 * t49;
    const double t52 = t12 + 0.10061550000000000000e2 * t13 + 0.101578e3;
    const double t53 = 0.1e1 / t52;
    const double t57 = log( t4 * t10 * t53 / 0.4e1 );
    const double t59 = t13 + 0.201231e2;
    const double t62 = atan( 0.11716852777089929792e1 / t59 );
    const double t64 = t27 + 0.743294e0;
    const double t65 = t64 * t64;
    const double t67 = log( t65 * t53 );
    const double t69 = 0.1554535e-1 * t57 + 0.61881802979060631482e0 * t62 + 0.26673100072733151594e-2 * t67;
    const double t71 = t69 * t43 * t47;


    eps = t50 + t71;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t44 = constants::m_cbrt_2;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t78 = t4 * t6;
    constexpr double t86 = t3 * t6;
    constexpr double t95 = t1 * t1;
    constexpr double t97 = 0.1e1 / t3;


    const double t7 = rho_a + rho_b;
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t4 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = sqrt( t11 );
    const double t15 = t12 + 0.65360000000000000000e1 * t13 + 0.427198e2;
    const double t16 = 0.1e1 / t15;
    const double t20 = log( t4 * t10 * t16 / 0.4e1 );
    const double t22 = t13 + 0.130720e2;
    const double t25 = atan( 0.44899888641287296627e-1 / t22 );
    const double t27 = t13 / 0.2e1;
    const double t28 = t27 + 0.409286e0;
    const double t29 = t28 * t28;
    const double t31 = log( t29 * t16 );
    const double t33 = 0.310907e-1 * t20 + 0.20521972937837502661e2 * t25 + 0.44313737677495382697e-2 * t31;
    const double t34 = rho_a - rho_b;
    const double t35 = 0.1e1 / t7;
    const double t36 = t34 * t35;
    const double t37 = 0.1e1 + t36;
    const double t38 = cbrt( t37 );
    const double t40 = 0.1e1 - t36;
    const double t41 = cbrt( t40 );
    const double t43 = t38 * t37 + t41 * t40 - 0.2e1;
    const double t47 = 0.1e1 / ( 0.2e1 * t44 - 0.2e1 );
    const double t49 = -t43 * t47 + 0.1e1;
    const double t50 = t33 * t49;
    const double t52 = t12 + 0.10061550000000000000e2 * t13 + 0.101578e3;
    const double t53 = 0.1e1 / t52;
    const double t57 = log( t4 * t10 * t53 / 0.4e1 );
    const double t59 = t13 + 0.201231e2;
    const double t62 = atan( 0.11716852777089929792e1 / t59 );
    const double t64 = t27 + 0.743294e0;
    const double t65 = t64 * t64;
    const double t67 = log( t65 * t53 );
    const double t69 = 0.1554535e-1 * t57 + 0.61881802979060631482e0 * t62 + 0.26673100072733151594e-2 * t67;
    const double t71 = t69 * t43 * t47;
    const double t73 = 0.1e1 / t8 / t7;
    const double t74 = t6 * t73;
    const double t79 = t15 * t15;
    const double t80 = 0.1e1 / t79;
    const double t81 = t9 * t80;
    const double t82 = t4 * t74;
    const double t83 = t82 / 0.12e2;
    const double t84 = 0.1e1 / t13;
    const double t85 = t84 * t1;
    const double t88 = t85 * t86 * t73;
    const double t90 = -t83 - 0.10893333333333333333e1 * t88;
    const double t98 = ( -t4 * t74 * t16 / 0.12e2 - t78 * t81 * t90 / 0.4e1 ) * t95 * t97;
    const double t99 = t5 * t8;
    const double t100 = t99 * t15;
    const double t103 = t22 * t22;
    const double t104 = 0.1e1 / t103;
    const double t106 = t104 * t84 * t1;
    const double t108 = 0.20160000000000000000e-2 * t104 + 0.1e1;
    const double t109 = 0.1e1 / t108;
    const double t114 = t28 * t16;
    const double t115 = t114 * t84;
    const double t118 = t29 * t80;
    const double t120 = -t115 * t82 / 0.6e1 - t118 * t90;
    const double t121 = 0.1e1 / t29;
    const double t122 = t120 * t121;
    const double t125 = 0.10363566666666666667e-1 * t98 * t100 + 0.15357238326806922974e0 * t106 * t86 * t73 * t109 + 0.44313737677495382697e-2 * t122 * t15;
    const double t126 = t125 * t49;
    const double t127 = t7 * t7;
    const double t128 = 0.1e1 / t127;
    const double t129 = t34 * t128;
    const double t130 = t35 - t129;
    const double t132 = -t130;
    const double t135 = 0.4e1 / 0.3e1 * t38 * t130 + 0.4e1 / 0.3e1 * t41 * t132;
    const double t137 = t33 * t135 * t47;
    const double t141 = t52 * t52;
    const double t142 = 0.1e1 / t141;
    const double t143 = t9 * t142;
    const double t145 = -t83 - 0.16769250000000000000e1 * t88;
    const double t151 = ( -t4 * t74 * t53 / 0.12e2 - t78 * t143 * t145 / 0.4e1 ) * t95 * t97;
    const double t152 = t99 * t52;
    const double t155 = t59 * t59;
    const double t156 = 0.1e1 / t155;
    const double t158 = t156 * t84 * t1;
    const double t160 = 0.13728463900000000000e1 * t156 + 0.1e1;
    const double t161 = 0.1e1 / t160;
    const double t166 = t64 * t53;
    const double t167 = t166 * t84;
    const double t170 = t65 * t142;
    const double t172 = -t167 * t82 / 0.6e1 - t170 * t145;
    const double t173 = 0.1e1 / t65;
    const double t174 = t172 * t173;
    const double t177 = 0.51817833333333333333e-2 * t151 * t152 + 0.12084332918108974175e0 * t158 * t86 * t73 * t161 + 0.26673100072733151594e-2 * t174 * t52;
    const double t179 = t177 * t43 * t47;
    const double t181 = t69 * t135 * t47;
    const double t184 = -t35 - t129;
    const double t186 = -t184;
    const double t189 = 0.4e1 / 0.3e1 * t38 * t184 + 0.4e1 / 0.3e1 * t41 * t186;
    const double t191 = t33 * t189 * t47;
    const double t193 = t69 * t189 * t47;


    eps = t50 + t71;
    vrho_a = t50 + t71 + t7 * ( t126 - t137 + t179 + t181 );
    vrho_b = t50 + t71 + t7 * ( t126 - t191 + t179 + t193 );

  }


};

struct BuiltinVWN_RPA : detail::BuiltinKernelImpl< BuiltinVWN_RPA > {

  BuiltinVWN_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN_RPA >(p) { }
  
  virtual ~BuiltinVWN_RPA() = default;

};



} // namespace ExchCXX