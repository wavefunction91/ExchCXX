#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPZ81 > :
  public lda_screening_interface< BuiltinPZ81 > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;

  static constexpr double dens_tol  = 1e-24;

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
  static constexpr double c_0 = 0.0020;
  static constexpr double c_1 = 0.0007;
  static constexpr double d_0 = -0.0116;
  static constexpr double d_1 = -0.0048;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_0;
    constexpr double t14 = beta1_0;
    constexpr double t19 = beta2_0 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_0;
    constexpr double t32 = c_0 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_0 * t1;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = log( t11 );


    eps = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_0 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_0;
    constexpr double t14 = beta1_0;
    constexpr double t19 = beta2_0 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_0;
    constexpr double t32 = c_0 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_0 * t1;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = log( t11 );
    const double t42 = t24 * t24;
    const double t44 = t13 / t42;
    const double t47 = t14 / t15 * t1;
    const double t49 = 0.1e1 / t7 / rho;
    const double t50 = t20 * t49;
    const double t54 = -t19 * t50 / 0.12e2 - t47 * t50 / 0.12e2;
    const double t56 = 0.1e1 / rho;
    const double t68 = piecewise_functor_3( t12, -t44 * t54, -t27 * t56 / 0.3e1 - t33 * t6 * t49 * t28 / 0.12e2 - t32 * t50 / 0.12e2 - t38 * t50 / 0.12e2 );


    eps = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_0 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 );
    vrho = rho * t68 + ( piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_0 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 ) );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_1;
    constexpr double t14 = beta1_1;
    constexpr double t19 = beta2_1 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_1;
    constexpr double t32 = c_1 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_1 * t1;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = log( t11 );


    eps = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_1 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double& eps, double& vrho ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t6 = t5 * t5;
    constexpr double t13 = gamma_1;
    constexpr double t14 = beta1_1;
    constexpr double t19 = beta2_1 * t1;
    constexpr double t20 = t3 * t6;
    constexpr double t27 = a_1;
    constexpr double t32 = c_1 * t1;
    constexpr double t33 = t32 * t3;
    constexpr double t38 = d_1 * t1;


    const double t7 = cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t9 = t6 * t8;
    const double t10 = t1 * t3 * t9;
    const double t11 = t10 / 0.4e1;
    const double t12 = 0.1e1 <= t11;
    const double t15 = sqrt( t10 );
    const double t21 = t20 * t8;
    const double t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
    const double t28 = log( t11 );
    const double t42 = t24 * t24;
    const double t44 = t13 / t42;
    const double t47 = t14 / t15 * t1;
    const double t49 = 0.1e1 / t7 / rho;
    const double t50 = t20 * t49;
    const double t54 = -t19 * t50 / 0.12e2 - t47 * t50 / 0.12e2;
    const double t56 = 0.1e1 / rho;
    const double t68 = piecewise_functor_3( t12, -t44 * t54, -t27 * t56 / 0.3e1 - t33 * t6 * t49 * t28 / 0.12e2 - t32 * t50 / 0.12e2 - t38 * t50 / 0.12e2 );


    eps = piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_1 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 );
    vrho = rho * t68 + ( piecewise_functor_3( t12, t13 / t24, t27 * t28 + b_1 + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1 ) );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t81 = constants::m_cbrt_2;
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
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t1 * t3 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = 0.1e1 <= t12;
    const double t16 = sqrt( t11 );
    const double t22 = t21 * t9;
    const double t25 = 0.1e1 + t15 * t16 / 0.2e1 + t20 * t22 / 0.4e1;
    const double t29 = log( t12 );
    const double t35 = t10 * t29;
    const double t43 = piecewise_functor_3( t13, t14 / t25, t28 * t29 + b_0 + t34 * t35 / 0.4e1 + t39 * t22 / 0.4e1 );
    const double t52 = 0.1e1 + t45 * t16 / 0.2e1 + t49 * t22 / 0.4e1;
    const double t68 = piecewise_functor_3( t13, t44 / t52, t55 * t29 + b_1 + t60 * t35 / 0.4e1 + t64 * t22 / 0.4e1 );
    const double t69 = t68 - t43;
    const double t70 = rho_a - rho_b;
    const double t71 = 0.1e1 / t7;
    const double t72 = t70 * t71;
    const double t73 = 0.1e1 + t72;
    const double t74 = cbrt( t73 );
    const double t76 = 0.1e1 - t72;
    const double t77 = cbrt( t76 );
    const double t79 = t74 * t73 + t77 * t76 - 0.2e1;
    const double t84 = 0.1e1 / ( 0.2e1 * t81 - 0.2e1 );
    const double t85 = t69 * t79 * t84;


    eps = t43 + t85;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t2 = constants::m_one_ov_pi;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t81 = constants::m_cbrt_2;
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
    const double t8 = cbrt( t7 );
    const double t9 = 0.1e1 / t8;
    const double t10 = t6 * t9;
    const double t11 = t1 * t3 * t10;
    const double t12 = t11 / 0.4e1;
    const double t13 = 0.1e1 <= t12;
    const double t16 = sqrt( t11 );
    const double t22 = t21 * t9;
    const double t25 = 0.1e1 + t15 * t16 / 0.2e1 + t20 * t22 / 0.4e1;
    const double t29 = log( t12 );
    const double t35 = t10 * t29;
    const double t43 = piecewise_functor_3( t13, t14 / t25, t28 * t29 + b_0 + t34 * t35 / 0.4e1 + t39 * t22 / 0.4e1 );
    const double t52 = 0.1e1 + t45 * t16 / 0.2e1 + t49 * t22 / 0.4e1;
    const double t68 = piecewise_functor_3( t13, t44 / t52, t55 * t29 + b_1 + t60 * t35 / 0.4e1 + t64 * t22 / 0.4e1 );
    const double t69 = t68 - t43;
    const double t70 = rho_a - rho_b;
    const double t71 = 0.1e1 / t7;
    const double t72 = t70 * t71;
    const double t73 = 0.1e1 + t72;
    const double t74 = cbrt( t73 );
    const double t76 = 0.1e1 - t72;
    const double t77 = cbrt( t76 );
    const double t79 = t74 * t73 + t77 * t76 - 0.2e1;
    const double t84 = 0.1e1 / ( 0.2e1 * t81 - 0.2e1 );
    const double t85 = t69 * t79 * t84;
    const double t86 = t25 * t25;
    const double t88 = t14 / t86;
    const double t89 = 0.1e1 / t16;
    const double t91 = t15 * t89 * t1;
    const double t93 = 0.1e1 / t8 / t7;
    const double t94 = t21 * t93;
    const double t98 = -t20 * t94 / 0.12e2 - t91 * t94 / 0.12e2;
    const double t103 = t6 * t93 * t29;
    const double t111 = piecewise_functor_3( t13, -t88 * t98, -t28 * t71 / 0.3e1 - t34 * t103 / 0.12e2 - t33 * t94 / 0.12e2 - t39 * t94 / 0.12e2 );
    const double t112 = t52 * t52;
    const double t114 = t44 / t112;
    const double t116 = t45 * t89 * t1;
    const double t120 = -t116 * t94 / 0.12e2 - t49 * t94 / 0.12e2;
    const double t131 = piecewise_functor_3( t13, -t114 * t120, -t55 * t71 / 0.3e1 - t60 * t103 / 0.12e2 - t59 * t94 / 0.12e2 - t64 * t94 / 0.12e2 );
    const double t132 = t131 - t111;
    const double t134 = t132 * t79 * t84;
    const double t135 = t7 * t7;
    const double t136 = 0.1e1 / t135;
    const double t137 = t70 * t136;
    const double t138 = t71 - t137;
    const double t140 = -t138;
    const double t143 = 0.4e1 / 0.3e1 * t74 * t138 + 0.4e1 / 0.3e1 * t77 * t140;
    const double t145 = t69 * t143 * t84;
    const double t148 = -t71 - t137;
    const double t150 = -t148;
    const double t153 = 0.4e1 / 0.3e1 * t74 * t148 + 0.4e1 / 0.3e1 * t77 * t150;
    const double t155 = t69 * t153 * t84;


    eps = t43 + t85;
    vrho_a = t43 + t85 + t7 * ( t111 + t134 + t145 );
    vrho_b = t43 + t85 + t7 * ( t111 + t134 + t155 );

  }


};

struct BuiltinPZ81 : detail::BuiltinKernelImpl< BuiltinPZ81 > {

  BuiltinPZ81( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPZ81 >(p) { }
  
  virtual ~BuiltinPZ81() = default;

};



} // namespace ExchCXX