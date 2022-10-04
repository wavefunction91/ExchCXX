#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinRevPBE_X > :
  public gga_screening_interface< BuiltinRevPBE_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;

  static constexpr double dens_tol  = 1e-32;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double kappa = 1.245;
  static constexpr double mu =  0.2195149727645171;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t12 = constants::m_cbrt_6;
    constexpr double t15 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t16 = t15 * t15;
    constexpr double t17 = 0.1e1 / t16;
    constexpr double t18 = mu * t12 * t17;


    const double t10 = safe_math::cbrt( rho );
    const double t20 = rho * rho;
    const double t21 = t10 * t10;
    const double t23 = 0.1e1 / t21 / t20;
    const double t27 = kappa + t18 * sigma * t9 * t23 / 0.24e2;
    const double t32 = 0.1e1 + kappa * ( 0.1e1 - kappa / t27 );
    const double t34 = t7 * t9 * t10 * t32;


    eps = -0.3e1 / 0.16e2 * t34;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t12 = constants::m_cbrt_6;
    constexpr double t15 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t16 = t15 * t15;
    constexpr double t17 = 0.1e1 / t16;
    constexpr double t18 = mu * t12 * t17;
    constexpr double t40 = t3 * t6;
    constexpr double t41 = t40 * t8;
    constexpr double t43 = kappa * kappa;


    const double t10 = safe_math::cbrt( rho );
    const double t20 = rho * rho;
    const double t21 = t10 * t10;
    const double t23 = 0.1e1 / t21 / t20;
    const double t27 = kappa + t18 * sigma * t9 * t23 / 0.24e2;
    const double t32 = 0.1e1 + kappa * ( 0.1e1 - kappa / t27 );
    const double t34 = t7 * t9 * t10 * t32;
    const double t42 = 0.1e1 / t10 / t20 * t1 * t41;
    const double t44 = t27 * t27;
    const double t46 = t43 / t44;
    const double t49 = t12 * t17 * sigma;
    const double t50 = t46 * mu * t49;
    const double t57 = t46 * t18;


    eps = -0.3e1 / 0.16e2 * t34;
    vrho = -t34 / 0.4e1 + t42 * t50 / 0.24e2;
    vsigma = -0.1e1 / t10 / rho * t1 * t41 * t57 / 0.64e2;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t9 = constants::m_cbrt_6;
    constexpr double t12 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t10 = mu * t9;
    constexpr double t13 = t12 * t12;
    constexpr double t14 = 0.1e1 / t13;


    const double t7 = safe_math::cbrt( rho );
    const double t16 = rho * rho;
    const double t17 = t7 * t7;
    const double t19 = 0.1e1 / t17 / t16;
    const double t23 = kappa + t10 * t14 * sigma * t19 / 0.24e2;
    const double t28 = 0.1e1 + kappa * ( 0.1e1 - kappa / t23 );
    const double t30 = t4 * t6 * t7 * t28;


    eps = -0.3e1 / 0.8e1 * t30;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t9 = constants::m_cbrt_6;
    constexpr double t12 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t10 = mu * t9;
    constexpr double t13 = t12 * t12;
    constexpr double t14 = 0.1e1 / t13;
    constexpr double t36 = t3 * t6;
    constexpr double t37 = kappa * kappa;
    constexpr double t38 = t36 * t37;


    const double t7 = safe_math::cbrt( rho );
    const double t16 = rho * rho;
    const double t17 = t7 * t7;
    const double t19 = 0.1e1 / t17 / t16;
    const double t23 = kappa + t10 * t14 * sigma * t19 / 0.24e2;
    const double t28 = 0.1e1 + kappa * ( 0.1e1 - kappa / t23 );
    const double t30 = t4 * t6 * t7 * t28;
    const double t35 = 0.1e1 / t7 / t16 * t1;
    const double t40 = t23 * t23;
    const double t41 = 0.1e1 / t40;
    const double t44 = t9 * t14 * sigma;
    const double t45 = t41 * mu * t44;
    const double t54 = t37 * t41 * t10 * t14;


    eps = -0.3e1 / 0.8e1 * t30;
    vrho = -t30 / 0.2e1 + t35 * t38 * t45 / 0.24e2;
    vsigma = -0.1e1 / t7 / rho * t1 * t36 * t54 / 0.64e2;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t18 = constants::m_cbrt_6;
    constexpr double t21 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t19 = mu * t18;
    constexpr double t22 = t21 * t21;
    constexpr double t23 = 0.1e1 / t22;


    const double t8 = rho_a - rho_b;
    const double t9 = rho_a + rho_b;
    const double t10 = 0.1e1 / t9;
    const double t11 = t8 * t10;
    const double t13 = 0.1e1 / 0.2e1 + t11 / 0.2e1;
    const double t14 = safe_math::cbrt( t13 );
    const double t15 = t14 * t13;
    const double t16 = safe_math::cbrt( t9 );
    const double t17 = t15 * t16;
    const double t24 = t23 * sigma_aa;
    const double t25 = rho_a * rho_a;
    const double t26 = safe_math::cbrt( rho_a );
    const double t27 = t26 * t26;
    const double t29 = 0.1e1 / t27 / t25;
    const double t33 = kappa + t19 * t24 * t29 / 0.24e2;
    const double t38 = 0.1e1 + kappa * ( 0.1e1 - kappa / t33 );
    const double t40 = t7 * t17 * t38;
    const double t42 = 0.1e1 / 0.2e1 - t11 / 0.2e1;
    const double t43 = safe_math::cbrt( t42 );
    const double t44 = t43 * t42;
    const double t45 = t44 * t16;
    const double t46 = t23 * sigma_bb;
    const double t47 = rho_b * rho_b;
    const double t48 = safe_math::cbrt( rho_b );
    const double t49 = t48 * t48;
    const double t51 = 0.1e1 / t49 / t47;
    const double t55 = kappa + t19 * t46 * t51 / 0.24e2;
    const double t60 = 0.1e1 + kappa * ( 0.1e1 - kappa / t55 );
    const double t62 = t7 * t45 * t60;


    eps = -0.3e1 / 0.8e1 * t40 - 0.3e1 / 0.8e1 * t62;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t18 = constants::m_cbrt_6;
    constexpr double t21 = constants::m_cbrt_pi_sq;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;
    constexpr double t19 = mu * t18;
    constexpr double t22 = t21 * t21;
    constexpr double t23 = 0.1e1 / t22;
    constexpr double t82 = kappa * kappa;
    constexpr double t136 = t3 * t6;
    constexpr double t141 = t18 * t23;


    const double t8 = rho_a - rho_b;
    const double t9 = rho_a + rho_b;
    const double t10 = 0.1e1 / t9;
    const double t11 = t8 * t10;
    const double t13 = 0.1e1 / 0.2e1 + t11 / 0.2e1;
    const double t14 = safe_math::cbrt( t13 );
    const double t15 = t14 * t13;
    const double t16 = safe_math::cbrt( t9 );
    const double t17 = t15 * t16;
    const double t24 = t23 * sigma_aa;
    const double t25 = rho_a * rho_a;
    const double t26 = safe_math::cbrt( rho_a );
    const double t27 = t26 * t26;
    const double t29 = 0.1e1 / t27 / t25;
    const double t33 = kappa + t19 * t24 * t29 / 0.24e2;
    const double t38 = 0.1e1 + kappa * ( 0.1e1 - kappa / t33 );
    const double t40 = t7 * t17 * t38;
    const double t42 = 0.1e1 / 0.2e1 - t11 / 0.2e1;
    const double t43 = safe_math::cbrt( t42 );
    const double t44 = t43 * t42;
    const double t45 = t44 * t16;
    const double t46 = t23 * sigma_bb;
    const double t47 = rho_b * rho_b;
    const double t48 = safe_math::cbrt( rho_b );
    const double t49 = t48 * t48;
    const double t51 = 0.1e1 / t49 / t47;
    const double t55 = kappa + t19 * t46 * t51 / 0.24e2;
    const double t60 = 0.1e1 + kappa * ( 0.1e1 - kappa / t55 );
    const double t62 = t7 * t45 * t60;
    const double t64 = 0.3e1 / 0.8e1 * t40;
    const double t65 = 0.3e1 / 0.8e1 * t62;
    const double t66 = t14 * t16;
    const double t67 = t9 * t9;
    const double t68 = 0.1e1 / t67;
    const double t69 = t8 * t68;
    const double t71 = t10 / 0.2e1 - t69 / 0.2e1;
    const double t72 = t38 * t71;
    const double t74 = t7 * t66 * t72;
    const double t75 = t74 / 0.2e1;
    const double t76 = t16 * t16;
    const double t77 = 0.1e1 / t76;
    const double t78 = t15 * t77;
    const double t80 = t7 * t78 * t38;
    const double t81 = t80 / 0.8e1;
    const double t84 = t7 * t17 * t82;
    const double t85 = t33 * t33;
    const double t86 = 0.1e1 / t85;
    const double t88 = t86 * mu * t18;
    const double t89 = t25 * rho_a;
    const double t91 = 0.1e1 / t27 / t89;
    const double t93 = t88 * t24 * t91;
    const double t94 = t84 * t93;
    const double t95 = t94 / 0.24e2;
    const double t96 = t43 * t16;
    const double t97 = -t71;
    const double t98 = t60 * t97;
    const double t100 = t7 * t96 * t98;
    const double t101 = t100 / 0.2e1;
    const double t102 = t44 * t77;
    const double t104 = t7 * t102 * t60;
    const double t105 = t104 / 0.8e1;
    const double t109 = -t10 / 0.2e1 - t69 / 0.2e1;
    const double t110 = t38 * t109;
    const double t112 = t7 * t66 * t110;
    const double t113 = t112 / 0.2e1;
    const double t114 = -t109;
    const double t115 = t60 * t114;
    const double t117 = t7 * t96 * t115;
    const double t118 = t117 / 0.2e1;
    const double t120 = t7 * t45 * t82;
    const double t121 = t55 * t55;
    const double t122 = 0.1e1 / t121;
    const double t124 = t122 * mu * t18;
    const double t125 = t47 * rho_b;
    const double t127 = 0.1e1 / t49 / t125;
    const double t129 = t124 * t46 * t127;
    const double t130 = t120 * t129;
    const double t131 = t130 / 0.24e2;
    const double t135 = t16 * t9 * t1;
    const double t138 = t135 * t136 * t15;
    const double t139 = t82 * t86;
    const double t140 = t139 * mu;
    const double t143 = t140 * t141 * t29;
    const double t147 = t135 * t136 * t44;
    const double t148 = t82 * t122;
    const double t149 = t148 * mu;
    const double t151 = t149 * t141 * t51;


    eps = -0.3e1 / 0.8e1 * t40 - 0.3e1 / 0.8e1 * t62;
    vrho_a = -t64 - t65 + t9 * ( -t75 - t81 + t95 - t101 - t105 );
    vrho_b = -t64 - t65 + t9 * ( -t113 - t81 - t118 - t105 + t131 );
    vsigma_aa = -t138 * t143 / 0.64e2;
    vsigma_ab = 0.0e0;
    vsigma_bb = -t147 * t151 / 0.64e2;

  }


};

struct BuiltinRevPBE_X : detail::BuiltinKernelImpl< BuiltinRevPBE_X > {

  BuiltinRevPBE_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinRevPBE_X >(p) { }
  
  virtual ~BuiltinRevPBE_X() = default;

};



} // namespace ExchCXX