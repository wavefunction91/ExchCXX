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
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-32;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 2.1544346900318956e-43;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double kappa = 1.245;
  static constexpr double mu =  0.2195149727645171;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );


    eps = 0.2e1 * t46;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t28 = t27 * t27;
    constexpr double t56 = kappa * kappa;
    constexpr double t78 = t20 * t25 * t28;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t37 = kappa + mu * t20 * t25 * sigma * t28 * t33 / 0.24e2;
    const double t42 = 0.1e1 + kappa * ( 0.1e1 - kappa / t37 );
    const double t46 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t17 * t18 * t42 );
    const double t52 = t30 * rho;
    const double t58 = t6 * t17 / t18 / t52 * t56;
    const double t59 = t37 * t37;
    const double t61 = 0.1e1 / t59 * mu;
    const double t64 = t25 * sigma * t28;
    const double t65 = t61 * t20 * t64;
    const double t69 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t31 * t42 / 0.8e1 + t58 * t65 / 0.24e2 );
    const double t79 = t61 * t78;
    const double t82 = piecewise_functor_3( t2, 0.0, -t6 * t17 / t18 / t30 * t56 * t79 / 0.64e2 );


    eps = 0.2e1 * t46;
    vrho = 0.2e1 * rho * t69 + 0.2e1 * t46;
    vsigma = 0.2e1 * rho * t82;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t81 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t77 );


    eps = t52 + t81;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = mu * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t100 = kappa * kappa;
    constexpr double t171 = t28 * t33;


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
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t43 = kappa + t29 * t34 * t39 / 0.24e2;
    const double t48 = 0.1e1 + kappa * ( 0.1e1 - kappa / t43 );
    const double t52 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t48 );
    const double t53 = rho_b <= dens_tol;
    const double t54 = -t16;
    const double t56 = piecewise_functor_5( t14, t11, t10, t15, t54 * t7 );
    const double t57 = 0.1e1 + t56;
    const double t58 = t57 <= zeta_tol;
    const double t59 = safe_math::cbrt( t57 );
    const double t61 = piecewise_functor_3( t58, t22, t59 * t57 );
    const double t62 = t61 * t26;
    const double t63 = t33 * sigma_bb;
    const double t64 = rho_b * rho_b;
    const double t65 = safe_math::cbrt( rho_b );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t72 = kappa + t29 * t63 * t68 / 0.24e2;
    const double t77 = 0.1e1 + kappa * ( 0.1e1 - kappa / t72 );
    const double t81 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t62 * t77 );
    const double t82 = t6 * t6;
    const double t83 = 0.1e1 / t82;
    const double t84 = t16 * t83;
    const double t86 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t84 );
    const double t89 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t86 );
    const double t90 = t89 * t26;
    const double t94 = t26 * t26;
    const double t95 = 0.1e1 / t94;
    const double t96 = t25 * t95;
    const double t99 = t5 * t96 * t48 / 0.8e1;
    const double t101 = t27 * t100;
    const double t102 = t5 * t101;
    const double t103 = t43 * t43;
    const double t105 = 0.1e1 / t103 * mu;
    const double t106 = t105 * t28;
    const double t107 = t35 * rho_a;
    const double t109 = 0.1e1 / t37 / t107;
    const double t111 = t106 * t34 * t109;
    const double t115 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t90 * t48 - t99 + t102 * t111 / 0.24e2 );
    const double t116 = t54 * t83;
    const double t118 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t116 );
    const double t121 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t118 );
    const double t122 = t121 * t26;
    const double t126 = t61 * t95;
    const double t129 = t5 * t126 * t77 / 0.8e1;
    const double t131 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t122 * t77 - t129 );
    const double t135 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t84 );
    const double t138 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t135 );
    const double t139 = t138 * t26;
    const double t144 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t139 * t48 - t99 );
    const double t146 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t116 );
    const double t149 = piecewise_functor_3( t58, 0.0, 0.4e1 / 0.3e1 * t59 * t146 );
    const double t150 = t149 * t26;
    const double t154 = t62 * t100;
    const double t155 = t5 * t154;
    const double t156 = t72 * t72;
    const double t158 = 0.1e1 / t156 * mu;
    const double t159 = t158 * t28;
    const double t160 = t64 * rho_b;
    const double t162 = 0.1e1 / t66 / t160;
    const double t164 = t159 * t63 * t162;
    const double t168 = piecewise_functor_3( t53, 0.0, -0.3e1 / 0.8e1 * t5 * t150 * t77 - t129 + t155 * t164 / 0.24e2 );
    const double t173 = t105 * t171 * t39;
    const double t176 = piecewise_functor_3( t1, 0.0, -t102 * t173 / 0.64e2 );
    const double t178 = t158 * t171 * t68;
    const double t181 = piecewise_functor_3( t53, 0.0, -t155 * t178 / 0.64e2 );


    eps = t52 + t81;
    vrho_a = t52 + t81 + t6 * ( t115 + t131 );
    vrho_b = t52 + t81 + t6 * ( t144 + t168 );
    vsigma_aa = t6 * t176;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t181;

  }


};

struct BuiltinRevPBE_X : detail::BuiltinKernelImpl< BuiltinRevPBE_X > {

  BuiltinRevPBE_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinRevPBE_X >(p) { }
  
  virtual ~BuiltinRevPBE_X() = default;

};



} // namespace ExchCXX
