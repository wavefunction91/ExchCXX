#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinSlaterExchange > :
  public lda_screening_interface< BuiltinSlaterExchange > {

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

  static constexpr double alpha = 1.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t15 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t10 * t11 );
    const double t16 = alpha * t15;


    eps = 0.2e1 * t16;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t6 = t3 / t4;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t8 = safe_math::cbrt( zeta_tol );
    const double t10 = piecewise_functor_3( 0.1e1 <= zeta_tol, t8 * zeta_tol, 1.0 );
    const double t11 = safe_math::cbrt( rho );
    const double t15 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t10 * t11 );
    const double t16 = alpha * t15;
    const double t17 = rho * alpha;
    const double t18 = t11 * t11;
    const double t23 = piecewise_functor_3( t2, 0.0, -t6 * t10 / t18 / 0.8e1 );


    eps = 0.2e1 * t16;
    vrho = 0.2e1 * t17 * t23 + 0.2e1 * t16;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t24 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t19 * t20 );
    const double t25 = alpha * t24;
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t39 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t35 * t20 );
    const double t40 = alpha * t39;


    eps = t25 + t40;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = 1.0/constants::m_cbrt_one_ov_pi;
    constexpr double t13 = constants::m_cbrt_2;
    constexpr double t5 = t2 / t3;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t8 = rho_a * t7;
    const double t10 = 0.2e1 * t8 <= zeta_tol;
    const double t11 = safe_math::cbrt( zeta_tol );
    const double t12 = t11 * zeta_tol;
    const double t14 = t13 * rho_a;
    const double t15 = safe_math::cbrt( t8 );
    const double t19 = piecewise_functor_3( t10, t12, 0.2e1 * t14 * t7 * t15 );
    const double t20 = safe_math::cbrt( t6 );
    const double t24 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t19 * t20 );
    const double t25 = alpha * t24;
    const double t26 = rho_b <= dens_tol;
    const double t27 = rho_b * t7;
    const double t29 = 0.2e1 * t27 <= zeta_tol;
    const double t30 = t13 * rho_b;
    const double t31 = safe_math::cbrt( t27 );
    const double t35 = piecewise_functor_3( t29, t12, 0.2e1 * t30 * t7 * t31 );
    const double t39 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t35 * t20 );
    const double t40 = alpha * t39;
    const double t41 = t13 * t7;
    const double t44 = t6 * t6;
    const double t45 = 0.1e1 / t44;
    const double t48 = 0.2e1 * t14 * t45 * t15;
    const double t49 = t15 * t15;
    const double t50 = 0.1e1 / t49;
    const double t51 = t7 * t50;
    const double t53 = -rho_a * t45 + t7;
    const double t58 = piecewise_functor_3( t10, 0.0, 0.2e1 * t41 * t15 - t48 + 0.2e1 / 0.3e1 * t14 * t51 * t53 );
    const double t62 = t20 * t20;
    const double t63 = 0.1e1 / t62;
    const double t66 = t5 * t19 * t63 / 0.8e1;
    const double t68 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t58 * t20 - t66 );
    const double t69 = alpha * t68;
    const double t72 = 0.2e1 * t30 * t45 * t31;
    const double t73 = rho_b * rho_b;
    const double t74 = t13 * t73;
    const double t75 = t44 * t6;
    const double t76 = 0.1e1 / t75;
    const double t77 = t31 * t31;
    const double t78 = 0.1e1 / t77;
    const double t79 = t76 * t78;
    const double t83 = piecewise_functor_3( t29, 0.0, -t72 - 0.2e1 / 0.3e1 * t74 * t79 );
    const double t89 = t5 * t35 * t63 / 0.8e1;
    const double t91 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t83 * t20 - t89 );
    const double t92 = alpha * t91;
    const double t95 = rho_a * rho_a;
    const double t96 = t13 * t95;
    const double t97 = t76 * t50;
    const double t101 = piecewise_functor_3( t10, 0.0, -t48 - 0.2e1 / 0.3e1 * t96 * t97 );
    const double t106 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t101 * t20 - t66 );
    const double t107 = alpha * t106;
    const double t110 = t7 * t78;
    const double t112 = -rho_b * t45 + t7;
    const double t117 = piecewise_functor_3( t29, 0.0, 0.2e1 * t41 * t31 - t72 + 0.2e1 / 0.3e1 * t30 * t110 * t112 );
    const double t122 = piecewise_functor_3( t26, 0.0, -0.3e1 / 0.8e1 * t5 * t117 * t20 - t89 );
    const double t123 = alpha * t122;


    eps = t25 + t40;
    vrho_a = t25 + t40 + t6 * ( t69 + t92 );
    vrho_b = t25 + t40 + t6 * ( t107 + t123 );

  }


};

struct BuiltinSlaterExchange : detail::BuiltinKernelImpl< BuiltinSlaterExchange > {

  BuiltinSlaterExchange( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSlaterExchange >(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

};



} // namespace ExchCXX