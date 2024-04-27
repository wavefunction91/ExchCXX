#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinEPC18_2 > :
  public lda_screening_interface< BuiltinEPC18_2 > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;

  static constexpr double a = 3.90; 
  static constexpr double b = 0.50; 
  static constexpr double c = 0.06; 

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;


    eps = piecewise_functor_3( t2, 0.0, -t9 * t18 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t7 = 0.1e1 + t6;
    const double t8 = t7 * t7;
    const double t9 = t8 * rho;
    const double t10 = b * t7;
    const double t13 = c * t8;
    const double t14 = rho * rho;
    const double t17 = -0.4e1 * t10 * rho + 0.16e2 * t13 * t14 + a;
    const double t18 = 0.1e1 / t17;
    const double t22 = t17 * t17;
    const double t23 = 0.1e1 / t22;
    const double t27 = 0.32e2 * t13 * rho - 0.4e1 * t10;
    const double t32 = piecewise_functor_3( t2, 0.0, t9 * t23 * t27 / 0.4e1 - t8 * t18 / 0.4e1 );


    eps = piecewise_functor_3( t2, 0.0, -t9 * t18 / 0.4e1 );
    vrho = rho * t32 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;


    eps = piecewise_functor_3( t3, 0.0, -t18 * t40 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    constexpr double t23 = constants::m_cbrt_2;
    constexpr double t24 = t23 * t23;


    const double t3 = rho_a <= dens_tol && rho_b <= dens_tol;
    const double t4 = rho_a + rho_b;
    const double t5 = 0.1e1 / t4;
    const double t8 = 0.2e1 * rho_a * t5 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t12 = 0.2e1 * rho_b * t5 <= zeta_tol;
    const double t13 = -t9;
    const double t14 = rho_a - rho_b;
    const double t16 = piecewise_functor_5( t8, t9, t12, t13, t14 * t5 );
    const double t17 = 0.1e1 + t16;
    const double t18 = t17 * t4;
    const double t19 = -t14;
    const double t21 = piecewise_functor_5( t12, t9, t8, t13, t19 * t5 );
    const double t22 = 0.1e1 + t21;
    const double t25 = safe_math::cbrt( t18 );
    const double t27 = t22 * t4;
    const double t28 = safe_math::cbrt( t27 );
    const double t31 = t24 * t25 / 0.2e1 + t24 * t28 / 0.2e1;
    const double t32 = t31 * t31;
    const double t33 = t32 * t31;
    const double t35 = t32 * t32;
    const double t38 = c * t35 * t32 - b * t33 + a;
    const double t39 = 0.1e1 / t38;
    const double t40 = t22 * t39;
    const double t43 = t4 * t4;
    const double t44 = 0.1e1 / t43;
    const double t45 = t14 * t44;
    const double t47 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t45 );
    const double t48 = t47 * t4;
    const double t50 = t17 * t22;
    const double t51 = t50 * t39;
    const double t52 = t19 * t44;
    const double t54 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t52 );
    const double t55 = t54 * t39;
    const double t57 = t38 * t38;
    const double t58 = 0.1e1 / t57;
    const double t59 = t22 * t58;
    const double t60 = b * t32;
    const double t61 = t25 * t25;
    const double t63 = t24 / t61;
    const double t64 = t48 + 0.1e1 + t16;
    const double t66 = t28 * t28;
    const double t68 = t24 / t66;
    const double t70 = t54 * t4 + t21 + 0.1e1;
    const double t73 = t63 * t64 / 0.6e1 + t68 * t70 / 0.6e1;
    const double t77 = c * t35 * t31;
    const double t80 = -0.3e1 * t60 * t73 + 0.6e1 * t77 * t73;
    const double t81 = t59 * t80;
    const double t85 = piecewise_functor_3( t3, 0.0, -t18 * t55 / 0.4e1 + t18 * t81 / 0.4e1 - t48 * t40 / 0.4e1 - t51 / 0.4e1 );
    const double t88 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t45 );
    const double t89 = t88 * t4;
    const double t92 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t52 );
    const double t93 = t92 * t39;
    const double t95 = t89 + 0.1e1 + t16;
    const double t98 = t92 * t4 + t21 + 0.1e1;
    const double t101 = t63 * t95 / 0.6e1 + t68 * t98 / 0.6e1;
    const double t106 = -0.3e1 * t60 * t101 + 0.6e1 * t77 * t101;
    const double t107 = t59 * t106;
    const double t111 = piecewise_functor_3( t3, 0.0, t18 * t107 / 0.4e1 - t18 * t93 / 0.4e1 - t89 * t40 / 0.4e1 - t51 / 0.4e1 );


    eps = piecewise_functor_3( t3, 0.0, -t18 * t40 / 0.4e1 );
    vrho_a = t4 * t85 + eps;
    vrho_b = t4 * t111 + eps;

  }


};

struct BuiltinEPC18_2 : detail::BuiltinKernelImpl< BuiltinEPC18_2 > {

  BuiltinEPC18_2( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinEPC18_2 >(p) { }
  
  virtual ~BuiltinEPC18_2() = default;

};



} // namespace ExchCXX