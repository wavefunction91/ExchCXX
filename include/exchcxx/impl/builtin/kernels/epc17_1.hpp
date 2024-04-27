#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinEPC17_1 > :
  public lda_screening_interface< BuiltinEPC17_1 > {

  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;

  static constexpr double dens_tol  = 1e-24;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.000000000000004e-32;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;

  static constexpr double a = 2.35; 
  static constexpr double b = 2.40; 
  static constexpr double c = 3.20; 

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;


    eps = piecewise_functor_3( t2, 0.0, -t9 * t19 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {



    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t3 = 0.1e1 <= zeta_tol;
    const double t4 = zeta_tol - 0.1e1;
    const double t6 = piecewise_functor_5( t3, t4, t3, -t4, 0.0 );
    const double t8 = square( 0.1e1 + t6 );
    const double t9 = t8 * rho;
    const double t10 = rho * rho;
    const double t11 = t8 * t10;
    const double t12 = safe_math::sqrt( t11 );
    const double t15 = c * t8;
    const double t18 = a - b * t12 / 0.2e1 + t15 * t10 / 0.4e1;
    const double t19 = 0.1e1 / t18;
    const double t23 = t18 * t18;
    const double t24 = 0.1e1 / t23;
    const double t26 = b / t12;
    const double t30 = t15 * rho / 0.2e1 - t26 * t9 / 0.2e1;
    const double t35 = piecewise_functor_3( t2, 0.0, t9 * t24 * t30 / 0.4e1 - t8 * t19 / 0.4e1 );


    eps = piecewise_functor_3( t2, 0.0, -t9 * t19 / 0.4e1 );
    vrho = rho * t35 + eps;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);


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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;


    eps = piecewise_functor_3( t3, 0.0, -t18 * t35 / 0.4e1 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {



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
    const double t23 = t4 * t4;
    const double t24 = t17 * t23;
    const double t25 = t24 * t22;
    const double t26 = safe_math::sqrt( t25 );
    const double t29 = c * t17;
    const double t30 = t23 * t22;
    const double t33 = a - b * t26 / 0.2e1 + t29 * t30 / 0.4e1;
    const double t34 = 0.1e1 / t33;
    const double t35 = t22 * t34;
    const double t38 = 0.1e1 / t23;
    const double t39 = t14 * t38;
    const double t41 = piecewise_functor_5( t8, 0.0, t12, 0.0, t5 - t39 );
    const double t42 = t41 * t4;
    const double t44 = t17 * t22;
    const double t45 = t44 * t34;
    const double t46 = t19 * t38;
    const double t48 = piecewise_functor_5( t12, 0.0, t8, 0.0, -t5 - t46 );
    const double t49 = t48 * t34;
    const double t51 = t33 * t33;
    const double t52 = 0.1e1 / t51;
    const double t53 = t22 * t52;
    const double t55 = b / t26;
    const double t56 = t41 * t23;
    const double t58 = t18 * t22;
    const double t59 = 0.2e1 * t58;
    const double t61 = t56 * t22 + t24 * t48 + t59;
    const double t64 = c * t41;
    const double t67 = t4 * t22;
    const double t69 = t29 * t67 / 0.2e1;
    const double t70 = t23 * t48;
    const double t73 = -t55 * t61 / 0.4e1 + t64 * t30 / 0.4e1 + t69 + t29 * t70 / 0.4e1;
    const double t74 = t53 * t73;
    const double t78 = piecewise_functor_3( t3, 0.0, -t18 * t49 / 0.4e1 + t18 * t74 / 0.4e1 - t42 * t35 / 0.4e1 - t45 / 0.4e1 );
    const double t81 = piecewise_functor_5( t8, 0.0, t12, 0.0, -t5 - t39 );
    const double t82 = t81 * t4;
    const double t85 = piecewise_functor_5( t12, 0.0, t8, 0.0, t5 - t46 );
    const double t86 = t85 * t34;
    const double t88 = t81 * t23;
    const double t91 = t88 * t22 + t24 * t85 + t59;
    const double t94 = c * t81;
    const double t97 = t23 * t85;
    const double t100 = -t55 * t91 / 0.4e1 + t94 * t30 / 0.4e1 + t69 + t29 * t97 / 0.4e1;
    const double t101 = t53 * t100;
    const double t105 = piecewise_functor_3( t3, 0.0, t18 * t101 / 0.4e1 - t18 * t86 / 0.4e1 - t82 * t35 / 0.4e1 - t45 / 0.4e1 );


    eps = piecewise_functor_3( t3, 0.0, -t18 * t35 / 0.4e1 );
    vrho_a = t4 * t78 + eps;
    vrho_b = t4 * t105 + eps;

  }


};

struct BuiltinEPC17_1 : detail::BuiltinKernelImpl< BuiltinEPC17_1 > {

  BuiltinEPC17_1( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinEPC17_1 >(p) { }
  
  virtual ~BuiltinEPC17_1() = default;

};



} // namespace ExchCXX