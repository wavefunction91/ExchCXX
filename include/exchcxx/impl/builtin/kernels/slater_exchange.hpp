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

  static constexpr double dens_tol  = 1e-24;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double alpha = 1.0;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;


    const double t11 = safe_math::cbrt( rho );
    const double t13 = t5 * t10 * t11;


    eps = -0.3e1 / 0.16e2 * t13;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;


    const double t11 = safe_math::cbrt( rho );
    const double t13 = t5 * t10 * t11;


    eps = -0.3e1 / 0.16e2 * t13;
    vrho = -t13 / 0.4e1;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t2 = alpha * t1;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;


    const double t8 = safe_math::cbrt( rho );
    const double t10 = t2 * t7 * t8;


    eps = -0.3e1 / 0.8e1 * t10;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double& eps, double& vrho ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t2 = alpha * t1;
    constexpr double t6 = t5 * t5;
    constexpr double t7 = t4 * t6;


    const double t8 = safe_math::cbrt( rho );
    const double t10 = t2 * t7 * t8;


    eps = -0.3e1 / 0.8e1 * t10;
    vrho = -t10 / 0.2e1;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double& eps ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;


    const double t11 = rho_a - rho_b;
    const double t12 = rho_a + rho_b;
    const double t13 = 0.1e1 / t12;
    const double t14 = t11 * t13;
    const double t15 = 0.1e1 + t14;
    const double t16 = safe_math::cbrt( t15 );
    const double t18 = 0.1e1 - t14;
    const double t19 = safe_math::cbrt( t18 );
    const double t21 = t16 * t15 + t19 * t18;
    const double t22 = safe_math::cbrt( t12 );
    const double t25 = t5 * t10 * t21 * t22;


    eps = -0.3e1 / 0.32e2 * t25;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double& eps, double& vrho_a, double& vrho_b ) {

    (void)(eps);
    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;
    constexpr double t31 = t4 * t7;


    const double t11 = rho_a - rho_b;
    const double t12 = rho_a + rho_b;
    const double t13 = 0.1e1 / t12;
    const double t14 = t11 * t13;
    const double t15 = 0.1e1 + t14;
    const double t16 = safe_math::cbrt( t15 );
    const double t18 = 0.1e1 - t14;
    const double t19 = safe_math::cbrt( t18 );
    const double t21 = t16 * t15 + t19 * t18;
    const double t22 = safe_math::cbrt( t12 );
    const double t25 = t5 * t10 * t21 * t22;
    const double t27 = t25 / 0.8e1;
    const double t30 = t22 * t12 * alpha * t1;
    const double t32 = t12 * t12;
    const double t33 = 0.1e1 / t32;
    const double t34 = t11 * t33;
    const double t35 = t13 - t34;
    const double t37 = -t35;
    const double t40 = 0.4e1 / 0.3e1 * t16 * t35 + 0.4e1 / 0.3e1 * t19 * t37;
    const double t45 = -t13 - t34;
    const double t47 = -t45;
    const double t52 = t31 * t9 * ( 0.4e1 / 0.3e1 * t16 * t45 + 0.4e1 / 0.3e1 * t19 * t47 );


    eps = -0.3e1 / 0.32e2 * t25;
    vrho_a = -t27 - 0.3e1 / 0.32e2 * t30 * t31 * t9 * t40;
    vrho_b = -t27 - 0.3e1 / 0.32e2 * t30 * t52;

  }


};

struct BuiltinSlaterExchange : detail::BuiltinKernelImpl< BuiltinSlaterExchange > {

  BuiltinSlaterExchange( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSlaterExchange >(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

};



} // namespace ExchCXX