#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>


namespace ExchCXX {

template <>
struct kernel_traits<BuiltinVWN_RPA> :
  public lda_screening_interface< BuiltinVWN_RPA > {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-24;

  BUILTIN_KERNEL_EVAL_RETURN 
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;

    const double t7 = std::pow( rho, constants::m_third );
    const double t8 = 1. / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = std::sqrt(t10);
    const double t14 = t10 / 4. + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 1. / t14;
    const double t19 = log(t4 * t9 * t15 / 4.);
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.130720e2;
    const double t24 = atan(0.44899888641287296627e-1 / t21);
    const double t25 = 0.20521972937837502661e2 * t24;
    const double t27 = t12 / 2. + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = log(t28 * t15);
    const double t31 = 0.44313737677495382697e-2 * t30;

    eps = t20 + t25 + t31;

  }



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vxc ) {

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_4;
    constexpr double t3 = constants::m_cbrt_one_ov_pi;
    constexpr double t4 = t1 * t3;
    constexpr double t6 = t5 * t5;

    const double t7 = std::pow( rho, constants::m_third );
    const double t8 = 1. / t7;
    const double t9 = t6 * t8;
    const double t10 = t4 * t9;
    const double t12 = std::sqrt(t10);
    const double t14 = t10 / 4. + 0.6536e1 * t12 + 0.427198e2;
    const double t15 = 1. / t14;
    const double t19 = log(t4 * t9 * t15 / 4.);
    const double t20 = 0.310907e-1 * t19;
    const double t21 = t12 + 0.130720e2;
    const double t24 = atan(0.44899888641287296627e-1 / t21);
    const double t25 = 0.20521972937837502661e2 * t24;
    const double t27 = t12 / 2. + 0.409286e0;
    const double t28 = t27 * t27;
    const double t30 = log(t28 * t15);
    const double t31 = 0.44313737677495382697e-2 * t30;

    eps = t20 + t25 + t31;

    const double t33 = 1. / t7 / rho;
    const double t34 = t6 * t33;
    const double t38 = t4 * t6;
    const double t39 = t14 * t14;
    const double t40 = 1. / t39;
    const double t41 = t8 * t40;
    const double t42 = t4 * t34;
    const double t44 = 1. / t12;
    const double t45 = t44 * t1;
    const double t46 = t3 * t6;
    const double t50 = -t42 / 12. - 0.10893333333333333333e1 * t45 * t46 * t33;
    const double t55 = t1 * t1;
    const double t57 = 1. / t3;
    const double t58 = (-t4 * t34 * t15 / 12. - t38 * t41 * t50 / 4.) * 
                       t55 * t57;
    const double t59 = t5 * t7;
    const double t60 = t59 * t14;
    const double t61 = t58 * t60;
    const double t63 = t21 * t21;
    const double t64 = 1. / t63;
    const double t66 = t64 * t44 * t1;
    const double t68 = 0.2016e-2 * t64 + 1.;
    const double t69 = 1. / t68;
    const double t72 = t66 * t46 * t33 * t69;
    const double t74 = t27 * t15;
    const double t75 = t74 * t44;
    const double t78 = t28 * t40;
    const double t80 = -t75 * t42 / 6. - t78 * t50;
    const double t81 = 1. / t28;
    const double t82 = t80 * t81;
    const double t83 = t82 * t14;
      
    vxc = t20 + t25 + t31 + 
      rho * (0.10363566666666666667e-1 * t61 + 0.15357238326806922974e0 * t72 + 
        0.44313737677495382697e-2 * t83);
      
  }

};



struct BuiltinVWN_RPA : detail::BuiltinKernelImpl< BuiltinVWN_RPA > {

  BuiltinVWN_RPA( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinVWN_RPA >(p) { }
  
  virtual ~BuiltinVWN_RPA() = default;

};

}
