#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {

template <>
struct kernel_traits<BuiltinB88> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-25;

  static constexpr double b88_beta  = 0.0042;
  static constexpr double b88_gamma = 6.0;


  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

#ifdef __CUDACC__
      rho   = fmax( rho, 0. );
      sigma = fmax( sigma, 1e-40 );
#else
      rho   = std::max( rho, 0. );
      sigma = std::max( sigma, 1e-40 );
#endif
      
      if( rho < dens_tol ) {
        eps = 0.;
        return;
      }

      constexpr double t1 = constants::m_cbrt_3;
      constexpr double t5 = constants::m_cbrt_4;
      constexpr double t8 = constants::m_cbrt_2;
      constexpr double t3 = constants::m_cbrt_one_ov_pi;
      constexpr double t6 = t5 * t5;
      constexpr double t7 = t1 * t3 * t6;
      constexpr double t9 = t8 * t8;
      constexpr double t12 = t1 * t1;
      constexpr double t13 = b88_beta * t12;
      constexpr double t14 = 1. / t3;
      constexpr double t15 = t14 * t5;
      constexpr double t16 = t13 * t15;
      constexpr double t22 = b88_gamma * b88_beta;

      const double t10 = std::pow( rho, constants::m_third );
      const double t11 = t9 * t10;
      const double t17 = sigma * t9;
      const double t18 = rho * rho;
      const double t19 = t10 * t10;
      const double t21 = 1. / t19 / t18;
      const double t23 = std::sqrt(sigma);
      const double t24 = t22 * t23;
      const double t25 = t10 * rho;
      const double t26 = 1. / t25;
      const double x1  = t23 * t8 * t26;
      const double t30 = log(t23 * t8 * t26 + sqrt(x1 * x1 + 1.));
      const double t31 = t8 * t26 * t30;
      const double t33 = 1. + t24 * t31;
      const double t34 = 1. / t33;
      const double t35 = t21 * t34;
      const double t39 = 1. + 2. / 9. * t16 * t17 * t35;
      const double t41 = t7 * t11 * t39;

      eps = -3. / 16. * t41;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

#ifdef __CUDACC__
      rho   = fmax( rho, 0. );
      sigma = fmax( sigma, 1e-40 );
#else
      rho   = std::max( rho, 0. );
      sigma = std::max( sigma, 1e-40 );
#endif
      
      if( rho < dens_tol ) {
        eps = 0.;
        vrho = 0.;
        vsigma = 0.;
        return;
      }

      constexpr double t1 = constants::m_cbrt_3;
      constexpr double t5 = constants::m_cbrt_4;
      constexpr double t8 = constants::m_cbrt_2;
      constexpr double t3 = constants::m_cbrt_one_ov_pi;
      constexpr double t6 = t5 * t5;
      constexpr double t7 = t1 * t3 * t6;
      constexpr double t9 = t8 * t8;
      constexpr double t12 = t1 * t1;
      constexpr double t13 = b88_beta * t12;
      constexpr double t14 = 1. / t3;
      constexpr double t15 = t14 * t5;
      constexpr double t16 = t13 * t15;
      constexpr double t22 = b88_gamma * b88_beta;

      constexpr double t46 = t6 * t9;
      constexpr double t81 = t5 * t9;

      const double t10 = std::pow( rho, constants::m_third );
      const double t11 = t9 * t10;
      const double t17 = sigma * t9;
      const double t18 = rho * rho;
      const double t19 = t10 * t10;
      const double t21 = 1. / t19 / t18;
      const double t23 = std::sqrt(sigma);
      const double t24 = t22 * t23;
      const double t25 = t10 * rho;
      const double t26 = 1. / t25;
      const double x1  = t23 * t8 * t26;
      const double t30 = log(t23 * t8 * t26 + sqrt(x1 * x1 + 1.));
      const double t31 = t8 * t26 * t30;
      const double t33 = 1. + t24 * t31;
      const double t34 = 1. / t33;
      const double t35 = t21 * t34;
      const double t39 = 1. + 2. / 9. * t16 * t17 * t35;
      const double t41 = t7 * t11 * t39;

      eps = -3. / 16. * t41;

      const double t45 = t25 * t1 * t3;
      const double t47 = t18 * rho;
      const double t49 = 1. / t19 / t47;
      const double t50 = t49 * t34;
      const double t54 = t33 * t33;
      const double t55 = 1. / t54;
      const double t56 = t21 * t55;
      const double t60 = t8 / t10 / t18 * t30;
      const double t62 = t22 * sigma;
      const double t63 = t9 * t49;
      const double t65 = t17 * t21 + 1.;
      const double t66 = std::sqrt(t65);
      const double t67 = 1. / t66;
      const double t68 = t63 * t67;
      const double t71 = -4. / 3. * t24 * t60 - 4. / 3. * t62 * t68;
      const double t76 = -16. / 27. * t16 * t17 * t50 - 
                         2. / 9. * t16 * t17 * t56 * t71;


      const double t80 = t13 * t14;
      const double t85 = t22 / t23;
      const double t87 = t9 * t21;
      const double t88 = t87 * t67;
      const double t91 = t22 * t88 / 2. + t85 * t31 / 2.;
      const double t97 = t46 * (-2. / 9. * t16 * t17 * t56 * t91 + 
                    2. / 9. * t80 * t81 * t35);

      vrho   = -t41 / 4. - 3. / 16. * t45 * t46 * t76;
      vsigma = -3. / 16. * t45 * t97;

  }
};


struct BuiltinB88 : detail::BuiltinKernelImpl< BuiltinB88 > {

  BuiltinB88( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinB88 >(p) { }
  
  virtual ~BuiltinB88() = default;

};


}
