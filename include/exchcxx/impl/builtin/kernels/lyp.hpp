#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {

template <>
struct kernel_traits<BuiltinLYP> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-32;


  static constexpr double lyp_A = 0.04918;
  static constexpr double lyp_B = 0.132;
  static constexpr double lyp_c = 0.2533;
  static constexpr double lyp_d = 0.349;


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

      constexpr double t26 = constants::m_cbrt_3;
      constexpr double t27 = t26 * t26;
      constexpr double t29 = constants::m_cbrt_pi_sq;
      constexpr double t30 = t29 * t29;

      double t7 = std::pow( rho, constants::m_third );;
      double t8 = 1. / t7;
      double t10 = lyp_d * t8 + 1.;
      double t11 = 1. / t10;
      double t13 = exp(-lyp_c * t8);
      double t14 = lyp_B * t13;
      double t15 = rho * rho;
      double t16 = t7 * t7;
      double t18 = 1. / t16 / t15;
      double t19 = sigma * t18;
      double t21 = lyp_d * t11 + lyp_c;
      double t22 = t21 * t8;
      double t24 = -1. / 72. - 7. / 72. * t22;

      double t34 = 5. / 2. - t22 / 18.;
      double t35 = t34 * sigma;
      double t38 = t22 - 11.;
      double t39 = t38 * sigma;
      double t43 = -t19 * t24 - 3. / 10. * t27 * t30 + t35 * t18 / 8. + 
                   t39 * t18 / 144. - 5. / 24. * t19;

      eps = lyp_A * (t14 * t11 * t43 - t11);

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

      constexpr double t26 = constants::m_cbrt_3;
      constexpr double t27 = t26 * t26;
      constexpr double t29 = constants::m_cbrt_pi_sq;
      constexpr double t30 = t29 * t29;
      constexpr double t55 = lyp_B * lyp_c;
      constexpr double t72 = lyp_d * lyp_d;

      double t7 = std::pow( rho, constants::m_third );;
      double t8 = 1. / t7;
      double t10 = lyp_d * t8 + 1.;
      double t11 = 1. / t10;
      double t13 = exp(-lyp_c * t8);
      double t14 = lyp_B * t13;
      double t15 = rho * rho;
      double t16 = t7 * t7;
      double t18 = 1. / t16 / t15;
      double t19 = sigma * t18;
      double t21 = lyp_d * t11 + lyp_c;
      double t22 = t21 * t8;
      double t24 = -1. / 72. - 7. / 72. * t22;
      double t34 = 5. / 2. - t22 / 18.;
      double t35 = t34 * sigma;
      double t38 = t22 - 11.;
      double t39 = t38 * sigma;
      double t43 = -t19 * t24 - 3. / 10. * t27 * t30 + t35 * t18 / 8. + 
                   t39 * t18 / 144. - 5. / 24. * t19;

      eps = lyp_A * (t14 * t11 * t43 - t11);

      double t47 = rho * lyp_A;
      double t48 = t10 * t10;
      double t49 = 1. / t48;
      double t50 = t49 * lyp_d;
      double t52 = 1. / t7 / rho;
      double t56 = t55 * t52;
      double t57 = t13 * t11;
      double t58 = t57 * t43;
      double t61 = t14 * t49;
      double t62 = t43 * lyp_d;
      double t66 = t15 * rho;
      double t68 = 1. / t16 / t66;
      double t69 = sigma * t68;
      double t73 = t72 * t49;
      double t75 = 1. / t16 / rho;
      double t78 = t21 * t52 - t73 * t75;
      double t79 = 7. / 216. * t78;
      double t81 = t78 / 54.;
      double t82 = t81 * sigma;
      double t88 = -t78 / 3.;
      double t89 = t88 * sigma;
      double t95 = 8. / 3. * t69 * t24 - t19 * t79 + t82 * t18 / 8. - 
                   t35 * t68 / 3. + t89 * t18 / 144. - t39 * t68 / 54. + 
                   5. / 9. * t69;
      double t98 = -t50 * t52 / 3. + t56 * t58 / 3. + t61 * t62 * t52 / 3. + 
                   t14 * t11 * t95;
      vrho = t47 * t98 + (lyp_A * (t14 * t11 * t43 - t11));

      double t100 = t47 * lyp_B;
      double t107 = -t18 * t24 + t34 * t18 / 8. + t38 * t18 / 144. - 
                    5. / 24. * t18;
      double t108 = t57 * t107;

      vsigma = t100 * t108;
  }

};





struct BuiltinLYP : detail::BuiltinKernelImpl< BuiltinLYP > {

  BuiltinLYP( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinLYP >(p) { }
  
  virtual ~BuiltinLYP() = default;

};

}
