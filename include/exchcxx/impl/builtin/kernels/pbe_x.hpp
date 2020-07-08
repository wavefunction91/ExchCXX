#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {

template <>
struct kernel_traits<BuiltinPBE_X> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-32;

  static constexpr double pbe_kappa = 0.8040;
  static constexpr double pbe_mu    = 0.2195149727645171;


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
      constexpr double t3 = constants::m_cbrt_one_ov_pi;
      constexpr double t4 = t1 * t3;
      constexpr double t5 = constants::m_cbrt_4;
      constexpr double t6 = t5 * t5;
      constexpr double t7 = t4 * t6;
      constexpr double t8 = constants::m_cbrt_2;
      constexpr double t9 = t8 * t8;

      constexpr double t12 = constants::m_cbrt_6;
      constexpr double t15 = constants::m_cbrt_pi_sq;
      constexpr double t16 = t15 * t15;
      constexpr double t17 = 1. / t16;
      constexpr double t18 = pbe_mu * t12 * t17;

      double t10 = std::pow(rho, constants::m_third); 
      double t20 = rho * rho;
      double t21 = t10 * t10;
      double t23 = 1. / t21 / t20;
      double t27 = pbe_kappa + t18 * sigma * t9 * t23 / 24.;
      double t32 = 1. + pbe_kappa * (1. - pbe_kappa / t27);
      double t34 = t7 * t9 * t10 * t32;
    
      eps = -3. / 16. * t34;

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
      constexpr double t3 = constants::m_cbrt_one_ov_pi;
      constexpr double t4 = t1 * t3;
      constexpr double t5 = constants::m_cbrt_4;
      constexpr double t6 = t5 * t5;
      constexpr double t7 = t4 * t6;
      constexpr double t8 = constants::m_cbrt_2;
      constexpr double t9 = t8 * t8;

      constexpr double t12 = constants::m_cbrt_6;
      constexpr double t15 = constants::m_cbrt_pi_sq;
      constexpr double t16 = t15 * t15;
      constexpr double t17 = 1. / t16;
      constexpr double t18 = pbe_mu * t12 * t17;

      constexpr double t40 = t3 * t6;
      constexpr double t41 = t40 * t8;
      constexpr double t43 = pbe_kappa * pbe_kappa;

      double t10 = std::pow(rho, constants::m_third); 
      double t20 = rho * rho;
      double t21 = t10 * t10;
      double t23 = 1. / t21 / t20;
      double t27 = pbe_kappa + t18 * sigma * t9 * t23 / 24.;
      double t32 = 1. + pbe_kappa * (1. - pbe_kappa / t27);
      double t34 = t7 * t9 * t10 * t32;
    
      eps = -3. / 16. * t34;

      double t42 = 1. / t10 / t20 * t1 * t41;
      double t44 = t27 * t27;
      double t46 = t43 / t44;
      double t49 = t12 * t17 * sigma;
      double t50 = t46 * pbe_mu * t49;

      double t57 = t46 * t18;

      vrho   = -t34 / 4. + t42 * t50 / 24.;
      vsigma = -1. / t10 / rho * t1 * t41 * t57 / 64.;
  }
};


struct BuiltinPBE_X : detail::BuiltinKernelImpl< BuiltinPBE_X > {

  BuiltinPBE_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_X >(p) { }
  
  virtual ~BuiltinPBE_X() = default;

};


}
