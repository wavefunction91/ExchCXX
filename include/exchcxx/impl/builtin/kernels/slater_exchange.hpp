#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>


namespace ExchCXX {

template <>
struct kernel_traits<BuiltinSlaterExchange> :
  public lda_screening_interface< BuiltinSlaterExchange > {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-24;

  BUILTIN_KERNEL_EVAL_RETURN 
    eval_exc_unpolar_impl( double rho, double& eps ) {

    constexpr double alpha = 1.;

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;

    double t11 = std::pow( rho, constants::m_third );
    double t13 = t5 * t10 * t11;

    eps = - 3. / 16. * t13;

  }



  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double& eps, double& vxc ) {


    constexpr double alpha = 1.;

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;

    double t11 = std::pow( rho, constants::m_third );
    double t13 = t5 * t10 * t11;

    eps = - 3. / 16. * t13;
    vxc = -t13 / 4.;

  }

};



struct BuiltinSlaterExchange : detail::BuiltinKernelImpl< BuiltinSlaterExchange > {

  BuiltinSlaterExchange( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSlaterExchange >(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

};

}
