#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

namespace ExchCXX {

template <typename KernelType>
struct lda_screening_interface {

  using traits = kernel_traits<KernelType>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double& eps ) {

    if( rho <= traits::dens_tol ) {
      eps = 0.;
    } else {
      traits::eval_exc_unpolar_impl( rho, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double& eps, double& vxc ) {

    if( rho <= traits::dens_tol ) {
      eps = 0.;
      vxc = 0.;
    } else {
      traits::eval_exc_vxc_unpolar_impl( rho, eps, vxc );
    }

  }

};


template <typename KernelType>
struct gga_screening_interface {

  using traits = kernel_traits<KernelType>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    if( rho <= traits::dens_tol ) {
      eps = 0.;
    } else {
      traits::eval_exc_unpolar_impl( rho, sigma, eps );
    }

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

    if( rho <= traits::dens_tol ) {
      eps    = 0.;
      vrho   = 0.;
      vsigma = 0.;
    } else {
      traits::eval_exc_vxc_unpolar_impl( rho, sigma, eps, vrho, vsigma );
    }

  }

};

}
