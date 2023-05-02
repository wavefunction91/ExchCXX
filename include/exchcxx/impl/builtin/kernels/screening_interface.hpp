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
    eval_exc_polar( double rho_a, double rho_b, double& eps ) {

    const double rho_s = rho_a + rho_b;
    const double rho_z = rho_a - rho_b;

    double zeta = 0;
    if( rho_s > 0 ) {
      zeta = rho_z / rho_s;
      zeta = safe_min(zeta,  1.);
      zeta = safe_max(zeta, -1.);
    }

    if( rho_s <= traits::dens_tol ) {
      eps = 0.;
    } else {
      traits::eval_exc_polar_impl( rho_a, rho_b, eps );
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


  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double& eps, double& vrho_a,
      double& vrho_b ) {

    const double rho_s = rho_a + rho_b;
    const double rho_z = rho_a - rho_b;

    double zeta = 0;
    if( rho_s > 0 ) {
      zeta = rho_z / rho_s;
      zeta = safe_min(zeta,  1.);
      zeta = safe_max(zeta, -1.);
    }

    if( rho_s <= traits::dens_tol ) {
      eps    = 0.;
      vrho_a = 0.;
      vrho_b = 0.;
    } else {
      traits::eval_exc_vxc_polar_impl( rho_a, rho_b, eps, vrho_a, vrho_b );
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
    eval_exc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps ) {

    const double rho_s = rho_a + rho_b;
    const double rho_z = rho_a - rho_b;

    double zeta = 0;
    if( rho_s > 0 ) {
      zeta = rho_z / rho_s;
      zeta = safe_min(zeta,  1.);
      zeta = safe_max(zeta, -1.);
    }

    if( rho_s <= traits::dens_tol ) {
      eps = 0.;
    } else {
      traits::eval_exc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab,
        sigma_bb, eps );
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

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar( double rho_a, double rho_b, double sigma_aa, 
      double sigma_ab, double sigma_bb, double& eps, double& vrho_a,
      double& vrho_b, double& vsigma_aa, double& vsigma_ab, 
      double& vsigma_bb ) {

    const double rho_s = rho_a + rho_b;
    const double rho_z = rho_a - rho_b;

    double zeta = 0;
    if( rho_s > 0 ) {
      zeta = rho_z / rho_s;
      zeta = safe_min(zeta,  1.);
      zeta = safe_max(zeta, -1.);
    }

    eps       = 0.;
    vrho_a    = 0.;
    vrho_b    = 0.;
    vsigma_aa = 0.;
    vsigma_ab = 0.;
    vsigma_bb = 0.;

    if( rho_s > traits::dens_tol ) {
      traits::eval_exc_vxc_polar_impl( rho_a, rho_b, sigma_aa, sigma_ab, 
        sigma_bb, eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb );
    }

  }

};

}
