/**
 * ExchCXX 
 *
 * Copyright (c) Microsoft Corporation.
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */


#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinPW86_X > :
  public gga_screening_interface< BuiltinPW86_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;


  static constexpr double aa = 1.296;
  static constexpr double bb = 14.0;
  static constexpr double cc = 0.2;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = aa * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t37 = t20 * t20;
    constexpr double t38 = bb * t37;
    constexpr double t40 = 0.1e1 / t23 / t22;
    constexpr double t41 = t38 * t40;
    constexpr double t51 = t22 * t22;
    constexpr double t53 = cc / t51;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t42 = sigma * sigma;
    const double t43 = t42 * t27;
    const double t44 = t30 * t30;
    const double t45 = t44 * rho;
    const double t47 = 0.1e1 / t18 / t45;
    const double t54 = t42 * sigma;
    const double t55 = t44 * t44;
    const double t56 = 0.1e1 / t55;
    const double t60 = 0.1e1 + t26 * t29 * t33 / 0.24e2 + t41 * t43 * t47 / 0.288e3 + t53 * t54 * t56 / 0.576e3;
    const double t61 = safe_math::pow( t60, 0.1e1 / 0.15e2 );
    const double t65 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t61 );


    eps = 0.2e1 * t65;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = aa * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t37 = t20 * t20;
    constexpr double t38 = bb * t37;
    constexpr double t40 = 0.1e1 / t23 / t22;
    constexpr double t41 = t38 * t40;
    constexpr double t51 = t22 * t22;
    constexpr double t53 = cc / t51;
    constexpr double t104 = t25 * t28;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t42 = sigma * sigma;
    const double t43 = t42 * t27;
    const double t44 = t30 * t30;
    const double t45 = t44 * rho;
    const double t47 = 0.1e1 / t18 / t45;
    const double t54 = t42 * sigma;
    const double t55 = t44 * t44;
    const double t56 = 0.1e1 / t55;
    const double t60 = 0.1e1 + t26 * t29 * t33 / 0.24e2 + t41 * t43 * t47 / 0.288e3 + t53 * t54 * t56 / 0.576e3;
    const double t61 = safe_math::pow( t60, 0.1e1 / 0.15e2 );
    const double t65 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t61 );
    const double t66 = 0.1e1 / t31;
    const double t71 = t6 * t17;
    const double t72 = t61 * t61;
    const double t73 = t72 * t72;
    const double t75 = t73 * t73;
    const double t76 = t75 * t73 * t72;
    const double t77 = 0.1e1 / t76;
    const double t78 = t18 * t77;
    const double t79 = t30 * rho;
    const double t81 = 0.1e1 / t31 / t79;
    const double t85 = t44 * t30;
    const double t87 = 0.1e1 / t18 / t85;
    const double t91 = t55 * rho;
    const double t92 = 0.1e1 / t91;
    const double t96 = -t26 * t29 * t81 / 0.9e1 - t41 * t43 * t87 / 0.54e2 - t53 * t54 * t92 / 0.72e2;
    const double t101 = piecewise_functor_3( t2, 0.0, -t6 * t17 * t66 * t61 / 0.8e1 - t71 * t78 * t96 / 0.4e2 );
    const double t108 = sigma * t27;
    const double t115 = t21 * t104 * t33 / 0.24e2 + t41 * t108 * t47 / 0.144e3 + t53 * t42 * t56 / 0.192e3;
    const double t119 = piecewise_functor_3( t2, 0.0, -t71 * t78 * t115 / 0.4e2 );


    eps = 0.2e1 * t65;
    vrho = 0.2e1 * rho * t101 + 0.2e1 * t65;
    vsigma = 0.2e1 * rho * t119;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double sigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = aa * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t37 = t20 * t20;
    constexpr double t38 = bb * t37;
    constexpr double t40 = 0.1e1 / t23 / t22;
    constexpr double t41 = t38 * t40;
    constexpr double t51 = t22 * t22;
    constexpr double t53 = cc / t51;
    constexpr double t104 = t25 * t28;
    constexpr double t191 = t40 * t27;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t42 = sigma * sigma;
    const double t43 = t42 * t27;
    const double t44 = t30 * t30;
    const double t45 = t44 * rho;
    const double t47 = 0.1e1 / t18 / t45;
    const double t54 = t42 * sigma;
    const double t55 = t44 * t44;
    const double t56 = 0.1e1 / t55;
    const double t60 = 0.1e1 + t26 * t29 * t33 / 0.24e2 + t41 * t43 * t47 / 0.288e3 + t53 * t54 * t56 / 0.576e3;
    const double t61 = safe_math::pow( t60, 0.1e1 / 0.15e2 );
    const double t66 = 0.1e1 / t31;
    const double t71 = t6 * t17;
    const double t72 = t61 * t61;
    const double t73 = t72 * t72;
    const double t75 = t73 * t73;
    const double t76 = t75 * t73 * t72;
    const double t77 = 0.1e1 / t76;
    const double t78 = t18 * t77;
    const double t79 = t30 * rho;
    const double t81 = 0.1e1 / t31 / t79;
    const double t85 = t44 * t30;
    const double t87 = 0.1e1 / t18 / t85;
    const double t91 = t55 * rho;
    const double t92 = 0.1e1 / t91;
    const double t96 = -t26 * t29 * t81 / 0.9e1 - t41 * t43 * t87 / 0.54e2 - t53 * t54 * t92 / 0.72e2;
    const double t101 = piecewise_functor_3( t2, 0.0, -t6 * t17 * t66 * t61 / 0.8e1 - t71 * t78 * t96 / 0.4e2 );
    const double t108 = sigma * t27;
    const double t115 = t21 * t104 * t33 / 0.24e2 + t41 * t108 * t47 / 0.144e3 + t53 * t42 * t56 / 0.192e3;
    const double t119 = piecewise_functor_3( t2, 0.0, -t71 * t78 * t115 / 0.4e2 );
    const double t123 = 0.1e1 / t31 / rho;
    const double t128 = t66 * t77;
    const double t133 = 0.1e1 / t76 / t60;
    const double t134 = t18 * t133;
    const double t135 = t96 * t96;
    const double t140 = 0.1e1 / t31 / t44;
    const double t144 = t44 * t79;
    const double t146 = 0.1e1 / t18 / t144;
    const double t151 = 0.1e1 / t55 / t30;
    const double t155 = 0.11e2 / 0.27e2 * t26 * t29 * t140 + 0.19e2 / 0.162e3 * t41 * t43 * t146 + t53 * t54 * t151 / 0.8e1;
    const double t160 = piecewise_functor_3( t2, 0.0, t6 * t17 * t123 * t61 / 0.12e2 - t71 * t128 * t96 / 0.6e2 + 0.7e1 / 0.3e3 * t71 * t134 * t135 - t71 * t78 * t155 / 0.4e2 );
    const double t166 = t115 * t96;
    const double t179 = -t21 * t104 * t81 / 0.9e1 - t41 * t108 * t87 / 0.27e2 - t53 * t42 * t92 / 0.24e2;
    const double t184 = piecewise_functor_3( t2, 0.0, -t71 * t128 * t115 / 0.12e3 + 0.7e1 / 0.3e3 * t71 * t134 * t166 - t71 * t78 * t179 / 0.4e2 );
    const double t187 = t115 * t115;
    const double t198 = t38 * t191 * t47 / 0.144e3 + t53 * sigma * t56 / 0.96e2;
    const double t203 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.3e3 * t71 * t134 * t187 - t71 * t78 * t198 / 0.4e2 );


    v2rho2 = 0.2e1 * rho * t160 + 0.4e1 * t101;
    v2rhosigma = 0.2e1 * rho * t184 + 0.2e1 * t119;
    v2sigma2 = 0.2e1 * rho * t203;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double sigma, double& vrho, double& vsigma, double& v2rho2, double& v2rhosigma, double& v2sigma2 ) {

    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t20 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t6 = t3 / t4;
    constexpr double t21 = aa * t20;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t37 = t20 * t20;
    constexpr double t38 = bb * t37;
    constexpr double t40 = 0.1e1 / t23 / t22;
    constexpr double t41 = t38 * t40;
    constexpr double t51 = t22 * t22;
    constexpr double t53 = cc / t51;
    constexpr double t104 = t25 * t28;
    constexpr double t191 = t40 * t27;


    const double t2 = rho / 0.2e1 <= dens_tol;
    const double t7 = 0.1e1 <= zeta_tol;
    const double t8 = zeta_tol - 0.1e1;
    const double t10 = piecewise_functor_5( t7, t8, t7, -t8, 0.0 );
    const double t11 = 0.1e1 + t10;
    const double t13 = safe_math::cbrt( zeta_tol );
    const double t15 = safe_math::cbrt( t11 );
    const double t17 = piecewise_functor_3( t11 <= zeta_tol, t13 * zeta_tol, t15 * t11 );
    const double t18 = safe_math::cbrt( rho );
    const double t19 = t17 * t18;
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t18 * t18;
    const double t33 = 0.1e1 / t31 / t30;
    const double t42 = sigma * sigma;
    const double t43 = t42 * t27;
    const double t44 = t30 * t30;
    const double t45 = t44 * rho;
    const double t47 = 0.1e1 / t18 / t45;
    const double t54 = t42 * sigma;
    const double t55 = t44 * t44;
    const double t56 = 0.1e1 / t55;
    const double t60 = 0.1e1 + t26 * t29 * t33 / 0.24e2 + t41 * t43 * t47 / 0.288e3 + t53 * t54 * t56 / 0.576e3;
    const double t61 = safe_math::pow( t60, 0.1e1 / 0.15e2 );
    const double t65 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t19 * t61 );
    const double t66 = 0.1e1 / t31;
    const double t71 = t6 * t17;
    const double t72 = t61 * t61;
    const double t73 = t72 * t72;
    const double t75 = t73 * t73;
    const double t76 = t75 * t73 * t72;
    const double t77 = 0.1e1 / t76;
    const double t78 = t18 * t77;
    const double t79 = t30 * rho;
    const double t81 = 0.1e1 / t31 / t79;
    const double t85 = t44 * t30;
    const double t87 = 0.1e1 / t18 / t85;
    const double t91 = t55 * rho;
    const double t92 = 0.1e1 / t91;
    const double t96 = -t26 * t29 * t81 / 0.9e1 - t41 * t43 * t87 / 0.54e2 - t53 * t54 * t92 / 0.72e2;
    const double t101 = piecewise_functor_3( t2, 0.0, -t6 * t17 * t66 * t61 / 0.8e1 - t71 * t78 * t96 / 0.4e2 );
    const double t108 = sigma * t27;
    const double t115 = t21 * t104 * t33 / 0.24e2 + t41 * t108 * t47 / 0.144e3 + t53 * t42 * t56 / 0.192e3;
    const double t119 = piecewise_functor_3( t2, 0.0, -t71 * t78 * t115 / 0.4e2 );
    const double t123 = 0.1e1 / t31 / rho;
    const double t128 = t66 * t77;
    const double t133 = 0.1e1 / t76 / t60;
    const double t134 = t18 * t133;
    const double t135 = t96 * t96;
    const double t140 = 0.1e1 / t31 / t44;
    const double t144 = t44 * t79;
    const double t146 = 0.1e1 / t18 / t144;
    const double t151 = 0.1e1 / t55 / t30;
    const double t155 = 0.11e2 / 0.27e2 * t26 * t29 * t140 + 0.19e2 / 0.162e3 * t41 * t43 * t146 + t53 * t54 * t151 / 0.8e1;
    const double t160 = piecewise_functor_3( t2, 0.0, t6 * t17 * t123 * t61 / 0.12e2 - t71 * t128 * t96 / 0.6e2 + 0.7e1 / 0.3e3 * t71 * t134 * t135 - t71 * t78 * t155 / 0.4e2 );
    const double t166 = t115 * t96;
    const double t179 = -t21 * t104 * t81 / 0.9e1 - t41 * t108 * t87 / 0.27e2 - t53 * t42 * t92 / 0.24e2;
    const double t184 = piecewise_functor_3( t2, 0.0, -t71 * t128 * t115 / 0.12e3 + 0.7e1 / 0.3e3 * t71 * t134 * t166 - t71 * t78 * t179 / 0.4e2 );
    const double t187 = t115 * t115;
    const double t198 = t38 * t191 * t47 / 0.144e3 + t53 * sigma * t56 / 0.96e2;
    const double t203 = piecewise_functor_3( t2, 0.0, 0.7e1 / 0.3e3 * t71 * t134 * t187 - t71 * t78 * t198 / 0.4e2 );


    vrho = 0.2e1 * rho * t101 + 0.2e1 * t65;
    vsigma = 0.2e1 * rho * t119;
    v2rho2 = 0.2e1 * rho * t160 + 0.4e1 * t101;
    v2rhosigma = 0.2e1 * rho * t184 + 0.2e1 * t119;
    v2sigma2 = 0.2e1 * rho * t203;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = aa * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t43 = t28 * t28;
    constexpr double t44 = bb * t43;
    constexpr double t46 = 0.1e1 / t31 / t30;
    constexpr double t56 = t30 * t30;
    constexpr double t58 = cc / t56;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t47 = sigma_aa * sigma_aa;
    const double t48 = t46 * t47;
    const double t49 = t35 * t35;
    const double t50 = t49 * rho_a;
    const double t52 = 0.1e1 / t36 / t50;
    const double t59 = t47 * sigma_aa;
    const double t60 = t49 * t49;
    const double t61 = 0.1e1 / t60;
    const double t65 = 0.1e1 + t29 * t34 * t39 / 0.24e2 + t44 * t48 * t52 / 0.576e3 + t58 * t59 * t61 / 0.2304e4;
    const double t66 = safe_math::pow( t65, 0.1e1 / 0.15e2 );
    const double t70 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t66 );
    const double t71 = rho_b <= dens_tol;
    const double t72 = -t16;
    const double t74 = piecewise_functor_5( t14, t11, t10, t15, t72 * t7 );
    const double t75 = 0.1e1 + t74;
    const double t76 = t75 <= zeta_tol;
    const double t77 = safe_math::cbrt( t75 );
    const double t79 = piecewise_functor_3( t76, t22, t77 * t75 );
    const double t80 = t79 * t26;
    const double t81 = t33 * sigma_bb;
    const double t82 = rho_b * rho_b;
    const double t83 = safe_math::cbrt( rho_b );
    const double t84 = t83 * t83;
    const double t86 = 0.1e1 / t84 / t82;
    const double t90 = sigma_bb * sigma_bb;
    const double t91 = t46 * t90;
    const double t92 = t82 * t82;
    const double t93 = t92 * rho_b;
    const double t95 = 0.1e1 / t83 / t93;
    const double t99 = t90 * sigma_bb;
    const double t100 = t92 * t92;
    const double t101 = 0.1e1 / t100;
    const double t105 = 0.1e1 + t29 * t81 * t86 / 0.24e2 + t44 * t91 * t95 / 0.576e3 + t58 * t99 * t101 / 0.2304e4;
    const double t106 = safe_math::pow( t105, 0.1e1 / 0.15e2 );
    const double t110 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t80 * t106 );


    eps = t70 + t110;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = aa * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t43 = t28 * t28;
    constexpr double t44 = bb * t43;
    constexpr double t46 = 0.1e1 / t31 / t30;
    constexpr double t56 = t30 * t30;
    constexpr double t58 = cc / t56;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t47 = sigma_aa * sigma_aa;
    const double t48 = t46 * t47;
    const double t49 = t35 * t35;
    const double t50 = t49 * rho_a;
    const double t52 = 0.1e1 / t36 / t50;
    const double t59 = t47 * sigma_aa;
    const double t60 = t49 * t49;
    const double t61 = 0.1e1 / t60;
    const double t65 = 0.1e1 + t29 * t34 * t39 / 0.24e2 + t44 * t48 * t52 / 0.576e3 + t58 * t59 * t61 / 0.2304e4;
    const double t66 = safe_math::pow( t65, 0.1e1 / 0.15e2 );
    const double t70 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t66 );
    const double t71 = rho_b <= dens_tol;
    const double t72 = -t16;
    const double t74 = piecewise_functor_5( t14, t11, t10, t15, t72 * t7 );
    const double t75 = 0.1e1 + t74;
    const double t76 = t75 <= zeta_tol;
    const double t77 = safe_math::cbrt( t75 );
    const double t79 = piecewise_functor_3( t76, t22, t77 * t75 );
    const double t80 = t79 * t26;
    const double t81 = t33 * sigma_bb;
    const double t82 = rho_b * rho_b;
    const double t83 = safe_math::cbrt( rho_b );
    const double t84 = t83 * t83;
    const double t86 = 0.1e1 / t84 / t82;
    const double t90 = sigma_bb * sigma_bb;
    const double t91 = t46 * t90;
    const double t92 = t82 * t82;
    const double t93 = t92 * rho_b;
    const double t95 = 0.1e1 / t83 / t93;
    const double t99 = t90 * sigma_bb;
    const double t100 = t92 * t92;
    const double t101 = 0.1e1 / t100;
    const double t105 = 0.1e1 + t29 * t81 * t86 / 0.24e2 + t44 * t91 * t95 / 0.576e3 + t58 * t99 * t101 / 0.2304e4;
    const double t106 = safe_math::pow( t105, 0.1e1 / 0.15e2 );
    const double t110 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t80 * t106 );
    const double t111 = t6 * t6;
    const double t112 = 0.1e1 / t111;
    const double t113 = t16 * t112;
    const double t115 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t113 );
    const double t118 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t115 );
    const double t123 = t26 * t26;
    const double t124 = 0.1e1 / t123;
    const double t128 = t5 * t25 * t124 * t66 / 0.8e1;
    const double t129 = t5 * t25;
    const double t130 = t66 * t66;
    const double t131 = t130 * t130;
    const double t133 = t131 * t131;
    const double t134 = t133 * t131 * t130;
    const double t135 = 0.1e1 / t134;
    const double t136 = t26 * t135;
    const double t137 = t35 * rho_a;
    const double t139 = 0.1e1 / t37 / t137;
    const double t143 = t49 * t35;
    const double t145 = 0.1e1 / t36 / t143;
    const double t149 = t60 * rho_a;
    const double t150 = 0.1e1 / t149;
    const double t154 = -t29 * t34 * t139 / 0.9e1 - t44 * t48 * t145 / 0.108e3 - t58 * t59 * t150 / 0.288e3;
    const double t155 = t136 * t154;
    const double t159 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t118 * t26 * t66 - t128 - t129 * t155 / 0.4e2 );
    const double t160 = t72 * t112;
    const double t162 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t160 );
    const double t165 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t162 );
    const double t173 = t5 * t79 * t124 * t106 / 0.8e1;
    const double t175 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t165 * t26 * t106 - t173 );
    const double t179 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t113 );
    const double t182 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t179 );
    const double t188 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t182 * t26 * t66 - t128 );
    const double t190 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t160 );
    const double t193 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t190 );
    const double t198 = t5 * t79;
    const double t199 = t106 * t106;
    const double t200 = t199 * t199;
    const double t202 = t200 * t200;
    const double t203 = t202 * t200 * t199;
    const double t204 = 0.1e1 / t203;
    const double t205 = t26 * t204;
    const double t206 = t82 * rho_b;
    const double t208 = 0.1e1 / t84 / t206;
    const double t212 = t92 * t82;
    const double t214 = 0.1e1 / t83 / t212;
    const double t218 = t100 * rho_b;
    const double t219 = 0.1e1 / t218;
    const double t223 = -t29 * t81 * t208 / 0.9e1 - t44 * t91 * t214 / 0.108e3 - t58 * t99 * t219 / 0.288e3;
    const double t224 = t205 * t223;
    const double t228 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t193 * t26 * t106 - t173 - t198 * t224 / 0.4e2 );
    const double t234 = t46 * sigma_aa;
    const double t241 = t29 * t33 * t39 / 0.24e2 + t44 * t234 * t52 / 0.288e3 + t58 * t47 * t61 / 0.768e3;
    const double t242 = t136 * t241;
    const double t245 = piecewise_functor_3( t1, 0.0, -t129 * t242 / 0.4e2 );
    const double t249 = t46 * sigma_bb;
    const double t256 = t29 * t33 * t86 / 0.24e2 + t44 * t249 * t95 / 0.288e3 + t58 * t90 * t101 / 0.768e3;
    const double t257 = t205 * t256;
    const double t260 = piecewise_functor_3( t71, 0.0, -t198 * t257 / 0.4e2 );


    eps = t70 + t110;
    vrho_a = t70 + t110 + t6 * ( t159 + t175 );
    vrho_b = t70 + t110 + t6 * ( t188 + t228 );
    vsigma_aa = t6 * t245;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t260;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = aa * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t43 = t28 * t28;
    constexpr double t44 = bb * t43;
    constexpr double t46 = 0.1e1 / t31 / t30;
    constexpr double t56 = t30 * t30;
    constexpr double t58 = cc / t56;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t47 = sigma_aa * sigma_aa;
    const double t48 = t46 * t47;
    const double t49 = t35 * t35;
    const double t50 = t49 * rho_a;
    const double t52 = 0.1e1 / t36 / t50;
    const double t59 = t47 * sigma_aa;
    const double t60 = t49 * t49;
    const double t61 = 0.1e1 / t60;
    const double t65 = 0.1e1 + t29 * t34 * t39 / 0.24e2 + t44 * t48 * t52 / 0.576e3 + t58 * t59 * t61 / 0.2304e4;
    const double t66 = safe_math::pow( t65, 0.1e1 / 0.15e2 );
    const double t71 = rho_b <= dens_tol;
    const double t72 = -t16;
    const double t74 = piecewise_functor_5( t14, t11, t10, t15, t72 * t7 );
    const double t75 = 0.1e1 + t74;
    const double t76 = t75 <= zeta_tol;
    const double t77 = safe_math::cbrt( t75 );
    const double t79 = piecewise_functor_3( t76, t22, t77 * t75 );
    const double t81 = t33 * sigma_bb;
    const double t82 = rho_b * rho_b;
    const double t83 = safe_math::cbrt( rho_b );
    const double t84 = t83 * t83;
    const double t86 = 0.1e1 / t84 / t82;
    const double t90 = sigma_bb * sigma_bb;
    const double t91 = t46 * t90;
    const double t92 = t82 * t82;
    const double t93 = t92 * rho_b;
    const double t95 = 0.1e1 / t83 / t93;
    const double t99 = t90 * sigma_bb;
    const double t100 = t92 * t92;
    const double t101 = 0.1e1 / t100;
    const double t105 = 0.1e1 + t29 * t81 * t86 / 0.24e2 + t44 * t91 * t95 / 0.576e3 + t58 * t99 * t101 / 0.2304e4;
    const double t106 = safe_math::pow( t105, 0.1e1 / 0.15e2 );
    const double t111 = t6 * t6;
    const double t112 = 0.1e1 / t111;
    const double t113 = t16 * t112;
    const double t115 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t113 );
    const double t118 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t115 );
    const double t123 = t26 * t26;
    const double t124 = 0.1e1 / t123;
    const double t128 = t5 * t25 * t124 * t66 / 0.8e1;
    const double t129 = t5 * t25;
    const double t130 = t66 * t66;
    const double t131 = t130 * t130;
    const double t133 = t131 * t131;
    const double t134 = t133 * t131 * t130;
    const double t135 = 0.1e1 / t134;
    const double t136 = t26 * t135;
    const double t137 = t35 * rho_a;
    const double t139 = 0.1e1 / t37 / t137;
    const double t143 = t49 * t35;
    const double t145 = 0.1e1 / t36 / t143;
    const double t149 = t60 * rho_a;
    const double t150 = 0.1e1 / t149;
    const double t154 = -t29 * t34 * t139 / 0.9e1 - t44 * t48 * t145 / 0.108e3 - t58 * t59 * t150 / 0.288e3;
    const double t155 = t136 * t154;
    const double t159 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t118 * t26 * t66 - t128 - t129 * t155 / 0.4e2 );
    const double t160 = t72 * t112;
    const double t162 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t160 );
    const double t165 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t162 );
    const double t173 = t5 * t79 * t124 * t106 / 0.8e1;
    const double t175 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t165 * t26 * t106 - t173 );
    const double t179 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t113 );
    const double t182 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t179 );
    const double t188 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t182 * t26 * t66 - t128 );
    const double t190 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t160 );
    const double t193 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t190 );
    const double t198 = t5 * t79;
    const double t199 = t106 * t106;
    const double t200 = t199 * t199;
    const double t202 = t200 * t200;
    const double t203 = t202 * t200 * t199;
    const double t204 = 0.1e1 / t203;
    const double t205 = t26 * t204;
    const double t206 = t82 * rho_b;
    const double t208 = 0.1e1 / t84 / t206;
    const double t212 = t92 * t82;
    const double t214 = 0.1e1 / t83 / t212;
    const double t218 = t100 * rho_b;
    const double t219 = 0.1e1 / t218;
    const double t223 = -t29 * t81 * t208 / 0.9e1 - t44 * t91 * t214 / 0.108e3 - t58 * t99 * t219 / 0.288e3;
    const double t224 = t205 * t223;
    const double t228 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t193 * t26 * t106 - t173 - t198 * t224 / 0.4e2 );
    const double t234 = t46 * sigma_aa;
    const double t241 = t29 * t33 * t39 / 0.24e2 + t44 * t234 * t52 / 0.288e3 + t58 * t47 * t61 / 0.768e3;
    const double t242 = t136 * t241;
    const double t245 = piecewise_functor_3( t1, 0.0, -t129 * t242 / 0.4e2 );
    const double t249 = t46 * sigma_bb;
    const double t256 = t29 * t33 * t86 / 0.24e2 + t44 * t249 * t95 / 0.288e3 + t58 * t90 * t101 / 0.768e3;
    const double t257 = t205 * t256;
    const double t260 = piecewise_functor_3( t71, 0.0, -t198 * t257 / 0.4e2 );
    const double t263 = t23 * t23;
    const double t264 = 0.1e1 / t263;
    const double t265 = t115 * t115;
    const double t268 = t111 * t6;
    const double t269 = 0.1e1 / t268;
    const double t270 = t16 * t269;
    const double t273 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t112 + 0.2e1 * t270 );
    const double t277 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t264 * t265 + 0.4e1 / 0.3e1 * t23 * t273 );
    const double t284 = t5 * t118 * t124 * t66;
    const double t286 = t5 * t118;
    const double t290 = 0.1e1 / t123 / t6;
    const double t294 = t5 * t25 * t290 * t66 / 0.12e2;
    const double t295 = t124 * t135;
    const double t296 = t295 * t154;
    const double t297 = t129 * t296;
    const double t300 = 0.1e1 / t134 / t65;
    const double t301 = t26 * t300;
    const double t302 = t154 * t154;
    const double t303 = t301 * t302;
    const double t307 = 0.1e1 / t37 / t49;
    const double t313 = 0.1e1 / t36 / t49 / t137;
    const double t318 = 0.1e1 / t60 / t35;
    const double t322 = 0.11e2 / 0.27e2 * t29 * t34 * t307 + 0.19e2 / 0.324e3 * t44 * t48 * t313 + t58 * t59 * t318 / 0.32e2;
    const double t323 = t136 * t322;
    const double t327 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t277 * t26 * t66 - t284 / 0.4e1 - t286 * t155 / 0.2e2 + t294 - t297 / 0.6e2 + 0.7e1 / 0.3e3 * t129 * t303 - t129 * t323 / 0.4e2 );
    const double t328 = t77 * t77;
    const double t329 = 0.1e1 / t328;
    const double t330 = t162 * t162;
    const double t333 = t72 * t269;
    const double t336 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t112 + 0.2e1 * t333 );
    const double t340 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.9e1 * t329 * t330 + 0.4e1 / 0.3e1 * t77 * t336 );
    const double t347 = t5 * t165 * t124 * t106;
    const double t352 = t5 * t79 * t290 * t106 / 0.12e2;
    const double t354 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t340 * t26 * t106 - t347 / 0.4e1 + t352 );
    const double t372 = t5 * t182 * t124 * t66;
    const double t396 = t5 * t193 * t124 * t106;
    const double t402 = t124 * t204;
    const double t403 = t402 * t223;
    const double t404 = t198 * t403;
    const double t412 = t179 * t179;
    const double t417 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t112 + 0.2e1 * t270 );
    const double t421 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t264 * t412 + 0.4e1 / 0.3e1 * t23 * t417 );
    const double t428 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t421 * t26 * t66 - t372 / 0.4e1 + t294 );
    const double t429 = t190 * t190;
    const double t434 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t112 + 0.2e1 * t333 );
    const double t438 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.9e1 * t329 * t429 + 0.4e1 / 0.3e1 * t77 * t434 );
    const double t444 = t5 * t193;
    const double t449 = 0.1e1 / t203 / t105;
    const double t450 = t26 * t449;
    const double t451 = t223 * t223;
    const double t452 = t450 * t451;
    const double t456 = 0.1e1 / t84 / t92;
    const double t462 = 0.1e1 / t83 / t92 / t206;
    const double t467 = 0.1e1 / t100 / t82;
    const double t471 = 0.11e2 / 0.27e2 * t29 * t81 * t456 + 0.19e2 / 0.324e3 * t44 * t91 * t462 + t58 * t99 * t467 / 0.32e2;
    const double t472 = t205 * t471;
    const double t476 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t438 * t26 * t106 - t396 / 0.4e1 - t444 * t224 / 0.2e2 + t352 - t404 / 0.6e2 + 0.7e1 / 0.3e3 * t198 * t452 - t198 * t472 / 0.4e2 );
    const double t481 = t295 * t241;
    const double t483 = t129 * t481 / 0.12e3;
    const double t484 = t241 * t154;
    const double t485 = t301 * t484;
    const double t497 = -t29 * t33 * t139 / 0.9e1 - t44 * t234 * t145 / 0.54e2 - t58 * t47 * t150 / 0.96e2;
    const double t498 = t136 * t497;
    const double t502 = piecewise_functor_3( t1, 0.0, -t286 * t242 / 0.4e2 - t483 + 0.7e1 / 0.3e3 * t129 * t485 - t129 * t498 / 0.4e2 );
    const double t506 = t402 * t256;
    const double t508 = t198 * t506 / 0.12e3;
    const double t519 = t256 * t223;
    const double t520 = t450 * t519;
    const double t532 = -t29 * t33 * t208 / 0.9e1 - t44 * t249 * t214 / 0.54e2 - t58 * t90 * t219 / 0.96e2;
    const double t533 = t205 * t532;
    const double t537 = piecewise_functor_3( t71, 0.0, -t444 * t257 / 0.4e2 - t508 + 0.7e1 / 0.3e3 * t198 * t520 - t198 * t533 / 0.4e2 );
    const double t539 = t241 * t241;
    const double t540 = t301 * t539;
    const double t549 = t44 * t46 * t52 / 0.288e3 + t58 * sigma_aa * t61 / 0.384e3;
    const double t550 = t136 * t549;
    const double t554 = piecewise_functor_3( t1, 0.0, 0.7e1 / 0.3e3 * t129 * t540 - t129 * t550 / 0.4e2 );
    const double t555 = t256 * t256;
    const double t556 = t450 * t555;
    const double t565 = t44 * t46 * t95 / 0.288e3 + t58 * sigma_bb * t101 / 0.384e3;
    const double t566 = t205 * t565;
    const double t570 = piecewise_functor_3( t71, 0.0, 0.7e1 / 0.3e3 * t198 * t556 - t198 * t566 / 0.4e2 );


    v2rho2_aa = 0.2e1 * t159 + 0.2e1 * t175 + t6 * ( t327 + t354 );
    v2rho2_bb = 0.2e1 * t188 + 0.2e1 * t228 + t6 * ( t428 + t476 );
    v2rhosigma_a_aa = t6 * t502 + t245;
    v2rhosigma_b_bb = t6 * t537 + t260;
    v2sigma2_aa_aa = t6 * t554;
    v2sigma2_bb_bb = t6 * t570;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;



  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb ) {

    (void)(sigma_ab);
    constexpr double t2 = constants::m_cbrt_3;
    constexpr double t3 = constants::m_cbrt_pi;
    constexpr double t28 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t5 = t2 / t3;
    constexpr double t29 = aa * t28;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t43 = t28 * t28;
    constexpr double t44 = bb * t43;
    constexpr double t46 = 0.1e1 / t31 / t30;
    constexpr double t56 = t30 * t30;
    constexpr double t58 = cc / t56;


    const double t1 = rho_a <= dens_tol;
    const double t6 = rho_a + rho_b;
    const double t7 = 0.1e1 / t6;
    const double t10 = 0.2e1 * rho_a * t7 <= zeta_tol;
    const double t11 = zeta_tol - 0.1e1;
    const double t14 = 0.2e1 * rho_b * t7 <= zeta_tol;
    const double t15 = -t11;
    const double t16 = rho_a - rho_b;
    const double t18 = piecewise_functor_5( t10, t11, t14, t15, t16 * t7 );
    const double t19 = 0.1e1 + t18;
    const double t20 = t19 <= zeta_tol;
    const double t21 = safe_math::cbrt( zeta_tol );
    const double t22 = t21 * zeta_tol;
    const double t23 = safe_math::cbrt( t19 );
    const double t25 = piecewise_functor_3( t20, t22, t23 * t19 );
    const double t26 = safe_math::cbrt( t6 );
    const double t27 = t25 * t26;
    const double t34 = t33 * sigma_aa;
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t47 = sigma_aa * sigma_aa;
    const double t48 = t46 * t47;
    const double t49 = t35 * t35;
    const double t50 = t49 * rho_a;
    const double t52 = 0.1e1 / t36 / t50;
    const double t59 = t47 * sigma_aa;
    const double t60 = t49 * t49;
    const double t61 = 0.1e1 / t60;
    const double t65 = 0.1e1 + t29 * t34 * t39 / 0.24e2 + t44 * t48 * t52 / 0.576e3 + t58 * t59 * t61 / 0.2304e4;
    const double t66 = safe_math::pow( t65, 0.1e1 / 0.15e2 );
    const double t70 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t27 * t66 );
    const double t71 = rho_b <= dens_tol;
    const double t72 = -t16;
    const double t74 = piecewise_functor_5( t14, t11, t10, t15, t72 * t7 );
    const double t75 = 0.1e1 + t74;
    const double t76 = t75 <= zeta_tol;
    const double t77 = safe_math::cbrt( t75 );
    const double t79 = piecewise_functor_3( t76, t22, t77 * t75 );
    const double t80 = t79 * t26;
    const double t81 = t33 * sigma_bb;
    const double t82 = rho_b * rho_b;
    const double t83 = safe_math::cbrt( rho_b );
    const double t84 = t83 * t83;
    const double t86 = 0.1e1 / t84 / t82;
    const double t90 = sigma_bb * sigma_bb;
    const double t91 = t46 * t90;
    const double t92 = t82 * t82;
    const double t93 = t92 * rho_b;
    const double t95 = 0.1e1 / t83 / t93;
    const double t99 = t90 * sigma_bb;
    const double t100 = t92 * t92;
    const double t101 = 0.1e1 / t100;
    const double t105 = 0.1e1 + t29 * t81 * t86 / 0.24e2 + t44 * t91 * t95 / 0.576e3 + t58 * t99 * t101 / 0.2304e4;
    const double t106 = safe_math::pow( t105, 0.1e1 / 0.15e2 );
    const double t110 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t80 * t106 );
    const double t111 = t6 * t6;
    const double t112 = 0.1e1 / t111;
    const double t113 = t16 * t112;
    const double t115 = piecewise_functor_5( t10, 0.0, t14, 0.0, t7 - t113 );
    const double t118 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t115 );
    const double t123 = t26 * t26;
    const double t124 = 0.1e1 / t123;
    const double t128 = t5 * t25 * t124 * t66 / 0.8e1;
    const double t129 = t5 * t25;
    const double t130 = t66 * t66;
    const double t131 = t130 * t130;
    const double t133 = t131 * t131;
    const double t134 = t133 * t131 * t130;
    const double t135 = 0.1e1 / t134;
    const double t136 = t26 * t135;
    const double t137 = t35 * rho_a;
    const double t139 = 0.1e1 / t37 / t137;
    const double t143 = t49 * t35;
    const double t145 = 0.1e1 / t36 / t143;
    const double t149 = t60 * rho_a;
    const double t150 = 0.1e1 / t149;
    const double t154 = -t29 * t34 * t139 / 0.9e1 - t44 * t48 * t145 / 0.108e3 - t58 * t59 * t150 / 0.288e3;
    const double t155 = t136 * t154;
    const double t159 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t118 * t26 * t66 - t128 - t129 * t155 / 0.4e2 );
    const double t160 = t72 * t112;
    const double t162 = piecewise_functor_5( t14, 0.0, t10, 0.0, -t7 - t160 );
    const double t165 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t162 );
    const double t173 = t5 * t79 * t124 * t106 / 0.8e1;
    const double t175 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t165 * t26 * t106 - t173 );
    const double t179 = piecewise_functor_5( t10, 0.0, t14, 0.0, -t7 - t113 );
    const double t182 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.3e1 * t23 * t179 );
    const double t188 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t182 * t26 * t66 - t128 );
    const double t190 = piecewise_functor_5( t14, 0.0, t10, 0.0, t7 - t160 );
    const double t193 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.3e1 * t77 * t190 );
    const double t198 = t5 * t79;
    const double t199 = t106 * t106;
    const double t200 = t199 * t199;
    const double t202 = t200 * t200;
    const double t203 = t202 * t200 * t199;
    const double t204 = 0.1e1 / t203;
    const double t205 = t26 * t204;
    const double t206 = t82 * rho_b;
    const double t208 = 0.1e1 / t84 / t206;
    const double t212 = t92 * t82;
    const double t214 = 0.1e1 / t83 / t212;
    const double t218 = t100 * rho_b;
    const double t219 = 0.1e1 / t218;
    const double t223 = -t29 * t81 * t208 / 0.9e1 - t44 * t91 * t214 / 0.108e3 - t58 * t99 * t219 / 0.288e3;
    const double t224 = t205 * t223;
    const double t228 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t193 * t26 * t106 - t173 - t198 * t224 / 0.4e2 );
    const double t234 = t46 * sigma_aa;
    const double t241 = t29 * t33 * t39 / 0.24e2 + t44 * t234 * t52 / 0.288e3 + t58 * t47 * t61 / 0.768e3;
    const double t242 = t136 * t241;
    const double t245 = piecewise_functor_3( t1, 0.0, -t129 * t242 / 0.4e2 );
    const double t249 = t46 * sigma_bb;
    const double t256 = t29 * t33 * t86 / 0.24e2 + t44 * t249 * t95 / 0.288e3 + t58 * t90 * t101 / 0.768e3;
    const double t257 = t205 * t256;
    const double t260 = piecewise_functor_3( t71, 0.0, -t198 * t257 / 0.4e2 );
    const double t263 = t23 * t23;
    const double t264 = 0.1e1 / t263;
    const double t265 = t115 * t115;
    const double t268 = t111 * t6;
    const double t269 = 0.1e1 / t268;
    const double t270 = t16 * t269;
    const double t273 = piecewise_functor_5( t10, 0.0, t14, 0.0, -0.2e1 * t112 + 0.2e1 * t270 );
    const double t277 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t264 * t265 + 0.4e1 / 0.3e1 * t23 * t273 );
    const double t284 = t5 * t118 * t124 * t66;
    const double t286 = t5 * t118;
    const double t290 = 0.1e1 / t123 / t6;
    const double t294 = t5 * t25 * t290 * t66 / 0.12e2;
    const double t295 = t124 * t135;
    const double t296 = t295 * t154;
    const double t297 = t129 * t296;
    const double t300 = 0.1e1 / t134 / t65;
    const double t301 = t26 * t300;
    const double t302 = t154 * t154;
    const double t303 = t301 * t302;
    const double t307 = 0.1e1 / t37 / t49;
    const double t313 = 0.1e1 / t36 / t49 / t137;
    const double t318 = 0.1e1 / t60 / t35;
    const double t322 = 0.11e2 / 0.27e2 * t29 * t34 * t307 + 0.19e2 / 0.324e3 * t44 * t48 * t313 + t58 * t59 * t318 / 0.32e2;
    const double t323 = t136 * t322;
    const double t327 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t277 * t26 * t66 - t284 / 0.4e1 - t286 * t155 / 0.2e2 + t294 - t297 / 0.6e2 + 0.7e1 / 0.3e3 * t129 * t303 - t129 * t323 / 0.4e2 );
    const double t328 = t77 * t77;
    const double t329 = 0.1e1 / t328;
    const double t330 = t162 * t162;
    const double t333 = t72 * t269;
    const double t336 = piecewise_functor_5( t14, 0.0, t10, 0.0, 0.2e1 * t112 + 0.2e1 * t333 );
    const double t340 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.9e1 * t329 * t330 + 0.4e1 / 0.3e1 * t77 * t336 );
    const double t347 = t5 * t165 * t124 * t106;
    const double t352 = t5 * t79 * t290 * t106 / 0.12e2;
    const double t354 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t340 * t26 * t106 - t347 / 0.4e1 + t352 );
    const double t372 = t5 * t182 * t124 * t66;
    const double t396 = t5 * t193 * t124 * t106;
    const double t402 = t124 * t204;
    const double t403 = t402 * t223;
    const double t404 = t198 * t403;
    const double t412 = t179 * t179;
    const double t417 = piecewise_functor_5( t10, 0.0, t14, 0.0, 0.2e1 * t112 + 0.2e1 * t270 );
    const double t421 = piecewise_functor_3( t20, 0.0, 0.4e1 / 0.9e1 * t264 * t412 + 0.4e1 / 0.3e1 * t23 * t417 );
    const double t428 = piecewise_functor_3( t1, 0.0, -0.3e1 / 0.8e1 * t5 * t421 * t26 * t66 - t372 / 0.4e1 + t294 );
    const double t429 = t190 * t190;
    const double t434 = piecewise_functor_5( t14, 0.0, t10, 0.0, -0.2e1 * t112 + 0.2e1 * t333 );
    const double t438 = piecewise_functor_3( t76, 0.0, 0.4e1 / 0.9e1 * t329 * t429 + 0.4e1 / 0.3e1 * t77 * t434 );
    const double t444 = t5 * t193;
    const double t449 = 0.1e1 / t203 / t105;
    const double t450 = t26 * t449;
    const double t451 = t223 * t223;
    const double t452 = t450 * t451;
    const double t456 = 0.1e1 / t84 / t92;
    const double t462 = 0.1e1 / t83 / t92 / t206;
    const double t467 = 0.1e1 / t100 / t82;
    const double t471 = 0.11e2 / 0.27e2 * t29 * t81 * t456 + 0.19e2 / 0.324e3 * t44 * t91 * t462 + t58 * t99 * t467 / 0.32e2;
    const double t472 = t205 * t471;
    const double t476 = piecewise_functor_3( t71, 0.0, -0.3e1 / 0.8e1 * t5 * t438 * t26 * t106 - t396 / 0.4e1 - t444 * t224 / 0.2e2 + t352 - t404 / 0.6e2 + 0.7e1 / 0.3e3 * t198 * t452 - t198 * t472 / 0.4e2 );
    const double t481 = t295 * t241;
    const double t483 = t129 * t481 / 0.12e3;
    const double t484 = t241 * t154;
    const double t485 = t301 * t484;
    const double t497 = -t29 * t33 * t139 / 0.9e1 - t44 * t234 * t145 / 0.54e2 - t58 * t47 * t150 / 0.96e2;
    const double t498 = t136 * t497;
    const double t502 = piecewise_functor_3( t1, 0.0, -t286 * t242 / 0.4e2 - t483 + 0.7e1 / 0.3e3 * t129 * t485 - t129 * t498 / 0.4e2 );
    const double t506 = t402 * t256;
    const double t508 = t198 * t506 / 0.12e3;
    const double t519 = t256 * t223;
    const double t520 = t450 * t519;
    const double t532 = -t29 * t33 * t208 / 0.9e1 - t44 * t249 * t214 / 0.54e2 - t58 * t90 * t219 / 0.96e2;
    const double t533 = t205 * t532;
    const double t537 = piecewise_functor_3( t71, 0.0, -t444 * t257 / 0.4e2 - t508 + 0.7e1 / 0.3e3 * t198 * t520 - t198 * t533 / 0.4e2 );
    const double t539 = t241 * t241;
    const double t540 = t301 * t539;
    const double t549 = t44 * t46 * t52 / 0.288e3 + t58 * sigma_aa * t61 / 0.384e3;
    const double t550 = t136 * t549;
    const double t554 = piecewise_functor_3( t1, 0.0, 0.7e1 / 0.3e3 * t129 * t540 - t129 * t550 / 0.4e2 );
    const double t555 = t256 * t256;
    const double t556 = t450 * t555;
    const double t565 = t44 * t46 * t95 / 0.288e3 + t58 * sigma_bb * t101 / 0.384e3;
    const double t566 = t205 * t565;
    const double t570 = piecewise_functor_3( t71, 0.0, 0.7e1 / 0.3e3 * t198 * t556 - t198 * t566 / 0.4e2 );


    vrho_a = t70 + t110 + t6 * ( t159 + t175 );
    vrho_b = t70 + t110 + t6 * ( t188 + t228 );
    vsigma_aa = t6 * t245;
    vsigma_ab = 0.e0;
    vsigma_bb = t6 * t260;
    v2rho2_aa = 0.2e1 * t159 + 0.2e1 * t175 + t6 * ( t327 + t354 );
    v2rho2_bb = 0.2e1 * t188 + 0.2e1 * t228 + t6 * ( t428 + t476 );
    v2rhosigma_a_aa = t6 * t502 + t245;
    v2rhosigma_b_bb = t6 * t537 + t260;
    v2sigma2_aa_aa = t6 * t554;
    v2sigma2_bb_bb = t6 * t570;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;



  }


};

struct BuiltinPW86_X : detail::BuiltinKernelImpl< BuiltinPW86_X > {

  BuiltinPW86_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPW86_X >(p) { }
  
  virtual ~BuiltinPW86_X() = default;

};



} // namespace ExchCXX
