/**
 * ExchCXX 
 *
 * Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). 
 *
 * Portions Copyright (c) Microsoft Corporation.
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
struct kernel_traits< BuiltinPKZB_X > :
  public mgga_screening_interface< BuiltinPKZB_X > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = true;
  static constexpr bool needs_laplacian = false;
  static constexpr bool is_kedf = false;
  static constexpr bool is_epc  = false;

  static constexpr double dens_tol  = 1e-15;
  static constexpr double zeta_tol  = 1e-15;
  static constexpr double sigma_tol  = 1.0000000000000027e-20;
  static constexpr double tau_tol = is_kedf ? 0.0 : 1e-20;




  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps ) {

    (void)(lapl);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t71 = piecewise_functor_3( t3, 0.0, -0.3e1 / 0.8e1 * t7 * t18 * t19 * t67 );


    eps = 0.2e1 * t71;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& eps, double& vrho, double& vsigma, double& vlapl, double& vtau ) {

    (void)(lapl);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;
    constexpr double t117 = t25 * t28;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t71 = piecewise_functor_3( t3, 0.0, -0.3e1 / 0.8e1 * t7 * t18 * t19 * t67 );
    const double t72 = 0.1e1 / t31;
    const double t77 = t4 * t18;
    const double t78 = t64 * t64;
    const double t79 = 0.1e1 / t78;
    const double t80 = t19 * t79;
    const double t81 = t30 * rho;
    const double t83 = 0.1e1 / t31 / t81;
    const double t84 = t29 * t83;
    const double t85 = t26 * t84;
    const double t91 = -0.5e1 / 0.12e2 * t26 * t37 * t33 + t85 / 0.108e3;
    const double t94 = t91 * t21;
    const double t95 = t94 * t25;
    const double t100 = t57 * t30;
    const double t102 = 0.1e1 / t19 / t100;
    const double t106 = -0.1e2 / 0.729e3 * t85 + 0.292e3 / 0.2025e4 * t44 * t91 - 0.73e2 / 0.972e4 * t95 * t34 + 0.73e2 / 0.3645e4 * t48 * t84 - 0.24436516267640330302e-2 * t54 * t56 * t102;
    const double t111 = piecewise_functor_3( t3, 0.0, -t7 * t18 * t72 * t67 / 0.8e1 - 0.16551095363746320496e0 * t77 * t80 * t106 );
    const double t118 = t117 * t33;
    const double t119 = t47 * t118;
    const double t123 = t54 * t27 * t60 * sigma;
    const double t125 = 0.5e1 / 0.972e3 * t26 * t28 * t33 - 0.146e3 / 0.18225e5 * t119 + 0.96852413827153753492e-3 * t123;
    const double t129 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t125 );
    const double t131 = t117 * t39;
    const double t140 = 0.73e2 / 0.2025e4 * t47 * t131 - 0.73e2 / 0.1944e5 * t54 * t27 / t19 / t57 * sigma;
    const double t144 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t140 );


    eps = 0.2e1 * t71;
    vrho = 0.2e1 * rho * t111 + 0.2e1 * t71;
    vsigma = 0.2e1 * rho * t129;
    vlapl = 0.e0;
    vtau = 0.2e1 * rho * t144;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, double& v2lapl2, double& v2lapltau, double& v2tau2 ) {

    (void)(lapl);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;
    constexpr double t117 = t25 * t28;
    constexpr double t252 = t53 * t27;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t72 = 0.1e1 / t31;
    const double t77 = t4 * t18;
    const double t78 = t64 * t64;
    const double t79 = 0.1e1 / t78;
    const double t80 = t19 * t79;
    const double t81 = t30 * rho;
    const double t83 = 0.1e1 / t31 / t81;
    const double t84 = t29 * t83;
    const double t85 = t26 * t84;
    const double t91 = -0.5e1 / 0.12e2 * t26 * t37 * t33 + t85 / 0.108e3;
    const double t94 = t91 * t21;
    const double t95 = t94 * t25;
    const double t100 = t57 * t30;
    const double t102 = 0.1e1 / t19 / t100;
    const double t106 = -0.1e2 / 0.729e3 * t85 + 0.292e3 / 0.2025e4 * t44 * t91 - 0.73e2 / 0.972e4 * t95 * t34 + 0.73e2 / 0.3645e4 * t48 * t84 - 0.24436516267640330302e-2 * t54 * t56 * t102;
    const double t111 = piecewise_functor_3( t3, 0.0, -t7 * t18 * t72 * t67 / 0.8e1 - 0.16551095363746320496e0 * t77 * t80 * t106 );
    const double t118 = t117 * t33;
    const double t119 = t47 * t118;
    const double t123 = t54 * t27 * t60 * sigma;
    const double t125 = 0.5e1 / 0.972e3 * t26 * t28 * t33 - 0.146e3 / 0.18225e5 * t119 + 0.96852413827153753492e-3 * t123;
    const double t129 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t125 );
    const double t131 = t117 * t39;
    const double t140 = 0.73e2 / 0.2025e4 * t47 * t131 - 0.73e2 / 0.1944e5 * t54 * t27 / t19 / t57 * sigma;
    const double t144 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t140 );
    const double t151 = t72 * t79;
    const double t156 = 0.1e1 / t78 / t64;
    const double t157 = t19 * t156;
    const double t158 = t106 * t106;
    const double t163 = 0.1e1 / t31 / t57;
    const double t164 = t29 * t163;
    const double t165 = t26 * t164;
    const double t167 = t91 * t91;
    const double t173 = 0.1e2 / 0.9e1 * t26 * t37 * t83 - 0.11e2 / 0.324e3 * t165;
    const double t176 = t173 * t21;
    const double t177 = t176 * t25;
    const double t184 = t57 * t81;
    const double t186 = 0.1e1 / t19 / t184;
    const double t190 = 0.11e3 / 0.2187e4 * t165 + 0.292e3 / 0.2025e4 * t167 + 0.292e3 / 0.2025e4 * t44 * t173 - 0.73e2 / 0.972e4 * t177 * t34 + 0.146e3 / 0.3645e4 * t95 * t84 - 0.803e3 / 0.10935e5 * t48 * t164 + 0.15476460302838875858e-1 * t54 * t56 * t186;
    const double t195 = piecewise_functor_3( t3, 0.0, t7 * t18 * t39 * t67 / 0.12e2 - 0.11034063575830880331e0 * t77 * t151 * t106 + 0.33102190727492640992e0 * t77 * t157 * t158 - 0.16551095363746320496e0 * t77 * t80 * t190 );
    const double t201 = t77 * t19;
    const double t202 = t156 * t125;
    const double t203 = t202 * t106;
    const double t209 = t94 * t118;
    const double t211 = t117 * t83;
    const double t212 = t47 * t211;
    const double t216 = t54 * t27 * t102 * sigma;
    const double t218 = -0.1e2 / 0.729e3 * t26 * t28 * t83 - 0.146e3 / 0.18225e5 * t209 + 0.1168e4 / 0.54675e5 * t212 - 0.51654620707815335196e-2 * t216;
    const double t223 = piecewise_functor_3( t3, 0.0, -0.55170317879154401653e-1 * t77 * t151 * t125 + 0.33102190727492640992e0 * t201 * t203 - 0.16551095363746320496e0 * t77 * t80 * t218 );
    const double t229 = t156 * t140;
    const double t230 = t229 * t106;
    const double t237 = 0.73e2 / 0.2025e4 * t94 * t131 - 0.73e2 / 0.1215e4 * t119 + 0.949e3 / 0.5832e5 * t123;
    const double t242 = piecewise_functor_3( t3, 0.0, -0.55170317879154401653e-1 * t77 * t151 * t140 + 0.33102190727492640992e0 * t201 * t230 - 0.16551095363746320496e0 * t77 * t80 * t237 );
    const double t245 = t125 * t125;
    const double t249 = 0.1e1 / t58;
    const double t253 = t79 * t51 * t252;
    const double t254 = t77 * t249 * t253;
    const double t257 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t77 * t157 * t245 - 0.16950901996748250202e-3 * t254 );
    const double t259 = t229 * t125;
    const double t262 = 0.1e1 / t57;
    const double t264 = t77 * t262 * t253;
    const double t267 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t201 * t259 + 0.66295196793057964127e-3 * t264 );
    const double t269 = t140 * t140;
    const double t273 = 0.1e1 / t81;
    const double t278 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t77 * t157 * t269 - 0.29832838556876083857e-2 * t77 * t273 * t253 );


    v2rho2 = 0.2e1 * rho * t195 + 0.4e1 * t111;
    v2rhosigma = 0.2e1 * rho * t223 + 0.2e1 * t129;
    v2rholapl = 0.e0;
    v2rhotau = 0.2e1 * rho * t242 + 0.2e1 * t144;
    v2sigma2 = 0.2e1 * rho * t257;
    v2sigmalapl = 0.e0;
    v2sigmatau = 0.2e1 * rho * t267;
    v2lapl2 = 0.e0;
    v2lapltau = 0.e0;
    v2tau2 = 0.2e1 * rho * t278;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_unpolar_impl( double rho, double sigma, double lapl, double tau, double& vrho, double& vsigma, double& vlapl, double& vtau, double& v2rho2, double& v2rhosigma, double& v2rholapl, double& v2rhotau, double& v2sigma2, double& v2sigmalapl, double& v2sigmatau, double& v2lapl2, double& v2lapltau, double& v2tau2 ) {

    (void)(lapl);
    constexpr double t4 = constants::m_cbrt_3;
    constexpr double t5 = constants::m_cbrt_pi;
    constexpr double t21 = constants::m_cbrt_6;
    constexpr double t22 = constants::m_pi_sq;
    constexpr double t23 = constants::m_cbrt_pi_sq;
    constexpr double t27 = constants::m_cbrt_2;
    constexpr double t7 = t4 / t5;
    constexpr double t24 = t23 * t23;
    constexpr double t25 = 0.1e1 / t24;
    constexpr double t26 = t21 * t25;
    constexpr double t28 = t27 * t27;
    constexpr double t51 = t21 * t21;
    constexpr double t53 = 0.1e1 / t23 / t22;
    constexpr double t54 = t51 * t53;
    constexpr double t117 = t25 * t28;
    constexpr double t252 = t53 * t27;


    const double t3 = rho / 0.2e1 <= dens_tol;
    const double t8 = 0.1e1 <= zeta_tol;
    const double t9 = zeta_tol - 0.1e1;
    const double t11 = piecewise_functor_5( t8, t9, t8, -t9, 0.0 );
    const double t12 = 0.1e1 + t11;
    const double t14 = safe_math::cbrt( zeta_tol );
    const double t16 = safe_math::cbrt( t12 );
    const double t18 = piecewise_functor_3( t12 <= zeta_tol, t14 * zeta_tol, t16 * t12 );
    const double t19 = safe_math::cbrt( rho );
    const double t29 = sigma * t28;
    const double t30 = rho * rho;
    const double t31 = t19 * t19;
    const double t33 = 0.1e1 / t31 / t30;
    const double t34 = t29 * t33;
    const double t35 = t26 * t34;
    const double t37 = tau * t28;
    const double t39 = 0.1e1 / t31 / rho;
    const double t44 = t26 * t37 * t39 / 0.4e1 - 0.9e1 / 0.2e2 - t35 / 0.288e3;
    const double t45 = t44 * t44;
    const double t47 = t44 * t21;
    const double t48 = t47 * t25;
    const double t55 = sigma * sigma;
    const double t56 = t55 * t27;
    const double t57 = t30 * t30;
    const double t58 = t57 * rho;
    const double t60 = 0.1e1 / t19 / t58;
    const double t64 = 0.804e0 + 0.5e1 / 0.972e3 * t35 + 0.146e3 / 0.2025e4 * t45 - 0.73e2 / 0.972e4 * t48 * t34 + 0.45818468001825619316e-3 * t54 * t56 * t60;
    const double t67 = 0.1804e1 - 0.646416e0 / t64;
    const double t71 = piecewise_functor_3( t3, 0.0, -0.3e1 / 0.8e1 * t7 * t18 * t19 * t67 );
    const double t72 = 0.1e1 / t31;
    const double t77 = t4 * t18;
    const double t78 = t64 * t64;
    const double t79 = 0.1e1 / t78;
    const double t80 = t19 * t79;
    const double t81 = t30 * rho;
    const double t83 = 0.1e1 / t31 / t81;
    const double t84 = t29 * t83;
    const double t85 = t26 * t84;
    const double t91 = -0.5e1 / 0.12e2 * t26 * t37 * t33 + t85 / 0.108e3;
    const double t94 = t91 * t21;
    const double t95 = t94 * t25;
    const double t100 = t57 * t30;
    const double t102 = 0.1e1 / t19 / t100;
    const double t106 = -0.1e2 / 0.729e3 * t85 + 0.292e3 / 0.2025e4 * t44 * t91 - 0.73e2 / 0.972e4 * t95 * t34 + 0.73e2 / 0.3645e4 * t48 * t84 - 0.24436516267640330302e-2 * t54 * t56 * t102;
    const double t111 = piecewise_functor_3( t3, 0.0, -t7 * t18 * t72 * t67 / 0.8e1 - 0.16551095363746320496e0 * t77 * t80 * t106 );
    const double t118 = t117 * t33;
    const double t119 = t47 * t118;
    const double t123 = t54 * t27 * t60 * sigma;
    const double t125 = 0.5e1 / 0.972e3 * t26 * t28 * t33 - 0.146e3 / 0.18225e5 * t119 + 0.96852413827153753492e-3 * t123;
    const double t129 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t125 );
    const double t131 = t117 * t39;
    const double t140 = 0.73e2 / 0.2025e4 * t47 * t131 - 0.73e2 / 0.1944e5 * t54 * t27 / t19 / t57 * sigma;
    const double t144 = piecewise_functor_3( t3, 0.0, -0.16551095363746320496e0 * t77 * t80 * t140 );
    const double t151 = t72 * t79;
    const double t156 = 0.1e1 / t78 / t64;
    const double t157 = t19 * t156;
    const double t158 = t106 * t106;
    const double t163 = 0.1e1 / t31 / t57;
    const double t164 = t29 * t163;
    const double t165 = t26 * t164;
    const double t167 = t91 * t91;
    const double t173 = 0.1e2 / 0.9e1 * t26 * t37 * t83 - 0.11e2 / 0.324e3 * t165;
    const double t176 = t173 * t21;
    const double t177 = t176 * t25;
    const double t184 = t57 * t81;
    const double t186 = 0.1e1 / t19 / t184;
    const double t190 = 0.11e3 / 0.2187e4 * t165 + 0.292e3 / 0.2025e4 * t167 + 0.292e3 / 0.2025e4 * t44 * t173 - 0.73e2 / 0.972e4 * t177 * t34 + 0.146e3 / 0.3645e4 * t95 * t84 - 0.803e3 / 0.10935e5 * t48 * t164 + 0.15476460302838875858e-1 * t54 * t56 * t186;
    const double t195 = piecewise_functor_3( t3, 0.0, t7 * t18 * t39 * t67 / 0.12e2 - 0.11034063575830880331e0 * t77 * t151 * t106 + 0.33102190727492640992e0 * t77 * t157 * t158 - 0.16551095363746320496e0 * t77 * t80 * t190 );
    const double t201 = t77 * t19;
    const double t202 = t156 * t125;
    const double t203 = t202 * t106;
    const double t209 = t94 * t118;
    const double t211 = t117 * t83;
    const double t212 = t47 * t211;
    const double t216 = t54 * t27 * t102 * sigma;
    const double t218 = -0.1e2 / 0.729e3 * t26 * t28 * t83 - 0.146e3 / 0.18225e5 * t209 + 0.1168e4 / 0.54675e5 * t212 - 0.51654620707815335196e-2 * t216;
    const double t223 = piecewise_functor_3( t3, 0.0, -0.55170317879154401653e-1 * t77 * t151 * t125 + 0.33102190727492640992e0 * t201 * t203 - 0.16551095363746320496e0 * t77 * t80 * t218 );
    const double t229 = t156 * t140;
    const double t230 = t229 * t106;
    const double t237 = 0.73e2 / 0.2025e4 * t94 * t131 - 0.73e2 / 0.1215e4 * t119 + 0.949e3 / 0.5832e5 * t123;
    const double t242 = piecewise_functor_3( t3, 0.0, -0.55170317879154401653e-1 * t77 * t151 * t140 + 0.33102190727492640992e0 * t201 * t230 - 0.16551095363746320496e0 * t77 * t80 * t237 );
    const double t245 = t125 * t125;
    const double t249 = 0.1e1 / t58;
    const double t253 = t79 * t51 * t252;
    const double t254 = t77 * t249 * t253;
    const double t257 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t77 * t157 * t245 - 0.16950901996748250202e-3 * t254 );
    const double t259 = t229 * t125;
    const double t262 = 0.1e1 / t57;
    const double t264 = t77 * t262 * t253;
    const double t267 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t201 * t259 + 0.66295196793057964127e-3 * t264 );
    const double t269 = t140 * t140;
    const double t273 = 0.1e1 / t81;
    const double t278 = piecewise_functor_3( t3, 0.0, 0.33102190727492640992e0 * t77 * t157 * t269 - 0.29832838556876083857e-2 * t77 * t273 * t253 );


    vrho = 0.2e1 * rho * t111 + 0.2e1 * t71;
    vsigma = 0.2e1 * rho * t129;
    vlapl = 0.e0;
    vtau = 0.2e1 * rho * t144;
    v2rho2 = 0.2e1 * rho * t195 + 0.4e1 * t111;
    v2rhosigma = 0.2e1 * rho * t223 + 0.2e1 * t129;
    v2rholapl = 0.e0;
    v2rhotau = 0.2e1 * rho * t242 + 0.2e1 * t144;
    v2sigma2 = 0.2e1 * rho * t257;
    v2sigmalapl = 0.e0;
    v2sigmatau = 0.2e1 * rho * t267;
    v2lapl2 = 0.e0;
    v2lapltau = 0.e0;
    v2tau2 = 0.2e1 * rho * t278;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t76 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t26 * t27 * t72 );
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t124 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t85 * t27 * t120 );


    eps = t76 + t124;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t76 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t26 * t27 * t72 );
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t124 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t85 * t27 * t120 );
    const double t125 = t7 * t7;
    const double t126 = 0.1e1 / t125;
    const double t127 = t17 * t126;
    const double t129 = piecewise_functor_5( t11, 0.0, t15, 0.0, t8 - t127 );
    const double t132 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t129 );
    const double t137 = t27 * t27;
    const double t138 = 0.1e1 / t137;
    const double t142 = t6 * t26 * t138 * t72 / 0.8e1;
    const double t143 = t3 * t26;
    const double t144 = t69 * t69;
    const double t145 = 0.1e1 / t144;
    const double t146 = t27 * t145;
    const double t147 = t35 * rho_a;
    const double t149 = 0.1e1 / t37 / t147;
    const double t151 = t34 * sigma_aa * t149;
    const double t157 = -0.5e1 / 0.12e2 * t34 * tau_a * t39 + t151 / 0.108e3;
    const double t160 = t157 * t29;
    const double t163 = t53 * t149;
    const double t166 = t62 * t35;
    const double t168 = 0.1e1 / t36 / t166;
    const double t172 = -0.1e2 / 0.729e3 * t151 + 0.292e3 / 0.2025e4 * t49 * t157 - 0.73e2 / 0.972e4 * t160 * t54 + 0.73e2 / 0.3645e4 * t52 * t163 - 0.12218258133820165151e-2 * t60 * t61 * t168;
    const double t173 = t146 * t172;
    const double t177 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t132 * t27 * t72 - t142 - 0.16551095363746320496e0 * t143 * t173 );
    const double t178 = t78 * t126;
    const double t180 = piecewise_functor_5( t15, 0.0, t11, 0.0, -t8 - t178 );
    const double t183 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t180 );
    const double t191 = t6 * t85 * t138 * t120 / 0.8e1;
    const double t193 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t183 * t27 * t120 - t191 );
    const double t197 = piecewise_functor_5( t11, 0.0, t15, 0.0, -t8 - t127 );
    const double t200 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t197 );
    const double t206 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t200 * t27 * t72 - t142 );
    const double t208 = piecewise_functor_5( t15, 0.0, t11, 0.0, t8 - t178 );
    const double t211 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t208 );
    const double t216 = t3 * t85;
    const double t217 = t117 * t117;
    const double t218 = 0.1e1 / t217;
    const double t219 = t27 * t218;
    const double t220 = t87 * rho_b;
    const double t222 = 0.1e1 / t89 / t220;
    const double t224 = t34 * sigma_bb * t222;
    const double t230 = -0.5e1 / 0.12e2 * t34 * tau_b * t91 + t224 / 0.108e3;
    const double t233 = t230 * t29;
    const double t236 = t105 * t222;
    const double t239 = t110 * t87;
    const double t241 = 0.1e1 / t88 / t239;
    const double t245 = -0.1e2 / 0.729e3 * t224 + 0.292e3 / 0.2025e4 * t101 * t230 - 0.73e2 / 0.972e4 * t233 * t106 + 0.73e2 / 0.3645e4 * t104 * t236 - 0.12218258133820165151e-2 * t60 * t109 * t241;
    const double t246 = t219 * t245;
    const double t250 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t211 * t27 * t120 - t191 - 0.16551095363746320496e0 * t216 * t246 );
    const double t255 = t33 * t39;
    const double t256 = t52 * t255;
    const double t259 = t60 * t65 * sigma_aa;
    const double t261 = 0.5e1 / 0.972e3 * t34 * t39 - 0.146e3 / 0.18225e5 * t256 + 0.48426206913576876746e-3 * t259;
    const double t262 = t146 * t261;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t262 );
    const double t268 = t33 * t91;
    const double t269 = t104 * t268;
    const double t272 = t60 * t113 * sigma_bb;
    const double t274 = 0.5e1 / 0.972e3 * t34 * t91 - 0.146e3 / 0.18225e5 * t269 + 0.48426206913576876746e-3 * t272;
    const double t275 = t219 * t274;
    const double t278 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t275 );
    const double t279 = t33 * t44;
    const double t283 = 0.1e1 / t36 / t62;
    const double t287 = 0.73e2 / 0.2025e4 * t52 * t279 - 0.73e2 / 0.3888e5 * t60 * t283 * sigma_aa;
    const double t288 = t146 * t287;
    const double t291 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t288 );
    const double t292 = t33 * t96;
    const double t296 = 0.1e1 / t88 / t110;
    const double t300 = 0.73e2 / 0.2025e4 * t104 * t292 - 0.73e2 / 0.3888e5 * t60 * t296 * sigma_bb;
    const double t301 = t219 * t300;
    const double t304 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t301 );


    eps = t76 + t124;
    vrho_a = t76 + t124 + t7 * ( t177 + t193 );
    vrho_b = t76 + t124 + t7 * ( t206 + t250 );
    vsigma_aa = t7 * t265;
    vsigma_ab = 0.e0;
    vsigma_bb = t7 * t278;
    vlapl_a = 0.e0;
    vlapl_b = 0.e0;
    vtau_a = t7 * t291;
    vtau_b = t7 * t304;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b, double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb, double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b, double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b, double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb, double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b, double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t125 = t7 * t7;
    const double t126 = 0.1e1 / t125;
    const double t127 = t17 * t126;
    const double t129 = piecewise_functor_5( t11, 0.0, t15, 0.0, t8 - t127 );
    const double t132 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t129 );
    const double t137 = t27 * t27;
    const double t138 = 0.1e1 / t137;
    const double t142 = t6 * t26 * t138 * t72 / 0.8e1;
    const double t143 = t3 * t26;
    const double t144 = t69 * t69;
    const double t145 = 0.1e1 / t144;
    const double t146 = t27 * t145;
    const double t147 = t35 * rho_a;
    const double t149 = 0.1e1 / t37 / t147;
    const double t151 = t34 * sigma_aa * t149;
    const double t157 = -0.5e1 / 0.12e2 * t34 * tau_a * t39 + t151 / 0.108e3;
    const double t160 = t157 * t29;
    const double t163 = t53 * t149;
    const double t166 = t62 * t35;
    const double t168 = 0.1e1 / t36 / t166;
    const double t172 = -0.1e2 / 0.729e3 * t151 + 0.292e3 / 0.2025e4 * t49 * t157 - 0.73e2 / 0.972e4 * t160 * t54 + 0.73e2 / 0.3645e4 * t52 * t163 - 0.12218258133820165151e-2 * t60 * t61 * t168;
    const double t173 = t146 * t172;
    const double t177 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t132 * t27 * t72 - t142 - 0.16551095363746320496e0 * t143 * t173 );
    const double t178 = t78 * t126;
    const double t180 = piecewise_functor_5( t15, 0.0, t11, 0.0, -t8 - t178 );
    const double t183 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t180 );
    const double t191 = t6 * t85 * t138 * t120 / 0.8e1;
    const double t193 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t183 * t27 * t120 - t191 );
    const double t197 = piecewise_functor_5( t11, 0.0, t15, 0.0, -t8 - t127 );
    const double t200 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t197 );
    const double t206 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t200 * t27 * t72 - t142 );
    const double t208 = piecewise_functor_5( t15, 0.0, t11, 0.0, t8 - t178 );
    const double t211 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t208 );
    const double t216 = t3 * t85;
    const double t217 = t117 * t117;
    const double t218 = 0.1e1 / t217;
    const double t219 = t27 * t218;
    const double t220 = t87 * rho_b;
    const double t222 = 0.1e1 / t89 / t220;
    const double t224 = t34 * sigma_bb * t222;
    const double t230 = -0.5e1 / 0.12e2 * t34 * tau_b * t91 + t224 / 0.108e3;
    const double t233 = t230 * t29;
    const double t236 = t105 * t222;
    const double t239 = t110 * t87;
    const double t241 = 0.1e1 / t88 / t239;
    const double t245 = -0.1e2 / 0.729e3 * t224 + 0.292e3 / 0.2025e4 * t101 * t230 - 0.73e2 / 0.972e4 * t233 * t106 + 0.73e2 / 0.3645e4 * t104 * t236 - 0.12218258133820165151e-2 * t60 * t109 * t241;
    const double t246 = t219 * t245;
    const double t250 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t211 * t27 * t120 - t191 - 0.16551095363746320496e0 * t216 * t246 );
    const double t255 = t33 * t39;
    const double t256 = t52 * t255;
    const double t259 = t60 * t65 * sigma_aa;
    const double t261 = 0.5e1 / 0.972e3 * t34 * t39 - 0.146e3 / 0.18225e5 * t256 + 0.48426206913576876746e-3 * t259;
    const double t262 = t146 * t261;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t262 );
    const double t268 = t33 * t91;
    const double t269 = t104 * t268;
    const double t272 = t60 * t113 * sigma_bb;
    const double t274 = 0.5e1 / 0.972e3 * t34 * t91 - 0.146e3 / 0.18225e5 * t269 + 0.48426206913576876746e-3 * t272;
    const double t275 = t219 * t274;
    const double t278 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t275 );
    const double t279 = t33 * t44;
    const double t283 = 0.1e1 / t36 / t62;
    const double t287 = 0.73e2 / 0.2025e4 * t52 * t279 - 0.73e2 / 0.3888e5 * t60 * t283 * sigma_aa;
    const double t288 = t146 * t287;
    const double t291 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t288 );
    const double t292 = t33 * t96;
    const double t296 = 0.1e1 / t88 / t110;
    const double t300 = 0.73e2 / 0.2025e4 * t104 * t292 - 0.73e2 / 0.3888e5 * t60 * t296 * sigma_bb;
    const double t301 = t219 * t300;
    const double t304 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t301 );
    const double t307 = t24 * t24;
    const double t308 = 0.1e1 / t307;
    const double t309 = t129 * t129;
    const double t312 = t125 * t7;
    const double t313 = 0.1e1 / t312;
    const double t314 = t17 * t313;
    const double t317 = piecewise_functor_5( t11, 0.0, t15, 0.0, -0.2e1 * t126 + 0.2e1 * t314 );
    const double t321 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.9e1 * t308 * t309 + 0.4e1 / 0.3e1 * t24 * t317 );
    const double t328 = t6 * t132 * t138 * t72;
    const double t330 = t3 * t132;
    const double t334 = 0.1e1 / t137 / t7;
    const double t338 = t6 * t26 * t334 * t72 / 0.12e2;
    const double t339 = t138 * t145;
    const double t340 = t339 * t172;
    const double t341 = t143 * t340;
    const double t344 = 0.1e1 / t144 / t69;
    const double t345 = t27 * t344;
    const double t346 = t172 * t172;
    const double t347 = t345 * t346;
    const double t351 = 0.1e1 / t37 / t62;
    const double t353 = t34 * sigma_aa * t351;
    const double t355 = t157 * t157;
    const double t361 = 0.1e2 / 0.9e1 * t34 * tau_a * t149 - 0.11e2 / 0.324e3 * t353;
    const double t364 = t361 * t29;
    const double t369 = t53 * t351;
    const double t372 = t62 * t147;
    const double t374 = 0.1e1 / t36 / t372;
    const double t378 = 0.11e3 / 0.2187e4 * t353 + 0.292e3 / 0.2025e4 * t355 + 0.292e3 / 0.2025e4 * t49 * t361 - 0.73e2 / 0.972e4 * t364 * t54 + 0.146e3 / 0.3645e4 * t160 * t163 - 0.803e3 / 0.10935e5 * t52 * t369 + 0.7738230151419437929e-2 * t60 * t61 * t374;
    const double t379 = t146 * t378;
    const double t383 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t321 * t27 * t72 - t328 / 0.4e1 - 0.33102190727492640992e0 * t330 * t173 + t338 - 0.11034063575830880331e0 * t341 + 0.33102190727492640992e0 * t143 * t347 - 0.16551095363746320496e0 * t143 * t379 );
    const double t384 = t83 * t83;
    const double t385 = 0.1e1 / t384;
    const double t386 = t180 * t180;
    const double t389 = t78 * t313;
    const double t392 = piecewise_functor_5( t15, 0.0, t11, 0.0, 0.2e1 * t126 + 0.2e1 * t389 );
    const double t396 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.9e1 * t385 * t386 + 0.4e1 / 0.3e1 * t83 * t392 );
    const double t403 = t6 * t183 * t138 * t120;
    const double t408 = t6 * t85 * t334 * t120 / 0.12e2;
    const double t410 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t396 * t27 * t120 - t403 / 0.4e1 + t408 );
    const double t428 = t6 * t200 * t138 * t72;
    const double t452 = t6 * t211 * t138 * t120;
    const double t458 = t138 * t218;
    const double t459 = t458 * t245;
    const double t460 = t216 * t459;
    const double t468 = t197 * t197;
    const double t473 = piecewise_functor_5( t11, 0.0, t15, 0.0, 0.2e1 * t126 + 0.2e1 * t314 );
    const double t477 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.9e1 * t308 * t468 + 0.4e1 / 0.3e1 * t24 * t473 );
    const double t484 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t477 * t27 * t72 - t428 / 0.4e1 + t338 );
    const double t485 = t208 * t208;
    const double t490 = piecewise_functor_5( t15, 0.0, t11, 0.0, -0.2e1 * t126 + 0.2e1 * t389 );
    const double t494 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.9e1 * t385 * t485 + 0.4e1 / 0.3e1 * t83 * t490 );
    const double t500 = t3 * t211;
    const double t505 = 0.1e1 / t217 / t117;
    const double t506 = t27 * t505;
    const double t507 = t245 * t245;
    const double t508 = t506 * t507;
    const double t512 = 0.1e1 / t89 / t110;
    const double t514 = t34 * sigma_bb * t512;
    const double t516 = t230 * t230;
    const double t522 = 0.1e2 / 0.9e1 * t34 * tau_b * t222 - 0.11e2 / 0.324e3 * t514;
    const double t525 = t522 * t29;
    const double t530 = t105 * t512;
    const double t533 = t110 * t220;
    const double t535 = 0.1e1 / t88 / t533;
    const double t539 = 0.11e3 / 0.2187e4 * t514 + 0.292e3 / 0.2025e4 * t516 + 0.292e3 / 0.2025e4 * t101 * t522 - 0.73e2 / 0.972e4 * t525 * t106 + 0.146e3 / 0.3645e4 * t233 * t236 - 0.803e3 / 0.10935e5 * t104 * t530 + 0.7738230151419437929e-2 * t60 * t109 * t535;
    const double t540 = t219 * t539;
    const double t544 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t494 * t27 * t120 - t452 / 0.4e1 - 0.33102190727492640992e0 * t500 * t246 + t408 - 0.11034063575830880331e0 * t460 + 0.33102190727492640992e0 * t216 * t508 - 0.16551095363746320496e0 * t216 * t540 );
    const double t549 = t339 * t261;
    const double t551 = 0.55170317879154401653e-1 * t143 * t549;
    const double t552 = t143 * t27;
    const double t553 = t344 * t261;
    const double t554 = t553 * t172;
    const double t559 = t160 * t255;
    const double t561 = t33 * t149;
    const double t562 = t52 * t561;
    const double t565 = t60 * t168 * sigma_aa;
    const double t567 = -0.1e2 / 0.729e3 * t34 * t149 - 0.146e3 / 0.18225e5 * t559 + 0.1168e4 / 0.54675e5 * t562 - 0.25827310353907667598e-2 * t565;
    const double t568 = t146 * t567;
    const double t572 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t330 * t262 - t551 + 0.33102190727492640992e0 * t552 * t554 - 0.16551095363746320496e0 * t143 * t568 );
    const double t576 = t458 * t274;
    const double t578 = 0.55170317879154401653e-1 * t216 * t576;
    const double t589 = t216 * t27;
    const double t590 = t505 * t274;
    const double t591 = t590 * t245;
    const double t596 = t233 * t268;
    const double t598 = t33 * t222;
    const double t599 = t104 * t598;
    const double t602 = t60 * t241 * sigma_bb;
    const double t604 = -0.1e2 / 0.729e3 * t34 * t222 - 0.146e3 / 0.18225e5 * t596 + 0.1168e4 / 0.54675e5 * t599 - 0.25827310353907667598e-2 * t602;
    const double t605 = t219 * t604;
    const double t609 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t500 * t275 - t578 + 0.33102190727492640992e0 * t589 * t591 - 0.16551095363746320496e0 * t216 * t605 );
    const double t613 = t339 * t287;
    const double t615 = 0.55170317879154401653e-1 * t143 * t613;
    const double t616 = t344 * t287;
    const double t617 = t616 * t172;
    const double t624 = 0.73e2 / 0.2025e4 * t160 * t279 - 0.73e2 / 0.1215e4 * t256 + 0.949e3 / 0.11664e6 * t259;
    const double t625 = t146 * t624;
    const double t629 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t330 * t288 - t615 + 0.33102190727492640992e0 * t552 * t617 - 0.16551095363746320496e0 * t143 * t625 );
    const double t633 = t458 * t300;
    const double t635 = 0.55170317879154401653e-1 * t216 * t633;
    const double t646 = t505 * t300;
    const double t647 = t646 * t245;
    const double t654 = 0.73e2 / 0.2025e4 * t233 * t292 - 0.73e2 / 0.1215e4 * t269 + 0.949e3 / 0.11664e6 * t272;
    const double t655 = t219 * t654;
    const double t659 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t500 * t301 - t635 + 0.33102190727492640992e0 * t589 * t647 - 0.16551095363746320496e0 * t216 * t655 );
    const double t661 = t261 * t261;
    const double t662 = t345 * t661;
    const double t665 = t145 * t57;
    const double t666 = t59 * t65;
    const double t667 = t665 * t666;
    const double t668 = t552 * t667;
    const double t671 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t143 * t662 - 0.84754509983741251008e-4 * t668 );
    const double t672 = t274 * t274;
    const double t673 = t506 * t672;
    const double t676 = t218 * t57;
    const double t677 = t59 * t113;
    const double t678 = t676 * t677;
    const double t679 = t589 * t678;
    const double t682 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t216 * t673 - 0.84754509983741251008e-4 * t679 );
    const double t683 = t616 * t261;
    const double t686 = t59 * t283;
    const double t687 = t665 * t686;
    const double t688 = t552 * t687;
    const double t691 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t552 * t683 + 0.33147598396528982063e-3 * t688 );
    const double t692 = t646 * t274;
    const double t695 = t59 * t296;
    const double t696 = t676 * t695;
    const double t697 = t589 * t696;
    const double t700 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t589 * t692 + 0.33147598396528982063e-3 * t697 );
    const double t701 = t287 * t287;
    const double t702 = t345 * t701;
    const double t706 = 0.1e1 / t36 / t147;
    const double t707 = t59 * t706;
    const double t708 = t665 * t707;
    const double t712 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t143 * t702 - 0.14916419278438041928e-2 * t552 * t708 );
    const double t713 = t300 * t300;
    const double t714 = t506 * t713;
    const double t718 = 0.1e1 / t88 / t220;
    const double t719 = t59 * t718;
    const double t720 = t676 * t719;
    const double t724 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t216 * t714 - 0.14916419278438041928e-2 * t589 * t720 );


    v2rho2_aa = 0.2e1 * t177 + 0.2e1 * t193 + t7 * ( t383 + t410 );
    v2rho2_bb = 0.2e1 * t206 + 0.2e1 * t250 + t7 * ( t484 + t544 );
    v2rhosigma_a_aa = t7 * t572 + t265;
    v2rhosigma_b_bb = t7 * t609 + t278;
    v2rholapl_a_a = 0.e0;
    v2rholapl_b_b = 0.e0;
    v2rhotau_a_a = t7 * t629 + t291;
    v2rhotau_b_b = t7 * t659 + t304;
    v2sigma2_aa_aa = t7 * t671;
    v2sigma2_bb_bb = t7 * t682;
    v2sigmalapl_aa_a = 0.e0;
    v2sigmalapl_bb_b = 0.e0;
    v2sigmatau_aa_a = t7 * t691;
    v2sigmatau_bb_b = t7 * t700;
    v2lapl2_aa = 0.e0;
    v2lapl2_bb = 0.e0;
    v2lapltau_a_a = 0.e0;
    v2lapltau_b_b = 0.e0;
    v2tau2_aa = t7 * t712;
    v2tau2_bb = t7 * t724;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2rholapl_a_b = 0.0;
    v2rholapl_b_a = 0.0;
    v2rhotau_a_b = 0.0;
    v2rhotau_b_a = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;
    v2sigmalapl_aa_b = 0.0;
    v2sigmalapl_ab_a = 0.0;
    v2sigmalapl_ab_b = 0.0;
    v2sigmalapl_bb_a = 0.0;
    v2sigmatau_aa_b = 0.0;
    v2sigmatau_ab_a = 0.0;
    v2sigmatau_ab_b = 0.0;
    v2sigmatau_bb_a = 0.0;
    v2lapl2_ab = 0.0;
    v2lapltau_a_b = 0.0;
    v2lapltau_b_a = 0.0;
    v2tau2_ab = 0.0;



  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_vxc_fxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double lapl_a, double lapl_b, double tau_a, double tau_b, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb, double& vlapl_a, double& vlapl_b, double& vtau_a, double& vtau_b, double& v2rho2_aa, double& v2rho2_ab, double& v2rho2_bb, double& v2rhosigma_a_aa, double& v2rhosigma_a_ab, double& v2rhosigma_a_bb, double& v2rhosigma_b_aa, double& v2rhosigma_b_ab, double& v2rhosigma_b_bb, double& v2rholapl_a_a, double& v2rholapl_a_b, double& v2rholapl_b_a, double& v2rholapl_b_b, double& v2rhotau_a_a, double& v2rhotau_a_b, double& v2rhotau_b_a, double& v2rhotau_b_b, double& v2sigma2_aa_aa, double& v2sigma2_aa_ab, double& v2sigma2_aa_bb, double& v2sigma2_ab_ab, double& v2sigma2_ab_bb, double& v2sigma2_bb_bb, double& v2sigmalapl_aa_a, double& v2sigmalapl_aa_b, double& v2sigmalapl_ab_a, double& v2sigmalapl_ab_b, double& v2sigmalapl_bb_a, double& v2sigmalapl_bb_b, double& v2sigmatau_aa_a, double& v2sigmatau_aa_b, double& v2sigmatau_ab_a, double& v2sigmatau_ab_b, double& v2sigmatau_bb_a, double& v2sigmatau_bb_b, double& v2lapl2_aa, double& v2lapl2_ab, double& v2lapl2_bb, double& v2lapltau_a_a, double& v2lapltau_a_b, double& v2lapltau_b_a, double& v2lapltau_b_b, double& v2tau2_aa, double& v2tau2_ab, double& v2tau2_bb ) {

    (void)(sigma_ab);
    (void)(lapl_a);
    (void)(lapl_b);
    constexpr double t3 = constants::m_cbrt_3;
    constexpr double t4 = constants::m_cbrt_pi;
    constexpr double t29 = constants::m_cbrt_6;
    constexpr double t30 = constants::m_pi_sq;
    constexpr double t31 = constants::m_cbrt_pi_sq;
    constexpr double t6 = t3 / t4;
    constexpr double t32 = t31 * t31;
    constexpr double t33 = 0.1e1 / t32;
    constexpr double t34 = t29 * t33;
    constexpr double t57 = t29 * t29;
    constexpr double t59 = 0.1e1 / t31 / t30;
    constexpr double t60 = t57 * t59;


    const double t2 = rho_a <= dens_tol;
    const double t7 = rho_a + rho_b;
    const double t8 = 0.1e1 / t7;
    const double t11 = 0.2e1 * rho_a * t8 <= zeta_tol;
    const double t12 = zeta_tol - 0.1e1;
    const double t15 = 0.2e1 * rho_b * t8 <= zeta_tol;
    const double t16 = -t12;
    const double t17 = rho_a - rho_b;
    const double t19 = piecewise_functor_5( t11, t12, t15, t16, t17 * t8 );
    const double t20 = 0.1e1 + t19;
    const double t21 = t20 <= zeta_tol;
    const double t22 = safe_math::cbrt( zeta_tol );
    const double t23 = t22 * zeta_tol;
    const double t24 = safe_math::cbrt( t20 );
    const double t26 = piecewise_functor_3( t21, t23, t24 * t20 );
    const double t27 = safe_math::cbrt( t7 );
    const double t35 = rho_a * rho_a;
    const double t36 = safe_math::cbrt( rho_a );
    const double t37 = t36 * t36;
    const double t39 = 0.1e1 / t37 / t35;
    const double t41 = t34 * sigma_aa * t39;
    const double t44 = 0.1e1 / t37 / rho_a;
    const double t49 = t34 * tau_a * t44 / 0.4e1 - 0.9e1 / 0.2e2 - t41 / 0.288e3;
    const double t50 = t49 * t49;
    const double t52 = t49 * t29;
    const double t53 = t33 * sigma_aa;
    const double t54 = t53 * t39;
    const double t61 = sigma_aa * sigma_aa;
    const double t62 = t35 * t35;
    const double t63 = t62 * rho_a;
    const double t65 = 0.1e1 / t36 / t63;
    const double t69 = 0.804e0 + 0.5e1 / 0.972e3 * t41 + 0.146e3 / 0.2025e4 * t50 - 0.73e2 / 0.972e4 * t52 * t54 + 0.22909234000912809658e-3 * t60 * t61 * t65;
    const double t72 = 0.1804e1 - 0.646416e0 / t69;
    const double t76 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t26 * t27 * t72 );
    const double t77 = rho_b <= dens_tol;
    const double t78 = -t17;
    const double t80 = piecewise_functor_5( t15, t12, t11, t16, t78 * t8 );
    const double t81 = 0.1e1 + t80;
    const double t82 = t81 <= zeta_tol;
    const double t83 = safe_math::cbrt( t81 );
    const double t85 = piecewise_functor_3( t82, t23, t83 * t81 );
    const double t87 = rho_b * rho_b;
    const double t88 = safe_math::cbrt( rho_b );
    const double t89 = t88 * t88;
    const double t91 = 0.1e1 / t89 / t87;
    const double t93 = t34 * sigma_bb * t91;
    const double t96 = 0.1e1 / t89 / rho_b;
    const double t101 = t34 * tau_b * t96 / 0.4e1 - 0.9e1 / 0.2e2 - t93 / 0.288e3;
    const double t102 = t101 * t101;
    const double t104 = t101 * t29;
    const double t105 = t33 * sigma_bb;
    const double t106 = t105 * t91;
    const double t109 = sigma_bb * sigma_bb;
    const double t110 = t87 * t87;
    const double t111 = t110 * rho_b;
    const double t113 = 0.1e1 / t88 / t111;
    const double t117 = 0.804e0 + 0.5e1 / 0.972e3 * t93 + 0.146e3 / 0.2025e4 * t102 - 0.73e2 / 0.972e4 * t104 * t106 + 0.22909234000912809658e-3 * t60 * t109 * t113;
    const double t120 = 0.1804e1 - 0.646416e0 / t117;
    const double t124 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t85 * t27 * t120 );
    const double t125 = t7 * t7;
    const double t126 = 0.1e1 / t125;
    const double t127 = t17 * t126;
    const double t129 = piecewise_functor_5( t11, 0.0, t15, 0.0, t8 - t127 );
    const double t132 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t129 );
    const double t137 = t27 * t27;
    const double t138 = 0.1e1 / t137;
    const double t142 = t6 * t26 * t138 * t72 / 0.8e1;
    const double t143 = t3 * t26;
    const double t144 = t69 * t69;
    const double t145 = 0.1e1 / t144;
    const double t146 = t27 * t145;
    const double t147 = t35 * rho_a;
    const double t149 = 0.1e1 / t37 / t147;
    const double t151 = t34 * sigma_aa * t149;
    const double t157 = -0.5e1 / 0.12e2 * t34 * tau_a * t39 + t151 / 0.108e3;
    const double t160 = t157 * t29;
    const double t163 = t53 * t149;
    const double t166 = t62 * t35;
    const double t168 = 0.1e1 / t36 / t166;
    const double t172 = -0.1e2 / 0.729e3 * t151 + 0.292e3 / 0.2025e4 * t49 * t157 - 0.73e2 / 0.972e4 * t160 * t54 + 0.73e2 / 0.3645e4 * t52 * t163 - 0.12218258133820165151e-2 * t60 * t61 * t168;
    const double t173 = t146 * t172;
    const double t177 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t132 * t27 * t72 - t142 - 0.16551095363746320496e0 * t143 * t173 );
    const double t178 = t78 * t126;
    const double t180 = piecewise_functor_5( t15, 0.0, t11, 0.0, -t8 - t178 );
    const double t183 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t180 );
    const double t191 = t6 * t85 * t138 * t120 / 0.8e1;
    const double t193 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t183 * t27 * t120 - t191 );
    const double t197 = piecewise_functor_5( t11, 0.0, t15, 0.0, -t8 - t127 );
    const double t200 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.3e1 * t24 * t197 );
    const double t206 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t200 * t27 * t72 - t142 );
    const double t208 = piecewise_functor_5( t15, 0.0, t11, 0.0, t8 - t178 );
    const double t211 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.3e1 * t83 * t208 );
    const double t216 = t3 * t85;
    const double t217 = t117 * t117;
    const double t218 = 0.1e1 / t217;
    const double t219 = t27 * t218;
    const double t220 = t87 * rho_b;
    const double t222 = 0.1e1 / t89 / t220;
    const double t224 = t34 * sigma_bb * t222;
    const double t230 = -0.5e1 / 0.12e2 * t34 * tau_b * t91 + t224 / 0.108e3;
    const double t233 = t230 * t29;
    const double t236 = t105 * t222;
    const double t239 = t110 * t87;
    const double t241 = 0.1e1 / t88 / t239;
    const double t245 = -0.1e2 / 0.729e3 * t224 + 0.292e3 / 0.2025e4 * t101 * t230 - 0.73e2 / 0.972e4 * t233 * t106 + 0.73e2 / 0.3645e4 * t104 * t236 - 0.12218258133820165151e-2 * t60 * t109 * t241;
    const double t246 = t219 * t245;
    const double t250 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t211 * t27 * t120 - t191 - 0.16551095363746320496e0 * t216 * t246 );
    const double t255 = t33 * t39;
    const double t256 = t52 * t255;
    const double t259 = t60 * t65 * sigma_aa;
    const double t261 = 0.5e1 / 0.972e3 * t34 * t39 - 0.146e3 / 0.18225e5 * t256 + 0.48426206913576876746e-3 * t259;
    const double t262 = t146 * t261;
    const double t265 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t262 );
    const double t268 = t33 * t91;
    const double t269 = t104 * t268;
    const double t272 = t60 * t113 * sigma_bb;
    const double t274 = 0.5e1 / 0.972e3 * t34 * t91 - 0.146e3 / 0.18225e5 * t269 + 0.48426206913576876746e-3 * t272;
    const double t275 = t219 * t274;
    const double t278 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t275 );
    const double t279 = t33 * t44;
    const double t283 = 0.1e1 / t36 / t62;
    const double t287 = 0.73e2 / 0.2025e4 * t52 * t279 - 0.73e2 / 0.3888e5 * t60 * t283 * sigma_aa;
    const double t288 = t146 * t287;
    const double t291 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t143 * t288 );
    const double t292 = t33 * t96;
    const double t296 = 0.1e1 / t88 / t110;
    const double t300 = 0.73e2 / 0.2025e4 * t104 * t292 - 0.73e2 / 0.3888e5 * t60 * t296 * sigma_bb;
    const double t301 = t219 * t300;
    const double t304 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t216 * t301 );
    const double t307 = t24 * t24;
    const double t308 = 0.1e1 / t307;
    const double t309 = t129 * t129;
    const double t312 = t125 * t7;
    const double t313 = 0.1e1 / t312;
    const double t314 = t17 * t313;
    const double t317 = piecewise_functor_5( t11, 0.0, t15, 0.0, -0.2e1 * t126 + 0.2e1 * t314 );
    const double t321 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.9e1 * t308 * t309 + 0.4e1 / 0.3e1 * t24 * t317 );
    const double t328 = t6 * t132 * t138 * t72;
    const double t330 = t3 * t132;
    const double t334 = 0.1e1 / t137 / t7;
    const double t338 = t6 * t26 * t334 * t72 / 0.12e2;
    const double t339 = t138 * t145;
    const double t340 = t339 * t172;
    const double t341 = t143 * t340;
    const double t344 = 0.1e1 / t144 / t69;
    const double t345 = t27 * t344;
    const double t346 = t172 * t172;
    const double t347 = t345 * t346;
    const double t351 = 0.1e1 / t37 / t62;
    const double t353 = t34 * sigma_aa * t351;
    const double t355 = t157 * t157;
    const double t361 = 0.1e2 / 0.9e1 * t34 * tau_a * t149 - 0.11e2 / 0.324e3 * t353;
    const double t364 = t361 * t29;
    const double t369 = t53 * t351;
    const double t372 = t62 * t147;
    const double t374 = 0.1e1 / t36 / t372;
    const double t378 = 0.11e3 / 0.2187e4 * t353 + 0.292e3 / 0.2025e4 * t355 + 0.292e3 / 0.2025e4 * t49 * t361 - 0.73e2 / 0.972e4 * t364 * t54 + 0.146e3 / 0.3645e4 * t160 * t163 - 0.803e3 / 0.10935e5 * t52 * t369 + 0.7738230151419437929e-2 * t60 * t61 * t374;
    const double t379 = t146 * t378;
    const double t383 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t321 * t27 * t72 - t328 / 0.4e1 - 0.33102190727492640992e0 * t330 * t173 + t338 - 0.11034063575830880331e0 * t341 + 0.33102190727492640992e0 * t143 * t347 - 0.16551095363746320496e0 * t143 * t379 );
    const double t384 = t83 * t83;
    const double t385 = 0.1e1 / t384;
    const double t386 = t180 * t180;
    const double t389 = t78 * t313;
    const double t392 = piecewise_functor_5( t15, 0.0, t11, 0.0, 0.2e1 * t126 + 0.2e1 * t389 );
    const double t396 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.9e1 * t385 * t386 + 0.4e1 / 0.3e1 * t83 * t392 );
    const double t403 = t6 * t183 * t138 * t120;
    const double t408 = t6 * t85 * t334 * t120 / 0.12e2;
    const double t410 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t396 * t27 * t120 - t403 / 0.4e1 + t408 );
    const double t428 = t6 * t200 * t138 * t72;
    const double t452 = t6 * t211 * t138 * t120;
    const double t458 = t138 * t218;
    const double t459 = t458 * t245;
    const double t460 = t216 * t459;
    const double t468 = t197 * t197;
    const double t473 = piecewise_functor_5( t11, 0.0, t15, 0.0, 0.2e1 * t126 + 0.2e1 * t314 );
    const double t477 = piecewise_functor_3( t21, 0.0, 0.4e1 / 0.9e1 * t308 * t468 + 0.4e1 / 0.3e1 * t24 * t473 );
    const double t484 = piecewise_functor_3( t2, 0.0, -0.3e1 / 0.8e1 * t6 * t477 * t27 * t72 - t428 / 0.4e1 + t338 );
    const double t485 = t208 * t208;
    const double t490 = piecewise_functor_5( t15, 0.0, t11, 0.0, -0.2e1 * t126 + 0.2e1 * t389 );
    const double t494 = piecewise_functor_3( t82, 0.0, 0.4e1 / 0.9e1 * t385 * t485 + 0.4e1 / 0.3e1 * t83 * t490 );
    const double t500 = t3 * t211;
    const double t505 = 0.1e1 / t217 / t117;
    const double t506 = t27 * t505;
    const double t507 = t245 * t245;
    const double t508 = t506 * t507;
    const double t512 = 0.1e1 / t89 / t110;
    const double t514 = t34 * sigma_bb * t512;
    const double t516 = t230 * t230;
    const double t522 = 0.1e2 / 0.9e1 * t34 * tau_b * t222 - 0.11e2 / 0.324e3 * t514;
    const double t525 = t522 * t29;
    const double t530 = t105 * t512;
    const double t533 = t110 * t220;
    const double t535 = 0.1e1 / t88 / t533;
    const double t539 = 0.11e3 / 0.2187e4 * t514 + 0.292e3 / 0.2025e4 * t516 + 0.292e3 / 0.2025e4 * t101 * t522 - 0.73e2 / 0.972e4 * t525 * t106 + 0.146e3 / 0.3645e4 * t233 * t236 - 0.803e3 / 0.10935e5 * t104 * t530 + 0.7738230151419437929e-2 * t60 * t109 * t535;
    const double t540 = t219 * t539;
    const double t544 = piecewise_functor_3( t77, 0.0, -0.3e1 / 0.8e1 * t6 * t494 * t27 * t120 - t452 / 0.4e1 - 0.33102190727492640992e0 * t500 * t246 + t408 - 0.11034063575830880331e0 * t460 + 0.33102190727492640992e0 * t216 * t508 - 0.16551095363746320496e0 * t216 * t540 );
    const double t549 = t339 * t261;
    const double t551 = 0.55170317879154401653e-1 * t143 * t549;
    const double t552 = t143 * t27;
    const double t553 = t344 * t261;
    const double t554 = t553 * t172;
    const double t559 = t160 * t255;
    const double t561 = t33 * t149;
    const double t562 = t52 * t561;
    const double t565 = t60 * t168 * sigma_aa;
    const double t567 = -0.1e2 / 0.729e3 * t34 * t149 - 0.146e3 / 0.18225e5 * t559 + 0.1168e4 / 0.54675e5 * t562 - 0.25827310353907667598e-2 * t565;
    const double t568 = t146 * t567;
    const double t572 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t330 * t262 - t551 + 0.33102190727492640992e0 * t552 * t554 - 0.16551095363746320496e0 * t143 * t568 );
    const double t576 = t458 * t274;
    const double t578 = 0.55170317879154401653e-1 * t216 * t576;
    const double t589 = t216 * t27;
    const double t590 = t505 * t274;
    const double t591 = t590 * t245;
    const double t596 = t233 * t268;
    const double t598 = t33 * t222;
    const double t599 = t104 * t598;
    const double t602 = t60 * t241 * sigma_bb;
    const double t604 = -0.1e2 / 0.729e3 * t34 * t222 - 0.146e3 / 0.18225e5 * t596 + 0.1168e4 / 0.54675e5 * t599 - 0.25827310353907667598e-2 * t602;
    const double t605 = t219 * t604;
    const double t609 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t500 * t275 - t578 + 0.33102190727492640992e0 * t589 * t591 - 0.16551095363746320496e0 * t216 * t605 );
    const double t613 = t339 * t287;
    const double t615 = 0.55170317879154401653e-1 * t143 * t613;
    const double t616 = t344 * t287;
    const double t617 = t616 * t172;
    const double t624 = 0.73e2 / 0.2025e4 * t160 * t279 - 0.73e2 / 0.1215e4 * t256 + 0.949e3 / 0.11664e6 * t259;
    const double t625 = t146 * t624;
    const double t629 = piecewise_functor_3( t2, 0.0, -0.16551095363746320496e0 * t330 * t288 - t615 + 0.33102190727492640992e0 * t552 * t617 - 0.16551095363746320496e0 * t143 * t625 );
    const double t633 = t458 * t300;
    const double t635 = 0.55170317879154401653e-1 * t216 * t633;
    const double t646 = t505 * t300;
    const double t647 = t646 * t245;
    const double t654 = 0.73e2 / 0.2025e4 * t233 * t292 - 0.73e2 / 0.1215e4 * t269 + 0.949e3 / 0.11664e6 * t272;
    const double t655 = t219 * t654;
    const double t659 = piecewise_functor_3( t77, 0.0, -0.16551095363746320496e0 * t500 * t301 - t635 + 0.33102190727492640992e0 * t589 * t647 - 0.16551095363746320496e0 * t216 * t655 );
    const double t661 = t261 * t261;
    const double t662 = t345 * t661;
    const double t665 = t145 * t57;
    const double t666 = t59 * t65;
    const double t667 = t665 * t666;
    const double t668 = t552 * t667;
    const double t671 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t143 * t662 - 0.84754509983741251008e-4 * t668 );
    const double t672 = t274 * t274;
    const double t673 = t506 * t672;
    const double t676 = t218 * t57;
    const double t677 = t59 * t113;
    const double t678 = t676 * t677;
    const double t679 = t589 * t678;
    const double t682 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t216 * t673 - 0.84754509983741251008e-4 * t679 );
    const double t683 = t616 * t261;
    const double t686 = t59 * t283;
    const double t687 = t665 * t686;
    const double t688 = t552 * t687;
    const double t691 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t552 * t683 + 0.33147598396528982063e-3 * t688 );
    const double t692 = t646 * t274;
    const double t695 = t59 * t296;
    const double t696 = t676 * t695;
    const double t697 = t589 * t696;
    const double t700 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t589 * t692 + 0.33147598396528982063e-3 * t697 );
    const double t701 = t287 * t287;
    const double t702 = t345 * t701;
    const double t706 = 0.1e1 / t36 / t147;
    const double t707 = t59 * t706;
    const double t708 = t665 * t707;
    const double t712 = piecewise_functor_3( t2, 0.0, 0.33102190727492640992e0 * t143 * t702 - 0.14916419278438041928e-2 * t552 * t708 );
    const double t713 = t300 * t300;
    const double t714 = t506 * t713;
    const double t718 = 0.1e1 / t88 / t220;
    const double t719 = t59 * t718;
    const double t720 = t676 * t719;
    const double t724 = piecewise_functor_3( t77, 0.0, 0.33102190727492640992e0 * t216 * t714 - 0.14916419278438041928e-2 * t589 * t720 );


    vrho_a = t76 + t124 + t7 * ( t177 + t193 );
    vrho_b = t76 + t124 + t7 * ( t206 + t250 );
    vsigma_aa = t7 * t265;
    vsigma_ab = 0.e0;
    vsigma_bb = t7 * t278;
    vlapl_a = 0.e0;
    vlapl_b = 0.e0;
    vtau_a = t7 * t291;
    vtau_b = t7 * t304;
    v2rho2_aa = 0.2e1 * t177 + 0.2e1 * t193 + t7 * ( t383 + t410 );
    v2rho2_bb = 0.2e1 * t206 + 0.2e1 * t250 + t7 * ( t484 + t544 );
    v2rhosigma_a_aa = t7 * t572 + t265;
    v2rhosigma_b_bb = t7 * t609 + t278;
    v2rholapl_a_a = 0.e0;
    v2rholapl_b_b = 0.e0;
    v2rhotau_a_a = t7 * t629 + t291;
    v2rhotau_b_b = t7 * t659 + t304;
    v2sigma2_aa_aa = t7 * t671;
    v2sigma2_bb_bb = t7 * t682;
    v2sigmalapl_aa_a = 0.e0;
    v2sigmalapl_bb_b = 0.e0;
    v2sigmatau_aa_a = t7 * t691;
    v2sigmatau_bb_b = t7 * t700;
    v2lapl2_aa = 0.e0;
    v2lapl2_bb = 0.e0;
    v2lapltau_a_a = 0.e0;
    v2lapltau_b_b = 0.e0;
    v2tau2_aa = t7 * t712;
    v2tau2_bb = t7 * t724;
    v2rho2_ab = 0.0;
    v2rhosigma_a_ab = 0.0;
    v2rhosigma_a_bb = 0.0;
    v2rhosigma_b_aa = 0.0;
    v2rhosigma_b_ab = 0.0;
    v2rholapl_a_b = 0.0;
    v2rholapl_b_a = 0.0;
    v2rhotau_a_b = 0.0;
    v2rhotau_b_a = 0.0;
    v2sigma2_aa_ab = 0.0;
    v2sigma2_aa_bb = 0.0;
    v2sigma2_ab_ab = 0.0;
    v2sigma2_ab_bb = 0.0;
    v2sigmalapl_aa_b = 0.0;
    v2sigmalapl_ab_a = 0.0;
    v2sigmalapl_ab_b = 0.0;
    v2sigmalapl_bb_a = 0.0;
    v2sigmatau_aa_b = 0.0;
    v2sigmatau_ab_a = 0.0;
    v2sigmatau_ab_b = 0.0;
    v2sigmatau_bb_a = 0.0;
    v2lapl2_ab = 0.0;
    v2lapltau_a_b = 0.0;
    v2lapltau_b_a = 0.0;
    v2tau2_ab = 0.0;



  }


};

struct BuiltinPKZB_X : detail::BuiltinKernelImpl< BuiltinPKZB_X > {

  BuiltinPKZB_X( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPKZB_X >(p) { }
  
  virtual ~BuiltinPKZB_X() = default;

};



} // namespace ExchCXX
