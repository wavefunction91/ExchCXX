#pragma once
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#include <exchcxx/impl/builtin/kernels/screening_interface.hpp>



namespace ExchCXX {

template <>
struct kernel_traits< BuiltinLYP > :
  public gga_screening_interface< BuiltinLYP > {

  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;

  static constexpr double dens_tol  = 1e-32;

  static constexpr bool is_hyb  = false;
  static constexpr double exx_coeff = 0.0;

  static constexpr double A = 0.04918;
  static constexpr double B = 0.132;
  static constexpr double c = 0.2533;
  static constexpr double d = 0.349;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar_impl( double rho, double sigma, double& eps ) {

    (void)(eps);
    constexpr double t26 = constants::m_cbrt_3;
    constexpr double t29 = constants::m_cbrt_pi_sq;
    constexpr double t27 = t26 * t26;
    constexpr double t30 = t29 * t29;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t10 = d * t8 + 0.1e1;
    const double t11 = 0.1e1 / t10;
    const double t13 = safe_math::exp( -c * t8 );
    const double t14 = B * t13;
    const double t15 = rho * rho;
    const double t16 = t7 * t7;
    const double t18 = 0.1e1 / t16 / t15;
    const double t19 = sigma * t18;
    const double t21 = d * t11 + c;
    const double t22 = t21 * t8;
    const double t24 = -0.1e1 / 0.72e2 - 0.7e1 / 0.72e2 * t22;
    const double t34 = 0.5e1 / 0.2e1 - t22 / 0.18e2;
    const double t35 = t34 * sigma;
    const double t38 = t22 - 0.11e2;
    const double t39 = t38 * sigma;
    const double t43 = -t19 * t24 - 0.3e1 / 0.10e2 * t27 * t30 + t35 * t18 / 0.8e1 + t39 * t18 / 0.144e3 - 0.5e1 / 0.24e2 * t19;


    eps = A * ( t14 * t11 * t43 - t11 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(eps);
    constexpr double t26 = constants::m_cbrt_3;
    constexpr double t29 = constants::m_cbrt_pi_sq;
    constexpr double t27 = t26 * t26;
    constexpr double t30 = t29 * t29;
    constexpr double t55 = B * c;
    constexpr double t72 = d * d;


    const double t7 = safe_math::cbrt( rho );
    const double t8 = 0.1e1 / t7;
    const double t10 = d * t8 + 0.1e1;
    const double t11 = 0.1e1 / t10;
    const double t13 = safe_math::exp( -c * t8 );
    const double t14 = B * t13;
    const double t15 = rho * rho;
    const double t16 = t7 * t7;
    const double t18 = 0.1e1 / t16 / t15;
    const double t19 = sigma * t18;
    const double t21 = d * t11 + c;
    const double t22 = t21 * t8;
    const double t24 = -0.1e1 / 0.72e2 - 0.7e1 / 0.72e2 * t22;
    const double t34 = 0.5e1 / 0.2e1 - t22 / 0.18e2;
    const double t35 = t34 * sigma;
    const double t38 = t22 - 0.11e2;
    const double t39 = t38 * sigma;
    const double t43 = -t19 * t24 - 0.3e1 / 0.10e2 * t27 * t30 + t35 * t18 / 0.8e1 + t39 * t18 / 0.144e3 - 0.5e1 / 0.24e2 * t19;
    const double t47 = rho * A;
    const double t48 = t10 * t10;
    const double t49 = 0.1e1 / t48;
    const double t50 = t49 * d;
    const double t52 = 0.1e1 / t7 / rho;
    const double t56 = t55 * t52;
    const double t57 = t13 * t11;
    const double t58 = t57 * t43;
    const double t61 = t14 * t49;
    const double t62 = t43 * d;
    const double t66 = t15 * rho;
    const double t68 = 0.1e1 / t16 / t66;
    const double t69 = sigma * t68;
    const double t73 = t72 * t49;
    const double t75 = 0.1e1 / t16 / rho;
    const double t78 = t21 * t52 - t73 * t75;
    const double t79 = 0.7e1 / 0.216e3 * t78;
    const double t81 = t78 / 0.54e2;
    const double t82 = t81 * sigma;
    const double t88 = -t78 / 0.3e1;
    const double t89 = t88 * sigma;
    const double t95 = 0.8e1 / 0.3e1 * t69 * t24 - t19 * t79 + t82 * t18 / 0.8e1 - t35 * t68 / 0.3e1 + t89 * t18 / 0.144e3 - t39 * t68 / 0.54e2 + 0.5e1 / 0.9e1 * t69;
    const double t98 = -t50 * t52 / 0.3e1 + t56 * t58 / 0.3e1 + t61 * t62 * t52 / 0.3e1 + t14 * t11 * t95;
    const double t100 = t47 * B;
    const double t107 = -t18 * t24 + t34 * t18 / 0.8e1 + t38 * t18 / 0.144e3 - 0.5e1 / 0.24e2 * t18;
    const double t108 = t57 * t107;


    eps = A * ( t14 * t11 * t43 - t11 );
    vrho = t47 * t98 + ( A * ( t14 * t11 * t43 - t11 ) );
    vsigma = t100 * t108;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_ferr_impl( double rho, double sigma, double& eps ) {

    (void)(rho);
    (void)(sigma);
    (void)(eps);




    eps = 0.0e0;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_ferr_impl( double rho, double sigma, double& eps, double& vrho, double& vsigma ) {

    (void)(rho);
    (void)(sigma);
    (void)(eps);




    eps = 0.0e0;
    vrho = 0.0e0;
    vsigma = 0.0e0;

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps ) {

    (void)(eps);
    constexpr double t38 = constants::m_cbrt_3;
    constexpr double t41 = constants::m_cbrt_pi_sq;
    constexpr double t60 = constants::m_cbrt_2;
    constexpr double t39 = t38 * t38;
    constexpr double t42 = t41 * t41;
    constexpr double t43 = t39 * t42;


    const double t7 = rho_a - rho_b;
    const double t8 = t7 * t7;
    const double t9 = rho_a + rho_b;
    const double t10 = t9 * t9;
    const double t11 = 0.1e1 / t10;
    const double t13 = -t8 * t11 + 0.1e1;
    const double t14 = safe_math::cbrt( t9 );
    const double t15 = 0.1e1 / t14;
    const double t17 = d * t15 + 0.1e1;
    const double t18 = 0.1e1 / t17;
    const double t21 = safe_math::exp( -c * t15 );
    const double t22 = B * t21;
    const double t24 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t25 = t14 * t14;
    const double t27 = 0.1e1 / t25 / t10;
    const double t28 = t24 * t27;
    const double t30 = d * t18 + c;
    const double t31 = t30 * t15;
    const double t33 = 0.47e2 - 0.7e1 * t31;
    const double t36 = t13 * t33 / 0.72e2 - 0.2e1 / 0.3e1;
    const double t44 = 0.1e1 / t9;
    const double t45 = t7 * t44;
    const double t46 = 0.1e1 + t45;
    const double t47 = t46 * t46;
    const double t48 = safe_math::cbrt( t46 );
    const double t49 = t48 * t48;
    const double t50 = t49 * t47;
    const double t51 = 0.1e1 - t45;
    const double t52 = t51 * t51;
    const double t53 = safe_math::cbrt( t51 );
    const double t54 = t53 * t53;
    const double t55 = t54 * t52;
    const double t56 = t50 + t55;
    const double t61 = t60 * t13;
    const double t63 = 0.5e1 / 0.2e1 - t31 / 0.18e2;
    const double t64 = rho_a * rho_a;
    const double t65 = safe_math::cbrt( rho_a );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t69 = sigma_aa * t68;
    const double t70 = t69 * t50;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = sigma_bb * t75;
    const double t77 = t76 * t55;
    const double t78 = t70 + t77;
    const double t79 = t63 * t78;
    const double t82 = t31 - 0.11e2;
    const double t84 = t49 * t47 * t46;
    const double t87 = t54 * t52 * t51;
    const double t89 = t69 * t84 + t76 * t87;
    const double t90 = t82 * t89;
    const double t95 = t47 * sigma_bb;
    const double t96 = t75 * t55;
    const double t99 = t52 * sigma_aa;
    const double t100 = t68 * t50;
    const double t106 = -t28 * t36 - 0.3e1 / 0.20e2 * t43 * t13 * t56 + t61 * t79 / 0.32e2 + t61 * t90 / 0.576e3 - t60 * ( 0.2e1 / 0.3e1 * t70 + 0.2e1 / 0.3e1 * t77 - t95 * t96 / 0.4e1 - t99 * t100 / 0.4e1 ) / 0.8e1;


    eps = A * ( t22 * t18 * t106 - t13 * t18 );

  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_polar_impl( double rho_a, double rho_b, double sigma_aa, double sigma_ab, double sigma_bb, double& eps, double& vrho_a, double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb ) {

    (void)(eps);
    constexpr double t38 = constants::m_cbrt_3;
    constexpr double t41 = constants::m_cbrt_pi_sq;
    constexpr double t60 = constants::m_cbrt_2;
    constexpr double t39 = t38 * t38;
    constexpr double t42 = t41 * t41;
    constexpr double t43 = t39 * t42;
    constexpr double t126 = B * c;
    constexpr double t143 = d * d;


    const double t7 = rho_a - rho_b;
    const double t8 = t7 * t7;
    const double t9 = rho_a + rho_b;
    const double t10 = t9 * t9;
    const double t11 = 0.1e1 / t10;
    const double t13 = -t8 * t11 + 0.1e1;
    const double t14 = safe_math::cbrt( t9 );
    const double t15 = 0.1e1 / t14;
    const double t17 = d * t15 + 0.1e1;
    const double t18 = 0.1e1 / t17;
    const double t21 = safe_math::exp( -c * t15 );
    const double t22 = B * t21;
    const double t24 = sigma_aa + 0.2e1 * sigma_ab + sigma_bb;
    const double t25 = t14 * t14;
    const double t27 = 0.1e1 / t25 / t10;
    const double t28 = t24 * t27;
    const double t30 = d * t18 + c;
    const double t31 = t30 * t15;
    const double t33 = 0.47e2 - 0.7e1 * t31;
    const double t36 = t13 * t33 / 0.72e2 - 0.2e1 / 0.3e1;
    const double t44 = 0.1e1 / t9;
    const double t45 = t7 * t44;
    const double t46 = 0.1e1 + t45;
    const double t47 = t46 * t46;
    const double t48 = safe_math::cbrt( t46 );
    const double t49 = t48 * t48;
    const double t50 = t49 * t47;
    const double t51 = 0.1e1 - t45;
    const double t52 = t51 * t51;
    const double t53 = safe_math::cbrt( t51 );
    const double t54 = t53 * t53;
    const double t55 = t54 * t52;
    const double t56 = t50 + t55;
    const double t61 = t60 * t13;
    const double t63 = 0.5e1 / 0.2e1 - t31 / 0.18e2;
    const double t64 = rho_a * rho_a;
    const double t65 = safe_math::cbrt( rho_a );
    const double t66 = t65 * t65;
    const double t68 = 0.1e1 / t66 / t64;
    const double t69 = sigma_aa * t68;
    const double t70 = t69 * t50;
    const double t71 = rho_b * rho_b;
    const double t72 = safe_math::cbrt( rho_b );
    const double t73 = t72 * t72;
    const double t75 = 0.1e1 / t73 / t71;
    const double t76 = sigma_bb * t75;
    const double t77 = t76 * t55;
    const double t78 = t70 + t77;
    const double t79 = t63 * t78;
    const double t82 = t31 - 0.11e2;
    const double t84 = t49 * t47 * t46;
    const double t87 = t54 * t52 * t51;
    const double t89 = t69 * t84 + t76 * t87;
    const double t90 = t82 * t89;
    const double t95 = t47 * sigma_bb;
    const double t96 = t75 * t55;
    const double t99 = t52 * sigma_aa;
    const double t100 = t68 * t50;
    const double t106 = -t28 * t36 - 0.3e1 / 0.20e2 * t43 * t13 * t56 + t61 * t79 / 0.32e2 + t61 * t90 / 0.576e3 - t60 * ( 0.2e1 / 0.3e1 * t70 + 0.2e1 / 0.3e1 * t77 - t95 * t96 / 0.4e1 - t99 * t100 / 0.4e1 ) / 0.8e1;
    const double t110 = t9 * A;
    const double t111 = t7 * t11;
    const double t112 = t10 * t9;
    const double t113 = 0.1e1 / t112;
    const double t114 = t8 * t113;
    const double t116 = -0.2e1 * t111 + 0.2e1 * t114;
    const double t118 = t17 * t17;
    const double t119 = 0.1e1 / t118;
    const double t120 = t13 * t119;
    const double t122 = 0.1e1 / t14 / t9;
    const double t123 = d * t122;
    const double t125 = t120 * t123 / 0.3e1;
    const double t127 = t126 * t122;
    const double t128 = t21 * t18;
    const double t129 = t128 * t106;
    const double t131 = t127 * t129 / 0.3e1;
    const double t132 = t22 * t119;
    const double t133 = t106 * d;
    const double t136 = t132 * t133 * t122 / 0.3e1;
    const double t138 = 0.1e1 / t25 / t112;
    const double t139 = t24 * t138;
    const double t141 = 0.8e1 / 0.3e1 * t139 * t36;
    const double t144 = t143 * t119;
    const double t146 = 0.1e1 / t25 / t9;
    const double t149 = t30 * t122 - t144 * t146;
    const double t150 = 0.7e1 / 0.3e1 * t149;
    const double t151 = t13 * t150;
    const double t153 = t116 * t33 / 0.72e2 + t151 / 0.72e2;
    const double t158 = t49 * t46;
    const double t159 = t44 - t111;
    const double t160 = t158 * t159;
    const double t161 = t54 * t51;
    const double t162 = -t159;
    const double t163 = t161 * t162;
    const double t165 = 0.8e1 / 0.3e1 * t160 + 0.8e1 / 0.3e1 * t163;
    const double t169 = t60 * t116;
    const double t172 = t149 / 0.54e2;
    const double t173 = t172 * t78;
    const double t175 = t61 * t173 / 0.32e2;
    const double t178 = 0.1e1 / t66 / t64 / rho_a;
    const double t179 = sigma_aa * t178;
    const double t180 = t179 * t50;
    const double t181 = t69 * t160;
    const double t182 = t76 * t163;
    const double t184 = -0.8e1 / 0.3e1 * t180 + 0.8e1 / 0.3e1 * t181 + 0.8e1 / 0.3e1 * t182;
    const double t185 = t63 * t184;
    const double t191 = -t149 / 0.3e1;
    const double t192 = t191 * t89;
    const double t194 = t61 * t192 / 0.576e3;
    const double t197 = t50 * t159;
    const double t200 = t55 * t162;
    const double t203 = -0.8e1 / 0.3e1 * t179 * t84 + 0.11e2 / 0.3e1 * t69 * t197 + 0.11e2 / 0.3e1 * t76 * t200;
    const double t204 = t82 * t203;
    const double t210 = t46 * sigma_bb;
    const double t214 = t75 * t161;
    const double t215 = t214 * t162;
    const double t218 = t51 * sigma_aa;
    const double t222 = t178 * t50;
    const double t225 = t68 * t158;
    const double t226 = t225 * t159;
    const double t232 = t141 - t28 * t153 - 0.3e1 / 0.20e2 * t43 * t116 * t56 - 0.3e1 / 0.20e2 * t43 * t13 * t165 + t169 * t79 / 0.32e2 + t175 + t61 * t185 / 0.32e2 + t169 * t90 / 0.576e3 + t194 + t61 * t204 / 0.576e3 - t60 * ( -0.16e2 / 0.9e1 * t180 + 0.16e2 / 0.9e1 * t181 + 0.16e2 / 0.9e1 * t182 - t210 * t96 * t159 / 0.2e1 - 0.2e1 / 0.3e1 * t95 * t215 - t218 * t100 * t162 / 0.2e1 + 0.2e1 / 0.3e1 * t99 * t222 - 0.2e1 / 0.3e1 * t99 * t226 ) / 0.8e1;
    const double t235 = t22 * t18 * t232 - t116 * t18 - t125 + t131 + t136;
    const double t238 = 0.2e1 * t111 + 0.2e1 * t114;
    const double t242 = t238 * t33 / 0.72e2 + t151 / 0.72e2;
    const double t247 = -t44 - t111;
    const double t248 = t158 * t247;
    const double t249 = -t247;
    const double t250 = t161 * t249;
    const double t252 = 0.8e1 / 0.3e1 * t248 + 0.8e1 / 0.3e1 * t250;
    const double t256 = t60 * t238;
    const double t259 = t69 * t248;
    const double t262 = 0.1e1 / t73 / t71 / rho_b;
    const double t263 = sigma_bb * t262;
    const double t264 = t263 * t55;
    const double t265 = t76 * t250;
    const double t267 = 0.8e1 / 0.3e1 * t259 - 0.8e1 / 0.3e1 * t264 + 0.8e1 / 0.3e1 * t265;
    const double t268 = t63 * t267;
    const double t273 = t50 * t247;
    const double t278 = t55 * t249;
    const double t281 = 0.11e2 / 0.3e1 * t69 * t273 - 0.8e1 / 0.3e1 * t263 * t87 + 0.11e2 / 0.3e1 * t76 * t278;
    const double t282 = t82 * t281;
    const double t288 = t96 * t247;
    const double t291 = t262 * t55;
    const double t294 = t214 * t249;
    const double t297 = t100 * t249;
    const double t300 = t225 * t247;
    const double t306 = t141 - t28 * t242 - 0.3e1 / 0.20e2 * t43 * t238 * t56 - 0.3e1 / 0.20e2 * t43 * t13 * t252 + t256 * t79 / 0.32e2 + t175 + t61 * t268 / 0.32e2 + t256 * t90 / 0.576e3 + t194 + t61 * t282 / 0.576e3 - t60 * ( 0.16e2 / 0.9e1 * t259 - 0.16e2 / 0.9e1 * t264 + 0.16e2 / 0.9e1 * t265 - t210 * t288 / 0.2e1 + 0.2e1 / 0.3e1 * t95 * t291 - 0.2e1 / 0.3e1 * t95 * t294 - t218 * t297 / 0.2e1 - 0.2e1 / 0.3e1 * t99 * t300 ) / 0.8e1;
    const double t309 = t22 * t18 * t306 - t238 * t18 - t125 + t131 + t136;
    const double t311 = t110 * B;
    const double t312 = t27 * t36;
    const double t313 = t63 * t68;
    const double t314 = t313 * t50;
    const double t317 = t82 * t68;
    const double t318 = t317 * t84;
    const double t322 = t52 * t68;
    const double t328 = -t312 + t61 * t314 / 0.32e2 + t61 * t318 / 0.576e3 - t60 * ( 0.2e1 / 0.3e1 * t100 - t322 * t50 / 0.4e1 ) / 0.8e1;
    const double t329 = t128 * t328;
    const double t330 = t146 * A;
    const double t331 = t330 * B;
    const double t332 = t128 * t36;
    const double t335 = t63 * t75;
    const double t336 = t335 * t55;
    const double t339 = t82 * t75;
    const double t340 = t339 * t87;
    const double t344 = t47 * t75;
    const double t350 = -t312 + t61 * t336 / 0.32e2 + t61 * t340 / 0.576e3 - t60 * ( 0.2e1 / 0.3e1 * t96 - t344 * t55 / 0.4e1 ) / 0.8e1;
    const double t351 = t128 * t350;


    eps = A * ( t22 * t18 * t106 - t13 * t18 );
    vrho_a = t110 * t235 + ( A * ( t22 * t18 * t106 - t13 * t18 ) );
    vrho_b = t110 * t309 + ( A * ( t22 * t18 * t106 - t13 * t18 ) );
    vsigma_aa = t311 * t329;
    vsigma_ab = -0.2e1 * t331 * t332;
    vsigma_bb = t311 * t351;

  }


};

struct BuiltinLYP : detail::BuiltinKernelImpl< BuiltinLYP > {

  BuiltinLYP( Spin p ) :
    detail::BuiltinKernelImpl< BuiltinLYP >(p) { }
  
  virtual ~BuiltinLYP() = default;

};



} // namespace ExchCXX