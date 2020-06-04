#pragma once

#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>


namespace ExchCXX {

template <>
struct kernel_traits<BuiltinPBE_C> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;
  static constexpr double dens_tol  = 1e-12;

  static constexpr double pbe_beta  = 0.06672455060314922;
  static constexpr double pbe_gamma = 0.031090690869654895034;
  static constexpr double pbe_B     = 1.;

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

      constexpr double t18 = t1 * t1;
      constexpr double t19 = t3 * t3;
      constexpr double t20 = t18 * t19;
      constexpr double t37 = constants::m_cbrt_2;
      constexpr double t44 = pbe_B * pbe_beta;
      constexpr double t45 = 1. / pbe_gamma;
      constexpr double t58 = t37 * t37;
      constexpr double t60 = 1. / t19;
      constexpr double t61 = t1 * t60;
      constexpr double t62 = t61 * t6;

      double t7 = std::pow(rho, constants::m_third);
      double t10 = t4 * t6 / t7;
      double t12 = 1. + 0.53425e-1 * t10;
      double t13 = std::sqrt(t10);
      double t16 = std::pow(t10,  3. / 2.);
      double t21 = t7 * t7;
      double t24 = t20 * t5 / t21;
      double t26 = 0.379785e1 * t13 + 0.8969 * t10 + 0.204775 * t16 + 
                   0.123235 * t24;
      double t29 = 1. + 0.16081979498692535067e2 / t26;
      double t30 = std::log(t29);
      double t31 = t12 * t30;
      double t32 = 0.621814e-1 * t31;
      double t33 = rho * rho;
      double t35 = 1. / t7 / t33;
      double t41 = t18 / t3 * t5;
      double t48 = std::exp(0.621814e-1 * t31 * t45);
      double t49 = t48 - 1.;
      double t50 = 1. / t49;
      double t51 = t45 * t50;
      double t52 = sigma * sigma;
      double t54 = t44 * t51 * t52;
      double t55 = t33 * t33;
      double t57 = 1. / t21 / t55;
      double t59 = t57 * t58;
      double t63 = t59 * t62;
      double t66 = sigma * t35 * t37 * t41 / 96. + t54 * t63 / 3072.;
      double t67 = pbe_beta * t66;
      double t68 = pbe_beta * t45;
      double t71 = t68 * t50 * t66 + 1.;
      double t72 = 1. / t71;
      double t73 = t45 * t72;
      double t75 = t67 * t73 + 1.;
      double t76 = std::log(t75);
      double t77 = pbe_gamma * t76;

      eps = t77 - t32;
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

      constexpr double t18 = t1 * t1;
      constexpr double t19 = t3 * t3;
      constexpr double t20 = t18 * t19;
      constexpr double t37 = constants::m_cbrt_2;
      constexpr double t44 = pbe_B * pbe_beta;
      constexpr double t45 = 1. / pbe_gamma;
      constexpr double t58 = t37 * t37;
      constexpr double t60 = 1. / t19;
      constexpr double t61 = t1 * t60;
      constexpr double t62 = t61 * t6;

      constexpr double t177 = pbe_beta * pbe_beta;
      constexpr double t179 = pbe_gamma * pbe_gamma;
      constexpr double t180 = 1. / t179;

      double t7 = std::pow(rho, constants::m_third);
      double t10 = t4 * t6 / t7;
      double t12 = 1. + 0.53425e-1 * t10;
      double t13 = std::sqrt(t10);
      double t16 = std::pow(t10,  3. / 2.);
      double t21 = t7 * t7;
      double t24 = t20 * t5 / t21;
      double t26 = 0.379785e1 * t13 + 0.8969 * t10 + 0.204775 * t16 + 
                   0.123235 * t24;
      double t29 = 1. + 0.16081979498692535067e2 / t26;
      double t30 = std::log(t29);
      double t31 = t12 * t30;
      double t32 = 0.621814e-1 * t31;
      double t33 = rho * rho;
      double t35 = 1. / t7 / t33;
      double t41 = t18 / t3 * t5;
      double t48 = std::exp(0.621814e-1 * t31 * t45);
      double t49 = t48 - 1.;
      double t50 = 1. / t49;
      double t51 = t45 * t50;
      double t52 = sigma * sigma;
      double t54 = t44 * t51 * t52;
      double t55 = t33 * t33;
      double t57 = 1. / t21 / t55;
      double t59 = t57 * t58;
      double t63 = t59 * t62;
      double t66 = sigma * t35 * t37 * t41 / 96. + t54 * t63 / 3072.;
      double t67 = pbe_beta * t66;
      double t68 = pbe_beta * t45;
      double t71 = t68 * t50 * t66 + 1.;
      double t72 = 1. / t71;
      double t73 = t45 * t72;
      double t75 = t67 * t73 + 1.;
      double t76 = std::log(t75);
      double t77 = pbe_gamma * t76;

      eps = t77 - t32;

      double t79 = 1. / t7 / rho;
      double t80 = t6 * t79;
      double t82 = t4 * t80 * t30;
      double t84 = t26 * t26;
      double t85 = 1. / t84;
      double t86 = t12 * t85;
      double t88 = 1. / t13 * t1;
      double t89 = t3 * t6;
      double t90 = t89 * t79;
      double t93 = t4 * t80;
      double t95 = sqrt(t10);
      double t96 = t95 * t1;
      double t104 = -0.632975 * t88 * t90 - 0.29896666666666666667 * t93 - 
                    0.1023875 * t96 * t90 - 
                    0.82156666666666666667e-1 * t20 * t5 / t21 / rho;
      double t105 = 1. / t29;
      double t106 = t104 * t105;
      double t107 = t86 * t106;
      double t109 = t33 * rho;
      double t111 = 1. / t7 / t109;
      double t116 = t44 * t45;
      double t117 = t49 * t49;
      double t118 = 1. / t117;
      double t119 = t118 * t52;
      double t121 = t116 * t119 * t57;
      double t122 = t58 * t1;
      double t123 = t122 * t60;
      double t124 = t4 * t6;
      double t132 = -0.11073470983333333333e-2 * t124 * t79 * t30 * t45 - 
                    t86 * t106 * t45;
      double t133 = t6 * t132;
      double t135 = t123 * t133 * t48;
      double t138 = t55 * rho;
      double t140 = 1. / t21 / t138;
      double t141 = t140 * t58;
      double t142 = t141 * t62;
      double t145 = -7. / 288. * sigma * t111 * t37 * t41 - 
                    t121 * t135 / 3072 - 7. / 4608 * t54 * t142;
      double t146 = pbe_beta * t145;
      double t148 = t71 * t71;
      double t149 = 1. / t148;
      double t150 = t45 * t149;
      double t151 = t68 * t118;
      double t152 = t66 * t132;
      double t157 = t68 * t50 * t145 - t151 * t152 * t48;
      double t158 = t150 * t157;
      double t160 = t146 * t73 - t67 * t158;
      double t162 = 1. / t75;
      double t163 = pbe_gamma * t160 * t162;

      vrho = -t32 + t77 + rho * (0.11073470983333333333e-2 * t82 + 
             t107 + t163);



      double t166 = rho * pbe_gamma;
      double t171 = t44 * t51 * sigma;
      double t174 = t35 * t37 * t41 / 96. + t171 * t63 / 1536.;
      double t175 = pbe_beta * t174;
      double t178 = t177 * t66;
      double t181 = t178 * t180;
      double t182 = t149 * t50;
      double t183 = t182 * t174;
      double t185 = t175 * t73 - t181 * t183;

      vsigma = t166 * t185 * t162;
  }

};


struct BuiltinPBE_C : detail::BuiltinKernelImpl< BuiltinPBE_C > {

  BuiltinPBE_C( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_C >(p) { }
  
  virtual ~BuiltinPBE_C() = default;

};


}
