#pragma once

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <cmath>

#include <exchcxx/impl/builtin/fwd.hpp>
#include <exchcxx/impl/builtin/constants.hpp>
#include <exchcxx/impl/builtin/kernel.hpp>
#include <exchcxx/impl/builtin/util.hpp>

#ifdef __CUDACC__

#define BUILTIN_KERNEL_EVAL_RETURN static inline constexpr void __host__ __device__

#else

#define BUILTIN_KERNEL_EVAL_RETURN static inline constexpr void 

#endif

namespace ExchCXX {


template <>
struct kernel_traits<BuiltinSlaterExchange> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = true;
  static constexpr bool is_gga  = false;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;

  BUILTIN_KERNEL_EVAL_RETURN 
    eval_exc_unpolar( double rho, double& eps ) {

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
    eval_exc_vxc_unpolar( double rho, double& eps, double& vxc ) {

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






template <>
struct kernel_traits<BuiltinLYP> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;


  static constexpr double lyp_A = 0.04918;
  static constexpr double lyp_B = 0.132;
  static constexpr double lyp_c = 0.2533;
  static constexpr double lyp_d = 0.349;


  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

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





template <>
struct kernel_traits<BuiltinPBE_X> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;

  static constexpr double pbe_kappa = 0.8040;
  static constexpr double pbe_mu    = 0.2195149727645171;


  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {


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


template <>
struct kernel_traits<BuiltinPBE_C> {

  static constexpr bool is_hyb  = false;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.;

  static constexpr double pbe_beta  = 0.06672455060314922;
  static constexpr double pbe_gamma = 0.031090690869654895034;
  static constexpr double pbe_B     = 1.;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

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



template <>
struct kernel_traits<BuiltinPBE0> {

  static constexpr bool is_hyb  = true;
  static constexpr bool is_lda  = false;
  static constexpr bool is_gga  = true;
  static constexpr bool is_mgga = false;
  static constexpr double exx_coeff = 0.25;

  using pbe_x_traits = kernel_traits<BuiltinPBE_X>;
  using pbe_c_traits = kernel_traits<BuiltinPBE_C>;

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_unpolar( double rho, double sigma, double& eps ) {

    pbe_x_traits::eval_exc_unpolar( rho, sigma, eps );
    double eps_x = eps;

    pbe_c_traits::eval_exc_unpolar( rho, sigma, eps );

    eps = (1. - exx_coeff) * eps_x + eps;
  }

  BUILTIN_KERNEL_EVAL_RETURN
    eval_exc_vxc_unpolar( double rho, double sigma, double& eps, double& vrho,
      double& vsigma ) {

    pbe_x_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );
    double eps_x    = eps;
    double vrho_x   = vrho;
    double vsigma_x = vsigma;


    pbe_c_traits::eval_exc_vxc_unpolar( rho, sigma, eps, vrho, vsigma );

    eps    = (1. - exx_coeff) * eps_x    + eps;
    vrho   = (1. - exx_coeff) * vrho_x   + vrho;
    vsigma = (1. - exx_coeff) * vsigma_x + vsigma;

  }

};


struct BuiltinSlaterExchange : detail::BuiltinKernelImpl< BuiltinSlaterExchange > {

  BuiltinSlaterExchange( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinSlaterExchange >(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

};

struct BuiltinLYP : detail::BuiltinKernelImpl< BuiltinLYP > {

  BuiltinLYP( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinLYP >(p) { }
  
  virtual ~BuiltinLYP() = default;

};

struct BuiltinPBE_X : detail::BuiltinKernelImpl< BuiltinPBE_X > {

  BuiltinPBE_X( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_X >(p) { }
  
  virtual ~BuiltinPBE_X() = default;

};

struct BuiltinPBE_C : detail::BuiltinKernelImpl< BuiltinPBE_C > {

  BuiltinPBE_C( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE_C >(p) { }
  
  virtual ~BuiltinPBE_C() = default;

};

struct BuiltinPBE0 : detail::BuiltinKernelImpl< BuiltinPBE0 > {

  BuiltinPBE0( XCKernel::Spin p ) :
    detail::BuiltinKernelImpl< BuiltinPBE0 >(p) { }
  
  virtual ~BuiltinPBE0() = default;

};

}

#include <exchcxx/impl/builtin/interface.hpp>
