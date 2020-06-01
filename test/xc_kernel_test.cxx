#include "catch2/catch.hpp"
#include <exchcxx/factory/xc_kernel.hpp>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>



using namespace ExchCXX;

constexpr std::array rho = {0.1, 0.2, 0.3, 0.4, 0.5};
constexpr std::array sigma = {0.2, 0.3, 0.4, 0.5, 0.6};

constexpr std::array rho_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
constexpr std::array sigma_polarized =
  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
   1.1, 1.2, 1.3, 1.5, 1.5 };

constexpr std::array exc_xc_lda_x_ref_unp = {
  -0.342808612301, -0.431911786723,
  -0.494415573788, -0.544174751790,
  -0.586194481348,
};

constexpr std::array vxc_xc_lda_x_ref_unp = {
  -0.457078149734, -0.575882382297,
  -0.659220765051, -0.725566335720,
  -0.781592641797
};

constexpr std::array exc_xc_lda_x_ref_pol = {
  -0.506753763434, -0.658748952120,
  -0.763800785778, -0.846274084184,
  -0.915314307811
};

constexpr std::array vxc_xc_lda_x_ref_pol = {
-0.575882382297, -0.725566335720,
-0.830566118415, -0.914156299468,
-0.984745021843, -1.046447735921,
-1.101623366705, -1.151764764594,
-1.197883627397, -1.240700981799
};



constexpr std::array exc_xc_gga_c_lyp_ref_unp = {
  -0.007040306272, -0.031424640440,
  -0.037479119388, -0.040429224120,
  -0.042290563929
};

constexpr std::array vrho_xc_gga_c_lyp_ref_unp = {
  -0.081854247031, -0.055198496086,
  -0.051617025994, -0.050995654065,
  -0.051084686930
};

constexpr std::array vsigma_xc_gga_c_lyp_ref_unp = {
  0.013598460611,
  0.004629650473,
  0.002429957976,
  0.001529632674,
  0.001065244937
};

constexpr std::array exc_xc_gga_c_lyp_ref_pol = {
  -0.031543975366, -0.043113613690,
  -0.046604883008, -0.048519647105,
  -0.049799110145
};

constexpr std::array vrho_xc_gga_c_lyp_ref_pol = {
  -0.089983823779, -0.035745759262,
  -0.062361000975, -0.045947114249,
  -0.059003605615, -0.049400798274,
  -0.058191535482, -0.051370405717,
  -0.058037798927, -0.052712179037
};

constexpr std::array vsigma_xc_gga_c_lyp_ref_pol = {
0.008447669161 , 0.006467154082,
-0.000638497084, 0.001421914705,
0.001031651601 , 0.000257537600,
0.000581699649 , 0.000435910598,
0.000202132738 , 0.000321427269,
0.000246773907 , 0.000146744820,
0.000206495563 , 0.000160996240,
0.000110118944
};


TEST_CASE("Hard coded", "[xc-builtin]") {

  const int npts = rho.size();
  std::vector<double> exc( npts );
  std::vector<double> vxc( npts );
  std::vector<double> vsigma( npts );

  double third = 1. / 3.;

  SECTION( "SLATER EXCHANGE" ) {

    XCKernel lda_x( XCKernel::Backend::builtin, 
      XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized ); 

    lda_x.eval_exc_vxc( npts, rho.data(), exc.data(), vxc.data() );
    for( auto i = 0; i < npts; ++i ) {
      CHECK( exc[i] == Approx(exc_xc_lda_x_ref_unp[i]) );
      CHECK( vxc[i] == Approx(vxc_xc_lda_x_ref_unp[i]) );
    }

  }

  SECTION( "LYP CORRELATION" ) {
    
    const double A = 0.04918;
    const double B = 0.132;
    const double c = 0.2533;
    const double d = 0.349;

    for( auto i = 0; i < npts; ++i ) {
      double t7 = std::pow( rho[i], third );;
      double t8 = 1. / t7;
      double t10 = d * t8 + 1.;
      double t11 = 1. / t10;
      double t13 = exp(-c * t8);
      double t14 = B * t13;
      double t15 = rho[i] * rho[i];
      double t16 = t7 * t7;
      double t18 = 1. / t16 / t15;
      double t19 = sigma[i] * t18;
      double t21 = d * t11 + c;
      double t22 = t21 * t8;
      double t24 = -1. / 72. - 7. / 72. * t22;
      double t26 = std::pow( 3., third );
      double t27 = t26 * t26;
      double t28 = M_PI * M_PI;
      double t29 = std::pow(t28,third);
      double t30 = t29 * t29;
      double t34 = 5. / 2. - t22 / 18.;
      double t35 = t34 * sigma[i];
      double t38 = t22 - 11.;
      double t39 = t38 * sigma[i];
      double t43 = -t19 * t24 - 3. / 10. * t27 * t30 + t35 * t18 / 8. + 
                   t39 * t18 / 144. - 5. / 24. * t19;

      exc[i] = A * (t14 * t11 * t43 - t11);

      double t47 = rho[i] * A;
      double t48 = t10 * t10;
      double t49 = 1. / t48;
      double t50 = t49 * d;
      double t52 = 1. / t7 / rho[i];
      double t55 = B * c;
      double t56 = t55 * t52;
      double t57 = t13 * t11;
      double t58 = t57 * t43;
      double t61 = t14 * t49;
      double t62 = t43 * d;
      double t66 = t15 * rho[i];
      double t68 = 1. / t16 / t66;
      double t69 = sigma[i] * t68;
      double t72 = d * d;
      double t73 = t72 * t49;
      double t75 = 1. / t16 / rho[i];
      double t78 = t21 * t52 - t73 * t75;
      double t79 = 7. / 216. * t78;
      double t81 = t78 / 54.;
      double t82 = t81 * sigma[i];
      double t88 = -t78 / 3.;
      double t89 = t88 * sigma[i];
      double t95 = 8. / 3. * t69 * t24 - t19 * t79 + t82 * t18 / 8. - 
                   t35 * t68 / 3. + t89 * t18 / 144. - t39 * t68 / 54. + 
                   5. / 9. * t69;
      double t98 = -t50 * t52 / 3. + t56 * t58 / 3. + t61 * t62 * t52 / 3. + 
                   t14 * t11 * t95;
      vxc[i] = t47 * t98 + (A * (t14 * t11 * t43 - t11));

      double t100 = t47 * B;
      double t107 = -t18 * t24 + t34 * t18 / 8. + t38 * t18 / 144. - 
                    5. / 24. * t18;
      double t108 = t57 * t107;

      vsigma[i] = t100 * t108;

      CHECK( exc[i]    == Approx(exc_xc_gga_c_lyp_ref_unp[i]) );
      CHECK( vxc[i]    == Approx(vrho_xc_gga_c_lyp_ref_unp[i]) );
      CHECK( vsigma[i] == Approx(vsigma_xc_gga_c_lyp_ref_unp[i]) );
    }
  }

  SECTION("PBE X") {

    XCKernel pbe_x( XCKernel::Kernel::PBE_X, XCKernel::Spin::Unpolarized );
    std::vector<double> exc_ref(npts), vrho_ref(npts), vsigma_ref(npts);

    pbe_x.eval_exc_vxc( npts, rho.data(), sigma.data(), exc_ref.data(), 
                        vrho_ref.data(), vsigma_ref.data() );


    const double kappa = 0.8040;
    const double mu    = 0.2195149727645171;

    for( auto i = 0; i < npts; ++i ) {

      double t1 = std::pow(3., third);
      double t3 = std::pow(1. / M_PI, third);
      double t4 = t1 * t3;
      double t5 = std::pow(4., third);
      double t6 = t5 * t5;
      double t7 = t4 * t6;
      double t8 = std::pow(2., third);
      double t9 = t8 * t8;
      double t10 = std::pow(rho[i], third); 
      double t12 = std::pow(6., third);;
      double t14 = M_PI * M_PI;
      double t15 = std::pow(t14,third);
      double t16 = t15 * t15;
      double t17 = 1. / t16;
      double t18 = mu * t12 * t17;
      double t20 = rho[i] * rho[i];
      double t21 = t10 * t10;
      double t23 = 1. / t21 / t20;
      double t27 = kappa + t18 * sigma[i] * t9 * t23 / 24.;
      double t32 = 1. + kappa * (1. - kappa / t27);
      double t34 = t7 * t9 * t10 * t32;

      exc[i] = -3. / 16. * t34;

      double t40 = t3 * t6;
      double t41 = t40 * t8;
      double t42 = 1. / t10 / t20 * t1 * t41;
      double t43 = kappa * kappa;
      double t44 = t27 * t27;
      double t46 = t43 / t44;
      double t49 = t12 * t17 * sigma[i];
      double t50 = t46 * mu * t49;

      double t57 = t46 * t18;

      vxc[i] = -t34 / 4. + t42 * t50 / 24.;
      vsigma[i] = -1. / t10 / rho[i] * t1 * t41 * t57 / 64.;

      CHECK( exc[i] == Approx(exc_ref[i]) );
      CHECK( vxc[i] == Approx(vrho_ref[i]) );
      CHECK( vsigma[i] == Approx(vsigma_ref[i]) );
    }
  }

  SECTION("PBE C") {

    XCKernel pbe_c( XCKernel::Kernel::PBE_C, XCKernel::Spin::Unpolarized );
    std::vector<double> exc_ref(npts), vrho_ref(npts), vsigma_ref(npts);

    pbe_c.eval_exc_vxc( npts, rho.data(), sigma.data(), exc_ref.data(), 
                        vrho_ref.data(), vsigma_ref.data() );


    const double beta  = 0.06672455060314922;
    const double gamma = 0.031090690869654895034;
    const double B     = 1.;

    for( auto i = 0; i < npts; ++i ) {

      double t1 = std::pow(3., third);
      double t2 = 1. / M_PI;
      double t3 = std::pow(t2, third);
      double t4 = t1 * t3;
      double t5 = std::pow(4., third);
      double t6 = t5 * t5;
      double t7 = std::pow(rho[i], third);
      double t10 = t4 * t6 / t7;
      double t12 = 1. + 0.53425e-1 * t10;
      double t13 = std::sqrt(t10);
      double t16 = std::pow(t10,3. / 2.);
      double t18 = t1 * t1;
      double t19 = t3 * t3;
      double t20 = t18 * t19;
      double t21 = t7 * t7;
      double t24 = t20 * t5 / t21;
      double t26 = 0.379785e1 * t13 + 0.8969 * t10 + 0.204775 * t16 + 0.123235 * t24;
      double t29 = 1. + 0.16081979498692535067e2 / t26;
      double t30 = std::log(t29);
      double t31 = t12 * t30;
      double t32 = 0.621814e-1 * t31;
      double t33 = rho[i] * rho[i];
      double t35 = 1. / t7 / t33;
      double t37 = std::pow(2., third);
      double t41 = t18 / t3 * t5;
      double t44 = B * beta;
      double t45 = 1. / gamma;
      double t48 = std::exp(0.621814e-1 * t31 * t45);
      double t49 = t48 - 1.;
      double t50 = 1. / t49;
      double t51 = t45 * t50;
      double t52 = sigma[i] * sigma[i];
      double t54 = t44 * t51 * t52;
      double t55 = t33 * t33;
      double t57 = 1. / t21 / t55;
      double t58 = t37 * t37;
      double t59 = t57 * t58;
      double t60 = 1. / t19;
      double t61 = t1 * t60;
      double t62 = t61 * t6;
      double t63 = t59 * t62;
      double t66 = sigma[i] * t35 * t37 * t41 / 96. + t54 * t63 / 3072.;
      double t67 = beta * t66;
      double t68 = beta * t45;
      double t71 = t68 * t50 * t66 + 1.;
      double t72 = 1. / t71;
      double t73 = t45 * t72;
      double t75 = t67 * t73 + 1.;
      double t76 = std::log(t75);
      double t77 = gamma * t76;

      exc[i] = t77 - t32;

      double t79 = 1. / t7 / rho[i];
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
                    0.82156666666666666667e-1 * t20 * t5 / t21 / rho[i];
      double t105 = 1. / t29;
      double t106 = t104 * t105;
      double t107 = t86 * t106;
      double t109 = t33 * rho[i];
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
      double t138 = t55 * rho[i];
      double t140 = 1. / t21 / t138;
      double t141 = t140 * t58;
      double t142 = t141 * t62;
      double t145 = -7. / 288. * sigma[i] * t111 * t37 * t41 - 
                    t121 * t135 / 3072 - 7. / 4608 * t54 * t142;
      double t146 = beta * t145;
      double t148 = t71 * t71;
      double t149 = 1. / t148;
      double t150 = t45 * t149;
      double t151 = t68 * t118;
      double t152 = t66 * t132;
      double t157 = t68 * t50 * t145 - t151 * t152 * t48;
      double t158 = t150 * t157;
      double t160 = t146 * t73 - t67 * t158;
      double t162 = 1. / t75;
      double t163 = gamma * t160 * t162;

      vxc[i] = -t32 + t77 + rho[i] * (0.11073470983333333333e-2 * t82 + 
               t107 + t163);

      double t166 = rho[i] * gamma;
      double t171 = t44 * t51 * sigma[i];
      double t174 = t35 * t37 * t41 / 96. + t171 * t63 / 1536.;
      double t175 = beta * t174;
      double t177 = beta * beta;
      double t178 = t177 * t66;
      double t179 = gamma * gamma;
      double t180 = 1. / t179;
      double t181 = t178 * t180;
      double t182 = t149 * t50;
      double t183 = t182 * t174;
      double t185 = t175 * t73 - t181 * t183;

      vsigma[i] = t166 * t185 * t162;

      CHECK( exc[i] == Approx(exc_ref[i]) );
      CHECK( vxc[i] == Approx(vrho_ref[i]) );
      CHECK( vsigma[i] == Approx(vsigma_ref[i]) );

    }
  }

}



TEST_CASE( "Unpolarized LDA Kernel Wrappers", "[xc-lda]" ) {

  const int npts = rho.size();
  std::vector<double> exc( npts );
  std::vector<double> vxc( npts );

  XCKernel lda( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );

  CHECK( lda.is_lda() );
  CHECK( not lda.is_polarized() );
  CHECK( not lda.is_gga() );
  CHECK( not lda.is_mgga() );
  CHECK( not lda.is_hyb() );

  CHECK( lda.hyb_exx() == Approx(0.0) );

  SECTION( "EXC only interface" ) {
    lda.eval_exc( npts, rho.data(), exc.data() );
    for( auto i = 0; i < npts; ++i )
      CHECK( exc[i] == Approx(exc_xc_lda_x_ref_unp[i]) );
  }
  
  SECTION( "EXC + VXC interface" ) {
    lda.eval_exc_vxc( npts, rho.data(), exc.data(), 
      vxc.data() );
    for( auto i = 0; i < npts; ++i ) {
      CHECK( exc[i] == Approx(exc_xc_lda_x_ref_unp[i]) );
      CHECK( vxc[i] == Approx(vxc_xc_lda_x_ref_unp[i]) );
    }
  }


}


TEST_CASE( "Polarized LDA Kernel Wrappers", "[xc-lda]" ) {

  const int npts = rho_polarized.size() / 2;
  std::vector<double> exc( npts );
  std::vector<double> vxc( 2*npts );

  XCKernel lda( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Polarized );

  CHECK( lda.is_lda() );
  CHECK( lda.is_polarized() );
  CHECK( not lda.is_gga() );
  CHECK( not lda.is_mgga() );
  CHECK( not lda.is_hyb() );

  CHECK( lda.hyb_exx() == Approx(0.0) );

  SECTION( "EXC only interface" ) {
    lda.eval_exc( npts, rho_polarized.data(), exc.data() );

    for( auto i = 0; i < npts; ++i ) 
      CHECK( exc[i] == Approx(exc_xc_lda_x_ref_pol[i]) );
  }
  
  SECTION( "EXC + VXC interface" ) {
    lda.eval_exc_vxc( npts, rho_polarized.data(), exc.data(), 
      vxc.data() );
    for( auto i = 0; i < npts; ++i ) {
      CHECK( exc[i] == Approx(exc_xc_lda_x_ref_pol[i]) );

      CHECK( vxc[2*i]   == Approx(vxc_xc_lda_x_ref_pol[2*i])   );
      CHECK( vxc[2*i+1] == Approx(vxc_xc_lda_x_ref_pol[2*i+1]) );
    }
  }

}






TEST_CASE( "Unpolarized GGA Kernel Wrappers", "[xc-gga]" ) {

  const int npts = rho.size();
  std::vector<double> exc( npts );
  std::vector<double> vrho( npts );
  std::vector<double> vsigma( npts );
  XCKernel lyp( XCKernel::Kernel::LYP, XCKernel::Spin::Unpolarized );

  CHECK( lyp.is_gga() );
  CHECK( not lyp.is_polarized() );
  CHECK( not lyp.is_lda() );
  CHECK( not lyp.is_mgga() );
  CHECK( not lyp.is_hyb() );

  CHECK( lyp.hyb_exx() == Approx(0.0) );

  SECTION( "EXC only interface" ) {
    lyp.eval_exc( npts, rho.data(), sigma.data(), exc.data() );
    for( auto i = 0; i < npts; ++i )
      CHECK( exc[i] == Approx(exc_xc_gga_c_lyp_ref_unp[i]) );
  }

  SECTION( "EXC + VXC interface" ) {
    lyp.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
      vrho.data(), vsigma.data() );
    for( auto i = 0; i < npts; ++i ) {
      CHECK( exc[i]    == Approx(exc_xc_gga_c_lyp_ref_unp[i])   );
      CHECK( vrho[i]   == Approx(vrho_xc_gga_c_lyp_ref_unp[i])  );
      CHECK( vsigma[i] == Approx(vsigma_xc_gga_c_lyp_ref_unp[i]));
    }
  }


}




TEST_CASE( "Polarized GGA Kernel Wrappers", "[xc-gga]" ) {

  const int npts = rho_polarized.size() / 2;
  std::vector<double> exc( npts );
  std::vector<double> vrho( 2*npts );
  std::vector<double> vsigma( 3*npts );
  XCKernel lyp( XCKernel::Kernel::LYP, XCKernel::Spin::Polarized );

  CHECK( lyp.is_gga() );
  CHECK( lyp.is_polarized() );
  CHECK( not lyp.is_lda() );
  CHECK( not lyp.is_mgga() );
  CHECK( not lyp.is_hyb() );

  CHECK( lyp.hyb_exx() == Approx(0.0) );

  SECTION( "EXC only interface" ) {
    lyp.eval_exc( npts, rho_polarized.data(), 
      sigma_polarized.data(), exc.data() );
    for( auto i = 0; i < npts; ++i )
      CHECK( exc[i] == Approx(exc_xc_gga_c_lyp_ref_pol[i]) );
  }

  SECTION( "EXC + VXC interface" ) {
    lyp.eval_exc_vxc( npts, rho_polarized.data(), 
      sigma_polarized.data(), exc.data(), vrho.data(), 
      vsigma.data() );
    for( auto i = 0; i < npts; ++i ) {
      CHECK( exc[i]        == Approx(exc_xc_gga_c_lyp_ref_pol[i])   );

      CHECK( vrho[2*i]     == Approx(vrho_xc_gga_c_lyp_ref_pol[2*i])  );
      CHECK( vrho[2*i+1]   == Approx(vrho_xc_gga_c_lyp_ref_pol[2*i+1])  );
      CHECK( vsigma[3*i]   == Approx(vsigma_xc_gga_c_lyp_ref_pol[3*i])  );
      CHECK( vsigma[3*i+1] == Approx(vsigma_xc_gga_c_lyp_ref_pol[3*i+1])  );
      CHECK( vsigma[3*i+2] == Approx(vsigma_xc_gga_c_lyp_ref_pol[3*i+2])  );
    }
  }


}



TEST_CASE( "Hybrid GGA Kernel Wrappers", "[xc-hyb-gga]" ) {

  XCKernel b3lyp( XCKernel::Kernel::B3LYP, XCKernel::Spin::Polarized );

  CHECK( b3lyp.is_gga() );
  CHECK( b3lyp.is_polarized() );
  CHECK( b3lyp.is_hyb() );
  CHECK( not b3lyp.is_lda() );
  CHECK( not b3lyp.is_mgga() );

  CHECK( b3lyp.hyb_exx() == Approx(0.2) );

}
