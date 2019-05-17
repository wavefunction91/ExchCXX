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






TEST_CASE( "Unpolarized LDA Kernel Wrappers", "[xc-lda]" ) {

  const int npts = rho.size();
  std::vector<double> exc( npts );
  std::vector<double> vxc( npts );

  XCKernel lda( "XC_LDA_X", XCKernel::Spin::Unpolarized );

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

  XCKernel lda( "XC_LDA_X", XCKernel::Spin::Polarized );

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




/*


TEST_CASE( "Unpolarized GGA Kernel Wrappers", "[xc-gga]" ) {

  const int npts = rho.size();
  std::vector<double> exc( npts );
  std::vector<double> vrho( npts );
  std::vector<double> vsigma( npts );
  XCKernel lyp( XC_GGA_C_LYP, XC_UNPOLARIZED );

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
  XCKernel lyp( XC_GGA_C_LYP, XC_POLARIZED );

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

  XCKernel b3lyp( XC_HYB_GGA_XC_B3LYP, XC_POLARIZED );

  CHECK( b3lyp.is_gga() );
  CHECK( b3lyp.is_polarized() );
  CHECK( b3lyp.is_hyb() );
  CHECK( not b3lyp.is_lda() );
  CHECK( not b3lyp.is_mgga() );

  CHECK( b3lyp.hyb_exx() == Approx(0.2) );

}
*/
