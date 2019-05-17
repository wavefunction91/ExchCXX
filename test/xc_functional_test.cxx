#include "catch2/catch.hpp"
#include <exchcxx/xc_functional.hpp>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <random>


template <typename T>
bool check_approx( const T& x, const T& y ) {
  if( std::isnan(x) )      return std::isnan(y);
  else if( std::isinf(x) ) return std::isinf(y);
  else                     return x == Approx(y);
}

using namespace ExchCXX;

TEST_CASE( "XC Functional Constructors", "[xc-misc]" ) {

  

  SECTION( "Pair vector copy constructor (copy)" ) {

    // Kernels
    XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );
    XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized );

    std::vector< std::pair<double, XCKernel> >
      kernels = {{ 1., slater}, {1., vwn3}};

    XCFunctional func(kernels);

  }

  SECTION( "Pair vector copy constructor (move)" ) {

    // Kernels
    XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );
    XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized );

    std::vector< std::pair<double, XCKernel> >
      kernels = {{ 1., std::move(slater)}, {1., std::move(vwn3)}};

    XCFunctional func(kernels);

  }

  SECTION( "Pair vector copy constructor (inplace)" ) {

    std::vector< std::pair<double, XCKernel> >
      kernels = {
        {1., XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized)}, 
        {1., XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)}
      };

    XCFunctional func(kernels);

  }


  SECTION( "Pair vector move constructor" ) {

    std::vector< std::pair<double, XCKernel> >
      kernels = {
        {1., XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized)},
        {1., XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)}
      };

    XCFunctional func(std::move(kernels));

  }

  SECTION( "Kernel vector copy constructor (copy)" ) {

    // Kernels
    XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );
    XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized );

    std::vector< XCKernel > kernels = {slater, vwn3};

    XCFunctional func(kernels);

  }

  SECTION( "Kernel vector copy constructor (move)" ) {

    // Kernels
    XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );
    XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized );

    std::vector< XCKernel > kernels = {std::move(slater), std::move(vwn3)};

    XCFunctional func(kernels);

  }

  SECTION( "Kernel vector copy constructor (inplace)" ) {

    std::vector< XCKernel >
      kernels = {
        XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized),
        XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)
      };

    XCFunctional func(kernels);

  }


  SECTION( "Kernel vector move constructor" ) {

    std::vector< XCKernel >
      kernels = {
        XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized),
        XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)
      };

    XCFunctional func(std::move(kernels));

  }

  SECTION( "Pair List constructor" ) {

    XCFunctional func(
      {
        {1., XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized)},
        {1., XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)}
      }
    );

  }

  SECTION( "Kernel List constructor" ) {

    XCFunctional func(
      {
        XCKernel(XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized),
        XCKernel(XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized)
      }
    );

  }
}

TEST_CASE( "Unpolarized LDA XC Functionals", "[xc-lda]" ) {

  const int npts = 100;

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  std::vector<double> rho( npts );

  std::generate( rho.begin(), rho.end(), 
    [&](){ return dist(gen); } );

  // Generate the reference data
  std::vector<double> exc_1( npts );
  std::vector<double> vxc_1( npts );
  std::vector<double> exc_2( npts );
  std::vector<double> vxc_2( npts );

  XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Unpolarized );
  XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Unpolarized );

  slater.eval_exc_vxc( npts, rho.data(), exc_1.data(), 
    vxc_1.data() );
  vwn3.eval_exc_vxc( npts, rho.data(), exc_2.data(), 
    vxc_2.data() );


  std::vector<double> exc( npts );
  std::vector<double> vxc( npts );

  SECTION( "Single kernel functionals" ) {

    std::vector<std::pair<double,XCKernel>> ks = {{1., slater}};

    XCFunctional func(ks);

    SECTION( "Exhibits the same meta data" ) {
      CHECK( func.is_lda()       == slater.is_lda()       );
      CHECK( func.is_polarized() == slater.is_polarized() );
      CHECK( func.is_gga()       == slater.is_gga()       );
      CHECK( func.is_mgga()      == slater.is_mgga()      );
      CHECK( func.is_hyb()       == slater.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( exc_1[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), exc.data(), 
        vxc.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( exc_1[i] ) );
        CHECK( vxc[i] == Approx( vxc_1[i] ) );
      }
    }

  }


  SECTION( "Two kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, slater},{beta, vwn3}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_lda()           );
      CHECK( not func.is_polarized() );
      CHECK( not func.is_gga()       );
      CHECK( not func.is_mgga()      );
      CHECK( not func.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), exc.data(), 
        vxc.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
        CHECK( vxc[i] == Approx( alpha*vxc_1[i] + beta*vxc_2[i] ) );
      }
    }

  }

}


TEST_CASE( "Polarized LDA XC Functionals", "[xc-lda]" ) {

  const int npts = 100;

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  std::vector<double> rho( 2*npts );

  std::generate( rho.begin(), rho.end(), 
    [&](){ return dist(gen); } );

  // Generate the reference data
  std::vector<double> exc_1( npts );
  std::vector<double> vxc_1( 2*npts );
  std::vector<double> exc_2( npts );
  std::vector<double> vxc_2( 2*npts );

  XCKernel slater( XCKernel::Kernel::SlaterExchange, XCKernel::Spin::Polarized );
  XCKernel vwn3  ( XCKernel::Kernel::VWN3,           XCKernel::Spin::Polarized );

  slater.eval_exc_vxc( npts, rho.data(), exc_1.data(), 
    vxc_1.data() );
  vwn3.eval_exc_vxc( npts, rho.data(), exc_2.data(), 
    vxc_2.data() );


  std::vector<double> exc( npts );
  std::vector<double> vxc( 2*npts );

  SECTION( "Single kernel functionals" ) {

    std::vector<std::pair<double,XCKernel>> ks = {{1., slater}};

    XCFunctional func(ks);

    SECTION( "Exhibits the same meta data" ) {
      CHECK( func.is_lda()       == slater.is_lda()       );
      CHECK( func.is_polarized() == slater.is_polarized() );
      CHECK( func.is_gga()       == slater.is_gga()       );
      CHECK( func.is_mgga()      == slater.is_mgga()      );
      CHECK( func.is_hyb()       == slater.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( exc_1[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), exc.data(), 
        vxc.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( exc_1[i] ) );
        CHECK( vxc[2*i] == Approx( vxc_1[2*i] ) );
        CHECK( vxc[2*i+1] == Approx( vxc_1[2*i+1] ) );
      }
    }

  }


  SECTION( "Two kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, slater},{beta, vwn3}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_lda()           );
      CHECK( func.is_polarized() );
      CHECK( not func.is_gga()       );
      CHECK( not func.is_mgga()      );
      CHECK( not func.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), exc.data(), 
        vxc.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
        CHECK( vxc[2*i] == Approx( alpha*vxc_1[2*i] + beta*vxc_2[2*i] ) );
        CHECK( vxc[2*i+1] == Approx( alpha*vxc_1[2*i+1] + beta*vxc_2[2*i+1] ) );
      }
    }

  }

}




TEST_CASE( "Unpolarized GGA XC Functionals", "[xc-gga]" ) {

  const int npts = 100;

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  std::vector<double> rho( npts );
  std::vector<double> sigma( npts );

  std::generate( rho.begin(), rho.end(), 
    [&](){ return dist(gen); } );
  std::generate( sigma.begin(), sigma.end(), 
    [&](){ return dist(gen); } );

  // Generate the reference data
  std::vector<double> exc_1( npts );
  std::vector<double> vrho_1( npts );
  std::vector<double> vsigma_1( npts );
  std::vector<double> exc_2( npts );
  std::vector<double> vrho_2( npts );
  std::vector<double> vsigma_2( npts );
  std::vector<double> exc_3( npts );
  std::vector<double> vxc_3( npts );

  XCKernel b88( XCKernel::Kernel::B88,  XCKernel::Spin::Unpolarized );
  XCKernel lyp( XCKernel::Kernel::LYP,  XCKernel::Spin::Unpolarized );
  XCKernel vwn( XCKernel::Kernel::VWN3, XCKernel::Spin::Unpolarized );

  b88.eval_exc_vxc( npts, rho.data(), sigma.data(), 
    exc_1.data(), vrho_1.data(), vsigma_1.data() );
  lyp.eval_exc_vxc( npts, rho.data(), sigma.data(), 
    exc_2.data(), vrho_2.data(), vsigma_2.data() );
  vwn.eval_exc_vxc( npts, rho.data(), exc_3.data(), 
    vxc_3.data());


  std::vector<double> exc( npts );
  std::vector<double> vrho( npts );
  std::vector<double> vsigma( npts );

  SECTION( "Single kernel functionals" ) {

    std::vector<std::pair<double,XCKernel>> ks = {{1., b88}};

    XCFunctional func(ks);

    SECTION( "Exhibits the same meta data" ) {
      CHECK( func.is_lda()       == b88.is_lda()       );
      CHECK( func.is_polarized() == b88.is_polarized() );
      CHECK( func.is_gga()       == b88.is_gga()       );
      CHECK( func.is_mgga()      == b88.is_mgga()      );
      CHECK( func.is_hyb()       == b88.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( exc_1[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( exc_1[i] ) );
        CHECK( vrho[i] == Approx( vrho_1[i] ) );
        CHECK( vsigma[i] == Approx( vsigma_1[i] ) );
      }
    }

  }


  SECTION( "Two GGA kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, b88},{beta, lyp}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_gga()           );
      CHECK( not func.is_lda()       );
      CHECK( not func.is_polarized() );
      CHECK( not func.is_mgga()      );
      CHECK( not func.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_2[i] ) );
        CHECK( vrho[i] == Approx( alpha*vrho_1[i] + beta*vrho_2[i] ) );
        CHECK( vsigma[i] == Approx( alpha*vsigma_1[i] + beta*vsigma_2[i] ) );
      }
    }

  }

  SECTION( "LDA + GGA kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, b88},{beta, vwn}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_gga()           );
      CHECK( not func.is_lda()       );
      CHECK( not func.is_polarized() );
      CHECK( not func.is_mgga()      );
      CHECK( not func.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_3[i] ) );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        CHECK( exc[i] == Approx( alpha*exc_1[i] + beta*exc_3[i] ) );
        CHECK( vrho[i] == Approx( alpha*vrho_1[i] + beta*vxc_3[i] ) );
        CHECK( vsigma[i] == Approx( alpha*vsigma_1[i] ) );
      }
    }

  }
}

TEST_CASE( "Polarized GGA XC Functionals", "[xc-gga]" ) {

  const int npts = 100;

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  std::vector<double> rho( 2*npts );
  std::vector<double> sigma( 3*npts );

  std::generate( rho.begin(), rho.end(), 
    [&](){ return dist(gen); } );
  std::generate( sigma.begin(), sigma.end(), 
    [&](){ return dist(gen); } );

  // Generate the reference data
  std::vector<double> exc_1( npts );
  std::vector<double> vrho_1( 2*npts );
  std::vector<double> vsigma_1( 3*npts );
  std::vector<double> exc_2( npts );
  std::vector<double> vrho_2( 2*npts );
  std::vector<double> vsigma_2( 3*npts );
  std::vector<double> exc_3( npts );
  std::vector<double> vxc_3( 2*npts );

  XCKernel b88( XCKernel::Kernel::B88,  XCKernel::Spin::Polarized );
  XCKernel lyp( XCKernel::Kernel::LYP,  XCKernel::Spin::Polarized );
  XCKernel vwn( XCKernel::Kernel::VWN3, XCKernel::Spin::Polarized );

  b88.eval_exc_vxc( npts, rho.data(), sigma.data(), 
    exc_1.data(), vrho_1.data(), vsigma_1.data() );
  lyp.eval_exc_vxc( npts, rho.data(), sigma.data(), 
    exc_2.data(), vrho_2.data(), vsigma_2.data() );
  vwn.eval_exc_vxc( npts, rho.data(), exc_3.data(), 
    vxc_3.data());


  std::vector<double> exc( npts );
  std::vector<double> vrho( 2*npts );
  std::vector<double> vsigma( 3*npts );

  SECTION( "Single kernel functionals" ) {

    std::vector<std::pair<double,XCKernel>> ks = {{1., b88}};

    XCFunctional func(ks);

    SECTION( "Exhibits the same meta data" ) {
      CHECK( func.is_lda()       == b88.is_lda()       );
      CHECK( func.is_polarized() == b88.is_polarized() );
      CHECK( func.is_gga()       == b88.is_gga()       );
      CHECK( func.is_mgga()      == b88.is_mgga()      );
      CHECK( func.is_hyb()       == b88.is_hyb()       );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i) {
        //CHECK( exc[i] == Approx( exc_1[i] ) );
        check_approx(exc[i],exc_1[i]);
      }
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        check_approx(exc[i],exc_1[i]);
        check_approx( vrho[2*i], vrho_1[2*i] );
        check_approx( vrho[2*i+1], vrho_1[2*i+1] );
        check_approx( vsigma[3*i], vsigma_1[3*i] );
        check_approx( vsigma[3*i+1], vsigma_1[3*i+1] );
        check_approx( vsigma[3*i+2], vsigma_1[3*i+2] );
      }
    }

  }


  SECTION( "Two GGA kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, b88},{beta, lyp}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_gga()       );
      CHECK( func.is_polarized() );
      CHECK( not func.is_lda()   );
      CHECK( not func.is_mgga()  );
      CHECK( not func.is_hyb()   );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        check_approx( exc[i], alpha*exc_1[i] + beta*exc_2[i] );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        check_approx( exc[i], alpha*exc_1[i] + beta*exc_2[i] );
        check_approx( vrho[2*i], alpha*vrho_1[2*i] + beta*vrho_2[2*i] );
        check_approx( vrho[2*i+1], alpha*vrho_1[2*i+1] + beta*vrho_2[2*i+1] );
        check_approx( vsigma[3*i], alpha*vsigma_1[3*i] + beta*vsigma_2[3*i] );
        check_approx( vsigma[3*i+1], alpha*vsigma_1[3*i+1] + beta*vsigma_2[3*i+1] );
        check_approx( vsigma[3*i+2], alpha*vsigma_1[3*i+2] + beta*vsigma_2[3*i+2] );
      }
    }

  }

  SECTION( "LDA + GGA kernel functionals" ) {

    const double alpha = dist(gen);
    const double beta  = dist(gen);
    std::vector<std::pair<double,XCKernel>> ks = {{alpha, b88},{beta, vwn}};

    XCFunctional func(ks);

    SECTION( "Exhibits the proper data" ) {
      CHECK( func.is_gga()       );
      CHECK( func.is_polarized() );
      CHECK( not func.is_lda()   );
      CHECK( not func.is_mgga()  );
      CHECK( not func.is_hyb()   );
    }


    SECTION( "EXC only interface" ) {
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      for(auto i = 0; i < npts; ++i)
        check_approx( exc[i], alpha*exc_1[i] + beta*exc_3[i] );
    }

    SECTION( "EXC + VXC interface" ) {
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(), 
        vrho.data(), vsigma.data() );

      for(auto i = 0; i < npts; ++i) {
        check_approx( exc[i], alpha*exc_1[i] + beta*exc_3[i] );
        check_approx( vrho[2*i], alpha*vrho_1[2*i] + beta*vxc_3[2*i] );
        check_approx( vrho[2*i+1], alpha*vrho_1[2*i+1] + beta*vxc_3[2*i+1] );
        check_approx( vsigma[3*i], alpha*vsigma_1[3*i] );
        check_approx( vsigma[3*i+1], alpha*vsigma_1[3*i+1] );
        check_approx( vsigma[3*i+2], alpha*vsigma_1[3*i+2] );
      }
    }

  }
}
