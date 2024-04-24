/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
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

#include "ut_common.hpp"

using namespace ExchCXX;

TEST_CASE( "XC Functional Constructors", "[xc-misc]" ) {

  

  SECTION( "Pair vector copy constructor (copy)" ) {

    // Kernels
    XCKernel slater( Kernel::SlaterExchange, Spin::Unpolarized );
    XCKernel vwn3  ( Kernel::VWN3,           Spin::Unpolarized );

    std::vector< std::pair<double, XCKernel> >
      kernels = {{ 1., slater}, {1., vwn3}};

    XCFunctional func(kernels);

  }

  SECTION( "Pair vector copy constructor (move)" ) {

    // Kernels
    XCKernel slater( Kernel::SlaterExchange, Spin::Unpolarized );
    XCKernel vwn3  ( Kernel::VWN3,           Spin::Unpolarized );

    std::vector< std::pair<double, XCKernel> >
      kernels = {{ 1., std::move(slater)}, {1., std::move(vwn3)}};

    XCFunctional func(kernels);

  }

  SECTION( "Pair vector copy constructor (inplace)" ) {

    std::vector< std::pair<double, XCKernel> >
      kernels = {
        {1., XCKernel(Kernel::SlaterExchange, Spin::Unpolarized)}, 
        {1., XCKernel(Kernel::VWN3,           Spin::Unpolarized)}
      };

    XCFunctional func(kernels);

  }


  SECTION( "Pair vector move constructor" ) {

    std::vector< std::pair<double, XCKernel> >
      kernels = {
        {1., XCKernel(Kernel::SlaterExchange, Spin::Unpolarized)},
        {1., XCKernel(Kernel::VWN3,           Spin::Unpolarized)}
      };

    XCFunctional func(std::move(kernels));

  }

  SECTION( "Kernel vector copy constructor (copy)" ) {

    // Kernels
    XCKernel slater( Kernel::SlaterExchange, Spin::Unpolarized );
    XCKernel vwn3  ( Kernel::VWN3,           Spin::Unpolarized );

    std::vector< XCKernel > kernels = {slater, vwn3};

    XCFunctional func(kernels);

  }

  SECTION( "Kernel vector copy constructor (move)" ) {

    // Kernels
    XCKernel slater( Kernel::SlaterExchange, Spin::Unpolarized );
    XCKernel vwn3  ( Kernel::VWN3,           Spin::Unpolarized );

    std::vector< XCKernel > kernels = {std::move(slater), std::move(vwn3)};

    XCFunctional func(kernels);

  }

  SECTION( "Kernel vector copy constructor (inplace)" ) {

    std::vector< XCKernel >
      kernels = {
        XCKernel(Kernel::SlaterExchange, Spin::Unpolarized),
        XCKernel(Kernel::VWN3,           Spin::Unpolarized)
      };

    XCFunctional func(kernels);

  }


  SECTION( "Kernel vector move constructor" ) {

    std::vector< XCKernel >
      kernels = {
        XCKernel(Kernel::SlaterExchange, Spin::Unpolarized),
        XCKernel(Kernel::VWN3,           Spin::Unpolarized)
      };

    XCFunctional func(std::move(kernels));

  }

  SECTION( "Pair List constructor" ) {

    XCFunctional func(
      {
        std::pair<double,XCKernel>{1., XCKernel(Kernel::SlaterExchange, Spin::Unpolarized)},
        std::pair<double,XCKernel>{1., XCKernel(Kernel::VWN3,           Spin::Unpolarized)}
      }
    );

  }

  SECTION( "Kernel List constructor" ) {

    XCFunctional func(
      {
        XCKernel(Kernel::SlaterExchange, Spin::Unpolarized),
        XCKernel(Kernel::VWN3,           Spin::Unpolarized)
      }
    );

  }
}


template <typename... Args>
void check_meta( Backend backend, Spin polar, Args&&... args ) {

  std::vector<XCKernel> kerns = { XCKernel(backend, args, polar)... };

  XCFunctional func( kerns );

  bool should_be_mgga = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.is_mgga(); } );
  bool should_be_gga = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.is_gga(); } ) 
    and not should_be_mgga;
  bool should_be_lda = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.is_lda(); } ) 
    and not should_be_mgga
    and not should_be_gga;

  bool should_be_polarized = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.is_polarized(); } ); 
  bool should_be_hyb = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.is_hyb(); } );
  bool should_need_lapl = std::any_of( kerns.begin(), kerns.end(),
    [](const auto& k) { return k.needs_laplacian(); } );

  CHECK( func.is_mgga()       == should_be_mgga       );
  CHECK( func.is_gga()        == should_be_gga        );
  CHECK( func.is_lda()        == should_be_lda        );
  CHECK( func.is_polarized()  == should_be_polarized  );
  CHECK( func.is_hyb()        == should_be_hyb        );
  CHECK( func.needs_laplacian() == should_need_lapl   );

  double total_exx = std::accumulate( kerns.begin(), kerns.end(), 0.,
    [](const auto &a, const auto &b){ return a + b.hyb_exx(); } );
  CHECK( func.hyb_exx() == Approx( total_exx ) );

}


TEST_CASE( "XC Functional Metadata", "[xc-meta]" ) {

  SECTION( "LDA" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::SlaterExchange 
    );
  }
  SECTION( "LDA + LDA" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::SlaterExchange, Kernel::VWN5 
    );
  }


  SECTION( "GGA" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B88 
    );
  }
  SECTION( "GGA + GGA" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B88, Kernel::LYP 
    );
  }


  SECTION( "Hyb" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B3LYP
    );
  }
  SECTION( "Hyb + Hyb" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B3LYP, Kernel::PBE0
    );
  }

  SECTION( "LDA + GGA" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B88, Kernel::VWN3 
    );
  }
  SECTION( "Hyb + Pure" ) {
    check_meta( Backend::libxc, Spin::Unpolarized, 
      Kernel::B3LYP, Kernel::LYP
    );
  }

  SECTION( "Polarized" ) {
    check_meta( Backend::libxc, Spin::Polarized, 
      Kernel::SlaterExchange, Kernel::VWN3
    );
  }

}


void gen_kern_vector( std::vector<std::pair<double,XCKernel>>& vec, Backend backend, Spin polar, double f, Kernel kern ) {
  vec.emplace_back( std::make_pair(f, XCKernel( backend, kern, polar ) ) );
}

template <typename... Args>
void gen_kern_vector( std::vector<std::pair<double,XCKernel>>& vec, Backend backend, Spin polar, double f, Kernel kern, Args&&... args ) {
  gen_kern_vector( vec, backend, polar, f, kern );
  gen_kern_vector( vec, backend, polar, std::forward<Args>(args)... ); 
}

template <typename... Args>
auto gen_kern_vector( Backend backend, Spin polar, Args&&... args ) {
  std::vector< std::pair<double, XCKernel> > vec;
  gen_kern_vector( vec, backend, polar, std::forward<Args>(args)... );
  return vec;
}






template <typename... Args>
void check_correctness( TestInterface interface, Backend backend, Spin polar, 
  Args&&... args ) {

  auto kern_pairs = gen_kern_vector( backend, polar, std::forward<Args>(args)... );
  XCFunctional func( kern_pairs );

  std::vector<double> coeffs;
  std::for_each( kern_pairs.begin(), kern_pairs.end(), 
    [&](const auto& a){ coeffs.emplace_back( a.first ); } 
  );

  size_t npts = 1024;

  std::vector<double>
    rho( func.rho_buffer_len(npts) ),
    sigma( func.sigma_buffer_len(npts) ),
    lapl( func.lapl_buffer_len(npts) ),
    tau( func.tau_buffer_len(npts) );

  // Randomly generate Rho / Sigma
  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);
  std::generate( rho.begin(), rho.end(), 
    [&](){ return dist(gen); } );
  std::generate( sigma.begin(), sigma.end(), 
    [&](){ return dist(gen); } );
  std::generate( lapl.begin(), lapl.end(), 
    [&](){ return dist(gen); } );
  std::generate( tau.begin(), tau.end(), 
    [&](){ return dist(gen); } );


  // Evaluate individual kernels
  std::vector< std::vector<double> >
    eps_refs, vrho_refs, vsigma_refs, vlapl_refs, vtau_refs;

  for( auto p : kern_pairs ) {

    XCKernel& kern = p.second;
    //double    f    = p.first;

    eps_refs.emplace_back( kern.exc_buffer_len(npts) );
    vrho_refs.emplace_back( kern.vrho_buffer_len(npts) );
    vsigma_refs.emplace_back( kern.vsigma_buffer_len(npts) );
    vlapl_refs.emplace_back( kern.vlapl_buffer_len(npts) );
    vtau_refs.emplace_back( kern.vtau_buffer_len(npts) );

    if( interface == TestInterface::EXC ) {

      if( kern.is_lda() )
        kern.eval_exc( npts, rho.data(), eps_refs.back().data() );
      else if( kern.is_gga() )
        kern.eval_exc( npts, rho.data(), sigma.data(), 
          eps_refs.back().data() );
      else if( kern.is_mgga() )
        kern.eval_exc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
          eps_refs.back().data() );

    } else if( interface == TestInterface::EXC_VXC ) {

      if( kern.is_lda() )
        kern.eval_exc_vxc( npts, rho.data(), eps_refs.back().data(), 
          vrho_refs.back().data() );
      else if( kern.is_gga() )
        kern.eval_exc_vxc( npts, rho.data(), sigma.data(), 
          eps_refs.back().data(), vrho_refs.back().data(),
          vsigma_refs.back().data() );
      else if( kern.is_mgga() )
        kern.eval_exc_vxc( npts, rho.data(), sigma.data(), 
          lapl.data(), tau.data(), eps_refs.back().data(), 
          vrho_refs.back().data(), vsigma_refs.back().data(),
          vlapl_refs.back().data(), vtau_refs.back().data() );

    }
    
  }

  
  // Evaluation combined functional
  std::vector<double>
    eps( func.exc_buffer_len( npts ) ),
    vrho( func.vrho_buffer_len( npts ) ),
    vsigma( func.vsigma_buffer_len( npts ) ),
    vlapl( func.vlapl_buffer_len( npts ) ),
    vtau( func.vtau_buffer_len( npts ) );

  if( func.is_lda() )
    func.eval_exc_vxc( npts, rho.data(), eps.data(), vrho.data() );
  else if( func.is_gga() )
    func.eval_exc_vxc( npts, rho.data(), sigma.data(), 
      eps.data(), vrho.data(), vsigma.data() );
  else if( func.is_mgga() )
    func.eval_exc_vxc( npts, rho.data(), sigma.data(), 
      lapl.data(), tau.data(), eps.data(), vrho.data(), 
      vsigma.data(), vlapl.data(), vtau.data() );


  for( auto i = 0ul; i < eps.size(); ++i ) {
    double ref_val = 0.;
    for( auto j = 0ul; j < kern_pairs.size(); ++j )
      ref_val += coeffs[j] * eps_refs[j][i];
    CHECK( eps[i] == Approx( ref_val ) );
  }

  if( interface != TestInterface::EXC ) {

    for( auto i = 0ul; i < vrho.size(); ++i ) {
      double ref_val = 0.;
      for( auto j = 0ul; j < kern_pairs.size(); ++j )
        ref_val += coeffs[j] * vrho_refs[j][i];
      CHECK( vrho[i] == Approx( ref_val ) );
    }

    if( func.is_gga() ) 
      for( auto i = 0ul; i < vsigma.size(); ++i ) {
        double ref_val = 0.;
        for( auto j = 0ul; j < kern_pairs.size(); ++j )
        if( kern_pairs[j].second.is_gga() )
          ref_val += coeffs[j] * vsigma_refs[j][i];

        CHECK( vsigma[i] == Approx( ref_val ) );
      }

    if( func.is_mgga() )  {
      for( auto i = 0ul; i < vlapl.size(); ++i ) {
        double ref_val = 0.;
        for( auto j = 0ul; j < kern_pairs.size(); ++j )
        if( kern_pairs[j].second.needs_laplacian() )
          ref_val += coeffs[j] * vlapl_refs[j][i];

        CHECK( vlapl[i] == Approx( ref_val ) );
      }
      for( auto i = 0ul; i < vtau.size(); ++i ) {
        double ref_val = 0.;
        for( auto j = 0ul; j < kern_pairs.size(); ++j )
        if( kern_pairs[j].second.is_mgga() )
          ref_val += coeffs[j] * vtau_refs[j][i];

        CHECK( vtau[i] == Approx( ref_val ) );
      }
    }

  }

}

TEST_CASE( "LDA XC Functionals", "[xc-lda]" ) {

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  SECTION("LDA Unpolarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SlaterExchange 
    );
  }

  SECTION("LDA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SlaterExchange 
    );
  }

  SECTION("LDA + LDA Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("LDA + LDA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("LDA Polarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange 
    );
  }

  SECTION("LDA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange 
    );
  }

  SECTION("LDA + LDA Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("LDA + LDA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::VWN5 
    );
  }

}


TEST_CASE( "GGA XC Functionals", "[xc-gga]" ) {

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  SECTION("GGA Unpolarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88 
    );
  }

  SECTION("GGA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88 
    );
  }

  SECTION("GGA + GGA Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::LYP 
    );
  }

  SECTION("GGA + GGA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::LYP 
    );
  }

  SECTION("GGA + LDA Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("GGA + LDA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::VWN5 
    );
  }







  SECTION("GGA Polarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88 
    );
  }

  SECTION("GGA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88 
    );
  }

  SECTION("GGA + GGA Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::LYP 
    );
  }

  SECTION("GGA + GGA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::LYP 
    );
  }

  SECTION("GGA + LDA Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("GGA + LDA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::B88,
      dist(gen), Kernel::VWN5 
    );
  }






  SECTION("LDA + GGA: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::LYP 
    );
  }
  SECTION("LDA + GGA: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::LYP 
    );
  }

}


TEST_CASE( "MGGA XC Functionals", "[xc-mgga]" ) {

  std::default_random_engine gen;
  std::normal_distribution<> dist(5.0,2.0);

  SECTION("MGGA-LAPL Unpolarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X 
    );
  }

  SECTION("MGGA-LAPL Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X 
    );
  }

  SECTION("MGGA-LAPL + MGGA-LAPL Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

  SECTION("MGGA-LAPL + MGGA-LAPL Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

  SECTION("MGGA-LAPL + LDA Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("MGGA-LAPL + LDA Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("MGGA-LAPL Unpolarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::R2SCANL_X 
    );
  }

  SECTION("MGGA-TAU Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SCAN_X 
    );
  }

  SECTION("MGGA-TAU + MGGA-TAU Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::SCAN_C 
    );
  }

  SECTION("MGGA-TAU + MGGA-TAU Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::SCAN_C 
    );
  }

  SECTION("MGGA-TAU + MGGA-LAPL Unpolarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

  SECTION("MGGA-TAU + MGGA-LAPL Unpolarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Unpolarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }






  SECTION("MGGA-LAPL Polarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X 
    );
  }

  SECTION("MGGA-LAPL Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X 
    );
  }

  SECTION("MGGA-LAPL + MGGA-LAPL Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

  SECTION("MGGA-LAPL + MGGA-LAPL Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

  SECTION("MGGA-LAPL + LDA Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("MGGA-LAPL + LDA Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::R2SCANL_X,
      dist(gen), Kernel::VWN5 
    );
  }

  SECTION("MGGA-TAU Polarized: EXC Only") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SCAN_X 
    );
  }

  SECTION("MGGA-TAU Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SCAN_X 
    );
  }

  SECTION("MGGA-TAU + MGGA-TAU Polarized: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::SCAN_C 
    );
  }

  SECTION("MGGA-TAU + MGGA-TAU Polarized: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SCAN_X,
      dist(gen), Kernel::SCAN_C 
    );
  }







  SECTION("LDA + MGGA-LAPL: EXC") {
    check_correctness( TestInterface::EXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::R2SCANL_C 
    );
  }
  SECTION("LDA + MGGA-LAPL: EXC + VXC") {
    check_correctness( TestInterface::EXC_VXC, Backend::libxc, Spin::Polarized, 
      dist(gen), Kernel::SlaterExchange,
      dist(gen), Kernel::R2SCANL_C 
    );
  }

}

TEST_CASE( "functional_map Test", "[xc-functional-map]") {

  SECTION("Conversion of String to Functional") {

    for (auto pair : string_functional_pairs) {
      CHECK(functional_map.value(pair.first) == pair.second);
    }

  }

  SECTION("Conversion of Functional to String") {

    for (auto pair : string_functional_pairs) {
      CHECK(functional_map.key(pair.second) == pair.first);
    }

  }

}
