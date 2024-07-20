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

TEST_CASE( "XCKernel Metadata Validity", "[xc-kernel]" ) {

  const int npts = 1024;

  auto lda_kernel_test = Kernel::SlaterExchange;
  auto gga_kernel_test = Kernel::LYP;
  auto hyb_kernel_test = Kernel::B3LYP;

  auto mgga_tau_kernel_test  = Kernel::SCAN_C;
  auto mgga_lapl_kernel_test = Kernel::R2SCANL_C;
  auto epc_lda_kernel_test   = Kernel::EPC17_2;

  Backend backend;

  SECTION( "Pure LDA Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel lda( backend, lda_kernel_test, Spin::Unpolarized );

    CHECK( lda.is_lda() );
    CHECK( not lda.is_polarized() );
    CHECK( not lda.is_gga() );
    CHECK( not lda.is_mgga() );
    CHECK( not lda.is_hyb() );
    CHECK( not lda.needs_laplacian() );

    CHECK( lda.rho_buffer_len( npts )    == npts );
    CHECK( lda.sigma_buffer_len( npts )  == 0    );
    CHECK( lda.lapl_buffer_len( npts )   == 0    );
    CHECK( lda.tau_buffer_len( npts )    == 0    );
    CHECK( lda.exc_buffer_len( npts )    == npts );
    CHECK( lda.vrho_buffer_len( npts )   == npts );
    CHECK( lda.vsigma_buffer_len( npts ) == 0    );
    CHECK( lda.vlapl_buffer_len( npts )  == 0    );
    CHECK( lda.vtau_buffer_len( npts )   == 0    );

  }

  SECTION( "Pure LDA Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel lda( backend, lda_kernel_test, Spin::Polarized );

    CHECK( lda.is_lda() );
    CHECK( lda.is_polarized() );
    CHECK( not lda.is_gga() );
    CHECK( not lda.is_mgga() );
    CHECK( not lda.is_hyb() );
    CHECK( not lda.needs_laplacian() );

    CHECK( lda.rho_buffer_len( npts )    == 2*npts );
    CHECK( lda.sigma_buffer_len( npts )  == 0      );
    CHECK( lda.lapl_buffer_len( npts )   == 0      );
    CHECK( lda.tau_buffer_len( npts )    == 0      );
    CHECK( lda.exc_buffer_len( npts )    == npts   );
    CHECK( lda.vrho_buffer_len( npts )   == 2*npts );
    CHECK( lda.vsigma_buffer_len( npts ) == 0      );
    CHECK( lda.vlapl_buffer_len( npts )  == 0      );
    CHECK( lda.vtau_buffer_len( npts )   == 0      );

  }


  SECTION( "Pure GGA Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel gga( backend, gga_kernel_test, Spin::Unpolarized );

    CHECK( gga.is_gga() );
    CHECK( not gga.is_polarized() );
    CHECK( not gga.is_lda() );
    CHECK( not gga.is_mgga() );
    CHECK( not gga.is_hyb() );
    CHECK( not gga.needs_laplacian() );

    CHECK( gga.rho_buffer_len( npts )    == npts );
    CHECK( gga.sigma_buffer_len( npts )  == npts );
    CHECK( gga.lapl_buffer_len( npts )   == 0    );
    CHECK( gga.tau_buffer_len( npts )    == 0    );
    CHECK( gga.exc_buffer_len( npts )    == npts );
    CHECK( gga.vrho_buffer_len( npts )   == npts );
    CHECK( gga.vsigma_buffer_len( npts ) == npts );
    CHECK( gga.vlapl_buffer_len( npts )  == 0    );
    CHECK( gga.vtau_buffer_len( npts )   == 0    );

  }

  SECTION( "Pure GGA Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel gga( backend, gga_kernel_test, Spin::Polarized );

    CHECK( gga.is_gga() );
    CHECK( gga.is_polarized() );
    CHECK( not gga.is_lda() );
    CHECK( not gga.is_mgga() );
    CHECK( not gga.is_hyb() );
    CHECK( not gga.needs_laplacian() );

    CHECK( gga.rho_buffer_len( npts )    == 2*npts );
    CHECK( gga.sigma_buffer_len( npts )  == 3*npts );
    CHECK( gga.lapl_buffer_len( npts )   == 0    );
    CHECK( gga.tau_buffer_len( npts )    == 0    );
    CHECK( gga.exc_buffer_len( npts )    == npts   );
    CHECK( gga.vrho_buffer_len( npts )   == 2*npts );
    CHECK( gga.vsigma_buffer_len( npts ) == 3*npts );
    CHECK( gga.vlapl_buffer_len( npts )  == 0    );
    CHECK( gga.vtau_buffer_len( npts )   == 0    );

  }

  SECTION( "Pure MGGA-TAU Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_tau_kernel_test, Spin::Unpolarized );

    CHECK( mgga.is_mgga() );
    CHECK( not mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
    CHECK( not mgga.is_hyb() );
    CHECK( not mgga.needs_laplacian() );

    CHECK( mgga.rho_buffer_len( npts )    == npts );
    CHECK( mgga.sigma_buffer_len( npts )  == npts );
    CHECK( mgga.lapl_buffer_len( npts )   == 0    );
    CHECK( mgga.tau_buffer_len( npts )    == npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts );
    CHECK( mgga.vrho_buffer_len( npts )   == npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == 0    );
    CHECK( mgga.vtau_buffer_len( npts )   == npts );

  }

  SECTION( "Pure MGGA-LAPL Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    //SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_lapl_kernel_test, Spin::Unpolarized );

    CHECK( mgga.is_mgga() );
    CHECK( mgga.needs_laplacian() );
    CHECK( not mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
    CHECK( not mgga.is_hyb() );

    CHECK( mgga.rho_buffer_len( npts )    == npts );
    CHECK( mgga.sigma_buffer_len( npts )  == npts );
    CHECK( mgga.lapl_buffer_len( npts )   == npts );
    CHECK( mgga.tau_buffer_len( npts )    == npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts );
    CHECK( mgga.vrho_buffer_len( npts )   == npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == npts );
    CHECK( mgga.vtau_buffer_len( npts )   == npts );

  }

  SECTION( "Pure MGGA-TAU Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_tau_kernel_test, Spin::Polarized );

    CHECK( mgga.is_mgga() );
    CHECK( mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
    CHECK( not mgga.is_hyb() );
    CHECK( not mgga.needs_laplacian() );

    CHECK( mgga.rho_buffer_len( npts )    == 2*npts );
    CHECK( mgga.sigma_buffer_len( npts )  == 3*npts );
    CHECK( mgga.lapl_buffer_len( npts )   == 0      );
    CHECK( mgga.tau_buffer_len( npts )    == 2*npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts   );
    CHECK( mgga.vrho_buffer_len( npts )   == 2*npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == 3*npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == 0      );
    CHECK( mgga.vtau_buffer_len( npts )   == 2*npts );

  }

  SECTION( "Pure MGGA-LAPL Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    //SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_lapl_kernel_test, Spin::Unpolarized );

    CHECK( mgga.is_mgga() );
    CHECK( mgga.needs_laplacian() );
    CHECK( not mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
    CHECK( not mgga.is_hyb() );

    CHECK( mgga.rho_buffer_len( npts )    == npts );
    CHECK( mgga.sigma_buffer_len( npts )  == npts );
    CHECK( mgga.lapl_buffer_len( npts )   == npts );
    CHECK( mgga.tau_buffer_len( npts )    == npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts );
    CHECK( mgga.vrho_buffer_len( npts )   == npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == npts );
    CHECK( mgga.vtau_buffer_len( npts )   == npts );

  }

  SECTION( "Pure MGGA-LAPL Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    //SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_lapl_kernel_test, Spin::Polarized );

    CHECK( mgga.is_mgga() );
    CHECK( mgga.needs_laplacian() );
    CHECK( mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
    CHECK( not mgga.is_hyb() );

    CHECK( mgga.rho_buffer_len( npts )    == 2*npts );
    CHECK( mgga.sigma_buffer_len( npts )  == 3*npts );
    CHECK( mgga.lapl_buffer_len( npts )   == 2*npts );
    CHECK( mgga.tau_buffer_len( npts )    == 2*npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts   );
    CHECK( mgga.vrho_buffer_len( npts )   == 2*npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == 3*npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == 2*npts );
    CHECK( mgga.vtau_buffer_len( npts )   == 2*npts );

  }

  SECTION( "Hybrid" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel hyb( backend, hyb_kernel_test, Spin::Unpolarized );
    CHECK( hyb.is_hyb() );

  }

  SECTION( "EPC LDA Polarized" ) {

    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel lda( backend, epc_lda_kernel_test, Spin::Polarized );

    CHECK( lda.is_lda() );
    CHECK( lda.is_epc() );
    CHECK( lda.is_polarized() );
    CHECK( not lda.is_gga() );
    CHECK( not lda.is_mgga() );
    CHECK( not lda.is_hyb() );
    CHECK( not lda.needs_laplacian() );

    CHECK( lda.rho_buffer_len( npts )    == 2*npts );
    CHECK( lda.sigma_buffer_len( npts )  == 0      );
    CHECK( lda.lapl_buffer_len( npts )   == 0      );
    CHECK( lda.tau_buffer_len( npts )    == 0      );
    CHECK( lda.exc_buffer_len( npts )    == npts   );
    CHECK( lda.vrho_buffer_len( npts )   == 2*npts );
    CHECK( lda.vsigma_buffer_len( npts ) == 0      );
    CHECK( lda.vlapl_buffer_len( npts )  == 0      );
    CHECK( lda.vtau_buffer_len( npts )   == 0      );

  }

}

TEST_CASE( "XCKernel Metadata Correctness", "[xc-kernel]" ) {

  Backend backend;
  SECTION( "LDA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : lda_kernels ) {
      XCKernel func( backend, kern, Spin::Unpolarized );
      auto exx = load_reference_exx( kern );

      CHECK( func.is_lda() );
      CHECK( exx == Approx( func.hyb_exx() ) );

      if( std::abs(exx) > 0 ) CHECK( func.is_hyb() );
    }

  }

  SECTION( "GGA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : gga_kernels ) {
      XCKernel func( backend, kern, Spin::Unpolarized );
      auto exx = load_reference_exx( kern );

      CHECK( func.is_gga() );
      CHECK( exx == Approx( func.hyb_exx() ) );

      if( std::abs(exx) > 0 ) CHECK( func.is_hyb() );
    }

  }

  SECTION( "MGGA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : mgga_kernels ) {
      if ( backend == Backend::builtin && kern == ExchCXX::Kernel::R2SCANL_X ) continue;
      if ( backend == Backend::builtin && kern == ExchCXX::Kernel::R2SCANL_C ) continue;
      XCKernel func( backend, kern, Spin::Unpolarized );
      auto exx = load_reference_exx( kern );

      CHECK( func.is_mgga() );
      CHECK( exx == Approx( func.hyb_exx() ) );

      if( std::abs(exx) > 0 ) CHECK( func.is_hyb() );
    }

  }

  SECTION( "EPC LDA Kernels" ) {

    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : epc_lda_kernels ) {
      XCKernel func( backend, kern, Spin::Polarized );
      auto exx = load_reference_exx( kern );

      CHECK( func.is_lda() );
      CHECK( func.is_epc() );
      CHECK( exx == Approx( func.hyb_exx() ) );

      if( std::abs(exx) > 0 ) CHECK( func.is_hyb() );
    }

  }
}









void kernel_test( TestInterface interface, Backend backend, Kernel kern,
  Spin polar ) {

  // Create the kernel
  XCKernel func( backend, kern, polar );

  const double alpha = 3.14;
  const double fill_val_e = 2.;
  const double fill_val_vr = 10.;
  const double fill_val_vs = 50.;
  const double fill_val_vl = 3.;
  const double fill_val_vt = 5.;

  const bool use_ref_values =
    (interface != TestInterface::EXC_INC) and
    (interface != TestInterface::EXC_VXC_INC);

  // LDA XC Kernels
  if( func.is_lda() ) {

    // Get reference values
    auto ref_vals = use_ref_values ?
      load_lda_reference_values( kern, polar ) :
      gen_lda_reference_values( backend,kern, polar );

    size_t  npts     = ref_vals.npts;
    const auto&   rho      = ref_vals.rho;
    const auto&   exc_ref  = ref_vals.exc;
    const auto&   vrho_ref = ref_vals.vrho;


    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ), fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ), fill_val_vr );

    if( interface == TestInterface::EXC ) {

      // Evaluate XC kernel
      func.eval_exc( npts, rho.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(exc_ref[i]) );

    }

    if( interface == TestInterface::EXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_inc( alpha, npts, rho.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );

    }

    if( interface == TestInterface::EXC_VXC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc( npts, rho.data(), exc.data(), vrho.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(exc_ref[i]) );
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(vrho_ref[i]) );

    }

    if( interface == TestInterface::EXC_VXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc_inc( alpha, npts, rho.data(), exc.data(), vrho.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );

    }

  // GGA XC Kernels
  } else if( func.is_gga() ) {

    // Get reference values
    auto ref_vals = use_ref_values ?
      load_gga_reference_values( kern, polar ) :
      gen_gga_reference_values( backend,kern, polar );

    size_t  npts       = ref_vals.npts;
    const auto&   sigma      = ref_vals.sigma;
    const auto&   rho        = ref_vals.rho;
    const auto&   exc_ref    = ref_vals.exc;
    const auto&   vrho_ref   = ref_vals.vrho;
    const auto&   vsigma_ref = ref_vals.vsigma;

    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ), fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ), fill_val_vr );
    std::vector<double> vsigma( func.vsigma_buffer_len( npts ), fill_val_vs );

    if( interface == TestInterface::EXC ) {

      // Evaluate XC kernel
      func.eval_exc( npts, rho.data(), sigma.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(exc_ref[i]) );

    }

    if( interface == TestInterface::EXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_inc( alpha, npts, rho.data(), sigma.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );

    }

    if( interface == TestInterface::EXC_VXC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc.data(),
        vrho.data(), vsigma.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(exc_ref[i]) );
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(vrho_ref[i]) );
      for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i )
        CHECK( vsigma[i] == Approx(vsigma_ref[i]) );

    }

    if( interface == TestInterface::EXC_VXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc_inc( alpha, npts, rho.data(), sigma.data(), exc.data(),
        vrho.data(), vsigma.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
      for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i )
        CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );

    }

  } else if( func.is_mgga() ) {

    // Get reference values
    auto ref_vals = use_ref_values ?
      load_mgga_reference_values( kern, polar, func.needs_laplacian() ) :
      gen_mgga_reference_values( backend,kern, polar );
    //auto ref_vals = 
    //  gen_mgga_reference_values( backend,kern, polar );

    size_t  npts       = ref_vals.npts;
    const auto&   rho        = ref_vals.rho;
    const auto&   sigma      = ref_vals.sigma;
    const auto&   lapl       = ref_vals.lapl;
    const auto&   tau        = ref_vals.tau;
    const auto&   exc_ref    = ref_vals.exc;
    const auto&   vrho_ref   = ref_vals.vrho;
    const auto&   vsigma_ref = ref_vals.vsigma;
    const auto&   vlapl_ref  = ref_vals.vlapl;
    const auto&   vtau_ref   = ref_vals.vtau;

    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ),       fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ),     fill_val_vr );
    std::vector<double> vsigma( func.vsigma_buffer_len( npts ), fill_val_vs );
    std::vector<double> vlapl( func.vlapl_buffer_len( npts ),   fill_val_vl );
    std::vector<double> vtau( func.vtau_buffer_len( npts ),     fill_val_vt );

    if( interface == TestInterface::EXC ) {

      // Evaluate XC kernel
      func.eval_exc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i ) {
        CHECK( exc[i] == Approx(exc_ref[i]) );
      }

    }

    if( interface == TestInterface::EXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_inc( alpha, npts, rho.data(), sigma.data(), lapl.data(), tau.data(), exc.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );

    }

    if( interface == TestInterface::EXC_VXC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(), exc.data(),
        vrho.data(), vsigma.data(), vlapl.data(), vtau.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i ) {
        CHECK( exc[i] == Approx(exc_ref[i]) );
      }
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i ) {
        CHECK( vrho[i] == Approx(vrho_ref[i]) );
      }
      for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i ) {
        CHECK( vsigma[i] == Approx(vsigma_ref[i]) );
      }
      for( auto i = 0ul; i < func.vlapl_buffer_len(npts); ++i ) {
        CHECK( vlapl[i] == Approx(vlapl_ref[i]) );
      }
      for( auto i = 0ul; i < func.vtau_buffer_len(npts); ++i ) {
#if XC_MAJOR_VERSION > 6
	if (func.needs_tau() )
#endif
          CHECK( vtau[i] == Approx(vtau_ref[i]) );
      }

    }

    if( interface == TestInterface::EXC_VXC_INC ) {

      // Evaluate XC kernel
      func.eval_exc_vxc_inc( alpha, npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        exc.data(), vrho.data(), vsigma.data(), vlapl.data(), vtau.data() );

      // Check correctness
      for( auto i = 0ul; i < func.exc_buffer_len(npts); ++i )
        CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
      for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i )
        CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
      for( auto i = 0ul; i < func.vlapl_buffer_len(npts); ++i )
        CHECK( vlapl[i] == Approx(fill_val_vl + alpha * vlapl_ref[i]) );
      for( auto i = 0ul; i < func.vtau_buffer_len(npts); ++i )
        CHECK( vtau[i] == Approx(fill_val_vt + alpha * vtau_ref[i]) );

    }

  }

}

TEST_CASE( "Libxc Correctness Check", "[xc-libxc]" ) {

  SECTION( "SlaterExchange Unpolarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::SlaterExchange,
      Spin::Unpolarized );
  }
  SECTION( "SlaterExchange Unpolarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::SlaterExchange,
      Spin::Unpolarized );
  }

  SECTION( "SlaterExchange Polarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::SlaterExchange,
      Spin::Polarized );
  }
  SECTION( "SlaterExchange Polarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::SlaterExchange,
      Spin::Polarized );
  }

  SECTION( "LYP Unpolarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::LYP,
      Spin::Unpolarized );
  }
  SECTION( "LYP Unpolarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::LYP,
      Spin::Unpolarized );
  }

  SECTION( "LYP Polarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::LYP,
      Spin::Polarized );
  }
  SECTION( "LYP Polarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::LYP,
      Spin::Polarized );
  }

  SECTION( "SCAN Unpolarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::SCAN_X,
      Spin::Unpolarized );
  }
  SECTION( "SCAN Unpolarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::SCAN_X,
      Spin::Unpolarized );
  }

  SECTION( "SCAN Polarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::SCAN_X,
      Spin::Polarized );
  }
  SECTION( "SCAN Polarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::SCAN_X,
      Spin::Polarized );
  }
  

  SECTION( "R2SCANL Unpolarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::R2SCANL_X,
      Spin::Unpolarized );
  }
  SECTION( "R2SCANL Unpolarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::R2SCANL_X,
      Spin::Unpolarized );
  }

  SECTION( "R2SCANL Polarized: EXC" ) {
    kernel_test( TestInterface::EXC, Backend::libxc, Kernel::R2SCANL_X,
      Spin::Polarized );
  }
  SECTION( "R2SCANL Polarized: EXC + VXC" ) {
    kernel_test( TestInterface::EXC_VXC, Backend::libxc, Kernel::R2SCANL_X,
      Spin::Polarized );
  }

}


void compare_libxc_builtin( TestInterface interface, EvalType evaltype,
  Kernel kern, Spin polar ) {

  size_t npts_lda, npts_gga, npts_mgga, npts_lapl;
  std::vector<double> ref_rho, ref_sigma, ref_lapl, ref_tau;
  std::tie(npts_lda, ref_rho  )  = load_reference_density( polar );
  std::tie(npts_gga, ref_sigma)  = load_reference_sigma  ( polar );
  std::tie(npts_lapl, ref_lapl)  = load_reference_lapl   ( polar );
  std::tie(npts_mgga, ref_tau)   = load_reference_tau    ( polar );

  REQUIRE( npts_lda == npts_gga );
  REQUIRE( npts_lda == npts_mgga );
  REQUIRE( npts_lda == npts_lapl );

  const int npts = npts_lda;

  XCKernel func_libxc  ( Backend::libxc,   kern, polar );
  XCKernel func_builtin( Backend::builtin, kern, polar );


  const int len_rho   = func_libxc.rho_buffer_len( npts );
  const int len_sigma = func_libxc.sigma_buffer_len( npts );
  const int len_lapl = func_libxc.lapl_buffer_len( npts );
  const int len_tau = func_libxc.tau_buffer_len( npts );

  std::vector<double> rho_small(len_rho, 1e-13);
  std::vector<double> sigma_small(len_sigma, 1e-14);
  std::vector<double> lapl_small(len_lapl, 1e-14);
  std::vector<double> tau_small(len_tau, 1e-14);

  std::vector<double> rho_zero(len_rho, 0.);
  std::vector<double> sigma_zero(len_sigma, 0.);
  std::vector<double> lapl_zero(len_lapl, 0.);
  std::vector<double> tau_zero(len_tau, 0.);

  std::vector<double> rho_use, sigma_use, lapl_use, tau_use;

  if( evaltype == EvalType::Regular ) {
    rho_use   = ref_rho;
    sigma_use = ref_sigma;
    lapl_use = ref_lapl;
    tau_use = ref_tau;
  }

  if( evaltype == EvalType::Small ) {
    rho_use   = rho_small;
    sigma_use = sigma_small;
    lapl_use = lapl_small;
    tau_use = tau_small;
  }

  if( evaltype == EvalType::Zero ) {
    rho_use   = rho_zero;
    sigma_use = sigma_zero;
    lapl_use = lapl_zero;
    tau_use = tau_zero;
  }



  std::vector<double> exc_libxc( func_builtin.exc_buffer_len(npts) );
  std::vector<double> vrho_libxc( func_builtin.vrho_buffer_len(npts) );
  std::vector<double> vsigma_libxc( func_builtin.vsigma_buffer_len(npts) );
  std::vector<double> vlapl_libxc( func_builtin.vlapl_buffer_len(npts) );
  std::vector<double> vtau_libxc( func_builtin.vtau_buffer_len(npts) );

  std::vector<double> exc_builtin( func_builtin.exc_buffer_len(npts) );
  std::vector<double> vrho_builtin( func_builtin.vrho_buffer_len(npts) );
  std::vector<double> vsigma_builtin( func_builtin.vsigma_buffer_len(npts) );
  std::vector<double> vlapl_builtin( func_builtin.vlapl_buffer_len(npts) );
  std::vector<double> vtau_builtin( func_builtin.vtau_buffer_len(npts) );

  if( func_libxc.is_lda() ) {

    if( interface == TestInterface::EXC ) {

      func_libxc.eval_exc( npts, rho_use.data(), exc_libxc.data() );
      func_builtin.eval_exc( npts, rho_use.data(), exc_builtin.data() );

    } else if( interface == TestInterface::EXC_VXC ) {

      func_libxc.eval_exc_vxc( npts, rho_use.data(), exc_libxc.data(),
        vrho_libxc.data() );
      func_builtin.eval_exc_vxc( npts, rho_use.data(), exc_builtin.data(),
        vrho_builtin.data() );

    }

  } else if( func_libxc.is_gga() ) {

    if( interface == TestInterface::EXC ) {

      func_libxc.eval_exc( npts, rho_use.data(), sigma_use.data(),
        exc_libxc.data() );
      func_builtin.eval_exc( npts, rho_use.data(), sigma_use.data(),
        exc_builtin.data() );

    } else if( interface == TestInterface::EXC_VXC ) {

      func_libxc.eval_exc_vxc( npts, rho_use.data(), sigma_use.data(),
        exc_libxc.data(), vrho_libxc.data(), vsigma_libxc.data() );
      func_builtin.eval_exc_vxc( npts, rho_use.data(), sigma_use.data(),
        exc_builtin.data(), vrho_builtin.data(), vsigma_builtin.data() );

    }

  } else if( func_libxc.is_mgga() ) {

    if( interface == TestInterface::EXC ) {

      func_libxc.eval_exc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), exc_libxc.data() );
      func_builtin.eval_exc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), exc_builtin.data() );

    } else if( interface == TestInterface::EXC_VXC ) {

      func_libxc.eval_exc_vxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), exc_libxc.data(), vrho_libxc.data(), vsigma_libxc.data(), vlapl_libxc.data(), vtau_libxc.data() );
      func_builtin.eval_exc_vxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), exc_builtin.data(), vrho_builtin.data(), vsigma_builtin.data(), vlapl_builtin.data(), vtau_builtin.data() );

    }

  }

  // Check correctness
  for( auto i = 0ul; i < func_libxc.exc_buffer_len(npts); ++i ) {
    INFO( "EXC Fails: Kernel is " << kern );
    CHECK( exc_builtin[i] == Approx(exc_libxc[i]) );
  }

  if( interface == TestInterface::EXC_VXC ) {
    for( auto i = 0ul; i < func_libxc.vrho_buffer_len(npts); ++i ) {
      INFO( "VRHO Fails: Kernel is " << kern );
      CHECK( vrho_builtin[i] == Approx(vrho_libxc[i]) );
    }
    for( auto i = 0ul; i < func_libxc.vsigma_buffer_len(npts); ++i ) {
      INFO( "VSIGMA Fails: Kernel is " << kern );
      CHECK( vsigma_builtin[i] == Approx(vsigma_libxc[i]) );
    }
  }

}

TEST_CASE( "Builtin Corectness Test", "[xc-builtin]" ) {

  SECTION( "Unpolarized Regular Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Regular,
        kern, Spin::Unpolarized );
    }    
  }

  SECTION( "Unpolarized Regular Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Regular,
        kern, Spin::Unpolarized );
    }    
  }

  SECTION( "Unpolarized Small Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Small,
        kern, Spin::Unpolarized );
    }
  }

  SECTION( "Unpolarized Small Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Small,
        kern, Spin::Unpolarized );
    }
  }

  SECTION( "Unpolarized Zero Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Zero,
        kern, Spin::Unpolarized );
    }    
  }

  SECTION( "Unpolarized Zero Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Zero,
        kern, Spin::Unpolarized );
    }    
  }






  SECTION( "Polarized Regular Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Regular,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Regular Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Regular,
        kern, Spin::Polarized );
    }    
  }

  SECTION( "Polarized Small Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Small,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Small Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Small,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Zero Eval : EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC, EvalType::Zero,
        kern, Spin::Polarized );
    }    
  }

  SECTION( "Polarized Zero Eval : EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::EXC_VXC, EvalType::Zero,
        kern, Spin::Polarized );
    }
  }

}


TEST_CASE( "Scale and Increment Interface", "[xc-inc]" ) {

  SECTION( "Builtin Unpolarized EXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      kernel_test( TestInterface::EXC_INC, Backend::builtin, kern,
        Spin::Unpolarized );
    }
  }

  SECTION( "Builtin Unpolarized EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_epc(kern)) continue;
      kernel_test( TestInterface::EXC_VXC_INC, Backend::builtin, kern,
        Spin::Unpolarized );
    }
  }

  SECTION( "Builtin Polarized EXC" ) {
    for( auto kern : builtin_supported_kernels )
      kernel_test( TestInterface::EXC_INC, Backend::builtin, kern,
        Spin::Polarized );
  }

  SECTION( "Builtin Polarized EXC + VXC" ) {
    for( auto kern : builtin_supported_kernels ) 
      kernel_test( TestInterface::EXC_VXC_INC, Backend::builtin, kern,
        Spin::Polarized );
  }
}

TEST_CASE( "kernel_map Test", "[xc-kernel-map]") {

  SECTION("Conversion of String to Kernel") {

    for (auto pair : string_kernal_pairs) {
      CHECK(kernel_map.value(pair.first) == pair.second);
    }

  }

  SECTION("Conversion of Kernel to String") {

    for (auto pair : string_kernal_pairs) {
      CHECK(kernel_map.key(pair.second) == pair.first);
    }

  }

}

#ifdef EXCHCXX_ENABLE_CUDA

template <typename T>
T* safe_cuda_malloc( size_t n ) {

  T* ptr = nullptr;;
  if( n ) {
    auto stat = cudaMalloc( (void**)&ptr, n*sizeof(T) );
    if( stat != cudaSuccess )
      throw std::runtime_error(cudaGetErrorString( stat ));
  }
  return ptr;

}

template <typename T>
void safe_cuda_cpy( T* dest, const T* src, size_t len ) {

  auto stat = cudaMemcpy( dest, src, len*sizeof(T), cudaMemcpyDefault );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

}

void cuda_free_all(){ }
template <typename T, typename... Args>
void cuda_free_all( T* ptr, Args&&... args ) {

  if( ptr ) {
    auto stat = cudaFree( ptr );
    if( stat != cudaSuccess )
      throw std::runtime_error(cudaGetErrorString( stat ));
  }

  cuda_free_all( std::forward<Args>(args)... );


}

void device_synchronize() {
  auto stat = cudaDeviceSynchronize();
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
}


void test_cuda_interface( TestInterface interface, EvalType evaltype,
  Backend backend, Kernel kern, Spin polar ) {

  size_t npts_lda, npts_gga, npts_mgga, npts_lapl;
  std::vector<double> ref_rho, ref_sigma, ref_lapl, ref_tau;
  std::tie(npts_lda, ref_rho  )  = load_reference_density( polar );
  std::tie(npts_gga, ref_sigma)  = load_reference_sigma  ( polar );
  std::tie(npts_lapl, ref_lapl)  = load_reference_lapl   ( polar );
  std::tie(npts_mgga, ref_tau)   = load_reference_tau    ( polar );

  REQUIRE( npts_lda == npts_gga );
  REQUIRE( npts_lda == npts_mgga );
  REQUIRE( npts_lda == npts_lapl );

  const int npts = npts_lda;

  XCKernel func( backend, kern, polar );

  size_t len_rho_buffer    = func.rho_buffer_len(npts);
  size_t len_sigma_buffer  = func.sigma_buffer_len(npts);
  size_t len_lapl_buffer   = func.lapl_buffer_len(npts);
  size_t len_tau_buffer    = func.tau_buffer_len(npts);
  size_t len_exc_buffer    = func.exc_buffer_len(npts);
  size_t len_vrho_buffer   = func.vrho_buffer_len(npts);
  size_t len_vsigma_buffer = func.vsigma_buffer_len(npts);
  size_t len_vlapl_buffer  = func.vlapl_buffer_len(npts);
  size_t len_vtau_buffer   = func.vtau_buffer_len(npts);


  std::vector<double> rho_small(len_rho_buffer, 1e-13);
  std::vector<double> sigma_small(len_sigma_buffer, 1e-14);
  std::vector<double> lapl_small(len_lapl_buffer, 1e-14);
  std::vector<double> tau_small(len_tau_buffer, 1e-14);

  std::vector<double> rho_zero(len_rho_buffer, 0.);
  std::vector<double> sigma_zero(len_sigma_buffer, 0.);
  std::vector<double> lapl_zero(len_lapl_buffer, 0.);
  std::vector<double> tau_zero(len_tau_buffer, 0.);

  std::vector<double> rho, sigma, lapl, tau;

  if( evaltype == EvalType::Regular ) {
    rho   = ref_rho;
    sigma = ref_sigma;
    lapl  = ref_lapl;
    tau   = ref_tau;
  }

  if( evaltype == EvalType::Small ) {
    rho   = rho_small;
    sigma = sigma_small;
    lapl  = lapl_small;
    tau   = tau_small;
  }

  if( evaltype == EvalType::Zero ) {
    rho   = rho_zero;
    sigma = sigma_zero;
    lapl  = lapl_zero;
    tau   = tau_zero;
  }

  // Get Reference Values
  std::vector<double>
    exc_ref( len_exc_buffer ),
    vrho_ref( len_vrho_buffer ),
    vsigma_ref( len_vsigma_buffer ),
    vlapl_ref( len_vlapl_buffer ),
    vtau_ref( len_vtau_buffer );

  if( interface == TestInterface::EXC or interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc( npts, rho.data(), exc_ref.data() );
    else if( func.is_gga() )
      func.eval_exc( npts, rho.data(), sigma.data(), exc_ref.data() );
    else if( func.is_mgga() )
      func.eval_exc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(), exc_ref.data() );

  } else if( interface == TestInterface::EXC_VXC or interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc( npts, rho.data(), exc_ref.data(), vrho_ref.data() );
    else if( func.is_gga() )
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc_ref.data(),
        vrho_ref.data(), vsigma_ref.data() );
    else if( func.is_mgga() )
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        exc_ref.data(), vrho_ref.data(), vsigma_ref.data(), vlapl_ref.data(), vtau_ref.data() );

  }






  // Allocate device memory
  double* rho_device    = safe_cuda_malloc<double>( len_rho_buffer    );
  double* sigma_device  = safe_cuda_malloc<double>( len_sigma_buffer  );
  double* lapl_device   = safe_cuda_malloc<double>( len_lapl_buffer  );
  double* tau_device    = safe_cuda_malloc<double>( len_tau_buffer  );
  double* exc_device    = safe_cuda_malloc<double>( len_exc_buffer    );
  double* vrho_device   = safe_cuda_malloc<double>( len_vrho_buffer   );
  double* vsigma_device = safe_cuda_malloc<double>( len_vsigma_buffer );
  double* vlapl_device  = safe_cuda_malloc<double>( len_vlapl_buffer );
  double* vtau_device   = safe_cuda_malloc<double>( len_vtau_buffer );

  // H2D Copy of rho / sigma
  safe_cuda_cpy( rho_device, rho.data(), len_rho_buffer );
  if( func.is_gga() or func.is_mgga() )
    safe_cuda_cpy( sigma_device, sigma.data(), len_sigma_buffer );
  if( func.is_mgga() )
    safe_cuda_cpy( tau_device, tau.data(), len_tau_buffer );
  if( func.needs_laplacian() )
    safe_cuda_cpy( lapl_device, lapl.data(), len_lapl_buffer );

  const double alpha = 3.14;
  const double fill_val_e = 2.;
  const double fill_val_vr = 10.;
  const double fill_val_vs = 50.;
  const double fill_val_vl = 20.;
  const double fill_val_vt = 30.;

  std::vector<double>
    exc( len_exc_buffer, fill_val_e ), vrho( len_vrho_buffer, fill_val_vr ),
    vsigma( len_vsigma_buffer, fill_val_vs ), vlapl(len_vlapl_buffer, fill_val_vl),
    vtau(len_vtau_buffer, fill_val_vt);

  // H2D copy of initial values, tests clobber / increment
  safe_cuda_cpy( exc_device, exc.data(), len_exc_buffer );
  safe_cuda_cpy( vrho_device, vrho.data(), len_vrho_buffer );
  if( func.is_gga() or func.is_mgga() )
    safe_cuda_cpy( vsigma_device, vsigma.data(), len_vsigma_buffer );
  if( func.is_mgga() )
    safe_cuda_cpy( vtau_device, vtau.data(), len_vtau_buffer );
  if( func.needs_laplacian() )
    safe_cuda_cpy( vlapl_device, vlapl.data(), len_vlapl_buffer );

  // Evaluate functional on device
  cudaStream_t stream = 0;
  if( interface == TestInterface::EXC ) {

    if( func.is_lda() )
      func.eval_exc_device( npts, rho_device, exc_device, stream );
    else if( func.is_gga() )
      func.eval_exc_device( npts, rho_device, sigma_device, exc_device,
        stream );
    else if( func.is_mgga() )
      func.eval_exc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        exc_device, stream );

  } else if( interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_inc_device( alpha, npts, rho_device, exc_device, stream );
    else if( func.is_gga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, exc_device,
        stream );
    else if( func.is_mgga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, lapl_device,
        tau_device, exc_device, stream );

  } else if( interface == TestInterface::EXC_VXC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_device( npts, rho_device, exc_device, vrho_device, stream );
    else if( func.is_gga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, exc_device,
        vrho_device, vsigma_device, stream );
    else if( func.is_mgga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        exc_device, vrho_device, vsigma_device, vlapl_device, vtau_device, stream );

  } else if( interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, exc_device,
        vrho_device, stream );
    else if( func.is_gga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        exc_device, vrho_device, vsigma_device, stream );
    else if( func.is_mgga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, exc_device, vrho_device, vsigma_device, 
        vlapl_device, vtau_device, stream );

  }

  device_synchronize();

  // D2H of results
  safe_cuda_cpy( exc.data(), exc_device, len_exc_buffer );
  safe_cuda_cpy( vrho.data(), vrho_device, len_vrho_buffer );
  if( func.is_gga() or func.is_mgga() )
    safe_cuda_cpy( vsigma.data(), vsigma_device, len_vsigma_buffer );
  if(func.is_mgga())
    safe_cuda_cpy( vtau.data(), vtau_device, len_vtau_buffer );
  if(func.needs_laplacian())
    safe_cuda_cpy( vlapl.data(), vlapl_device, len_vlapl_buffer );

  // Check correctness
  if( interface == TestInterface::EXC_INC or interface == TestInterface::EXC_VXC_INC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
  } else {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(exc_ref[i]) );
  }

  if( interface == TestInterface::EXC_VXC_INC ) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
    for( auto i = 0ul; i < len_vlapl_buffer; ++i )
      CHECK( vlapl[i] == Approx(fill_val_vl + alpha * vlapl_ref[i]) );
    for( auto i = 0ul; i < len_vtau_buffer; ++i )
      CHECK( vtau[i] == Approx(fill_val_vt + alpha * vtau_ref[i]) );

  } else if(interface == TestInterface::EXC_VXC)  {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i ) {
      INFO( "Kernel is " << kern );
      CHECK( vsigma[i] == Approx(vsigma_ref[i]) );
    }
    for( auto i = 0ul; i < len_vlapl_buffer; ++i ) {
      INFO( "Kernel is " << kern );
      CHECK( vlapl[i] == Approx(vlapl_ref[i]).margin(std::numeric_limits<double>::epsilon()) );
    }
    for( auto i = 0ul; i < len_vtau_buffer; ++i ) {
      INFO( "Kernel is " << kern << std::scientific << " " << vtau[i] << " " << vtau_ref[i] );
      CHECK( vtau[i] == Approx(vtau_ref[i]).margin(std::numeric_limits<double>::epsilon()) );
    }

  }

  cuda_free_all( rho_device, sigma_device, exc_device, vrho_device, vsigma_device, lapl_device, tau_device,
    vlapl_device, vtau_device );
}



TEST_CASE( "CUDA Interfaces", "[xc-device]" ) {

  SECTION( "Libxc Functionals" ) {

    SECTION( "LDA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }








    SECTION( "LDA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

  }

  SECTION( "Builtin Functionals" ) {

    SECTION("EXC Regular: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC, EvalType::Regular,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + VXC Regular: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_VXC, EvalType::Regular,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + INC Regular: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_INC, EvalType::Regular,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + VXC + INC Regular: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_VXC_INC, EvalType::Regular,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC Small: Unpolarized") {
      for (auto kern : builtin_supported_kernels) {
        if (is_unstable_small(kern) || !supports_unpolarized(kern))
          continue;
        test_cuda_interface(TestInterface::EXC, EvalType::Small,
                            Backend::builtin, kern, Spin::Unpolarized);
      }
    }

    SECTION("EXC + VXC Small: Unpolarized") {
      for (auto kern : builtin_supported_kernels) {
        if (is_unstable_small(kern) || !supports_unpolarized(kern))
          continue;
        test_cuda_interface(TestInterface::EXC_VXC, EvalType::Small,
                            Backend::builtin, kern, Spin::Unpolarized);
      }
    }

    SECTION("EXC + INC Small: Unpolarized") {
      for (auto kern : builtin_supported_kernels) {
        if (is_unstable_small(kern) || !supports_unpolarized(kern))
          continue;
        test_cuda_interface(TestInterface::EXC_INC, EvalType::Small,
                            Backend::builtin, kern, Spin::Unpolarized);
      }
    }

    SECTION("EXC + VXC + INC Small: Unpolarized") {
      for (auto kern : builtin_supported_kernels) {
        if (is_unstable_small(kern) || !supports_unpolarized(kern))
          continue;
        test_cuda_interface(TestInterface::EXC_VXC_INC, EvalType::Small,
                            Backend::builtin, kern, Spin::Unpolarized);
      }
    }

    SECTION("EXC Zero: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC, EvalType::Zero,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + VXC Zero: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_VXC, EvalType::Zero,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + INC Zero: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_INC, EvalType::Zero,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC + VXC + INC Zero: Unpolarized") {
      for (auto kern : builtin_supported_kernels)
        if (supports_unpolarized(kern))
          test_cuda_interface(TestInterface::EXC_VXC_INC, EvalType::Zero,
                              Backend::builtin, kern, Spin::Unpolarized);
    }

    SECTION("EXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + VXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + VXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_cuda_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_cuda_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

  }


}

#endif



#ifdef EXCHCXX_ENABLE_HIP

template <typename T>
T* safe_hip_malloc( size_t n ) {

  T* ptr = nullptr;;
  if( n ) {
    auto stat = hipMalloc( (void**)&ptr, n*sizeof(T) );
    if( stat != hipSuccess )
      throw std::runtime_error(hipGetErrorString( stat ));
  }
  return ptr;

}

template <typename T>
void safe_hip_cpy( T* dest, const T* src, size_t len ) {

  auto stat = hipMemcpy( dest, src, len*sizeof(T), hipMemcpyDefault );
  if( stat != hipSuccess )
    throw std::runtime_error(hipGetErrorString( stat ));

}

void hip_free_all(){ }
template <typename T, typename... Args>
void hip_free_all( T* ptr, Args&&... args ) {

  if( ptr ) {
    auto stat = hipFree( ptr );
    if( stat != hipSuccess )
      throw std::runtime_error(hipGetErrorString( stat ));
  }

  hip_free_all( std::forward<Args>(args)... );


}

void device_synchronize() {
  auto stat = hipDeviceSynchronize();
  if( stat != hipSuccess )
    throw std::runtime_error(hipGetErrorString( stat ));
}


void test_hip_interface( TestInterface interface, EvalType evaltype,
  Backend backend, Kernel kern, Spin polar ) {

  size_t npts_lda, npts_gga;
  std::vector<double> ref_rho, ref_sigma;
  std::tie(npts_lda, ref_rho  )  = load_reference_density( polar );
  std::tie(npts_gga, ref_sigma)  = load_reference_sigma  ( polar );

  REQUIRE( npts_lda == npts_gga );

  const int npts = npts_lda;

  XCKernel func( backend, kern, polar );

  size_t len_rho_buffer    = func.rho_buffer_len(npts);
  size_t len_sigma_buffer  = func.sigma_buffer_len(npts);
  size_t len_exc_buffer    = func.exc_buffer_len(npts);
  size_t len_vrho_buffer   = func.vrho_buffer_len(npts);
  size_t len_vsigma_buffer = func.vsigma_buffer_len(npts);


  std::vector<double> rho_small(len_rho_buffer, 1e-13);
  std::vector<double> sigma_small(len_sigma_buffer, 1e-14);

  std::vector<double> rho_zero(len_rho_buffer, 0.);
  std::vector<double> sigma_zero(len_sigma_buffer, 0.);

  std::vector<double> rho, sigma;

  if( evaltype == EvalType::Regular ) {
    rho   = ref_rho;
    sigma = ref_sigma;
  }

  if( evaltype == EvalType::Small ) {
    rho   = rho_small;
    sigma = sigma_small;
  }

  if( evaltype == EvalType::Zero ) {
    rho   = rho_zero;
    sigma = sigma_zero;
  }

  // Get Reference Values
  std::vector<double>
    exc_ref( len_exc_buffer ),
    vrho_ref( len_vrho_buffer ),
    vsigma_ref( len_vsigma_buffer );

  if( interface == TestInterface::EXC or interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc( npts, rho.data(), exc_ref.data() );
    else if( func.is_gga() )
      func.eval_exc( npts, rho.data(), sigma.data(), exc_ref.data() );

  } else if( interface == TestInterface::EXC_VXC or interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc( npts, rho.data(), exc_ref.data(), vrho_ref.data() );
    else if( func.is_gga() )
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc_ref.data(),
        vrho_ref.data(), vsigma_ref.data() );

  }






  // Allocate device memory
  double* rho_device    = safe_hip_malloc<double>( len_rho_buffer    );
  double* sigma_device  = safe_hip_malloc<double>( len_sigma_buffer  );
  double* exc_device    = safe_hip_malloc<double>( len_exc_buffer    );
  double* vrho_device   = safe_hip_malloc<double>( len_vrho_buffer   );
  double* vsigma_device = safe_hip_malloc<double>( len_vsigma_buffer );

  // H2D Copy of rho / sigma
  safe_hip_cpy( rho_device, rho.data(), len_rho_buffer );
  if( func.is_gga() )
    safe_hip_cpy( sigma_device, sigma.data(), len_sigma_buffer );

  const double alpha = 3.14;
  const double fill_val_e = 2.;
  const double fill_val_vr = 10.;
  const double fill_val_vs = 50.;

  std::vector<double>
    exc( len_exc_buffer, fill_val_e ), vrho( len_vrho_buffer, fill_val_vr ),
    vsigma( len_vsigma_buffer, fill_val_vs );

  // H2D copy of initial values, tests clobber / increment
  safe_hip_cpy( exc_device, exc.data(), len_exc_buffer );
  safe_hip_cpy( vrho_device, vrho.data(), len_vrho_buffer );
  if( func.is_gga() )
    safe_hip_cpy( vsigma_device, vsigma.data(), len_vsigma_buffer );

  // Evaluate functional on device
  hipStream_t stream = 0;
  if( interface == TestInterface::EXC ) {

    if( func.is_lda() )
      func.eval_exc_device( npts, rho_device, exc_device, stream );
    else if( func.is_gga() )
      func.eval_exc_device( npts, rho_device, sigma_device, exc_device,
        stream );

  } else if( interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_inc_device( alpha, npts, rho_device, exc_device, stream );
    else if( func.is_gga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, exc_device,
        stream );

  } else if( interface == TestInterface::EXC_VXC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_device( npts, rho_device, exc_device, vrho_device, stream );
    else if( func.is_gga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, exc_device,
        vrho_device, vsigma_device, stream );

  } else if( interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, exc_device,
        vrho_device, stream );
    else if( func.is_gga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        exc_device, vrho_device, vsigma_device, stream );

  }

  device_synchronize();

  // D2H of results
  safe_hip_cpy( exc.data(), exc_device, len_exc_buffer );
  safe_hip_cpy( vrho.data(), vrho_device, len_vrho_buffer );
  if(func.is_gga())
    safe_hip_cpy( vsigma.data(), vsigma_device, len_vsigma_buffer );

  // Check correctness
  if( interface == TestInterface::EXC_INC or interface == TestInterface::EXC_VXC_INC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
  } else {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(exc_ref[i]) );
  }

  if( interface == TestInterface::EXC_VXC_INC ) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );

  } else if(interface == TestInterface::EXC_VXC)  {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i ) {
      INFO( "Kernel is " << kern );
      CHECK( vsigma[i] == Approx(vsigma_ref[i]) );
    }

  }

  hip_free_all( rho_device, sigma_device, exc_device, vrho_device, vsigma_device );
}



TEST_CASE( "HIP Interfaces", "[xc-device]" ) {

  SECTION( "Libxc Functionals" ) {

    SECTION( "LDA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }








    SECTION( "LDA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

  }

  SECTION( "Builtin Functionals" ) {

    SECTION("EXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_hip_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }


  }


}

#endif








#ifdef EXCHCXX_ENABLE_SYCL

template <typename T>
T* safe_sycl_malloc( size_t n, sycl::queue& q ) {
  if( n ) {
    T* ptr = sycl::malloc_device<T>(n, q);
    return ptr;
  } else return nullptr;
}

template <typename T>
void safe_sycl_cpy( T* dest, const T* src, size_t len, sycl::queue& q ) {

    q.memcpy( (void*)dest, (const void*)src, len*sizeof(T) );

}

void sycl_free_all(sycl::queue&){ }
template <typename T, typename... Args>
void sycl_free_all( sycl::queue& q, T* ptr, Args&&... args ) {

  if( ptr ) {
    sycl::free( (void*)ptr, q );
  }

  sycl_free_all( q, std::forward<Args>(args)... );

}

void device_synchronize( sycl::queue& q ) {
q.wait_and_throw();
}


void test_sycl_interface( TestInterface interface, EvalType evaltype,
                          Backend backend, Kernel kern, Spin polar, sycl::queue& q ) {

  auto [npts_lda, ref_rho]   = load_reference_density( polar );
  auto [npts_gga, ref_sigma] = load_reference_sigma  ( polar );

  REQUIRE( npts_lda == npts_gga );

  const int npts = npts_lda;

  XCKernel func( backend, kern, polar );

  size_t len_rho_buffer    = func.rho_buffer_len(npts);
  size_t len_sigma_buffer  = func.sigma_buffer_len(npts);
  size_t len_exc_buffer    = func.exc_buffer_len(npts);
  size_t len_vrho_buffer   = func.vrho_buffer_len(npts);
  size_t len_vsigma_buffer = func.vsigma_buffer_len(npts);


  std::vector<double> rho_small(len_rho_buffer, 1e-13);
  std::vector<double> sigma_small(len_sigma_buffer, 1e-14);

  std::vector<double> rho_zero(len_rho_buffer, 0.);
  std::vector<double> sigma_zero(len_sigma_buffer, 0.);

  std::vector<double> rho, sigma;

  if( evaltype == EvalType::Regular ) {
    rho   = ref_rho;
    sigma = ref_sigma;
  }

  if( evaltype == EvalType::Small ) {
    rho   = rho_small;
    sigma = sigma_small;
  }

  if( evaltype == EvalType::Zero ) {
    rho   = rho_zero;
    sigma = sigma_zero;
  }

  // Get Reference Values
  std::vector<double>
    exc_ref( len_exc_buffer ),
    vrho_ref( len_vrho_buffer ),
    vsigma_ref( len_vsigma_buffer );

  if( interface == TestInterface::EXC or interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc( npts, rho.data(), exc_ref.data() );
    else if( func.is_gga() )
      func.eval_exc( npts, rho.data(), sigma.data(), exc_ref.data() );

  } else if( interface == TestInterface::EXC_VXC or interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc( npts, rho.data(), exc_ref.data(), vrho_ref.data() );
    else if( func.is_gga() )
      func.eval_exc_vxc( npts, rho.data(), sigma.data(), exc_ref.data(),
        vrho_ref.data(), vsigma_ref.data() );

  }






  // Allocate device memory
  double* rho_device    = safe_sycl_malloc<double>( len_rho_buffer   , q );
  double* sigma_device  = safe_sycl_malloc<double>( len_sigma_buffer , q );
  double* exc_device    = safe_sycl_malloc<double>( len_exc_buffer   , q );
  double* vrho_device   = safe_sycl_malloc<double>( len_vrho_buffer  , q );
  double* vsigma_device = safe_sycl_malloc<double>( len_vsigma_buffer, q );

  // H2D Copy of rho / sigma
  safe_sycl_cpy( rho_device, rho.data(), len_rho_buffer, q );
  if( func.is_gga() )
    safe_sycl_cpy( sigma_device, sigma.data(), len_sigma_buffer, q );

  const double alpha = 3.14;
  const double fill_val_e = 2.;
  const double fill_val_vr = 10.;
  const double fill_val_vs = 50.;

  std::vector<double>
    exc( len_exc_buffer, fill_val_e ), vrho( len_vrho_buffer, fill_val_vr ),
    vsigma( len_vsigma_buffer, fill_val_vs );

  // H2D copy of initial values, tests clobber / increment
  safe_sycl_cpy( exc_device, exc.data(), len_exc_buffer, q );
  safe_sycl_cpy( vrho_device, vrho.data(), len_vrho_buffer, q );
  if( func.is_gga() )
    safe_sycl_cpy( vsigma_device, vsigma.data(), len_vsigma_buffer, q );

  q.wait();

  // Evaluate functional on device
  if( interface == TestInterface::EXC ) {

    if( func.is_lda() )
      func.eval_exc_device( npts, rho_device, exc_device, &q );
    else if( func.is_gga() )
      func.eval_exc_device( npts, rho_device, sigma_device, exc_device,
        &q );

  } else if( interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_inc_device( alpha, npts, rho_device, exc_device, &q );
    else if( func.is_gga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, exc_device,
        &q );

  } else if( interface == TestInterface::EXC_VXC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_device( npts, rho_device, exc_device, vrho_device, &q );
    else if( func.is_gga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, exc_device,
        vrho_device, vsigma_device, &q );

  } else if( interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, exc_device,
        vrho_device, &q );
    else if( func.is_gga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        exc_device, vrho_device, vsigma_device, &q );

  }

  device_synchronize( q );

  // D2H of results
  safe_sycl_cpy( exc.data(), exc_device, len_exc_buffer, q );
  safe_sycl_cpy( vrho.data(), vrho_device, len_vrho_buffer, q );
  if(func.is_gga())
    safe_sycl_cpy( vsigma.data(), vsigma_device, len_vsigma_buffer, q );

  device_synchronize( q );
  // Check correctness
  if( interface == TestInterface::EXC_INC or interface == TestInterface::EXC_VXC_INC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
  } else {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(exc_ref[i]) );
  }

  if( interface == TestInterface::EXC_VXC_INC ) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );

  } else if(interface == TestInterface::EXC_VXC)  {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i ) {
      INFO( "Kernel is " << kern );
      CHECK( vsigma[i] == Approx(vsigma_ref[i]) );
    }

  }

  device_synchronize( q );
  sycl_free_all( q, rho_device, sigma_device, exc_device, vrho_device,
    	 vsigma_device );

  device_synchronize( q );
}


#if 0
struct SYCLTestFeature {
  sycl::queue q;
  SYCLTestFeature() :
    q( sycl::gpu_selector_v,
       sycl::property_list{sycl::property::queue::in_order{}} ) { }
};
#else
struct SYCLTestFeature {
  static sycl::queue q;

  SYCLTestFeature() {}
};

sycl::queue SYCLTestFeature::q(
       sycl::gpu_selector_v,
       sycl::property_list{sycl::property::queue::in_order{}} );
#endif

TEST_CASE_METHOD( SYCLTestFeature, "SYCL Interfaces", "[xc-device]" ) {

  //std::cout << "Running on "
  //          << q.get_device().get_info<sycl::info::device::name>()
  //          << "\n";

  SECTION( "Libxc Functionals" ) {

    SECTION( "LDA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "LDA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "LDA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized, q );
    }








    SECTION( "LDA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "LDA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "LDA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized, q );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized, q );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized, q );
    }

  }

  SECTION( "Builtin Functionals" ) {

    SECTION("EXC Regular: Unpolarized") {
      //std::cout << "EXC Regular: Unpolarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
                             Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC Regular: Unpolarized") {
      //std::cout << "EXC + VXC Regular: Unpolarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + INC Regular: Unpolarized") {
      //std::cout << "EXC + INC Regular: Unpolarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC + INC Regular: Unpolarized") {
      //std::cout << "EXC + VXC + INC Regular: Unpolarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC + VXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized, q );
    }

    SECTION("EXC Regular: Polarized") {
      //std::cout << "EXC Regular: Polarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Regular,
                             Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC Regular: Polarized") {
      //std::cout << "EXC + VXC Regular: Polarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + INC Regular: Polarized") {
      //std::cout << "EXC + INC Regular: Polarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC + INC Regular: Polarized") {
      //std::cout << "EXC + VXC + INC Regular: Polarized" << std::endl;
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized, q );
    }

    SECTION("EXC + VXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_sycl_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized, q );
    }

  }


}

#endif
