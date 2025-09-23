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

#include "exchcxx/enums/kernels.hpp"
#include "ut_common.hpp"

using namespace ExchCXX;

TEST_CASE( "XCKernel Metadata Validity", "[xc-kernel]" ) {

  const int npts = 1024;

  auto lda_kernel_test = Kernel::SlaterExchange;
  auto gga_kernel_test = Kernel::LYP;

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

    CHECK( lda.v2rho2_buffer_len( npts )   == npts );
    CHECK( lda.v2rhosigma_buffer_len( npts ) == 0  );
    CHECK( lda.v2rholapl_buffer_len( npts ) == 0   );
    CHECK( lda.v2rhotau_buffer_len( npts ) == 0    );
    CHECK( lda.v2sigma2_buffer_len( npts ) == 0    );
    CHECK( lda.v2sigmalapl_buffer_len( npts ) == 0 );
    CHECK( lda.v2sigmatau_buffer_len( npts ) == 0  );
    CHECK( lda.v2lapl2_buffer_len( npts ) == 0     );
    CHECK( lda.v2lapltau_buffer_len( npts ) == 0   );
    CHECK( lda.v2tau2_buffer_len( npts ) == 0      );

  }

  SECTION( "Pure LDA Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel lda( backend, lda_kernel_test, Spin::Polarized );

    CHECK( lda.is_lda() );
    CHECK( lda.is_polarized() );
    CHECK( not lda.is_gga() );
    CHECK( not lda.is_mgga() );
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

    CHECK( lda.v2rho2_buffer_len( npts )   == 3*npts );
    CHECK( lda.v2rhosigma_buffer_len( npts ) == 0    );
    CHECK( lda.v2rholapl_buffer_len( npts ) == 0     );
    CHECK( lda.v2rhotau_buffer_len( npts ) == 0      );
    CHECK( lda.v2sigma2_buffer_len( npts ) == 0      );
    CHECK( lda.v2sigmalapl_buffer_len( npts ) == 0   );
    CHECK( lda.v2sigmatau_buffer_len( npts ) == 0    );
    CHECK( lda.v2lapl2_buffer_len( npts ) == 0       );
    CHECK( lda.v2lapltau_buffer_len( npts ) == 0     );
    CHECK( lda.v2tau2_buffer_len( npts ) == 0        );
  }


  SECTION( "Pure GGA Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel gga( backend, gga_kernel_test, Spin::Unpolarized );

    CHECK( gga.is_gga() );
    CHECK( not gga.is_polarized() );
    CHECK( not gga.is_lda() );
    CHECK( not gga.is_mgga() );
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

    CHECK( gga.v2rho2_buffer_len( npts )     == npts );
    CHECK( gga.v2rhosigma_buffer_len( npts ) == npts );
    CHECK( gga.v2rholapl_buffer_len( npts )  == 0    );
    CHECK( gga.v2rhotau_buffer_len( npts )   == 0    );
    CHECK( gga.v2sigma2_buffer_len( npts )   == npts );
    CHECK( gga.v2sigmalapl_buffer_len( npts ) == 0   );
    CHECK( gga.v2sigmatau_buffer_len( npts ) == 0    );
    CHECK( gga.v2lapl2_buffer_len( npts ) == 0       );
    CHECK( gga.v2lapltau_buffer_len( npts ) == 0     );
    CHECK( gga.v2tau2_buffer_len( npts )     == 0    );
  }

  SECTION( "Pure GGA Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel gga( backend, gga_kernel_test, Spin::Polarized );

    CHECK( gga.is_gga() );
    CHECK( gga.is_polarized() );
    CHECK( not gga.is_lda() );
    CHECK( not gga.is_mgga() );
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

    CHECK( gga.v2rho2_buffer_len( npts )     == 3*npts );
    CHECK( gga.v2rhosigma_buffer_len( npts ) == 6*npts );
    CHECK( gga.v2rholapl_buffer_len( npts )  == 0      );
    CHECK( gga.v2rhotau_buffer_len( npts )   == 0      );
    CHECK( gga.v2sigma2_buffer_len( npts )   == 6*npts );
    CHECK( gga.v2sigmalapl_buffer_len( npts ) == 0     );
    CHECK( gga.v2sigmatau_buffer_len( npts ) == 0      );
    CHECK( gga.v2lapl2_buffer_len( npts ) == 0         );
    CHECK( gga.v2lapltau_buffer_len( npts ) == 0       );
    CHECK( gga.v2tau2_buffer_len( npts )     == 0      );

  }

  SECTION( "Pure MGGA-TAU Unpolarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_tau_kernel_test, Spin::Unpolarized );

    CHECK( mgga.is_mgga() );
    CHECK( not mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
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

    CHECK( mgga.v2rho2_buffer_len( npts )     == npts );
    CHECK( mgga.v2rhosigma_buffer_len( npts ) == npts );
    CHECK( mgga.v2rholapl_buffer_len( npts )  == 0    );
    CHECK( mgga.v2rhotau_buffer_len( npts )   == npts );
    CHECK( mgga.v2sigma2_buffer_len( npts )   == npts );
    CHECK( mgga.v2sigmalapl_buffer_len( npts ) == 0   );
    CHECK( mgga.v2sigmatau_buffer_len( npts ) == npts );
    CHECK( mgga.v2lapl2_buffer_len( npts ) == 0       );
    CHECK( mgga.v2lapltau_buffer_len( npts ) == 0     );
    CHECK( mgga.v2tau2_buffer_len( npts )     == npts );
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

    CHECK( mgga.rho_buffer_len( npts )    == npts );
    CHECK( mgga.sigma_buffer_len( npts )  == npts );
    CHECK( mgga.lapl_buffer_len( npts )   == npts );
    CHECK( mgga.tau_buffer_len( npts )    == npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts );
    CHECK( mgga.vrho_buffer_len( npts )   == npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == npts );
    CHECK( mgga.vtau_buffer_len( npts )   == npts );

    CHECK( mgga.v2rho2_buffer_len( npts )     == npts );
    CHECK( mgga.v2rhosigma_buffer_len( npts ) == npts );
    CHECK( mgga.v2rholapl_buffer_len( npts )  == npts );
    CHECK( mgga.v2rhotau_buffer_len( npts )   == npts );
    CHECK( mgga.v2sigma2_buffer_len( npts )   == npts );
    CHECK( mgga.v2sigmalapl_buffer_len( npts ) == npts );
    CHECK( mgga.v2sigmatau_buffer_len( npts ) == npts );
    CHECK( mgga.v2lapl2_buffer_len( npts ) == npts     );
    CHECK( mgga.v2lapltau_buffer_len( npts ) == npts   );
    CHECK( mgga.v2tau2_buffer_len( npts )     == npts );
  }

  SECTION( "Pure MGGA-TAU Polarized" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel mgga( backend, mgga_tau_kernel_test, Spin::Polarized );

    CHECK( mgga.is_mgga() );
    CHECK( mgga.is_polarized() );
    CHECK( not mgga.is_lda() );
    CHECK( not mgga.is_gga() );
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

    CHECK( mgga.v2rho2_buffer_len( npts )     == 3*npts );
    CHECK( mgga.v2rhosigma_buffer_len( npts ) == 6*npts );
    CHECK( mgga.v2rholapl_buffer_len( npts )  == 0      );
    CHECK( mgga.v2rhotau_buffer_len( npts )   == 4*npts );
    CHECK( mgga.v2sigma2_buffer_len( npts )   == 6*npts );
    CHECK( mgga.v2sigmalapl_buffer_len( npts ) == 0     );
    CHECK( mgga.v2sigmatau_buffer_len( npts ) == 6*npts );
    CHECK( mgga.v2lapl2_buffer_len( npts ) == 0         );
    CHECK( mgga.v2lapltau_buffer_len( npts ) == 0       );
    CHECK( mgga.v2tau2_buffer_len( npts )     == 3*npts );
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

    CHECK( mgga.rho_buffer_len( npts )    == 2*npts );
    CHECK( mgga.sigma_buffer_len( npts )  == 3*npts );
    CHECK( mgga.lapl_buffer_len( npts )   == 2*npts );
    CHECK( mgga.tau_buffer_len( npts )    == 2*npts );
    CHECK( mgga.exc_buffer_len( npts )    == npts   );
    CHECK( mgga.vrho_buffer_len( npts )   == 2*npts );
    CHECK( mgga.vsigma_buffer_len( npts ) == 3*npts );
    CHECK( mgga.vlapl_buffer_len( npts )  == 2*npts );
    CHECK( mgga.vtau_buffer_len( npts )   == 2*npts );

    CHECK( mgga.v2rho2_buffer_len( npts )     == 3*npts );
    CHECK( mgga.v2rhosigma_buffer_len( npts ) == 6*npts );
    CHECK( mgga.v2rholapl_buffer_len( npts )  == 4*npts );
    CHECK( mgga.v2rhotau_buffer_len( npts )   == 4*npts );
    CHECK( mgga.v2sigma2_buffer_len( npts )   == 6*npts );
    CHECK( mgga.v2sigmalapl_buffer_len( npts ) == 6*npts );
    CHECK( mgga.v2sigmatau_buffer_len( npts ) == 6*npts );
    CHECK( mgga.v2lapl2_buffer_len( npts ) == 3*npts     );
    CHECK( mgga.v2lapltau_buffer_len( npts ) == 4*npts   );
    CHECK( mgga.v2tau2_buffer_len( npts )     == 3*npts );
  }



  SECTION( "EPC LDA Polarized" ) {

    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    XCKernel lda( backend, epc_lda_kernel_test, Spin::Polarized );

    CHECK( lda.is_lda() );
    CHECK( lda.is_epc() );
    CHECK( lda.is_polarized() );
    CHECK( not lda.is_gga() );
    CHECK( not lda.is_mgga() );
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

    CHECK( lda.v2rho2_buffer_len( npts )   == 3*npts );
    CHECK( lda.v2rhosigma_buffer_len( npts ) == 0    );
    CHECK( lda.v2rholapl_buffer_len( npts ) == 0     );
    CHECK( lda.v2rhotau_buffer_len( npts ) == 0      );
    CHECK( lda.v2sigma2_buffer_len( npts ) == 0      );
    CHECK( lda.v2sigmalapl_buffer_len( npts ) == 0   );
    CHECK( lda.v2sigmatau_buffer_len( npts ) == 0    );
    CHECK( lda.v2lapl2_buffer_len( npts ) == 0       );
    CHECK( lda.v2lapltau_buffer_len( npts ) == 0     );
    CHECK( lda.v2tau2_buffer_len( npts ) == 0        );
  }

}

TEST_CASE( "XCKernel Metadata Correctness", "[xc-kernel]" ) {

  Backend backend;
  SECTION( "LDA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : lda_kernels ) {
      XCKernel func( backend, kern, Spin::Unpolarized );
      CHECK( func.is_lda() );

    }

  }

  SECTION( "GGA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : gga_kernels ) {
      XCKernel func( backend, kern, Spin::Unpolarized );
      CHECK( func.is_gga() );
    }

  }

  SECTION( "MGGA Kernels" ) {

    SECTION( "Libxc Backend" )   { backend = Backend::libxc; }
    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : mgga_kernels ) {
      if ( backend == Backend::builtin && kern == ExchCXX::Kernel::R2SCANL_X ) continue;
      if ( backend == Backend::builtin && kern == ExchCXX::Kernel::R2SCANL_C ) continue;
      XCKernel func( backend, kern, Spin::Unpolarized );
      CHECK( func.is_mgga() );
    }

  }

  SECTION( "EPC LDA Kernels" ) {

    SECTION( "Builtin Backend" ) { backend = Backend::builtin; }

    for( const auto& kern : epc_lda_kernels ) {
      XCKernel func( backend, kern, Spin::Polarized );
      CHECK( func.is_lda() );
      CHECK( func.is_epc() );
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

  const double fill_val_v2rho2 = 10.;
  const double fill_val_v2rhosigma = 11.;
  const double fill_val_v2rholapl = 12.;
  const double fill_val_v2rhotau = 13.;
  const double fill_val_v2sigma2 = 14.;
  const double fill_val_v2sigmalapl = 15.;
  const double fill_val_v2sigmatau = 16.;
  const double fill_val_v2lapl2 = 17.;
  const double fill_val_v2lapltau = 18.;
  const double fill_val_v2tau2 = 19.;

  const bool use_ref_values =
    (interface != TestInterface::EXC_INC) and
    (interface != TestInterface::EXC_VXC_INC) and
    (interface != TestInterface::FXC_INC) and
    (interface != TestInterface::VXC_FXC_INC);

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
    const auto&   v2rho2_ref = ref_vals.v2rho2;


    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ), fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ), fill_val_vr );
    std::vector<double> v2rho2( func.v2rho2_buffer_len( npts ), fill_val_v2rho2 );

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

    if( interface == TestInterface::FXC_INC ) {
      // Evaluate XC kernel
      func.eval_fxc_inc( alpha, npts, rho.data(), v2rho2.data() );
      // Check correctness
      for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
        CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );
    }

    if( interface == TestInterface::VXC_FXC_INC ) {

      // Evaluate XC kernel
      func.eval_vxc_fxc_inc( alpha, npts, rho.data(), vrho.data(), v2rho2.data() );

      // Check correctness
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
      for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
        CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );

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
    const auto&   v2rho2_ref = ref_vals.v2rho2;
    const auto&   v2rhosigma_ref = ref_vals.v2rhosigma;
    const auto&   v2sigam2_ref = ref_vals.v2sigma2;

    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ), fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ), fill_val_vr );
    std::vector<double> vsigma( func.vsigma_buffer_len( npts ), fill_val_vs );

    std::vector<double> v2rho2( func.v2rho2_buffer_len( npts ), fill_val_v2rho2 );
    std::vector<double> v2rhosigma( func.v2rhosigma_buffer_len( npts ), fill_val_v2rhosigma );
    std::vector<double> v2sigma2( func.v2sigma2_buffer_len( npts ), fill_val_v2sigma2 );

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

    if( interface == TestInterface::FXC_INC ) {

      // Evaluate XC kernel
      func.eval_fxc_inc( alpha, npts, rho.data(), sigma.data(), v2rho2.data(),
        v2rhosigma.data(), v2sigma2.data() );

      // Check correctness
      for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
        CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );
      for( auto i = 0ul; i < func.v2rhosigma_buffer_len(npts); ++i )
        CHECK( v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) );
      for( auto i = 0ul; i < func.v2sigma2_buffer_len(npts); ++i )
        CHECK( v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigam2_ref[i]) );

    }

    if( interface == TestInterface::VXC_FXC_INC ) {

      // Evaluate XC kernel
      func.eval_vxc_fxc_inc( alpha, npts, rho.data(), sigma.data(), vrho.data(),
        vsigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data() );

      // Check correctness
      for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
        CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
      for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i )
        CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
      for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
        CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );
      for( auto i = 0ul; i < func.v2rhosigma_buffer_len(npts); ++i )
        CHECK( v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) );
      for( auto i = 0ul; i < func.v2sigma2_buffer_len(npts); ++i )
        CHECK( v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigam2_ref[i]) );

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

    const auto&   v2rho2_ref = ref_vals.v2rho2;
    const auto&   v2rhosigma_ref = ref_vals.v2rhosigma;
    const auto&   v2rholapl_ref = ref_vals.v2rholapl;
    const auto&   v2rhotau_ref = ref_vals.v2rhotau;
    const auto&   v2sigma2_ref = ref_vals.v2sigma2;
    const auto&   v2sigmalapl_ref = ref_vals.v2sigmalapl;
    const auto&   v2sigmatau_ref = ref_vals.v2sigmatau;
    const auto&   v2lapl2_ref = ref_vals.v2lapl2;
    const auto&   v2lapltau_ref = ref_vals.v2lapltau;
    const auto&   v2tau2_ref = ref_vals.v2tau2;

    // Allocate buffers
    std::vector<double> exc( func.exc_buffer_len( npts ),       fill_val_e );
    std::vector<double> vrho( func.vrho_buffer_len( npts ),     fill_val_vr );
    std::vector<double> vsigma( func.vsigma_buffer_len( npts ), fill_val_vs );
    std::vector<double> vlapl( func.vlapl_buffer_len( npts ),   fill_val_vl );
    std::vector<double> vtau( func.vtau_buffer_len( npts ),     fill_val_vt );

    std::vector<double> v2rho2( func.v2rho2_buffer_len( npts ), fill_val_v2rho2 );
    std::vector<double> v2rhosigma( func.v2rhosigma_buffer_len( npts ), fill_val_v2rhosigma );
    std::vector<double> v2rholapl( func.v2rholapl_buffer_len( npts ), fill_val_v2rholapl );
    std::vector<double> v2rhotau( func.v2rhotau_buffer_len( npts ), fill_val_v2rhotau );
    std::vector<double> v2sigma2( func.v2sigma2_buffer_len( npts ), fill_val_v2sigma2 );
    std::vector<double> v2sigmalapl( func.v2sigmalapl_buffer_len( npts ), fill_val_v2sigmalapl );
    std::vector<double> v2sigmatau( func.v2sigmatau_buffer_len( npts ), fill_val_v2sigmatau );
    std::vector<double> v2lapl2( func.v2lapl2_buffer_len( npts ), fill_val_v2lapl2 );
    std::vector<double> v2lapltau( func.v2lapltau_buffer_len( npts ), fill_val_v2lapltau );
    std::vector<double> v2tau2( func.v2tau2_buffer_len( npts ), fill_val_v2tau2 );

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

    if( interface == TestInterface::FXC_INC ) {

      // Evaluate XC kernel
      func.eval_fxc_inc( alpha, npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        v2rho2.data(), v2rhosigma.data(), v2rholapl.data(), v2rhotau.data(),
        v2sigma2.data(), v2sigmalapl.data(), v2sigmatau.data(), v2lapl2.data(),
        v2lapltau.data(), v2tau2.data() );

      // Check correctness
      for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
        CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );
      for( auto i = 0ul; i < func.v2rhosigma_buffer_len(npts); ++i )
        CHECK( v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) );
      for( auto i = 0ul; i < func.v2rholapl_buffer_len(npts); ++i )
        CHECK( v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]) );
      for( auto i = 0ul; i < func.v2rhotau_buffer_len(npts); ++i )
        CHECK( v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]) );
      for( auto i = 0ul; i < func.v2sigma2_buffer_len(npts); ++i )
        CHECK( v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]) );
      for( auto i = 0ul; i < func.v2sigmalapl_buffer_len(npts); ++i )
        CHECK( v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]) );
      for( auto i = 0ul; i < func.v2sigmatau_buffer_len(npts); ++i )
        CHECK( v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]) );
      for( auto i = 0ul; i < func.v2lapl2_buffer_len(npts); ++i )
        CHECK( v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]) );
      for( auto i = 0ul; i < func.v2lapltau_buffer_len(npts); ++i )
        CHECK( v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]) );
      for( auto i = 0ul; i < func.v2tau2_buffer_len(npts); ++i )
        CHECK( v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]) );

  }

  if( interface == TestInterface::VXC_FXC_INC ) {

    // Evaluate XC kernel
    func.eval_vxc_fxc_inc( alpha, npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
      vrho.data(), vsigma.data(), vlapl.data(), vtau.data(),
      v2rho2.data(), v2rhosigma.data(), v2rholapl.data(), v2rhotau.data(),
      v2sigma2.data(), v2sigmalapl.data(), v2sigmatau.data(), v2lapl2.data(),
      v2lapltau.data(), v2tau2.data() );
    // Check correctness
    for( auto i = 0ul; i < func.vrho_buffer_len(npts); ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < func.vsigma_buffer_len(npts); ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
    for( auto i = 0ul; i < func.vlapl_buffer_len(npts); ++i )
      CHECK( vlapl[i] == Approx(fill_val_vl + alpha * vlapl_ref[i]) );
    for( auto i = 0ul; i < func.vtau_buffer_len(npts); ++i )
        CHECK( vtau[i] == Approx(fill_val_vt + alpha * vtau_ref[i]) );

    for( auto i = 0ul; i < func.v2rho2_buffer_len(npts); ++i )
      CHECK( v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) );
    for( auto i = 0ul; i < func.v2rhosigma_buffer_len(npts); ++i )
      CHECK( v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) );
    for( auto i = 0ul; i < func.v2rholapl_buffer_len(npts); ++i )
      CHECK( v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]) );
    for( auto i = 0ul; i < func.v2rhotau_buffer_len(npts); ++i )
      CHECK( v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]) );
    for( auto i = 0ul; i < func.v2sigma2_buffer_len(npts); ++i )
      CHECK( v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]) );
    for( auto i = 0ul; i < func.v2sigmalapl_buffer_len(npts); ++i )
      CHECK( v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]) );
    for( auto i = 0ul; i < func.v2sigmatau_buffer_len(npts); ++i )
      CHECK( v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]) );
    for( auto i = 0ul; i < func.v2lapl2_buffer_len(npts); ++i )
      CHECK( v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]) );
    for( auto i = 0ul; i < func.v2lapltau_buffer_len(npts); ++i )
      CHECK( v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]) );
    for( auto i = 0ul; i < func.v2tau2_buffer_len(npts); ++i )
      CHECK( v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]) );

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

  const int len_v2rho2 = func_libxc.v2rho2_buffer_len( npts );
  const int len_v2rhosigma = func_libxc.v2rhosigma_buffer_len( npts );
  const int len_v2rholapl = func_libxc.v2rholapl_buffer_len( npts );
  const int len_v2rhotau = func_libxc.v2rhotau_buffer_len( npts );
  const int len_v2sigma2 = func_libxc.v2sigma2_buffer_len( npts );
  const int len_v2sigmalapl = func_libxc.v2sigmalapl_buffer_len( npts );
  const int len_v2sigmatau  = func_libxc.v2sigmatau_buffer_len( npts );
  const int len_v2lapl2     = func_libxc.v2lapl2_buffer_len( npts );
  const int len_v2lapltau = func_libxc.v2lapltau_buffer_len( npts );
  const int len_v2tau2 = func_libxc.v2tau2_buffer_len( npts );

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

  std::vector<double> v2rho2_libxc      ( len_v2rho2 );
  std::vector<double> v2rhosigma_libxc  ( len_v2rhosigma );
  std::vector<double> v2rholapl_libxc   ( len_v2rholapl );
  std::vector<double> v2rhotau_libxc    ( len_v2rhotau );
  std::vector<double> v2sigma2_libxc    ( len_v2sigma2 );
  std::vector<double> v2sigmalapl_libxc ( len_v2sigmalapl );
  std::vector<double> v2sigmatau_libxc  ( len_v2sigmatau );
  std::vector<double> v2lapl2_libxc     ( len_v2lapl2 );
  std::vector<double> v2lapltau_libxc   ( len_v2lapltau );
  std::vector<double> v2tau2_libxc      ( len_v2tau2 );

  std::vector<double> exc_builtin( func_builtin.exc_buffer_len(npts) );
  std::vector<double> vrho_builtin( func_builtin.vrho_buffer_len(npts) );
  std::vector<double> vsigma_builtin( func_builtin.vsigma_buffer_len(npts) );
  std::vector<double> vlapl_builtin( func_builtin.vlapl_buffer_len(npts) );
  std::vector<double> vtau_builtin( func_builtin.vtau_buffer_len(npts) );

  std::vector<double> v2rho2_builtin      ( func_builtin.v2rho2_buffer_len(npts) );
  std::vector<double> v2rhosigma_builtin  ( func_builtin.v2rhosigma_buffer_len(npts) );
  std::vector<double> v2rholapl_builtin   ( func_builtin.v2rholapl_buffer_len(npts) );
  std::vector<double> v2rhotau_builtin    ( func_builtin.v2rhotau_buffer_len(npts) );
  std::vector<double> v2sigma2_builtin    ( func_builtin.v2sigma2_buffer_len(npts) );
  std::vector<double> v2sigmalapl_builtin ( func_builtin.v2sigmalapl_buffer_len(npts) );
  std::vector<double> v2sigmatau_builtin  ( func_builtin.v2sigmatau_buffer_len(npts) );
  std::vector<double> v2lapl2_builtin     ( func_builtin.v2lapl2_buffer_len(npts) );
  std::vector<double> v2lapltau_builtin   ( func_builtin.v2lapltau_buffer_len(npts) );
  std::vector<double> v2tau2_builtin      ( func_builtin.v2tau2_buffer_len(npts) );

  if( func_libxc.is_lda() ) {

    if( interface == TestInterface::EXC ) {

      func_libxc.eval_exc( npts, rho_use.data(), exc_libxc.data() );
      func_builtin.eval_exc( npts, rho_use.data(), exc_builtin.data() );

    } else if( interface == TestInterface::EXC_VXC ) {

      func_libxc.eval_exc_vxc( npts, rho_use.data(), exc_libxc.data(),
        vrho_libxc.data() );
      func_builtin.eval_exc_vxc( npts, rho_use.data(), exc_builtin.data(),
        vrho_builtin.data() );

    } else if( interface == TestInterface::FXC ) {

      func_libxc.eval_fxc( npts, rho_use.data(), v2rho2_libxc.data() );
      func_builtin.eval_fxc( npts, rho_use.data(), v2rho2_builtin.data() );

    } else if( interface == TestInterface::VXC_FXC ) {

      func_libxc.eval_vxc_fxc( npts, rho_use.data(), vrho_libxc.data(), v2rho2_libxc.data() );
      func_builtin.eval_vxc_fxc( npts, rho_use.data(), vrho_builtin.data(), v2rho2_builtin.data() );

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

    } else if( interface == TestInterface::FXC ) {

      func_libxc.eval_fxc( npts, rho_use.data(), sigma_use.data(),
        v2rho2_libxc.data(), v2rhosigma_libxc.data(), v2sigma2_libxc.data() );
      func_builtin.eval_fxc( npts, rho_use.data(), sigma_use.data(),
        v2rho2_builtin.data(), v2rhosigma_builtin.data(), v2sigma2_builtin.data() );

    } else if( interface == TestInterface::VXC_FXC ) {

      func_libxc.eval_vxc_fxc( npts, rho_use.data(), sigma_use.data(), vrho_libxc.data(), vsigma_libxc.data(),
        v2rho2_libxc.data(), v2rhosigma_libxc.data(), v2sigma2_libxc.data() );
      func_builtin.eval_vxc_fxc( npts, rho_use.data(), sigma_use.data(), vrho_builtin.data(), vsigma_builtin.data(),
        v2rho2_builtin.data(), v2rhosigma_builtin.data(), v2sigma2_builtin.data() );

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

    } else if( interface == TestInterface::FXC ) {

      func_libxc.eval_fxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(),
        v2rho2_libxc.data(), v2rhosigma_libxc.data(), v2rholapl_libxc.data(),
        v2rhotau_libxc.data(), v2sigma2_libxc.data(), v2sigmalapl_libxc.data(),
        v2sigmatau_libxc.data(), v2lapl2_libxc.data(), v2lapltau_libxc.data(),
        v2tau2_libxc.data() );
      func_builtin.eval_fxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(),
        v2rho2_builtin.data(), v2rhosigma_builtin.data(), v2rholapl_builtin.data(),
        v2rhotau_builtin.data(), v2sigma2_builtin.data(), v2sigmalapl_builtin.data(),
        v2sigmatau_builtin.data(), v2lapl2_builtin.data(), v2lapltau_builtin.data(),
        v2tau2_builtin.data() );

    } else if( interface == TestInterface::VXC_FXC ) {

      func_libxc.eval_vxc_fxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), vrho_libxc.data(), vsigma_libxc.data(),
        vlapl_libxc.data(), vtau_libxc.data(),
        v2rho2_libxc.data(), v2rhosigma_libxc.data(), v2rholapl_libxc.data(),
        v2rhotau_libxc.data(), v2sigma2_libxc.data(), v2sigmalapl_libxc.data(),
        v2sigmatau_libxc.data(), v2lapl2_libxc.data(), v2lapltau_libxc.data(),
        v2tau2_libxc.data() );
      func_builtin.eval_vxc_fxc( npts, rho_use.data(), sigma_use.data(),
        lapl_use.data(), tau_use.data(), vrho_builtin.data(), vsigma_builtin.data(),
        vlapl_builtin.data(), vtau_builtin.data(),
        v2rho2_builtin.data(), v2rhosigma_builtin.data(), v2rholapl_builtin.data(),
        v2rhotau_builtin.data(), v2sigma2_builtin.data(), v2sigmalapl_builtin.data(),
        v2sigmatau_builtin.data(), v2lapl2_builtin.data(), v2lapltau_builtin.data(),
        v2tau2_builtin.data() );

    }

  }

  // Check correctness
  if ( interface == TestInterface::EXC || interface == TestInterface::EXC_VXC ) {
    for( auto i = 0ul; i < func_libxc.exc_buffer_len(npts); ++i ) {
      INFO( "EXC Fails: Kernel is " << kern );
      CHECK( exc_builtin[i] == Approx(exc_libxc[i]) );
    }
  }

  if( interface == TestInterface::EXC_VXC || interface == TestInterface::VXC_FXC ) {
    for( auto i = 0ul; i < func_libxc.vrho_buffer_len(npts); ++i ) {
      INFO( "VRHO Fails: Kernel is " << kern );
      CHECK( vrho_builtin[i] == Approx(vrho_libxc[i]) );
    }
    for( auto i = 0ul; i < func_libxc.vsigma_buffer_len(npts); ++i ) {
      INFO( "VSIGMA Fails: Kernel is " << kern );
      CHECK( vsigma_builtin[i] == Approx(vsigma_libxc[i]) );
    }
    for( auto i = 0ul; i < func_libxc.vlapl_buffer_len(npts); ++i ) {
      INFO( "VLAPL Fails: Kernel is " << kern );
      CHECK( vlapl_builtin[i] == Approx(vlapl_libxc[i]) );
    }
    for( auto i = 0ul; i < func_libxc.vtau_buffer_len(npts); ++i ) {
      INFO( "VTAU Fails: Kernel is " << kern );
      CHECK( vtau_builtin[i] == Approx(vtau_libxc[i]) );
    }
  }

  if( interface == TestInterface::FXC || interface == TestInterface::VXC_FXC ) {
    for( auto i = 0ul; i < len_v2rho2; ++i ) {
      INFO( "V2RHO2 Fails: Kernel is " << kern << ", builtin = " << v2rho2_builtin[i] << ", libxc = " << v2rho2_libxc[i] );
      bool is_close = (v2rho2_builtin[i] == Approx(v2rho2_libxc[i]) || v2rho2_builtin[i] == Approx(v2rho2_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhosigma; ++i ) {
      INFO( "V2RHOSIGMA Fails: Kernel is " << kern << ", builtin = " << v2rhosigma_builtin[i] << ", libxc = " << v2rhosigma_libxc[i] );
      bool is_close = (v2rhosigma_builtin[i] == Approx(v2rhosigma_libxc[i]) || v2rhosigma_builtin[i] == Approx(v2rhosigma_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rholapl; ++i ) {
      INFO( "V2RHOLAPL Fails: Kernel is " << kern << ", builtin = " << v2rholapl_builtin[i] << ", libxc = " << v2rholapl_libxc[i] );
      bool is_close = (v2rholapl_builtin[i] == Approx(v2rholapl_libxc[i]) || v2rholapl_builtin[i] == Approx(v2rholapl_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhotau; ++i ) {
      INFO( "V2RHOTAU Fails: Kernel is " << kern << ", builtin = " << v2rhotau_builtin[i] << ", libxc = " << v2rhotau_libxc[i] );
      bool is_close = (v2rhotau_builtin[i] == Approx(v2rhotau_libxc[i]) || v2rhotau_builtin[i] == Approx(v2rhotau_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigma2; ++i ) {
      INFO( "V2SIGMA2 Fails: Kernel is " << kern << ", builtin = " << v2sigma2_builtin[i] << ", libxc = " << v2sigma2_libxc[i] );
      bool is_close = (v2sigma2_builtin[i] == Approx(v2sigma2_libxc[i]) || v2sigma2_builtin[i] == Approx(v2sigma2_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmalapl; ++i ) {
      INFO( "V2SIGMALAPL Fails: Kernel is " << kern << ", builtin = " << v2sigmalapl_builtin[i] << ", libxc = " << v2sigmalapl_libxc[i] );
      bool is_close = (v2sigmalapl_builtin[i] == Approx(v2sigmalapl_libxc[i]) || v2sigmalapl_builtin[i] == Approx(v2sigmalapl_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmatau; ++i ) {
      INFO( "V2SIGMATAU Fails: Kernel is " << kern << ", builtin = " << v2sigmatau_builtin[i] << ", libxc = " << v2sigmatau_libxc[i] );
      bool is_close = (v2sigmatau_builtin[i] == Approx(v2sigmatau_libxc[i]) || v2sigmatau_builtin[i] == Approx(v2sigmatau_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapl2; ++i ) {
      INFO( "V2LAPL2 Fails: Kernel is " << kern << ", builtin = " << v2lapl2_builtin[i] << ", libxc = " << v2lapl2_libxc[i] );
      bool is_close = (v2lapl2_builtin[i] == Approx(v2lapl2_libxc[i]) || v2lapl2_builtin[i] == Approx(v2lapl2_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapltau; ++i ) {
      INFO( "V2LAPLTAU Fails: Kernel is " << kern << ", builtin = " << v2lapltau_builtin[i] << ", libxc = " << v2lapltau_libxc[i] );
      bool is_close = (v2lapltau_builtin[i] == Approx(v2lapltau_libxc[i]) || v2lapltau_builtin[i] == Approx(v2lapltau_libxc[i]).margin(1e-12));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2tau2; ++i ) {
      INFO( "V2TAU2 Fails: Kernel is " << kern << ", builtin = " << v2tau2_builtin[i] << ", libxc = " << v2tau2_libxc[i] );
      bool is_close = (v2tau2_builtin[i] == Approx(v2tau2_libxc[i]) || v2tau2_builtin[i] == Approx(v2tau2_libxc[i]).margin(1e-12));
      CHECK( is_close );
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

  SECTION( "Unpolarized Regular Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Regular,
        kern, Spin::Unpolarized );
    }
  }

  SECTION( "Unpolarized Regular Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Regular,
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

  SECTION( "Unpolarized Small Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Small,
        kern, Spin::Unpolarized );
    }
  }

  SECTION( "Unpolarized Small Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_unstable_small(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Small,
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

  SECTION( "Unpolarized Zero Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Zero,
        kern, Spin::Unpolarized );
    }
  }

  SECTION( "Unpolarized Zero Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Zero,
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

  SECTION( "Polarized Regular Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Regular,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Regular Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Regular,
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

  SECTION( "Polarized Small Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_unstable_small(kern)) continue;
      if(is_unstable_small_polarized_fxc_exchange_due_to_libxc_bug(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Small,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Small Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_unstable_small(kern)) continue;
      if(is_unstable_small_polarized_fxc_exchange_due_to_libxc_bug(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Small,
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

  SECTION( "Polarized Zero Eval : FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::FXC, EvalType::Zero,
        kern, Spin::Polarized );
    }
  }

  SECTION( "Polarized Zero Eval : VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      compare_libxc_builtin( TestInterface::VXC_FXC, EvalType::Zero,
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

  SECTION( "Builtin Unpolarized FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      kernel_test( TestInterface::FXC_INC, Backend::builtin, kern,
        Spin::Unpolarized );
    }
  }

  SECTION( "Builtin Unpolarized VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      if(is_epc(kern)) continue;
      kernel_test( TestInterface::VXC_FXC_INC, Backend::builtin, kern,
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

  SECTION( "Builtin Polarized FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      kernel_test( TestInterface::FXC_INC, Backend::builtin, kern,
        Spin::Polarized );
    }
  }

  SECTION( "Builtin Polarized VXC + FXC" ) {
    for( auto kern : builtin_supported_kernels ) {
      if(is_deorbitalized(kern)) continue;
      kernel_test( TestInterface::VXC_FXC_INC, Backend::builtin, kern,
        Spin::Polarized );
    }
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

#if defined(EXCHCXX_ENABLE_CUDA) || defined(EXCHCXX_ENABLE_HIP)

#ifdef EXCHCXX_ENABLE_HIP
#define cudaFree hipFree
#define cudaGetErrorString hipGetErrorString
#define cudaMalloc hipMalloc
#define cudaStream_t hipStream_t
#define cudaSuccess hipSuccess
#define cudaMemcpy hipMemcpy
#define cudaMemcpyDefault hipMemcpyDefault
#define cudaDeviceSynchronize hipDeviceSynchronize
#endif

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


void test_device_interface( TestInterface interface, EvalType evaltype,
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

  const int npts = 1;//npts_lda;

  if (polar == Spin::Unpolarized && !supports_unpolarized(kern)){
    CHECK_THROWS( XCKernel( backend, kern, polar ) );
    return;
  }
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

  size_t len_v2rho2      = func.v2rho2_buffer_len(npts);
  size_t len_v2rhosigma  = func.v2rhosigma_buffer_len(npts);
  size_t len_v2rholapl   = func.v2rholapl_buffer_len(npts);
  size_t len_v2rhotau   = func.v2rhotau_buffer_len(npts);
  size_t len_v2sigma2    = func.v2sigma2_buffer_len(npts);
  size_t len_v2sigmalapl = func.v2sigmalapl_buffer_len(npts);
  size_t len_v2sigmatau  = func.v2sigmatau_buffer_len(npts);
  size_t len_v2lapl2     = func.v2lapl2_buffer_len(npts);
  size_t len_v2lapltau   = func.v2lapltau_buffer_len(npts);
  size_t len_v2tau2      = func.v2tau2_buffer_len(npts);

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

  std::vector<double>
    v2rho2_ref( len_v2rho2 ),
    v2rhosigma_ref( len_v2rhosigma ),
    v2rholapl_ref( len_v2rholapl ),
    v2rhotau_ref( len_v2rhotau ),
    v2sigma2_ref( len_v2sigma2 ),
    v2sigmalapl_ref( len_v2sigmalapl ),
    v2sigmatau_ref( len_v2sigmatau ),
    v2lapl2_ref( len_v2lapl2 ),
    v2lapltau_ref( len_v2lapltau ),
    v2tau2_ref( len_v2tau2 );

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

  } else if( interface == TestInterface::FXC or interface == TestInterface::FXC_INC ) {

    if( func.is_lda() )
      func.eval_fxc( npts, rho.data(), v2rho2_ref.data() );
    else if( func.is_gga() )
      func.eval_fxc( npts, rho.data(), sigma.data(), v2rho2_ref.data(),
        v2rhosigma_ref.data(), v2sigma2_ref.data() );
    else if( func.is_mgga() )
      func.eval_fxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        v2rho2_ref.data(), v2rhosigma_ref.data(), v2rholapl_ref.data(),
        v2rhotau_ref.data(), v2sigma2_ref.data(), v2sigmalapl_ref.data(),
        v2sigmatau_ref.data(), v2lapl2_ref.data(), v2lapltau_ref.data(),
        v2tau2_ref.data() );
  } else if( interface == TestInterface::VXC_FXC or interface == TestInterface::VXC_FXC_INC ) {
    if( func.is_lda() )
      func.eval_vxc_fxc( npts, rho.data(), vrho_ref.data(), v2rho2_ref.data() );
    else if( func.is_gga() )
      func.eval_vxc_fxc( npts, rho.data(), sigma.data(), vrho_ref.data(),
        vsigma_ref.data(), v2rho2_ref.data(), v2rhosigma_ref.data(),
        v2sigma2_ref.data() );
    else if( func.is_mgga() )
      func.eval_vxc_fxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        vrho_ref.data(), vsigma_ref.data(), vlapl_ref.data(), vtau_ref.data(),
        v2rho2_ref.data(), v2rhosigma_ref.data(), v2rholapl_ref.data(),
        v2rhotau_ref.data(), v2sigma2_ref.data(), v2sigmalapl_ref.data(),
        v2sigmatau_ref.data(), v2lapl2_ref.data(), v2lapltau_ref.data(),
        v2tau2_ref.data() );
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

  double* v2rho2_device      = safe_cuda_malloc<double>( len_v2rho2      );
  double* v2rhosigma_device  = safe_cuda_malloc<double>( len_v2rhosigma  );
  double* v2rholapl_device   = safe_cuda_malloc<double>( len_v2rholapl   );
  double* v2rhotau_device    = safe_cuda_malloc<double>( len_v2rhotau    );
  double* v2sigma2_device    = safe_cuda_malloc<double>( len_v2sigma2    );
  double* v2sigmalapl_device = safe_cuda_malloc<double>( len_v2sigmalapl );
  double* v2sigmatau_device  = safe_cuda_malloc<double>( len_v2sigmatau  );
  double* v2lapl2_device     = safe_cuda_malloc<double>( len_v2lapl2     );
  double* v2lapltau_device   = safe_cuda_malloc<double>( len_v2lapltau   );
  double* v2tau2_device      = safe_cuda_malloc<double>( len_v2tau2      );

  // H2D Copy of rho / sigma
  safe_cuda_cpy( rho_device, rho.data(), len_rho_buffer );
  if( func.is_gga() or func.is_mgga() )
    safe_cuda_cpy( sigma_device, sigma.data(), len_sigma_buffer );
  if( func.is_mgga() )
    safe_cuda_cpy( tau_device, tau.data(), len_tau_buffer );
  if( func.needs_laplacian() )
    safe_cuda_cpy( lapl_device, lapl.data(), len_lapl_buffer );

  const double alpha = 3.14;
  const double fill_val_e = 0.1;
  const double fill_val_vr = 1.;
  const double fill_val_vs = 2.;
  const double fill_val_vl = 3.;
  const double fill_val_vt = 4.;
  const double fill_val_v2rho2 = 10.;
  const double fill_val_v2rhosigma = 11.;
  const double fill_val_v2rholapl = 12.;
  const double fill_val_v2rhotau = 13.;
  const double fill_val_v2sigma2 = 14.;
  const double fill_val_v2sigmalapl = 15.;
  const double fill_val_v2sigmatau = 16.;
  const double fill_val_v2lapl2 = 17.;
  const double fill_val_v2lapltau = 18.;
  const double fill_val_v2tau2 = 19.;

  std::vector<double>
    exc( len_exc_buffer, fill_val_e ), vrho( len_vrho_buffer, fill_val_vr ),
    vsigma( len_vsigma_buffer, fill_val_vs ), vlapl(len_vlapl_buffer, fill_val_vl),
    vtau(len_vtau_buffer, fill_val_vt);
  std::vector<double>
    v2rho2( len_v2rho2, fill_val_v2rho2 ),
    v2rhosigma( len_v2rhosigma, fill_val_v2rhosigma ),
    v2rholapl( len_v2rholapl, fill_val_v2rholapl ),
    v2rhotau( len_v2rhotau, fill_val_v2rhotau ),
    v2sigma2( len_v2sigma2, fill_val_v2sigma2 ),
    v2sigmalapl( len_v2sigmalapl, fill_val_v2sigmalapl ),
    v2sigmatau( len_v2sigmatau, fill_val_v2sigmatau ),
    v2lapl2( len_v2lapl2, fill_val_v2lapl2 ),
    v2lapltau( len_v2lapltau, fill_val_v2lapltau ),
    v2tau2( len_v2tau2, fill_val_v2tau2 );

  // H2D copy of initial values, tests clobber / increment
  safe_cuda_cpy( exc_device, exc.data(), len_exc_buffer );
  safe_cuda_cpy( vrho_device, vrho.data(), len_vrho_buffer );
  safe_cuda_cpy( v2rho2_device, v2rho2.data(), len_v2rho2 );
  if( func.is_gga() or func.is_mgga() ){
    safe_cuda_cpy( vsigma_device, vsigma.data(), len_vsigma_buffer );
    safe_cuda_cpy( v2rhosigma_device, v2rhosigma.data(), len_v2rhosigma );
    safe_cuda_cpy( v2sigma2_device, v2sigma2.data(), len_v2sigma2 );
  }
  if( func.is_mgga() ){
    safe_cuda_cpy( vtau_device, vtau.data(), len_vtau_buffer );
    safe_cuda_cpy( v2rhotau_device, v2rhotau.data(), len_v2rhotau );
    safe_cuda_cpy( v2sigmatau_device, v2sigmatau.data(), len_v2sigmatau );
    safe_cuda_cpy( v2tau2_device, v2tau2.data(), len_v2tau2 );
  }
  if( func.needs_laplacian() ){
    safe_cuda_cpy( vlapl_device, vlapl.data(), len_vlapl_buffer );
    safe_cuda_cpy( v2rholapl_device, v2rholapl.data(), len_v2rholapl );
    safe_cuda_cpy( v2sigmalapl_device, v2sigmalapl.data(), len_v2sigmalapl );
    safe_cuda_cpy( v2lapl2_device, v2lapl2.data(), len_v2lapl2 );
    safe_cuda_cpy( v2lapltau_device, v2lapltau.data(), len_v2lapltau );
  }


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

  } else if( interface == TestInterface::FXC ) {

    if( func.is_lda() )
      func.eval_fxc_device( npts, rho_device, v2rho2_device, stream );
    else if( func.is_gga() )
      func.eval_fxc_device( npts, rho_device, sigma_device, v2rho2_device,
        v2rhosigma_device, v2sigma2_device, stream );
    else if( func.is_mgga() )
      func.eval_fxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        v2rho2_device, v2rhosigma_device, v2rholapl_device, v2rhotau_device,
        v2sigma2_device, v2sigmalapl_device, v2sigmatau_device,
        v2lapl2_device, v2lapltau_device, v2tau2_device, stream );
  } else if( interface == TestInterface::FXC_INC ) {
    if( func.is_lda() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, v2rho2_device, stream );
    else if( func.is_gga() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        v2rho2_device, v2rhosigma_device, v2sigma2_device, stream );
    else if( func.is_mgga() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, v2rho2_device, v2rhosigma_device,
        v2rholapl_device, v2rhotau_device, v2sigma2_device,
        v2sigmalapl_device, v2sigmatau_device, v2lapl2_device,
        v2lapltau_device, v2tau2_device, stream );
  } else if( interface == TestInterface::VXC_FXC ) {

    if( func.is_lda() )
      func.eval_vxc_fxc_device( npts, rho_device, vrho_device, v2rho2_device, stream );
    else if( func.is_gga() )
      func.eval_vxc_fxc_device( npts, rho_device, sigma_device, vrho_device,
        vsigma_device, v2rho2_device, v2rhosigma_device, v2sigma2_device, stream );
    else if( func.is_mgga() )
      func.eval_vxc_fxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        vrho_device, vsigma_device, vlapl_device, vtau_device,
        v2rho2_device, v2rhosigma_device, v2rholapl_device,
        v2rhotau_device, v2sigma2_device, v2sigmalapl_device,
        v2sigmatau_device, v2lapl2_device, v2lapltau_device,
        v2tau2_device, stream );
  } else if( interface == TestInterface::VXC_FXC_INC ) {
    if( func.is_lda() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, vrho_device,
        v2rho2_device, stream );
    else if( func.is_gga() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        vrho_device, vsigma_device, v2rho2_device, v2rhosigma_device,
        v2sigma2_device, stream );
    else if( func.is_mgga() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, vrho_device, vsigma_device,
        vlapl_device, vtau_device, v2rho2_device, v2rhosigma_device,
        v2rholapl_device, v2rhotau_device, v2sigma2_device,
        v2sigmalapl_device, v2sigmatau_device, v2lapl2_device,
        v2lapltau_device, v2tau2_device, stream );
  }

  device_synchronize();

  // D2H of results
  safe_cuda_cpy( exc.data(), exc_device, len_exc_buffer );
  safe_cuda_cpy( vrho.data(), vrho_device, len_vrho_buffer );
  safe_cuda_cpy( v2rho2.data(), v2rho2_device, len_v2rho2 );
  if( func.is_gga() or func.is_mgga() ){
    safe_cuda_cpy( vsigma.data(), vsigma_device, len_vsigma_buffer );
    safe_cuda_cpy( v2rhosigma.data(), v2rhosigma_device, len_v2rhosigma );
    safe_cuda_cpy( v2sigma2.data(), v2sigma2_device, len_v2sigma2 );
  }
  if( func.is_mgga() ){
    safe_cuda_cpy( vtau.data(), vtau_device, len_vtau_buffer );
    safe_cuda_cpy( v2rhotau.data(), v2rhotau_device, len_v2rhotau );
    safe_cuda_cpy( v2sigmatau.data(), v2sigmatau_device, len_v2sigmatau );
    safe_cuda_cpy( v2tau2.data(), v2tau2_device, len_v2tau2 );
  }
  if( func.needs_laplacian() ){
    safe_cuda_cpy( vlapl.data(), vlapl_device, len_vlapl_buffer );
    safe_cuda_cpy( v2rholapl.data(), v2rholapl_device, len_v2rholapl );
    safe_cuda_cpy( v2sigmalapl.data(), v2sigmalapl_device, len_v2sigmalapl );
    safe_cuda_cpy( v2lapl2.data(), v2lapl2_device, len_v2lapl2 );
    safe_cuda_cpy( v2lapltau.data(), v2lapltau_device, len_v2lapltau );
  }

  // Check correctness
  if( interface == TestInterface::EXC_INC or interface == TestInterface::EXC_VXC_INC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
  } else if( interface == TestInterface::EXC or interface == TestInterface::EXC_VXC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(exc_ref[i]) );
  }

  if( interface == TestInterface::EXC_VXC_INC or interface == TestInterface::VXC_FXC_INC ) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
    for( auto i = 0ul; i < len_vlapl_buffer; ++i )
      CHECK( vlapl[i] == Approx(fill_val_vl + alpha * vlapl_ref[i]) );
    for( auto i = 0ul; i < len_vtau_buffer; ++i )
      CHECK( vtau[i] == Approx(fill_val_vt + alpha * vtau_ref[i]) );

  } else if(interface == TestInterface::EXC_VXC or interface == TestInterface::VXC_FXC) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i ){
      INFO( "Kernel is " << kern );
      CHECK( vrho[i] == Approx(vrho_ref[i]) );
    }
    for( auto i = 0ul; i < len_vsigma_buffer; ++i ) {
      INFO( "vsigma Fails: Kernel is " << kern << ", builtin device = " << vsigma[i] << ", builtin = " << vsigma_ref[i] );
      bool is_close = (vsigma[i] == Approx(vsigma_ref[i]) || vsigma[i] == Approx(vsigma_ref[i]).margin(1e-13));
      CHECK( is_close );
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

  if( interface == TestInterface::FXC or interface == TestInterface::VXC_FXC ) {
    for( auto i = 0ul; i < len_v2rho2; ++i ) {
      INFO( "V2RHO2 Fails: Kernel is " << kern << ", builtin device = " << v2rho2[i] << ", builtin = " << v2rho2_ref[i] );
      bool is_close = (v2rho2[i] == Approx(v2rho2_ref[i]) || v2rho2[i] == Approx(v2rho2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhosigma; ++i ) {
      INFO( "V2RHOSIGMA Fails: Kernel is " << kern << ", builtin device = " << v2rhosigma[i] << ", builtin = " << v2rhosigma_ref[i] );
      bool is_close = (v2rhosigma[i] == Approx(v2rhosigma_ref[i]) || v2rhosigma[i] == Approx(v2rhosigma_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rholapl; ++i ) {
      INFO( "V2RHOLAPL Fails: Kernel is " << kern << ", builtin device = " << v2rholapl[i] << ", builtin = " << v2rholapl_ref[i] );
      bool is_close = (v2rholapl[i] == Approx(v2rholapl_ref[i]) || v2rholapl[i] == Approx(v2rholapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhotau; ++i ) {
      INFO( "V2RHOTAU Fails: Kernel is " << kern << ", builtin device = " << v2rhotau[i] << ", builtin = " << v2rhotau_ref[i] );
      bool is_close = (v2rhotau[i] == Approx(v2rhotau_ref[i]) || v2rhotau[i] == Approx(v2rhotau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigma2; ++i ) {
      INFO( "V2SIGMA2 Fails: Kernel is " << kern << ", builtin device = " << v2sigma2[i] << ", builtin = " << v2sigma2_ref[i] );
      bool is_close = (v2sigma2[i] == Approx(v2sigma2_ref[i]) || v2sigma2[i] == Approx(v2sigma2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmalapl; ++i ) {
      INFO( "V2SIGMALAPL Fails: Kernel is " << kern << ", builtin device = " << v2sigmalapl[i] << ", builtin = " << v2sigmalapl_ref[i] );
      bool is_close = (v2sigmalapl[i] == Approx(v2sigmalapl_ref[i]) || v2sigmalapl[i] == Approx(v2sigmalapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmatau; ++i ) {
      INFO( "V2SIGMATAU Fails: Kernel is " << kern << ", builtin device = " << v2sigmatau[i] << ", builtin = " << v2sigmatau_ref[i] );
      bool is_close = (v2sigmatau[i] == Approx(v2sigmatau_ref[i]) || v2sigmatau[i] == Approx(v2sigmatau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapl2; ++i ) {
      INFO( "V2LAPL2 Fails: Kernel is " << kern << ", builtin device = " << v2lapl2[i] << ", builtin = " << v2lapl2_ref[i] );
      bool is_close = (v2lapl2[i] == Approx(v2lapl2_ref[i]) || v2lapl2[i] == Approx(v2lapl2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapltau; ++i ) {
      INFO( "V2LAPLTAU Fails: Kernel is " << kern << ", builtin device = " << v2lapltau[i] << ", builtin = " << v2lapltau_ref[i] );
      bool is_close = (v2lapltau[i] == Approx(v2lapltau_ref[i]) || v2lapltau[i] == Approx(v2lapltau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2tau2; ++i ) {
      INFO( "V2TAU2 Fails: Kernel is " << kern << ", builtin device = " << v2tau2[i] << ", builtin = " << v2tau2_ref[i] );
      bool is_close = (v2tau2[i] == Approx(v2tau2_ref[i]) || v2tau2[i] == Approx(v2tau2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
  } else if( interface == TestInterface::FXC_INC or interface == TestInterface::VXC_FXC_INC ) {
    for( auto i = 0ul; i < len_v2rho2; ++i ) {
      INFO( "V2RHO2 Fails: Kernel is " << kern << ", builtin device = " << v2rho2[i] << ", builtin = " << v2rho2_ref[i] );
      bool is_close = (v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) || v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhosigma; ++i ) {
      INFO( "V2RHOSIGMA Fails: Kernel is " << kern << ", builtin device = " << v2rhosigma[i] << ", builtin = " << v2rhosigma_ref[i] );
      bool is_close = (v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) || v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rholapl; ++i ) {
      INFO( "V2RHOLAPL Fails: Kernel is " << kern << ", builtin device = " << v2rholapl[i] << ", builtin = " << v2rholapl_ref[i] );
      bool is_close = (v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]) || v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhotau; ++i ) {
      INFO( "V2RHOTAU Fails: Kernel is " << kern << ", builtin device = " << v2rhotau[i] << ", builtin = " << v2rhotau_ref[i] );
      bool is_close = (v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]) || v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigma2; ++i ) {
      INFO( "V2SIGMA2 Fails: Kernel is " << kern << ", builtin device = " << v2sigma2[i] << ", builtin = " << v2sigma2_ref[i] );
      bool is_close = (v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]) || v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmalapl; ++i ) {
      INFO( "V2SIGMALAPL Fails: Kernel is " << kern << ", builtin device = " << v2sigmalapl[i] << ", builtin = " << v2sigmalapl_ref[i] );
      bool is_close = (v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]) || v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmatau; ++i ) {
      INFO( "V2SIGMATAU Fails: Kernel is " << kern << ", builtin device = " << v2sigmatau[i] << ", builtin = " << v2sigmatau_ref[i] );
      bool is_close = (v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]) || v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapl2; ++i ) {
      INFO( "V2LAPL2 Fails: Kernel is " << kern << ", builtin device = " << v2lapl2[i] << ", builtin = " << v2lapl2_ref[i] );
      bool is_close = (v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]) || v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapltau; ++i ) {
      INFO( "V2LAPLTAU Fails: Kernel is " << kern << ", builtin device = " << v2lapltau[i] << ", builtin = " << v2lapltau_ref[i] );
      bool is_close = (v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]) || v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2tau2; ++i ) {
      INFO( "V2TAU2 Fails: Kernel is " << kern << ", builtin device = " << v2tau2[i] << ", builtin = " << v2tau2_ref[i] );
      bool is_close = (v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]) || v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
  }
  // Free device memory
  cuda_free_all( rho_device, sigma_device, exc_device, vrho_device, vsigma_device, lapl_device, tau_device,
    vlapl_device, vtau_device,
    v2rho2_device, v2rhosigma_device, v2rholapl_device, v2rhotau_device,
    v2sigma2_device, v2sigmalapl_device, v2sigmatau_device,
    v2lapl2_device, v2lapltau_device, v2tau2_device );

}

#endif // EXCHCXX_ENABLE_CUDA/HIP



#ifdef EXCHCXX_ENABLE_SYCL

inline sycl::queue q{ sycl::default_selector_v,
    sycl::property_list{sycl::property::queue::in_order{}} };

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


void test_device_interface( TestInterface interface, EvalType evaltype,
                          Backend backend, Kernel kern, Spin polar) {

  size_t npts_lda, npts_gga, npts_mgga, npts_lapl;
  std::vector<double> ref_rho, ref_sigma, ref_lapl, ref_tau;
  std::tie(npts_lda, ref_rho  )  = load_reference_density( polar );
  std::tie(npts_gga, ref_sigma)  = load_reference_sigma  ( polar );
  std::tie(npts_lapl, ref_lapl)  = load_reference_lapl   ( polar );
  std::tie(npts_mgga, ref_tau)   = load_reference_tau    ( polar );

  REQUIRE( npts_lda == npts_gga );
  REQUIRE( npts_lda == npts_mgga );
  REQUIRE( npts_lda == npts_lapl );

  const int npts = 1;//npts_lda;

  if (polar == Spin::Unpolarized && !supports_unpolarized(kern)){
    CHECK_THROWS( XCKernel( backend, kern, polar ) );
    return;
  }
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

  size_t len_v2rho2      = func.v2rho2_buffer_len(npts);
  size_t len_v2rhosigma  = func.v2rhosigma_buffer_len(npts);
  size_t len_v2rholapl   = func.v2rholapl_buffer_len(npts);
  size_t len_v2rhotau   = func.v2rhotau_buffer_len(npts);
  size_t len_v2sigma2    = func.v2sigma2_buffer_len(npts);
  size_t len_v2sigmalapl = func.v2sigmalapl_buffer_len(npts);
  size_t len_v2sigmatau  = func.v2sigmatau_buffer_len(npts);
  size_t len_v2lapl2     = func.v2lapl2_buffer_len(npts);
  size_t len_v2lapltau   = func.v2lapltau_buffer_len(npts);
  size_t len_v2tau2      = func.v2tau2_buffer_len(npts);

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

  std::vector<double>
    v2rho2_ref( len_v2rho2 ),
    v2rhosigma_ref( len_v2rhosigma ),
    v2rholapl_ref( len_v2rholapl ),
    v2rhotau_ref( len_v2rhotau ),
    v2sigma2_ref( len_v2sigma2 ),
    v2sigmalapl_ref( len_v2sigmalapl ),
    v2sigmatau_ref( len_v2sigmatau ),
    v2lapl2_ref( len_v2lapl2 ),
    v2lapltau_ref( len_v2lapltau ),
    v2tau2_ref( len_v2tau2 );

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

  } else if( interface == TestInterface::FXC or interface == TestInterface::FXC_INC ) {

    if( func.is_lda() )
      func.eval_fxc( npts, rho.data(), v2rho2_ref.data() );
    else if( func.is_gga() )
      func.eval_fxc( npts, rho.data(), sigma.data(), v2rho2_ref.data(),
        v2rhosigma_ref.data(), v2sigma2_ref.data() );
    else if( func.is_mgga() )
      func.eval_fxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        v2rho2_ref.data(), v2rhosigma_ref.data(), v2rholapl_ref.data(),
        v2rhotau_ref.data(), v2sigma2_ref.data(), v2sigmalapl_ref.data(),
        v2sigmatau_ref.data(), v2lapl2_ref.data(), v2lapltau_ref.data(),
        v2tau2_ref.data() );
  } else if( interface == TestInterface::VXC_FXC or interface == TestInterface::VXC_FXC_INC ) {
    if( func.is_lda() )
      func.eval_vxc_fxc( npts, rho.data(), vrho_ref.data(), v2rho2_ref.data() );
    else if( func.is_gga() )
      func.eval_vxc_fxc( npts, rho.data(), sigma.data(), vrho_ref.data(),
        vsigma_ref.data(), v2rho2_ref.data(), v2rhosigma_ref.data(),
        v2sigma2_ref.data() );
    else if( func.is_mgga() )
      func.eval_vxc_fxc( npts, rho.data(), sigma.data(), lapl.data(), tau.data(),
        vrho_ref.data(), vsigma_ref.data(), vlapl_ref.data(), vtau_ref.data(),
        v2rho2_ref.data(), v2rhosigma_ref.data(), v2rholapl_ref.data(),
        v2rhotau_ref.data(), v2sigma2_ref.data(), v2sigmalapl_ref.data(),
        v2sigmatau_ref.data(), v2lapl2_ref.data(), v2lapltau_ref.data(),
        v2tau2_ref.data() );
  }







  // Allocate device memory
  double* rho_device    = safe_sycl_malloc<double>( len_rho_buffer    , q);
  double* sigma_device  = safe_sycl_malloc<double>( len_sigma_buffer  , q);
  double* lapl_device   = safe_sycl_malloc<double>( len_lapl_buffer  , q);
  double* tau_device    = safe_sycl_malloc<double>( len_tau_buffer  , q);
  double* exc_device    = safe_sycl_malloc<double>( len_exc_buffer    , q);
  double* vrho_device   = safe_sycl_malloc<double>( len_vrho_buffer   , q);
  double* vsigma_device = safe_sycl_malloc<double>( len_vsigma_buffer , q);
  double* vlapl_device  = safe_sycl_malloc<double>( len_vlapl_buffer , q);
  double* vtau_device   = safe_sycl_malloc<double>( len_vtau_buffer , q);

  double* v2rho2_device      = safe_sycl_malloc<double>( len_v2rho2      , q);
  double* v2rhosigma_device  = safe_sycl_malloc<double>( len_v2rhosigma  , q);
  double* v2rholapl_device   = safe_sycl_malloc<double>( len_v2rholapl   , q);
  double* v2rhotau_device    = safe_sycl_malloc<double>( len_v2rhotau    , q);
  double* v2sigma2_device    = safe_sycl_malloc<double>( len_v2sigma2    , q);
  double* v2sigmalapl_device = safe_sycl_malloc<double>( len_v2sigmalapl , q);
  double* v2sigmatau_device  = safe_sycl_malloc<double>( len_v2sigmatau  , q);
  double* v2lapl2_device     = safe_sycl_malloc<double>( len_v2lapl2     , q);
  double* v2lapltau_device   = safe_sycl_malloc<double>( len_v2lapltau   , q);
  double* v2tau2_device      = safe_sycl_malloc<double>( len_v2tau2      , q);

  // H2D Copy of rho / sigma
  safe_sycl_cpy( rho_device, rho.data(), len_rho_buffer, q);
  if( func.is_gga() or func.is_mgga() )
    safe_sycl_cpy( sigma_device, sigma.data(), len_sigma_buffer, q);
  if( func.is_mgga() )
    safe_sycl_cpy( tau_device, tau.data(), len_tau_buffer, q);
  if( func.needs_laplacian() )
    safe_sycl_cpy( lapl_device, lapl.data(), len_lapl_buffer, q);

  const double alpha = 3.14;
  const double fill_val_e = 0.1;
  const double fill_val_vr = 1.;
  const double fill_val_vs = 2.;
  const double fill_val_vl = 3.;
  const double fill_val_vt = 4.;
  const double fill_val_v2rho2 = 10.;
  const double fill_val_v2rhosigma = 11.;
  const double fill_val_v2rholapl = 12.;
  const double fill_val_v2rhotau = 13.;
  const double fill_val_v2sigma2 = 14.;
  const double fill_val_v2sigmalapl = 15.;
  const double fill_val_v2sigmatau = 16.;
  const double fill_val_v2lapl2 = 17.;
  const double fill_val_v2lapltau = 18.;
  const double fill_val_v2tau2 = 19.;

  std::vector<double>
    exc( len_exc_buffer, fill_val_e ), vrho( len_vrho_buffer, fill_val_vr ),
    vsigma( len_vsigma_buffer, fill_val_vs ), vlapl(len_vlapl_buffer, fill_val_vl),
    vtau(len_vtau_buffer, fill_val_vt);
  std::vector<double>
    v2rho2( len_v2rho2, fill_val_v2rho2 ),
    v2rhosigma( len_v2rhosigma, fill_val_v2rhosigma ),
    v2rholapl( len_v2rholapl, fill_val_v2rholapl ),
    v2rhotau( len_v2rhotau, fill_val_v2rhotau ),
    v2sigma2( len_v2sigma2, fill_val_v2sigma2 ),
    v2sigmalapl( len_v2sigmalapl, fill_val_v2sigmalapl ),
    v2sigmatau( len_v2sigmatau, fill_val_v2sigmatau ),
    v2lapl2( len_v2lapl2, fill_val_v2lapl2 ),
    v2lapltau( len_v2lapltau, fill_val_v2lapltau ),
    v2tau2( len_v2tau2, fill_val_v2tau2 );

  // H2D copy of initial values, tests clobber / increment
  safe_sycl_cpy( exc_device, exc.data(), len_exc_buffer, q);
  safe_sycl_cpy( vrho_device, vrho.data(), len_vrho_buffer, q);
  safe_sycl_cpy( v2rho2_device, v2rho2.data(), len_v2rho2, q);
  if( func.is_gga() or func.is_mgga() ){
    safe_sycl_cpy( vsigma_device, vsigma.data(), len_vsigma_buffer, q);
    safe_sycl_cpy( v2rhosigma_device, v2rhosigma.data(), len_v2rhosigma, q);
    safe_sycl_cpy( v2sigma2_device, v2sigma2.data(), len_v2sigma2, q);
  }
  if( func.is_mgga() ){
    safe_sycl_cpy( vtau_device, vtau.data(), len_vtau_buffer, q);
    safe_sycl_cpy( v2rhotau_device, v2rhotau.data(), len_v2rhotau, q);
    safe_sycl_cpy( v2sigmatau_device, v2sigmatau.data(), len_v2sigmatau, q);
    safe_sycl_cpy( v2tau2_device, v2tau2.data(), len_v2tau2, q);
  }
  if( func.needs_laplacian() ){
    safe_sycl_cpy( vlapl_device, vlapl.data(), len_vlapl_buffer, q);
    safe_sycl_cpy( v2rholapl_device, v2rholapl.data(), len_v2rholapl, q);
    safe_sycl_cpy( v2sigmalapl_device, v2sigmalapl.data(), len_v2sigmalapl, q);
    safe_sycl_cpy( v2lapl2_device, v2lapl2.data(), len_v2lapl2, q);
    safe_sycl_cpy( v2lapltau_device, v2lapltau.data(), len_v2lapltau, q);
  }


  // Evaluate functional on device
  if( interface == TestInterface::EXC ) {

    if( func.is_lda() )
      func.eval_exc_device( npts, rho_device, exc_device, &q );
    else if( func.is_gga() )
      func.eval_exc_device( npts, rho_device, sigma_device, exc_device,
        &q );
    else if( func.is_mgga() )
      func.eval_exc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        exc_device, &q );

  } else if( interface == TestInterface::EXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_inc_device( alpha, npts, rho_device, exc_device, &q );
    else if( func.is_gga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, exc_device,
        &q );
    else if( func.is_mgga() )
      func.eval_exc_inc_device( alpha, npts, rho_device, sigma_device, lapl_device,
        tau_device, exc_device, &q );

  } else if( interface == TestInterface::EXC_VXC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_device( npts, rho_device, exc_device, vrho_device, &q );
    else if( func.is_gga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, exc_device,
        vrho_device, vsigma_device, &q );
    else if( func.is_mgga() )
      func.eval_exc_vxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        exc_device, vrho_device, vsigma_device, vlapl_device, vtau_device, &q );

  } else if( interface == TestInterface::EXC_VXC_INC ) {

    if( func.is_lda() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, exc_device,
        vrho_device, &q );
    else if( func.is_gga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        exc_device, vrho_device, vsigma_device, &q );
    else if( func.is_mgga() )
      func.eval_exc_vxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, exc_device, vrho_device, vsigma_device,
        vlapl_device, vtau_device, &q );

  } else if( interface == TestInterface::FXC ) {

    if( func.is_lda() )
      func.eval_fxc_device( npts, rho_device, v2rho2_device, &q );
    else if( func.is_gga() )
      func.eval_fxc_device( npts, rho_device, sigma_device, v2rho2_device,
        v2rhosigma_device, v2sigma2_device, &q );
    else if( func.is_mgga() )
      func.eval_fxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        v2rho2_device, v2rhosigma_device, v2rholapl_device, v2rhotau_device,
        v2sigma2_device, v2sigmalapl_device, v2sigmatau_device,
        v2lapl2_device, v2lapltau_device, v2tau2_device, &q );
  } else if( interface == TestInterface::FXC_INC ) {
    if( func.is_lda() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, v2rho2_device, &q );
    else if( func.is_gga() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        v2rho2_device, v2rhosigma_device, v2sigma2_device, &q );
    else if( func.is_mgga() )
      func.eval_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, v2rho2_device, v2rhosigma_device,
        v2rholapl_device, v2rhotau_device, v2sigma2_device,
        v2sigmalapl_device, v2sigmatau_device, v2lapl2_device,
        v2lapltau_device, v2tau2_device, &q );
  } else if( interface == TestInterface::VXC_FXC ) {

    if( func.is_lda() )
      func.eval_vxc_fxc_device( npts, rho_device, vrho_device, v2rho2_device, &q );
    else if( func.is_gga() )
      func.eval_vxc_fxc_device( npts, rho_device, sigma_device, vrho_device,
        vsigma_device, v2rho2_device, v2rhosigma_device, v2sigma2_device, &q );
    else if( func.is_mgga() )
      func.eval_vxc_fxc_device( npts, rho_device, sigma_device, lapl_device, tau_device,
        vrho_device, vsigma_device, vlapl_device, vtau_device,
        v2rho2_device, v2rhosigma_device, v2rholapl_device,
        v2rhotau_device, v2sigma2_device, v2sigmalapl_device,
        v2sigmatau_device, v2lapl2_device, v2lapltau_device,
        v2tau2_device, &q );
  } else if( interface == TestInterface::VXC_FXC_INC ) {
    if( func.is_lda() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, vrho_device,
        v2rho2_device, &q );
    else if( func.is_gga() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        vrho_device, vsigma_device, v2rho2_device, v2rhosigma_device,
        v2sigma2_device, &q );
    else if( func.is_mgga() )
      func.eval_vxc_fxc_inc_device( alpha, npts, rho_device, sigma_device,
        lapl_device, tau_device, vrho_device, vsigma_device,
        vlapl_device, vtau_device, v2rho2_device, v2rhosigma_device,
        v2rholapl_device, v2rhotau_device, v2sigma2_device,
        v2sigmalapl_device, v2sigmatau_device, v2lapl2_device,
        v2lapltau_device, v2tau2_device, &q );
  }

  device_synchronize( q );

  // D2H of results
  safe_sycl_cpy( exc.data(), exc_device, len_exc_buffer, q);
  safe_sycl_cpy( vrho.data(), vrho_device, len_vrho_buffer, q);
  safe_sycl_cpy( v2rho2.data(), v2rho2_device, len_v2rho2, q);
  if( func.is_gga() or func.is_mgga() ){
    safe_sycl_cpy( vsigma.data(), vsigma_device, len_vsigma_buffer, q);
    safe_sycl_cpy( v2rhosigma.data(), v2rhosigma_device, len_v2rhosigma, q);
    safe_sycl_cpy( v2sigma2.data(), v2sigma2_device, len_v2sigma2, q);
  }
  if( func.is_mgga() ){
    safe_sycl_cpy( vtau.data(), vtau_device, len_vtau_buffer, q);
    safe_sycl_cpy( v2rhotau.data(), v2rhotau_device, len_v2rhotau, q);
    safe_sycl_cpy( v2sigmatau.data(), v2sigmatau_device, len_v2sigmatau, q);
    safe_sycl_cpy( v2tau2.data(), v2tau2_device, len_v2tau2, q);
  }
  if( func.needs_laplacian() ){
    safe_sycl_cpy( vlapl.data(), vlapl_device, len_vlapl_buffer, q);
    safe_sycl_cpy( v2rholapl.data(), v2rholapl_device, len_v2rholapl, q);
    safe_sycl_cpy( v2sigmalapl.data(), v2sigmalapl_device, len_v2sigmalapl, q);
    safe_sycl_cpy( v2lapl2.data(), v2lapl2_device, len_v2lapl2, q);
    safe_sycl_cpy( v2lapltau.data(), v2lapltau_device, len_v2lapltau, q);
  }

  // Check correctness
  if( interface == TestInterface::EXC_INC or interface == TestInterface::EXC_VXC_INC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(fill_val_e + alpha * exc_ref[i]) );
  } else if( interface == TestInterface::EXC or interface == TestInterface::EXC_VXC ) {
    for( auto i = 0ul; i < len_exc_buffer; ++i )
      CHECK( exc[i] == Approx(exc_ref[i]) );
  }

  if( interface == TestInterface::EXC_VXC_INC or interface == TestInterface::VXC_FXC_INC ) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i )
      CHECK( vrho[i] == Approx(fill_val_vr + alpha * vrho_ref[i]) );
    for( auto i = 0ul; i < len_vsigma_buffer; ++i )
      CHECK( vsigma[i] == Approx(fill_val_vs + alpha * vsigma_ref[i]) );
    for( auto i = 0ul; i < len_vlapl_buffer; ++i )
      CHECK( vlapl[i] == Approx(fill_val_vl + alpha * vlapl_ref[i]) );
    for( auto i = 0ul; i < len_vtau_buffer; ++i )
      CHECK( vtau[i] == Approx(fill_val_vt + alpha * vtau_ref[i]) );

  } else if(interface == TestInterface::EXC_VXC or interface == TestInterface::VXC_FXC) {

    for( auto i = 0ul; i < len_vrho_buffer; ++i ){
      INFO( "Kernel is " << kern );
      CHECK( vrho[i] == Approx(vrho_ref[i]) );
    }
    for( auto i = 0ul; i < len_vsigma_buffer; ++i ) {
      INFO( "vsigma Fails: Kernel is " << kern << ", builtin device = " << vsigma[i] << ", builtin = " << vsigma_ref[i] );
      bool is_close = (vsigma[i] == Approx(vsigma_ref[i]) || vsigma[i] == Approx(vsigma_ref[i]).margin(1e-13));
      CHECK( is_close );
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

  if( interface == TestInterface::FXC or interface == TestInterface::VXC_FXC ) {
    for( auto i = 0ul; i < len_v2rho2; ++i ) {
      INFO( "V2RHO2 Fails: Kernel is " << kern << ", builtin device = " << v2rho2[i] << ", builtin = " << v2rho2_ref[i] );
      bool is_close = (v2rho2[i] == Approx(v2rho2_ref[i]) || v2rho2[i] == Approx(v2rho2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhosigma; ++i ) {
      INFO( "V2RHOSIGMA Fails: Kernel is " << kern << ", builtin device = " << v2rhosigma[i] << ", builtin = " << v2rhosigma_ref[i] );
      bool is_close = (v2rhosigma[i] == Approx(v2rhosigma_ref[i]) || v2rhosigma[i] == Approx(v2rhosigma_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rholapl; ++i ) {
      INFO( "V2RHOLAPL Fails: Kernel is " << kern << ", builtin device = " << v2rholapl[i] << ", builtin = " << v2rholapl_ref[i] );
      bool is_close = (v2rholapl[i] == Approx(v2rholapl_ref[i]) || v2rholapl[i] == Approx(v2rholapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhotau; ++i ) {
      INFO( "V2RHOTAU Fails: Kernel is " << kern << ", builtin device = " << v2rhotau[i] << ", builtin = " << v2rhotau_ref[i] );
      bool is_close = (v2rhotau[i] == Approx(v2rhotau_ref[i]) || v2rhotau[i] == Approx(v2rhotau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigma2; ++i ) {
      INFO( "V2SIGMA2 Fails: Kernel is " << kern << ", builtin device = " << v2sigma2[i] << ", builtin = " << v2sigma2_ref[i] );
      bool is_close = (v2sigma2[i] == Approx(v2sigma2_ref[i]) || v2sigma2[i] == Approx(v2sigma2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmalapl; ++i ) {
      INFO( "V2SIGMALAPL Fails: Kernel is " << kern << ", builtin device = " << v2sigmalapl[i] << ", builtin = " << v2sigmalapl_ref[i] );
      bool is_close = (v2sigmalapl[i] == Approx(v2sigmalapl_ref[i]) || v2sigmalapl[i] == Approx(v2sigmalapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmatau; ++i ) {
      INFO( "V2SIGMATAU Fails: Kernel is " << kern << ", builtin device = " << v2sigmatau[i] << ", builtin = " << v2sigmatau_ref[i] );
      bool is_close = (v2sigmatau[i] == Approx(v2sigmatau_ref[i]) || v2sigmatau[i] == Approx(v2sigmatau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapl2; ++i ) {
      INFO( "V2LAPL2 Fails: Kernel is " << kern << ", builtin device = " << v2lapl2[i] << ", builtin = " << v2lapl2_ref[i] );
      bool is_close = (v2lapl2[i] == Approx(v2lapl2_ref[i]) || v2lapl2[i] == Approx(v2lapl2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapltau; ++i ) {
      INFO( "V2LAPLTAU Fails: Kernel is " << kern << ", builtin device = " << v2lapltau[i] << ", builtin = " << v2lapltau_ref[i] );
      bool is_close = (v2lapltau[i] == Approx(v2lapltau_ref[i]) || v2lapltau[i] == Approx(v2lapltau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2tau2; ++i ) {
      INFO( "V2TAU2 Fails: Kernel is " << kern << ", builtin device = " << v2tau2[i] << ", builtin = " << v2tau2_ref[i] );
      bool is_close = (v2tau2[i] == Approx(v2tau2_ref[i]) || v2tau2[i] == Approx(v2tau2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
  } else if( interface == TestInterface::FXC_INC or interface == TestInterface::VXC_FXC_INC ) {
    for( auto i = 0ul; i < len_v2rho2; ++i ) {
      INFO( "V2RHO2 Fails: Kernel is " << kern << ", builtin device = " << v2rho2[i] << ", builtin = " << v2rho2_ref[i] );
      bool is_close = (v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]) || v2rho2[i] == Approx(fill_val_v2rho2 + alpha * v2rho2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhosigma; ++i ) {
      INFO( "V2RHOSIGMA Fails: Kernel is " << kern << ", builtin device = " << v2rhosigma[i] << ", builtin = " << v2rhosigma_ref[i] );
      bool is_close = (v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]) || v2rhosigma[i] == Approx(fill_val_v2rhosigma + alpha * v2rhosigma_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rholapl; ++i ) {
      INFO( "V2RHOLAPL Fails: Kernel is " << kern << ", builtin device = " << v2rholapl[i] << ", builtin = " << v2rholapl_ref[i] );
      bool is_close = (v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]) || v2rholapl[i] == Approx(fill_val_v2rholapl + alpha * v2rholapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2rhotau; ++i ) {
      INFO( "V2RHOTAU Fails: Kernel is " << kern << ", builtin device = " << v2rhotau[i] << ", builtin = " << v2rhotau_ref[i] );
      bool is_close = (v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]) || v2rhotau[i] == Approx(fill_val_v2rhotau + alpha * v2rhotau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigma2; ++i ) {
      INFO( "V2SIGMA2 Fails: Kernel is " << kern << ", builtin device = " << v2sigma2[i] << ", builtin = " << v2sigma2_ref[i] );
      bool is_close = (v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]) || v2sigma2[i] == Approx(fill_val_v2sigma2 + alpha * v2sigma2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmalapl; ++i ) {
      INFO( "V2SIGMALAPL Fails: Kernel is " << kern << ", builtin device = " << v2sigmalapl[i] << ", builtin = " << v2sigmalapl_ref[i] );
      bool is_close = (v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]) || v2sigmalapl[i] == Approx(fill_val_v2sigmalapl + alpha * v2sigmalapl_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2sigmatau; ++i ) {
      INFO( "V2SIGMATAU Fails: Kernel is " << kern << ", builtin device = " << v2sigmatau[i] << ", builtin = " << v2sigmatau_ref[i] );
      bool is_close = (v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]) || v2sigmatau[i] == Approx(fill_val_v2sigmatau + alpha * v2sigmatau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapl2; ++i ) {
      INFO( "V2LAPL2 Fails: Kernel is " << kern << ", builtin device = " << v2lapl2[i] << ", builtin = " << v2lapl2_ref[i] );
      bool is_close = (v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]) || v2lapl2[i] == Approx(fill_val_v2lapl2 + alpha * v2lapl2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2lapltau; ++i ) {
      INFO( "V2LAPLTAU Fails: Kernel is " << kern << ", builtin device = " << v2lapltau[i] << ", builtin = " << v2lapltau_ref[i] );
      bool is_close = (v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]) || v2lapltau[i] == Approx(fill_val_v2lapltau + alpha * v2lapltau_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
    for( auto i = 0ul; i < len_v2tau2; ++i ) {
      INFO( "V2TAU2 Fails: Kernel is " << kern << ", builtin device = " << v2tau2[i] << ", builtin = " << v2tau2_ref[i] );
      bool is_close = (v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]) || v2tau2[i] == Approx(fill_val_v2tau2 + alpha * v2tau2_ref[i]).margin(1e-11));
      CHECK( is_close );
    }
  }
  // Free device memory
  sycl_free_all( q, rho_device, sigma_device, exc_device, vrho_device, vsigma_device, lapl_device, tau_device,
    vlapl_device, vtau_device,
    v2rho2_device, v2rhosigma_device, v2rholapl_device, v2rhotau_device,
    v2sigma2_device, v2sigmalapl_device, v2sigmatau_device,
    v2lapl2_device, v2lapltau_device, v2tau2_device );

}

#endif // EXCHCXX_ENABLE_SYCL

#ifdef EXCHCXX_ENABLE_DEVICE
TEST_CASE( "GPU Interfaces", "[xc-device]" ) {

  SECTION( "Libxc Functionals" ) {


    SECTION( "LDA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: FXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Regular Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: FXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Regular Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: FXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Regular Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }



    SECTION( "LDA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: FXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Small Eval Unpolarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: FXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Small Eval Unpolarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC + VXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: FXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;

        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Small Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "LDA Functionals: FXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Zero Eval Unpolarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "GGA Functionals: FXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Zero Eval Unpolarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
    }

    SECTION( "MGGA Functionals: FXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Zero Eval Unpolarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Unpolarized );
      }
    }






    SECTION( "LDA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: FXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Regular Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: FXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Regular Eval Polarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: FXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Regular Eval Polarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: FXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Small Eval Polarized" ) {
      for( auto kern : lda_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }


    SECTION( "GGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: FXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Small Eval Polarized" ) {
      for( auto kern : gga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC + VXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: FXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Small Eval Polarized" ) {
      for( auto kern : mgga_kernels ) {
        if(is_unstable_small(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }


    SECTION( "LDA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "LDA Functionals: FXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "LDA Functionals: VXC + FXC Zero Eval Polarized" ) {
      for( auto kern : lda_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "GGA Functionals: FXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "GGA Functionals: VXC + FXC Zero Eval Polarized" ) {
      for( auto kern : gga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: EXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: EXC + VXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
    }

    SECTION( "MGGA Functionals: FXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

    SECTION( "MGGA Functionals: VXC + FXC Zero Eval Polarized" ) {
      for( auto kern : mgga_kernels ){
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::libxc, kern, Spin::Polarized );
      }
    }

  }

  SECTION( "Builtin Functionals" ) {

    SECTION("EXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("FXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("FXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC + INC Regular: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC + VXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("FXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC + VXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("FXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC + INC Small: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("FXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("EXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("EXC + VXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
    }

    SECTION("FXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }

    SECTION("VXC + FXC + INC Zero: Unpolarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Unpolarized );
      }
    }



    SECTION("EXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("FXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("FXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC + INC Regular: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Regular,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + VXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("FXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + VXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small(kern)) continue;
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("FXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC + INC Small: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_unstable_small_2nd_deriv_device(kern)) continue;
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Small,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("FXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("EXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("EXC + VXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels )
        test_device_interface( TestInterface::EXC_VXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
    }

    SECTION("FXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::FXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

    SECTION("VXC + FXC + INC Zero: Polarized") {
      for( auto kern : builtin_supported_kernels ) {
        if(is_deorbitalized(kern)) continue;
        test_device_interface( TestInterface::VXC_FXC_INC, EvalType::Zero,
          Backend::builtin, kern, Spin::Polarized );
      }
    }

  }


}

#endif // EXCHCXX_ENABLE_DEVICE
