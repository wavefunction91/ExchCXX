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

#include <exchcxx/xc_functional.hpp>

namespace ExchCXX {


BidirectionalMap<std::string, Functional> functional_map{
    {{"SVWN3", Functional::SVWN3},
    {"SVWN5", Functional::SVWN5},
    {"BLYP", Functional::BLYP},
    {"B3LYP", Functional::B3LYP},
    {"PBE", Functional::PBE},
    {"SCAN", Functional::SCAN},
    {"R2SCAN", Functional::R2SCAN},
    {"R2SCANL", Functional::R2SCANL},
    {"M062X", Functional::M062X},
    {"PKZB", Functional::PKZB},
    {"revPBE", Functional::revPBE},
    {"PBE0", Functional::PBE0},
    {"EPC17_1", Functional::EPC17_1},
    {"EPC17_2", Functional::EPC17_2},
    {"EPC18_1", Functional::EPC18_1},
    {"EPC18_2", Functional::EPC18_2}}};

std::ostream &operator<<(std::ostream &out, Functional functional) {
  out << functional_map.key(functional);
  return out;
}

std::vector< XCKernel > functional_factory( 
  const Backend        backend,
  const Functional func,
  const Spin          polar
) {

  std::vector< XCKernel > kerns;

  if( func == Functional::SVWN3 )
    kerns = { 
      XCKernel( backend, Kernel::SlaterExchange, polar ),
      XCKernel( backend, Kernel::VWN3,           polar )
    };
  else if( func == Functional::SVWN5 )
    kerns = { 
      XCKernel( backend, Kernel::SlaterExchange, polar ),
      XCKernel( backend, Kernel::VWN5,           polar )
    };
  else if( func == Functional::BLYP )
    kerns = { 
      XCKernel( backend, Kernel::B88, polar ),
      XCKernel( backend, Kernel::LYP, polar )
    };
  else if( func == Functional::B3LYP )
    kerns = { 
      XCKernel( backend, Kernel::B3LYP, polar )
    };
  else if( func == Functional::PBE )
    kerns = {
        XCKernel( backend, Kernel::PBE_X, polar ),
        XCKernel( backend, Kernel::PBE_C, polar )
    };
  else if( func == Functional::SCAN )
    kerns = {
        XCKernel( backend, Kernel::SCAN_X, polar ),
        XCKernel( backend, Kernel::SCAN_C, polar )
    };
  else if( func == Functional::R2SCAN )
    kerns = {
        XCKernel( backend, Kernel::R2SCAN_X, polar ),
        XCKernel( backend, Kernel::R2SCAN_C, polar )
    };
  else if( func == Functional::M062X )
    kerns = {
        XCKernel( backend, Kernel::M062X_X, polar ),
        XCKernel( backend, Kernel::M062X_C, polar )
    };
  else if( func == Functional::R2SCANL )
    kerns = {
        XCKernel( backend, Kernel::R2SCANL_X, polar ),
        XCKernel( backend, Kernel::R2SCANL_C, polar )
    };
  else if( func == Functional::revPBE )
    kerns = {
        XCKernel( backend, Kernel::revPBE_X, polar ),
        XCKernel( backend, Kernel::PBE_C, polar )
    };
  else if( func == Functional::PBE0 )
    kerns = { 
      XCKernel( backend, Kernel::PBE0, polar )
    };
  else if( func == Functional::EPC17_1 )
    kerns = { 
      XCKernel( backend, Kernel::EPC17_1, polar )
    };
  else if( func == Functional::EPC17_2 )
    kerns = { 
      XCKernel( backend, Kernel::EPC17_2, polar )
    };
  else if( func == Functional::EPC18_1 )
    kerns = { 
      XCKernel( backend, Kernel::EPC18_1, polar )
    };
  else if( func == Functional::EPC18_2 )
    kerns = { 
      XCKernel( backend, Kernel::EPC18_2, polar )
    };
  else {
    EXCHCXX_BOOL_CHECK( "Functional NYI Through Builtin Backend", false );
  }


  return kerns;
}

XCFunctional::XCFunctional() = default;

XCFunctional::XCFunctional( const std::vector< XCKernel > &ks ) {

  for(const auto& k : ks )
    kernels_.push_back( { 1., k } );

}

XCFunctional::XCFunctional( const std::initializer_list< value_type >& list ) : kernels_{ list } { }

XCFunctional::XCFunctional( const std::vector<value_type>& ks ) :
  kernels_(ks) { }
XCFunctional::XCFunctional( std::vector<value_type>&& ks ) :
  kernels_(std::move(ks)) { }


XCFunctional::XCFunctional( 
  const Backend        backend, 
  const Functional func,
  const Spin           polar
) :
  XCFunctional(functional_factory(backend,func,polar)) { }



XCFunctional& XCFunctional::operator=( const XCFunctional& ) = default;
XCFunctional& XCFunctional::operator=( XCFunctional&&      ) noexcept = default;

XCFunctional::XCFunctional( const XCFunctional& )       = default;
XCFunctional::XCFunctional( XCFunctional&& )  noexcept  = default;





void _scal( size_t len, double alpha, double* v ) {
  for( auto k = 0ul; k < len; ++k ) v[k] *= alpha;
}
void _addscal( size_t len, double alpha, double* v, const double* w ) {
  for( auto k = 0ul; k < len; ++k ) v[k] += alpha * w[k];
}



LDA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_inc( 
        kernels_[i].first, N, rho, eps
      );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;
      kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );
  
    }
  }

}


LDA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);
  const size_t len_vxc_buffer = vrho_buffer_len(N);

  std::vector<double> eps_scr, vxc_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vxc_scr.resize( len_vxc_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vxc, len_vxc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_vxc_inc(
        kernels_[i].first, N, rho, eps, vxc
      );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;
      double* vxc_eval = i ? vxc_scr.data() : vxc;
      kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vxc_eval);

      if( i ) {

        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
        _addscal( len_vxc_buffer, kernels_[i].first, vxc, vxc_eval );

      } else {

        _scal( len_exc_buffer, kernels_[i].first, eps );
        _scal( len_vxc_buffer, kernels_[i].first, vxc );

      }

    }

  }

}



// GGA Interfaces

GGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, eps
        );
      else
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, eps
        );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
      else
        kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );

    }
  
  }

}


GGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, eps, vrho
        );

    } else {

      double* eps_eval    = i ? eps_scr.data()    : eps;
      double* vrho_eval   = i ? vrho_scr.data()   : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, 
          vsigma_eval );
      else
        kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

      if( i ) {

        _addscal( len_exc_buffer,    kernels_[i].first, eps,    eps_eval  );
        _addscal( len_vrho_buffer,   kernels_[i].first, vrho,   vrho_eval );
   
        if( kernels_[i].second.is_gga() )
          _addscal( len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval );

      } else {

        _scal( len_exc_buffer,    kernels_[i].first, eps  );
        _scal( len_vrho_buffer,   kernels_[i].first, vrho );

        if( kernels_[i].second.is_gga() )
          _scal( len_vsigma_buffer, kernels_[i].first, vsigma );

      }

    }
  
  }

}




// mGGA Interfaces

MGGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, eps
        );
      else
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, eps
        );

    } else { 

      double* eps_eval = i ? eps_scr.data() : eps;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc(N, rho, sigma, lapl, tau, eps_eval);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
      else
        kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );

    }
  
  }

}


MGGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_vlapl_buffer  = vlapl_buffer_len(N);
  const size_t len_vtau_buffer   = vtau_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr, vlapl_scr, vtau_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
    vlapl_scr.resize( len_vlapl_buffer );
    vtau_scr.resize( len_vtau_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );
  std::fill_n( vlapl, len_vlapl_buffer, 0. );
  std::fill_n( vtau, len_vtau_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps, 
          vrho, vsigma, vlapl, vtau 
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, eps, vrho
        );
    
    } else {

      double* eps_eval    = i ? eps_scr.data()    : eps;
      double* vrho_eval   = i ? vrho_scr.data()   : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;
      double* vlapl_eval  = i ? vlapl_scr.data()  : vlapl;
      double* vtau_eval   = i ? vtau_scr.data()   : vtau;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, lapl, tau, eps_eval, 
          vrho_eval, vsigma_eval, vlapl_eval, vtau_eval );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, 
          vsigma_eval );
      else
        kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

      if( i ) {

        _addscal( len_exc_buffer,    kernels_[i].first, eps,    eps_eval  );
        _addscal( len_vrho_buffer,   kernels_[i].first, vrho,   vrho_eval );
   
        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          _addscal( len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval );

        if( kernels_[i].second.needs_laplacian() ) 
          _addscal( len_vlapl_buffer, kernels_[i].first, vlapl, vlapl_eval );

        if( kernels_[i].second.is_mgga() ) 
          _addscal( len_vtau_buffer,  kernels_[i].first, vtau,  vtau_eval  );

      } else {

        _scal( len_exc_buffer,    kernels_[i].first, eps  );
        _scal( len_vrho_buffer,   kernels_[i].first, vrho );

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          _scal( len_vsigma_buffer, kernels_[i].first, vsigma );

        if( kernels_[i].second.needs_laplacian() ) 
          _scal( len_vlapl_buffer, kernels_[i].first, vlapl );

        if( kernels_[i].second.is_mgga() ) 
          _scal( len_vtau_buffer,  kernels_[i].first, vtau  );

      }

    }
  
  }

}




} // ExchCXX


