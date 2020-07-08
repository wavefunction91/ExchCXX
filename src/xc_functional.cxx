#include <exchcxx/xc_functional.hpp>

namespace ExchCXX {


std::vector< XCKernel > functional_factory( 
  const Backend        backend,
  const XCFunctional::Functional func,
  const Spin          polar
) {

  std::vector< XCKernel > kerns;

  using Kernel     = Kernel;
  using Functional = XCFunctional::Functional;

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
  else if( func == Functional::PBE0 )
    kerns = { 
      XCKernel( backend, Kernel::PBE0, polar )
    };
  else {
    assert( false && "FUNCTIONAL NYS" );
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
  const XCFunctional::Functional func,
  const Spin           polar
) :
  XCFunctional(functional_factory(backend,func,polar)) { }




XCFunctional& XCFunctional::operator=( const XCFunctional& ) = default;
XCFunctional& XCFunctional::operator=( XCFunctional&&      ) noexcept = default;

XCFunctional::XCFunctional( const XCFunctional& )       = default;
XCFunctional::XCFunctional( XCFunctional&& )  noexcept  = default;


LDA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;
    kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i )
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] *= kernels_[i].first;
  
  }

}


LDA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);
  const size_t len_vxc_buffer = vrho_buffer_len(N);

  std::vector<double> eps_scr, vxc_scr;
  if( kernels_.size() > 1 ) {
    eps_scr.resize( len_exc_buffer );
    vxc_scr.resize( len_vxc_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vxc, len_vxc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;
    double* vxc_eval = i ? vxc_scr.data() : vxc;
    kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vxc_eval);

    if( i ) {

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] += kernels_[i].first * eps_eval[k];
      for( auto k = 0ul; k < len_vxc_buffer; ++k ) 
        vxc[k] += kernels_[i].first * vxc_eval[k];

    } else {

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] *= kernels_[i].first;
      for( auto k = 0ul; k < len_vxc_buffer; ++k ) 
        vxc[k] *= kernels_[i].first;

    }

  }

}



// GGA Interfaces

GGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i )
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] *= kernels_[i].first;
  
  }

}


GGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr;
  if( kernels_.size() > 1 ) {
    eps_scr.resize( len_exc_buffer );
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr.data()    : eps;
    double* vrho_eval   = i ? vrho_scr.data()   : vrho;
    double* vsigma_eval = i ? vsigma_scr.data() : vsigma;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, 
        vsigma_eval );
    else
      kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

    if( i ) {

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] += kernels_[i].first * eps_eval[k];
      for( auto k = 0ul; k < len_vrho_buffer; ++k ) 
        vrho[k] += kernels_[i].first * vrho_eval[k];

      if( kernels_[i].second.is_gga() )
        for( auto k = 0ul; k < len_vsigma_buffer; ++k ) 
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

    } else {

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] *= kernels_[i].first;
      for( auto k = 0ul; k < len_vrho_buffer; ++k ) 
        vrho[k] *= kernels_[i].first;

      if( kernels_[i].second.is_gga() )
        for( auto k = 0ul; k < len_vsigma_buffer; ++k ) 
          vsigma[k] *= kernels_[i].first;

    }
  
  }

}




// mGGA Interfaces

MGGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_mgga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc(N, rho, sigma, lapl, tau, eps_eval);
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i )
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else
      for( auto k = 0ul; k < len_exc_buffer; ++k )
        eps[k] *= kernels_[i].first;
  
  }

}


MGGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_vlapl_buffer  = vlapl_buffer_len(N);
  const size_t len_vtau_buffer   = vtau_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr, vlapl_scr, vtau_scr;
  if( kernels_.size() > 1 ) {
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

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] += kernels_[i].first * eps_eval[k];
      for( auto k = 0ul; k < len_vrho_buffer; ++k ) 
        vrho[k] += kernels_[i].first * vrho_eval[k];

      if( kernels_[i].second.is_mgga() or kernels_[i].second.is_gga() )
        for( auto k = 0ul; k < len_vsigma_buffer; ++k ) 
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

      if( kernels_[i].second.is_mgga() ) {
        for( auto k = 0ul; k < len_vlapl_buffer; ++k ) 
          vlapl[k] += kernels_[i].first * vlapl_eval[k];
        for( auto k = 0ul; k < len_vtau_buffer; ++k ) 
          vtau[k] += kernels_[i].first * vtau_eval[k];
      }

    } else {

      for( auto k = 0ul; k < len_exc_buffer; ++k ) 
        eps[k] *= kernels_[i].first;
      for( auto k = 0ul; k < len_vrho_buffer; ++k ) 
        vrho[k] *= kernels_[i].first;

      if( kernels_[i].second.is_mgga() or kernels_[i].second.is_gga() )
        for( auto k = 0ul; k < len_vsigma_buffer; ++k ) 
          vsigma[k] *= kernels_[i].first;

      if( kernels_[i].second.is_mgga() ) {
        for( auto k = 0ul; k < len_vlapl_buffer; ++k ) 
          vlapl[k] *= kernels_[i].first;
        for( auto k = 0ul; k < len_vtau_buffer; ++k ) 
          vtau[k] *= kernels_[i].first;
      }

    }
  
  }

}




}; // ExchCXX


