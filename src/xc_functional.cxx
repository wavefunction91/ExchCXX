#include <exchcxx/xc_functional.hpp>

namespace ExchCXX {


std::vector< XCKernel > functional_factory( 
  const XCKernel::Backend        backend,
  const XCFunctional::Functional func,
  const XCKernel::Spin          polar
) {

  std::vector< XCKernel > kerns;

  using Kernel     = XCKernel::Kernel;
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



bool XCFunctional::sanity_check() const {

  // Must have one kernel
  if( not kernels_.size() ) return false;

  // Polarization is all or nothing
  int polar_one = kernels_[0].second.is_polarized();
  bool polar_all = std::all_of(
    kernels_.begin(), kernels_.end(),
    [&](const auto& a){ 
      return (int)a.second.is_polarized() == polar_one; 
    }
  ); 

  if( not polar_all ) return false;

  // If we made it, kernel is sane
  return true;
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
  const XCKernel::Backend        backend, 
  const XCFunctional::Functional func,
  const XCKernel::Spin           polar
) :
  XCFunctional(functional_factory(backend,func,polar)) { }




XCFunctional& XCFunctional::operator=( const XCFunctional& ) = default;
XCFunctional& XCFunctional::operator=( XCFunctional&&      ) noexcept = default;

XCFunctional::XCFunctional( const XCFunctional& )       = default;
XCFunctional::XCFunctional( XCFunctional&& )  noexcept  = default;


LDA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_lda() );

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;
    kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i ) 
      for( auto k = 0ul; k < N; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else
      for( auto k = 0ul; k < N; ++k )
        eps[k] = kernels_[i].first * eps[k];
  
  }

}


LDA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_lda() );

  std::vector<double> eps_scr, vxc_scr;
  if( kernels_.size() > 1 ) {
    eps_scr.resize( N );
    vxc_scr.resize( is_polarized() ? 2*N : N );
  }


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;
    double* vxc_eval = i ? vxc_scr.data() : vxc;
    kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vxc_eval);

    if( i ) {
      for( auto k = 0ul; k < N; ++k ) {
        eps[k] += kernels_[i].first * eps_eval[k];
        vxc[k] += kernels_[i].first * vxc_eval[k];
      }
      if( is_polarized() )
      for( auto k = N; k < 2*N; ++k ) 
        vxc[k] += kernels_[i].first * vxc_eval[k];
    } else {
      for( auto k = 0ul; k < N; ++k ) {
        eps[k] = kernels_[i].first * eps[k];
        vxc[k] = kernels_[i].first * vxc[k];
      }
      if( is_polarized() )
      for( auto k = N; k < 2*N; ++k ) 
        vxc[k] = kernels_[i].first * vxc[k];
    }
  
  }

}



// GGA Interfaces

GGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i ) 
      for( auto k = 0ul; k < N; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else 
      for( auto k = 0ul; k < N; ++k )
        eps[k] = kernels_[i].first * eps_eval[k];
  
  }

}


GGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  std::vector<double> eps_scr, vrho_scr, vsigma_scr;
  if( kernels_.size() > 1 ) {
    eps_scr.resize( N );
    vrho_scr.resize( is_polarized() ? 2*N : N );
    vsigma_scr.resize( is_polarized() ? 3*N : N );
  }


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr.data()    : eps;
    double* vrho_eval   = i ? vrho_scr.data()   : vrho;
    double* vsigma_eval = i ? vsigma_scr.data() : vsigma;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval );
    else
      kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

    if( i ) {
      for( auto k = 0ul; k < N; ++k ) {

        eps[k] += kernels_[i].first * eps_eval[k];
        vrho[k] += kernels_[i].first * vrho_eval[k];

        if( kernels_[i].second.is_gga() )
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

      }
      if( is_polarized() ) {

        for( auto k = N; k < 2*N; ++k ) 
          vrho[k] += kernels_[i].first * vrho_eval[k];

        if( kernels_[i].second.is_gga() )
        for( auto k = N; k < 3*N; ++k ) 
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

      }
    } else {
      for( auto k = 0ul; k < N; ++k ) {

        eps[k] = kernels_[i].first * eps[k];
        vrho[k] = kernels_[i].first * vrho[k];

        if( kernels_[i].second.is_gga() )
          vsigma[k] = kernels_[i].first * vsigma[k];

      }
      if( is_polarized() ) {

        for( auto k = N; k < 2*N; ++k ) 
          vrho[k] = kernels_[i].first * vrho[k];

        if( kernels_[i].second.is_gga() )
        for( auto k = N; k < 3*N; ++k ) 
          vsigma[k] = kernels_[i].first * vsigma[k];

      }
    }
  
  }

}




// mGGA Interfaces

MGGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  assert( is_mgga() );

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 ) 
    eps_scr.resize( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr.data() : eps;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc(N, rho, sigma, lapl, tau, eps_eval);
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc(N, rho, eps_eval);

    if( i ) 
      for( auto k = 0ul; k < N; ++k )
        eps[k] += kernels_[i].first * eps_eval[k];
    else 
      for( auto k = 0ul; k < N; ++k )
        eps[k] = kernels_[i].first * eps_eval[k];
  
  }

}


MGGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  assert( is_gga() );

  std::vector<double> eps_scr, vrho_scr, vsigma_scr, vlapl_scr, vtau_scr;
  if( kernels_.size() > 1 ) {
    eps_scr.resize( N );
    vrho_scr.resize( is_polarized() ? 2*N : N );
    vsigma_scr.resize( is_polarized() ? 3*N : N );
    vlapl_scr.resize( is_polarized() ? 2*N : N );
    vtau_scr.resize( is_polarized() ? 2*N : N );
  }


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr.data()    : eps;
    double* vrho_eval   = i ? vrho_scr.data()   : vrho;
    double* vsigma_eval = i ? vsigma_scr.data() : vsigma;
    double* vlapl_eval  = i ? vlapl_scr.data()  : vlapl;
    double* vtau_eval   = i ? vtau_scr.data()   : vtau;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_vxc(N, rho, sigma, lapl, tau, eps_eval, vrho_eval, vsigma_eval, vlapl_eval, vtau_eval );
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval );
    else
      kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

    if( i ) {
      for( auto k = 0ul; k < N; ++k ) {

        eps[k] += kernels_[i].first * eps_eval[k];
        vrho[k] += kernels_[i].first * vrho_eval[k];

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

        if( kernels_[i].second.is_mgga() ) {
          vlapl[k] += kernels_[i].first * vlapl_eval[k];
          vtau[k]  += kernels_[i].first * vtau_eval[k];
        }

      }
      if( is_polarized() ) {

        for( auto k = N; k < 2*N; ++k ) 
          vrho[k] += kernels_[i].first * vrho_eval[k];

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
        for( auto k = N; k < 3*N; ++k ) 
          vsigma[k] += kernels_[i].first * vsigma_eval[k];

        if( kernels_[i].second.is_mgga() )
        for( auto k = N; k < 2*N; ++k ) {
          vlapl[k] += kernels_[i].first * vlapl_eval[k];
          vtau[k]  += kernels_[i].first * vtau_eval[k];
        }

      }
    } else {
      for( auto k = 0ul; k < N; ++k ) {

        eps[k] = kernels_[i].first * eps[k];
        vrho[k] = kernels_[i].first * vrho[k];

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          vsigma[k] = kernels_[i].first * vsigma[k];

        if( kernels_[i].second.is_mgga() ) {
          vlapl[k] = kernels_[i].first * vlapl[k];
          vtau[k]  = kernels_[i].first * vtau[k];
        }

      }
      if( is_polarized() ) {

        for( auto k = N; k < 2*N; ++k ) 
          vrho[k] = kernels_[i].first * vrho[k];

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
        for( auto k = N; k < 3*N; ++k ) 
          vsigma[k] = kernels_[i].first * vsigma[k];

        if( kernels_[i].second.is_mgga() )
        for( auto k = N; k < 2*N; ++k ) { 
          vlapl[k] = kernels_[i].first * vlapl[k];
          vtau[k]  = kernels_[i].first * vtau[k];
        }

      }
    }
  
  }

}




}; // ExchCXX


