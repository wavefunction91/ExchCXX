#include <exchcxx/xc_functional.hpp>

namespace ExchCXX {

bool XCFunctional::sanity_check() {

  // Must have one kernel
  if( not kernels_.size() ) return false;

  // Polarization is all or nothing
  bool polar_one = kernels_[0].second.is_polarized();
  bool polar_all = std::any_of(
    kernels_.begin(), kernels_.end(),
    [&](const auto& a){ 
      return a.second.is_polarized() != polar_one; 
    }
  ); 

  if( not polar_all ) return false;

  // If we made it, kernel is sane
  return true;
}

void XCFunctional::eval_exc( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

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


void XCFunctional::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc
) const {

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

void XCFunctional::eval_exc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps 
) const {

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


void XCFunctional::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps, 
  double*       vrho,
  double*       vsigma
) const {

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


}; // ExchCXX


