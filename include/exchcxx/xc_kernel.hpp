#ifndef __INCLUDED_XC_KERNEL_HPP__
#define __INCLUDED_XC_KERNEL_HPP__

#include <xc.h> // Libxc

// Standard Libs
#include <vector>
#include <cassert>
#include <algorithm>

namespace ExchCXX {

class XCKernel {

protected:

  const int    polar_;  ///< Spin polarization
  xc_func_type kernel_; ///< Libxc kernel definition


  inline auto xc_info() const noexcept { return kernel_.info; };

public:

  XCKernel() = delete;
  
  XCKernel( const int kern, const int spin_polar );


  XCKernel( const XCKernel& );

  XCKernel( XCKernel&&      ) = default;

  // Destroy interal Libxc data
  ~XCKernel(){ xc_func_end( &kernel_ ); }



  inline bool is_lda() const noexcept {
    return kernel_.info->family == XC_FAMILY_LDA;
  }

  inline bool is_gga() const noexcept {
    return 
      (kernel_.info->family == XC_FAMILY_GGA    ) or
      (kernel_.info->family == XC_FAMILY_HYB_GGA);
  }

  inline bool is_mgga() const noexcept {
    return 
      (kernel_.info->family == XC_FAMILY_MGGA    ) or
      (kernel_.info->family == XC_FAMILY_HYB_MGGA);
  }

  inline bool is_hyb() const noexcept {
    return
      (kernel_.info->family == XC_FAMILY_HYB_GGA ) or
      (kernel_.info->family == XC_FAMILY_HYB_MGGA);
  }

  inline bool is_polarized() const noexcept {
    return polar_ == XC_POLARIZED;
  }


  inline double hyb_exx() const noexcept {
    return xc_hyb_exx_coef( &kernel_ );
  }



  // LDA interfaces
  void eval_exc( 
    const int     N, 
    const double* rho, 
    double*       eps 
  ) const; 
  


  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const; 

  // TODO: LDA kxc interfaces

  // GGA interface
  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const; 

  
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const;

  // TODO: GGA kxc interfaces  

};




class XCFunctional {

private:

  std::vector< std::pair<double, XCKernel> > kernels_;

public:

  XCFunctional() = delete;

/*
  template <typename... Args>
  XCFunctional( Args&&... args ) :
    kernels_( std::forward<Args>(args)... ){ }  
*/



  inline bool is_lda() const {
    return std::all_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_lda(); }
    );
  }

  inline bool is_gga() const {
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_gga(); }
    );
  }

  inline bool is_polarized() const {
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_polarized(); }
    );
  }


  // LDA Interfaces

  void eval_exc( 
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
        eps[k] += eps_eval[k];
    
    }

  }


  void eval_exc_vxc( 
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
          eps[k] += eps_eval[k];
          vxc[k] += vxc_eval[k];
        }
        if( is_polarized() )
        for( auto k = N; k < 2*N; ++k ) 
          vxc[k] += vxc_eval[k];
      }
    
    }

  }



  // GGA Interfaces

  void eval_exc( 
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
        eps[k] += eps_eval[k];
    
    }

  }


  void eval_exc_vxc( 
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

          eps[k] += eps_eval[k];
          vrho[k] += vrho_eval[k];

          if( kernels_[i].second.is_gga() )
            vsigma[k] += vsigma_eval[k];

        }
        if( is_polarized() ) {

          for( auto k = N; k < 2*N; ++k ) 
            vrho[k] += vrho_eval[k];

          if( kernels_[i].second.is_gga() )
          for( auto k = N; k < 3*N; ++k ) 
            vsigma[k] += vsigma_eval[k];

        }
      }
    
    }

  }
};





};

#endif
