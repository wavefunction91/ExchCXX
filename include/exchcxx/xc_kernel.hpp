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

  int          polar_;  ///< Spin polarization
  xc_func_type kernel_; ///< Libxc kernel definition

  bool initialized_ = false;
  void throw_if_uninitialized() const { assert( initialized_ ); }


  inline auto xc_info() const { 
    throw_if_uninitialized();
    return kernel_.info; 
  };

  XCKernel( 
    xc_func_type  kern, 
    const int     spin_polar, 
    const bool    init
  ) : polar_(spin_polar), kernel_(kern), initialized_(init){ }

public:

  XCKernel() = delete;
  
  XCKernel( const int kern, const int spin_polar );
  XCKernel( const XCKernel& )                    ;
  XCKernel( XCKernel&&      )            noexcept;
  XCKernel& operator=( const XCKernel& )         ;
  XCKernel& operator=( XCKernel&&      ) noexcept;

  // Destroy interal Libxc data
  ~XCKernel();



  inline bool is_lda() const {
    throw_if_uninitialized();
    return kernel_.info->family == XC_FAMILY_LDA;
  }

  inline bool is_gga() const {
    throw_if_uninitialized();
    return 
      (kernel_.info->family == XC_FAMILY_GGA    ) or
      (kernel_.info->family == XC_FAMILY_HYB_GGA);
  }

  inline bool is_mgga() const {
    throw_if_uninitialized();
    return 
      (kernel_.info->family == XC_FAMILY_MGGA    ) or
      (kernel_.info->family == XC_FAMILY_HYB_MGGA);
  }

  inline bool is_hyb() const {
    throw_if_uninitialized();
    return
      (kernel_.info->family == XC_FAMILY_HYB_GGA ) or
      (kernel_.info->family == XC_FAMILY_HYB_MGGA);
  }

  inline bool is_polarized() const {
    throw_if_uninitialized();
    return polar_ == XC_POLARIZED;
  }


  inline double hyb_exx() const {
    throw_if_uninitialized();
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

  // TODO: LDA fxc/kxc interfaces  
    
    
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

  // TODO: GGA fxc/kxc interfaces  

  // mGGA interface
  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const; 

  
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       vlapl,
    double*       vtau
  ) const;

  // TODO: mGGA fxc/kxc interfaces  
    
};









};

#endif
