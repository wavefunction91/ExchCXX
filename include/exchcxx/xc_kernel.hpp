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









};

#endif
