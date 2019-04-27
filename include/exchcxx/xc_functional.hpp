#ifndef __INCLUDED_XC_FUNCTIONAL_HPP__
#define __INCLUDED_XC_FUNCTIONAL_HPP__

#include <exchcxx/xc_kernel.hpp>

#include <numeric>

namespace ExchCXX {

class XCFunctional {

  using value_type = std::pair<double, XCKernel>;
private:

  std::vector< value_type > kernels_;

  bool sanity_check() const ; 
  void throw_if_not_sane() const { assert( sanity_check() ); }

public:

  XCFunctional();
  XCFunctional( const std::vector< XCKernel >& );
  XCFunctional( const std::initializer_list< value_type >& list );
  XCFunctional( const decltype(kernels_)& ks );
  XCFunctional( decltype(kernels_)&& ks );

  XCFunctional( const XCFunctional& )                    ;
  XCFunctional( XCFunctional&& )                 noexcept;
  XCFunctional& operator=( const XCFunctional& )         ;
  XCFunctional& operator=( XCFunctional&&      ) noexcept;



  inline bool is_lda() const {
    throw_if_not_sane();
    return std::all_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_lda(); }
    );
  }

  inline bool is_gga() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_gga(); }
    ) and not is_mgga();
  }

  inline bool is_mgga() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_mgga(); }
    );
  }

  inline bool is_polarized() const {
    throw_if_not_sane();
    // Polarization is all or nothing
    return kernels_[0].second.is_polarized(); 
  }

  inline bool is_hyb() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_hyb(); }
    );
  }

  inline double hyb_exx() const {
    throw_if_not_sane();
    return std::accumulate( 
      kernels_.begin(), kernels_.end(), 0.,
      [](const auto& x, const auto &y) { 
        return x + y.second.hyb_exx(); 
      }
    );

  }

  // LDA Interfaces

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



  // GGA Interfaces

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

};

}; // namespace ExchCXX

#endif
