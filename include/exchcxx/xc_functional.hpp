#ifndef __INCLUDED_XC_FUNCTIONAL_HPP__
#define __INCLUDED_XC_FUNCTIONAL_HPP__

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/util/exchcxx_macros.hpp>

#include <numeric>

namespace ExchCXX {

class XCFunctional {

  using value_type = std::pair<double, XCKernel>;
private:

  std::vector< value_type > kernels_;

  inline bool sanity_check() const {

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


  void throw_if_not_sane() const { assert( sanity_check() ); }

public:

  enum Functional {
    SVWN3,
    SVWN5,
    BLYP,
    B3LYP,
    PBE0,
  };

  XCFunctional();
  XCFunctional( const std::vector< XCKernel >& );
  XCFunctional( const std::initializer_list< value_type >& list );
  XCFunctional( const decltype(kernels_)& ks );
  XCFunctional( decltype(kernels_)&& ks );

  XCFunctional( const XCKernel::Backend, const Functional, const XCKernel::Spin );
  XCFunctional( const Functional func, const XCKernel::Spin polar) :
    XCFunctional( XCKernel::Backend::libxc, func, polar) { };

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
  LDA_EXC_GENERATOR(     eval_exc     ) const;
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // GGA Interfaces
  GGA_EXC_GENERATOR(     eval_exc     ) const;
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // mGGA interface
  MGGA_EXC_GENERATOR(     eval_exc     ) const;
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;


  // Device code
#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA Interfaces
  LDA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

  // GGA Interfaces
  GGA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

  // mGGA interface
  MGGA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

#endif

};

} // namespace ExchCXX

#endif
