#include "libxc_common.hpp"
namespace ExchCXX {


/*
template< XCKernel::Kernel Kern >
struct libxc_supports_kernel : public std::false_type { };

#define add_libxc_support(KERN) \
template<> \
struct libxc_supports_kernel< XCKernel::Kernel::KERN > : public std::true_type { };

add_libxc_support( SlaterExchange )
add_libxc_support( VWN3           )
add_libxc_support( VWN5           )
*/

std::unordered_map< XCKernel::Kernel, int > libxc_kernel_map {
  // LDA Functionals
  { XCKernel::Kernel::SlaterExchange, XC_LDA_X            },
  { XCKernel::Kernel::VWN3,           XC_LDA_C_VWN_3      },
  { XCKernel::Kernel::VWN5,           XC_LDA_C_VWN_RPA    },

  // GGA Functionals
  { XCKernel::Kernel::PBE_X,          XC_GGA_X_PBE        },
  { XCKernel::Kernel::PBE_C,          XC_GGA_C_PBE        },
  { XCKernel::Kernel::B88,            XC_GGA_X_B88        },
  { XCKernel::Kernel::LYP,            XC_GGA_C_LYP        },

  // Hybrid GGA Functionals
  { XCKernel::Kernel::B3LYP,          XC_HYB_GGA_XC_B3LYP },
  { XCKernel::Kernel::PBE0,           XC_HYB_GGA_XC_PBEH  },
};


inline bool libxc_supports_kernel( const XCKernel::Kernel kern ) {
  return libxc_kernel_map.find(kern) != libxc_kernel_map.end();
}





XCKernel libxc_kernel_factory(const XCKernel::Kernel kern, 
  const XCKernel::Spin spin_polar ) {

  assert( libxc_supports_kernel(kern) );

  return XCKernel( 
    std::make_unique< detail::LibxcKernelImpl >( kern, spin_polar ) );

}


namespace detail {

LibxcKernelImpl::LibxcKernelImpl(
  const XCKernel::Kernel    kern,
  const XCKernel::Spin  spin_polar
) : LibxcKernelImpl( libxc_kernel_map[kern], libxc_polar_map[spin_polar] ) { }  

LibxcKernelImpl::LibxcKernelImpl( 
  xc_func_type  kern, 
  const int     spin_polar, 
  const bool    init
) : polar_(spin_polar), kernel_(kern), initialized_(init){ }

LibxcKernelImpl::LibxcKernelImpl( 
  const int kern, 
  const int spin_polar
) : polar_(spin_polar) {

  // Initialize XC Kernel using Libxc
  int info = xc_func_init( &kernel_, kern, spin_polar );

  assert( info == 0 );

  initialized_ = true;
} 

LibxcKernelImpl::LibxcKernelImpl( const LibxcKernelImpl& other ) :
  LibxcKernelImpl( other.xc_info()->number, other.polar_ ){ };



LibxcKernelImpl::~LibxcKernelImpl() noexcept {
  if( initialized_ ) xc_func_end( &kernel_ );
}



std::unique_ptr< XCKernelImpl > LibxcKernelImpl::clone_() const {

  return std::make_unique< LibxcKernelImpl >( *this );

}



bool LibxcKernelImpl::is_lda_() const noexcept {
  return kernel_.info->family == XC_FAMILY_LDA;
}

bool LibxcKernelImpl::is_gga_() const noexcept {
  return 
    (kernel_.info->family == XC_FAMILY_GGA    ) or
    (kernel_.info->family == XC_FAMILY_HYB_GGA);
}

bool LibxcKernelImpl::is_mgga_() const noexcept {
  return 
    (kernel_.info->family == XC_FAMILY_MGGA    ) or
    (kernel_.info->family == XC_FAMILY_HYB_MGGA);
}

bool LibxcKernelImpl::is_hyb_() const noexcept {
  return
    (kernel_.info->family == XC_FAMILY_HYB_GGA ) or
    (kernel_.info->family == XC_FAMILY_HYB_MGGA);
}

bool LibxcKernelImpl::is_polarized_() const noexcept {
  return polar_ == XC_POLARIZED;
}


double LibxcKernelImpl::hyb_exx_() const noexcept {
  return xc_hyb_exx_coef( &kernel_ );
}


// LDA interfaces
LDA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc( &kernel_, N, rho, eps );

}


LDA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc_vxc( &kernel_, N, rho, eps, vxc );

}


// GGA interface
GGA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc( &kernel_, N, rho, sigma, eps );

}


GGA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc_vxc( &kernel_, N, rho, sigma, eps, vrho, vsigma );

}

  
// mGGA interface
MGGA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  assert( is_mgga() );
  xc_mgga_exc( &kernel_, N, rho, sigma, lapl, tau, eps );

}


MGGA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  assert( is_mgga() );
  xc_mgga_exc_vxc( &kernel_, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau );

}


};
};
