#include "libxc_common.hpp"

#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/util/unused.hpp>
namespace ExchCXX {


/*
template< Kernel Kern >
struct libxc_supports_kernel : public std::false_type { };

#define add_libxc_support(KERN) \
template<> \
struct libxc_supports_kernel< Kernel::KERN > : public std::true_type { };

add_libxc_support( SlaterExchange )
add_libxc_support( VWN3           )
add_libxc_support( VWN5           )
*/

std::unordered_map< Kernel, int > libxc_kernel_map {
  // LDA Functionals
  { Kernel::SlaterExchange, XC_LDA_X            },
  { Kernel::VWN3,           XC_LDA_C_VWN_3      },
  { Kernel::VWN5,           XC_LDA_C_VWN_RPA    },
  { Kernel::PZ81,           XC_LDA_C_PZ         },
  { Kernel::PZ81_MOD,       XC_LDA_C_PZ_MOD     },
  { Kernel::PW91_LDA,       XC_LDA_C_PW         },
  { Kernel::PW91_LDA_MOD,   XC_LDA_C_PW_MOD     },
  { Kernel::PW91_LDA_RPA,   XC_LDA_C_PW_RPA     },

  // GGA Functionals
  { Kernel::PBE_X,          XC_GGA_X_PBE        },
  { Kernel::PBE_C,          XC_GGA_C_PBE        },
  { Kernel::revPBE_X,       XC_GGA_X_PBE_R      },
  { Kernel::B88,            XC_GGA_X_B88        },
  { Kernel::LYP,            XC_GGA_C_LYP        },

  // Hybrid GGA Functionals
  { Kernel::B3LYP,          XC_HYB_GGA_XC_B3LYP },
  { Kernel::PBE0,           XC_HYB_GGA_XC_PBEH  },
};


inline bool libxc_supports_kernel( const Kernel kern ) {
  return libxc_kernel_map.find(kern) != libxc_kernel_map.end();
}





XCKernel libxc_kernel_factory(const Kernel kern, 
  const Spin spin_polar ) {

  assert( libxc_supports_kernel(kern) );

  return XCKernel( 
    std::make_unique< detail::LibxcKernelImpl >( kern, spin_polar ) );

}


namespace detail {

LibxcKernelImpl::LibxcKernelImpl(
  const Kernel    kern,
  const Spin  spin_polar
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
  if( info )
    throw std::runtime_error("Libxc Kernel Init Failed");

  initialized_ = true;
} 

LibxcKernelImpl::LibxcKernelImpl( const LibxcKernelImpl& other ) :
  LibxcKernelImpl( other.xc_info()->number, other.polar_ ){ }



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

bool LibxcKernelImpl::supports_inc_interface_() const noexcept {
  return false;
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


UNUSED_INC_INTERFACE_GENERATOR( LDA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( LDA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )


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

UNUSED_INC_INTERFACE_GENERATOR( GGA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( GGA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )

  
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


UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )


}
}
