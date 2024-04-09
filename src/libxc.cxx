/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#include "libxc_common.hpp"

#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/util/unused.hpp>
#include <exchcxx/exceptions/exchcxx_exception.hpp>
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

  // MGGA Functionals
  { Kernel::SCAN_C,         XC_MGGA_X_SCAN      },
  { Kernel::SCAN_X,         XC_MGGA_C_SCAN      },
  { Kernel::R2SCAN_C,       XC_MGGA_X_R2SCAN    },
  { Kernel::R2SCAN_X,       XC_MGGA_C_R2SCAN    },
  { Kernel::R2SCANL_C,      XC_MGGA_X_R2SCANL   },
  { Kernel::R2SCANL_X,      XC_MGGA_C_R2SCANL   },
  { Kernel::FT98_X,         XC_MGGA_X_FT98      },

  // KEDFs
  { Kernel::PC07_K,         XC_MGGA_K_PC07      },
  { Kernel::PC07OPT_K,      XC_MGGA_K_PC07_OPT  },

  // Hybrid GGA Functionals
  { Kernel::B3LYP,          XC_HYB_GGA_XC_B3LYP },
  { Kernel::PBE0,           XC_HYB_GGA_XC_PBEH  },
};


inline bool libxc_supports_kernel( const Kernel kern ) {
  return libxc_kernel_map.find(kern) != libxc_kernel_map.end();
}


inline bool libxc_supports_functional( const std::string xc_name ) {
  return xc_functional_get_number(xc_name.c_str()) != -1;
}


XCKernel libxc_kernel_factory(const Kernel kern, 
  const Spin spin_polar ) {

  EXCHCXX_BOOL_CHECK( "KERNEL NYI FOR Libxc BACKEND", libxc_supports_kernel(kern) )

  return XCKernel( 
    std::make_unique< detail::LibxcKernelImpl >( kern, spin_polar ) );

}

XCKernel libxc_kernel_factory( const std::string xc_name,
  const Spin spin_polar ) {

  EXCHCXX_BOOL_CHECK( "LibXC FUNCTIONAL NOT FOUND", libxc_supports_functional(xc_name) )

  return XCKernel(
      std::make_unique< detail::LibxcKernelImpl >( xc_name, spin_polar ) );

}

namespace detail {

LibxcKernelImpl::LibxcKernelImpl(
  const Kernel    kern,
  const Spin  spin_polar
) : LibxcKernelImpl( libxc_kernel_map[kern], libxc_polar_map[spin_polar] ) { }  

LibxcKernelImpl::LibxcKernelImpl(
  const std::string xc_name,
  const Spin  spin_polar
) : LibxcKernelImpl( xc_functional_get_number(xc_name.c_str()), libxc_polar_map[spin_polar] ) { }  

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
    (kernel_.info->family == XC_FAMILY_GGA    ) 
#if XC_MAJOR_VERSION < 7
    or (kernel_.info->family == XC_FAMILY_HYB_GGA)
#endif
    ;
}

bool LibxcKernelImpl::is_mgga_() const noexcept {
  return 
    (kernel_.info->family == XC_FAMILY_MGGA    )
#if XC_MAJOR_VERSION < 7
    or (kernel_.info->family == XC_FAMILY_HYB_MGGA)
#endif
    ;
}

bool LibxcKernelImpl::is_hyb_() const noexcept {
#if XC_MAJOR_VERSION < 7
  return
    (kernel_.info->family == XC_FAMILY_HYB_GGA ) or
    (kernel_.info->family == XC_FAMILY_HYB_MGGA);
#else
  return xc_hyb_type(&kernel_) == XC_HYB_HYBRID;
#endif
}

bool LibxcKernelImpl::needs_laplacian_() const noexcept {
  return kernel_.info->flags & XC_FLAGS_NEEDS_LAPLACIAN;
}

bool LibxcKernelImpl::is_polarized_() const noexcept {
  return polar_ == XC_POLARIZED;
}


double LibxcKernelImpl::hyb_exx_() const noexcept {
  return is_hyb_() ? xc_hyb_exx_coef(&kernel_) : 0.0;
}

bool LibxcKernelImpl::supports_inc_interface_() const noexcept {
  return false;
}


// LDA interfaces
LDA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );
  xc_lda_exc( &kernel_, N, rho, eps );

}


LDA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );
  xc_lda_exc_vxc( &kernel_, N, rho, eps, vxc );

}


UNUSED_INC_INTERFACE_GENERATOR( LDA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( LDA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )


// GGA interface
GGA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );
  xc_gga_exc( &kernel_, N, rho, sigma, eps );

}


GGA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );
  xc_gga_exc_vxc( &kernel_, N, rho, sigma, eps, vrho, vsigma );

}

UNUSED_INC_INTERFACE_GENERATOR( GGA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( GGA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )

  
// mGGA interface
MGGA_EXC_GENERATOR( LibxcKernelImpl::eval_exc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );
  xc_mgga_exc( &kernel_, N, rho, sigma, lapl, tau, eps );

}


MGGA_EXC_VXC_GENERATOR( LibxcKernelImpl::eval_exc_vxc_ ) const {

  throw_if_uninitialized();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );
  xc_mgga_exc_vxc( &kernel_, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau );

}


UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC, LibxcKernelImpl::eval_exc_inc_,     
                                const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, LibxcKernelImpl::eval_exc_vxc_inc_, 
                                const )


}
}
