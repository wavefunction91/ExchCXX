#include <exchcxx/impl/libxc.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

#include <unordered_map>

namespace ExchCXX {


std::unordered_map< XCKernel::Spin, int > libxc_polar_map {
  { XCKernel::Spin::Polarized,   XC_POLARIZED   },
  { XCKernel::Spin::Unpolarized, XC_UNPOLARIZED }
};


/*
template< XCKernel::Kernel Kern >
struct libxc_supports_kernel : public std::false_type { };

template<>
struct libxc_supports_kernel< XCKernel::Kernel::SlaterExchange > : public std::true_type { };
*/



XCKernel libxc_kernel_factory(const std::string& kname, 
  const XCKernel::Spin spin_polar ) {

  return XCKernel( 
    std::make_unique< detail::LibxcKernelImpl >( kname, spin_polar ) );

}


namespace detail {

LibxcKernelImpl::LibxcKernelImpl(
  const std::string&    kname,
  const XCKernel::Spin  spin_polar
) : LibxcKernelImpl( XC_LDA_X, libxc_polar_map[spin_polar] ) { }  

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

/*
LibxcXCKernel::XCKernel( XCKernel&& other ) noexcept :
  XCKernel( other.kernel_, other.polar_, other.initialized_) { 
  other.initialized_ = false; // Avoid double destruction
};

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept {
  kernel_     = other.kernel_;
  polar_      = other.polar_;
  initialized_= other.initialized_;

  other.initialized_ = false; // Avoid double destruction

  return *this;
}

XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = std::move( XCKernel(other) );
}
*/


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
void LibxcKernelImpl::eval_exc_( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc( &kernel_, N, rho, eps );

}


void LibxcKernelImpl::eval_exc_vxc_( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc 
) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc_vxc( &kernel_, N, rho, eps, vxc );

}

// TODO: LDA kxc interfaces

// GGA interface
void LibxcKernelImpl::eval_exc_( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps
) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc( &kernel_, N, rho, sigma, eps );

}


void LibxcKernelImpl::eval_exc_vxc_( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  double*       vrho,
  double*       vsigma
) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc_vxc( &kernel_, N, rho, sigma, eps, vrho, vsigma );

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void LibxcKernelImpl::eval_exc_( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps
) const {

  throw_if_uninitialized();
  assert( is_mgga() );
  xc_mgga_exc( &kernel_, N, rho, sigma, lapl, tau, eps );

}


void LibxcKernelImpl::eval_exc_vxc_( 
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
) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_mgga_exc_vxc( &kernel_, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau );

}

};
};
