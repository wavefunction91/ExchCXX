#include <exchcxx/impl/builtin/kernel.hpp>
#include <exchcxx/impl/builtin/interface.hpp>

#include <exchcxx/impl/builtin.hpp>

namespace ExchCXX {
namespace detail  {





std::unique_ptr<BuiltinKernel> 
  gen_from_kern( XCKernel::Kernel kern, XCKernel::Spin polar ) {

  if( kern == XCKernel::Kernel::SlaterExchange )
    return std::make_unique<BuiltinSlaterExchange>( polar );
  else if( kern == XCKernel::Kernel::LYP )
    return std::make_unique<BuiltinLYP>( polar );
  else if( kern == XCKernel::Kernel::PBE_X )
    return std::make_unique<BuiltinPBE_X>( polar );
  else if( kern == XCKernel::Kernel::PBE_C )
    return std::make_unique<BuiltinPBE_C>( polar );
  else
    throw std::runtime_error("Specified kernel does not have a builtin implementation");


}

BuiltinKernelInterface::~BuiltinKernelInterface() noexcept = default;

BuiltinKernelInterface::BuiltinKernelInterface( XCKernel::Kernel kern, 
  XCKernel::Spin polar ) : 
    whatami_(kern), impl_(gen_from_kern( kern, polar )) { }

BuiltinKernelInterface::BuiltinKernelInterface( 
  const BuiltinKernelInterface& other 
) : BuiltinKernelInterface( other.whatami_, other.impl_->polar() ) { }

std::unique_ptr<XCKernelImpl> BuiltinKernelInterface::clone_() const {
  return std::make_unique<BuiltinKernelInterface>( *this );
}







bool BuiltinKernelInterface::is_lda_()       const noexcept {
  return impl_->is_lda();
};
bool BuiltinKernelInterface::is_gga_()       const noexcept {
  return impl_->is_gga();
}
bool BuiltinKernelInterface::is_mgga_()      const noexcept {
  return impl_->is_mgga();
}
bool BuiltinKernelInterface::is_hyb_()       const noexcept {
  return impl_->is_hyb();
}
bool BuiltinKernelInterface::is_polarized_() const noexcept {
  return impl_->is_polarized();
}
double BuiltinKernelInterface::hyb_exx_()    const noexcept {
  return impl_->hyb_exx();
}



#define FORWARD_FOR_BUILTIN( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const );
#define FORWARD_FOR_BUILTIN_DEVICE( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS_DEVICE( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const );







FORWARD_FOR_BUILTIN( LDA,  EXC,     eval_exc     );
FORWARD_FOR_BUILTIN( LDA,  EXC_VXC, eval_exc_vxc );
FORWARD_FOR_BUILTIN( GGA,  EXC,     eval_exc     );
FORWARD_FOR_BUILTIN( GGA,  EXC_VXC, eval_exc_vxc );
FORWARD_FOR_BUILTIN( MGGA, EXC,     eval_exc     );
FORWARD_FOR_BUILTIN( MGGA, EXC_VXC, eval_exc_vxc );


#ifdef EXCHCXX_ENABLE_DEVICE

FORWARD_FOR_BUILTIN( LDA,  EXC,     eval_exc_device     );
FORWARD_FOR_BUILTIN( LDA,  EXC_VXC, eval_exc_vxc_device );
FORWARD_FOR_BUILTIN( GGA,  EXC,     eval_exc_device     );
FORWARD_FOR_BUILTIN( GGA,  EXC_VXC, eval_exc_vxc_device );
FORWARD_FOR_BUILTIN( MGGA, EXC,     eval_exc_device     );
FORWARD_FOR_BUILTIN( MGGA, EXC_VXC, eval_exc_vxc_device );

FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC,     eval_exc_device_async     );
FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC_VXC, eval_exc_vxc_device_async );
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC,     eval_exc_device_async     );
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC_VXC, eval_exc_vxc_device_async );
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC,     eval_exc_device_async     );
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_device_async );


#endif















}
}

