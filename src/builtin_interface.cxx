#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/interface.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>


namespace ExchCXX {
namespace detail  {





std::unique_ptr<BuiltinKernel> 
  gen_from_kern( Kernel kern, Spin polar ) {

  if( kern == Kernel::SlaterExchange )
    return std::make_unique<BuiltinSlaterExchange>( polar );
  else if( kern == Kernel::VWN3 )
    return std::make_unique<BuiltinVWN3>( polar );
  else if( kern == Kernel::VWN5 )
    return std::make_unique<BuiltinVWN_RPA>( polar );
  else if( kern == Kernel::PW91_LDA )
    return std::make_unique<BuiltinPW91_LDA>( polar );
  else if( kern == Kernel::PW91_LDA_MOD )
    return std::make_unique<BuiltinPW91_LDA_MOD>( polar );
  else if( kern == Kernel::PW91_LDA_RPA )
    return std::make_unique<BuiltinPW91_LDA_RPA>( polar );
  else if( kern == Kernel::PZ81 )
    return std::make_unique<BuiltinPZ81>( polar );
  else if( kern == Kernel::PZ81_MOD )
    return std::make_unique<BuiltinPZ81_MOD>( polar );


  else if( kern == Kernel::B88 )
    return std::make_unique<BuiltinB88>( polar );
  else if( kern == Kernel::LYP )
    return std::make_unique<BuiltinLYP>( polar );
  else if( kern == Kernel::PBE_X )
    return std::make_unique<BuiltinPBE_X>( polar );
  else if( kern == Kernel::revPBE_X )
    return std::make_unique<BuiltinRevPBE_X>( polar );
  else if( kern == Kernel::PBE_C )
    return std::make_unique<BuiltinPBE_C>( polar );


  else if( kern == Kernel::PBE0 )
    return std::make_unique<BuiltinPBE0>( polar );
  else if( kern == Kernel::B3LYP )
    return std::make_unique<BuiltinB3LYP>( polar );


  else
    throw std::runtime_error("Specified kernel does not have a builtin implementation");


}

BuiltinKernelInterface::~BuiltinKernelInterface() noexcept = default;

BuiltinKernelInterface::BuiltinKernelInterface( Kernel kern, 
  Spin polar ) : 
    whatami_(kern), impl_(gen_from_kern( kern, polar )) { }

BuiltinKernelInterface::BuiltinKernelInterface( 
  const BuiltinKernelInterface& other 
) : BuiltinKernelInterface( other.whatami_, other.impl_->polar() ) { }

std::unique_ptr<XCKernelImpl> BuiltinKernelInterface::clone_() const {
  return std::make_unique<BuiltinKernelInterface>( *this );
}







bool BuiltinKernelInterface::is_lda_()       const noexcept {
  return impl_->is_lda();
}
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

bool BuiltinKernelInterface::supports_inc_interface_() const noexcept {
  return true;
}



#define FORWARD_FOR_BUILTIN( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_DEVICE( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS_DEVICE( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_INC( APPROX, TYPE, func ) \
  FORWARD_XC_INC_ARGS( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_INC_DEVICE( APPROX, TYPE, func ) \
  FORWARD_XC_INC_ARGS_DEVICE( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )







FORWARD_FOR_BUILTIN( LDA,  EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( LDA,  EXC_VXC, eval_exc_vxc )
FORWARD_FOR_BUILTIN( GGA,  EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( GGA,  EXC_VXC, eval_exc_vxc )
FORWARD_FOR_BUILTIN( MGGA, EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( MGGA, EXC_VXC, eval_exc_vxc )

FORWARD_FOR_BUILTIN_INC( LDA,  EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( LDA,  EXC_VXC, eval_exc_vxc_inc )
FORWARD_FOR_BUILTIN_INC( GGA,  EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( GGA,  EXC_VXC, eval_exc_vxc_inc )
FORWARD_FOR_BUILTIN_INC( MGGA, EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( MGGA, EXC_VXC, eval_exc_vxc_inc )


#ifdef EXCHCXX_ENABLE_DEVICE

FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC_VXC, eval_exc_vxc_device )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC_VXC, eval_exc_vxc_device )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_device )

FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  EXC_VXC, eval_exc_vxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  EXC_VXC, eval_exc_vxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_inc_device )

#endif















}
}

