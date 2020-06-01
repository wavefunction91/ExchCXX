#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/builtin.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

#include <iostream>

namespace ExchCXX {

XCKernel builtin_kernel_factory( XCKernel::Kernel kernel, XCKernel::Spin polar) {

    return XCKernel(std::make_unique< detail::BuiltinKernelInterface >( kernel, polar ));
  //return XCKernel(
  //  std::make_unique< detail::BuiltinKernelInterface >( kernel, polar );
  //);

}


namespace detail {


std::unique_ptr<BuiltinKernel> gen_from_kern( XCKernel::Kernel kern, XCKernel::Spin polar ) {

  if( kern == XCKernel::Kernel::SlaterExchange )
    return std::make_unique<BuiltinSlaterExchange>( polar );
  else
    throw std::runtime_error("Specified kernel does not have a builtin implementation");


}

BuiltinKernelInterface::~BuiltinKernelInterface() noexcept {

}

BuiltinKernelInterface::BuiltinKernelInterface( XCKernel::Kernel kern, 
  XCKernel::Spin polar ) : 
    whatami_(kern), impl_(gen_from_kern( kern, polar )) { }

BuiltinKernelInterface::BuiltinKernelInterface( const BuiltinKernelInterface& other ) :
  BuiltinKernelInterface( other.whatami_, other.impl_->polar() ) { }

std::unique_ptr<XCKernelImpl> BuiltinKernelInterface::clone_() const {
  return std::make_unique<BuiltinKernelInterface>( *this );
}

}

};
