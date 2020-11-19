#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX {

KerMap str2ker;

std::ostream& operator<<( std::ostream& out, Kernel kern ) {
  out << str2ker.to_str(kern);
  return out;
}

XCKernel::XCKernel( 
  const Backend backend, 
  const Kernel kern, 
  const Spin polar) : 
XCKernel( std::move(kernel_factory( backend, kern, polar )) ) { }
  

XCKernel::XCKernel( impl_ptr&& ptr ) :
  pimpl_(std::move(ptr)) { }

XCKernel::XCKernel( const XCKernel& other ) :
  pimpl_(other.pimpl_->clone()) { }

XCKernel::XCKernel( XCKernel&& other ) noexcept = default;

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept =default;


XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = std::move( XCKernel(other) );
}

XCKernel::~XCKernel() noexcept = default;

}
