#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX {

BidirectionalMap<std::string, Kernel> kernel_map{
    {{"SlaterExchange", Kernel::SlaterExchange},
     {"PBE_X", Kernel::PBE_X},
     {"PBE_C", Kernel::PBE_C},
     {"LYP", Kernel::LYP},
     {"B3LYP", Kernel::B3LYP},
     {"PBE0", Kernel::PBE0},
     {"VWN3", Kernel::VWN3},
     {"VWN5", Kernel::VWN5},
     {"PZ81", Kernel::PZ81},
     {"PZ81_MOD", Kernel::PZ81_MOD},
     {"PW91_LDA", Kernel::PW91_LDA},
     {"PW91_LDA_MOD", Kernel::PW91_LDA_MOD},
     {"PW91_LDA_RPA", Kernel::PW91_LDA_RPA},
     {"B88", Kernel::B88}}};

std::ostream& operator<<( std::ostream& out, Kernel kern ) {
  out << kernel_map.key(kern);
  return out;
}

XCKernel::XCKernel( 
  const Backend backend, 
  const Kernel kern, 
  const Spin polar) : 
XCKernel( kernel_factory( backend, kern, polar ) ) { }
  

XCKernel::XCKernel( impl_ptr&& ptr ) :
  pimpl_(std::move(ptr)) { }

XCKernel::XCKernel( const XCKernel& other ) :
  pimpl_(other.pimpl_->clone()) { }

XCKernel::XCKernel( XCKernel&& other ) noexcept = default;

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept =default;


XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = XCKernel(other);
}

XCKernel::~XCKernel() noexcept = default;

}
