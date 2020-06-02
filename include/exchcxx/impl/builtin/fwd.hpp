#pragma once

namespace ExchCXX {

class BuiltinKernel;

template <typename KernelType>
struct kernel_traits;




namespace detail {
  template <typename KernelType, typename = void>
  struct BuiltinKernelImpl;

  class BuiltinKernelInterface;
}

}
