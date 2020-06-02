#pragma once

namespace ExchCXX {

class BuiltinKernel;

struct BuiltinSlaterExchange;
struct BuiltinLYP;
struct BuiltinPBE_X;
struct BuiltinPBE_C;

struct BuiltinPBE0;


template <typename KernelType>
struct kernel_traits;




namespace detail {
  template <typename KernelType, typename = void>
  struct BuiltinKernelImpl;

  class BuiltinKernelInterface;
}

}
