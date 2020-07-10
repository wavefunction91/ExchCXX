#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif


namespace ExchCXX {

class BuiltinKernel;

struct BuiltinSlaterExchange;
struct BuiltinVWN3;
struct BuiltinVWN_RPA;

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
