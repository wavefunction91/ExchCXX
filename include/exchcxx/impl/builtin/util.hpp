#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

namespace ExchCXX {

inline static void disabled_lda_interface() {
  throw std::runtime_error("LDA Interface is disabled for the specified kernel");
}

inline static void disabled_gga_interface() {
  throw std::runtime_error("GGA Interface is disabled for the specified kernel");
}

inline static void disabled_mgga_interface() {
  throw std::runtime_error("MGGA Interface is disabled for the specified kernel");
}

inline static void disabled_inc_interface() {
  throw std::runtime_error("Scale and Increment Interface is disabled for the specified kernel / backend");
}

#ifdef EXCHCXX_ENABLE_DEVICE 
inline static void disabled_lda_device_interface() {
  throw std::runtime_error("LDA Device Interface is disabled for the specified kernel");
}

inline static void disabled_gga_device_interface() {
  throw std::runtime_error("GGA Device Interface is disabled for the specified kernel");
}

inline static void disabled_mgga_device_interface() {
  throw std::runtime_error("MGGA Device Interface is disabled for the specified kernel");
}

inline static void disabled_inc_device_interface() {
  throw std::runtime_error("Scale and Increment Interface is disabled for the specified kernel / backend");
}
#endif

#ifdef __CUDACC__

#define BUILTIN_KERNEL_EVAL_RETURN static inline constexpr void __host__ __device__

#else

#define BUILTIN_KERNEL_EVAL_RETURN static inline constexpr void 

#endif

}
