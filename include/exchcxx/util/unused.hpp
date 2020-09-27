#include <utility>
#include <exchcxx/util/exchcxx_macros.hpp>

namespace ExchCXX {

inline static void unused( ) { }

template <typename T, typename... Args>
void unused( const T& t, Args&&... args ) {
  (void)(t);
  unused( std::forward<Args>(args)... );
}


template <typename... Args>
inline static void disabled_LDA_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("LDA Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_GGA_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("GGA Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_MGGA_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("MGGA Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_INC_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("INC Interface is disabled for the specified kernel");
}

#ifdef EXCHCXX_ENABLE_DEVICE 
template <typename... Args>
inline static void disabled_LDA_device_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("LDA Device Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_GGA_device_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("GGA Device Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_MGGA_device_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("MGGA Device Interface is disabled for the specified kernel");
}

template <typename... Args>
inline static void disabled_INC_device_interface(Args&&... args) {
  unused( std::forward<Args>(args)... );
  throw std::runtime_error("INC Device Interface is disabled for the specified kernel");
}
#endif

}

#define UNUSED_INTERFACE_GENERATOR( APPROX, TYPE, func, qualifiers )      \
  FORWARD_XC_ARGS( APPROX, TYPE, func, disabled_ ## APPROX ## _interface, \
                   qualifiers )

#define UNUSED_INC_INTERFACE_GENERATOR( APPROX, TYPE, func, qualifiers )      \
  FORWARD_XC_INC_ARGS( APPROX, TYPE, func, disabled_INC_interface, \
                       qualifiers )

#ifdef EXCHCXX_ENABLE_DEVICE 

#define UNUSED_DEVICE_INTERFACE_GENERATOR( APPROX, TYPE, func, qualifiers )      \
  FORWARD_XC_ARGS_DEVICE( APPROX, TYPE, func,                                    \
                          disabled_ ## APPROX ## _device_interface,              \
                          qualifiers )

#define UNUSED_DEVICE_INC_INTERFACE_GENERATOR( APPROX, TYPE, func, qualifiers )  \
  FORWARD_XC_INC_ARGS_DEVICE( APPROX, TYPE, func, disabled_INC_device_interface, \
                              qualifiers )

#endif
