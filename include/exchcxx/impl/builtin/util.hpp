#pragma once
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


}
