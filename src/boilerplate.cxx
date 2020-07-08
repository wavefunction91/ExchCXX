#include <exchcxx/boilerplate.hpp>

namespace ExchCXX {

static bool initialized_ = false;

void initialize(Spin polar) {
  initialized_ = true;
}

void finalize() {
  initialized_ = false;
}

bool is_initialized(){ return initialized_; }

}
