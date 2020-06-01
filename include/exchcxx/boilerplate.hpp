#pragma once

#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

void initialize(XCKernel::Spin);
void finalize();

bool is_initialized();

}
