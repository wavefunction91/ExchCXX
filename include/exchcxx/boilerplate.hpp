#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

void initialize(XCKernel::Spin);
void finalize();

bool is_initialized();

}
