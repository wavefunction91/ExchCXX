#pragma once
#include <cuda_runtime.h>

namespace ExchCXX {
namespace device {

  struct cuda_stream_t {
    cudaStream_t* stream;
  };

}
};
