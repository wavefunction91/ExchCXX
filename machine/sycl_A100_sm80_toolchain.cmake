set( CMAKE_C_COMPILER   clang )
set( CMAKE_CXX_COMPILER clang++ )

set( CMAKE_CXX_FLAGS " -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend '--cuda-gpu-arch=sm_80' " )

set( EXCHCXX_ENABLE_SYCL ON CACHE BOOL "" )
