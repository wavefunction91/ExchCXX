set( CMAKE_C_COMPILER   clang-15 )
set( CMAKE_CXX_COMPILER clang-15 )

set( CMAKE_CXX_FLAGS " -fsycl -fsycl-unnamed-lambda -fsycl-device-code-split=per_kernel -fsycl-targets=amdgcn-amd-amdhsa -Xsycl-target-backend --offload-arch=gfx906 " )

set( EXCHCXX_ENABLE_SYCL ON CACHE BOOL "" )
