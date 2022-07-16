set( CMAKE_C_COMPILER   icx  )
set( CMAKE_CXX_COMPILER icpx )

set( CMAKE_CXX_FLAGS " -fsycl -fno-sycl-id-queries-fit-in-int -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.1.0,12.4.0\" " )

set( EXCHCXX_ENABLE_SYCL ON CACHE BOOL "" )
