file(GLOB EXCHCXX_CUDA_SOURCES "cuda/*.cu")
if( EXCHCXX_ENABLE_LIBXC )
  list(APPEND EXCHCXX_CUDA_SOURCES cuda/libxc_device.cxx)
endif()

find_package( CUDAToolkit REQUIRED )

target_sources( exchcxx PRIVATE ${EXCHCXX_CUDA_SOURCES} )
target_link_libraries( exchcxx PUBLIC CUDA::cudart )
target_compile_features( exchcxx PRIVATE cuda_std_14 )
target_compile_options( exchcxx
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
)
