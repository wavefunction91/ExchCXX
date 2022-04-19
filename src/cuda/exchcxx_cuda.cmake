set( EXCHCXX_CUDA_SOURCES 
  cuda/xc_functional_device.cu 
  cuda/libxc_device.cxx 
  cuda/builtin.cu 
)

find_package( CUDAToolkit REQUIRED )
#add_library( exchcxx_device OBJECT ${EXCHCXX_CUDA_SOURCES} )
#target_compile_features( exchcxx_device PUBLIC cuda_std_14 cxx_std_14 )
#target_compile_definitions( exchcxx_device PUBLIC "EXCHCXX_HAS_CONFIG_H=1" )
#target_link_libraries( exchcxx_device PUBLIC Libxc::xc CUDA::cudart )
#
#target_include_directories( exchcxx_device
#  PRIVATE
#    $<INSTALL_INTERFACE:include>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
#)
#
#target_compile_options( exchcxx_device 
#  PRIVATE
#    $<$<COMPILE_LANGUAGE:CXX>: -Wall -Wextra -Wpedantic -Wnon-virtual-dtor>
#    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
#)

target_sources( exchcxx PRIVATE ${EXCHCXX_CUDA_SOURCES} )
target_link_libraries( exchcxx PUBLIC CUDA::cudart )
target_compile_features( exchcxx PRIVATE cuda_std_14 )
target_compile_options( exchcxx
  PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>: -Xcudafe --diag_suppress=partial_override -Xptxas -v > 
)
