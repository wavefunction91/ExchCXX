cmake_minimum_required( VERSION 3.17 FATAL_ERROR) # Require CMake 3.17+

# Set up project definition + version information
project( ExchCXX VERSION 0.0.1 LANGUAGES C CXX )

# ExchCXX Options
option( EXCHCXX_ENABLE_TESTS "Enable Unit Tests" ON         )
option( EXCHCXX_ENABLE_DEVICE "Enable Device Code (GPU)" ON )

if( EXCHCXX_ENABLE_DEVICE )
  enable_language( CUDA )
endif()

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()


## Find LibXC
find_package( Libxc REQUIRED )
        
add_subdirectory( src  )

if( EXCHCXX_ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif()
