cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( exchcxx_cmake_subproject LANGUAGES C CXX )


if( NOT GIT_REVISION )
  # Get current Git Revision
  find_package( Git REQUIRED )
  execute_process( 
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    OUTPUT_VARIABLE GIT_REVISION
    ERROR_QUIET
  )
  string( STRIP "${GIT_REVISION}" GIT_REVISION )
endif()

message( STATUS "Pulling GIT_REVISION = ${GIT_REVISION}" )

include( FetchContent )
FetchContent_Declare( exchcxx
  GIT_REPOSITORY https://github.com/wavefunction91/exchcxx.git
  GIT_TAG        ${GIT_REVISION}
)
FetchContent_MakeAvailable( exchcxx )

add_executable( exchcxx_link_tester exchcxx_link_tester.cxx )
target_link_libraries( exchcxx_link_tester PUBLIC ExchCXX::ExchCXX )
