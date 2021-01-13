cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( exchcxx_cmake_discovery LANGUAGES C CXX )

find_package( ExchCXX REQUIRED )
add_executable( exchcxx_link_tester exchcxx_link_tester.cxx )
target_link_libraries( exchcxx_link_tester PUBLIC ExchCXX::ExchCXX )
