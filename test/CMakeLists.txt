
find_package( Catch2 CONFIG QUIET )
if( NOT Catch2_FOUND )

  FetchContent_Declare(
    catch2_download
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  )
  if( NOT catch2_download_POPULATED )
    FetchContent_Populate( catch2_download )
    set(CATCH_BUILD_TESTING OFF CACHE BOOL "Build SelfTest project" FORCE)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "Install documentation alongside library" FORCE)
    set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "Install contrib alongside library" FORCE)
    add_subdirectory( ${catch2_download_SOURCE_DIR} ${catch2_download_BINARY_DIR} )
  endif()

  #add_library( catch2 ALIAS Catch2 )
  set( CATCH2_LIBS Catch2 )

else()

  #add_library( catch2 ALIAS Catch2::Catch2 )
  set( CATCH2_LIBS Catch2::Catch2 )

endif()

# Global Catch2 executable
add_library( catch2_main STATIC ut_main.cxx)
target_link_libraries( catch2_main PUBLIC ${CATCH2_LIBS} )

#add_executable( boilerplate_test boilerplate_test.cxx )
#target_link_libraries( boilerplate_test PUBLIC exchcxx catch2_main )
#target_compile_features( boilerplate_test PRIVATE cxx_std_14 )

add_executable( xc_kernel_test xc_kernel_test.cxx reference_values.cxx )
target_link_libraries( xc_kernel_test PUBLIC exchcxx catch2_main )
target_compile_features( xc_kernel_test PRIVATE cxx_std_14 )

add_executable( xc_functional_test xc_functional_test.cxx )
target_link_libraries( xc_functional_test PUBLIC exchcxx catch2_main )
target_compile_features( xc_functional_test PRIVATE cxx_std_14 )

add_test( NAME BOILERPLATE   COMMAND boilerplate_test   )
add_test( NAME XC_KERNEL     COMMAND xc_kernel_test     )
add_test( NAME XC_FUNCTIONAL COMMAND xc_functional_test )
