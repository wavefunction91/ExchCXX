find_package( Catch2 REQUIRED )

# Global Catch2 executable
add_library( catch2_main STATIC ut_main.cxx)
target_link_libraries( catch2_main PUBLIC Catch2::Catch2 )

add_executable( xc_kernel_test xc_kernel_test.cxx )
target_link_libraries( xc_kernel_test PUBLIC exchcxx catch2_main )
target_compile_features( xc_kernel_test PRIVATE cxx_std_17 )

add_executable( xc_functional_test xc_functional_test.cxx )
target_link_libraries( xc_functional_test PUBLIC exchcxx catch2_main )
target_compile_features( xc_functional_test PRIVATE cxx_std_17 )
