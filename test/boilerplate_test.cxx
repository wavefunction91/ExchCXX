#include "catch2/catch.hpp"
#include <exchcxx/boilerplate.hpp>

TEST_CASE( "ExchCXX Initialize / Finalize" ) {

  // Default is not initialized
  REQUIRE( !ExchCXX::is_initialized() );

  // Make sure the initialization is toggled
  ExchCXX::initialize(ExchCXX::XCKernel::Spin::Unpolarized);
  REQUIRE( ExchCXX::is_initialized() );

  // Make sure the initialiazation is reversable
  ExchCXX::finalize();
  REQUIRE( !ExchCXX::is_initialized() );

}
