#define CATCH_CONFIG_MAIN
#include <spheroidal/grid_functions.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "Test Name", "[Will Pass]" )
{
    REQUIRE( 4 / 2 == 2 );
}

TEST_CASE( "Test Name 2", "[Will Fail]" )
{
    REQUIRE( 5 / 2 == 2 );
}