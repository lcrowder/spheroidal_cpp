#define CATCH_CONFIG_MAIN
#include <spheroidal/grid_functions.hpp>
#include <catch2/catch.hpp>

TEST_CASE( "Test Name", "[Will Pass]" )
{
    REQUIRE( 4 / 2 == 2 );

}

TEST_CASE( "Test Name 2", "[Will Fail]" )
{
    REQUIRE( 5 / 2 == 2 );

}