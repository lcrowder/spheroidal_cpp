#define CATCH_CONFIG_MAIN
#include <spheroidal/grid_functions.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "DING", "[Will Pass]" )
{
    REQUIRE( 6 / 2 == 3 );

}

TEST_CASE( "DONG", "[Will Fail]" )
{
    REQUIRE( 6 / 2 == 2 );

}