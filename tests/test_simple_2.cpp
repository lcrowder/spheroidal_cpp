#include <spheroidal/grid_functions.hpp>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "Another Test", "[Will Pass]" )
{
    REQUIRE( 6 / 2 == 3 );

}

TEST_CASE( "Another Test 2", "[Will Fail]" )
{
    REQUIRE( 6 / 2 == 2 );

}