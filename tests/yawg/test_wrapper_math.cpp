#define CATCH_CONFIG_MAIN
#include <yawg/core.h>
#include <yawg/utils.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

// Use Catch2 to test the imag, real, abs, abs2, and arg methods of gsl::complex
TEST_CASE("gsl::complex methods", "[gsl::complex]")
{
    gsl::complex c1(1.0, 2.0);
    REQUIRE(c1.real() == 1.0);
    REQUIRE(c1.imag() == 2.0);
    REQUIRE_THAT(c1.abs(), Catch::Matchers::WithinAbs(std::sqrt(5.0), 1e-12));
    REQUIRE_THAT(c1.abs2(), Catch::Matchers::WithinAbs(5.0, 1e-12));
    REQUIRE_THAT(c1.arg(), Catch::Matchers::WithinAbs(std::atan2(2.0, 1.0), 1e-12));
}

// Use Catch2 to test the += -= *= /= operators of gsl::complex
TEST_CASE("gsl::complex assignment operators", "[gsl::complex]")
{
    gsl::complex c1(1.0, 2.0), c2(3.0, 4.0), c3(5.0, 6.0), c4(7.0, 8.0);

    // Test +=
    SECTION("Test += operator")
    {
        c1 += c2;
        REQUIRE(c1.dat[0] == 4.0);
        REQUIRE(c1.dat[1] == 6.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test -=
    SECTION("Test -= operator")
    {
        c1 -= c2;
        REQUIRE(c1.dat[0] == -2.0);
        REQUIRE(c1.dat[1] == -2.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test *=
    SECTION("Test *= operator")
    {
        c1 *= c2;
        REQUIRE(c1.dat[0] == -5.0);
        REQUIRE(c1.dat[1] == 10.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test /=
    SECTION("Test /= operator")
    {
        c1 /= c2;
        REQUIRE_THAT(c1.dat[0], Catch::Matchers::WithinAbs(0.44, 1e-12));
        REQUIRE_THAT(c1.dat[1], Catch::Matchers::WithinAbs(0.08, 1e-12));
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test =
    SECTION("Test = operator")
    {
        c1 = c2;
        REQUIRE(c1.dat[0] == 3.0);
        REQUIRE(c1.dat[1] == 4.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test nested assignment
    SECTION("Test nested = operator")
    {
        c1 = c2 = c3;
        REQUIRE(c1.dat[0] == 5.0);
        REQUIRE(c1.dat[1] == 6.0);
        REQUIRE(c2.dat[0] == 5.0);
        REQUIRE(c2.dat[1] == 6.0);
        REQUIRE(c3.dat[0] == 5.0);
        REQUIRE(c3.dat[1] == 6.0);
    }

    // Test literal operator
    SECTION("Test literal operator")
    {
        using namespace gsl::complex_literals;
        c1 = 1.0;
        REQUIRE(c1.dat[0] == 1.0);
        REQUIRE(c1.dat[1] == 0.0);

        c2 = 2.0_i;
        REQUIRE(c2.dat[0] == 0.0);
        REQUIRE(c2.dat[1] == 2.0);

        c3 = 3.0 + 4.0_i;
        REQUIRE(c3.dat[0] == 3.0);
        REQUIRE(c3.dat[1] == 4.0);

        c4 = 5.0 - 6.0_i;
        REQUIRE(c4.dat[0] == 5.0);
        REQUIRE(c4.dat[1] == -6.0);
    }
}

// Use Catch2 to test the +, -, *, /, and == operators of gsl::complex
TEST_CASE("gsl::complex operators", "[gsl::complex]")
{
    gsl::complex c1(1.0, 2.0), c2(3.0, 4.0);

    // Test +
    SECTION("Test + operator")
    {
        gsl::complex c3 = c1 + c2;
        REQUIRE(c3.dat[0] == 4.0);
        REQUIRE(c3.dat[1] == 6.0);
        REQUIRE(c1.dat[0] == 1.0);
        REQUIRE(c1.dat[1] == 2.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test -
    SECTION("Test - operator")
    {
        gsl::complex c3 = c1 - c2;
        REQUIRE(c3.dat[0] == -2.0);
        REQUIRE(c3.dat[1] == -2.0);
        REQUIRE(c1.dat[0] == 1.0);
        REQUIRE(c1.dat[1] == 2.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test *
    SECTION("Test * operator")
    {
        gsl::complex c3 = c1 * c2;
        REQUIRE(c3.dat[0] == -5.0);
        REQUIRE(c3.dat[1] == 10.0);
        REQUIRE(c1.dat[0] == 1.0);
        REQUIRE(c1.dat[1] == 2.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test /
    SECTION("Test / operator")
    {
        gsl::complex c3 = c1 / c2;
        REQUIRE_THAT(c3.dat[0], Catch::Matchers::WithinAbs(0.44, 1e-12));
        REQUIRE_THAT(c3.dat[1], Catch::Matchers::WithinAbs(0.08, 1e-12));
        REQUIRE(c1.dat[0] == 1.0);
        REQUIRE(c1.dat[1] == 2.0);
        REQUIRE(c2.dat[0] == 3.0);
        REQUIRE(c2.dat[1] == 4.0);
    }

    // Test ==s
    SECTION("Test == operator")
    {
        using namespace gsl::complex_literals;
        gsl::complex c3(1.0, 2.0);
        REQUIRE(c1 == c3);
        REQUIRE_FALSE(c1 == c2);

        REQUIRE((1.0 + 2.0_i) == gsl::complex(1.0, 2.0));
        REQUIRE_FALSE((1.0 + 2.0_i) == gsl::complex(1.0, 3.0));
        REQUIRE(1.0 == gsl::complex(1.0, 0.0));
        REQUIRE_FALSE(1.0 == gsl::complex(1.0, 1.0));
        REQUIRE(1.0_i == gsl::complex(0.0, 1.0));
        REQUIRE_FALSE(1.0_i == gsl::complex(1.0, 1.0));
    }
}

//! Use Catch2 to test addition of gsl::vectors
TEST_CASE("gsl::vector addition", "[gsl::vector]")
{
    gsl::vector v1( gsl::arange(1.0, 4.0) );
    gsl::vector v2( gsl::arange(4.0, 7.0) );

    SECTION("Test lvalue + lvalue")
    {
        gsl::vector v3 = v1 + v2;
        REQUIRE( v3.get(0) == 5.0 );
        REQUIRE( v3.get(1) == 7.0 );
        REQUIRE( v3.get(2) == 9.0 );
    }

    SECTION("Test lvalue + rvalue")
    {
        gsl::vector v3 = v1 + gsl::vector( gsl::arange(4.0, 7.0) );
        REQUIRE( v3.get(0) == 5.0 );
        REQUIRE( v3.get(1) == 7.0 );
        REQUIRE( v3.get(2) == 9.0 );
    }

    SECTION("Test rvalue + lvalue")
    {
        gsl::vector v3 = gsl::vector( gsl::arange(1.0, 4.0) ) + v2;
        REQUIRE( v3.get(0) == 5.0 );
        REQUIRE( v3.get(1) == 7.0 );
        REQUIRE( v3.get(2) == 9.0 );
    }

    SECTION("Test rvalue + rvalue")
    {
        gsl::vector v3 = gsl::vector( gsl::arange(1.0, 4.0) ) + gsl::vector( gsl::arange(4.0, 7.0) );
        REQUIRE( v3.get(0) == 5.0 );
        REQUIRE( v3.get(1) == 7.0 );
        REQUIRE( v3.get(2) == 9.0 );
    }
}

//! Use Catch2 to test addition of gsl::vectors and gsl::cvectors
TEST_CASE("gsl::vector and gsl::cvector addition", "[gsl::vector][gsl::cvector]")
{
    using namespace gsl::complex_literals;

    gsl::vector v1( gsl::arange(1.0, 4.0) );
    gsl::cvector v2( 1.0_i * gsl::arange(4.0, 7.0) );

    SECTION("Test addition")
    {
        gsl::cvector v3 = v1 + v2;
        REQUIRE( v3.get(0) == gsl::complex(1.0, 4.0) );
        REQUIRE( v3.get(1) == gsl::complex(2.0, 5.0) );
        REQUIRE( v3.get(2) == gsl::complex(3.0, 6.0) );
    }
}

//! Use Catch2 to test scalar multiplication of gsl::vectors and gsl::cvectors
TEST_CASE("gsl::vector and gsl::cvector scalar multiplication", "[gsl::vector][gsl::cvector]")
{
    using namespace gsl::complex_literals;

    gsl::vector v1( 2.0 * gsl::arange(1.0, 4.0) );
    REQUIRE( v1.get(0) == 2.0 );
    REQUIRE( v1.get(1) == 4.0 );
    REQUIRE( v1.get(2) == 6.0 );   

    gsl::cvector v2( 2.0 * gsl::cvector( gsl::arange(1.0, 4.0) ) ); 
    REQUIRE( v2.get(0) == gsl::complex(2.0, 0.0) );
    REQUIRE( v2.get(1) == gsl::complex(4.0, 0.0) );
    REQUIRE( v2.get(2) == gsl::complex(6.0, 0.0) );

    gsl::cvector v3( 1.0_i * gsl::cvector( gsl::arange(1.0, 4.0) ) );
    REQUIRE( v3.get(0) == gsl::complex(0.0, 1.0) );
    REQUIRE( v3.get(1) == gsl::complex(0.0, 2.0) );
    REQUIRE( v3.get(2) == gsl::complex(0.0, 3.0) );

    gsl::cvector v4( 1.0_i * gsl::arange(1.0, 4.0) );
    REQUIRE( v4.get(0) == gsl::complex(0.0, 1.0) );
    REQUIRE( v4.get(1) == gsl::complex(0.0, 2.0) );
    REQUIRE( v4.get(2) == gsl::complex(0.0, 3.0) );
}

//! Use Catch2 to test +=, -=, *=, /= of gsl::vectors
TEST_CASE("gsl::vector +=, -=, *=, /=", "[gsl::vector]")
{
    gsl::vector v1( gsl::arange(1.0, 4.0) );
    gsl::vector v2( gsl::arange(4.0, 7.0) );
    double x = 1.5;

    SECTION("Test += lvalue")
    {
        v1 += v2;
        REQUIRE( v1.get(0) == 5.0 );
        REQUIRE( v1.get(1) == 7.0 );
        REQUIRE( v1.get(2) == 9.0 );
    }

    SECTION("Test += rvalue")
    {
        v1 += gsl::vector( gsl::arange(4.0, 7.0) );
        REQUIRE( v1.get(0) == 5.0 );
        REQUIRE( v1.get(1) == 7.0 );
        REQUIRE( v1.get(2) == 9.0 );
    }

    SECTION("Test -= lvalue")
    {
        v1 -= v2;
        REQUIRE( v1.get(0) == -3.0 );
        REQUIRE( v1.get(1) == -3.0 );
        REQUIRE( v1.get(2) == -3.0 );
    }

    SECTION("Test -= rvalue")
    {
        v1 -= gsl::vector( gsl::arange(4.0, 7.0) );
        REQUIRE( v1.get(0) == -3.0 );
        REQUIRE( v1.get(1) == -3.0 );
        REQUIRE( v1.get(2) == -3.0 );
    }

    SECTION("Test *=")
    {
        v1 *= x;
        REQUIRE( v1.get(0) == 1.5 );
        REQUIRE( v1.get(1) == 3.0 );
        REQUIRE( v1.get(2) == 4.5 );
    }

    SECTION("Test /=")
    {
        v1 /= x;
        REQUIRE_THAT(v1.get(0), Catch::Matchers::WithinAbs(2.0/3.0, 1e-12));
        REQUIRE_THAT(v1.get(1), Catch::Matchers::WithinAbs(4.0/3.0, 1e-12));
        REQUIRE_THAT(v1.get(2), Catch::Matchers::WithinAbs(6.0/3.0, 1e-12));
    }
}