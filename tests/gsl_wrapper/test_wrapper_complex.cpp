#define CATCH_CONFIG_MAIN
#include <gsl_wrapper/core.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

// Use Catch2 to test the various constructors for gsl::complex
TEST_CASE("gsl::complex constructors", "[gsl::complex]")
{
    // Default constructor
    gsl::complex c1;
    REQUIRE(c1.get_gsl_data().dat[0] == 0.0);
    REQUIRE(c1.get_gsl_data().dat[1] == 0.0);

    // Constructor with real and imaginary parts
    gsl::complex c2(1.0, 2.0);
    REQUIRE(c2.get_gsl_data().dat[0] == 1.0);
    REQUIRE(c2.get_gsl_data().dat[1] == 2.0);

    // Constructor with one parameter
    gsl::complex c3(1.0);
    REQUIRE(c3.get_gsl_data().dat[0] == 1.0);
    REQUIRE(c3.get_gsl_data().dat[1] == 0.0);

    // Copy constructor
    gsl::complex c4(c3);
    REQUIRE(c3.get_gsl_data().dat[0] == 1.0);
    REQUIRE(c3.get_gsl_data().dat[1] == 0.0);
    REQUIRE(c4.get_gsl_data().dat[0] == 1.0);
    REQUIRE(c4.get_gsl_data().dat[1] == 0.0);

    // Test that copy is deep
    c3.set(2.0, 3.0);
    REQUIRE(c3.get_gsl_data().dat[0] == 2.0);
    REQUIRE(c3.get_gsl_data().dat[1] == 3.0);
    REQUIRE(c4.get_gsl_data().dat[0] == 1.0);
    REQUIRE(c4.get_gsl_data().dat[1] == 0.0);
}

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
        REQUIRE(c1.get_gsl_data().dat[0] == 4.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 6.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test -=
    SECTION("Test -= operator")
    {
        c1 -= c2;
        REQUIRE(c1.get_gsl_data().dat[0] == -2.0);
        REQUIRE(c1.get_gsl_data().dat[1] == -2.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test *=
    SECTION("Test *= operator")
    {
        c1 *= c2;
        REQUIRE(c1.get_gsl_data().dat[0] == -5.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 10.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test /=
    SECTION("Test /= operator")
    {
        c1 /= c2;
        REQUIRE_THAT(c1.get_gsl_data().dat[0], Catch::Matchers::WithinAbs(0.44, 1e-12));
        REQUIRE_THAT(c1.get_gsl_data().dat[1], Catch::Matchers::WithinAbs(0.08, 1e-12));
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test =
    SECTION("Test = operator")
    {
        c1 = c2;
        REQUIRE(c1.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 4.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test nested assignment
    SECTION("Test nested = operator")
    {
        c1 = c2 = c3;
        REQUIRE(c1.get_gsl_data().dat[0] == 5.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 6.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 5.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 6.0);
        REQUIRE(c3.get_gsl_data().dat[0] == 5.0);
        REQUIRE(c3.get_gsl_data().dat[1] == 6.0);
    }

    // Test literal operator
    SECTION("Test literal operator")
    {
        using namespace gsl::complex_literals;
        c1 = 1.0;
        REQUIRE(c1.get_gsl_data().dat[0] == 1.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 0.0);

        c2 = 2.0_i;
        REQUIRE(c2.get_gsl_data().dat[0] == 0.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 2.0);

        c3 = 3.0 + 4.0_i;
        REQUIRE(c3.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c3.get_gsl_data().dat[1] == 4.0);

        c4 = 5.0 - 6.0_i;
        REQUIRE(c4.get_gsl_data().dat[0] == 5.0);
        REQUIRE(c4.get_gsl_data().dat[1] == -6.0);
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
        REQUIRE(c3.get_gsl_data().dat[0] == 4.0);
        REQUIRE(c3.get_gsl_data().dat[1] == 6.0);
        REQUIRE(c1.get_gsl_data().dat[0] == 1.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 2.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test -
    SECTION("Test - operator")
    {
        gsl::complex c3 = c1 - c2;
        REQUIRE(c3.get_gsl_data().dat[0] == -2.0);
        REQUIRE(c3.get_gsl_data().dat[1] == -2.0);
        REQUIRE(c1.get_gsl_data().dat[0] == 1.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 2.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test *
    SECTION("Test * operator")
    {
        gsl::complex c3 = c1 * c2;
        REQUIRE(c3.get_gsl_data().dat[0] == -5.0);
        REQUIRE(c3.get_gsl_data().dat[1] == 10.0);
        REQUIRE(c1.get_gsl_data().dat[0] == 1.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 2.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
    }

    // Test /
    SECTION("Test / operator")
    {
        gsl::complex c3 = c1 / c2;
        REQUIRE_THAT(c3.get_gsl_data().dat[0], Catch::Matchers::WithinAbs(0.44, 1e-12));
        REQUIRE_THAT(c3.get_gsl_data().dat[1], Catch::Matchers::WithinAbs(0.08, 1e-12));
        REQUIRE(c1.get_gsl_data().dat[0] == 1.0);
        REQUIRE(c1.get_gsl_data().dat[1] == 2.0);
        REQUIRE(c2.get_gsl_data().dat[0] == 3.0);
        REQUIRE(c2.get_gsl_data().dat[1] == 4.0);
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

// Use Catch2 to test the various constructors for gsl::cvector
TEST_CASE("gsl::cvector constructors", "[gsl::cvector]")
{
    // Default constructor
    gsl::cvector v1;
    REQUIRE(v1.size() == 0);
    REQUIRE(v1.get_gsl_ptr() == nullptr);

    // Constructor with size
    gsl::cvector v2(10);
    REQUIRE(v2.size() == 10);
    REQUIRE(v2.get_gsl_ptr() != nullptr);

    // Copy constructor
    gsl::cvector v3(v2);
    REQUIRE(v3.size() == 10);
    REQUIRE(v3.get_gsl_ptr() != nullptr);
    REQUIRE(v3.get_gsl_ptr() != v2.get_gsl_ptr());

    // Move constructor
    gsl::cvector v4(std::move(v3));
    REQUIRE(v4.size() == 10);
    REQUIRE(v4.get_gsl_ptr() != nullptr);
    REQUIRE(v3.size() == 0);
    REQUIRE(v3.get_gsl_ptr() == nullptr);
}

// Use Catch2 to test conversion between gsl::vector and gsl::cvector
TEST_CASE("gsl::vector <-> gsl::cvector conversion", "[gsl::vector][gsl::cvector]")
{
    gsl::vector v1(10);
    for (size_t i = 0; i < 10; ++i)
        v1.set(i, (double)i);

    gsl::cvector v2(10);
    for (size_t i = 0; i < 10; ++i)
        v2.set(i, gsl::complex((double)i, (double)(i + 1)));

    SECTION("gsl::vector -> gsl::cvector")
    {
        v2 = v1;
        REQUIRE(v2.get_gsl_ptr() != nullptr);
        REQUIRE(v2.size() == 10);
        for (size_t i = 0; i < 10; ++i)
            REQUIRE(v2.get(i) == gsl::complex((double)i, 0.0));
    }

    SECTION("gsl::cvector -> gsl::vector")
    {
        v1 = v2;
        REQUIRE(v1.get_gsl_ptr() != nullptr);
        REQUIRE(v1.size() == 10);
        for (size_t i = 0; i < 10; ++i)
            REQUIRE(v1.get(i) == (double)i);
    }
}

// Use Catch2 to test the assignment operators for gsl::cvector
TEST_CASE("gsl::cvector assignment operators", "[gsl::cvector]")
{
    gsl::cvector v1(3), v2(3), v3(5);

    // Copy assignment
    v1 = v2;
    REQUIRE(v1.size() == 3);
    REQUIRE(v1.get_gsl_ptr() != nullptr);
    REQUIRE(v1.get_gsl_ptr() != v2.get_gsl_ptr());

    v2 = v3;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get_gsl_ptr() != nullptr);
    REQUIRE(v2.get_gsl_ptr() != v3.get_gsl_ptr());

    // Test self-assignment
    v2 = v2;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get_gsl_ptr() != nullptr);

    // Test Move self-assignment
    v2 = std::move(v2);
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get_gsl_ptr() != nullptr);
}

// Use Catch2 to test the resize and clear method of gsl::cvector
TEST_CASE("gsl::cvector resize", "[gsl::cvector]")
{
    gsl::cvector v(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get_gsl_ptr() != nullptr);

    v.resize(5);
    REQUIRE(v.size() == 5);
    REQUIRE(v.get_gsl_ptr() != nullptr);

    v.resize(0);
    REQUIRE(v.size() == 0);
    REQUIRE(v.get_gsl_ptr() == nullptr);

    v.resize(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get_gsl_ptr() != nullptr);

    v.clear();
    REQUIRE(v.size() == 0);
    REQUIRE(v.get_gsl_ptr() == nullptr);
}

// Use Catch2 to test the assignment and access of array elements
TEST_CASE("gsl::cvector element access", "[gsl::cvector]")
{
    using namespace std::literals::complex_literals;
    auto z = 1.0i;

    gsl::cvector v(10);
    for (size_t i = 0; i < v.size(); i++)
        // v(i) = gsl_complex_rect(i, i); // Desired usage
        v.set(i, gsl_complex_rect(i, i)); // Current usage

    for (size_t i = 0; i < v.size(); i++)
    {
        // REQUIRE(v(i).real() == i); // Desired usage
        // REQUIRE(v(i).imag() == i);
        REQUIRE(v.get(i).real() == i); // Current usage
        REQUIRE(v.get(i).imag() == i); // Current usage
    }
}

// Use Catch2 to test that gsl::cvector::print doesn't crash
TEST_CASE("gsl::cvector print", "[gsl::cvector]")
{
    gsl::cvector v(10);
    v.print();
}

// Use Catch2 to test the various constructors for gsl::cmatrix
TEST_CASE("gsl::cmatrix constructors", "[gsl::cmatrix]")
{
    // Default constructor
    gsl::cmatrix m1;
    REQUIRE(m1.nrows() == 0);
    REQUIRE(m1.ncols() == 0);
    REQUIRE(m1.get_gsl_ptr() == nullptr);

    // Constructor with size
    gsl::cmatrix m2(10, 10);
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 10);
    REQUIRE(m2.get_gsl_ptr() != nullptr);

    // Copy constructor
    gsl::cmatrix m3(m2);
    REQUIRE(m3.nrows() == 10);
    REQUIRE(m3.ncols() == 10);
    REQUIRE(m3.get_gsl_ptr() != nullptr);
    REQUIRE(m3.get_gsl_ptr() != m2.get_gsl_ptr());

    // Move constructor
    gsl::cmatrix m4(std::move(m3));
    REQUIRE(m4.nrows() == 10);
    REQUIRE(m4.ncols() == 10);
    REQUIRE(m4.get_gsl_ptr() != nullptr);
    REQUIRE(m3.nrows() == 0);
    REQUIRE(m3.ncols() == 0);
    REQUIRE(m3.get_gsl_ptr() == nullptr);
}

// Use Catch2 to test conversion between gsl::matrix and gsl::cmatrix
TEST_CASE("gsl::matrix <-> gsl::cmatrix conversion", "[gsl::matrix][gsl::cmatrix]")
{
    gsl::matrix m1(10, 10);
    for (size_t i = 0; i < m1.nrows(); i++)
        for (size_t j = 0; j < m1.ncols(); j++)
            m1.set(i, j, i + j);

    gsl::cmatrix m2(m1);
    for (size_t i = 0; i < m2.nrows(); i++)
        for (size_t j = 0; j < m2.ncols(); j++)
            m2.set(i, j, gsl::complex(i + j, i * j));

    SECTION( "gsl::cmatrix -> gsl::matrix" )
    {
        m1 = m2;
        for (size_t i = 0; i < m1.nrows(); i++)
            for (size_t j = 0; j < m1.ncols(); j++)
                REQUIRE(m1.get(i, j) == i + j);    
    }

    SECTION( "gsl::matrix -> gsl::cmatrix" )
    {
        m2 = m1;
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                REQUIRE(m2.get(i, j) == gsl::complex(i + j, 0.0));
    }
}

// Use Catch2 to test the assignment operators for gsl::cmatrix
TEST_CASE("gsl::cmatrix assignment operators", "[gsl::cmatrix]")
{
    gsl::cmatrix m1(3, 3), m2(3, 3), m3(5, 5);

    // Copy assignment
    m1 = m2;
    REQUIRE(m1.nrows() == 3);
    REQUIRE(m1.ncols() == 3);
    REQUIRE(m1.get_gsl_ptr() != nullptr);
    REQUIRE(m1.get_gsl_ptr() != m2.get_gsl_ptr());

    m2 = m3;
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);
    REQUIRE(m2.get_gsl_ptr() != m3.get_gsl_ptr());

    // Test self-assignment
    m2 = m2;
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);

    // Test Move self-assignment
    m2 = std::move(m2);
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);
}

// Use Catch2 to test the resize and clear method of gsl::cmatrix
TEST_CASE("gsl::cmatrix resize", "[gsl::cmatrix]")
{
    gsl::cmatrix m(10, 10);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 10);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.resize(5, 5);
    REQUIRE(m.nrows() == 5);
    REQUIRE(m.ncols() == 5);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.resize(0, 0);
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
    REQUIRE(m.get_gsl_ptr() == nullptr);

    m.resize(10, 10);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 10);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.clear();
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
    REQUIRE(m.get_gsl_ptr() == nullptr);
}

// Use Catch2 to test element access of gsl::cmatrix
TEST_CASE("gsl::cmatrix element access", "[gsl::cmatrix]")
{
    gsl::cmatrix m(10, 10);

    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m.set(i, j, gsl_complex_rect(i, j));
    // m(i, j) = gsl::complex(i, j);

    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
        {
            REQUIRE(m.get(i, j).real() == i);
            REQUIRE(m.get(i, j).imag() == j);
            // REQUIRE(m(i, j).real() == i);
            // REQUIRE(m(i, j).imag() == j);
        }
}

// Use Catch2 to test that gsl::cmatrix::print doesn't crash
TEST_CASE("gsl::cmatrix print", "[gsl::cmatrix]")
{
    gsl::cmatrix m(10, 10);
    m.print();
}