#define CATCH_CONFIG_MAIN
#include <gsl_wrapper/core.h>
#include <catch2/catch_test_macros.hpp>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

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
        //v(i) = gsl_complex_rect(i, i); // Desired usage
        v.set( i, gsl_complex_rect(i, i)); // Current usage

    for (size_t i = 0; i < v.size(); i++)
    {
        //REQUIRE(v(i).real() == i); // Desired usage
        //REQUIRE(v(i).imag() == i);
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
    using namespace gsl::literals;
    gsl::cmatrix m(10, 10);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m.set( i, j, gsl_complex_rect(i, j) );
            //m(i, j) = gsl::complex(i, j);
            
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
        {
            REQUIRE(m.get(i, j).real() == i);
            REQUIRE(m.get(i, j).imag() == j);
            //REQUIRE(m(i, j).real() == i);
            //REQUIRE(m(i, j).imag() == j);
        }
}

// Use Catch2 to test that gsl::cmatrix::print doesn't crash
TEST_CASE("gsl::cmatrix print", "[gsl::cmatrix]")
{
    gsl::cmatrix m(10, 10);
    m.print();
}