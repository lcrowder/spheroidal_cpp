#define CATCH_CONFIG_MAIN
#include <spheroidal/gsl_wrapper.h>
#include <catch2/catch_test_macros.hpp>

// Use Catch2 to test the various constructors for gsl::vector
TEST_CASE("gsl::vector constructors", "[gsl::vector]")
{
    // Default constructor
    gsl::vector v1;
    REQUIRE(v1.size() == 0);
    REQUIRE(v1.get_gsl_ptr() == nullptr);

    // Constructor with size
    gsl::vector v2(10);
    REQUIRE(v2.size() == 10);
    REQUIRE(v2.get_gsl_ptr() != nullptr);

    // Copy constructor
    gsl::vector v3(v2);
    REQUIRE(v3.size() == 10);
    REQUIRE(v3.get_gsl_ptr() != nullptr);
    REQUIRE(v3.get_gsl_ptr() != v2.get_gsl_ptr());

    // Copy assignment
    v2 = v3;
    REQUIRE(v2.size() == 10);
    REQUIRE(v2.get_gsl_ptr() != nullptr);
    REQUIRE(v2.get_gsl_ptr() != v3.get_gsl_ptr());

    // Test self-assignment
    v2 = v2;
    REQUIRE(v2.size() == 10);
    REQUIRE(v2.get_gsl_ptr() != nullptr);
}

// Use Catch2 to test the resize method of gsl::vector
TEST_CASE("gsl::vector resize", "[gsl::vector]")
{
    gsl::vector v(10);
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
TEST_CASE("gsl::vector element access", "[gsl::vector]")
{
    gsl::vector v(10);
    for (size_t i = 0; i < v.size(); i++)
        v(i) = i;

    for (size_t i = 0; i < v.size(); i++)
        REQUIRE(v(i) == i);
}

// Use Catch2 to test that gsl::vector::print doesn't crash
TEST_CASE("gsl::vector print", "[gsl::vector]")
{
    gsl::vector v(10);
    v.print();
}

// Use Catch2 to test the various constructors for gsl::matrix
TEST_CASE("gsl::matrix constructors", "[gsl::matrix]")
{
    // Default constructor
    gsl::matrix m1;
    REQUIRE(m1.nrows() == 0);
    REQUIRE(m1.ncols() == 0);
    REQUIRE(m1.get_gsl_ptr() == nullptr);

    // Constructor with size
    gsl::matrix m2(10, 5);
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);

    // Copy constructor
    gsl::matrix m3(m2);
    REQUIRE(m3.nrows() == 10);
    REQUIRE(m3.ncols() == 5);
    REQUIRE(m3.get_gsl_ptr() != nullptr);
    REQUIRE(m3.get_gsl_ptr() != m2.get_gsl_ptr());

    // Copy assignment
    m2 = m3;
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);
    REQUIRE(m2.get_gsl_ptr() != m3.get_gsl_ptr());

    // Test self-assignment
    m2 = m2;
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get_gsl_ptr() != nullptr);
}

// Use Catch2 to test the element accessors of gsl::matrix
TEST_CASE("gsl::matrix element access", "[gsl::matrix]")
{
    gsl::matrix m(10, 5);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = i + j;

    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            REQUIRE(m(i, j) == i + j);
}

// Use Catch2 to test that gsl::matrix::print doesn't crash
TEST_CASE("gsl::matrix print", "[gsl::matrix]")
{
    gsl::matrix m(3, 4);
    m.print();
}

// Use Catch2 to test the resize and clear methods of gsl::matrix
TEST_CASE("gsl::matrix resize", "[gsl::matrix]")
{
    gsl::matrix m(10, 5);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 5);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.resize(5, 10);
    REQUIRE(m.nrows() == 5);
    REQUIRE(m.ncols() == 10);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.resize(0, 0);
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
    REQUIRE(m.get_gsl_ptr() == nullptr);

    m.resize(10, 5);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 5);
    REQUIRE(m.get_gsl_ptr() != nullptr);

    m.clear();
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
    REQUIRE(m.get_gsl_ptr() == nullptr);
}