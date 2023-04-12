#define CATCH_CONFIG_MAIN
#include <gsl_wrapper/core.h>
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

    // Move constructor
    gsl::vector v4(std::move(v3));
    REQUIRE(v4.size() == 10);
    REQUIRE(v4.get_gsl_ptr() != nullptr);
    REQUIRE(v3.size() == 0);
    REQUIRE(v3.get_gsl_ptr() == nullptr);
}

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
TEST_CASE("gsl::vector -> gsl::cvector conversion", "[gsl::vector][gsl::cvector]")
{
    gsl::vector v1(10);
    for (size_t i = 0; i < 10; ++i)
        v1.set(i, (double)i);

    gsl::cvector v2(10);
    for (size_t i = 0; i < 10; ++i)
        v2.set(i, gsl::complex((double)i, (double)(i + 1)));

    v2 = v1;
    REQUIRE(v2.get_gsl_ptr() != nullptr);
    REQUIRE(v2.size() == 10);
    for (size_t i = 0; i < 10; ++i)
        REQUIRE(v2.get(i) == gsl::complex((double)i, 0.0));
}

// Use Catch2 to test conversion between gsl::matrix and gsl::cmatrix
TEST_CASE("gsl::matrix -> gsl::cmatrix conversion", "[gsl::matrix][gsl::cmatrix]")
{
    gsl::matrix m1(10, 10);
    for (size_t i = 0; i < m1.nrows(); i++)
        for (size_t j = 0; j < m1.ncols(); j++)
            m1.set(i, j, i + j);

    gsl::cmatrix m2(m1);
    for (size_t i = 0; i < m2.nrows(); i++)
        for (size_t j = 0; j < m2.ncols(); j++)
            m2.set(i, j, gsl::complex(i + j, i * j));

    SECTION("gsl::matrix -> gsl::cmatrix")
    {
        m2 = m1;
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                REQUIRE(m2.get(i, j) == gsl::complex(i + j, 0.0));
    }
}

TEST_CASE("gsl::vector assignment operators", "[gsl::vector]")
{
    gsl::vector v1(3), v2(3), v3(5);

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
TEST_CASE("gsl::vector element access", "[gsl::vector]")
{
    gsl::vector v(10);
    for (size_t i = 0; i < v.size(); i++)
        v(i) = i;

    for (size_t i = 0; i < v.size(); i++)
        REQUIRE(v(i) == i);
}

// Use Catch2 to test the assignment and access of array elements
TEST_CASE("gsl::cvector element access", "[gsl::cvector]")
{
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

// Use Catch2 to test that print statements don't crash
TEST_CASE("gsl::vector print", "[gsl::vector][gsl::cector][gsl::matrix][gsl::cmatrix]")
{
    gsl::vector v(5);
    gsl::cvector cv(5);
    gsl::matrix m(2, 3);
    gsl::cmatrix cm(2, 3);

    v.print();
    cv.print();
    m.print();
    cm.print();
}

// Use Catch2 to test the various constructors for gsl::matrix
TEST_CASE("gsl::matrix constructors", "[gsl::matrix]")
{
    // Default constructor
    gsl::matrix m0;
    REQUIRE(m0.nrows() == 0);
    REQUIRE(m0.ncols() == 0);
    REQUIRE(m0.get_gsl_ptr() == nullptr);

    // Default constructor
    gsl::matrix m1(0, 0);
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

    // Move constructor
    gsl::matrix m4(11, 6);
    gsl::matrix m5(std::move(m4));
    REQUIRE(m5.nrows() == 11);
    REQUIRE(m5.ncols() == 6);
    REQUIRE(m5.get_gsl_ptr() != nullptr);
    REQUIRE(m4.nrows() == 0);
    REQUIRE(m4.ncols() == 0);
    REQUIRE(m4.get_gsl_ptr() == nullptr);

    // Conversion constructor
    gsl::matrix m6(12, 7);
    gsl::matrix m7(m6.get_gsl_ptr());
    REQUIRE(m7.nrows() == 12);
    REQUIRE(m7.ncols() == 7);
    REQUIRE(m7.get_gsl_ptr() != nullptr);
    REQUIRE(m7.get_gsl_ptr() != m6.get_gsl_ptr());
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

// Use Catch2 to test the assignment and move assignment of gsl::matrix
TEST_CASE("gsl::matrix assignment", "[gsl::matrix]")
{
    gsl::matrix m1(10, 5);
    gsl::matrix m2(11, 6);
    gsl::matrix m3(12, 7);

    // Copy assignment
    m1 = m2;
    REQUIRE(m1.nrows() == 11);
    REQUIRE(m1.ncols() == 6);
    REQUIRE(m1.get_gsl_ptr() != nullptr);
    REQUIRE(m1.get_gsl_ptr() != m2.get_gsl_ptr());

    // Move assignment
    m2 = std::move(m3);
    REQUIRE(m2.nrows() == 12);
    REQUIRE(m2.ncols() == 7);
    REQUIRE(m2.get_gsl_ptr() != nullptr);
    REQUIRE(m3.nrows() == 0);
    REQUIRE(m3.ncols() == 0);
    REQUIRE(m3.get_gsl_ptr() == nullptr);

    // Test self assignment
    m1 = m1;
    REQUIRE(m1.nrows() == 11);
    REQUIRE(m1.ncols() == 6);
    REQUIRE(m1.get_gsl_ptr() != nullptr);

    // Test self move assignment
    m1 = std::move(m1);
    REQUIRE(m1.nrows() == 11);
    REQUIRE(m1.ncols() == 6);
    REQUIRE(m1.get_gsl_ptr() != nullptr);
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