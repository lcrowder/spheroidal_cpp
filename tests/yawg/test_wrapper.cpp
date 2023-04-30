#define CATCH_CONFIG_MAIN
#include <yawg/core.h>
#include <gsl/gsl_math.h>
#include <unistd.h>
#include <catch2/catch_test_macros.hpp>

// Use Catch2 to test the various constructors for gsl::vector
TEST_CASE("gsl::vector constructors", "[gsl::vector]")
{
    // Default constructor
    gsl::vector v1;
    REQUIRE(v1.size() == 0);

    // Constructor with size
    gsl::vector v2(10);
    REQUIRE(v2.size() == 10);

    // Copy constructor
    gsl::vector v3(v2);
    REQUIRE(v3.size() == 10);
    REQUIRE(v3.get() != v2.get());

    // Move constructor
    gsl::vector v4(std::move(v3));
    REQUIRE(v3.get() == nullptr);
    REQUIRE(v4.size() == 10);
    REQUIRE(v3.size() == 0);
}

// Use Catch2 to test the various constructors for gsl::complex
TEST_CASE("gsl::complex constructors", "[gsl::complex]")
{
    // Default constructor
    gsl::complex c1;
    REQUIRE(c1.dat[0] == 0.0);
    REQUIRE(c1.dat[1] == 0.0);

    // Constructor with real and imaginary parts
    gsl::complex c2(1.0, 2.0);
    REQUIRE(c2.dat[0] == 1.0);
    REQUIRE(c2.dat[1] == 2.0);

    // Constructor with one parameter
    gsl::complex c3(1.0);
    REQUIRE(c3.dat[0] == 1.0);
    REQUIRE(c3.dat[1] == 0.0);

    // Copy constructor
    gsl::complex c4(c3);
    REQUIRE(c3.dat[0] == 1.0);
    REQUIRE(c3.dat[1] == 0.0);
    REQUIRE(c4.dat[0] == 1.0);
    REQUIRE(c4.dat[1] == 0.0);

    // Test that copy is deep
    c3.set(2.0, 3.0);
    REQUIRE(c3.dat[0] == 2.0);
    REQUIRE(c3.dat[1] == 3.0);
    REQUIRE(c4.dat[0] == 1.0);
    REQUIRE(c4.dat[1] == 0.0);
}

// Use Catch2 to test the various constructors for gsl::cvector
TEST_CASE("gsl::cvector constructors", "[gsl::cvector]")
{
    // Default constructor
    gsl::cvector v1;
    REQUIRE(v1.size() == 0);

    // Constructor with size
    gsl::cvector v2(10);
    REQUIRE(v2.size() == 10);

    // Copy constructor
    gsl::cvector v3(v2);
    REQUIRE(v3.size() == 10);
    REQUIRE(v3.get() != v2.get());

    // Move constructor
    gsl::cvector v4(std::move(v3));
    REQUIRE(v4.size() == 10);
    REQUIRE(v3.get() == nullptr);
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
    REQUIRE(v2.get() != nullptr);
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
    REQUIRE(v1.get() != nullptr);
    REQUIRE(v1.get() != v2.get());

    v2 = v3;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);
    REQUIRE(v2.get() != v3.get());

    // Test self-assignment
    v2 = v2;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);

    // Test Move self-assignment
    v2 = std::move(v2);
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);
}

// Use Catch2 to test the assignment operators for gsl::cvector
TEST_CASE("gsl::cvector assignment operators", "[gsl::cvector]")
{
    gsl::cvector v1(3), v2(3), v3(5);

    // Copy assignment
    v1 = v2;
    REQUIRE(v1.size() == 3);
    REQUIRE(v1.get() != nullptr);
    REQUIRE(v1.get() != v2.get());

    v2 = v3;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);
    REQUIRE(v2.get() != v3.get());

    // Test self-assignment
    v2 = v2;
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);

    // Test Move self-assignment
    v2 = std::move(v2);
    REQUIRE(v2.size() == 5);
    REQUIRE(v2.get() != nullptr);
}

// Use Catch2 to test the resize method of gsl::vector
TEST_CASE("gsl::vector resize", "[gsl::vector]")
{
    gsl::vector v(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get() != nullptr);

    v.resize(5);
    REQUIRE(v.size() == 5);
    REQUIRE(v.get() != nullptr);

    v.resize(0);
    REQUIRE(v.size() == 0);
    REQUIRE(v.get() == nullptr);

    v.resize(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get() != nullptr);

    v.clear();
    REQUIRE(v.size() == 0);
    REQUIRE(v.get() == nullptr);
}

// Use Catch2 to test the resize and clear method of gsl::cvector
TEST_CASE("gsl::cvector resize", "[gsl::cvector]")
{
    gsl::cvector v(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get() != nullptr);

    v.resize(5);
    REQUIRE(v.size() == 5);
    REQUIRE(v.get() != nullptr);

    v.resize(0);
    REQUIRE(v.size() == 0);
    REQUIRE(v.get() == nullptr);

    v.resize(10);
    REQUIRE(v.size() == 10);
    REQUIRE(v.get() != nullptr);

    v.clear();
    REQUIRE(v.size() == 0);
    REQUIRE(v.get() == nullptr);
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
        v(i) = gsl_complex_rect(i, i); // Desired usage

    for (size_t i = 0; i < v.size(); i++)
    {
        REQUIRE(v(i).real() == i); // Desired usage
        REQUIRE(v(i).imag() == i);
    }
}

// Use Catch2 to test that print statements don't crash
TEST_CASE("gsl::vector print", "[gsl::vector][gsl::cector][gsl::matrix][gsl::cmatrix]")
{
    gsl::vector v(5);
    for (size_t i = 0; i < v.size(); i++)
        v(i) = i;

    gsl::cvector cv(5);
    for (size_t i = 0; i < cv.size(); i++)
        cv(1) = gsl_complex_rect(i, i);

    gsl::matrix m(2, 3);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = i + j;

    gsl::cmatrix cm(2, 3);
    for (size_t i = 0; i < cm.nrows(); i++)
        for (size_t j = 0; j < cm.ncols(); j++)
            cm(i, j) = gsl_complex_rect(i + j, i + j);

    v.print();
    cv.print();
    m.print();
    cm.print();
}

// Use Catch2 to test the save_csv and load_csv methods of gsl::matrix
TEST_CASE("gsl::matrix save_csv and load_csv", "[gsl::matrix]")
{
    gsl::matrix m(2, 3);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = i + j;

    // Open a file for writing, saving the matrix
    FILE *fp = fopen("./data/test_matrix.csv", "w");
    if (fp != nullptr)
    {
        REQUIRE(fp != nullptr);
        m.print_csv(fp);
        fclose(fp);

        // Open a file for reading, loading the matrix
        fp = fopen("./data/test_matrix.csv", "r");
        REQUIRE(fp != nullptr);
        gsl::matrix m2(fp);
        fclose(fp);

        // Check that the matrices are the same
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                REQUIRE(m(i, j) == m2(i, j));
    }
}

// Use Catch2 to test the save_csv and load_csv methods of gsl::cmatrix
TEST_CASE("gsl::cmatrix save_csv and load_csv", "[gsl::cmatrix]")
{
    gsl::vector v(4);
    v(0) = -1.0;
    v(1) = -5.0;
    v(2) = 5.0;
    v(3) = 1.0;

    gsl::cmatrix m(4, 4);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = gsl_complex_rect(v(i), v(j));

    // Open a file for writing, saving the matrix
    FILE *fp = fopen("./data/test_cmatrix.csv", "w");
    if (fp != nullptr)
    {
        m.print_csv(fp);
        fclose(fp);

        // Open a file for reading, loading the matrix
        fp = fopen("./data/test_cmatrix.csv", "r");
        REQUIRE(fp != nullptr);
        gsl::cmatrix m2(fp);
        m2.print();
        fclose(fp);

        // Check that the matrices are the same
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
            {
                REQUIRE(m(i, j).real() == m2(i, j).real());
                REQUIRE(m(i, j).imag() == m2(i, j).imag());
            }
    }
}

// Use Catch2 to test the various constructors for gsl::matrix
TEST_CASE("gsl::matrix constructors", "[gsl::matrix]")
{
    // Default constructor
    gsl::matrix m0;
    REQUIRE(m0.nrows() == 0);
    REQUIRE(m0.ncols() == 0);

    // Default constructor
    gsl::matrix m1(0, 0);
    REQUIRE(m1.nrows() == 0);
    REQUIRE(m1.ncols() == 0);

    // Constructor with size
    gsl::matrix m2(10, 5);
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 5);

    // Copy constructor
    gsl::matrix m3(m2);
    REQUIRE(m3.nrows() == 10);
    REQUIRE(m3.ncols() == 5);
    REQUIRE(m3.get() != m2.get());

    // Move constructor
    gsl::matrix m4(11, 6);
    gsl::matrix m5(std::move(m4));
    REQUIRE(m5.nrows() == 11);
    REQUIRE(m5.ncols() == 6);
    REQUIRE(m5.get() != nullptr);
    REQUIRE(m4.nrows() == 0);
    REQUIRE(m4.ncols() == 0);
    REQUIRE(m4.get() == nullptr);
}

// Use Catch2 to test the various constructors for gsl::cmatrix
TEST_CASE("gsl::cmatrix constructors", "[gsl::cmatrix]")
{
    // Default constructor
    gsl::cmatrix m1;
    REQUIRE(m1.nrows() == 0);
    REQUIRE(m1.ncols() == 0);

    // Constructor with size
    gsl::cmatrix m2(10, 10);
    REQUIRE(m2.nrows() == 10);
    REQUIRE(m2.ncols() == 10);

    // Copy constructor
    gsl::cmatrix m3(m2);
    REQUIRE(m3.nrows() == 10);
    REQUIRE(m3.ncols() == 10);
    REQUIRE(m3.get() != m2.get());

    // Move constructor
    gsl::cmatrix m4(std::move(m3));
    REQUIRE(m4.nrows() == 10);
    REQUIRE(m4.ncols() == 10);
    REQUIRE(m4.get() != nullptr);
    REQUIRE(m3.nrows() == 0);
    REQUIRE(m3.ncols() == 0);
    REQUIRE(m3.get() == nullptr);
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
    REQUIRE(m1.get() != nullptr);
    REQUIRE(m1.get() != m2.get());

    // Move assignment
    m2 = std::move(m3);
    REQUIRE(m2.nrows() == 12);
    REQUIRE(m2.ncols() == 7);
    REQUIRE(m2.get() != nullptr);
    REQUIRE(m3.nrows() == 0);
    REQUIRE(m3.ncols() == 0);
    REQUIRE(m3.get() == nullptr);

    // Test self assignment
    m1 = m1;
    REQUIRE(m1.nrows() == 11);
    REQUIRE(m1.ncols() == 6);
    REQUIRE(m1.get() != nullptr);

    // Test self move assignment
    m1 = std::move(m1);
    REQUIRE(m1.nrows() == 11);
    REQUIRE(m1.ncols() == 6);
    REQUIRE(m1.get() != nullptr);
}

// Use Catch2 to test the assignment operators for gsl::cmatrix
TEST_CASE("gsl::cmatrix assignment operators", "[gsl::cmatrix]")
{
    gsl::cmatrix m1(3, 3), m2(3, 3), m3(5, 5);

    // Copy assignment
    m1 = m2;
    REQUIRE(m1.nrows() == 3);
    REQUIRE(m1.ncols() == 3);
    REQUIRE(m1.get() != nullptr);
    REQUIRE(m1.get() != m2.get());

    m2 = m3;
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get() != nullptr);
    REQUIRE(m2.get() != m3.get());

    // Test self-assignment
    m2 = m2;
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get() != nullptr);

    // Test Move self-assignment
    m2 = std::move(m2);
    REQUIRE(m2.nrows() == 5);
    REQUIRE(m2.ncols() == 5);
    REQUIRE(m2.get() != nullptr);
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
            m(i, j) = gsl::complex(i, j);

    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
        {
            REQUIRE(m(i, j).real() == i);
            REQUIRE(m(i, j).imag() == j);
        }

    // Test nested complex assignment
    gsl::cmatrix m1(3, 3);
    gsl::cmatrix m2(3, 3);
    m1(0, 0) = gsl::complex(10, 20);
    m2(0, 0) = gsl::complex(100, 200);

    m1(0, 0) = m2(0, 0) = gsl::complex(1, 2);

    REQUIRE(m2(0, 0) == gsl::complex(1, 2));
    REQUIRE(m1(0, 0) == gsl::complex(1, 2));
}

// Use Catch2 to test the resize and clear methods of gsl::matrix
TEST_CASE("gsl::matrix resize", "[gsl::matrix]")
{
    gsl::matrix m(10, 5);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 5);

    m.resize(5, 10);
    REQUIRE(m.nrows() == 5);
    REQUIRE(m.ncols() == 10);

    m.resize(0, 0);
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);

    m.resize(10, 5);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 5);

    m.clear();
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
}

// Use Catch2 to test the resize and clear method of gsl::cmatrix
TEST_CASE("gsl::cmatrix resize", "[gsl::cmatrix]")
{
    gsl::cmatrix m(10, 10);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 10);

    m.resize(5, 5);
    REQUIRE(m.nrows() == 5);
    REQUIRE(m.ncols() == 5);

    m.resize(0, 0);
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);

    m.resize(10, 10);
    REQUIRE(m.nrows() == 10);
    REQUIRE(m.ncols() == 10);

    m.clear();
    REQUIRE(m.nrows() == 0);
    REQUIRE(m.ncols() == 0);
}

// Use Catch2 to test the reshape method of gsl::matrix
TEST_CASE("gsl::matrix reshape", "[gsl::matrix]")
{
    gsl::matrix m(2, 3);
    for( size_t i = 0; i < m.nrows(); i++ )
        for( size_t j = 0; j < m.ncols(); j++ )
            m(i, j) = i + j;

    gsl::matrix m1 = m.reshape(3, 2);

    for( size_t t = 0; t < m.nrows() * m.ncols(); t++ )
        REQUIRE(m(t / m.ncols(), t % m.ncols()) == m1(t / m1.ncols(), t % m1.ncols()));
}

// Use Catch2 to test the reshape method of gsl::cmatrix
TEST_CASE( "gsl::cmatrix reshape", "[gsl::cmatrix]" )
{
    gsl::cmatrix m(2, 3);
    for( size_t i = 0; i < m.nrows(); i++ )
        for( size_t j = 0; j < m.ncols(); j++ )
            m(i, j) = gsl::complex(i, j);

    gsl::cmatrix m1 = m.reshape(3, 2);

    for( size_t t = 0; t < m.nrows() * m.ncols(); t++ )
        REQUIRE(m(t / m.ncols(), t % m.ncols()) == m1(t / m1.ncols(), t % m1.ncols()));
}