#define CATCH_CONFIG_MAIN
#include <spheroidal/gsl_wrapper.h>
#include <spheroidal/gsl_utils.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gsl/gsl_math.h>
#include <cmath>

// Use Catch2 to test the linspace function
TEST_CASE("linspace", "[linspace]")
{
    // Test linspace with 0 elements
    // Should return an empty vector
    gsl::vector v0 = gsl::linspace(0.0, 1.0, 0);
    REQUIRE(v0.size() == 0);
    REQUIRE(v0.get_gsl_ptr() == nullptr);

    // Test linspace with 1 element
    // Should return a vector with a single element equal to the end point
    gsl::vector v1 = gsl::linspace(0.0, 1.0, 1);
    REQUIRE(v1.size() == 1);
    REQUIRE(v1.get_gsl_ptr() != nullptr);
    REQUIRE(v1(0) == 1.0);

    // Test linspace with 2 elements
    gsl::vector v2 = gsl::linspace(0.0, 1.0, 2);
    REQUIRE(v2.size() == 2);
    REQUIRE(v2.get_gsl_ptr() != nullptr);
    REQUIRE(v2(0) == 0.0);
    REQUIRE(v2(1) == 1.0);

    // Test linspace with 3 elements
    gsl::vector v3 = gsl::linspace(0.0, 1.0, 3);
    REQUIRE(v3.size() == 3);
    REQUIRE(v3.get_gsl_ptr() != nullptr);
    REQUIRE(v3(0) == 0.0);
    REQUIRE(v3(1) == 0.5);
    REQUIRE(v3(2) == 1.0);

    // Test linspace with 4 elements
    gsl::vector v4 = gsl::linspace(0.0, 1.0, 4);
    REQUIRE(v4.size() == 4);
    REQUIRE(v4.get_gsl_ptr() != nullptr);
    REQUIRE(v4(0) == 0.0);
    REQUIRE(v4(1) == 0.3333333333333333);
    REQUIRE(v4(2) == 0.6666666666666666);
    REQUIRE(v4(3) == 1.0);

    // Test linspace with 5 elements
    gsl::vector v5 = gsl::linspace(0.0, 1.0, 5);
    REQUIRE(v5.size() == 5);
    REQUIRE(v5.get_gsl_ptr() != nullptr);
    REQUIRE(v5(0) == 0.0);
    REQUIRE(v5(1) == 0.25);
    REQUIRE(v5(2) == 0.5);
    REQUIRE(v5(3) == 0.75);
    REQUIRE(v5(4) == 1.0);
}

// Use Catch2 to test the leggauss functions
TEST_CASE("leggauss", "[leggauss]")
{
    // Test leggauss with 0 elements
    // Should return an empty vector
    gsl::vector x0, w0;
    gsl::leggauss(0, x0, w0);
    REQUIRE(x0.size() == 0);
    REQUIRE(x0.get_gsl_ptr() == nullptr);
    REQUIRE(w0.size() == 0);
    REQUIRE(w0.get_gsl_ptr() == nullptr);

    // Test leggauss with 1 element
    gsl::vector x1, w1;
    gsl::leggauss(1, x1, w1);
    REQUIRE(x1.size() == 1);
    REQUIRE(x1.get_gsl_ptr() != nullptr);
    REQUIRE(x1(0) == 0.0);
    REQUIRE(w1.size() == 1);
    REQUIRE(w1.get_gsl_ptr() != nullptr);
    REQUIRE(w1(0) == 2.0);

    // Test leggauss with 2 elements
    gsl::vector x2, w2;
    gsl::leggauss(2, x2, w2);
    REQUIRE(x2.size() == 2);
    REQUIRE(x2.get_gsl_ptr() != nullptr);
    REQUIRE(x2(0) == -0.5773502691896257);
    REQUIRE(x2(1) == 0.5773502691896257);
    REQUIRE(w2.size() == 2);
    REQUIRE(w2.get_gsl_ptr() != nullptr);
    REQUIRE(w2(0) == 1.0);
    REQUIRE(w2(1) == 1.0);

    // Test leggauss with 3 elements
    gsl::vector x3, w3;
    gsl::leggauss(3, x3, w3);
    REQUIRE(x3.size() == 3);
    REQUIRE(x3.get_gsl_ptr() != nullptr);
    REQUIRE(x3(0) == -0.7745966692414834);
    REQUIRE(x3(1) == 0.0);
    REQUIRE(x3(2) == 0.7745966692414834);
    REQUIRE(w3.size() == 3);
    REQUIRE(w3.get_gsl_ptr() != nullptr);
    REQUIRE(w3(0) == 0.5555555555555556);
    REQUIRE(w3(1) == 0.8888888888888888);
    REQUIRE(w3(2) == 0.5555555555555556);
}

// Use Catch2 to test the leggauss functions
TEST_CASE("leggauss return", "[leggauss]")
{
    // Test leggauss with 0 elements
    // Should return an empty vector
    gsl::vector x0 = gsl::leggauss(0);
    REQUIRE(x0.size() == 0);
    REQUIRE(x0.get_gsl_ptr() == nullptr);

    // Test leggauss with 1 element
    gsl::vector x1 = gsl::leggauss(1);
    REQUIRE(x1.size() == 1);
    REQUIRE(x1.get_gsl_ptr() != nullptr);
    REQUIRE(x1(0) == 0.0);

    // Test leggauss with 2 elements
    gsl::vector x2 = gsl::leggauss(2);
    REQUIRE(x2.size() == 2);
    REQUIRE(x2.get_gsl_ptr() != nullptr);
    REQUIRE(x2(0) == -0.5773502691896257);
    REQUIRE(x2(1) == 0.5773502691896257);

    // Test leggauss with 3 elements
    gsl::vector x3 = gsl::leggauss(3);
    REQUIRE(x3.size() == 3);
    REQUIRE(x3.get_gsl_ptr() != nullptr);
    REQUIRE(x3(0) == -0.7745966692414834);
    REQUIRE(x3(1) == 0.0);
    REQUIRE(x3(2) == 0.7745966692414834);
}

// Use Catch2 to test the gsl::meshgrid function
TEST_CASE("meshgrid", "[meshgrid]")
{
    // Test meshgrid with 0 elements
    // Should return an empty vector
    gsl::vector x, y;
    gsl::matrix X, Y;
    gsl::meshgrid(x, y, X, Y);
    REQUIRE(X.get_gsl_ptr() == nullptr);
    REQUIRE(Y.get_gsl_ptr() == nullptr);

    // Test meshgrid with size 1 vectors
    x = gsl::vector(1);
    y = gsl::vector(1);

    x(0) = 5.0;
    y(0) = 10.0;

    gsl::meshgrid(x, y, X, Y);
    REQUIRE(X(0, 0) == 5.0);
    REQUIRE(Y(0, 0) == 10.0);

    // Test meshgrid with larger vectors created from linspaces
    x = gsl::linspace(-1.0, 1.0, 4);
    y = gsl::linspace(0.0, 2.0, 5);

    gsl::meshgrid(x, y, X, Y);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            REQUIRE(X(i, j) == x(i));
            REQUIRE(Y(i, j) == y(j));
        }
    }
}

// Use Catch2 to test the gsl::eye function
TEST_CASE("eye", "[eye]")
{
    // Test eye with 0 elements
    // Should return an empty vector
    gsl::matrix I0 = gsl::eye(0);
    REQUIRE(I0.get_gsl_ptr() == nullptr);

    // Test eye with 1 element
    gsl::matrix I1 = gsl::eye(1);
    REQUIRE(I1.get_gsl_ptr() != nullptr);
    REQUIRE(I1(0, 0) == 1.0);

    // Test eye with 3 elements
    gsl::matrix I3 = gsl::eye(3);
    REQUIRE(I3.get_gsl_ptr() != nullptr);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (i == j)
                REQUIRE(I3(i, j) == 1.0);
            else
                REQUIRE(I3(i, j) == 0.0);
}

// Use Catch2 to test the gsl::arrayfun function for vectors
TEST_CASE("arrayfun vector", "[arrayfun]")
{
    // Test applyfun with 0 elements
    // Should return an empty vector
    gsl::vector x0;
    gsl::arrayfun([](double x)
                  { return x * x; },
                  x0);
    REQUIRE(x0.size() == 0);
    REQUIRE(x0.get_gsl_ptr() == nullptr);

    // Test applyfun with 1 element
    gsl::vector t1 = gsl::linspace(0, 5.0, 1);
    gsl::vector x1 = gsl::arrayfun([](double x)
                                   { return x * x; },
                                   std::move(t1));
    REQUIRE(x1.size() == 1);
    REQUIRE(x1.get_gsl_ptr() != nullptr);
    REQUIRE(t1.get_gsl_ptr() == nullptr);
    REQUIRE(x1(0) == 25.0);

    // Test applyfun with 2 elements and cmath function
    gsl::vector x2(2);
    x2(0) = 5.0;
    x2(1) = 10.0;
    x2 = gsl::arrayfun([](double x)
                       { return std::cos(x); },
                       x2);
    REQUIRE(x2.size() == 2);
    REQUIRE(x2.get_gsl_ptr() != nullptr);
    REQUIRE(x2(0) == std::cos(5.0));
    REQUIRE(x2(1) == std::cos(10.0));

    // Test applyfun with 3 elements and gsl_math.h function
    gsl::vector x3(3);
    x3(0) = 2.0;
    x3(1) = 3.0;
    x3(2) = 4.0;
    x3 = gsl::arrayfun([](double x)
                       { return tan(x); },
                       x3);
    REQUIRE(x3.size() == 3);
    REQUIRE(x3.get_gsl_ptr() != nullptr);
    REQUIRE(x3(0) == tan(2.0));
    REQUIRE(x3(1) == tan(3.0));
    REQUIRE(x3(2) == tan(4.0));
}

// Use Catch2 to test the gsl::arrayfun function for matrices
TEST_CASE("arrayfun matrix", "[arrayfun]")
{
    // Test applyfun with 0 elements
    // Should return an empty vector
    gsl::matrix x0;
    gsl::arrayfun([](double x)
                  { return x * x; },
                  x0);
    REQUIRE(x0.nrows() == 0);
    REQUIRE(x0.ncols() == 0);
    REQUIRE(x0.get_gsl_ptr() == nullptr);

    // Test applyfun with 1 element
    gsl::matrix t1 = gsl::eye(1);
    gsl::matrix x1 = gsl::arrayfun([](double x)
                                   { return x + 1; },
                                   std::move(t1));
    REQUIRE(x1.nrows() == 1);
    REQUIRE(x1.ncols() == 1);
    REQUIRE(t1.get_gsl_ptr() == nullptr);
    REQUIRE(x1.get_gsl_ptr() != nullptr);
    REQUIRE(x1(0, 0) == 2.0);

    // Test applyfun with 2 elements and cmath function
    gsl::matrix x2(2, 2);
    x2(0, 0) = 5.0;
    x2(0, 1) = 10.0;
    x2(1, 0) = 15.0;
    x2(1, 1) = 20.0;
    x2 = gsl::arrayfun([](double x)
                       { return std::cos(x); },
                       x2);
    REQUIRE(x2.nrows() == 2);
    REQUIRE(x2.ncols() == 2);
    REQUIRE(x2.get_gsl_ptr() != nullptr);
    REQUIRE(x2(0, 0) == std::cos(5.0));
    REQUIRE(x2(0, 1) == std::cos(10.0));
    REQUIRE(x2(1, 0) == std::cos(15.0));
    REQUIRE(x2(1, 1) == std::cos(20.0));

    // Test applyfun with 3 elements and gsl_math.h function
    gsl::matrix x3(3, 3);
    x3(0, 0) = 2.0;
    x3(0, 1) = 3.0;
    x3(0, 2) = 4.0;
    x3(1, 0) = 5.0;
    x3(1, 1) = 6.0;
    x3(1, 2) = 7.0;
    x3(2, 0) = 8.0;
    x3(2, 1) = 9.0;
    x3(2, 2) = 10.0;
    x3 = gsl::arrayfun([](double x)
                       { return tan(x); },
                       x3);
    REQUIRE(x3.nrows() == 3);
    REQUIRE(x3.ncols() == 3);
    REQUIRE(x3.get_gsl_ptr() != nullptr);
    REQUIRE(x3(0, 0) == tan(2.0));
    REQUIRE(x3(0, 1) == tan(3.0));
    REQUIRE(x3(0, 2) == tan(4.0));
    REQUIRE(x3(1, 0) == tan(5.0));
    REQUIRE(x3(1, 1) == tan(6.0));
    REQUIRE(x3(1, 2) == tan(7.0));
    REQUIRE(x3(2, 0) == tan(8.0));
    REQUIRE(x3(2, 1) == tan(9.0));
    REQUIRE(x3(2, 2) == tan(10.0));
}
