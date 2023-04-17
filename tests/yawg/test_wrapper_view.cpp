#define CATCH_CONFIG_MAIN
#include <yawg/core.h>
#include <catch2/catch_test_macros.hpp>

// Use Catch2 to test the getting and setting of vector views
TEST_CASE("gsl::vector view methods", "[gsl::vector_view]")
{
    gsl::vector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = i;

    gsl::vector v2(5);
    for (int i = 0; i < 5; ++i)
        v2(i) = 10 * i;

    // Test the view method
    SECTION("Test view method")
    {
        // Create view to first five elements
        gsl::vector_view v1_view = v1.subvector(0, 5);

        // Verify that the view and vector point to the same data
        REQUIRE(v1_view.get_gsl_ptr()->data == v1.get_gsl_ptr()->data);

        // Set the view to the values of v2
        v1_view = v2;
        for (int i = 0; i < 5; ++i)
            REQUIRE(v1(i) == 10 * i);

        // Verify that the view and vector still point to the same data
        REQUIRE(v1_view.get_gsl_ptr()->data == v1.get_gsl_ptr()->data);

        // "Dereference" the view to get a new vector
        gsl::vector v1_view_deref = v1_view;
        for (int i = 0; i < 5; ++i)
            REQUIRE(v1_view_deref(i) == 10 * i);

        // After dereferencing, the view and vector should point to different data
        REQUIRE(v1_view_deref.get_gsl_ptr()->data != v1.get_gsl_ptr()->data);
        REQUIRE(v1_view_deref.get_gsl_ptr()->data != v1_view.get_gsl_ptr()->data);

        // Do it again with different constructor
        gsl::vector u1(v1_view);
        for (int i = 0; i < 5; ++i)
            REQUIRE(u1(i) == 10 * i);
        REQUIRE(u1.get_gsl_ptr()->data != v1.get_gsl_ptr()->data);
        REQUIRE(u1.get_gsl_ptr()->data != v1_view.get_gsl_ptr()->data);
    }
}

// Use Catch2 to test the getting and setting of cvector views
TEST_CASE("gsl::cvector view methods", "[gsl::cvector_view]")
{
    gsl::cvector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = gsl::complex(i, 2 * i);

    gsl::cvector v2(5);
    for (int i = 0; i < 5; ++i)
        v2(i) = gsl::complex(10 * i, 20 * i);

    // Test the view method
    SECTION("Test view method")
    {
        // Create view to first five elements
        gsl::cvector_view v1_view = v1.subvector(0, 5);

        // Verify that the view and vector point to the same data
        REQUIRE(v1_view.get_gsl_ptr()->data == v1.get_gsl_ptr()->data);

        // Set the view to the values of v2
        v1_view = v2;
        for (int i = 0; i < 5; ++i)
            REQUIRE(v1(i) == gsl::complex(10 * i, 20 * i));

        // Verify that the view and vector still point to the same data
        REQUIRE(v1_view.get_gsl_ptr()->data == v1.get_gsl_ptr()->data);

        // "Dereference" the view to get a new vector
        gsl::cvector v1_view_deref = v1_view;
        for (int i = 0; i < 5; ++i)
            REQUIRE(v1_view_deref(i) == gsl::complex(10 * i, 20 * i));

        // After dereferencing, the view and vector should point to different data
        REQUIRE(v1_view_deref.get_gsl_ptr()->data != v1.get_gsl_ptr()->data);
        REQUIRE(v1_view_deref.get_gsl_ptr()->data != v1_view.get_gsl_ptr()->data);

        // Do it again with different constructor
        gsl::cvector u1(v1_view);
        for (int i = 0; i < 5; ++i)
            REQUIRE(u1(i) == gsl::complex(10 * i, 20 * i));
        REQUIRE(u1.get_gsl_ptr()->data != v1.get_gsl_ptr()->data);
        REQUIRE(u1.get_gsl_ptr()->data != v1_view.get_gsl_ptr()->data);
    }
}

// Use Catch2 to test the getting and setting of matrix views
TEST_CASE("gsl::matrix view methods", "[gsl::matrix_view]")
{
    gsl::matrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = i + j;

    gsl::matrix m2(5, 5);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            m2(i, j) = 10 * (i + j);

    // Test the view method
    SECTION("Test view method")
    {
        // Create view to first five rows and columns
        gsl::matrix_view m1_view = m1.submatrix(0, 0, 5, 5);

        // Verify that the view and matrix point to the same data
        REQUIRE(m1_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);

        // Set the view to the values of m2
        m1_view = m2;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1(i, j) == 10 * (i + j));

        // Verify that the view and matrix still point to the same data
        REQUIRE(m1_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);

        // "Dereference" the view to get a new matrix
        gsl::matrix m1_view_deref = m1_view;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1_view_deref(i, j) == 10 * (i + j));

        // After dereferencing, the view and matrix should point to different data
        REQUIRE(m1_view_deref.get_gsl_ptr()->data != m1.get_gsl_ptr()->data);
        REQUIRE(m1_view_deref.get_gsl_ptr()->data != m1_view.get_gsl_ptr()->data);

        // Do it again with different constructor
        gsl::matrix u1(m1_view);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(u1(i, j) == 10 * (i + j));
        REQUIRE(u1.get_gsl_ptr()->data != m1.get_gsl_ptr()->data);
        REQUIRE(u1.get_gsl_ptr()->data != m1_view.get_gsl_ptr()->data);
    }
}

// Use Catch2 to test the getting and setting of cmatrix views
TEST_CASE("gsl::cmatrix view methods", "[gsl::cmatrix_view]")
{
    gsl::cmatrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = gsl::complex(i + j, 2 * (i + j));

    gsl::cmatrix m2(5, 5);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            m2(i, j) = gsl::complex(10 * (i + j), 20 * (i + j));

    // Test the view method
    SECTION("Test view method")
    {
        // Create view to first five rows and columns
        gsl::cmatrix_view m1_view = m1.submatrix(0, 0, 5, 5);

        // Verify that the view and matrix point to the same data
        REQUIRE(m1_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);

        // Set the view to the values of m2
        m1_view = m2;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1(i, j) == gsl::complex(10 * (i + j), 20 * (i + j)));

        // Verify that the view and matrix still point to the same data
        REQUIRE(m1_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);

        // "Dereference" the view to get a new matrix
        gsl::cmatrix m1_view_deref = m1_view;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1_view_deref(i, j) == gsl::complex(10 * (i + j), 20 * (i + j)));

        // After dereferencing, the view and matrix should point to different data
        REQUIRE(m1_view_deref.get_gsl_ptr()->data != m1.get_gsl_ptr()->data);
        REQUIRE(m1_view_deref.get_gsl_ptr()->data != m1_view.get_gsl_ptr()->data);

        // Do it again with different constructor
        gsl::cmatrix u1(m1_view);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(u1(i, j) == gsl::complex(10 * (i + j), 20 * (i + j)));
        REQUIRE(u1.get_gsl_ptr()->data != m1.get_gsl_ptr()->data);
        REQUIRE(u1.get_gsl_ptr()->data != m1_view.get_gsl_ptr()->data);
    }
}

// Use Catch2 to test the getting and setting of matrix columns and rows
TEST_CASE("gsl::matrix column and row view methods", "[gsl::matrix][gsl::row_view][gsl::column_view]")
{
    gsl::matrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = i + j;

    gsl::vector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = i;

    // Test the column method
    SECTION("Test column method")
    {
        // Get the first and second column
        gsl::column_view v0_view = m1.column(0);
        gsl::column_view v1_view = m1.column(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 1));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(i, 0) == i + 1);

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 1));
    }

    // Test the row method
    SECTION("Test row method")
    {
        // Get the first and second row
        gsl::row_view v0_view = m1.row(0);
        gsl::row_view v1_view = m1.row(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 10));

        // Set the first row to the second row
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(0, i) == i + 1);

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 10));

        // Verify the stride of each row view is 1
        REQUIRE(v0_view.get_gsl_ptr()->stride == 1);
        REQUIRE(v1_view.get_gsl_ptr()->stride == 1);
    }
}

// Use Catch2 to test the getting and setting of cmatrix columns and rows
TEST_CASE("gsl::cmatrix column and row view methods", "[gsl::cmatrix][gsl::crow_view][gsl::ccolumn_view]")
{
    gsl::cmatrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = gsl::complex(i + j, 2 * (i + j));

    gsl::cvector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = gsl::complex(i, 2 * i);

    // Test the column method
    SECTION("Test column method")
    {
        // Get the first and second column
        gsl::ccolumn_view v0_view = m1.column(0);
        gsl::ccolumn_view v1_view = m1.column(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 2));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(i, 0) == gsl::complex(i + 1, 2 * (i + 1)));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 2));
    }

    // Test the row method
    SECTION("Test row method")
    {
        // Get the first and second row
        gsl::crow_view v0_view = m1.row(0);
        gsl::crow_view v1_view = m1.row(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 20));

        // Set the first row to the second row
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(0, i) == gsl::complex(i + 1, 2 * (i + 1)));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get_gsl_ptr()->data == m1.get_gsl_ptr()->data);
        REQUIRE(v1_view.get_gsl_ptr()->data == (m1.get_gsl_ptr()->data + 20));

        // Verify the stride of each row view is 1
        REQUIRE(v0_view.get_gsl_ptr()->stride == 1);
        REQUIRE(v1_view.get_gsl_ptr()->stride == 1);
    }
}

// Use Catch2 to catch the edge cases of a vector view
TEST_CASE("gsl::vector view edge cases", "[gsl::vector][gsl::vector_view]")
{
    gsl::vector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = i;

    gsl::vector v2(10);
    for (int i = 0; i < 10; ++i)
        v2(i) = i + 1;

    // Test nested assignment of views
    SECTION("Test nested assignment of views")
    {
        gsl::vector_view v0_view = v1.subvector(0, 10);
        gsl::vector_view v1_view = v1.subvector(0, 10);
        v0_view = v1_view = v2;

        for (int i = 0; i < 10; ++i)
        {
            REQUIRE(v1(i) == i + 1);
            REQUIRE(v2(i) == i + 1);
        }
    }

    // Test overlapping vector view assignment
    SECTION("Test overlapping vector view assignment")
    {
        gsl::vector u1(5);
        for (int i = 0; i < 5; ++i)
            u1(i) = 10 * i;

        gsl::vector_view v0_view = v1.subvector(0, 5);
        gsl::vector_view v1_view = v1.subvector(3, 5);
        v0_view = v1_view = u1;

        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(v1(i) == 10 * i);
        }

        for (int i = 5; i < 7; ++i)
        {
            REQUIRE(v1(i) == 10 * (i - 3));
        }
    }
}

// Use Catch2 to catch the edge cases of a cvector view
TEST_CASE("gsl::cvector view edge cases", "[gsl::cvector][gsl::cvector_view]")
{
    gsl::cvector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = gsl::complex(i, 2 * i);

    gsl::cvector v2(10);
    for (int i = 0; i < 10; ++i)
        v2(i) = gsl::complex(i + 1, 2 * (i + 1));

    // Test nested assignment of views
    SECTION("Test nested assignment of views")
    {
        gsl::cvector_view v0_view = v1.subvector(0, 10);
        gsl::cvector_view v1_view = v1.subvector(0, 10);
        v0_view = v1_view = v2;

        for (int i = 0; i < 10; ++i)
        {
            REQUIRE(v1(i) == gsl::complex(i + 1, 2 * (i + 1)));
            REQUIRE(v2(i) == gsl::complex(i + 1, 2 * (i + 1)));
        }
    }

    // Test overlapping vector view assignment
    SECTION("Test overlapping vector view assignment")
    {
        gsl::cvector u1(5);
        for (int i = 0; i < 5; ++i)
            u1(i) = gsl::complex(10 * i, 20 * i);

        gsl::cvector_view v0_view = v1.subvector(0, 5);
        gsl::cvector_view v1_view = v1.subvector(3, 5);
        v0_view = v1_view = u1;

        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(v1(i) == gsl::complex(10 * i, 20 * i));
        }

        for (int i = 5; i < 7; ++i)
        {
            REQUIRE(v1(i) == gsl::complex(10 * (i - 3), 20 * (i - 3)));
        }
    }
}

// Use Catch2 to catch the edge cases of a matrix view
TEST_CASE("gsl::matrix view edge cases", "[gsl::matrix][gsl::matrix_view]")
{
    gsl::matrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = i + j;

    gsl::matrix m2(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m2(i, j) = i + j + 1;

    // Test nested assignment of views
    SECTION("Test nested assignment of views")
    {
        gsl::matrix_view m0_view = m1.submatrix(0, 0, 10, 10);
        gsl::matrix_view m1_view = m1.submatrix(0, 0, 10, 10);
        m0_view = m1_view = m2;

        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 10; ++j)
            {
                REQUIRE(m1(i, j) == i + j + 1);
                REQUIRE(m2(i, j) == i + j + 1);
            }
    }

    // Test overlapping matrix view assignment
    SECTION("Test overlapping matrix view assignment")
    {
        gsl::matrix u1(5, 5);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                u1(i, j) = 10 * i + j;

        gsl::matrix_view m0_view = m1.submatrix(0, 0, 5, 5);
        gsl::matrix_view m1_view = m1.submatrix(3, 3, 5, 5);
        m0_view = m1_view = u1;

        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1(i, j) == 10 * i + j);

        for (int i = 5; i < 7; ++i)
            for (int j = 5; j < 7; ++j)
                REQUIRE(m1(i, j) == 10 * (i - 3) + (j - 3));
    }
}

// Use Catch2 to catch the edge cases of a catrix view
TEST_CASE("gsl::cmatrix view edge cases", "[gsl::cmatrix][gsl::cmatrix_view]")
{
    gsl::cmatrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = gsl::complex(i + j, 2 * (i + j));

    gsl::cmatrix m2(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m2(i, j) = gsl::complex(i + j + 1, 2 * (i + j + 1));

    // Test nested assignment of views
    SECTION("Test nested assignment of views")
    {
        gsl::cmatrix_view m0_view = m1.submatrix(0, 0, 10, 10);
        gsl::cmatrix_view m1_view = m1.submatrix(0, 0, 10, 10);
        m0_view = m1_view = m2;

        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 10; ++j)
            {
                REQUIRE(m1(i, j) == gsl::complex(i + j + 1, 2 * (i + j + 1)));
                REQUIRE(m2(i, j) == gsl::complex(i + j + 1, 2 * (i + j + 1)));
            }
    }

    // Test overlapping matrix view assignment
    SECTION("Test overlapping matrix view assignment")
    {
        gsl::cmatrix u1(5, 5);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                u1(i, j) = gsl::complex(10 * i + j, 20 * i + j);

        gsl::cmatrix_view m0_view = m1.submatrix(0, 0, 5, 5);
        gsl::cmatrix_view m1_view = m1.submatrix(3, 3, 5, 5);
        m0_view = m1_view = u1;

        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                REQUIRE(m1(i, j) == gsl::complex(10 * i + j, 20 * i + j));

        for (int i = 5; i < 7; ++i)
            for (int j = 5; j < 7; ++j)
                REQUIRE(m1(i, j) == gsl::complex(10 * (i - 3) + (j - 3), 20 * (i - 3) + (j - 3)));
    }
}
