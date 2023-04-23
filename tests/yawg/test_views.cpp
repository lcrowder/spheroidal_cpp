#define CATCH_CONFIG_MAIN
#include <yawg/core.h>
#include <gsl/gsl_math.h>
#include <catch2/catch_test_macros.hpp>

// Use Catch2 to test various operations on a vector_view,
//  such as assignment, addition, norm.
TEST_CASE("gsl::vector_view operations", "[gsl::vector_view]")
{
    gsl::vector v(10);
    for (size_t i = 0; i < v.size(); i++)
        v(i) = i;

    gsl::vector_view vv = v.subvector(2, 5);

    // SECTION to test element access of the subvector.
    // Changes should be reflected in the original vector
    SECTION("subvector element access")
    {
        // Add 2 to each element of the subvector
        for (size_t i = 0; i < vv.size(); i++)
            vv(i) += 2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == i + 2);
            else
                REQUIRE(v(i) == i);
    }

    // SECTION to test the size, norm operators on the view vv
    SECTION("subvector size and norm")
    {
        REQUIRE(vv.size() == 5);
        REQUIRE(vv.norm() == sqrt(2 * 2 + 3 * 3 + 4 * 4 + 5 * 5 + 6 * 6));
    }

    // SECTION to test the assignment operator of the subvector
    SECTION("subvector assignment from vector")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        vv = v2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2));
            else
                REQUIRE(v(i) == i);
    }

    SECTION("Subvector assignment from vector view")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::vector_view vv2 = v2.view();

        vv = vv2;

        // Check that the original vector and view have been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
            {
                REQUIRE(v(i) == 100 * (i - 2));
                REQUIRE(vv(i - 2) == 100 * (i - 2));
            }
            else
                REQUIRE(v(i) == i);
    }

    SECTION("vector assignment from vector_view")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::vector_view vv2 = v2.view();

        v = vv2;

        for (size_t i = 0; i < v.size(); i++)
            v(i) = 100 * i;
    }

    // Test vector_view += vector
    SECTION("Subvector += vector")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        vv += v2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2) + i);
            else
                REQUIRE(v(i) == i);
    }

    // Test vector_view += vector_view
    SECTION("Subvector += vector_view")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::vector_view vv2 = v2.view();

        vv += vv2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2) + i);
            else
                REQUIRE(v(i) == i);
    }

    // Test creation of views from views
    SECTION("Views of views")
    {
        gsl::vector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::vector_view vv2 = v2.view();

        gsl::vector_view vv3 = vv2.view();

        vv = vv3;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2));
            else
                REQUIRE(v(i) == i);
    }

    // Test nested assignment of views
    SECTION("Nested assignment of views")
    {
        gsl::vector u(5);
        for (int i = 0; i < 5; ++i)
            u(i) = 10 * i;

        gsl::vector_view vv2 = v.subvector(0, 5);
        gsl::vector_view vv3 = v.subvector(3, 5);
        vv2 = vv3 = u;

        for (int i = 0; i < 5; ++i)
            REQUIRE(v(i) == 10 * i);

        for (int i = 5; i < 7; ++i)
            REQUIRE(v(i) == 10 * (i - 3));
    }
}

// Use Catch2 to test various operations on a cvector
//  such as assignment, addition, norm.
TEST_CASE("gsl::cvector operations", "[gsl::cvector]")
{
    gsl::cvector v(10);
    for (size_t i = 0; i < v.size(); i++)
        v(i) = i;

    gsl::cvector_view vv = v.subvector(2, 5);

    // SECTION to test element access of the subvector.
    // Changes should be reflected in the original vector
    SECTION("subvector element access")
    {
        // Add 2 to each element of the subvector
        for (size_t i = 0; i < vv.size(); i++)
            vv(i) += 2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == i + 2);
            else
                REQUIRE(v(i) == i);
    }

    // SECTION to test the size, norm operators on the view vv
    SECTION("subvector size and norm")
    {
        REQUIRE(vv.size() == 5);
        REQUIRE(vv.norm() == sqrt(2 * 2 + 3 * 3 + 4 * 4 + 5 * 5 + 6 * 6));
    }

    // SECTION to test the assignment operator of the subvector
    SECTION("subvector assignment from vector")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        vv = v2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2));
            else
                REQUIRE(v(i) == i);
    }

    SECTION("Subvector assignment from vector view")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::cvector_view vv2 = v2.view();

        vv = vv2;

        // Check that the original vector and view have been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
            {
                REQUIRE(v(i) == 100 * (i - 2));
                REQUIRE(vv(i - 2) == 100 * (i - 2));
            }
            else
                REQUIRE(v(i) == i);
    }

    SECTION("vector assignment from vector_view")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::cvector_view vv2 = v2.view();

        v = vv2;

        for (size_t i = 0; i < v.size(); i++)
            v(i) = 100 * i;
    }

    // Test vector_view += vector
    SECTION("Subvector += vector")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        vv += v2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2) + i);
            else
                REQUIRE(v(i) == i);
    }

    // Test vector_view += vector_view
    SECTION("Subvector += vector_view")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::cvector_view vv2 = v2.view();

        vv += vv2;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2) + i);
            else
                REQUIRE(v(i) == i);
    }

    // Test creation of views from views
    SECTION("Views of views")
    {
        gsl::cvector v2(5);
        for (size_t i = 0; i < v2.size(); i++)
            v2(i) = 100 * i;

        gsl::cvector_view vv2 = v2.view();

        gsl::cvector_view vv3 = vv2.view();

        vv = vv3;

        // Check that the original vector has been modified
        for (size_t i = 0; i < v.size(); i++)
            if (i >= 2 && i < 7)
                REQUIRE(v(i) == 100 * (i - 2));
            else
                REQUIRE(v(i) == i);
    }

    // Test nested assignment of views
    SECTION("Nested assignment of views")
    {
        gsl::cvector u(5);
        for (size_t i = 0; i < u.size(); i++)
            u(i) = 10 * i;

        gsl::cvector_view vv2 = v.subvector(0, 5);
        gsl::cvector_view vv3 = v.subvector(3, 5);

        vv2 = vv3 = u;

        for (int i = 0; i < 5; ++i)
            REQUIRE(v(i) == 10 * i);

        for (int i = 5; i < 7; ++i)
            REQUIRE(v(i) == 10 * (i - 3));
    }
}

// Use Catch2 to test various operations on a matrix_view,
//  such as assignment, addition, norm.
TEST_CASE("matrix_view", "[matrix_view]")
{
    gsl::matrix m(5, 5);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = i * m.ncols() + j;

    gsl::matrix_view mv = m.submatrix(2, 2, 3, 3);

    SECTION("submatrix assignment from matrix")
    {
        gsl::matrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        mv = m2;

        // Check that the original matrix has been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                    REQUIRE(m(i, j) == 100 * (i - 2) + (j - 2));
                else
                    REQUIRE(m(i, j) == i * m.ncols() + j);
    }

    SECTION("submatrix assignment from matrix view")
    {
        gsl::matrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        gsl::matrix_view mv2 = m2.view();

        mv = mv2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                {
                    REQUIRE(m(i, j) == 100 * (i - 2) + (j - 2));
                    REQUIRE(mv(i - 2, j - 2) == 100 * (i - 2) + (j - 2));
                }
                else
                    REQUIRE(m(i, j) == i * m.ncols() + j);
    }

    SECTION("matrix assignment from matrix_view")
    {
        gsl::matrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        gsl::matrix_view mv2 = m2.view();

        m = mv2;

        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                m(i, j) = 100 * i + j;
    }

    SECTION("submatrix += matrix")
    {
        gsl::matrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        mv += m2;

        // Check that the original matrix has been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                    REQUIRE(m(i, j) == 100 * (i - 2) + (j - 2) + i * m.ncols() + j);
                else
                    REQUIRE(m(i, j) == i * m.ncols() + j);
    }

    SECTION("submatrix += matrix_view")
    {
        gsl::matrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        gsl::matrix_view mv2 = m2.view();

        mv += mv2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                {
                    REQUIRE(m(i, j) == 100 * (i - 2) + (j - 2) + i * m.ncols() + j);
                    REQUIRE(mv(i - 2, j - 2) == 100 * (i - 2) + (j - 2) + i * m.ncols() + j);
                }
                else
                    REQUIRE(m(i, j) == i * m.ncols() + j);
    }

    SECTION("Views of views")
    {
        gsl::matrix_view mv2 = m.submatrix(0, 0, 3, 3);
        gsl::matrix_view mv3 = mv2.submatrix(0, 0, 2, 2);

        for (size_t i = 0; i < mv3.nrows(); i++)
            for (size_t j = 0; j < mv3.ncols(); j++)
                mv3(i, j) = 100 * i + j;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 0 && i < 2 && j >= 0 && j < 2)
                    REQUIRE(m(i, j) == 100 * i + j);
                else
                    REQUIRE(m(i, j) == i * m.ncols() + j);
    }

    SECTION("Nested assignment of views")
    {
        gsl::matrix m2(2, 2);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = 100 * i + j;

        gsl::matrix_view mv2 = m.submatrix(0, 0, 2, 2);
        gsl::matrix_view mv3 = m.submatrix(2, 2, 2, 2);

        mv2 = mv3 = m2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < 2; i++)
            for (size_t j = 0; j < 2; j++)
                REQUIRE(m(i, j) == 100 * i + j);

        for (size_t i = 0; i < mv3.ncols(); i++)
            for (size_t j = 0; j < mv3.ncols(); j++)
                REQUIRE(mv3(i, j) == 100 * i + j);
    }
}

TEST_CASE("cmatrix_view", "[cmatrix_view]")
{
    gsl::cmatrix m(5, 5);
    for (size_t i = 0; i < m.nrows(); i++)
        for (size_t j = 0; j < m.ncols(); j++)
            m(i, j) = gsl::complex(i * m.ncols() + j, 0);

    gsl::cmatrix_view mv = m.submatrix(2, 2, 3, 3);

    SECTION("submatrix assignment from complex matrix")
    {
        gsl::cmatrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);

        mv = m2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                    REQUIRE(m(i, j) == gsl::complex(100 * (i - 2) + (j - 2), 0));
                else
                    REQUIRE(m(i, j) == gsl::complex(i * m.ncols() + j, 0));
    }

    SECTION("submatrix assignment from complex matrix view")
    {
        gsl::cmatrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);

        gsl::cmatrix_view mv2 = m2.view();

        mv = mv2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                {
                    REQUIRE(m(i, j) == gsl::complex(100 * (i - 2) + (j - 2), 0));
                    REQUIRE(mv(i - 2, j - 2) == gsl::complex(100 * (i - 2) + (j - 2), 0));
                }
                else
                    REQUIRE(m(i, j) == gsl::complex(i * m.ncols() + j, 0));
    }

    SECTION("matrix assignment from complex matrix view")
    {
        gsl::cmatrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);

        gsl::cmatrix_view mv2 = m2.view();

        m = mv2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                m(i, j) = 100 * i + j;
    }

    SECTION("submatrix += matrix")
    {
        gsl::cmatrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);


        mv += m2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                    REQUIRE(m(i, j) == gsl::complex(100 * (i-2) + (j-2), 0) + gsl::complex(i * m.ncols() + j, 0));
                else
                    REQUIRE(m(i, j) == gsl::complex(i * m.ncols() + j, 0));
    }

    SECTION("submatrix += matrix view")
    {
        gsl::cmatrix m2(3, 3);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);

        gsl::cmatrix_view mv2 = m2.view();

        mv += mv2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 2 && i < 5 && j >= 2 && j < 5)
                {
                    REQUIRE(m(i, j) == gsl::complex(100 * (i-2) + j-2 , 0) + gsl::complex(i * m.ncols() + j, 0));
                    REQUIRE(mv(i-2, j-2) == gsl::complex(100 * (i-2) + j-2, 0) + gsl::complex(i * m.ncols() + j, 0));    
                }
                else
                    REQUIRE(m(i, j) == gsl::complex(i * m.ncols() + j, 0));
    }

    SECTION("Views of views")
    {
        gsl::cmatrix_view mv2 = m.submatrix(0, 0, 3, 3);
        gsl::cmatrix_view mv3 = mv2.submatrix(0, 0, 2, 2);

        for (size_t i = 0; i < mv3.nrows(); i++)
            for (size_t j = 0; j < mv3.ncols(); j++)
                mv3(i, j) = gsl::complex(100 * i + j, 0);

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < m.nrows(); i++)
            for (size_t j = 0; j < m.ncols(); j++)
                if (i >= 0 && i < 2 && j >= 0 && j < 2)
                    REQUIRE(m(i, j) == gsl::complex(100 * i + j, 0));
                else
                    REQUIRE(m(i, j) == gsl::complex(i * m.ncols() + j, 0));
    }

    SECTION("Nested assignment of views")
    {
        gsl::cmatrix m2(2, 2);
        for (size_t i = 0; i < m2.nrows(); i++)
            for (size_t j = 0; j < m2.ncols(); j++)
                m2(i, j) = gsl::complex(100 * i + j, 0);

        gsl::cmatrix_view mv2 = m.submatrix(0, 0, 2, 2);
        gsl::cmatrix_view mv3 = m.submatrix(2, 2, 2, 2);

        mv2 = mv3 = m2;

        // Check that the original matrix and view have been modified
        for (size_t i = 0; i < 2; i++)
            for (size_t j = 0; j < 2; j++)
                REQUIRE(m(i, j) == gsl::complex(100 * i + j, 0));

        for (size_t i = 0; i < mv3.ncols(); i++)
            for (size_t j = 0; j < mv3.ncols(); j++)
                REQUIRE(mv3(i, j) == gsl::complex(100 * i + j, 0));
    }
}

TEST_CASE("gsl::matrix column and row view methods", "[gsl::matrix][gsl::row_view][gsl::column_view]")
{
    gsl::matrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = i * 10 + j;

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
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 1));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(i, 0) == m1(i, 1));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 1));
    }

    // Test the row method
    SECTION("Test row method")
    {
        // Get the first and second column
        gsl::row_view v0_view = m1.row(0);
        gsl::row_view v1_view = m1.row(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 10));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(0,i) == m1(1, i));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 10));

        REQUIRE(v0_view.get()->stride == 1);
        REQUIRE(v0_view.get()->stride == 1);

        // Reshape v0_view into a matrix and check the rows and columns
        gsl::matrix_view m2 = v0_view.reshape(2, 5);
        REQUIRE(m2.nrows() == 2);
        REQUIRE(m2.ncols() == 5);
    }
}

TEST_CASE("gsl::cmatrix column and row view methods", "[gsl::cmatrix][gsl::crow_view][gsl::ccolumn_view]")
{
    gsl::cmatrix m1(10, 10);
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            m1(i, j) = i * 10 + j;

    gsl::cvector v1(10);
    for (int i = 0; i < 10; ++i)
        v1(i) = i;

    // Test the column method
    SECTION("Test column method")
    {
        // Get the first and second column
        gsl::ccolumn_view v0_view = m1.column(0);
        gsl::ccolumn_view v1_view = m1.column(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 2));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(i, 0) == m1(i, 1));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 2));
    }

    // Test the row method
    SECTION("Test row method")
    {
        // Get the first and second column
        gsl::crow_view v0_view = m1.row(0);
        gsl::crow_view v1_view = m1.row(1);

        // Verify that the views point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 20));

        // Set the first column to the second colum
        v0_view = v1_view;
        for (int i = 0; i < 10; ++i)
            REQUIRE(m1(0,i) == m1(1, i));

        // Verify that the views still point to the matrix data
        REQUIRE(v0_view.get()->data == m1.get()->data);
        REQUIRE(v1_view.get()->data == (m1.get()->data + 20));

        REQUIRE(v0_view.get()->stride == 1);
        REQUIRE(v0_view.get()->stride == 1);

        // Reshape v0_view into a matrix and check the rows and columns
        gsl::cmatrix_view m2 = v0_view.reshape(2, 5);
        REQUIRE(m2.nrows() == 2);
        REQUIRE(m2.ncols() == 5);
    }
}