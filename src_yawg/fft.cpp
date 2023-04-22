#include <yawg/core.h>
#include <yawg/fft.h>
#include <gsl/gsl_fft_complex.h>

/*! \brief Compute in-place fft of a gsl::(c)vector
 * Computes the in-place complex Fast Fourier Transform of the data in x.
 * \param x The cvector to be transformed.
 *
 * \note More efficient implementations of this function (and the others in this file)
 * would use gsl_fft_real_* functions for real arguments, or cache the wavetable
 * and workspace for repeated calls in different parts of execution.
 *
 * \return The transformed vector, reclaiming the memory of the original vector.
 */
gsl::cvector gsl::fft(gsl::cvector &&x)
{
    const gsl_vector_complex *xp = x.get();

    // evil floating point bit level hacking
    if ((xp->size & (xp->size - 1)) == 0)
        gsl_fft_complex_radix2_forward(xp->data, xp->stride, xp->size);
    else
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(xp->size);
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(xp->size);
        gsl_fft_complex_forward(xp->data, xp->stride, xp->size, wavetable, workspace);
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

/*! \brief Compute fft of a gsl::(c)vector
 * Computes the complex Fast Fourier Transform of the data in x.
 * \param x The cvector to be transformed.
 *
 * \return The transformed vector
 */
gsl::cvector gsl::fft(const gsl::cvector &x)
{
    gsl::cvector y(x);
    return gsl::fft(std::move(y));
}

/*! \brief Compute in-place inverse fft of a gsl::(c)vector
 * Computes the in-place complex Inverse Fast Fourier Transform of the data in x.
 * \param x The cvector to be transformed.
 *
 * \return The transformed vector, reclaiming the memory of the original vector.
 */
gsl::cvector gsl::ifft(gsl::cvector &&x)
{
    const gsl_vector_complex *xp = x.get();

    // evil floating point bit level hacking
    if ((xp->size & (xp->size - 1)) == 0)
        gsl_fft_complex_radix2_inverse(xp->data, xp->stride, xp->size);
    else
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(xp->size);
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(xp->size);
        gsl_fft_complex_inverse(xp->data, xp->stride, xp->size, wavetable, workspace);
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

/*! \brief Compute inverse fft of a gsl::(c)vector
 * Computes the complex inverse Fast Fourier Transform of the data in x.
 * \param x The cvector to be transformed.
 *
 * \return The transformed vector
 */
gsl::cvector gsl::ifft(const gsl::cvector &x)
{
    gsl::cvector y(x);
    return gsl::ifft(std::move(y));
}

/*! \brief Compute in-place 1D fft of each column/row of a gsl::(c)matrix
 * Computes the 1D complex Fast Fourier Transform of each row (dim=1) or column (dim=2)
 * of the data in x.
 * \param x The cmatrix to be transformed.
 * \param dim The dimension to transform. 1 for columns, 2 for rows.
 *
 * \return The transformed matrix, reclaiming memory of the original matrix
 */
gsl::cmatrix gsl::fft(gsl::cmatrix &&x, int dim)
{
    const gsl_matrix_complex *xp = x.get();

    // If operating on rows and length of row is power of 2, use radix2 method
    if (dim == 2 && ((x.ncols() & (x.ncols() - 1)) == 0))
    {
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get(), i);
            gsl_fft_complex_radix2_forward((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size);
        }
    }
    // If operating on columns and length of column is power of 2, use radix2 method
    else if (dim == 1 && ((x.nrows() & (x.nrows() - 1)) == 0))
    {
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get(), j);
            gsl_fft_complex_radix2_forward((&col.vector)->data, (&col.vector)->stride, (&col.vector)->size);
        }
    }
    // If operating on columns and length of column is not power of 2, use general method
    else if (dim == 1)
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(x.nrows());
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(x.nrows());
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get(), j);
            gsl_fft_complex_forward((&col.vector)->data, (&col.vector)->stride, (&col.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }
    // If operating on rows (assumed to be the case) and nrows is not power of 2, use general method
    else
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(x.ncols());
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(x.ncols());
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get(), i);
            gsl_fft_complex_forward((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

/*! \brief Compute 1D fft of each column/row of a gsl::(c)matrix
 * Computes the 1D complex Fast Fourier Transform of each row (dim=1) or column (dim=2)
 * of the data in x.
 * \param x The cmatrix to be transformed.
 * \param dim The dimension to transform. 1 for columns, 2 for rows.
 *
 * \return The transformed matrix
 */
gsl::cmatrix gsl::fft(const gsl::cmatrix &x, int dim)
{
    gsl::cmatrix y(x);
    return gsl::fft(std::move(y), dim);
}

/*! \brief Compute in-place 1D inverse fft of each column/row of a gsl::(c)matrix
 * Computes the 1D complex inverse Fast Fourier Transform of each row (dim=1) or column (dim=2)
 * of the data in x.
 * \param x The cmatrix to be transformed.
 * \param dim The dimension to transform. 1 for columns, 2 for rows.
 *
 * \return The transformed matrix, reclaiming memory of the original matrix
 */
gsl::cmatrix gsl::ifft(gsl::cmatrix &&x, int dim)
{
    gsl_matrix_complex *xp = x.get();

    // If operating on rows and length of row is power of 2, use radix2 method
    if (dim == 2 && ((x.ncols() & (x.ncols() - 1)) == 0))
    {
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get(), i);
            gsl_fft_complex_radix2_inverse((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size);
        }
    }
    // If operating on columns and length of column is power of 2, use radix2 method
    else if (dim == 1 && ((x.nrows() & (x.nrows() - 1)) == 0))
    {
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get(), j);
            gsl_fft_complex_radix2_inverse((&col.vector)->data, (&col.vector)->stride, (&col.vector)->size);
        }
    }
    // If operating on columns and length of column is not power of 2, use general method
    else if (dim == 1)
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(x.nrows());
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(x.nrows());
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get(), j);
            gsl_fft_complex_inverse((&col.vector)->data, (&col.vector)->stride, (&col.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }
    // If operating on rows (assumed to be the case) and nrows is not power of 2, use general method
    else
    {
        gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(x.ncols());
        gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(x.ncols());
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get(), i);
            gsl_fft_complex_inverse((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

/*! \brief Compute 1D inverse fft of each column/row of a gsl::(c)matrix
 * Computes the 1D complex inverse Fast Fourier Transform of each row (dim=1) or column (dim=2)
 * of the data in x.
 * \param x The cmatrix to be transformed.
 * \param dim The dimension to transform. 1 for columns, 2 for rows.
 *
 * \return The transformed matrix
 */
gsl::cmatrix gsl::ifft(const gsl::cmatrix &x, int dim)
{
    gsl::cmatrix y(x);
    return gsl::ifft(std::move(y), dim);
}