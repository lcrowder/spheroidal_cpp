#include <gsl_wrapper/core.h>
#include <gsl_wrapper/fft.h>
#include <gsl/gsl_fft_complex.h>

/*! \brief Compute fft of a gsl::(c)vector */
//! \note If x has only real data, using gsl_fft_real_* methods would be faster.
//!   A more efficient implementation would cache the wavetable and workspace
gsl::cvector gsl::fft(gsl::cvector &&x)
{
    gsl_vector_complex *xp = x.get_gsl_ptr();

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

gsl::cvector gsl::fft(const gsl::cvector &x)
{
    gsl::cvector y(x);
    return gsl::fft(std::move(y));
}

/*! \brief Compute ifft of a gsl::(c)vector */
gsl::cvector gsl::ifft(gsl::cvector &&x)
{
    gsl_vector_complex *xp = x.get_gsl_ptr();

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

gsl::cvector gsl::ifft(const gsl::cvector &x)
{
    gsl::cvector y(x);
    return gsl::ifft(std::move(y));
}

/*! \brief Compute fft of each row (dim=1) / column (dim=2) of a gsl::(c)matrix */
gsl::cmatrix gsl::fft(gsl::cmatrix &&x, int dim)
{
    gsl_matrix_complex *xp = x.get_gsl_ptr();

    // If operating on rows and length of row is power of 2, use radix2 method
    if (dim == 2 && ((x.ncols() & (x.ncols() - 1)) == 0))
    {
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get_gsl_ptr(), i);
            gsl_fft_complex_radix2_forward((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size);
        }
    }
    // If operating on columns and length of column is power of 2, use radix2 method
    else if (dim == 1 && ((x.nrows() & (x.nrows() - 1)) == 0))
    {
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get_gsl_ptr(), j);
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
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get_gsl_ptr(), j);
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
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get_gsl_ptr(), i);
            gsl_fft_complex_forward((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

gsl::cmatrix gsl::fft(const gsl::cmatrix &x, int dim)
{
    gsl::cmatrix y(x);
    return gsl::fft(std::move(y), dim);
}

gsl::cmatrix gsl::ifft( gsl::cmatrix&& x, int dim )
{
    gsl_matrix_complex *xp = x.get_gsl_ptr();

    // If operating on rows and length of row is power of 2, use radix2 method
    if (dim == 2 && ((x.ncols() & (x.ncols() - 1)) == 0))
    {
        for (int i = 0; i < x.nrows(); ++i)
        {
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get_gsl_ptr(), i);
            gsl_fft_complex_radix2_inverse((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size);
        }
    }
    // If operating on columns and length of column is power of 2, use radix2 method
    else if (dim == 1 && ((x.nrows() & (x.nrows() - 1)) == 0))
    {
        for (int j = 0; j < x.ncols(); ++j)
        {
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get_gsl_ptr(), j);
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
            gsl_vector_complex_view col = gsl_matrix_complex_column(x.get_gsl_ptr(), j);
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
            gsl_vector_complex_view row = gsl_matrix_complex_row(x.get_gsl_ptr(), i);
            gsl_fft_complex_inverse((&row.vector)->data, (&row.vector)->stride, (&row.vector)->size, wavetable, workspace);
        }
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
    }

    return std::move(x);
}

gsl::cmatrix gsl::ifft( const gsl::cmatrix& x, int dim )
{
    gsl::cmatrix y(x);
    return gsl::ifft(std::move(y), dim);
}