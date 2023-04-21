#define CATCH_CONFIG_MAIN
#include <yawg/core.h>
#include <yawg/utils.hpp>
#include <yawg/fft.h>
#include <yawg/legendre.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Use Catch2 to test gsl::fft and gsl::ifft
TEST_CASE("gsl::fft", "[gsl::fft]")
{
    using namespace gsl::complex_literals;
    gsl::cvector x = gsl::linspace(0, 1, 4) + 1.0_i * gsl::linspace(1, 0, 4);

    SECTION("Test fft")
    {
        gsl::cvector y = gsl::fft(x);

        REQUIRE(x.get_gsl_ptr() != nullptr);
        REQUIRE(x.get_gsl_ptr() != y.get_gsl_ptr());
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(2.0, 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(2.0, 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(4.0 / 3.0, 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(-2.0 / 3.0, 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(2.0 / 3.0, 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(-4.0 / 3.0, 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(0.0, 1e-15));
    }

    SECTION("Test ifft")
    {
        gsl::cvector y = gsl::ifft(x);

        REQUIRE(x.get_gsl_ptr() != nullptr);
        REQUIRE(x.get_gsl_ptr() != y.get_gsl_ptr());
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(0.5, 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(0.5, 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(-1.0 / 3.0, 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(0.0, 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(-1.0 / 6.0, 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(1.0 / 6.0, 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-15));
    }

    SECTION("Test that x = ifft(fft(x))")
    {
        gsl::cvector y = gsl::ifft(gsl::fft(x));
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(x.get(0).real(), 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(x.get(0).imag(), 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(x.get(1).real(), 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(x.get(1).imag(), 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(x.get(2).real(), 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(x.get(2).imag(), 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(x.get(3).real(), 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(x.get(3).imag(), 1e-15));
    }

    SECTION("Test that x = fft(ifft(x))")
    {
        gsl::cvector y = gsl::fft(gsl::ifft(x));
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(x.get(0).real(), 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(x.get(0).imag(), 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(x.get(1).real(), 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(x.get(1).imag(), 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(x.get(2).real(), 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(x.get(2).imag(), 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(x.get(3).real(), 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(x.get(3).imag(), 1e-15));
    }
}

// Use Catch2 to test gsl::fft and gsl::ifft when the input size is not a power of 2
TEST_CASE("gsl::fft with non-power-of-2 input size", "[gsl::fft]")
{
    using namespace gsl::complex_literals;

    gsl::cvector x = gsl::linspace(0, 1, 5) + 1.0_i * gsl::linspace(1, 0, 5);

    SECTION("Test fft")
    {
        gsl::cvector y = gsl::fft(x);

        REQUIRE(x.get_gsl_ptr() != nullptr);
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(2.5, 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(2.5, 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(0.235238700294483, 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(1.485238700294483, 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(-0.421925189854434, 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(0.828074810145566, 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(-0.828074810145566, 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(0.421925189854434, 1e-15));

        REQUIRE_THAT(y.get(4).real(), Catch::Matchers::WithinAbs(-1.485238700294483, 1e-15));
        REQUIRE_THAT(y.get(4).imag(), Catch::Matchers::WithinAbs(-0.235238700294483, 1e-15));
    }

    SECTION("Test ifft")
    {
        gsl::cvector y = gsl::ifft(x);

        REQUIRE(x.get_gsl_ptr() != nullptr);
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(0.5, 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(0.5, 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(-0.297047740058897, 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(-0.047047740058897, 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(-0.165614962029113, 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(0.084385037970887, 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(-0.084385037970887, 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(0.165614962029113, 1e-15));

        REQUIRE_THAT(y.get(4).real(), Catch::Matchers::WithinAbs(0.047047740058897, 1e-15));
        REQUIRE_THAT(y.get(4).imag(), Catch::Matchers::WithinAbs(0.297047740058897, 1e-15));
    }

    SECTION("Test that x = ifft(fft(x))")
    {
        gsl::cvector y = gsl::ifft(gsl::fft(x));
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(x.get(0).real(), 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(x.get(0).imag(), 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(x.get(1).real(), 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(x.get(1).imag(), 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(x.get(2).real(), 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(x.get(2).imag(), 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(x.get(3).real(), 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(x.get(3).imag(), 1e-15));

        REQUIRE_THAT(y.get(4).real(), Catch::Matchers::WithinAbs(x.get(4).real(), 1e-15));
        REQUIRE_THAT(y.get(4).imag(), Catch::Matchers::WithinAbs(x.get(4).imag(), 1e-15));
    }

    SECTION("Test that x = fft(ifft(x))")
    {
        gsl::cvector y = gsl::fft(gsl::ifft(x));
        REQUIRE_THAT(y.get(0).real(), Catch::Matchers::WithinAbs(x.get(0).real(), 1e-15));
        REQUIRE_THAT(y.get(0).imag(), Catch::Matchers::WithinAbs(x.get(0).imag(), 1e-15));

        REQUIRE_THAT(y.get(1).real(), Catch::Matchers::WithinAbs(x.get(1).real(), 1e-15));
        REQUIRE_THAT(y.get(1).imag(), Catch::Matchers::WithinAbs(x.get(1).imag(), 1e-15));

        REQUIRE_THAT(y.get(2).real(), Catch::Matchers::WithinAbs(x.get(2).real(), 1e-15));
        REQUIRE_THAT(y.get(2).imag(), Catch::Matchers::WithinAbs(x.get(2).imag(), 1e-15));

        REQUIRE_THAT(y.get(3).real(), Catch::Matchers::WithinAbs(x.get(3).real(), 1e-15));
        REQUIRE_THAT(y.get(3).imag(), Catch::Matchers::WithinAbs(x.get(3).imag(), 1e-15));

        REQUIRE_THAT(y.get(4).real(), Catch::Matchers::WithinAbs(x.get(4).real(), 1e-15));
        REQUIRE_THAT(y.get(4).imag(), Catch::Matchers::WithinAbs(x.get(4).imag(), 1e-15));
    }
}

// Use Catch2 to test gsl::fft( cmatrix )
TEST_CASE("gsl::fft on 3x4 matrix", "[gsl::fft]")
{
    gsl::cmatrix x(3, 4);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            x.set(i, j, gsl::complex(0.5 * i + 0.25 * j, 0.25 * i + 0.5 * j));

    SECTION("Test fft on each length 3 column (not power of 2)")
    {
        gsl::cmatrix y = gsl::fft(x, 1);
        REQUIRE_THAT(y.get(0, 0).real(), Catch::Matchers::WithinAbs(1.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 0).imag(), Catch::Matchers::WithinAbs(0.7500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 1).real(), Catch::Matchers::WithinAbs(2.2500000000, 1e-10));
        REQUIRE_THAT(y.get(0, 1).imag(), Catch::Matchers::WithinAbs(2.2500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 2).real(), Catch::Matchers::WithinAbs(3.0000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 2).imag(), Catch::Matchers::WithinAbs(3.7500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 3).real(), Catch::Matchers::WithinAbs(3.7500000000, 1e-10));
        REQUIRE_THAT(y.get(0, 3).imag(), Catch::Matchers::WithinAbs(5.2500000000, 1e-10));

        REQUIRE_THAT(y.get(1, 0).real(), Catch::Matchers::WithinAbs(-0.9665063509, 1e-10));
        REQUIRE_THAT(y.get(1, 0).imag(), Catch::Matchers::WithinAbs(0.0580127019, 1e-10));

        REQUIRE_THAT(y.get(1, 1).real(), Catch::Matchers::WithinAbs(-0.9665063509, 1e-10));
        REQUIRE_THAT(y.get(1, 1).imag(), Catch::Matchers::WithinAbs(0.0580127019, 1e-10));

        REQUIRE_THAT(y.get(1, 2).real(), Catch::Matchers::WithinAbs(-0.9665063509, 1e-10));
        REQUIRE_THAT(y.get(1, 2).imag(), Catch::Matchers::WithinAbs(0.0580127019, 1e-10));

        REQUIRE_THAT(y.get(1, 3).real(), Catch::Matchers::WithinAbs(-0.9665063509, 1e-10));
        REQUIRE_THAT(y.get(1, 3).imag(), Catch::Matchers::WithinAbs(0.0580127019, 1e-10));

        REQUIRE_THAT(y.get(2, 0).real(), Catch::Matchers::WithinAbs(-0.5334936491, 1e-10));
        REQUIRE_THAT(y.get(2, 0).imag(), Catch::Matchers::WithinAbs(-0.8080127019, 1e-10));

        REQUIRE_THAT(y.get(2, 1).real(), Catch::Matchers::WithinAbs(-0.5334936491, 1e-10));
        REQUIRE_THAT(y.get(2, 1).imag(), Catch::Matchers::WithinAbs(-0.8080127019, 1e-10));

        REQUIRE_THAT(y.get(2, 2).real(), Catch::Matchers::WithinAbs(-0.5334936491, 1e-10));
        REQUIRE_THAT(y.get(2, 2).imag(), Catch::Matchers::WithinAbs(-0.8080127019, 1e-10));

        REQUIRE_THAT(y.get(2, 3).real(), Catch::Matchers::WithinAbs(-0.5334936491, 1e-10));
        REQUIRE_THAT(y.get(2, 3).imag(), Catch::Matchers::WithinAbs(-0.8080127019, 1e-10));
    }

    SECTION("Test fft on each length 4 row (power of 2)")
    {
        gsl::cmatrix y = gsl::fft(x, 2);
        REQUIRE_THAT(y.get(0, 0).real(), Catch::Matchers::WithinAbs(1.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 0).imag(), Catch::Matchers::WithinAbs(3.0000000000, 1e-10));

        REQUIRE_THAT(y.get(0, 1).real(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 1).imag(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));

        REQUIRE_THAT(y.get(0, 2).real(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 2).imag(), Catch::Matchers::WithinAbs(-1.0000000000, 1e-10));

        REQUIRE_THAT(y.get(0, 3).real(), Catch::Matchers::WithinAbs(0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 3).imag(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));

        REQUIRE_THAT(y.get(1, 0).real(), Catch::Matchers::WithinAbs(3.5000000000, 1e-10));
        REQUIRE_THAT(y.get(1, 0).imag(), Catch::Matchers::WithinAbs(4.0000000000, 1e-10));

        REQUIRE_THAT(y.get(1, 1).real(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));
        REQUIRE_THAT(y.get(1, 1).imag(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));

        REQUIRE_THAT(y.get(1, 2).real(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(1, 2).imag(), Catch::Matchers::WithinAbs(-1.0000000000, 1e-10));

        REQUIRE_THAT(y.get(1, 3).real(), Catch::Matchers::WithinAbs(0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(1, 3).imag(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));

        REQUIRE_THAT(y.get(2, 0).real(), Catch::Matchers::WithinAbs(5.5000000000, 1e-10));
        REQUIRE_THAT(y.get(2, 0).imag(), Catch::Matchers::WithinAbs(5.0000000000, 1e-10));

        REQUIRE_THAT(y.get(2, 1).real(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));
        REQUIRE_THAT(y.get(2, 1).imag(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));

        REQUIRE_THAT(y.get(2, 2).real(), Catch::Matchers::WithinAbs(-0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(2, 2).imag(), Catch::Matchers::WithinAbs(-1.0000000000, 1e-10));

        REQUIRE_THAT(y.get(2, 3).real(), Catch::Matchers::WithinAbs(0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(2, 3).imag(), Catch::Matchers::WithinAbs(-1.5000000000, 1e-10));
    }
}

// Use Catch2 to test gsl::ifft( cmatrix ) function
TEST_CASE("gsl::ifft on 3x4 matrix", "[gsl::ifft]")
{
    gsl::cmatrix x(3, 4);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            x.set(i, j, gsl::complex(0.5 * i + 0.25 * j, 0.25 * i + 0.5 * j));

    SECTION("Test fft on each length 3 column (not power of 2)")
    {
        gsl::cmatrix y = gsl::ifft(x, 1);
        REQUIRE_THAT(y.get(0, 0).real(), Catch::Matchers::WithinAbs(0.5000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 0).imag(), Catch::Matchers::WithinAbs(0.2500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 1).real(), Catch::Matchers::WithinAbs(0.7500000000, 1e-10));
        REQUIRE_THAT(y.get(0, 1).imag(), Catch::Matchers::WithinAbs(0.7500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 2).real(), Catch::Matchers::WithinAbs(1.0000000000, 1e-10));
        REQUIRE_THAT(y.get(0, 2).imag(), Catch::Matchers::WithinAbs(1.2500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 3).real(), Catch::Matchers::WithinAbs(1.2500000000, 1e-10));
        REQUIRE_THAT(y.get(0, 3).imag(), Catch::Matchers::WithinAbs(1.7500000000, 1e-10));

        REQUIRE_THAT(y.get(1, 0).real(), Catch::Matchers::WithinAbs(-0.1778312164, 1e-10));
        REQUIRE_THAT(y.get(1, 0).imag(), Catch::Matchers::WithinAbs(-0.2693375673, 1e-10));

        REQUIRE_THAT(y.get(1, 1).real(), Catch::Matchers::WithinAbs(-0.1778312164, 1e-10));
        REQUIRE_THAT(y.get(1, 1).imag(), Catch::Matchers::WithinAbs(-0.2693375673, 1e-10));

        REQUIRE_THAT(y.get(1, 2).real(), Catch::Matchers::WithinAbs(-0.1778312164, 1e-10));
        REQUIRE_THAT(y.get(1, 2).imag(), Catch::Matchers::WithinAbs(-0.2693375673, 1e-10));

        REQUIRE_THAT(y.get(1, 3).real(), Catch::Matchers::WithinAbs(-0.1778312164, 1e-10));
        REQUIRE_THAT(y.get(1, 3).imag(), Catch::Matchers::WithinAbs(-0.2693375673, 1e-10));

        REQUIRE_THAT(y.get(2, 0).real(), Catch::Matchers::WithinAbs(-0.3221687836, 1e-10));
        REQUIRE_THAT(y.get(2, 0).imag(), Catch::Matchers::WithinAbs(0.0193375673, 1e-10));

        REQUIRE_THAT(y.get(2, 1).real(), Catch::Matchers::WithinAbs(-0.3221687836, 1e-10));
        REQUIRE_THAT(y.get(2, 1).imag(), Catch::Matchers::WithinAbs(0.0193375673, 1e-10));

        REQUIRE_THAT(y.get(2, 2).real(), Catch::Matchers::WithinAbs(-0.3221687836, 1e-10));
        REQUIRE_THAT(y.get(2, 2).imag(), Catch::Matchers::WithinAbs(0.0193375673, 1e-10));

        REQUIRE_THAT(y.get(2, 3).real(), Catch::Matchers::WithinAbs(-0.3221687836, 1e-10));
        REQUIRE_THAT(y.get(2, 3).imag(), Catch::Matchers::WithinAbs(0.0193375673, 1e-10));
    }

    SECTION("Test fft on each length 4 row (power of 2)")
    {
        gsl::cmatrix y = gsl::ifft(x, 2);
        REQUIRE_THAT(y.get(0, 0).real(), Catch::Matchers::WithinAbs(0.3750000000, 1e-10));
        REQUIRE_THAT(y.get(0, 0).imag(), Catch::Matchers::WithinAbs(0.7500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 1).real(), Catch::Matchers::WithinAbs(0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(0, 1).imag(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));

        REQUIRE_THAT(y.get(0, 2).real(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(0, 2).imag(), Catch::Matchers::WithinAbs(-0.2500000000, 1e-10));

        REQUIRE_THAT(y.get(0, 3).real(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));
        REQUIRE_THAT(y.get(0, 3).imag(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));

        REQUIRE_THAT(y.get(1, 0).real(), Catch::Matchers::WithinAbs(0.8750000000, 1e-10));
        REQUIRE_THAT(y.get(1, 0).imag(), Catch::Matchers::WithinAbs(1.0000000000, 1e-10));

        REQUIRE_THAT(y.get(1, 1).real(), Catch::Matchers::WithinAbs(0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(1, 1).imag(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));

        REQUIRE_THAT(y.get(1, 2).real(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(1, 2).imag(), Catch::Matchers::WithinAbs(-0.2500000000, 1e-10));

        REQUIRE_THAT(y.get(1, 3).real(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));
        REQUIRE_THAT(y.get(1, 3).imag(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));

        REQUIRE_THAT(y.get(2, 0).real(), Catch::Matchers::WithinAbs(1.3750000000, 1e-10));
        REQUIRE_THAT(y.get(2, 0).imag(), Catch::Matchers::WithinAbs(1.2500000000, 1e-10));

        REQUIRE_THAT(y.get(2, 1).real(), Catch::Matchers::WithinAbs(0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(2, 1).imag(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));

        REQUIRE_THAT(y.get(2, 2).real(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));
        REQUIRE_THAT(y.get(2, 2).imag(), Catch::Matchers::WithinAbs(-0.2500000000, 1e-10));

        REQUIRE_THAT(y.get(2, 3).real(), Catch::Matchers::WithinAbs(-0.3750000000, 1e-10));
        REQUIRE_THAT(y.get(2, 3).imag(), Catch::Matchers::WithinAbs(-0.1250000000, 1e-10));
    }
}

// Use Catch2 to test legendre_P
TEST_CASE("legendre_P", "[legendre_P]")
{
    gsl::vector x = gsl::linspace(-0.9, 0.9, 10);
    gsl::vector y = gsl::legendre_P(2, x);
    
    y.print();
}