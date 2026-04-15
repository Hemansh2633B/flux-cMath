#import "hyperbolic.fx";

namespace cmath
{
    def isfinite(Complex z) -> bool
    {
        return isfinite(z.real) & isfinite(z.imag);
    };

    def isinf(Complex z) -> bool
    {
        return isinf(z.real) | isinf(z.imag);
    };

    def isnan(Complex z) -> bool
    {
        return isnan(z.real) | isnan(z.imag);
    };

    def isreal(Complex z) -> bool
    {
        return fabs(z.imag) == 0.0;
    };

    def isimag(Complex z) -> bool
    {
        return fabs(z.real) == 0.0 & fabs(z.imag) != 0.0;
    };

    def iszero(Complex z) -> bool
    {
        return z.real == 0.0 & z.imag == 0.0;
    };

    def isclose(Complex a, Complex b, double rel_tol, double abs_tol) -> bool
    {
        if (isnan(a) | isnan(b))
        {
            return false;
        };

        if (isinf(a) | isinf(b))
        {
            return a.real == b.real & a.imag == b.imag;
        };

        double diff = abs_c(sub(a, b));
        double mag_a = abs_c(a);
        double mag_b = abs_c(b);
        double larger = (mag_a > mag_b) ? mag_a : mag_b;
        double scaled = rel_tol * larger;
        double tol = (scaled > abs_tol) ? scaled : abs_tol;
        return diff <= tol;
    };

    def isclose(Complex a, Complex b) -> bool
    {
        return isclose(a, b, 0.000000001, 0.0);
    };

    def isunit(Complex z, double tol) -> bool
    {
        return fabs(abs_c(z) - 1.0) <= tol;
    };

    def quadrant(Complex z) -> int
    {
        if (z.real > 0.0 & z.imag > 0.0)
        {
            return 1;
        };
        if (z.real < 0.0 & z.imag > 0.0)
        {
            return 2;
        };
        if (z.real < 0.0 & z.imag < 0.0)
        {
            return 3;
        };
        if (z.real > 0.0 & z.imag < 0.0)
        {
            return 4;
        };
        return 0;
    };
};
