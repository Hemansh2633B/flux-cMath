#import "polar.fx";

namespace cmath
{
    struct RootsArray
    {
        Complex[16] roots;
        int count;
    };

    def exp(Complex z) -> Complex
    {
        double scale_val = exp(z.real);
        return Complex(scale_val * cos(z.imag), scale_val * sin(z.imag));
    };

    def log(Complex z) -> Complex
    {
        if (z.real == 0.0 & z.imag == 0.0)
        {
            return Complex(0.0 - INF, 0.0);
        };
        double radius = hypot(z.real, z.imag);
        double angle = atan2(z.imag, z.real);
        return Complex(log(radius), angle);
    };

    def log(Complex z, double base) -> Complex
    {
        Complex numerator = log(z);
        Complex denominator = log(Complex(base, 0.0));
        return div(numerator, denominator);
    };

    def log10(Complex z) -> Complex
    {
        return div(log(z), log(Complex(10.0, 0.0)));
    };

    def log2(Complex z) -> Complex
    {
        return div(log(z), log(Complex(2.0, 0.0)));
    };

    def sqrt(Complex z) -> Complex
    {
        if (z.real == 0.0 & z.imag == 0.0)
        {
            return ZERO;
        };

        double radius = hypot(z.real, z.imag);
        double w = sqrt((radius + fabs(z.real)) / 2.0);

        if (w == 0.0)
        {
            return Complex(0.0, sqrt(fabs(z.imag)));
        };

        if (z.real >= 0.0 & z.imag >= 0.0)
        {
            return Complex(w, z.imag / (2.0 * w));
        };

        if (z.real >= 0.0 & z.imag < 0.0)
        {
            return Complex(w, z.imag / (2.0 * w));
        };

        if (z.real < 0.0 & z.imag >= 0.0)
        {
            return Complex(fabs(z.imag) / (2.0 * w), w);
        };

        return Complex(fabs(z.imag) / (2.0 * w), 0.0 - w);
    };

    def cbrt(Complex z) -> Complex
    {
        double radius = hypot(z.real, z.imag);
        double phi = atan2(z.imag, z.real);
        double cr = pow(radius, 1.0 / 3.0);
        double ang = phi / 3.0;
        return Complex(cr * cos(ang), cr * sin(ang));
    };

    def root_n(Complex z, int n) -> Complex
    {
        if (n == 0)
        {
            throw(CmathError("root_n: n cannot be zero\0"));
        };
        double radius = hypot(z.real, z.imag);
        double phi = atan2(z.imag, z.real);
        double dn = (double)n;
        double rn = pow(radius, 1.0 / dn);
        double ang = phi / dn;
        return Complex(rn * cos(ang), rn * sin(ang));
    };

    def roots_of_unity(int n) -> RootsArray
    {
        if (n < 1)
        {
            throw(CmathError("roots_of_unity: n must be at least 1\0"));
        };
        if (n > 16)
        {
            throw(CmathError("roots_of_unity: n exceeds 16\0"));
        };

        RootsArray result;
        result.count = n;
        for (int k = 0; k < n; k++)
        {
            double angle = TAU * (double)k / (double)n;
            result.roots[k] = Complex(cos(angle), sin(angle));
        };
        return result;
    };

    def pow_c(Complex base, Complex exponent) -> Complex
    {
        if (base.real == 0.0 & base.imag == 0.0)
        {
            if (exponent.real == 0.0 & exponent.imag == 0.0)
            {
                return ONE;
            };
            if (exponent.imag == 0.0 & exponent.real > 0.0)
            {
                return ZERO;
            };
            throw(CmathError("pow_c: undefined for zero base and non-positive exponent\0"));
        };
        return exp(mul(exponent, log(base)));
    };

    def pow_cd(Complex base, double exponent) -> Complex
    {
        if (base.real == 0.0 & base.imag == 0.0)
        {
            if (exponent == 0.0)
            {
                return ONE;
            };
            if (exponent > 0.0)
            {
                return ZERO;
            };
            throw(CmathError("pow_cd: undefined for zero base and non-positive exponent\0"));
        };
        return exp(scale(log(base), exponent));
    };

    def pow_dc(double base, Complex exponent) -> Complex
    {
        if (base == 0.0)
        {
            if (exponent.real == 0.0 & exponent.imag == 0.0)
            {
                return ONE;
            };
            if (exponent.imag == 0.0 & exponent.real > 0.0)
            {
                return ZERO;
            };
            throw(CmathError("pow_dc: undefined for zero base and non-positive exponent\0"));
        };
        return exp(mul(exponent, log(Complex(base, 0.0))));
    };
};
