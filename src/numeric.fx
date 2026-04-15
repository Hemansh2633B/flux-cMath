#import "mobius.fx";

namespace cmath
{
    def deriv(byte* f_ptr, Complex z, double h) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        Complex hz = Complex(h, 0.0);
        Complex fp = f(add(z, hz));
        Complex fm = f(sub(z, hz));
        return scale(sub(fp, fm), 1.0 / (2.0 * h));
    };

    def deriv2(byte* f_ptr, Complex z, double h) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        Complex hz = Complex(h, 0.0);
        Complex fp = f(add(z, hz));
        Complex f0 = f(z);
        Complex fm = f(sub(z, hz));
        Complex num = add(sub(fp, scale(f0, 2.0)), fm);
        return scale(num, 1.0 / (h * h));
    };

    def newton_root(byte* f_ptr, byte* df_ptr,
                    Complex z0, int max_iter, double tol) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        def{}* df(Complex) -> Complex = (byte*)df_ptr;
        Complex z = z0;
        for (int i = 0; i < max_iter; i++)
        {
            Complex fz = f(z);
            Complex dfz = df(z);
            if (iszero(dfz))
            {
                throw(CmathError("newton: zero derivative\0"));
            };
            Complex dz = div(fz, dfz);
            z = sub(z, dz);
            if (abs_c(dz) < tol)
            {
                return z;
            };
        };
        throw(CmathError("newton: did not converge\0"));
        return z;
    };

    def newton_root_auto(byte* f_ptr, Complex z0, int max_iter, double tol) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        Complex z = z0;
        double h = 0.000001;

        for (int i = 0; i < max_iter; i++)
        {
            Complex fz = f(z);
            Complex dfz = deriv(f_ptr, z, h);
            if (iszero(dfz))
            {
                throw(CmathError("newton: zero derivative\0"));
            };
            Complex dz = div(fz, dfz);
            z = sub(z, dz);
            if (abs_c(dz) < tol)
            {
                return z;
            };
        };

        throw(CmathError("newton: did not converge\0"));
        return z;
    };

    def halley_root(byte* f_ptr, Complex z0, int max_iter, double tol) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        Complex z = z0;
        double h = 0.000001;

        for (int i = 0; i < max_iter; i++)
        {
            Complex fz = f(z);
            Complex dfz = deriv(f_ptr, z, h);
            if (iszero(dfz))
            {
                throw(CmathError("halley: zero derivative\0"));
            };

            Complex d2fz = deriv2(f_ptr, z, h);
            Complex num_h = mul(scale(fz, 2.0), dfz);
            Complex den_h = sub(scale(mul(dfz, dfz), 2.0), mul(fz, d2fz));
            Complex dz;

            if (abs_c(den_h) < tol)
            {
                dz = div(fz, dfz);
            }
            else
            {
                dz = div(num_h, den_h);
            };

            z = sub(z, dz);
            if (abs_c(dz) < tol)
            {
                return z;
            };
        };

        throw(CmathError("halley: did not converge\0"));
        return z;
    };

    def integrate_rect(byte* f_ptr, double a, double b, int n) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        if (n <= 0)
        {
            throw(CmathError("integrate_rect: n must be positive\0"));
        };

        double h = (b - a) / (double)n;
        Complex sum = ZERO;

        for (int i = 0; i < n; i++)
        {
            double x = a + (double)i * h;
            sum = add(sum, f(Complex(x, 0.0)));
        };

        return scale(sum, h);
    };

    def integrate_trapz(byte* f_ptr, double a, double b, int n) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        if (n <= 0)
        {
            throw(CmathError("integrate_trapz: n must be positive\0"));
        };

        double h = (b - a) / (double)n;
        Complex sum = add(scale(f(Complex(a, 0.0)), 0.5), scale(f(Complex(b, 0.0)), 0.5));

        for (int i = 1; i < n; i++)
        {
            double x = a + (double)i * h;
            sum = add(sum, f(Complex(x, 0.0)));
        };

        return scale(sum, h);
    };

    def integrate_simp(byte* f_ptr, double a, double b, int n) -> Complex
    {
        def{}* f(Complex) -> Complex = (byte*)f_ptr;
        if (n % 2 != 0)
        {
            n = n + 1;
        };
        double h = (b - a) / (double)n;
        Complex sum = add(f(Complex(a, 0.0)), f(Complex(b, 0.0)));
        for (int i = 1; i < n; i++)
        {
            double x = a + (double)i * h;
            double w = (i % 2 == 0) ? 2.0 : 4.0;
            sum = add(sum, scale(f(Complex(x, 0.0)), w));
        };
        return scale(sum, h / 3.0);
    };
};
