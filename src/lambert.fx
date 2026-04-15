#import "zeta.fx";

namespace cmath
{
    def lambertw(Complex z) -> Complex
    {
        return lambertw(z, 0);
    };

    def lambertw(Complex z, int branch) -> Complex
    {
        if (iszero(z))
        {
            return ZERO;
        };

        Complex branch_point = Complex(0.0 - 1.0 / E, 0.0);
        if (isclose(z, branch_point, 0.00000000000001, 0.00000000000001))
        {
            return Complex(-1.0, 0.0);
        };

        Complex w;
        if (branch == 0)
        {
            double magnitude = abs_c(z);
            if (magnitude <= 0.5)
            {
                w = z;
            }
            elif (z.real > 0.0)
            {
                w = log(z);
            }
            else
            {
                w = Complex(log(magnitude), arg(z));
            };
        }
        else
        {
            w = add(log(z), Complex(0.0, TAU * (double)branch));
        };

        for (int i = 0; i < 6; i++)
        {
            Complex ew = exp(w);
            Complex wew = mul(w, ew);
            Complex num = sub(wew, z);
            Complex wp1 = add(w, ONE);
            Complex wp2 = add(w, Complex(2.0, 0.0));
            Complex fprime = mul(ew, wp1);

            Complex denom;
            if (iszero(wp1))
            {
                denom = fprime;
            }
            else
            {
                Complex corr = div(mul(num, wp2), scale(wp1, 2.0));
                denom = sub(fprime, corr);
            };

            if (iszero(denom))
            {
                if (iszero(fprime))
                {
                    break;
                };
                denom = fprime;
            };

            Complex step = div(num, denom);
            w = sub(w, step);
        };

        return w;
    };

    def lambertw_real(double x) -> double
    {
        double min_x = 0.0 - 1.0 / E;
        if (x < min_x)
        {
            throw(CmathError("lambertw_real: domain x >= -1/e\0"));
        };

        if (x == 0.0)
        {
            return 0.0;
        };

        if (x == min_x)
        {
            return -1.0;
        };

        double w;
        if (fabs(x) < 1.0)
        {
            w = x;
        }
        else
        {
            w = log(x);
        };

        for (int i = 0; i < 30; i++)
        {
            double ew = exp(w);
            double f = w * ew - x;
            double fp = ew * (w + 1.0);
            if (fp == 0.0)
            {
                break;
            };
            double step = f / fp;
            w = w - step;
            if (fabs(step) < 0.00000000000001)
            {
                break;
            };
        };

        return w;
    };
};
