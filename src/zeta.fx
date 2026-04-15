#import "bessel.fx";

namespace cmath
{
    def eta(Complex s) -> Complex
    {
        Complex sum = ZERO;
        double sign_val = 1.0;
        for (int n = 1; n <= 64; n++)
        {
            Complex term = pow_dc((double)n, neg(s));
            sum = add(sum, scale(term, sign_val));
            sign_val = 0.0 - sign_val;
        };
        return sum;
    };

    def zeta(Complex s) -> Complex
    {
        if (isclose(s, ONE))
        {
            throw(CmathError("zeta: pole at s=1\0"));
        };
        Complex e = eta(s);
        Complex denom = sub(ONE, pow_dc(2.0, sub(ONE, s)));
        return div(e, denom);
    };

    def hurwitz_zeta(Complex s, double a) -> Complex
    {
        if (a <= 0.0)
        {
            throw(CmathError("hurwitz_zeta: a must be positive\0"));
        };
        Complex sum = ZERO;
        for (int n = 0; n < 64; n++)
        {
            double base = (double)n + a;
            Complex term = pow_dc(base, neg(s));
            sum = add(sum, term);
        };
        return sum;
    };

    def ei(double x) -> double
    {
        if (x == 0.0)
        {
            return 0.0 - INF;
        };

        if (x > 40.0)
        {
            double sum = 1.0;
            double term = 1.0;
            for (int k = 1; k <= 20; k++)
            {
                term = term * (double)k / x;
                sum = sum + term;
            };
            return exp(x) * sum / x;
        };

        double ax = fabs(x);
        double sum_series = 0.0;
        double term_pow = x;
        double fact = 1.0;

        for (int n = 1; n <= 80; n++)
        {
            fact = fact * (double)n;
            sum_series = sum_series + term_pow / ((double)n * fact);
            term_pow = term_pow * x;
        };

        return GAMMA_EM + log(ax) + sum_series;
    };

    def ei_c(Complex z) -> Complex
    {
        if (iszero(z))
        {
            return Complex(0.0 - INF, 0.0);
        };

        Complex sum = ZERO;
        Complex z_pow = z;
        double fact = 1.0;

        for (int n = 1; n <= 64; n++)
        {
            fact = fact * (double)n;
            Complex contrib = scale(z_pow, 1.0 / ((double)n * fact));
            sum = add(sum, contrib);
            z_pow = mul(z_pow, z);
        };

        Complex lead = add(Complex(GAMMA_EM, 0.0), log(z));
        return add(lead, sum);
    };

    def li(Complex z) -> Complex
    {
        return ei_c(log(z));
    };
};
