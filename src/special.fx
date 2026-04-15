#import "classify.fx";

namespace cmath
{
    def gamma_c(Complex z) -> Complex
    {
        if (z.imag == 0.0)
        {
            if (z.real <= 0.0)
            {
                double zr_floor = floor(z.real);
                if (z.real == zr_floor)
                {
                    throw(CmathError("gamma_c: pole at non-positive integer\0"));
                };
            };
        };

        double[9] p = [
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            0.0000099843695780195716,
            0.00000015056327351493116
        ];

        if (z.real < 0.5)
        {
            Complex pi_z = scale(z, PI);
            Complex sin_term = sin(pi_z);
            Complex g = gamma_c(sub(ONE, z));
            Complex den = mul(sin_term, g);
            return div(Complex(PI, 0.0), den);
        };

        Complex w = sub(z, ONE);
        Complex x = Complex(p[0], 0.0);
        for (int k = 1; k <= 8; k++)
        {
            Complex den = add(w, Complex((double)k, 0.0));
            Complex term = div(Complex(p[k], 0.0), den);
            x = add(x, term);
        };

        Complex t = add(w, Complex(7.5, 0.0));
        Complex root_tau = Complex(sqrt(TAU), 0.0);
        Complex t_pow = pow_c(t, add(w, Complex(0.5, 0.0)));
        Complex e_term = exp(neg(t));
        return mul(mul(root_tau, x), mul(t_pow, e_term));
    };

    def lgamma_c(Complex z) -> Complex
    {
        return log(gamma_c(z));
    };

    def factorial_c(int n) -> double
    {
        if (n < 0)
        {
            throw(CmathError("factorial_c: n must be non-negative\0"));
        };
        return tgamma((double)(n + 1));
    };

    def beta_c(Complex a, Complex b) -> Complex
    {
        Complex sum_lg = add(lgamma_c(a), lgamma_c(b));
        Complex lg_ab = lgamma_c(add(a, b));
        return exp(sub(sum_lg, lg_ab));
    };

    def digamma_c(Complex z) -> Complex
    {
        Complex w = z;
        Complex acc = ZERO;

        while (w.real < 6.0)
        {
            acc = sub(acc, recip(w));
            w = add(w, ONE);
        };

        Complex log_w = log(w);
        Complex inv = recip(w);
        Complex inv2 = mul(inv, inv);
        Complex inv4 = mul(inv2, inv2);
        Complex inv6 = mul(inv4, inv2);

        Complex result = sub(log_w, scale(inv, 0.5));
        result = sub(result, scale(inv2, 1.0 / 12.0));
        result = add(result, scale(inv4, 1.0 / 120.0));
        result = sub(result, scale(inv6, 1.0 / 252.0));
        return add(result, acc);
    };
};
