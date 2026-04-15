#import "erf.fx";

namespace cmath
{
    def besselj0(double x) -> double
    {
        return j0(x);
    };

    def besselj1(double x) -> double
    {
        return j1(x);
    };

    def besselj(int n, double x) -> double
    {
        return jn(n, x);
    };

    def besselj_nu(Complex nu, Complex z) -> Complex
    {
        Complex half_z = scale(z, 0.5);
        Complex term = pow_c(half_z, nu);
        Complex z2neg = neg(mul(half_z, half_z));
        Complex sum = ZERO;
        double fact_k = 1.0;

        for (int k = 0; k < 40; k++)
        {
            Complex gamma_term = gamma_c(add(nu, Complex((double)k + 1.0, 0.0)));
            Complex denom = mul(Complex(fact_k, 0.0), gamma_term);
            Complex contrib = div(term, denom);
            sum = add(sum, contrib);
            term = mul(term, z2neg);
            fact_k = fact_k * (double)(k + 1);
        };

        return sum;
    };

    def besselj_c(int n, Complex z) -> Complex
    {
        if (n < 0)
        {
            int m = 0 - n;
            Complex jm = besselj_c(m, z);
            if ((m % 2) == 0)
            {
                return jm;
            };
            return neg(jm);
        };

        if (iszero(z))
        {
            if (n == 0)
            {
                return ONE;
            };
            return ZERO;
        };

        Complex half_z = scale(z, 0.5);
        Complex pow_hz = pow_c(half_z, Complex((double)n, 0.0));
        Complex sum = ZERO;
        double fact_k = 1.0;
        Complex z2neg = neg(mul(half_z, half_z));
        Complex term = pow_hz;
        double gam_nk = tgamma((double)(n + 1));

        for (int k = 0; k < 40; k++)
        {
            Complex contrib = scale(term, 1.0 / (fact_k * gam_nk));
            sum = add(sum, contrib);
            term = mul(term, z2neg);
            fact_k = fact_k * (double)(k + 1);
            gam_nk = gam_nk * (double)(n + k + 1);
        };

        return sum;
    };

    def bessely0(double x) -> double
    {
        return y0(x);
    };

    def bessely1(double x) -> double
    {
        return y1(x);
    };

    def bessely(int n, double x) -> double
    {
        return yn(n, x);
    };

    def bessely_c(int n, Complex z) -> Complex
    {
        if (z.imag == 0.0)
        {
            return Complex(yn(n, z.real), 0.0);
        };

        double eps = 0.0000001;
        double nu_real = (double)n + eps;
        Complex nu = Complex(nu_real, 0.0);
        Complex j_pos = besselj_nu(nu, z);
        Complex j_neg = besselj_nu(neg(nu), z);
        double angle = PI * nu_real;
        Complex num = sub(scale(j_pos, cos(angle)), j_neg);
        double den = sin(angle);
        return scale(num, 1.0 / den);
    };

    def besseli_c(int n, Complex z) -> Complex
    {
        int m = n;
        if (m < 0)
        {
            m = 0 - m;
        };
        Complex iz = mul(I, z);
        Complex pref = pow_c(Complex(0.0, -1.0), Complex((double)m, 0.0));
        return mul(pref, besselj_c(m, iz));
    };

    def besselk0_series(Complex z) -> Complex
    {
        Complex half_z = scale(z, 0.5);
        Complex z2_over4 = mul(half_z, half_z);
        Complex pow_term = ONE;
        double fact = 1.0;
        Complex sum = ZERO;

        for (int k = 0; k < 30; k++)
        {
            Complex psi = digamma_c(Complex((double)k + 1.0, 0.0));
            Complex contrib = scale(mul(psi, pow_term), 1.0 / (fact * fact));
            sum = add(sum, contrib);
            pow_term = mul(pow_term, z2_over4);
            fact = fact * (double)(k + 1);
        };

        Complex lead = neg(mul(log(half_z), besseli_c(0, z)));
        return add(lead, sum);
    };

    def besselk1_series(Complex z) -> Complex
    {
        Complex half_z = scale(z, 0.5);
        Complex z2_over4 = mul(half_z, half_z);
        Complex pow_term = half_z;
        double fact_k = 1.0;
        double fact_kp1 = 1.0;
        Complex sum = ZERO;

        for (int k = 0; k < 30; k++)
        {
            Complex psi1 = digamma_c(Complex((double)k + 1.0, 0.0));
            Complex psi2 = digamma_c(Complex((double)k + 2.0, 0.0));
            Complex psi_sum = add(psi1, psi2);
            Complex contrib = scale(mul(psi_sum, pow_term), 1.0 / (fact_k * fact_kp1));
            sum = add(sum, contrib);
            pow_term = mul(pow_term, z2_over4);
            fact_k = fact_k * (double)(k + 1);
            fact_kp1 = fact_kp1 * (double)(k + 2);
        };

        Complex term0 = recip(z);
        Complex term1 = mul(log(half_z), besseli_c(1, z));
        Complex term2 = scale(sum, -0.5);
        return add(add(term0, term1), term2);
    };

    def besselk_c(int n, Complex z) -> Complex
    {
        if (iszero(z))
        {
            throw(CmathError("besselk_c: singular at zero\0"));
        };

        int m = n;
        if (m < 0)
        {
            m = 0 - m;
        };

        if (m == 0)
        {
            return besselk0_series(z);
        };

        if (m == 1)
        {
            return besselk1_series(z);
        };

        Complex k_prev = besselk0_series(z);
        Complex k_curr = besselk1_series(z);

        for (int order = 1; order < m; order++)
        {
            Complex coeff = div(Complex(2.0 * (double)order, 0.0), z);
            Complex k_next = add(mul(coeff, k_curr), k_prev);
            k_prev = k_curr;
            k_curr = k_next;
        };

        return k_curr;
    };
};
