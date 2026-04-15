#import "lambert.fx";

namespace cmath
{
    struct FresnelResult
    {
        double s;
        double c;
    };

    def fresnel(double x) -> FresnelResult
    {
        FresnelResult result;

        if (x == 0.0)
        {
            result.s = 0.0;
            result.c = 0.0;
            return result;
        };

        bool neg_x = x < 0.0;
        double ax = neg_x ? (0.0 - x) : x;

        if (ax < 1.6)
        {
            double pi_half = PI * 0.5;
            double s_sum = 0.0;
            double c_sum = 0.0;
            double sign_val = 1.0;

            for (int n = 0; n < 20; n++)
            {
                double dn = (double)n;

                double c_pow_pi = pow(pi_half, 2.0 * dn);
                double c_pow_x = pow(ax, 4.0 * dn + 1.0);
                double c_fact = tgamma(2.0 * dn + 1.0);
                double c_denom = c_fact * (4.0 * dn + 1.0);
                c_sum = c_sum + sign_val * c_pow_pi * c_pow_x / c_denom;

                double s_pow_pi = pow(pi_half, 2.0 * dn + 1.0);
                double s_pow_x = pow(ax, 4.0 * dn + 3.0);
                double s_fact = tgamma(2.0 * dn + 2.0);
                double s_denom = s_fact * (4.0 * dn + 3.0);
                s_sum = s_sum + sign_val * s_pow_pi * s_pow_x / s_denom;

                sign_val = 0.0 - sign_val;
            };

            result.s = s_sum;
            result.c = c_sum;
        }
        else
        {
            double pix = PI * ax;
            double inv_pix = 1.0 / pix;
            double inv_pix2 = inv_pix * inv_pix;
            double inv_pix4 = inv_pix2 * inv_pix2;
            double f = inv_pix * (1.0 - 2.0 * inv_pix2 + 24.0 * inv_pix4);
            double g_pref = 1.0 / (PI * PI * ax * ax * ax);
            double g = g_pref * (1.0 - 6.0 * inv_pix2 + 120.0 * inv_pix4);
            double theta = 0.5 * PI * ax * ax;

            result.s = 0.5 - f * cos(theta) - g * sin(theta);
            result.c = 0.5 + f * sin(theta) - g * cos(theta);
        };

        if (neg_x)
        {
            result.s = 0.0 - result.s;
            result.c = 0.0 - result.c;
        };

        return result;
    };

    def fresnel_s(double x) -> double
    {
        return fresnel(x).s;
    };

    def fresnel_c(double x) -> double
    {
        return fresnel(x).c;
    };

    def fresnel_c_cplx(Complex z) -> Complex
    {
        Complex arg_plus = mul(scale(sub(ONE, I), 0.5 * SQRT_PI), z);
        Complex arg_minus = mul(scale(add(ONE, I), 0.5 * SQRT_PI), z);
        Complex f_plus = mul(scale(add(ONE, I), 0.5), erf_c(arg_plus));
        Complex f_minus = mul(scale(sub(ONE, I), 0.5), erf_c(arg_minus));
        return scale(add(f_plus, f_minus), 0.5);
    };

    def fresnel_s_cplx(Complex z) -> Complex
    {
        Complex arg_plus = mul(scale(sub(ONE, I), 0.5 * SQRT_PI), z);
        Complex arg_minus = mul(scale(add(ONE, I), 0.5 * SQRT_PI), z);
        Complex f_plus = mul(scale(add(ONE, I), 0.5), erf_c(arg_plus));
        Complex f_minus = mul(scale(sub(ONE, I), 0.5), erf_c(arg_minus));
        Complex diff = sub(f_plus, f_minus);
        return mul(Complex(0.0, -0.5), diff);
    };
};
