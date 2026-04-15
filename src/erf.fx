#import "special.fx";

namespace cmath
{
    def erf_c(Complex z) -> Complex
    {
        if (z.imag == 0.0)
        {
            return Complex(erf(z.real), 0.0);
        };

        if (z.real < 0.0)
        {
            return neg(erf_c(neg(z)));
        };

        double magnitude = abs_c(z);
        if (magnitude < 3.0)
        {
            Complex sum = ZERO;
            Complex zp = z;
            Complex z2neg = neg(mul(z, z));
            double fact = 1.0;
            for (int n = 0; n < 30; n++)
            {
                double denom = fact * (double)(2 * n + 1);
                Complex contrib = scale(zp, 1.0 / denom);
                sum = add(sum, contrib);
                zp = mul(zp, z2neg);
                fact = fact * (double)(n + 1);
            };
            return scale(sum, 2.0 / SQRT_PI);
        };

        Complex z2 = mul(z, z);
        Complex inv_z = recip(z);
        Complex inv_z2 = mul(inv_z, inv_z);
        Complex inv_z4 = mul(inv_z2, inv_z2);
        Complex inv_z6 = mul(inv_z4, inv_z2);

        Complex poly = ONE;
        poly = sub(poly, scale(inv_z2, 0.5));
        poly = add(poly, scale(inv_z4, 0.75));
        poly = sub(poly, scale(inv_z6, 1.875));

        Complex tail = mul(exp(neg(z2)), mul(inv_z, poly));
        tail = scale(tail, 1.0 / SQRT_PI);
        return sub(ONE, tail);
    };

    def erfc_c(Complex z) -> Complex
    {
        if (z.imag == 0.0)
        {
            return Complex(erfc(z.real), 0.0);
        };
        return sub(ONE, erf_c(z));
    };

    def erfi_c(Complex z) -> Complex
    {
        Complex iz = mul(I, z);
        Complex erf_iz = erf_c(iz);
        return mul(Complex(0.0, -1.0), erf_iz);
    };

    def erfcx_c(Complex z) -> Complex
    {
        return mul(exp(mul(z, z)), erfc_c(z));
    };
};
