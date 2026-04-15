#import "matrix.fx";

namespace cmath
{
    struct Mobius
    {
        Complex a;
        Complex b;
        Complex c;
        Complex d;
    };

    def mobius_det(Mobius m) -> Complex
    {
        return sub(mul(m.a, m.d), mul(m.b, m.c));
    };

    def mobius_trace(Mobius m) -> Complex
    {
        return add(m.a, m.d);
    };

    def cross_ratio(Complex z1, Complex z2, Complex z3, Complex z4) -> Complex
    {
        Complex num = mul(sub(z1, z3), sub(z2, z4));
        Complex den = mul(sub(z1, z4), sub(z2, z3));
        return div(num, den);
    };

    def mobius_apply(Mobius m, Complex z) -> Complex
    {
        Complex num = add(mul(m.a, z), m.b);
        Complex den = add(mul(m.c, z), m.d);
        if (iszero(den))
        {
            throw(CmathError("mobius: pole\0"));
        };
        return div(num, den);
    };

    def mobius_compose(Mobius m1, Mobius m2) -> Mobius
    {
        Mobius out;
        out.a = add(mul(m1.a, m2.a), mul(m1.b, m2.c));
        out.b = add(mul(m1.a, m2.b), mul(m1.b, m2.d));
        out.c = add(mul(m1.c, m2.a), mul(m1.d, m2.c));
        out.d = add(mul(m1.c, m2.b), mul(m1.d, m2.d));
        return out;
    };

    def mobius_inverse(Mobius m) -> Mobius
    {
        Complex det = mobius_det(m);
        if (isclose(det, ZERO, 0.000000000001, 0.000000000001))
        {
            throw(CmathError("mobius_inverse: singular transform\0"));
        };
        Complex inv_det = recip(det);
        Mobius out;
        out.a = mul(m.d, inv_det);
        out.b = mul(neg(m.b), inv_det);
        out.c = mul(neg(m.c), inv_det);
        out.d = mul(m.a, inv_det);
        return out;
    };

    def mobius_fixed_points(Mobius m) -> Complex[2]
    {
        Complex[2] roots;
        Complex B = sub(m.d, m.a);

        if (iszero(m.c))
        {
            if (iszero(B))
            {
                roots[0] = NANJ;
                roots[1] = NANJ;
                return roots;
            };
            Complex z0 = div(m.b, B);
            roots[0] = z0;
            roots[1] = z0;
            return roots;
        };

        Complex disc = add(mul(B, B), scale(mul(m.c, m.b), 4.0));
        Complex sq = sqrt(disc);
        Complex two_c = scale(m.c, 2.0);
        Complex minus_B = neg(B);

        roots[0] = div(add(minus_B, sq), two_c);
        roots[1] = div(sub(minus_B, sq), two_c);
        return roots;
    };

    def mobius_to_canonical(Complex z1, Complex z2, Complex z3) -> Mobius
    {
        Mobius m;
        Complex alpha = sub(z2, z3);
        Complex gamma = sub(z2, z1);
        m.a = alpha;
        m.b = neg(mul(z1, alpha));
        m.c = gamma;
        m.d = neg(mul(z3, gamma));
        return m;
    };

    def mobius_from_three_points(Complex z1, Complex z2, Complex z3,
                                 Complex w1, Complex w2, Complex w3) -> Mobius
    {
        Mobius mz = mobius_to_canonical(z1, z2, z3);
        Mobius mw = mobius_to_canonical(w1, w2, w3);
        Mobius mw_inv = mobius_inverse(mw);
        return mobius_compose(mw_inv, mz);
    };

    def is_elliptic(Mobius m) -> bool
    {
        Complex det = mobius_det(m);
        if (isclose(det, ZERO, 0.000000000001, 0.000000000001))
        {
            return false;
        };
        Complex tr = mobius_trace(m);
        Complex ratio = div(mul(tr, tr), det);
        if (fabs(ratio.imag) > 0.000000001)
        {
            return false;
        };
        return ratio.real >= 0.0 & ratio.real < 4.0;
    };

    def is_hyperbolic(Mobius m) -> bool
    {
        Complex det = mobius_det(m);
        if (isclose(det, ZERO, 0.000000000001, 0.000000000001))
        {
            return false;
        };
        Complex tr = mobius_trace(m);
        Complex ratio = div(mul(tr, tr), det);
        if (fabs(ratio.imag) > 0.000000001)
        {
            return false;
        };
        return ratio.real > 4.0;
    };

    def is_parabolic(Mobius m) -> bool
    {
        Complex det = mobius_det(m);
        if (isclose(det, ZERO, 0.000000000001, 0.000000000001))
        {
            return false;
        };
        Complex tr = mobius_trace(m);
        Complex ratio = div(mul(tr, tr), det);
        return isclose(ratio, Complex(4.0, 0.0), 0.000000001, 0.000000001);
    };

    def is_loxodromic(Mobius m) -> bool
    {
        if (is_parabolic(m))
        {
            return false;
        };
        if (is_elliptic(m))
        {
            return false;
        };
        if (is_hyperbolic(m))
        {
            return false;
        };
        return true;
    };
};
