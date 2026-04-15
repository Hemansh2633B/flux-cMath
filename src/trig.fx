#import "power.fx";

namespace cmath
{
    def sin(Complex z) -> Complex
    {
        double a = z.real;
        double b = z.imag;
        return Complex(sin(a) * cosh(b), cos(a) * sinh(b));
    };

    def cos(Complex z) -> Complex
    {
        double a = z.real;
        double b = z.imag;
        return Complex(cos(a) * cosh(b), 0.0 - sin(a) * sinh(b));
    };

    def tan(Complex z) -> Complex
    {
        return div(sin(z), cos(z));
    };

    def asin(Complex z) -> Complex
    {
        Complex iz = Complex(0.0 - z.imag, z.real);
        Complex z2 = Complex(1.0 - z.real * z.real + z.imag * z.imag, 0.0 - 2.0 * z.real * z.imag);
        Complex sq = sqrt(z2);
        Complex inner = add(iz, sq);
        Complex l = log(inner);
        return Complex(l.imag, 0.0 - l.real);
    };

    def acos(Complex z) -> Complex
    {
        Complex z2 = Complex(1.0 - z.real * z.real + z.imag * z.imag, 0.0 - 2.0 * z.real * z.imag);
        Complex sq = sqrt(z2);
        Complex isq = Complex(0.0 - sq.imag, sq.real);
        Complex inner = add(z, isq);
        Complex l = log(inner);
        return Complex(l.imag, 0.0 - l.real);
    };

    def atan(Complex z) -> Complex
    {
        Complex iz = Complex(0.0 - z.imag, z.real);
        Complex one = ONE;
        Complex lnum = log(sub(one, iz));
        Complex lden = log(add(one, iz));
        Complex diff = sub(lnum, lden);
        return Complex(0.0 - diff.imag * 0.5, diff.real * 0.5);
    };

    def sec(Complex z) -> Complex
    {
        return recip(cos(z));
    };

    def csc(Complex z) -> Complex
    {
        return recip(sin(z));
    };

    def cot(Complex z) -> Complex
    {
        return div(cos(z), sin(z));
    };

    def asec(Complex z) -> Complex
    {
        return acos(recip(z));
    };

    def acsc(Complex z) -> Complex
    {
        return asin(recip(z));
    };

    def acot(Complex z) -> Complex
    {
        return atan(recip(z));
    };

    def versin(Complex z) -> Complex
    {
        return sub(ONE, cos(z));
    };

    def coversin(Complex z) -> Complex
    {
        return sub(ONE, sin(z));
    };

    def haversin(Complex z) -> Complex
    {
        return scale(versin(z), 0.5);
    };

    def exsec(Complex z) -> Complex
    {
        return sub(sec(z), ONE);
    };
};
