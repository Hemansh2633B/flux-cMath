#import "trig.fx";

namespace cmath
{
    def sinh(Complex z) -> Complex
    {
        double a = z.real;
        double b = z.imag;
        return Complex(sinh(a) * cos(b), cosh(a) * sin(b));
    };

    def cosh(Complex z) -> Complex
    {
        double a = z.real;
        double b = z.imag;
        return Complex(cosh(a) * cos(b), sinh(a) * sin(b));
    };

    def tanh(Complex z) -> Complex
    {
        return div(sinh(z), cosh(z));
    };

    def asinh(Complex z) -> Complex
    {
        return log(add(z, sqrt(add(mul(z, z), ONE))));
    };

    def acosh(Complex z) -> Complex
    {
        return log(add(z, sqrt(sub(mul(z, z), ONE))));
    };

    def atanh(Complex z) -> Complex
    {
        Complex num = sub(log(add(ONE, z)), log(sub(ONE, z)));
        return scale(num, 0.5);
    };

    def sech(Complex z) -> Complex
    {
        return recip(cosh(z));
    };

    def csch(Complex z) -> Complex
    {
        return recip(sinh(z));
    };

    def coth(Complex z) -> Complex
    {
        return div(cosh(z), sinh(z));
    };

    def asech(Complex z) -> Complex
    {
        return acosh(recip(z));
    };

    def acsch(Complex z) -> Complex
    {
        return asinh(recip(z));
    };

    def acoth(Complex z) -> Complex
    {
        return atanh(recip(z));
    };
};
