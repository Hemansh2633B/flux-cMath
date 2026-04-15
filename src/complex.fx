#import "constants.fx";

namespace cmath
{
    object Complex
    {
        double real;
        double imag;

        def __init(double r, double i) -> this
        {
            this.real = r;
            this.imag = i;
            return this;
        };

        def __init(double r) -> this
        {
            this.real = r;
            this.imag = 0.0;
            return this;
        };

        def __exit() -> void
        {
            return;
        };
    };

    def add(Complex a, Complex b) -> Complex
    {
        return Complex(a.real + b.real, a.imag + b.imag);
    };

    def sub(Complex a, Complex b) -> Complex
    {
        return Complex(a.real - b.real, a.imag - b.imag);
    };

    def mul(Complex a, Complex b) -> Complex
    {
        double real_part = a.real * b.real - a.imag * b.imag;
        double imag_part = a.real * b.imag + a.imag * b.real;
        return Complex(real_part, imag_part);
    };

    def mul(double s, Complex z) -> Complex
    {
        return Complex(s * z.real, s * z.imag);
    };

    def mul(Complex z, double s) -> Complex
    {
        return Complex(z.real * s, z.imag * s);
    };

    def div(Complex a, Complex b) -> Complex
    {
        if (b.real == 0.0 & b.imag == 0.0)
        {
            throw(CmathError("cmath: division by zero\0"));
        };
        double denom = b.real * b.real + b.imag * b.imag;
        double real_num = a.real * b.real + a.imag * b.imag;
        double imag_num = a.imag * b.real - a.real * b.imag;
        return Complex(real_num / denom, imag_num / denom);
    };

    def neg(Complex z) -> Complex
    {
        return Complex(0.0 - z.real, 0.0 - z.imag);
    };

    def conj(Complex z) -> Complex
    {
        return Complex(z.real, 0.0 - z.imag);
    };

    def abs_c(Complex z) -> double
    {
        return hypot(z.real, z.imag);
    };

    def norm(Complex z) -> double
    {
        return z.real * z.real + z.imag * z.imag;
    };

    def arg(Complex z) -> double
    {
        return atan2(z.imag, z.real);
    };

    def sign(Complex z) -> Complex
    {
        double magnitude = abs_c(z);
        if (magnitude == 0.0)
        {
            return ZERO;
        };
        return Complex(z.real / magnitude, z.imag / magnitude);
    };

    def proj(Complex z) -> Complex
    {
        if (isinf(z.real) | isinf(z.imag))
        {
            return Complex(INF, copysign(0.0, z.imag));
        };
        return z;
    };

    def scale(Complex z, double s) -> Complex
    {
        return Complex(z.real * s, z.imag * s);
    };

    def recip(Complex z) -> Complex
    {
        double magnitude2 = norm(z);
        if (magnitude2 == 0.0)
        {
            throw(CmathError("cmath: division by zero\0"));
        };
        Complex conj_z = conj(z);
        return Complex(conj_z.real / magnitude2, conj_z.imag / magnitude2);
    };
};
