extern
{
    def !!
        sin(double) -> double,
        cos(double) -> double,
        tan(double) -> double,
        asin(double) -> double,
        acos(double) -> double,
        atan(double) -> double,
        atan2(double, double) -> double,
        sinh(double) -> double,
        cosh(double) -> double,
        tanh(double) -> double,
        asinh(double) -> double,
        acosh(double) -> double,
        atanh(double) -> double,
        exp(double) -> double,
        log(double) -> double,
        log10(double) -> double,
        log2(double) -> double,
        sqrt(double) -> double,
        pow(double, double) -> double,
        fabs(double) -> double,
        ceil(double) -> double,
        floor(double) -> double,
        fmod(double, double) -> double,
        hypot(double, double) -> double,
        copysign(double, double) -> double,
        tgamma(double) -> double,
        lgamma(double) -> double,
        erf(double) -> double,
        erfc(double) -> double,
        j0(double) -> double,
        j1(double) -> double,
        jn(int, double) -> double,
        y0(double) -> double,
        y1(double) -> double,
        yn(int, double) -> double,
        isfinite(double) -> bool,
        isinf(double) -> bool,
        isnan(double) -> bool;
};

namespace cmath
{
    const double PI = 3.141592653589793;
    const double E = 2.718281828459045;
    const double TAU = 6.283185307179586;
    const double LN2 = 0.6931471805599453;
    const double LN10 = 2.302585092994046;
    const double SQRT2 = 1.4142135623730951;
    const double SQRT3 = 1.7320508075688772;
    const double PHI = 1.618033988749895;
    const double SQRT_PI = 1.7724538509055159;
    const double GAMMA_EM = 0.5772156649015329;
    const double SQRT2_INV = 0.7071067811865476;
    const double INF = 1.0 / 0.0;
    const double NAN = 0.0 / 0.0;
    const Complex INFJ = Complex(0.0, INF);
    const Complex NANJ = Complex(0.0, NAN);
    const Complex I = Complex(0.0, 1.0);
    const Complex ZERO = Complex(0.0, 0.0);
    const Complex ONE = Complex(1.0, 0.0);

    object CmathError
    {
        byte[] msg;

        def __init(byte[] m) -> this
        {
            this.msg = m;
            return this;
        };

        def __exit() -> void
        {
            return;
        };
    };
};
