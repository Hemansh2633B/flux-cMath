#import "fresnel.fx";

namespace cmath
{
    def scale_pi(Complex z) -> Complex
    {
        return scale(z, PI);
    };

    def dft_twiddle(int k, int N) -> Complex
    {
        if (N == 0)
        {
            throw(CmathError("dft_twiddle: N cannot be zero\0"));
        };
        double angle = 0.0 - TAU * (double)k / (double)N;
        return rect(1.0, angle);
    };

    def idft_twiddle(int k, int N) -> Complex
    {
        if (N == 0)
        {
            throw(CmathError("idft_twiddle: N cannot be zero\0"));
        };
        double angle = TAU * (double)k / (double)N;
        return rect(1.0, angle);
    };

    def sinc(double x) -> double
    {
        if (x == 0.0)
        {
            return 1.0;
        };
        double px = PI * x;
        return sin(px) / px;
    };

    def sinc_c(Complex z) -> Complex
    {
        Complex pz = scale_pi(z);
        if (iszero(pz))
        {
            return ONE;
        };
        return div(sin(pz), pz);
    };

    def rect_window(int n, int N) -> double
    {
        if (n >= 0 & n < N)
        {
            return 1.0;
        };
        return 0.0;
    };

    def hann_window(int n, int N) -> double
    {
        if (N <= 1)
        {
            return 1.0;
        };
        double dn = (double)n;
        double dN = (double)(N - 1);
        return 0.5 * (1.0 - cos(TAU * dn / dN));
    };

    def hamming_window(int n, int N) -> double
    {
        if (N <= 1)
        {
            return 1.0;
        };
        double dn = (double)n;
        double dN = (double)(N - 1);
        return 0.54 - 0.46 * cos(TAU * dn / dN);
    };

    def blackman_window(int n, int N) -> double
    {
        if (N <= 1)
        {
            return 1.0;
        };
        double dn = (double)n;
        double dN = (double)(N - 1);
        double x = TAU * dn / dN;
        return 0.42 - 0.5 * cos(x) + 0.08 * cos(2.0 * x);
    };

    def phase_shift(Complex z, double phi) -> Complex
    {
        return mul(z, rect(1.0, phi));
    };

    def instantaneous_freq(Complex z1, Complex z2, double dt) -> double
    {
        if (dt == 0.0)
        {
            throw(CmathError("instantaneous_freq: dt cannot be zero\0"));
        };
        double dphi = arg(z2) - arg(z1);
        while (dphi > PI)
        {
            dphi = dphi - TAU;
        };
        while (dphi < 0.0 - PI)
        {
            dphi = dphi + TAU;
        };
        return dphi / (TAU * dt);
    };

    def group_delay(Complex H1, Complex H2, double dw) -> double
    {
        if (dw == 0.0)
        {
            throw(CmathError("group_delay: dw cannot be zero\0"));
        };
        double dphi = arg(H2) - arg(H1);
        while (dphi > PI)
        {
            dphi = dphi - TAU;
        };
        while (dphi < 0.0 - PI)
        {
            dphi = dphi + TAU;
        };
        return -1.0 * dphi / dw;
    };
};
