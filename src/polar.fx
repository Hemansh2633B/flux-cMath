#import "complex.fx";

namespace cmath
{
    struct PolarResult
    {
        double r;
        double phi;
    };

    def phase(Complex z) -> double
    {
        return atan2(z.imag, z.real);
    };

    def polar(Complex z) -> PolarResult
    {
        PolarResult result;
        result.r = hypot(z.real, z.imag);
        result.phi = atan2(z.imag, z.real);
        return result;
    };

    def rect(double r, double phi) -> Complex
    {
        return Complex(r * cos(phi), r * sin(phi));
    };

    def from_polar(PolarResult p) -> Complex
    {
        return rect(p.r, p.phi);
    };
};
