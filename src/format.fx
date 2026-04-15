#import "numeric.fx";
#import "standard.fx";

using standard::io::console;
using standard::strings;

namespace cmath
{
    def to_string(Complex z) -> string
    {
        if (z.imag == 0.0)
        {
            return f"{z.real}\0";
        };
        if (z.real == 0.0)
        {
            return f"{z.imag}i\0";
        };
        if (z.imag >= 0.0)
        {
            return f"{z.real}+{z.imag}i\0";
        };
        return f"{z.real}{z.imag}i\0";
    };

    def to_string_polar(Complex z) -> string
    {
        PolarResult p = polar(z);
        return f"{p.r}∠{p.phi}\0";
    };

    def to_string_exp(Complex z) -> string
    {
        PolarResult p = polar(z);
        return f"{p.r}·e^(i{p.phi})\0";
    };

    def to_string_rect(double r, double phi) -> string
    {
        return f"({r}, {phi})\0";
    };

    def print_c(Complex z) -> void
    {
        println(to_string(z));
        return;
    };

    def print_polar(Complex z) -> void
    {
        println(to_string_polar(z));
        return;
    };

    def print_matrix(Complex2x2 A) -> void
    {
        string a = to_string(A.a);
        string b = to_string(A.b);
        string c = to_string(A.c);
        string d = to_string(A.d);
        println(f"[{a}, {b}]\0");
        println(f"[{c}, {d}]\0");
        return;
    };
};
