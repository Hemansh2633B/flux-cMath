#import "signal.fx";

namespace cmath
{
    struct Complex2x2
    {
        Complex a;
        Complex b;
        Complex c;
        Complex d;
    };

    const Complex2x2 PAULI_X = {Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 0.0), Complex(0.0, 0.0)};
    const Complex2x2 PAULI_Y = {Complex(0.0, 0.0), Complex(0.0, -1.0), Complex(0.0, 1.0), Complex(0.0, 0.0)};
    const Complex2x2 PAULI_Z = {Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(-1.0, 0.0)};
    const Complex2x2 HADAMARD = {Complex(SQRT2_INV, 0.0), Complex(SQRT2_INV, 0.0), Complex(SQRT2_INV, 0.0), Complex(-SQRT2_INV, 0.0)};
    const Complex2x2 IDENTITY = {ONE, ZERO, ZERO, ONE};

    def mat_add(Complex2x2 A, Complex2x2 B) -> Complex2x2
    {
        Complex2x2 out = {add(A.a, B.a), add(A.b, B.b), add(A.c, B.c), add(A.d, B.d)};
        return out;
    };

    def mat_sub(Complex2x2 A, Complex2x2 B) -> Complex2x2
    {
        Complex2x2 out = {sub(A.a, B.a), sub(A.b, B.b), sub(A.c, B.c), sub(A.d, B.d)};
        return out;
    };

    def mat_mul(Complex2x2 A, Complex2x2 B) -> Complex2x2
    {
        Complex r00 = add(mul(A.a, B.a), mul(A.b, B.c));
        Complex r01 = add(mul(A.a, B.b), mul(A.b, B.d));
        Complex r10 = add(mul(A.c, B.a), mul(A.d, B.c));
        Complex r11 = add(mul(A.c, B.b), mul(A.d, B.d));
        Complex2x2 out = {r00, r01, r10, r11};
        return out;
    };

    def mat_scale(Complex2x2 A, Complex s) -> Complex2x2
    {
        Complex2x2 out = {mul(A.a, s), mul(A.b, s), mul(A.c, s), mul(A.d, s)};
        return out;
    };

    def mat_det(Complex2x2 A) -> Complex
    {
        return sub(mul(A.a, A.d), mul(A.b, A.c));
    };

    def mat_trace(Complex2x2 A) -> Complex
    {
        return add(A.a, A.d);
    };

    def mat_inv(Complex2x2 A) -> Complex2x2
    {
        Complex det = mat_det(A);
        if (isclose(det, ZERO, 0.000000000001, 0.000000000001))
        {
            throw(CmathError("mat_inv: singular matrix\0"));
        };
        Complex inv_det = recip(det);
        Complex2x2 adj = {A.d, neg(A.b), neg(A.c), A.a};
        return mat_scale(adj, inv_det);
    };

    def mat_conj(Complex2x2 A) -> Complex2x2
    {
        Complex2x2 out = {conj(A.a), conj(A.b), conj(A.c), conj(A.d)};
        return out;
    };

    def mat_transpose(Complex2x2 A) -> Complex2x2
    {
        Complex2x2 out = {A.a, A.c, A.b, A.d};
        return out;
    };

    def mat_adjoint(Complex2x2 A) -> Complex2x2
    {
        return mat_conj(mat_transpose(A));
    };

    def mat_vec_mul(Complex2x2 A, Complex v0, Complex v1) -> Complex[2]
    {
        Complex[2] out;
        out[0] = add(mul(A.a, v0), mul(A.b, v1));
        out[1] = add(mul(A.c, v0), mul(A.d, v1));
        return out;
    };

    def mat_pow(Complex2x2 A, int n) -> Complex2x2
    {
        if (n == 0)
        {
            return IDENTITY;
        };

        Complex2x2 base = A;
        int exp_n = n;

        if (n < 0)
        {
            base = mat_inv(A);
            exp_n = 0 - n;
        };

        Complex2x2 result = IDENTITY;

        while (exp_n > 0)
        {
            if ((exp_n % 2) == 1)
            {
                result = mat_mul(result, base);
            };
            base = mat_mul(base, base);
            exp_n = exp_n / 2;
        };

        return result;
    };

    def mat_exp(Complex2x2 A) -> Complex2x2
    {
        Complex tr = add(A.a, A.d);
        Complex half = scale(tr, 0.5);
        Complex2x2 A_shifted = {sub(A.a, half), A.b, A.c, sub(A.d, half)};
        Complex det = mat_det(A_shifted);
        Complex q = sqrt(det);
        Complex exp_h = exp(half);
        Complex sh;
        Complex ch;

        if (iszero(q))
        {
            ch = ONE;
            sh = ONE;
        }
        else
        {
            ch = cosh(q);
            sh = div(sinh(q), q);
        };

        Complex2x2 i_term = {ch, ZERO, ZERO, ch};
        Complex2x2 a_term = mat_scale(A_shifted, sh);
        Complex2x2 inner = mat_add(i_term, a_term);
        return mat_scale(inner, exp_h);
    };

    def mat_is_unitary(Complex2x2 A, double tol) -> bool
    {
        Complex2x2 prod = mat_mul(A, mat_adjoint(A));
        bool ok_a = isclose(prod.a, ONE, tol, tol);
        bool ok_b = isclose(prod.b, ZERO, tol, tol);
        bool ok_c = isclose(prod.c, ZERO, tol, tol);
        bool ok_d = isclose(prod.d, ONE, tol, tol);
        return ok_a & ok_b & ok_c & ok_d;
    };

    def mat_is_hermitian(Complex2x2 A, double tol) -> bool
    {
        Complex2x2 adj = mat_adjoint(A);
        bool ok_a = isclose(A.a, adj.a, tol, tol);
        bool ok_b = isclose(A.b, adj.b, tol, tol);
        bool ok_c = isclose(A.c, adj.c, tol, tol);
        bool ok_d = isclose(A.d, adj.d, tol, tol);
        return ok_a & ok_b & ok_c & ok_d;
    };

    def mat_eigenvalues(Complex2x2 A) -> Complex[2]
    {
        Complex tr = mat_trace(A);
        Complex det = mat_det(A);
        Complex disc = sqrt(sub(mul(tr, tr), scale(det, 4.0)));
        Complex[2] vals;
        vals[0] = scale(add(tr, disc), 0.5);
        vals[1] = scale(sub(tr, disc), 0.5);
        return vals;
    };
};
