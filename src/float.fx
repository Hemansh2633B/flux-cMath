module math_float

fn isfinite(x: number) -> bool {
    return x != INF && x != -INF;
}

fn isnan(x: number) -> bool {
    return x != x;
}

fn isinf(x: number) -> bool {
    return x == INF || x == -INF;
}

fn isclose(a: number, b: number, tol: number) -> bool {
    diff := a - b;
    if (diff < 0) diff = -diff;
    return diff < tol;
}