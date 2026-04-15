module math_core

fn abs(x: number) -> number {
    if (x < 0) return -x;
    return x;
}

fn floor(x: number) -> number {
    i := int(x);
    if (i > x) return i - 1;
    return i;
}

fn ceil(x: number) -> number {
    i := int(x);
    if (i < x) return i + 1;
    return i;
}

fn trunc(x: number) -> number {
    return int(x);
}