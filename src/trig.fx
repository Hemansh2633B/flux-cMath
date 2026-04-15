module math_trig

fn sin(x: number) -> number {
    return x - (x*x*x)/6 + (x*x*x*x*x)/120;
}

fn cos(x: number) -> number {
    return 1 - (x*x)/2 + (x*x*x*x)/24;
}

fn tan(x: number) -> number {
    return sin(x) / cos(x);
}

fn atan(x: number) -> number {
    return x - (x*x*x)/3 + (x*x*x*x*x)/5;
}