module math_algebra

fn pow(base: number, exp: number) -> number {
    result := 1;
    i := 0;
    while (i < exp) {
        result = result * base;
        i = i + 1;
    }
    return result;
}

fn sqrt(x: number) -> number {
    guess := x / 2;
    i := 0;

    while (i < 20) {
        guess = (guess + x / guess) / 2;
        i = i + 1;
    }

    return guess;
}

fn exp(x: number) -> number {
    return 1 + x + (x*x)/2 + (x*x*x)/6;
}

fn log(x: number) -> number {
    // basic approximation
    return (x - 1) - ((x - 1)*(x - 1))/2;
}