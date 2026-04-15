module math_number

fn gcd(a: number, b: number) -> number {
    while (b != 0) {
        temp := b;
        b = a % b;
        a = temp;
    }
    return a;
}

fn lcm(a: number, b: number) -> number {
    return (a * b) / gcd(a, b);
}

fn factorial(n: number) -> number {
    result := 1;
    i := 2;
    while (i <= n) {
        result = result * i;
        i = i + 1;
    }
    return result;
}