# PP chia doi theo cong thuc sai so hau nghiem

import math

# NHAP HAM f(x) O DAY
def f(x):
    return math.log(x) - 1


def bisection(a, b, epsilon):
    steps = 0
    c = a
    border = "-" * 63

    print(f"{'k':>3}{'a_k':>15}{'b_k':>15}{'x_k+1':>15}{'Dau f(x_k+1)':>15}")
    print(border)

    while b - a >= epsilon:
        c = (a + b) / 2
        sign = "+" if f(c) > 0 else "-"
        print(f"{steps:>3}{a:>15.8f}{b:>15.8f}{c:>15.8f}{sign:^15}")
        if f(c) == 0:
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
        steps += 1

    print(border)
    print(f"Nghiem gan dung: {c:.10g}")


def main():
    a, b, epsilon = map(float, input("Nhap a, b, epsilon: ").split())
    bisection(a, b, epsilon)


if __name__ == "__main__":
    main()
