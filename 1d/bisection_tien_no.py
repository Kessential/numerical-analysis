# PP chia doi theo cong thuc sai so tien nghiem

import math

# NHAP HAM f(x) O DAY
def f(x):
    return math.log(x) - 1


def bisection(a, b, epsilon):
    # 1. Tinh truoc so lan lap toi thieu (tien nghiem)
    n = math.ceil(math.log2((b - a) / epsilon))

    print(f"So buoc lap tinh truoc (n) = {n}")

    c = a
    border = "-" * 63

    print(f"{'k':>3}{'a_k':>15}{'b_k':>15}{'x_k+1':>15}{'Dau f(x_k+1)':>15}")
    print(border)

    # 2. Vong lap chay chinh xac n lan
    for k in range(0, n):
        c = (a + b) / 2
        sign = "+" if f(c) > 0 else "-"
        print(f"{k:>3}{a:>15.8f}{b:>15.8f}{c:>15.8f}{sign:^15}")

        if f(c) == 0:
            break
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c

    print(border)
    print(f"Nghiem gan dung sau {n} buoc: {c:.10g}")


def main():
    a, b, epsilon = map(float, input("Nhap a, b, epsilon: ").split())
    bisection(a, b, epsilon)


if __name__ == "__main__":
    main()
