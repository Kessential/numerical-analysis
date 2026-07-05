# Phuong phap tiep tuyen (Newton - Raphson)
# theo cong thuc sai so 2 lan xap xi lien tiep
# | x_k - x* | <= ( M2 / 2.m1 ) * | x_k - x_k-1 |^2

import sympy as sp

_x = sp.symbols("x")

# NHAP HAM f(_x) o day
_f_expr = sp.exp(_x) - 4 * sp.cos(_x / 2) + 1 # <- Ham (f_x) voi bien _x

_df_expr = sp.diff(_f_expr, _x)
_d2f_expr = sp.diff(_df_expr, _x)

f = sp.lambdify(_x, _f_expr, "math")
df = sp.lambdify(_x, _df_expr, "math")
d2f = sp.lambdify(_x, _d2f_expr, "math")


def solve(a, b, epsilon, max_iter=100):
    x_old = a if f(a) * d2f(a) > 0 else b

    M2 = max(abs(d2f(a)), abs(d2f(b)))
    m1 = min(abs(df(a)), abs(df(b)))

    print(f"Diem bat dau x0 = {x_old:.8f}")
    print("-" * 52)
    print(f"{'k':>3}{'x_k':>15}{'x_k+1':>15}{'Sai so':>15}")

    k = 0
    while True:
        x_new = x_old - f(x_old) / df(x_old)
        diff = abs(x_new - x_old)
        error = (M2 / (2 * m1)) * diff * diff

        print(f"{k:>3}{x_old:>15.8f}{x_new:>15.8f}{error:>15.8f}")

        if error <= epsilon:
            print("-" * 52)
            print(f"Nghiem xap xi tim duoc: {x_new:.8f}")
            break

        x_old = x_new
        k += 1

        if k > max_iter:
            print(f"Thuat toan chay hon {max_iter} lan lap!")
            break


def main():
    default_max_iter = 100
    data = input("Nhap a, b, epsilon va (tuy chon) so lan lap toi da: ").split()
    a, b, epsilon = map(float, data[:3])

    if len(data) > 3:
        max_iter = int(data[3])
    else:
        try:
            extra = input(f"So lan lap toi da (Enter = {default_max_iter}): ").strip()
        except EOFError:
            extra = ""
        max_iter = int(extra) if extra else default_max_iter

    solve(a, b, epsilon, max_iter)


if __name__ == "__main__":
    main()
