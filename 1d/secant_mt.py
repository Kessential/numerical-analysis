# Phuong phap day cung
# dua theo cong thuc sai so muc tieu
# | x_k - x* | <= | f(x_k) | / m1

import sympy as sp

_x = sp.symbols("x")

# CHINH HAM f(x) O DAY:
_f_expr = _x**5 - 17 # <- Ham f(x), bien x go thanh _x

_df_expr = sp.diff(_f_expr, _x)
_d2f_expr = sp.diff(_df_expr, _x)

f = sp.lambdify(_x, _f_expr, "math")
df = sp.lambdify(_x, _df_expr, "math")
d2f = sp.lambdify(_x, _d2f_expr, "math")


def solve(a, b, epsilon, max_iter=100):
    if f(a) * d2f(a) > 0:
        d, x_old = a, b
    else:
        d, x_old = b, a

    m1 = min(abs(df(a)), abs(df(b)))
    k = 0

    print(f"Diem Fourier d = {d:.10f}")
    print(f"Diem bat dau x0 = {x_old:.10f}")
    print("-" * 52)
    print(f"{'k':>3}{'x_k':>15}{'x_k+1':>15}{'Sai so':>15}")

    while True:
        x_new = x_old - (f(x_old) * (x_old - d)) / (f(x_old) - f(d))
        error = abs(f(x_new)) / m1

        print(f"{k:>3}{x_old:>15.10f}{x_new:>15.10f}{error:>15.10f}")

        if error <= epsilon:
            print("-" * 52)
            print(f"Nghiem xap xi tim duoc: {x_new:.10f}")
            break

        x_old = x_new
        k += 1

        if k > max_iter:
            print(f"Lap qua {max_iter} lan!")
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
