import math

import sympy as sp

_x = sp.symbols("x")

# NHAP HAM f(_x) o day
_f_expr = _x**3 - 2 * _x - 5  # <- Ham f(_x) voi bien _x

_df_expr = sp.diff(_f_expr, _x)
_d2f_expr = sp.diff(_df_expr, _x)

f = sp.lambdify(_x, _f_expr, "math")
df = sp.lambdify(_x, _df_expr, "math")
d2f = sp.lambdify(_x, _d2f_expr, "math")


def combined(a, b, epsilon=1e-5, max_iter=100):
    # Buoc 1: Kiem tra dieu kien co nghiem
    if f(a) * f(b) >= 0:
        print("Loi: [a, b] khong chua nghiem")
        return math.nan

    # Buoc 2: Chon diem xuat phat theo dieu kien Fourier f(x) * f''(x) > 0
    if f(a) * d2f(a) > 0:
        x_tt, x_dc = a, b
    elif f(b) * d2f(b) > 0:
        x_tt, x_dc = b, a
    else:
        print("Loi: f(x)*f''(x) > 0 tai hai dau mut")
        return math.nan

    x_dc_new = x_dc

    border = "-" * 49
    print(border)
    print(f"{'k':>3}{'x_{2k+2}':>15}{'x_{2k+3}':>15}{'Sai so':>15}")
    print(border)

    # Buoc 3: Vong lap
    for i in range(max_iter):
        # 1. Tinh diem Tiep tuyen moi
        df_x_tt = df(x_tt)
        if abs(df_x_tt) < 1e-9:
            print("Loi: Dao ham bang 0, dung thuat toan.")
            return math.nan
        x_tt_new = x_tt - f(x_tt) / df_x_tt

        # 2. Tinh diem Day cung moi
        mau_so = f(x_dc) - f(x_tt)
        if abs(mau_so) < 1e-9:
            print("Loi: Mau so day cung bang 0, dung thuat toan.")
            return math.nan
        tu_so = f(x_dc) * (x_dc - x_tt)
        x_dc_new = x_dc - tu_so / mau_so

        # 3. Danh gia sai so
        sai_so = abs(x_dc_new - x_dc)

        print(f"{i:>3}{x_tt_new:>15.8f}{x_dc_new:>15.8f}{sai_so:>15.8f}")

        if sai_so <= epsilon:
            print(border)
            print(f"=> Hoi tu sau {i + 1} buoc lap.")
            return x_dc_new

        x_tt = x_tt_new
        x_dc = x_dc_new

    print("Canh bao: Dat den so vong lap toi da ma chua hoi tu.")
    return x_dc_new


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

    x = combined(a, b, epsilon, max_iter)

    if not math.isnan(x):
        print(f"\nNghiem gan dung cuoi cung: {x:.8f}")


if __name__ == "__main__":
    main()
