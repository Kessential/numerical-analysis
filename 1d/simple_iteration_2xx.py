# Phuong phap lap don
# theo 2 xap xi lien tiep

import sympy as sp

_x = sp.symbols("x")

# NHAP HAM phi(x) O DAY
_phi_expr = 1 / sp.sqrt(_x + 3)  # <- ham phi(x), bien x go thanh _x

_dphi_expr = sp.diff(_phi_expr, _x)

phi = sp.lambdify(_x, _phi_expr, "math")
dphi = sp.lambdify(_x, _dphi_expr, "math")


def estimate_q(a, b, samples=2000):
    """Uoc luong q = sup|phi'(x)| tren [a,b] bang lay mau."""
    return max(abs(dphi(a + i * (b - a) / (samples - 1))) for i in range(samples))


def check_convergence(a, b, samples=2000):
    q = estimate_q(a, b, samples)
    vals = [phi(a + i * (b - a) / (samples - 1)) for i in range(samples)]
    phi_min, phi_max = min(vals), max(vals)
    maps_into_ab = phi_min >= a - 1e-9 and phi_max <= b + 1e-9
    return q, maps_into_ab, phi_min, phi_max


def solve(a, b, epsilon, max_iter=1000):
    q, maps_into_ab, phi_min, phi_max = check_convergence(a, b)

    print(f"He so co uoc luong q = sup|phi'(x)| tren [a,b] = {q:.10f}")
    if not maps_into_ab:
        print(f"Canh bao: phi([a,b]) = [{phi_min:.6f}, {phi_max:.6f}] khong nam trong [a,b]!")
    if q >= 1:
        print("Canh bao: q >= 1, dieu kien anh xa co khong thoa man, dung lai.")
        return

    x_0 = (a + b) / 2
    print(f"Chon x0 = {x_0:.10f} (trung diem [a,b])\n")

    epsilon_0 = (1 - q) * epsilon / q
    x_old = x_0
    n = 0

    print(f"{'k':>3}{'x_k':>15}{'x_k+1':>15}{'Sai so':>15}")
    print("-" * 49)

    while True:
        x_new = phi(x_old)
        error = abs(x_old - x_new)

        print(f"{n:>3}{x_old:>15.10f}{x_new:>15.10f}{error:>15.10f}")

        if error < epsilon_0:
            print("-" * 49)
            print(f"Nghiem gan dung cua PT la: {x_new:.10f}")
            break

        x_old = x_new
        n += 1

        if n > max_iter:
            print(f"Thuat toan ko hoi tu sau {max_iter} lan lap")
            break


def main():
    default_max_iter = 1000
    data = input("Nhap a, b, epsilon va (tuy chon) so lan lap toi da: ").split()
    a, b, epsilon = map(float, data[:3])

    if len(data) > 3:
        max_iter = int(float(data[3]))
    else:
        try:
            extra = input(f"So lan lap toi da (Enter = {default_max_iter}): ").strip()
        except EOFError:
            extra = ""
        max_iter = int(extra) if extra else default_max_iter

    solve(a, b, epsilon, max_iter)


if __name__ == "__main__":
    main()
