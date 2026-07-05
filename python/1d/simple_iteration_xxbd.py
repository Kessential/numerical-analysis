# Phuong phap lap don
# theo xap xi ban dau

import math
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


def solve(a, b, epsilon):
    q, maps_into_ab, phi_min, phi_max = check_convergence(a, b)

    print(f"He so co uoc luong q = sup|phi'(x)| tren [a,b] = {q:.10f}")
    if not maps_into_ab:
        print(f"Canh bao: phi([a,b]) = [{phi_min:.6f}, {phi_max:.6f}] khong nam trong [a,b]!")
    if q >= 1:
        print("Canh bao: q >= 1, dieu kien anh xa co khong thoa man, dung lai.")
        return

    x_0 = (a + b) / 2
    print(f"Chon x0 = {x_0:.10f} (trung diem [a,b])")

    x_1 = phi(x_0)
    d = abs(x_1 - x_0)

    tu = math.log((epsilon * (1 - q)) / d)
    mau = math.log(q)

    n_target = math.ceil(tu / mau)

    print(f"|x1 - x0| = {d:.8f}")
    print(f"So buoc lap du kien de dat sai so: {n_target}")
    border = "-" * 54
    print(border)
    print(f"{'k':>3}{'x_k':>15}{'x_k+1':>15}{'|x_k+1 - x_k|':>18}")
    print(border)

    # Chay vong lap dung n_target lan
    x_current = x_0
    for k in range(0, n_target):
        x_next = phi(x_current)
        diff = abs(x_next - x_current)
        print(f"{k:>3}{x_current:>15.8f}{x_next:>15.8f}{diff:>18.8f}")
        x_current = x_next

    print(border)
    print(f"Nghiem tim duoc sau {n_target} buoc: {x_current:.8f}")


def main():
    a, b, epsilon = map(float, input("Nhap a, b, epsilon: ").split())
    solve(a, b, epsilon)


if __name__ == "__main__":
    main()
