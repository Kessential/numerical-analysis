import numpy as np
import sympy as sp

# ==== Sua truc tiep cac muc nay theo de bai (so chieu tu suy ra tu do dai Phi_exprs_str) ====
Phi_exprs_str = [
    "1 - cos(x1*x2*x3)",
    "1 - (1 - x1)**(1/4) - 0.05*x3**2 + 0.15*x3",
    "x1**2 + 0.1*x2**2 - 0.01*x2 + 1",
]
X0 = [0.0, 0.1, 0.3]
epsilon = 1e-6
max_iter = 100
D = [(-0.1, 0.1), (-0.1, 0.3), (0.5, 1.1)]  # mien de uoc luong he so co q tu Jacobi cua Phi
# ==============================================================================================

n = len(Phi_exprs_str)
X_syms = sp.symbols(f"x1:{n + 1}")
_symbol_map = {str(s): s for s in X_syms}

Phi_exprs = [sp.sympify(e, locals=_symbol_map) for e in Phi_exprs_str]
Phi_matrix = sp.Matrix(Phi_exprs)
J_matrix = Phi_matrix.jacobian(X_syms)  # Jacobi cua Phi tu suy ra, khong go tay

Phi_func = sp.lambdify(X_syms, Phi_matrix, "numpy")
J_func = sp.lambdify(X_syms, J_matrix, "numpy")


def estimate_q(samples_per_dim=25):
    """Uoc luong q_row, q_col, K bang cach lay sup |dPhi_i/dx_j| tren luoi mau trong D."""
    grids = [np.linspace(a, b, samples_per_dim) for a, b in D]
    mesh = np.meshgrid(*grids, indexing="ij")
    points = np.stack([m.ravel() for m in mesh], axis=-1)

    sup_abs_J = np.zeros((n, n))
    for X in points:
        Jx = np.abs(np.array(J_func(*X), dtype=float))
        sup_abs_J = np.maximum(sup_abs_J, Jx)

    q_row = np.max(np.sum(sup_abs_J, axis=1))
    q_col = np.max(np.sum(sup_abs_J, axis=0))
    K = np.max(sup_abs_J) * n

    return q_row, q_col, K


def main():
    q_row, q_col, K = estimate_q()
    print(f"He so co uoc luong tren D: q_row={q_row:.6f}, q_col={q_col:.6f}, K={K:.6f}")

    # Chuan do sai so phai khop voi chuan sinh ra q: q_row <-> chuan vo cung (max),
    # q_col <-> chuan 1 (tong). K chan tren ca hai nen dung chuan vo cung cho an toan.
    candidates = [
        ("q_row", q_row, "chuan vo cung ||.||_inf (max tri tuyet doi)", lambda d: np.max(np.abs(d))),
        ("q_col", q_col, "chuan 1 ||.||_1 (tong tri tuyet doi)", lambda d: np.sum(np.abs(d))),
        ("K", K, "chuan vo cung ||.||_inf (max tri tuyet doi)", lambda d: np.max(np.abs(d))),
    ]

    valid = [c for c in candidates if c[1] < 1]
    if not valid:
        print("Canh bao: khong cach nao cho q<1 => chua chac hoi tu tren D da chon!")
        return

    name, q, norm_label, norm_func = min(valid, key=lambda c: c[1])
    print(f"=> Chon q = {name} = {q:.6f} (< 1), do sai so bang {norm_label}\n")

    tol = epsilon * (1 - q) / q

    X = np.array(X0, dtype=float)
    header = f"{'k':>3}" + "".join(f"{'x' + str(i + 1) + ',k':>15}" for i in range(n)) + f"{'Sai so':>15}"
    border = "-" * len(header)

    print(border)
    print(header)
    print(border)

    for k in range(1, max_iter + 1):
        X_new = np.array(Phi_func(*X), dtype=float).flatten()
        diff = norm_func(X_new - X)

        print(f"{k:>3}" + "".join(f"{v:>15.8f}" for v in X_new) + f"{diff:>15.8f}")

        if diff <= tol:
            print(border)
            print(f"[+] Hoi tu sau {k} buoc lap.")
            actual_error = (q / (1 - q)) * diff
            print(f"[-] Sai so hau nghiem thuc te: {actual_error:.8f} <= {epsilon}")
            print("[-] Nghiem gan dung la:")
            for i, v in enumerate(X_new):
                print(f"    x{i + 1} = {v:.8f}")
            return

        X = X_new

    print(border)
    print(f"[!] Khong hoi tu sau {max_iter} buoc lap.")


if __name__ == "__main__":
    main()
