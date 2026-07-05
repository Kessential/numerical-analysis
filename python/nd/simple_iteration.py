import numpy as np
import sympy as sp

# ==== Sua truc tiep cac muc nay theo de bai (so chieu tu suy ra tu do dai Phi_exprs_str) ====
Phi_exprs_str = [
    "(cos(x2*x3) + 0.5) / 3",
    "(1/25) * sqrt(x1**2 + 0.3125) - 0.03",
    "-(1/20) * exp(-x1*x2) - (10*pi - 3)/60",
]
X0 = [0.0, 0.0, 0.0]
epsilon = 1e-6
max_iter = 100
D = [(-1, 1), (-1, 1), (-1, 1)]  # mien de uoc luong he so co q tu Jacobi cua Phi
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

    q = min(q_row, q_col, K)
    if q >= 1:
        print("Canh bao: khong cach nao cho q<1 => chua chac hoi tu tren D da chon!")
        return
    print(f"=> Chon q = {q:.6f} (< 1) de danh gia sai so\n")

    tol = epsilon * (1 - q) / q

    X = np.array(X0, dtype=float)
    header = f"{'k':>3}" + "".join(f"{'x' + str(i + 1) + ',k':>15}" for i in range(n)) + f"{'Sai so':>15}"
    border = "-" * len(header)

    print(border)
    print(header)
    print(border)

    for k in range(1, max_iter + 1):
        X_new = np.array(Phi_func(*X), dtype=float).flatten()
        diff = np.max(np.abs(X_new - X))

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
