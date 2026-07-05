import numpy as np
import sympy as sp

# ==== Sua truc tiep 2 muc nay theo de bai (so chieu tu suy ra tu do dai) ====
F_exprs_str = [
    "10*x1 - 2*x2**2 + x2 - 2*x3 - 5",
    "8*x2**2 + 4*x3**2 - 9",
    "8*x2*x3 + 4",
]
# X0=(0,0,0) lam J(X0) suy bien; X0=(0.1,0.1,-0.1) nhu newton.py lai qua xa nghiem
# khien Jacobi co dinh sai lech va Newton modified phan ky (mien hoi tu hep hon Newton
# thuong), nen chon diem gan nghiem hon de hoi tu
X0 = [0.6, 0.9, -0.4]
n_iterations = 5
# =============================================================================

n = len(F_exprs_str)
X_syms = sp.symbols(f"x1:{n + 1}")
_symbol_map = {str(s): s for s in X_syms}

F_exprs = [sp.sympify(e, locals=_symbol_map) for e in F_exprs_str]
F_matrix = sp.Matrix(F_exprs)
J_matrix = F_matrix.jacobian(X_syms)  # Jacobi tu suy ra tu F, khong go tay

F_func = sp.lambdify(X_syms, F_matrix, "numpy")
J_func = sp.lambdify(X_syms, J_matrix, "numpy")


def main():
    X = np.array(X0, dtype=float)

    # Newton modified: chi tinh va nghich dao J tai X0 mot lan, dung lai cho moi buoc
    J0 = np.array(J_func(*X), dtype=float)
    try:
        J0_inv = np.linalg.inv(J0)
    except np.linalg.LinAlgError:
        print("Loi: Ma tran Jacobi tai X0 bi suy bien.")
        return

    header = f"{'k':>3}" + "".join(f"{'x' + str(i + 1) + ',k':>15}" for i in range(n))
    border = "-" * len(header)

    print(header)
    print(border)
    print(f"{0:>3}" + "".join(f"{v:>15.6f}" for v in X))

    for k in range(1, n_iterations + 1):
        Fx = np.array(F_func(*X), dtype=float).flatten()
        X = X - J0_inv @ Fx
        print(f"{k:>3}" + "".join(f"{v:>15.6f}" for v in X))

    print(border)
    print(f"Xap xi nghiem X_{n_iterations}:")
    for i, val in enumerate(X):
        print(f"  x{i + 1} = {val:.6f}")


if __name__ == "__main__":
    main()
