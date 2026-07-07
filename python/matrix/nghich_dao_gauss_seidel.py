import numpy as np
from collections import deque


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def check_dominance(A):
    diag = np.abs(np.diag(A))
    row_sum = np.sum(np.abs(A), axis=1) - diag
    col_sum = np.sum(np.abs(A), axis=0) - diag

    row_dominant = np.all(diag > row_sum)
    col_dominant = np.all(diag > col_sum)
    return row_dominant, col_dominant


def forward_solve_step(A, B_col, x_prev):
    n = A.shape[0]
    x_next = np.zeros_like(x_prev)
    for i in range(n):
        s = B_col[i] - A[i, :i] @ x_next[:i] - A[i, i + 1 :] @ x_prev[i + 1 :]
        x_next[i] = s / A[i, i]
    return x_next


def gauss_seidel_iterate(A, B, eps, max_iter=1000):
    n = A.shape[0]
    diag = np.diag(A)

    if np.any(np.abs(diag) < 1e-9):
        print("Loi: co phan tu duong cheo a[i][i] = 0 -> khong the lap Gauss-Seidel!")
        return None, 0

    row_dominant, col_dominant = check_dominance(A)

    if row_dominant:
        print("=> Ma tran A cheo troi hang: dung chuan vo cung ||.||_inf.")
        norm_ord = np.inf
        s = 0.0
        off_upper = np.sum(np.abs(np.triu(A, 1)), axis=1)
        off_lower = np.sum(np.abs(np.tril(A, -1)), axis=1)
        q = np.max(off_upper / (np.abs(diag) - off_lower))
    elif col_dominant:
        print("=> Ma tran A cheo troi cot: dung chuan 1 ||.||_1.")
        norm_ord = 1
        off_lower_col = np.sum(np.abs(np.tril(A, -1)), axis=0)
        off_upper_col = np.sum(np.abs(np.triu(A, 1)), axis=0)
        s = np.max(off_lower_col / np.abs(diag))
        q = np.max(off_upper_col / (np.abs(diag) - off_lower_col))
    else:
        print("Canh bao: A khong cheo troi hang cung khong cheo troi cot -> khong dam bao hoi tu!")
        norm_ord = np.inf
        s = 0.0
        off_upper = np.sum(np.abs(np.triu(A, 1)), axis=1)
        off_lower = np.sum(np.abs(np.tril(A, -1)), axis=1)
        q = np.max(off_upper / (np.abs(diag) - off_lower))

    print(f"\ns = {s:.6f}, q = {q:.6f} -> " + ("KHONG dam bao hoi tu (q >= 1)!" if q >= 1 else "dam bao hoi tu (q < 1)."))

    D_A = np.diag(diag)
    L_A = -np.tril(A, -1)
    U_A = -np.triu(A, 1)

    print("\n--- MA TRAN D_A ---")
    print_matrix(D_A)
    print("--- MA TRAN L_A ---")
    print_matrix(L_A)
    print("--- MA TRAN U_A ---")
    print_matrix(U_A)

    p = B.shape[1]
    first_two = {}
    last_two = deque(maxlen=2)
    X_prev = np.zeros((n, p))
    for it in range(1, max_iter + 1):
        X_next = np.zeros((n, p))
        for c in range(p):
            X_next[:, c] = forward_solve_step(A, B[:, c], X_prev[:, c])

        diff = np.max(np.abs(X_next - X_prev)) if norm_ord == np.inf else np.sum(np.abs(X_next - X_prev))
        post_err = q / ((1 - s) * (1 - q)) * diff if q < 1 else float("inf")

        if it <= 2:
            first_two[it] = (X_next.copy(), post_err)
        last_two.append((it, X_next.copy(), post_err))

        if post_err < eps or it == max_iter:
            break

        X_prev = X_next

    last_its = {item[0] for item in last_two}
    for it in (1, 2):
        if it in first_two and it not in last_its:
            X_it, err_it = first_two[it]
            print(f"\n--- LAN LAP {it}: sai so hau nghiem uoc luong = {err_it:.6g} ---")
            print_matrix(X_it)

    for i, (it, X_it, err_it) in enumerate(last_two):
        label = "CUOI" if i == len(last_two) - 1 else "AP CHOT"
        print(f"\n--- LAN LAP {it} ({label}): sai so hau nghiem uoc luong = {err_it:.6g} ---")
        print_matrix(X_it)

    last_it, X_last, _ = last_two[-1]
    return X_last, last_it


def main():
    filename = "test.txt"
    try:
        A = np.loadtxt(filename, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename}'. Kiem tra lai duong dan!")
        return

    if A.size == 0:
        print("Loi: File rong hoac khong chua du lieu hop le!")
        return

    n, cols = A.shape
    if n != cols:
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the tim nghich dao bang PP lap Gauss-Seidel!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    E = np.eye(n)
    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-4): "))

    X, it = gauss_seidel_iterate(A, E, eps)
    if X is None:
        return

    print(f"\n=> Dung sau {it} lan lap.")
    print("--- MA TRAN NGHICH DAO XAP XI A^-1 ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
