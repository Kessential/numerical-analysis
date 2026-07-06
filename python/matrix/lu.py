import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def lu_decompose(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.eye(n)

    show = sorted(set(range(min(2, n))) | set(range(max(0, n - 2), n)))

    for t in range(n):
        for i in range(t, n):
            L[i, t] = A[i, t] - L[i, :t] @ U[:t, t]

        if abs(L[t, t]) < 1e-9:
            print(f"Loi: l[{t + 1}][{t + 1}] = 0 -> Ma tran A suy bien, khong the phan tach LU!")
            return None, None

        for k in range(t + 1, n):
            U[t, k] = (A[t, k] - L[t, :t] @ U[:t, k]) / L[t, t]

        if t in show:
            print(f"--- Lan lap t = {t + 1} ---")
            print(f"L[:,{t + 1}] =")
            print_matrix(L[:, t].reshape(-1, 1))
            print(f"U[{t + 1},:] =")
            print_matrix(U[t, :].reshape(1, -1))

    return L, U


def forward_substitution(L, b):
    n = L.shape[0]
    y = np.zeros_like(b, dtype=float)
    for i in range(n):
        y[i] = (b[i] - L[i, :i] @ y[:i]) / L[i, i]
        print(f"y[{i + 1}] = {y[i]:.4f}")
    return y


def backward_substitution(U, y):
    n = U.shape[0]
    x = np.zeros_like(y, dtype=float)
    for i in range(n - 1, -1, -1):
        x[i] = y[i] - U[i, i + 1 :] @ x[i + 1 :]
        print(f"x[{i + 1}] = {x[i]:.4f}")
    return x


def main():
    filename = "test.txt"
    try:
        Aug = np.loadtxt(filename, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename}'. Kiem tra lai duong dan!")
        return

    if Aug.size == 0:
        print("Loi: File rong hoac khong chua du lieu hop le!")
        return

    m, total_cols = Aug.shape
    print(f"Ma tran mo rong co tong cong {total_cols} cot.")
    cols_b = int(input("Nhap so cot cua ma tran B (vi du p = 1, 2...): "))
    cols_a = total_cols - cols_b

    if cols_a != m:
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {m} x {cols_a}), khong the phan tach LU!")
        return

    A = Aug[:, :cols_a]
    B = Aug[:, cols_a:]

    print(f"\n--- MA TRAN A ({m} x {cols_a}) ---")
    print_matrix(A)
    print(f"--- MA TRAN B ({m} x {cols_b}) ---")
    print_matrix(B)

    L, U = lu_decompose(A)
    if L is None:
        return

    print("--- MA TRAN L ---")
    print_matrix(L)
    print("--- MA TRAN U ---")
    print_matrix(U)

    X = np.zeros((cols_a, cols_b))
    for c in range(cols_b):
        print(f"\n=== Giai cot thu {c + 1} cua B ===")
        print("The xuoi (Ly = b):")
        y = forward_substitution(L, B[:, c])
        print("The nguoc (Ux = y):")
        X[:, c] = backward_substitution(U, y)

    print("--- MA TRAN NGHIEM X ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
