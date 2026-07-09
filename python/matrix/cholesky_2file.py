import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def cholesky_decompose(M):
    n = M.shape[0]
    U = np.zeros((n, n))

    show = sorted(set(range(min(2, n))) | set(range(max(0, n - 2), n)))

    for i in range(n):
        s = M[i, i] - np.sum(U[:i, i] ** 2)

        if s <= 1e-9:
            print(f"Loi: u[{i + 1}][{i + 1}]^2 = {s:.4f} <= 0 -> Ma tran khong thoa dieu kien cua PP Choleski!")
            return None

        U[i, i] = np.sqrt(s)

        for k in range(i + 1, n):
            U[i, k] = (M[i, k] - np.sum(U[:i, i] * U[:i, k])) / U[i, i]

        if i in show:
            print(f"--- Lan lap i = {i + 1} ---")
            print(f"U[{i + 1},:] =")
            print_matrix(U[i, :].reshape(1, -1))

    return U


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
        x[i] = (y[i] - U[i, i + 1 :] @ x[i + 1 :]) / U[i, i]
        print(f"x[{i + 1}] = {x[i]:.4f}")
    return x


def main():
    filename_a = "testA.txt"
    filename_b = "testB.txt"
    try:
        A = np.loadtxt(filename_a, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename_a}'. Kiem tra lai duong dan!")
        return

    try:
        B = np.loadtxt(filename_b, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename_b}'. Kiem tra lai duong dan!")
        return

    if A.size == 0 or B.size == 0:
        print("Loi: File rong hoac khong chua du lieu hop le!")
        return

    m, cols_a = A.shape
    m_b, cols_b = B.shape

    if m_b != m:
        print(f"Loi: Ma tran A ({m} hang) va ma tran B ({m_b} hang) khong khop so hang!")
        return

    if cols_a != m:
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {m} x {cols_a}), khong the ap dung PP Choleski!")
        return

    print(f"\n--- MA TRAN A ({m} x {cols_a}) ---")
    print_matrix(A)
    print(f"--- MA TRAN B ({m} x {cols_b}) ---")
    print_matrix(B)

    if np.allclose(A, A.T, atol=1e-9):
        print("=> Ma tran A doi xung: ap dung truc tiep PP Choleski cho A.\n")
        M, D = A, B
    else:
        print("=> Ma tran A khong doi xung: chuyen ve giai M x = d voi M = A^T.A, d = A^T.b\n")
        M, D = A.T @ A, A.T @ B

    U = cholesky_decompose(M)
    if U is None:
        return

    print("--- MA TRAN U (A = U^T.U) ---")
    print_matrix(U)

    X = np.zeros((cols_a, cols_b))
    for c in range(cols_b):
        print(f"\n=== Giai cot thu {c + 1} cua B ===")
        print("The xuoi (U^T.y = d):")
        y = forward_substitution(U.T, D[:, c])
        print("The nguoc (U.x = y):")
        X[:, c] = backward_substitution(U, y)

    print("\n--- MA TRAN NGHIEM X ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
