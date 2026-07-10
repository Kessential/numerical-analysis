import numpy as np
from collections import deque


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def newton_inverse(A, eps, max_iter=1000):
    n = A.shape[0]
    E = np.eye(n)

    norm1 = np.linalg.norm(A, 1)
    norm_inf = np.linalg.norm(A, np.inf)
    X = A.T / (norm1 * norm_inf)

    print(f"||A||_1 = {norm1:.6f}, ||A||_inf = {norm_inf:.6f}")
    print("--- MA TRAN BAN DAU X_0 = A^T / (||A||_1 * ||A||_inf) ---")
    print_matrix(X)

    G = E - A @ X
    # Dieu kien hoi tu ly thuyet ||G_0||_2 < 1 dung chuan pho (spectral norm), nhung
    # o day dung ||.||_inf de nhat quan voi cac PP lap khac trong repo va de tinh
    # (khong can SVD) - day chi la uoc luong thuc dung, KHONG phai dieu kien du/can
    # chat che: co truong hop ||G_0||_inf >= 1 nhung ||G_0||_2 < 1 (van hoi tu), hoac
    # nguoc lai. Neu nghi ngo, doi chieu ket qua cuoi voi cac PP khac (Gauss-Jordan...).
    q = np.linalg.norm(G, np.inf)
    print(f"G_0 = E - A.X_0, ||G_0||_inf = q = {q:.6f} -> " + ("CANH BAO (q >= 1, chi la uoc luong theo chuan vo cung, khong hoan toan chac chan)!" if q >= 1 else "uoc luong hoi tu (q < 1)."))

    first_two = {}
    last_two = deque(maxlen=2)

    for it in range(1, max_iter + 1):
        X = X @ (2 * E - A @ X)
        G = E - A @ X
        err = np.linalg.norm(G, np.inf)

        if it <= 2:
            first_two[it] = (X.copy(), err)
        last_two.append((it, X.copy(), err))

        if err < eps or it == max_iter:
            break

    last_its = {item[0] for item in last_two}
    for it in (1, 2):
        if it in first_two and it not in last_its:
            X_it, err_it = first_two[it]
            print(f"\n--- LAN LAP {it}: ||E - A.X_n||_inf = {err_it:.6g} ---")
            print_matrix(X_it)

    for i, (it, X_it, err_it) in enumerate(last_two):
        label = "CUOI" if i == len(last_two) - 1 else "AP CHOT"
        print(f"\n--- LAN LAP {it} ({label}): ||E - A.X_n||_inf = {err_it:.6g} ---")
        print_matrix(X_it)

    last_it, X_last, last_err = last_two[-1]
    if last_err >= eps:
        print(
            f"\nCANH BAO: sau {last_it} lan lap (toi da) van chua dat sai so yeu cau "
            f"(||E-A.X_n||_inf = {last_err:.6g} >= eps = {eps:.6g}). Ket qua co the khong dang tin cay -> "
            "A gan suy bien / dieu kien xau, hoac can tang max_iter."
        )

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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the ap dung PP lap Newton!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-4): "))

    X, it = newton_inverse(A, eps)

    print(f"\n=> Dung sau {it} lan lap.")
    print("--- MA TRAN NGHICH DAO XAP XI A^-1 ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
