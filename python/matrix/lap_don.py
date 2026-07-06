import numpy as np
from collections import deque

def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def simple_iteration(C, D, x0, eps, max_iter=1000):
    q_inf = np.linalg.norm(C, np.inf)
    q_1 = np.linalg.norm(C, 1)
    print(f"\n||C||_inf = {q_inf:.6f}, ||C||_1 = {q_1:.6f}")

    if q_inf < 1:
        norm_ord, q = np.inf, q_inf
        print(f"=> Dung chuan vo cung: q = {q:.6f} < 1 -> dam bao hoi tu.")
    elif q_1 < 1:
        norm_ord, q = 1, q_1
        print(f"=> Chuan vo cung khong hoi tu (q >= 1), chuyen sang chuan 1: q = {q:.6f} < 1 -> dam bao hoi tu.")
    else:
        norm_ord, q = np.inf, q_inf
        print(f"Canh bao: ca hai chuan deu cho q >= 1 -> khong dam bao hoi tu! (dung tam chuan vo cung, q = {q:.6f})")

    first_two = {}
    last_two = deque(maxlen=2)
    x_prev = x0.copy()
    for n in range(1, max_iter + 1):
        x_next = C @ x_prev + D
        diff = np.max(np.abs(x_next - x_prev)) if norm_ord == np.inf else np.sum(np.abs(x_next - x_prev))
        post_err = q / (1 - q) * diff if q < 1 else float("inf")

        if n <= 2:
            first_two[n] = (x_next.copy(), post_err)
        last_two.append((n, x_next.copy(), post_err))

        if post_err < eps or n == max_iter:
            break

        x_prev = x_next

    last_ns = {item[0] for item in last_two}
    for n in (1, 2):
        if n in first_two and n not in last_ns:
            X_n, err_n = first_two[n]
            print(f"\n--- LAN LAP {n}: sai so hau nghiem uoc luong = {err_n:.6g} ---")
            print_matrix(X_n)

    for i, (n, X_n, err_n) in enumerate(last_two):
        label = "CUOI" if i == len(last_two) - 1 else "AP CHOT"
        print(f"\n--- LAN LAP {n} ({label}): sai so hau nghiem uoc luong = {err_n:.6g} ---")
        print_matrix(X_n)

    last_n, X_last, _ = last_two[-1]
    return X_last, last_n


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
    print("Quy uoc: he can giai da duoc dua ve dang x = Cx + D, file la ma tran mo rong [C | D].")
    cols_d = int(input("Nhap so cot cua ma tran D (vi du p = 1, 2...): "))
    cols_c = total_cols - cols_d

    if cols_c != m:
        print(f"Loi: Ma tran C phai vuong (dang doc duoc la {m} x {cols_c}), khong the lap!")
        return

    C = Aug[:, :cols_c]
    D = Aug[:, cols_c:]

    print(f"\n--- MA TRAN C ({m} x {cols_c}) ---")
    print_matrix(C)
    print(f"--- MA TRAN D ({m} x {cols_d}) ---")
    print_matrix(D)

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-4): "))

    X, n = simple_iteration(C, D, D.copy(), eps)

    print(f"\n=> Dung sau {n} lan lap.")


if __name__ == "__main__":
    main()
