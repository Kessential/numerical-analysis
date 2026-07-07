import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the tim nghich dao!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    # Ghep [A | E]
    Aug = np.hstack([A, np.eye(n)])

    print("--- LICH SU CHON PHAN TU KHOA (PIVOT) ---")

    # Quy trinh thuan: khu Gauss, doi hang khi pivot = 0
    for i in range(n):
        if abs(Aug[i, i]) < 1e-9:
            found = False
            for t in range(i + 1, n):
                if abs(Aug[t, i]) > 1e-9:
                    print(f"[!] Phat hien a[{i + 1}][{i + 1}] = 0. Doi cho hang {i + 1} va hang {t + 1}")
                    Aug[[i, t]] = Aug[[t, i]]
                    found = True
                    break
            if not found:
                print(f"Loi: khong tim duoc pivot khac 0 tai cot {i + 1} -> A suy bien, khong kha nghich!")
                return

        print(f"Lan lap {i + 1}: Chon pivot a[{i + 1}][{i + 1}] = {Aug[i, i]:.4f}")

        if i < n - 1:
            factors = Aug[i + 1 :, i] / Aug[i, i]
            Aug[i + 1 :, i:] -= np.outer(factors, Aug[i, i:])

        if i == 0:
            print("\n--- MA TRAN SAU LAN LAP 1 ---")
            print_matrix(Aug)

    print("\n--- MA TRAN SAU QUY TRINH THUAN [U | ...] ---")
    print_matrix(Aug)

    # Quy trinh nghich: the nguoc tung cot cua E de duoc tung cot cua A^-1
    X = np.zeros((n, n))
    for r in range(n - 1, -1, -1):
        s = Aug[r, r + 1 : n] @ X[r + 1 : n, :]
        X[r, :] = (Aug[r, n:] - s) / Aug[r, r]

    print("--- MA TRAN NGHICH DAO A^-1 ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
