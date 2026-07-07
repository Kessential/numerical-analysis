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

    ind = [-1] * n
    row_used = [False] * n
    col_used = [False] * n
    step = 1

    print("--- LICH SU CHON PHAN TU KHOA (PIVOT) ---")

    for _ in range(n):
        max_val = 0.0
        p, q = -1, -1
        found_priority_1 = False

        # Uu tien 1: pivot co |gia tri| = 1
        for r in range(n):
            if row_used[r]:
                continue
            for c in range(n):
                if not col_used[c] and abs(abs(Aug[r, c]) - 1) < 1e-9:
                    p, q = r, c
                    max_val = abs(Aug[r, c])
                    found_priority_1 = True
                    break
            if found_priority_1:
                break

        # Uu tien 2: pivot co |gia tri| lon nhat
        if not found_priority_1:
            for r in range(n):
                if row_used[r]:
                    continue
                for c in range(n):
                    if not col_used[c] and abs(Aug[r, c]) > max_val:
                        max_val = abs(Aug[r, c])
                        p, q = r, c

        if max_val < 1e-9:
            print(f"Loi: khong tim duoc pivot khac 0 o lan lap {step} -> A suy bien, khong kha nghich!")
            return

        row_used[p] = True
        col_used[q] = True
        ind[p] = q

        print(f"Lan lap {step}: Chon pivot a[{p + 1}][{q + 1}] = {Aug[p, q]:.4f}")

        # Chuan hoa hang chua pivot
        Aug[p, :] /= Aug[p, q]

        # Khu ca ben tren va ben duoi
        factors = Aug[:, q].copy()
        factors[p] = 0.0
        Aug -= np.outer(factors, Aug[p, :])

        if step == 1:
            print("\n--- MA TRAN SAU LAN LAP 1 ---")
            print_matrix(Aug)

        step += 1

    print("\n--- MA TRAN SAU KHI KHU GAUSS-JORDAN (chua sap xep lai theo dung thu tu an) ---")
    print_matrix(Aug)

    # Sap xep lai hang r ve dung vi tri an ind[r] de duoc [E | A^-1]
    inv_A = np.zeros((n, n))
    for r in range(n):
        inv_A[ind[r], :] = Aug[r, n:]

    print("--- MA TRAN NGHICH DAO A^-1 ---")
    print_matrix(inv_A)


if __name__ == "__main__":
    main()
