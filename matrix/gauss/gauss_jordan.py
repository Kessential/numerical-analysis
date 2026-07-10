import numpy as np

def print_matrix(mat):
    for row in mat:
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def main():
    filename = "test_vo_so_nghiem.txt"
    try:
        Aug = np.loadtxt(filename, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: {filename}'. Kiem tra lai duong dan!")
        return

    if Aug.size == 0:
        print("Loi: File rong hoac khong chua du lieu hop le!")
        return

    m, total_cols = Aug.shape

    print(f"Ma tran mo rong co tong cong {total_cols} cot.")
    cols_b = int(input("Nhap so cot cua ma tran B (vi du p = 1, 2...): "))
    cols_a = total_cols - cols_b

    print(f"=> Ma tran A co kich thuoc {m} x {cols_a}.")

    # 3. Thuat toan Gauss-Jordan (khu toan dien)
    ind = [-1] * m
    row_used = [False] * m
    col_used = [False] * cols_a
    step = 1

    print("--- LICH SU CHON PHAN TU KHOA (PIVOT) ---")

    for _ in range(min(m, cols_a)):
        max_val = 0.0
        p, q = -1, -1
        found_priority_1 = False

        # Uu tien 1: pivot co |gia tri| = 1
        for r in range(m):
            if row_used[r]:
                continue
            for c in range(cols_a):
                if not col_used[c] and abs(abs(Aug[r, c]) - 1) < 1e-9:
                    p, q = r, c
                    max_val = abs(Aug[r, c])
                    found_priority_1 = True
                    break
            if found_priority_1:
                break

        # Uu tien 2: pivot co |gia tri| lon nhat
        if not found_priority_1:
            for r in range(m):
                if row_used[r]:
                    continue
                for c in range(cols_a):
                    if not col_used[c] and abs(Aug[r, c]) > max_val:
                        max_val = abs(Aug[r, c])
                        p, q = r, c

        if max_val < 1e-9:
            break

        row_used[p] = True
        col_used[q] = True
        ind[p] = q

        print(f"Lan lap {step}: Chon pivot a[{p + 1}][{q + 1}] = {Aug[p, q]:.8f}")

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

    print("\n--- MA TRAN KET THUC KHOA GAUSS-JORDAN (RREF) ---")
    print_matrix(Aug)

    # 4. Kiem tra nghiem tu ma tran rut gon
    is_consistent = True
    for r in range(m):
        all_zero_in_a = not np.any(np.abs(Aug[r, :cols_a]) > 1e-9)
        if all_zero_in_a and np.any(np.abs(Aug[r, cols_a : cols_a + cols_b]) > 1e-9):
            is_consistent = False
            break

    if not is_consistent:
        print("\n============================================")
        print(" KET LUAN: HE PHUONG TRINH VO NGHIEM!")
        print("============================================")
        return

    is_basic = [False] * cols_a
    rank = 0
    for r in range(m):
        if ind[r] != -1 and ind[r] < cols_a:
            is_basic[ind[r]] = True
            rank += 1

    if rank == cols_a:
        print("\n============================================")
        print(" KET LUAN: HE CO NGHIEM DUY NHAT")
        print("============================================")

        # Gauss-Jordan: khong can the nguoc, nghiem nam thang o cot B
        X = np.zeros((cols_a, cols_b))
        for r in range(m):
            pivot_col = ind[r]
            if pivot_col != -1 and pivot_col < cols_a:
                X[pivot_col, :] = Aug[r, cols_a : cols_a + cols_b]

        print("--- MA TRAN NGHIEM X ---")
        print_matrix(X)

    else:
        print("\n============================================")
        print(" KET LUAN: HE CO VO SO NGHIEM")
        print("============================================")

        free_vars = [c for c in range(cols_a) if not is_basic[c]]
        print(f"So an tu do: {len(free_vars)} (Gom cac an: " + ", ".join(f"x{f + 1}" for f in free_vars) + ")\n")

        # Buoc A: vecto co so - khong can the nguoc
        V = np.zeros((len(free_vars), cols_a))
        for vi, f in enumerate(free_vars):
            V[vi, f] = 1.0
            for r in range(m):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    V[vi, p_col] = -Aug[r, f]

        # Buoc B: nghiem rieng (X0) cho tung cot cua B
        for cb in range(cols_b):
            print(f">>> XET MA TRAN B COT THU {cb + 1}:")
            X0 = np.zeros(cols_a)

            for r in range(m):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    X0[p_col] = Aug[r, cols_a + cb]

            print("Nghiem tong quat co dang: X = X0" + "".join(f" + t{vi + 1}*V{vi + 1}" for vi in range(len(free_vars))))
            print("\nBang toa do cac vector:")

            header = "An       X0(Rieng)" + "".join(f"     V{vi + 1}(t{vi + 1})" for vi in range(len(free_vars)))
            print(header)
            print("-" * 49)

            for jx in range(cols_a):
                x0_val = 0.0 if abs(X0[jx]) < 1e-9 else X0[jx]
                line = f"x{jx + 1}{x0_val:19.8f}"
                for vi in range(len(free_vars)):
                    v_val = 0.0 if abs(V[vi, jx]) < 1e-9 else V[vi, jx]
                    line += f"{v_val:17.8f}"
                print(line)
            print("-" * 49 + "\n")


if __name__ == "__main__":
    main()
