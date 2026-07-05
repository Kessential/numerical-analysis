import numpy as np


def print_matrix(mat):
    for row in mat:
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def main():
    filename = "test.txt"
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

    # 3. Thuat toan Gauss - quy trinh thuan
    ind = [-1] * m
    i, j = 0, 0
    step = 1

    print("--- LICH SU CHON PHAN TU KHOA (PIVOT) ---")

    while i < m and j < cols_a:
        if abs(Aug[i, j]) < 1e-9:
            found = False
            for t in range(i + 1, m):
                if abs(Aug[t, j]) > 1e-9:
                    Aug[[i, t]] = Aug[[t, i]]
                    found = True
                    break
            if not found:
                j += 1
                continue

        ind[i] = j
        print(f"Lan lap {step}: Chon pivot a[{i + 1}][{j + 1}] = {Aug[i, j]:.4f}")

        if i == m - 1:
            break

        # Khu cac phan tu ben duoi (chay tiep den het totalCols de tru ca phan B)
        factors = Aug[i + 1 :, j] / Aug[i, j]
        Aug[i + 1 :, j:] -= np.outer(factors, Aug[i, j:])

        if step == 1:
            print("\n--- MA TRAN SAU LAN LAP 1 ---")
            print_matrix(Aug)

        if j == cols_a - 1:
            break

        i += 1
        j += 1
        step += 1

    print("\n--- MA TRAN SAU LAN LAP CUOI (KET THUC QUY TRINH THUAN) ---")
    print_matrix(Aug)

    # 4. Kiem tra va quy trinh nghich (ho tro vo so nghiem)
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

        X = np.zeros((cols_a, cols_b))
        for r in range(m - 1, -1, -1):
            pivot_col = ind[r]
            if pivot_col != -1 and pivot_col < cols_a:
                s = Aug[r, pivot_col + 1 : cols_a] @ X[pivot_col + 1 : cols_a, :]
                X[pivot_col, :] = (Aug[r, cols_a : cols_a + cols_b] - s) / Aug[r, pivot_col]

        print("--- MA TRAN NGHIEM X ---")
        print_matrix(X)

    else:
        print("\n============================================")
        print(" KET LUAN: HE CO VO SO NGHIEM")
        print("============================================")

        free_vars = [c for c in range(cols_a) if not is_basic[c]]
        print(f"So an tu do: {len(free_vars)} (Gom cac an: " + ", ".join(f"x{f + 1}" for f in free_vars) + ")\n")

        # Buoc A: tinh cac vecto co so (V_i) cua khong gian nghiem thuan nhat
        V = np.zeros((len(free_vars), cols_a))
        for vi, f in enumerate(free_vars):
            V[vi, f] = 1.0
            for r in range(m - 1, -1, -1):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    s = Aug[r, p_col + 1 : cols_a] @ V[vi, p_col + 1 : cols_a]
                    V[vi, p_col] = -s / Aug[r, p_col]

        # Buoc B: tinh nghiem rieng (X0) cho tung cot cua B va in ket qua
        for cb in range(cols_b):
            print(f">>> XET MA TRAN B COT THU {cb + 1}:")
            X0 = np.zeros(cols_a)

            for r in range(m - 1, -1, -1):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    s = Aug[r, p_col + 1 : cols_a] @ X0[p_col + 1 : cols_a]
                    X0[p_col] = (Aug[r, cols_a + cb] - s) / Aug[r, p_col]

            print("Nghiem tong quat co dang: X = X0" + "".join(f" + t{vi + 1}*V{vi + 1}" for vi in range(len(free_vars))))
            print("\nBang toa do cac vector:")

            header = "An       X0(Rieng)" + "".join(f"     V{vi + 1}(t{vi + 1})" for vi in range(len(free_vars)))
            print(header)
            print("-" * 49)

            for jx in range(cols_a):
                x0_val = 0.0 if abs(X0[jx]) < 1e-9 else X0[jx]
                line = f"x{jx + 1}{x0_val:15.4f}"
                for vi in range(len(free_vars)):
                    v_val = 0.0 if abs(V[vi, jx]) < 1e-9 else V[vi, jx]
                    line += f"{v_val:13.4f}"
                print(line)
            print("-" * 49 + "\n")


if __name__ == "__main__":
    main()
