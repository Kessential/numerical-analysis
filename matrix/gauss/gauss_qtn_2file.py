import numpy as np


def print_matrix(mat):
    for row in mat:
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def quy_trinh_nghich_doc_lap(Aug, m, cols_a, cols_b):
    # 1. Tu dong quet tim vi tri pivot (do khong chay QTT nen phai tu tim)
    ind = [-1] * m
    for i in range(m):
        for j in range(cols_a):
            if abs(Aug[i, j]) > 1e-9:
                ind[i] = j
                break

    # 2. Kiem tra he vo nghiem
    is_consistent = True
    for r in range(m):
        if ind[r] == -1:
            if np.any(np.abs(Aug[r, cols_a : cols_a + cols_b]) > 1e-9):
                is_consistent = False
                break

    if not is_consistent:
        print("============================================")
        print(" KET LUAN: HE PHUONG TRINH VO NGHIEM!")
        print("============================================")
        return

    # 3. Tinh hang va phan loai nghiem
    is_basic = [False] * cols_a
    rank = 0
    for r in range(m):
        if ind[r] != -1 and ind[r] < cols_a:
            is_basic[ind[r]] = True
            rank += 1

    # 4. Giai va in ket qua theo phan nhanh
    if rank == cols_a:
        print("============================================")
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
        print("============================================")
        print(" KET LUAN: HE CO VO SO NGHIEM")
        print("============================================")

        free_vars = [j for j in range(cols_a) if not is_basic[j]]

        print(f"So an tu do: {len(free_vars)} (Gom cac an: " + ", ".join(f"x{f + 1}" for f in free_vars) + ")\n")

        # Vecto co so
        V = np.zeros((len(free_vars), cols_a))
        for i, f in enumerate(free_vars):
            V[i, f] = 1.0
            for r in range(m - 1, -1, -1):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    s = Aug[r, p_col + 1 : cols_a] @ V[i, p_col + 1 : cols_a]
                    V[i, p_col] = -s / Aug[r, p_col]

        # Nghiem rieng
        for cb in range(cols_b):
            print(f">>> XET MA TRAN B COT THU {cb + 1}:")
            X0 = np.zeros(cols_a)

            for r in range(m - 1, -1, -1):
                p_col = ind[r]
                if p_col != -1 and p_col < cols_a:
                    s = Aug[r, p_col + 1 : cols_a] @ X0[p_col + 1 : cols_a]
                    X0[p_col] = (Aug[r, cols_a + cb] - s) / Aug[r, p_col]

            print("Bang toa do cac vector:")
            header = "An       X0(Rieng)" + "".join(f"     V{i + 1}(t{i + 1})" for i in range(len(free_vars)))
            print(header)
            print("-" * 49)

            for j in range(cols_a):
                x0_val = 0.0 if abs(X0[j]) < 1e-9 else X0[j]
                line = f"x{j + 1}{x0_val:19.8f}"
                for i in range(len(free_vars)):
                    v_val = 0.0 if abs(V[i, j]) < 1e-9 else V[i, j]
                    line += f"{v_val:17.8f}"
                print(line)
            print("-" * 49 + "\n")


def main():
    filename_a = "testA1.txt"
    filename_b = "testB1.txt"
    try:
        A = np.loadtxt(filename_a, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename_a}'!")
        return

    try:
        B = np.loadtxt(filename_b, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename_b}'!")
        return

    if A.size == 0 or B.size == 0:
        print("Loi: File rong!")
        return

    m, cols_a = A.shape
    m_b, cols_b = B.shape

    if m_b != m:
        print(f"Loi: Ma tran A ({m} hang) va ma tran B ({m_b} hang) khong khop so hang!")
        return

    Aug = np.hstack([A, B])
    print(f"Ma tran A: {m} x {cols_a}, ma tran B: {m} x {cols_b}.")

    print("\n--- MA TRAN DAU VAO ---")
    print_matrix(Aug)

    quy_trinh_nghich_doc_lap(Aug, m, cols_a, cols_b)


if __name__ == "__main__":
    main()
