import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def all_principal_minors_invertible(A):
    # Buoc 0 (doc): kiem tra truoc moi A_k (k=1..n, goc tren trai) co kha nghich
    # khong, bang cach tinh dinh thuc tung A_k.
    n = A.shape[0]
    for k in range(1, n + 1):
        if abs(np.linalg.det(A[:k, :k])) < 1e-9:
            return False
    return True


def bordering_inverse(M):
    # PP vien quanh thuan tuy: doi hoi moi ma tran con chinh M_k (k=1..n) kha
    # nghich. main() da dam bao dieu nay tu Buoc 0 truoc khi goi ham nay (chon
    # M=A neu moi A_k kha nghich, nguoc lai M=A^T.A vi luon xac dinh duong).
    n = M.shape[0]

    if abs(M[0, 0]) < 1e-9:
        print("Loi: m[1][1] = 0 -> khong the bat dau PP vien quanh!")
        return None

    inv = np.array([[1.0 / M[0, 0]]])

    show = sorted(set(range(min(2, n))) | set(range(max(0, n - 2), n)))

    print("--- Buoc k = 1 ---")
    print("M_1^-1 =")
    print_matrix(inv)

    skipped_msg_done = False
    for k in range(2, n + 1):
        alpha_col = M[: k - 1, k - 1]  # alpha_{k-1,1}: cot cuoi, tru phan tu goc
        alpha_row = M[k - 1, : k - 1]  # alpha_{1,k-1}: hang cuoi, tru phan tu goc
        a_kk = M[k - 1, k - 1]

        s = a_kk - alpha_row @ inv @ alpha_col
        if abs(s) < 1e-9:
            print(f"Loi: m[{k}][{k}] - alpha.M^-1.alpha = 0 tai buoc k={k} -> M_{k} khong kha nghich!")
            return None

        b_kk = 1.0 / s
        inv_alpha_col = inv @ alpha_col
        alpha_row_inv = alpha_row @ inv

        beta_col = -b_kk * inv_alpha_col  # beta_{k-1,1}
        beta_row = -b_kk * alpha_row_inv  # beta_{1,k-1}
        B = inv + b_kk * np.outer(inv_alpha_col, alpha_row_inv)  # xem nghich_dao.md muc 4 (da sua dau)

        new_inv = np.zeros((k, k))
        new_inv[: k - 1, : k - 1] = B
        new_inv[: k - 1, k - 1] = beta_col
        new_inv[k - 1, : k - 1] = beta_row
        new_inv[k - 1, k - 1] = b_kk

        inv = new_inv

        if k - 1 in show:
            print(f"--- Buoc k = {k} ---")
            print(f"s = m_{k}{k} - alpha_(1,{k-1}).M_{k-1}^-1.alpha_({k-1},1) = {s:.8f}")
            print(f"b_{k}{k} = 1/s = {b_kk:.8f}")
            print(f"beta_({k-1},1) (cot) =")
            print_matrix(beta_col.reshape(-1, 1))
            print(f"beta_(1,{k-1}) (hang) =")
            print_matrix(beta_row.reshape(1, -1))
            print(f"B_{k-1} =")
            print_matrix(B)
            print(f"M_{k}^-1 =")
            print_matrix(inv)
        elif not skipped_msg_done:
            print("... (bo qua cac buoc giua, chi hien 2 buoc dau va 2 buoc cuoi) ...")
            skipped_msg_done = True

    return inv


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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the ap dung PP vien quanh!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    # Buoc 0: kiem tra truoc moi A_k co kha nghich khong, chon M va co assym.
    if all_principal_minors_invertible(A):
        print(">>> Moi A_k (k=1..n) deu kha nghich -> ap dung PP vien quanh truc tiep tren A.")
        M, assym = A, False
    else:
        # A kha nghich bat ky khong dam bao moi ma tran con chinh A_k cung kha
        # nghich (vd A=[[0,1],[1,2]]: det A=-1 nhung a11=0). Theo slide: M=A^T.A
        # luon xac dinh duong khi A kha nghich (x^t.M.x = ||Ax||^2 >= 0, "=" <=>
        # Ax=0 <=> x=0), nen moi dinh thuc con chinh cua M deu duong (Sylvester)
        # -> PP vien quanh tren M luon thuc hien duoc. Sau do A^-1 = M^-1.A^T.
        print(">>> Khong phai moi A_k deu kha nghich -> chuyen sang dung M = A^T.A:")
        M, assym = A.T @ A, True
        print("--- MA TRAN M = A^T.A ---")
        print_matrix(M)

    inv_M = bordering_inverse(M)
    if inv_M is None:
        print("Loi: A khong kha nghich!")
        return

    inv_A = inv_M @ A.T if assym else inv_M

    print("--- MA TRAN NGHICH DAO A^-1 ---")
    print_matrix(inv_A)


if __name__ == "__main__":
    main()
