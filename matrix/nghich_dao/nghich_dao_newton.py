import numpy as np
from collections import deque


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def _largest_eigenvalue_symmetric(M, eps=1e-9, max_iter=500):
    """PP luy thua tim tri rieng troi (tri tuyet doi lon nhat) cua ma tran doi xung M."""
    n = M.shape[0]
    x = np.ones(n)
    x = x / x[np.argmax(np.abs(x))]
    lam_prev = None
    lam = 0.0
    for _ in range(max_iter):
        y = M @ x
        scale = y[np.argmax(np.abs(y))]
        if abs(scale) < 1e-300:
            return 0.0
        x = y / scale
        lam = scale
        if lam_prev is not None and abs(lam - lam_prev) < eps * max(1.0, abs(lam)):
            return lam
        lam_prev = lam
    return lam


def spectral_norm(M, eps=1e-9, max_iter=500):
    """Chuan pho ||M||_2 = sqrt(tri rieng troi cua M^T.M) (PP luy thua)."""
    lam = _largest_eigenvalue_symmetric(M.T @ M, eps=eps, max_iter=max_iter)
    return np.sqrt(max(lam, 0.0))


def _norm(M, norm_choice):
    return spectral_norm(M) if norm_choice == "2" else np.linalg.norm(M, np.inf)


def newton_inverse(A, eps, max_iter=1000, norm_choice="inf"):
    n = A.shape[0]
    E = np.eye(n)

    norm1 = np.linalg.norm(A, 1)
    norm_inf = np.linalg.norm(A, np.inf)
    X = A.T / (norm1 * norm_inf)
    norm_label = "||.||_2 (chuan pho, PP luy thua tren G^T.G)" if norm_choice == "2" else "||.||_inf"

    print(f"||A||_1 = {norm1:.6f}, ||A||_inf = {norm_inf:.6f}")
    print("--- MA TRAN BAN DAU X_0 = A^T / (||A||_1 * ||A||_inf) ---")
    print_matrix(X)
    print(f"Dieu kien hoi tu duoc kiem tra bang chuan: {norm_label}.")

    G = E - A @ X
    # Dieu kien hoi tu ly thuyet dung chuan pho (spectral norm) ||G_0||_2 < 1.
    # Voi norm_choice="inf" (mac dinh, nhat quan voi cac PP lap khac trong repo,
    # khong can tinh tri ky di) day chi la uoc luong thay the, KHONG phai dieu kien
    # du/can chat che: co truong hop ||G_0||_inf >= 1 nhung ||G_0||_2 < 1 (van hoi
    # tu), hoac nguoc lai. Voi norm_choice="2" thi day la dieu kien chuan, chinh xac.
    q = _norm(G, norm_choice)
    print(f"G_0 = E - A.X_0, {norm_label} = q = {q:.6f} -> " + ("CANH BAO (q >= 1)!" if q >= 1 else "uoc luong hoi tu (q < 1)."))

    first_two = {}
    last_two = deque(maxlen=2)

    for it in range(1, max_iter + 1):
        X = X @ (2 * E - A @ X)
        G = E - A @ X
        err = _norm(G, norm_choice)

        if it <= 2:
            first_two[it] = (X.copy(), err)
        last_two.append((it, X.copy(), err))

        if err < eps or it == max_iter:
            break

    last_its = {item[0] for item in last_two}
    for it in (1, 2):
        if it in first_two and it not in last_its:
            X_it, err_it = first_two[it]
            print(f"\n--- LAN LAP {it}: {norm_label} cua (E - A.X_n) = {err_it:.6g} ---")
            print_matrix(X_it)

    for i, (it, X_it, err_it) in enumerate(last_two):
        label = "CUOI" if i == len(last_two) - 1 else "AP CHOT"
        print(f"\n--- LAN LAP {it} ({label}): {norm_label} cua (E - A.X_n) = {err_it:.6g} ---")
        print_matrix(X_it)

    last_it, X_last, last_err = last_two[-1]
    if last_err >= eps:
        print(
            f"\nCANH BAO: sau {last_it} lan lap (toi da) van chua dat sai so yeu cau "
            f"({norm_label} cua (E-A.X_n) = {last_err:.6g} >= eps = {eps:.6g}). Ket qua co the khong dang tin cay -> "
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
    chuan = input("Chon chuan kiem tra hoi tu (1 = vo cung ||.||_inf [mac dinh], 2 = chuan pho ||.||_2): ").strip()
    norm_choice = "2" if chuan == "2" else "inf"

    X, it = newton_inverse(A, eps, norm_choice=norm_choice)

    print(f"\n=> Dung sau {it} lan lap.")
    print("--- MA TRAN NGHICH DAO XAP XI A^-1 ---")
    print_matrix(X)


if __name__ == "__main__":
    main()
