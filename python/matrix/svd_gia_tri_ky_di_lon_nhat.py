import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            v = 0.0 if abs(val) < 1e-9 else val
            print(f"{v:10.4f}", end=" ")
        print()
    print()


def normalize_inf(v):
    scale = v[np.argmax(np.abs(v))]
    if abs(scale) < 1e-300:
        return v, 0.0
    return v / scale, scale


def power_method_symmetric(M, eps=1e-9, max_iter=500):
    """PP luy thua ap dung cho ma tran doi xung ban xac dinh duong M (vd M=A^T A):
    gia tri rieng luon thuc va >=0 nen chi co the roi vao TH1 (xem algo/svd.md
    muc 1) - khong can xu ly TH2/TH3 nhu ban tong quat o luy_thua.py.
    Tra ve (lambda, v, so lan lap, converged)."""
    n = M.shape[0]
    x = np.ones(n)
    x, _ = normalize_inf(x)

    prev_lam = None
    for k in range(1, max_iter + 1):
        y = M @ x
        mask = np.abs(x) > 1e-8
        if not np.any(mask):
            break
        ratios = y[mask] / x[mask]
        spread = np.max(ratios) - np.min(ratios) if ratios.size > 1 else 0.0
        lam = np.mean(ratios)

        if spread < eps * max(1.0, abs(lam)):
            if prev_lam is not None and abs(lam - prev_lam) < eps * max(1.0, abs(lam)):
                v, _ = normalize_inf(y)
                return lam, v, k, True
            prev_lam = lam
        else:
            prev_lam = None

        x, scale = normalize_inf(y)
        if scale == 0.0:
            return 0.0, x, k, True

    y = M @ x
    lam = float(x @ y)
    return lam, x, max_iter, False


def largest_singular_triplet(A, eps=1e-9, max_iter=500):
    """Tim sigma_1 lon nhat va cap vecto ky di v_1 (phai), u_1 (trai) tuong
    ung, bang PP luy thua thuan tren M=A^T A - KHONG xuong thang
    (algo/svd.md muc 1). Day la buoc co so cho PP xuong thang (muc 2).
    Tra ve None neu A = 0 (khong co gia tri ky di khac 0)."""
    M = A.T @ A
    lam1, v1, it, converged = power_method_symmetric(M, eps=eps, max_iter=max_iter)
    if not converged:
        print(f"Canh bao: PP luy thua khong hoi tu sau {max_iter} lan lap.")

    if lam1 <= eps:
        print("A = 0 (hoac gan bang 0): khong co gia tri ky di khac 0.")
        return None

    v1 = v1 / np.linalg.norm(v1)
    sigma1 = np.sqrt(lam1)
    u1 = A @ v1 / sigma1
    return sigma1, v1, u1, it


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

    m, n = A.shape
    print(f"--- MA TRAN A ({m} x {n}) ---")
    print_matrix(A)

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-9): "))

    result = largest_singular_triplet(A, eps=eps)
    if result is None:
        return
    sigma1, v1, u1, it = result

    print(f"\n=== Gia tri ky di lon nhat (hoi tu sau {it} lan lap) ===")
    print(f"sigma_1 = {sigma1:.6f}")
    print("v_1 (vector ky di phai) =")
    print_matrix(v1)
    print("u_1 (vector ky di trai) =")
    print_matrix(u1)

    print("--- (Doi chieu) numpy.linalg.svd ---")
    U_np, S_np, Vt_np = np.linalg.svd(A, full_matrices=False)
    print(f"sigma_1 (numpy) = {S_np[0]:.6f}")


if __name__ == "__main__":
    main()
