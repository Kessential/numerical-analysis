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
    """Ban sao doc lap (xem svd_gia_tri_ky_di.py / algo/svd.md muc 1).
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


def condition_number(A, eps=1e-9, max_iter=500):
    """cond(A) = sigma_max/sigma_min (algo/svd.md muc 4):
    sigma_max^2 = gtr troi cua M=A^T A (PP luy thua thuan);
    sigma_min^2 = 1/(gtr troi cua M^-1) (PP luy thua nghich dao)."""
    M = A.T @ A
    lam_max, _, it_max, conv_max = power_method_symmetric(M, eps=eps, max_iter=max_iter)
    if not conv_max:
        print(f"Canh bao: PP luy thua tim sigma_max khong hoi tu sau {max_iter} lan lap.")

    M_inv = np.linalg.inv(M)
    mu_max, _, it_min, conv_min = power_method_symmetric(M_inv, eps=eps, max_iter=max_iter)
    if not conv_min:
        print(f"Canh bao: PP luy thua nghich dao tim sigma_min khong hoi tu sau {max_iter} lan lap.")
    lam_min = 1.0 / mu_max

    sigma_max = np.sqrt(lam_max)
    sigma_min = np.sqrt(lam_min)
    return sigma_max, sigma_min, sigma_max / sigma_min, it_max, it_min


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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}) de tinh so dieu kien theo cach nay!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    if abs(np.linalg.det(A)) < 1e-12:
        print("A suy bien (det ~ 0): cond(A) = +inf, khong tinh tiep duoc.")
        return

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-9): "))

    sigma_max, sigma_min, cond, it_max, it_min = condition_number(A, eps=eps)

    print(f"\nsigma_max = {sigma_max:.6f} (PP luy thua tren A^T A, {it_max} lan lap)")
    print(f"sigma_min = {sigma_min:.6f} (PP luy thua nghich dao tren (A^T A)^-1, {it_min} lan lap)")
    print(f"cond(A) = sigma_max/sigma_min = {cond:.6f}")

    print(f"--- (Doi chieu) numpy.linalg.cond(A) = {np.linalg.cond(A):.6f} ---")


if __name__ == "__main__":
    main()
