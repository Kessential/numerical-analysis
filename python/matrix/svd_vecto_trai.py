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


def deflate(M, v1):
    v1 = np.asarray(v1, dtype=float)
    s = int(np.argmax(np.abs(v1)))
    v1 = v1 / v1[s]
    a_s = M[s, :].copy()
    B = M - np.outer(v1, a_s)
    return B, s, a_s, v1


def recover_eigenvector(u, lam1, lam_target, v1, a_s):
    return (lam1 - lam_target) * u - (a_s @ u) * v1


def singular_values_right_vectors(A, k, eps=1e-9, max_iter=500):
    """Xem svd_gia_tri_ky_di.py (Phan 1, algo/svd.md muc 1)."""
    M = A.T @ A
    M_cur = M
    chain = []
    sigmas, vs = [], []

    for step in range(1, k + 1):
        lam, v_cur, _, converged = power_method_symmetric(M_cur, eps=eps, max_iter=max_iter)
        if not converged:
            print(f"Canh bao: PP luy thua khong hoi tu sau {max_iter} lan lap o buoc {step}"
                  " (2 gia tri ky di ke nhau qua gan nhau) - dung lai o day.")
            break
        if lam <= eps:
            break

        v = v_cur
        for lam_j, v_j, a_j in reversed(chain):
            v = recover_eigenvector(v, lam_j, lam, v_j, a_j)
        v, _ = normalize_inf(v)
        v = v / np.linalg.norm(v)

        sigmas.append(np.sqrt(lam))
        vs.append(v)

        if step < k:
            M_next, s, a_s, v1n = deflate(M_cur, v_cur)
            chain.append((lam, v1n, a_s))
            M_cur = M_next

    return sigmas, vs


def left_singular_vectors(A, sigmas, vs):
    """u_i = A v_i / sigma_i (algo/svd.md muc 2)."""
    return [A @ v / sig for sig, v in zip(sigmas, vs)]


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
    r_max = int(input(f"Nhap so gia tri ky di can tim (toi da {min(m, n)}): "))

    sigmas, vs = singular_values_right_vectors(A, r_max, eps=eps)
    us = left_singular_vectors(A, sigmas, vs)

    print(f"\n=== Vector ky di trai u_i = A.v_i / sigma_i ({len(us)} vector) ===")
    for i, (sig, u) in enumerate(zip(sigmas, us), start=1):
        print(f"sigma_{i} = {sig:.6f}, u_{i} =")
        print_matrix(u)

    if len(us) > 1:
        U = np.column_stack(us)
        gram = U.T @ U
        print("--- Kiem tra truc chuan: U^T U (phai xap xi I) ---")
        print_matrix(gram)

    print("--- (Doi chieu) numpy.linalg.svd, cot dau cua U (co the lech dau +-1 so voi tren) ---")
    U_np, S_np, _ = np.linalg.svd(A, full_matrices=False)
    print_matrix(U_np[:, : len(us)])


if __name__ == "__main__":
    main()
