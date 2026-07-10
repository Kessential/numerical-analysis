import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            v = 0.0 if abs(val) < 1e-9 else val
            print(f"{v:14.8f}", end=" ")
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


def full_svd(A, eps=1e-9, max_iter=500):
    """Khai trien ky di day du A = U Sigma V^T (algo/svd.md muc 3): chay
    Phan 1 (sigma_i, v_i) cho toi khi het hang cua A, roi suy u_i = A v_i/sigma_i."""
    r_max = min(A.shape)
    sigmas, vs = singular_values_right_vectors(A, r_max, eps=eps, max_iter=max_iter)
    us = [A @ v / sig for sig, v in zip(sigmas, vs)]

    U = np.column_stack(us) if us else np.zeros((A.shape[0], 0))
    S = np.array(sigmas)
    V = np.column_stack(vs) if vs else np.zeros((A.shape[1], 0))
    return U, S, V


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

    U, S, V = full_svd(A, eps=eps)
    r = len(S)
    print(f"\n=== Hang xap xi (so gia tri ky di khac 0 tim duoc): r = {r} ===")

    print("U =")
    print_matrix(U)
    print("Sigma (duong cheo) =")
    print_matrix(S.reshape(1, -1))
    print("V =")
    print_matrix(V)

    A_rec = U @ np.diag(S) @ V.T
    print("A tai tao lai U.Sigma.V^T =")
    print_matrix(A_rec)
    err = np.linalg.norm(A - A_rec, 2)
    # sai so tich luy qua r buoc xuong thang lien tiep, nen dung nguong long
    # hon eps cua tung buoc PP luy thua rieng le (khong phai eps truc tiep)
    tol = max(eps * 100, 1e-8)
    if err <= tol:
        print(f"||A - U.Sigma.V^T||_2 = {err:.3e} <= {tol:.1e}: khai trien dung.")
    else:
        print(f"||A - U.Sigma.V^T||_2 = {err:.3e} > {tol:.1e}: co buoc PP luy thua o Phan 2 chua hoi tu.")

    print("--- (Doi chieu) numpy.linalg.svd ---")
    U_np, S_np, Vt_np = np.linalg.svd(A, full_matrices=False)
    print_matrix(S_np[:r].reshape(1, -1))


if __name__ == "__main__":
    main()
