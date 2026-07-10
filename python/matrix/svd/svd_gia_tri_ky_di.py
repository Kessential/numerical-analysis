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
    """PP luy thua ap dung cho ma tran doi xung ban xac dinh duong M (vd M=A^T A):
    gia tri rieng luon thuc va >=0 nen chi co the roi vao TH1 (xem algo/svd.md
    muc 1) - khong can xu ly TH2/TH3 nhu ban tong quat o luy_thua.py.
    Tra ve (lambda, v, so lan lap, converged). Neu het max_iter ma chua hoi tu
    (2 gia tri rieng ke nhau qua gan nhau -> hoi tu cham), converged=False va
    khong nen tin gia tri lambda tra ve (chi la uoc luong Rayleigh tam thoi)."""
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
    """Cach chon 2: B = M - v1.a_s, s la toa do lon nhat cua v1 (da chuan hoa
    de v1[s]=1), a_s la hang s cua M. Xem algo/luy_thua_xuong_thang.md."""
    v1 = np.asarray(v1, dtype=float)
    s = int(np.argmax(np.abs(v1)))
    v1 = v1 / v1[s]
    a_s = M[s, :].copy()
    B = M - np.outer(v1, a_s)
    return B, s, a_s, v1


def recover_eigenvector(u, lam1, lam_target, v1, a_s):
    """v = (lambda1-lambda_target).u - (a_s.u).v1 (xem xuong_thang.py)."""
    return (lam1 - lam_target) * u - (a_s @ u) * v1


def singular_values_right_vectors(A, k, eps=1e-9, max_iter=500):
    """Tim k cap (sigma_i, v_i) lon nhat bang PP luy thua + xuong thang tren
    M=A^T A (algo/svd.md muc 1). Dung lai khi lambda_k <= eps (het hang cua A)."""
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

    print(f"\n=== Tim duoc {len(sigmas)} gia tri ky di ===")
    for i, (sig, v) in enumerate(zip(sigmas, vs), start=1):
        print(f"sigma_{i} = {sig:.6f}, v_{i} =")
        print_matrix(v)

    print("--- (Doi chieu) numpy.linalg.svd ---")
    _, S_np, _ = np.linalg.svd(A)
    print_matrix(S_np[: len(sigmas)].reshape(1, -1))


if __name__ == "__main__":
    main()
