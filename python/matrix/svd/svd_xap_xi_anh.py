import numpy as np
from PIL import Image


def normalize_inf(v):
    scale = v[np.argmax(np.abs(v))]
    if abs(scale) < 1e-300:
        return v, 0.0
    return v / scale, scale


def power_method_symmetric(M, eps=1e-8, max_iter=3000):
    """Ban sao doc lap (xem svd_gia_tri_ky_di.py / algo/svd.md muc 1).
    Tra ve (lambda, v, so lan lap, converged). Pho gia tri ky di cua anh tu
    nhien thuong giam rat cham (sigma_{k+1}/sigma_k gan 1) nen can max_iter
    lon hon nhieu so voi vi du ma tran nho o cac file svd_*.py khac."""
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


def top_k_singular_triplets(A, k, eps=1e-8, max_iter=3000):
    """Chi lay k cap ky di lon nhat (algo/svd.md muc 5): dung Gram matrix
    co, giua A^T A (n x n) va A A^T (m x m), de xuong thang nhanh hon."""
    m, n = A.shape
    use_ATA = n <= m
    M = (A.T @ A) if use_ATA else (A @ A.T)

    M_cur = M
    chain = []
    sigmas, vecs = [], []

    for step in range(1, k + 1):
        lam, v_cur, _, converged = power_method_symmetric(M_cur, eps=eps, max_iter=max_iter)
        if not converged:
            print(f"Canh bao: PP luy thua khong hoi tu sau {max_iter} lan lap o buoc {step}"
                  " (2 gia tri ky di ke nhau qua gan nhau) - dung lai o day, chi lay duoc"
                  f" {step - 1} cap ky di dau tien.")
            break
        if lam <= eps:
            break

        v = v_cur
        for lam_j, v_j, a_j in reversed(chain):
            v = recover_eigenvector(v, lam_j, lam, v_j, a_j)
        v, _ = normalize_inf(v)
        v = v / np.linalg.norm(v)

        sigmas.append(np.sqrt(lam))
        vecs.append(v)

        if step < k:
            M_next, s, a_s, v1n = deflate(M_cur, v_cur)
            chain.append((lam, v1n, a_s))
            M_cur = M_next

    if use_ATA:
        vs = vecs
        us = [A @ v / sig for sig, v in zip(sigmas, vs)]
    else:
        us = vecs
        vs = [A.T @ u / sig for sig, u in zip(sigmas, us)]

    return sigmas, us, vs


def rank_k_approximation(A, k, eps=1e-8, max_iter=3000):
    sigmas, us, vs = top_k_singular_triplets(A, k, eps=eps, max_iter=max_iter)
    A_k = np.zeros_like(A)
    for sig, u, v in zip(sigmas, us, vs):
        A_k += sig * np.outer(u, v)
    return A_k, sigmas


def main():
    default_path = "../slide/image.png"
    path = input(f"Duong dan anh (Enter de dung mac dinh {default_path}): ").strip() or default_path

    try:
        img = Image.open(path).convert("L")
    except OSError:
        print(f"Loi: Khong mo duoc anh '{path}'.")
        return

    max_side = 200
    if max(img.size) > max_side:
        ratio = max_side / max(img.size)
        img = img.resize((max(1, int(img.width * ratio)), max(1, int(img.height * ratio))))

    A = np.asarray(img, dtype=float)
    m, n = A.shape
    print(f"--- Anh xam ({m} x {n} pixel), hang toi da = {min(m, n)} ---")

    k = int(input(f"Nhap hang xap xi k can dung (vi du 10, toi da {min(m, n)}): "))

    A_k, sigmas = rank_k_approximation(A, k)
    r = len(sigmas)

    err_rel = np.linalg.norm(A - A_k, "fro") / np.linalg.norm(A, "fro")
    storage_ratio = r * (m + n + 1) / (m * n)

    print(f"\nDa dung r = {r} cap ky di lon nhat (yeu cau k = {k}).")
    print(f"sigma_1..sigma_r = {[round(s, 3) for s in sigmas]}")
    print(f"Sai so tuong doi ||A-A_k||_F / ||A||_F = {err_rel:.4%}")
    print(f"Ti le luu tru r(m+n+1)/(mn) = {storage_ratio:.4%}")

    A_k_img = np.clip(A_k, 0, 255).astype(np.uint8)
    out_path = "svd_xap_xi_anh_output.png"
    Image.fromarray(A_k_img).save(out_path)
    print(f"Da luu anh xap xi hang-{r} vao '{out_path}'.")


if __name__ == "__main__":
    main()
