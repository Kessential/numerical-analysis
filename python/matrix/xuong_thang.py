import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            if np.iscomplexobj(val):
                re = 0.0 if abs(val.real) < 1e-9 else val.real
                im = 0.0 if abs(val.imag) < 1e-9 else val.imag
                s = f"{re:.8f}" + ("" if im == 0 else f"{'+' if im >= 0 else '-'}{abs(im):.8f}i")
                print(f"{s:>24}", end=" ")
            else:
                v = 0.0 if abs(val) < 1e-9 else val
                print(f"{v:14.8f}", end=" ")
        print()
    print()


def normalize_inf(v):
    scale = v[np.argmax(np.abs(v))]
    if abs(scale) < 1e-300:
        return v, 0.0
    return v / scale, scale


def power_method(A, x0=None, eps=1e-6, max_iter=300):
    """Ban sao doc lap cua PP luy thua (xem luy_thua.py / algo/luy_thua_xuong_thang.md),
    dung lam buoc con cho PP xuong thang: tim gtr troi (TH1) cua ma tran da xuong thang."""
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.asarray(x0, dtype=float).copy()
    x, _ = normalize_inf(x)

    prev_lam1 = None
    prev_pq = None

    for k in range(1, max_iter + 1):
        y1 = A @ x
        y2 = A @ y1

        mask = np.abs(x) > 1e-8
        if np.any(mask):
            ratios = y1[mask] / x[mask]
            spread = np.max(ratios) - np.min(ratios) if ratios.size > 1 else 0.0
            lam1_est = np.mean(ratios)
        else:
            spread, lam1_est = np.inf, None

        if lam1_est is not None and spread < eps * max(1.0, abs(lam1_est)):
            if prev_lam1 is not None and abs(lam1_est - prev_lam1) < eps * max(1.0, abs(lam1_est)):
                v1, _ = normalize_inf(y1)
                return {"case": 1, "iterations": k, "eigenvalues": [lam1_est], "eigenvectors": [v1]}
            prev_lam1 = lam1_est
        else:
            prev_lam1 = None

        order = np.argsort(-np.abs(x))
        pq = None
        for a in range(n):
            for b in range(a + 1, n):
                r, s = order[a], order[b]
                M = np.array([[-y1[r], x[r]], [-y1[s], x[s]]])
                if abs(np.linalg.det(M)) > 1e-9:
                    rhs = np.array([-y2[r], -y2[s]])
                    pq = np.linalg.solve(M, rhs)
                    break
            if pq is not None:
                break

        if pq is not None:
            p, q = pq
            if prev_pq is not None and abs(p - prev_pq[0]) < eps * max(1.0, abs(p)) and abs(q - prev_pq[1]) < eps * max(1.0, abs(q)):
                delta = p * p - 4 * q
                if delta >= 0:
                    sq = np.sqrt(delta)
                    lam_a, lam_b = (p + sq) / 2, (p - sq) / 2
                    # 2 nghiem thuc: chi la TH2 (lambda1=-lambda2) neu p=lambda1+lambda2~0,
                    # nguoc lai day chi la 2 gtr troi thuc khong dong module tinh duoc som
                    # hon phep thu ti so don gian -> tra ve dung dang TH1.
                    if abs(p) < eps * max(1.0, abs(lam_a), abs(lam_b)):
                        v_b = normalize_inf((y1 - lam_a * x).astype(complex))[0]
                        v_a = normalize_inf((y1 - lam_b * x).astype(complex))[0]
                        return {"case": 2, "iterations": k, "eigenvalues": [lam_a, lam_b], "eigenvectors": [v_a, v_b]}
                    lam1, lam_other = (lam_a, lam_b) if abs(lam_a) >= abs(lam_b) else (lam_b, lam_a)
                    v1 = normalize_inf(y1 - lam_other * x)[0]
                    return {"case": 1, "iterations": k, "eigenvalues": [lam1], "eigenvectors": [v1]}
                sq = np.sqrt(-delta)
                lam_a, lam_b = complex(p / 2, sq / 2), complex(p / 2, -sq / 2)
                v_b = normalize_inf((y1 - lam_a * x).astype(complex))[0]
                v_a = normalize_inf((y1 - lam_b * x).astype(complex))[0]
                return {"case": 3, "iterations": k, "eigenvalues": [lam_a, lam_b], "eigenvectors": [v_a, v_b]}
            prev_pq = (p, q)
        else:
            prev_pq = None

        x, scale = normalize_inf(y1)
        if scale == 0.0:
            break

    return {"case": 0, "iterations": max_iter, "eigenvalues": [], "eigenvectors": []}


def deflate(A, lam1, v1):
    """Cach chon 2 cua slide: B = A - v1.a_s, voi a_s la hang s cua A (s la
    chi so thanh phan lon nhat cua v1, da chuan hoa de v1[s] = 1). Tra ve
    (B, s, a_s). Dam bao B v1 ~ 0."""
    v1 = np.asarray(v1)
    s = int(np.argmax(np.abs(v1)))
    v1 = v1 / v1[s]
    a_s = A[s, :].astype(v1.dtype)
    B = A - np.outer(v1, a_s)
    return B, s, a_s, v1


def recover_eigenvector(u, lam1, lam_target, v1, a_s):
    """v = (lambda1-lambda_target).u - (a_s.u).v1.

    Slide ghi cong thuc tong quat v_i=(lambda1-lambda_i)u_i+lambda1(x^Tu_i)v_1 (dau +),
    nhung thay truc tiep vao A v_i = (B+lambda1.v1.x^T)v_i va dung Bu_i=lambda_i.u_i,
    B.v1=0, x^Tv1=1 se thay dau dung phai la TRU thi moi co A v_i = lambda_i v_i (da
    kiem tra bang so, xem python/algo/luy_thua_xuong_thang.md muc Luu y) - slide bi
    sai dau o cong thuc nay."""
    return (lam1 - lam_target) * u - (a_s @ u) * v1


def main():
    filename = "test.txt"
    try:
        A0 = np.loadtxt(filename, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: '{filename}'. Kiem tra lai duong dan!")
        return

    if A0.size == 0:
        print("Loi: File rong hoac khong chua du lieu hop le!")
        return

    n, cols = A0.shape
    if n != cols:
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the ap dung PP xuong thang!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A0)

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-6): "))
    so_gtr = int(input("Nhap so gia tri rieng troi can tim lien tiep (vi du 2, 3...): "))

    print("\n--- Buoc 1: PP luy thua tim gtr troi dau tien tren A ---")
    result = power_method(A0, eps=eps)
    if result["case"] == 0:
        print("Khong hoi tu: khong tim duoc gtr troi dau tien, dung thuat toan.")
        return
    if result["case"] != 1:
        print("A co ngay 1 cap gtr troi (TH2/TH3) o buoc dau tien:")
        for lam, v in zip(result["eigenvalues"], result["eigenvectors"]):
            print(f"lambda = {lam}")
            print_matrix(v)
        print("PP xuong thang (cach chon 2) o day chi xuong thang duoc tu 1 gtr troi DON (TH1);")
        print("xuong thang cho 1 cap gtr troi (TH2/TH3) can suy rong (rank-2), khong cai dat o day.")
        return

    lam1 = result["eigenvalues"][0]
    v1 = result["eigenvectors"][0]
    print(f"lambda_1 = {lam1:.6f}, hoi tu sau {result['iterations']} lan lap.")
    print_matrix(v1)

    # chain[i] = (lam_i, v_i trong khong gian ma tran dang xuong thang, a_s_i, ma_tran_i)
    A_cur = A0
    chain = []
    B, s, a_s, v1n = deflate(A_cur, lam1, v1)
    chain.append((lam1, v1n, a_s))
    print(f"\n--- Xuong thang lan 1 (s = {s + 1}): B = A - v1.a_s ---")
    print_matrix(B)
    print(f"Kiem tra B.v1 ~ 0: max|B.v1| = {np.max(np.abs(B @ v1n)):.3e}")

    eigen_goc = [(lam1, v1n)]
    A_cur = B

    for step in range(2, so_gtr + 1):
        print(f"\n--- Buoc {step}: PP luy thua tim gtr troi tiep theo tren ma tran da xuong thang ---")
        res = power_method(A_cur, eps=eps)
        if res["case"] == 0:
            print("Khong hoi tu: dung lai o day.")
            break
        if res["case"] != 1:
            print("Gap cap gtr troi (TH2/TH3) o buoc nay, khong xuong thang tiep duoc (can rank-2):")
            for lam, u in zip(res["eigenvalues"], res["eigenvectors"]):
                v = u
                for lam_k, v_k, a_s_k in reversed(chain):
                    v = recover_eigenvector(v, lam_k, lam, v_k, a_s_k)
                v, _ = normalize_inf(v)
                print(f"lambda = {lam}")
                print_matrix(v)
            break

        lam_k, u_k = res["eigenvalues"][0], res["eigenvectors"][0]
        print(f"lambda_{step} = {lam_k:.6f} (trong khong gian ma tran da xuong thang), hoi tu sau {res['iterations']} lan lap.")

        v_goc = u_k
        for lam_j, v_j, a_s_j in reversed(chain):
            v_goc = recover_eigenvector(v_goc, lam_j, lam_k, v_j, a_s_j)
        v_goc, _ = normalize_inf(v_goc)
        print(f"-> Khoi phuc ve khong gian goc cua A: vecto rieng ung lambda_{step} =")
        print_matrix(v_goc)
        eigen_goc.append((lam_k, v_goc))

        if step < so_gtr:
            B_next, s_next, a_s_next, u_kn = deflate(A_cur, lam_k, u_k)
            chain.append((lam_k, u_kn, a_s_next))
            print(f"Xuong thang lan {step} (s = {s_next + 1}): B_{step} = B_{step - 1} - u_{step}.a_s")
            print_matrix(B_next)
            A_cur = B_next

    print("\n=== TONG HOP CAC GIA TRI RIENG TROI TIM DUOC (theo thu tu) ===")
    for lam, v in eigen_goc:
        print(f"lambda = {lam:.6f}")
        print_matrix(v)

    print("--- (Doi chieu) Gia tri rieng tinh boi numpy.linalg.eigvals, sap theo |lambda| giam dan ---")
    eig_np = np.linalg.eigvals(A0)
    order = np.argsort(-np.abs(eig_np))
    print_matrix(eig_np[order].reshape(1, -1))


if __name__ == "__main__":
    main()
