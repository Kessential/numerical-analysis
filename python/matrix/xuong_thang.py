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
    history = []

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

        th1_stable = (
            lam1_est is not None
            and spread < eps * max(1.0, abs(lam1_est))
            and prev_lam1 is not None
            and abs(lam1_est - prev_lam1) < eps * max(1.0, abs(lam1_est))
        )

        # Luon tinh truoc he p,q (dung de phan biet TH2/TH3) ngay ca khi TH1 co
        # ve da on dinh: neu p,q cua chinh vong lap nay cung cho thay dau hieu
        # 1 cap gia tri rieng doi nhau/phuc lien hop (xem Bug 1/2, bug_report.md),
        # thi KHONG duoc chot TH1 voi do on dinh tinh co cua ti so y1[i]/x[i].
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

        # Nguong "canh bao" phai LONG hon nhieu so voi eps dung de chot TH2/TH3
        # (o duoi): p hoi tu ve 0 CHAM hon ti so y1[i]/x[i] on dinh, nen tai
        # dung luc TH1 "co ve" da on dinh (2 vong lien tiep), p thuong van con
        # cach 0 mot khoang > eps - dung eps truc tiep se khong bat duoc dau
        # hieu nay (xem Bug 1/2, bug_report.md).
        guard_thresh = np.sqrt(eps)
        pq_indicates_pair = False
        if pq is not None:
            p_chk, q_chk = pq
            delta_chk = p_chk * p_chk - 4 * q_chk
            if delta_chk >= 0:
                sq_chk = np.sqrt(delta_chk)
                lam_a_chk, lam_b_chk = (p_chk + sq_chk) / 2, (p_chk - sq_chk) / 2
                if abs(p_chk) < guard_thresh * max(1.0, abs(lam_a_chk), abs(lam_b_chk)):
                    pq_indicates_pair = True
            else:
                pq_indicates_pair = True

        history.append((k, lam1_est, x.copy()))

        if th1_stable and not pq_indicates_pair:
            v1, _ = normalize_inf(y1)
            return {"case": 1, "iterations": k, "eigenvalues": [lam1_est], "eigenvectors": [v1], "history": history}

        if lam1_est is not None and spread < eps * max(1.0, abs(lam1_est)):
            prev_lam1 = lam1_est
        else:
            prev_lam1 = None

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
                    # Dung guard_thresh (long hon eps, xem Bug 1/2 bug_report.md) - p hoi
                    # tu ve 0 cham hon do on dinh cua chinh he p,q.
                    if abs(p) < guard_thresh * max(1.0, abs(lam_a), abs(lam_b)):
                        v_b = normalize_inf((y1 - lam_a * x).astype(complex))[0]
                        v_a = normalize_inf((y1 - lam_b * x).astype(complex))[0]
                        return {"case": 2, "iterations": k, "eigenvalues": [lam_a, lam_b], "eigenvectors": [v_a, v_b], "history": history}
                    lam1, lam_other = (lam_a, lam_b) if abs(lam_a) >= abs(lam_b) else (lam_b, lam_a)
                    v1 = normalize_inf(y1 - lam_other * x)[0]
                    return {"case": 1, "iterations": k, "eigenvalues": [lam1], "eigenvectors": [v1], "history": history}
                sq = np.sqrt(-delta)
                lam_a, lam_b = complex(p / 2, sq / 2), complex(p / 2, -sq / 2)
                v_b = normalize_inf((y1 - lam_a * x).astype(complex))[0]
                v_a = normalize_inf((y1 - lam_b * x).astype(complex))[0]
                return {"case": 3, "iterations": k, "eigenvalues": [lam_a, lam_b], "eigenvectors": [v_a, v_b], "history": history}
            prev_pq = (p, q)
        else:
            prev_pq = None

        x, scale = normalize_inf(y1)
        if scale == 0.0:
            break

    return {"case": 0, "iterations": max_iter, "eigenvalues": [], "eigenvectors": [], "history": history}


def print_history(history):
    """In 2 lan lap dau va 2 lan lap cuoi cua qua trinh PP luy thua (show truncation,
    cung pattern voi cholesky.py/lu.py)."""
    show = sorted(set(range(min(2, len(history)))) | set(range(max(0, len(history) - 2), len(history))))
    for idx in show:
        k, lam1_est, x_k = history[idx]
        lam_str = f"{lam1_est:.6f}" if lam1_est is not None else "(chua uoc luong duoc)"
        print(f"--- Lan lap {k}: uoc luong lambda1 (TH1, ti so y1[i]/x[i]) = {lam_str} ---")
        print("x^(k) (vector lap, da chuan hoa) =")
        print_matrix(x_k)


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
    print_history(result["history"])
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

    # Du phong: du chi yeu cau so_gtr=1, van chay them Buoc 2 (khong in day du,
    # khong tinh vao ket qua) de xac nhan lambda_1 that su la gtr troi DON nhat -
    # neu khong, PP luy thua o Buoc 1 co the da nham TH2/TH3 thanh TH1 (Bug 2,
    # bug_report.md) va cham lang bo sot 1 gtr rieng khac cung module.
    steps_to_run = max(so_gtr, 2)
    for step in range(2, steps_to_run + 1):
        official = step <= so_gtr
        if official:
            print(f"\n--- Buoc {step}: PP luy thua tim gtr troi tiep theo tren ma tran da xuong thang ---")
        else:
            print(f"\n--- Buoc {step} (kiem tra du phong, ngoai pham vi so_gtr={so_gtr} da yeu cau) ---")

        res = power_method(A_cur, eps=eps)
        print_history(res["history"])

        if res["case"] == 0:
            if official:
                print("Khong hoi tu: dung lai o day.")
            else:
                print("CANH BAO: kiem tra du phong khong hoi tu -> khong the xac nhan lambda_1 la gtr troi DON"
                      " nhat; ket qua o Buoc 1 co the chua day du.")
            break
        if res["case"] != 1:
            if official:
                print("Gap cap gtr troi (TH2/TH3) o buoc nay, khong xuong thang tiep duoc (can rank-2):")
                for lam, u in zip(res["eigenvalues"], res["eigenvectors"]):
                    v = u
                    for lam_k, v_k, a_s_k in reversed(chain):
                        v = recover_eigenvector(v, lam_k, lam, v_k, a_s_k)
                    v, _ = normalize_inf(v)
                    print(f"lambda = {lam}")
                    print_matrix(v)
                    eigen_goc.append((lam, v))
            else:
                print("CANH BAO: kiem tra du phong phat hien them 1 cap gtr rieng (TH2/TH3) cung module voi"
                      " lambda_1 ma yeu cau so_gtr=1 se bo sot -> lambda_1 o Buoc 1 CO THE KHONG PHAI la gtr"
                      " troi DON nhat.")
            break

        lam_k, u_k = res["eigenvalues"][0], res["eigenvectors"][0]

        if not official:
            # Du phong xac nhan lambda_1 dung la TH1 don, khong co gi bat thuong -> khong can lam gi them.
            break

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
