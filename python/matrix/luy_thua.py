import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            if isinstance(val, complex) or np.iscomplexobj(val):
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
    """PP luy thua tim gtr troi, tu dong nhan dien ca 3 truong hop:
    TH1 |l1|>|l2| (mot gtr troi thuc); TH2/TH3 |l1|=|l2|>|l3| (dung chung
    1 cong thuc p,q: delta>=0 la TH2 - hai so thuc doi nhau, delta<0 la
    TH3 - cap phuc lien hop). Xem python/algo/luy_thua_xuong_thang.md."""
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
        # 1 cap gia tri rieng doi nhau/phuc lien hop (xem Bug 1, bug_report.md),
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
        # cach 0 mot khoang > eps (vi du ~1e-3 lan lon hon) - dung eps truc
        # tiep se khong bat duoc dau hieu nay (xem Bug 1, bug_report.md).
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
                    # 2 nghiem thuc: chi la TH2 (lambda1=-lambda2) neu p=lambda1+lambda2~0.
                    # Neu p khac 0 ro ret, day chi la 2 gtr troi thuc KHONG dong module ma
                    # phuong phap tong quat nay tinh duoc som hon phep thu ti so don gian cua
                    # TH1 (thuong xay ra khi |lambda2/lambda1| khong qua nho) -> tra ve dung
                    # dang TH1 voi gtr co module lon hon, tranh gan nham "TH2".
                    # Dung guard_thresh (long hon eps, xem Bug 1 bug_report.md) thay vi eps
                    # truc tiep: p hoi tu ve 0 cham hon do on dinh cua chinh he p,q, nen voi
                    # nguong eps qua chat, TH2 that de bi tra ve nham TH1.
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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the ap dung PP luy thua!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    eps = float(input("Nhap sai so cho phep epsilon (vi du 1e-6): "))

    result = power_method(A, eps=eps)

    if result["case"] == 0:
        print(f"\nKhong hoi tu sau {result['iterations']} lan lap: co the co > 2 gtr trien cung module,")
        print("hoac |l1|,|l2| qua gan |l3| nen hoi tu rat cham.")
        return

    history = result["history"]
    show = sorted(set(range(min(2, len(history)))) | set(range(max(0, len(history) - 2), len(history))))
    for idx in show:
        k, lam1_est, x_k = history[idx]
        lam_str = f"{lam1_est:.6f}" if lam1_est is not None else "(chua uoc luong duoc)"
        print(f"\n--- Lan lap {k}: uoc luong lambda1 (TH1, ti so y1[i]/x[i]) = {lam_str} ---")
        print("x^(k) (vector lap, da chuan hoa) =")
        print_matrix(x_k)

    case = result["case"]
    it = result["iterations"]
    if case == 1:
        print(f"\n=== TH1: |lambda1| > |lambda2| (hoi tu sau {it} lan lap) ===")
    elif case == 2:
        print(f"\n=== TH2: |lambda1| = |lambda2| > |lambda3|, lambda1 = -lambda2 (hoi tu sau {it} lan lap) ===")
    else:
        print(f"\n=== TH3: |lambda1| = |lambda2| > |lambda3|, lambda1 = conj(lambda2) (hoi tu sau {it} lan lap) ===")

    for lam, v in zip(result["eigenvalues"], result["eigenvectors"]):
        if isinstance(lam, complex) and abs(lam.imag) > 1e-9:
            print(f"lambda = {lam.real:.6f} {'+' if lam.imag >= 0 else '-'} {abs(lam.imag):.6f}i")
        else:
            lam_r = lam.real if isinstance(lam, complex) else lam
            print(f"lambda = {lam_r:.6f}")
        print_matrix(v)

    print("--- (Doi chieu) Gia tri rieng tinh boi numpy.linalg.eigvals, sap theo |lambda| giam dan ---")
    eig_np = np.linalg.eigvals(A)
    order = np.argsort(-np.abs(eig_np))
    print_matrix(eig_np[order].reshape(1, -1))


if __name__ == "__main__":
    main()
