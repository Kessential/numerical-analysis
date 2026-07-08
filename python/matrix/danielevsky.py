import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:10.4f}", end=" ")
        print()
    print()


def print_poly(coeffs):
    deg = len(coeffs) - 1
    parts = []
    for i, c in enumerate(coeffs):
        power = deg - i
        c = 0.0 if abs(c) < 1e-9 else c
        if c == 0:
            continue
        sign = "-" if c < 0 else ("+" if parts else "")
        mag = abs(c)
        if power == 0:
            term = f"{mag:.4f}"
        elif power == 1:
            term = f"{mag:.4f}*lambda"
        else:
            term = f"{mag:.4f}*lambda^{power}"
        parts.append(f"{sign} {term}".strip())
    print(" ".join(parts) if parts else "0")


def swap_matrix(n, i, j):
    C = np.eye(n)
    C[[i, j], :] = C[[j, i], :]
    return C


def build_M(A, k, m, n):
    row = A[k, :].copy()
    pivot = row[m]
    M = np.eye(n)
    M[m, :] = row
    Minv = np.eye(n)
    Minv[m, :] = -row / pivot
    Minv[m, m] = 1.0 / pivot
    return M, Minv


def solve_coupling(Ak, F, B):
    """Giai K (k x m) tu he Sylvester Ak.K - K.F = B, voi F la khoi Frobenius
    (dong 0 chua he so -p_i, duoi duong cheo phu la 1). Dung de day B ve 0
    (buoc S_q trong slide) ma khong lam thay doi F va Ak."""
    k = Ak.shape[0]
    m = F.shape[0]
    p = -F[0, :]
    I_k = np.eye(k)

    Rs = [I_k.copy()]
    ds = [np.zeros(k)]
    for j in range(m - 1):
        Rs.append(Ak @ Rs[-1] + p[j] * I_k)
        ds.append(Ak @ ds[-1] + B[:, j])

    lhs = Ak @ Rs[-1] + p[m - 1] * I_k
    rhs = Ak @ ds[-1] + B[:, m - 1]
    K0 = np.linalg.solve(lhs, rhs)

    K = np.zeros((k, m))
    for j in range(m):
        K[:, j] = Rs[j] @ K0 - ds[j]
    return K


def danielevsky_reduce(A):
    """Bien doi dong dang A -> F = P^-1 A P, F dang khoi (block-diagonal) cac
    khoi Frobenius. Tra ve (F, P, blocks) voi blocks = [(start, size), ...]
    (chi so 0-based) mo ta vi tri/co cua tung khoi Frobenius tren duong cheo."""
    n = A.shape[0]
    A = A.astype(float).copy()
    P = np.eye(n)
    blocks = []

    r = n
    while r > 1:
        k = r - 1
        split_at = None
        while k >= 1:
            m = k - 1
            if abs(A[k, m]) < 1e-9:
                j_found = None
                for j in range(m):
                    if abs(A[k, j]) > 1e-9:
                        j_found = j
                        break
                if j_found is not None:
                    print(f"Hang {k + 1}: a[{k + 1}][{m + 1}] = 0 -> hoan doi hang/cot {j_found + 1} <-> {m + 1}.")
                    C = swap_matrix(n, j_found, m)
                    A = C @ A @ C
                    P = P @ C
                else:
                    print(f"Hang {k + 1}: toan bo a[{k + 1}][1..{m + 1}] = 0 -> tach khoi Frobenius co {r - k}, bat dau tu hang {k + 1}.")
                    split_at = k
                    break
            pivot = A[k, m]
            print(f"Hang {k + 1}: pivot a[{k + 1}][{m + 1}] = {pivot:.4f} != 0 -> dung M dua hang {k + 1} ve don vi tai cot {m + 1}.")
            M, Minv = build_M(A, k, m, n)
            A = M @ A @ Minv
            P = P @ Minv
            k -= 1

        if split_at is not None:
            k = split_at
            size = r - k
            if k > 0:
                Ak = A[:k, :k]
                F = A[k:r, k:r]
                Bc = A[:k, k:r]
                if np.max(np.abs(Bc)) > 1e-9:
                    print(f"Triet tieu khoi lien ket B ({k}x{size}) giua phan 1..{k} va khoi Frobenius {k + 1}..{r}:")
                    print_matrix(Bc)
                    try:
                        K = solve_coupling(Ak, F, Bc)
                    except np.linalg.LinAlgError:
                        print(
                            f"CANH BAO: he Sylvester Ak.K - K.F = B suy bien (phan 1..{k} va khoi Frobenius "
                            f"{k + 1}..{r} co chung tri rieng - ma tran co tri rieng boi/derogatory). "
                            "Khong triet tieu duoc B, bo qua buoc nay. Da thuc dac trung va tri rieng cua "
                            "cac khoi van dung, nhung vecto rieng ung voi khoi nay co the KHONG chinh xac."
                        )
                    else:
                        S = np.eye(n)
                        S[:k, k:r] = K
                        Sinv = np.eye(n)
                        Sinv[:k, k:r] = -K
                        A = S @ A @ Sinv
                        P = P @ Sinv
                        print("-> Sau khi triet tieu, B =")
                        print_matrix(A[:k, k:r])
            blocks.append((k, size))
            r = k
        else:
            blocks.append((0, r))
            r = 0

    if r == 1:
        blocks.append((0, 1))

    return A, P, blocks


def block_polynomial(F, start, size):
    row = F[start, start : start + size]
    return np.concatenate(([1.0], -row))


def eigen_from_blocks(F, P, blocks):
    n = F.shape[0]
    results = []
    for start, size in blocks:
        coeffs = block_polynomial(F, start, size)
        print(f"--- Khoi Frobenius o vi tri {start + 1}, co {size}: da thuc dac trung ---")
        print_poly(coeffs)

        roots = np.roots(coeffs) if size > 1 else np.array([-coeffs[1]])
        for lam in roots:
            v_local = np.array([lam ** (size - 1 - i) for i in range(size)], dtype=complex)
            v_full = np.zeros(n, dtype=complex)
            v_full[start : start + size] = v_local
            v = P @ v_full
            scale = v[np.argmax(np.abs(v))]
            if abs(scale) > 1e-12:
                v = v / scale
            results.append((lam, v))
    return results


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
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {n} x {cols}), khong the ap dung PP Danielevsky!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    F, P, blocks = danielevsky_reduce(A)

    print("--- DANG CHUAN FROBENIUS F = P^-1.A.P ---")
    print_matrix(F)
    print("--- MA TRAN CHUYEN CO SO P ---")
    print_matrix(P)

    results = eigen_from_blocks(F, P, blocks)

    print("=== GIA TRI RIENG VA VECTO RIENG ===")
    for lam, v in results:
        if abs(lam.imag) < 1e-6:
            print(f"lambda = {lam.real:.4f}")
        else:
            sign = "+" if lam.imag >= 0 else "-"
            print(f"lambda = {lam.real:.4f} {sign} {abs(lam.imag):.4f}i")
        v_show = v.real if np.all(np.abs(v.imag) < 1e-6) else v
        print_matrix(v_show)

    print("--- (Doi chieu) Gia tri rieng tinh boi numpy.linalg.eigvals ---")
    eig_np = np.linalg.eigvals(A)
    order = np.argsort(-eig_np.real)
    print_matrix(eig_np[order].reshape(1, -1))


if __name__ == "__main__":
    main()
