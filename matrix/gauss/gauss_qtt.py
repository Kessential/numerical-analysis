import numpy as np


def print_matrix(mat):
    for row in mat:
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def quy_trinh_thuan(A):
    m, n = A.shape
    i, j = 0, 0
    ind = [-1] * m
    step = 1

    print("\n--- LICH SU CHON KHOA (PIVOT) ---")

    while i < m and j < n:
        if abs(A[i, j]) > 1e-9:
            ind[i] = j

            print(f"Lan lap {step}: Chon pivot a[{i + 1}][{j + 1}] = {A[i, j]:.8f}")

            # Khu cac phan tu ben duoi pivot
            factors = A[i + 1 :, j] / A[i, j]
            A[i + 1 :, j:] -= np.outer(factors, A[i, j:])

            if step == 1:
                print("\n--- MA TRAN SAU LAN LAP 1 ---")
                print_matrix(A)

            i += 1
            j += 1
            step += 1
        else:
            found_nonzero = False
            for t in range(i + 1, m):
                if abs(A[t, j]) > 1e-9:
                    print(f"[!] Phat hien a[{i + 1}][{j + 1}] = 0. Doi cho hang {i + 1} va hang {t + 1}")
                    A[[i, t]] = A[[t, i]]
                    found_nonzero = True
                    break

            if not found_nonzero:
                j += 1

    print("\n--- MANG LUU VI TRI PIVOT ---")
    print("ind = [", " ".join(str(p + 1 if p != -1 else 0) for p in ind), "]")

    print("\n--- MA TRAN BAC THANG (SAU QUY TRINH THUAN) ---")
    print_matrix(A)


def main():
    filename = "test.txt"
    try:
        A = np.loadtxt(filename, dtype=float, ndmin=2)
    except OSError:
        print(f"Loi: Khong the mo duoc file: {filename}")
        return

    if A.size == 0:
        print("Loi: File rong!")
        return

    m, n = A.shape
    print(f"=> Ma tran doc duoc co kich thuoc {m} x {n}.")

    quy_trinh_thuan(A)


if __name__ == "__main__":
    main()
