import numpy as np


def print_matrix(mat):
    for row in np.atleast_2d(mat):
        for val in row:
            val = 0.0 if abs(val) < 1e-9 else val
            print(f"{val:14.8f}", end=" ")
        print()
    print()


def minor(M, i, j):
    return np.delete(np.delete(M, i, axis=0), j, axis=1)


def determinant(M):
    n = M.shape[0]
    if n == 1:
        return M[0, 0]
    if n == 2:
        return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]

    det = 0.0
    for j in range(n):
        det += ((-1) ** j) * M[0, j] * determinant(minor(M, 0, j))
    return det


def cofactor_matrix(A):
    n = A.shape[0]
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            C[i, j] = ((-1) ** (i + j)) * determinant(minor(A, i, j))
        print(f"--- Da tinh xong hang {i + 1} cua ma tran phu hop C ---")
        print_matrix(C[i, :].reshape(1, -1))
    return C


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
    if m != n:
        print(f"Loi: Ma tran A phai vuong (dang doc duoc la {m} x {n}), khong the tim ma tran phu hop!")
        return

    print(f"--- MA TRAN A ({n} x {n}) ---")
    print_matrix(A)

    detA = determinant(A)
    print(f"det(A) = {detA:.8f}\n")

    if abs(detA) < 1e-9:
        print("Loi: det(A) = 0 -> Ma tran A khong kha nghich!")
        return

    print("--- MA TRAN PHU HOP C, c[i][j] = (-1)^(i+j) * det(minor(A,i,j)) ---")
    C = cofactor_matrix(A)
    print("--- MA TRAN PHU HOP C (day du) ---")
    print_matrix(C)

    inv_A = C.T / detA
    print("--- MA TRAN NGHICH DAO A^-1 = C^T / det(A) ---")
    print_matrix(inv_A)


if __name__ == "__main__":
    main()
