# Thuật toán PP phân tách LU và PP Choleski giải Ax = b

Tổng hợp lại từ slide `lu-cholesky.md` thành thuật toán chi tiết để cài đặt (code tại `lu.py` và `cholesky.py`).

## 1. PP phân tách LU (Doolittle, chọn $u_{ii} = 1$)

**Input:**

* Ma trận vuông $A$ cỡ $n \times n$.
* Ma trận $B$ cỡ $n \times p$ (vế phải, có thể nhiều cột).

**Output:**

* Ma trận $L$ (tam giác dưới), $U$ (tam giác trên, đường chéo bằng 1) sao cho $A = LU$.
* Ma trận nghiệm $X$ cỡ $n \times p$ thỏa $AX = B$, hoặc thông báo A suy biến nếu không phân tách được.

**Bước 1 — Kiểm tra đầu vào:** A phải vuông ($n \times n$). Nếu không, dừng và báo lỗi.

**Bước 2 — Phân tách $A = LU$:** với $t = \overline{1, n}$ (lần lượt từng cột):

* Bước 2.1: tính cột $t$ của $L$:

$$
l_{it} = a_{it} - \sum_{k=1}^{t-1} l_{ik} u_{kt}, \quad i = \overline{t, n}
$$

* Bước 2.2: nếu $l_{tt} = 0$ thì **A suy biến** (hoặc cần đổi hàng) — thông báo và dừng thuật toán.
* Bước 2.3: đặt $u_{tt} = 1$ và tính hàng $t$ của $U$:

$$
u_{tk} = \frac{a_{tk} - \sum_{j=1}^{t-1} l_{tj} u_{jk}}{l_{tt}}, \quad k = \overline{t+1, n}
$$

**Bước 3 — Giải hệ** $AX = B \Leftrightarrow LUX = B$: với mỗi cột thứ $j$ $(j = \overline{1, p})$ của $B$, đặt $b^{(j)}_i := B_{ij}$ ($i = \overline{1,n}$, tức vector cột $j$ của B), giải hệ

$$
\begin{cases} Ly^{(j)} = b^{(j)} \\ Ux^{(j)} = y^{(j)} \end{cases}
$$

để được nghiệm $x^{(j)}$, rồi gán $X_{ij} := x^{(j)}_i$ ($i = \overline{1,n}$) là cột thứ $j$ của ma trận nghiệm X. Cụ thể mỗi lần giải gồm:

* Bước 3.1 (thế xuôi, vì $L$ tam giác dưới):

$$
y^{(j)}_i = \frac{b^{(j)}_i - \sum_{k=1}^{i-1} l_{ik} y^{(j)}_k}{l_{ii}}, \quad i = \overline{1, n}
$$

* Bước 3.2 (thế ngược, vì $U$ tam giác trên, đường chéo bằng 1):

$$
x^{(j)}_i = y^{(j)}_i - \sum_{k=i+1}^{n} u_{ik} x^{(j)}_k, \quad i = \overline{n, 1}
$$

Lặp Bước 3.1–3.2 cho $j = \overline{1, p}$ để thu được toàn bộ ma trận nghiệm $X$.

---

## 2. PP phân tách Choleski (chọn $U = L^T \Rightarrow A = U^T U$)

**Input:**

* Ma trận vuông $A$ cỡ $n \times n$.
* Ma trận $B$ cỡ $n \times p$ (vế phải, có thể nhiều cột).

**Output:**

* Ma trận $U$ (tam giác trên) sao cho $M = U^T U$ (M định nghĩa ở Bước 2).
* Ma trận nghiệm $X$ cỡ $n \times p$ thỏa $AX = B$, hoặc thông báo suy biến nếu không thỏa điều kiện PP Choleski.

**Bước 1 — Kiểm tra đầu vào:** A phải vuông ($n \times n$). Nếu không, dừng và báo lỗi.

**Bước 2 — Đưa về hệ đối xứng:**

* Nếu A đối xứng: đặt $M := A$, $D := B$.
* Nếu A không đối xứng: đặt $M := A^T A$, $D := A^T B$ (tương đương $AX=B \Leftrightarrow A^TAX = A^TB$). Không cần kiểm tra riêng $\det A \ne 0$ vì A suy biến sẽ khiến $u_{ii} = 0$ ở Bước 3.

**Bước 3 — Phân tách $M = U^T U$ (U tam giác trên):** với $i = \overline{1, n}$:

* Bước 3.1: tính

$$
s = m_{ii} - \sum_{j=1}^{i-1} u_{ji}^2
$$

* Bước 3.2: xét các trường hợp đặc biệt của $s$:
  * Nếu $s = 0$: **M suy biến**, thông báo và dừng — không thỏa điều kiện PP Choleski.
  * Nếu $s < 0$: về lý thuyết cần chuyển sang tính trên tập số phức $\mathbb{C}$ ($u_{ii} = \sqrt{s}$ là căn của số âm). Thuật toán (cài đặt) này bỏ qua trường hợp này: coi như **M không thỏa điều kiện PP Choleski**, thông báo và dừng.
  * Ngược lại ($s > 0$): $u_{ii} = \sqrt{s}$ (số thực dương).
* Bước 3.3: tính hàng $i$ còn lại của U:

$$
u_{ik} = \frac{m_{ik} - \sum_{j=1}^{i-1} u_{ji} u_{jk}}{u_{ii}}, \quad k = \overline{i+1, n}
$$

**Bước 4 — Giải hệ** $MX = D \Leftrightarrow U^T U X = D$: với mỗi cột thứ $j$ $(j = \overline{1, p})$ của $D$, đặt $d^{(j)}_i := D_{ij}$ ($i = \overline{1,n}$, tức vector cột $j$ của D), giải hệ

$$
\begin{cases} U^T y^{(j)} = d^{(j)} \\ Ux^{(j)} = y^{(j)} \end{cases}
$$

để được nghiệm $x^{(j)}$, rồi gán $X_{ij} := x^{(j)}_i$ ($i = \overline{1,n}$) là cột thứ $j$ của ma trận nghiệm X. Cụ thể mỗi lần giải gồm:

* Bước 4.1 (thế xuôi, vì $U^T$ tam giác dưới):

$$
y^{(j)}_i = \frac{d^{(j)}_i - \sum_{k=1}^{i-1} u_{ki} y^{(j)}_k}{u_{ii}}, \quad i = \overline{1, n}
$$

* Bước 4.2 (thế ngược, vì $U$ tam giác trên):

$$
x^{(j)}_i = \frac{y^{(j)}_i - \sum_{k=i+1}^{n} u_{ik} x^{(j)}_k}{u_{ii}}, \quad i = \overline{n, 1}
$$

Lặp Bước 4.1–4.2 cho $j = \overline{1, p}$ để thu được toàn bộ ma trận nghiệm $X$.
