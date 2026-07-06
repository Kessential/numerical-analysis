# PP lặp đơn và PP lặp Jacobi giải Ax = B

## 1. PP lặp đơn (giải $x = Cx + D$ bằng lặp)

**Input:**

* $C_{m \times m}$, $D_{m \times p}$, sao cho $x = Cx + D$.
* Giá trị ban đầu $x_0 \in \mathbb{R}^{m \times p}$ (mặc định $x_0 = D$).
* Sai số cho phép $\varepsilon > 0$, số lần lặp tối đa $N$.

**Output:** nghiệm xấp xỉ $x_n$ có sai số hậu nghiệm nhỏ hơn $\varepsilon$, hoặc cảnh báo nếu $\|C\| \ge 1$ (không đảm bảo hội tụ).

**Bước 1 — Kiểm tra điều kiện hội tụ:** tính $q = \|C\|$ theo một chuẩn ma trận cố định (ví dụ $\|\cdot\|_\infty$ — max tổng trị tuyệt đối theo hàng). Nếu $q \ge 1$: in cảnh báo "không đảm bảo hội tụ".

**Bước 2 — (Tùy chọn) Ước lượng trước số lần lặp cần thiết:** với $x_1 = Cx_0 + D$, công thức sai số tiên nghiệm

$$
\|x_n - x^*\| \le \frac{q^n}{1-q}\|x_1-x_0\| < \varepsilon \;\Rightarrow\; n \ge \frac{\ln\big(\varepsilon(1-q)/\|x_1-x_0\|\big)}{\ln q}
$$

Chỉ mang tính tham khảo, có thể bỏ qua nếu dừng lặp trực tiếp theo sai số hậu nghiệm ở Bước 3.

**Bước 3 — Lặp:** với $n = 1, 2, \dots$

$$
x_n = Cx_{n-1} + D
$$

tính sai số hậu nghiệm

$$
e_n = \frac{q}{1-q}\|x_n - x_{n-1}\|
$$

($\|\cdot\|$ ở đây là chuẩn đã chọn ở Bước 1.) Nếu $e_n < \varepsilon$ hoặc $n = N$: dừng, trả về $x_n$ là nghiệm xấp xỉ. Ngược lại tăng $n$, quay lại Bước 3.

---

## 2. PP lặp Jacobi (giải $Ax = B$)

**Input:**

* Ma trận vuông $A_{m \times m}$, $B_{m \times p}$.
* Sai số cho phép $\varepsilon > 0$, số lần lặp tối đa $N$.

**Output:** nghiệm xấp xỉ $X$, hoặc thông báo nếu $A$ không chéo trội (không đảm bảo hội tụ) hay có $a_{ii} = 0$ (không dựng được $T$).

**Bước 1 — Kiểm tra đầu vào:** Nếu $a_{ii} = 0, \ i = \overline{1, m}$: dừng và báo lỗi.

**Bước 2 — Kiểm tra tính chéo trội của $A$, xác định chuẩn $\|\cdot\|$ và hệ số $\lambda$ dùng để đánh giá sai số:**

* **Chéo trội hàng** nếu với mọi $i$: $\displaystyle |a_{ii}| > \sum_{j \ne i} |a_{ij}|$. Sử dụng chuẩn $\|\cdot\|_\infty$ và $\lambda = 1$.
* **Chéo trội cột** nếu với mọi $i$: $\displaystyle |a_{ii}| > \sum_{j \ne i} |a_{ji}|$. Sử dụng chuẩn $\|\cdot\|_1$ và $\displaystyle \lambda = \frac{\max_i |a_{ii}|}{\min_i |a_{ii}|}$.
* Nếu không thỏa cả hai: in cảnh báo "không đảm bảo hội tụ", coi như chéo trội hàng (dùng $C$, chuẩn $\|\cdot\|_\infty$, $\lambda=1$).

**Bước 3 — Kiểm tra điều kiện hội tụ:** đặt $T = \text{diag}(1/a_{11}, \dots, 1/a_{mm})$, tính

$$
q = \begin{cases}
    \|I - TA\|_\infty & \text{nếu $A$ chéo trội hàng (hoặc mặc định)} \\
    \|I - AT\|_1 & \text{nếu $A$ chéo trội cột}
\end{cases}
$$

Nếu $q \ge 1$: in cảnh báo "không đảm bảo hội tụ" nhưng vẫn tiếp tục lặp tối đa $N$ lần.

**Bước 4 — Xây dựng ma trận lặp:**

$$
C = I - TA, \qquad D = TB
$$

Dùng $C$ này để lặp cho **cả hai** trường hợp (vì $C=I-TA$ và $C_1=I-AT$ đồng dạng qua $T$: $TC_1 = CT$, nên lặp trên $x$ bằng $C$ tương đương lặp trên $y=T^{-1}x$ bằng $C_1$). Ở trường hợp chéo trội cột, $C_1=I-AT$ chỉ dùng để tính $q$ ở Bước 3 (vì $\|C\|_1$ nói chung không nhỏ hơn 1) — **không** dùng $C_1$ để lặp trực tiếp trên $x$.

**Bước 5 — Lặp:** chọn $X_0 = D$ (tương ứng $X_{-1} = 0$); với $n = 1, 2, \dots$

$$
X_n = CX_{n-1} + D
$$

tính sai số hậu nghiệm

$$
e_n = \lambda \frac{q}{1-q}\|X_n - X_{n-1}\|
$$

($\|\cdot\|$ lấy theo chuẩn đã chọn ở Bước 2.) Nếu $e_n < \varepsilon$ hoặc $n = N$: dừng, trả về $X_n$ là nghiệm xấp xỉ của $Ax = B$. Ngược lại tăng $n$, quay lại Bước 5.