## PP lặp Gauss-Seidel

**Input:**

* Ma trận vuông $A_{m \times m}$, $B_{m \times p}$.
* Sai số cho phép $\varepsilon > 0$, số lần lặp tối đa $N$.

**Output:** nghiệm xấp xỉ $X$, hoặc thông báo nếu $A$ không chéo trội (không đảm bảo hội tụ) hay có $a_{ii} = 0$ (không lặp được).

**Bước 1 — Kiểm tra đầu vào:** Nếu $a_{ii} = 0, \ i = \overline{1, m}$: dừng và báo lỗi.

**Bước 2 — Phân rã** $A = D_A - L_A - U_A$, với 
* $D_A = \text{diag}(a_{11}, \dots, a_{mm})$ 
* $(L_A)_{ij} = -a_{ij}$ với $j < i$ là phần dưới đường chéo
* $(U_A)_{ij} = -a_{ij}$ với $j > i$ là phần trên đường chéo. Khi đó $D_A - L_A$ là phần tam giác dưới của $A$.

**Bước 3 — Kiểm tra tính chéo trội của $A$, xác định chuẩn $\|\cdot\|$ và các hệ số $s, q$ dùng để đánh giá sai số:**

* **Chéo trội hàng** nếu với mọi $i$: $\displaystyle |a_{ii}| > \sum_{j \ne i} |a_{ij}|$. Dùng chuẩn $\|\cdot\|_\infty$ và

$$
s = 0, \qquad q = \max_{i=\overline{1,m}} \frac{\sum_{j>i} |a_{ij}|}{|a_{ii}| - \sum_{j<i} |a_{ij}|}
$$

* **Chéo trội cột** nếu với mọi $j$: $\displaystyle |a_{jj}| > \sum_{i \ne j} |a_{ij}|$. Dùng chuẩn $\|\cdot\|_1$ và

$$
s = \max_{j=\overline{1,m}} \sum_{i>j} \left| \frac{a_{ij}}{a_{jj}} \right|, \qquad q = \max_{j=\overline{1,m}} \frac{\sum_{i<j} |a_{ij}|}{|a_{jj}| - \sum_{i>j} |a_{ij}|}
$$

* Nếu không thỏa cả hai: in cảnh báo "không đảm bảo hội tụ", coi như chéo trội hàng (dùng chuẩn $\|\cdot\|_\infty$, $s=0$, $q$ theo công thức chéo trội hàng).

**Bước 4 — Kiểm tra điều kiện hội tụ:** Nếu $q \ge 1$: in cảnh báo "không đảm bảo hội tụ" nhưng vẫn tiếp tục lặp tối đa $N$ lần.

**Bước 5 — Lặp:** chọn $X_0 = 0$ (hoặc giá trị ban đầu bất kỳ); với $n = 1, 2, \dots$, giải hệ tam giác dưới

$$
(D_A - L_A) X_n = U_A X_{n-1} + B
$$

thế xuôi theo từng dòng $i = \overline{1, m}$:

$$
x_i^{(n)} = \frac{b_i - \sum_{j<i} a_{ij} x_j^{(n)} - \sum_{j>i} a_{ij} x_j^{(n-1)}}{a_{ii}}
$$

Tính sai số hậu nghiệm (chuẩn $\|\cdot\|$ đã chọn ở Bước 3):

$$
e_n = \frac{q}{(1-s)(1-q)} \|X_n - X_{n-1}\|
$$

Nếu $e_n < \varepsilon$ hoặc $n = N$: dừng, trả về $X_n$ là nghiệm xấp xỉ của $Ax = B$. Ngược lại tăng $n$, quay lại Bước 5.
