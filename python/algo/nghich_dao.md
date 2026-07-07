# Các phương pháp tìm ma trận nghịch đảo

## 1. PP dùng ma trận phụ hợp

**Input:** ma trận vuông $A_{n \times n}$.

**Output:** $A^{-1}$, hoặc thông báo nếu $\det A = 0$ (A suy biến, không kha nghịch).

**Bước 1 — Kiểm tra đầu vào:** A phải vuông. Nếu không, dừng và báo lỗi.

**Bước 2 — Tính định thức bằng khai triển Laplace đệ quy theo hàng 1:** với minor $A_{ij}$ là ma trận có được từ A bằng cách xóa hàng $i$, cột $j$:

$$
\det M = \begin{cases}
M_{11} & \text{nếu } M \text{ cỡ } 1\times1 \\
M_{11}M_{22}-M_{12}M_{21} & \text{nếu } M \text{ cỡ } 2\times2 \\
\displaystyle\sum_{j=1}^{k} (-1)^{1+j} M_{1j} \det(\text{minor}(M,1,j)) & \text{nếu } M \text{ cỡ } k\times k,\ k>2
\end{cases}
$$

Áp dụng cho $M = A$ để được $\det A$.

**Bước 3 — Kiểm tra suy biến:** nếu $|\det A| \approx 0$: dừng, báo A không khả nghịch.

**Bước 4 — Tính ma trận phụ hợp** $C = [c_{ij}]_{n\times n}$:

$$
c_{ij} = (-1)^{i+j} \det(\text{minor}(A,i,j)), \quad i,j = \overline{1,n}
$$

(mỗi $c_{ij}$ dùng lại công thức đệ quy ở Bước 2 cho minor cỡ $(n-1)\times(n-1)$.)

**Bước 5 — Kết luận:**

$$
A^{-1} = \frac{1}{\det A} C^{T}
$$

**Lưu ý:** khai triển Laplace đệ quy có độ phức tạp $O(n!)$ nên phương pháp này chỉ nên dùng minh họa với $n$ nhỏ (khoảng $\le 6$), không phù hợp tính toán số lớn.

---

## 2. PP Gauss và Gauss-Jordan

Slide gộp chung PP Gauss và PP Gauss-Jordan vì cả hai đều dựa trên cùng một công thức $[A\mid E] \Leftrightarrow [E\mid A^{-1}]$. Tuy nhiên đây là hai thuật toán khử khác nhau (giống cặp file `gauss.py`/`gauss_jordan.py` đã có để giải $Ax=b$), nên cài đặt thành 2 file riêng: `nghich_dao_gauss.py` (khử xuôi + thế ngược) và `nghich_dao_gauss_jordan.py` (khử toàn phần).

### 2a. PP Gauss (khử xuôi + thế ngược)

**Input:** ma trận vuông $A_{n \times n}$.

**Output:** $A^{-1}$, hoặc thông báo nếu A suy biến (không tìm được pivot khác 0 ở cột nào đó).

**Bước 1 — Kiểm tra đầu vào:** A phải vuông. Nếu không, dừng và báo lỗi.

**Bước 2 — Ghép ma trận mở rộng** $[A \mid E_n]$.

**Bước 3 — Quy trình thuận (khử Gauss, đổi hàng khi cần):** với $i = \overline{1,n}$:

* Nếu $a_{ii} = 0$: tìm hàng $t > i$ có $a_{ti} \ne 0$ để đổi chỗ. Nếu không có: dừng, báo **A suy biến, không khả nghịch**.
* Khử các phần tử bên dưới pivot: $\text{Aug}[t,\,i{:}] \mathrel{-}= \dfrac{a_{ti}}{a_{ii}}\,\text{Aug}[i,\,i{:}], \quad \forall t>i$.

Sau bước này, phần $A$ của $\text{Aug}$ trở thành ma trận tam giác trên $U$.

**Bước 4 — Quy trình nghịch (thế ngược từng cột của $E$):** với $r = \overline{n,1}$ (giảm dần):

$$
X_{r,:} = \frac{\text{Aug}[r,\, n{+}1{:}2n] - \sum_{k>r} u_{rk}X_{k,:}}{u_{rr}}
$$

Kết quả $X = A^{-1}$.

### 2b. PP Gauss-Jordan (khử toàn phần)

**Input:** ma trận vuông $A_{n \times n}$.

**Output:** $A^{-1}$, hoặc thông báo nếu không tìm được đủ $n$ phần tử khóa (A suy biến).

**Bước 1 — Kiểm tra đầu vào:** A phải vuông. Nếu không, dừng và báo lỗi.

**Bước 2 — Ghép ma trận mở rộng** $[A \mid E_n]$.

**Bước 3 — Khử toàn diện (Gauss-Jordan, có chọn khóa từng bước):** lặp $n$ lần, mỗi lần chọn 1 phần tử khóa (pivot) chưa dùng theo hàng/cột (ưu tiên $|a_{rc}|=1$, nếu không có thì chọn $|a_{rc}|$ lớn nhất trong các hàng/cột còn lại):

* Nếu không còn phần tử khóa nào khác 0: dừng, báo **A suy biến, không khả nghịch**.
* Chuẩn hóa hàng khóa $p$: $\text{Aug}[p,:] \mathrel{/}= \text{Aug}[p,q]$.
* Khử tất cả các hàng khác (cả trên lẫn dưới) theo cột khóa $q$: $\text{Aug}[r,:] \mathrel{-}= \text{Aug}[r,q]\cdot\text{Aug}[p,:] \; \forall r \ne p$.
* Ghi nhớ $\text{ind}[p] = q$ (hàng $p$ tương ứng ẩn/cột thứ $q$).

**Bước 4 — Kết luận:** vì A vuông và khả nghịch nên sau $n$ bước, mỗi hàng $r$ đã khử xong ứng với $\text{ind}[r]$ là một hoán vị của $\{1,\dots,n\}$. Xếp lại các hàng theo đúng thứ tự ẩn để được nghiệm:

$$
A^{-1}[\text{ind}[r], :] = \text{Aug}[r,\; n{+}1 : 2n], \quad r = \overline{1,n}
$$

(đây chính là công thức $[A\mid E] \Leftrightarrow [E\mid A^{-1}]$ của slide, chỉ khác là có sắp xếp lại hàng do dùng khử toàn phần/full pivoting để ổn định số học.)

---

## 3. PP Choleski, Jacobi, Gauss-Seidel tìm ma trận nghịch đảo

Ba phương pháp này **không có thuật toán mới**: chúng chính là thuật toán giải $AX=B$ đã có ở `lu_cholesky.md` (mục 2), `lapdon_jacobi.md` (mục 2) và `gauss_seidel.md`, chỉ thay $B := E_n$ (ma trận đơn vị) thay vì đọc từ file. Code tương ứng (`nghich_dao_cholesky.py`, `nghich_dao_jacobi.py`, `nghich_dao_gauss_seidel.py`) chỉ khác bản gốc ở chỗ tự dựng $E_n$ và không hỏi số cột của B.

---

## 4. PP viền quanh (bordering method)

**Input:** ma trận vuông $A_{n \times n}$ khả nghịch.

**Output:** $A^{-1}$, hoặc thông báo A không khả nghịch.

**Bước 0 — Kiểm tra điều kiện và chọn $M$:** 

* Nếu **không phải** mọi $A_k$ ($k=\overline{1,n}$) đều khả nghịch: đặt $M := A^TA$ và cờ $\text{assym} := \text{true}$.
* Trái lại: đặt $M := A$ và cờ $\text{assym} := \text{false}$.

**Bước 1 — Khởi tạo** ($k=1$): Nếu $m_{11}=0$, thông báo A suy biến và dừng thuật toán. Trái lại, gán $M_1^{-1}=[1/m_{11}]$.

**Bước 2 — Với $k=2,\dots,n$, viền thêm 1 hàng + 1 cột:** đặt $\alpha_{k-1,1}=M[1{:}k{-}1,k]$, $\alpha_{1,k-1}=M[k,1{:}k{-}1]$, $m_{kk}=M[k,k]$. Tính trước

$$
s = m_{kk} - \alpha_{1,k-1}M_{k-1}^{-1}\alpha_{k-1,1}
$$

Nếu $s\approx0$: thông báo và dừng hẳn. Ngược lại:

$$
b_{kk} = \frac{1}{s}
$$

$$
\beta_{k-1,1} = -b_{kk}\,M_{k-1}^{-1}\alpha_{k-1,1}
$$

$$
\beta_{1,k-1} = -b_{kk}\,\alpha_{1,k-1}M_{k-1}^{-1}
$$

$$
B_{k-1} = M_{k-1}^{-1} + b_{kk}\left(M_{k-1}^{-1}\alpha_{k-1,1}\right)\left(\alpha_{1,k-1}M_{k-1}^{-1}\right)
$$

$$
M_k^{-1} = \begin{bmatrix} B_{k-1} & \beta_{k-1,1} \\ \beta_{1,k-1} & b_{kk} \end{bmatrix}
$$

**Bước 3 — Kết luận:** lặp Bước 2 đến $k=n$ được $M_n^{-1}=M^{-1}$. Nếu $\text{assym}=\text{true}$: $A^{-1}=M^{-1}A^T$; trái lại $A^{-1}=M^{-1}$.

---

## 5. PP lặp Newton

**Input:** ma trận vuông $A_{n\times n}$ khả nghịch, sai số $\varepsilon>0$, số lần lặp tối đa $N$.

**Output:** $A^{-1}$ xấp xỉ, hoặc cảnh báo nếu chưa đạt sai số yêu cầu sau $N$ lần lặp.

**Bước 1 — Chọn giá trị lặp ban đầu:**

$$
X_0 = \frac{A^{T}}{\|A\|_1 \|A\|_\infty}, \qquad G_0 = E - AX_0
$$

Cách chọn này đảm bảo $q_2 := \|G_0\|_2 < 1$ (chuẩn phổ) với mọi $A$ khả nghịch — điều kiện hội tụ **lý thuyết** của phương pháp.

**Bước 2 — Ước lượng thực hành:** tính $\|.\|_2$ tốn kém (cần SVD), nên dùng $q := \|G_0\|_\infty$ làm ước lượng thay thế, không phải điều kiện chặt như $q_2$ (có thể $q\ge1$ dù $q_2<1$ vẫn hội tụ, hoặc ngược lại). Nếu $q\ge1$: in cảnh báo "ước lượng $\|.\|_\infty$ không đảm bảo hội tụ" nhưng vẫn tiếp tục lặp tối đa $N$ lần.

**Bước 3 — Lặp:** với $n=1,2,\dots$

$$
X_n = X_{n-1}(2E - AX_{n-1}), \qquad G_n = E - AX_n
$$

Sai số tiên nghiệm theo $q_2$: $\|X_n-A^{-1}\|_2 \le \dfrac{q_2^{2^n}}{1-q_2}\|X_0\|_2$ (chỉ tham khảo, không dùng để dừng vòng lặp vì $q_2$ tốn kém tính). Tiêu chí dừng thực hành: nếu $\|G_n\|_\infty < \varepsilon$ hoặc $n=N$: dừng, trả về $X_n$.

**Bước 4 — Cảnh báo nếu không hội tụ:** nếu dừng do $n=N$ mà $\|G_N\|_\infty \ge \varepsilon$: cảnh báo rõ kết quả $X_N$ không đáng tin cậy (A gần suy biến/điều kiện xấu, hoặc cần tăng $N$), thay vì âm thầm trả về.
