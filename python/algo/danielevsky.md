# Thuật toán Danielevsky: đa thức đặc trưng, giá trị riêng, véctơ riêng

Tổng hợp lại từ slide `danielevsky.md` thành thuật toán chi tiết để cài đặt (code tại `danielevsky.py`).

## 1. Đưa A về dạng khối Frobenius (tìm đa thức đặc trưng)

**Input:** ma trận vuông $A_{n\times n}$.

**Output:** ma trận $F=P^{-1}AP$ dạng khối chéo gồm các khối Frobenius $F_1^{(m_1)},\dots,F_s^{(m_s)}$ ($\sum m_i=n$), và ma trận chuyển cơ sở $P$.

**Bước 1 — Khởi tạo:** $A^{(1)} := A$, $P := I_n$, $r := n$ (cỡ của khối đang xét).

**Bước 2 — Trong khi $r>1$, với $k=r,r-1,\dots,2$** (đặt $m=k-1$ là cột đích của phần tử trục):

* Bước 2.1 — Nếu $a_{k,m}=0$ và tìm được $s<m$ với $a_{ks}\ne0$: hoán đổi đồng thời hàng $s\leftrightarrow m$ và cột $s\leftrightarrow m$ bằng ma trận hoán vị $C_{s,m}$: $A\leftarrow C_{s,m}AC_{s,m}$, $P\leftarrow PC_{s,m}$ (khi đó $a_{k,m}^{new}=a_{ks}\ne0$).
* Bước 2.2 — Nếu $a_{k,m}=0$ và mọi $j<m$ đều có $a_{kj}=0$: sang Bước 3, không thực hiện Bước 2.3 cho dòng $k$ này.
* Bước 2.3 — Dựng ma trận $M$: dòng $m$ của $M$ bằng dòng $k$ hiện tại của $A$, các dòng còn lại là dòng đơn vị:

$$
M_{m,j}=a_{kj},\quad j=\overline{1,n}
$$

  Nghịch đảo $M^{-1}$ cũng chỉ khác đơn vị ở dòng $m$:

$$
M^{-1}_{m,j}=\begin{cases}\dfrac{1}{a_{k,m}} & j=m \\[4pt] -\dfrac{a_{kj}}{a_{k,m}} & j\ne m\end{cases}
$$

  Cập nhật $A\leftarrow MAM^{-1}$, $P\leftarrow PM^{-1}$. Giảm $k$ đi 1, quay lại Bước 2. Nếu $k$ giảm xuống dưới $2$ (đã xử lý hết dòng $2,\dots,r$): ghi nhận khối $(vị\ trí=1,\ cỡ=r)$, đặt $r:=0$.

**Bước 3 — Tách khối tại dòng $k$:** ghi nhận khối Frobenius $(vị\ trí=k,\ cỡ=r-k+1)$ (dòng $k$ làm dòng hệ số, các dòng $k+1,\dots,r$ đã là đơn vị). Đặt $B$ là khối con của $A$ gồm các dòng $\overline{1,k-1}$, cột $\overline{k,r}$; $F$ là khối con của $A$ gồm các dòng và cột $\overline{k,r}$. Sang Bước 4, sau đó đặt $r:=k-1$ và quay lại Bước 2.

**Bước 4 — Triệt tiêu khối liên kết $B$** (cột của $B,F$ đánh số $1,\dots,m$ với $m=r-k+1$; gọi $b^{(q)}\in\mathbb R^{k-1}$ là cột $q$ của $B$, tức $b^{(q)}_i=b_{iq}$):

* Bước 4.1 — Với $q=1,\dots,m-1$ (đã có $q-1$ cột đầu của $B$ bằng 0): dựng $S_q$ = ma trận đơn vị, với

$$
(S_q)_{i,\,k+q} = -b^{(q)}_i,\qquad i=\overline{1,k-1}
$$

  (và $(S_q^{-1})_{i,\,k+q}=b^{(q)}_i$, các phần tử khác như ma trận đơn vị). Cập nhật $A\leftarrow S_qAS_q^{-1}$, $P\leftarrow PS_q^{-1}$. Sau bước này cột $q$ của $B$ bằng 0, cột $q+1$ trở thành $A_{k-1}b^{(q)}+b^{(q+1)}$.

* Bước 4.2 — Sau khi hết $q=1,\dots,m-1$: nếu cột $m$ của $B$ khác 0: hoán đổi hàng/cột $(k-1)\leftrightarrow r$ bằng $C_{k-1,r}$: $A\leftarrow C_{k-1,r}AC_{k-1,r}$, $P\leftarrow PC_{k-1,r}$, rồi quay lại Bước 2 xử lý lại toàn bộ khối cỡ $r$ từ đầu (dòng $r$).

* Bước 4.3 — Nếu cột $m$ của $B$ bằng 0: khối $\begin{bmatrix}A_{k-1}&0\\0&F\end{bmatrix}$ đã tách hẳn, hoàn tất Bước 4.

**Bước 5 — Kết luận:** khi $r\le1$: nếu $r=1$, ghi nhận thêm khối $(vị\ trí=1,\ cỡ=1)$. Kết quả $F=\operatorname{diag}(F_1^{(m_1)},\dots,F_s^{(m_s)})$, với mỗi khối cho đa thức đặc trưng thành phần đọc trực tiếp từ dòng đầu $(-p_1,\dots,-p_{m_i})$:

$$
\det(F_i^{(m_i)}-\lambda I) = (-1)^{m_i}\big[\lambda^{m_i}+p_1\lambda^{m_i-1}+\dots+p_{m_i}\big]
$$

---

## 2. Giá trị riêng

**Input:** các khối Frobenius $F_i^{(m_i)}$ thu được ở mục 1.

**Output:** toàn bộ $n$ giá trị riêng của $A$ (kể cả phức, tính cả bội).

**Bước 1 — Với mỗi khối $i$:** giải đa thức bậc $m_i$ bằng phương pháp tìm nghiệm đa thức số:

$$
\lambda^{m_i}+p_1^{(i)}\lambda^{m_i-1}+\dots+p_{m_i}^{(i)}=0
$$

**Bước 2 — Gộp kết quả:** hợp nghiệm của mọi khối lại được $n$ giá trị riêng của $A$.

---

## 3. Véctơ riêng

**Input:** ma trận $P$ và các khối Frobenius $F_i^{(m_i)}$ cùng vị trí $\text{start}_i$ của chúng (từ mục 1), một giá trị riêng $\lambda$ là nghiệm của khối $i$ (từ mục 2).

**Output:** véctơ riêng $x\in\mathbb C^n$ của $A$ ứng với $\lambda$.

**Bước 1 — Véctơ riêng trong hệ khối:**

$$
u=\begin{bmatrix}\lambda^{m_i-1}\\ \vdots\\ \lambda\\ 1\end{bmatrix}
$$

**Bước 2 — Ghép vào véctơ đầy đủ:** đặt $v\in\mathbb C^n$ toàn 0, gán $v_j=u_{j-\text{start}_i+1}$ với $j=\overline{\text{start}_i,\ \text{start}_i+m_i-1}$.

**Bước 3 — Ánh xạ ngược về hệ gốc:** $x=Pv$.

**Bước 4 — (Tuỳ chọn) Chuẩn hoá:** chia $x$ cho thành phần có trị tuyệt đối lớn nhất.

## Lưu ý

* Hệ $A_{k-1}K-KF=B$ ở Bước 4 (mục 1) chỉ có nghiệm duy nhất khi $A_{k-1}$ và khối Frobenius $F$ vừa tách không chung trị riêng. Nếu $A$ có trị riêng bội bị chia vào cả hai phần (ma trận derogatory — hạn chế của bản thân PP Danielevsky), hệ suy biến: `danielevsky.py` bắt lỗi này, in cảnh báo và bỏ qua việc triệt tiêu $B$ — đa thức đặc trưng/trị riêng từng khối vẫn đúng, nhưng véctơ riêng ứng với khối đó có thể sai.
