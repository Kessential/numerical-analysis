# Khai triển kỳ dị (SVD): giá trị kỳ dị, khai triển, số điều kiện, PP xuống thang, xấp xỉ ảnh

Tổng hợp lại từ slide `svd.md` thành thuật toán chi tiết để cài đặt (code tại `svd_gia_tri_ky_di_lon_nhat.py`, `svd_gia_tri_ky_di.py`, `svd_vecto_trai.py`, `svd_khai_trien.py`, `svd_so_dieu_kien.py`, `svd_xap_xi_anh.py`).

## 1. Tìm giá trị kỳ dị lớn nhất và vector kỳ dị tương ứng — PP lũy thừa

**Input:** $A \in \mathbb{R}^{m\times n}$, sai số $\varepsilon$, số lần lặp tối đa $N$.

**Output:** $(\sigma_1, v_1, u_1)$ — giá trị kỳ dị lớn nhất và cặp véc-tơ kỳ dị tương ứng.

**Bước 1 — PP lũy thừa trên $M:=A^TA$:** nếu $\|A\|_2\approx 0$ ($A$ suy biến hoàn toàn về 0): dừng, báo không có giá trị kỳ dị khác 0. Ngược lại, áp dụng PP lũy thừa (mục 1, `luy_thua_xuong_thang.md`) cho $M$, thu được $(\lambda_1, v_1)$ (luôn rơi vào TH1).

**Bước 2 — Chuẩn hoá:** $v_1 \leftarrow v_1/\|v_1\|_2$ (PP lũy thừa chỉ chuẩn hoá $v_1$ theo toạ độ lớn nhất, chưa phải chuẩn 2).

**Bước 3 — Tính giá trị kỳ dị:** $\sigma_1 := \sqrt{\lambda_1}$.

**Bước 4 — Tính vector kỳ dị trái:** $u_1 := Av_1/\sigma_1$.

---

## 2. Tìm các giá trị kỳ dị còn lại và vector kỳ dị phải — PP xuống thang

**Input:** $A \in \mathbb{R}^{m\times n}$, số giá trị kỳ dị cần tìm $r'\le \text{rank}A$, sai số $\varepsilon$, số lần lặp tối đa $N$.

**Output:** $(\sigma_1,v_1),\dots,(\sigma_{r'},v_{r'})$ với $\sigma_1\ge\sigma_2\ge\dots>0$.

**Bước 1 — Khởi tạo:** $M^{(1)} := A^TA$.

**Bước 2 — Với $k=1,2,\dots,r'$:**

* Bước 2.1 — Tìm giá trị kỳ dị trội của $M^{(k)}$ bằng PP lũy thừa (mục 1, `luy_thua_xuong_thang.md`): thu được $(\lambda_k, v_k^{(k)})$ (luôn rơi vào TH1; với $k=1$, bước này trùng với Phần 1). Nếu $\lambda_k \le \varepsilon$: dừng thuật toán — đã hết các giá trị kỳ dị khác 0 ($k-1 = \text{rank}A$).
* Bước 2.2 — Tính giá trị kỳ dị: $\sigma_k := \sqrt{\lambda_k}$.
* Bước 2.3 — Xuống thang (cách chọn 2): $s_k:=$ toạ độ $i$ làm $|v_{k,i}^{(k)}|$ lớn nhất, chuẩn hoá $v_k^{(k)} \leftarrow v_k^{(k)}/v_{k,s_k}^{(k)}$, $a_k:=$ hàng $s_k$ của $M^{(k)}$, $M^{(k+1)} := M^{(k)} - v_k^{(k)}a_k$.
* Bước 2.4 — Khôi phục $v_k$ về không gian gốc bằng công thức (mục 2, `luy_thua_xuong_thang.md`), áp dụng lần lượt $j=k-1,\dots,1$:

$$
v_k^{(j)} := (\lambda_j-\lambda_k)v_k^{(j+1)} - (a_j\cdot v_k^{(j+1)})v_j^{(j)}, \qquad v_k := v_k^{(1)}/\|v_k^{(1)}\|_2
$$

---

## 3. Vector kỳ dị trái (tổng quát)

**Input:** $A$, các cặp $(\sigma_i,v_i)_{i=1}^{r'}$ từ Phần 2.

**Output:** $u_1,\dots,u_{r'}$ trực chuẩn.

**Bước 1 — Tính:** $u_i := \dfrac{Av_i}{\sigma_i}, \quad i=\overline{1,r'}$.

---

## 4. Khai triển kỳ dị đầy đủ

**Input:** $A \in \mathbb{R}^{m\times n}$, sai số $\varepsilon$, số lần lặp tối đa $N$ (dùng cho PP lũy thừa ở Phần 2).

**Output:** $U,\Sigma,V$ và $A$ tái tạo lại (kiểm tra sai số $\|A-U\Sigma V^T\|_2$).

**Bước 1 — Tìm toàn bộ phổ kỳ dị:** chạy Phần 2 với $r'=\min(m,n)$ (dừng tự nhiên khi $\lambda_k\le\varepsilon$, tức $r'$ giảm còn $r=\text{rank}A$), rồi Phần 3, thu được $(\sigma_i,v_i,u_i)_{i=1}^r$.

**Bước 2 — Ghép ma trận:** $U:=[u_1\ \dots\ u_r]$, $\Sigma:=\text{diag}(\sigma_1,\dots,\sigma_r)$, $V:=[v_1\ \dots\ v_r]$.

**Bước 3 — Tái tạo và kiểm tra:** $A_{rec}:=U\Sigma V^T$, $\delta:=\max(100\varepsilon, 10^{-8})$ (nới hơn $\varepsilon$ của từng bước lũy thừa, vì sai số tích luỹ qua $r$ bước xuống thang liên tiếp ở Phần 2). Nếu $\|A-A_{rec}\|_2 \le \delta$: kết luận khai triển đúng, trả về $U,\Sigma,V$. Ngược lại: báo lỗi (PP lũy thừa ở một bước nào đó của Phần 2 chưa hội tụ).

---

## 5. Số điều kiện

**Input:** $A$ vuông khả nghịch, sai số $\varepsilon$, số lần lặp tối đa $N$.

**Output:** $\text{cond}\ A = \|A\|_2\,\|A^{-1}\|_2 = \sigma_1/\sigma_n$, hoặc thông báo nếu $A$ suy biến (không khả nghịch, $\text{cond}\ A=\infty$).

**Bước 1 — Kiểm tra đầu vào:** nếu $\det A \approx 0$ ($A$ suy biến, $\sigma_n=0$): dừng, báo $\text{cond}\ A=\infty$.

**Bước 2 — Tìm $\sigma_{\max}$:** chạy Phần 1 trên $A$, được $\sigma_1$ và $M:=A^TA$ (dùng lại ở Bước 3).

**Bước 3 — Tìm $\sigma_{\min}$:** tính $M^{-1}$ (nghịch đảo — dùng phương pháp bất kỳ trong `algo/nghich_dao.md`). PP lũy thừa trên $M^{-1} \to \mu_1 \Rightarrow \sigma_n := 1/\sqrt{\mu_1}$.

**Bước 4 — Kết luận:** $\text{cond}A := \sigma_1/\sigma_n$.

---

## 6. Bài toán xấp xỉ ảnh (nén ảnh bằng SVD hạng thấp)

**Input:** ảnh xám $A_{m\times n}$, hạng xấp xỉ $k$, sai số $\varepsilon$, số lần lặp tối đa $N$ (dùng cho PP lũy thừa ở Phần 2).

**Output:** $A_k$ (ảnh nén), tỉ lệ nén $\dfrac{k(m+n+1)}{mn}$, sai số $\dfrac{\|A-A_k\|_F}{\|A\|_F}$.

**Bước 1 — Chọn ma trận Gram nhỏ hơn:** dùng $M=A^TA$ ($n\times n$) nếu $n\le m$, ngược lại $M=AA^T$ ($m\times m$), để giảm chi phí xuống thang.

**Bước 2 — Lấy $k$ cặp kỳ dị lớn nhất:** chạy Phần 2 cho $M$ đã chọn, $k$ bước, được $k$ cặp $(\sigma_i, \text{vtr riêng } w_i \text{ của } M)$. Suy ra vế còn lại tương tự công thức Phần 3: $u_i=Av_i/\sigma_i$ ($v_i:=w_i$) nếu $M=A^TA$, hoặc $v_i=A^Tu_i/\sigma_i$ ($u_i:=w_i$) nếu $M=AA^T$.

**Bước 3 — Dựng ảnh xấp xỉ:** $A_k:=\sum_{i=1}^k\sigma_iu_iv_i^T$.

**Bước 4 — Đánh giá:** tỉ lệ lưu trữ $k(m+n+1)/(mn)$, sai số tương đối $\|A-A_k\|_F/\|A\|_F$.