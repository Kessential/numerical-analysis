# PP GAUSS VÀ GAUSS-JORDAN
## GIẢI PT Ax = b

**Hà Thị Ngọc Yến**
**Hà nội, 2026**

---

## 3.1 Bài toán

* $Ax = b$, $A$ cỡ $m \times n$, $b$ cỡ $m \times 1$, $x$ cỡ $n \times 1$
* Mô hình kinh tế Leotief:
  * Mô hình định giá: $X = CX$
    $C(i,j)$ = tỷ lệ cung cấp sp của ngành $j$ cho ngành $i$ = nhu cầu của ngành $i$ đối với ngành $j$
  * Mô hình định lượng: $X = CX + D$

---

## Bài toán

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{bmatrix}
\begin{bmatrix}
x_{11} & \cdots & x_{1k} \\
x_{21} & \cdots & x_{2k} \\
\vdots & \vdots & \vdots \\
x_{n1} & \cdots & x_{nk}
\end{bmatrix}
=
\begin{bmatrix}
b_{11} & \cdots & b_{1k} \\
b_{21} & \cdots & b_{2k} \\
\vdots & \vdots & \vdots \\
b_{m1} & \cdots & b_{mk}
\end{bmatrix}
$$

---

## Phương pháp Gauss

**Ý tưởng:**

* **Quy trình thuận (QTT):** Dùng phép khử dần ẩn khỏi các phương trình đề đưa ma trận bổ sung về dạng bậc thang.
* **Quy trình nghịch (QTN):** Dùng phép thế từ hệ bậc thang để tìm dần giá trị các ẩn.

---

## PP Gauss – QTT

* **B1:** Khởi tạo $i=1;\ j=1;\ ind=[0,0,...,0]_{1 \times m}$
* **B2:** Kiểm tra nếu $a_{ij} \neq 0 \Rightarrow ind[i]=j \Rightarrow$ B3, trái lại thì sang B6
* **B3:** Nếu $i = m$ thì kết thúc QTT
  Nếu không, thì sang B4
* **B4:** Cho $k$ chạy từ $i+1$ đến $m$, thực hiện biến đổi
  $$L_k - \frac{a_{kj}}{a_{ij}} \times L_i \Rightarrow L_k$$

---

## PP Gauss - QTT

* **B5:** Nếu $j = n$ thì QTT kết thúc, trái lại $i = i+1;\ j = j+1 \Rightarrow$ B2
* **B6:** Cho $t = i+1$.
* **B7:** Kiểm tra nếu $a_{tj} \neq 0$ thì đổi chỗ 2 hàng $t$ và $i$, và $ind[i] = j \Rightarrow$ B3, trái lại sang B8
* **B8:** Nếu $t = m;\ j = n$ thì QTT kết thúc
  Nếu $t = m;\ j < n$ thì $j = j+1 \Rightarrow$ B2,
  Nếu $t < m$ thì $t = t+1 \Rightarrow$ B7

---

## PP Gauss - QTN

Dành cho các bạn tự viết

---

## Phương pháp Gauss - Jordan

**Ý tưởng:**

* Hạn chế sai số tính toán khi gặp các phép chia cho số gần 0 bằng cách chọn phần tử khử thích hợp
* Dùng phép khử ẩn thứ $k$ (tương ứng với cột có chứa phần tử khử) khỏi tất cả các hàng không chứa phần tử khử

---

## Chọn phần tử khử

* **Ưu tiên 1:**
  Chọn $a_{pq}^{(k)} = 1; 2; 4; 5...$ để các phép chia cho $a_{pq}^{(k)}$ không có sai số hoặc sai số nhỏ.
* **Ưu tiên 2:**
  Chọn $a_{pq}^{(k)}$ sao cho $\left| a_{pq}^{(k)} \right| = \max_{i,j} \left| a_{ij}^{(k-1)} \right|$

*Chú ý: Phần tử khử thứ $k$ được chọn từ các hàng và cột không chứa các phần tử khử đã chọn trước đó.*

---

## Quá trình khử

Cho $t$ chạy từ $1$ đến $m$, $t \neq p$, thực hiện phép biến đổi

$$L_t - \frac{a_{tq}}{a_{pq}} L_p \Rightarrow L_t$$

---

## Quá trình khử 

*(Phiên bản áp dụng cho giải hệ $Ax = b$ từ File 2)*

$$
\left[
\begin{array}{cccc|c}
a_{11} & a_{12} & \cdots & a_{1n} & b_1 \\
a_{21} & a_{22} & \cdots & a_{2n} & b_2 \\
\vdots & \vdots & \cdots & \vdots & \vdots \\
a_{p1} & a_{p2} & \cdots & a_{pn} & b_p \\
\vdots & \vdots & \cdots & \vdots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn} & b_m
\end{array}
\right]
\Leftrightarrow
\left[
\begin{array}{cccccc|c}
a_{11}^1 & a_{12}^1 & \cdots & 0 & \cdots & a_{1n}^1 & b_1^1 \\
a_{21}^1 & a_{22}^1 & \cdots & 0 & \cdots & a_{2n}^1 & b_2^1 \\
\vdots & \vdots & \cdots & \vdots & \cdots & \vdots & \vdots \\
a_{p1} & a_{p2} & \cdots & a_{pq} & \cdots & a_{pn} & b_p \\
\vdots & \vdots & \cdots & \vdots & \vdots & \vdots & \vdots \\
a_{m1}^1 & a_{m2}^1 & \cdots & 0 & \cdots & a_{mn}^1 & b_m^1
\end{array}
\right]
$$

> **Ghi chú (Ở File 1):** Cột $b$ được thay bằng ma trận $B$ có $r$ cột ($b_{11} \dots b_{1r}$, $b_{21} \dots b_{2r}$, v.v.), khi đó các giá trị $b^1_i$ ở vế phải cũng tương ứng mở rộng thành $b^1_{1} \dots b^1_{1r}$, v.v.

---

## Bảng giá trị ví dụ minh họa (Slide 12)

$$
\begin{bmatrix}
2 & 4 & 5 & -6 & 7 & 3 \\
2 & 3 & 5 & 2 & 1 & 4 \\
5 & 8.5 & 13 & -3 & 8.5 & 9 \\
-3 & -4 & -9 & -11 & 4.3 & -8
\end{bmatrix}
$$