# PP PHÂN TÁCH LU VÀ PP CHOLESKI GIẢI PT Ax=b.

**Hà Thị Ngọc Yến**
**Hà nội, 2026**

---

## Bài toán

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nn}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
\vdots \\
x_n
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_n
\end{bmatrix}
$$

---

## PP phân tách LU

**Ý tưởng:**

* Khi giải HĐSTT, nếu ma trận liên kết là ma trận tam giác thì việc giải nó là đơn giản, tính toán ít vì chỉ cần thực hiện QTN (Quy trình nghịch) của Gauss là đủ.
* Tách ma trận liên kết thành tích của hai ma trận tam giác. 

---

## PP phân tách LU

Tách Ma trận A thành tích của ma trận tam giác dưới (**L**ower) và ma trận tam giác trên (**U**pper)

$$A = LU$$

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nn}
\end{bmatrix}
=
\begin{bmatrix}
l_{11} & & & \\
l_{21} & l_{22} & & \\
\vdots & \vdots & \ddots & \\
l_{n1} & l_{n2} & \cdots & l_{nn}
\end{bmatrix}
\times
\begin{bmatrix}
u_{11} & u_{12} & \cdots & u_{1n} \\
& u_{22} & \cdots & u_{2n} \\
& & \ddots & \vdots \\
& & & u_{nn}
\end{bmatrix}
$$

---

## PP phân tách LU

Phương trình trên tương ứng với một hệ có

* $n \times n \qquad \text{phương trình}$
* $n \times n + n \quad \text{ẩn}$

Có thể chọn trước $n$ giá trị, các giá trị còn lại tìm được bằng cách giải phương trình.

---

## PP phân tách LU

Chọn $u_{ii} = 1, i = \overline{1,n}$.

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nn}
\end{bmatrix}
=
\begin{bmatrix}
l_{11} & & & \\
l_{21} & l_{22} & & \\
\vdots & \vdots & \ddots & \\
l_{n1} & l_{n2} & \cdots & l_{nn}
\end{bmatrix}
\times
\begin{bmatrix}
1 & u_{12} & \cdots & u_{1n} \\
& 1 & \cdots & u_{2n} \\
& & \ddots & \vdots \\
& & & 1
\end{bmatrix}
$$

$$
\Leftrightarrow
\begin{cases}
l_{11} = a_{11} \\
l_{11}u_{1k} = a_{1k}, \quad k = \overline{2,n} \\
l_{i1} = a_{i1}, \quad i = \overline{2,n} \\
l_{i1}u_{1k} + l_{i2}u_{2k} + \cdots + l_{ii}u_{ik} = a_{ik}, \quad k = \overline{2,n}, \ i = \overline{2,n}
\end{cases}
$$

---

## PP phân tách LU

$$
\begin{cases}
l_{i1} = a_{i1}, \quad i = \overline{1,n} \\
u_{1k} = \frac{a_{1k}}{a_{11}}, \quad k = \overline{2,n}
\end{cases}
$$

$$
\begin{cases}
l_{i2} = a_{i2} - l_{i1}u_{12} = a_{i2} - \frac{a_{i1}a_{12}}{a_{11}}, \quad i \ge 2 \\
u_{2k} = \frac{a_{2k} - l_{21}u_{1k}}{l_{22}}, \quad k \ge 2
\end{cases}
$$

$$
\begin{cases}
l_{i3} = a_{i3} - l_{i1}u_{13} - l_{i2}u_{23}, \quad i \ge 3 \\
u_{3k} = \frac{a_{3k} - (l_{31}u_{1k} + l_{32}u_{2k})}{l_{33}}, \quad k \ge 3
\end{cases}
$$

$$
\begin{cases}
l_{it} = a_{it} - l_{i1}u_{1t} - l_{i2}u_{2t} - \cdots - l_{i,t-1}u_{t-1,t}, \quad i \ge t \\
u_{tk} = \frac{a_{tk} - (l_{t1}u_{1k} + l_{t2}u_{2k} + \cdots + l_{t,t-1}u_{t-1,k})}{l_{tt}}, \quad k \ge t
\end{cases}
$$

---

## PP phân tách Choleski

Chọn $\quad U = L^T \Rightarrow A = U^T U$
Điều kiện: A là ma trận đối xứng không suy biến

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{12} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{1n} & a_{2n} & \cdots & a_{nn}
\end{bmatrix}
=
\begin{bmatrix}
u_{11} & & & \\
u_{12} & u_{22} & & \\
\vdots & \vdots & \ddots & \\
u_{1n} & u_{2n} & \cdots & u_{nn}
\end{bmatrix}
\times
\begin{bmatrix}
u_{11} & u_{12} & \cdots & u_{1n} \\
& u_{22} & \cdots & u_{2n} \\
& & \ddots & \vdots \\
& & & u_{nn}
\end{bmatrix}
$$

---

## PP phân tách Choleski

Chọn $\quad U = L^T \Rightarrow A = U^T U$
Điều kiện: A là ma trận đối xứng không suy biến

$$
\Leftrightarrow
\begin{cases}
u_{11}^2 = a_{11} \\
u_{11}u_{1k} = a_{1k}, \quad k = \overline{2,n} \\
u_{12}^2 + u_{22}^2 = a_{22} \\
u_{12}u_{1k} + u_{22}u_{2k} = a_{2k}, \quad k = \overline{3,n} \\
u_{1i}^2 + u_{2i}^2 + \cdots + u_{ii}^2 = a_{ii}, \quad i = \overline{3,n} \\
u_{1i}u_{1k} + u_{2i}u_{2k} + \cdots + u_{ii}u_{ik} = a_{ik}, \quad k = \overline{i+1,n}, \ i = \overline{3,n-1}
\end{cases}
$$

---

## PP phân tách Choleski

Chọn $\quad U = L^T \Rightarrow A = U^T U$
Điều kiện: A là ma trận đối xứng không suy biến

$$
\Leftrightarrow
\begin{cases}
u_{11} = \sqrt{a_{11}} \\
u_{1k} = \frac{a_{1k}}{u_{11}}, \quad k = \overline{2,n} \\
u_{22} = \sqrt{a_{22} - u_{12}^2} \\
u_{2k} = \frac{a_{2k} - u_{12}u_{1k}}{u_{22}}, \quad k = \overline{3,n} \\
...........
\end{cases}
$$

---

## PP phân tách Choleski

Chọn $\quad U = L^T \Rightarrow A = U^T U$
Điều kiện: A là ma trận đối xứng không suy biến

$$
\Rightarrow
\begin{cases}
u_{ii} = \sqrt{a_{ii} - \left( u_{1i}^2 + u_{2i}^2 + \cdots + u_{i-1,i}^2 \right)} \\
\\
u_{ik} = \frac{a_{ik} - \left( u_{1i}u_{1k} + u_{2i}u_{2k} + \cdots + u_{i-1,i}u_{i-1,k} \right)}{u_{ii}}
\end{cases}
$$
$k = \overline{i+1,n}, \ i = \overline{2,n-1}$

$$
u_{nn} = \sqrt{a_{nn} - \left( u_{1n}^2 + u_{2n}^2 + \cdots + u_{n-1,n}^2 \right)}
$$

---

## PP phân tách Choleski

Chú ý: Khi viết thuật toán, cần chú ý một số trường hợp đặc biệt:

1. $u_{ii} = 0$ cần thông báo ma trận A suy biến, không thỏa mãn điều kiện thực hiện phương pháp.

2. $u_{ii} < 0$ cần chuyển sang tính toán với số phức và thực hiện các phép toán như trong công thức nhưng trên tập số phức.

---

## Phương pháp Choleski

1. Trường hợp A đối xứng, cần sử dụng quy trình ngược của PP Gauss cho hệ (1):

$$
Ax = b \Leftrightarrow U^T U x = b \Leftrightarrow
\begin{cases}
Ux = y \\
U^T y = b
\end{cases} \quad (1)
$$

2. Trường hợp A không đối xứng:

$$ \exists A^{-1} $$
$$ Ax = b \Leftrightarrow A^T Ax = A^T b, \quad M := A^T A, \ d := A^T b $$

$$
Mx = d \Leftrightarrow U^T U x = d \Leftrightarrow
\begin{cases}
Ux = y \\
U^T y = d
\end{cases}
$$

---

## Phương pháp Choleski

Chú ý:

$$ x^T A^T Ax > 0 \quad \forall x \neq 0 \Leftrightarrow \det A \neq 0 $$

Do đó, khi lập trình, không cần kiểm tra điều kiện không suy biến của A vì nếu A suy biến sẽ có giá trị $u_{ii} = 0$.

Thuật toán: dành cho các bạn tự viết