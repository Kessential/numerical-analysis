# PP LẶP GAUSS - SEIDEL

**Hà Thị Ngọc Yến**
Hà Nội, 2025

---

## Ý tưởng phương pháp

- Cải tiến phương pháp lặp đơn khi gặp phương trình lặp hội tụ với tốc độ chậm
- Thông tin mới được sử dụng càng sớm càng tốt

$$ x^{(n)} = Cx^{(n-1)} + D, \quad x^{(0)} \in \mathbb{R}^{m \times p} $$

---

## Ý tưởng phương pháp

- Phương pháp lặp đơn đã chứng tỏ dãy lặp
$$ x^{(n)} = Cx^{(n-1)} + D, \quad x^{(0)} \in \mathbb{R}^{m \times p} $$
khi $\|B\| < 1$ thì lần xấp xỉ sau tốt hơn lần xấp xỉ trước đó

$$ \|x^{(n)} - x^*\| < \|x^{(n-1)} - x^*\| $$

---

## Ý tưởng phương pháp

Phương trình lặp Jacobi

$$
\begin{aligned}
x_{1i}^{(n+1)} &= & & c_{12}x_{2i}^{(n)} &+ \cdots &+ c_{1m}x_{mi}^{(n)} &+ d_{1i} \\
x_{2i}^{(n+1)} &= \boxed{c_{21}x_{1i}^{(n)}} & & &+ \cdots &+ c_{2m}x_{mi}^{(n)} &+ d_{2i} \\
x_{3i}^{(n+1)} &= \boxed{c_{31}x_{1i}^{(n)}} &+ \boxed{c_{32}x_{2i}^{(n)}} & &+ \cdots &+ c_{3m}x_{mi}^{(n)} &+ d_{3i} \\
\vdots \quad &= \quad \vdots & \vdots \quad \quad & \vdots & \vdots \quad & \quad \vdots \\
x_{mi}^{(n+1)} &= \boxed{c_{m1}x_{1i}^{(n)}} &+ \boxed{c_{m2}x_{2i}^{(n)}} &+ \cdots & & &
\end{aligned}
$$

---

## Ý tưởng phương pháp

Thay ngay khi có giá trị mới

$$
\begin{aligned}
x_{1i}^{(n+1)} &= & & c_{12}x_{2i}^{(n)} &+ \cdots &+ c_{1m}x_{mi}^{(n)} &+ d_{1i} \\
x_{2i}^{(n+1)} &= \boxed{c_{21}x_{1i}^{(n+1)}} & & &+ \cdots &+ c_{2m}x_{mi}^{(n)} &+ d_{2i} \\
x_{3i}^{(n+1)} &= \boxed{c_{31}x_{1i}^{(n)}} &+ \boxed{c_{32}x_{2i}^{(n)}} & &+ \cdots &+ c_{3m}x_{mi}^{(n)} &+ d_{3i} \\
\vdots \quad &= \quad \vdots & \vdots \quad \quad & \vdots & \vdots \quad & \quad \vdots \\
x_{mi}^{(n+1)} &= \boxed{c_{m1}x_{1i}^{(n)}} &+ \boxed{c_{m2}x_{2i}^{(n)}} &+ \cdots & & &
\end{aligned}
$$

---

## Ý tưởng phương pháp

Thay ngay khi có giá trị mới

$$
\begin{aligned}
x_{1i}^{(n+1)} &= & & c_{12}x_{2i}^{(n)} &+ \cdots &+ c_{1m}x_{mi}^{(n)} &+ d_{1i} \\
x_{2i}^{(n+1)} &= \boxed{c_{21}x_{1i}^{(n+1)}} & & &+ \cdots &+ c_{2m}x_{mi}^{(n)} &+ d_{2i} \\
x_{3i}^{(n+1)} &= \boxed{c_{31}x_{1i}^{(n+1)}} &+ \boxed{c_{32}x_{2i}^{(n+1)}} & &+ \cdots &+ c_{3m}x_{mi}^{(n)} &+ d_{3i} \\
\vdots \quad &= \quad \vdots & \vdots \quad \quad & \vdots & \vdots \quad & \quad \vdots \\
x_{mi}^{(n+1)} &= \boxed{c_{m1}x_{1i}^{(n)}} &+ \boxed{c_{m2}x_{2i}^{(n)}} &+ \cdots &+ & & d_{mi}
\end{aligned}
$$

---

## Ý tưởng phương pháp

Phương trình lặp Seidel

$$
\begin{aligned}
x_{1i}^{(n+1)} &= & & c_{12}x_{2i}^{(n)} &+ \cdots &+ c_{1m}x_{mi}^{(n)} &+ d_{1i} \\
x_{2i}^{(n+1)} &= \boxed{c_{21}x_{1i}^{(n+1)}} & & &+ \cdots &+ c_{2m}x_{mi}^{(n)} &+ d_{2i} \\
x_{3i}^{(n+1)} &= \boxed{c_{31}x_{1i}^{(n+1)}} &+ \boxed{c_{32}x_{2i}^{(n+1)}} & &+ \cdots &+ c_{3m}x_{mi}^{(n)} &+ d_{3i} \\
\vdots \quad &= \quad \vdots & \vdots \quad \quad & \vdots & \vdots \quad & \quad \vdots \\
x_{mi}^{(n+1)} &= \boxed{c_{m1}x_{1i}^{(n+1)}} &+ \boxed{c_{m2}x_{2i}^{(n+1)}} &+ \cdots &+ & & d_{mi}
\end{aligned}
$$

---

## Ý tưởng phương pháp

Thay ngay khi có giá trị mới

$$ x^{(n+1)} = Ux^{(n)} + Lx^{(n+1)} + d $$

$$
L = \begin{bmatrix} 0 & & 0 \\ \vdots & \ddots & \\ b_{ij} & \dots & 0 \end{bmatrix} ; \quad U = \begin{bmatrix} 0 & \dots & b_{ij} \\ & \ddots & \vdots \\ 0 & & 0 \end{bmatrix}
$$

$$ x^{(n+1)} = (I - L)^{-1}Ux^{(n)} + d $$

---

## SỰ HỘI TỤ

Định lý 1: Giả sử $B$ là ma trận vuông cấp $m$, khi đó, các mệnh đề sau là tương đương:

1. Ma trận hội tụ, tức là $\lim_{n \to \infty} B^n = 0$
2. $\lim_{n \to \infty} \|B\|^n = 0$
3. Mọi giá trị riêng đều có modun nhỏ hơn 1 
$$ \rho(B) = \max_i |\lambda_i| < 1. $$

---

## SỰ HỘI TỤ

Định lý 2: Cho ma trận vuông $B \in \mathbb{R}^{m \times m}$.
Khi đó, với mỗi $\varepsilon > 0$ tồn tại một chuẩn trên $\mathbb{R}^m$
sao cho $$ \rho(B) \le \|B\| \le \rho(B) + \varepsilon $$

Hệ quả 3: Nếu $\|B\| < 1$ với một chuẩn nào đó thì $B$ là ma trận hội tụ

---

## SỰ HỘI TỤ

Viết lại pt lặp Gauss-Seidel

$$ A = D_A - L_A - U_A \Rightarrow I - TA = TL_A + TU_A, \quad T = (D_A)^{-1} $$
$$ (D_A - L_A)x^{(n+1)} = U_Ax^{(n)} + b $$
$$ x^{(n+1)} = (D_A - L_A)^{-1}U_Ax^{(n)} + b $$
$$ M := (D_A - L_A)^{-1}U_A $$

---

## SỰ HỘI TỤ

$$ \det(\lambda I - M) := \det[\lambda I - (D_A - L_A)^{-1}U_A] = 0 $$
$$ \Leftrightarrow \det[\lambda(D_A - L_A) - U_A] = 0 $$

Bổ đề: Nếu $|\lambda| \ge 1$ thì ma trận $\lambda(D_A - L_A) - U_A$ khả nghịch

---

## SỰ HỘI TỤ

Định lý 4:

Nếu $A$ là ma trận chéo trội thì

$$ \rho(M) < 1, \quad M := (D_A - L_A)^{-1}U_A $$

---

## Sai số

$$ s_{row} = 0, \qquad q_{row} = \max_{i=1,m} \frac{\sum_{j>i} |a_{ij}|}{|a_{ii}| - \sum_{j<i} |a_{ij}|} $$

$$ s_{col} = \max_{j=1,m} \sum_{i>j} \left| \frac{a_{ij}}{a_{jj}} \right|, \qquad q_{col} = \max_{j=1,m} \frac{\sum_{i<j} |a_{ij}|}{|a_{jj}| - \sum_{i>j} |a_{ij}|} $$

$$ \|X^{(k)} - X^*\| \le \frac{q^k}{(1-s)(1-q)} \|X^{(1)} - X^{(0)}\| $$

$$ \|X^{(k)} - X^*\| \le \frac{q}{(1-s)(1-q)} \|X^{(k)} - X^{(k-1)}\| $$

---

## Bảng Dữ Liệu

|    |    |     |   |    |
|---:|---:|----:|---|---:|
| 10 |  5 |   7 |   | 11 |
|  2 | 15 |   3 |   | 12 |
| -3 |  1 | -30 |   | 19 |

<br>

| 0 | 1 | 2 | 3 | 4 | 5 |
|--:|--:|--:|--:|--:|--:|
| 1 | -0.1 | 1.215355556 | 1.231443995 | 1.22043199 | 1.218805285 |
| 1 | 0.613333333 | 0.75853037 | 0.781724376 | 0.783359785 | 0.783345537 |
| 1 | -0.602888889 | -0.729584543 | -0.730420254 | -0.729264539 | -0.729102344 |
| <br> | | | | | |
| 1 | -0.1 | 1.356666667 | 1.048777778 | 1.245722222 | 1.195661852 |
| 1 | 0.466666667 | 0.953333333 | 0.740666667 | 0.807607407 | 0.776608148 |
| 1 | -0.7 | -0.607777778 | -0.737222222 | -0.713522222 | -0.730985309 |

---

## Ví dụ

|   |    |    |    |   |
|--:|---:|---:|---:|--:|
| 15| -2 | -4 |  2 | 3 |
|  3| 15 | -8 |  3 | 5 |
|  2|  4 | 15 | -5 | 7 |
|  3| -2 |  1 |  8 |-4 |