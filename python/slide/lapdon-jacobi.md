# PP LẶP ĐƠN – LẶP JACOBI

**Hà Thị Ngọc Yến**
*Hà Nội, 2026*

---

## Ý tưởng phương pháp

- Đưa về phương trình tương đương:
  $$Ax = B \in \mathbb{R}^{m \times p} \Leftrightarrow x = Cx + D$$

- Lập dãy:
  $$x_n = Cx_{n-1} + D, \quad x_0 \in \mathbb{R}^{m \times p}$$

- Nếu dãy hội tụ thì giới hạn là nghiệm của phương trình.

---

## Sự hội tụ của PP lặp đơn

- Nếu $||C|| < 1$ thì $x_n = Cx_{n-1} + D, \quad x_0 \in \mathbb{R}^{m \times p}$ hội tụ tới nghiệm đúng duy nhất của phương trình theo đánh giá:
  
  $$||x_n - x^*|| \le \frac{||B||^n}{1 - ||B||} ||x_1 - x_0||$$
  
  $$||x_n - x^*|| \le \frac{||B||}{1 - ||B||} ||x_n - x_{n-1}||$$

---

## Các bước cm sự hội tụ của PP

- Dãy $\{x_n\}$ là dãy Cauchy nên hội tụ.
- Giới hạn của dãy là nghiệm duy nhất của phương trình.
- Cm hai công thức sai số.

---

## Phương pháp lặp Jacobi

- **Ma trận chéo trội hàng:**
  $$|a_{ii}| > \sum_{\substack{j=1 \\ j \neq i}}^m |a_{ij}|$$

- **Ma trận chéo trội cột:**
  $$|a_{ii}| > \sum_{\substack{j=1 \\ j \neq i}}^m |a_{ji}|$$

---

## PP lặp Jacobi (A là ma trận chéo trội hàng)

$$T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{mm}}\right);$$

$$Ax = B \Leftrightarrow x = (I - TA)x + TB,$$

$$C = I - TA, \quad D = Tb$$

$$x^{(0)} \in \mathbb{R}^{m \times p}, \quad x^{(n+1)} = Cx^{(n)} + D.$$

---

## PP lặp Jacobi (A là ma trận chéo trội hàng)

$$B = \begin{bmatrix} 
0 & \frac{-a_{12}}{a_{11}} & \dots & \frac{-a_{1m}}{a_{11}} \\ 
\frac{-a_{21}}{a_{22}} & 0 & \dots & \frac{-a_{2m}}{a_{22}} \\ 
\dots & \dots & \dots & \dots \\ 
\frac{-a_{m1}}{a_{mm}} & \frac{-a_{m2}}{a_{mm}} & \dots & 0 
\end{bmatrix}; \quad 
D = \begin{bmatrix} 
\frac{b_{11}}{a_{11}} & \frac{b_{1p}}{a_{11}} \\ 
\frac{b_{21}}{a_{22}} & \dots & \frac{b_{2p}}{a_{22}} \\ 
\vdots & & \vdots \\ 
\frac{b_{m1}}{a_{mm}} & \frac{b_{mp}}{a_{mm}} 
\end{bmatrix}$$

$$||B||_\infty = \max_{i=\overline{1,m}} \left\{ \frac{1}{|a_{ii}|} \sum_{\substack{j=1 \\ j \neq i}}^m |a_{ij}| \right\} < 1$$

---

## PP lặp Jacobi (A là ma trận chéo trội cột)

$$T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{mm}}\right); \quad x = Ty$$

$$Ax = b \Leftrightarrow ATy = B \Leftrightarrow y = (I - AT)y + B,$$

$$C_1 = I - AT$$

$$y^{(0)} \in \mathbb{R}^{m \times p}, \quad y^{(n+1)} = C_1 y^{(n)} + B.$$

---

## PP lặp Jacobi (A là ma trận chéo trội cột)

$$C_1 = \begin{bmatrix} 
0 & \frac{-a_{12}}{a_{22}} & \dots & \frac{-a_{1m}}{a_{mm}} \\ 
\frac{-a_{21}}{a_{11}} & 0 & \dots & \frac{-a_{2m}}{a_{mm}} \\ 
\dots & \dots & \dots & \dots \\ 
\frac{-a_{m1}}{a_{11}} & \frac{-a_{m2}}{a_{22}} & \dots & 0 
\end{bmatrix}; \quad 
B = \begin{bmatrix} 
b_{11} & b_{1p} \\ 
b_{21} & \dots & b_{2p} \\ 
\vdots & & \vdots \\ 
b_{m1} & b_{mp} 
\end{bmatrix}$$

$$||B||_1 = \max_{j=\overline{1,m}} \left\{ \frac{1}{|a_{jj}|} \sum_{\substack{i=1 \\ i \neq j}}^m |a_{ij}| \right\} < 1$$

---

## PP lặp Jacobi (A là ma trận chéo trội cột)

$$y^{(n+1)} = C_1 y^{(n)} + B$$

$$\Leftrightarrow Ty^{(n+1)} = T(I - AT)T^{-1}Ty^{(n)} + TB$$

$$\Leftrightarrow x^{(n+1)} = (I - TA)x^{(n)} + TB$$

$$\Leftrightarrow x^{(n+1)} = Cx^{(n)} + TB$$

---

## PP lặp Jacobi (A là ma trận chéo trội cột)

- Liên hệ về chuẩn qua phép đổi biến:

  $$||x|| = ||Ty|| \le ||T|| \cdot ||y|| = \frac{||y||}{\min |a_{ii}|} \quad \forall x \in \mathbb{R}^{m \times p}$$

  $$||y|| = ||T^{-1}x|| \le ||T^{-1}|| \cdot ||x|| = \max |a_{ii}| \cdot ||x|| \quad \forall y \in \mathbb{R}^{m \times p}$$

---

## PP lặp Jacobi (A là ma trận chéo trội cột)

- Hệ quả là:

  $$||x^{(n)} - x^*||_1 \le \lambda \frac{||C_1||_1}{1 - ||B_1||_1} ||x^{(n)} - x^{(n-1)}||_1$$

  $$||x^{(n)} - x^*||_1 \le \lambda \frac{||C_1||_1^n}{1 - ||B_1||_1} ||x^{(1)} - x^{(0)}||_1$$

  $$\lambda = \frac{\max |a_{ii}|}{\min |a_{ii}|}.$$