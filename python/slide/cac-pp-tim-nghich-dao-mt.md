Dưới đây là nội dung của bài thuyết trình được chuyển đổi sang định dạng Markdown. Các công thức toán học được viết bằng mã LaTeX để đảm bảo hiển thị chính xác.

---

# Các phương pháp tìm ma trận nghịch đảo

Hà nội, 2026

---

## CÁC PHƯƠNG PHÁP TÌM ĐÚNG MA TRẬN NGHỊCH ĐẢO

* Phương pháp dùng ma trận phần phụ
* Phương pháp Gauss và Gauss-Jordan
* Phương pháp Choleski
* Phương pháp viền quanh

---

## PP dùng ma trận phần phụ

* Định lý:

$$ C = [c_{ij}]_{n \times n} = \left[ (-1)^{i+j} \det A_{ij} \right] $$

$$ A C^t = (\det A).E $$

---

## PP Gauss và Gauss-Jordan

* Đưa về giải phương trình $\quad AX = E$

bằng phương pháp G và GJ

$$ [A | E] \Leftrightarrow [E | A^{-1}] $$

---

## PP Choleski

* Phân tách ma trận đối xứng $A$ (hoặc $A^t A$) thành tích hai ma trận tam giác là chuyển vị của nhau:

$$ A = S^t S \qquad (A^t A = S^t S) $$

$$ S = [s_{ij}]_{n \times n}, \quad s_{ij} = 0 \quad \forall i > j $$

---

## PP Choleski

* Giải n phương trình

$$ S^t S x^{(i)} = e_i \quad \lor \quad S^t S x^{(i)} = A^t e_i $$

* Ghép các cột nghiệm $x^{(i)}$ tạo thành ma trận nghịch đảo cần tìm

* C2: Tạo gói tìm nghịch đảo của ma trận tam giác

---

## PP viền quanh

* Phân khối ma trận:

$$ A_n = \begin{bmatrix} A_{n-1} & \alpha_{n-1, 1} \\ \alpha_{1, n-1} & a_{n,n} \end{bmatrix}, \qquad A_n^{-1} = \begin{bmatrix} B_{n-1} & \beta_{n-1, 1} \\ \beta_{1, n-1} & b_{n,n} \end{bmatrix} $$

---

## PP viền quanh

* Nhân hai ma trận

$$ A.A^{-1} = E \iff \begin{cases} A_{n-1} B_{n-1} + \alpha_{n-1, 1} \beta_{1, n-1} = E_{n-1} & (1) \\ A_{n-1} \beta_{n-1, 1} + \alpha_{n-1, 1} b_{n,n} = 0_{n-1, 1} & (2) \\ \alpha_{1, n-1} B_{n-1} + a_{n,n} \beta_{1, n-1} = 0_{1, n-1} & (3) \\ \alpha_{1, n-1} \beta_{n-1, 1} + a_{n,n} b_{n,n} = 1 & (4) \end{cases} $$

---

## PP viền quanh

* Lấy (2) nhân với $\alpha_{1,n-1} A_{n-1}^{-1}$ từ bên trái rồi trừ đi (4), ta tìm được

$$ b_{n,n} = \frac{1}{a_{n,n} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1}} \qquad (5) $$

* Thay vào (2) ta có

$$ \beta_{n-1,1} = \frac{-A_{n-1}^{-1} \alpha_{n-1,1}}{a_{n,n} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1}} \qquad (6) $$

---

## PP viền quanh

* Làm tương tự với cặp phương trình (1) và (3) ta thu được:

$$ \beta_{1,n-1} = \frac{-\alpha_{1,n-1} A_{n-1}^{-1}}{a_{n,n} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1}} \qquad (7) $$

$$ B_{n-1} = A_{n-1}^{-1} \left( E_{n-1} - \frac{\alpha_{n-1,1} \alpha_{1,n-1} A_{n-1}^{-1}}{a_{n,n} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1}} \right) \qquad (8) $$

---

## PP viền quanh

* Điều kiện:

$$ \exists A_{n-1}^{-1} $$

$$ a_{n,n} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1} \neq 0 $$

---

## PP viền quanh

* Điều kiện:

$$ \begin{bmatrix} A_{n-1}^{-1} & 0_{n-1,1} \\ 0_{1,n-1} & 1 \end{bmatrix} \times \begin{bmatrix} A_{n-1} & \alpha_{n-1,1} \\ \alpha_{1,n-1} & a_{nn} \end{bmatrix} = \begin{bmatrix} E_{n-1} & A_{n-1}^{-1}\alpha_{n-1,1} \\ \alpha_{1,n-1} & a_{nn} \end{bmatrix} $$

$$ \det A_{n-1}^{-1} \det A_n = a_{nn} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1} $$

* Do đó,

$$ \det A \neq 0 \iff a_{nn} - \alpha_{1,n-1} A_{n-1}^{-1} \alpha_{n-1,1} \neq 0 $$

---

## PP viền quanh

* Điều kiện: $\quad \exists A_k^{-1} \quad \forall k = \overline{1, n}.$
* A khả nghịch bất kỳ:

$$ x^t(A^t A)x = \langle x, A^t Ax \rangle = \langle Ax, Ax \rangle \ge 0 $$

$$ "=" \iff Ax = 0 \iff x = 0. $$

* Ma trận $M = A^t A$ thỏa mãn điều kiện thực hiện phương pháp. Khi đó

$$ A^{-1} = M^{-1} A^t $$

---

## CÁC PP LẶP

* Phương pháp lặp Jacobi (A chéo trội)

$$ X = (I - TA)X + T $$

* Phương pháp lặp Gauss - Seildel

$$ A = D - L - U $$

$$ X = (I - T(D - L - U))X + T = (TL + TU)X + T $$

$$ X = (I - TL)^{-1} TUX + T = (D - L)^{-1} UX + T $$

---

## PP LẶP NEWTON

* Ý tưởng: Tìm nghịch đảo của số thực bằng pp tiếp tuyến

$$ ax = 1 \iff f(x) = \frac{1}{x} - a = 0, \quad a > 0 $$

$$ f'(x) = \frac{-1}{x^2} < 0; \quad f''(x) = \frac{2}{x^3} > 0 \quad \forall x > 0 $$

$$ x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} = x_n - \frac{\frac{1}{x_n} - a}{\frac{-1}{x_n^2}} = x_n - x_n(ax_n - 1) $$

---

## PP LẶP NEWTON

* Ý tưởng: mô phỏng pp tìm nghịch đảo của số thực

$$ x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} = x_n - \frac{\frac{1}{x_n} - a}{\frac{-1}{x_n^2}} = x_n - x_n(ax_n - 1) $$

$$ X_{n+1} = X_n - X_n(AX_n - E) $$

---

## PP LẶP NEWTON

* Sự hội tụ của phương pháp:

$$ X_{n+1} = X_n + X_n(E - AX_n) $$

$$ G_n = E - AX_n \implies X_{n+1} = X_n + X_n G_n = X_n(E + G_n) $$

$$ \begin{aligned} G_{n+1} &= E - AX_{n+1} = E - AX_n(E + G_n) \\ &= E - AX_n - AX_n G_n = G_n - AX_n G_n \\ &= G_n^2 = G_{n-1}^{2^2} = \dots = G_{n-k}^{2^{k+1}} = G_0^{2^{n+1}} \end{aligned} $$

$$ \implies \|G_{n+1}\| \le \|G_0\|^{2^{n+1}} \xrightarrow{n \to \infty} 0 \iff \|G_0\| < 1 $$

---

## PP LẶP NEWTON

* Sai số:

$$ G_n = AX_n - E = A(X_n - A^{-1}) \implies X_n - A^{-1} = A^{-1} G_n $$

$$ \|X_n - A^{-1}\| \le \|A^{-1}\| \|G_n\| $$

$$ (AX_0)^{-1} = (I - G_0)^{-1} = \sum_{k=0}^{\infty} G_0^k $$

$$ \iff A^{-1} = X_0(I - G_0)^{-1} = X_0 \sum_{k=0}^{\infty} G_0^k $$

$$ \implies \|A^{-1}\| \le \frac{\|X_0\|}{1 - q} \implies \boxed{ \|X_n - A^{-1}\| \le \frac{q^{2^n}}{1 - q} \|X_0\| }, \quad q = \|G_0\| < 1 $$