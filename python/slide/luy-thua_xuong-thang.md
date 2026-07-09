
## PP lũy thừa tìm giá trị riêng trội
## PP xuống thang tìm giá trị riêng trội tiếp theo

**Hà Nội, 2026**

---

## Phương pháp lặp tìm giá trị riêng trội

## Giá trị riêng trội

*   Giả sử ma trận $A$ vuông cỡ $n$ thực có các giá trị riêng khác nhau xếp theo thứ tự:
    $$|\lambda_1| \ge |\lambda_2| \ge \dots \ge |\lambda_s|$$
*   Các véc tơ riêng ứng với các giá trị riêng:
    $$Av_i = \lambda_i v_i, \quad i = \overline{1, s}.$$
*   Khi đó:
    *   $|\lambda_1| > |\lambda_2| \Rightarrow \lambda_1$ là giá trị riêng trội
    *   $|\lambda_1| = |\lambda_2| > |\lambda_3| \Rightarrow \lambda_1, \lambda_2$ là các giá trị riêng trội
    *   $\dots$

---

## PP lũy thừa tìm GTR

*   Giả sử
    $$x = a_1 v_1 + a_2 v_2 + \dots + a_s v_s$$
    $$\Rightarrow Ax = a_1 A v_1 + a_2 A v_2 + \dots + a_s A v_s$$
    $$= a_1 \lambda_1 v_1 + a_2 \lambda_2 v_2 + \dots + a_s \lambda_s v_s$$
    $$\Rightarrow A^2 x = a_1 \lambda_1^2 v_1 + a_2 \lambda_2^2 v_2 + \dots + a_s \lambda_s^2 v_s$$
    $$\Rightarrow A^k x = a_1 \lambda_1^k v_1 + a_2 \lambda_2^k v_2 + \dots + a_s \lambda_s^k v_s$$
    $$\Rightarrow \frac{A^k x}{\lambda_1^k} = a_1 v_1 + a_2 \frac{\lambda_2^k}{\lambda_1^k} v_2 + \dots + a_s \frac{\lambda_s^k}{\lambda_1^k} v_s$$

---

## PP lũy thừa tìm GTR

*   Trường hợp $|\lambda_1| > |\lambda_2|$

    $$\frac{A^k x}{\lambda_1^k} = a_1 v_1 + a_2 \frac{\lambda_2^k}{\lambda_1^k} v_2 + \dots + a_s \frac{\lambda_s^k}{\lambda_1^k} v_s$$
    
    $$\lim_{k \to \infty} \frac{A^k x}{\lambda_1^k} = \lim_{k \to \infty} \left[ a_1 v_1 + a_2 \frac{\lambda_2^k}{\lambda_1^k} v_2 + \dots + a_s \frac{\lambda_s^k}{\lambda_1^k} v_s \right] = a_1 v_1$$
    
    $$\Rightarrow \frac{A^k x}{\lambda_1^k} \approx a_1 v_1 \Rightarrow \frac{A^{k+1} x}{\lambda_1^k} = A \left( \frac{A^k x}{\lambda_1^k} \right) \approx \lambda_1 (a_1 v_1)$$
    
    $$\Rightarrow \lambda_1 \approx \frac{(A^{k+1} x)_i}{(A^k x)_i} \quad \forall i = \overline{1, n}.$$

---

## PP lũy thừa tìm GTR

*   Trường hợp $|\lambda_1| = |\lambda_2| > |\lambda_3|, \lambda_1 = -\lambda_2$

    $$\frac{A^k x}{\lambda_1^k} = a_1 v_1 + a_2 \frac{\lambda_2^k}{\lambda_1^k} v_2 + \dots + a_s \frac{\lambda_s^k}{\lambda_1^k} v_s$$
    
    $$\lim_{n \to \infty} \frac{A^{2n} x}{\lambda_1^{2n}} = \lim_{n \to \infty} \left[ a_1 v_1 + a_2 (-1)^{2n} v_2 + a_3 \frac{\lambda_3^{2n}}{\lambda_1^{2n}} v_3 + \dots + a_s \frac{\lambda_s^{2n}}{\lambda_1^{2n}} v_s \right]$$
    $$= a_1 v_1 + a_2 v_2$$
    
    $$\Rightarrow \frac{A^{2n} x}{\lambda_1^{2n}} \approx a_1 v_1 + a_2 v_2 \Rightarrow A^2 \left( \frac{A^{2n} x}{\lambda_1^{2n}} \right) = \lambda_1^2 \left( \frac{A^{2n} x}{\lambda_1^{2n}} \right)$$
    
    $$\Rightarrow \lambda_1^2 \approx \frac{(A^{2n+2} x)_i}{(A^{2n} x)_i}$$

---

## PP lũy thừa tìm GTR

*   Trường hợp $|\lambda_1| = |\lambda_2| > |\lambda_3|, \lambda_1 = \overline{\lambda_2}.$

    $$\frac{A^k x}{\lambda_1^k} = a_1 v_1 + a_2 \frac{\lambda_2^k}{\lambda_1^k} v_2 + \dots + a_s \frac{\lambda_s^k}{\lambda_1^k} v_s, \quad \lambda_{1,2} = \alpha \pm i\beta$$
    
    $$\lim_{n \to \infty} \left[ a_3 \frac{\lambda_3^n}{\lambda_1^n} v_3 + \dots + a_s \frac{\lambda_s^n}{\lambda_1^n} v_s \right] = 0$$
    
    $$\Rightarrow \frac{A^n x}{\lambda_1^n} \approx a_1 v_1 + \frac{\lambda_2^n}{\lambda_1^n} a_2 v_2 \Rightarrow A^n x = \lambda_1^n a_1 v_1 + \lambda_2^n a_2 v_2$$
    
    $$\Rightarrow A^{n+2} x - (\lambda_1 + \lambda_2) A^{n+1} x + \lambda_1 \lambda_2 A^n x = 0.$$

---

## PP lũy thừa tìm GTR

*   Trường hợp $|\lambda_1| = |\lambda_2| > |\lambda_3|, \lambda_1 = \overline{\lambda_2}.$
    
    $$p := \lambda_1 + \lambda_2; \quad q = \lambda_1 \lambda_2; \quad (1, -p, q) \neq (0, 0, 0)$$
    
    $$\Rightarrow \begin{cases} (A^{n+2} x)_i - p(A^{n+1} x)_i + q(A^n x)_i = 0 \quad \forall i = \overline{1, m} \\ t^2 - pt + q = 0 \end{cases}$$
    
    $$\Rightarrow \begin{vmatrix} \lambda^2 & \lambda & 1 \\ (A^{n+2} x)_r & (A^{n+1} x)_r & (A^n x)_r \\ (A^{n+2} x)_s & (A^{n+1} x)_s & (A^n x)_s \end{vmatrix} = 0 \Rightarrow \lambda_{1,2}$$

---

## PP lũy thừa tìm GTR, VTR

*   Trường hợp $|\lambda_1| = |\lambda_2| > |\lambda_3|, \lambda_1 = \overline{\lambda_2}.$

    $$A^n x \approx \lambda_1^n a_1 v_1 + \lambda_2^n a_2 v_2$$
    
    $$\Rightarrow \begin{cases} A^{n+1} x - \lambda_1 A^n x \approx a_2 \lambda_2^n (\lambda_2 - \lambda_1) v_2 \\ A^{n+1} x - \lambda_2 A^n x \approx a_1 \lambda_1^n (\lambda_1 - \lambda_2) v_1 \end{cases}$$
    
    $$\Rightarrow \begin{cases} A \left( A^{n+1} x - \lambda_1 A^n x \right) = \lambda_2 \left( A^{n+1} x - \lambda_1 A^n x \right) \\ A \left( A^{n+1} x - \lambda_2 A^n x \right) = \lambda_1 \left( A^{n+1} x - \lambda_2 A^n x \right) \end{cases}$$

---

## PP XUỐNG THANG

*   Cho ma trận vuông cấp n: $A_{n \times n}$
*   Ký hiệu các giá trị riêng và véc-tơ riêng t.ư của $A$: 
    $$\lambda_1, \lambda_2, \dots, \lambda_n, \quad A v_i = \lambda_i v_i.$$
*   Véc-tơ $x$ thỏa mãn: $x^T v_1 = 1$
*   Ma trận xuống thang: $B = A - \lambda_1 v_1 x^T$
*   Giá trị riêng và véc-tơ riêng t.ư. của $B$: 
    $$0, \lambda_2, \dots, \lambda_n, \quad B u_i = \lambda_i u_i.$$
*   Liên hệ giữa các vtr của $A$ và $B$:
    $$v_i = (\lambda_1 - \lambda_i) u_i + \lambda_1 (x^T u_i) v_1$$

---

## PP XUỐNG THANG

*   Cách chọn 1: 
*   $$A v_1 = \lambda_1 v_1 \quad A^T w_1 = \lambda_1 w_1, \quad x = \frac{w_1}{w_1^T v_1} \Rightarrow x^T v_1 = 1.$$
    
    $$B = A - \frac{\lambda_1}{w_1^T v_1} v_1 w_1^T$$
    
    $$\Rightarrow B v_1 = A v_1 - \lambda_1 v_1 = 0$$
    
    $$B v_k = A v_k - \frac{\lambda_1}{w_1^T v_1} v_1 w_1^T v_k = A v_k = \lambda_k v_k$$
    
    $$\left( A^T w_1 = \lambda_1 w_1 \Rightarrow w_1^T A v_k = \lambda_1 w_1^T v_k \Rightarrow w_1^T v_k = 0 \quad \forall \lambda_k \neq \lambda_1 \right)$$

---

## PP XUỐNG THANG

*   Cách chọn 2:
    $$x = \frac{1}{\lambda_1} \begin{bmatrix} a_{s1} & \dots & a_{s,s-1} & a_{s,s} & a_{s,s+1} & \dots & a_{s,n} \end{bmatrix}^T$$
    
    $$x^T v_1 = 1, \quad v_{1,s} = 1.$$
    
    $$B = A - v_1 \begin{bmatrix} a_{s1} & \dots & a_{s,s-1} & a_{s,s} & a_{s,s+1} & \dots & a_{s,n} \end{bmatrix} = \Theta A$$
    
    $$\Theta(s, v_1) = I - \begin{bmatrix} 0 & \dots & 0 & v_1 & 0 & \dots & 0 \end{bmatrix}$$
    
    $$= \begin{bmatrix} e_1 & \dots & e_{s-1} & e_s - v_1 & e_{s+1} & \dots & e_n \end{bmatrix}$$

---

## PP XUỐNG THANG

*   Xét
    $$x \neq 0 \Rightarrow \exists s : x_s \neq 0 \Rightarrow x_s = 1.$$
*   Đặt 
    $$\Theta(s, x) = \begin{bmatrix} e_1 & \dots & e_{s-1} & e_s - x & e_{s+1} & \dots & e_m \end{bmatrix}$$
    $$\Rightarrow \Theta(s, x) z = z - z_s x \quad \forall z \in \mathbb{R}^m$$
    $$\Rightarrow \Theta(s, x) x = x - x_s x = 0$$
*   Ma trận mới:
    $$A^{(2)} = \Theta(s, v_1) A$$
    $$A^{(2)} v_1 = \Theta(s, v_1) \lambda_1 v_1 = 0$$
    $$A^{(2)} \Theta(s, v_1) v_i = \lambda_i \Theta(s, v_1) v_i, \quad i > 1.$$

---

## PP XUỐNG THANG

*   GT: $\quad A v_i = \lambda_i v_i, \quad i = \overline{1, r}$
*   Sau $k - 1$ bước xuống thang, đã tìm được 
    $$(\lambda_i, v_i), \quad i = \overline{1, k - 1}, \quad A^{(k)}$$
*   Tìm gtr, vtr tiếp theo từ ma trận $A^{(k)}$
    
    $$A^{(k)} = \Theta(v_{k-1}, s_{k-1}) A^{(k-1)}$$
    
    $$A^{(k)} v_{k-1} = \Theta(v_{k-1}, s_{k-1}) \lambda_{k-1} v_{k-1} = 0$$
    
    $$A^{(k)} \Theta(v_{k-1}, s_{k-1}) v_i = \Theta(v_{k-1}, s_{k-1}) 0. v_{k-1} = 0, \quad i < k - 1$$
    
    $$A^{(k)} \Theta(v_{k-1}, s_{k-1}) v_i = \lambda_i \Theta(v_{k-1}, s_{k-1}) v_i, \quad i > k - 1.$$