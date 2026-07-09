# Phương pháp Danielevsky tìm đa thức đặc trưng

## Hà Nội, 2026

## Giá trị riêng – Véctơ riêng

$$ \exists x \neq 0, \quad Ax = \lambda x $$

$$ \iff \det(A - \lambda I) = 0 $$

$$ \iff (-1)^n \left[ \lambda^n + p_1\lambda^{n-1} + p_2\lambda^{n-2} + \dots + p_n \right] = 0 $$

---

## Giá trị riêng, véctơ riêng

*   Hai ma trận đồng dạng có cùng bộ giá trị riêng

$$ B = P^{-1}AP $$

$$ \Rightarrow B - \lambda I = P^{-1}AP - \lambda P^{-1}P = P^{-1}(A - \lambda I)P $$

$$ \Rightarrow \det(B - \lambda I) = \det(A - \lambda I) $$

---

## Giá trị riêng, véctơ riêng

*   Véc tơ riêng của 2 ma trận đồng dạng

$$ Bx = \lambda x $$

$$ \iff P^{-1}APx = \lambda x $$

$$ \iff A(Px) = \lambda (Px) $$

---

## Câu hỏi 1

*   Ma trận nào nhận đa thức 

$$ p(\lambda) = (-1)^n \left[ \lambda^n + p_1\lambda^{n-1} + \dots + p_{n-1}\lambda + p_n \right] $$

là đa thức đặc trưng?

---

## Khối Jordan

$$ J_{r,\lambda} = \begin{bmatrix} 
-\lambda & 0 & \cdots & 0 & 0 \\
1 & -\lambda & \cdots & 0 & 0 \\
0 & 1 & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & -\lambda 
\end{bmatrix}; \text{ or } J^t_{r,\lambda} $$

$$ \det(J_{r,\lambda}) = (-1)^r \lambda^r $$

---

## Khối Jordan

*   Từ khối Jordan cỡ r, bỏ đi cột k và hàng 1, ta thu được ma trận có định thức là

$$ \det \begin{bmatrix} 
1 & -\lambda & \cdots & 0 & 0 & \\
0 & 1 & \cdots & 0 & 0 & \\
\vdots & \vdots & \ddots & \vdots & \vdots & 0 \\
0 & 0 & \cdots & 1 & -\lambda & \\
0 & 0 & \cdots & 0 & 1 & \\
& & 0 & & & J_{r-k,\lambda} 
\end{bmatrix} = (-1)^{r-k} \lambda^{r-k} $$

---

## Khối Frobenius dạng 1

$$ C^{(r)} = \begin{bmatrix} 
-p_1 & -p_2 & \cdots & -p_{r-1} & -p_r \\
1 & 0 & \cdots & 0 & 0 \\
0 & 1 & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & 0 
\end{bmatrix}; \quad C^{(1)} = -p_1 $$

---

## Khối Frobenius

$$ \left| C^{(2)} - \lambda I_2 \right| = (-p_1 - \lambda)(-\lambda) + (-1)^2 p_2 = (-1)^2 \left[ \lambda^2 + p_1\lambda + p_2 \right] $$

$$ \left| C^{(r)} - \lambda I_r \right| = -\lambda \left| C^{(r-1)} - \lambda I_{r-1} \right| + (-1)^r p_r $$

$$ = -\lambda(-1)^{r-1} \left[ \lambda^{r-1} + p_1\lambda^{r-2} + \dots + p_{r-1} \right] + (-1)^r p_r $$

$$ = (-1)^r \left[ \lambda^r + p_1\lambda^{r-1} + \dots + p_{r-1}\lambda + p_r \right] $$

---

## Khối Frobenius dạng 2

*   Khối Frobenius

$$ C = \begin{bmatrix} 
0 & 1 & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & 0 \\
0 & 0 & \cdots & 0 & 1 \\
-p_1 & -p_2 & \cdots & -p_{r-1} & -p_r 
\end{bmatrix}; $$

---

## Dạng chuẩn Frobenius của ma trận

*   Dạng chuẩn Frobenius:

$$ P^{-1}AP = \begin{bmatrix} 
F_s & 0 & \cdots & 0 \\
0 & F_{s-1} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & F_1 
\end{bmatrix} = F $$

$$ F_i \in \mathbb{R}^{m_i \times m_i}, \quad \sum_{i=1}^{s} m_i = n $$

---

## Dạng chuẩn Frobenius của ma trận

trong đó

$$ F_i = \begin{bmatrix} 
-p_1 & -p_2 & \cdots & -p_{m_i-1} & -p_{m_i} \\
1 & 0 & \cdots & 0 & 0 \\
0 & 1 & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & 0 
\end{bmatrix} $$

---

## Dạng chuẩn Frobenius của ma trận

*   Đa thức đặc trưng của ma trận dạng chuẩn Frobenius:

$$ \det(A - \lambda I_n) = \det(F_1^{(r_1)} - \lambda I_{r_1}) \dots \det(F_k^{(r_k)} - \lambda I_{r_k}) $$

---

## Câu hỏi 2

Ma trận vuông A bất kỳ có thể có dạng chuẩn Frobenius hay không?

---

## Phương pháp Danielevsky

*   **Ý tưởng:** Dùng các ma trận chuyển cơ sở đưa ma trận A về ma trận tương đương dạng chuẩn Frobenius.

---

## Phương pháp Danielevsky

*   Đưa hàng thứ n của ma trận A về dạng

$$ \begin{matrix} 0 & \cdots & 0 & 1 & 0 \end{matrix} $$

Hoặc

$$ \begin{matrix} 0 & \cdots & 0 & 0 & a_{nn} \end{matrix} $$

Qui ước: $a_k = \begin{bmatrix} a_{k1} & a_{k2} & \cdots & a_{kn} \end{bmatrix}^t$

---

## Phương pháp Danielevsky
Đưa hàng thứ n của ma trận A về dạng

*   **Trường hợp 1:** $\quad a_{n,n-1} = 0, \; a_{nk} \neq 0, \; k < n-1$

Đưa giá trị $a_{nk}$ vào vị trí hàng n cột n-1

$$ C_{k \leftrightarrow n-1} = \begin{bmatrix} e_1 & \cdots & e_{k-1} & e_{n-1} & e_{k+1} & \cdots & e_{n-2} & e_k & e_n \end{bmatrix} $$
$$ \Rightarrow C_{k \leftrightarrow n-1}^{-1} = C_{k \leftrightarrow n-1} $$
$$ A^{(1)} := C_{k \leftrightarrow n-1}A^{(1)}C_{k \leftrightarrow n-1} = \left[ a_{ij}^{(1)} \right] $$
$$ P = PC_{k \leftrightarrow n-1} $$
$$ a_{n,n-1}^{(1)} = e_n^t A e_k = a_{nk} \neq 0 $$

---

## Phương pháp Danielevsky
Đưa hàng thứ n của ma trận A về dạng F

*   **Trường hợp 2:** $\quad a_{n,n-1} \neq 0$

$$ M_1 = \begin{bmatrix} e_1 & \cdots & e_{n-2} & a_n & e_n \end{bmatrix}^t $$

$$ = \begin{bmatrix} 
1 & \cdots & 0 & 0 & 0 \\
\vdots & \ddots & \vdots & \vdots & \vdots \\
0 & \cdots & 1 & 0 & 0 \\
a_{n1} & \cdots & a_{n,n-2} & a_{n,n-1} & a_{nn} \\
0 & \cdots & 0 & 0 & 1 
\end{bmatrix} $$

---

## Phương pháp Danielevsky
Đưa hàng thứ n của ma trận A về dạng F

*   **Trường hợp 2:** $\quad a_{n,n-1} \neq 0$

$$ \det M_1 = a_{n,n-1} \neq 0 $$

$$ M_1^{-1} = \begin{bmatrix} 
1 & \cdots & 0 & 0 & 0 \\
\vdots & \ddots & \vdots & \vdots & \vdots \\
0 & \cdots & 1 & 0 & 0 \\
-\frac{a_{n1}}{a_{n,n-1}} & \cdots & -\frac{a_{n,n-2}}{a_{n,n-1}} & \frac{1}{a_{n,n-1}} & -\frac{a_{nn}}{a_{n,n-1}} \\
0 & \cdots & 0 & 0 & 1 
\end{bmatrix} $$

$$ = \begin{bmatrix} e_1 - \frac{a_{n1}}{a_{n,n-1}}e_{n-1} & \cdots & e_{n-2} - \frac{a_{n,n-2}}{a_{n,n-1}}e_{n-1} & \frac{e_{n-1}}{a_{n,n-1}} & e_n - \frac{a_{nn}}{a_{n,n-1}}e_{n-1} \end{bmatrix}^t $$

---

## Phương pháp Danielevsky
Đưa hàng thứ n của ma trận A về dạng F

*   **Trường hợp 2:** $\quad a_{n,n-1} \neq 0$

$$ A^{(1)} = A; \quad A^{(2)} = M_1A^{(1)}M_1^{-1} = \left[ a_{ij}^{(2)} \right]; \quad P = PM_1^{-1} $$

$$ a_{ij}^{(2)} = e_i^t A \left( e_j - \frac{a_{nj}}{a_{n,n-1}}e_{n-1} \right) = a_{ij} - \frac{a_{nj}a_{i,n-1}}{a_{n,n-1}}, \quad i \neq n-1, \; j \neq n-1 $$

$$ a_{i,n-1}^{(2)} = e_i^t A \left( \frac{1}{a_{n,n-1}}e_{n-1} \right), \quad j \neq n-1 $$

$$ \Rightarrow a_{nj}^{(2)} = \begin{cases} 
0, & j \neq n-1 \\ 
1, & j = n-1 
\end{cases} $$

---

## Phương pháp Danielevsky
Đưa hàng thứ n của ma trận A về dạng

*   **Trường hợp 3:** $\quad a_{nk} = 0, \; \forall k \leq n-1$

Dòng cuối cùng có dạng khối Frobenius cấp 1

$$ A_n = \begin{bmatrix} A_{n-1} & \square \\ \theta & a_{nn} \end{bmatrix}, \quad \theta = \begin{bmatrix} 0 & \cdots & 0 \end{bmatrix} $$

$$ \det(A_n - \lambda I_n) = (a_{nn} - \lambda)\det(A_{n-1} - \lambda I_{n-1}) $$

---

## Phương pháp Danielevsky

Gặp *TH1* đưa về *TH2*.  
Gặp *TH2*: thực hiện tiếp tục cho hàng n-1:

*   Đưa hàng thứ n - 1 của ma trận A về dạng

$$ \begin{bmatrix} 0 & \cdots & 1 & 0 & 0 \end{bmatrix}_n \quad (TH1, 2) $$

hoặc

$$ \begin{bmatrix} 0 & \cdots & 0 & a_{n-1,n-1} & a_{n-1,n} \end{bmatrix}_n \quad (TH3) $$

---

## Phương pháp Danielevsky
Đưa hàng thứ k của ma trận A về dạng F

*   **TH1:** $\quad a_{k,k-1} = 0, \; a_{ks} \neq 0, \; s < k-1$
*   Đưa giá trị $a_{ks}$ vào vị trí hàng k cột k-1

$$ C_{s,k} = \begin{bmatrix} e_1 & \cdots & e_{s-1} & e_k & e_{s+1} & \cdots & e_{k-1} & e_s & e_{k+1} & \cdots & e_n \end{bmatrix} $$
$$ \Rightarrow C_{s,k}^{-1} = C_{s,k} $$
$$ A^{(n-k)} := C_{s,k}A^{(n-k)}C_{s,k} = \left[ a_{ij}^{(n-k)} \right] $$
$$ P = C_{s,k} P $$
$$ a_{k,k-1}^{(n-k)} = e_k^t A^{(n-k)} e_s = a_{ks} \neq 0 $$

---

## Phương pháp Danielevsky
Đưa hàng thứ k của ma trận A về dạng F

*   **Tổng quát, TH2:** $\quad a_{k,k-1} \neq 0$

$$ M_{n-k+1} = \begin{bmatrix} & & & \overset{(k-1)}{} & & \\ e_1 & \cdots & e_{k-2} & a_k & e_k & \cdots & e_n \end{bmatrix}^t $$

$$ M_{n-k+1}^{-1} = \begin{bmatrix} e_1 & \cdots & e_{k-2} & -\frac{a_k^{(n-k)} - \left( a_{k,k-1}^{(n-k)} + 1 \right)e_{k-1}}{a_{k,k-1}} & e_k & \cdots & e_n \end{bmatrix}^t $$

---

## Phương pháp Danielevsky

Giả sử khi xét đến hàng n-m+1, ta gặp *TH3*, tức là

$$ A^{(2)} = M_1 A M_1^{-1} $$
$$ A^{(3)} = M_2 A^{(2)} M_2^{-1} = M_2 M_1 A (M_2 M_1)^{-1} $$
$$ \cdots $$
$$ A^{(m)} = M_{m-1}\dots M_1 A (M_{m-1}\dots M_1)^{-1} $$
$$ P = M_1^{-1} M_2^{-1} \dots M_{m-1}^{-1} $$

---

## Phương pháp Danielevsky

trong đó

$$ A^{(m)} = \begin{bmatrix} 
A_{n-m}^{(m)} & B_{(n-m) \times m} \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix} $$

$F_1^{(m)}$ là khối Frobenius cỡ m.

---

## Phương pháp Danielevsky

Dùng các phép biến đổi đồng dạng đưa ma trận

$$ A^{(m)} = \begin{bmatrix} 
A_{n-m}^{(m)} & B_{(n-m) \times m} \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix} 
\longrightarrow 
\begin{bmatrix} 
A_{n-m}^{(m)} & 0_{(n-m) \times m} \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix} $$

---

## Phương pháp Danielevsky

Đưa cột thứ 1 của B về 0:

$$ S_1 = \begin{bmatrix} 
E_{n-m} & \begin{matrix} 0 & -b_{1,1} & 0 & \cdots & 0 \\ \vdots & \vdots & \vdots & \vdots & \vdots \\ 0 & -b_{n-m,1} & 0 & \cdots & 0 \end{matrix} \\ 
0 & E_m 
\end{bmatrix} \Rightarrow S_1^{-1} = \begin{bmatrix} 
E_{n-m} & \begin{matrix} 0 & b_{1,1} & 0 & \cdots & 0 \\ \vdots & \vdots & \vdots & \vdots & \vdots \\ 0 & b_{n-m,1} & 0 & \cdots & 0 \end{matrix} \\ 
0 & E_m 
\end{bmatrix} $$

$$ A^{(m)} = S_1 A^{(m)} S_1^{-1} = \begin{bmatrix} 
A_{n-m}^{(m)} & A_{n-m}^{(m)} [+] + B + [-]F_1^{(m)} \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix}, \quad P = PS_1^{-1} $$

---

## Phương pháp Danielevsky

Đưa cột thứ 1 của B về 0:

$$ A_{n-m}^{(m)} [+] + B + [-]F_1^{(m)} $$

$$ = A_{n-m}^{(m)} \begin{bmatrix} 
0 & b_{1,1} & 0 & \cdots & 0 \\ 
\vdots & \vdots & \vdots & \vdots & \vdots \\ 
0 & b_{n-m,1} & 0 & \cdots & 0 
\end{bmatrix} + B + \begin{bmatrix} 
0 & -b_{1,1} & 0 & \cdots & 0 \\ 
\vdots & \vdots & \vdots & \vdots & \vdots \\ 
0 & -b_{n-m,1} & 0 & \cdots & 0 
\end{bmatrix} \begin{bmatrix} 
-p_{11} & -p_{12} & \cdots \\ 
1 & 0 & \cdots \\ 
0 & 1 & \cdots \\ 
\vdots & \vdots & 
\end{bmatrix} $$

$$ = \begin{bmatrix} 0 & \cdots \\ \vdots & \cdots \\ 0 & \cdots \end{bmatrix} + \begin{bmatrix} b_{1,1} & \cdots \\ \vdots & \cdots \\ b_{n-m,1} & \cdots \end{bmatrix} + \begin{bmatrix} -b_{1,1} & \cdots \\ \vdots & \cdots \\ -b_{n-m,1} & \cdots \end{bmatrix} = \begin{bmatrix} 0 & \cdots \\ \vdots & \cdots \\ 0 & \cdots \end{bmatrix} $$

---

## Phương pháp Danielevsky

Đưa cột q<m của B về 0 khi có q-1 cột đầu là 0:

$$ S_q = \begin{bmatrix} 
E_{n-m} & 0_{(n-m) \times q} & \begin{matrix} -b_{1,q} \\ \vdots \\ -b_{n-m,q} \end{matrix} & 0_{(n-m) \times (n-m-q-1)} \\ 
0 & & E_m & 
\end{bmatrix} $$

$$ A^{(m)} = S_q A^{(m)} S_q^{-1} = \begin{bmatrix} 
A_{n-m}^{(m)} & A_{n-m}^{(m)} [+] + B + [-]F_1^{(m)} \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix} $$

$$ P = PS_q^{-1} $$

---

## Phương pháp Danielevsky

Đưa cột q của B về 0 khi có q-1 cột đầu là 0:

$$ A_{n-m}^{(m)} \begin{bmatrix} 
0_{(n-m) \times q} & \begin{matrix} b_{1,q} \\ \vdots \\ b_{n-m,q} \end{matrix} & 0_{(n-m) \times (n-m-q-1)} 
\end{bmatrix} = \begin{bmatrix} 
0_{(n-m) \times q} & \cdots 
\end{bmatrix} $$

---

## Phương pháp Danielevsky

Đưa cột q của B về 0 khi có q-1 cột đầu là 0:

$$ B = \begin{bmatrix} 
0_{(n-m) \times (q-1)} & \begin{matrix} b_{1,q} \\ \vdots \\ b_{n-m,1} \end{matrix} & \cdots 
\end{bmatrix} $$

$$ \begin{bmatrix} 
0_{(n-m) \times q} & \begin{matrix} -b_{1,q} \\ \vdots \\ -b_{n-m,q} \end{matrix} & \cdots 
\end{bmatrix} \begin{bmatrix} 
-p_{11} & \cdots & -p_{1,q} & \cdots \\ 
1_{2,1} & \vdots & \vdots & \cdots \\ 
& \vdots & 1_{q+1,q} & \cdots \\ 
& \vdots & \vdots & \cdots 
\end{bmatrix} = \begin{bmatrix} 
0_{(n-m) \times (q-1)} & \begin{matrix} -b_{1,q} \\ \vdots \\ -b_{n-m,q} \end{matrix} & \cdots 
\end{bmatrix} $$

---

## Phương pháp Danielevsky

*   Nếu cột cuối của B khác 0 thì dùng ma trận hoán vị cột n và cột n-m, rồi thực hiện lại thuật toán từ hàng n.

$$ A^{(m)} = C_{n-k \leftrightarrow n}A^{(m)}C_{n-k \leftrightarrow n} $$
$$ P = PC_{n-k \leftrightarrow n} $$

---

## Phương pháp Danielevsky

Sau hữu hạn lần biến đổi ta thu được

$$ A^{(m)} = \begin{bmatrix} 
A_{n-m}^{(m)} & 0 \\ 
0_{m \times (n-m)} & F_1^{(m)} 
\end{bmatrix} $$

Lặp lại quy trình với ma trận $A_{n-m}^{(m)}$

---

## Phương pháp Danielevsky

*   Giải phương trình đa thức đặc trưng
*   Tìm véctơ riêng của từng khối Frobenius của dạng chuẩn F, ghép lại thành véctơ riêng của ma trận dạng chuẩn F, kết hợp với ma trận của phép biến đổi đồng dạng, ta thu được véctơ riêng của ma trận A ban đầu.

---

## Phương pháp Danielevsky

*   Véctơ riêng của khối

$$ u_{(1)} = \begin{bmatrix} \lambda^{m-1} \\ \vdots \\ \lambda \\ 1 \end{bmatrix}, \quad u_{(2)} = \begin{bmatrix} 1 \\ \lambda \\ \vdots \\ \lambda^{m-1} \end{bmatrix} $$

*   Véctơ riêng của ma trận

$$ P \left[ 0_{m_s}^T, \dots, 0_{m_{i+1}}^T, u_{(1)}^T, 0_{m_{i-1}}^T, \dots, 0_{m_1}^T \right]^T $$