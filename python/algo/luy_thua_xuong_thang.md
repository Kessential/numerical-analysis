# PP lũy thừa tìm giá trị riêng trội, PP xuống thang tìm giá trị riêng trội tiếp theo

## 1. PP lũy thừa tìm giá trị riêng trội (đủ 3 trường hợp)

**Input:** ma trận vuông $A_{n\times n}$ thực, véctơ khởi đầu $x_0\ne0$ (mặc định véctơ toàn 1), sai số $\varepsilon>0$, số lần lặp tối đa $N$.

**Output:** tự động phát hiện 1 trong 2 dạng:
* **TH1** ($|\lambda_1|>|\lambda_2|$): $\lambda_1$ thực và véctơ riêng $v_1$.
* **TH2/TH3** ($|\lambda_1|=|\lambda_2|>|\lambda_3|$): cặp $\lambda_{1,2}$ (TH2: thực đối nhau, TH3: phức liên hợp) và 2 véctơ riêng tương ứng.

**Bước 1 — Chuẩn hoá:** $x^{(0)} = x_0/\|x_0\|_\infty$.

**Bước 2 — Với $k=0,1,\dots,N$:** đặt $y:=Ax^{(k)}$, $z:=Ay$ (vai trò $A^{k+1}x_0$, $A^{k+2}x_0$).

* **Bước 2.1 — TH1:** gọi $I=\{i:x_i^{(k)}\ne0\}$, $m=|I|$; tính $r_i=y_i/x_i^{(k)}$ với mọi $i\in I$. Nếu các $r_i$ xấp xỉ bằng nhau và giá trị
$$
\lambda_1 := \frac{1}{m}\sum_{i\in I} r_i
$$
ổn định qua 2 vòng lặp liên tiếp: $v_1:=y$. Dừng thuật toán, trả về $(\lambda_1,v_1)$.

* **Bước 2.2 — TH2/TH3 (nếu Bước 2.1 chưa thoả):** chọn 2 toạ độ $r\ne s$ ($|x_i^{(k)}|$ lớn nhất, để hệ dưới đây không suy biến), giải $p,q$ từ hệ
$$
z_r-py_r+qx_r^{(k)}=0,\qquad z_s-py_s+qx_s^{(k)}=0
$$
Nếu $(p,q)$ ổn định qua 2 vòng lặp liên tiếp: $\Delta=p^2-4q$, $\lambda_{a, b}=(p\pm\sqrt\Delta)/2$, véctơ riêng tương ứng $v_{(\lambda_b)}\approx y-\lambda_ax^{(k)}$, $v_{(\lambda_a)}\approx y-\lambda_bx^{(k)}$.
  * $\Delta\ge0$, $p\approx0$: **TH2** ($\lambda_1=-\lambda_2$, đặt $\lambda_1:=\lambda_a,\lambda_2:=\lambda_b$). Dừng, trả về $(\lambda_a,\lambda_b)$ và $(v_{(\lambda_a)},v_{(\lambda_b)})$ đã tính ở trên.
  * $\Delta\ge0$, $p\not\approx0$: Dừng, trả về dạng **TH1**: $\lambda_1:=\lambda_a$ nếu $|\lambda_a|\ge|\lambda_b|$ (ngược lại $\lambda_1:=\lambda_b$), $v_1:=v_{(\lambda_1)}$ đã tính ở trên.
  * $\Delta<0$: **TH3** (phức liên hợp, $\lambda_{a,b}=\alpha\pm i\beta$). Dừng, trả về $(\lambda_a,\lambda_b)$ và $(v_{(\lambda_a)},v_{(\lambda_b)})$ đã tính ở trên.

* **Bước 2.3:** nếu chưa hội tụ nhánh nào: $x^{(k+1)}\leftarrow y/\|y\|_\infty$, quay lại Bước 2 với $k+1$.

**Bước 3 — Không hội tụ:** hết $N$ vòng lặp: cảnh báo (có thể >2 gtr trội cùng module, hoặc $|\lambda_1|,|\lambda_2|$ quá gần $|\lambda_3|$).

---

## 2. PP xuống thang tìm giá trị riêng trội tiếp theo

Giả thiết đã có $(\lambda_1,v_1)$ từ TH1 (mục 1). Dùng **Cách chọn 2** của slide (không cần véctơ riêng trái $w_1$).

**Input:** $A_{n\times n}$, $(\lambda_1,v_1)$.

**Output:** $(\lambda_2,v_2)$ — giá trị riêng trội tiếp theo và véctơ riêng tương ứng.

**Bước 1 — Chuẩn hoá $v_1$:** $s:=$ toạ độ $i$ làm $|v_{1,i}|$ lớn nhất, rồi $v_1\leftarrow v_1/v_{1,s}$.

**Bước 2 — Ma trận xuống thang:** với $a_s=$ hàng $s$ của $A$: $B=A-v_1a_s$.

**Bước 3 — Tìm gtr trội của $B$:** áp dụng PP lũy thừa (mục 1) cho $B$. Phổ của $B$ là $\{0,\lambda_2,\dots,\lambda_n\}$, nên nếu hội tụ TH1: kết quả $(\lambda_2,u_2)$.

**Bước 4 — Khôi phục véctơ riêng của $A$:**
$$
v_2 = (\lambda_1-\lambda_2)\,u_2 - (a_s\cdot u_2)\,v_1
$$

**Bước 5 — Tổng quát hoá: tìm dãy $\lambda_1,\lambda_2,\lambda_3,\dots$:** đặt $A^{(1)}:=A$, $v_1^{(1)}:=v_1$, $a_1:=a_s$ (khớp Bước 1–2). Với $k=1,2,\dots$ (đã có $A^{(k)}$):
$$
(\lambda_k,\,v_k^{(k)}) := \text{PP lũy thừa TH1 cho } A^{(k)}
$$
$$
s_k := \text{toạ độ } i \text{ làm } |v_{k,i}^{(k)}| \text{ lớn nhất}, \qquad v_k^{(k)} \leftarrow v_k^{(k)}/v_{k,s_k}^{(k)}, \qquad a_k := \text{hàng } s_k \text{ của } A^{(k)}
$$
$$
A^{(k+1)} := A^{(k)} - v_k^{(k)}a_k
$$
(với $k=1$: $A^{(2)}=B$ — đúng Bước 2; với $k=2$: $v_2^{(2)}=u_2$, $a_2=b_{s_2}$, $A^{(3)}=C$.)

Khôi phục $v_k^{(k)}$ về không gian gốc của $A$ bằng công thức Bước 4, áp dụng lần lượt $j=k-1,k-2,\dots,1$:
$$
v_k^{(j)} := (\lambda_j-\lambda_k)\,v_k^{(j+1)} - (a_j\cdot v_k^{(j+1)})\,v_j^{(j)}
$$
Kết quả cuối $v_k:=v_k^{(1)}$ là véctơ riêng của $A$ ứng với $\lambda_k$.

Nếu PP lũy thừa cho $A^{(k)}$ hội tụ ra **cặp** (TH2/TH3) thay vì TH1: dừng lại — xuống thang cho 1 cặp giá trị riêng cần suy rộng sang rank-2, không cài đặt ở đây.
