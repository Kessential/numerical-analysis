Vấn đề gốc: hệ $p,q$ ở Bước 2.2 (giải từ $z_r-py_r+qx_r^{(k)}=0$, $z_s-py_s+qx_s^{(k)}=0$) về mặt toán học luôn hội tụ về 2 giá trị riêng có module lớn nhất của $A$ — bất kể 2 giá trị đó có bằng module nhau hay không. Lý do: quan hệ truy hồi $A^{k+2}x-pA^{k+1}x+qA^kx\approx0$ với $p=\lambda_1+\lambda_2,\ q=\lambda_1\lambda_2$ chỉ cần $|\lambda_3|$ nhỏ hơn hẳn $|\lambda_2|$ để bỏ qua phần dư — không cần giả thiết $|\lambda_1|=|\lambda_2|$ ở đâu cả trong phép suy diễn.

Hệ quả: ngay cả khi thực chất đang ở TH1 thật sự ($|\lambda_1|>|\lambda_2|$ rõ rệt, không hề bằng nhau), hệ $p,q$ này vẫn "ăn gian" tìm ra đúng $\lambda_1,\lambda_2$ — chỉ là nó không biết/không cần biết rằng $|\lambda_1|\ne|\lambda_2|$.

Tốc độ hội tụ khác nhau:
- Bước 2.1 (tỉ số đơn giản) hội tụ với tốc độ $|\lambda_2/\lambda_1|^k$ (sai số do lẫn thành phần $v_2$).
- Bước 2.2 (hệ $p,q$) hội tụ với tốc độ $|\lambda_3/\lambda_2|^k$ (sai số do lẫn thành phần $v_3$).

Hai tốc độ này độc lập nhau. Nếu $|\lambda_3/\lambda_2| < |\lambda_2/\lambda_1|$ (khoảng cách $\lambda_2\to\lambda_3$ tỉ đối lớn hơn khoảng cách $\lambda_1\to\lambda_2$), thì Bước 2.2 ổn định trước Bước 2.1 — dù đang ở TH1. Ví dụ đã kiểm bằng số: $\lambda=(5,-4,3,-2,1)$ có $|\lambda_2/\lambda_1|=4/5=0.8$ nhưng $|\lambda_3/\lambda_2|=3/4=0.75<0.8$, nên hệ $p,q$ hội tụ nhanh hơn (ra đúng $\lambda_a=5,\lambda_b=-4$) trước khi Bước 2.1 kịp ổn định.

Vì sao không được gán luôn là TH2: TH2 nghĩa là $\lambda_1=-\lambda_2$ (tức $p=\lambda_1+\lambda_2=0$). Trong ví dụ trên $p=5+(-4)=1\ne0$ — hai giá trị này không đối nhau, chỉ tình cờ là 2 giá trị module lớn nhất. Nếu cứ thấy $\Delta\ge0$ (2 nghiệm thực) mà kết luận ngay "TH2" thì sẽ sai (đây chính là lỗi tôi từng bắt được lúc code báo nhầm TH2 cho ma trận này).

Cách xử lý: kiểm tra thêm $p\approx0$ hay không.
- $p\approx0$: đúng là $\lambda_1=-\lambda_2$ → TH2 thật.
- $p\not\approx0$: chỉ là "vé số" — hệ $p,q$ tìm ra top-2 eigenvalue sớm hơn dự kiến, nhưng bản chất vẫn là TH1. Lúc này lấy nghiệm có $|\cdot|$ lớn hơn làm $\lambda_1$ (nghiệm nhỏ hơn — $\lambda_2$ thật — bị bỏ, không trả về, vì kết quả TH1 chỉ cần 1 cặp $(\lambda_1,v_1)$).