## Checklist kiểm tra file

### 📁 matrix/

#### cholesky/

- [x] cholesky_2file.py
- [x] cholesky.py

vấn đề: chưa hiện ma trận M và D khi A ko đối xứng và kết quả
kiểm tra lại với trường hợp B nhiều cột

#### jacobi/

- [x] jacobi_2file.py
- [x] jacobi.py

chưa tính lambda cụ thể;
trường hợp phân kỳ hiện sai số lần 1 và 2 là vô cùng; a[i][i] chỗ tính nghịch đảo = 0 tại i nào, chĩ rõ ra

#### lu/

- [x] lu_2file.py
- [x] lu.py

#### nghich_dao/

- [x] nghich_dao_cholesky.py
- [x] nghich_dao_gauss_jordan.py
- [x] nghich_dao_gauss_seidel.py
- [x] nghich_dao_gauss.py
- [x] nghich_dao_jacobi.py
- [x] nghich_dao_newton.py
- [x] nghich_dao_phu_hop.py
- [x] nghich_dao_vien_quanh.py

test newton với trường hợp chuẩn 2

lỗi UI khi phân kỳ

viền quanh hiển thị các chỉ số của quá trình lặp của 2 lần đầu và cuối, kiểm tra cả ma trận

#### svd/

- [x] svd_gia_tri_ky_di_lon_nhat.py
- [x] svd_gia_tri_ky_di.py
- [x] svd_khai_trien.py
- [x] svd_so_dieu_kien.py
- [x] svd_vecto_trai.py
- [x] svd_xap_xi_anh.py

svd_gia_tri_ky_di_lon_nhat.py: bug đáng chú ý — vẫn in cảnh báo "khong hoi tu" nhưng sau đó vẫn in tiêu đề "=== ... (hoi tu sau 500 lan lap) ===" kèm
  σ₁=5.521358 (sai ~10% so với giá trị thật 5.000000) — không có guard chặn việc dùng kết quả chưa hội tụ. -> bug

svd_gia_tri_ky_di.py: in cảnh báo không hội tụ sau 500 lần lặp ở bước 1, trả về 0 trị kỳ dị. -> bug ko hiển thị giá trị kì dị với trường hợp đặc biệt; chỉ hiện v_i mà ko hiện u_i; ma trận toàn 0 ko báo lỗi luôn

svd_khai_trien.py: có vấn đề khi sử dụng ma trận dưới; epsilon quá chặt; bug khi ma trận toàn 0

svd_so_dieu_kien.py: ko kiểm tra khi sigma_min = 0 -> cond(A) ra số rất lớn với trường hợp ma trận ở dưới (có phải sai ko???)

svd_vecto_trai.py: phải đâu (thấy giống với file svd_gia_tri_ky_di.py) nên có thể gộp vào cũng đc

Ma trận đặc biệt test các vấn đề ở trên
$$
\begin{pmatrix}
0.771027  & 2.082821  & -0.087476  & 3.846622 \\
0.756989  & -0.097691 &  -0.092846 & -1.157548 \\
-0.656637 & -1.62074  & 4.200717   & 1.628382 \\
0.898679  & 1.510563  & -1.029806  & 1.351133 \\
\end{pmatrix}
$$

#### matrix/ (file gốc)

- [x] danielevsky.py
- [x] lap_don.py
- [x] luy_thua.py
- [x] xuong_thang.py

danielevsky hiện cả ma trận M nữa, kiểm tra phần chuẩn hoá

trường hợp sylvester có 2 trị riêng ra vector riêng bị sai hệ Sylvester bị suy biến

pp luỹ thừa nhận diện sai trường hợp, chưa in đủ quá trình của các lần lặp

xuống thang với trường hợp in 3 gtr trội mà chỉ hiện 1 gtr trội, 2 gtr còn lại trị tuyệt đối bằng nhau

## Báo cáo lỗi AI test
  Bug 1 (nghiêm trọng) — luy_thua.py: TH2 bị nhận nhầm thành TH1

  Vị trí: matrix/luy_thua.py, hàm power_method, dòng 38-97.

  Nguyên nhân: Vòng lặp kiểm tra điều kiện hội tụ TH1 (dòng 50-56, dựa trên tỉ số y1[i]/x[i] ổn định 2 vòng liên tiếp) trước khi kiểm tra hệ phương trình
  TH2/TH3 (dòng 71-95). Khi trị riêng trội thứ 3 (λ3) có |λ3| gần với |λ1|=|λ2| (thường tỉ số |λ3|/|λ1| trong khoảng 0.5–0.95), vector lặp có xác suất đáng
  kể khiến tỉ số từng phần tử "tình cờ" gần bằng nhau ở đúng 2 vòng liên tiếp dù bản chất là TH2 (cặp λ, -λ) — code trả về case=1 (TH1) ngay, bỏ sót hẳn trị
  riêng thứ hai cùng module.

  Tần suất: quét 300–600 ma trận TH2 ngẫu nhiên (λ3/λ1 ∈ [0.3, 0.95]) → tỷ lệ nhận nhầm ~3–11% tùy vùng tỉ số. Chiều ngược lại (TH1 thật hoặc TH3 thật bị
  nhận nhầm) không xảy ra (0/200 mỗi loại) — lỗi chỉ một chiều.

  Repro: test_luythua_th2_nham_th1_1.txt (σ thật 4, -4, 2.6 → báo TH1, lambda=4.000001), test_luythua_th2_nham_th1_2.txt (σ thật 4,-4,3.4 → báo TH1,
  lambda=-4.000000), test_luythua_th2_nham_th1_3.txt (σ thật 2.113,-2.113,1.103 → báo TH1, lambda=-2.113100). Chạy với eps=1e-6, đối chiếu bằng
  numpy.linalg.eigvals in cuối mỗi output.

  Hướng sửa gợi ý: trước khi return TH1 ở dòng 53, kiểm tra thêm xem hệ p,q (TH2/TH3) có đồng thời hội tụ về nghiệm gần đối nhau hay không — hoặc đảo thứ tự
  ưu tiên, chạy đủ vài vòng dự phòng trước khi chốt case.

  Bug 2 (nghiêm trọng hơn) — xuong_thang.py: cùng lỗi, và im lặng hoàn toàn khi so_gtr=1

  Vị trí: matrix/xuong_thang.py, hàm power_method, dòng 26-99 — bản copy gần như y hệt hàm ở Bug 1.

  Biểu hiện:
  - Với so_gtr≥2: bug xảy ra ở bước 1 (~11% ma trận TH2 test được), in lambda_1 = ... hoi tu sau N lan lap không cảnh báo gì; nhờ cơ chế xuống thang
  (deflation) không phụ thuộc nhãn case, bước 2 thường vẫn tìm ra trị riêng còn lại đúng — kết quả cuối cùng thường vẫn đúng (33/33 trường hợp test khớp
  numpy).
  - Với so_gtr=1: vòng lặp bước 2 trở đi (dòng 179) không chạy, nên không có bất kỳ cảnh báo nào — chương trình in trị riêng bị nhận nhầm như thể là trị
  riêng trội duy nhất, hoàn toàn im lặng về việc thực chất có 1 trị riêng khác cùng module. Đây là tình huống tệ nhất vì người dùng không có cách nào biết
  kết quả không đầy đủ.

  Repro: test_xuongthang_th2_nham_th1.txt, chạy với eps=1e-6:
  - so_gtr=2 → bước 1 báo nhầm lambda_1=4.000001 (không cảnh báo), bước 2 tự sửa lại đúng lambda_2=-4.000000.
  - so_gtr=1 → chỉ in lambda_1=4.000001 như kết quả cuối cùng, không có dấu hiệu gì cho biết còn thiếu -4.

  Hướng sửa gợi ý: cùng gốc với Bug 1 (sửa power_method dùng chung); ngoài ra nên cân nhắc luôn chạy thêm 1 bước kiểm tra TH2/TH3 dự phòng bất kể so_gtr
  được yêu cầu bao nhiêu, để ít nhất cảnh báo cho người dùng.

  Bug 3 (trung bình) — svd_gia_tri_ky_di_lon_nhat.py: header tự mâu thuẫn với cảnh báo không hội tụ

  Vị trí: matrix/svd/svd_gia_tri_ky_di_lon_nhat.py, dòng 63-64 (cảnh báo) và dòng 99 (header) trong main().

  Mô tả: Khi power_method_symmetric không hội tụ sau max_iter=500 (2 giá trị kỳ di lớn nhất quá gần nhau), code in cảnh báo "Canh bao: PP luy thua khong hoi
  tu sau 500 lan lap." (dòng 63-64) nhưng không return/dừng — vẫn tiếp tục in "=== Gia tri ky di lon nhat (hoi tu sau {it} lan lap) ===" (dòng 99) ngay sau
  đó, tự mâu thuẫn với cảnh báo vừa in. Giá trị sigma_1 trả về khi đó chỉ là ước lượng Rayleigh tạm thời (dòng 51-52 trong power_method_symmetric), không
  đáng tin.

  Số liệu cụ thể: test với test_svd_gan_bang_nhau.txt (σ thật = 5, 4.9999, 1) → in cảnh báo "khong hoi tu" NHƯNG NGAY SAU ĐÓ vẫn ghi "hoi tu sau 500 lan
  lap", trả về sigma_1 = 5.521358 — sai ~10.4% so với giá trị thật 5.000000 (đối chiếu numpy.linalg.svd).

  Hướng sửa gợi ý: đổi câu chữ header thành có điều kiện (vd. "khong hoi tu, ket qua uoc luong tam thoi") khi converged=False, thay vì luôn ghi "hoi tu".

  Bug 4 (nhẹ, chỉ là thông báo gây hiểu lầm) — svd_khai_trien.py: đổ lỗi sai nguyên nhân sai số tái tạo

  Vị trí: matrix/svd/svd_khai_trien.py, dòng 145-152.

  Mô tả: Khi ma trận A hạng thiếu (một số σ ~0 bị bỏ qua đúng theo thiết kế ở nhánh lam<=eps, không phải do không hội tụ), sai số tái tạo ||A-U.Sigma.V^T||
  có thể vượt ngưỡng tol=max(eps*100,1e-8) chỉ vì tích lũy sai số làm tròn — nhưng message lại luôn ghi cứng "co buoc PP luy thua o Phan 2 chua hoi tu" dù
  thực chất không có bước nào không hội tụ.

  Repro: test_svd_rank_deficient.txt, eps=1e-9 → dừng đúng ở r=2 (hạng thật), nhưng sai số 5.824e-7 > tol=1e-7 → in nhầm thông báo "chua hoi tu".