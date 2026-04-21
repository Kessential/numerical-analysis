#include <iostream>
#include <cmath>
#include <functional>
#include <iomanip> 

using namespace std;

// Hàm thực hiện thuật toán kết hợp
double combined(
    function<double(double)> f, 
    function<double(double)> df, 
    function<double(double)> d2f, 
    double a, double b, 
    double epsilon = 1e-5, 
    int max_iter = 100) 
{
    // Bước 1: Kiểm tra điều kiện có nghiệm
    if (f(a) * f(b) >= 0) {
        cout << "Loi: [a, b] khong chua nghiem\n";
        return NAN; // Trả về giá trị Not-A-Number báo lỗi
    }

    double x_tt, x_dc;

    // Bước 2: Chọn điểm xuất phát theo điều kiện Fourier f(x) * f''(x) > 0
    if (f(a) * d2f(a) > 0) {
        x_tt = a; 
        x_dc = b;   
    } else if (f(b) * d2f(b) > 0) {
        x_tt = b;
        x_dc = a;
    } else {
        cout << "Loi: f(x)*f''(x) > 0 tai hai dau mut\n";
        return NAN;
    }

    double x_dc_new = x_dc;

    // In tiêu đề bảng
    cout << fixed << setprecision(8);
    cout << string(60, '-') << "\n";
    cout << setw(3) << "k" 
         << setw(15) << "x_{2k+2}" 
         << setw(15) << "x_{2k+3}" 
         << setw(15) << "Sai so\n"; 
    cout << string(60, '-') << "\n";

    // Bước 3: Vòng lặp
    for (int i = 0; i < max_iter; ++i) {
        
        // 1. Tính điểm Tiếp tuyến mới
        double df_x_tt = df(x_tt);
        if (abs(df_x_tt) < 1e-9) {
            cout << "Loi: Dao ham bang 0, dung thuat toan.\n";
            return NAN;
        }
        double x_tt_new = x_tt - f(x_tt) / df_x_tt;

        // 2. Tính điểm Dây cung mới
        double mau_so = f(x_dc) - f(x_tt);
        if (abs(mau_so) < 1e-9) {
            cout << "Loi: Mau so day cung bang 0, dung thuat toan.\n";
            return NAN;
        }
        double tu_so = f(x_dc) * (x_dc - x_tt);
        x_dc_new = x_dc - tu_so / mau_so;

        // 3. Đánh giá sai số
        double sai_so = abs(x_dc_new - x_dc);

        // In kết quả của lần lặp hiện tại
        cout << setw(3) << i + 1 
             << setw(15) << x_tt_new 
             << setw(15) << x_dc_new 
             << setw(15) << sai_so << "\n";

        // Kiểm tra điều kiện dừng
        if (sai_so <= epsilon) {
            cout << string(60, '-') << "\n";
            cout << "=> Hoi tu sau " << i + 1 << " buoc lap.\n";
            return x_dc_new;
        }

        // Cập nhật giá trị
        x_tt = x_tt_new;
        x_dc = x_dc_new;
    }

    cout << "Canh bao: Dat den so vong lap toi da ma chua hoi tu.\n";
    return x_dc_new;
}

int main() {
    // Giai phuong trinh x^3 - 2x - 5 = 0 tren khoang [2, 3]
    auto f = [](double x) { return pow(x, 3) - 2 * x - 5; };
    auto df = [](double x) { return 3 * pow(x, 2) - 2; };
    auto d2f = [](double x) { return 6 * x; };

    double a = 2.0;
    double b = 3.0;
    double epsilon = 1e-5;

    double x = combined(f, df, d2f, a, b, epsilon);
    
    if (!isnan(x)) {
        cout << "\nNghiem gan dung cuoi cung: " << x << "\n";
    }

    return 0;
}