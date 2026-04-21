#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Hàm in ma trận
void printMatrix(const vector<vector<double>>& mat) {
    for (const auto& row : mat) {
        for (double val : row) {
            cout << setw(10) << fixed << setprecision(4) << (abs(val) < 1e-9 ? 0.0 : val) << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

int main() {
    // 1. ĐỌC FILE DỮ LIỆU
    string filename = "test.txt" ; // Đặt tên file của bạn ở đây
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Loi: Khong the mo duoc file: " << filename << "'. Kiem tra lai duong dan!\n";
        return 1;
    }

    vector<vector<double>> Aug;
    string line;

    // Đọc từng dòng trong file
    while (getline(file, line)) {
        if (line.empty()) continue; // Bỏ qua nếu có dòng trống

        stringstream ss(line);
        vector<double> row;
        double val;

        // Tách các số cách nhau bởi khoảng trắng trên dòng đó
        while (ss >> val) {
            row.push_back(val);
        }

        if (!row.empty()) {
            Aug.push_back(row);
        }
    }
    file.close();

    // 2. LẤY KÍCH THƯỚC TỰ ĐỘNG
    if (Aug.empty()) {
        cerr << "Lỗi: File rỗng hoặc không chứa dữ liệu hợp lệ!\n";
        return 1;
    }

    int m = Aug.size();          // Số hàng = số dòng đọc được
    int n = Aug[0].size();       // Số cột = số phần tử của dòng đầu tiên
    
    // Giả định ma trận A là ma trận vuông m x m (phổ biến nhất khi giải Ax=B)
    int colsA = m;               
    int colsB = n - colsA;       // Phần còn lại là số cột của B

    cout << "Da doc ma tran kich thuoc " << m << " x " << n << " tu file.\n";
    cout << "--- MA TRAN BAN DAU [A|B] ---\n";
    printMatrix(Aug);

    // ==========================================
    // 3. THUẬT TOÁN GAUSS - QUY TRÌNH THUẬN 
    // ==========================================
    vector<int> ind(m, -1);
    int i = 0, j = 0;
    int step = 1;

    // Đã xóa phần in MA TRAN KHOI TAO ở đây theo ý bạn
    cout << "--- LICH SU CHON PHAN TU KHOA (PIVOT) ---\n";
    
    while (i < m && j < n) {
        if (abs(Aug[i][j]) < 1e-9) {
            bool found = false;
            for (int t = i + 1; t < m; ++t) {
                if (abs(Aug[t][j]) > 1e-9) {
                    swap(Aug[i], Aug[t]);
                    found = true;
                    break;
                }
            }
            if (!found) { 
                j++;
                continue; 
            }
        }

        ind[i] = j;
        cout << "Lan lap " << step << ": Chon pivot a[" << i+1 << "][" << j+1 << "] = " 
             << fixed << setprecision(4) << Aug[i][j] << "\n";

        if (i == m - 1) {
            // Nếu là hàng cuối, in luôn ra nếu đây vô tình là lần lặp đầu tiên (ma trận 1 hàng)
            if (step == 1) {
                cout << "\n--- MA TRAN SAU LAN LAP DAU TIEN ---\n";
                printMatrix(Aug);
            }
            break; 
        }

        // Khử các phần tử bên dưới pivot
        for (int k = i + 1; k < m; ++k) {
            double factor = Aug[k][j] / Aug[i][j];
            for (int c = j; c < n; ++c) {
                Aug[k][c] -= factor * Aug[i][c];
            }
        }

        // ---> IN MA TRẬN SAU LẦN LẶP ĐẦU TIÊN TẠI ĐÂY <---
        if (step == 1) {
            cout << "\n--- MA TRAN SAU LAN LAP DAU TIEN ---\n";
            printMatrix(Aug);
            cout << "-----------------------------------------\n";
        }

        step++; // Tăng biến đếm sau khi đã hoàn tất mọi việc của lần lặp hiện tại

        if (j == n - 1) break; 
        i++;
        j++;
    }

    cout << "\n--- MA TRAN SAU LAN LAP CUOI (KET THUC QUY TRINH THUAN) ---\n";
    printMatrix(Aug);
    // ==========================================
    // 4. THUẬT TOÁN GAUSS - QUY TRÌNH NGHỊCH
    // ==========================================
    vector<vector<double>> X(colsA, vector<double>(colsB, 0.0));

    for (int r = m - 1; r >= 0; --r) {
        int pivot_col = ind[r];
        if (pivot_col != -1 && pivot_col < colsA) {
            for (int cb = 0; cb < colsB; ++cb) {
                double sum = 0;
                for (int k = pivot_col + 1; k < colsA; ++k) {
                    sum += Aug[r][k] * X[k][cb];
                }
                X[pivot_col][cb] = (Aug[r][colsA + cb] - sum) / Aug[r][pivot_col];
            }
        }
    }

    cout << "--- MA TRAN NGHIEM X ---\n";
    printMatrix(X);

    return 0;
}