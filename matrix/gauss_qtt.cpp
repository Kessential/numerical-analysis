#pragma GCC optimize("O2")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

// ==========================================
// HÀM TIỆN ÍCH: IN MA TRẬN
// ==========================================
void printMatrix(const vector<vector<double>>& mat) {
    for (const auto& row : mat) {
        for (double val : row) {
            cout << setw(10) << fixed << setprecision(4) << (abs(val) < 1e-9 ? 0.0 : val) << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

// ==========================================
// HÀM: CHỈ CHẠY QUY TRÌNH THUẬN (Đưa về dạng bậc thang)
// ==========================================
void quyTrinhThuan(vector<vector<double>>& A, int m, int n) {
    int i = 0, j = 0;
    int step = 1;
    
    cout << "\n--- LICH SU CHON PHAN TU KHOA (PIVOT) ---\n";
    
    while (i < m && j < n) { 
        // 1. Tìm phần tử trội (Partial Pivoting) để giảm sai số
        int max_row = i;
        for (int t = i + 1; t < m; ++t) {
            if (abs(A[t][j]) > abs(A[max_row][j])) {
                max_row = t;
            }
        }

        // 2. Đổi chỗ hàng nếu phần tử lớn nhất không nằm ở hàng hiện tại
        if (max_row != i) {
            swap(A[i], A[max_row]);
        }

        // 3. Nếu cả cột toàn số 0 (hoặc xấp xỉ 0) thì bỏ qua cột này
        if (abs(A[i][j]) < 1e-9) {
            j++;
            continue;
        }

        // In ra pivot được chọn
        cout << "Lan lap " << step << ": Chon pivot a[" << i+1 << "][" << j+1 << "] = " 
             << fixed << setprecision(4) << A[i][j] << "\n";

        // Nếu đã đến hàng cuối cùng thì dừng
        if (i == m - 1) break; 

        // 4. Khử các phần tử bên dưới thành số 0
        for (int k = i + 1; k < m; ++k) {
            double factor = A[k][j] / A[i][j];
            for (int c = j; c < n; ++c) { 
                A[k][c] -= factor * A[i][c];
            }
        }

        if (step == 1) {
            cout << "\n--- MA TRAN SAU LAN LAP 1 ---\n";
            printMatrix(A);
        }

        i++;
        j++;
        step++; 
    }

    cout << "\n--- MA TRAN BAC THANG (SAU QUY TRINH THUAN) ---\n";
    printMatrix(A);
}

// ==========================================
// HÀM MAIN
// ==========================================
int main() {
    string filename = "test.txt";

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Loi: Khong the mo duoc file!\n";
        return 1;
    }

    vector<vector<double>> A;
    string line;
    while (getline(file, line)) {
        if (line.empty()) continue; 
        stringstream ss(line);
        vector<double> row;
        double val;
        while (ss >> val) {
            row.push_back(val);
        }
        if (!row.empty()) {
            A.push_back(row);
        }
    }
    file.close();

    if (A.empty()) {
        cerr << "Loi: File rong!\n";
        return 1;
    }

    int m = A.size();          // Số hàng
    int n = A[0].size();       // Số cột
    
    cout << "=> Ma tran doc duoc co kich thuoc " << m << " x " << n << ".\n";

    // Chỉ gọi duy nhất 1 hàm
    quyTrinhThuan(A, m, n);

    return 0;
}