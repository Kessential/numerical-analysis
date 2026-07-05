#include <iostream>
#include <iomanip>
#include <cmath>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen; // Sử dụng namespace của Eigen để code ngắn gọn hơn

const double PI = acos(-1);
const double EXP = exp(1);
const double TERM = (10 * PI - 3) / 3.0;

// 1. Khai báo hệ hàm F(X) trả về một Vector2d
Vector3d F(const Vector3d& X) {
    Vector3d f;
    double x1 = X(0), x2 = X(1), x3 = X(2);
    f(0) = 10 * x1 - 2 * x2 * x2 + x2 - 2 * x3 - 5;
    f(1) = 8 * x2 * x2 + 4 * x3 * x3 - 9;
    f(2) = 8 * x2 * x3 + 4;
    return f;
}

// 2. Khai báo ma trận Jacobi J(X) trả về một Matrix2d
Matrix3d J(const Vector3d& X) {
    Matrix3d jacobi;
    double x1 = X(0), x2 = X(1), x3 = X(2);
    jacobi(0, 0) = 10.0;
    jacobi(0, 1) = -4 * x2 + 1;
    jacobi(0, 2) = -2.0;

    jacobi(1, 0) = 0;
    jacobi(1, 1) = 16 * x2;
    jacobi(1, 2) = 8 * x3;

    jacobi(2, 0) = 0;
    jacobi(2, 1) = 8 * x3;
    jacobi(2, 2) = 8 * x2;
    return jacobi;
}

int main() {
    // Xấp xỉ ban đầu X_0 = (0, 0)^T
    Vector3d X(0, 0, 0); 
    int n_iterations = 5; // Tính đến X_5
    
    cout << fixed << setprecision(6);
    cout << setw(3) << "k"
         << setw(15) << "x_1,k"
         << setw(15) << "x_2,k"
         << setw(15) << "x_3,k\n";

    cout << string(50, '-') + "\n";

    cout << setw(3) << 0
         << setw(15) << X(0)
         << setw(15) << X(1)
         << setw(15) << X(2) << "\n";

        // 3. Vòng lặp Newton
    for (int k = 1; k <= n_iterations; ++k) {
        
        // Giải hệ J(X) * Delta_X = -F(X) 
        // Với ma trận 2x2, dùng .inverse() là nhanh và đủ chính xác.
        // Hoặc có thể dùng solver chuẩn: Vector2d deltaX = J(X).colPivHouseholderQr().solve(-F(X));
        Vector3d deltaX = J(X).inverse() * (-F(X)); 
        
        // Cập nhật nghiệm X_{k+1} = X_k + Delta_X
        X = X + deltaX;

        cout << setw(3) << k
             << setw(15) << X(0)
             << setw(15) << X(1) 
             << setw(15) << X(2) << "\n";
    }

    cout << "-----------------------------------\n";
    cout << "Xap xi nghiem X_5: \n" << X << "\n";

    return 0;
}