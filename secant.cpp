#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double f(double x) { return exp(x) - cos(2 * x); }
double df(double x) { return exp(x) + 2 * sin(2 * x); }
double d2f(double x) { return exp (x) + 4 * cos(2 * x); }

void solve(double a, double b, double epsilon) {
  double d, x_old;

  // 1. Chọn điểm Fourier d và điểm bắt đầu x0 [cite: 16, 17, 18]
  if (f(a) * d2f(a) > 0) {
    d = a;
    x_old = b;
  } else {
    d = b;
    x_old = a;
  }

  // 2. Tính m1 = min |f'(x)| trên [a, b]
  double m1 = min(abs(df(a)), abs(df(b)));

  cout << fixed << setprecision(10);
  cout << "Diem Fourier d = " << d << "\n";
  cout << "Diem bat dau x0 = " << x_old << "\n";
  cout << "----------------------------------------------------" << "\n";
  cout << setw(3) << "k" << setw(15) << "x_k" << setw(15) << "x_k+1" << setw(15)
       << "Sai so"
       << "\n";

  int k = 1;
  while (true) {
    // 3. Công thức lặp dây cung [cite: 47, 106]
    double x_new = x_old - (f(x_old) * (x_old - d)) / (f(x_old) - f(d));

    // 4. Tính sai số mục tiêu
    double error = abs(f(x_new)) / m1;

    cout << setw(3) << k << setw(15) << x_old << setw(15) << x_new << setw(15)
         << error << "\n";

    // 5. Kiểm tra điều kiện dừng [cite: 59]
    if (error <= epsilon) {
      cout << "----------------------------------------------------\n";
      cout << "Nghiem xap xi tim duoc: " << x_new << "\n";
      break;
    }
    x_old = x_new;
    k++;
    if (k > 100)
      break;
  }
}

int main() {
  double a, b, epsilon;
  if (!(cin >> a >> b >> epsilon))
    return 0;
  solve(a, b, epsilon);
  return 0;
}
