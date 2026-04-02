// Phương pháp lặp đơn
// theo 2 xấp xỉ liên tiếp

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double phi(double x) { return 17 / (pow(x, 4) + 25); }

void solve(double x_0, double epsilon, double q) {
  if (q >= 1 || q <= 0) {
    cout << "He so co q phai nho hon 1!\n";
    return;
  }

  double epsilon_0 = (1 - q) * epsilon / q;
  double x_old = x_0, x_new;
  int n = 1;

  cout << fixed << setprecision(8);
  cout << setw(3) << "k"
       << setw(15) << "x_k"
       << setw(15) << "x_k+1"
       << setw(15) << "Sai so\n";
  cout << "-------------------------------------------------\n";

  while (true) {
    x_new = phi(x_old);

    double error = abs(x_old - x_new);

    cout << setw(3) << n 
         << setw(15) << x_old 
         << setw(15) << x_new 
         << setw(15) << error << "\n";

    if (error < epsilon_0) {
      cout << "-------------------------------------------------\n";
      cout << "Nghiem gan dung cua PT la: " << x_new << "\n";
      break;
    }

    x_old = x_new;
    ++n;

    if (n > 1000) {
      cout << "Thuat toan ko hoi tu sau 1000 lan lap\n";
      break;
    }
  }
}

int main() {
  double x_0, epsilon, q;
  cin >> x_0 >> epsilon >> q;
  solve(x_0, epsilon, q);
}
