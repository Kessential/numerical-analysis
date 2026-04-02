// Phương pháp lặp đơn
// theo xấp xỉ ban đầu

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double phi(double x) { return 17 / (pow(x, 4) + 25); }

void solve(double x_0, double epsilon, double q) {
  if (q >= 1 || q <= 0) {
    cout << "He so co q phai nho hon 1!\n";
    return;
  }

  double x_1 = phi(x_0);
  double d = abs(x_1 - x_0);

  double tu = log((epsilon * (1 - q)) / d);
  double mau = log(q);

  int n_target = ceil(tu / mau);

  cout << fixed << setprecision(8);
  cout << "|x1 - x0| = " << d << endl;
  cout << "So buoc lap du kien de dat sai so: " << n_target << endl;
  cout << "-------------------------------------------------\n";
  cout << setw(3) << "k" << setw(15) << "x_k" << "\n";

  // Buoc 3: Chay vong lap dung n_target lan
  double x_current = x_0;
  for (int k = 1; k <= n_target; ++k) {
    x_current = phi(x_current);
    cout << setw(3) << k << setw(15) << x_current << "\n";
  }

  cout << "-------------------------------------------------\n";
  cout << "Nghiem tim duoc sau " << n_target << " buoc: " << x_current << "\n";
}

int main() {
  double x_0, epsilon, q;
  cin >> x_0 >> epsilon >> q;
  solve(x_0, epsilon, q);
}
