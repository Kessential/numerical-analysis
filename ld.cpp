#include <bits/stdc++.h>

using namespace std;

double phi(double x) { return 17 / (pow(x, 4) + 25); }

void solve(double x_0, double epsilon, double q) {
  if (q >= 1)
    return;
  double epsilon_0 = (1 - q) * epsilon / q;
  double x_old = x_0, x_new;
  int n = 0;
  cout << fixed << setprecision(8) << setw(3) << "k" << setw(15) << "x_k"
       << setw(15) << "x_k+1" << setw(15) << "Sai so" << "\n";
  while (true) {
    ++n;
    x_new = phi(x_old);
    double error = abs(x_old - x_new);
    cout << setw(3) << n << setw(15) << x_old << setw(15) << x_new << setw(15)
         << error << "\n";
    if (error < epsilon_0) {
      cout << "Nghiem gan dung cua PT la: " << x_new << "\n";
      break;
    }
    x_old = x_new;
    if (n > 1000)
      break;
  }
}
int main() {
  double x_0, epsilon, q;
  cin >> x_0 >> epsilon >> q;
  solve(x_0, epsilon, q);
}
