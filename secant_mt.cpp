// Phương pháp dây cung
// dựa theo công thức sai số mục tiêu
// | x_k - x* | <= | f(x_k) | / m1

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double mypow(double x, int n) {
  double res = 1;
  while (n--) res *= x;
  return res;
}

double f(double x) { return mypow(x, 5) - 3 * mypow(x, 3) + 2 * x * x - x + 5;}
double df(double x) { return 5*mypow(x, 4) - 9*x*x + 4*x - 1;}
double d2f(double x) { return 20* mypow(x, 3) - 18*x + 4; }

void solve(double a, double b, double epsilon) {
  double d, x_old;

  if (f(a) * d2f(a) > 0) {
    d = a;
    x_old = b;
  } else {
    d = b;
    x_old = a;
  }

  double m1 = min(abs(df(a)), abs(df(b)));
  int k = 1;

  cout << fixed << setprecision(10);
  cout << "Diem Fourier d = " << d << "\n";
  cout << "Diem bat dau x0 = " << x_old << "\n";
  cout << "----------------------------------------------------" << "\n";
  cout << setw(3) << "k"
       << setw(15) << "x_k"
       << setw(15) << "x_k+1"
       << setw(15) << "Sai so\n";

  while (true) {
    double x_new = x_old - (f(x_old) * (x_old - d)) / (f(x_old) - f(d));
    double error = abs(f(x_new)) / m1;

    cout << setw(3) << k
         << setw(15) << x_old
         << setw(15) << x_new
         << setw(15) << error << "\n";

    if (error <= epsilon) {
      cout << "----------------------------------------------------\n";
      cout << "Nghiem xap xi tim duoc: " << x_new << "\n";
      break;
    }

    x_old = x_new;
    k++;

    if (k > 100) {
      cout << "Lap qua 100 lan!\n";
      break;
    }
  }
}

int main() {
  double a, b, epsilon;
  if (!(cin >> a >> b >> epsilon))
    return 0;
  solve(a, b, epsilon);
  return 0;
}
