// Phương pháp tiếp tuyến (Newton - Raphson) 
// theo công thức sai số mục tiêu
// | x_k - x* | <= | f(x_k) | / m1

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double f(double x) { return exp(x) - 4 * cos(x / 2) + 1; }
double df(double x) { return exp(x) + 2 * sin(x / 2); }
double d2f(double x) { return exp(x) + cos(x / 2); }

void solve(double a, double b, double epsilon) {
  double x_old;

  if (f(a) * d2f(a) > 0) {
    x_old = a;
  } else {
    x_old = b;
  }

  double m1 = min(abs(df(a)), abs(df(b)));

  cout << fixed << setprecision(10);
  cout << "Diem bat dau x0 = " << x_old << "\n";
  cout << "----------------------------------------------------" << "\n";
  cout << setw(3) << "k" 
       << setw(15) << "x_k"
       << setw(15) << "x_k+1"
       << setw(15) << "Sai so\n";

  int k = 1;

  while (true) {
    double x_new = x_old - f(x_old) / df(x_old);
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
      cout << "Thuat toan chay hon 100 lan lap!\n";
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
