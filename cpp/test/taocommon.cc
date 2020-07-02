#include <aocommon/matrix2x2.h>
#include <iostream>

int main() {
  double unit[4] = {1.0, 0.0, 0.0, 1.0};
  double e1, e2;
  Matrix2x2::EigenValues(unit, e1, e2);
  std::cout << "EigenValue 1 " << e1 << "\nEigenvalue 2 " << e2 << std::endl;
  return 0;
}
