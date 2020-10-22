#include <cassert>
#include <iostream>

extern "C" {
double double_me_(double *);

// Wrongly prototyped
double double_me_f_(double *);

// Left to the compiler
double double_me_auto_(double *);

void logical_test_(int *);
}

int main() {
  double i = 5E100;
  std::cout << "[CFORTRAN]: Check i = " << i << " *2 = " << double_me_(&i)
            << "(FORTRAN FUNCTION DOUBLE PRECISION)" << std::endl;
  assert(i * 2 == double_me_(&i));

  std::cout << "[CFORTRAN]: Check i = " << i << " *2 = " << double_me_f_(&i)
            << " (Should be wrong due to incorrect prototype)" << std::endl;
  assert(i * 2 != double_me_f_(&i));

  // std::cout << "[CFORTRAN]: Check i = " << i << " *2 = " << double_me_f_(&i)
  //           << " (Hopefully correct, assuming default REAL type is double)"
  //           << std::endl;
  // assert(i * 2 == double_me_auto_(&i));

  std::cout << "Should be true" << std::endl;
  int a = true, b = false;
  logical_test_(&a);
  std::cout << "Should be false" << std::endl;
  logical_test_(&b);
}
