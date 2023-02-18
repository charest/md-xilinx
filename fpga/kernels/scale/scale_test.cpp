#include "scale.hpp"

#include <iostream>

constexpr int size = 4096;
constexpr int fact = 2;

int main() {
    
  real_t inout[size], res[size];

  for (int i = 0; i < size; ++i){
    inout[i] = i;
    res[i] = fact*inout[i];
  }

  scale_kernel(inout, fact, size);

  for (int i = 0; i < size; ++i){
    if (res[i] != inout[i])
      return EXIT_FAILURE;
    }

  std::cout << "Test passed.\n";
  return EXIT_SUCCESS;
}
