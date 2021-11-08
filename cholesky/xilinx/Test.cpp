#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "Cholesky.h"

void Reference(float A[]) {
  A[0] = sqrt(A[0]);
  for (int i = 1; i < N; i++) {
    for (int j = 0; j < i; j++) {
      float dot = 0;
      for (int k = 0; k < j; k++) {
        dot += A[i * N + k] * A[i * N + k];
      }
      A[i * N + j] -= dot;
      A[i * N + j] /= A[j * N + j];
    }
    float dot = 0;
    for (int k = 0; k < i; k++) {
      dot += A[i * N + k] * A[i * N + k];
    }
    A[i * N + i] -= dot;
    A[i * N + i] = sqrt(A[i * N + i]);
  }
}

int main() {
  std::vector<float> A(N * N);
  std::vector<float> A_ref(N * N);

  std::default_random_engine rng;
  std::uniform_real_distribution<float> dist;
  std::for_each(A.begin(), A.end(), [&](float &i) { i = dist(rng); });
  A_ref = A;

  // Run simulation
  Cholesky(A.data());

  // Reference implementation for comparing the result
  Reference(A_ref.data());

  // Verify correctness
  for (int i = 0; i < N * N; i++) {
    const auto diff = std::abs(A_ref[i] - A[i]);
    if (diff >= 1e-3) {
      std::cout << "Mismatch at (" << i / N << ", " << i % N << "): "
                << A[i]
                << " (should be " << A_ref[i] << ").\n";
      return 1;
    }
  }
  std::cout << "Test ran successfully.\n";

  return 0;
}
