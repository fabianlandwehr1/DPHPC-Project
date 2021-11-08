#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "Durbin.h"

//  def kernel(r):
//
//      y = np.empty_like(r)
//      alpha = -r[0]
//      beta = 1.0
//      y[0] = -r[0]
//
//      for k in range(1, r.shape[0]):
//          beta *= 1.0 - alpha * alpha
//          alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
//          y[:k] += alpha * np.flip(y[:k])
//          y[k] = alpha
//
//      return y

void Reference(float const r[], float y[]) {

  float alpha = -r[0];
  float beta = 1.0;
  y[0] = -r[0];

  // for k in range(1, r.shape[0]):
  for (int k = 1; k < N; k++) {

    // beta *= 1.0 - alpha * alpha
    beta *= 1.0 - alpha * alpha;

    // alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
    float dot = 0;
    for (int i = 0; i < k; ++i) {
      dot += r[k - i - 1] * y[i];
    }
    alpha = -(r[k] + dot) / beta;

    // y[:k] += alpha * np.flip(y[:k])
    float y_op[k] = {0};
    for (int i = 0; i < k; ++i) {
      y_op[i] = alpha * y[k - i - 1];
    }
    for (int i = 0; i < k; ++i) {
      y[i] += y_op[i];
    }

    // y[k] = alpha
    y[k] = alpha;

  }

}

int main() {
  std::vector<float> r(N);
  std::vector<float> y(N);
  std::vector<float> y_ref(N);

  for (int i = 0; i < N; ++i) {
    r[i] = N + 1 - i;
  }

  // Run simulation
  Durbin(r.data(), y.data());

  // Reference implementation for comparing the result
  Reference(r.data(), y_ref.data());

  // Verify correctness
  for (int i = 0; i < N; i++) {
    const auto diff = std::abs(y_ref[i] - y[i]);
    if (diff >= 1e-3) {
      std::cout << "Mismatch at (" << i << "): "
                << y[i]
                << " (should be " << y_ref[i] << ").\n";
      return 1;
    }
  }
  std::cout << "Test ran successfully.\n";

  return 0;
}
