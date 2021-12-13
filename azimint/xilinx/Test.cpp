#include <iostream>
#include <limits>
#include <cmath>
#include <vector>

#include "Azimint.h"

//  def kernel(data, radius, npt):
//    rmax = radius.max()
//    res = np.zeros(npt, dtype=np.float64)
//    for i in range(npt):
//        r1 = rmax * i / npt
//        r2 = rmax * (i + 1) / npt
//        mask_r12 = np.logical_and((r1 <= radius), (radius < r2))
//        values_r12 = data[mask_r12]
//        res[i] = values_r12.mean()
//    return res

void Reference(double const data[], double const radius[], double res[]) {

  //    rmax = radius.max()
  double rmax = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < N; i++) {
    rmax = radius[i] > rmax ? radius[i] : rmax;
  }

  //    for i in range(npt):
  for (int i = 0; i < npt; i++) {

    //        r1 = rmax * i / npt
    //        r2 = rmax * (i + 1) / npt
    double r1 = rmax * i / npt;
    double r2 = rmax * (i + 1) / npt;

    //        mask_r12 = np.logical_and((r1 <= radius), (radius < r2))
    //        values_r12 = data[mask_r12]
    //        res[i] = values_r12.mean()
    double sum = 0;
    int num = 0;
    for (int j = 0; j < N; ++j) {
      if (r1 <= radius[j] && radius[j] < r2) {
        sum += data[j];
        num++;
      }
    }
    res[i] = num > 0 ? sum / num : 0;
  }
}

int main() {
  std::vector<double> data(N);
  std::vector<double> radius(N);
  std::vector<double> res(N);
  std::vector<double> res_ref(N);

  // Initialize inputs
  uint32_t seed = 1;
  for (int i = 0; i < N; i++) {
    data[i] = (seed = seed * 0x510f506fu + 0xf59e9d39u) / (float) 0xFFFFFFFFu;
    radius[i] = (seed = seed * 0x510f506fu + 0xf59e9d39u) / (float) 0xFFFFFFFFu;
  }

  // Run simulation
  Azimint(data.data(), radius.data(), res.data());

  // Reference implementation for comparing the result
  Reference(data.data(), radius.data(), res_ref.data());

  for (int i = 0; i < 10; i++) {
    std::cout << res_ref[i] << "\t";
  }
  std::cout << std::endl;

  // Verify correctness
  for (int i = 0; i < N; i++) {
    const auto diff = std::abs(res_ref[i] - res[i]);
    if (diff >= 1e-3 || std::isnan(res_ref[i]) != std::isnan(res[i])) {
      std::cout << "Mismatch at (" << i << "): "
                << res[i]
                << " (should be " << res_ref[i] << ").\n";
      return 1;
    }
  }
  std::cout << "Test ran successfully.\n";

  return 0;
}
