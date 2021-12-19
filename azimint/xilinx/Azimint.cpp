#include "Azimint.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

using hlslib::Stream;

void ProcessingElement(float const *data, float const *radius, float *res) {

  float rmax = -std::numeric_limits<float>::infinity();
  for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE II=1
    rmax = radius[i] > rmax ? radius[i] : rmax;
  }

  for (int i = 0; i < npt; i++) {

    float r1 = rmax * i / npt;
    float r2 = rmax * (i + 1) / npt;

    float sum = 0;
    int num = 0;
    for (int j = 0; j < N; ++j) {
      #pragma HLS PIPELINE II=3
      if (r1 <= radius[j] && radius[j] < r2) {
        sum += data[j];
        num++;
      }
    }
    res[i] = num > 0 ? sum / num : 0;
  }

}

void Azimint(float const *data, float const *radius, float *res) {

  #pragma HLS INTERFACE m_axi port=data offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=radius offset=slave bundle=gmem1
  #pragma HLS INTERFACE m_axi port=res offset=slave bundle=gmem2
  #pragma HLS INTERFACE s_axilite port=data bundle=control
  #pragma HLS INTERFACE s_axilite port=radius bundle=control
  #pragma HLS INTERFACE s_axilite port=res bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  ProcessingElement(data, radius, res);

}
