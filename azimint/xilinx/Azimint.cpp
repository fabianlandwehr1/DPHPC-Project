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

  float r1[npt], r2[npt], sum[npt], num[npt];
  for (int i = 0; i < npt; i++) {
    #pragma HLS PIPELINE II=1
    r1[i] = rmax * i / npt;
    r2[i] = rmax * (i + 1) / npt;
    sum[i] = 0;
    num[i] = 0;
  }

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < npt; ++i) {
      #pragma HLS PIPELINE II=1
      if (r1[i] <= radius[j] && radius[j] < r2[i]) {
        sum[i] += data[j];
        num[i]++;
      }
    }
  }

  for (int i = 0; i < npt; i++) {
#pragma HLS PIPELINE II=1
    res[i] = num[i] > 0 ? sum[i] / num[i] : 0;
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
