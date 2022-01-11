#include "Azimint.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"
#include "hlslib/xilinx/DataPack.h"

using hlslib::Stream;
using hlslib::DataPack;

void ReadMemory(DataPack<float, W> const *in, Stream<DataPack<float, W>> stream[D]) {
  for (int j = 0; j < N / W; j++) {
    #pragma HLS PIPELINE II = 1
    auto v = in[j];
    for (int k = 0; k < D; k++) {
      #pragma HLS UNROLL
      stream[k].Push(v);
    }
  }
}

void WriteMemory(Stream<float> stream[D], float *out) {
  for (int i = 0; i < npt; i++) {
    #pragma HLS PIPELINE II=1
    out[i] = stream[i / d].Pop();
  }
}

void sum(Stream<DataPack<float, W>> &radius, Stream<DataPack<float, W>> &data, float rmax, int start, Stream<float> &res) {

  // Calculate radius bounds
  float r1[d], r2[d], sum[d];
  int num[d];
  for (int i = 0; i < d; i++) {
    #pragma HLS PIPELINE II=1
    r1[i] = rmax * (start + i) / npt;
    r2[i] = rmax * (start + i + 1) / npt;
    sum[i] = 0;
    num[i] = 0;
  }

  for (int j = 0; j < N / W; j++) {

    #pragma HLS PIPELINE II=9

    // Receive next input
    auto rad = radius.Pop();
    auto dat = data.Pop();

    // Process input (Depth = 10)
    for (int i = 0; i < d; i++) {
      #pragma HLS UNROLL

      // ==================
      // = AUTO GENERATED =
      // ==================

      bool cond1 = r1[i] <= rad[0] && rad[0] < r2[i];
      float s1 = cond1 ? dat[0] : 0.0f;
      int n1 = cond1 ? 1 : 0;

      bool cond2 = r1[i] <= rad[1] && rad[1] < r2[i];
      float s2 = cond2 ? dat[1] : 0.0f;
      int n2 = cond2 ? 1 : 0;

      bool cond3 = r1[i] <= rad[2] && rad[2] < r2[i];
      float s3 = cond3 ? dat[2] : 0.0f;
      int n3 = cond3 ? 1 : 0;

      bool cond4 = r1[i] <= rad[3] && rad[3] < r2[i];
      float s4 = cond4 ? dat[3] : 0.0f;
      int n4 = cond4 ? 1 : 0;

      float s5 = s1 + s2;
      float s6 = s3 + s4;
      float s7 = s5 + s6;

      int n5 = n1 + n2;
      int n6 = n3 + n4;
      int n7 = n5 + n6;

      sum[i] += s7;
      num[i] += n7;

      // ==================

    }

  }

  // Push out results
  for (int i = 0; i < d; i++) {
    #pragma HLS PIPELINE II=1
    res.Push(num[i] > 0 ? sum[i] / num[i] : 0);
  }

}

void sumAll(float const *radius, float const *data, float rmax, float *res) {
  #pragma HLS DATAFLOW

  HLSLIB_DATAFLOW_INIT();

  Stream<DataPack<float, W>> stream_radius[D];
  Stream<DataPack<float, W>> stream_data[D];
  Stream<float> stream_res[D];

  HLSLIB_DATAFLOW_FUNCTION(ReadMemory, reinterpret_cast<const DataPack<float, W> *>(radius), stream_radius);
  HLSLIB_DATAFLOW_FUNCTION(ReadMemory, reinterpret_cast<const DataPack<float, W> *>(data), stream_data);

  for (int k = 0; k < D; k++) {
    #pragma HLS UNROLL
    HLSLIB_DATAFLOW_FUNCTION(sum, stream_radius[k], stream_data[k], rmax, k * d, stream_res[k]);
  }

  HLSLIB_DATAFLOW_FUNCTION(WriteMemory, stream_res, res);

  HLSLIB_DATAFLOW_FINALIZE();
}

float max(float const *radius) {

  constexpr int w = 32;
  float rmax[w];

  auto radius_v = reinterpret_cast<const DataPack<float, w> *>(radius);

  for (int k = 0; k < w; k++) {
    #pragma HLS UNROLL
    rmax[k] = -std::numeric_limits<float>::infinity();
  }

  for (int i = 0; i < N / w; i++) {
    #pragma HLS PIPELINE II=1
    auto rad = radius_v[i];
    for (int k = 0; k < w; k++) {
      #pragma HLS UNROLL
      rmax[k] = rad[k] > rmax[k] ? rad[k] : rmax[k];
    }
  }

  for (int k = 1; k < w; k++) {
    #pragma HLS UNROLL
    rmax[0] = rmax[0] > rmax[k] ? rmax[0] : rmax[k];
  }

  return rmax[0];

}

void Azimint(float const *data, float const *radius, float *res) {

  #pragma HLS INTERFACE m_axi port=data offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=radius offset=slave bundle=gmem1
  #pragma HLS INTERFACE m_axi port=res offset=slave bundle=gmem2
  #pragma HLS INTERFACE s_axilite port=data bundle=control
  #pragma HLS INTERFACE s_axilite port=radius bundle=control
  #pragma HLS INTERFACE s_axilite port=res bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  sumAll(radius, data, max(radius), res);


}
