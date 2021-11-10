#include "Durbin.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

using hlslib::Stream;

void ProcessingElement(float const r[N], float y[N]) {

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
      #pragma HLS PIPELINE II=3
      dot += r[k - i - 1] * y[i];
    }
    alpha = -(r[k] + dot) / beta;

    // y[:k] += alpha * np.flip(y[:k])
    Stream<float> y_in;
    Stream<float> y_out;
    for (int i = 0; i < k; i++) {
      #pragma HLS PIPELINE II=1
      y_in.Push(y[i]);
    }
    for (int i = 0; i < k; i++) {
      #pragma HLS PIPELINE II=1
      y_out.Push(y_in.Pop() + alpha * y[k - i - 1]);
    }
    for (int i = 0; i < k; i++) {
      #pragma HLS PIPELINE II=1
      y[i] = y_out.Pop();
    }

    // y[k] = alpha
    y[k] = alpha;

  }

}

void Durbin(float const *r, float *y) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  ProcessingElement(r, y);

}
