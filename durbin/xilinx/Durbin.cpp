#include "Durbin.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

using hlslib::Stream;

// Reads into the head of the pipeline
void ReadMemory(float const *in, Stream<float> &stream, int k) {
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    stream.Push(in[i]);
  }
}

// Update each element of y
void UpdateY(Stream<float> &in, Stream<float> &out, float *y, float alpha, int k) {
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out.Push(in.Pop() + alpha * y[k - i - 1]);
  }
  out.Push(alpha);
}

// Writes from the tail of the pipeline
void WriteMemory(Stream<float> &stream, float *out, int k) {
  for (int i = 0; i <= k; i++) {
    #pragma HLS PIPELINE II=1
    out[i] = stream.Pop();
  }
}

void ProcessingElement(float const r[N], float y[N]) {

  float buf[N];
  float *y_tmp = buf;

  float alpha = -r[0];
  float beta = 1.0;
  y[0] = -r[0];

  // for k in range(1, r.shape[0]):
  for (int k = 1; k < N; k++) {

    // beta *= 1.0 - alpha * alpha
    beta *= 1.0 - alpha * alpha;

    // alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
    float dot = 0;
    for (int i = 0; i < k; i++) {
      #pragma HLS PIPELINE II=3
      dot += r[k - i - 1] * y[i];
    }
    alpha = -(r[k] + dot) / beta;

    // y[:k] += alpha * np.flip(y[:k])
    // y[k] = alpha
    HLSLIB_DATAFLOW_INIT();
    Stream<float> y_in;
    Stream<float> y_out;
    HLSLIB_DATAFLOW_FUNCTION(ReadMemory, y, y_in, k);
    HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_in, y_out, y, alpha, k);
    HLSLIB_DATAFLOW_FUNCTION(WriteMemory, y_out, y_tmp, k);
    HLSLIB_DATAFLOW_FINALIZE();

    float *tmp = y;
    y = y_tmp;
    y_tmp = tmp;

  }

  if (N % 2 == 0) {

    // Switch y arrays again
    float *tmp = y;
    y = y_tmp;
    y_tmp = tmp;

    // Copy y_tmp to y
    for (int i = 0; i < N; i++) {
      y[i] = y_tmp[i];
    }
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
