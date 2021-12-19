
void conv2d(float const *input, float const *weights, float const *bias, float *output);

// Top-level entry function, not relevant for this example
void Example1(float const *input, float const *weights, float const *bias, float *output) {
  #pragma HLS INTERFACE m_axi port=input bundle=gmem0 offset=slave
  #pragma HLS INTERFACE m_axi port=weights bundle=gmem1 offset=slave
  #pragma HLS INTERFACE m_axi port=bias bundle=gmem2 offset=slave
  #pragma HLS INTERFACE m_axi port=ouput bundle=gmem3 offset=slave
  #pragma HLS INTERFACE s_axilite port=input
  #pragma HLS INTERFACE s_axilite port=weights
  #pragma HLS INTERFACE s_axilite port=bias
  #pragma HLS INTERFACE s_axilite port=ouput
  #pragma HLS INTERFACE s_axilite port=return

  conv2d(input,weights,bias,output)
}
