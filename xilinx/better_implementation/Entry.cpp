void gramschmidt_naive(float *A, float *Q, float *R);

// Top-level entry function, not relevant for this example
void Example1(float *A, float *Q, float *R) {
  #pragma HLS INTERFACE m_axi port=A bundle=gmem0 offset=slave
  #pragma HLS INTERFACE m_axi port=Q bundle=gmem1 offset=slave
  #pragma HLS INTERFACE m_axi port=R bundle=gmem2 offset=slave
  #pragma HLS INTERFACE s_axilite port=A
  #pragma HLS INTERFACE s_axilite port=Q
  #pragma HLS INTERFACE s_axilite port=R
  #pragma HLS INTERFACE s_axilite port=return
  gramschmidt_naive(A,Q,R);
}
