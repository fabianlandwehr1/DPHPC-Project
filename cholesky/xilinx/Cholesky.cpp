#include "Cholesky.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

using hlslib::Stream;

void ProcessingElement(float A[]) {

  A[0] = A[0]*A[0];

}

void Cholesky(float A[]) {

  #pragma HLS INTERFACE m_axi port=A offset=slave bundle=gmem0
  #pragma HLS INTERFACE s_axilite port=A bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  #pragma HLS DATAFLOW
  ProcessingElement(A);

}
