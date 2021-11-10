#include "Durbin.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

using hlslib::Stream;
void reverseArray(float y[N], float revY[N], int k)
{
  for(int i=0;i<k;i++)
  {
    #pragma HLS UNROLL
    revY[i] = y[k-1-i];
  }
}

void ProcessingElement(float const r[N], float y[N]) {
  
  float beta = 1.0;
  float alpha = -r[0];
  y[0] = -r[0];
  float reverseArr[N];

  for (int i=1;i<N;i++)
  {

    beta *= (1.0 - alpha * alpha);
    float dotProdAccumulator = 0.0;

    for(int j=0;j<i;j++)
    {
      #pragma HLS PIPELINE II=3
      dotProdAccumulator += y[j] * r[i-1 - j];
    } 
    
    alpha = -1 * (r[i] + dotProdAccumulator) / beta;

    reverseArray(y, reverseArr, i);
    for(int j=0;j<i;j++)
    {
      #pragma HLS UNROLL
      y[j] += alpha * reverseArr[j];
    } 

    y[i] = alpha;
  }

}

void Durbin(float const *r, float *y) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  #pragma HLS DATAFLOW
  ProcessingElement(r, y);

}
