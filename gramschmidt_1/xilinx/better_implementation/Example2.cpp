#include <stdio.h>
#include <math.h> 
#include "hls_math.h"
constexpr int N = 64;

void gramschmidt_naive(float *A, float *Q, float *R) {
  for(int i=0; i<N; i++){
    float squared_norms[N];
    float buffer[N];
    for(int j=0; j<N; j++){
      #pragma HLS PIPELINE II = 1
      buffer[j] = A[j*N+i];
    } 
    for(int j=0; j<N; j++){
      float val = buffer[j];
      for(int k=i; k<N; k++){
        #pragma HLS PIPELINE II = 1
        float val2 = A[j*N+k];
        squared_norms[k] += val*val2; //row-wise access instead column-wise, so we have no data dependencies in the inner loop
        #pragma HLS DEPENDENCE variable=squared_norms false
      }// = <a_i, a_j>
    }
    float norm = hls::rsqrt(squared_norms[i]); //1/|a_i|
    for(int k=0; k<N; k++){
      #pragma HLS PIPELINE II = 1
      Q[k*N+i] = buffer[k]*norm;
    } 
    float squared_norm = squared_norms[i];

    for(int k=0; k<N; k++){
      float factor = buffer[k]/squared_norm;
      for(int j=i+1; j<N; j++){
        #pragma HLS PIPELINE II = 1
        A[k*N+j] -= squared_norms[j]*factor;
        #pragma HLS DEPENDENCE variable=A false
      }
    }

    for(int j=i; j<N; j++){
      #pragma HLS PIPELINE II = 1
      R[i*N+j] = squared_norms[j]*norm;
    }

  }
}

