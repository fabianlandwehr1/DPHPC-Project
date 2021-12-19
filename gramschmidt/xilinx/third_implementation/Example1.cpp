#include <stdio.h>
#include <math.h> 
#include "hls_math.h"
constexpr int N = 64;

void gramschmidt_naive(float *A, float *Q, float *R) {
  float p[N];
  float s[N];
  for(int i=0; i<N; i++){
    float val = A[i*N];
    for(int j=0; j<N; j++){
      #pragma HLS PIPELINE II = 1
      p[j] += val*A[i*N+j];
      #pragma HLS DEPENDENCE variable=p false
    }
  }
  float squared_norm = p[0];
  float r_norm = hls::rsqrt(squared_norm);
  for(int i=1; i<N; i++){
    #pragma HLS PIPELINE II = 1
    s[i] = A[i]/squared_norm;
  }
  for(int i=0; i<N; i++){
    #pragma HLS PIPELINE II = 1
    R[i*N] = p[i]*r_norm;
  }
  for(int i=0; i<N-1; i++){
    float q[N];
    for(int j=0; j<N; j++){
      #pragma HLS PIPELINE II = 1
      q[j] = A[i*N+j]*r_norm;
    }
    for(int j=0; j<N; j++){
      #pragma HLS PIPELINE II = 1
      Q[i*N+j] = q[j];
    }
    for(int k=0; k<N; k++){
      float val = A[k*N+i];
      for(int j=i+1; j<N; j++){
        #pragma HLS PIPELINE II = 1
        A[k*N+j] -= s[j]*val;
        #pragma HLS DEPENDENCE variable=A false
      }
    }
    for(int k=0; k<N; k++){
      float val = A[k*N+i+1];
      for(int j=i+1; j<N; j++){
        #pragma HLS PIPELINE II = 1
        p[j] += A[k*N+j]*val; //p[j] = <a_i+1. a_j>
        #pragma HLS DEPENDENCE variable=p false
      }
    }
    squared_norm = p[i+1];
    r_norm = hls:: rsqrt(squared_norm);
    for(int j=i+2; j<N; j++){
      #pragma HLS PIPELINE II = 1
      s[j] = p[j]/squared_norm;
    }
    for(int j=i+1; j<N; j++){
      #pragma HLS PIPELINE II = 1
      R[(i+1)*N+j] = p[j]*r_norm;
    }
  }
  for(int i=0; i<N; i++){
    #pragma HLS PIPELINE II = 2
    Q[(N-1)*N+i] = A[(N-1)*N+i]*r_norm;
  }
  
}

