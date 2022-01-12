#include <stdio.h>
#include <math.h> 
#include "hls_math.h"
constexpr int C_in = 3;
constexpr int C_out = 16;
constexpr int H = 32;
constexpr int K = 2;
constexpr int W = 32;
constexpr int N = 8; 
//input of size N*H*W*C_in, weights of size K*K*C_in*C_out, bias of size C_out, output of size N*H_out*W_out*C_out
//output[i,:,:,j] = conv(input[i,:,:,:]*weights[:,:,:,j])+bias[j]
//output[i,j,p,q] = sum(input[i,j:j+k,p:p+k,:]*weights[:,:,:,q])+bias[q]
void ProcessingElement(float const *input, float const *weights, float const *bias, float *output){
  int H_out = H-K+1;
  int W_out = W-K+1;
  int i_o1 = H_out*W_out*C_out;  //offsets for indexing. i: input; w: weights; o: output 
  int i_o2 = W_out*C_out;
  int i_o3 = C_out;
  int w_o1 = K*C_in*C_out;
  int w_o2 = C_in*C_out;
  int w_o3 = C_out;
  int o_o1 = H_out*W_out*C_out;
  int o_o2 = W_out*C_out;
  int o_o3 = C_out;
  for(int i=0; i<N; i++){
    for(int j=0; j<H_out; j++){
      for(int p=0; p<W_out; p++){
        for(int q=0; q<C_out; q++){
          float val = bias[q]
          for(int k1=0; k1<K; k1++){
            for(int k2=0; k2<K; k2++){
              for(int c=0; c<C_in; c++){
                val += input[i*i_o1+(j+k1)*i_o2+(p+k2)*i_o3+c]*weights[k1*w_o1+k2*w_o2+c*w_o3+q];
              }
            }
          }
          output[i*o_o1+j*o_o2+p*o_o3+q] = val;
        }
      }
    }
  }
  
}

void conv2d1(float const *input, float const *weights, float const *bias, float *output) {
  #pragma HLS INTERFACE m_axi port=input bundle=gmem0 offset=slave
  #pragma HLS INTERFACE m_axi port=weights bundle=gmem1 offset=slave
  #pragma HLS INTERFACE m_axi port=bias bundle=gmem2 offset=slave
  #pragma HLS INTERFACE m_axi port=ouput bundle=gmem3 offset=slave
  #pragma HLS INTERFACE s_axilite port=input
  #pragma HLS INTERFACE s_axilite port=weights
  #pragma HLS INTERFACE s_axilite port=bias
  #pragma HLS INTERFACE s_axilite port=ouput
  #pragma HLS INTERFACE s_axilite port=return

  ProcessingElement(input,weights,bias,output)
}
