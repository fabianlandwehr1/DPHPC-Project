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
void UpdateY(Stream<float> &in, Stream<float> &out, Stream<float> &alpha_in, Stream<float> &alpha_out, Stream<float> &reverse_helper, int k) {
  float alpha = alpha_in.Pop();
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out.Push(in.Pop() + alpha * reverse_helper.Pop());
  }
  // Pushing the last element y[k] to the next stream
  out.Push(in.Pop());
  alpha_out.Push(alpha);
}

void ReversePopulate(Stream<float> &in, Stream<float>& out, Stream<float> & reverse, int k)
{ 
  if(k == 0)
    return;
  
  float f = in.Pop();
  out.Push(f);
  ReversePopulate(in, out, reverse, k-1);
  reverse.Push(f);
}

// Writes from the tail of the pipeline
void WriteMemory(Stream<float> &stream, float *out, int k) {
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out[i] = stream.Pop();
  }
}


void ProcessingElement(float const r[N], Stream<float>& y_in, Stream<float>& y_out, Stream<float>& beta_in, Stream<float>& beta_out, Stream<float>& alpha_in, Stream<float>& alpha_out, Stream<float>& reverse_helper, int N, int k) {

    float alpha = alpha_in.Pop();
    float beta = beta_in.Pop();

    // beta *= 1.0 - alpha * alpha
    beta *= 1.0 - alpha * alpha;

    // alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
    float dot = 0;
    std::cout<<"Reached this stage"<<std::endl;
    for (int i = 0; i < k; i++) {
        #pragma HLS PIPELINE II=1
        // #pragma HLS DEPENDENCE variable=dot false
        float y = y_in.Pop();
        dot += r[k - i - 1] * y;
        y_out.Push(y);
    }
    std::cout<<"Completed for loop"<<std::endl;
    alpha = -(r[k] + dot) / beta;

    alpha_out.Push(alpha);    
    // y[k] = alpha
    y_out.Push(alpha);

    beta_out.Push(beta);
}

void InitYAlphaBeta(Stream<float> &y, Stream<float> &alpha, Stream<float> &beta, float initVal) {
  y.Push(initVal);
  alpha.Push(initVal);
  beta.Push(1.0);
}

void MoveToEnd(Stream<float> &stream)
{
  stream.Push(stream.Pop());
}
void DurbinMod(float const *r, float *y_out) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y_out offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y_out bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

    #pragma HLS DATAFLOW

    // Stream<float> y[N];
    // Stream<float> y_unupdated[N-1];
    // Stream<float> beta;
    // Stream<float> alpha;
    // Stream<float> reverse_helper("reverse_helper");
    
    Stream<float> y0("y0");
    Stream<float> y1("y1");
    Stream<float> y2("y2");
    Stream<float> y_unupdated0("y_unupdated0");
    Stream<float> beta0("beta0");
    Stream<float> alpha0("alpha0");
    Stream<float> reverse_helper("reverse_helper");
    
    // y[0].Push(-r[0]);
    // alpha[0].Push(-r[0]);
    // beta[0].Push(1.0);

    HLSLIB_DATAFLOW_INIT();
    
    HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y0, alpha0, beta0, -r[0]);
    // for k in range(1, r.shape[0])
    // for(int k=1;k<N;k++)
    // {
    //     #pragma HLS UNROLL
    //     // std::cout<<"Unrolling k "<< k<<std::endl;
    //   HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r, y[k-1], y_unupdated[k-1], beta[k-1], beta[k], alpha[k-1], alpha_intermediate[k-1], reverse_helper, N, k);
    //   HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated[k-1], y_unupdated[k-1], reverse_helper, k);
    //   HLSLIB_DATAFLOW_FUNCTION(MoveToEnd, y_unupdated[k-1]);
    //   HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_unupdated[k-1], y[k], alpha_intermediate[k-1], alpha[k], reverse_helper, k);
    // }

      std::cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<std::endl;

      HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r, y0, y_unupdated0, beta0, beta0, alpha0, alpha0, reverse_helper, N, 1);
      HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated0, y_unupdated0, reverse_helper, 1);
      HLSLIB_DATAFLOW_FUNCTION(MoveToEnd, y_unupdated0);
      HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_unupdated0, y1, alpha0, alpha0, reverse_helper, 1);
    
      HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r, y1, y_unupdated0, beta0, beta0, alpha0, alpha0, reverse_helper, N, 2);
      HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated0, y_unupdated0, reverse_helper, 2);
      HLSLIB_DATAFLOW_FUNCTION(MoveToEnd, y_unupdated0);
      HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_unupdated0, y2, alpha0, alpha0, reverse_helper, 2);
    

    HLSLIB_DATAFLOW_FUNCTION(WriteMemory, y2, y_out, N);

    HLSLIB_DATAFLOW_FINALIZE();
}
