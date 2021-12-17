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
void UpdateY(Stream<float> &in, Stream<float> &out, Stream<float> &alpha_in, Stream<float> &reverse_helper, int k) {
  float alpha = alpha_in.Pop();
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out.Push(in.Pop() + alpha * reverse_helper.Pop());
  }
  // Pushing the last element y[k] to the next stream
  out.Push(in.Pop());
}

void ReversePopulate(Stream<float> &in, Stream<float>& out, Stream<float> & reverse, int k_, int k)
{ 
  if(k == 0)
    return;
  float y[N];
  for(int i=0;i<k_;i++)
  {
    y[i] = in.Pop();
  }
  for(int i=0;i<k_;i++)
  {
    out.Push(y[i]);
    reverse.Push(y[k_-1-i]);
  }
  out.Push(in.Pop());
  // float f = in.Pop();
  // out.Push(f);
  // ReversePopulate(in, out, reverse, N, k-1);
  // reverse.Push(f);

  // if (N == k)
  // {
  //   out.Push(in.Pop());
  // }
}

// Writes from the tail of the pipeline
void WriteMemory(Stream<float> &stream, float *out, int k) {
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out[i] = stream.Pop();
  }
}


void ProcessingElement(Stream<float>& r, Stream<float>& y_in, Stream<float>& y_out, Stream<float>& beta_in, Stream<float>& beta_out, Stream<float>& alpha_in, Stream<float>& alpha_out, Stream<float>& alpha_out_real, int N, int k) {

    float alpha = alpha_in.Pop();
    float beta = beta_in.Pop();

    // beta *= 1.0 - alpha * alpha
    beta *= 1.0 - alpha * alpha;

    // alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
    float dot = 0;

    for (int i = 0; i < k; i++) {
        #pragma HLS PIPELINE II=1
        #pragma HLS DEPENDENCE variable=dot false
        float y = y_in.Pop();
        dot += r.Pop() * y;
        y_out.Push(y);
    }
    alpha = -(r.Pop() + dot) / beta;

    // y[k] = alpha
    y_out.Push(alpha);

    beta_out.Push(beta);
    alpha_out_real.Push(alpha);    
    alpha_out.Push(alpha);    
}

void InitYAlphaBeta(Stream<float> &y, Stream<float> &alpha, Stream<float> &beta, Stream<float>& initVal) {
  float r = initVal.Pop();
  y.Push(-r);
  alpha.Push(-r);
  beta.Push(1.0);
}


// void DebugPrint(Stream<float> &stream_in,  std::string name, int k)
// {
//   std::string print  = name + ":";
//   int n = stream_in.Size();

//   for(int i=0;i<k;i++)
//   {
//     float y = stream_in.Pop();
//     print += std::to_string(y) + " ";
//   }
//   print += "\n";
//   // std::cout<<print;
// }

void InitR(Stream<float> r_mod[N], const float r[N])
{
  r_mod[0].Push(r[0]);
    for(int i=1;i<N;i++)
    {
      for(int j=0;j<i;j++)
      {
        r_mod[i].Push(r[i-1-j]);
      }
      r_mod[i].Push(r[i]);
    } 
}

void DurbinMod(float const *r, float *y_out) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y_out offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y_out bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

    #pragma HLS DATAFLOW

    Stream<float, N> y[N];
    Stream<float, N> r_mod[N];
    Stream<float, N> y_unupdated[N-1];
    Stream<float, N> beta[N];
    Stream<float, N> alpha[N];
    Stream<float, N> alpha_interim[N-1];
    Stream<float, N> reverse_helper[N-1];
    Stream<float, N> y_reverse_supported[N-1];
    
    // Stream<float> y0("y0");
    // Stream<float> y1("y1");
    // Stream<float> y_reverse_supported0("y_reverse_supported0");
    // Stream<float> y_unupdated0("y_unupdated0");
    // Stream<float> y2("y2");
    // Stream<float> y_reverse_supported1("y_reverse_supported1");
    // Stream<float> y_unupdated1("y_unupdated1");
    // Stream<float> beta0("beta0");
    // Stream<float> alpha0("alpha0");
    // Stream<float> reverse_helper0("reverse_helper");
    
    // y[0].Push(-r[0]);
    // alpha[0].Push(-r[0]);
    // beta[0].Push(1.0);

    

    HLSLIB_DATAFLOW_INIT();
    

    HLSLIB_DATAFLOW_FUNCTION(InitR, r_mod, r);
    
    HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y[0], alpha[0], beta[0], r_mod[0]);
    // for k in range(1, r.shape[0])
    for(int k=1;k<N;k++)
    {
      // #pragma HLS DEPENDENCE variable=r false
      // #pragma HLS UNROLL
        // std::cout<<"Unrolling k "<< k<<std::endl;
      HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r_mod[k], y[k-1], y_unupdated[k-1], beta[k-1], beta[k], alpha[k-1], alpha_interim[k-1], alpha[k], N, k);
      HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated[k-1], y_reverse_supported[k-1], reverse_helper[k-1], k, k);
      HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_reverse_supported[k-1], y[k], alpha_interim[k-1], reverse_helper[k-1], k);
    }
      
    HLSLIB_DATAFLOW_FUNCTION(WriteMemory, y[N-1], y_out, N);

    HLSLIB_DATAFLOW_FINALIZE();
}