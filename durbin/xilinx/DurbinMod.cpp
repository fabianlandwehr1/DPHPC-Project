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
    reverse.Push(y[k_-1-i]);
    out.Push(y[i]);
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
        y_out.Push(y);
        dot += r.Pop() * y;
    }
    alpha = -(r.Pop() + dot) / beta;

    // y[k] = alpha
    y_out.Push(alpha);

    beta_out.Push(beta);
    alpha_out_real.Push(alpha);    
    alpha_out.Push(alpha);    
}

void InitYAlphaBeta(Stream<float> &y, Stream<float> &alpha, Stream<float> &beta, float initVal) {
  float r = initVal;
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
  // r_mod[0].Push(r[0]);
    for(int i=1;i<N;i++)
    {
      for(int j=0;j<i;j++)
      {
        r_mod[i%num_streams].Push(r[i-1-j]);
      }
      r_mod[i%num_streams].Push(r[i]);
    } 
}

void PEDriverFunc(Stream<Vec_t> &r_inp, Stream)

void DurbinMod(float const *r, float *y_out) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y_out offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y_out bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

    #pragma HLS DATAFLOW

    Stream<float, N> y[num_streams];
    Stream<float, N> r_mod[num_streams];
    Stream<float, N> y_unupdated[num_streams];
    Stream<float, N> beta[num_streams];
    Stream<float, N> alpha[num_streams];
    Stream<float, N> alpha_interim[num_streams];
    Stream<float, N> reverse_helper[num_streams];
    Stream<float, N> y_reverse_supported[num_streams];
    
    Stream<float> y0("y0");
    Stream<float> y1("y1");
    // Stream<float> y_reverse_supported0("y_reverse_supported0");
    // Stream<float> y_unupdated0("y_unupdated0");
    // Stream<float> y2("y2");
    // Stream<float> y_reverse_supported1("y_reverse_supported1");
    // Stream<float> y_unupdated1("y_unupdated1");
    Stream<float> beta0("beta0");
    Stream<float> alpha0("alpha0");
    Stream<float> reverse_helper0("reverse_helper");
    
    // y[0].Push(-r[0]);
    // alpha[0].Push(-r[0]);
    // beta[0].Push(1.0);

    

    InitYAlphaBeta(y[0], alpha[0], beta[0], r[0]);
    HLSLIB_DATAFLOW_INIT();
    

    HLSLIB_DATAFLOW_FUNCTION(InitR, r_mod, r);
    
    // for k in range(1, r.shape[0])
    for(int k=1;k<N;k++)
    {
      // #pragma HLS DEPENDENCE variable=r false
      // #pragma HLS UNROLL
        // std::cout<<"Unrolling k "<< k<<std::endl;
      HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r_mod[(k)%num_streams], y[(k-1)%num_streams], y_unupdated[(k-1)%num_streams], beta[(k-1)%num_streams], beta[(k)%num_streams], alpha[(k-1)%num_streams], alpha_interim[(k-1)%num_streams], alpha[(k)%num_streams], N, k);
      HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated[(k-1)%num_streams], y_reverse_supported[(k-1)%num_streams], reverse_helper[(k-1)%num_streams], k, k);
      HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_reverse_supported[(k-1)%num_streams], y[(k)%num_streams], alpha_interim[(k-1)%num_streams], reverse_helper[(k-1)%num_streams], k);
    }
      

    HLSLIB_DATAFLOW_FINALIZE();
    WriteMemory(y[(N-1)%num_streams], y_out, N);
}
