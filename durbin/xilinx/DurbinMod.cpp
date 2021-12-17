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
void UpdateY(Stream<float> &in, Stream<float> &out, Stream<float> &alpha_in, Stream<float> &reverse_helper, Stream<float> &alpha_out, int k) {
  float alpha = alpha_in.Pop();
  alpha_out.Push(alpha);
  // std::cout<<"UY k " + std::to_string(k) + " " + std::to_string(in.Size()) + "\n";
  for (int i = 0; i < k; i++) {
    #pragma HLS PIPELINE II=1
    out.Push(in.Pop() + alpha * reverse_helper.Pop());
  }
  // Pushing the last element y[k] to the next stream
  out.Push(in.Pop());
  // std::cout<<"Exiting UY k " + std::to_string(k) + " " + std::to_string(out.Size()) + "\n";
}

void ReversePopulate(Stream<float> &in, Stream<float>& out, Stream<float> & reverse, Stream<float> &alpha_in, Stream<float> &alpha_out, int k_, int k)
{ 
  float alpha = alpha_in.Pop();
  alpha_out.Push(alpha);
  // std::cout<<"RP k " + std::to_string(k_) + " " + std::to_string(in.Size()) + "\n";
  if(k == 0)
    return;
  float y[N];
  for(int i=0;i<k_;i++)
  {
    // std::cout<<"RP k " + std::to_string(k_) + " i " + std::to_string(i) + "\n";
    y[i] = in.Pop();
  }
  for(int i=0;i<k_;i++)
  {
    reverse.Push(y[k_-1-i]);
    out.Push(y[i]);
  }
  out.Push(in.Pop());
  // std::cout<<"Exiting RP k " + std::to_string(k_) + " " + std::to_string(out.Size()) + " " + std::to_string(reverse.Size())+ "\n";
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


void ProcessingElement(Stream<float>& r, Stream<float>& y_in, Stream<float>& y_out, Stream<float>& beta_in, Stream<float>& beta_out, Stream<float>& alpha_in, Stream<float>& alpha_out, int k) {
    
    float alpha = alpha_in.Pop();
    float beta = beta_in.Pop();

    // beta *= 1.0 - alpha * alpha
    beta *= 1.0 - alpha * alpha;

    // alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
    float dot = 0;

    for (int i = 0; i < k; i++) {
      // std::cout<<"In PE k " + std::to_string(k) + " i " + std::to_string(i) +  " " +  std::to_string(y_in.Size()) + " " +  std::to_string(r.Size()) +"\n";
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
    // std::cout<<"exiting PE" + std::to_string(k) + "\n";
    alpha_out.Push(alpha);    

}

void InitYAlphaBeta(Stream<float> &y, Stream<float> &y_inp, Stream<float> &alpha, Stream<float> &alpha_inp, Stream<float> &beta, Stream<float> &beta_inp, int size) {
  // std::cout<<"In InitYAB k " + std::to_string(size) + " "+ std::to_string(y_inp.Size()) + "\n";
  
  for(int i=0;i<size;i++)
  {
    y.Push(y_inp.Pop());
  }

  alpha.Push(alpha_inp.Pop());
  beta.Push(beta_inp.Pop());
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

void InitR(Stream<float> r_mod[num_streams], const float r[N], int d)
{
  // r_mod[0].Push(r[0]);
  
  
    for(int i=0;i<num_streams;i++)
    {
      for(int j=0;j<i + d;j++)
      {
        r_mod[i].Push(r[i + d-1-j]);
      }
      r_mod[i].Push(r[i+d]);
    } 
}
void PEDriverFuncFirstIter(const float r[], Stream<float>& y_inp, Stream<float>& y_out,Stream<float>& alpha_inp, Stream<float>& alpha_out, Stream<float>& beta_inp, Stream<float>& beta_out, int d)
{
  int d_offset = D * d;
  int initial_arg_2_init = d_offset == 0? 1: d_offset;
  int final_arg_2_init = d_offset + num_streams;
  
  // #pragma HLS DATAFLOW

  Stream<float, N> y[num_streams+1];
  Stream<float, N> r_mod[num_streams+1];
  Stream<float, N> y_unupdated[num_streams+1];
  Stream<float, N> beta[num_streams+1];
  Stream<float, N> alpha[num_streams+1];
  Stream<float, N> alpha_interim[num_streams+1];
  Stream<float, N> alpha_interim_sync[num_streams+1];
  Stream<float, N> reverse_helper[num_streams+1];
  Stream<float, N> y_reverse_supported[num_streams+1];

  HLSLIB_DATAFLOW_INIT();

  HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y[0] , y_inp, alpha[0], alpha_inp, beta[0], beta_inp, initial_arg_2_init);
  HLSLIB_DATAFLOW_FUNCTION(InitR, r_mod, r, d_offset);
  
  // for k in range(1, r.shape[0])
  int k;
  for(int loopIdx = 0;loopIdx<numItersFirstIter-1;loopIdx++)
  {
    #pragma HLS UNROLL
    k = loopIdx+1;
    // #pragma HLS DEPENDENCE variable=r false
      // std::cout<<"Unrolling k "<< k<<std::endl;
    int newK = d == 0? k: k + d_offset-1;
    HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r_mod[(d == 0)?(k):(k-1)], y[(k-1)], y_unupdated[(k-1)], beta[(k-1)], beta[(k)], alpha[(k-1)], alpha_interim[(k-1)], newK);
    HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated[(k-1)], y_reverse_supported[(k-1)], reverse_helper[(k-1)], alpha_interim[(k-1)], alpha_interim_sync[(k-1)], newK, newK);
    HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_reverse_supported[(k-1)], y[(k)], alpha_interim_sync[(k-1)], reverse_helper[(k-1)], alpha[(k)], newK);
  }
   
  HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y_out, y[k], alpha_out, alpha[k], beta_out, beta[k], final_arg_2_init);
  HLSLIB_DATAFLOW_FINALIZE();
}


void PEDriverFuncMain(const float r[], Stream<float>& y_inp, Stream<float>& y_out,Stream<float>& alpha_inp, Stream<float>& alpha_out, Stream<float>& beta_inp, Stream<float>& beta_out, int d)
{
  int numIters = N/D + 1;
  int d_offset = D * d;
  int initial_arg_2_init = d_offset == 0? 1: d_offset;
  int final_arg_2_init = d_offset + num_streams;
  
  // #pragma HLS DATAFLOW

  Stream<float, N> y[num_streams+1];
  Stream<float, N> r_mod[num_streams+1];
  Stream<float, N> y_unupdated[num_streams+1];
  Stream<float, N> beta[num_streams+1];
  Stream<float, N> alpha[num_streams+1];
  Stream<float, N> alpha_interim[num_streams+1];
  Stream<float, N> alpha_interim_sync[num_streams+1];
  Stream<float, N> reverse_helper[num_streams+1];
  Stream<float, N> y_reverse_supported[num_streams+1];

  HLSLIB_DATAFLOW_INIT();

  HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y[0] , y_inp, alpha[0], alpha_inp, beta[0], beta_inp, initial_arg_2_init);
  HLSLIB_DATAFLOW_FUNCTION(InitR, r_mod, r, d_offset);
  
  // for k in range(1, r.shape[0])
  int k;
  for(int loopIdx = 0;loopIdx<numItersMain-1;loopIdx++)
  {
    #pragma HLS UNROLL

    k = loopIdx + 1;
    // #pragma HLS DEPENDENCE variable=r false
      // std::cout<<"Unrolling k "<< k<<std::endl;
    int newK = d == 0? k: k + d_offset-1;
    HLSLIB_DATAFLOW_FUNCTION(ProcessingElement, r_mod[(d == 0)?(k):(k-1)], y[(k-1)], y_unupdated[(k-1)], beta[(k-1)], beta[(k)], alpha[(k-1)], alpha_interim[(k-1)], newK);
    HLSLIB_DATAFLOW_FUNCTION(ReversePopulate, y_unupdated[(k-1)], y_reverse_supported[(k-1)], reverse_helper[(k-1)], alpha_interim[(k-1)], alpha_interim_sync[(k-1)], newK, newK);
    HLSLIB_DATAFLOW_FUNCTION(UpdateY, y_reverse_supported[(k-1)], y[(k)], alpha_interim_sync[(k-1)], reverse_helper[(k-1)], alpha[(k)], newK);
  }
   
  HLSLIB_DATAFLOW_FUNCTION(InitYAlphaBeta, y_out, y[k], alpha_out, alpha[k], beta_out, beta[k], final_arg_2_init);
  HLSLIB_DATAFLOW_FINALIZE();
}

void DurbinMod(float const *r, float *y_out) {

  #pragma HLS INTERFACE m_axi port=r offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=y_out offset=slave bundle=gmem1
  #pragma HLS INTERFACE s_axilite port=r bundle=control
  #pragma HLS INTERFACE s_axilite port=y_out bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control
  
  Stream<float, N> y[D + 1];
  Stream<float, N> alpha[D + 1];
  Stream<float, N> beta[D + 1];
  float y_init = -r[0];
  y[0].Push(y_init);
  alpha[0].Push(y_init);
  beta[0].Push(1.0);

  PEDriverFuncFirstIter(r, y[0], y[1], alpha[0], alpha[1], beta[0], beta[1], 0);
  for(int i=1;i<D;i++)
  {
    // std::cout<<"Callign PEDriverFunc "<<i<<std::endl;
    PEDriverFuncMain(r, y[i], y[i+1], alpha[i], alpha[i+1], beta[i], beta[i+1], i);
  }
  // std::cout<<"reched final write"<<std::endl;
  WriteMemory(y[(D)], y_out, N);
    
}
