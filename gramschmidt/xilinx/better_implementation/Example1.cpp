#include <stdio.h>
#include <math.h> 
constexpr int N = 64;

void gramschmidt_naive(float *A, float *Q, float *R) {

  for(int i=0; i<N; i++){
    float squared_norms[N];
    for(int j=i; j<N; j++){
      float squared_norm = 0.0;
      for(int k=0; k<N; k++){
        float val1 = A[k*N+i];
        float val2 = A[k*N+j];
        squared_norm = squared_norm + val1*val2;
      }
      squared_norms[j] = squared_norm; //norm = <a_i, a_j>
    }
    float norm = squared_norms[i]; //|a_i|
    for(int k=0; k<N; k++){
      Q[k*N+i] = A[k*N+i]/norm;
    } 
    for(int j=i+1; j<N; j++){
      float factor = squared_norms[j]/squared_norms[i];
      for(int k=0; k<N; k++){
        A[k*N+j] -= factor*A[k*N+i];
      }
      R[i*N+j] = factor*norm;
    }
  }
}

