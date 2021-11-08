#include <stdio.h>
#include <math.h> 
constexpr int N = 64;

void gramschmidt_naive(float *A, float *Q, float *R) {

  for(int k=0; k<N; k++){
    float dot_product = 0.0;
    for(int j=0; j<N; j++){
      float val = A[j*N+k];
      dot_product += val*val;   //calculate dot product of the i-th column. Can be optimized. Data Dependencies + column-major accesses
      
    }
    
    R[k*N+k] = sqrt(dot_product);
    
    for(int j=0; j<N; j++){
      Q[j*N+k] = A[j*N+k]/R[k*N+k];
    }
    for(int j=k;j<N; j++){ //Note that the report is not able to caputure the value of k at compile time. Therefore there will be min and max.
      dot_product = 0.0   //BUT WE CAN STILL CALCULATE THE CORRECT LOOP LETENCY BY LINEAR INTERPOLATING.
      for(int i=0; i<N; i++){
        float val = Q[i*N+k];
        float val2 = A[i*N+j];
        dot_product += val*val2;
      }
      R[k*N+j] = dot_product;
      float val = R[k*N+j];
      for(int i=0; i<N; i++){
        A[i*N+j] -= Q[i*N+k]*val;
      }
      
    }

  }

}

