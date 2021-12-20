#include "CavityFlow.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"
#include "hlslib/xilinx/DataPack.h"

using hlslib::Stream;

using Vec_t = hlslib::DataPack<float, W>;


void CopyArray5Times(float * inp, float* out0, float* out1, float* out2, float* out3, float* out4, int size)
{
  for(int i=0;i<size;i++)
  {
    #pragma HLS PIPELINE II=1
    out0[i] = inp[i];
    out1[i] = inp[i];
    out2[i] = inp[i];
    out3[i] = inp[i];
    out4[i] = inp[i];
  }
}

void CopyStream5Times(Stream<float>& inp, Stream<float>& out0, Stream<float>& out1, Stream<float>& out2, Stream<float>& out3, Stream<float>& out4, int size)
{
  for(int i=0;i<size;i++)
  {
    #pragma HLS PIPELINE II=1
    float val = inp.Pop();
    out0.Push(val);
    out1.Push(val);
    out2.Push(val);
    out3.Push(val);
    out4.Push(val);
  }
}

void SplitStream(Stream<float>& inp, Stream<float>& out1, Stream<float>& out2, int size)
{
  for(int i=0;i<size;i++)
  {
    #pragma HLS PIPELINE II=1
    float val = inp.Pop();
    out1.Push(val);
    out2.Push(val);
  }
}
// def build_up_b_fpga(b, rho, dt, u, v, dx, dy):

//     b[1:-1,
//       1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
//                                 (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
//                       ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 *
//                       ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
//                        (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
//                       ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

void build_up_b_fpga(float b[], 
  float rho, float dt, 
  float u0[],  float u1[],  float u2[],  float u3[],  float u4[], 
  float v0[],  float v1[],  float v2[],  float v3[],  float v4[], 
  float dx, float dy)
{
  
  
  for(int i=1;i<NY-1;i++)
  {
	
  	for(int j=1;j<NX-1;j++)
    { // Can potentially be squzeezed into one single loop (may not provide any benefits)
      #pragma HLS PIPELINE II=1
      // Note the difference in indexing and denominators of the two\
        sets of variables...

			float u_var = u0[(i) * NX + (j)];
			float u_left_var = u1[(i) * NX + (j-1)];
			float u_right_var = u2[(i) * NX + (j+1)];
			float u_below_var = u3[(i+1) * NX + j];
			float u_above_var = u4[(i-1) * NX + (j)];

			float v_var = v0[(i) * NX + (j)];
			float v_left_var = v1[(i) * NX + (j-1)];
			float v_right_var = v2[(i) * NX + (j+1)];
			float v_below_var = v3[(i+1) * NX + j];
			float v_above_var = v4[(i-1) * NX + (j)];
      // Don't use multiple memory accesses in same loop. (for each interface i.e)

      // (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)
      float var_1 = (u_right_var - u_left_var) / (2 * dx); 
      // (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)
      float var_2 = (v_below_var - v_above_var) / (2 * dy); // Parallelise/vectorize
      // (u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy)
      float var_3 = (u_below_var - u_above_var) / (2 * dy); 
      // (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)
      float var_4 = (v_right_var - v_left_var) / (2 * dx); // Parallelise/vectorize

      // (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
      // (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
      // ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 *
      // ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
      // (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
      // ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)
      float varB = (1/dt * (var_1 + var_2)) - var_1 * var_1 - 2*var_3*var_4 - var_2 * var_2;

      // b[1:-1,1:-1] = rho * varB
      b[i * NX + j] = rho * varB;
    }

  }
}

// def pressure_poisson_fpga(nit, p, dx, dy, b):
//     pn = np.empty_like(p)
//     pn = p.copy()

//     for q in range(nit):
//         pn = p.copy()
//         p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
//                           (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
//                          (2 * (dx**2 + dy**2)) - dx**2 * dy**2 /
//                          (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

//         p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
//         p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
//         p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
//         p[-1, :] = 0  # p = 0 at y = 2


void pressure_poisson_fpga(
  float p0[],  float p1[],  float p2[],  float p3[],  float p4[], 
  float dx, float dy, float b[])
{
  // pn = np.empty_like(p)
  float pn0[NX * NY];
  float pn1[NX * NY];
  float pn2[NX * NY];
  float pn3[NX * NY];
  float pn4[NX * NY];

  // dx**2
  float dx_2 = dx * dx;
  // dy**2
  float dy_2 = dy * dy;
  
  // for q in range(nit):
  for (int q = 0;q < NIT;q++)
  { 
    //pn = p.copy()
    CopyArray5Times(p0, pn0, pn1, pn2, pn3, pn4, NX * NY);

    // iterate [1:-1, :]
    for(int i=1;i<NY-1;i++)
    {

      // iterate [1:-1, 1:-1]
      for(int j=1;j<NX-1;j++)
      {
        // Due to the decision logic the worst case II will be 2, 
        // but this shouldn't be much of a problem for the general calculation
        #pragma HLS PIPELINE II=2
        
        // (pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2
        float varA = (pn0[i * NX + j+1] + pn2[i * NX + j-1]) * dy_2;
        // (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2)
        float varB = (pn1[(i+1) * NX + j] + pn3[(i-1) * NX + j]) * dx_2;

        // ((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
        // (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
        // (2 * (dx**2 + dy**2))
        float varC = (varA + varB) / (2 * (dx_2 + dy_2));

        // dx**2 * dy**2 /
        // (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        float varD = dx_2 * dy_2/(2 * (dx_2 + dy_2)) * b[i*NX + j];

        // (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
        // (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
        //(2 * (dx**2 + dy**2)) - dx**2 * dy**2 /
        //(2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        float p_val = varC - varD;
        p0[i*NX + j] = p_val;
				p1[i*NX + j] = p_val;
				p2[i*NX + j] = p_val;
				p3[i*NX + j] = p_val;
				p4[i*NX + j] = p_val;
        
        // p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        if (i == 1)
        {
          p0[j] = p_val;
          p1[j] = p_val;
          p2[j] = p_val;
          p3[j] = p_val;
          p4[j] = p_val;
        }
        
        // p[-1, :] = 0  # p = 0 at y = 2
        if (i == (NY - 2))
        {
          p0[(i+1) * NX + j] = 0;
          p1[(i+1) * NX + j] = 0;
          p2[(i+1) * NX + j] = 0;
          p3[(i+1) * NX + j] = 0;
          p4[(i+1) * NX + j] = 0;
        }

        // p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
        if (j == (1))
        {
          p0[(i) * NX + 0] = p_val;
          p1[(i) * NX + 0] = p_val;
          p2[(i) * NX + 0] = p_val;
          p3[(i) * NX + 0] = p_val;
          p4[(i) * NX + 0] = p_val;
        }

        // p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
        if (j == (NX - 2))
        {
          p0[(i) * NX + (NX - 1)] = p_val;
          p1[(i) * NX + (NX - 1)] = p_val;
          p2[(i) * NX + (NX - 1)] = p_val;
          p3[(i) * NX + (NX - 1)] = p_val;
          p4[(i) * NX + (NX - 1)] = p_val;
        }
      }
    }
  }
}

// def cavity_flow(nx, ny, nt, nit, u, v, dt, dx, dy, p, rho, nu):
//     un = np.empty_like(u)
//     vn = np.empty_like(v)
//     b = np.zeros((ny, nx))

//     for n in range(nt):
//         un = u.copy()
//         vn = v.copy()

//         build_up_b_fpga(b, rho, dt, u, v, dx, dy)
//         pressure_poisson_fpga(nit, p, dx, dy, b)

//         u[1:-1,
//           1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
//                    (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
//                    vn[1:-1, 1:-1] * dt / dy *
//                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) - dt / (2 * rho * dx) *
//                    (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu *
//                    (dt / dx**2 *
//                     (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
//                     dt / dy**2 *
//                     (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

//         v[1:-1,
//           1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
//                    (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
//                    vn[1:-1, 1:-1] * dt / dy *
//                    (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) *
//                    (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
//                    (dt / dx**2 *
//                     (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
//                     dt / dy**2 *
//                     (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

//         u[0, :] = 0
//         u[:, 0] = 0
//         u[:, -1] = 0
//         u[-1, :] = 1  # set velocity on cavity lid equal to 1
//         v[0, :] = 0
//         v[-1, :] = 0
//         v[:, 0] = 0
//         v[:, -1] = 0
//         return u, v


void ProcessingElement(
  Stream<float>& u0_inp, Stream<float>& u1_inp, Stream<float>& u2_inp, Stream<float>& u3_inp, Stream<float>& u4_inp, 
  Stream<float>& u0_out, Stream<float>& u1_out, Stream<float>& u2_out, Stream<float>& u3_out, Stream<float>& u4_out, 
  Stream<float>& v0_inp, Stream<float>& v1_inp, Stream<float>& v2_inp, Stream<float>& v3_inp, Stream<float>& v4_inp, 
  Stream<float>& v0_out, Stream<float>& v1_out, Stream<float>& v2_out, Stream<float>& v3_out, Stream<float>& v4_out, 
  float p0[], float p1[], float p2[], float p3[], float p4[],
  float b[], float dt, float dx, float dy
)
{   

    //un = np.empty_like(u)
    float un0[NY * NX];
    float un1[NY * NX];
    float un2[NY * NX];
    float un3[NY * NX];
    float un4[NY * NX];
    
    //vn = np.empty_like(v)
    float vn0[NY * NX];
    float vn1[NY * NX];
    float vn2[NY * NX];
    float vn3[NY * NX];
    float vn4[NY * NX];

    float u0[NY*NX];
    float u1[NY*NX];
    float u2[NY*NX];
    float u3[NY*NX];
    float u4[NY*NX];

    CopyArray5Times(u, u0, u1, u2, u3, u4, NX * NY);

    float v0[NY*NX];
    float v1[NY*NX];
    float v2[NY*NX];
    float v3[NY*NX];
    float v4[NY*NX];

  CopyArray5Times(v, v0, v1, v2, v3, v4, NX * NY);

    
    CopyArray5Times(u0, un0, un1, un2, un3, un4, NX * NY);
    CopyArray5Times(v0, vn0, vn1, vn2, vn3, vn4, NX * NY);

    // build_up_b_fpga(b, rho, dt, u, v, dx, dy)
    build_up_b_fpga(b,RHO, dt, u0, u1, u2, u3, u4, v0, v1, v2, v3, v4, dx, dy);
    // pressure_poisson_fpga(nit, p, dx, dy, b)
    pressure_poisson_fpga(p0, p1, p2, p3, p4, dx, dy, b);

    for(int i=0;i<NY;i++)
    {
		  
      for(int j=0;j<NX;j++)
      {
        #pragma HLS PIPELINE II=1
        float u_val, v_val;

        if ((i == 0) || (j == 0) || (i == (NY - 1)) || (j == (NX -1)))
        {
          u_val = 0;
          v_val = 0;

          if (i == (NY - 1))
          {
            u_val = 1;
          }
        }
        else{

          float p_var = p0[(i)*NX + j];
          float p_left_var = p1[(i)*NX + j-1];
          float p_right_var = p2[(i)*NX + (j+1)];
          float p_below_var = p3[(i+1) * NX + j];
          float p_above_var = p4[(i-1) * NX + j];

          float un_var = un0[(i)*NX + j];
          float un_left_var = un1[(i)*NX + j-1];
          float un_right_var = un2[(i)*NX + (j+1)];
          float un_below_var = un3[(i+1) * NX + j];
          float un_above_var = un4[(i-1) * NX + j];
          
          float vn_var = vn0[(i)*NX + j];
          float vn_left_var = vn1[(i)*NX + j-1];
          float vn_right_var = vn2[(i)*NX + (j+1)];
          float vn_below_var = vn3[(i+1) * NX + j];
          float vn_above_var = vn4[(i-1) * NX + j];



          // u[1:-1,1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
          // (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
          // vn[1:-1, 1:-1] * dt / dy *
          // (un[1:-1, 1:-1] - un[0:-2, 1:-1]) - dt / (2 * rho * dx) *
          // (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu *
          // (dt / dx**2 *
          // (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
          // dt / dy**2 *
          // (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))
          u_val = (un_var - un_var * dt/dx *\
              (un_var - un_left_var) -\
              vn_var * dt/dy *\
              (un_var - un_above_var) - dt/(2 * RHO * dx) *\
              (p_right_var - p_left_var) + NU *\
              (dt/(dx * dx) * \
              (un_right_var - 2* un_var + un_left_var) + \
              dt / (dy*dy) *\
              (un_below_var -2 * un_var + un_above_var)));
          
          
          
          
          //  v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
          // (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
          // vn[1:-1, 1:-1] * dt / dy *
          // (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) *
          // (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
          // (dt / dx**2 *
          // (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
          // dt / dy**2 *
          // (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
          v_val = (vn_var - un_var*dt/dx *\
                          (vn_var - vn_left_var) - 
                          vn_var * dt/dy * \
                          (vn_var - vn_above_var) - dt/(2*RHO*dy) *\
                          (p_below_var - p_above_var) + NU*\
                          (dt / (dx * dx) * \
                          (vn_right_var - 2*vn_var + vn_left_var) + \
                          dt / (dy*dy) * \
                          (vn_below_var - 2*vn_var + vn_above_var))
                        );
        }
        u0[i*NX + j] = u_val;
        u1[i*NX + j] = u_val;
        u2[i*NX + j] = u_val;
        u3[i*NX + j] = u_val;
        u4[i*NX + j] = u_val;

        v0[i*NX + j] = v_val;
        v1[i*NX + j] = v_val;
        v2[i*NX + j] = v_val;
        v3[i*NX + j] = v_val;
        v4[i*NX + j] = v_val;

      }
    }
}

void CavityFlow(float *u, float *v, float* p) {

  #pragma HLS INTERFACE m_axi port=u offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=v offset=slave bundle=gmem1
  #pragma HLS INTERFACE m_axi port=p offset=slave bundle=gmem2
  #pragma HLS INTERFACE s_axilite port=u bundle=control
  #pragma HLS INTERFACE s_axilite port=v bundle=control
  #pragma HLS INTERFACE s_axilite port=p bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control

  float u0[NY*NX];
  float u1[NY*NX];
  float u2[NY*NX];
  float u3[NY*NX];
  float u4[NY*NX];

  CopyArray5Times(u, u0, u1, u2, u3, u4, NX * NY);

  float v0[NY*NX];
  float v1[NY*NX];
  float v2[NY*NX];
  float v3[NY*NX];
  float v4[NY*NX];

  CopyArray5Times(v, v0, v1, v2, v3, v4, NX * NY);

  
  // np.zeros((ny, nx))
  float b[NY*NX] = {0};

  float p0[NX*NY];
  float p1[NX*NY];
  float p2[NX*NY];
  float p3[NX*NY];
  float p4[NX*NY];

  CopyArray5Times(p, p0, p1, p2, p3, p4, NX * NY);

  float dx = 2.0 / (NX - 1);
  float dy = 2.0 / (NY - 1);
  float dt = 0.1 / ((NX - 1) * (NY - 1));

  // for n in range(nt)
  for(int n = 0;n < NT;n++)
  {
    ProcessingElement(u0, u1, u2, u3, u4,
                      v0, v1, v2, v3, v4,
                      p0, p1, p2, p3, p4,
                      b, dt, dx, dy);    
  }

  for (int i=0;i<NX*NY;i++)
  {
    #pragma HLS PIPELINE II=1
    u[i] = u0[i];
    v[i] = v0[i];
  }

}
