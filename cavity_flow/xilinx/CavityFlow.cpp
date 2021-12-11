#include "CavityFlow.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"
#include "hlslib/xilinx/DataPack.h"

using hlslib::Stream;

using Vec_t = hlslib::DataPack<float, W>;

// def build_up_b_fpga(b, rho, dt, u, v, dx, dy):

//     b[1:-1,
//       1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
//                                 (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
//                       ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 *
//                       ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
//                        (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
//                       ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

void build_up_b_fpga(float b[], float rho, float dt, float u[], float v[], float dx, float dy)
{
  float u_above[NX];
  float u_centre_h[NX];
  
	float u_start_col[NY];
	float u_end_col[NY];

  for(int j = 0;j<NX;j++)
  {
    #pragma HLS PIPELINE II=1
    u_above[j] = u[ 0 * NX + j];
  }

  for(int j=0;j<NX;j++)
  {
    #pragma HLS PIPELINE II=1
    u_centre_h[j] = u[1 * NX + j];
  }

	for(int i=0;i<NY;i++)
	{
		#pragma HLS PIPELINE II=1
		u_start_col[i] = u[i*NX + 0];
	}
  
	for(int i=0;i<NY;i++)
	{
		#pragma HLS PIPELINE II=1
		u_end_col[i] = u[i*NX + NX-1];
	}

  float v_above[NX];
  float v_centre_h[NX];
  
	float v_start_col[NY];
	float v_end_col[NY];

  for(int j = 0;j<NX;j++)
  {
    #pragma HLS PIPELINE II=1
    v_above[j] = v[ 0 * NX + j];
  }

  for(int j=0;j<NX;j++)
  {
    #pragma HLS PIPELINE II=1
    v_centre_h[j] = v[1 * NX + j];
  }

	for(int i=0;i<NY;i++)
	{
		#pragma HLS PIPELINE II=1
		v_start_col[i] = v[i*NX + 0];
	}
  
	for(int i=0;i<NY;i++)
	{
		#pragma HLS PIPELINE II=1
		v_end_col[i] = v[i*NX + NX-1];
	}
	

  for(int i=1;i<NY-1;i++)
  {
		float u_below_var_for_next_iter = u_start_col[i+1];
		float v_below_var_for_next_iter = v_start_col[i+1];

		for(int j=1;j<NX-1;j++)
    { // Can potentially be squzeezed into one single loop (may not provide any benefits)
      #pragma HLS PIPELINE II=1
      // Note the difference in indexing and denominators of the two\
        sets of variables...

			float u_var = u_centre_h[j];
			float u_left_var = u_centre_h[j-1];
			float u_right_var = u_centre_h[(j+1)];
			float u_below_var = u[(i+1) * NX + j];
			float u_above_var = u_above[j];

			float v_var = v_centre_h[j];
			float v_left_var = v_centre_h[j-1];
			float v_right_var = v_centre_h[(j+1)];
			float v_below_var = v[(i+1) * NX + j];
			float v_above_var = v_above[j];
      // Don't use multiple memory accesses in same loop. (for each interface i.e)

      // (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)
      float var_1 = (u_right_var - u_left_var) / (2 * dx); 
      // (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)
      float var_2 = (v_below_var - v_above_var) / (2 * dy); // Parallelise/vectorize
      // (u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy)
      float var_3 = (u_below_var - u_above_var) / (2 * dy); 
      // (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)
      float var_4 = (v_right_var - v_left_var) / (2 * dx); // Parallelise/vectorize
      
			u_above[j] = u_var;
			u_centre_h[j-1] = u_below_var_for_next_iter;
			u_below_var_for_next_iter = u_below_var;

			v_above[j] = v_var;
			v_centre_h[j-1] = v_below_var_for_next_iter;
			v_below_var_for_next_iter = v_below_var;

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

		u_centre_h[NX-2] = u_below_var_for_next_iter;
		u_centre_h[NX-1] = u_end_col[i+1];

		v_centre_h[NX-2] = v_below_var_for_next_iter;
		v_centre_h[NX-1] = v_end_col[i+1];
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


void pressure_poisson_fpga(float p[], float dx, float dy, float b[])
{
  // pn = np.empty_like(p)
  float pn[NX * NY];
  float p_tmp[NX * NY];
  
	for(int i=0;i<NY*NX;i++)
	{
		#pragma HLS PIPELINE II=1
		p_tmp[i] = p[i];
	}
  // for q in range(nit):
  for (int q = 0;q < NIT;q++)
  { 
    //pn = p.copy()
    for(int i=0;i<NY*NX;i++)
    {
      #pragma HLS PIPELINE II=1
      pn[i] = p[i];
    }


    // iterate [1:-1, :]
    for(int i=1;i<NY-1;i++)
    {

      // iterate [1:-1, 1:-1]
      for(int j=1;j<NX-1;j++)
      {
        #pragma HLS PIPELINE II=1
        // dx**2
        float dx_2 = dx * dx;
        // dy**2
        float dy_2 = dy * dy;
        
        // (pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2
        float varA = (pn[i * NX + j+1] + pn[i * NX + j-1]) * dy_2;
        // (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2)
        float varB = (pn[(i+1) * NX + j] + pn[(i-1) * NX + j]) * dx_2;

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
        p_tmp[i * NX + j] = varC - varD;
				p[i*NX + j] = varC - varD;
      }
    }

    // p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      p[i * NX + NX-1] = p_tmp[i*NX + NX-2];
      #pragma HLS DEPENDENCE variable=p false
    } 
    // p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      p[0 * NX + i] = p_tmp[1*NX + i];
      #pragma HLS DEPENDENCE variable=p false
    } 
    // p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      p[i * NX + 0] = p_tmp[i*NX + 1];
      #pragma HLS DEPENDENCE variable=p false
    } 
    // p[-1, :] = 0  # p = 0 at y = 2
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      p[(NY-1) * NX + i] = 0;
      #pragma HLS DEPENDENCE variable=p false
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


void CavityFlowFunc(float u[], float v[], float p[]) {
  //un = np.empty_like(u)
  float un[NY * NX];
  //vn = np.empty_like(v)
  float vn[NY * NX];
  // np.zeros((ny, nx))
  float b[NY*NX] = {0};

  float dx = 2.0 / (NX - 1);
  float dy = 2.0 / (NY - 1);
  float dt = 0.1 / ((NX - 1) * (NY - 1));

  // for n in range(nt)
  for(int n = 0;n < NT;n++)
  {
    // un = u.copy()
    // vn = v.copy()
    for(int i = 0;i<NY*NX;i++)
    {
      #pragma HLS PIPELINE II=1
      un[i] = u[i];
      vn[i] = v[i];
    }

    // build_up_b_fpga(b, rho, dt, u, v, dx, dy)
    build_up_b_fpga(b, RHO, dt, u, v, dx, dy);
    // pressure_poisson_fpga(nit, p, dx, dy, b)
    pressure_poisson_fpga(p, dx, dy, b);

    for(int i=1;i<NY-1;i++)
    {
      for(int j=1;j<NX-1;j++)
      {
        #pragma HLS LOOP_FLATTEN
        #pragma HLS PIPELINE II=1
        // u[1:-1,1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
        // (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
        // vn[1:-1, 1:-1] * dt / dy *
        // (un[1:-1, 1:-1] - un[0:-2, 1:-1]) - dt / (2 * rho * dx) *
        // (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu *
        // (dt / dx**2 *
        // (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
        // dt / dy**2 *
        // (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))
        u[i*NX + j] = (un[i * NX + j] - un[i * NX + j] * dt/dx *\
            (un[i * NX + j] - un[i * NX + j-1]) -\
            vn[i * NX + j] * dt/dy *\
            (un[i * NX + j] - un[(i-1) * NX + j]) - dt/(2 * RHO * dx) *\
            (p[i * NX + j+1] - p[i * NX + j-1]) + NU *\
            (dt/(dx * dx) * \
            (un[i * NX + j+1] - 2* un[i * NX + j] + un[i * NX + j-1]) + \
            dt / (dy*dy) *\
            (un[(i + 1)* NX + j] -2 * un[i * NX + j] + un[(i-1) * NX + j])));
        
        
        
        //  v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
        // (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
        // vn[1:-1, 1:-1] * dt / dy *
        // (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) *
        // (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
        // (dt / dx**2 *
        // (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
        // dt / dy**2 *
        // (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        v[i*NX + j] = (vn[i*NX + j] - un[i*NX + j]*dt/dx *\
                        (vn[i*NX + j] - vn[i*NX + j-1]) - 
                        vn[i*NX + j] * dt/dy * \
                        (vn[i*NX + j] - vn[(i-1)*NX + j]) - dt/(2*RHO*dy) *\
                        (p[(i+1)*NX + j] - p[(i-1)*NX + j]) + NU*\
                        (dt / (dx * dx) * \
                        (vn[i*NX + j+1] - 2*vn[i*NX + j] + vn[i*NX + j-1]) + \
                        dt / (dy*dy) * \
                        (vn[(i+1)*NX + j] - 2*vn[i*NX + j] + vn[(i-1)*NX + j]))
                      );

      }
    }

    // u[0, :] = 0
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      u[0 * NX + i] = 0;
    }
    // u[:, 0] = 0
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      u[i * NX + 0] = 0;
    }
    // u[:, -1] = 0
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      u[i * NX + NX-1] = 0;
    }
    // u[-1, :] = 1
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      u[(NY-1) * NX + i] = 1;
    }

    // v[0, :] = 0
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      v[0 * NX + i] = 0;
    }
    // v[-1, :] = 0
    for(int i=0;i<NX;i++)
    {
      #pragma HLS PIPELINE II=1
      v[(NY-1) * NX + i] = 0;
    }
    // v[:, 0] = 0
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      v[i * NX + 0] = 0;
    }
    // v[:, -1] = 0
    for(int i=0;i<NY;i++)
    {
      #pragma HLS PIPELINE II=1
      v[i * NX + NX-1] = 0;
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

  CavityFlowFunc(u, v, p);
}
