#include <iostream>
#include <thrust/device_vector.h>
#include "rism3d.h"

double RISM3D :: cal_exdft () {
  double excp_fd(double, double);
  __global__ void sum(double *, double2 *);
  __global__ void sum2(double *, double2 *, double2 *);
  __global__ void cal_dwork(double2 *, const double2 *, 
			    const double *, const double *, 
			    const int *, int, int);
  __global__ void rho0wfkhk(double2 *, const double2 *, const double *,
			    const int *, double);
  __global__ void rho0(double2 *, double);
  __global__ void rho0guvfmbex(double2 *, const double2 *, 
			       double, double, double);
  __global__ void wfka(double2 *, const double *, const int *);
  __global__ void cal_eda(double *, const double2 *, const double2 *, double);

  int ng = ce -> ngrid;

  double xmu1 = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    sum <<< g, b, b.x * sizeof(double) >>> (ds, dhuv + (iv * ng));
    thrust::device_ptr<double> ds_ptr(ds);
    double tmp = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
    xmu1 -= tmp * sv -> rhov[iv];
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), - 1);
  }

  double2 * dwork;
  cudaMalloc(&dwork, ng * sv -> natv * sizeof(double2));

  for (int iv = 0; iv < sv -> natv; ++iv) {
    cal_dwork <<< g, b >>> (dwork + (iv * ng), dhuv, 
			    sv -> dc + (iv * sv -> natv * nga), sv -> drho,
			    dindga, nga, sv -> natv);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), 1);
    fft -> execute(dwork + (iv * ng), 1);
  }

  double xmu2 = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    sum <<< g, b, b.x * sizeof(double) >>> (ds, dwork + (iv * ng));
    thrust::device_ptr<double> ds_ptr(ds);
    double tmp = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
    xmu2 += tmp * sv -> rhov[iv];
  }

  double xmu3 = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    sum2 <<< g, b, b.x * sizeof(double) >>> (ds, dwork + (iv * ng), 
					     dhuv + (iv * ng));
    thrust::device_ptr<double> ds_ptr(ds);
    double tmp = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
    xmu3 += tmp * sv -> rhov[iv];
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    rho0wfkhk <<< g, b >>> (dwork + (iv * ng), dhuv + (iv * ng), 
			    sv -> dw + (iv * nga), dindga, sv -> rhov[0]);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), 1);
    fft -> execute(dwork + (iv * ng), 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    rho0 <<< g, b >>> (dwork + (iv * ng), sv -> rhov[0]);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    double bexcp_fd = excp_fd(sv -> pfhs[iv], sv -> rhov[0]);
    double a = sv -> pfhs[iv] / sv -> rhov[0];
    rho0guvfmbex <<< g, b >>> (dwork + (iv * ng), dhuv + (iv * ng),
			       a, bexcp_fd, sv -> rhov[0]);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dwork + (iv * ng), - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    wfka <<< g, b >>> (dwork + (iv * ng), sv -> dw + (iv * nga), dindga);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dwork + (iv * ng), 1);
  }

  double eda = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double bexcp_fd = excp_fd(sv -> pfhs[iv], sv -> rhov[0]);
    double excp0_st = sv -> rhov[0] * bexcp_fd * sv -> wfk0[iv];
    cal_eda <<< g, b, b.x * sizeof(double) >>> (ds, dwork + (iv * ng), 
						dhuv + (iv * ng), excp0_st);
    thrust::device_ptr<double> ds_ptr(ds);
    double tmp = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
    eda += tmp * sv -> rhov[iv];
  }

  cudaFree(dwork);
  return (xmu1 + xmu2 + 0.5 * xmu3 - eda) * ce -> dv;
} 

__global__ void sum(double * ds, double2 * dhuv) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;

  sdata[threadIdx.x] = dhuv[ip].x;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    volatile double *smem = sdata;
    smem[threadIdx.x] += smem[threadIdx.x + 32];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 16];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 8];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 4];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 2];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 1];
  }
  if (threadIdx.x == 0) ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
}

__global__ void sum2(double * ds, double2 * dwork, double2 * dhuv) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;

  sdata[threadIdx.x] = dwork[ip].x * dhuv[ip].x;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    volatile double *smem = sdata;
    smem[threadIdx.x] += smem[threadIdx.x + 32];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 16];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 8];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 4];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 2];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 1];
  }
  if (threadIdx.x == 0) ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
}

__global__ void cal_dwork(double2 * dwork, 
			  const double2 * __restrict__ dhuv, 
			  const double * __restrict__ dc, 
			  const double * __restrict__ drho, 
			  const int * __restrict__ dindga, 
			  int nga, int natv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  unsigned int ip2 = dindga[ip];
  unsigned int ngr = blockDim.x * gridDim.x * gridDim.y;

  double wr = 0.0;
  double wi = 0.0;
  for (unsigned int iv = 0; iv < natv; ++iv) {
    unsigned int i = ip + iv * ngr;
    unsigned int i2 = ip2 + iv * nga;
    wr += dhuv[i].x * dc[i2] * drho[iv];
    wi += dhuv[i].y * dc[i2] * drho[iv];
  }
  dwork[ip].x = wr;
  dwork[ip].y = wi;
}

__global__ void rho0wfkhk(double2 * dwork, 
			  const double2 * __restrict__ dhuv,
			  const double * __restrict__ dw,
			  const int * __restrict__ dindga,
			  double rho0) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  unsigned int ip2 = dindga[ip];

  dwork[ip].x = rho0 * dw[ip2] * dhuv[ip].x;
  dwork[ip].y = rho0 * dw[ip2] * dhuv[ip].y; 
}

__global__ void rho0(double2 * dwork, double rho0) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;

  dwork[ip].x += rho0;
}

__global__ void rho0guvfmbex(double2 * dwork, 
			     const double2 * __restrict__ dhuv, 
			     double a, double bexcp_fd, double rho0) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  
  double dens = dwork[ip].x;
  double pfhs = a * dens;
  double fd = (4.0 - 2.0 * pfhs) / ((1.0 - pfhs) * (1.0 - pfhs) * (1.0 - pfhs))
    * (pfhs / dens);

  dwork[ip].x = rho0 * ((dhuv[ip].x + 1.0) * fd - bexcp_fd);
  dwork[ip].y = 0.0;
}

__global__ void wfka(double2 * dwork, 
		     const double * __restrict__ dw,
		     const int * __restrict__ dindga) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  unsigned int ip2 = dindga[ip];

  dwork[ip].x *= dw[ip2];
  dwork[ip].y *= dw[ip2];
}

__global__ void cal_eda(double * ds, const double2 * __restrict__ dwork, 
			const double2 * __restrict__ dhuv, double excp0_st) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;

  sdata[threadIdx.x] = (dwork[ip].x + excp0_st) * (dhuv[ip].x + 1.0) 
    - excp0_st;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    volatile double *smem = sdata;
    smem[threadIdx.x] += smem[threadIdx.x + 32];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 16];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 8];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 4];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 2];
    __syncwarp();
    smem[threadIdx.x] += smem[threadIdx.x + 1];
  }
  if (threadIdx.x == 0) ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
}
