#include "solute.h"

void Solute :: init(int n) {
  num = n;
  q = new double[num];
  sig = new double[num];
  eps = new double[num];
  r = new double[num * 4];
}


void Solute :: setup_cuda() {
  cudaMalloc(&dq, num * sizeof(double));
  cudaMalloc(&dr, num * sizeof(double4));
  cudaMemcpyAsync(dq, q, num * sizeof(double), cudaMemcpyDefault);
  cudaMemcpyAsync(dr, r, num * sizeof(double4), cudaMemcpyDefault);
}
