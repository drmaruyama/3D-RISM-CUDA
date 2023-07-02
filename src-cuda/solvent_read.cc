#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "solvent.h"
#include "physical.h"

void Solvent :: read(string fsolvent) {
  void alloc3D (vector <vector <double * > > &, int, int, int);

  ifstream in_file;
  in_file.open(fsolvent.c_str());

  if (!in_file) {
    std::cerr << "Solvent file could not be opened!" << std::endl;
    exit (1);
  }

  double * sig;
  double * eps;
  double * q;
  double * den;
  int * sol;
  string dummy;  
  double dr;
  int num;
  auto i = 0;
  std::string line;
  while (std::getline(in_file, line)) {
    std::istringstream iss(line);
    if (i == 1) {
      iss  >> dummy >> num >> natv >> ntab >> dr >> dummy;
      alloc3D (xvv, natv, natv, ntab);
      sig = new double[num];
      eps = new double[num];
      q = new double[num];
      den = new double[num];
      sol = new int[num];
    }
    if (i == 2) {
      iss  >> dummy >> dummy >> temper >> xt;
    }
    if (i > 3 && i <= num + 3) {
      int i2 = i - 4;
      iss >> dummy >> dummy >> dummy >> sol[i2]
          >> sig[i2] >> eps[i2] >> q[i2]
          >> dummy >> dummy >> dummy >> den[i2];
    }
    if (i > num + 3 && i < num + 4 + ntab) {
      int i2 = i - num - 4;
      for (int iv2 = 0; iv2 < natv; ++iv2) {
	for (int iv1 = 0; iv1 < natv; ++iv1) {
	  in_file >> xvv[iv2][iv1][i2];
	}
      }
    }
    ++i;
  }

  sigv = new double[natv];
  epsv = new double[natv];
  qv = new double[natv];
  rhov = new double[natv];

  for (int i = 0; i < natv; ++i) {
    rhov[i] = 0.0;
  }
  for (int i = 0; i < num; ++i) {
    int i2 = abs(sol[i]) - 1;
    sigv[i2] = sig[i];
    epsv[i2] = eps[i];
    qv[i2] = q[i];
    rhov[i2] += den[i];
  }
  for (int i = 0; i < natv; ++i) {
    rhov[i] *= (avogadoro * 1.0e-27);
  }

  double dk = M_PI / (dr * ntab);
  ttab = new double[ntab];
  for (int i = 0; i < ntab; ++i) {
    ttab[i] = dk * (i + 1);
  }

  in_file.close();

  delete[] sig;
  delete[] eps;
  delete[] q;
  delete[] den;
  delete[] sol;
}
