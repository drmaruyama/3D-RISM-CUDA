#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_xmu(double * & xmu, double pmv) {
    
  ofstream out_file;
  out_file.open((fname + extxmu).c_str());

  double ibeta = avogadoro * boltzmann * sv -> temper / kcal2J;

  double xmua = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    xmua += xmu[iv];
  }

  double gf = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    gf += xmu[sv -> natv + iv];
  }

  out_file << "Solvation Free Energy = " << ibeta * xmua
	   << " (kcal/mol)" << endl;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  " << iv << " -----> " << ibeta * xmu[iv] << endl;
  }
  out_file << endl;
  out_file << "SFE (GF) = " << ibeta * gf
	   << " (kcal/mol)" << endl;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  " << iv << " -----> " << ibeta * xmu[sv -> natv + iv] 
	     << endl;
  }
  out_file << endl;
  out_file << "PMV  = " << pmv << " (cc/mol)" << endl;

  out_file.close();
} 
